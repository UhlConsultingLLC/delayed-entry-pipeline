"""The GENIE BPC adapter: a BPC release into the standard tables, with its capability profile.

    from delayed_entry.adapters import genie_bpc
    result = genie_bpc.read(clinical_dir, genomic_dir, settings)
    result.tables, result.profile, result.input_files

Built from the CRC BPC project's cohort code (copied, not imported). The cohort's own files live in
cohorts/<cohort>/: the raw file list with column types, the release list for the check after download,
and the cohort's extra columns. What is specific to a study (which labels mean this cancer, what to do
with a withheld survival time) comes from the settings' source.genie_bpc section.
"""
from __future__ import annotations

import dataclasses
import re
from pathlib import Path

import pandas as pd
import yaml

from delayed_entry.adapters import AdapterResult
from delayed_entry.adapters.genie_bpc import raw
from delayed_entry.privacy import show
from delayed_entry.profile import CapabilityProfile, profile_from_dict
from delayed_entry.settings import Settings
from delayed_entry.tables import StandardTables, conform, table_fields

HERE = Path(__file__).parent
STANDARD_SQL = HERE / "standard.sql"
LABELS = HERE / "labels.yml"
PROFILE = HERE / "profile.yml"
RELEASES = HERE / "synapse_releases.yml"
COVERAGE_FILES = ("genomic_information*", "*.bed")    # a panel coverage file or panel sizes
FIRST_ORDER_SQL = {"earliest_day": "t.days_to_entry_event ASC NULLS LAST",
                   "largest_panel": "n_genes DESC NULLS LAST",
                   "lowest_test_order": "t.test_order ASC NULLS LAST"}


def profile() -> CapabilityProfile:
    return profile_from_dict(yaml.safe_load(PROFILE.read_text(encoding="utf-8")), adapter="genie_bpc")


def with_release(prof: CapabilityProfile, cohort: str) -> CapabilityProfile:
    """The profile with its attribution's dataset version filled in: the release the settings name, and its
    Synapse ID from synapse_releases.yml. A version the profile already gives is kept."""
    if prof.attribution is None or prof.attribution.get("version"):
        return prof
    release = (yaml.safe_load(RELEASES.read_text(encoding="utf-8")) or {}).get(cohort, {}).get("release")
    version = f"GENIE BPC {cohort}" + (f" (Synapse {release})" if release else "")
    return dataclasses.replace(prof, attribution={**prof.attribution, "version": version})


def read_steps(sql_file: Path = STANDARD_SQL) -> list[tuple[str, str]]:
    """The SQL file as [(step name, statement)], each step starting with a line '-- step: name'."""
    parts = re.split(r"^-- step: *(\S+) *$", sql_file.read_text(encoding="utf-8"), flags=re.M)
    return [(name, sql.strip().rstrip(";")) for name, sql in zip(parts[1::2], parts[2::2])]


def labels_table() -> list[tuple[str, str, str]]:
    spec = yaml.safe_load(LABELS.read_text(encoding="utf-8"))
    return [(test, label, reading) for test, mapping in spec.items() for label, reading in mapping.items()]


def coverage_files(genomic_dir: Path) -> list[str]:
    return sorted(p.name for pattern in COVERAGE_FILES for p in Path(genomic_dir).glob(pattern))


def _frame(con, table: str, standard: str) -> pd.DataFrame:
    df = con.execute(f"SELECT * FROM {table}").df()
    keep = [f.name for f in table_fields(standard) if f.name in df.columns]
    extras = [c for c in df.columns if c not in keep] if standard == "patients" else []
    return df[keep + extras]


EMPTY_PATHOLOGY_RESULT = ("CREATE OR REPLACE TABLE pathology_result (path_row BIGINT, record_id VARCHAR, test VARCHAR, "
                          "slot INTEGER, label VARCHAR, dob_result_days DOUBLE)")


def cohort_spec(cohort: str) -> dict:
    """The cohort's own settings in cohorts/<cohort>/cohort.yml: extra patient columns, the results read from
    pathology reports, and files used besides those in raw_files.yml."""
    spec = yaml.safe_load((raw.cohort_dir(cohort) / "cohort.yml").read_text(encoding="utf-8")) or {}
    unknown = set(spec.get("pathology_results") or []) - {"msi_mmr"}
    if unknown:
        raise ValueError(f"cohort {cohort}: reading {sorted(unknown)} from pathology reports isn't written")
    return spec


def build(con, settings: Settings, cohort: str, prof: CapabilityProfile) -> None:
    """Run standard.sql on a connection holding the raw tables. A cohort whose pathology reports carry no
    MSI/MMR results (cohort.yml) gets an empty pathology_result table instead of that step."""
    order = settings.read("entry.first_test_order")
    first_order = ", ".join([FIRST_ORDER_SQL[k] for k in order] + ["t.test_id ASC"])
    spec = cohort_spec(cohort)
    extras = "".join(f",\n    {col} AS {name}" for name, col in (spec.get("extra_patient_fields") or {}).items())
    msi_mmr = "msi_mmr" in (spec.get("pathology_results") or [])
    regimen_types = raw.raw_files(cohort).get("raw_regimen", ("", "", {}))[2]
    regimen = "r.regimen_drugs" if "regimen_drugs" in regimen_types else "CAST(NULL AS VARCHAR)"
    params = {
        "censor": settings.read("source.genie_bpc.withheld_survival") == "censor_at_age_limit",
        "age_limit": float(prof.age_top_code),
        "days_per_year": float(settings.read("time.days_per_year")),
        "labels": list(settings.read("source.genie_bpc.cancer_type_labels")),
    }
    con.execute(raw.EMPTY_FUSION)       # hand-made raw tables (tests) may have no fusion table
    con.execute("CREATE OR REPLACE TABLE bpc_label (test VARCHAR, label VARCHAR, reading VARCHAR)")
    con.executemany("INSERT INTO bpc_label VALUES (?, ?, ?)", labels_table())
    for name, sql in read_steps():
        if name == "pathology_result" and not msi_mmr:
            con.execute(EMPTY_PATHOLOGY_RESULT)
            continue
        sql = sql.replace("{first_order}", first_order).replace("{extras}", extras).replace("{regimen}", regimen)
        used = {k: v for k, v in params.items() if f"${k}" in sql}
        con.execute(sql, used) if used else con.execute(sql)
        if name == "pathology_result":
            unknown = con.execute("""SELECT DISTINCT r.test, r.label FROM pathology_result AS r
                                     LEFT JOIN bpc_label AS l ON l.test = r.test AND l.label = r.label
                                     WHERE l.reading IS NULL ORDER BY ALL""").fetchall()
            if unknown:
                raise ValueError(f"result labels the adapter doesn't know; add them to labels.yml: {unknown}")


def tables_from_connection(con, settings: Settings, cohort: str | None = None,
                           prof: CapabilityProfile | None = None) -> tuple[StandardTables, dict[str, int]]:
    """The standard tables from a connection that holds the raw tables (raw.load, or hand-made ones in a
    test), and how many values in each field failed to convert to its type."""
    cohort = cohort or settings.read("source.genie_bpc.cohort")
    prof = prof or profile()
    build(con, settings, cohort, prof)
    tests = pd.concat([_frame(con, "std_tests_ngs", "tests"), _frame(con, "std_tests_pathology", "tests")],
                      ignore_index=True)
    panels = con.execute("SELECT panel AS panel_id, gene, CAST(NULL AS INTEGER) AS codon_from, "
                         "CAST(NULL AS INTEGER) AS codon_to FROM raw_panel_gene").df()
    tables = StandardTables(
        patients=_frame(con, "std_patients", "patients"),
        tests=tests,
        results=_frame(con, "std_results", "results"),
        panels=panels,
        treatments=_frame(con, "std_treatments", "treatments"),
        panel_sizes=None,
    )
    return conform(tables)


def read(clinical_dir: str | Path, genomic_dir: str | Path | None, settings: Settings) -> AdapterResult:
    """Read a GENIE BPC release into the standard tables."""
    cohort = settings.read("source.genie_bpc.cohort")
    genomic_dir = Path(genomic_dir or clinical_dir)
    found = coverage_files(genomic_dir)
    if found:
        raise NotImplementedError(f"the genomic folder has a panel coverage or size file ({', '.join(found)}); "
                                  "reading it isn't written yet, and it would change the capability profile")
    prof = with_release(profile(), cohort)
    con, files = raw.load(cohort, clinical_dir, genomic_dir)
    tables, failures = tables_from_connection(con, settings, cohort, prof)
    min_count = max(settings.read("privacy.min_count"), prof.privacy_min_count or 0)
    no_cancer = raw.scalar(con, "SELECT count(*) FROM raw_patient WHERE record_id NOT IN (SELECT record_id FROM index_cancer)")
    unlinked = raw.scalar(con, "SELECT count(*) FROM raw_mutation WHERE Tumor_Sample_Barcode NOT IN "
                               "(SELECT test_id FROM std_tests_ngs)")
    unlinked_fusions = raw.scalar(con, "SELECT count(*) FROM raw_fusion WHERE Tumor_Sample_Barcode NOT IN "
                                       "(SELECT test_id FROM std_tests_ngs)")
    undated = raw.scalar(con, "SELECT count(*) FROM std_tests_pathology WHERE days_to_entry_event IS NULL")
    notes = [f"cohort {cohort}: {len(files)} files read",
             f"{show(no_cancer, min_count)} patients in the release have no BPC Project cancer: left out at the "
             "flow diagram's first step",
             f"{show(unlinked, min_count)} mutation calls are on samples that aren't the index cancer's sequencing "
             "reports: not used",
             f"{show(unlinked_fusions, min_count)} fusion rows are on samples that aren't the index cancer's sequencing "
             "reports: not used",
             f"{show(undated, min_count)} MSI/MMR results have a withheld date: not known to come before the entry test"]
    notes += [f"{k}: {show(n, min_count)} values didn't convert to the field's type and were read as blank"
              for k, n in failures.items()]
    link = raw.run_checks(con)
    notes += [f"link check {k}: {show(n, min_count)} rows" for k, n in link.items() if n]
    con.close()
    return AdapterResult(tables=tables, profile=prof, input_files=files,
                         source_steps=[("Has a BPC Project cancer", int(no_cancer))], read_notes=notes)
