"""The check to run right after downloading a GENIE BPC release, before anything is analyzed.

    python -m delayed_entry.adapters.genie_bpc.check_download <settings> <clinical folder> <genomic folder> <work folder>
           [--settings-lock <settings_lock.json>] [--rehearsal]

Copied from the CRC BPC project (crcbpc.check_download) and generalized: the release list comes from the
cohort's spec (cohorts/<cohort>/expected_release.yml, with the check limits under `checks`), and every
label the study's rules search for (stage, genes, variant classes, protein changes, this cancer's labels
and site codes) comes from the settings. It asks one question: are these the files the settings and the
adapter expect? It checks the folders, the files (present, and byte for byte what Synapse holds), that
every file can be read, that every column the adapter uses is there with values of the right kind, that
coded columns hold only the documented values, that the units are what the adapter assumes, and that the
files link to each other. Then it builds the standard tables once, to see that they build. It reports
no survival estimate and no count of deaths.

Each check ends as one of:
    pass     as expected
    FAIL     the analysis would break or quietly go wrong: stop here
    review   the settings already say what happens to these rows, so the run can go on, but read the count
    skipped  doesn't apply, or can't run until a failure above is fixed

The report is safe to show to someone else: file names, check names and counts, with any count below the
privacy minimum shown as "<minimum". It never prints a patient, a sample or a line of a file, and an error
from a file reader is cut down to its type and column name, because the full message quotes the line.
"""
from __future__ import annotations

import argparse
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path

import duckdb
import pandas as pd
import yaml

from delayed_entry.adapters import genie_bpc
from delayed_entry.adapters.genie_bpc import genomic_information, raw
from delayed_entry.cohort import drug_matches
from delayed_entry.lock import settings_lock_matches
from delayed_entry.privacy import check_not_cloud
from delayed_entry.settings import Settings, load_settings

PACKAGE_ROOT = Path(__file__).resolve().parents[4]      # the delayed-entry-pipeline folder
REPORT_FILE = "check_download_report.txt"
MADE_UP_MARKER = "truth.csv"                            # only made-up data hold this file
PASS, FAIL, REVIEW, SKIPPED = "pass", "FAIL", "review", "skipped"

# What each link check in link_checks.sql means if it finds rows: FAIL, or review with the answer.
LINK_CHECKS = {
    "patient_id_unique": (FAIL, ""),
    "patient_has_cancer": (REVIEW, "left out at the flow diagram's first step"),
    "cancer_has_patient": (FAIL, ""),
    "cancer_key_unique": (FAIL, ""),
    "report_has_cancer": (REVIEW, "not used: only reports for the index cancer are read"),
    "report_sample_unique": (FAIL, ""),
    "report_number_unique": (FAIL, ""),
    "report_has_panel": (REVIEW, "a first report with no panel file fails the genotype rule"),
    "mutation_has_report": (REVIEW, "calls for a sample outside the report file are not used"),
    "call_gene_on_panel": (REVIEW, "the data guide warns gene names can differ between sources; see the gene checks"),
    "fusion_has_report": (REVIEW, "fusions for a sample outside the report file are not used"),
    "other_cancer_has_patient": (REVIEW, "not used"),
    "other_cancer_not_index": (FAIL, ""),
    "other_cancer_flag_blank": (REVIEW, "only the index file defines BPC Project cancers; see the check on site codes"),
    "pathology_has_patient": (REVIEW, "not used"),
    "pathology_key_unique": (REVIEW, "a repeated report can't change an any-result rule"),
    "regimen_number_unique": (REVIEW, "only the index cancer matters; see the check on index-cancer regimens"),
    "report_seq_year_is_a_year": (FAIL, ""),
}
# Why the checks on other primaries of this cancer are skipped when no analysis reads them.
NO_OTHER_PRIMARY_RULE = ("no analysis that runs leaves out patients with another primary of this cancer "
                         "(other_primary_same_cancer), so these labels decide nothing")
# The standard field each settings rule may name, and the raw column it comes from.
RAW_COLUMN = {"stage_at_diagnosis": ("raw_cancer", "stage_dx"), "histology": ("raw_cancer", "ca_hist_adeno_squamous"),
              "primary_site_code": ("raw_cancer", "ca_d_site"), "sex": ("raw_patient", "naaccr_sex_code"),
              "center": ("raw_patient", "institution")}


@dataclass(frozen=True)
class Check:
    section: str
    name: str
    status: str
    detail: str = ""


class Shower:
    """Counts as the report prints them: below the privacy minimum as "<minimum"."""

    def __init__(self, min_count: int):
        self.min_count = min_count

    def __call__(self, n: int) -> str:
        return f"<{self.min_count}" if 1 <= n < self.min_count else f"{n:,}"


def shape(text: str) -> str:
    """A value's form with its content removed: digits become 9 and letters A, so ">32485" is ">99999"."""
    return re.sub(r"[A-Za-z]", "A", re.sub(r"[0-9]", "9", text))[:24]


def safe_error(e: Exception) -> str:
    """What may be said about an error from reading a file. The full message quotes the line it broke on."""
    column = re.search(r'converting column "?([\w.]+)"?', str(e))
    return type(e).__name__ + (f" in column {column.group(1)}" if column else "")


def plain(text: str) -> str:
    return text.encode("ascii", "replace").decode("ascii")


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def folder_of(row_folder: str, clinical: Path, genomic: Path) -> Path:
    return clinical if row_folder == "clinical" else genomic


# ---- What the settings search for ----
def searched_labels(settings: Settings) -> dict:
    """Every label the settings' rules match by exact text."""
    genes, classes, variants, fusion_genes, field_values = set(), set(), set(), set(), []
    for g in settings.read("genotype.groups"):
        for m in g.get("any_of", []):
            if m["kind"] == "fusion":
                fusion_genes |= set(m["genes"])     # matched in the fusion file
            if m["kind"] != "mutation":
                continue            # only mutation rules are matched in the mutation file
            genes |= set(m["genes"])
            classes |= set(m.get("variant_classes", []))
            variants |= set(m.get("variants", []))
    for r in settings.read("cohort.rules"):
        for c in r["all_of"]:
            if c.get("field") in RAW_COLUMN and c["op"] in ("in", "starts_with_any"):
                field_values.append((c["field"], c.get("normalize", "none"), [str(v) for v in c["values"]]))
    other_primary = any(a["run"] and a["kind"] == "biomarker" and "other_primary_same_cancer" in a["exclude_if"]
                        for a in settings.read("analyses"))
    return {"genes": sorted(genes), "classes": sorted(classes),
            "variants": sorted(variants), "fusion_genes": sorted(fusion_genes), "fields": field_values,
            "fusion_panels": settings.read("genotype.coverage.fusion_panels"),
            "drugs": settings.read("genotype.reference_drugs"), "other_primary_used": other_primary,
            "cancer_labels": list(settings.read("source.genie_bpc.cancer_type_labels")),
            "site_pattern": settings.read("source.genie_bpc.cancer_site_pattern")}


# ---- A. Folders and safety ----
def check_folders(clinical: Path, genomic: Path, work: Path, rehearsal: bool, settings_path: Path,
                  settings_lock: Path | None) -> list[Check]:
    s = "A. Folders and safety"
    out = []
    for label, p in (("clinical", clinical), ("genomic", genomic)):
        out.append(Check(s, f"The {label} folder exists", PASS if p.is_dir() else FAIL, "" if p.is_dir() else p.name))
    for label, p in (("clinical", clinical), ("genomic", genomic), ("work", work)):
        try:
            check_not_cloud(p.resolve())
            out.append(Check(s, f"The {label} folder is not cloud-synced", PASS))
        except ValueError:
            out.append(Check(s, f"The {label} folder is not cloud-synced", FAIL, "row-level data must stay on local disk"))
    inside = PACKAGE_ROOT in work.resolve().parents
    out.append(Check(s, "The work folder is outside the package folder", FAIL if inside else PASS))
    if rehearsal:
        for name in ("The data folders are outside the package folder", "The settings were locked before the download",
                     "No made-up data in the folders"):
            out.append(Check(s, name, SKIPPED, "rehearsal"))
        return out
    inside = [p.name for p in (clinical, genomic) if PACKAGE_ROOT in p.resolve().parents]
    out.append(Check(s, "The data folders are outside the package folder", FAIL if inside else PASS, ", ".join(inside)))
    if settings_lock is None or not settings_lock.exists():
        out.append(Check(s, "The settings were locked before the download", FAIL,
                         "no settings lock: lock the settings (delayed_entry.lock.lock_settings) before downloading"))
    else:
        ok = settings_lock_matches(settings_lock, settings_path)
        out.append(Check(s, "The settings were locked before the download", PASS if ok else FAIL,
                         "" if ok else "the settings file has changed since it was locked"))
    made_up = [p.name for p in (clinical, genomic) if (p / MADE_UP_MARKER).exists()]
    out.append(Check(s, "No made-up data in the folders", FAIL if made_up else PASS,
                     f"{MADE_UP_MARKER} is there" if made_up else ""))
    return out


# ---- B. The files themselves ----
def release_files(cohort: str, expected: dict) -> list[dict] | None:
    """The release's files with their sizes and checksums: listed in the release list itself, or in the file
    it names (written by scripts/synapse_release_list.py). None if that file hasn't been made yet."""
    files = expected["files"]
    if isinstance(files, list):
        return files
    path = raw.cohort_dir(cohort) / files
    return yaml.safe_load(path.read_text(encoding="utf-8"))["files"] if path.exists() else None


def check_genomic_information(genomic: Path, rule: dict | None) -> list[Check]:
    """A center's BED file isn't read. A genomic_information file is read for panel sizes when the settings give the
    rule (source.genie_bpc.panel_sizes): it must then be there, once, with the columns the rule reads."""
    s = "B. The files themselves"
    beds = sorted(p.name for p in genomic.glob("*.bed")) if genomic.is_dir() else []
    infos = genomic_information.find(genomic) if genomic.is_dir() else []
    out = [Check(s, "No center BED file (reading one isn't written; a genomic_information file can give panel sizes)",
                 FAIL if beds else PASS, ", ".join(beds))]
    name = "The genomic_information file, for panel sizes"
    if len(infos) > 1:
        out.append(Check(s, name, FAIL, f"more than one: {', '.join(p.name for p in infos)}"))
    elif rule is None:
        out.append(Check(s, name, REVIEW if infos else PASS,
                         f"{infos[0].name} is there but not read: the settings give no rule (source.genie_bpc.panel_sizes)"
                         if infos else "not needed: the settings work out no panel sizes"))
    elif not infos:
        out.append(Check(s, name, FAIL, "the settings say how to work out panel sizes from it "
                                        "(source.genie_bpc.panel_sizes), but the genomic folder has none"))
    else:
        bad = genomic_information.problems(genomic_information.read(infos[0]))
        out.append(Check(s, name, FAIL if bad else PASS, "; ".join(bad) or f"{infos[0].name}, read by the settings' rule"))
    return out


def check_files(clinical: Path, genomic: Path, cohort: str, expected: dict, rehearsal: bool,
                size_rule: dict | None = None) -> list[Check]:
    s = "B. The files themselves"
    out = []
    needed = {name: folder_of(folder, clinical, genomic) for name, folder, _ in raw.raw_files(cohort).values()}
    missing = [n for n, f in needed.items() if not (f / n).exists()]
    out.append(Check(s, "Every file the adapter reads is there", FAIL if missing else PASS, ", ".join(missing)))
    panels = sorted(genomic.glob(raw.PANEL_FILES)) if genomic.is_dir() else []
    out.append(Check(s, "The gene-panel files are there", PASS if panels else FAIL, f"{len(panels)} files"))
    out += check_genomic_information(genomic, size_rule)
    if rehearsal:
        out.append(Check(s, "Size and checksum of each file match Synapse", SKIPPED, "rehearsal"))
        return out
    files = release_files(cohort, expected)
    if files is None:
        out.append(Check(s, "The release list holds each file's size and checksum", FAIL,
                         f"no {expected['files']} yet: run scripts/synapse_release_list.py {cohort}"))
        return out
    for used, label, bad in ((True, "file the adapter reads", FAIL), (False, "other file in the release", REVIEW)):
        rows = [f for f in files if f["used"] is used]
        absent, wrong_size, wrong_md5 = [], [], []
        for f in rows:
            p = folder_of(f["folder"], clinical, genomic) / f["name"]
            if not p.exists():
                absent.append(f["name"])
            elif p.stat().st_size != f["bytes"]:
                wrong_size.append(f["name"])
            elif md5(p) != f["md5"]:
                wrong_md5.append(f["name"])
        for name, found in (("is there", absent), ("has the size Synapse reports", wrong_size),
                            ("has the checksum Synapse reports", wrong_md5)):
            out.append(Check(s, f"Every {label} {name} ({len(rows)} files)", bad if found else PASS, ", ".join(found)))
    for label, p in (("clinical", clinical), ("genomic", genomic)):
        found = bool(list(p.glob("manifest_*.csv"))) if p.is_dir() else False
        out.append(Check(s, f"The {label} download left its checksum manifest", PASS if found else FAIL,
                         "" if found else "needed for the provenance log"))
    return out


# ---- C. Reading the files ----
def text_table(table: str) -> str:
    return f"text_{table}"


def read_as_text(con: duckdb.DuckDBPyConnection, table: str, path: Path) -> None:
    """Read a file with every column as text, so a bad value can be counted instead of stopping the read."""
    delim = "\t" if path.suffix == ".txt" else ","
    con.execute(f"CREATE OR REPLACE TABLE {text_table(table)} AS SELECT * FROM read_csv({raw.sql_text(path.as_posix())}, "
                f"header = true, delim = {raw.sql_text(delim)}, all_varchar = true, nullstr = '', sample_size = -1)")


def columns_of(con, table: str) -> set[str]:
    return {r[0] for r in con.execute(f"SELECT column_name FROM (DESCRIBE {text_table(table)})").fetchall()}


def check_reading(con, clinical: Path, genomic: Path, cohort: str, show: Shower) -> tuple[list[Check], set[str]]:
    s = "C. Reading the files"
    not_utf8, unreadable, empty, missing_cols, bad_values, na_text = [], [], [], [], [], []
    readable = set()
    for table, (name, folder, types) in raw.raw_files(cohort).items():
        path = folder_of(folder, clinical, genomic) / name
        if not path.exists():
            continue
        try:
            path.read_bytes().decode("utf-8")
        except UnicodeDecodeError:
            not_utf8.append(name)
            continue
        try:
            read_as_text(con, table, path)
        except duckdb.Error as e:
            unreadable.append(f"{name} ({safe_error(e)})")
            continue
        readable.add(table)
        t = text_table(table)
        if raw.scalar(con, f"SELECT count(*) FROM {t}") == 0:
            empty.append(name)
        have = columns_of(con, table)
        absent = [c for c in types if c not in have]
        if absent:
            missing_cols.append(f"{name}: {', '.join(absent)}")
        for col, kind in types.items():
            if col not in have:
                continue
            q = f'"{col}"'
            n_na = raw.scalar(con, f"SELECT count(*) FROM {t} WHERE trim({q}) = 'NA'")
            if n_na:
                na_text.append(f"{col}: {show(n_na)}")
            if kind == "VARCHAR":
                continue
            bad = f"nullif(trim({q}), '') IS NOT NULL AND trim({q}) <> 'NA' AND TRY_CAST(trim({q}) AS {kind}) IS NULL"
            n = raw.scalar(con, f"SELECT count(*) FROM {t} WHERE {bad}")
            if n:
                forms = [plain(shape(r[0])) for r in con.execute(f"SELECT DISTINCT {q} FROM {t} WHERE {bad} LIMIT 5").fetchall()]
                bad_values.append(f"{col} should be {kind}: {show(n)} values like {sorted(set(forms))}")
    out = [Check(s, "Every file is valid UTF-8 text", FAIL if not_utf8 else PASS, ", ".join(not_utf8)),
           Check(s, "Every file can be read", FAIL if unreadable else PASS, "; ".join(unreadable)),
           Check(s, "No file is empty", FAIL if empty else PASS, ", ".join(empty)),
           Check(s, "Every column the adapter uses is there", FAIL if missing_cols else PASS, "; ".join(missing_cols)),
           Check(s, "Every value converts to the type the adapter reads it as", FAIL if bad_values else PASS,
                 "; ".join(bad_values)),
           Check(s, "No column uses the text NA for a missing value", REVIEW if na_text else PASS,
                 "read as missing, as genieBPC does: " + "; ".join(na_text) if na_text else "")]
    return out, readable


# ---- D. Gene panels ----
def check_panels(genomic: Path, expected: dict, genes: list[str], rehearsal: bool,
                 fusion_panels: dict | None = None) -> list[Check]:
    """The panel files. For a gene in `fusion_panels` (the settings' panels that detect its fusions), a panel
    reads it only if it is on that list too."""
    s = "D. Gene panels"
    fusion_panels = {g: ps for g, ps in (fusion_panels or {}).items() if g in genes}
    bad, lacking, full = [], [], []
    names = []
    for path in sorted(genomic.glob(raw.PANEL_FILES)) if genomic.is_dir() else []:
        try:
            panel, panel_genes = raw.read_panel_file(path)
        except (KeyError, ValueError, UnicodeDecodeError):
            bad.append(path.name)
            continue
        names.append(panel)
        if path.name != f"data_gene_panel_{panel}.txt" or len(panel_genes) < 2:
            bad.append(path.name)
        elif set(genes) <= set(panel_genes) and all(panel in ps for ps in fusion_panels.values()):
            full.append(panel)
        else:
            lacking.append(panel)
    listed = ", ".join(genes)
    out = [Check(s, "Each panel file names its panel and lists its genes", FAIL if bad else PASS, ", ".join(bad)),
           Check(s, f"At least one panel reads every gene the groups need ({listed})", PASS if full else FAIL,
                 f"{len(full)} panels do" + (" (for fusions, by the settings' list)" if fusion_panels else ""))]
    if fusion_panels:
        unknown = sorted({p for ps in fusion_panels.values() for p in ps} - set(names))
        out.append(Check(s, "Every panel the settings list as detecting fusions has a panel file", FAIL if unknown else PASS,
                         "not a panel here, so its patients would be unassignable (check the spelling): "
                         + ", ".join(unknown) if unknown else ""))
    out.append(Check(s, f"Every panel reads {listed}", REVIEW if lacking else PASS,
                     "tumors on these fail the genotype rule: " + ", ".join(lacking) if lacking else ""))
    if rehearsal:
        out.append(Check(s, "The panels are the ones in the release", SKIPPED, "rehearsal"))
    else:
        differ = sorted(set(names) ^ set(expected["panels"]))
        out.append(Check(s, "The panels are the ones in the release", FAIL if differ else PASS, ", ".join(differ)))
    return out


# ---- E. Documented values ----
def check_values(con, readable: set[str], expected: dict, show: Shower) -> list[Check]:
    s = "E. Coded values"
    blocking, review, checked = [], [], 0
    have = {t: columns_of(con, t) for t in readable}
    for table, columns in expected["values"].items():
        if table not in readable:
            continue
        for col, spec in columns.items():
            if col not in have[table]:
                continue
            checked += 1
            allowed = ", ".join(raw.sql_text(v.strip()) for v in spec["allowed"])
            rows = con.execute(f"""SELECT trim("{col}"), count(*) FROM {text_table(table)}
                                   WHERE nullif(trim("{col}"), '') IS NOT NULL AND trim("{col}") <> 'NA'
                                     AND trim("{col}") NOT IN ({allowed}) GROUP BY 1 ORDER BY 1""").fetchall()
            if rows:
                listed = "; ".join(f"{plain(str(v))[:60]!r} {show(n)}" for v, n in rows[:5])
                (blocking if spec["stops"] else review).append(f"{table}.{col}: {listed}")
    out = [Check(s, f"Coded columns the rules read hold only documented values ({checked} columns)",
                 FAIL if blocking else PASS, " | ".join(blocking))]
    out.append(Check(s, "Descriptive columns hold only documented values", REVIEW if review else PASS, " | ".join(review)))
    return out


def has_value(con, table: str, column: str, values, normalize: str = "none") -> list[str]:
    """Which of `values` never appear in the column, compared the way the rule compares them."""
    if column not in columns_of(con, table):
        return list(values)
    col = f'"{column}"'
    if normalize == "trim":
        col = f"trim({col})"
    elif normalize == "letters_digits_upper":
        col = f"regexp_replace(upper({col}), '[^A-Z0-9]', '', 'g')"
    found = {r[0] for r in con.execute(f"SELECT DISTINCT {col} FROM {text_table(table)}").fetchall()}
    return [v for v in values if v not in found]


def check_labels(con, readable: set[str], labels: dict, rehearsal: bool, show: Shower) -> list[Check]:
    """Labels the rules search for by exact text: if the spelling differed, a rule would quietly match nobody."""
    s = "E. Coded values"
    out = []
    if {"raw_cancer", "raw_mutation", "raw_pathology", "raw_patient"} - readable:
        return [Check(s, "The labels the rules search for are in the files", SKIPPED, "a file above can't be read")]
    for field, normalize, values in labels["fields"]:
        table, column = RAW_COLUMN[field]
        gone = has_value(con, table, column, values, normalize)
        out.append(Check(s, f"The labels a cohort rule matches on {field} appear as the rule compares them",
                         FAIL if gone else PASS, ", ".join(gone)))
    gone = has_value(con, "raw_mutation", "Hugo_Symbol", labels["genes"], "trim")
    out.append(Check(s, "The genes the groups name appear in the mutation file", FAIL if gone else PASS, ", ".join(gone)))
    if labels["classes"]:
        gone = has_value(con, "raw_mutation", "Variant_Classification", labels["classes"], "trim")
        status = FAIL if len(gone) == len(labels["classes"]) else (REVIEW if gone else PASS)
        out.append(Check(s, "The variant classes the groups count appear in the mutation file", status,
                         "not in the file, or spelled differently (a rule would skip them): " + ", ".join(gone) if gone else ""))
    if labels["variants"]:
        gone = has_value(con, "raw_mutation", "HGVSp_Short", labels["variants"], "trim")
        out.append(Check(s, "The protein changes the groups name are written in that form", FAIL if gone else PASS,
                         ", ".join(gone)))
    out += check_fusion_labels(con, readable, labels, show)
    columns = [c for c in columns_of(con, "raw_pathology") if re.fullmatch(r"path_ca_type\d+", c)]
    found = set()
    for c in columns:
        found |= {r[0] for r in con.execute(f'SELECT DISTINCT trim("{c}") FROM {text_table("raw_pathology")}').fetchall()}
    out.append(Check(s, "A specimen label of this cancer appears in the pathology file",
                     PASS if found & set(labels["cancer_labels"]) else FAIL))
    name = "Other primaries of this cancer in the non-index file use the labels the settings match"
    if not labels["other_primary_used"]:
        out.append(Check(s, name, SKIPPED, NO_OTHER_PRIMARY_RULE))
    elif "raw_cancer_other" not in readable:
        out.append(Check(s, name, SKIPPED, "the non-index cancer file can't be read"))
    elif "ca_type" not in columns_of(con, "raw_cancer_other"):
        out.append(Check(s, name, FAIL, "the file has no ca_type column"))
    else:
        types = {r[0] for r in con.execute(f'SELECT DISTINCT trim("ca_type") FROM {text_table("raw_cancer_other")}').fetchall()}
        # a near spelling: a label sharing the first four letters of a word in this cancer's labels ("colo", "rect")
        stems = {w.lower()[:4] for lab in labels["cancer_labels"] for w in re.findall(r"[A-Za-z]{4,}", lab)} - {"canc"}
        near = sorted(t for t in types if t and t not in labels["cancer_labels"] and any(w in t.lower() for w in stems))
        out.append(Check(s, name, REVIEW if near else PASS,
                         "spelled otherwise, so the exclusion would keep these patients: " + ", ".join(near) if near else ""))
    out.append(check_sites(con, readable, labels, rehearsal, show))
    return out


def gene_in_text(gene: str) -> str:
    """A regular expression finding the gene's name in a fusion's description as a whole name ("ALK" in
    "EML4-ALK fusion", not in "ALKBH5")."""
    return f"(^|[^A-Za-z0-9]){re.escape(gene)}([^A-Za-z0-9]|$)"


def check_fusion_labels(con, readable: set[str], labels: dict, show: Shower) -> list[Check]:
    """The genes the groups count as fusions, in the fusion file. A rule finds a fusion by its gene rows, so a
    fusion whose description names a gene that has no row of its own on that sample would be missed."""
    s = "E. Coded values"
    genes = labels["fusion_genes"]
    if not genes:
        return []
    if "raw_fusion" not in readable:
        return [Check(s, "The genes the groups count as fusions appear in the fusion file", FAIL,
                      "the groups count fusions, but the cohort's files include no fusion file that can be read")]
    gone = has_value(con, "raw_fusion", "Hugo_Symbol", genes, "trim")
    status = FAIL if len(gone) == len(genes) else (REVIEW if gone else PASS)
    out = [Check(s, "The genes the groups count as fusions appear in the fusion file", status,
                 "not in the file, or spelled differently: " + ", ".join(gone) if gone else "")]
    t = text_table("raw_fusion")
    missed = []
    for g in genes:
        n = raw.scalar(con, f"""
            SELECT count(DISTINCT f.Tumor_Sample_Barcode) FROM {t} AS f
            WHERE regexp_matches(coalesce(f.Fusion, ''), {raw.sql_text(gene_in_text(g))})
              AND NOT EXISTS (SELECT 1 FROM {t} AS o WHERE o.Tumor_Sample_Barcode = f.Tumor_Sample_Barcode
                                                       AND trim(o.Hugo_Symbol) = {raw.sql_text(g)})""")
        if n:
            missed.append(f"{g}: {show(n)} samples")
    out.append(Check(s, "Every fusion named for a gene the groups count has a row for that gene",
                     REVIEW if missed else PASS,
                     "described with the gene but listed under the partner only, so a rule on the gene misses them: "
                     + "; ".join(missed) if missed else ""))
    return out


def check_drug_labels(con, readable: set[str], labels: dict, reads_drugs: bool) -> list[Check]:
    """The drugs the settings list for the count in the reference group, found in the regimen file the way the
    screen finds them. A name written otherwise would quietly count nobody."""
    s = "E. Coded values"
    drugs = labels["drugs"]
    if not drugs:
        return []
    name = "The drugs the settings list appear in the regimen file as the screen finds them"
    if not reads_drugs:
        return [Check(s, name, REVIEW, "the cohort's spec doesn't read regimen_drugs (raw_files.yml), so the screen "
                                       "can't count them: stated as a limitation")]
    if "raw_regimen" not in readable or "regimen_drugs" not in columns_of(con, "raw_regimen"):
        return [Check(s, name, SKIPPED, "the regimen file or its regimen_drugs column can't be read (see above)")]
    found = pd.Series([r[0] for r in con.execute(f"SELECT DISTINCT regimen_drugs FROM {text_table('raw_regimen')} "
                                                 "WHERE regimen_drugs IS NOT NULL").fetchall()], dtype="string")
    listed = list(dict.fromkeys(d for names in drugs.values() for d in names))
    gone = [d for d in listed if not drug_matches(found, [d]).any()]
    return [Check(s, name, REVIEW if gone else PASS,
                  "not in the file, or written otherwise, so the count misses them: " + ", ".join(gone) if gone else "")]


def check_sites(con, readable: set[str], labels: dict, rehearsal: bool, show: Shower) -> Check:
    """A biomarker analysis may leave out patients with another primary of this cancer, found by its type. A cancer
    whose type is blank or "Other" would be missed, so check every non-index cancer at this cancer's site codes is
    typed as this cancer."""
    s = "E. Coded values"
    name = "Non-index cancers at this cancer's site codes are typed as this cancer"
    if not labels["other_primary_used"]:
        return Check(s, name, SKIPPED, NO_OTHER_PRIMARY_RULE)
    if "raw_cancer_other" not in readable:
        return Check(s, name, SKIPPED, "the non-index cancer file can't be read")
    if not {"ca_type", "ca_d_site"} <= columns_of(con, "raw_cancer_other"):
        if rehearsal:
            return Check(s, name, SKIPPED, "the made-up file has no site codes")
        return Check(s, name, FAIL, "the file has no ca_type or ca_d_site column")
    listed = ", ".join(raw.sql_text(t) for t in labels["cancer_labels"])
    n = con.execute(f"""SELECT count(*) FROM {text_table('raw_cancer_other')}
                        WHERE regexp_matches(trim(ca_d_site), {raw.sql_text(labels['site_pattern'])})
                          AND coalesce(trim(ca_type), '') NOT IN ({listed})""").fetchone()[0]
    return Check(s, name, FAIL if n else PASS, f"{show(n)} cancers: the exclusion would keep these patients" if n else "")


# ---- F. Links, units and a trial build ----
def unit_check(s: str, name: str, con, query: str, limit: float, show: Shower, too_many: str = FAIL) -> Check:
    """`query` returns (rows that disagree, rows compared). More than `limit` disagreeing is `too_many`."""
    bad, total = con.execute(query).fetchall()[0]
    if total == 0:
        return Check(s, name, FAIL, "no rows to compare")
    return Check(s, name, too_many if bad / total > limit else (REVIEW if bad else PASS),
                 f"{show(bad)} of {show(total)} rows disagree" if bad else "")


def check_links(con, limits: dict, genes: list[str], show: Shower) -> list[Check]:
    s = "F. Links between the files"
    counts = raw.run_checks(con)
    unlisted = sorted(set(counts) - set(LINK_CHECKS))
    out = [Check(s, "Every link check has a stated meaning", FAIL if unlisted else PASS, ", ".join(unlisted))]
    for name, n in counts.items():
        status, why = LINK_CHECKS.get(name, (FAIL, ""))
        out.append(Check(s, name, PASS if n == 0 else status, f"{show(n)} rows" + (f": {why}" if why else "") if n else ""))
    with_calls, reports = con.execute("""
        SELECT count(DISTINCT r.cpt_genie_sample_id) FILTER (WHERE m.Tumor_Sample_Barcode IS NOT NULL),
               count(DISTINCT r.cpt_genie_sample_id)
        FROM raw_report AS r LEFT JOIN raw_mutation AS m ON m.Tumor_Sample_Barcode = r.cpt_genie_sample_id""").fetchall()[0]
    linked = reports > 0 and with_calls / reports >= limits["min_call_link"]
    out.append(Check(s, f"Most report samples have mutation calls (at least {limits['min_call_link']:.0%})",
                     PASS if linked else FAIL,
                     "" if linked else "the sample IDs may not match, which would make every tumor look wild type"))
    gene_list = ", ".join(raw.sql_text(g) for g in genes)
    n = raw.scalar(con, f"""
        SELECT count(*) FROM raw_mutation AS m
        JOIN raw_report AS r ON r.cpt_genie_sample_id = m.Tumor_Sample_Barcode
        LEFT JOIN raw_panel_gene AS g ON g.panel = r.cpt_seq_assay_id AND g.gene = m.Hugo_Symbol
        WHERE m.Hugo_Symbol IN ({gene_list}) AND g.gene IS NULL""")
    out.append(Check(s, "Calls in the groups' genes are on panels that list those genes", REVIEW if n else PASS,
                     f"{show(n)} calls: not used, since the panel can't confirm them" if n else ""))
    n = raw.scalar(con, """
        SELECT count(*) FROM (
            SELECT r.record_id FROM raw_regimen AS r
            JOIN raw_cancer AS c ON c.record_id = r.record_id AND c.ca_seq = r.ca_seq
            GROUP BY r.record_id, r.ca_seq, r.regimen_number_within_cancer HAVING count(*) > 1)""")
    out.append(Check(s, "An index cancer's regimens have different numbers", FAIL if n else PASS, f"{show(n)} cases" if n else ""))
    n = raw.scalar(con, """
        SELECT count(*) FROM raw_cancer AS c
        LEFT JOIN raw_report AS r ON r.record_id = c.record_id AND r.ca_seq = c.ca_seq WHERE r.record_id IS NULL""")
    out.append(Check(s, "Every index cancer has a sequencing report (the data guide says so)", REVIEW if n else PASS,
                     f"{show(n)} cancers: they fail the has-test rule" if n else ""))
    return out


def check_units(con, limits: dict, show: Shower) -> list[Check]:
    s = "G. Ranges and units"
    lim = limits["max_unit_mismatch"]
    lo, hi = limits["age_range"]
    by0, by1 = limits["birth_year_range"]
    return [
        unit_check(s, f"Report months match report days (days / 30.4375, within {limits['report_months_tolerance']} month)", con, f"""
            SELECT count(*) FILTER (WHERE abs(dx_cpt_rep_mos - dx_cpt_rep_days / 30.4375) > {limits['report_months_tolerance']}), count(*)
            FROM raw_report WHERE dx_cpt_rep_mos IS NOT NULL AND dx_cpt_rep_days IS NOT NULL""", lim, show),
        unit_check(s, f"Age at diagnosis is in years (matches days from birth, within {limits['age_years_tolerance']} years)", con, f"""
            SELECT count(*) FILTER (WHERE abs(age_dx - dob_ca_dx_days / 365.25) > {limits['age_years_tolerance']}), count(*)
            FROM raw_cancer WHERE age_dx IS NOT NULL AND dob_ca_dx_days IS NOT NULL""", lim, show),
        unit_check(s, "Report day from birth = diagnosis day from birth + report day from diagnosis "
                      f"(within {limits['report_days_tolerance']} days)", con, f"""
            SELECT count(*) FILTER (WHERE abs(r.dob_cpt_report_days - c.dob_ca_dx_days - r.dx_cpt_rep_days) > {limits['report_days_tolerance']}), count(*)
            FROM raw_report AS r JOIN raw_cancer AS c ON c.record_id = r.record_id AND c.ca_seq = r.ca_seq
            WHERE r.dob_cpt_report_days IS NOT NULL AND c.dob_ca_dx_days IS NOT NULL AND r.dx_cpt_rep_days IS NOT NULL""",
                   lim, show, too_many=REVIEW),
        unit_check(s, f"Age at diagnosis is between {lo} and {hi}", con, f"""
            SELECT count(*) FILTER (WHERE age_dx < {lo} OR age_dx > {hi}), count(*) FROM raw_cancer WHERE age_dx IS NOT NULL""",
                   lim, show),
        unit_check(s, f"Year of birth is between {by0} and {by1}", con, f"""
            SELECT count(*) FILTER (WHERE birth_year < {by0} OR birth_year > {by1}), count(*) FROM raw_patient
            WHERE birth_year IS NOT NULL""", lim, show),
        unit_check(s, "Survival times are not negative", con, """
            SELECT count(*) FILTER (WHERE tt_os_dx_days < 0), count(*) FROM raw_cancer WHERE tt_os_dx_days IS NOT NULL""",
                   lim, show),
    ]


def check_build(clinical: Path, genomic: Path, settings: Settings) -> list[Check]:
    """Build the standard tables once, to see that they build. Nothing is counted."""
    s = "H. A trial build"
    try:
        genie_bpc.read(clinical, genomic, settings)
        return [Check(s, "The standard tables build, and every result label is known", PASS)]
    except Exception as e:      # any failure here is what this check exists to find
        return [Check(s, "The standard tables build, and every result label is known", FAIL, safe_error(e))]


# ---- All of it ----
def run_checks(clinical: Path, genomic: Path, work: Path, settings_path: Path, settings_lock: Path | None = None,
               rehearsal: bool = False, expected_yml: Path | None = None) -> list[Check]:
    """Every check, in order. `expected_yml` replaces the cohort's release list (tests use a made-up one)."""
    settings = load_settings(settings_path)
    cohort = settings.read("source.genie_bpc.cohort")
    own = yaml.safe_load((raw.cohort_dir(cohort) / "expected_release.yml").read_text(encoding="utf-8"))
    expected = yaml.safe_load(expected_yml.read_text(encoding="utf-8")) if expected_yml else own
    limits = expected.get("checks") or own["checks"]
    show = Shower(settings.read("privacy.min_count"))
    labels = searched_labels(settings)
    checks = check_folders(clinical, genomic, work, rehearsal, settings_path, settings_lock)
    size_rule = settings.read("source.genie_bpc.panel_sizes") if settings.has("source.genie_bpc.panel_sizes") else None
    checks += check_files(clinical, genomic, cohort, expected, rehearsal, size_rule)
    if not (clinical.is_dir() and genomic.is_dir()):
        return checks
    scratch = duckdb.connect()
    reading, readable = check_reading(scratch, clinical, genomic, cohort, show)
    checks += reading
    checks += check_panels(genomic, expected, list(settings.read("genotype.coverage.required")), rehearsal,
                           labels["fusion_panels"])
    checks += check_values(scratch, readable, expected, show)
    checks += check_labels(scratch, readable, labels, rehearsal, show)
    regimen_types = raw.raw_files(cohort).get("raw_regimen", ("", "", {}))[2]
    checks += check_drug_labels(scratch, readable, labels, "regimen_drugs" in regimen_types)
    scratch.close()
    s = "F. Links between the files"
    if any(c.status == FAIL for c in checks if c.section.startswith(("B.", "C."))):
        checks.append(Check(s, "The links, units and trial build", SKIPPED, "they can't run until the failures above are fixed"))
        return checks
    try:
        con, _ = raw.load(cohort, clinical, genomic)
    except Exception as e:      # any failure here is what this check exists to find
        checks.append(Check(s, "The adapter can load every file", FAIL, safe_error(e)))
        return checks
    checks.append(Check(s, "The adapter can load every file", PASS))
    checks += check_links(con, limits, labels["genes"], show)
    checks += check_units(con, limits, show)
    con.close()
    checks += check_build(clinical, genomic, settings)
    return checks


def verdict(checks: list[Check]) -> str:
    failed, review = (sum(c.status == s for c in checks) for s in (FAIL, REVIEW))
    if failed:
        return f"STOP: {failed} checks failed. Don't go on until each one is understood and fixed."
    if review:
        return f"PASSED, with {review} items to read. Read each one and note it in the run log, then go on."
    return "PASSED: every check passed. Go on to the next step."


def report(checks: list[Check], title: str = "Check of the downloaded files") -> str:
    lines = [title, "=" * len(title)]
    section = None
    for c in checks:
        if c.section != section:
            section = c.section
            lines += ["", section]
        lines.append(f"  {c.status:<8} {plain(c.name)}" + (f"  [{plain(c.detail)}]" if c.detail else ""))
    tally = {s: sum(c.status == s for c in checks) for s in (PASS, FAIL, REVIEW, SKIPPED)}
    lines += ["", f"{tally[PASS]} passed, {tally[FAIL]} failed, {tally[REVIEW]} to read, {tally[SKIPPED]} skipped", "",
              "RESULT: " + verdict(checks), ""]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Check a downloaded GENIE BPC release before anything is analyzed.")
    ap.add_argument("settings", type=Path, help="the study's settings file")
    ap.add_argument("clinical", type=Path, help="folder with the clinical files")
    ap.add_argument("genomic", type=Path, help="folder with the genomic files")
    ap.add_argument("work", type=Path, help="local folder for the report (not cloud-synced)")
    ap.add_argument("--settings-lock", type=Path, default=None, help="the settings lock made before the download")
    ap.add_argument("--rehearsal", action="store_true", help="made-up data: skip the checks only a real release can pass")
    a = ap.parse_args(argv)
    work = check_not_cloud(a.work.resolve())
    work.mkdir(parents=True, exist_ok=True)
    checks = run_checks(a.clinical.resolve(), a.genomic.resolve(), work, a.settings.resolve(), a.settings_lock, a.rehearsal)
    text = report(checks)
    (work / REPORT_FILE).write_bytes(text.encode("utf-8"))
    print(text)
    if any(c.status == FAIL for c in checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
