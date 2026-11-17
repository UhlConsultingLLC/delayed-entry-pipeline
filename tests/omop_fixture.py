"""A made-up OMOP CDM database, built from the made-up GENIE BPC fixture's standard tables, and its mapping.

Every concept made up for it is a local concept (2,000,000,000 and above, the range OMOP leaves to local
codes). Each diagnosis is placed on July 1 of the year of diagnosis (2012 where the year is unknown), and every
other date is that date plus its interval, rounded down to a whole day. Decoys the adapter must pass over are
added: a later diagnosis row and an earlier one of another condition, a stage farther from diagnosis, a variant on
a procedure that isn't a sequencing test, and a drug outside the drug concept set.
"""
from __future__ import annotations

import datetime as dt
import math
import re
from pathlib import Path

import pandas as pd
import yaml

from delayed_entry.adapters.omop import AMINO_ACIDS
from delayed_entry.tables import StandardTables

DIAGNOSIS, OTHER_CONDITION = 2000000001, 2000000002
STAGE = {"Stage I": 2000000101, "Stage I-III NOS": 2000000102, "Stage II": 2000000103, "Stage III": 2000000104,
         "Stage IV": 2000000105}
MSI_TEST, MMR_TEST = 2000000201, 2000000202
READING = {"MSI-H": 2000000301, "MSI-L": 2000000302, "MSS": 2000000303, "dMMR": 2000000304, "pMMR": 2000000305,
           "indeterminate": 2000000306}
FLAGS = {"censored_at_age_limit": 2000000401, "multiple_cohort_cancers": 2000000402,
         "other_primary_same_cancer": 2000000403, "tumor_location": 2000000404}
SEQUENCING, OTHER_PROCEDURE = 2000000501, 2000000502
CANCER_DRUG, OTHER_DRUG = 2000000601, 2000000602
MORPHOLOGY = {"Adenocarcinoma": "8140/3", "Carcinoma": "8010/3", "Other histologies/mixed tumor": "8255/3",
              "Squamous cell": "8070/3"}
ONE_TO_THREE = {v: k for k, v in AMINO_ACIDS.items()}
UNKNOWN_YEAR = 2012


def long_form(change: str, i: int) -> str:
    """p.R1450* as p.Arg1450Ter; every third one also with a reference sequence and parentheses."""
    body = re.sub(r"[A-Z*]", lambda m: ONE_TO_THREE[m.group(0)], change[2:])
    return f"NP_000000.1:p.({body})" if i % 3 == 0 else f"p.{body}"


def _date(anchor: dt.date, days) -> str:
    if days is None or pd.isna(days):
        return ""
    return (anchor + dt.timedelta(days=math.floor(float(days)))).isoformat()


def build(t: StandardTables, folder: Path) -> dict:
    """Write the CDM tables to folder/cdm and a panels file to folder; return {GENIE test_id: OMOP test_id}."""
    cdm = folder / "cdm"
    cdm.mkdir(parents=True, exist_ok=True)
    p = t.patients.reset_index(drop=True)
    anchors, person, cond, meas, obs, op, death = {}, [], [], [], [], [], []
    sites = {c: i + 1 for i, c in enumerate(sorted(p["center"].dropna().unique()))}
    for i, r in p.iterrows():
        pid = 1000 + i
        year = UNKNOWN_YEAR if pd.isna(r["year_of_diagnosis"]) else int(r["year_of_diagnosis"])
        anchor = dt.date(year, 7, 1)
        anchors[r["patient_id"]] = (pid, anchor)
        person.append({"person_id": pid, "person_source_value": r["patient_id"],
                       "year_of_birth": "" if pd.isna(r["age_at_diagnosis"]) else year - int(r["age_at_diagnosis"]),
                       "gender_source_value": r["sex"], "race_source_value": r["race_ethnicity"],
                       "care_site_id": sites.get(r["center"], "")})
        source = f"{MORPHOLOGY.get(r['histology'], '')}-{r['primary_site_code']}"
        cond.append({"condition_occurrence_id": 3 * i + 1, "person_id": pid, "condition_concept_id": DIAGNOSIS,
                     "condition_start_date": anchor.isoformat(), "condition_source_value": source})
        if i % 50 == 0:      # decoys: the same cancer recorded again later, another condition earlier
            cond.append({"condition_occurrence_id": 3 * i + 2, "person_id": pid, "condition_concept_id": DIAGNOSIS,
                         "condition_start_date": _date(anchor, 90), "condition_source_value": "8000/3-C80.9"})
            cond.append({"condition_occurrence_id": 3 * i + 3, "person_id": pid, "condition_concept_id": OTHER_CONDITION,
                         "condition_start_date": _date(anchor, -400), "condition_source_value": "other"})
        meas.append({"measurement_id": 10 * i + 1, "person_id": pid, "measurement_concept_id": STAGE[r["stage_at_diagnosis"]],
                     "measurement_date": anchor.isoformat(), "value_as_concept_id": ""})
        if i % 40 == 0:      # decoy: a stage recorded long after diagnosis
            meas.append({"measurement_id": 10 * i + 2, "person_id": pid, "measurement_concept_id": STAGE["Stage IV"],
                         "measurement_date": _date(anchor, 400), "value_as_concept_id": ""})
        for j, (name, concept) in enumerate(FLAGS.items()):
            value = r[name]
            obs.append({"observation_id": 10 * i + j, "person_id": pid, "observation_concept_id": concept,
                        "observation_date": anchor.isoformat(),
                        "value_as_string": str(bool(value)).lower() if isinstance(value, (bool,)) or str(value) in ("True", "False")
                        else value})
        days = r["days_to_death_or_last_followup"]
        if pd.notna(days):
            op.append({"observation_period_id": pid, "person_id": pid,
                       "observation_period_start_date": _date(anchor, -365),
                       "observation_period_end_date": _date(anchor, days)})
            if r["vital_status"] == 1:
                death.append({"person_id": pid, "death_date": _date(anchor, days)})
    tests = t.tests
    ids, proc, variants = {}, [], []
    ngs = tests[tests["modality"] == "ngs"].sort_values(["patient_id", "test_order"], kind="mergesort")
    for k, (_, r) in enumerate(ngs.iterrows()):
        pid, anchor = anchors[r["patient_id"]]
        proc_id = 500000 + k
        ids[r["test_id"]] = f"procedure:{proc_id}"
        proc.append({"procedure_occurrence_id": proc_id, "person_id": pid, "procedure_concept_id": SEQUENCING,
                     "procedure_date": _date(anchor, r["days_to_entry_event"]), "procedure_source_value": r["panel_id"]})
    proc.append({"procedure_occurrence_id": 999999, "person_id": 1000, "procedure_concept_id": OTHER_PROCEDURE,
                 "procedure_date": "2012-01-01", "procedure_source_value": "not sequencing"})
    mutations = t.results[t.results["kind"] == "mutation"].reset_index(drop=True)
    for k, r in mutations.iterrows():
        variants.append({"variant_occurrence_id": k + 1, "procedure_occurrence_id": int(ids[r["test_id"]].split(":")[1]),
                         "target_gene1_symbol": r["gene"], "hgvs_p": long_form(r["variant"], k),
                         "sequence_alteration": r["variant_class"]})
    variants.append({"variant_occurrence_id": len(variants) + 1, "procedure_occurrence_id": 999999,
                     "target_gene1_symbol": "KRAS", "hgvs_p": "p.Gly12Asp", "sequence_alteration": "Missense_Mutation"})
    path_tests = tests[tests["modality"] == "pathology"].reset_index(drop=True)
    readings = t.results[t.results["kind"] == "msi_mmr"].set_index("test_id")["value"]
    for k, r in path_tests.iterrows():
        pid, anchor = anchors[r["patient_id"]]
        mid = 9000000 + k
        ids[r["test_id"]] = f"measurement:{mid}"
        meas.append({"measurement_id": mid, "person_id": pid,
                     "measurement_concept_id": MMR_TEST if ":MMR" in r["test_id"] else MSI_TEST,
                     "measurement_date": _date(anchor, r["days_to_entry_event"]),
                     "value_as_concept_id": READING[readings[r["test_id"]]]})
    drugs = []
    tr = t.treatments if t.treatments is not None else pd.DataFrame()
    for k, (_, r) in enumerate(tr.iterrows()):
        pid, anchor = anchors[r["patient_id"]]
        drugs.append({"drug_exposure_id": k + 1, "person_id": pid, "drug_concept_id": CANCER_DRUG,
                      "drug_exposure_start_date": _date(anchor, r["days_to_start"]),
                      "drug_source_value": f"drug {int(r['regimen_order'])}"})
    drugs.append({"drug_exposure_id": len(drugs) + 1, "person_id": 1000, "drug_concept_id": OTHER_DRUG,
                  "drug_exposure_start_date": "2012-07-02", "drug_source_value": "not a cancer drug"})
    care = [{"care_site_id": i, "care_site_name": c} for c, i in sites.items()]
    for name, rows in (("person", person), ("care_site", care), ("condition_occurrence", cond), ("measurement", meas),
                       ("observation", obs), ("observation_period", op), ("death", death),
                       ("procedure_occurrence", proc), ("variant_occurrence", variants), ("drug_exposure", drugs)):
        pd.DataFrame(rows).to_csv(cdm / f"{name}.csv", index=False)
    t.panels.to_csv(folder / "panels.csv", index=False)
    return ids


def mapping(**changes) -> dict:
    """The fixture's mapping, confirmed. changes replace top-level keys."""
    labels = lambda d: {v: k for k, v in d.items()}
    m = {
        "mapping_format": 1, "kind": "omop_cdm", "status": "confirmed", "confirmed_by": "A. Reviewer",
        "confirmed_on": "2026-11-17",
        "cdm": {"folder": "cdm", "date_format": "%Y-%m-%d"},
        "concept_sets": {
            "colorectal cancer": {"concepts": [DIAGNOSIS], "include_descendants": False},
            "stage": {"concepts": sorted(STAGE.values()), "include_descendants": False},
            "msi or mmr": {"concepts": [MSI_TEST, MMR_TEST], "include_descendants": False},
            "sequencing": {"concepts": [SEQUENCING], "include_descendants": False},
            "cancer drugs": {"concepts": [CANCER_DRUG], "include_descendants": False},
            **{name: {"concepts": [c], "include_descendants": False} for name, c in FLAGS.items()},
        },
        "patient_id": "person_source_value",
        "diagnosis": "colorectal cancer",
        "patients": {
            "age_at_diagnosis": {"computed": True},
            "sex": {"person": "gender_source_value"},
            "stage_at_diagnosis": {"measurement": "stage", "value": "measurement_concept_id", "labels": labels(STAGE)},
            "histology": {"index_condition": "condition_source_value", "extract": r"^(\d{4}/\d)-",
                          "labels": labels(MORPHOLOGY)},
            "primary_site_code": {"index_condition": "condition_source_value", "extract": r"-(C\d\d\.\d)$"},
            "year_of_diagnosis": {"computed": True},
            "center": {"care_site": "care_site_name"},
            **{name: {"observation": name, "value": "value_as_string"} for name in FLAGS if name != "tumor_location"},
        },
        "extras": {"tumor_location": {"observation": "tumor_location", "value": "value_as_string"},
                   "race_ethnicity": {"person": "race_source_value"}},
        "sequencing": {
            "procedures": "sequencing", "modality": "ngs", "entry_event": "report",
            "panel_id": {"procedure": "procedure_source_value"},
            "variants": {"table": "variant_occurrence", "test_link": "procedure_occurrence_id",
                         "gene": "target_gene1_symbol", "protein_change": "hgvs_p", "partner_gene": None,
                         "fusion_description": None, "variant_class": "sequence_alteration",
                         "variant_class_labels": None, "copy_number": None},
        },
        "biomarker_readings": [{"measurement": "msi or mmr", "kind": "msi_mmr", "gene": None,
                                "value": "value_as_concept_id", "labels": labels(READING),
                                "modality": "pathology", "entry_event": "report"}],
        "panels_file": "panels.csv",
        "panel_sizes_file": None,
        "treatments": {"drugs": "cancer drugs", "name": "drug_source_value"},
        "profile": {"entry_event": "report", "panel_coverage": "genes_only", "panel_sizes": "not_given",
                    "death_capture": "linked", "time_format": "real_dates", "privacy_min_count": 11,
                    "histology_vocabulary": "the source's categories, read from ICD-O-3 morphology codes",
                    "age_top_code": None, "known_quirks": [], "sources": {"entry_event": "made up"}},
    }
    m.update(changes)
    return m


def write_mapping(folder: Path, m: dict | None = None, name: str = "mapping.yaml") -> Path:
    path = folder / name
    path.write_bytes(yaml.safe_dump(m if m is not None else mapping(), sort_keys=False, allow_unicode=True).encode("utf-8"))
    return path


def omop_settings(folder: Path, resamples: int | None = None) -> Path:
    """The CRC settings read through the OMOP adapter: source.adapter omop, and no GENIE BPC choices."""
    from helpers import CRC_SETTINGS
    data = yaml.safe_load(CRC_SETTINGS.read_text(encoding="utf-8"))
    data["source"]["adapter"]["value"] = "omop"
    del data["source"]["genie_bpc"]
    if resamples is not None:
        data["uncertainty"]["bootstrap_resamples"]["value"] = resamples
    path = folder / "omop_settings.yaml"
    path.write_bytes(yaml.safe_dump(data, sort_keys=False, allow_unicode=True).encode("utf-8"))
    return path
