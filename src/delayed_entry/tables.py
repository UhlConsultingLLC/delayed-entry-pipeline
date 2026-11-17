"""The standard tables every adapter produces (design plan section 3).

Five tables, plus one optional extra:
    patients     one row per patient: outcome, cohort and adjustment fields
    tests        one row per test: its entry event and when it happened
    results      one row per result of a test: genotype calls and biomarker readings
    panels       one row per panel and gene: what a panel reads
    treatments   optional: one row per regimen, for the treatment clock
    panel_sizes  optional: each panel's coding size, for a mutation-burden fallback

All times are days from diagnosis, whatever the source gives (real dates, shifted dates or intervals).
A negative time comes before diagnosis.

Each field is marked, for each analysis component, as required, optional or needed for one check only,
so the readiness screen can say exactly which analyses a missing field blocks. The cohort rules and the
adjustment variables name further fields in the settings; the screen adds those as required.

The outcome fields (vital status and days to death or last follow-up) are marked as outcomes. The
screen never sees them: it is given `outcome_blind(tables)`, where they are replaced by one column saying
only whether each patient's outcome is recorded.
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields as dc_fields

import pandas as pd

TEXT, NUMBER, INTEGER, BOOLEAN = "text", "number", "integer", "boolean"
REQUIRED, OPTIONAL, CHECK_ONLY = "required", "optional", "check-only"

# Analysis components: the parts of the analysis a field can be needed for.
COMPONENTS = {
    "cohort": "the cohort rules and the flow diagram",
    "survival": "survival on either clock",
    "corrected_clock": "the corrected clock's entry time",
    "genotype": "the genotype groups",
    "wild_type_confirmation": "confirming wild type from panel coverage",
    "adjustment": "the adjusted and weighted models",
    "treatment_clock": "the treatment-clock analysis",
    "biomarker": "the biomarker analyses",
    "time_at_risk_check": "the check for a test that comes back after follow-up ended",
    "flags": "the flags counted in the cohort",
}


@dataclass(frozen=True)
class Field:
    table: str
    name: str
    kind: str
    description: str
    used_by: dict = field(default_factory=dict)     # component: REQUIRED / OPTIONAL / CHECK_ONLY
    outcome: bool = False
    key: bool = False                               # an identity column: the table can't be used without it


def _f(table, name, kind, description, outcome=False, key=False, **used_by):
    return Field(table, name, kind, description, dict(used_by), outcome, key)


FIELDS: list[Field] = [
    # ---- patients ----
    _f("patients", "patient_id", TEXT, "the patient's ID in the source", key=True, cohort=REQUIRED, survival=REQUIRED),
    _f("patients", "vital_status", INTEGER, "1 = died, 0 = alive at last follow-up", outcome=True,
       survival=REQUIRED),
    _f("patients", "days_to_death_or_last_followup", NUMBER, "days from diagnosis to death or last follow-up",
       outcome=True, survival=REQUIRED),
    _f("patients", "age_at_diagnosis", NUMBER, "age in years at diagnosis", cohort=OPTIONAL, adjustment=OPTIONAL),
    _f("patients", "sex", TEXT, "sex as the source records it", adjustment=OPTIONAL),
    _f("patients", "stage_at_diagnosis", TEXT, "stage at diagnosis, in the source's labels", cohort=OPTIONAL),
    _f("patients", "histology", TEXT, "histology: an ICD-O-3 code or the source's category (the profile says which)",
       cohort=OPTIONAL),
    _f("patients", "primary_site_code", TEXT, "primary site, ICD-O-3 topography code", cohort=OPTIONAL),
    _f("patients", "year_of_diagnosis", INTEGER, "calendar year of diagnosis (exact or estimated)", adjustment=OPTIONAL),
    _f("patients", "center", TEXT, "the center or institution", adjustment=OPTIONAL),
    _f("patients", "censored_at_age_limit", BOOLEAN,
       "survival withheld by the source's age rule, so follow-up ends at the age limit (a source quirk)",
       flags=CHECK_ONLY),
    _f("patients", "multiple_cohort_cancers", BOOLEAN, "more than one cancer of the study's kind in the source",
       biomarker=CHECK_ONLY),
    _f("patients", "other_primary_same_cancer", BOOLEAN,
       "another primary of the same cancer type, whose reports could be mistaken for this one's", biomarker=CHECK_ONLY),
    # ---- tests ----
    _f("tests", "patient_id", TEXT, "the patient the test belongs to", key=True, corrected_clock=REQUIRED,
       genotype=REQUIRED, biomarker=REQUIRED),
    _f("tests", "test_id", TEXT, "the test's ID; results link to it", key=True, corrected_clock=REQUIRED,
       genotype=REQUIRED, biomarker=REQUIRED),
    _f("tests", "modality", TEXT, "what kind of test: for example ngs (sequencing) or pathology",
       corrected_clock=REQUIRED, biomarker=REQUIRED),
    _f("tests", "entry_event", TEXT, "which event the time is: report, sample_collection or order",
       corrected_clock=REQUIRED),
    _f("tests", "days_to_entry_event", NUMBER, "days from diagnosis to the event", corrected_clock=REQUIRED,
       biomarker=OPTIONAL),
    _f("tests", "days_from_birth", NUMBER,
       "days from birth to the event, if the source dates events from birth; when given, a result is judged to come "
       "before the entry test on this clock, as the source's own dates are", biomarker=OPTIONAL),
    _f("tests", "panel_id", TEXT, "the panel the test used", wild_type_confirmation=REQUIRED, corrected_clock=OPTIONAL),
    _f("tests", "test_order", INTEGER, "the source's own sequence number for the patient's tests (a tie-break)",
       corrected_clock=OPTIONAL),
    _f("tests", "after_follow_up_end", BOOLEAN,
       "the event came on or after the end of the patient's follow-up, as the source records it",
       time_at_risk_check=CHECK_ONLY),
    _f("tests", "linked_to_other_cancers", BOOLEAN, "the test is linked to more than one of the patient's cancers",
       flags=CHECK_ONLY),
    # ---- results ----
    _f("results", "test_id", TEXT, "the test the result came from", key=True, genotype=REQUIRED, biomarker=REQUIRED),
    _f("results", "kind", TEXT, "mutation, fusion, copy_number, protein or msi_mmr", key=True, genotype=REQUIRED,
       biomarker=REQUIRED),
    _f("results", "gene", TEXT, "the gene (HGNC symbol)", genotype=REQUIRED),
    _f("results", "variant", TEXT, "the variant: for a mutation, the protein change in HGVSp short form, for example "
       "p.G12D; for a fusion, the source's description of it", genotype=REQUIRED),
    _f("results", "variant_class", TEXT, "the variant class, for example Missense_Mutation; a biomarker's "
       "mutation-burden fallback counts mutations by it", genotype=OPTIONAL, biomarker=CHECK_ONLY),
    _f("results", "value", TEXT, "a reading, for example a microsatellite status", biomarker=REQUIRED),
    # ---- panels ----
    _f("panels", "panel_id", TEXT, "the panel", key=True, wild_type_confirmation=REQUIRED),
    _f("panels", "gene", TEXT, "a gene the panel reads", key=True, wild_type_confirmation=REQUIRED),
    _f("panels", "codon_from", INTEGER, "first codon of a region read (blank: the whole gene is listed only)",
       wild_type_confirmation=OPTIONAL),
    _f("panels", "codon_to", INTEGER, "last codon of that region", wild_type_confirmation=OPTIONAL),
    # ---- panel_sizes (optional table) ----
    _f("panel_sizes", "panel_id", TEXT, "the panel", key=True, biomarker=CHECK_ONLY),
    _f("panel_sizes", "coding_megabases", NUMBER, "the panel's coding size in megabases", biomarker=CHECK_ONLY),
    # ---- treatments (optional table) ----
    _f("treatments", "patient_id", TEXT, "the patient", key=True, treatment_clock=REQUIRED),
    _f("treatments", "regimen", TEXT, "the regimen's drugs, as text", treatment_clock=OPTIONAL),
    _f("treatments", "regimen_order", INTEGER, "the source's sequence number for the patient's regimens",
       treatment_clock=OPTIONAL),
    _f("treatments", "line", INTEGER, "line of therapy, if the source gives it", treatment_clock=OPTIONAL),
    _f("treatments", "days_to_start", NUMBER, "days from diagnosis to the regimen's start", treatment_clock=REQUIRED),
]

TABLES = ("patients", "tests", "results", "panels", "panel_sizes", "treatments")
OPTIONAL_TABLES = ("panel_sizes", "treatments")
ENTRY_EVENTS = ("report", "sample_collection", "order")
RESULT_KINDS = ("mutation", "fusion", "copy_number", "protein", "msi_mmr")
OUTCOME_FIELDS = tuple(f.name for f in FIELDS if f.outcome)
OUTCOME_RECORDED = "outcome_recorded"     # what the screen sees in place of the outcome fields


def field_spec(table: str, name: str) -> Field | None:
    for f in FIELDS:
        if f.table == table and f.name == name:
            return f
    return None


def table_fields(table: str) -> list[Field]:
    return [f for f in FIELDS if f.table == table]


@dataclass
class StandardTables:
    """What an adapter produces. Columns not in FIELDS are kept as extras (for example a covariate
    only one cancer has, such as tumor location); the settings may name them."""
    patients: pd.DataFrame
    tests: pd.DataFrame
    results: pd.DataFrame
    panels: pd.DataFrame
    treatments: pd.DataFrame | None = None
    panel_sizes: pd.DataFrame | None = None

    def table(self, name: str) -> pd.DataFrame | None:
        return getattr(self, name)

    def columns(self) -> dict[str, list[str]]:
        return {t: ([] if self.table(t) is None else list(self.table(t).columns)) for t in TABLES}


@dataclass(frozen=True)
class OutcomeBlindTables:
    """The standard tables with every outcome value removed. The patients table has, in place of the
    outcome fields, `outcome_recorded`: True if both outcome fields hold a value, False if either is blank,
    and blank if the source has no outcome columns at all. This is all the readiness screen is given."""
    patients: pd.DataFrame
    tests: pd.DataFrame
    results: pd.DataFrame
    panels: pd.DataFrame
    treatments: pd.DataFrame | None
    panel_sizes: pd.DataFrame | None
    outcome_fields_present: tuple[str, ...]

    def table(self, name: str) -> pd.DataFrame | None:
        return getattr(self, name)


def _drop_outcomes(df: pd.DataFrame | None) -> pd.DataFrame | None:
    if df is None:
        return None
    return df.drop(columns=[c for c in OUTCOME_FIELDS if c in df.columns]).copy()


def outcome_blind(t: StandardTables) -> OutcomeBlindTables:
    """The tables the screen may see. Outcome values never pass through: only whether they are recorded."""
    p = t.patients
    present = tuple(c for c in OUTCOME_FIELDS if c in p.columns)
    blind = _drop_outcomes(p)
    if len(present) == len(OUTCOME_FIELDS):
        recorded = p[list(OUTCOME_FIELDS)].notna().all(axis=1)
        blind[OUTCOME_RECORDED] = recorded.astype("boolean").to_numpy()
    else:
        blind[OUTCOME_RECORDED] = pd.array([pd.NA] * len(p), dtype="boolean")
    return OutcomeBlindTables(blind, _drop_outcomes(t.tests), _drop_outcomes(t.results), _drop_outcomes(t.panels),
                              _drop_outcomes(t.treatments), _drop_outcomes(t.panel_sizes), present)


# ---- Types ----
def _convert(series: pd.Series, kind: str) -> pd.Series:
    """The column converted to its standard type. A value that can't convert becomes blank (counted by bad_values)."""
    if kind == TEXT:
        out = series.astype("string")
        return out.mask(out.str.len() == 0)
    if kind in (NUMBER, INTEGER):
        if pd.api.types.is_bool_dtype(series):
            return pd.Series(pd.NA, index=series.index, dtype="Float64" if kind == NUMBER else "Int64")
        num = pd.to_numeric(series, errors="coerce")
        if kind == NUMBER:
            return num.astype("Float64")
        whole = num.where(num.isna() | (num == num.round()))
        return whole.astype("Int64")
    if kind == BOOLEAN:
        if pd.api.types.is_bool_dtype(series):
            return series.astype("boolean")
        text = series.astype("string").str.strip().str.lower()
        mapped = text.map({"true": True, "1": True, "1.0": True, "yes": True,
                           "false": False, "0": False, "0.0": False, "no": False})
        return mapped.astype("boolean")
    raise ValueError(f"unknown kind {kind}")


def bad_values(series: pd.Series, kind: str) -> int:
    """How many filled values don't convert to the field's type."""
    filled = series.notna() & ~(series.astype("string").str.strip().fillna("") == "")
    return int((filled & _convert(series, kind).isna()).sum())


def conform(t: StandardTables) -> tuple[StandardTables, dict[str, int]]:
    """Every standard column converted to its type, and how many values in each failed to convert
    ({"table.field": n}, only fields with failures). Extras are left as they are."""
    failures, out = {}, {}
    for name in TABLES:
        df = t.table(name)
        if df is None:
            out[name] = None
            continue
        df = df.copy()
        for f in table_fields(name):
            if f.name in df.columns:
                n = bad_values(df[f.name], f.kind)
                if n:
                    failures[f"{name}.{f.name}"] = n
                df[f.name] = _convert(df[f.name], f.kind)
        out[name] = df.reset_index(drop=True)
    return StandardTables(**out), failures


def empty_table(name: str) -> pd.DataFrame:
    """A table with the standard columns and no rows."""
    return pd.DataFrame({f.name: pd.Series(dtype="object") for f in table_fields(name)})


def table_names() -> list[str]:
    return [f.name for f in dc_fields(StandardTables)]
