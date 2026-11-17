"""The standard tables (design plan section 3)."""
import pandas as pd

from delayed_entry import tables as tb
from delayed_entry.tables import (CHECK_ONLY, COMPONENTS, FIELDS, OPTIONAL, REQUIRED, StandardTables, conform,
                                  outcome_blind)
from helpers import make_tables


def test_the_five_tables_and_their_required_fields():
    names = {f.table for f in FIELDS}
    assert names == {"patients", "tests", "results", "panels", "treatments", "panel_sizes"}
    required = {(f.table, f.name) for f in FIELDS if REQUIRED in f.used_by.values()}
    assert {("patients", "patient_id"), ("patients", "vital_status"), ("patients", "days_to_death_or_last_followup"),
            ("tests", "patient_id"), ("tests", "test_id"), ("tests", "days_to_entry_event"), ("tests", "entry_event"),
            ("results", "test_id")} <= required


def test_every_field_is_marked_for_at_least_one_component():
    for f in FIELDS:
        assert f.used_by, f.name
        assert set(f.used_by) <= set(COMPONENTS), f.name
        assert set(f.used_by.values()) <= {REQUIRED, OPTIONAL, CHECK_ONLY}, f.name


def test_the_outcome_fields_are_the_two_survival_fields():
    assert tb.OUTCOME_FIELDS == ("vital_status", "days_to_death_or_last_followup")


def test_outcome_blind_tables_say_only_whether_an_outcome_is_recorded():
    t = make_tables({"Wild type": 3})
    t.patients.loc[1, "vital_status"] = pd.NA
    blind = outcome_blind(t)
    assert "vital_status" not in blind.patients and "days_to_death_or_last_followup" not in blind.patients
    assert blind.patients["outcome_recorded"].tolist() == [True, False, True]
    assert blind.outcome_fields_present == tb.OUTCOME_FIELDS


def test_without_outcome_columns_the_record_is_blank():
    t = make_tables({"Wild type": 2})
    t = StandardTables(t.patients.drop(columns=list(tb.OUTCOME_FIELDS)), t.tests, t.results, t.panels)
    blind = outcome_blind(t)
    assert blind.patients["outcome_recorded"].isna().all() and blind.outcome_fields_present == ()


def test_conform_converts_types_and_counts_failures():
    p = pd.DataFrame({"patient_id": ["a", "b", "c"], "age_at_diagnosis": ["61", "x", ""],
                      "vital_status": ["1", "0", "1.5"], "censored_at_age_limit": ["Yes", "no", None],
                      "tumor_location": ["Left", "Right", None]})
    t, failures = conform(StandardTables(p, tb.empty_table("tests"), tb.empty_table("results"), tb.empty_table("panels")))
    assert failures == {"patients.age_at_diagnosis": 1, "patients.vital_status": 1}
    assert t.patients["age_at_diagnosis"].tolist()[0] == 61 and pd.isna(t.patients["age_at_diagnosis"][1])
    assert t.patients["censored_at_age_limit"].tolist()[:2] == [True, False]
    assert t.patients["tumor_location"].tolist()[:2] == ["Left", "Right"]       # extras kept as they are
