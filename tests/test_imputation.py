"""Multiple imputation: when the complete-case rule calls for it, the settings must say how, the screen says so
before any outcome is read, and the run carries it into the spec R reads. Made-up data only (the R part is tested in
tests/test_run.py and tests/R/test_r_units.R)."""
import pytest
import yaml

from delayed_entry.analysis import RunStop, build_analysis
from delayed_entry.run import build_spec
from delayed_entry.screen import BLOCK, DEVIATION, NOTE, screen
from delayed_entry.settings import SCHEMA_BY_PATH, load_settings
from delayed_entry.tables import outcome_blind
from helpers import CRC_SETTINGS, bpc_profile, crc_settings, make_tables

GROUPS = {"RAS-mutant": 40, "BRAF V600E": 35, "Wild type": 40}
MODEL = SCHEMA_BY_PATH["missing_data.imputation_model"].example


def tables(missing_sex: int):
    """The made-up tables, years spread over the three bands, some patients' sex unknown, and enough early deaths
    that entry stays at the report."""
    t = make_tables(GROUPS)
    t.patients["year_of_diagnosis"] = [(2010, 2014, 2017)[i % 3] for i in range(len(t.patients))]
    t.patients.loc[t.patients.index[:missing_sex], "sex"] = None
    early = t.patients["patient_id"].str.startswith("RAS-") & (t.patients["patient_id"] < "RAS-0020")
    t.patients.loc[early, "days_to_death_or_last_followup"] = 500.0
    return t


def with_model(tmp_path, **changes):
    """The CRC settings with an imputation model (the template's example) and any other values changed."""
    data = yaml.safe_load(CRC_SETTINGS.read_text(encoding="utf-8"))
    for path, value in changes.items():
        node = data
        for part in path.split("."):
            node = node[part]
        node["value"] = value
    data["missing_data"]["imputation_model"] = {"value": MODEL, "made_by": "a test", "source": "made up"}
    p = tmp_path / "with_model.yaml"
    p.write_bytes(yaml.safe_dump(data, sort_keys=False).encode("utf-8"))
    return load_settings(p)


def items(rep, topic):
    return [i for i in rep.items if i.topic == topic]


# ---- The screen: before any outcome is read ----
def test_too_much_missing_with_no_imputation_model_blocks_the_plan():
    rep = screen(outcome_blind(tables(10)), crc_settings(), bpc_profile())
    (item,) = items(rep, "missing data")
    assert item.severity == BLOCK and item.analysis == "primary" and "missing_data.imputation_model" in item.effect
    assert "10" not in item.finding                               # the count itself is never shown


def test_with_an_imputation_model_the_plan_states_it_as_a_note(tmp_path):
    """Imputation the settings set in advance is the plan, not a departure from it: a note, so the primary
    analysis can still run as planned."""
    rep = screen(outcome_blind(tables(10)), with_model(tmp_path), bpc_profile())
    (item,) = items(rep, "missing data")
    assert item.severity == NOTE and "multiple imputation (20 imputations)" in item.effect
    assert any(i.finding == "the primary analysis can run as planned" for i in items(rep, "runs as standard"))
    assert not [i for i in rep.items if i.analysis == "primary" and i.severity == DEVIATION]


def test_a_complete_case_study_gets_no_missing_data_item():
    assert items(screen(outcome_blind(tables(3)), crc_settings(), bpc_profile()), "missing data") == []


# ---- The run's analysis cohort ----
def test_the_run_stops_when_imputation_is_needed_and_the_settings_dont_say_how():
    with pytest.raises(RunStop, match="missing_data.imputation_model"):
        build_analysis(tables(10), crc_settings(), bpc_profile())


def test_with_an_imputation_model_the_run_goes_on_and_the_spec_carries_it(tmp_path):
    s = with_model(tmp_path, **{"missing_data.imputations": 7})
    data = build_analysis(tables(10), s, bpc_profile())
    assert data.covariate_plan["missing_data_rule"] == "multiple imputation"
    assert data.covariates["sex"].isna().sum() == 10            # left missing for R to impute
    spec = build_spec(data, s)
    assert spec["missing_data"]["rule"] == "multiple imputation" and spec["missing_data"]["imputations"] == 7
    assert spec["missing_data"]["model"]["outcome_terms"] == ["event", "cumulative_hazard_corrected", "entry_time"]


def test_a_complete_case_spec_says_so_and_holds_no_model():
    data = build_analysis(tables(3), crc_settings(), bpc_profile())
    spec = build_spec(data, crc_settings())
    assert spec["missing_data"] == {"rule": "complete case", "imputations": 20, "model": None}


def test_year_bands_are_counted_where_the_year_is_known_when_imputing(tmp_path):
    s = with_model(tmp_path)
    t = tables(10)
    t.patients.loc[t.patients.index[20:30], "year_of_diagnosis"] = None
    data = build_analysis(t, s, bpc_profile())
    assert set(data.covariate_plan["year_bands"].values()) == {"2012 or earlier", "2013-2015", "2016-2018"}
