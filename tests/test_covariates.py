"""The adjustment variables: merging levels under new names, keeping levels, year bands and the missing-data rule.
Made-up values only."""
import pandas as pd
import pytest

from delayed_entry.covariates import (band_of, code_adjustment, merge_bands, merged_label, missing_data_rule,
                                      model_bands, apply_model_bands)
from delayed_entry.settings import SettingsError
from helpers import crc_settings, settings_with

SMOKING = {"field": "smoking", "form": "categorical", "reference": "Never",
           "merge": {"Never": ["never"], "Former": ["former, over a year", "former, unknown when"],
                     "Current": ["current", "former, within a year"]}}
CRC_VARIABLES = crc_settings().read("adjustment.variables")


def with_variables(tmp_path, variables):
    return settings_with(tmp_path, {"adjustment.variables": variables})


# ---- The settings: what a merge may say ----
def test_a_merge_loads(tmp_path):
    s = with_variables(tmp_path, CRC_VARIABLES + [SMOKING])
    assert s.read("adjustment.variables")[-1]["merge"]["Current"] == ["current", "former, within a year"]


@pytest.mark.parametrize("change, problem", [
    ({"merge": {"A": ["x"], "B": ["x"]}}, "merged into both 'A' and 'B'"),
    ({"keep": ["Never", "Sometimes"]}, "'Sometimes', which is not one of the merged levels"),
    ({"reference": "Unknown"}, "the reference 'Unknown' is not one of its levels"),
    ({"form": "linear"}, "merge is for form categorical"),
])
def test_a_merge_that_cant_work_names_its_problem(tmp_path, change, problem):
    with pytest.raises(SettingsError, match=problem.replace("(", r"\(").replace(")", r"\)")):
        with_variables(tmp_path, CRC_VARIABLES + [{**SMOKING, **change}])


def test_a_spline_needs_three_rising_knots(tmp_path):
    bad = [{"field": "age_at_diagnosis", "form": "spline", "knots_at_percentiles": [50, 10, 90]}] + CRC_VARIABLES[1:]
    with pytest.raises(SettingsError, match="must rise"):
        with_variables(tmp_path, bad)
    bad[0]["knots_at_percentiles"] = [10, 90]
    with pytest.raises(SettingsError, match="at least three"):
        with_variables(tmp_path, bad)


def test_a_field_listed_twice_is_refused(tmp_path):
    with pytest.raises(SettingsError, match="listed twice"):
        with_variables(tmp_path, CRC_VARIABLES + [CRC_VARIABLES[1]])


# ---- Coding ----
def patients(**columns):
    n = len(next(iter(columns.values())))
    base = {"age_at_diagnosis": [60.0] * n, "sex": ["Female"] * n, "center": ["MSK"] * n,
            "tumor_location": ["Left colon"] * n, "year_of_diagnosis": [2015] * n}
    return pd.DataFrame({**base, **columns}, index=[f"P{i}" for i in range(n)])


def test_merged_levels_take_their_new_names_and_others_are_missing(tmp_path):
    s = with_variables(tmp_path, CRC_VARIABLES + [SMOKING])
    p = patients(smoking=["never", " current ", "former, within a year", "former, unknown when", "unknown", None])
    coded = code_adjustment(p, s)
    assert list(coded.frame["smoking"].astype(object).where(coded.frame["smoking"].notna(), None)) == \
        ["Never", "Current", "Current", "Former", None, None]
    term = next(t for t in coded.terms if t.field == "smoking")
    assert term.levels == ("Never", "Former", "Current")        # the reference first, then the settings' order


def test_kept_levels_and_the_reference_come_from_the_settings():
    coded = code_adjustment(patients(center=["DFCI", "MSK", "VICC", "UHN"]), crc_settings())
    assert list(coded.frame["center"].astype(object).where(coded.frame["center"].notna(), None)) == \
        ["DFCI", "MSK", "VICC", None]
    assert next(t for t in coded.terms if t.field == "center").levels == ("MSK", "DFCI", "VICC")
    assert next(t for t in coded.terms if t.field == "sex").levels == ("Male", "Female")


def test_years_fall_in_the_settings_bands():
    bands = crc_settings().read("adjustment.year_bands")
    got = band_of(pd.Series([2005, 2012, 2013, 2018, 2019, None]), bands)
    assert list(got.astype(object).where(got.notna(), None)) == \
        ["2012 or earlier", "2012 or earlier", "2013-2015", "2016-2018", None, None]


# ---- Merging small year bands (the CRC BPC plan's rule, made general) ----
BANDS = [{"label": "early", "to": 2012}, {"label": "middle", "from": 2013, "to": 2015}, {"label": "late", "from": 2016}]


@pytest.mark.parametrize("counts, expected", [
    ({"early": 40, "middle": 50, "late": 60}, {"early": "early", "middle": "middle", "late": "late"}),
    ({"early": 12, "middle": 50, "late": 60}, {"early": "2015 or earlier", "middle": "2015 or earlier", "late": "late"}),
    ({"early": 40, "middle": 50, "late": 7}, {"early": "early", "middle": "2013 or later", "late": "2013 or later"}),
    ({"early": 40, "middle": 9, "late": 60}, {"early": "2015 or earlier", "middle": "2015 or earlier", "late": "late"}),
    ({"early": 70, "middle": 9, "late": 40}, {"early": "early", "middle": "2013 or later", "late": "2013 or later"}),
    ({"early": 12, "middle": 25, "late": 60}, {"early": "2015 or earlier", "middle": "2015 or earlier", "late": "late"}),
    ({"early": 40, "middle": 9, "late": 40}, {"early": "2015 or earlier", "middle": "2015 or earlier", "late": "late"}),
])
def test_small_bands_merge_by_the_plans_rule(counts, expected):
    assert merge_bands(counts, BANDS, 30) == expected


def test_a_band_still_small_after_its_end_joins_merges_with_its_other_neighbor():
    four = [{"label": "a", "to": 2010}, {"label": "b", "from": 2011, "to": 2012}, {"label": "c", "from": 2013, "to": 2014},
            {"label": "d", "from": 2015}]
    got = merge_bands({"a": 5, "b": 20, "c": 50, "d": 60}, four, 30)
    assert got == {"a": "2014 or earlier", "b": "2014 or earlier", "c": "2014 or earlier", "d": "d"}


def test_bands_that_cant_be_merged_into_two_stop():
    with pytest.raises(ValueError, match="every year band is small"):
        merge_bands({"early": 5, "middle": 5, "late": 5}, BANDS, 30)


def test_four_bands_merge_the_smallest_neighbor():
    four = [{"label": "a", "to": 2010}, {"label": "b", "from": 2011, "to": 2012}, {"label": "c", "from": 2013, "to": 2014},
            {"label": "d", "from": 2015, "to": 2016}]
    got = merge_bands({"a": 50, "b": 10, "c": 40, "d": 60}, four, 30)
    assert got == {"a": "a", "b": "2011-2014", "c": "2011-2014", "d": "d"}


def test_merged_labels():
    assert merged_label([{"label": "x", "to": 2012}, {"label": "y", "from": 2013, "to": 2015}]) == "2015 or earlier"
    assert merged_label([{"label": "x", "from": 2013, "to": 2015}, {"label": "y", "from": 2016}]) == "2013 or later"
    assert merged_label([{"label": "x", "from": 2013, "to": 2015}, {"label": "y", "from": 2016, "to": 2018}]) == "2013-2018"


def test_the_model_bands_come_from_the_patients_in_the_model():
    s = crc_settings()
    years = [2010] * 5 + [2014] * 40 + [2017] * 40
    coded = code_adjustment(patients(year_of_diagnosis=years), s)
    mapping = model_bands(coded, s, coded.frame.index)
    assert mapping["2012 or earlier"] == mapping["2013-2015"] == "2015 or earlier"
    banded = apply_model_bands(coded, mapping)
    assert next(t for t in banded.terms if t.form == "bands").levels == ("2015 or earlier", "2016-2018")


# ---- The missing-data rule, with the settings' share ----
def test_complete_case_up_to_the_settings_share(tmp_path):
    p = patients(sex=["Female"] * 95 + [None] * 5)
    assert missing_data_rule(code_adjustment(p, crc_settings()), crc_settings())[0] == "complete case"
    stricter = settings_with(tmp_path, {"missing_data.complete_case_max_share": 0.04})
    assert missing_data_rule(code_adjustment(p, stricter), stricter)[0] == "multiple imputation"
