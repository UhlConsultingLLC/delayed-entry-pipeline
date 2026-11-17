"""The run's analysis cohort: the checks that need outcomes, who each analysis keeps, and the spec R reads.
Made-up tables only (helpers.make_tables): every patient stage IV, one report on day 100, died on day 900."""
import dataclasses

import pandas as pd
import pytest

from delayed_entry.analysis import (APPLIED, PASS, WARNING, RunStop, build_analysis, compare_with_plan, delays,
                                    keep_plan_decisions, treatment_clock)
from delayed_entry.run import build_spec, labels
from delayed_entry.screen import BURDEN_HIGH, BURDEN_LOW
from delayed_entry.tables import conform
from helpers import bpc_profile, crc_settings, make_tables, settings_with

GROUPS = {"RAS-mutant": 40, "BRAF V600E": 35, "Wild type": 40}


def tables(groups=GROUPS, **kw):
    """The made-up tables, with years of diagnosis spread over the settings' three bands."""
    t = make_tables(groups, **kw)
    t.patients["year_of_diagnosis"] = [(2010, 2014, 2017)[i % 3] for i in range(len(t.patients))]
    return t


def test_bands_too_small_to_model_stop_the_run():
    with pytest.raises(RunStop, match="every year band is small"):
        build_analysis(early_deaths(make_tables(GROUPS)), crc_settings(), bpc_profile())   # every year in one band


def set_outcome(t, ids, *, exit_day=None, status=None):
    p = t.patients
    rows = p["patient_id"].isin(ids)
    if exit_day is not None:
        p.loc[rows, "days_to_death_or_last_followup"] = exit_day
    if status is not None:
        p.loc[rows, "vital_status"] = status
    return t


def ids(group, n):
    return [f"{group[:3]}-{i:04d}" for i in range(n)]


def check(data, name):
    return next(c for c in data.checks if c.check == name)


def early_deaths(t, n=20):
    """Give n patients a death 400 days after their report: enough for the survival-after-entry check to pass."""
    return set_outcome(t, ids("Wil", 0) + [f"RAS-{i:04d}" for i in range(n)], exit_day=500.0)


# ---- Time at risk from the times themselves ----
def test_follow_up_ending_by_entry_excludes_at_the_rule_that_holds_the_time_check():
    t = early_deaths(tables())
    set_outcome(t, ["Wil-0003"], exit_day=80.0, status=0)          # censored before the report came back
    data = build_analysis(t, crc_settings(), bpc_profile())
    c = data.blind_cohort
    assert c.patients.loc["Wil-0003", "excluded_at"] == "entry_before_exit"
    assert "Wil-0003" not in set(data.cohort["patient_id"])
    assert check(data, "time at risk from the times themselves").verdict == APPLIED
    step = data.flow[data.flow["step"] == "entry_before_exit"].iloc[0]
    assert step["excluded"] == 1                                    # made-up data: an exact count is fine here


def test_nothing_to_exclude_by_the_times_passes():
    data = build_analysis(early_deaths(tables()), crc_settings(), bpc_profile())
    assert check(data, "time at risk from the times themselves").verdict == PASS
    assert data.cohort["in_analysis"].all()


# ---- Survival after the entry event ----
def test_enough_early_deaths_keep_entry_at_the_event():
    data = build_analysis(early_deaths(tables()), crc_settings(), bpc_profile())
    assert check(data, "survival after the entry event").verdict == PASS
    assert (data.cohort["entry_day"] == 100.0).all() and data.d2["entry_shift_days"] == 0


def test_few_early_deaths_move_entry_by_the_settings_window():
    t = tables()                                                    # every death 800 days after the report
    set_outcome(t, ids("Wild type", 3), exit_day=820.0, status=0)   # censored: no time left once entry moves
    data = build_analysis(t, crc_settings(), bpc_profile())
    assert check(data, "survival after the entry event").verdict == APPLIED
    assert (data.cohort["entry_day"] == 100.0 + 730).all()
    left = data.cohort[~data.cohort["in_analysis"]]
    assert sorted(left["patient_id"]) == ids("Wild type", 3)
    assert len(data.cohort) == sum(GROUPS.values())                 # still eligible, out of the analysis


def test_the_window_and_the_action_come_from_the_settings(tmp_path):
    s = settings_with(tmp_path, {"survival_after_entry.if_fewer": "warn_only"})
    data = build_analysis(tables(), s, bpc_profile())
    assert check(data, "survival after the entry event").verdict == WARNING
    assert (data.cohort["entry_day"] == 100.0).all()
    s = settings_with(tmp_path, {"survival_after_entry.window_days": 900}, "wider.yaml")
    data = build_analysis(early_deaths(tables()), s, bpc_profile())
    assert data.d2["early_deaths"] == len(data.cohort)              # every death within 900 days of the report


# ---- The smallest-group rule, again ----
def test_a_group_that_drops_below_the_minimum_is_described_and_the_change_is_stated():
    t = early_deaths(tables({"RAS-mutant": 40, "BRAF V600E": 31, "Wild type": 40}))
    set_outcome(t, ids("BRAF V600E", 2), exit_day=50.0, status=0)
    data = build_analysis(t, crc_settings(), bpc_profile())
    d1 = data.d1.set_index(["analysis", "comparison"])["decision"]
    assert d1[("primary", "BRAF V600E vs Wild type")] == "descriptive only"
    assert data.analysis("primary").described == {"BRAF V600E": "smallest group"}
    plan = {"primary": {"RAS-mutant vs Wild type": "compare", "BRAF V600E vs Wild type": "compare"}}
    out = compare_with_plan(data, plan)
    assert out[0].verdict == APPLIED and "BRAF V600E vs Wild type is now descriptive only" in out[0].detail


def test_a_group_the_settings_block_stops_the_run():
    t = early_deaths(tables({"RAS-mutant": 31, "BRAF V600E": 35, "Wild type": 40}))
    set_outcome(t, ids("RAS-mutant", 40)[25:27], exit_day=50.0, status=0)
    with pytest.raises(RunStop, match="blocks"):
        build_analysis(t, crc_settings(), bpc_profile())


def test_a_comparison_the_plan_describes_stays_described():
    data = build_analysis(early_deaths(tables()), crc_settings(), bpc_profile())
    plan = {"primary": {"RAS-mutant vs Wild type": "compare", "BRAF V600E vs Wild type": "descriptive only"}}
    keep_plan_decisions(data, plan)
    d1 = data.d1.set_index(["analysis", "comparison"])["decision"]
    assert d1[("primary", "BRAF V600E vs Wild type")] == "descriptive only"
    assert data.analysis("primary").described == {"BRAF V600E": "smallest group"}


def test_an_analysis_the_plan_doesnt_run_must_not_run():
    data = build_analysis(early_deaths(tables()), crc_settings(), bpc_profile())
    with pytest.raises(RunStop, match="doesn't run it"):
        compare_with_plan(data, {"SA1": {"all comparisons": "not run"}})


# ---- Who each analysis keeps ----
def test_short_waits_are_kept_by_the_settings_cutoff(tmp_path):
    data = build_analysis(early_deaths(tables()), crc_settings(), bpc_profile())
    assert data.secondary["keep_SA1"].all()                         # day 100 is within 6 months
    s = settings_with(tmp_path, {"analyses": [{**a, "keep": {"within_months": 3}} if a["id"] == "SA1" else a
                                              for a in crc_settings().read("analyses")]})
    data = build_analysis(early_deaths(tables()), s, bpc_profile())
    assert not data.secondary["keep_SA1"].any()


def test_the_treatment_clock_counts_from_the_first_line():
    cohort = pd.DataFrame({"patient_id": ["A", "B", "C"], "genotype": ["Wild type"] * 3, "entry_day": [100.0, 10.0, 100.0],
                           "exit_day": [900.0, 900.0, 900.0], "event": [1.0, 0.0, 1.0], "in_analysis": [True, True, True]})
    tr = pd.DataFrame({"patient_id": ["A", "B", "B"], "regimen": None, "regimen_order": [1, 1, 2], "line": None,
                       "days_to_start": [30.0, 50.0, 70.0]})
    clock = treatment_clock(cohort, tr).set_index("patient_id")
    assert clock.loc["A", "entry_day"] == 70.0 and clock.loc["A", "exit_day"] == 870.0
    assert clock.loc["B", "entry_day"] == 0.0                       # the report came before the start
    assert not clock.loc["C", "in_analysis"]                        # no first line


def test_a_biomarker_mostly_unknown_uses_its_fallback_and_its_variants_dont_run():
    data = build_analysis(early_deaths(tables()), crc_settings(), bpc_profile())
    s1 = data.analysis("S1")
    assert s1.described["MSI-H/dMMR vs MSS/pMMR"] == "biomarker fallback"
    assert data.analysis("SA8").status == "not run" and data.analysis("SA9").status == "not run"


def test_a_biomarker_with_results_is_compared_and_its_variants_run():
    msi = {g: (n // 2, n - n // 2, 0) for g, n in GROUPS.items()}
    data = build_analysis(early_deaths(tables(msi=msi)), crc_settings(), bpc_profile())
    assert data.analysis("S1").described == {}
    assert data.analysis("SA8").status == "run" and data.analysis("SA9").status == "run"
    assert data.biomarker_counts["S1"]["Unknown"] == 0


def with_burden(t, high_ids=(), size=0.5):
    """Panel P8 given `size` megabases. Each patient in high_ids gets 9 more missense rows on the entry test: with
    the call, 10 over 0.5 megabases, which is S1's 20 per megabase. Everyone else has 2 per megabase."""
    rows = [{"test_id": f"{pid}-T01", "kind": "mutation", "gene": "GENE4", "variant": f"p.X{i}",
             "variant_class": "Missense_Mutation", "value": None} for pid in high_ids for i in range(9)]
    t = dataclasses.replace(t, results=pd.concat([t.results, pd.DataFrame(rows)], ignore_index=True),
                            panel_sizes=pd.DataFrame({"panel_id": ["P8"], "coding_megabases": [size]}))
    return conform(t)[0]


def test_a_biomarker_mostly_unknown_with_panel_sizes_compares_mutation_burden():
    high = ids("RAS-mutant", 20) + ids("Wild type", 15)                 # 35 high, 80 low
    data = build_analysis(with_burden(early_deaths(tables()), high), crc_settings(), bpc_profile(panel_sizes="given"))
    s1 = data.analysis("S1")
    assert s1.status == "run" and s1.grouping.levels == (BURDEN_LOW, BURDEN_HIGH) and s1.described == {}
    assert s1.adjusted and "hypermutation proxy" in s1.label and "MSI-H/dMMR vs MSS/pMMR" in s1.label
    assert data.biomarker_counts["S1"] == {BURDEN_LOW: 80, BURDEN_HIGH: 35, "Unknown": 0, "not in S1": 0}
    rows = data.d1[data.d1["analysis"] == "S1"]
    assert rows["comparison"].tolist() == [f"{BURDEN_HIGH} vs {BURDEN_LOW}"] and rows["decision"].tolist() == ["compare"]
    assert data.analysis("SA8").status == "not run" and data.analysis("SA9").status == "not run"
    spec = next(x for x in build_spec(data, crc_settings())["secondary"] if x["id"] == "S1")
    assert spec["grouping"]["levels"] == [BURDEN_LOW, BURDEN_HIGH] and spec["grouping"]["reference"] == BURDEN_LOW
    member = data.secondary.set_index("patient_id")[s1.member_column]
    assert member["RAS-0000"] == BURDEN_HIGH and member["BRA-0000"] == BURDEN_LOW


def test_without_sizes_in_the_profile_the_fallback_stays_descriptive():
    data = build_analysis(with_burden(early_deaths(tables()), ids("RAS-mutant", 35)), crc_settings(), bpc_profile())
    assert data.analysis("S1").described["MSI-H/dMMR vs MSS/pMMR"] == "biomarker fallback"


def test_a_short_burden_group_is_described_and_one_that_cant_be_counted_doesnt_run():
    data = build_analysis(with_burden(early_deaths(tables()), ids("RAS-mutant", 5)), crc_settings(),
                          bpc_profile(panel_sizes="given"))
    assert data.analysis("S1").described == {BURDEN_LOW: "smallest group", BURDEN_HIGH: "smallest group"}
    t = with_burden(early_deaths(tables()), ids("RAS-mutant", 35))
    t = dataclasses.replace(t, panel_sizes=pd.DataFrame({"panel_id": ["P8", "P8"], "coding_megabases": [0.5, 0.6]}))
    data = build_analysis(t, crc_settings(), bpc_profile(panel_sizes="given"))
    assert data.analysis("S1").status == "not run" and "two different coding sizes" in data.analysis("S1").why
    assert "S1" not in data.biomarker_counts and data.analysis("SA9").status == "not run"


def test_too_much_missing_needs_imputation_which_stops_the_run():
    t = early_deaths(tables())
    t.patients.loc[t.patients.index[:10], "sex"] = None             # 10 of 115: over the settings' 5%
    with pytest.raises(RunStop, match="multiple imputation"):
        build_analysis(t, crc_settings(), bpc_profile())


def test_the_wait_is_described_in_months():
    data = build_analysis(early_deaths(tables()), crc_settings(), bpc_profile())
    d = delays(data, crc_settings(), ["genotype", "center", "year_band"])
    row = d[d["by"] == "All patients"].iloc[0]
    assert row["median"] == pytest.approx(100 / (365.25 / 12))
    assert set(d["by"]) == {"All patients", "Genotype group", "Center", "Year band"}


# ---- The spec R reads: every value from the settings ----
def test_the_spec_holds_the_settings_values(tmp_path):
    s = settings_with(tmp_path, {"uncertainty.bootstrap_resamples": 77, "materiality.rmst_months": 1.5,
                                 "stability.risk_set": 15})
    data = build_analysis(early_deaths(tables()), s, bpc_profile())
    spec = build_spec(data, s)
    assert spec["resamples"] == 77 and spec["bars"]["rmst_diff"] == 1.5 and spec["stable_risk_set"] == 15
    assert spec["genotype"]["levels"][0] == "Wild type"             # the reference first
    assert spec["adjustment"]["knots"] == {"age_at_diagnosis": [10, 50, 90]}
    assert "never reaches 15" in labels(s)["unstable"]
    assert {a["id"] for a in spec["secondary"]} == {"S1", "S2", "SA1", "SA2", "SA3", "SA4", "SA7"}
