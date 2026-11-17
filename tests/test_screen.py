"""The readiness screen (design plan section 6), on made-up data only."""
import dataclasses
import re

import pandas as pd
import pytest

from delayed_entry.screen import (BLOCK, BLOCKED, COMPARE, DESCRIBE, DEVIATION, DEVIATIONS, MISSING, NOT_RUN, NOTE,
                                  PRIMARY, READY, screen)
from delayed_entry.tables import OUTCOME_FIELDS, outcome_blind
from helpers import bpc_profile, crc_settings, make_tables, scrambled, settings_with

S = crc_settings()
P = bpc_profile()
BALANCED = {"RAS-mutant": 60, "BRAF V600E": 40, "Wild type": 70}


def run(t, settings=S, profile=P):
    return screen(outcome_blind(t), settings, profile)


# ---- It never reads an outcome ----
@pytest.mark.parametrize("seed", [1, 2, 3])
def test_scrambling_every_outcome_value_leaves_the_report_unchanged(fixture_result, seed):
    r = fixture_result
    before = screen(outcome_blind(r.tables), S, r.profile, r.source_steps)
    after = screen(outcome_blind(scrambled(r.tables, seed)), S, r.profile, r.source_steps)
    assert after.to_dict() == before.to_dict()
    assert after.to_markdown() == before.to_markdown()


def test_scrambling_changes_the_outcomes_themselves(fixture_result):
    """So the test above is not empty: the scrambled tables really do differ."""
    t = fixture_result.tables
    s = scrambled(t, 1)
    assert not s.patients["days_to_death_or_last_followup"].equals(t.patients["days_to_death_or_last_followup"])
    assert s.patients["days_to_death_or_last_followup"].isna().equals(t.patients["days_to_death_or_last_followup"].isna())


def test_the_screen_refuses_tables_with_outcomes(fixture_result):
    with pytest.raises(TypeError, match="outcome_blind"):
        screen(fixture_result.tables, S, P)
    blind = outcome_blind(fixture_result.tables)
    sneaked = dataclasses.replace(blind, patients=blind.patients.assign(vital_status=1))
    with pytest.raises(TypeError, match="outcome columns reached the screen"):
        screen(sneaked, S, P)


def test_the_outcome_blind_tables_hold_no_outcome_value(fixture_result):
    blind = outcome_blind(fixture_result.tables)
    for name in ("patients", "tests", "results", "panels", "treatments"):
        assert not set(OUTCOME_FIELDS) & set(blind.table(name).columns)
    assert set(blind.patients["outcome_recorded"].dropna().unique()) <= {True, False}


# ---- Thresholds come from the settings file ----
def test_a_changed_threshold_in_the_settings_file_changes_the_result(tmp_path, fixture_result):
    """The same data and the same code: only the settings file differs."""
    blind, prof = outcome_blind(fixture_result.tables), fixture_result.profile
    registered = screen(blind, S, prof)
    assert registered.decisions[PRIMARY]["BRAF V600E vs Wild type"] == COMPARE
    assert registered.decisions["SA1"]["BRAF V600E vs Wild type"] == DESCRIBE
    stricter = screen(blind, settings_with(tmp_path, {"smallest_group.patients": 60}, "strict.yaml"), prof)
    assert stricter.decisions[PRIMARY]["BRAF V600E vs Wild type"] == DESCRIBE
    looser = screen(blind, settings_with(tmp_path, {"smallest_group.patients": 5}, "loose.yaml"), prof)
    assert looser.decisions["SA1"]["BRAF V600E vs Wild type"] == COMPARE
    assert looser.decisions["S2"]["BRAF V600E vs Wild type"] == COMPARE


def test_a_changed_short_wait_cutoff_changes_who_sa1_keeps(tmp_path):
    t = make_tables(BALANCED, report_day=200.0)                 # every report at day 200
    assert run(t).decisions["SA1"]["RAS-mutant vs Wild type"] == DESCRIBE        # 6 months is 182.6 days: nobody
    analyses = S.read("analyses")
    next(a for a in analyses if a["id"] == "SA1")["keep"] = {"within_months": 7}
    assert run(t, settings_with(tmp_path, {"analyses": analyses})).decisions["SA1"]["RAS-mutant vs Wild type"] == COMPARE


def test_every_setting_the_screen_read_is_listed_with_its_origin(fixture_result):
    rep = screen(outcome_blind(fixture_result.tables), S, fixture_result.profile)
    listed = {s["choice"]: s for s in rep.settings_used}
    for choice in ("smallest_group.patients", "biomarker.unknown_max_share", "privacy.min_count", "cohort.rules",
                   "genotype.groups", "entry.first_test_order", "survival_after_entry.window_days"):
        assert listed[choice]["origin"].startswith("settings file, line"), choice
        assert listed[choice]["made_by"] and listed[choice]["source"]
    assert rep.defaults_used == []


def test_every_default_the_screen_applies_is_listed(tmp_path, fixture_result):
    s = settings_with(tmp_path, {"smallest_group.patients": "default", "entry.first_test_order": "default"})
    rep = screen(outcome_blind(fixture_result.tables), s, fixture_result.profile)
    used = {d["choice"]: d for d in rep.defaults_used}
    assert set(used) == {"smallest_group.patients", "entry.first_test_order"}
    assert all("CRC BPC plan" in d["source"] for d in used.values())
    assert "## Defaults applied (2)" in rep.to_markdown()


# ---- The smallest-group rule (D1), from counts alone ----
def test_balanced_groups_are_all_compared():
    rep = run(make_tables(BALANCED))
    assert rep.decisions[PRIMARY] == {"RAS-mutant vs Wild type": COMPARE, "BRAF V600E vs Wild type": COMPARE}
    assert rep.find(analysis=PRIMARY, topic="runs as standard")


def test_a_short_braf_group_is_described_not_compared():
    rep = run(make_tables({**BALANCED, "BRAF V600E": 25}))
    assert rep.decisions[PRIMARY]["BRAF V600E vs Wild type"] == DESCRIBE
    assert rep.decisions[PRIMARY]["RAS-mutant vs Wild type"] == COMPARE
    item = rep.find(analysis=PRIMARY, topic="smallest group")[0]
    assert item.severity == DEVIATION and "BRAF V600E has fewer than 30" in item.finding


@pytest.mark.parametrize("short", ["Wild type", "RAS-mutant"])
def test_a_short_group_the_settings_mark_block_blocks(short):
    rep = run(make_tables({**BALANCED, short: 20}))
    assert BLOCKED in rep.decisions[PRIMARY].values()
    item = rep.find(analysis=PRIMARY, topic="smallest group", severity=BLOCK)[0]
    assert "a person decides" in item.effect
    assert not rep.find(analysis=PRIMARY, topic="runs as standard")


def test_what_a_short_group_does_comes_from_the_settings(tmp_path):
    s = settings_with(tmp_path, {"smallest_group.if_short": {"RAS-mutant": "describe", "BRAF V600E": "block",
                                                             "Wild type": "block"}})
    rep = run(make_tables({**BALANCED, "BRAF V600E": 25}), s)
    assert rep.decisions[PRIMARY]["BRAF V600E vs Wild type"] == BLOCKED


def test_biomarker_groups_get_the_rule_too():
    t = make_tables(BALANCED, msi={"RAS-mutant": (5, 50, 5), "BRAF V600E": (10, 25, 5), "Wild type": (5, 60, 5)})
    rep = run(t)
    assert rep.decisions["S1"] == {"MSI-H/dMMR vs MSS/pMMR": DESCRIBE}        # 20 MSI-H: fewer than 30
    assert rep.decisions["SA8"] == {"MSI-H/dMMR vs MSS/pMMR": DESCRIBE}


def test_treatment_clock_counts_patients_with_a_first_line_regimen():
    t = make_tables(BALANCED, treated={"RAS-mutant": 60, "BRAF V600E": 20, "Wild type": 70})
    assert run(t).decisions["S2"] == {"RAS-mutant vs Wild type": COMPARE, "BRAF V600E vs Wild type": DESCRIBE}


def test_no_treatments_table_blocks_the_treatment_clock():
    t = dataclasses.replace(make_tables(BALANCED), treatments=None)
    rep = run(t)
    assert rep.decisions["S2"] == {"all comparisons": NOT_RUN}
    assert rep.find(analysis="S2", severity=BLOCK)


def test_an_analysis_switched_off_is_not_run(tmp_path):
    analyses = S.read("analyses")
    next(a for a in analyses if a["id"] == "SA7")["run"] = False
    rep = run(make_tables(BALANCED), settings_with(tmp_path, {"analyses": analyses}))
    assert rep.decisions["SA7"] == {"all comparisons": "switched off"}


# ---- Small counts ----
def shown_numbers(rep) -> list[int]:
    nums = [v for g in rep.group_sizes.values() for v in g.values()]
    nums += [str(r[k]) for r in (rep.flow or []) for k in ("excluded", "remaining")]
    nums += list(rep.missingness.values())
    nums += [i.finding for i in rep.items + rep.warnings]
    return [int(n) for text in nums for n in re.findall(r"(?<![\w.<-])(\d+)(?![\w.%])", str(text))]


def test_small_counts_are_masked_and_flagged():
    rep = run(make_tables({**BALANCED, "BRAF V600E": 6}))
    assert rep.group_sizes[PRIMARY]["BRAF V600E"] == "<11"
    assert rep.find(analysis=PRIMARY, topic="small count")
    assert not [n for n in shown_numbers(rep) if 1 <= n <= 10]


def test_the_fixture_report_shows_no_small_count(fixture_result):
    rep = screen(outcome_blind(fixture_result.tables), S, fixture_result.profile, fixture_result.source_steps)
    assert not [n for n in shown_numbers(rep) if 1 <= n <= 10]


def test_a_sensitivity_count_close_to_the_primary_one_is_not_shown():
    """Shown beside the primary count, it would give away the few patients the subset leaves out."""
    t = make_tables(BALANCED)
    t.patients.loc[t.patients["patient_id"] == "BRA-0000", "histology"] = "Carcinoma"   # SA7 drops one
    rep = run(t)
    assert rep.group_sizes["SA7"]["BRAF V600E"].startswith("at least 30 (differs from its reference")
    assert rep.group_sizes["SA7"]["RAS-mutant"] == "60"


# ---- The capability profile ----
def test_sample_collection_is_a_deviation_when_the_settings_accept_it(tmp_path):
    t = make_tables(BALANCED)
    t.tests["entry_event"] = "sample_collection"
    prof = bpc_profile(entry_event="sample_collection")
    refused = run(t, S, prof)
    assert refused.find(analysis=PRIMARY, topic="entry event", severity=BLOCK)
    s = settings_with(tmp_path, {"entry.acceptable_substitutes": ["sample_collection"]})
    item = run(t, s, prof).find(analysis=PRIMARY, topic="entry event")[0]
    assert item.severity == DEVIATION and "entry at sample collection" in item.effect


def test_no_entry_event_at_all_is_a_loud_block():
    item = run(make_tables(BALANCED), S, bpc_profile(entry_event="none")).find(topic="entry event")[0]
    assert item.severity == BLOCK and "NO ENTRY EVENT" in item.finding and "naive clock only" in item.effect


def test_no_panel_coverage_leaves_wild_type_unconfirmed():
    t = dataclasses.replace(make_tables(BALANCED), panels=make_tables(BALANCED).panels.iloc[0:0])
    rep = run(t, S, bpc_profile(panel_coverage="none"))
    item = rep.find(topic="panel coverage")[0]
    assert item.severity == DEVIATION and "unconfirmed" in item.effect
    assert rep.decisions[PRIMARY]["RAS-mutant vs Wild type"] == COMPARE      # groups still assigned


def test_gene_level_coverage_as_a_deviation_when_the_settings_say_so(tmp_path):
    s = settings_with(tmp_path, {"genotype.coverage.when_only_genes_listed": "deviation"})
    assert run(make_tables(BALANCED), s).find(topic="panel coverage")[0].severity == DEVIATION
    assert run(make_tables(BALANCED)).find(topic="panel coverage")[0].severity == NOTE


@pytest.mark.parametrize("capture", ["partial", "unknown"])
def test_incomplete_death_capture_is_a_stated_limit(capture):
    item = run(make_tables(BALANCED), S, bpc_profile(death_capture=capture)).find(topic="death capture")[0]
    assert item.list == DEVIATIONS and item.severity == DEVIATION


def test_the_privacy_minimum_can_never_be_below_the_agreement(tmp_path):
    rep = run(make_tables(BALANCED), settings_with(tmp_path, {"privacy.min_count": 5}))
    assert rep.find(topic="privacy minimum", severity=BLOCK)
    raised = run(make_tables(BALANCED), settings_with(tmp_path, {"privacy.min_count": 20}, "raised.yaml"))
    assert not raised.find(topic="privacy minimum")


# ---- Data checks and missing fields ----
def test_a_field_a_cohort_rule_needs_is_reported_missing():
    t = make_tables(BALANCED)
    t = dataclasses.replace(t, patients=t.patients.drop(columns=["primary_site_code"]))
    rep = run(t)
    item = rep.find(topic="cohort rule")[0]
    assert item.list == MISSING and item.severity == BLOCK and "primary_site_code" in item.finding


def test_an_absent_adjustment_variable_is_reported():
    t = make_tables(BALANCED)
    t = dataclasses.replace(t, patients=t.patients.drop(columns=["tumor_location"]))
    item = run(t).find(topic="needed for the adjusted models")[0]
    assert "tumor_location" in item.finding


def test_no_record_of_a_report_after_follow_up_is_a_note():
    t = make_tables(BALANCED)
    t = dataclasses.replace(t, tests=t.tests.drop(columns=["after_follow_up_end"]))
    item = run(t).find(topic="time at risk")[0]
    assert item.list == MISSING and item.severity == NOTE and "upper bounds" in item.effect


def test_no_outcome_columns_blocks_every_survival_analysis():
    t = make_tables(BALANCED)
    t = dataclasses.replace(t, patients=t.patients.drop(columns=list(OUTCOME_FIELDS)))
    rep = run(t)
    assert rep.find(topic="outcome", severity=BLOCK)


def test_duplicate_patients_block():
    t = make_tables(BALANCED)
    t = dataclasses.replace(t, patients=pd.concat([t.patients, t.patients.iloc[:1]], ignore_index=True))
    assert run(t).find(topic="duplicates", severity=BLOCK)


def test_a_report_before_diagnosis_enters_at_day_0_as_a_note():
    t = make_tables(BALANCED, report_day=-15.0)
    item = run(t).find(topic="entry before diagnosis")[0]
    assert item.severity == NOTE and item.effect == "they enter at day 0"


def test_missing_values_are_counted_for_every_field():
    rep = run(make_tables(BALANCED))
    assert "patients.sex" in rep.missingness and "tests.panel_id" in rep.missingness


# ---- Warnings, never overrides ----
def test_a_smallest_group_below_the_privacy_minimum_is_warned(tmp_path):
    rep = run(make_tables(BALANCED), settings_with(tmp_path, {"smallest_group.patients": 8}))
    assert [w for w in rep.warnings if w.topic == "smallest group"]
    assert rep.decisions[PRIMARY]["BRAF V600E vs Wild type"] == COMPARE      # warned, not overridden


def test_a_survival_check_window_past_the_horizon_is_warned(tmp_path):
    rep = run(make_tables(BALANCED), settings_with(tmp_path, {"survival_after_entry.window_days": 1200}))
    assert [w for w in rep.warnings if w.topic == "survival-after-entry check"]


def test_an_age_rule_above_the_top_code_is_warned(tmp_path):
    rules = S.read("cohort.rules")
    rules[2]["all_of"][0]["value"] = 95
    rep = run(make_tables(BALANCED), settings_with(tmp_path, {"cohort.rules": rules}))
    assert [w for w in rep.warnings if w.topic == "cohort rule"]


def test_the_registered_settings_raise_no_warning():
    assert run(make_tables(BALANCED)).warnings == []


# ---- What needs outcomes is listed, not done ----
def test_checks_that_need_outcomes_are_listed_with_their_settings(tmp_path):
    rep = run(make_tables(BALANCED), settings_with(tmp_path, {"survival_after_entry.window_days": 365}))
    first = rep.deferred[0]
    assert first["check"] == "survival after the entry event" and "within 365 days" in first["rule"]
    assert {d["check"] for d in rep.deferred} >= {"time at risk from the times themselves", "the horizon against follow-up"}


# ---- The three lists ----
def test_a_clean_dataset_is_ready_for_every_analysis():
    t = make_tables(BALANCED, msi={g: (35, 20, 5) for g in BALANCED},
                    treated={g: n for g, n in BALANCED.items()})
    rep = run(t)
    ready = {i.analysis for i in rep.list(READY) if i.topic == "runs as standard"}
    assert ready == {PRIMARY, "S1", "S2", "S3", "SA1", "SA2", "SA3", "SA4", "SA5", "SA6", "SA7", "SA8", "SA9"}


def test_the_report_is_the_same_every_time(fixture_result):
    a = screen(outcome_blind(fixture_result.tables), S, fixture_result.profile)
    b = screen(outcome_blind(fixture_result.tables), S, fixture_result.profile)
    assert a.to_yaml() == b.to_yaml()


def test_the_fixture_screen_matches_the_crc_projects_design_checks(fixture_result):
    """On this made-up fixture the CRC project's own D1 compares everything except BRAF V600E in S2 and SA1."""
    rep = screen(outcome_blind(fixture_result.tables), S, fixture_result.profile, fixture_result.source_steps)
    described = {(a, c) for a, comps in rep.decisions.items() for c, d in comps.items() if d != COMPARE}
    assert described == {("S2", "BRAF V600E vs Wild type"), ("SA1", "BRAF V600E vs Wild type")}
    assert rep.group_sizes[PRIMARY] == {"RAS-mutant": "324", "BRAF V600E": "38", "Wild type": "296"}
