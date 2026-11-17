"""The readiness screen: cases an independent review of the first build found. Made-up data only."""
import dataclasses
import re

import pandas as pd

from delayed_entry.screen import BLOCK, MISSING, NOT_RUN, PRIMARY, READY, screen
from delayed_entry.tables import outcome_blind
from helpers import bpc_profile, crc_settings, make_tables, settings_with

S = crc_settings()
P = bpc_profile()
BALANCED = {"RAS-mutant": 60, "BRAF V600E": 40, "Wild type": 70}


def run(t, settings=S, profile=P):
    return screen(outcome_blind(t), settings, profile)


def exact_numbers(text: str) -> list[int]:
    return [int(n) for n in re.findall(r"(?<![\w.<-])(\d+)(?![\w.%])", text)]


def add_result(t, pid, value, day=50.0):
    """Give a patient one more MSI/MMR result, dated `day`."""
    test = {"patient_id": pid, "test_id": f"{pid}-X{value}", "modality": "pathology", "entry_event": "report",
            "days_to_entry_event": day}
    res = {"test_id": f"{pid}-X{value}", "kind": "msi_mmr", "value": value}
    return dataclasses.replace(t, tests=pd.concat([t.tests, pd.DataFrame([test])], ignore_index=True),
                               results=pd.concat([t.results, pd.DataFrame([res])], ignore_index=True))


# ---- A hidden count can't be worked out by subtraction ----
def test_one_small_primary_group_turns_every_size_into_a_range():
    """Otherwise the flow's last count minus the two shown groups gives the hidden one."""
    rep = run(make_tables({"RAS-mutant": 200, "BRAF V600E": 7, "Wild type": 150}))
    sizes = rep.group_sizes[PRIMARY]
    assert sizes == {"RAS-mutant": "at least 30", "BRAF V600E": "<11", "Wild type": "at least 30"}
    for analysis, groups in rep.group_sizes.items():         # and so everywhere, once the primary is ranges
        assert not [n for v in groups.values() for n in exact_numbers(v) if n not in (0, 11, 30)], analysis


def test_a_biomarker_whose_exclusions_are_few_shows_ranges():
    t = make_tables(BALANCED, msi={g: (20, 20, 0) for g in BALANCED})
    for pid in ("RAS-0000", "RAS-0001", "RAS-0002"):                  # three left out of S1
        t.patients.loc[t.patients["patient_id"] == pid, "multiple_cohort_cancers"] = True
    sizes = run(t).group_sizes["S1"]
    assert all(not v.isdigit() for v in sizes.values())


def test_a_variant_close_to_its_biomarker_is_not_shown_exactly():
    t = make_tables(BALANCED, msi={g: (20, 20, 0) for g in BALANCED})
    for pid in ("RAS-0000", "RAS-0001", "RAS-0002"):                  # three discordant: MSI-H and MSS
        t = add_result(t, pid, "MSS")
    rep = run(t)
    assert rep.group_sizes["S1"]["MSI-H/dMMR"] == "60"
    assert rep.group_sizes["SA8"]["MSI-H/dMMR"].startswith("at least 30 (differs from its reference")


def test_missing_values_show_no_table_total():
    t = make_tables(BALANCED)
    t.patients.loc[:150, "sex"] = None
    rep = run(t)
    assert " of " not in " ".join(rep.missingness.values())
    assert rep.missingness["patients.sex"] == "151"
    t.patients["sex"] = None
    t.patients.loc[0, "sex"] = "Female"
    assert run(t).missingness["patients.sex"] == "all but fewer than 11"


def test_a_settings_minimum_below_the_agreement_never_weakens_the_masking(tmp_path):
    s = settings_with(tmp_path, {"privacy.min_count": 5})
    rep = run(make_tables({**BALANCED, "BRAF V600E": 7}), s)
    assert rep.min_count == 11 and rep.group_sizes[PRIMARY]["BRAF V600E"] == "<11"
    assert rep.find(topic="privacy minimum", severity=BLOCK)


# ---- Absent fields are reported, never a crash ----
def test_no_entry_times_at_all_runs_the_naive_clock_only():
    t = make_tables(BALANCED)
    t = dataclasses.replace(t, tests=t.tests.drop(columns=["days_to_entry_event"]))
    rep = run(t, S, bpc_profile(entry_event="none"))
    assert rep.find(topic="entry event", severity=BLOCK)
    assert rep.decisions[PRIMARY]["RAS-mutant vs Wild type"] == "compare"     # groups still compared, naive clock
    assert rep.decisions["SA1"] == {"all comparisons": NOT_RUN}
    assert rep.decisions["S1"] == {"all comparisons": NOT_RUN}


def test_no_modality_is_a_block_not_a_crash():
    t = make_tables(BALANCED)
    rep = run(dataclasses.replace(t, tests=t.tests.drop(columns=["modality"])))
    assert [i for i in rep.items if i.list == MISSING and i.severity == BLOCK and "tests.modality" in i.finding]
    assert rep.decisions == {}


def test_no_result_values_stop_the_biomarker_and_its_variants():
    t = make_tables(BALANCED)
    rep = run(dataclasses.replace(t, results=t.results.drop(columns=["value"])))
    assert rep.decisions["S1"] == rep.decisions["SA8"] == rep.decisions["SA9"] == {"all comparisons": NOT_RUN}
    ready = {i.analysis for i in rep.items if i.list == READY and i.topic == "runs as standard"}
    assert not ready & {"S1", "SA8", "SA9"}


def test_an_absent_key_field_stops_after_the_structure():
    t = make_tables(BALANCED)
    rep = run(dataclasses.replace(t, results=t.results.drop(columns=["kind"])))
    assert rep.find(topic="standard tables", severity=BLOCK) and rep.decisions == {}


def test_a_numeric_extra_compared_with_text_matches_nothing(tmp_path):
    t = make_tables(BALANCED)
    t.patients["score"] = 3.0
    analyses = S.read("analyses")
    next(a for a in analyses if a["id"] == "SA7")["keep"] = {"field": "score", "op": "in", "values": ["high"],
                                                             "missing": "fail"}
    rep = run(t, settings_with(tmp_path, {"analyses": analyses}))
    assert rep.group_sizes["SA7"] == {"RAS-mutant": "0", "BRAF V600E": "0", "Wild type": "0"}
