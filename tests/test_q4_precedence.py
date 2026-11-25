"""Question 4, option A (decided 2026-11-24): a tumor whose first test matches two or more groups goes to the first
listed, in the draft order: EGFR common, ALK fusion, KRAS (one group, KRAS mutant, since question 5 joined the two
drafted KRAS groups), then Other recognized driver. On made-up tumors: each pair of main groups, in both orders of
the calls. Every call here is made up."""
import pytest

from delayed_entry.settings import load_settings
from test_nsclc_spec import NSCLC_SETTINGS
from test_q3_egfr_rules import case_id, group_of

COMMON, ALK, KRAS, FIFTH, REF = ("EGFR common", "ALK fusion", "KRAS mutant", "Other recognized driver",
                                 "No recognized driver")

DEL19 = ("mutation", "EGFR", "In_Frame_Del", "p.E746_A750del")
L858R = ("mutation", "EGFR", "Missense_Mutation", "p.L858R")
EML4_ALK = ("fusion", "ALK", None, "EML4-ALK fusion")
KRAS_G12C = ("mutation", "KRAS", "Missense_Mutation", "p.G12C")
G12D = ("mutation", "KRAS", "Missense_Mutation", "p.G12D")
G12V = ("mutation", "KRAS", "Missense_Mutation", "p.G12V")
G13D = ("mutation", "KRAS", "Missense_Mutation", "p.G13D")
Q61H = ("mutation", "KRAS", "Missense_Mutation", "p.Q61H")
V600E = ("mutation", "BRAF", "Missense_Mutation", "p.V600E")
CD74_ROS1 = ("fusion", "ROS1", None, "CD74-ROS1 fusion")
T790M = ("mutation", "EGFR", "Missense_Mutation", "p.T790M")


@pytest.fixture(scope="module")
def q4():
    return load_settings(NSCLC_SETTINGS)


def test_question_4_is_recorded_and_the_order_holds(q4):
    assert [g["name"] for g in q4.read("genotype.groups")] == [COMMON, ALK, KRAS, FIFTH, REF]
    groups = q4.provenance(["genotype.groups"])[0]
    assert groups["made_by"] == "NSCLC study author"      # every part decided since the fusion rows (2026-11-24)
    assert "Question 4, decided 2026-11-24 (option A)" in groups["source"]
    assert "Question 4 decided 2026-11-24 (option A)" in groups["note"]
    assert "questions 4-5" not in groups["note"] and "TO DECIDE: question 5" not in groups["note"]
    analyses = q4.provenance(["analyses"])[0]
    assert "SA3" not in {a["id"] for a in q4.read("analyses")}
    assert "question 4, decided 2026-11-24, option A" in analyses["source"]
    assert "SA3 stays out (question 4, decided 2026-11-24, option A)" in analyses["note"]


@pytest.mark.parametrize("calls, group", [
    ([DEL19], COMMON), ([L858R], COMMON), ([EML4_ALK], ALK), ([KRAS_G12C], KRAS),
    ([G12D], KRAS), ([G13D], KRAS), ([Q61H], KRAS), ([V600E], FIFTH), ([CD74_ROS1], FIFTH),
], ids=case_id)
def test_each_made_up_call_alone_lands_in_its_own_group(q4, calls, group):
    assert group_of(q4, calls) == group


@pytest.mark.parametrize("reverse", [False, True], ids=["as listed", "reversed"])
@pytest.mark.parametrize("calls, group", [
    ([DEL19, EML4_ALK], COMMON),
    ([L858R, EML4_ALK], COMMON),
    ([DEL19, KRAS_G12C], COMMON),
    ([L858R, Q61H], COMMON),
    ([DEL19, G12D], COMMON),
    ([EML4_ALK, KRAS_G12C], ALK),
    ([EML4_ALK, G12D], ALK),
], ids=case_id)
def test_each_pair_of_main_groups_goes_to_the_first_listed(q4, calls, group, reverse):
    assert group_of(q4, calls[::-1] if reverse else calls) == group


@pytest.mark.parametrize("calls, group", [
    ([L858R, EML4_ALK, KRAS_G12C], COMMON),     # three main groups
    ([EML4_ALK, KRAS_G12C, G13D], ALK),
    ([EML4_ALK, CD74_ROS1], ALK),               # a main group with the fifth
    ([KRAS_G12C, CD74_ROS1], KRAS),
    ([G12D, V600E], KRAS),
    ([L858R, V600E], COMMON),
    ([DEL19, T790M, KRAS_G12C], COMMON),        # a same-gene and a cross-gene overlap together
], ids=case_id)
def test_more_than_two_matches_and_main_groups_over_the_fifth(q4, calls, group):
    assert group_of(q4, calls) == group


@pytest.mark.parametrize("calls, group", [
    ([KRAS_G12C, G13D], KRAS),
    ([KRAS_G12C, Q61H], KRAS),
    ([KRAS_G12C, G12V], KRAS),                  # G12F (c.34_35GG>TT) written as two single-base records
    ([KRAS_G12C, G12D], KRAS),                  # G12Y (c.34_35GG>TA) written the same way
    ([("mutation", "KRAS", "Missense_Mutation", "p.G12F")], KRAS),       # G12F written as one record
], ids=case_id)
def test_two_kras_calls_stay_in_the_one_kras_group(q4, calls, group):
    """Since question 5 there is one KRAS group, so a G12F or G12Y split into two single-base records is a KRAS
    mutant tumor either way."""
    assert group_of(q4, calls) == group
