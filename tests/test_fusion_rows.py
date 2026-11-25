"""The fusion rows (option C) and question 2's ERBB2 and MET scope (none of the optional changes), both decided
2026-11-24: an ALK, ROS1 or RET fusion row counts when its description names a second gene, whatever the frame, and a
row describing a rearrangement inside the gene ("GENE-intragenic" at MSK and VICC, "GENE-INTRAGENIC" at DFCI, as in
the CRC release) doesn't. On made-up tumors, through the settings' rules and precedence. Every row here is made
up."""
import pytest

from delayed_entry.settings import load_settings
from test_nsclc_spec import NSCLC_SETTINGS
from test_q3_egfr_rules import case_id, group_of

COMMON, ALK, FIFTH, REF = "EGFR common", "ALK fusion", "Other recognized driver", "No recognized driver"
ALK_RULE = {"kind": "fusion", "genes": ["ALK"], "except_variants": ["ALK-intragenic", "ALK-INTRAGENIC"]}
ROS1_RET_RULE = {"kind": "fusion", "genes": ["ROS1", "RET"],
                 "except_variants": ["ROS1-intragenic", "ROS1-INTRAGENIC", "RET-intragenic", "RET-INTRAGENIC"]}


def fusion(gene, text):
    return ("fusion", gene, None, text)


@pytest.fixture(scope="module")
def decided():
    return load_settings(NSCLC_SETTINGS)


def test_the_fusion_rules_except_the_intragenic_rows_by_their_exact_text(decided):
    groups = {g["name"]: g.get("any_of", []) for g in decided.read("genotype.groups")}
    assert groups[ALK] == [ALK_RULE]
    assert groups[FIFTH][0] == ROS1_RET_RULE
    assert [m for rules in groups.values() for m in rules if m["kind"] == "fusion"] == [ALK_RULE, ROS1_RET_RULE]


def test_the_fusion_rows_and_the_scope_are_recorded_and_the_groups_decided(decided):
    groups = decided.provenance(["genotype.groups"])[0]
    assert groups["made_by"] == "NSCLC study author"
    assert "Fusion rows, decided 2026-11-24 (option C" in groups["source"]
    assert "ERBB2 and MET scope, decided 2026-11-24: none of the optional changes" in groups["source"]
    assert "TO DECIDE" not in groups["note"] and "RET rearrangements inside the gene" in groups["note"]
    assert "every ROS1 or RET fusion row counts" not in groups["note"]


@pytest.mark.parametrize("calls, group", [
    # a second gene named, in either order and either center's style: the fusion counts
    ([fusion("ALK", "EML4-ALK fusion")], ALK),
    ([fusion("ALK", "ALK-EML4")], ALK),                                    # DFCI style: the row's own gene first
    ([fusion("ALK", "EML4-ALK")], ALK),
    ([fusion("ALK", "ALK-EML4"), fusion("EML4", "EML4-ALK")], ALK),       # a DFCI pair: a row for each gene
    ([fusion("ROS1", "CD74-ROS1 fusion")], FIFTH),
    ([fusion("ROS1", "ROS1-CD74")], FIFTH),
    ([fusion("RET", "KIF5B-RET fusion")], FIFTH),
    ([fusion("RET", "RET-CCDC6")], FIFTH),
    # a rearrangement inside the gene, no second gene named: not counted, in both spellings
    ([fusion("ALK", "ALK-intragenic")], REF),
    ([fusion("ALK", "ALK-INTRAGENIC")], REF),
    ([fusion("ALK", " ALK-intragenic ")], REF),                            # the variant text is trimmed
    ([fusion("ROS1", "ROS1-intragenic")], REF),
    ([fusion("ROS1", "ROS1-INTRAGENIC")], REF),
    ([fusion("RET", "RET-intragenic")], REF),
    ([fusion("RET", "RET-INTRAGENIC")], REF),
    # an intragenic row takes no one out of the group another call puts them in
    ([fusion("ALK", "ALK-intragenic"), fusion("ROS1", "CD74-ROS1 fusion")], FIFTH),
    ([fusion("ALK", "ALK-INTRAGENIC"), fusion("ALK", "EML4-ALK fusion")], ALK),
    ([fusion("RET", "RET-INTRAGENIC"), ("mutation", "EGFR", "Missense_Mutation", "p.L858R")], COMMON),
    # only the exact texts are excepted: another spelling counts unless it is dropped, as a deviation, after the
    # release's descriptions are listed, in the check after download
    ([fusion("ALK", "ALK-Intragenic")], ALK),
], ids=case_id)
def test_an_intragenic_row_is_not_counted_and_a_partner_named_row_is(decided, calls, group):
    assert group_of(decided, calls) == group
