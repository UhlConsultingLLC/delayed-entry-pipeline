"""Fusion groups and the panels that detect them, groups described only, and the count of targeted drugs in the
reference group. Made-up data only: the CRC-shaped hand-made tables, with a made-up fusion gene (GENE6) and
made-up drug names."""
import dataclasses

import pandas as pd
import pytest

from delayed_entry.cohort import build_cohort, drug_matches, panel_coverage
from delayed_entry.screen import BLOCKED, COMPARE, DESCRIBE, DEVIATION, DEVIATIONS, MISSING, NOTE, PRIMARY, screen
from delayed_entry.settings import SettingsError
from delayed_entry.tables import conform, outcome_blind
from helpers import bpc_profile, crc_settings, make_tables, settings_with

S = crc_settings()
P = bpc_profile()
BALANCED = {"RAS-mutant": 60, "BRAF V600E": 40, "Wild type": 70}


def run(t, settings=S):
    return screen(outcome_blind(t), settings, P)


def with_fusions(t, group: str, n: int, gene: str = "GENE6"):
    """The first n patients of the group get a made-up fusion of the gene on their sequencing test."""
    tests = t.tests[(t.tests["modality"] == "ngs") & t.tests["patient_id"].str.startswith(group[:3])].head(n)
    rows = pd.DataFrame({"test_id": tests["test_id"], "kind": "fusion", "gene": gene,
                         "variant": f"PART1-{gene} fusion", "variant_class": None, "value": None}).astype("string")
    return conform(dataclasses.replace(t, results=pd.concat([t.results, rows], ignore_index=True)))[0]


def fusion_settings(tmp_path, changes: dict | None = None):
    """The CRC settings with a GENE6 fusion group before the reference, GENE6 required of every panel."""
    groups = S.read("genotype.groups")
    groups.insert(2, {"name": "GENE6 fusion", "any_of": [{"kind": "fusion", "genes": ["GENE6"]}]})
    return settings_with(tmp_path, {"genotype.groups": groups,
                                    "genotype.coverage.required": S.read("genotype.coverage.required") | {"GENE6": "whole_gene"},
                                    "smallest_group.if_short": S.read("smallest_group.if_short") | {"GENE6 fusion": "describe"},
                                    **(changes or {})})


# ---- The new choices in the CRC settings keep the registered plan as it was ----
def test_the_crc_settings_count_no_fusions_describe_no_group_and_count_no_drugs():
    assert S.read("genotype.coverage.fusion_panels") == {}
    assert S.read("genotype.described_only") == []
    assert S.read("genotype.reference_drugs") == {}


# ---- Fusion groups and fusion detection ----
def test_a_fusion_puts_its_patient_in_the_fusion_group(tmp_path):
    t = with_fusions(make_tables(BALANCED), "Wild type", 12)
    c = build_cohort(outcome_blind(t), fusion_settings(tmp_path))
    counts = c.eligible["genotype"].value_counts()
    assert counts["GENE6 fusion"] == 12 and counts["Wild type"] == 58


def test_a_fusion_gene_with_no_panel_list_rests_on_the_gene_list(tmp_path):
    t = with_fusions(make_tables(BALANCED), "Wild type", 12)
    rep = run(t, fusion_settings(tmp_path))                     # the CRC settings accept gene-level coverage
    item = rep.find(topic="fusion detection")[0]
    assert item.severity == NOTE and "name no panels that detect them" in item.finding
    s = fusion_settings(tmp_path, {"genotype.coverage.when_only_genes_listed": "deviation"})
    item = run(t, s).find(topic="fusion detection")[0]
    assert item.severity == DEVIATION and "stated as a deviation" in item.effect


def test_the_panel_list_decides_fusion_detection(tmp_path):
    t = with_fusions(make_tables(BALANCED), "Wild type", 12)
    rep = run(t, fusion_settings(tmp_path, {"genotype.coverage.fusion_panels": {"GENE6": ["P8"]}}))
    item = rep.find(topic="fusion detection")[0]
    assert item.severity == NOTE and "taken from the settings' list of panels" in item.finding
    assert not [w for w in rep.warnings if w.topic == "fusion detection"]
    assert rep.decisions[PRIMARY]["GENE6 fusion vs Wild type"] == DESCRIBE      # 12: short, and described


def test_a_panel_not_on_the_list_leaves_its_patients_unassignable(tmp_path):
    t = with_fusions(make_tables(BALANCED), "Wild type", 12)
    s = fusion_settings(tmp_path, {"genotype.coverage.fusion_panels": {"GENE6": ["P9"]}})
    c = build_cohort(outcome_blind(t), s)
    assert c.patients["genotype"].isna().all()                  # every made-up patient was sequenced on P8
    rep = run(t, s)
    warning = next(w for w in rep.warnings if w.topic == "fusion detection")
    assert "P9" in warning.finding and "misspelled" in warning.effect
    assert BLOCKED in rep.decisions[PRIMARY].values()


def test_panel_coverage_needs_both_the_gene_and_the_list():
    panels = pd.DataFrame({"panel_id": ["A", "A", "B", "B", "C"], "gene": ["G1", "G2", "G1", "G2", "G1"],
                           "codon_from": None, "codon_to": None})
    cov = panel_coverage(panels, {"G1": "whole_gene", "G2": "whole_gene"}, "accept_gene_level", {"G2": ["A", "C"]})
    assert cov["covered"].to_dict() == {"A": True, "B": False, "C": False}    # C is listed but doesn't read G2


def test_a_fusion_panel_list_for_a_gene_not_required_is_refused(tmp_path):
    with pytest.raises(SettingsError, match="'GENE6' is not a required gene"):
        settings_with(tmp_path, {"genotype.coverage.fusion_panels": {"GENE6": ["P8"]}})


# ---- Groups described only ----
def test_a_group_described_only_is_never_compared_whatever_its_size(tmp_path):
    s = settings_with(tmp_path, {"genotype.described_only": ["BRAF V600E"],
                                 "smallest_group.if_short": {"RAS-mutant": "block", "Wild type": "block"}})
    rep = run(make_tables(BALANCED), s)                         # 40 BRAF V600E: large enough to compare
    assert rep.decisions[PRIMARY] == {"RAS-mutant vs Wild type": COMPARE, "BRAF V600E vs Wild type": DESCRIBE}
    assert rep.decisions["SA1"]["BRAF V600E vs Wild type"] == DESCRIBE
    item = rep.find(analysis=PRIMARY, topic="described only")[0]
    assert item.severity == NOTE and item.settings == ("genotype.described_only",)
    assert rep.find(analysis=PRIMARY, topic="runs as standard")         # a planned description, not a deviation


@pytest.mark.parametrize("described, short, message", [
    (["BRAF"], None, "'BRAF' is not one of the groups"),
    (["Wild type"], None, "'Wild type' is the reference"),
    (["RAS-mutant", "BRAF V600E"], {"Wild type": "block"}, "at least one group besides the reference must be compared"),
    (["BRAF V600E"], None, "must name each group that is compared exactly once: RAS-mutant, Wild type"),
])
def test_what_described_only_may_name(tmp_path, described, short, message):
    changes = {"genotype.described_only": described}
    if short is not None:
        changes["smallest_group.if_short"] = short
    with pytest.raises(SettingsError, match=message):
        settings_with(tmp_path, changes)


# ---- Targeted drugs in the reference group ----
def with_regimens(t, group: str, k: int, text: str, rest: str = "Drug R"):
    """Every treated patient gets the regimen `rest`; the first k of the group get `text` instead."""
    tr = t.treatments.copy()
    tr["regimen"] = rest
    rows = tr.index[tr["patient_id"].str.startswith(group[:3])][:k]
    tr.loc[rows, "regimen"] = text
    return conform(dataclasses.replace(t, treatments=tr))[0]


TREATED = {"RAS-mutant": 60, "BRAF V600E": 40, "Wild type": 70}
DRUGS = {"genotype.reference_drugs": {"GENE6 fusion": ["Drug Q"]}}


@pytest.mark.parametrize("k, shown", [(15, "15"), (3, "<11"), (65, "all but fewer than 11"), (0, "none")])
def test_reference_patients_given_a_listed_drug_are_counted_and_masked(tmp_path, k, shown):
    t = with_regimens(make_tables(BALANCED, treated=TREATED), "Wild type", k, "Drug P, Drug Q Hydrochloride")
    item = run(t, settings_with(tmp_path, DRUGS)).find(analysis=PRIMARY, topic="targeted drugs in the reference")[0]
    assert item.list == DEVIATIONS and item.severity == NOTE
    assert item.finding.startswith(f"{shown} of the Wild type group's eligible patients received a drug the settings "
                                   "list for GENE6 fusion (Drug Q)")
    assert "they stay in the group" in item.effect


def test_counts_that_differ_by_a_small_number_are_shown_as_ranges(tmp_path):
    """Drug P's patients include Drug Q's, so two exact counts would give away the few given Drug P alone."""
    t = with_regimens(make_tables(BALANCED, treated=TREATED), "Wild type", 40, "Drug Q, Drug P")
    tr = t.treatments.copy()
    rows = tr.index[tr["patient_id"].str.startswith("Wil")][40:43]
    tr.loc[rows, "regimen"] = "Drug P"
    t = conform(dataclasses.replace(t, treatments=tr))[0]
    s = settings_with(tmp_path, {"genotype.reference_drugs": {"GENE6 fusion": ["Drug Q"], "GENE7 fusion": ["Drug P"]}})
    found = [i.finding for i in run(t, s).find(topic="targeted drugs in the reference")]
    assert [f.split(" of the ")[0] for f in found] == ["at least 30", "at least 30"]     # 40 and 43, not shown


def test_only_the_reference_group_is_counted(tmp_path):
    t = with_regimens(make_tables(BALANCED, treated=TREATED), "RAS-mutant", 40, "Drug Q")
    item = run(t, settings_with(tmp_path, DRUGS)).find(topic="targeted drugs in the reference")[0]
    assert item.finding.startswith("none of the Wild type group's")


def test_no_drug_names_in_the_treatments_table_is_a_stated_limit(tmp_path):
    rep = run(make_tables(BALANCED, treated=TREATED), settings_with(tmp_path, DRUGS))     # regimen left blank
    item = rep.find(topic="targeted drugs in the reference")[0]
    assert item.list == MISSING and item.severity == NOTE and "can't be counted" in item.effect


def test_no_drugs_listed_means_no_count():
    rep = run(with_regimens(make_tables(BALANCED, treated=TREATED), "Wild type", 15, "Drug Q"))
    assert not rep.find(topic="targeted drugs in the reference")


def test_a_drug_is_found_as_a_whole_name_ignoring_case():
    regimens = pd.Series(["Drug P, Drug Q", "drug q hydrochloride", "Pro-Drug Q", "Drug QX", None, "Drug P"],
                         dtype="string")
    assert drug_matches(regimens, ["Drug Q"]).tolist() == [True, True, False, False, False, False]
    assert drug_matches(regimens, ["Drug Q", "Drug P"]).tolist() == [True, True, False, False, False, True]
