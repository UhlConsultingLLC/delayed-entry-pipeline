"""Question 3, option C (decided 2026-11-23): EGFR common stays as drafted, and Other recognized driver gains three
EGFR rules (missense changes at six codons minus 19 single-base changes, in-frame insertions in exons 19 and 20,
and in-frame deletions touching codon 709 or 710). On made-up tumors and on the made-up NSCLC-shaped files. Every
call here is made up."""
import shutil

import pandas as pd
import pytest

from delayed_entry.adapters import genie_bpc
from delayed_entry.cohort import assign_genotype, build_cohort
from delayed_entry.settings import load_settings
from delayed_entry.tables import OutcomeBlindTables, outcome_blind
from test_nsclc_spec import NSCLC_SETTINGS, nsclc_files, nsclc, read, write  # noqa: F401

COMMON, FIFTH, REF = "EGFR common", "Other recognized driver", "No recognized driver"

# Table 12's single-base EGFR changes (Lindeman 2013), without L858R, which EGFR common takes first
KEPT = ["p.E709K", "p.E709A", "p.E709G", "p.E709V", "p.E709D", "p.E709Q", "p.G719A", "p.G719C", "p.G719S",
        "p.G719D", "p.S768I", "p.T790M", "p.L858M", "p.L861Q", "p.L861R"]
EXCEPTED = ["p.G719R", "p.G719V",
            "p.S768C", "p.S768G", "p.S768N", "p.S768R", "p.S768T",
            "p.T790A", "p.T790K", "p.T790P", "p.T790R", "p.T790S",
            "p.L858P", "p.L858Q", "p.L858R", "p.L858V",
            "p.L861M", "p.L861P", "p.L861V"]
EGFR_RULES = [
    {"kind": "mutation", "genes": ["EGFR"], "variant_classes": ["Missense_Mutation"],
     "codons": [709, 719, 768, 790, 858, 861], "except_variants": EXCEPTED},
    {"kind": "mutation", "genes": ["EGFR"], "variant_classes": ["In_Frame_Ins"], "codon_ranges": [[729, 761], [762, 823]]},
    {"kind": "mutation", "genes": ["EGFR"], "variant_classes": ["In_Frame_Del"], "codon_ranges": [[709, 710]]},
]


def egfr(variant_class, variant):
    return ("mutation", "EGFR", variant_class, variant)


def case_id(value):
    """A readable test id: the tumor's calls, or the group."""
    return " + ".join(f"{gene} {variant}" for _, gene, _, variant in value) if isinstance(value, list) else value


@pytest.fixture(scope="module")
def q3():
    return load_settings(NSCLC_SETTINGS)


def group_of(settings, calls):
    """The group one made-up tumor lands in by the settings' rules and precedence (assign_genotype's group before
    the coverage rule, which is a separate step)."""
    first = pd.DataFrame({"test_id": ["T1"], "panel_id": ["MADE-UP-PANEL"]}, index=pd.Index(["P1"], name="patient_id"))
    results = pd.DataFrame([("T1", *c) for c in calls], columns=["test_id", "kind", "gene", "variant_class", "variant"])
    panels = pd.DataFrame({"panel_id": "MADE-UP-PANEL", "gene": ["EGFR", "KRAS", "ALK"]})
    blind = OutcomeBlindTables(patients=pd.DataFrame({"patient_id": ["P1"]}), tests=first.reset_index(),
                               results=results, panels=panels, treatments=None, panel_sizes=None,
                               outcome_fields_present=())
    return assign_genotype(blind, settings, first).loc["P1", "group_by_rules"]


def test_the_three_egfr_rules_are_written_exactly_after_question_2s_rules(q3):
    groups = {g["name"]: g for g in q3.read("genotype.groups")}
    fifth = groups[FIFTH]["any_of"]
    assert fifth[-3:] == EGFR_RULES
    assert [m for m in fifth if "EGFR" in m["genes"]] == EGFR_RULES          # no other EGFR rule in the group
    assert len(fifth) == 10                                                    # Q2's seven rules, then these three
    assert len(EXCEPTED) == len(set(EXCEPTED)) == 19 and not set(KEPT) & set(EXCEPTED)
    assert groups[COMMON]["any_of"] == [                                       # EGFR common as drafted
        {"kind": "mutation", "genes": ["EGFR"], "variant_classes": ["In_Frame_Del"], "codon_ranges": [[729, 761]]},
        {"kind": "mutation", "genes": ["EGFR"], "variants": ["p.L858R"]}]
    assert q3.read("genotype.coverage.required")["EGFR"] == list(range(729, 762)) + [858]   # coverage unchanged


def test_question_3_is_recorded_in_the_decided_groups(q3):
    groups = q3.provenance(["genotype.groups"])[0]
    assert groups["made_by"] == "NSCLC study author"      # every part decided since the fusion rows (2026-11-24)
    assert "Question 3, decided 2026-11-23 (option C)" in groups["source"]
    assert "Question 3 decided 2026-11-23 (option C)" in groups["note"] and "TO DECIDE" not in groups["note"]
    assert "questions 3-5" not in groups["note"] and "NM_005228 / ENST00000275493 (EGFR)" in groups["note"]
    coverage = q3.provenance(["genotype.coverage.required"])[0]
    assert "follows question 3" not in coverage["note"] and "Question 3 (decided 2026-11-23" in coverage["note"]


@pytest.mark.parametrize("variant", KEPT)
def test_each_kept_single_base_change_lands_in_the_fifth_group(q3, variant):
    assert group_of(q3, [egfr("Missense_Mutation", variant)]) == FIFTH


@pytest.mark.parametrize("variant", EXCEPTED)
def test_each_excepted_change_alone_stays_in_the_reference_but_l858r_is_egfr_common(q3, variant):
    assert group_of(q3, [egfr("Missense_Mutation", variant)]) == (COMMON if variant == "p.L858R" else REF)


@pytest.mark.parametrize("calls, group", [
    # exon 19 deletions stay EGFR common
    ([egfr("In_Frame_Del", "p.E746_A750del")], COMMON),
    # exon 20 insertions, in the three spellings seen in public data, all classed In_Frame_Ins
    ([egfr("In_Frame_Ins", "p.A767_V769dup")], FIFTH),
    ([egfr("In_Frame_Ins", "p.D770_N771insSVD")], FIFTH),
    ([egfr("In_Frame_Ins", "p.N771delinsGY")], FIFTH),
    # an exon 19 insertion
    ([egfr("In_Frame_Ins", "p.I744_K745insKIPVAI")], FIFTH),
    # insertions at an exon boundary count (exons 19 and 20 share the group; a rare exon 18 or 21 one is counted too)
    ([egfr("In_Frame_Ins", "p.D761_E762insEAFQ")], FIFTH),
    ([egfr("In_Frame_Ins", "p.K728_G729insX")], FIFTH),
    ([egfr("In_Frame_Ins", "p.K823_G824insX")], FIFTH),
    # exon 18 and exon 21 insertions away from the boundary codons stay in the reference
    ([egfr("In_Frame_Ins", "p.E709_T710insX")], REF),
    ([egfr("In_Frame_Ins", "p.L862_G863insL")], REF),
    # exon 18 deletions touching codon 709 or 710 (E709_T710delinsD) count; other exon 18 deletions don't
    ([egfr("In_Frame_Del", "p.E709_T710delinsD")], FIFTH),
    ([egfr("In_Frame_Del", "p.E709_T710del")], FIFTH),
    ([egfr("In_Frame_Del", "p.T710_L718del")], FIFTH),               # a longer one reaching codon 710
    ([egfr("In_Frame_Del", "p.S720_A722del")], REF),
    # changes that swap codons for the same number of new ones are classed missense: they count only at the codons
    ([egfr("Missense_Mutation", "p.E746_L747delinsIP")], REF),
    ([egfr("Missense_Mutation", "p.H773_V774delinsLM")], REF),
    ([egfr("Missense_Mutation", "p.S768_V769delinsIL")], FIFTH),
    ([egfr("Missense_Mutation", "p.L858_A859delinsRS")], FIFTH),       # L858R written as two codons: not EGFR common
    # two-base changes in one codon aren't in the except list, so the missense rule keeps them
    ([egfr("Missense_Mutation", "p.T790L")], FIFTH),
    ([egfr("Missense_Mutation", "p.G719F")], FIFTH),
    # complex exon 19 changes split by net length
    ([egfr("In_Frame_Ins", "p.L747_E749delinsPKAN")], FIFTH),           # a net gain
    ([egfr("In_Frame_Del", "p.L747_P753delinsS")], COMMON),             # a net loss
    # what C leaves in the reference
    ([egfr("Missense_Mutation", "p.L747P")], REF),                      # an exon 19 point change
    ([egfr("In_Frame_Del", "p.V769_D770del")], REF),                    # an exon 20 in-frame deletion
    ([egfr("Missense_Mutation", "p.R776H")], REF),                      # a codon the rule doesn't name
    ([("fusion", "EGFR", None, "EGFR-RAD51 fusion")], REF),              # an EGFR fusion row
    ([("mutation", "SOS1", "Missense_Mutation", "p.G719A")], REF),       # the label in another gene
    # two calls in one tumor: precedence
    ([egfr("Missense_Mutation", "p.T790M")], FIFTH),                    # de novo T790M
    ([egfr("Missense_Mutation", "p.L858R"), egfr("Missense_Mutation", "p.T790M")], COMMON),
    ([egfr("In_Frame_Del", "p.E746_A750del"), egfr("Missense_Mutation", "p.T790M")], COMMON),
    ([egfr("Missense_Mutation", "p.L858R"), egfr("Missense_Mutation", "p.G719A")], COMMON),
    ([egfr("Missense_Mutation", "p.G719A"), egfr("Missense_Mutation", "p.S768I")], FIFTH),
    ([("mutation", "KRAS", "Missense_Mutation", "p.G12C"), egfr("Missense_Mutation", "p.G719A")], "KRAS mutant"),
    ([("mutation", "KRAS", "Missense_Mutation", "p.G12D"), egfr("In_Frame_Ins", "p.D770_N771insSVD")], "KRAS mutant"),
], ids=case_id)
def test_made_up_egfr_tumors_land_where_option_c_says(q3, calls, group):
    assert group_of(q3, calls) == group


@pytest.fixture(scope="module")
def egfr_cohorts(nsclc_files, nsclc, tmp_path_factory):
    """The made-up NSCLC-shaped files before and after made-up EGFR calls are added: three to patients in the
    reference and T790M to one EGFR common patient."""
    before = build_cohort(outcome_blind(genie_bpc.read(nsclc_files, None, nsclc).tables), nsclc).patients
    ref = sorted(before.index[before["genotype"] == REF])[:3]
    common = sorted(before.index[before["genotype"] == COMMON])[0]
    calls = {ref[0]: ("Missense_Mutation", "p.G719A", FIFTH), ref[1]: ("In_Frame_Ins", "p.D770_N771insSVD", FIFTH),
             ref[2]: ("Missense_Mutation", "p.E746_L747delinsIP", REF), common: ("Missense_Mutation", "p.T790M", COMMON)}
    d = tmp_path_factory.mktemp("nsclc_egfr")
    for p in nsclc_files.iterdir():
        shutil.copy(p, d / p.name)
    reports = read(d / "cancer_panel_test_level_dataset.csv")
    extra = [dict(Hugo_Symbol="EGFR", Variant_Classification=c, HGVSp_Short=v, Tumor_Sample_Barcode=s)
             for pid, (c, v, _) in calls.items() for s in reports.loc[reports["record_id"] == pid, "cpt_genie_sample_id"]]
    maf = read(d / "data_mutations_extended.txt", "\t")
    write(pd.concat([maf, pd.DataFrame(extra)], ignore_index=True), d / "data_mutations_extended.txt", "\t")
    after = build_cohort(outcome_blind(genie_bpc.read(d, None, nsclc).tables), nsclc).patients
    return before, after, calls


def test_made_up_egfr_calls_in_the_files_move_only_their_patients(egfr_cohorts):
    before, after, calls = egfr_cohorts
    for pid, (_, variant, group) in calls.items():
        assert after.loc[pid, "genotype"] == group, (pid, variant)
    others = before.index.difference(list(calls))
    pd.testing.assert_series_equal(after.loc[others, "genotype"], before.loc[others, "genotype"])
