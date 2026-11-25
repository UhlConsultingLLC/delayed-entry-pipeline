"""Question 5, option D with K (decided 2026-11-24): one KRAS group, KRAS mutant (p.G12C, and KRAS missense or
in-frame calls at codons 12, 13 and 61), with KRAS G12C kept as a sensitivity analysis, SA13; the codons stay 12,
13 and 61, and the RAS changes left in the reference are stated. On made-up tumors and on the made-up NSCLC-shaped
files. Every call here is made up."""
import pandas as pd
import pytest

from delayed_entry.analysis import keep_mask
from delayed_entry.settings import load_settings
from test_alk_sensitivity import planned
from test_nsclc_spec import NSCLC_SETTINGS, nsclc, nsclc_files  # noqa: F401 (fixtures)
from test_q3_egfr_rules import case_id, group_of

KRAS, COMMON, REF = "KRAS mutant", "EGFR common", "No recognized driver"
CLASSES = ["Missense_Mutation", "In_Frame_Del", "In_Frame_Ins"]
KRAS_RULES = [
    {"kind": "mutation", "genes": ["KRAS"], "variants": ["p.G12C"]},
    {"kind": "mutation", "genes": ["KRAS"], "variant_classes": CLASSES, "codons": [12, 13, 61]},
]
SA13_RESULT = {"kind": "mutation", "genes": ["KRAS"], "variant_classes": CLASSES, "codons": [12, 13, 61],
               "except_variants": ["p.G12C"]}


def kras(variant, variant_class="Missense_Mutation"):
    return ("mutation", "KRAS", variant_class, variant)


@pytest.fixture(scope="module")
def q5():
    return load_settings(NSCLC_SETTINGS)


def sa13(settings):
    return next(a for a in settings.read("analyses") if a["id"] == "SA13")


def test_one_kras_group_with_two_rules_and_no_except_entry(q5):
    groups = {g["name"]: g for g in q5.read("genotype.groups")}
    assert groups[KRAS]["any_of"] == KRAS_RULES
    assert "KRAS G12C" not in groups and "KRAS other" not in groups
    assert q5.read("smallest_group.if_short")[KRAS] == "describe"
    assert q5.read("genotype.coverage.required")["KRAS"] == [12, 13, 61]       # K: the codons are kept


def test_question_5_is_recorded_and_the_reference_is_stated(q5):
    groups = q5.provenance(["genotype.groups"])[0]
    assert groups["made_by"] == "NSCLC study author"      # every part decided since the fusion rows (2026-11-24)
    assert "Question 5, decided 2026-11-24 (option D, with the codons kept)" in groups["source"]
    assert "Question 5 decided 2026-11-24 (option D, with the codons kept)" in groups["note"]
    assert "TO DECIDE" not in groups["note"] and "Fusion rows decided 2026-11-24 (option C)" in groups["note"]
    assert "NM_033360 / ENST00000256078 (KRAS" in groups["note"]
    assert "KRAS changes outside codons 12, 13 and 61" in groups["note"] and "NRAS and HRAS" in groups["note"]
    assert "five groups" not in q5.provenance(["genotype.reference"])[0]["source"]
    analyses = q5.provenance(["analyses"])[0]
    assert "SA13 (KRAS G12C only) is added by question 5 (decided 2026-11-24, option D)" in analyses["note"]


def test_sa13_is_kras_g12c_only_and_its_except_entry_is_pinned(q5):
    a = sa13(q5)
    assert (a["run"], a["kind"], a["label"], a["if_short"]) == (True, "subset", "KRAS G12C only", "describe")
    assert a["keep"] == {"exclude_from_group": KRAS, "if_result": SA13_RESULT}


@pytest.mark.parametrize("calls, group", [
    ([kras("p.G12C")], KRAS), ([kras("p.G12D")], KRAS), ([kras("p.G13D")], KRAS), ([kras("p.Q61H")], KRAS),
    ([kras("p.G12F")], KRAS), ([kras("p.G12_G13delinsVC")], KRAS), ([kras("p.A11_G12dup", "In_Frame_Ins")], KRAS),
    ([kras("p.G12C", None)], KRAS),                                     # a G12C with a blank class still counts
    # left in the reference, and stated
    ([kras("p.A146T")], REF), ([kras("p.K117N")], REF), ([kras("p.A59T")], REF), ([kras("p.L19F")], REF),
    ([kras("p.Q22K")], REF), ([kras("p.G10_A11dup", "In_Frame_Ins")], REF),
    ([("mutation", "NRAS", "Missense_Mutation", "p.Q61K")], REF),
    ([("mutation", "HRAS", "Missense_Mutation", "p.G12V")], REF),
], ids=case_id)
def test_made_up_ras_tumors_land_where_option_d_with_k_says(q5, calls, group):
    assert group_of(q5, calls) == group


@pytest.mark.parametrize("calls, kept", [
    ([kras("p.G12C")], True),
    ([kras("p.G12C", None)], True),                                     # a blank class
    ([kras("p.G12C"), kras("p.A146T")], True),                          # a second call away from 12, 13 and 61
    ([kras("p.G12C"), kras("p.G12=", "Silent")], True),                 # a silent second call
    ([kras("p.G12C"), kras("p.G13D")], False),
    ([kras("p.G12C"), kras("p.G12V")], False),                          # a G12F written as two records
    ([kras("p.G12C"), kras("p.G12_G13delinsVC")], False),
    ([kras("p.G12D")], False),
    ([kras("p.Q61H")], False),
], ids=case_id)
def test_sa13_keeps_a_kras_mutant_tumor_only_without_another_call_at_12_13_or_61(q5, calls, kept):
    first = pd.DataFrame({"first_test_id": ["T1"], "genotype": [KRAS]}, index=pd.Index(["P1"], name="patient_id"))
    results = pd.DataFrame([("T1", *c) for c in calls], columns=["test_id", "kind", "gene", "variant_class", "variant"])
    assert keep_mask(first, None, results, sa13(q5)["keep"], q5).tolist() == [kept]


def test_sa13_keeps_every_other_group_whole(q5):
    first = pd.DataFrame({"first_test_id": ["T1", "T2"], "genotype": [COMMON, REF]},
                         index=pd.Index(["P1", "P2"], name="patient_id"))
    results = pd.DataFrame([("T1", *kras("p.G12D")), ("T1", "mutation", "EGFR", "Missense_Mutation", "p.L858R")],
                           columns=["test_id", "kind", "gene", "variant_class", "variant"])
    assert keep_mask(first, None, results, sa13(q5)["keep"], q5).tolist() == [True, True]


def test_the_run_plans_sa13_and_leaves_out_only_kras_mutant_patients(nsclc_files, nsclc):
    """On the made-up NSCLC-shaped files, with their panels standing in for the release's."""
    data = planned(nsclc_files, nsclc)
    a = data.analysis("SA13")
    assert a.status == "run"
    sec = data.secondary.set_index("patient_id")
    groups = data.in_analysis().set_index("patient_id")["genotype"]
    left_out = groups.loc[sec.index[~sec[a.member_column].astype(bool)]]
    kept_kras = groups.loc[sec.index[sec[a.member_column].astype(bool)]]
    assert len(left_out) and (left_out == KRAS).all()
    assert (kept_kras == KRAS).any()
