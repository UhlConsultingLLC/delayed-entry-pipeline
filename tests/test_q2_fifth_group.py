"""Question 2, option C (decided 2026-11-23): the fifth group, Other recognized driver, on made-up rows and on the
made-up NSCLC-shaped files. Every call here is made up."""
import pandas as pd
import pytest

from delayed_entry.adapters import genie_bpc
from delayed_entry.cohort import build_cohort, result_matches
from delayed_entry.settings import SettingsError, load_settings
from delayed_entry.tables import outcome_blind
from test_nsclc_spec import NSCLC_SETTINGS, ROS1_PATIENT, RET_PATIENT, OTHER_DRIVER_CALLS, nsclc_files, nsclc  # noqa: F401
import yaml

FIFTH, REF = "Other recognized driver", "No recognized driver"


def test_five_groups_in_order_with_the_new_reference_and_no_sa4():
    s = load_settings(NSCLC_SETTINGS)
    assert [g["name"] for g in s.read("genotype.groups")] == ["EGFR common", "ALK fusion", "KRAS mutant", FIFTH,
                                                             REF]
    assert s.read("genotype.reference") == REF
    assert "SA4" not in {a["id"] for a in s.read("analyses")}
    assert s.read("genotype.described_only") == ["Other recognized driver"]            # question 13, option B
    assert "Other recognized driver" not in s.read("smallest_group.if_short")


def fifth_hits(rows):
    s = load_settings(NSCLC_SETTINGS)
    rules = next(g for g in s.read("genotype.groups") if g["name"] == FIFTH)["any_of"]
    res = pd.DataFrame(rows, columns=["kind", "gene", "variant_class", "variant"])
    hit = pd.Series(False, index=res.index)
    for m in rules:
        hit |= result_matches(res, m, "trimmed")
    return hit.tolist()


@pytest.mark.parametrize("row, counts", [
    (("mutation", "MET", "Splice_Site", "p.X1010_splice"), True),
    (("mutation", "MET", "Splice_Site", "p.X963_splice"), True),
    (("mutation", "MET", "Splice_Site", "p.X1004_splice"), True),
    (("mutation", "MET", "Missense_Mutation", "p.D1010N"), True),
    (("mutation", "MET", "Missense_Mutation", "p.D1010A"), False),
    (("mutation", "MET", "Missense_Mutation", "p.Y1003S"), False),
    (("mutation", "MET", "Intron", None), False),
    (("mutation", "MET", "Splice_Site", "p.X1028_splice"), False),
    (("mutation", "ERBB2", "In_Frame_Ins", "p.Y772_A775dup"), True),
    (("mutation", "ERBB2", "In_Frame_Ins", "p.G776delinsVC"), True),
    (("mutation", "ERBB2", "Missense_Mutation", "p.S310F"), True),
    (("mutation", "ERBB2", "Missense_Mutation", "p.S310A"), False),
    (("mutation", "ERBB2", "Missense_Mutation", "p.V777L"), True),
    (("mutation", "ERBB2", "Missense_Mutation", "p.V777M"), False),
    (("mutation", "ERBB2", "Missense_Mutation", "p.I655V"), False),
    (("mutation", "ERBB2", "Amplification", None), False),
    (("mutation", "BRAF", "Missense_Mutation", "p.V600E"), True),
    (("mutation", "BRAF", "Missense_Mutation", "p.V600K"), False),
    (("fusion", "RET", None, "KIF5B-RET fusion"), True),
    (("fusion", "ROS1", None, "CD74-ROS1 fusion"), True),
    (("fusion", "NTRK1", None, "TPM3-NTRK1 fusion"), False),
])
def test_the_fifth_groups_rules_on_made_up_rows(row, counts):
    assert fifth_hits([row]) == [counts]


def test_made_up_other_drivers_land_in_the_fifth_group(nsclc_files, nsclc):
    c = build_cohort(outcome_blind(genie_bpc.read(nsclc_files, None, nsclc).tables), nsclc).patients
    for pid in [ROS1_PATIENT, RET_PATIENT, *OTHER_DRIVER_CALLS]:
        assert c.loc[pid, "genotype"] == FIFTH, pid


def test_naming_the_fifth_group_in_both_or_neither_place_is_refused(tmp_path):
    data = yaml.safe_load(NSCLC_SETTINGS.read_text(encoding="utf-8"))
    short = data["smallest_group"]["if_short"]["value"]
    both, neither = dict(short), {k: v for k, v in short.items() if k != FIFTH}
    both[FIFTH] = "describe"
    for described, if_short in (([FIFTH], both), ([], neither)):
        data["genotype"]["described_only"]["value"] = described
        data["smallest_group"]["if_short"]["value"] = if_short
        p = tmp_path / "q13.yaml"
        p.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
        with pytest.raises(SettingsError, match="if_short"):
            load_settings(p)
