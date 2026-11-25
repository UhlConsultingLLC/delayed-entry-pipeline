"""The fusion_detecting_panel subset rule and the NSCLC ALK sensitivity analysis (SA12), on made-up data: the
NSCLC-shaped files the NSCLC tests build from the made-up GENIE BPC fixture.

In 2014-2017 ALK fusions were often found by FISH or IHC rather than DNA panels, and a DNA panel detects them
only if it tiles ALK's introns, so an ALK-positive tumor can sit unseen in the "No recognized driver" reference. SA12
keeps the reference only where the entry test's panel is listed as detecting ALK fusions; the primary analysis
keeps everyone."""
import pandas as pd
import pytest
import yaml

from delayed_entry.adapters import genie_bpc
from delayed_entry.adapters.genie_bpc import check_download as cd
from delayed_entry.cohort import fusion_panel_keep
from delayed_entry.screen import NOT_RUN, screen
from delayed_entry.settings import SettingsError, load_settings
from delayed_entry.tables import outcome_blind
from test_nsclc_spec import MADE_UP_FUSION_PANELS, NSCLC_SETTINGS, nsclc_files  # noqa: F401 (a fixture)

REF = "No recognized driver"


def settings_with_sa12(tmp_path, panels=None, groups="reference_only", fusion_panels=None):
    """The NSCLC settings made over into a study whose primary keeps every patient (no list of fusion-detecting
    panels, unless one is given) and which runs SA12 with the given panels: made-up values. The NSCLC study itself
    chose otherwise (question 6, option A), so its own SA12 is switched off."""
    data = yaml.safe_load(NSCLC_SETTINGS.read_text(encoding="utf-8"))
    sa12 = next(a for a in data["analyses"]["value"] if a["id"] == "SA12")
    sa12["run"] = True
    sa12["keep"] = {"fusion_detecting_panel": {"gene": "ALK", "panels": list(panels or [])}, "groups": groups}
    data["genotype"]["coverage"]["fusion_panels"]["value"] = {} if fusion_panels is None else fusion_panels
    path = tmp_path / "nsclc_sa12.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False, allow_unicode=True), encoding="utf-8")
    return load_settings(path)


# ---- The rule itself ----
def made_up_patients():
    return pd.DataFrame({"genotype": [REF, REF, "ALK fusion", "EGFR common", REF, None],
                         "first_panel_id": ["P-ALK", "P-OTHER", "P-ALK", "P-OTHER", None, "P-OTHER"]})


def test_reference_only_restricts_the_reference_and_keeps_every_other_group():
    keep = {"fusion_detecting_panel": {"gene": "ALK", "panels": ["P-ALK"]}, "groups": "reference_only"}
    assert fusion_panel_keep(made_up_patients(), keep, REF).tolist() == [True, False, True, True, False, True]


def test_all_groups_restricts_every_group():
    keep = {"fusion_detecting_panel": {"gene": "ALK", "panels": ["P-ALK"]}, "groups": "all_groups"}
    assert fusion_panel_keep(made_up_patients(), keep, REF).tolist() == [True, False, True, False, False, False]


# ---- The settings ----
def test_the_nsclc_settings_switch_sa12_off_under_option_a():
    s = load_settings(NSCLC_SETTINGS)
    sa12 = next(a for a in s.read("analyses") if a["id"] == "SA12")
    assert sa12["run"] is False and sa12["keep"]["groups"] == "reference_only"
    # every assignable patient's first panel is on the list already, so SA12 would repeat the primary
    assert sa12["keep"]["fusion_detecting_panel"]["panels"] == s.read("genotype.coverage.fusion_panels")["ALK"]


def test_one_list_of_alk_panels_must_serve_both_rules(tmp_path):
    with pytest.raises(SettingsError, match="SA12: its panels for ALK differ from genotype.coverage.fusion_panels"):
        settings_with_sa12(tmp_path, ["DFCI-PANEL-1"], fusion_panels={"ALK": ["DFCI-PANEL-1", "MSK-PANEL-1"]})
    settings_with_sa12(tmp_path, ["DFCI-PANEL-1"], fusion_panels={"ALK": ["DFCI-PANEL-1"]})      # the same list: fine


def test_a_gene_whose_fusions_no_group_counts_is_refused(tmp_path):
    data = yaml.safe_load(NSCLC_SETTINGS.read_text(encoding="utf-8"))
    next(a for a in data["analyses"]["value"] if a["id"] == "SA12")["keep"]["fusion_detecting_panel"]["gene"] = "BRAF"
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    with pytest.raises(SettingsError, match="no genotype group counts BRAF's fusions"):
        load_settings(path)


# ---- The screen and the check after download, on the NSCLC-shaped files ----
def screened(nsclc_files, s):
    res = genie_bpc.read(nsclc_files, None, s)
    return res, screen(outcome_blind(res.tables), s, res.profile, res.source_steps, res.input_files)


def test_with_no_panels_listed_sa12_does_not_run_and_says_why(nsclc_files, tmp_path):
    _, rep = screened(nsclc_files, settings_with_sa12(tmp_path, []))
    assert rep.decisions["SA12"] == {"all comparisons": NOT_RUN}
    assert rep.find(topic="fusion-detecting panels")


def test_sa12_keeps_the_reference_only_on_the_listed_panels(nsclc_files, tmp_path):
    s = settings_with_sa12(tmp_path, MADE_UP_FUSION_PANELS)          # every panel but VICC-PANEL-2
    res, rep = screened(nsclc_files, s)
    assert "SA12" in rep.decisions and rep.decisions["SA12"] != {"all comparisons": NOT_RUN}
    assert rep.decisions["primary"]                                   # the primary is unchanged by SA12
    from delayed_entry.cohort import build_cohort
    c = build_cohort(outcome_blind(res.tables), s)
    e = c.eligible
    kept = e[fusion_panel_keep(e, {"fusion_detecting_panel": {"gene": "ALK", "panels": MADE_UP_FUSION_PANELS},
                                   "groups": "reference_only"}, REF)]
    dropped = e.loc[~e.index.isin(kept.index)]
    assert len(dropped) and (dropped["genotype"] == REF).all() and (dropped["first_panel_id"] == "VICC-PANEL-2").all()
    assert (kept["genotype"] != REF).sum() == (e["genotype"] != REF).sum()       # the other groups keep everyone


def test_the_check_after_download_names_an_sa12_panel_with_no_panel_file(nsclc_files, tmp_path):
    s = settings_with_sa12(tmp_path, MADE_UP_FUSION_PANELS + ["DFCI-PANEL-9"])
    checks = {c.name: c for c in cd.run_checks(nsclc_files, nsclc_files, tmp_path, s.path, rehearsal=True)}
    named = checks["Every panel SA12 lists as detecting fusions has a panel file"]
    assert named.status == cd.FAIL and named.detail.endswith("DFCI-PANEL-9")


# ---- The run ----
def planned(nsclc_files, s):
    from delayed_entry.analysis import build_analysis
    res = genie_bpc.read(nsclc_files, None, s)
    return build_analysis(res.tables, s, res.profile, res.source_steps)


def test_the_run_leaves_sa12_out_while_no_panel_is_listed(nsclc_files, tmp_path):
    a = planned(nsclc_files, settings_with_sa12(tmp_path, [])).analysis("SA12")
    assert a.status == NOT_RUN and "no panel is listed as detecting ALK fusions" in a.why


def test_the_run_keeps_the_reference_only_on_the_listed_panels(nsclc_files, tmp_path):
    data = planned(nsclc_files, settings_with_sa12(tmp_path, MADE_UP_FUSION_PANELS))
    a = data.analysis("SA12")
    assert a.status == "run"
    sec = data.secondary.set_index("patient_id")
    groups = data.in_analysis().set_index("patient_id")["genotype"]
    left_out = groups.loc[sec.index[~sec[a.member_column].astype(bool)]]
    assert len(left_out) and (left_out == REF).all()
