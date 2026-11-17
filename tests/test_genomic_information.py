"""Panel sizes from a GENIE genomic_information file, on made-up regions written at test time: the merging, the
settings' rule, the adapter (with the made-up GENIE BPC fixture), the burden fallback it enables, and the check
after download."""
import shutil

import pandas as pd
import pytest
import yaml

from delayed_entry.adapters import genie_bpc
from delayed_entry.adapters.genie_bpc import check_download as cd
from delayed_entry.adapters.genie_bpc.genomic_information import COLUMNS, merged_length, panel_sizes
from delayed_entry.screen import BURDEN_HIGH, BURDEN_LOW, screen
from delayed_entry.settings import SettingsError, load_settings
from delayed_entry.tables import outcome_blind
from helpers import CRC_SETTINGS, FIXTURE

RULE = {"from": "genomic_information", "feature_types": ["exon"], "include_in_panel": "only_marked"}
PANELS = ["DFCI-PANEL-1", "DFCI-PANEL-2", "MSK-PANEL-1", "MSK-PANEL-2", "VICC-PANEL-1", "VICC-PANEL-2"]


def regions(rows) -> pd.DataFrame:
    """Made-up rows: (panel, chromosome, start, end, feature type, includeInPanel)."""
    return pd.DataFrame([dict(zip(COLUMNS, map(str, r))) for r in rows])


def test_overlapping_and_repeated_regions_are_counted_once():
    assert merged_length([0, 50, 50, 300], [100, 150, 150, 400]) == 250        # 0-150 and 300-400
    assert merged_length([0, 100], [100, 200]) == 200                           # touching regions
    assert merged_length([], []) == 0


def test_the_rule_says_which_regions_count():
    df = regions([("A", "1", 0, 100_000, "exon", "True"), ("A", "1", 50_000, 150_000, "exon", "True"),
                  ("A", "1", 50_000, 150_000, "exon", "True"),                  # listed again: one per exon it touches
                  ("A", "2", 0, 50_000, "exon", "True"), ("A", "1", 500_000, 600_000, "intron", "True"),
                  ("A", "3", 0, 400_000, "exon", "False"), ("B", "7", 0, 120_000, "exon", "TRUE"),
                  ("B", "7", "x", 10, "exon", "True")])
    sizes, unmeasured = panel_sizes(df, RULE)
    assert dict(zip(sizes["panel_id"], sizes["coding_megabases"])) == {"A": 0.2, "B": 0.12} and unmeasured == 1
    wider = dict(panel_sizes(df, {**RULE, "feature_types": ["exon", "intron"], "include_in_panel": "all"})[0].values)
    assert wider["A"] == pytest.approx(0.2 + 0.1 + 0.4)


# ---- The adapter ----
def settings_file(tmp_path, rule=RULE, **values):
    data = yaml.safe_load(CRC_SETTINGS.read_text(encoding="utf-8"))
    if rule is not None:
        data["source"]["genie_bpc"]["panel_sizes"] = {"value": rule, "made_by": "a test", "source": "made up"}
    for path, v in values.items():
        node = data
        for part in path.split("."):
            node = node[part]
        node["value"] = v
    path = tmp_path / "settings.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return load_settings(path)


def release(tmp_path, megabases=0.15, panels=PANELS, name="genomic_information.txt"):
    data = shutil.copytree(FIXTURE, tmp_path / "data")
    if panels:
        size = int(megabases * 1_000_000)
        rows = [(p, "1", 0, size, "exon", "True") for p in panels]
        rows += [(p, "1", size // 2, size, "exon", "True") for p in panels]            # an overlap, counted once
        rows += [(p, "2", 0, 9_000_000, "intron", "True") for p in panels]             # not exonic: left out
        regions(rows).to_csv(data / name, sep="\t", index=False)
    return data


def test_the_adapter_reads_panel_sizes_by_the_settings_rule(tmp_path):
    res = genie_bpc.read(release(tmp_path), None, settings_file(tmp_path))
    sizes = dict(zip(res.tables.panel_sizes["panel_id"], res.tables.panel_sizes["coding_megabases"]))
    assert sizes == {p: 0.15 for p in PANELS}
    assert res.profile.panel_sizes == "given" and "genomic_information.txt" in res.profile.sources["panel_sizes"]
    assert any(p.name == "genomic_information.txt" for p in res.input_files)
    assert any("panel sizes for" in n for n in res.read_notes)


def test_with_sizes_the_burden_fallback_compares_high_with_low(tmp_path):
    s = settings_file(tmp_path, **{"biomarker.unknown_max_share": 0.25})      # the fixture's S1 is 26% unknown
    res = genie_bpc.read(release(tmp_path), None, s)
    rep = screen(outcome_blind(res.tables), s, res.profile, res.source_steps, res.input_files)
    assert rep.decisions["S1"] == {f"{BURDEN_HIGH} vs {BURDEN_LOW}": "compare"}


def test_a_file_the_settings_dont_read_is_noted_and_left(tmp_path):
    res = genie_bpc.read(release(tmp_path), None, settings_file(tmp_path, rule=None))
    assert res.tables.panel_sizes is None and res.profile.panel_sizes == "not_given"
    assert any("not read: the settings give no rule" in n for n in res.read_notes)


@pytest.mark.parametrize("make, problem", [
    (lambda d: (d / "genomic_information.txt").unlink(), "the genomic folder has none"),
    (lambda d: shutil.copy(d / "genomic_information.txt", d / "genomic_information_2.txt"), "more than one"),
    (lambda d: (d / "genomic_information.txt").write_text("SEQ_ASSAY_ID\tChromosome\nA\t1\n", encoding="utf-8"),
     "no column Start_Position"),
])
def test_a_file_that_cant_give_sizes_stops_the_adapter(tmp_path, make, problem):
    data = release(tmp_path)
    make(data)
    with pytest.raises(ValueError, match=problem):
        genie_bpc.read(data, None, settings_file(tmp_path))


def test_a_panel_the_file_doesnt_list_has_no_size(tmp_path):
    res = genie_bpc.read(release(tmp_path, panels=PANELS[:2]), None, settings_file(tmp_path))
    assert set(res.tables.panel_sizes["panel_id"]) == set(PANELS[:2])


# ---- The settings and the check after download ----
@pytest.mark.parametrize("rule, problem", [
    ({**RULE, "feature_types": ["coding"]}, "panel_sizes"),
    ({**RULE, "include_in_panel": "some"}, "panel_sizes"),
])
def test_a_malformed_rule_is_refused(tmp_path, rule, problem):
    with pytest.raises(SettingsError, match=problem):
        settings_file(tmp_path, rule=rule)


def test_the_check_after_download_judges_the_file_against_the_rule(tmp_path):
    data = release(tmp_path)
    status = lambda checks: {c.name: c.status for c in checks}["The genomic_information file, for panel sizes"]
    assert status(cd.check_genomic_information(data, RULE)) == cd.PASS
    assert status(cd.check_genomic_information(data, None)) == cd.REVIEW
    empty = tmp_path / "empty"
    empty.mkdir()
    assert status(cd.check_genomic_information(empty, RULE)) == cd.FAIL
    assert status(cd.check_genomic_information(empty, None)) == cd.PASS
