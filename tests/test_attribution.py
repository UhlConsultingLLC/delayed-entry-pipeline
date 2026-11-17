"""The attribution a source's terms require: an optional capability-profile item, filled in for GENIE BPC with the
release the settings name, printed in the run report, and warned about while a part of it is blank. Made-up values."""
import dataclasses

import pandas as pd
import pytest

from delayed_entry.adapters import genie_bpc
from delayed_entry.profile import ProfileError, profile_from_dict
from delayed_entry.run_report import build_report
from helpers import FIXTURE, crc_settings

BASE = {"entry_event": "report", "panel_coverage": "none", "panel_sizes": "not_given", "death_capture": "linked",
        "time_format": "intervals", "privacy_min_count": 11, "histology_vocabulary": "made up"}
FULL = {"citation": "A. Consortium. A made-up paper. J Made Up. 2020;1:1-2.", "version": "release 9",
        "acknowledgment": "We thank the made-up consortium."}


def test_a_profile_without_attribution_is_as_before():
    p = profile_from_dict(BASE, adapter="flat_file")
    assert p.attribution is None and "attribution" not in p.to_dict() and p.attribution_gaps() == []


def test_an_attribution_is_kept_and_its_blank_parts_named():
    p = profile_from_dict({**BASE, "attribution": {**FULL, "acknowledgment": None}}, adapter="flat_file")
    assert p.to_dict()["attribution"]["citation"] == FULL["citation"]
    assert p.attribution_gaps() == ["acknowledgment"]
    assert profile_from_dict({**BASE, "attribution": FULL}, adapter="flat_file").attribution_gaps() == []


@pytest.mark.parametrize("attribution, problem", [
    ({"citation": "x"}, "attribution must be null, or"),
    ({**FULL, "citation": ""}, "attribution.citation must be given"),
    ({**FULL, "version": 2020}, "attribution.version must be text or null"),
    ("cite us", "attribution must be null, or"),
])
def test_a_malformed_attribution_is_refused(attribution, problem):
    with pytest.raises(ProfileError, match=problem.replace("(", r"\(").replace(".", r"\.")):
        profile_from_dict({**BASE, "attribution": attribution}, adapter="flat_file")


def test_genie_bpc_carries_the_consortium_citation_and_the_release_read():
    bare = genie_bpc.profile()
    assert "Cancer Discov. 2017;7(8):818-831" in bare.attribution["citation"] and bare.attribution["version"] is None
    res = genie_bpc.read(FIXTURE, None, crc_settings())
    assert res.profile.attribution["version"] == "GENIE BPC crc_v2.0 (Synapse syn30991637)"
    assert "version" not in res.profile.attribution_gaps()


def test_a_version_the_profile_gives_is_kept():
    given = dataclasses.replace(genie_bpc.profile(), attribution={**genie_bpc.profile().attribution, "version": "mine"})
    assert genie_bpc.with_release(given, "crc_v2.0").attribution["version"] == "mine"


def info(profile: dict) -> dict:
    """The least a run report needs, all made up."""
    return {"study": "A made-up study", "registered_plan": "not registered", "settings_file": "s.yaml",
            "settings_sha256": "0" * 64, "plan_version": 1, "locked_at": "a time", "locked_by": "A. Reviewer",
            "plan_files": {}, "tool": "delayed-entry-pipeline", "run_at": "a time", "r_versions": {}, "resamples": 10,
            "seed": 1, "inputs": [], "profile": profile,
            "screen": {"ready": [], "missing": [], "deviations_required": [], "warnings": []},
            "checks": [], "analyses": [], "min_count": 11, "horizon": 3, "interval_level": 0.95, "files": {}}


def report(profile: dict) -> str:
    empty = pd.DataFrame(columns=["analysis", "pass", "quantity", "group", "estimate", "lower", "upper", "label"])
    return build_report(info(profile), {"estimates.csv": empty, "bias.csv": empty}, [], []).markdown()


def test_the_run_report_prints_the_attribution_and_says_what_is_blank():
    md = report({"adapter": "made_up", "attribution": {**FULL, "acknowledgment": None}})
    section = md.split("## Data source and attribution")[1].split("## ")[0]
    assert FULL["citation"] in section and "Dataset version: release 9" in section
    assert "Acknowledgment: NOT FILLED IN" in section
    assert "| attribution |" not in md                     # not repeated in the profile's table
    done = report({"adapter": "made_up", "attribution": FULL})
    assert FULL["acknowledgment"] in done and "NOT FILLED IN" not in done


def test_a_source_without_attribution_gets_no_section():
    assert "Data source and attribution" not in report({"adapter": "made_up"})
