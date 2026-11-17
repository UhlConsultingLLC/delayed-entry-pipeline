"""Locking the plan (design plan section 7)."""
import json
from datetime import datetime, timezone

import pytest
import yaml

from delayed_entry.lock import (LockError, lock_plan, lock_settings, require_lock, settings_lock_matches, verify_lock)
from delayed_entry.screen import screen
from delayed_entry.tables import outcome_blind
from helpers import CRC_SETTINGS, bpc_profile, crc_settings, make_tables, settings_with

WHEN = datetime(2026, 11, 13, 12, 0, tzinfo=timezone.utc)
GROUPS = {"RAS-mutant": 60, "BRAF V600E": 25, "Wild type": 70}


def made_up_inputs(folder):
    folder.mkdir(exist_ok=True)
    inputs = []
    for name in ("a.csv", "b.csv"):
        (folder / name).write_text(f"made-up,{name}\n", encoding="utf-8")
        inputs.append(folder / name)
    return inputs


@pytest.fixture
def setup(tmp_path):
    s = crc_settings()
    inputs = made_up_inputs(tmp_path / "data")
    rep = screen(outcome_blind(make_tables(GROUPS)), s, bpc_profile(), inputs=inputs)
    return s, rep, inputs, tmp_path / "plan"


def test_the_lock_writes_the_plan_twice_and_a_manifest(setup):
    s, rep, inputs, out = setup
    res = lock_plan(rep, s, inputs, out, locked_by="A. Reviewer", now=WHEN)
    assert res.version == 1 and {p.name for p in res.folder.iterdir()} == {"plan.md", "plan.yaml", "settings.yaml", "manifest.json"}
    m = json.loads(res.manifest.read_text(encoding="utf-8"))
    assert m["locked_at_utc"] == "2026-11-13T12:00:00+00:00" and m["locked_by"] == "A. Reviewer"
    assert set(m["plan_files"]) == {"plan.md", "plan.yaml", "settings.yaml"}
    assert all(len(d) == 64 for d in m["plan_files"].values())
    assert m["settings"]["sha256"] == s.sha256 and m["settings"]["resolved_sha256"] == s.resolved_sha256
    assert len(m["defaults"]["sha256"]) == 64 and len(m["inputs"]) == 2
    assert (res.folder / "settings.yaml").read_text(encoding="utf-8") == CRC_SETTINGS.read_text(encoding="utf-8")


def test_the_plan_holds_the_settings_the_profile_and_the_three_lists(setup):
    s, rep, inputs, out = setup
    res = lock_plan(rep, s, inputs, out, locked_by="A. Reviewer", read_notes=["a note from reading"], now=WHEN)
    plan = yaml.safe_load((res.folder / "plan.yaml").read_text(encoding="utf-8"))
    assert {p["choice"] for p in plan["settings"]} == {p["choice"] for p in s.provenance()}
    assert plan["capability_profile"]["entry_event"] == "report"
    screen_part = plan["readiness_screen"]
    assert {"ready", "missing", "deviations_required"} <= set(screen_part)
    assert screen_part["decisions"]["primary"]["BRAF V600E vs Wild type"] == "descriptive only"
    assert {i["file"] for i in screen_part["inputs"]} == {"a.csv", "b.csv"}
    assert plan["reading_the_source"] == ["a note from reading"]
    md = (res.folder / "plan.md").read_text(encoding="utf-8")
    assert "Reviewed and locked by: A. Reviewer." in md and "## Settings" in md and "### Deviations required" in md


def test_nothing_changed_means_the_lock_holds(setup):
    s, rep, inputs, out = setup
    res = lock_plan(rep, s, inputs, out, locked_by="A. Reviewer")
    assert verify_lock(res.manifest, s, inputs) == []
    require_lock(res.manifest, s, inputs)


@pytest.mark.parametrize("newline", [b"\n", b"\r\n"], ids=["unix", "windows"])
def test_the_lock_holds_whatever_line_endings_the_settings_file_has(tmp_path, newline):
    """The settings' hash is the file's own bytes, so a file saved with Windows line endings verifies, and the
    plan's copy of it is the same bytes."""
    from delayed_entry.hashing import sha256_file
    from delayed_entry.settings import load_settings
    path = tmp_path / "settings.yaml"
    path.write_bytes(CRC_SETTINGS.read_bytes().replace(b"\r\n", b"\n").replace(b"\n", newline))
    s = load_settings(path)
    assert s.sha256 == sha256_file(path)
    assert s.resolved_sha256 == crc_settings().resolved_sha256         # the same values either way
    inputs = made_up_inputs(tmp_path / "data")
    rep = screen(outcome_blind(make_tables(GROUPS)), s, bpc_profile(), inputs=inputs)
    res = lock_plan(rep, s, inputs, tmp_path / "plan", locked_by="A. Reviewer")
    assert verify_lock(res.manifest) == [] and verify_lock(res.manifest, load_settings(path), inputs) == []
    require_lock(res.manifest, load_settings(path), inputs)
    assert (res.folder / "settings.yaml").read_bytes() == path.read_bytes()
    lock_settings(path, tmp_path / "settings_lock", "A. Reviewer")
    assert settings_lock_matches(tmp_path / "settings_lock" / "settings_lock.json", path)


def test_changing_only_the_line_endings_after_the_lock_is_a_change(tmp_path):
    """The lock fingerprints the file as it is: re-saving it with other line endings breaks the lock."""
    from delayed_entry.settings import load_settings
    path = tmp_path / "settings.yaml"
    path.write_bytes(CRC_SETTINGS.read_bytes().replace(b"\r\n", b"\n"))
    inputs = made_up_inputs(tmp_path / "data")
    s = load_settings(path)
    res = lock_plan(screen(outcome_blind(make_tables(GROUPS)), s, bpc_profile(), inputs=inputs), s, inputs,
                    tmp_path / "plan", locked_by="A. Reviewer")
    path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
    assert verify_lock(res.manifest) == ["the settings file has changed since the lock"]


def test_a_changed_input_breaks_the_lock(setup):
    s, rep, inputs, out = setup
    res = lock_plan(rep, s, inputs, out, locked_by="A. Reviewer")
    inputs[0].write_text("made-up,2\n", encoding="utf-8")
    assert verify_lock(res.manifest) == ["input a.csv has changed since the lock"]
    with pytest.raises(LockError):
        require_lock(res.manifest, s, inputs)


def test_a_changed_plan_file_breaks_the_lock(setup):
    s, rep, inputs, out = setup
    res = lock_plan(rep, s, inputs, out, locked_by="A. Reviewer")
    p = res.folder / "plan.md"
    p.write_text(p.read_text(encoding="utf-8") + "\nedited\n", encoding="utf-8")
    assert verify_lock(res.manifest) == ["plan.md has changed since the lock"]


def test_an_edited_settings_file_breaks_the_lock_even_unasked(tmp_path):
    """verify_lock hashes the settings file again where the manifest recorded it."""
    path = tmp_path / "study.yaml"
    path.write_bytes(CRC_SETTINGS.read_bytes())
    from delayed_entry.settings import load_settings
    s = load_settings(path)
    inputs = made_up_inputs(tmp_path / "data")
    res = lock_plan(screen(outcome_blind(make_tables(GROUPS)), s, bpc_profile(), inputs=inputs), s, inputs,
                    tmp_path / "plan", locked_by="A. Reviewer")
    path.write_text(path.read_text(encoding="utf-8").replace("value: 730", "value: 731"), encoding="utf-8")
    assert verify_lock(res.manifest) == ["the settings file has changed since the lock"]


def test_a_changed_defaults_file_breaks_the_lock(setup, tmp_path):
    s, rep, inputs, out = setup
    res = lock_plan(rep, s, inputs, out, locked_by="A. Reviewer")
    m = json.loads(res.manifest.read_text(encoding="utf-8"))
    copy = tmp_path / "defaults.yml"
    copy.write_bytes(open(m["defaults"]["path"], "rb").read())
    m["defaults"]["path"] = str(copy)                 # point the record at a copy, then change the copy
    res.manifest.write_text(json.dumps(m), encoding="utf-8")
    copy.write_text(copy.read_text(encoding="utf-8") + "\n# changed\n", encoding="utf-8")
    assert verify_lock(res.manifest) == ["the defaults file has changed since the lock"]


def test_other_settings_or_other_inputs_break_the_lock(setup, tmp_path):
    s, rep, inputs, out = setup
    res = lock_plan(rep, s, inputs, out, locked_by="A. Reviewer")
    other = settings_with(tmp_path, {"smallest_group.patients": 31})
    assert verify_lock(res.manifest, settings=other) == ["the settings are not the ones locked"]
    assert verify_lock(res.manifest, inputs=inputs[:1])[0].startswith("the inputs are not the ones locked")


def test_a_person_must_review_the_plan(setup):
    s, rep, inputs, out = setup
    with pytest.raises(LockError, match="a person reviews the plan"):
        lock_plan(rep, s, inputs, out, locked_by=" ")


def test_the_report_must_come_from_the_same_settings(setup, tmp_path):
    _, rep, inputs, out = setup
    with pytest.raises(LockError, match="other settings"):
        lock_plan(rep, settings_with(tmp_path, {"smallest_group.patients": 31}), inputs, out, locked_by="A. Reviewer")


def test_the_report_must_come_from_the_same_inputs(setup, tmp_path):
    s, rep, _, out = setup
    others = made_up_inputs(tmp_path / "other")
    others[0].write_text("different\n", encoding="utf-8")
    with pytest.raises(LockError, match="other input files"):
        lock_plan(rep, s, others, out, locked_by="A. Reviewer")
    unrecorded = screen(outcome_blind(make_tables(GROUPS)), s, bpc_profile())
    with pytest.raises(LockError, match="doesn't record the files"):
        lock_plan(unrecorded, s, others, out, locked_by="A. Reviewer")


def test_a_plan_without_inputs_is_refused(setup):
    s, rep, _, out = setup
    with pytest.raises(LockError, match="give the input files"):
        lock_plan(rep, s, [], out, locked_by="A. Reviewer")


def test_a_privacy_minimum_below_the_agreement_is_never_locked(tmp_path):
    s = settings_with(tmp_path, {"privacy.min_count": 5})
    inputs = made_up_inputs(tmp_path / "data")
    rep = screen(outcome_blind(make_tables(GROUPS)), s, bpc_profile(), inputs=inputs)
    with pytest.raises(LockError, match="privacy minimum"):
        lock_plan(rep, s, inputs, tmp_path / "plan", locked_by="A. Reviewer")


def test_a_locked_plan_changes_only_by_a_new_version_with_a_reason(setup):
    s, rep, inputs, out = setup
    lock_plan(rep, s, inputs, out, locked_by="A. Reviewer", now=WHEN)
    with pytest.raises(LockError, match="needs a reason"):
        lock_plan(rep, s, inputs, out, locked_by="A. Reviewer")
    res = lock_plan(rep, s, inputs, out, locked_by="B. Reviewer", reason="a corrected input file", now=WHEN)
    assert res.version == 2 and (out / "plan_v1" / "manifest.json").exists()
    log = (out / "deviations_log.md").read_text(encoding="utf-8")
    assert "| 2 | 2026-11-13T12:00:00+00:00 | B. Reviewer | a corrected input file |" in log


def test_a_failed_lock_leaves_nothing_behind(setup):
    s, rep, inputs, out = setup
    inputs[1].unlink()
    with pytest.raises(LockError, match="input file not found"):
        lock_plan(rep, s, inputs, out, locked_by="A. Reviewer")
    assert not out.exists() or list(out.iterdir()) == []
    inputs[1].write_text("made-up,b.csv\n", encoding="utf-8")
    assert lock_plan(rep, s, inputs, out, locked_by="A. Reviewer").version == 1     # no reason asked: nothing was locked
    assert not (out / "deviations_log.md").exists()


def test_the_run_must_name_its_settings_and_inputs(setup):
    s, rep, inputs, out = setup
    res = lock_plan(rep, s, inputs, out, locked_by="A. Reviewer")
    with pytest.raises(TypeError):
        require_lock(res.manifest)                                  # noqa: the run can't skip them


def test_the_settings_can_be_locked_before_any_download(tmp_path):
    copy = tmp_path / "crc.yaml"
    copy.write_bytes(CRC_SETTINGS.read_bytes())
    lock = lock_settings(copy, tmp_path / "lock", locked_by="A. Reviewer", now=WHEN)
    assert settings_lock_matches(lock, copy)
    with pytest.raises(LockError, match="already exists"):
        lock_settings(copy, tmp_path / "lock", locked_by="A. Reviewer")
    copy.write_text(copy.read_text(encoding="utf-8").replace("value: 730", "value: 731"), encoding="utf-8")
    assert not settings_lock_matches(lock, copy)
