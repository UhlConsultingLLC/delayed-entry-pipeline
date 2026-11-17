"""The download script's guard: patient-level data only after the study's settings are locked. No network:
nothing here downloads anything."""
import importlib.util

import pytest

from delayed_entry.lock import lock_settings
from helpers import CRC_SETTINGS, ROOT

spec = importlib.util.spec_from_file_location("download_synapse", ROOT / "scripts" / "download_synapse.py")
dl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dl)


@pytest.fixture
def crc_lock(tmp_path):
    return lock_settings(CRC_SETTINGS, tmp_path / "lock", locked_by="A. Reviewer")


def test_every_release_lists_its_three_parts_and_which_hold_patient_data():
    for name, r in dl.releases().items():
        assert set(r["parts"]) == {"documentation", "clinical", "genomic"}, name
        assert [p for p, s in r["parts"].items() if s["patient_data"]] == ["clinical", "genomic"], name
        assert all(s["folder"].startswith("syn") for s in r["parts"].values())


def test_the_crc_release_matches_the_folders_the_adapter_checks():
    import yaml
    expected = yaml.safe_load((ROOT / "src/delayed_entry/adapters/genie_bpc/cohorts/crc_v2.0/expected_release.yml")
                              .read_text(encoding="utf-8"))
    parts = dl.releases()["crc_v2.0"]["parts"]
    assert {k: parts[k]["folder"] for k in ("clinical", "genomic")} == expected["folders"]


def test_documentation_needs_no_lock():
    assert dl.check_allowed("nsclc_v2.0", "documentation", None, None) == "syn29288688"


@pytest.mark.parametrize("part", ["clinical", "genomic"])
def test_patient_data_is_refused_without_locked_settings(part):
    with pytest.raises(dl.Refused, match="patient-level data"):
        dl.check_allowed("nsclc_v2.0", part, None, None)


def test_patient_data_is_allowed_with_the_releases_own_locked_settings(crc_lock):
    assert dl.check_allowed("crc_v2.0", "clinical", CRC_SETTINGS, crc_lock) == "syn39802567"


def test_settings_for_another_release_unlock_nothing(crc_lock):
    with pytest.raises(dl.Refused, match="not for nsclc_v2.0"):
        dl.check_allowed("nsclc_v2.0", "clinical", CRC_SETTINGS, crc_lock)


def test_settings_changed_since_the_lock_are_refused(tmp_path):
    copy = tmp_path / "crc.yaml"
    copy.write_bytes(CRC_SETTINGS.read_bytes())
    lock = lock_settings(copy, tmp_path / "lock", locked_by="A. Reviewer")
    copy.write_text(copy.read_text(encoding="utf-8").replace("value: 730", "value: 731"), encoding="utf-8")
    with pytest.raises(dl.Refused, match="changed since they were locked"):
        dl.check_allowed("crc_v2.0", "genomic", copy, lock)


def test_settings_that_do_not_load_are_refused(tmp_path, crc_lock):
    bad = tmp_path / "bad.yaml"
    bad.write_text("settings_format: 1\n", encoding="utf-8")
    with pytest.raises(dl.Refused, match="doesn't load"):
        dl.check_allowed("crc_v2.0", "clinical", bad, crc_lock)


def test_unknown_releases_and_parts_are_named():
    with pytest.raises(dl.Refused, match="unknown release 'brca_v1.0'"):
        dl.check_allowed("brca_v1.0", "documentation", None, None)
    with pytest.raises(dl.Refused, match="unknown part 'imaging'"):
        dl.check_allowed("nsclc_v2.0", "imaging", None, None)


def test_the_data_folder_must_be_local_and_outside_the_package(tmp_path):
    assert dl.destination(tmp_path / "bpc_nsclc_2.0", "documentation") == (tmp_path / "bpc_nsclc_2.0" / "documentation").resolve()
    with pytest.raises(dl.Refused, match="cloud-synced"):
        dl.destination(tmp_path / "OneDrive" / "data", "documentation")
    with pytest.raises(dl.Refused, match="outside this package"):
        dl.destination(ROOT / "data", "documentation")
