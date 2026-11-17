"""The command line, and checks on the project folder itself."""

import pytest

from delayed_entry.__main__ import main
from delayed_entry.settings import blank_template
from helpers import CRC_SETTINGS, FIXTURE, ROOT

IGNORE = {".git", ".venv", "venv", "__pycache__", ".pytest_cache", "build", "dist"}
DATA_SUFFIXES = {".csv", ".tsv", ".txt", ".parquet", ".rds", ".maf", ".seg", ".xlsx", ".sas7bdat", ".feather",
                 ".h5ad", ".duckdb", ".vcf"}
TEXT_SUFFIXES = {".py", ".md", ".yml", ".yaml", ".sql", ".toml", ".txt", ".cfg"}


def project_files():
    for p in ROOT.rglob("*"):
        parts = p.relative_to(ROOT).parts
        # an install's build folder (*.egg-info, made by pip install -e .) is not the project's own
        if p.is_file() and not (set(parts) & IGNORE) and not any(x.endswith(".egg-info") for x in parts):
            yield p


def test_the_command_line_writes_the_template(tmp_path, capsys):
    main(["template", str(tmp_path / "t.yaml")])
    assert (tmp_path / "t.yaml").read_text(encoding="utf-8") == blank_template()


def test_the_command_line_checks_settings(tmp_path, capsys):
    main(["check-settings", str(CRC_SETTINGS)])
    assert "settings OK: 58 choices, 0 taken from the defaults file" in capsys.readouterr().out
    main(["template", str(tmp_path / "blank.yaml")])
    capsys.readouterr()
    with pytest.raises(SystemExit) as stop:
        main(["check-settings", str(tmp_path / "blank.yaml")])
    assert stop.value.code == 1 and "59 problem(s)" in capsys.readouterr().out   # 58 choices and the optional one


def test_the_command_line_locks_the_settings(tmp_path, capsys):
    main(["lock-settings", str(CRC_SETTINGS), str(tmp_path / "lock"), "--locked-by", "A. Reviewer"])
    assert (tmp_path / "lock" / "settings_lock.json").exists()
    main(["template", str(tmp_path / "blank.yaml")])
    with pytest.raises(Exception):
        main(["lock-settings", str(tmp_path / "blank.yaml"), str(tmp_path / "lock2"), "--locked-by", "A. Reviewer"])
    assert not (tmp_path / "lock2").exists()


def test_the_command_line_screens_locks_and_verifies(tmp_path, capsys):
    work = tmp_path / "work"
    main(["screen", str(CRC_SETTINGS), str(work), "--clinical", str(FIXTURE)])
    out = capsys.readouterr().out
    assert "SA1: BRAF V600E vs Wild type: descriptive only" in out
    assert (work / "readiness_report.md").exists() and (work / "readiness_report.yaml").exists()
    main(["lock", str(CRC_SETTINGS), str(work), "--clinical", str(FIXTURE), "--locked-by", "A. Reviewer"])
    assert "locked: version 1" in capsys.readouterr().out
    main(["verify-lock", str(work / "plan_v1" / "manifest.json")])
    assert "lock matches" in capsys.readouterr().out


def test_the_command_line_refuses_a_cloud_synced_work_folder(tmp_path):
    with pytest.raises(ValueError, match="cloud-synced"):
        main(["screen", str(CRC_SETTINGS), str(tmp_path / "Dropbox" / "work"), "--clinical", str(FIXTURE)])


def test_no_data_files_outside_the_test_fixtures():
    fixtures = ROOT / "tests" / "fixtures"
    stray = [p for p in project_files() if p.suffix.lower() in DATA_SUFFIXES and fixtures not in p.parents]
    assert stray == []


def test_the_fixtures_are_made_up():
    readme = (FIXTURE / "README.md").read_text(encoding="utf-8")
    assert "Made-up" in readme and "No real patients" in readme


def test_no_large_files():
    assert [p for p in project_files() if p.stat().st_size > 2_000_000] == []


def test_the_package_data_is_declared():
    toml = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert '"*.yml"' in toml and '"**/*.yml"' in toml and '"**/*.sql"' in toml
    pkg = ROOT / "src" / "delayed_entry"
    for f in ("defaults.yml", "adapters/flatfile_synonyms.yml", "adapters/genie_bpc/profile.yml",
              "adapters/genie_bpc/labels.yml", "adapters/genie_bpc/standard.sql", "adapters/genie_bpc/link_checks.sql",
              "adapters/genie_bpc/cohorts/crc_v2.0/raw_files.yml"):
        assert (pkg / f).exists(), f
