"""Guards on what this repository may hold: code, settings, documentation and made-up test data, never real data,
a local machine's paths or a credential. They read the files git tracks, so they run in a git checkout (and on
every push, on GitHub) and are skipped elsewhere."""
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = "tests/fixtures/"
MAX_BYTES = 1_000_000
DATA_SUFFIXES = {".csv", ".tsv", ".txt", ".parquet", ".xlsx", ".xls", ".duckdb", ".db", ".sqlite", ".feather",
                 ".rds", ".rdata", ".sas7bdat", ".xpt", ".dta", ".sav", ".maf", ".vcf", ".bam", ".gz", ".zip", ".pdf", ".docx"}
LOCAL_ONLY = ("RELEASE_AUDIT.md", "docs/local/")
LOCAL_PATH = re.compile(r"\b[A-Za-z]:[\\/]+Users[\\/]+|(?<![\w.~])/(?:Users|home)/[A-Za-z0-9._-]+/")
CREDENTIAL = re.compile("|".join([
    r"gh[pousr]_[A-Za-z0-9]{36,}",                         # GitHub tokens
    r"github_pat_[A-Za-z0-9_]{40,}",
    r"AKIA[0-9A-Z]{16}",                                   # AWS access key IDs
    r"-----BEGIN [A-Z ]*PRIVATE KEY-----",
    r"eyJ[A-Za-z0-9_-]{10,}\.eyJ[A-Za-z0-9_-]{10,}",       # a JSON web token, such as a Synapse access token
]))


def tracked() -> list[str]:
    if shutil.which("git") is None or not (ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    out = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, check=True).stdout
    return [p for p in out.decode("utf-8").split("\0") if p]


def texts():
    for p in tracked():
        try:
            yield p, (ROOT / p).read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue


def test_no_data_file_outside_the_made_up_fixtures():
    assert [p for p in tracked() if Path(p).suffix.lower() in DATA_SUFFIXES and not p.startswith(FIXTURES)] == []


def test_the_fixtures_say_they_are_made_up():
    folders = {p.split("/")[2] for p in tracked() if p.startswith(FIXTURES) and p.count("/") >= 3}
    for f in folders:
        readme = ROOT / FIXTURES / f / "README.md"
        assert readme.exists() and "Made-up" in readme.read_text(encoding="utf-8"), f


def test_no_file_is_over_1_mb():
    assert [p for p in tracked() if (ROOT / p).stat().st_size > MAX_BYTES] == []


def test_no_file_holds_a_local_user_folder():
    assert [p for p, text in texts() if LOCAL_PATH.search(text)] == []


def test_no_file_holds_a_credential():
    assert [p for p, text in texts() if CREDENTIAL.search(text)] == []


@pytest.mark.parametrize("text, pattern", [
    ("C:" + "\\Users\\someone\\data.csv", LOCAL_PATH),
    ("C:" + "/Us" + "ers/someone/data.csv", LOCAL_PATH),
    ("/ho" + "me/someone/data.csv", LOCAL_PATH),
    ("gh" + "p_" + "a1" * 18, CREDENTIAL),
    ("AK" + "IA" + "ABCDEFGH12345678", CREDENTIAL),
    ("-----BEGIN " + "RSA PRIVATE KEY-----", CREDENTIAL),
    ("ey" + "J" + "a" * 12 + ".ey" + "J" + "b" * 12, CREDENTIAL),
])
def test_the_patterns_catch_what_they_guard_against(text, pattern):
    """Made-up examples, built here so that no file holds one."""
    assert pattern.search(text)


@pytest.mark.parametrize("text", ["~/.synapseConfig", "tmp_path / 'OneDrive'", "C:" + "\\R\\bin\\x64\\Rterm.exe"])
def test_ordinary_text_is_not_flagged(text):
    assert not LOCAL_PATH.search(text) and not CREDENTIAL.search(text)


def test_the_local_only_files_are_not_tracked():
    assert [p for p in tracked() if p.startswith(LOCAL_ONLY)] == []
