"""Small counts and local folders. Ported from the CRC BPC project's tests/test_privacy.py and
test_paths_and_repo.py, with the minimum count passed in."""

import pandas as pd
import pytest

from delayed_entry.privacy import check_not_cloud, is_small, show, suppress_small_cells


def test_small_counts_are_masked_and_dependents_blanked():
    t = pd.DataFrame({"group": ["a", "b", "c"], "n": [5, 11, 0], "median": [10.0, 12.0, None]})
    out = suppress_small_cells(t, ["n"], ["median"], min_count=11)
    assert out["n"].tolist() == ["<11", 11, 0]
    assert out["median"].tolist()[:2] == ["", 12.0]


def test_the_boundary_is_the_minimum_count():
    assert [is_small(n, 11) for n in (0, 1, 10, 11)] == [False, True, True, False]
    assert [is_small(n, 20) for n in (11, 19, 20)] == [True, True, False]
    assert show(7, 11) == "<11" and show(0, 11) == "0" and show(25, 20) == "25"


def test_the_original_table_is_not_modified():
    t = pd.DataFrame({"n": [3]})
    suppress_small_cells(t, ["n"], min_count=11)
    assert t["n"].tolist() == [3]


def test_cloud_synced_folders_are_refused(tmp_path):
    for name in ("OneDrive", "Dropbox", "Google Drive", "iCloud"):
        with pytest.raises(ValueError, match="cloud-synced"):
            check_not_cloud(tmp_path / name / "data")
    assert check_not_cloud(tmp_path / "data") == tmp_path / "data"
