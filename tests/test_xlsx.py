"""The standard-library .xlsx reader."""
import pandas as pd
import pytest

from delayed_entry.adapters.xlsx import read_xlsx, sheet_names
from helpers import write_xlsx


@pytest.fixture
def book(tmp_path):
    return write_xlsx(tmp_path / "book.xlsx", {
        "Patients": [["id", "age", "dead", "dx_date", "note"],
                     ["a", 61, True, 42005, "first"],        # 42005 is 2015-01-01
                     ["b", 70.5, False, None, None],
                     ["c", None, None, 42400.5, "x & y"]],
        "Second": [["only"], ["value"]]},
        date_cells={("Patients", 1, 3), ("Patients", 3, 3)})


def test_the_first_sheet_is_read_with_its_header(book):
    df = read_xlsx(book)
    assert list(df.columns) == ["id", "age", "dead", "dx_date", "note"]
    assert df["id"].tolist() == ["a", "b", "c"]


def test_numbers_booleans_and_blanks(book):
    df = read_xlsx(book)
    assert df.loc[0, "age"] == 61 and df.loc[1, "age"] == 70.5 and pd.isna(df.loc[2, "age"])
    assert df.loc[0, "dead"] is True and df.loc[1, "dead"] is False


def test_dates_come_back_as_text(book):
    df = read_xlsx(book)
    assert df.loc[0, "dx_date"] == "2015-01-01" and df.loc[2, "dx_date"] == "2016-01-31 12:00:00"
    assert pd.isna(df.loc[1, "dx_date"])


def test_shared_strings_keep_their_text(book):
    assert read_xlsx(book).loc[2, "note"] == "x & y"


def test_a_sheet_by_name(book):
    assert sheet_names(book) == ["Patients", "Second"]
    assert read_xlsx(book, "Second")["only"].tolist() == ["value"]
    with pytest.raises(ValueError, match="no sheet 'Third'"):
        read_xlsx(book, "Third")
