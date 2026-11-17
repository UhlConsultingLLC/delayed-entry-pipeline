"""SAS and Stata files through the flat-file adapter, on made-up data: a Stata patients file with value labels and
dates, SAS transport (.xpt) files whose dates are numbers of days since 1960-01-01, and numeric patient IDs that must
still join across files. The files are written to a temporary folder; none is kept in the repository."""
import datetime as dt

import pandas as pd
import pytest
import yaml

from delayed_entry.adapters import flatfile
from delayed_entry.adapters.flatfile import MappingError, read_source
from delayed_entry.tables import table_fields
from helpers import crc_settings
from sas_xport import write_xport

S = crc_settings()
SAS_EPOCH = dt.date(1960, 1, 1)
NOTE = "a made-up note, here only to make each record longer than 80 bytes"
PROFILE = {"entry_event": "report", "panel_coverage": "genes_only", "panel_sizes": "not_given",
           "death_capture": "linked", "time_format": "real_dates", "privacy_min_count": 11,
           "histology_vocabulary": "made up", "age_top_code": None, "known_quirks": [], "sources": {}}
DX = {1001: dt.date(2015, 1, 10), 1002: dt.date(2016, 3, 1), 1003: dt.date(2017, 6, 30), 1004: dt.date(2018, 2, 14)}
LAST = {1001: dt.date(2016, 1, 10), 1002: dt.date(2019, 3, 1), 1003: dt.date(2017, 12, 31), 1004: dt.date(2020, 2, 14)}
REPORTS = [(1001, "S-1", dt.date(2015, 2, 9), "PANEL-1", 1), (1001, "S-2", dt.date(2015, 8, 1), "PANEL-2", 2),
           (1002, "S-3", dt.date(2016, 2, 20), "PANEL-1", 1), (1003, "S-4", None, "PANEL-2", 1),
           (1004, "S-5", dt.date(2018, 4, 1), "PANEL-1", 1)]


def sas_days(d: dt.date | None):
    return None if d is None else (d - SAS_EPOCH).days


@pytest.fixture
def files(tmp_path):
    ids = list(DX)
    pd.DataFrame({
        "pid": ids, "age": [61, 70, 55, 48],
        "sex": pd.Categorical(["Male", "Female", "Female", "Male"]),
        "stage": pd.Categorical(["Stage IV", "Stage IV", "Stage III", "Stage IV"]),
        "dx": pd.to_datetime([DX[i] for i in ids]), "last": pd.to_datetime([LAST[i] for i in ids]),
        "dead": [1, 0, 1, 0],
    }).to_stata(tmp_path / "patients.dta", write_index=False, convert_dates={"dx": "td", "last": "td"})
    write_xport(tmp_path / "tests.xpt", {
        "PT_ID": [r[0] for r in REPORTS], "SAMPLE": [r[1] for r in REPORTS],
        "REPORT": [sas_days(r[2]) for r in REPORTS], "PANEL": [r[3] for r in REPORTS],
        "ORDER": [r[4] for r in REPORTS], "NOTE": [NOTE] * len(REPORTS)}, name="TESTS")
    write_xport(tmp_path / "results.xpt", {
        "SAMPLE": ["S-1", "S-3", "S-5"], "GENE": ["KRAS", "BRAF", "TP53"],
        "HGVSP": ["p.G12D", "p.V600E", "p.R175H"], "CLASS": ["Missense_Mutation"] * 3, "NOTE": [NOTE] * 3},
        name="RESULTS")
    pd.DataFrame({"panel_id": ["PANEL-1", "PANEL-2"], "gene": ["KRAS", "BRAF"]}).to_csv(tmp_path / "panels.csv", index=False)
    return tmp_path


def mapping(folder, **changes):
    absent = {"absent": True}
    m = {"mapping_format": 1, "status": "confirmed", "confirmed_by": "A. Data Manager", "confirmed_on": "2026-11-17",
         "missing_values": [""], "diagnosis_date": {"column": "dx"},
         "files": {"patients": "patients.dta", "tests": "tests.xpt", "results": "results.xpt", "panels": "panels.csv"},
         "profile": dict(PROFILE),
         "tables": {
             "patients": {f.name: absent for f in table_fields("patients")} | {
                 "patient_id": {"column": "pid"}, "vital_status": {"column": "dead"},
                 "days_to_death_or_last_followup": {"column": "last", "transform": {"days_since_diagnosis": None}},
                 "age_at_diagnosis": {"column": "age"}, "sex": {"column": "sex"}, "stage_at_diagnosis": {"column": "stage"}},
             "tests": {f.name: absent for f in table_fields("tests")} | {
                 "patient_id": {"column": "PT_ID"}, "test_id": {"column": "SAMPLE"}, "modality": {"constant": "ngs"},
                 "entry_event": {"constant": "report"}, "panel_id": {"column": "PANEL"}, "test_order": {"column": "ORDER"},
                 "days_to_entry_event": {"column": "REPORT",
                                         "transform": {"days_since_diagnosis": {"days_from": "1960-01-01"}}}},
             "results": {f.name: absent for f in table_fields("results")} | {
                 "test_id": {"column": "SAMPLE"}, "kind": {"constant": "mutation"}, "gene": {"column": "GENE"},
                 "variant": {"column": "HGVSP"}, "variant_class": {"column": "CLASS"}},
             "panels": {f.name: absent for f in table_fields("panels")} | {
                 "panel_id": {"column": "panel_id"}, "gene": {"column": "gene"}}}}
    m.update(changes)
    path = folder / "mapping.yaml"
    path.write_text(yaml.safe_dump(m, sort_keys=False), encoding="utf-8")
    return path


def test_stata_and_sas_files_read_into_the_standard_tables(files):
    res = flatfile.read(mapping(files), S)
    p = res.tables.patients.set_index("patient_id")
    assert sorted(p.index) == ["1001", "1002", "1003", "1004"]          # whole numbers, not "1001.0"
    assert p.loc["1001", "sex"] == "Male" and p.loc["1003", "stage_at_diagnosis"] == "Stage III"   # Stata's labels
    assert p.loc["1002", "vital_status"] == 0 and p.loc["1001", "age_at_diagnosis"] == 61
    for pid in DX:
        assert p.loc[str(pid), "days_to_death_or_last_followup"] == (LAST[pid] - DX[pid]).days
    t = res.tables.tests.set_index("test_id")
    for pid, sample, when, panel, order in REPORTS:          # SAS dates: days since 1960-01-01, less the diagnosis
        got = t.loc[sample, "days_to_entry_event"]
        assert pd.isna(got) if when is None else got == (when - DX[pid]).days
        assert t.loc[sample, "patient_id"] == str(pid) and t.loc[sample, "panel_id"] == panel
        assert t.loc[sample, "test_order"] == order
    r = res.tables.results.set_index("test_id")
    assert r.loc["S-3", "variant"] == "p.V600E" and r.loc["S-1", "gene"] == "KRAS"
    assert {p.name for p in res.input_files} >= {"patients.dta", "tests.xpt", "results.xpt"}


def test_stata_codes_can_be_read_in_place_of_labels(files):
    labels = read_source(files / "patients.dta")
    codes = read_source(files / "patients.dta", value_labels=False)
    assert set(labels["sex"]) == {"Male", "Female"} and set(codes["sex"]) <= {0, 1}
    m = mapping(files, files={"patients": {"path": "patients.dta", "value_labels": False}, "tests": "tests.xpt",
                              "results": "results.xpt", "panels": "panels.csv"})
    sex = flatfile.read(m, S).tables.patients.set_index("patient_id")["sex"]
    assert set(sex) <= {"0", "1"}


def test_stata_labels_that_are_not_unique_name_the_way_out(tmp_path):
    path = tmp_path / "twice.dta"
    try:
        pd.DataFrame({"x": [1, 2]}).to_stata(path, write_index=False, value_labels={"x": {1: "same", 2: "same"}})
    except (ValueError, TypeError):
        pytest.skip("this pandas won't write a repeated value label")
    with pytest.raises(MappingError, match="value_labels: false"):
        read_source(path)
    assert list(read_source(path, value_labels=False)["x"]) == [1, 2]


def test_a_sas_transport_files_text_is_decoded_as_the_file_says(tmp_path):
    latin = write_xport(tmp_path / "latin.xpt", {"NAME": ["Müller"], "NOTE": [NOTE + " " + NOTE]})
    assert read_source(latin)["NAME"].tolist() == ["Müller"]                    # ISO-8859-1 by default
    utf8 = write_xport(tmp_path / "utf8.xpt", {"NAME": ["Müller"], "NOTE": [NOTE + " " + NOTE]}, encoding="utf-8")
    assert read_source(utf8, encoding="utf-8")["NAME"].tolist() == ["Müller"]
    assert read_source(utf8)["NAME"].tolist() != ["Müller"]                     # read as Latin-1, it is garbled


def test_a_date_held_as_seconds_since_1960_is_read_with_seconds_from():
    seconds = pd.Series([(dt.datetime(2015, 2, 9, 13, 30) - dt.datetime(1960, 1, 1)).total_seconds(), None])
    got = flatfile.to_dates(seconds, {"seconds_from": "1960-01-01"})
    assert got.iloc[0] == pd.Timestamp("2015-02-09 13:30") and pd.isna(got.iloc[1])


def test_the_draft_keeps_a_files_options(files):
    draft = flatfile.write_draft({"patients": {"path": "patients.dta", "value_labels": False},
                                  "tests": {"path": "tests.xpt", "encoding": "latin-1"}}, files / "draft.yaml")
    m = yaml.safe_load(draft.read_text(encoding="utf-8"))
    assert m["files"]["patients"] == {"path": "patients.dta", "value_labels": False}
    assert m["files"]["tests"] == {"path": "tests.xpt", "encoding": "latin-1"}
    assert m["tables"]["tests"]["test_id"]["column"] in (None, "SAMPLE")


@pytest.mark.parametrize("change, problem", [
    ({"files": {"patients": {"path": "patients.dta", "value_labels": "no"}}}, "value_labels must be true or false"),
    ({"files": {"patients": {"path": "patients.dta", "labels": True}}}, "give a path, or"),
    ({"diagnosis_date": {"column": "dx", "format": "%Y-%m-%d", "days_from": "1960-01-01"}}, "give one of format"),
    ({"diagnosis_date": {"column": "dx", "days_from": "soon"}}, "days_from must be a date"),
])
def test_a_badly_given_option_is_named(files, change, problem):
    path = mapping(files)
    m = yaml.safe_load(path.read_text(encoding="utf-8"))
    for key, value in change.items():
        if key == "files":
            m["files"].update(value)
        else:
            m[key] = value
    path.write_text(yaml.safe_dump(m, sort_keys=False), encoding="utf-8")
    with pytest.raises(MappingError, match=problem):
        flatfile.load_mapping(path)


def test_a_days_since_diagnosis_with_two_forms_is_refused(files):
    path = mapping(files)
    m = yaml.safe_load(path.read_text(encoding="utf-8"))
    m["tables"]["tests"]["days_to_entry_event"]["transform"]["days_since_diagnosis"]["format"] = "%Y-%m-%d"
    path.write_text(yaml.safe_dump(m, sort_keys=False), encoding="utf-8")
    with pytest.raises(MappingError, match="give one of format, days_from or seconds_from"):
        flatfile.load_mapping(path)
