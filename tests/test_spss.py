"""SPSS (.sav) files through the flat-file adapter, on made-up files written at test time with pyreadstat (the
spss extra). Where pyreadstat isn't installed, or can't be loaded (Windows' Smart App Control can block its compiled
part), the tests that need it are skipped; the message a person then gets is tested everywhere."""
import datetime as dt
import sys

import pandas as pd
import pytest
import yaml

from delayed_entry.adapters import flatfile
from delayed_entry.adapters.flatfile import MappingError, read_source
from delayed_entry.tables import table_fields
from helpers import crc_settings

try:
    import pyreadstat
except ImportError:            # not installed, or blocked from loading
    pyreadstat = None
needs_pyreadstat = pytest.mark.skipif(pyreadstat is None, reason="pyreadstat isn't installed or can't be loaded here")

S = crc_settings()
DX = {1001: dt.date(2015, 1, 10), 1002: dt.date(2016, 3, 1), 1003: dt.date(2017, 6, 30)}
LAST = {1001: dt.date(2016, 1, 10), 1002: dt.date(2019, 3, 1), 1003: dt.date(2017, 12, 31)}


def test_without_pyreadstat_an_spss_file_says_what_it_needs(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "pyreadstat", None)          # as if it can't be imported
    (tmp_path / "patients.sav").write_bytes(b"")
    with pytest.raises(MappingError, match=r"needs the pyreadstat package .*delayed-entry-pipeline\[spss\]"):
        read_source(tmp_path / "patients.sav")


def write_patients(path):
    ids = list(DX)
    df = pd.DataFrame({"pid": [float(i) for i in ids], "sex": [1.0, 2.0, 2.0], "name": ["Müller", "Okafor", "Lee"],
                       "dx": [DX[i] for i in ids], "last": [LAST[i] for i in ids], "dead": [1.0, 0.0, 1.0]})
    pyreadstat.write_sav(df, str(path), variable_value_labels={"sex": {1.0: "Male", 2.0: "Female"}})
    return path


@needs_pyreadstat
def test_labels_dates_text_and_whole_numbers_are_read(tmp_path):
    df = read_source(write_patients(tmp_path / "patients.sav"))
    assert df["sex"].tolist() == ["Male", "Female", "Female"]               # SPSS's value labels
    assert str(df["pid"].dtype) == "Int64" and df["pid"].tolist() == [1001, 1002, 1003]
    assert df["name"].tolist() == ["Müller", "Okafor", "Lee"]
    assert pd.Timestamp(df.loc[0, "dx"]) == pd.Timestamp(DX[1001])


@needs_pyreadstat
def test_codes_can_be_read_in_place_of_labels(tmp_path):
    codes = read_source(write_patients(tmp_path / "patients.sav"), value_labels=False)
    assert codes["sex"].tolist() == [1, 2, 2]


@needs_pyreadstat
def test_an_spss_file_reads_into_the_standard_tables(tmp_path):
    write_patients(tmp_path / "patients.sav")
    pd.DataFrame({"pid": ["1001", "1002", "1003"], "sid": ["S-1", "S-2", "S-3"], "days": ["30", "40", "50"]}).to_csv(
        tmp_path / "tests.csv", index=False)
    pd.DataFrame({"sid": ["S-1"], "gene": ["GENE1"], "p": ["p.G12D"]}).to_csv(tmp_path / "results.csv", index=False)
    pd.DataFrame({"panel": ["P"], "gene": ["GENE1"]}).to_csv(tmp_path / "panels.csv", index=False)
    absent = {"absent": True}
    m = {"mapping_format": 1, "status": "confirmed", "confirmed_by": "A. Data Manager", "confirmed_on": "2026-11-18",
         "missing_values": [""], "diagnosis_date": {"column": "dx"},
         "files": {"patients": "patients.sav", "tests": "tests.csv", "results": "results.csv", "panels": "panels.csv"},
         "profile": {"entry_event": "report", "panel_coverage": "genes_only", "panel_sizes": "not_given",
                     "death_capture": "linked", "time_format": "real_dates", "privacy_min_count": 11,
                     "histology_vocabulary": "made up", "age_top_code": None, "known_quirks": [], "sources": {}},
         "tables": {
             "patients": {f.name: absent for f in table_fields("patients")} | {
                 "patient_id": {"column": "pid"}, "vital_status": {"column": "dead"}, "sex": {"column": "sex"},
                 "days_to_death_or_last_followup": {"column": "last", "transform": {"days_since_diagnosis": None}}},
             "tests": {f.name: absent for f in table_fields("tests")} | {
                 "patient_id": {"column": "pid"}, "test_id": {"column": "sid"}, "modality": {"constant": "ngs"},
                 "entry_event": {"constant": "report"}, "days_to_entry_event": {"column": "days"}},
             "results": {f.name: absent for f in table_fields("results")} | {
                 "test_id": {"column": "sid"}, "kind": {"constant": "mutation"}, "gene": {"column": "gene"},
                 "variant": {"column": "p"}},
             "panels": {f.name: absent for f in table_fields("panels")} | {
                 "panel_id": {"column": "panel"}, "gene": {"column": "gene"}}}}
    path = tmp_path / "mapping.yaml"
    path.write_text(yaml.safe_dump(m, sort_keys=False), encoding="utf-8")
    p = flatfile.read(path, S).tables.patients.set_index("patient_id")
    assert sorted(p.index) == ["1001", "1002", "1003"] and p.loc["1001", "sex"] == "Male"
    for pid in DX:
        assert p.loc[str(pid), "days_to_death_or_last_followup"] == (LAST[pid] - DX[pid]).days
    assert p.loc["1002", "vital_status"] == 0


@needs_pyreadstat
def test_a_file_that_isnt_spss_is_named(tmp_path):
    (tmp_path / "broken.sav").write_bytes(b"not an SPSS file")
    with pytest.raises(MappingError, match="can't be read as an SPSS file"):
        read_source(tmp_path / "broken.sav")
