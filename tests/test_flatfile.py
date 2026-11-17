"""The generic flat-file adapter: CSV/TSV, Excel and Parquet, with a mapping a person confirms."""
import duckdb
import pandas as pd
import pytest
import yaml

from delayed_entry.adapters import flatfile
from delayed_entry.adapters.flatfile import MappingError, UnconfirmedMappingError
from delayed_entry.screen import screen
from delayed_entry.tables import TABLES, outcome_blind, table_fields
from helpers import crc_settings, make_tables, write_xlsx

S = crc_settings()
PROFILE = {"entry_event": "report", "panel_coverage": "genes_only", "panel_sizes": "not_given",
           "death_capture": "linked", "time_format": "intervals", "privacy_min_count": 11,
           "histology_vocabulary": "registry categories", "age_top_code": None, "known_quirks": [],
           "sources": {"entry_event": "made-up data"}}


def write_table(df: pd.DataFrame, path):
    """Write a made-up table in the format its suffix names."""
    df = df.astype(object).where(df.notna(), None)
    if path.suffix == ".csv":
        df.to_csv(path, index=False)
    elif path.suffix == ".tsv":
        df.to_csv(path, index=False, sep="\t")
    elif path.suffix == ".xlsx":
        write_xlsx(path, {"Sheet1": [list(df.columns)] + df.values.tolist()})
    else:
        con = duckdb.connect()
        con.register("t", df.astype(str).where(df.notna(), None))
        con.execute(f"COPY t TO '{path.as_posix()}' (FORMAT parquet)")
        con.close()
    return path


def confirm(draft_path, **edits):
    """What a person does: check the draft, fill in the profile and confirm it."""
    m = yaml.safe_load(draft_path.read_text(encoding="utf-8"))
    m.update(status="confirmed", confirmed_by="A. Data Manager", confirmed_on="2026-11-13", profile=dict(PROFILE))
    for t, fields in (m.get("tables") or {}).items():
        for name, entry in fields.items():
            if entry.get("column") is None:
                fields[name] = {"absent": True}
    for key, value in edits.items():
        m[key] = value
    draft_path.write_text(yaml.safe_dump(m, sort_keys=False), encoding="utf-8")
    return draft_path


@pytest.fixture
def source(tmp_path):
    """Made-up standard tables written as four kinds of file, with source column names that differ."""
    t = make_tables({"RAS-mutant": 40, "BRAF V600E": 35, "Wild type": 50}, msi={"Wild type": (35, 10, 5)},
                    treated={"RAS-mutant": 40, "BRAF V600E": 35, "Wild type": 50})
    patients = t.patients.rename(columns={"patient_id": "Patient ID", "age_at_diagnosis": "age_dx",
                                          "stage_at_diagnosis": "Stage", "vital_status": "dead",
                                          "days_to_death_or_last_followup": "os_days"})
    files = {"patients": write_table(patients, tmp_path / "patients.csv"),
             "tests": write_table(t.tests.rename(columns={"test_id": "sample_id"}), tmp_path / "tests.tsv"),
             "results": write_table(t.results.rename(columns={"gene": "Hugo_Symbol", "variant": "HGVSp_Short"}),
                                    tmp_path / "results.parquet"),
             "panels": write_table(t.panels, tmp_path / "panels.xlsx"),
             "treatments": write_table(t.treatments, tmp_path / "treatments.csv")}
    return t, {k: v.name for k, v in files.items()}, tmp_path


def test_the_draft_matches_columns_by_name(source):
    _, files, d = source
    draft = flatfile.write_draft(files, d / "mapping.yaml")
    m = yaml.safe_load(draft.read_text(encoding="utf-8"))
    assert m["status"] == "draft"
    assert m["tables"]["patients"]["patient_id"] == {"column": "Patient ID"}
    assert m["tables"]["patients"]["age_at_diagnosis"] == {"column": "age_dx"}
    assert m["tables"]["results"]["variant"] == {"column": "HGVSp_Short"}
    assert m["tables"]["tests"]["test_id"] == {"column": "sample_id"}
    assert all(v is None for k, v in m["profile"].items() if k not in ("known_quirks", "sources"))


def test_a_draft_is_never_used(source):
    _, files, d = source
    draft = flatfile.write_draft(files, d / "mapping.yaml")
    with pytest.raises(UnconfirmedMappingError, match="never used"):
        flatfile.read(draft, S)


@pytest.mark.parametrize("edit", [{"confirmed_by": None}, {"confirmed_on": None}, {"status": "draft"}])
def test_a_mapping_needs_a_person_a_date_and_the_confirmed_status(source, edit):
    _, files, d = source
    path = confirm(flatfile.write_draft(files, d / "mapping.yaml"), **edit)
    with pytest.raises(UnconfirmedMappingError):
        flatfile.load_mapping(path)


def test_a_field_left_undecided_is_refused(source):
    _, files, d = source
    path = confirm(flatfile.write_draft(files, d / "mapping.yaml"))
    m = yaml.safe_load(path.read_text(encoding="utf-8"))
    m["tables"]["patients"]["sex"] = {"column": None}
    path.write_text(yaml.safe_dump(m), encoding="utf-8")
    with pytest.raises(MappingError, match=r"tables.patients.sex: not decided"):
        flatfile.load_mapping(path)


def test_the_profile_must_be_filled_in(source):
    _, files, d = source
    path = confirm(flatfile.write_draft(files, d / "mapping.yaml"))
    m = yaml.safe_load(path.read_text(encoding="utf-8"))
    m["profile"]["death_capture"] = None
    path.write_text(yaml.safe_dump(m), encoding="utf-8")
    with pytest.raises(MappingError, match="death_capture must be one of"):
        flatfile.load_mapping(path)


def test_a_confirmed_mapping_reads_every_format_into_the_standard_tables(source):
    t, files, d = source
    res = flatfile.read(confirm(flatfile.write_draft(files, d / "mapping.yaml")), S)
    assert res.profile.adapter == "flat_file" and len(res.input_files) == 6
    for name in ("patients", "tests", "results", "panels", "treatments"):
        got, want = res.tables.table(name), t.table(name)
        assert len(got) == len(want), name
    p = res.tables.patients.set_index("patient_id")
    assert p.loc["RAS-0000", "age_at_diagnosis"] == 60 and p.loc["RAS-0000", "stage_at_diagnosis"] == "Stage IV"
    assert res.tables.results["variant"].dropna().isin(["p.G12D", "p.V600E", "p.R175H"]).all()


def test_the_flat_file_route_screens_the_same_as_the_tables_themselves(source):
    t, files, d = source
    res = flatfile.read(confirm(flatfile.write_draft(files, d / "mapping.yaml")), S)
    direct = screen(outcome_blind(t), S, res.profile)
    via_files = screen(outcome_blind(res.tables), S, res.profile)
    assert via_files.decisions == direct.decisions and via_files.group_sizes == direct.group_sizes


def test_constants_and_transforms(tmp_path):
    pd.DataFrame({"id": ["a", "b", "c"], "dx": ["2015-01-01", "2016-03-01", "2017-06-30"],
                  "status": ["Dead", "Alive", "Dead"], "os_months": ["10", "20.5", "x"]}).to_csv(tmp_path / "pt.csv", index=False)
    pd.DataFrame({"id": ["a", "b", "c"], "sid": ["a1", "b1", "c1"],
                  "report_date": ["2015-02-01", "2016-02-20", ""]}).to_csv(tmp_path / "tests.csv", index=False)
    pd.DataFrame({"sid": ["a1"], "gene": ["GENE1"], "p": ["p.G12D"]}).to_csv(tmp_path / "res.csv", index=False)
    pd.DataFrame({"panel": ["P"], "gene": ["GENE1"]}).to_csv(tmp_path / "panels.csv", index=False)
    absent = {"absent": True}
    m = {"mapping_format": 1, "status": "confirmed", "confirmed_by": "A. Data Manager", "confirmed_on": "2026-11-13",
         "missing_values": ["", "NA"], "diagnosis_date": {"column": "dx", "format": "%Y-%m-%d"},
         "files": {"patients": "pt.csv", "tests": "tests.csv", "results": "res.csv", "panels": "panels.csv"},
         "profile": dict(PROFILE, time_format="real_dates"),
         "tables": {
             "patients": {f.name: absent for f in table_fields("patients")} | {
                 "patient_id": {"column": "id"},
                 "vital_status": {"column": "status", "transform": {"map": {"Dead": 1, "Alive": 0}}},
                 "days_to_death_or_last_followup": {"column": "os_months", "transform": {"multiply": 30.4375}}},
             "tests": {f.name: absent for f in table_fields("tests")} | {
                 "patient_id": {"column": "id"}, "test_id": {"column": "sid"}, "modality": {"constant": "ngs"},
                 "entry_event": {"constant": "report"},
                 "days_to_entry_event": {"column": "report_date", "transform": {"days_since_diagnosis": {"format": "%Y-%m-%d"}}}},
             "results": {f.name: absent for f in table_fields("results")} | {
                 "test_id": {"column": "sid"}, "kind": {"constant": "mutation"}, "gene": {"column": "gene"},
                 "variant": {"column": "p"}},
             "panels": {f.name: absent for f in table_fields("panels")} | {
                 "panel_id": {"column": "panel"}, "gene": {"column": "gene"}}}}
    path = tmp_path / "mapping.yaml"
    path.write_text(yaml.safe_dump(m, sort_keys=False), encoding="utf-8")
    res = flatfile.read(path, S)
    p = res.tables.patients.set_index("patient_id")
    assert p["vital_status"].tolist() == [1, 0, 1]
    assert p.loc["b", "days_to_death_or_last_followup"] == pytest.approx(20.5 * 30.4375)
    assert pd.isna(p.loc["c", "days_to_death_or_last_followup"])
    tests = res.tables.tests.set_index("test_id")
    assert tests.loc["a1", "days_to_entry_event"] == 31 and tests.loc["b1", "days_to_entry_event"] == -10
    assert pd.isna(tests.loc["c1", "days_to_entry_event"]) and (tests["modality"] == "ngs").all()
    assert any("os_months" in n or "days_to_death_or_last_followup" in n for n in res.read_notes)


def test_a_mapped_column_that_is_not_in_the_file_is_named(source):
    _, files, d = source
    path = confirm(flatfile.write_draft(files, d / "mapping.yaml"))
    m = yaml.safe_load(path.read_text(encoding="utf-8"))
    m["tables"]["patients"]["sex"] = {"column": "Gender"}
    path.write_text(yaml.safe_dump(m), encoding="utf-8")
    with pytest.raises(MappingError, match="no column 'Gender'"):
        flatfile.read(path, S)


def test_an_unreadable_kind_of_file_is_refused(tmp_path):
    with pytest.raises(MappingError, match="can't read this kind of file"):
        flatfile.read_source(tmp_path / "data.sas7bdat")


def test_every_standard_table_appears_in_a_draft(source):
    _, files, d = source
    m = yaml.safe_load(flatfile.write_draft(files, d / "mapping.yaml").read_text(encoding="utf-8"))
    assert set(m["files"]) == set(TABLES)


# ---- Mapping checks an independent review asked for ----
def edit_mapping(path, change):
    m = yaml.safe_load(path.read_text(encoding="utf-8"))
    change(m)
    path.write_text(yaml.safe_dump(m, sort_keys=False), encoding="utf-8")
    return path


def test_a_misspelled_table_name_is_named(source):
    _, files, d = source
    path = confirm(flatfile.write_draft(files, d / "mapping.yaml"))
    edit_mapping(path, lambda m: m["files"].update(treatment=m["files"].pop("treatments")))
    with pytest.raises(MappingError, match="'treatment' is not a standard table"):
        flatfile.load_mapping(path)


def test_an_unquoted_yes_or_no_map_key_is_refused(source):
    _, files, d = source
    path = confirm(flatfile.write_draft(files, d / "mapping.yaml"))
    path.write_text(path.read_text(encoding="utf-8").replace(
        "    sex:\n      column: sex\n", "    sex:\n      column: sex\n      transform:\n        map: {Yes: Female, No: Male}\n"),
        encoding="utf-8")
    with pytest.raises(MappingError, match="put it in quotes"):
        flatfile.load_mapping(path)


@pytest.mark.parametrize("transform, message", [({"map": ["a", "b"]}, "map must be a mapping"),
                                                ({"multiply": "thirty"}, "multiply needs a number"),
                                                ({"days_since_diagnosis": {"format": "%Y"}}, "needs diagnosis_date")])
def test_a_badly_shaped_transform_is_refused(source, transform, message):
    _, files, d = source
    path = confirm(flatfile.write_draft(files, d / "mapping.yaml"))
    edit_mapping(path, lambda m: m["tables"]["patients"]["age_at_diagnosis"].update(transform=transform))
    with pytest.raises(MappingError, match=message):
        flatfile.load_mapping(path)


def test_repeated_patient_ids_with_diagnosis_dates_are_named_without_a_small_count(tmp_path):
    pd.DataFrame({"id": ["a", "a", "b"], "dx": ["2015-01-01"] * 3}).to_csv(tmp_path / "pt.csv", index=False)
    for name in ("tests", "res", "panels"):
        pd.DataFrame({"x": ["1"]}).to_csv(tmp_path / f"{name}.csv", index=False)
    absent = {"absent": True}
    m = {"mapping_format": 1, "status": "confirmed", "confirmed_by": "A. Data Manager", "confirmed_on": "2026-11-13",
         "diagnosis_date": {"column": "dx", "format": "%Y-%m-%d"},
         "files": {"patients": "pt.csv", "tests": "tests.csv", "results": "res.csv", "panels": "panels.csv"},
         "profile": dict(PROFILE),
         "tables": {t: {f.name: ({"column": "x"} if f.key else absent) for f in table_fields(t)}
                    for t in ("tests", "results", "panels")}
         | {"patients": {f.name: absent for f in table_fields("patients")} | {"patient_id": {"column": "id"}}}}
    path = tmp_path / "mapping.yaml"
    path.write_text(yaml.safe_dump(m, sort_keys=False), encoding="utf-8")
    with pytest.raises(MappingError, match="<11 patient IDs appear more than once"):
        flatfile.read(path, S)
