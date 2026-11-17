"""The OMOP CDM adapter, on made-up OMOP data only: the made-up GENIE BPC fixture written as an OMOP database
(omop_fixture.py), read back, screened and run.

The tests that run R need it configured (see test_run.py); without it they are skipped."""
import math
import os
import shutil
from pathlib import Path

import duckdb
import pandas as pd
import pytest
import yaml

from delayed_entry.__main__ import main
from delayed_entry.adapters import genie_bpc, omop
from delayed_entry.adapters.omop import MAPPED, MappingError, UnconfirmedMappingError, short_protein_change
from delayed_entry.lock import lock_plan
from delayed_entry.profile import ITEMS as PROFILE_ITEMS
from delayed_entry.run import ROptions, run
from delayed_entry.screen import screen
from delayed_entry.settings import load_settings
from delayed_entry.tables import outcome_blind
from helpers import FIXTURE, crc_settings
import omop_fixture as fx

R = os.environ.get("DELAYED_ENTRY_R")
needs_r = pytest.mark.skipif(not R, reason="R is not configured (set DELAYED_ENTRY_R)")


@pytest.fixture(scope="module")
def genie():
    return genie_bpc.read(FIXTURE, None, crc_settings())


@pytest.fixture(scope="module")
def source(tmp_path_factory, genie):
    """The fixture's OMOP database, its mapping and settings, and the GENIE test IDs' new names."""
    base = tmp_path_factory.mktemp("omop")
    ids = fx.build(genie.tables, base)
    return base, ids, fx.write_mapping(base), load_settings(fx.omop_settings(base))


@pytest.fixture(scope="module")
def read_back(source):
    base, ids, mapping, settings = source
    return omop.read(mapping, settings)


def copy_of(source, tmp_path) -> Path:
    shutil.copytree(source[0], tmp_path / "db")
    return tmp_path / "db"


def text(s: pd.Series) -> pd.Series:
    return s.astype("string").fillna("<NA>")


# ---- Protein changes ----
@pytest.mark.parametrize("given, short", [
    ("p.Gly12Asp", "p.G12D"), ("NP_004976.2:p.(Arg1450Ter)", "p.R1450*"), ("p.Glu1309AspfsTer4", "p.E1309Dfs*4"),
    ("p.Val600Glu", "p.V600E"), ("p.G12D", "p.G12D"), ("c.35G>A", "c.35G>A"),
])
def test_protein_changes_are_read_in_short_form(given, short):
    assert short_protein_change(given) == short
    assert pd.isna(short_protein_change(None))


# ---- The mapping ----
def test_the_template_is_a_draft_that_lists_every_choice_and_is_never_used(tmp_path):
    path = omop.write_template(tmp_path / "m.yaml")
    body = path.read_text(encoding="utf-8")
    assert all(f"  {name}: null" in body for name in MAPPED) and all(f"  {item}: null" in body for item in PROFILE_ITEMS)
    with pytest.raises(UnconfirmedMappingError, match="never used"):
        omop.load_mapping(path)


def test_the_command_line_writes_the_template(tmp_path, capsys):
    main(["omop-template", str(tmp_path / "m.yaml")])
    assert (tmp_path / "m.yaml").exists() and "m.yaml" in capsys.readouterr().out


def broken(**changes):
    m = fx.mapping()
    for path, value in changes.items():
        node = m
        *parts, last = path.split("__")
        for p in parts:
            node = node[p]
        if value is KeyError:
            del node[last]
        else:
            node[last] = value
    return m


@pytest.mark.parametrize("changes, problem", [
    ({"patients__sex": KeyError}, "patients.sex: not decided"),
    ({"patients__sex": {"computed": True}}, "only age_at_diagnosis, year_of_diagnosis can be computed"),
    ({"patients__sex": {"person": "gender_source_value", "index_condition": "x"}}, "give exactly one of"),
    ({"patients__stage_at_diagnosis": {"measurement": "no such set", "value": "measurement_concept_id"}},
     "'no such set' is not a concept set"),
    ({"patients__histology": {"index_condition": "condition_source_value", "extract": "no group"}},
     "exactly one group"),
    ({"patients__sex": {"person": "gender_source_value", "labels": {True: "Male"}}}, "read as true or false"),
    ({"patients__vital_status": {"person": "x"}}, "the adapter works it out"),
    ({"extras__center": {"person": "x"}}, "extras.center: a standard field's name"),
    ({"diagnosis": "no such set"}, "diagnosis: must name a concept set"),
    ({"concept_sets__stage": {"concepts": ["Stage I"], "include_descendants": False}}, "concept IDs, whole numbers"),
    ({"sequencing__entry_event": "sometime"}, "sequencing.entry_event"),
    ({"biomarker_readings": [{"measurement": "msi or mmr", "kind": "mutation"}]}, "biomarker_readings[0]"),
    ({"treatments": {"drugs": "cancer drugs", "name": "brand"}}, "treatments: null, or"),
    ({"profile__time_format": "intervals"}, "holds dates"),
    ({"profile__time_format": "shifted_dates"}, "shifted dates don't give the calendar year"),
    ({"profile__entry_event": "order"}, "must say what sequencing.entry_event says"),
    ({"profile__panel_coverage": "none"}, "none exactly when there is no panels_file"),
    ({"profile__panel_sizes": "given"}, "not_given exactly when there is no panel_sizes_file"),
])
def test_a_mapping_problem_is_named(tmp_path, changes, problem):
    path = fx.write_mapping(tmp_path, broken(**changes))
    with pytest.raises(MappingError) as e:
        omop.load_mapping(path)
    assert problem in str(e.value)


# ---- Reading: the fixture comes back as it went in ----
def test_the_patients_come_back_with_their_fields_and_outcomes(genie, read_back):
    g = genie.tables.patients.set_index("patient_id")
    o = read_back.tables.patients.set_index("patient_id").reindex(g.index)
    assert set(read_back.tables.patients["patient_id"]) == set(g.index)
    for c in g.columns:
        a, b = g[c], o[c]
        if c == "days_to_death_or_last_followup":        # dates are whole days: compared as numbers, rounded down
            a = pd.to_numeric(a).map(lambda v: float(math.floor(v)) if pd.notna(v) else v).astype("Float64")
            b = pd.to_numeric(b).astype("Float64")
            assert a.isna().equals(b.isna()) and (a[a.notna()] == b[b.notna()]).all(), c
            continue
        if c == "year_of_diagnosis":                      # a date always has a year; compare where it was known
            a, b = a[a.notna()], b[a.notna()]
        assert (text(a) == text(b)).all(), c


def test_the_tests_results_treatments_and_panels_come_back(genie, source, read_back):
    _, ids, _, _ = source
    g, o = genie.tables, read_back.tables
    gt = g.tests.assign(test_id=g.tests["test_id"].map(ids)).set_index("test_id")
    ot = o.tests.set_index("test_id").reindex(gt.index)
    assert len(o.tests) == len(gt)
    for c in ("patient_id", "modality", "entry_event", "days_to_entry_event", "panel_id"):
        assert (text(gt[c]) == text(ot[c])).all(), c
    key = lambda df: sorted(map(tuple, df[["test_id", "kind", "gene", "variant", "variant_class", "value"]]
                                .astype("string").fillna("<NA>").values.tolist()))
    assert key(g.results.assign(test_id=g.results["test_id"].map(ids))) == key(o.results)
    days = lambda t: t.groupby("patient_id")["days_to_start"].apply(lambda s: sorted(s.dropna().tolist())).to_dict()
    assert days(g.treatments) == days(o.treatments)
    assert (o.treatments.groupby("patient_id")["regimen_order"].apply(lambda s: s.tolist() == list(range(1, len(s) + 1)))).all()
    assert g.panels.astype("string").equals(o.panels.astype("string"))


def test_after_follow_up_end_is_worked_out_from_the_dates(genie, source, read_back):
    """Where both dates are known it agrees with the source's own flag; where either is missing it is blank."""
    _, ids, _, _ = source
    g = genie.tables.tests.assign(test_id=genie.tables.tests["test_id"].map(ids)).set_index("test_id")
    o = read_back.tables.tests.set_index("test_id").reindex(g.index)
    both = g["after_follow_up_end"].notna() & o["after_follow_up_end"].notna()
    assert both.any() and (g.loc[both, "after_follow_up_end"] == o.loc[both, "after_follow_up_end"]).all()
    assert o.loc[o["days_to_entry_event"].isna(), "after_follow_up_end"].isna().all()


def test_what_reading_found_is_counted_and_masked(read_back, source):
    base = source[0]
    notes = "\n".join(read_back.read_notes)
    assert "<11 variant rows are not on a patient's sequencing test" in notes      # one decoy
    assert "have no date" in notes and "neither a death nor an observation period" in notes
    assert read_back.source_steps == [("Has a condition in the concept set 'colorectal cancer'", 0)]
    files = {p.name for p in read_back.input_files}
    assert {"mapping.yaml", "person.csv", "death.csv", "variant_occurrence.csv", "panels.csv"} <= files
    assert all(Path(p).is_relative_to(base) for p in read_back.input_files)
    assert read_back.profile.adapter == "omop"


def test_a_person_without_the_diagnosis_is_left_out_at_the_first_step(source, tmp_path):
    db = copy_of(source, tmp_path)
    cond = pd.read_csv(db / "cdm" / "condition_occurrence.csv", dtype=str)
    first = cond["person_id"].iloc[0]
    cond[cond["person_id"] != first].to_csv(db / "cdm" / "condition_occurrence.csv", index=False)
    res = omop.read(db / "mapping.yaml", source[3])
    assert res.source_steps == [("Has a condition in the concept set 'colorectal cancer'", 1)]
    assert len(res.tables.patients) == len(pd.read_csv(db / "cdm" / "person.csv")) - 1


@pytest.mark.parametrize("form", ["parquet", "duckdb"])
def test_parquet_files_and_a_duckdb_file_read_the_same(source, read_back, tmp_path, form):
    """Typed columns (whole-number IDs, real dates) read as the text files do."""
    db = copy_of(source, tmp_path)
    con = duckdb.connect(str(db / "omop.duckdb") if form == "duckdb" else ":memory:")
    if form == "duckdb":
        con.execute("CREATE SCHEMA cdm")
    for csv in sorted((db / "cdm").glob("*.csv")):
        q = f"SELECT * FROM read_csv_auto('{csv.as_posix()}')"
        if form == "duckdb":
            con.execute(f"CREATE TABLE cdm.{csv.stem} AS {q}")
        else:
            con.execute(f"COPY ({q}) TO '{csv.with_suffix('.parquet').as_posix()}' (FORMAT parquet)")
            csv.unlink()
    con.close()
    cdm = {"duckdb": "omop.duckdb", "schema": "cdm"} if form == "duckdb" else {"folder": "cdm", "date_format": "%Y-%m-%d"}
    res = omop.read(fx.write_mapping(db, fx.mapping(cdm=cdm)), source[3])
    for name in ("patients", "tests", "results", "treatments"):
        a = read_back.tables.table(name).astype("string").fillna("<NA>")
        b = res.tables.table(name).astype("string").fillna("<NA>")
        assert a.equals(b), name
    if form == "duckdb":
        assert {p.name for p in res.input_files} == {"mapping.yaml", "omop.duckdb", "panels.csv"}


def test_fusions_and_copy_numbers_are_read_from_the_variant_table(source, tmp_path):
    db = copy_of(source, tmp_path)
    path = db / "cdm" / "variant_occurrence.csv"
    v = pd.read_csv(path, dtype=str)
    test = v["procedure_occurrence_id"].iloc[0]
    v["target_gene2_symbol"], v["copy_number"] = pd.NA, pd.NA
    extra = pd.DataFrame([{"variant_occurrence_id": "90001", "procedure_occurrence_id": test, "target_gene1_symbol": "EML4",
                           "target_gene2_symbol": "ALK"},
                          {"variant_occurrence_id": "90002", "procedure_occurrence_id": test, "target_gene1_symbol": "ERBB2",
                           "copy_number": "8"}])
    pd.concat([v, extra], ignore_index=True).to_csv(path, index=False)
    m = fx.mapping()
    m["sequencing"]["variants"].update(partner_gene="target_gene2_symbol", copy_number="copy_number")
    r = omop.read(fx.write_mapping(db, m), source[3]).tables.results
    on = r[r["test_id"] == f"procedure:{test}"]
    fusion = on[on["kind"] == "fusion"]
    assert sorted(fusion["gene"]) == ["ALK", "EML4"] and set(fusion["variant"]) == {"EML4-ALK"}
    cn = on[on["kind"] == "copy_number"]
    assert cn["gene"].tolist() == ["ERBB2"] and cn["value"].tolist() == ["8"]


def test_a_concept_set_can_take_in_its_descendants(source, tmp_path):
    db = copy_of(source, tmp_path)
    parent = 2000000900
    m = fx.mapping()
    m["concept_sets"]["colorectal cancer"] = {"concepts": [parent], "include_descendants": True}
    with pytest.raises(MappingError, match="no concept_ancestor table"):
        omop.read(fx.write_mapping(db, m), source[3])
    pd.DataFrame([{"ancestor_concept_id": parent, "descendant_concept_id": fx.DIAGNOSIS}]).to_csv(
        db / "cdm" / "concept_ancestor.csv", index=False)
    res = omop.read(db / "mapping.yaml", source[3])
    assert res.source_steps[0][1] == 0 and len(res.tables.patients) > 0


def test_labels_must_list_every_value_found(source, tmp_path):
    db = copy_of(source, tmp_path)
    m = fx.mapping()
    del m["patients"]["histology"]["labels"]["8255/3"]
    with pytest.raises(MappingError, match="patients.histology: values the labels don't list: 8255/3"):
        omop.read(fx.write_mapping(db, m), source[3])
    m["patients"]["histology"]["labels"]["8255/3"] = None                 # listed as blank
    res = omop.read(fx.write_mapping(db, m), source[3])
    assert "Other histologies/mixed tumor" not in set(res.tables.patients["histology"].dropna())


def test_same_day_drugs_form_one_regimen_named_by_concept(source, tmp_path):
    db = copy_of(source, tmp_path)
    path = db / "cdm" / "drug_exposure.csv"
    d = pd.read_csv(path, dtype=str)
    first = d.iloc[0]
    twin = first.copy()
    twin["drug_exposure_id"], twin["drug_concept_id"] = "99999", str(fx.CANCER_DRUG + 10)
    pd.concat([d, twin.to_frame().T], ignore_index=True).to_csv(path, index=False)
    pd.DataFrame([{"concept_id": fx.CANCER_DRUG, "concept_name": "oxaliplatin"},
                  {"concept_id": fx.CANCER_DRUG + 10, "concept_name": "fluorouracil"}]).to_csv(db / "cdm" / "concept.csv", index=False)
    m = fx.mapping(treatments={"drugs": "cancer drugs", "name": "concept_name"})
    m["concept_sets"]["cancer drugs"]["concepts"].append(fx.CANCER_DRUG + 10)
    tr = omop.read(fx.write_mapping(db, m), source[3]).tables.treatments
    person = pd.read_csv(db / "cdm" / "person.csv", dtype=str).set_index("person_id")["person_source_value"]
    mine = tr[tr["patient_id"] == person[first["person_id"]]]
    assert mine["regimen"].iloc[0] == "fluorouracil, oxaliplatin"
    assert mine["regimen_order"].tolist() == list(range(1, len(mine) + 1))


# ---- Screening and running ----
def test_the_omop_route_screens_the_same_as_genie_bpc(genie, read_back, source):
    g = screen(outcome_blind(genie.tables), crc_settings(), genie.profile, genie.source_steps, genie.input_files)
    o = screen(outcome_blind(read_back.tables), source[3], read_back.profile, read_back.source_steps,
               read_back.input_files)
    assert o.decisions == g.decisions and o.decisions


@needs_r
def test_a_whole_run_on_the_omop_database(source, tmp_path):
    db = copy_of(source, tmp_path)
    settings_path = fx.omop_settings(db, resamples=10)
    settings = load_settings(settings_path)
    res = omop.read(db / "mapping.yaml", settings)
    lock_plan(screen(outcome_blind(res.tables), settings, res.profile, res.source_steps, res.input_files), settings,
              res.input_files, tmp_path / "work", "A. Reviewer", read_notes=res.read_notes)
    summary = run(settings_path, tmp_path / "work", tmp_path / "out", mapping=db / "mapping.yaml",
                  r=ROptions(R, os.environ.get("DELAYED_ENTRY_R_LIBRARY")), log=lambda *_: None)
    assert all(v["verdict"] in ("pass", "warning", "applied") for v in summary["checks"])
    r_steps = [c for c in summary["checks"] if c["stage"] == "R"]
    assert r_steps and all(c["verdict"] == "pass" for c in r_steps)
    assert {"counts.csv", "estimates.csv", "bias.csv", "run_report.md"} <= set(summary["files"])
    bias = pd.read_csv(tmp_path / "out" / "bias.csv", keep_default_na=False)
    assert {"primary", "3a", "3b", "S1", "S2", "SA1"} <= set(bias["analysis"])
    report = (tmp_path / "out" / "run_report.md").read_text(encoding="utf-8")
    assert "omop" in report.lower()
    # the made-up mapping's profile gives no attribution: no section, no check
    assert "Data source and attribution" not in report and not [c for c in summary["checks"] if c["stage"] == "source"]
