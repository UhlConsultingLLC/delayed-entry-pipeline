"""Loading the GENIE BPC files (adapters/genie_bpc/raw.py) and the link checks. Ported from the CRC BPC
project's tests/test_db.py, on the same made-up fixture."""
import csv
import shutil

import duckdb
import pytest

from delayed_entry.adapters.genie_bpc import raw
from helpers import FIXTURE

COHORT = "crc_v2.0"
CHECKS = {"patient_id_unique", "patient_has_cancer", "cancer_has_patient", "cancer_key_unique", "report_has_cancer",
          "report_sample_unique", "report_number_unique", "report_has_panel", "mutation_has_report",
          "call_gene_on_panel", "fusion_has_report", "other_cancer_has_patient", "other_cancer_not_index", "other_cancer_flag_blank",
          "pathology_has_patient", "pathology_key_unique", "regimen_number_unique", "report_seq_year_is_a_year"}


@pytest.fixture(scope="module")
def con():
    return raw.load(COHORT, FIXTURE)[0]


def file_rows(path):
    delim = "\t" if path.suffix == ".txt" else ","
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter=delim))


def rewrite(path, rows):
    delim = "\t" if path.suffix == ".txt" else ","
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter=delim, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def test_every_row_of_every_file_is_loaded(con):
    for table, (name, _, _) in raw.raw_files(COHORT).items():
        assert raw.scalar(con, f"SELECT count(*) FROM {table}") == len(file_rows(FIXTURE / name)), table
    genes_in_files = sum(len(raw.read_panel_file(p)[1]) for p in FIXTURE.glob(raw.PANEL_FILES))
    assert raw.scalar(con, "SELECT count(*) FROM raw_panel_gene") == genes_in_files


def test_load_lists_every_file_it_read():
    _, files = raw.load(COHORT, FIXTURE)
    assert {p.name for p in files} == {n for n, _, _ in raw.raw_files(COHORT).values()} | \
           {p.name for p in FIXTURE.glob(raw.PANEL_FILES)}


def test_the_made_up_truth_file_is_never_loaded(con):
    assert "raw_truth" not in {r[0] for r in con.execute("SHOW TABLES").fetchall()}


def test_typed_columns_keep_their_type_and_blanks_become_null(con):
    types = dict(con.execute("SELECT column_name, column_type FROM (DESCRIBE raw_patient)").fetchall())
    assert types["redacted"] == "VARCHAR"        # left to guess, DuckDB reads Yes/No as true/false
    assert types["birth_year"] == "INTEGER"
    assert {r[0] for r in con.execute("SELECT DISTINCT redacted FROM raw_patient").fetchall()} == {"Yes", "No"}
    blanks = sum(r["dx_cpt_rep_days"] == "" for r in file_rows(FIXTURE / "cancer_panel_test_level_dataset.csv"))
    assert blanks > 0
    assert raw.scalar(con, "SELECT count(*) FROM raw_report WHERE dx_cpt_rep_days IS NULL") == blanks


def test_panel_files_become_one_row_per_gene():
    panel, genes = raw.read_panel_file(FIXTURE / "data_gene_panel_VICC-PANEL-1.txt")
    assert panel == "VICC-PANEL-1"
    assert genes == ["APC", "BRAF", "KRAS", "PIK3CA", "SMAD4", "TP53"]   # no NRAS on this panel


def test_every_link_check_passes_on_the_fixture(con):
    results = raw.run_checks(con)
    assert set(results) == CHECKS
    assert all(n == 0 for n in results.values()), results


def test_link_checks_catch_broken_links(tmp_path):
    """Break eight things on purpose; each check must find exactly its own problem."""
    broken = shutil.copytree(FIXTURE, tmp_path / "broken")
    reports = file_rows(broken / "cancer_panel_test_level_dataset.csv")
    reports[0]["ca_seq"] = "5"                                          # a cancer this patient doesn't have
    reports.append({**reports[1], "cpt_genie_sample_id": reports[1]["cpt_genie_sample_id"] + "X"})  # a test number used twice
    reports[2]["cpt_seq_date"] = "2016-05"                             # a year of sequencing that isn't a year
    rewrite(broken / "cancer_panel_test_level_dataset.csv", reports)
    calls = file_rows(broken / "data_mutations_extended.txt")
    calls.append({**calls[0], "Tumor_Sample_Barcode": "GENIE-DFCI-S99999-T01"})   # a sample with no report
    calls.append({**calls[1], "Hugo_Symbol": "ZZZ1"})                  # a gene no panel reads
    rewrite(broken / "data_mutations_extended.txt", calls)
    other = file_rows(broken / "cancer_level_dataset_non_index.csv")
    other[0]["ca_seq"] = "1"                                            # the same number as the index cancer
    rewrite(broken / "cancer_level_dataset_non_index.csv", other)
    path = file_rows(broken / "pathology_report_level_dataset.csv")
    path.append(dict(path[0]))                                          # one pathology report listed twice
    rewrite(broken / "pathology_report_level_dataset.csv", path)
    regimens = file_rows(broken / "regimen_cancer_level_dataset.csv")
    regimens.append(dict(regimens[0]))                                  # one regimen listed twice
    rewrite(broken / "regimen_cancer_level_dataset.csv", regimens)

    results = raw.run_checks(raw.load(COHORT, broken)[0])
    expected = dict.fromkeys(CHECKS, 0) | {"report_has_cancer": 1, "report_number_unique": 1,
                                           "mutation_has_report": 1, "call_gene_on_panel": 1,
                                           "other_cancer_not_index": 1, "pathology_key_unique": 1,
                                           "regimen_number_unique": 1, "report_seq_year_is_a_year": 1}
    assert results == expected


@pytest.mark.parametrize("flag, not_index, blank", [("", 0, 1), ("Unknown", 1, 0), ("Yes", 1, 0), ("No", 0, 0)])
def test_the_non_index_flag_checks(tmp_path, flag, not_index, blank):
    """A blank index flag is its own review check; any other flag but "No" still fails."""
    broken = shutil.copytree(FIXTURE, tmp_path / "flag")
    other = file_rows(broken / "cancer_level_dataset_non_index.csv")
    other[0]["redcap_ca_index"] = flag
    rewrite(broken / "cancer_level_dataset_non_index.csv", other)
    results = raw.run_checks(raw.load(COHORT, broken)[0])
    assert (results["other_cancer_not_index"], results["other_cancer_flag_blank"]) == (not_index, blank)


def test_a_renamed_column_stops_the_load(tmp_path):
    broken = shutil.copytree(FIXTURE, tmp_path / "renamed")
    path = broken / "cancer_panel_test_level_dataset.csv"
    path.write_text(path.read_text(encoding="utf-8").replace("dx_cpt_rep_days", "dx_cpt_report_days", 1),
                    encoding="utf-8", newline="\n")
    with pytest.raises(duckdb.Error, match="dx_cpt_rep_days"):
        raw.load(COHORT, broken)


def test_clinical_and_genomic_files_can_sit_in_separate_folders(tmp_path):
    clinical, genomic = tmp_path / "clinical", tmp_path / "genomic"
    clinical.mkdir()
    genomic.mkdir()
    for p in FIXTURE.iterdir():
        if p.name not in ("truth.csv", "README.md"):
            shutil.copy(p, (genomic if p.name.startswith("data_") else clinical) / p.name)
    con, _ = raw.load(COHORT, clinical, genomic)
    assert all(n == 0 for n in raw.run_checks(con).values())
    with pytest.raises(FileNotFoundError):
        raw.load(COHORT, clinical)            # the genomic files aren't there


def test_the_text_na_is_a_missing_value_like_a_blank(tmp_path):
    """genieBPC's reader treats a blank and the text NA alike, so the loader must too."""
    shutil.copytree(FIXTURE, tmp_path / "data")
    path = tmp_path / "data" / "cancer_level_dataset_index.csv"
    rows = file_rows(path)
    rows[0]["tt_os_dx_days"] = "NA"                 # a typed number
    rows[0]["crc_type"] = "NA"                      # a typed text
    rewrite(path, rows)
    con, _ = raw.load(COHORT, tmp_path / "data")
    first = con.execute("SELECT tt_os_dx_days, crc_type FROM raw_cancer WHERE record_id = ? AND ca_seq = ?",
                        [rows[0]["record_id"], int(rows[0]["ca_seq"])]).fetchall()[0]
    assert first == (None, None)


def test_a_guessed_column_that_turns_to_text_deep_in_a_file_does_not_stop_the_load(tmp_path):
    shutil.copytree(FIXTURE, tmp_path / "data")
    path = tmp_path / "data" / "data_mutations_extended.txt"
    rows = [dict(r, guess_me="1") for r in file_rows(path)]
    rows = (rows * (30000 // len(rows) + 1))[:30000]
    rows.append(dict(rows[0], guess_me="not a number"))
    rewrite(path, rows)
    con, _ = raw.load(COHORT, tmp_path / "data")
    assert raw.scalar(con, "SELECT count(*) FROM raw_mutation") == 30001


def test_an_unknown_cohort_is_named():
    with pytest.raises(ValueError, match="no GENIE BPC cohort 'nsclc_v9'"):
        raw.raw_files("nsclc_v9")
