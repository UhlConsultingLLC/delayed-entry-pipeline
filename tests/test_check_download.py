"""The check after downloading a GENIE BPC release. Ported from the CRC BPC project's
tests/test_check_download.py, on a made-up release built from the fixture. Two of that project's checks
(a blank stage IV histology, and cancers out of order in the file) served its genieBPC cross-check,
which this pipeline doesn't have, so their tests are not carried over."""
import csv
import hashlib
import shutil
from pathlib import Path

import pytest
import yaml

from delayed_entry.adapters.genie_bpc import check_download as cd
from delayed_entry.adapters.genie_bpc import raw
from delayed_entry.lock import lock_settings
from helpers import CRC_SETTINGS, FIXTURE, ROOT

GENOMIC_FILES = ("data_mutations_extended.txt",)
COHORT = "crc_v2.0"
CRC_LABELS = ("Colon Cancer", "Colon/Rectum Cancer", "Rectum and Rectosigmoid Cancer")


class Release:
    def __init__(self, root: Path):
        self.clinical, self.genomic, self.work = root / "clinical", root / "genomic", root / "work"
        for d in (self.clinical, self.genomic, self.work):
            d.mkdir()
        for p in FIXTURE.iterdir():
            if p.name in ("README.md", "truth.csv"):
                continue
            to = self.genomic if p.name in GENOMIC_FILES or p.name.startswith("data_gene_panel_") else self.clinical
            shutil.copy(p, to / p.name)
        for d in (self.clinical, self.genomic):
            (d / f"manifest_{d.name}_2026-11-10.csv").write_text("synapse_id\n", encoding="utf-8")
        self.lock = lock_settings(CRC_SETTINGS, root / "lock", locked_by="A. Reviewer")
        self.yml = root / "expected.yml"
        # The real non-index file has site codes and the made-up one doesn't: give each cancer one that fits its type.
        other = self.clinical / "cancer_level_dataset_non_index.csv"
        with open(other, encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f))
        for r in rows:
            r["ca_d_site"] = "C18.7" if r["ca_type"] in CRC_LABELS else "C50.9"
        with open(other, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        self.refresh()

    def refresh(self):
        real = yaml.safe_load((raw.cohort_dir(COHORT) / "expected_release.yml").read_text(encoding="utf-8"))
        files = []
        for folder, d in (("clinical", self.clinical), ("genomic", self.genomic)):
            for p in sorted(d.iterdir()):
                if not p.name.startswith("manifest_"):
                    files.append(dict(folder=folder, name=p.name, bytes=p.stat().st_size, md5=cd.md5(p), used=True))
        panels = [p.stem[len("data_gene_panel_"):] for p in sorted(self.genomic.glob(raw.PANEL_FILES))]
        self.yml.write_text(yaml.safe_dump(dict(files=files, panels=panels, values=real["values"])), encoding="utf-8")

    def run(self, **kw):
        kw.setdefault("settings_lock", self.lock)
        return cd.run_checks(self.clinical, self.genomic, self.work, CRC_SETTINGS, expected_yml=self.yml, **kw)

    def check(self, name: str):
        return next(c for c in self.run() if c.name.startswith(name))

    def status(self, name: str) -> str:
        return self.check(name).status

    def detail(self, name: str) -> str:
        return self.check(name).detail

    def path(self, name: str) -> Path:
        return (self.genomic if name.startswith("data_") else self.clinical) / name

    def set_cell(self, name: str, column: str, value: str, row: int = 0, every: bool = False):
        p = self.path(name)
        delim = "\t" if p.suffix == ".txt" else ","
        with open(p, encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f, delimiter=delim))
        for i, r in enumerate(rows):
            if every or i == row:
                r[column] = value
        with open(p, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter=delim, lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        self.refresh()

    def edit_text(self, name: str, old: str, new: str, count: int = -1):
        p = self.path(name)
        p.write_text(p.read_text(encoding="utf-8").replace(old, new, count), encoding="utf-8", newline="")
        self.refresh()


PATIENT, CANCER, REPORT, MUTATION, PATHOLOGY = ("patient_level_dataset.csv", "cancer_level_dataset_index.csv",
                                                "cancer_panel_test_level_dataset.csv", "data_mutations_extended.txt",
                                                "pathology_report_level_dataset.csv")
OTHER = "cancer_level_dataset_non_index.csv"
SITES = "Non-index cancers at this cancer's site codes"


@pytest.fixture
def release(tmp_path):
    return Release(tmp_path)


# ---- A good release ----
def test_a_good_release_has_no_failure(release):
    checks = release.run()
    assert [c for c in checks if c.status == cd.FAIL] == []
    assert cd.verdict(checks).startswith("PASSED")


def test_the_rehearsal_on_the_made_up_data_has_no_failure(tmp_path):
    checks = cd.run_checks(FIXTURE, FIXTURE, tmp_path, CRC_SETTINGS, rehearsal=True)
    assert [c for c in checks if c.status == cd.FAIL] == []
    skipped = {c.name for c in checks if c.status == cd.SKIPPED}
    assert "The settings were locked before the download" in skipped
    assert "Size and checksum of each file match Synapse" in skipped
    site = next(c for c in checks if c.name.startswith(SITES))
    assert (site.status, site.detail) == (cd.SKIPPED, "the made-up file has no site codes")


def test_the_list_of_the_real_release_covers_every_file_the_adapter_reads():
    real = yaml.safe_load((raw.cohort_dir(COHORT) / "expected_release.yml").read_text(encoding="utf-8"))
    used = {f["name"]: f for f in real["files"] if f["used"]}
    for name, folder, _ in raw.raw_files(COHORT).values():
        assert used[name]["folder"] == folder
    assert {n for n in used if n.startswith("data_gene_panel_")} == {f"data_gene_panel_{p}.txt" for p in real["panels"]}
    assert all(len(f["md5"]) == 32 and f["bytes"] > 0 for f in real["files"])
    assert set(real["checks"]) == {"min_call_link", "max_unit_mismatch", "report_months_tolerance",
                                   "age_years_tolerance", "report_days_tolerance", "age_range", "birth_year_range"}


def test_every_link_check_has_a_stated_meaning():
    assert set(cd.LINK_CHECKS) == set(raw.read_checks())


def test_the_labels_searched_for_come_from_the_settings():
    from helpers import crc_settings
    labels = cd.searched_labels(crc_settings())
    assert labels["genes"] == ["BRAF", "KRAS", "NRAS"] and labels["variants"] == ["p.V600E"]
    assert labels["cancer_labels"] == list(CRC_LABELS)
    assert ("stage_at_diagnosis", "trim", ["Stage IV"]) in labels["fields"]


# ---- This cancer's labels and site codes ----
def test_a_label_of_this_cancer_spelled_otherwise_in_the_non_index_file_is_flagged(release):
    name = "Other primaries of this cancer in the non-index file"
    assert release.status(name) == cd.PASS
    release.set_cell(OTHER, "ca_type", "Colorectal Cancer")
    assert release.status(name) == cd.REVIEW


def test_a_site_of_this_cancer_with_a_blank_or_other_type_fails(release):
    assert release.status(SITES) == cd.PASS
    release.set_cell(OTHER, "ca_d_site", "C20.9")
    release.set_cell(OTHER, "ca_type", "")
    assert release.status(SITES) == cd.FAIL
    assert release.detail(SITES) == "<11 cancers: the exclusion would keep these patients"
    release.set_cell(OTHER, "ca_type", "Other")
    assert release.status(SITES) == cd.FAIL


def test_the_appendix_is_not_a_site_of_this_cancer(release):
    release.set_cell(OTHER, "ca_d_site", "C18.1")
    release.set_cell(OTHER, "ca_type", "Appendix Cancer")
    assert release.status(SITES) == cd.PASS


def test_a_blank_index_flag_in_the_non_index_file_is_read_not_failed(release):
    release.set_cell(OTHER, "redcap_ca_index", "")
    assert release.status("other_cancer_flag_blank") == cd.REVIEW
    assert release.status("other_cancer_not_index") == cd.PASS
    release.set_cell(OTHER, "redcap_ca_index", "Yes")
    assert release.status("other_cancer_not_index") == cd.FAIL


# ---- Folders and safety ----
def test_settings_not_locked_before_the_download_stop_it(release):
    assert release.status("The settings were locked before the download") == cd.PASS
    checks = release.run(settings_lock=None)
    assert next(c.status for c in checks if c.name == "The settings were locked before the download") == cd.FAIL


def test_settings_changed_since_the_lock_stop_it(release, tmp_path):
    copy = tmp_path / "copy.yaml"
    copy.write_bytes(CRC_SETTINGS.read_bytes())
    lock = lock_settings(copy, tmp_path / "lock2", locked_by="A. Reviewer")
    copy.write_text(copy.read_text(encoding="utf-8").replace("value: 730", "value: 731"), encoding="utf-8")
    checks = cd.run_checks(release.clinical, release.genomic, release.work, copy, lock, expected_yml=release.yml)
    assert next(c for c in checks if c.name == "The settings were locked before the download").status == cd.FAIL


def test_made_up_data_in_the_real_folders_stops_it(release):
    shutil.copy(FIXTURE / "truth.csv", release.clinical / "truth.csv")
    assert release.status("No made-up data") == cd.FAIL


def test_a_cloud_synced_folder_stops_it(release, tmp_path):
    cloud = tmp_path / "OneDrive" / "clinical"
    shutil.copytree(release.clinical, cloud)
    checks = cd.run_checks(cloud, release.genomic, release.work, CRC_SETTINGS, release.lock, expected_yml=release.yml)
    assert next(c.status for c in checks if c.name == "The clinical folder is not cloud-synced") == cd.FAIL


def test_a_work_folder_inside_the_package_stops_it(release):
    checks = cd.check_folders(release.clinical, release.genomic, ROOT / "work", False, CRC_SETTINGS, release.lock)
    assert next(c.status for c in checks if c.name == "The work folder is outside the package folder") == cd.FAIL


def test_a_missing_folder_stops_it_without_crashing(release):
    shutil.rmtree(release.genomic)
    checks = release.run()
    assert next(c.status for c in checks if c.name == "The genomic folder exists") == cd.FAIL


# ---- The files themselves ----
def test_a_missing_file_stops_it(release):
    release.path(PATIENT).unlink()
    assert release.status("Every file the adapter reads is there") == cd.FAIL


def test_a_missing_file_the_adapter_does_not_read_is_only_reported(release):
    y = yaml.safe_load(release.yml.read_text(encoding="utf-8"))
    y["files"].append(dict(folder="clinical", name="imaging_level_dataset.csv", bytes=1, md5="0" * 32, used=False))
    release.yml.write_text(yaml.safe_dump(y), encoding="utf-8")
    assert release.status("Every other file in the release is there") == cd.REVIEW
    assert cd.verdict(release.run()).startswith("PASSED, with")


def test_a_truncated_download_fails_the_size_check(release):
    p = release.path(CANCER)
    p.write_bytes(p.read_bytes()[:-500])
    assert release.status("Every file the adapter reads has the size") == cd.FAIL


def test_a_changed_file_of_the_same_size_fails_the_checksum(release):
    p = release.path(CANCER)
    p.write_bytes(p.read_bytes().replace(b"Stage IV", b"Stage XX", 1))
    assert release.status("Every file the adapter reads has the checksum") == cd.FAIL


def test_a_failed_file_check_holds_back_the_deeper_checks(release):
    p = release.path(CANCER)
    p.write_bytes(p.read_bytes()[:-500])
    assert release.status("The links, units") == cd.SKIPPED


def test_a_missing_manifest_stops_it(release):
    next(release.clinical.glob("manifest_*.csv")).unlink()
    assert release.status("The clinical download left its checksum manifest") == cd.FAIL


def test_a_panel_coverage_file_stops_it(release):
    (release.genomic / "genomic_information.txt").write_text("x\n", encoding="utf-8")
    release.refresh()
    assert release.status("No panel coverage or size file") == cd.FAIL


# ---- Reading the files ----
def test_a_renamed_column_is_named(release):
    release.edit_text(CANCER, "stage_dx", "stage_at_dx", 1)
    assert release.status("Every column the adapter uses") == cd.FAIL
    assert "stage_dx" in release.detail("Every column the adapter uses")


def test_a_withheld_value_written_as_text_is_caught_and_shown_only_by_its_form(release):
    release.set_cell(REPORT, "dx_cpt_rep_days", ">32485")
    assert release.status("Every value converts") == cd.FAIL
    detail = release.detail("Every value converts")
    assert "dx_cpt_rep_days" in detail and ">99999" in detail and "32485" not in detail and "<11" in detail


def test_the_text_na_is_read_as_missing_and_reported(release):
    release.set_cell(CANCER, "tt_os_dx_days", "NA")
    assert release.status("No column uses the text NA") == cd.REVIEW
    assert release.status("Every value converts") == cd.PASS
    assert release.status("The adapter can load every file") == cd.PASS


def test_a_file_that_is_not_utf8_is_caught_without_printing_it(release):
    p = release.path(PATHOLOGY)
    p.write_bytes(p.read_bytes().replace(b"Colon Cancer", b"Colon Canc\xe9r", 1))
    release.refresh()
    assert release.status("Every file is valid UTF-8") == cd.FAIL


def test_a_header_only_file_is_caught(release):
    p = release.path(PATIENT)
    p.write_text(p.read_text(encoding="utf-8").split("\n")[0] + "\n", encoding="utf-8", newline="")
    release.refresh()
    assert release.status("No file is empty") == cd.FAIL


# ---- Gene panels ----
def test_genes_separated_by_spaces_are_caught(release):
    p = next(release.genomic.glob(raw.PANEL_FILES))
    p.write_text(p.read_text(encoding="utf-8").replace("\t", " "), encoding="utf-8", newline="")
    release.refresh()
    assert release.status("Each panel file names") == cd.FAIL


def test_a_panel_name_that_is_not_the_file_name_is_caught(release):
    p = next(release.genomic.glob(raw.PANEL_FILES))
    p.write_text(p.read_text(encoding="utf-8").replace("stable_id: ", "stable_id: X", 1), encoding="utf-8", newline="")
    release.refresh()
    assert release.status("Each panel file names") == cd.FAIL


def test_no_panel_reading_every_gene_the_groups_need_stops_it(release):
    for p in release.genomic.glob(raw.PANEL_FILES):
        p.write_text(p.read_text(encoding="utf-8").replace("KRAS", "KRASX"), encoding="utf-8", newline="")
    release.refresh()
    assert release.status("At least one panel reads every gene") == cd.FAIL


def test_a_panel_that_is_not_in_the_release_is_named(release):
    y = yaml.safe_load(release.yml.read_text(encoding="utf-8"))
    y["panels"].append("MADE-UP-9")
    release.yml.write_text(yaml.safe_dump(y), encoding="utf-8")
    assert release.status("The panels are the ones") == cd.FAIL


# ---- Coded values and the labels the rules search for ----
def test_a_stage_label_not_in_the_documentation_stops_it(release):
    release.set_cell(CANCER, "stage_dx", "Stage IVA")
    assert release.status("Coded columns the rules read") == cd.FAIL
    assert "Stage IVA" in release.detail("Coded columns the rules read")


def test_a_trailing_space_on_a_rule_label_passes_because_the_rule_trims(release):
    """The CRC project also needed the exact spelling for its genieBPC cross-check; here the settings' rule trims."""
    release.set_cell(CANCER, "stage_dx", "Stage IV ", every=True)
    assert release.status("The labels a cohort rule matches on stage_at_diagnosis") == cd.PASS


def test_a_rule_label_absent_from_the_file_stops_it(release):
    release.edit_text(CANCER, "Stage IV", "Stage 4")
    assert release.status("The labels a cohort rule matches on stage_at_diagnosis") == cd.FAIL


def test_a_survival_status_other_than_0_or_1_stops_it(release):
    release.set_cell(CANCER, "os_dx_status", "2")
    assert release.status("Coded columns the rules read") == cd.FAIL


def test_an_unknown_label_in_a_descriptive_column_is_only_reported(release):
    release.set_cell(PATIENT, "race_ethnicity", "Something new")
    assert release.status("Descriptive columns") == cd.REVIEW
    assert release.status("Coded columns the rules read") == cd.PASS


def test_an_msi_label_not_documented_stops_it(release):
    release.set_cell(PATHOLOGY, "msi_result", "MSI high")
    assert release.status("Coded columns the rules read") == cd.FAIL


def test_a_differently_written_protein_change_stops_it(release):
    release.edit_text(MUTATION, "p.V600E", "V600E")
    assert release.status("The protein changes the groups name") == cd.FAIL


def test_a_renamed_gene_stops_it(release):
    release.edit_text(MUTATION, "\nNRAS\t", "\nNRAS1\t")
    assert release.status("The genes the groups name appear") == cd.FAIL


# ---- Links, units and the trial build ----
def test_a_repeated_patient_stops_it(release):
    p = release.path(PATIENT)
    lines = p.read_text(encoding="utf-8").split("\n")
    p.write_text("\n".join(lines[:2] + [lines[1]] + lines[2:]), encoding="utf-8", newline="")
    release.refresh()
    assert release.status("patient_id_unique") == cd.FAIL


def test_sample_ids_that_do_not_match_the_mutation_file_stop_it(release):
    """The worst silent failure: with no calls linked, every tumor would be read as wild type."""
    release.edit_text(MUTATION, "GENIE-", "GNIE-")
    assert release.status("Most report samples have mutation calls") == cd.FAIL


def test_a_repeated_regimen_number_on_an_index_cancer_stops_it(release):
    p = release.path("regimen_cancer_level_dataset.csv")
    lines = p.read_text(encoding="utf-8").split("\n")
    p.write_text("\n".join(lines[:2] + [lines[1]] + lines[2:]), encoding="utf-8", newline="")
    release.refresh()
    assert release.status("An index cancer's regimens") == cd.FAIL


def rescale(release, name, column, factor):
    p = release.path(name)
    with open(p, encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        if r[column]:
            r[column] = str(float(r[column]) * factor)
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    release.refresh()


def test_report_times_in_months_are_caught(release):
    rescale(release, REPORT, "dx_cpt_rep_days", 1 / 30.4375)
    assert release.status("Report months match report days") == cd.FAIL


def test_age_in_days_is_caught(release):
    rescale(release, CANCER, "age_dx", 365.25)
    assert release.status("Age at diagnosis is in years") == cd.FAIL


def test_an_msi_label_is_also_stopped_by_the_trial_build(release):
    """If the value check were ever loosened, the adapter itself would still refuse the label."""
    p = release.path(PATHOLOGY)
    with open(p, encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    row = next(r for r in rows if r["msi_result"] and any(r.get(f"path_ca_type{i}") == "Colon Cancer" for i in range(1, 46)))
    row["msi_result"] = "MSI high"
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    release.refresh()
    y = yaml.safe_load(release.yml.read_text(encoding="utf-8"))
    del y["values"]["raw_pathology"]
    release.yml.write_text(yaml.safe_dump(y), encoding="utf-8")
    assert release.status("The standard tables build") == cd.FAIL


# ---- What the report may say ----
def test_small_counts_are_hidden():
    show = cd.Shower(11)
    assert [show(n) for n in (0, 1, 10, 11, 1234)] == ["0", "<11", "<11", "11", "1,234"]


def test_an_error_is_cut_down_to_its_type_and_column():
    e = ValueError('CSV Error on Line: 752\nOriginal Line: CRC,MADE-UP-S00457,MSK\nError when converting column "tt_os_dx_days".')
    assert cd.safe_error(e) == "ValueError in column tt_os_dx_days"


def test_the_report_of_a_broken_release_names_no_patient(release):
    release.set_cell(REPORT, "dx_cpt_rep_days", ">32485")
    release.set_cell(CANCER, "stage_dx", "Stage IVA")
    text = cd.report(release.run())
    assert "GENIE-" not in text and "32485" not in text
    assert text.isascii()


def test_a_failed_check_exits_with_an_error_and_leaves_the_report(release, monkeypatch):
    real_run = cd.run_checks
    monkeypatch.setattr(cd, "run_checks", lambda *a, **k: real_run(*a, **k, expected_yml=release.yml))
    release.path(PATIENT).unlink()
    with pytest.raises(SystemExit) as stop:
        cd.main([str(CRC_SETTINGS), str(release.clinical), str(release.genomic), str(release.work),
                 "--settings-lock", str(release.lock)])
    assert stop.value.code == 1
    assert "RESULT: STOP" in (release.work / cd.REPORT_FILE).read_text(encoding="utf-8")


def test_md5_is_the_standard_digest(tmp_path):
    p = tmp_path / "a.txt"
    p.write_bytes(b"hello")
    assert cd.md5(p) == hashlib.md5(b"hello").hexdigest()
