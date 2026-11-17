"""The adapter's NSCLC v2.0-public spec, on made-up files shaped like that release.

There is no made-up NSCLC release yet, so the made-up CRC files are reshaped the way the NSCLC files differ
(variable synopsis, 2026-11-13): no crc_type but ca_lung_cigarette; 25 specimen columns, no MSI/MMR; lung
labels; panels that read EGFR and ALK; the made-up BRAF V600E calls turned into EGFR L858R calls; and a made-up
fusion file in cBioPortal's format: EML4-ALK on three patients' samples (a row for each gene), a CD74-ROS1 on
one, and one EML4-ALK listed under EML4 only."""
import importlib.util
import shutil

import duckdb
import pandas as pd
import pytest

from delayed_entry.adapters import genie_bpc
from delayed_entry.adapters.genie_bpc import check_download as cd
from delayed_entry.adapters.genie_bpc import raw
from delayed_entry.cohort import build_cohort
from delayed_entry.screen import COMPARE, DESCRIBE, screen
from delayed_entry.settings import load_settings
from delayed_entry.tables import outcome_blind
from helpers import FIXTURE, ROOT

NSCLC_SETTINGS = ROOT / "settings" / "nsclc_bpc.yaml"
COLORECTAL = ["Colon Cancer", "Colon/Rectum Cancer", "Rectum and Rectosigmoid Cancer"]
SMOKING = ["Never used", "Current user", "Former user (quit <1 year)", "Former user (quit >1 year)",
           "Former user (unknown time)", "Unknown"]
ALK_PATIENTS = ["GENIE-DFCI-S00015", "GENIE-DFCI-S00016", "GENIE-DFCI-S00076"]    # made-up stage IV, no driver
PARTNER_ONLY = "GENIE-DFCI-S00087"      # its EML4-ALK has an EML4 row only
ROS1_PATIENT = "GENIE-DFCI-S00104"
FUSION_COLUMNS = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "Tumor_Sample_Barcode", "Fusion", "DNA_support",
                  "RNA_support", "Method", "Frame"]


def fusion_file(reports: pd.DataFrame) -> pd.DataFrame:
    """Made-up fusions on every sample of the chosen patients."""
    def rows(pid, genes, fusion):
        samples = reports.loc[reports["record_id"] == pid, "cpt_genie_sample_id"]
        return [dict(Hugo_Symbol=g, Entrez_Gene_Id="", Center="DFCI", Tumor_Sample_Barcode=s, Fusion=fusion,
                     DNA_support="yes", RNA_support="unknown", Method="NGS", Frame="in frame")
                for s in samples for g in genes]
    out = [r for pid in ALK_PATIENTS for r in rows(pid, ["EML4", "ALK"], "EML4-ALK fusion")]
    out += rows(PARTNER_ONLY, ["EML4"], "EML4-ALK fusion")
    out += rows(ROS1_PATIENT, ["CD74", "ROS1"], "CD74-ROS1 fusion")
    return pd.DataFrame(out, columns=FUSION_COLUMNS)


def read(path, sep=","):
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)


def write(df, path, sep=","):
    df.to_csv(path, sep=sep, index=False, lineterminator="\n")


@pytest.fixture(scope="module")
def nsclc_files(tmp_path_factory):
    d = tmp_path_factory.mktemp("nsclc")
    for p in FIXTURE.iterdir():
        if p.name not in ("truth.csv", "README.md"):
            shutil.copy(p, d / p.name)
    ca = read(d / "cancer_level_dataset_index.csv").drop(columns=["crc_type"])
    ca["ca_lung_cigarette"] = [SMOKING[i % len(SMOKING)] for i in range(len(ca))]
    write(ca, d / "cancer_level_dataset_index.csv")
    path = read(d / "pathology_report_level_dataset.csv")
    path = path.drop(columns=[c for c in path.columns if c.startswith(("msi_", "mmr_"))]
                     + [f"path_ca_type{i}" for i in range(26, 46)])
    for i in range(1, 26):
        path[f"path_ca_type{i}"] = path[f"path_ca_type{i}"].replace(COLORECTAL, "Non Small Cell Lung Cancer")
    write(path, d / "pathology_report_level_dataset.csv")
    other = read(d / "cancer_level_dataset_non_index.csv")
    other["ca_type"] = other["ca_type"].replace(COLORECTAL, "Non Small Cell Lung Cancer")
    write(other, d / "cancer_level_dataset_non_index.csv")
    for p in d.glob(raw.PANEL_FILES):
        p.write_text(p.read_text(encoding="utf-8").rstrip("\n") + "\tEGFR\tALK\n", encoding="utf-8")
    maf = read(d / "data_mutations_extended.txt", "\t")
    v600e = (maf["Hugo_Symbol"] == "BRAF") & (maf["HGVSp_Short"].str.strip() == "p.V600E")
    maf.loc[v600e, ["Hugo_Symbol", "HGVSp_Short"]] = ["EGFR", "p.L858R"]
    write(maf, d / "data_mutations_extended.txt", "\t")
    write(fusion_file(read(d / "cancer_panel_test_level_dataset.csv")), d / "data_fusions.txt", "\t")
    return d


@pytest.fixture(scope="module")
def nsclc():
    return load_settings(NSCLC_SETTINGS)


def test_the_nsclc_spec_lists_the_release_files_with_their_columns():
    files = raw.raw_files("nsclc_v2.0")
    assert {name for name, _, _ in files.values()} == {
        "patient_level_dataset.csv", "cancer_level_dataset_index.csv", "cancer_panel_test_level_dataset.csv",
        "data_mutations_extended.txt", "data_fusions.txt", "cancer_level_dataset_non_index.csv",
        "pathology_report_level_dataset.csv", "regimen_cancer_level_dataset.csv"}
    assert set(files["raw_fusion"][2]) == {"Hugo_Symbol", "Tumor_Sample_Barcode", "Fusion"}
    assert files["raw_fusion"][1] == "genomic" and "regimen_drugs" in files["raw_regimen"][2]
    assert "ca_lung_cigarette" in files["raw_cancer"][2] and "crc_type" not in files["raw_cancer"][2]
    assert not [c for c in files["raw_pathology"][2] if c.startswith(("msi", "mmr"))]
    assert sum(c.startswith("path_ca_type") for c in files["raw_pathology"][2]) == 25


def test_nsclc_shaped_files_load_and_link(nsclc_files):
    con, _ = raw.load("nsclc_v2.0", nsclc_files)
    assert all(n == 0 for n in raw.run_checks(con).values())


def test_the_crc_spec_does_not_fit_nsclc_files(nsclc_files):
    with pytest.raises(duckdb.Error, match="crc_type"):
        raw.load("crc_v2.0", nsclc_files)


def test_the_standard_tables_carry_smoking_and_no_msi_mmr(nsclc_files, nsclc):
    res = genie_bpc.read(nsclc_files, None, nsclc)
    p = res.tables.patients
    assert "smoking_at_diagnosis" in p and "tumor_location" not in p
    assert set(p["smoking_at_diagnosis"].dropna()) <= set(SMOKING)
    assert set(res.tables.results["kind"]) == {"mutation", "fusion"}
    assert (res.tables.tests["modality"] == "ngs").all()


def test_fusions_become_results_of_kind_fusion_with_their_description(nsclc_files, nsclc):
    res = genie_bpc.read(nsclc_files, None, nsclc)
    fusions = res.tables.results[res.tables.results["kind"] == "fusion"]
    assert set(fusions["gene"]) == {"EML4", "ALK", "CD74", "ROS1"}
    assert set(fusions.loc[fusions["gene"] == "ALK", "variant"]) == {"EML4-ALK fusion"}
    assert fusions["variant_class"].isna().all() and fusions["value"].isna().all()
    assert res.tables.treatments["regimen"].notna().any()          # regimen_drugs is read


def test_a_two_gene_fusion_puts_its_patients_in_the_alk_group(nsclc_files, nsclc):
    res = genie_bpc.read(nsclc_files, None, nsclc)
    c = build_cohort(outcome_blind(res.tables), nsclc)
    assert (c.patients.loc[ALK_PATIENTS, "genotype"] == "ALK fusion").all()
    assert c.patients.loc[PARTNER_ONLY, "genotype"] != "ALK fusion"       # no ALK row: the download check flags it


def test_the_nsclc_settings_screen_nsclc_shaped_files(nsclc_files, nsclc):
    res = genie_bpc.read(nsclc_files, None, nsclc)
    rep = screen(outcome_blind(res.tables), nsclc, res.profile, res.source_steps)
    primary = rep.decisions["primary"]
    assert primary["EGFR common vs No listed driver"] == COMPARE
    assert primary["ALK fusion vs No listed driver"] == DESCRIBE           # a few made-up fusions: short
    assert rep.group_sizes["primary"]["ALK fusion"] == "<11"
    assert rep.find(topic="panel coverage")[0].severity == "deviation"     # genes listed only, as the settings say
    fusion = rep.find(topic="fusion detection")[0]                           # no panel list yet (review question 6)
    assert fusion.severity == "deviation" and "ALK" in fusion.finding


def nsclc_settings_file(tmp_path, changes: dict):
    """The NSCLC draft with some values changed, written to a new file (made-up values, for the tests only)."""
    data = cd.yaml.safe_load(NSCLC_SETTINGS.read_text(encoding="utf-8"))
    for path, value in changes.items():
        node = data
        for part in path.split("."):
            node = node[part]
        node["value"] = value
    p = tmp_path / "nsclc_changed.yaml"
    p.write_text(cd.yaml.safe_dump(data, sort_keys=False, allow_unicode=True), encoding="utf-8")
    return p


MADE_UP_FUSION_PANELS = ["DFCI-PANEL-1", "DFCI-PANEL-2", "MSK-PANEL-1", "MSK-PANEL-2", "VICC-PANEL-1"]   # not VICC-PANEL-2


def test_the_fusion_panel_list_decides_who_is_assignable(nsclc_files, tmp_path):
    s = load_settings(nsclc_settings_file(tmp_path, {"genotype.coverage.fusion_panels": {"ALK": MADE_UP_FUSION_PANELS}}))
    c = build_cohort(outcome_blind(genie_bpc.read(nsclc_files, None, s).tables), s)
    unlisted = c.patients["first_panel_id"] == "VICC-PANEL-2"
    assert unlisted.any() and c.patients.loc[unlisted, "genotype"].isna().all()
    assert c.patients.loc[~unlisted & c.patients["first_panel_id"].notna(), "genotype"].notna().all()
    assert (c.patients.loc[ALK_PATIENTS, "genotype"] == "ALK fusion").all()          # sequenced at DFCI


def test_the_download_check_reads_the_fusion_panel_list(nsclc_files, tmp_path):
    path = nsclc_settings_file(tmp_path, {"genotype.coverage.fusion_panels": {"ALK": MADE_UP_FUSION_PANELS + ["DFCI-PANEL-9"]}})
    checks = {c.name: c for c in cd.run_checks(nsclc_files, nsclc_files, tmp_path, path, rehearsal=True)}
    named = checks["Every panel the settings list as detecting fusions has a panel file"]
    assert named.status == cd.FAIL and named.detail.endswith("DFCI-PANEL-9")
    every = checks["Every panel reads EGFR, KRAS, ALK"]
    assert every.status == cd.REVIEW and every.detail.endswith("VICC-PANEL-2")


def test_the_checks_on_other_primaries_are_skipped_when_no_analysis_reads_them(nsclc_files, tmp_path):
    checks = {c.name: c for c in cd.run_checks(nsclc_files, nsclc_files, tmp_path, NSCLC_SETTINGS, rehearsal=True)}
    for name in ("Non-index cancers at this cancer's site codes are typed as this cancer",
                 "Other primaries of this cancer in the non-index file use the labels the settings match"):
        assert checks[name].status == cd.SKIPPED and "other_primary_same_cancer" in checks[name].detail


def test_the_download_check_finds_the_listed_drugs_as_the_screen_does(nsclc_files, tmp_path):
    path = nsclc_settings_file(tmp_path, {"genotype.reference_drugs": {"ALK": ["capecitabine", "Drug Never Given"]}})
    checks = {c.name: c for c in cd.run_checks(nsclc_files, nsclc_files, tmp_path, path, rehearsal=True)}
    drugs = checks["The drugs the settings list appear in the regimen file as the screen finds them"]
    assert drugs.status == cd.REVIEW and drugs.detail.endswith(": Drug Never Given")
    s = load_settings(path)
    res = genie_bpc.read(nsclc_files, None, s)
    item = screen(outcome_blind(res.tables), s, res.profile).find(topic="targeted drugs in the reference")[0]
    assert "of the No listed driver group's eligible patients received a drug the settings list for ALK" in item.finding


def test_the_rehearsal_download_check_passes_on_nsclc_shaped_files(nsclc_files, tmp_path):
    checks = cd.run_checks(nsclc_files, nsclc_files, tmp_path, NSCLC_SETTINGS, rehearsal=True)
    assert [c.name for c in checks if c.status == cd.FAIL] == []
    names = {c.name: c.status for c in checks}
    assert names["The genes the groups name appear in the mutation file"] == cd.PASS          # EGFR, KRAS; not ALK
    assert names["The genes the groups count as fusions appear in the fusion file"] == cd.PASS
    partner = next(c for c in checks if c.name == "Every fusion named for a gene the groups count has a row for that gene")
    assert partner.status == cd.REVIEW and partner.detail.endswith("ALK: <11 samples")
    assert names["A specimen label of this cancer appears in the pathology file"] == cd.PASS


def test_without_its_release_list_the_real_check_stops(nsclc_files, tmp_path):
    expected = cd.yaml.safe_load((raw.cohort_dir("nsclc_v2.0") / "expected_release.yml").read_text(encoding="utf-8"))
    expected["files"] = "not_made_yet.yml"
    checks = cd.check_files(nsclc_files, nsclc_files, "nsclc_v2.0", expected, rehearsal=False)
    stop = next(c for c in checks if c.name == "The release list holds each file's size and checksum")
    assert stop.status == cd.FAIL and "synapse_release_list.py nsclc_v2.0" in stop.detail


def test_the_nsclc_release_list_covers_every_file_the_adapter_reads():
    """The list made from Synapse's metadata (scripts/synapse_release_list.py, 2026-11-13)."""
    expected = cd.yaml.safe_load((raw.cohort_dir("nsclc_v2.0") / "expected_release.yml").read_text(encoding="utf-8"))
    files = cd.release_files("nsclc_v2.0", expected)
    assert files is not None and len(files) == 56
    used = {f["name"]: f for f in files if f["used"]}
    for name, folder, _ in raw.raw_files("nsclc_v2.0").values():
        assert used[name]["folder"] == folder
    assert {n for n in used if n.startswith("data_gene_panel_")} == {f"data_gene_panel_{p}.txt" for p in expected["panels"]}
    assert "data_fusions.txt" in used and len(used) == 19
    assert all(len(f["md5"]) == 32 and f["bytes"] > 0 for f in files)
    assert {f["folder"] for f in files} == {"clinical", "genomic"}


def test_the_nsclc_release_lists_its_folders_and_panels():
    expected = cd.yaml.safe_load((raw.cohort_dir("nsclc_v2.0") / "expected_release.yml").read_text(encoding="utf-8"))
    assert expected["folders"] == {"clinical": "syn30358089", "genomic": "syn30358098"}
    assert len(expected["panels"]) == 11 and "UHN-50-V2" in expected["panels"]
    assert expected["values"]["raw_patient"]["institution"]["allowed"] == ["DFCI", "MSK", "UHN", "VICC"]


# ---- The release-list script, without Synapse ----
spec = importlib.util.spec_from_file_location("synapse_release_list", ROOT / "scripts" / "synapse_release_list.py")
srl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(srl)


def test_used_marks_reproduce_the_crc_release_list():
    crc = cd.yaml.safe_load((raw.cohort_dir("crc_v2.0") / "expected_release.yml").read_text(encoding="utf-8"))
    assert all(srl.is_used(f["name"], "crc_v2.0") is f["used"] for f in crc["files"])


@pytest.mark.parametrize("name, used", [("data_fusions.txt", True), ("data_gene_panel_UHN-48-V1.txt", True),
                                        ("patient_level_dataset.csv", True), ("imaging_level_dataset.csv", False),
                                        ("case_lists/cases_all.txt", False), ("tmb.tsv", False)])
def test_used_marks_for_nsclc(name, used):
    assert srl.is_used(name, "nsclc_v2.0") is used


def test_the_release_list_file_is_what_the_check_reads(tmp_path):
    entry = srl.release_file_entry("genomic", "data_fusions.txt", "syn30358105", 2, 1234, "0" * 32, "nsclc_v2.0")
    path = srl.write_release_files([entry], "nsclc_v2.0", tmp_path / "release_files.yml")
    assert cd.yaml.safe_load(path.read_text(encoding="utf-8"))["files"] == [entry]
    assert entry["used"] is True


def test_size_and_md5_are_read_from_either_place():
    class Legacy:
        _file_handle = {"contentSize": 10, "contentMd5": "a" * 32}

    class Newer:
        contentSize, md5 = 20, "b" * 32

    assert srl.size_and_md5(Legacy()) == (10, "a" * 32)
    assert srl.size_and_md5(Newer()) == (20, "b" * 32)
