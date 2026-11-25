"""The adapter's NSCLC v2.0-public spec, on made-up files shaped like that release.

There is no made-up NSCLC release yet, so the made-up CRC files are reshaped the way the NSCLC files differ
(variable synopsis, 2026-11-13): no crc_type but ca_lung_cigarette; 25 specimen columns, no MSI/MMR; lung
labels; panels that read EGFR and ALK; the made-up BRAF V600E calls turned into EGFR L858R calls; and a made-up
fusion file in cBioPortal's format, written DFCI-style (a row for each fused gene, each with its own gene first:
"ALK-EML4" under ALK, "EML4-ALK" under EML4): EML4-ALK on three patients' samples, a CD74-ROS1 on one, a
KIF5B-RET on one, and one EML4-ALK listed under EML4 only. For question 2's fifth group, two made-up
patients gain a BRAF p.V600E or a MET p.X1010_splice call, and the panels list MET."""
import importlib.util
import re
import shutil

import duckdb
import pandas as pd
import pytest

from delayed_entry.adapters import genie_bpc
from delayed_entry.adapters.genie_bpc import check_download as cd
from delayed_entry.adapters.genie_bpc import raw
from delayed_entry.cohort import build_cohort
from delayed_entry.screen import COMPARE, DESCRIBE, NOTE, screen
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
RET_PATIENT = "GENIE-DFCI-S00019"      # made-up
OTHER_DRIVER_CALLS = {"GENIE-DFCI-S00005": [("BRAF", "Missense_Mutation", "p.V600E")],     # made-up
                      "GENIE-DFCI-S00014": [("MET", "Splice_Site", "p.X1010_splice")]}
FUSION_COLUMNS = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "Tumor_Sample_Barcode", "Fusion", "DNA_support",
                  "RNA_support", "Method", "Frame"]


def dfci_rows(reports: pd.DataFrame, pid, pair, genes) -> list[dict]:
    """Made-up DFCI-style rows of one fusion on every sample of a patient: a row for each of `genes` (the fused
    genes the panel reads), each with its own gene first, and the frame left blank, as DFCI writes them."""
    samples = reports.loc[reports["record_id"] == pid, "cpt_genie_sample_id"]
    return [dict(Hugo_Symbol=g, Entrez_Gene_Id="", Center="DFCI", Tumor_Sample_Barcode=s,
                 Fusion=f"{g}-{pair[1] if g == pair[0] else pair[0]}", DNA_support="yes", RNA_support="unknown",
                 Method="NGS", Frame="")
            for s in samples for g in genes]


def fusion_file(reports: pd.DataFrame) -> pd.DataFrame:
    """Made-up fusions on every sample of the chosen patients, DFCI-style."""
    out = [r for pid in ALK_PATIENTS for r in dfci_rows(reports, pid, ("EML4", "ALK"), ["EML4", "ALK"])]
    out += dfci_rows(reports, PARTNER_ONLY, ("EML4", "ALK"), ["EML4"])
    out += dfci_rows(reports, ROS1_PATIENT, ("CD74", "ROS1"), ["CD74", "ROS1"])
    out += dfci_rows(reports, RET_PATIENT, ("KIF5B", "RET"), ["KIF5B", "RET"])
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
        p.write_text(p.read_text(encoding="utf-8").rstrip("\n") + "\tEGFR\tALK\tMET\n", encoding="utf-8")
    maf = read(d / "data_mutations_extended.txt", "\t")
    v600e = (maf["Hugo_Symbol"] == "BRAF") & (maf["HGVSp_Short"].str.strip() == "p.V600E")
    maf.loc[v600e, ["Hugo_Symbol", "HGVSp_Short"]] = ["EGFR", "p.L858R"]
    reports = read(d / "cancer_panel_test_level_dataset.csv")
    extra = [dict(Hugo_Symbol=g, Variant_Classification=c, HGVSp_Short=v, Tumor_Sample_Barcode=s)
             for pid, calls in OTHER_DRIVER_CALLS.items() for g, c, v in calls
             for s in reports.loc[reports["record_id"] == pid, "cpt_genie_sample_id"]]
    maf = pd.concat([maf, pd.DataFrame(extra)], ignore_index=True)
    write(maf, d / "data_mutations_extended.txt", "\t")
    write(fusion_file(read(d / "cancer_panel_test_level_dataset.csv")), d / "data_fusions.txt", "\t")
    return d


@pytest.fixture(scope="module")
def nsclc(tmp_path_factory):
    """The NSCLC settings as decided, with the made-up files' panels standing in for the release's."""
    return load_settings(nsclc_settings_file(tmp_path_factory.mktemp("nsclc"), {}))


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
    assert set(fusions["gene"]) == {"EML4", "ALK", "CD74", "ROS1", "KIF5B", "RET"}
    assert set(fusions.loc[fusions["gene"] == "ALK", "variant"]) == {"ALK-EML4"}          # DFCI-style: own gene first
    assert set(fusions.loc[fusions["gene"] == "EML4", "variant"]) == {"EML4-ALK"}
    assert fusions["variant_class"].isna().all() and fusions["value"].isna().all()
    assert res.tables.treatments["regimen"].notna().any()          # regimen_drugs is read


def test_a_two_gene_fusion_puts_its_patients_in_the_alk_group(nsclc_files, nsclc):
    res = genie_bpc.read(nsclc_files, None, nsclc)
    c = build_cohort(outcome_blind(res.tables), nsclc)
    assert (c.patients.loc[ALK_PATIENTS, "genotype"] == "ALK fusion").all()
    assert c.patients.loc[PARTNER_ONLY, "genotype"] != "ALK fusion"       # no ALK row: the download check flags it


def test_a_call_whose_gene_the_panel_doesnt_list_still_places_its_patient(nsclc_files, nsclc, tmp_path):
    """What the download check says of such calls: the cohort uses them. The coverage rule reads only EGFR, KRAS
    and ALK, so a made-up BRAF V600E on panels that no longer list BRAF still places its patient."""
    from test_check_download import drop_gene_from_panels
    files = tmp_path / "files"
    shutil.copytree(nsclc_files, files)
    drop_gene_from_panels(files, "BRAF")
    c = build_cohort(outcome_blind(genie_bpc.read(files, None, nsclc).tables), nsclc)
    assert c.patients.loc["GENIE-DFCI-S00005", "genotype"] == "Other recognized driver"
    assert OTHER_DRIVER_CALLS["GENIE-DFCI-S00005"][0][0] == "BRAF"


def test_the_nsclc_settings_screen_nsclc_shaped_files(nsclc_files, nsclc):
    res = genie_bpc.read(nsclc_files, None, nsclc)
    rep = screen(outcome_blind(res.tables), nsclc, res.profile, res.source_steps)
    primary = rep.decisions["primary"]
    assert primary["EGFR common vs No recognized driver"] == COMPARE
    assert primary["ALK fusion vs No recognized driver"] == DESCRIBE           # a few made-up fusions: short
    assert rep.group_sizes["primary"]["ALK fusion"] == "<11"
    assert rep.find(topic="panel coverage")[0].severity == "note"     # genes listed only, checked at gene level
    assert primary["Other recognized driver vs No recognized driver"] == DESCRIBE      # question 13, option B
    described = rep.find(analysis="primary", topic="described only")
    assert len(described) == 1 and "Other recognized driver" in described[0].finding
    fusion = rep.find(topic="fusion detection")[0]                           # from the list of panels (question 6)
    assert fusion.severity == "note" and "ALK" in fusion.finding and "list of panels" in fusion.finding
    assert not [w for w in rep.warnings if w.topic == "fusion detection"]    # every listed panel is in the files


def test_with_no_fusion_panel_list_a_panel_listing_the_gene_counts_as_a_deviation(nsclc_files, tmp_path):
    """Under "deviation", the draft's value; as decided (accept_gene_level, 2026-11-24) the same finding is a note."""
    s = load_settings(nsclc_settings_file(tmp_path, {"genotype.coverage.fusion_panels": {},
                                                     "genotype.coverage.when_only_genes_listed": "deviation"}))
    res = genie_bpc.read(nsclc_files, None, s)
    fusion = screen(outcome_blind(res.tables), s, res.profile, res.source_steps).find(topic="fusion detection")[0]
    assert fusion.severity == "deviation" and "name no panels that detect them" in fusion.finding


MADE_UP_FUSION_PANELS = ["DFCI-PANEL-1", "DFCI-PANEL-2", "MSK-PANEL-1", "MSK-PANEL-2", "VICC-PANEL-1"]   # not VICC-PANEL-2


def nsclc_settings_file(tmp_path, changes: dict):
    """The NSCLC settings with some values changed, written to a new file (made-up values, for the tests only).
    The made-up files' panels stand in for the release's: every made-up panel but VICC-PANEL-2 detects ALK fusions,
    as the release's capture panels do and its hotspot panels don't. SA12's panels follow the ALK list, since the
    settings require one list for both."""
    data = cd.yaml.safe_load(NSCLC_SETTINGS.read_text(encoding="utf-8"))
    data["genotype"]["coverage"]["fusion_panels"]["value"] = {"ALK": list(MADE_UP_FUSION_PANELS)}
    for path, value in changes.items():
        node = data
        for part in path.split("."):
            node = node[part]
        node["value"] = value
    sa12 = next(a for a in data["analyses"]["value"] if a["id"] == "SA12")
    sa12["keep"]["fusion_detecting_panel"]["panels"] = list(
        (data["genotype"]["coverage"]["fusion_panels"]["value"] or {}).get("ALK", []))
    p = tmp_path / "nsclc_changed.yaml"
    p.write_text(cd.yaml.safe_dump(data, sort_keys=False, allow_unicode=True), encoding="utf-8")
    return p


def test_the_decided_alk_panels_are_the_release_capture_panels():
    """Question 6, option A: the eight capture panels detect ALK fusions; the three hotspot panels, UHN's two among
    them, don't. So no UHN patient is assignable: SA11 and SA12 would repeat the primary, and center loses UHN."""
    s = load_settings(NSCLC_SETTINGS)
    expected = cd.yaml.safe_load((raw.cohort_dir("nsclc_v2.0") / "expected_release.yml").read_text(encoding="utf-8"))
    alk = s.read("genotype.coverage.fusion_panels")["ALK"]
    assert set(alk) < set(expected["panels"]) and len(alk) == len(set(alk))
    assert set(expected["panels"]) - set(alk) == {"UHN-48-V1", "UHN-50-V2", "VICC-01-SOLIDTUMOR"}
    assert s.read("genotype.coverage.applies_to") == "all_groups"
    runs = {a["id"]: a["run"] for a in s.read("analyses")}
    assert runs["SA11"] is False and runs["SA12"] is False
    center = next(v for v in s.read("adjustment.variables") if v["field"] == "center")
    assert center["keep"] == ["DFCI", "MSK", "VICC"] and center["reference"] == "MSK"


def test_the_decided_population_is_nonsquamous():
    """Question 1, option B: the histology rule keeps the nonsquamous labels and passes a blank, so the release's
    other labels are the ones left out; SA10 would repeat the primary, SA7 still narrows it."""
    s = load_settings(NSCLC_SETTINGS)
    nonsquamous = ["Adenocarcinoma", "Carcinoma", "Other histologies/mixed tumor"]
    rule = next(r for r in s.read("cohort.rules") if r["id"] == "histology")
    (cond,) = rule["all_of"]
    assert cond == {"field": "histology", "op": "in", "values": nonsquamous, "normalize": "trim", "missing": "pass"}
    expected = cd.yaml.safe_load((raw.cohort_dir("nsclc_v2.0") / "expected_release.yml").read_text(encoding="utf-8"))
    allowed = expected["values"]["raw_cancer"]["ca_hist_adeno_squamous"]["allowed"]
    assert set(allowed) - set(nonsquamous) == {"Squamous cell", "Small cell carcinoma", "Sarcoma"}
    runs = {a["id"]: a["run"] for a in s.read("analyses")}
    assert runs["SA10"] is False and runs["SA7"] is True
    histology = next(v for v in s.read("adjustment.variables") if v["field"] == "histology")
    assert histology["reference"] == "Adenocarcinoma"            # question 11: the three levels merged into two
    assert [x for xs in histology["merge"].values() for x in xs] == nonsquamous


def test_squamous_patients_leave_at_the_histology_rule(nsclc_files, nsclc):
    c = build_cohort(outcome_blind(genie_bpc.read(nsclc_files, None, nsclc).tables), nsclc)
    for pid in ("GENIE-DFCI-S01323", "GENIE-MSK-S00153"):           # made-up stage IV, Squamous cell
        assert c.patients.loc[pid, "excluded_at"] == "histology"


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
    path = nsclc_settings_file(tmp_path, {})
    checks = {c.name: c for c in cd.run_checks(nsclc_files, nsclc_files, tmp_path, path, rehearsal=True)}
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
    assert "of the No recognized driver group's eligible patients received a drug the settings list for ALK" in item.finding


def test_the_rehearsal_download_check_passes_on_nsclc_shaped_files(nsclc_files, tmp_path):
    checks = cd.run_checks(nsclc_files, nsclc_files, tmp_path, nsclc_settings_file(tmp_path, {}), rehearsal=True)
    assert [c.name for c in checks if c.status == cd.FAIL] == []
    names = {c.name: c.status for c in checks}
    assert names["The genes the groups name appear in the mutation file"] == cd.PASS          # EGFR, KRAS; not ALK
    assert names["The genes the groups count as fusions appear in the fusion file"] == cd.PASS
    partner = next(c for c in checks if c.name == "Every fusion named for a gene the groups count has a row for that gene")
    assert partner.status == cd.REVIEW and partner.detail.endswith("ALK: <11 samples")
    assert names["A specimen label of this cancer appears in the pathology file"] == cd.PASS
    assert names["The labels a cohort rule matches on histology appear as the rule compares them"] == cd.PASS


# ---- 6a (o): the label check is gene-aware ----
def test_a_met_splice_label_moved_to_another_gene_stops_the_download_check(nsclc_files, tmp_path):
    """MET's p.X1010_splice guard holds only if the label is in a MET call: the made-up MET call moved to EGFR (a
    gene the groups name) FAILs, though the label is still in the file; a MET call of another kind keeps MET in it."""
    files = tmp_path / "files"
    shutil.copytree(nsclc_files, files)
    maf = read(files / "data_mutations_extended.txt", "\t")
    met = maf["HGVSp_Short"] == "p.X1010_splice"
    extra = maf[met].assign(Variant_Classification="Missense_Mutation", HGVSp_Short="p.R988C")
    maf.loc[met, "Hugo_Symbol"] = "EGFR"
    write(pd.concat([maf, extra], ignore_index=True), files / "data_mutations_extended.txt", "\t")
    checks = {c.name: c for c in cd.run_checks(files, files, tmp_path, nsclc_settings_file(tmp_path, {}), rehearsal=True)}
    label = checks["The protein changes the groups name are written in that form, in the genes their rules name"]
    assert (label.status, label.detail) == (cd.FAIL, "p.X1010_splice in MET")
    assert checks["The genes the groups name appear in the mutation file"].status == cd.PASS
    unchanged = {c.name: c for c in cd.run_checks(nsclc_files, nsclc_files, tmp_path, nsclc_settings_file(tmp_path, {}),
                                                  rehearsal=True)}
    assert unchanged[label.name].status == cd.PASS


# ---- 6a (q): the transcript check ----
TRANSCRIPTS = {"MET": "ENST00000397752", "EGFR": "ENST00000275493", "ERBB2": "ENST00000269571",
               "KRAS": "ENST00000256078", "BRAF": "ENST00000288602"}
RULE_GENES = ["BRAF", "EGFR", "ERBB2", "KRAS", "MET"]


def test_the_nsclc_spec_expects_the_transcripts_genie_annotates_on_and_the_crc_spec_none():
    nsclc = cd.yaml.safe_load((raw.cohort_dir("nsclc_v2.0") / "expected_release.yml").read_text(encoding="utf-8"))
    crc = cd.yaml.safe_load((raw.cohort_dir("crc_v2.0") / "expected_release.yml").read_text(encoding="utf-8"))
    assert nsclc["transcripts"] == TRANSCRIPTS and "transcripts" not in crc
    assert cd.searched_labels(load_settings(NSCLC_SETTINGS))["genes"] == RULE_GENES      # every rule gene has one


def with_transcripts(nsclc_files, dst, edit=None):
    """A copy of the made-up mutation file with GENIE's Transcript_ID and HGVSc columns, every call on its gene's
    expected transcript (other genes on a made-up ID); `edit(maf)` then changes some rows."""
    maf = read(nsclc_files / "data_mutations_extended.txt", "\t")
    maf["Transcript_ID"] = maf["Hugo_Symbol"].map(TRANSCRIPTS).fillna("ENST00000999999")
    maf["HGVSc"] = maf["Transcript_ID"] + ".1:c.1A>G"
    if edit:
        edit(maf)
    write(maf, dst, "\t")
    return dst


def transcript_check(path, genes=RULE_GENES, transcripts=TRANSCRIPTS, rehearsal=False):
    con = duckdb.connect()
    cd.read_as_text(con, "raw_mutation", path)
    (check,) = cd.check_transcripts(con, {"raw_mutation"}, {"genes": genes}, transcripts, rehearsal, cd.Shower(11))
    return check


def test_calls_on_the_expected_transcripts_pass_with_or_without_a_version(nsclc_files, tmp_path):
    check = transcript_check(with_transcripts(nsclc_files, tmp_path / "maf.txt"))
    assert check.status == cd.PASS and "rows, all on ENST00000275493" in check.detail and "MET: <11 rows" in check.detail

    def versioned(maf):
        maf["Transcript_ID"] = maf["Transcript_ID"] + ".2"
    assert transcript_check(with_transcripts(nsclc_files, tmp_path / "v.txt", versioned)).status == cd.PASS


def test_a_few_calls_on_another_transcript_are_a_review_and_most_a_failure(nsclc_files, tmp_path):
    def two_egfr(maf):
        rows = maf.index[maf["Hugo_Symbol"] == "EGFR"][:2]
        maf.loc[rows, "Transcript_ID"] = "ENST00000455089"
    check = transcript_check(with_transcripts(nsclc_files, tmp_path / "a.txt", two_egfr))
    assert check.status == cd.REVIEW and "<11 on another transcript (ENST00000455089)" in check.detail
    assert check.detail.endswith("a deviation decision, made before any outcome is seen")

    def met_long(maf):           # MET on the longer transcript, whose codons run 18 higher
        maf.loc[maf["Hugo_Symbol"] == "MET", ["Transcript_ID", "HGVSc"]] = ["ENST00000318493", "ENST00000318493.6:c.1A>G"]
    check = transcript_check(with_transcripts(nsclc_files, tmp_path / "b.txt", met_long))
    assert check.status == cd.FAIL and "MET: <11 rows, <11 on another transcript (ENST00000318493)" in check.detail


def test_hgvsc_on_another_transcript_or_no_transcript_is_a_review(nsclc_files, tmp_path):
    def hgvsc_other(maf):
        maf.loc[maf.index[maf["Hugo_Symbol"] == "KRAS"][:1], "HGVSc"] = "ENST00000311936.8:c.35G>T"
    check = transcript_check(with_transcripts(nsclc_files, tmp_path / "a.txt", hgvsc_other))
    assert check.status == cd.REVIEW and "(ENST00000311936)" in check.detail

    def blank(maf):
        maf.loc[maf.index[maf["Hugo_Symbol"] == "KRAS"][:1], ["Transcript_ID", "HGVSc"]] = ["", "c.35G>T"]
    check = transcript_check(with_transcripts(nsclc_files, tmp_path / "b.txt", blank))
    assert check.status == cd.REVIEW and "<11 with no transcript" in check.detail


def test_a_gene_with_no_expected_transcript_or_a_file_with_no_transcripts(nsclc_files, tmp_path):
    path = with_transcripts(nsclc_files, tmp_path / "maf.txt")
    check = transcript_check(path, ["EGFR", "NRAS"], {"EGFR": TRANSCRIPTS["EGFR"]})
    assert check.status == cd.REVIEW and "NRAS: no transcript expected in the cohort's spec" in check.detail
    plain = nsclc_files / "data_mutations_extended.txt"           # the made-up file has no Transcript_ID column
    assert transcript_check(plain, rehearsal=True).status == cd.SKIPPED
    real = transcript_check(plain)
    assert real.status == cd.REVIEW and "count the rows by transcript by hand" in real.detail
    con = duckdb.connect()
    assert cd.check_transcripts(con, set(), {"genes": ["EGFR"]}, None, False, cd.Shower(11))[0].status == cd.SKIPPED


# ---- 6a (p): each group rule's matches, report only ----
def test_the_download_check_counts_each_group_rule_and_never_fails_on_it(nsclc_files, tmp_path):
    path = nsclc_settings_file(tmp_path, {})
    checks = cd.run_checks(nsclc_files, nsclc_files, tmp_path, path, rehearsal=True)
    rules = {c.name: c for c in checks if c.section.startswith("I. ")}
    groups = load_settings(path).read("genotype.groups")
    assert len(rules) == sum(len(g.get("any_of", [])) for g in groups) == 15
    assert {c.status for c in rules.values()} <= {cd.PASS, cd.REVIEW}
    l858r = rules["EGFR common, rule 2: mutation in EGFR; p.L858R"]
    assert l858r.status == cd.PASS and re.fullmatch(r"[\d,]+ calls in [\d,]+ patients", l858r.detail)
    alk = rules["ALK fusion, rule 1: fusion in ALK; less its listed exceptions"]
    assert alk.status == cd.PASS and alk.detail == "<11 calls in <11 patients"
    assert not [c for c in rules.values() if "intragenic" in c.name + c.detail]       # except-list text isn't shown
    met = rules["Other recognized driver, rule 3: mutation in MET; Splice_Site; codons 963-1010"]
    assert met.detail == ("<11 calls in <11 patients; first codons matched: 1010; its genes' calls of these classes "
                          "at other codons: 0")
    kras = rules["KRAS mutant, rule 2: mutation in KRAS; Missense_Mutation, In_Frame_Del, In_Frame_Ins; codons 12, 13, 61"]
    assert kras.status == cd.PASS and "; codons matched: 12" in kras.detail and "at other codons: " in kras.detail
    erbb2 = rules["Other recognized driver, rule 6: mutation in ERBB2; In_Frame_Ins; codons 770-831"]
    assert erbb2.status == cd.REVIEW
    assert erbb2.detail.startswith("0 calls in 0 patients") and "no call matches this rule" in erbb2.detail
    assert cd.verdict(checks).startswith("PASSED")


def test_the_rule_counts_wait_for_the_trial_build_and_mask_small_codon_counts(nsclc):
    (check,) = cd.check_rule_counts(None, nsclc, cd.Shower(11))
    assert check.status == cd.SKIPPED
    assert cd.counts_unless_small({12: 30, 13: 15}, cd.Shower(11)) == "12 30, 13 15"
    assert cd.counts_unless_small({12: 30, 13: 2}, cd.Shower(11)) == "12, 13"     # 13's count would follow


# ---- 6a (r): fusion genes no panel is checked for ----
UNCHECKED = "but the coverage rule doesn't require"


def test_the_readiness_report_notes_fusion_genes_the_coverage_rule_does_not_require(nsclc_files, tmp_path):
    """From the settings alone: the groups count ROS1 and RET fusions, and genotype.coverage.required names EGFR, KRAS
    and ALK, so one note says no panel is checked for detecting them. A made-up change that requires RET leaves ROS1
    alone in the note; requiring both removes it."""
    def note(changes):
        s = load_settings(nsclc_settings_file(tmp_path, changes))
        res = genie_bpc.read(nsclc_files, None, s)
        return [i for i in screen(outcome_blind(res.tables), s, res.profile).find(topic="fusion detection")
                if UNCHECKED in i.finding]
    (item,) = note({})
    assert item.severity == NOTE and item.analysis == "primary" and item.settings == ("genotype.coverage.required",)
    assert item.finding == ("the groups count RET, ROS1 fusions, but the coverage rule doesn't require these genes, so no "
                            "panel is checked for detecting them")
    assert item.effect == ("a tumor whose first panel can't see such a fusion is placed by its other calls, in No "
                           "recognized driver if it has none; stated in the plan")
    required = load_settings(NSCLC_SETTINGS).read("genotype.coverage.required")
    (item,) = note({"genotype.coverage.required": {**required, "RET": "whole_gene"}})
    assert item.finding.startswith("the groups count ROS1 fusions, but the coverage rule doesn't require this gene")
    assert note({"genotype.coverage.required": {**required, "RET": "whole_gene", "ROS1": "whole_gene"}}) == []


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


# ---- Calls no rule can read, and fusions in genes no fusion rule names (6a (g) and (i)) ----
def with_unreadable_calls(nsclc_files, nsclc, folder):
    """A copy of the made-up files where four made-up reference patients gain one call each that no rule reads:
    an EGFR deletion with its class written In_Frame_DEL, an EGFR call with a blank class, an EGFR deletion with
    a blank protein label, and an EGFR fusion row. Returns the folder and the four patients, in that order."""
    c = build_cohort(outcome_blind(genie_bpc.read(nsclc_files, None, nsclc).tables), nsclc)
    chosen = sorted(c.eligible.index[c.eligible["genotype"] == "No recognized driver"])[:4]
    shutil.copytree(nsclc_files, folder)
    reports = read(folder / "cancer_panel_test_level_dataset.csv")
    samples = lambda pid: reports.loc[reports["record_id"] == pid, "cpt_genie_sample_id"]
    calls = [("EGFR", "In_Frame_DEL", "p.E746_A750del"), ("EGFR", "", "p.L747_P753delinsS"), ("EGFR", "In_Frame_Del", "")]
    rows = [dict(Hugo_Symbol=g, Variant_Classification=cl, HGVSp_Short=v, Tumor_Sample_Barcode=s)
            for pid, (g, cl, v) in zip(chosen, calls) for s in samples(pid)]
    maf = read(folder / "data_mutations_extended.txt", "\t")
    write(pd.concat([maf, pd.DataFrame(rows)], ignore_index=True), folder / "data_mutations_extended.txt", "\t")
    fus = [dict(Hugo_Symbol="EGFR", Entrez_Gene_Id="", Center="DFCI", Tumor_Sample_Barcode=s, Fusion="EGFR-SEPT14 fusion",
                DNA_support="yes", RNA_support="unknown", Method="NGS", Frame="in frame") for s in samples(chosen[3])]
    old = read(folder / "data_fusions.txt", "\t")
    write(pd.concat([old, pd.DataFrame(fus, columns=FUSION_COLUMNS)], ignore_index=True), folder / "data_fusions.txt", "\t")
    return folder, chosen


def test_driver_calls_no_rule_can_read_are_counted_and_change_no_group(nsclc_files, nsclc, tmp_path):
    files, chosen = with_unreadable_calls(nsclc_files, nsclc, tmp_path / "files")
    res = genie_bpc.read(files, None, nsclc)
    c = build_cohort(outcome_blind(res.tables), nsclc)
    assert (c.patients.loc[chosen, "genotype"] == "No recognized driver").all()
    items = screen(outcome_blind(res.tables), nsclc, res.profile, res.source_steps).find(
        topic="driver calls the rules can't read")
    findings = [i.finding for i in items]
    assert len(findings) == 3 and all(f.startswith("<11 eligible patients' first test has a call") for f in findings)
    assert any("written otherwise than the rules write it (In_Frame_DEL)" in f for f in findings)
    assert any("a blank variant class" in f for f in findings)
    assert any("a protein label that is blank or has no codon number" in f for f in findings)
    assert all("matches no rule" in i.effect for i in items)


def test_fusion_rows_in_a_gene_no_fusion_rule_names_are_counted_by_text(nsclc_files, nsclc, tmp_path):
    files, chosen = with_unreadable_calls(nsclc_files, nsclc, tmp_path / "files")
    res = genie_bpc.read(files, None, nsclc)
    assert build_cohort(outcome_blind(res.tables), nsclc).patients.loc[chosen[3], "genotype"] == "No recognized driver"
    items = screen(outcome_blind(res.tables), nsclc, res.profile, res.source_steps).find(topic="fusions no rule reads")
    assert [i.finding for i in items] == [
        "<11 eligible patients' first test has a fusion row in EGFR, which no group's rule reads "
        "(fusion texts: EGFR-SEPT14 fusion)"]


def test_the_made_up_files_as_they_are_have_no_unreadable_driver_calls(nsclc_files, nsclc):
    res = genie_bpc.read(nsclc_files, None, nsclc)
    rep = screen(outcome_blind(res.tables), nsclc, res.profile, res.source_steps)
    assert not rep.find(topic="driver calls the rules can't read") and not rep.find(topic="fusions no rule reads")


# ---- The overlapping-groups note, split by gene (6a (h)) ----
def test_the_screen_splits_overlaps_into_different_genes_and_one_gene(nsclc_files, nsclc, tmp_path):
    """Two made-up reference patients gain two drivers each: EGFR L858R with T790M (one gene, read by EGFR common
    and Other recognized driver) and EGFR L858R with KRAS G12C (two genes). Each part is counted on its own."""
    from delayed_entry.cohort import overlap_counts
    before = overlap_counts(build_cohort(outcome_blind(genie_bpc.read(nsclc_files, None, nsclc).tables), nsclc).eligible,
                            nsclc)
    c0 = build_cohort(outcome_blind(genie_bpc.read(nsclc_files, None, nsclc).tables), nsclc)
    one, two = sorted(c0.eligible.index[c0.eligible["genotype"] == "No recognized driver"])[-2:]
    files = tmp_path / "files"
    shutil.copytree(nsclc_files, files)
    reports = read(files / "cancer_panel_test_level_dataset.csv")
    calls = {one: [("EGFR", "p.L858R"), ("EGFR", "p.T790M")], two: [("EGFR", "p.L858R"), ("KRAS", "p.G12C")]}
    rows = [dict(Hugo_Symbol=g, Variant_Classification="Missense_Mutation", HGVSp_Short=v, Tumor_Sample_Barcode=smp)
            for pid, cs in calls.items() for g, v in cs for smp in reports.loc[reports["record_id"] == pid, "cpt_genie_sample_id"]]
    maf = read(files / "data_mutations_extended.txt", "	")
    write(pd.concat([maf, pd.DataFrame(rows)], ignore_index=True), files / "data_mutations_extended.txt", "	")
    res = genie_bpc.read(files, None, nsclc)
    c = build_cohort(outcome_blind(res.tables), nsclc)
    assert c.patients.loc[[one, two], "match_genes"].tolist() == [1, 2]
    after = overlap_counts(c.eligible, nsclc)
    assert after["same_gene"] == before["same_gene"] + 1 and after["different_genes"] == before["different_genes"] + 1
    notes = [i.finding for i in screen(outcome_blind(res.tables), nsclc, res.profile, res.source_steps)
             .find(topic="overlapping groups")]
    shown = lambda k: f"<11" if 0 < k < 11 else str(k)
    assert f"{shown(after['different_genes'])} eligible patients match more than one group's rules through calls in "            "different genes" in notes
    assert f"{shown(after['same_gene'])} eligible patients match more than one group's rules through calls in one "            "gene only" in notes
    total_line = [n for n in notes if n.endswith("in all")]
    small = any(0 < after[p] < 11 for p in ("same_gene", "different_genes"))
    assert total_line == ([] if small else [f"{after['total']} eligible patients match more than one group's rules in all"])



def test_an_alk_row_the_rule_excepts_is_counted_among_the_fusions_no_rule_reads(nsclc_files, nsclc, tmp_path):
    """The fusion rows (option C): an intragenic ALK row is excepted by its text, so its made-up patient stays in the
    reference, and the screen lists it with the rows no rule reads."""
    c = build_cohort(outcome_blind(genie_bpc.read(nsclc_files, None, nsclc).tables), nsclc)
    pid = sorted(c.eligible.index[c.eligible["genotype"] == "No recognized driver"])[0]
    files = tmp_path / "files"
    shutil.copytree(nsclc_files, files)
    reports = read(files / "cancer_panel_test_level_dataset.csv")
    rows = [dict(Hugo_Symbol="ALK", Entrez_Gene_Id="", Center="DFCI", Tumor_Sample_Barcode=s, Fusion="ALK-INTRAGENIC",
                 DNA_support="yes", RNA_support="unknown", Method="NGS", Frame="")
            for s in reports.loc[reports["record_id"] == pid, "cpt_genie_sample_id"]]
    old = read(files / "data_fusions.txt", "\t")
    write(pd.concat([old, pd.DataFrame(rows, columns=FUSION_COLUMNS)], ignore_index=True), files / "data_fusions.txt", "\t")
    res = genie_bpc.read(files, None, nsclc)
    assert build_cohort(outcome_blind(res.tables), nsclc).patients.loc[pid, "genotype"] == "No recognized driver"
    items = screen(outcome_blind(res.tables), nsclc, res.profile, res.source_steps).find(topic="fusions no rule reads")
    assert [i.finding for i in items] == [
        "<11 eligible patients' first test has a fusion row in ALK, which no group's rule reads (fusion texts: ALK-INTRAGENIC)"]
