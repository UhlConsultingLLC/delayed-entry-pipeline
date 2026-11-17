"""MAF and annotated VCF files through the flat-file adapter, on made-up text written at test time (none is kept in
the repository): VEP's CSQ and SnpEff's ANN annotations, the two reading rules (which annotation, which records),
alleles of a multi-allelic record, one VCF per sample read through a pattern, and the hgvsp_short transform."""
import gzip

import pandas as pd
import pytest
import yaml

from delayed_entry.adapters import flatfile
from delayed_entry.adapters.flatfile import MappingError, read_source
from delayed_entry.adapters.hgvs import short_protein_change
from delayed_entry.tables import table_fields
from helpers import crc_settings

S = crc_settings()
CSQ = ("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: "
       "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|HGVSc|HGVSp|CANONICAL\">")
ANN = ("##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | "
       "Annotation_Impact | Gene_Name | Gene_ID | Feature_ID | HGVS.c | HGVS.p' \">")
HEAD = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"


def csq(allele, consequence, gene, protein, canonical="YES", transcript="ENST0001"):
    return "|".join([allele, consequence, "MODERATE", gene, "ENSG0001", transcript, "c.1A>T", protein, canonical])


# Made-up records. KRAS lists a non-canonical transcript first, with its own numbering; TP53 fails the filter;
# GENE4 is multi-allelic (S-1 carries A, S-2 carries C); GENE5 is a deletion, which VEP names "-".
RECORDS = [
    ("12", "100", "C", "T", "PASS", [csq("T", "missense_variant", "KRAS", "ENSP2.1:p.Gly39Asp", "", "ENST0002"),
                                     csq("T", "missense_variant", "KRAS", "ENSP1.1:p.Gly12Asp")], ["0/1", "0/0"]),
    ("7", "200", "A", "T", "PASS", [csq("T", "missense_variant", "BRAF", "ENSP3.1:p.Val600Glu")], ["0/0", "0/1"]),
    ("17", "300", "C", "T", "LowQual", [csq("T", "missense_variant", "TP53", "ENSP4.1:p.Arg175His")], ["0/1", "0/1"]),
    ("4", "400", "G", "A,C", "PASS", [csq("A", "missense_variant", "GENE4", "ENSP5.1:p.Ala10Thr"),
                                      csq("C", "missense_variant", "GENE4", "ENSP5.1:p.Ala10Pro")], ["0/1", "0/2"]),
    ("13", "500", "AT", "A", "PASS", [csq("-", "frameshift_variant", "GENE5", "ENSP6.1:p.Ile5LysfsTer3")], ["1/1", "./."]),
]


def vcf_text(records=RECORDS, samples=("S-1", "S-2"), header_lines=(CSQ,), field="CSQ") -> str:
    lines = ["##fileformat=VCFv4.2", *header_lines, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    lines.append(HEAD + ("\tFORMAT\t" + "\t".join(samples) if samples else ""))
    for chrom, pos, ref, alt, flt, anns, gts in records:
        row = [chrom, pos, ".", ref, alt, "50", flt, f"{field}=" + ",".join(anns)]
        lines.append("\t".join(row + (["GT", *gts[:len(samples)]] if samples else [])))
    return "\n".join(lines) + "\n"


def write(path, text):
    path.write_text(text, encoding="utf-8")
    return path


def rows(df):
    return sorted(zip(df["SAMPLE"], df["SYMBOL"], df["HGVSp"]))


# ---- The VCF reader ----
def test_canonical_annotations_of_passing_records_for_each_carrier(tmp_path):
    df = read_source(write(tmp_path / "a.vcf", vcf_text()), annotation={"pick": "canonical", "filters": "pass"})
    assert rows(df) == [("S-1", "GENE4", "ENSP5.1:p.Ala10Thr"), ("S-1", "GENE5", "ENSP6.1:p.Ile5LysfsTer3"),
                        ("S-1", "KRAS", "ENSP1.1:p.Gly12Asp"),
                        ("S-2", "BRAF", "ENSP3.1:p.Val600Glu"), ("S-2", "GENE4", "ENSP5.1:p.Ala10Pro")]
    gene4 = df[df["SYMBOL"] == "GENE4"].set_index("SAMPLE")
    assert gene4.loc["S-1", "ALT"] == "A" and gene4.loc["S-2", "ALT"] == "C" and gene4.loc["S-2", "GT"] == "0/2"


@pytest.mark.parametrize("pick, kras", [("first", ["ENSP2.1:p.Gly39Asp"]),
                                        ("all", ["ENSP1.1:p.Gly12Asp", "ENSP2.1:p.Gly39Asp"])])
def test_the_pick_rule_decides_which_annotations_count(tmp_path, pick, kras):
    df = read_source(write(tmp_path / "a.vcf", vcf_text()), annotation={"pick": pick, "filters": "pass"})
    assert sorted(df.loc[df["SYMBOL"] == "KRAS", "HGVSp"]) == kras


def test_the_filters_rule_keeps_or_drops_records_that_failed_a_filter(tmp_path):
    path = write(tmp_path / "a.vcf", vcf_text())
    kept = read_source(path, annotation={"pick": "canonical", "filters": "all"})
    assert sorted(kept.loc[kept["SYMBOL"] == "TP53", "SAMPLE"]) == ["S-1", "S-2"]
    assert "TP53" not in set(read_source(path, annotation={"pick": "canonical", "filters": "pass"})["SYMBOL"])


def test_snpeff_annotations_and_a_gzipped_file(tmp_path):
    rec = [("12", "100", "C", "T", "PASS", ["T|missense_variant|MODERATE|KRAS|ENSG0001|ENST0001|c.35G>A|p.Gly12Asp"],
            ["0/1"])]
    text = vcf_text(rec, samples=("S-9",), header_lines=(ANN,), field="ANN")
    with gzip.open(tmp_path / "s.vcf.gz", "wt", encoding="utf-8") as fh:
        fh.write(text)
    df = read_source(tmp_path / "s.vcf.gz", annotation={"pick": "first", "filters": "pass"})
    assert df[["SAMPLE", "Gene_Name", "HGVS.p", "Annotation"]].values.tolist() == [["S-9", "KRAS", "p.Gly12Asp",
                                                                                     "missense_variant"]]
    with pytest.raises(MappingError, match="pick canonical needs a CANONICAL subfield"):
        read_source(tmp_path / "s.vcf.gz", annotation={"pick": "canonical", "filters": "pass"})


@pytest.mark.parametrize("text, annotation, problem", [
    (vcf_text(header_lines=()), {"pick": "all", "filters": "all"}, "annotate it first"),
    (vcf_text(header_lines=(CSQ, ANN)), {"pick": "all", "filters": "all"}, "both CSQ and ANN"),
    (vcf_text(), None, "a VCF needs its annotation's reading rules"),
    ("chrom\tpos\n1\t2\n", {"pick": "all", "filters": "all"}, "a record comes before the #CHROM header line"),
    ("#chrom\tpos\n1\t2\n", {"pick": "all", "filters": "all"}, "no VCF header line"),
])
def test_a_vcf_that_cant_be_read_says_why(tmp_path, text, annotation, problem):
    with pytest.raises(MappingError, match=problem):
        read_source(write(tmp_path / "a.vcf", text), annotation=annotation)


def test_with_both_annotations_the_field_says_which(tmp_path):
    df = read_source(write(tmp_path / "a.vcf", vcf_text(header_lines=(CSQ, ANN))),
                     annotation={"pick": "canonical", "filters": "pass", "field": "CSQ"})
    assert "KRAS" in set(df["SYMBOL"])


def test_a_maf_is_read_after_its_comment_lines(tmp_path):
    maf = write(tmp_path / "calls.maf", "#version 2.4\n#a made-up comment\n"
                "Hugo_Symbol\tVariant_Classification\tTumor_Sample_Barcode\tHGVSp_Short\n"
                "KRAS\tMissense_Mutation\tS-1\tp.G12D\nBRAF\tMissense_Mutation\tS-2\tp.V600E\n")
    df = read_source(maf)
    assert list(df.columns)[:2] == ["Hugo_Symbol", "Variant_Classification"] and df["HGVSp_Short"].tolist() == [
        "p.G12D", "p.V600E"]
    draft = yaml.safe_load(flatfile.write_draft({"results": "calls.maf"}, tmp_path / "d.yaml").read_text(encoding="utf-8"))
    assert draft["tables"]["results"]["test_id"] == {"column": "Tumor_Sample_Barcode"}
    assert draft["tables"]["results"]["gene"] == {"column": "Hugo_Symbol"}
    assert draft["tables"]["results"]["variant"] == {"column": "HGVSp_Short"}


@pytest.mark.parametrize("given, short", [("ENSP1.1:p.Gly12Asp", "p.G12D"), ("p.Leu10%3D", "p.L10="),
                                          ("p.Ile5LysfsTer3", "p.I5Kfs*3"), ("", None)])
def test_vep_protein_changes_in_short_form(given, short):
    got = short_protein_change(given)
    assert pd.isna(got) if short is None else got == short


# ---- Through a confirmed mapping ----
def results_mapping(tmp_path, results_file, test_id_column="SAMPLE", **changes):
    pd.DataFrame({"pid": ["P1", "P2"], "dx": ["2015-01-01", "2016-01-01"], "dead": ["0", "1"],
                  "days": ["900", "400"]}).to_csv(tmp_path / "patients.csv", index=False)
    pd.DataFrame({"pid": ["P1", "P2"], "sid": ["S-1", "S-2"], "days": ["30", "45"],
                  "panel": ["P", "P"]}).to_csv(tmp_path / "tests.csv", index=False)
    pd.DataFrame({"panel_id": ["P"] * 5, "gene": ["KRAS", "BRAF", "TP53", "GENE4", "GENE5"]}).to_csv(
        tmp_path / "panels.csv", index=False)
    absent = {"absent": True}
    m = {"mapping_format": 1, "status": "confirmed", "confirmed_by": "A. Data Manager", "confirmed_on": "2026-11-17",
         "missing_values": [""],
         "files": {"patients": "patients.csv", "tests": "tests.csv", "panels": "panels.csv", "results": results_file},
         "profile": {"entry_event": "report", "panel_coverage": "genes_only", "panel_sizes": "not_given",
                     "death_capture": "linked", "time_format": "intervals", "privacy_min_count": 11,
                     "histology_vocabulary": "made up", "age_top_code": None, "known_quirks": [], "sources": {}},
         "tables": {
             "patients": {f.name: absent for f in table_fields("patients")} | {
                 "patient_id": {"column": "pid"}, "vital_status": {"column": "dead"},
                 "days_to_death_or_last_followup": {"column": "days"}},
             "tests": {f.name: absent for f in table_fields("tests")} | {
                 "patient_id": {"column": "pid"}, "test_id": {"column": "sid"}, "modality": {"constant": "ngs"},
                 "entry_event": {"constant": "report"}, "days_to_entry_event": {"column": "days"},
                 "panel_id": {"column": "panel"}},
             "results": {f.name: absent for f in table_fields("results")} | {
                 "test_id": {"column": test_id_column}, "kind": {"constant": "mutation"}, "gene": {"column": "SYMBOL"},
                 "variant": {"column": "HGVSp", "transform": {"hgvsp_short": True}},
                 "variant_class": {"column": "Consequence", "transform": {"map": {
                     "missense_variant": "Missense_Mutation", "frameshift_variant": "Frame_Shift_Del"}}}},
             "panels": {f.name: absent for f in table_fields("panels")} | {
                 "panel_id": {"column": "panel_id"}, "gene": {"column": "gene"}}}}
    m.update(changes)
    path = tmp_path / "mapping.yaml"
    path.write_text(yaml.safe_dump(m, sort_keys=False), encoding="utf-8")
    return path


def test_a_vcf_reads_into_the_results_table(tmp_path):
    write(tmp_path / "calls.vcf", vcf_text())
    res = flatfile.read(results_mapping(tmp_path, {"path": "calls.vcf",
                                                  "annotation": {"pick": "canonical", "filters": "pass"}}), S)
    got = sorted(map(tuple, res.tables.results[["test_id", "gene", "variant", "variant_class"]].astype(str).values))
    assert got == [("S-1", "GENE4", "p.A10T", "Missense_Mutation"), ("S-1", "GENE5", "p.I5Kfs*3", "Frame_Shift_Del"),
                   ("S-1", "KRAS", "p.G12D", "Missense_Mutation"), ("S-2", "BRAF", "p.V600E", "Missense_Mutation"),
                   ("S-2", "GENE4", "p.A10P", "Missense_Mutation")]


def test_one_vcf_per_sample_is_read_through_a_pattern(tmp_path):
    (tmp_path / "vcf").mkdir()
    for sample, record in (("S-1", RECORDS[0]), ("S-2", RECORDS[1])):          # sites only: no sample columns
        write(tmp_path / "vcf" / f"{sample}.vcf", vcf_text([record], samples=()))
    res = flatfile.read(results_mapping(tmp_path, {"path": "vcf/*.vcf",
                                                  "annotation": {"pick": "canonical", "filters": "pass"}},
                                       test_id_column="source_file"), S)
    assert sorted(map(tuple, res.tables.results[["test_id", "variant"]].astype(str).values)) == [
        ("S-1", "p.G12D"), ("S-2", "p.V600E")]
    assert {p.name for p in res.input_files} >= {"S-1.vcf", "S-2.vcf"}


def test_a_pattern_that_matches_nothing_is_named(tmp_path):
    with pytest.raises(MappingError, match="no file matches vcf/"):
        flatfile.read(results_mapping(tmp_path, {"path": "vcf/*.vcf",
                                                "annotation": {"pick": "all", "filters": "all"}}), S)


@pytest.mark.parametrize("results, problem", [
    ("calls.vcf", "a VCF needs its annotation's reading rules"),
    ({"path": "calls.vcf", "annotation": {"pick": "best", "filters": "pass"}}, "annotation.pick must be one of"),
    ({"path": "calls.vcf", "annotation": {"pick": "all", "filters": "most"}}, "annotation.filters must be one of"),
    ({"path": "calls.vcf", "annotation": {"pick": "all", "filters": "all", "field": "EFF"}}, "annotation.field"),
    ({"path": "panels.csv", "annotation": {"pick": "all", "filters": "all"}}, "annotation applies only to a VCF"),
])
def test_a_vcfs_reading_rules_must_be_given(tmp_path, results, problem):
    write(tmp_path / "calls.vcf", vcf_text())
    with pytest.raises(MappingError, match=problem):
        flatfile.load_mapping(results_mapping(tmp_path, results))


def test_hgvsp_short_takes_true(tmp_path):
    path = results_mapping(tmp_path, {"path": "calls.vcf", "annotation": {"pick": "all", "filters": "all"}})
    m = yaml.safe_load(path.read_text(encoding="utf-8"))
    m["tables"]["results"]["variant"]["transform"] = {"hgvsp_short": "yes"}
    path.write_text(yaml.safe_dump(m, sort_keys=False), encoding="utf-8")
    with pytest.raises(MappingError, match="hgvsp_short takes true"):
        flatfile.load_mapping(path)


def test_the_draft_asks_for_the_vcf_rules_and_suggests_the_annotation_columns(tmp_path):
    write(tmp_path / "calls.vcf", vcf_text())
    draft = flatfile.write_draft({"results": "calls.vcf"}, tmp_path / "d.yaml")
    m = yaml.safe_load(draft.read_text(encoding="utf-8"))
    assert m["files"]["results"] == {"path": "calls.vcf", "annotation": {"pick": None, "filters": None}}
    r = m["tables"]["results"]
    assert (r["test_id"], r["gene"], r["variant"], r["variant_class"]) == (
        {"column": "SAMPLE"}, {"column": "SYMBOL"}, {"column": "HGVSp"}, {"column": "Consequence"})
    assert "pick: canonical / first / all" in draft.read_text(encoding="utf-8")
