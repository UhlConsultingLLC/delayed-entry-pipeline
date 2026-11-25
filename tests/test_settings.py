"""The settings: schema, loader, the CRC file and the blank template."""
import ast
import re
from pathlib import Path

import pandas as pd
import pytest
import yaml

from delayed_entry import settings as st
from delayed_entry.covariates import band_of, merge_bands
from delayed_entry.settings import SettingsError, blank_template, load_settings, loads_settings, unfilled_choices
from helpers import CRC_SETTINGS, ROOT, crc_settings, settings_with

TEMPLATE = ROOT / "settings" / "template_blank.yaml"


# ---- The CRC settings carry the registered plan's values (design plan section 4's table) ----
REGISTERED = {
    "entry.event": "report",
    "time.horizon_years": 3,
    "summary.survival_at_years": [1, 2, 3, 4, 5],
    "summary.median_limit_years": 5,
    "materiality.rmst_months": 2,
    "materiality.median_months": 3,
    "materiality.hr_ratio_change": 0.10,
    "materiality.rmst_gap_change_months": 2,
    "labels.clearly_material": "estimate_and_whole_interval",
    "labels.material": "estimate_only",
    "smallest_group.patients": 30,
    "biomarker.unknown_max_share": 0.30,
    "survival_after_entry.window_days": 730,
    "survival_after_entry.min_deaths": 11,
    "survival_after_entry.if_fewer": "move_entry_by_window",
    "stability.risk_set": 20,
    "uncertainty.bootstrap_resamples": 1000,
    "uncertainty.interval_level": 0.95,
    "weighting.trim_percentiles": [1, 99],
    "weighting.balance_target_smd": 0.1,
    "ph_check.p_threshold": 0.01,
    "ph_check.split_months": 12,
    "assumption_check.estimators": ["MB"],
    "privacy.min_count": 11,
    "genotype.reference": "Wild type",
    "genotype.coverage.applies_to": "all_groups",
}


@pytest.mark.parametrize("choice", list(REGISTERED))
def test_the_crc_settings_carry_the_registered_values(choice):
    assert crc_settings().read(choice) == REGISTERED[choice]


def test_the_crc_cohort_rules_and_groups_are_the_registered_ones():
    s = crc_settings()
    assert [r["id"] for r in s.read("cohort.rules")] == ["stage_iv", "histology", "age_18", "not_appendix", "has_report",
                                                       "entry_known", "exit_known", "entry_before_exit", "genotype_assignable"]
    groups = s.read("genotype.groups")
    assert [g["name"] for g in groups] == ["RAS-mutant", "BRAF V600E", "Wild type"]
    assert groups[0]["any_of"][0]["codons"] == [12, 13, 59, 61, 117, 146]
    assert s.read("genotype.coverage.required") == {"KRAS": [12, 13, 59, 61, 117, 146],
                                                    "NRAS": [12, 13, 59, 61, 117, 146], "BRAF": [600]}
    assert s.read("smallest_group.if_short") == {"RAS-mutant": "block", "BRAF V600E": "describe", "Wild type": "block"}
    assert {a["id"] for a in s.read("analyses")} == {"S1", "S2", "S3", "SA1", "SA2", "SA3", "SA4", "SA5", "SA6", "SA7",
                                                     "SA8", "SA9"}


def test_every_crc_choice_says_who_made_it_and_its_source():
    for p in crc_settings().provenance():
        assert p["made_by"].strip() and p["source"].strip(), p["choice"]
        assert p["origin"].startswith("settings file, line ")


def test_the_crc_settings_use_no_default():
    assert crc_settings().defaults_used() == []


def test_reading_a_choice_is_recorded():
    s = crc_settings()
    s.read("privacy.min_count")
    s.read("time.horizon_years")
    assert s.reads() == ["privacy.min_count", "time.horizon_years"]


def test_the_nsclc_settings_load_and_every_choice_is_decided():
    """Every numbered question, the coverage setting, C1, the fusion rows, question 2's scope, item C2, the
    choices never put as questions and the registered plan are decided (the last on 2026-11-24). The registered
    plan points to the provenance log, where the registration's link goes, so the link stays out of the file
    that is registered and locked. The log is published only after the lock, so it isn't read here."""
    s = load_settings(ROOT / "settings" / "nsclc_bpc.yaml")
    assert s.read("source.genie_bpc.cohort") == "nsclc_v2.0"
    plan = s.read("study.registered_plan")
    assert "docs/provenance/nsclc_bpc_v2.0.md" in plan and "osf.io" not in plan
    assert not {p["choice"]: p["made_by"] for p in s.provenance() if p["made_by"] != "NSCLC study author"}
    assert not [a for a in s.read("analyses") if a["kind"].startswith("biomarker")]   # question 10: PD-L1 out
    assert "oncologist" not in (ROOT / "settings" / "nsclc_bpc.yaml").read_text(encoding="utf-8")
    assert [g["name"] for g in s.read("genotype.groups")][-1] == s.read("genotype.reference") == "No recognized driver"


def test_the_nsclc_answers_to_questions_7_9_and_13_and_the_coverage_setting():
    """Questions 7 (A), 9 (A) and 13 (B) and the coverage setting (A), decided 2026-11-24. Question 14 (B) changes
    no setting: its two lines are in docs/nsclc_settings_review.md."""
    s = load_settings(ROOT / "settings" / "nsclc_bpc.yaml")
    assert s.read("time.horizon_years") == 3
    assert s.read("source.genie_bpc.cancer_type_labels") == ["Non Small Cell Lung Cancer", "Lung Cancer, NOS"]
    assert s.read("genotype.coverage.when_only_genes_listed") == "accept_gene_level"
    assert s.read("genotype.described_only") == ["Other recognized driver"]
    assert s.read("smallest_group.patients") == 30
    assert s.read("smallest_group.if_short") == {"EGFR common": "block", "ALK fusion": "describe",
                                                "KRAS mutant": "describe", "No recognized driver": "block"}
    still_open = " ".join(p["note"] for p in s.provenance() if "TO DECIDE" in p.get("note", ""))
    assert not re.search(r"question (7|9|13)\b", still_open) and "numbered questions" not in still_open


READING_RULE = ('for the HR-ratio and RMST-gap comparison rows only "clearly material" counts as material bias; '
                '"material" is reported with its interval as "past the bar; chance not excluded"; a group\'s own '
                'overstatement is read as in the CRC plan, section 11')


def test_the_nsclc_answers_to_questions_8_and_11_and_item_c1():
    """Question 8 (option B': the CRC bars and labels, with a stated reading), question 11 and item C1 (four rows),
    decided 2026-11-24. The reading rule is in both labels' source, since the locked plan prints each choice's source,
    not its note."""
    s = load_settings(ROOT / "settings" / "nsclc_bpc.yaml")
    assert {k: s.read(f"materiality.{k}") for k in ("rmst_months", "median_months", "hr_ratio_change",
                                                    "rmst_gap_change_months")} == \
        {"rmst_months": 2, "median_months": 3, "hr_ratio_change": 0.10, "rmst_gap_change_months": 2}
    assert s.read("labels.clearly_material") == "estimate_and_whole_interval"
    assert s.read("labels.material") == "estimate_only"
    for label in ("labels.clearly_material", "labels.material"):
        source = s.provenance([label])[0]["source"]
        assert READING_RULE in source and "design difference from the CRC study" in source, label
    assert ("(Ellis LM et al., J Clin Oncol 2014;32:1277-1280, as reprinted in Cortinovis D et al., J Thorac Dis "
            "2015;7:E541-E547, Table 1, and adapted in Fenchel K et al.") in \
        s.provenance(["materiality.median_months"])[0]["source"]
    variables = {v["field"]: v for v in s.read("adjustment.variables")}
    assert list(variables) == ["age_at_diagnosis", "sex", "center", "smoking_at_diagnosis", "histology",
                               "year_of_diagnosis"]                                    # race described only
    smoking = variables["smoking_at_diagnosis"]
    assert smoking["reference"] == "Never" and smoking["merge"] == {
        "Never": ["Never used"], "Former": ["Former user (quit >1 year)", "Former user (unknown time)"],
        "Current": ["Current user", "Former user (quit <1 year)"]}
    release = yaml.safe_load((ROOT / "src" / "delayed_entry" / "adapters" / "genie_bpc" / "cohorts" / "nsclc_v2.0"
                              / "expected_release.yml").read_text(encoding="utf-8"))
    allowed = release["values"]["raw_cancer"]["ca_lung_cigarette"]["allowed"]
    merged = [x for xs in smoking["merge"].values() for x in xs]
    assert sorted(merged + ["Unknown"]) == sorted(allowed)      # spelled as the release does; Unknown is missing
    assert variables["histology"] == {"field": "histology", "form": "categorical", "reference": "Adenocarcinoma",
                                      "merge": {"Adenocarcinoma": ["Adenocarcinoma"],
                                                "Carcinoma or other": ["Carcinoma", "Other histologies/mixed tumor"]}}
    assert s.read("missing_data.imputation_model") == {
        "engine": "chained_equations", "outcome_terms": ["event", "cumulative_hazard_corrected", "entry_time"],
        "include_group": True,
        "methods": {"linear": "pmm", "spline": "pmm", "categorical": "polyreg", "bands": "polyreg"},
        "pmm_donors": 5, "iterations": 5, "pooling": "rubins_rules", "bootstrap": "impute_each_resample_once",
        "point_estimates": "average_across_imputations"}
    assert s.read("genotype.reference_drugs") == {
        "EGFR": ["Gefitinib", "Afatinib", "Osimertinib", "Icotinib", "Rociletinib"],
        "EGFR (erlotinib)": ["Erlotinib"],
        "ALK or ROS1": ["Crizotinib", "Ceritinib", "Alectinib", "Brigatinib", "Lorlatinib", "Ensartinib",
                        "Entrectinib"],
        "BRAF V600E": ["Dabrafenib", "Trametinib", "Vemurafenib"]}
    drugs = s.provenance(["genotype.reference_drugs"])[0]["source"]
    assert "until October 18, 2016" in drugs and "Approval dates (Drugs@FDA approval letters)" in drugs
    assert "NDA 021743, supplement 25, approved October 18, 2016" in drugs and "to confirm" not in drugs
    assert not [p["choice"] for p in s.provenance() if "TO DECIDE" in p.get("note", "") + p["source"]]


C2_KEPT = {"summary.survival_at_years": [1, 2, 3, 4, 5], "summary.median_limit_years": 5,
           "survival_after_entry.window_days": 730, "survival_after_entry.min_deaths": 11,
           "survival_after_entry.if_fewer": "move_entry_by_window", "stability.risk_set": 20,
           "uncertainty.bootstrap_resamples": 1000, "uncertainty.seed": 20261113, "uncertainty.interval_level": 0.95,
           "missing_data.complete_case_max_share": 0.05, "missing_data.imputations": 20,
           "weighting.variables": "same_as_adjustment", "weighting.trim_percentiles": [1, 99],
           "weighting.balance_target_smd": 0.1, "ph_check.split_months": 12,
           "assumption_check.estimators": ["MB", "IPW1", "IPW2"], "privacy.min_count": 11}
NEVER_ASKED = ("study.name", "study.cancer", "source.adapter", "source.genie_bpc.cohort",
               "source.genie_bpc.cancer_site_pattern", "source.genie_bpc.withheld_survival", "cohort.rules",
               "entry.event", "entry.acceptable_substitutes", "entry.modalities", "entry.first_test_order",
               "entry.before_diagnosis", "genotype.results_from", "genotype.variant_text",
               "genotype.coverage.required", "time.days_per_year", "adjustment.band_min_patients",
               "ph_check.p_threshold", "uncertainty.when_a_resample_has_no_value", "analyses")


def test_the_nsclc_answers_to_item_c2_and_the_choices_never_put_as_questions():
    """Item C2 (the choices carried from the CRC plan) and the choices never put as numbered questions, decided
    2026-11-24: every value kept, except the last year band, now open ("2016 or later"), and study.cancer, which names
    the three centers (text only)."""
    s = load_settings(ROOT / "settings" / "nsclc_bpc.yaml")
    assert {c: s.read(c) for c in C2_KEPT} == C2_KEPT
    bands = s.read("adjustment.year_bands")
    assert bands == [{"label": "2013 or earlier", "to": 2013}, {"label": "2014-2015", "from": 2014, "to": 2015},
                     {"label": "2016 or later", "from": 2016}]
    # the year of diagnosis is estimated, so a 2017 diagnosis can be coded 2018: it lands in a band, never in none
    assert band_of(pd.Series([2009, 2013, 2014, 2015, 2016, 2017, 2018]), bands).tolist() == [
        "2013 or earlier", "2013 or earlier", "2014-2015", "2014-2015", "2016 or later", "2016 or later",
        "2016 or later"]
    assert merge_bands({"2013 or earlier": 12, "2014-2015": 400, "2016 or later": 500}, bands, 30) == {
        "2013 or earlier": "2015 or earlier", "2014-2015": "2015 or earlier", "2016 or later": "2016 or later"}
    assert s.read("study.cancer") == ("nonsquamous non-small cell lung cancer, stage IV at diagnosis, first "
                                      "sequenced at DFCI, MSK or VICC (reports 2014-2017)")
    assert next(a for a in s.read("analyses") if a["id"] == "SA1")["keep"] == {"within_months": 6}
    assert [r["id"] for r in s.read("cohort.rules")] == ["stage_iv", "histology", "age_18", "has_report", "entry_known",
                                                       "exit_known", "entry_before_exit", "genotype_assignable"]
    assert s.read("genotype.coverage.required")["ALK"] == "whole_gene"
    sources = {p["choice"]: p["source"] for p in s.provenance()}
    notes = {p["choice"]: p.get("note", "") for p in s.provenance()}
    for choice in C2_KEPT:
        assert sources[choice].lower().startswith("item c2 (the dry run), decided 2026-11-24 (kept): "), choice
    for choice in NEVER_ASKED:
        assert "confirmed 2026-11-24 " in sources[choice].lower(), choice
    assert sources["adjustment.year_bands"].startswith("Item C2 (the dry run), decided 2026-11-24 (the last band open)")
    assert "item C2 (the dry run), kept" in sources["analyses"] and "182.6 days" in sources["analyses"]
    assert "If that exclusion is applied as written" in notes["source.genie_bpc.withheld_survival"]
    assert "skipped since question 9" in notes["source.genie_bpc.cancer_site_pattern"]
    assert not [c for c in sources if "pending" in sources[c] + notes[c]]


OPEN_CITATION = ("to check (C4)", "to confirm (C4)", "to check, C4", "to confirm in v2.0-public",
                 "as reported by a secondary source")
ERBB2_EXCEPT = ["p.S310A", "p.S310C", "p.S310P", "p.S310T",
                "p.V659A", "p.V659F", "p.V659G", "p.V659I", "p.V659L",
                "p.G660A", "p.G660C", "p.G660R", "p.G660S", "p.G660V",
                "p.L755F", "p.L755M", "p.L755V", "p.L755W",
                "p.D769A", "p.D769E", "p.D769G", "p.D769V",
                "p.V777A", "p.V777E", "p.V777G", "p.V777M",
                "p.V842A", "p.V842E", "p.V842G", "p.V842L"]


def test_the_nsclc_checks_before_the_lock_leave_no_citation_open_and_change_no_rule():
    """The checks before the lock (transcripts, the MET splice guard, ROS1/RET detection and the citations; the dry
    run's C5), decided 2026-11-24: sources and notes only. No source or note still calls a citation to check or to confirm; the ERBB2
    except lists stay as they were (the excepted changes a cited source supports are stated, not counted), and the
    MET rules still read Splice_Site, not Splice_Region (stated)."""
    s = load_settings(ROOT / "settings" / "nsclc_bpc.yaml")
    texts = {p["choice"]: p["source"] + " " + p.get("note", "") for p in s.provenance()}
    assert not [(choice, m) for choice, text in texts.items() for m in OPEN_CITATION if m in text]
    other = next(g for g in s.read("genotype.groups") if g["name"] == "Other recognized driver")
    erbb2 = [(r["variant_classes"], r.get("codons"), r.get("codon_ranges"), r.get("except_variants"))
             for r in other["any_of"] if r.get("genes") == ["ERBB2"]]
    assert erbb2 == [(["In_Frame_Ins"], None, [[770, 831]], None),
                     (["Missense_Mutation"], [310, 659, 660, 755, 769, 777, 842], None, ERBB2_EXCEPT)]
    met = [r for r in other["any_of"] if r.get("genes") == ["MET"]]
    assert [r["variant_classes"] for r in met] == [["Splice_Site"], ["Splice_Site"], ["Missense_Mutation"]]
    note = s.provenance(["genotype.groups"])[0]["note"]
    assert "among them G660R, V777M, S310P, L755M and L755W" in note and "classes Splice_Region or Intron" in note
    assert "NM_000245 / ENST00000397752 (MET)" in note and "NM_004448 / ENST00000269571 (ERBB2)" in note
    coverage = s.provenance(["genotype.coverage.required"])[0]["note"]
    assert "RET intron 7 isn't on MSK-IMPACT341" in coverage and "ROS1 introns 32-34 are targeted" in coverage


def test_the_nsclc_groups_note_matches_the_download_check_after_the_build():
    """The two settings-note fixes after the build (the dry run's C6), decided 2026-11-24: the download check applies
    except lists (item p's per-rule count) but never checks them against the release, so only the intragenic texts
    are checked by eye; BRAF's transcript is named with the others. MET's three rules, with their except list, are
    pinned whole, since no other test pins that list."""
    s = load_settings(ROOT / "settings" / "nsclc_bpc.yaml")
    note = s.provenance(["genotype.groups"])[0]["note"]
    assert "never reads except lists" not in note and "never checks except lists against the release" in note
    assert "so the intragenic texts are checked by eye against the listed descriptions" in note
    assert "NM_004333 / ENST00000288602 (BRAF)" in note and "every MET, EGFR, ERBB2, KRAS and BRAF row" in note
    assert "V600E is numbered p.V640E" in note and "except-list entry would go unnoticed" in note
    other = next(g for g in s.read("genotype.groups") if g["name"] == "Other recognized driver")
    met = [r for r in other["any_of"] if "MET" in r["genes"]]
    assert met == [
        {"kind": "mutation", "genes": ["MET"], "variant_classes": ["Splice_Site"], "codon_ranges": [[963, 1010]]},
        {"kind": "mutation", "genes": ["MET"], "variant_classes": ["Splice_Site"], "variants": ["p.X1010_splice"]},
        {"kind": "mutation", "genes": ["MET"], "variant_classes": ["Missense_Mutation"], "codons": [1010],
         "except_variants": ["p.D1010A", "p.D1010E", "p.D1010G", "p.D1010V"]}]


# ---- The blank template ----
def test_the_shipped_template_matches_the_schema():
    assert TEMPLATE.read_text(encoding="utf-8") == blank_template()


def test_the_blank_template_makes_every_choice_explicit():
    with pytest.raises(SettingsError) as e:
        load_settings(TEMPLATE)
    assert set(unfilled_choices(e.value)) == {c.path for c in st.SCHEMA}
    assert all(re.search(r"line \d+", line) for line in str(e.value).splitlines()[1:])


def test_the_template_shows_each_default_and_its_source():
    text = blank_template()
    for choice, d in st.load_defaults().items():
        assert d["source"] in text, choice


def test_a_filled_in_template_loads(tmp_path):
    """Filling every blank of the template with the CRC values gives the CRC settings. An optional choice the CRC
    file leaves out is deleted from the template, as its comment says."""
    crc = yaml.safe_load(CRC_SETTINGS.read_text(encoding="utf-8"))
    text = blank_template()
    data = yaml.safe_load(text)
    for c in st.SCHEMA:
        src, dst = crc, data
        parts = c.path.split(".")
        for part in parts[:-1]:
            src, dst = src[part], dst[part]
        if c.optional and parts[-1] not in src:
            del dst[parts[-1]]
            continue
        dst[parts[-1]].update(value=src[parts[-1]]["value"], made_by="a test", source="the CRC file")
    s = loads_settings(yaml.safe_dump(data, sort_keys=False))
    assert all(s.read(c.path) == crc_settings().read(c.path) for c in st.SCHEMA if crc_settings().has(c.path))


def test_an_optional_choice_may_be_left_out_but_not_left_blank(tmp_path):
    choice = "missing_data.imputation_model"
    assert st.SCHEMA_BY_PATH[choice].optional and not crc_settings().has(choice)
    model = st.SCHEMA_BY_PATH[choice].example
    s = settings_with_entry(tmp_path, {"value": model, "made_by": "a test", "source": "made up"})
    assert s.read(choice)["methods"]["linear"] == "pmm"
    with pytest.raises(SettingsError, match="'missing_data.imputation_model' is not filled in"):
        settings_with_entry(tmp_path, {"value": None, "made_by": "a test", "source": "made up"})
    with pytest.raises(SettingsError, match="must be one of: rubins_rules"):
        settings_with_entry(tmp_path, {"value": {**model, "pooling": "average"}, "made_by": "a test", "source": "made up"})


def settings_with_entry(tmp_path, entry):
    data = yaml.safe_load(CRC_SETTINGS.read_text(encoding="utf-8"))
    data["missing_data"]["imputation_model"] = entry
    p = tmp_path / "with_imputation.yaml"
    p.write_bytes(yaml.safe_dump(data, sort_keys=False).encode("utf-8"))
    return load_settings(p)


# ---- What the loader refuses ----
def crc_text() -> str:
    return CRC_SETTINGS.read_text(encoding="utf-8")


def line_of(text: str, needle: str) -> int:
    return next(i for i, line in enumerate(text.splitlines(), start=1) if needle in line)


def test_a_blank_value_names_its_line():
    text = crc_text().replace("    value: 730\n", "    value:\n")
    with pytest.raises(SettingsError, match=r"line \d+: 'survival_after_entry.window_days' is not filled in"):
        loads_settings(text)


def test_a_missing_choice_is_named():
    data = yaml.safe_load(crc_text())
    del data["stability"]
    with pytest.raises(SettingsError, match="'stability.risk_set' is missing"):
        loads_settings(yaml.safe_dump(data))


def test_an_unknown_key_is_named_with_its_line():
    text = crc_text().replace("smallest_group:\n  patients:", "smallest_group:\n  patiens:\n    value: 3\n  patients:")
    with pytest.raises(SettingsError) as e:
        loads_settings(text)
    assert f"line {line_of(text, 'patiens:')}: 'smallest_group.patiens' is not a known choice" in str(e.value)


def test_a_misspelled_entry_key_is_caught():
    text = crc_text().replace("    value: 20\n    made_by", "    value: 20\n    made_bye: someone\n    made_by", 1)
    with pytest.raises(SettingsError, match="made_bye"):
        loads_settings(text)


def test_a_value_of_the_wrong_kind_names_its_line():
    text = crc_text().replace("    value: 0.30\n", "    value: thirty percent\n", 1)
    with pytest.raises(SettingsError) as e:
        loads_settings(text)
    assert f"line {line_of(text, 'thirty percent')}: 'biomarker.unknown_max_share' must be a number" in str(e.value)


def test_a_share_above_one_is_refused():
    with pytest.raises(SettingsError, match="must be at most 1"):
        loads_settings(crc_text().replace("    value: 0.30\n", "    value: 30\n", 1))


def test_a_problem_inside_a_list_points_at_the_item(tmp_path):
    text = crc_text().replace("      kind: subset\n      keep: {within_months: 6}", "      kind: subset\n      keep: {within_months: -6}")
    with pytest.raises(SettingsError) as e:
        loads_settings(text)
    assert f"line {line_of(text, '- id: SA1')}: 'analyses' [3]" in str(e.value)


def test_who_made_a_choice_is_required():
    text = crc_text().replace("    value: 730\n    made_by: CRC BPC study author\n", "    value: 730\n")
    with pytest.raises(SettingsError, match="'survival_after_entry.window_days' needs 'made_by'"):
        loads_settings(text)


def test_the_source_of_a_choice_is_required():
    text = re.sub(r"(    value: 730\n    made_by: CRC BPC study author\n)    source: [^\n]*\n", r"\1", crc_text())
    with pytest.raises(SettingsError, match="'survival_after_entry.window_days' needs 'source'"):
        loads_settings(text)


def test_the_reference_must_be_a_group(tmp_path):
    with pytest.raises(SettingsError, match="'Wildtype' is not one of the groups"):
        settings_with(tmp_path, {"genotype.reference": "Wildtype"})


def test_every_structural_check_must_be_placed_in_the_rules(tmp_path):
    rules = [r for r in crc_settings().read("cohort.rules") if r["id"] != "entry_before_exit"]
    with pytest.raises(SettingsError, match="'entry_before_follow_up_end' must appear exactly once"):
        settings_with(tmp_path, {"cohort.rules": rules})


def test_a_rule_must_say_what_a_blank_means(tmp_path):
    rules = crc_settings().read("cohort.rules")
    del rules[0]["all_of"][0]["missing"]
    with pytest.raises(SettingsError, match="say what a blank means"):
        settings_with(tmp_path, {"cohort.rules": rules})


@pytest.mark.parametrize("ranges, message", [
    ([[478, 439]], r"runs backward: write \[439, 478\]"),
    ([[439, 460, 478]], r"must be two codons, \[first, last\], not 3"),
    ([[439]], "must list at least 2"),
])
def test_a_codon_range_must_be_first_then_last(tmp_path, ranges, message):
    """A reversed range would silently match nothing, and a three-number one stopped the screen."""
    analyses = crc_settings().read("analyses")
    rule = next(a["keep"]["if_result"] for a in analyses if "codon_ranges" in a.get("keep", {}).get("if_result", {}))
    rule["codon_ranges"] = ranges
    with pytest.raises(SettingsError, match=message):
        settings_with(tmp_path, {"analyses": analyses})


def test_a_one_codon_range_is_accepted(tmp_path):
    analyses = crc_settings().read("analyses")
    rule = next(a["keep"]["if_result"] for a in analyses if "codon_ranges" in a.get("keep", {}).get("if_result", {}))
    rule["codon_ranges"] = [[600, 600]]
    assert settings_with(tmp_path, {"analyses": analyses}).read("analyses")


def test_if_short_must_name_each_group(tmp_path):
    with pytest.raises(SettingsError, match="must name each group that is compared exactly once"):
        settings_with(tmp_path, {"smallest_group.if_short": {"RAS-mutant": "block", "Wild type": "block"}})


def test_a_variant_must_belong_to_a_biomarker_analysis(tmp_path):
    analyses = crc_settings().read("analyses")
    next(a for a in analyses if a["id"] == "SA8")["of"] = "SA1"
    with pytest.raises(SettingsError, match="SA8: 'of' must name a biomarker analysis"):
        settings_with(tmp_path, {"analyses": analyses})


def test_genie_bpc_options_are_refused_for_another_adapter(tmp_path):
    with pytest.raises(SettingsError, match="applies only when source.adapter is genie_bpc"):
        settings_with(tmp_path, {"source.adapter": "flat_file"})


def test_bad_yaml_names_its_line():
    with pytest.raises(SettingsError, match=r"line \d+: "):
        loads_settings("settings_format: 1\nstudy: [unclosed\n")


def test_a_key_written_twice_is_refused():
    """Plain YAML would quietly keep the second value."""
    text = crc_text().replace("    value: 730\n", "    value: 730\n    value: 5\n", 1)
    with pytest.raises(SettingsError, match=r"line \d+: the key 'value' appears twice"):
        loads_settings(text)
    with pytest.raises(SettingsError, match="the key 'privacy' appears twice"):
        loads_settings(crc_text() + "privacy:\n  min_count:\n    value: 5\n")


def test_a_value_yaml_reads_as_true_or_false_must_be_quoted(tmp_path):
    rules = crc_settings().read("cohort.rules")
    rules[1]["all_of"][0]["values"] = [True]          # what an unquoted Yes becomes
    with pytest.raises(SettingsError, match="put it in quotes"):
        settings_with(tmp_path, {"cohort.rules": rules})


def test_a_value_of_the_wrong_kind_for_its_field_is_refused(tmp_path):
    rules = crc_settings().read("cohort.rules")
    rules.insert(0, {"id": "year", "label": "Known year", "all_of": [
        {"field": "year_of_diagnosis", "op": "not_in", "values": ["unknown"], "missing": "fail"}]})
    with pytest.raises(SettingsError, match="holds numbers, so its values must be numbers"):
        settings_with(tmp_path, {"cohort.rules": rules})
    rules[0]["all_of"][0] = {"field": "histology", "op": "ge", "value": 3, "missing": "fail"}
    with pytest.raises(SettingsError, match="can't be compared with 'ge'"):
        settings_with(tmp_path, {"cohort.rules": rules}, "two.yaml")


@pytest.mark.parametrize("where", ["rule", "subset", "adjustment", "stratum"])
def test_an_outcome_can_never_choose_who_is_analyzed(tmp_path, where):
    s = crc_settings()
    rules, analyses, adjust = s.read("cohort.rules"), s.read("analyses"), s.read("adjustment.variables")
    cond = {"field": "days_to_death_or_last_followup", "op": "ge", "value": 30, "missing": "fail"}
    if where == "rule":
        rules[0]["all_of"].append(cond)
    elif where == "subset":
        next(a for a in analyses if a["id"] == "SA7")["keep"] = cond
    elif where == "adjustment":
        adjust.append({"field": "vital_status", "form": "categorical"})
    else:
        next(a for a in analyses if a["id"] == "SA2")["stratify_by"] = "vital_status"
    with pytest.raises(SettingsError, match="outcome"):
        settings_with(tmp_path, {"cohort.rules": rules, "analyses": analyses, "adjustment.variables": adjust})


def test_an_overlap_must_name_groups_with_rules(tmp_path):
    analyses = crc_settings().read("analyses")
    next(a for a in analyses if a["id"] == "SA3")["keep"] = {"exclude_overlap": ["RAS-mutant", "Wild type"]}
    with pytest.raises(SettingsError, match="'Wild type' is not a genotype group with its own rules"):
        settings_with(tmp_path, {"analyses": analyses})


def test_a_variant_cannot_run_without_its_parent(tmp_path):
    analyses = crc_settings().read("analyses")
    next(a for a in analyses if a["id"] == "S1")["run"] = False
    with pytest.raises(SettingsError, match="SA8 runs, but S1, which it varies, is switched off"):
        settings_with(tmp_path, {"analyses": analyses})


def test_the_fallback_may_stop_only_the_biomarkers_own_variants(tmp_path):
    analyses = crc_settings().read("analyses")
    next(a for a in analyses if a["id"] == "S1")["not_run_under_fallback"] = ["SA8", "SA1"]
    with pytest.raises(SettingsError, match="'SA1' is not one of its own variants"):
        settings_with(tmp_path, {"analyses": analyses})


def test_a_biomarkers_reference_must_be_one_of_its_groups(tmp_path):
    analyses = crc_settings().read("analyses")
    next(a for a in analyses if a["id"] == "S1")["reference"] = "MSS"
    with pytest.raises(SettingsError, match=r"S1: the reference must be one of its groups \(MSI-H/dMMR, MSS/pMMR\)"):
        settings_with(tmp_path, {"analyses": analyses})


def test_an_outcome_stratum_is_named_once_and_no_reference_is_asked_for(tmp_path):
    analyses = crc_settings().read("analyses")
    next(a for a in analyses if a["id"] == "SA2")["stratify_by"] = "vital_status"
    with pytest.raises(SettingsError) as e:
        settings_with(tmp_path, {"analyses": analyses})
    assert "SA2: an outcome can't be a stratum" in str(e.value) and "reference" not in str(e.value)


def test_the_resolved_hash_changes_with_the_defaults_file(tmp_path):
    s1 = settings_with(tmp_path, {"stability.risk_set": "default"})
    other = {**st.load_defaults(), "stability.risk_set": {"value": 25, "source": "a changed defaults file"}}
    s2 = st.load_settings(tmp_path / "changed.yaml", defaults=other)
    assert s1.sha256 == s2.sha256 and s1.resolved_sha256 != s2.resolved_sha256


# ---- Defaults ----
def test_a_default_is_taken_from_the_defaults_file_and_listed(tmp_path):
    s = settings_with(tmp_path, {"stability.risk_set": "default", "uncertainty.interval_level": "default"})
    assert s.read("stability.risk_set") == 20
    used = {d["choice"]: d for d in s.defaults_used()}
    assert set(used) == {"stability.risk_set", "uncertainty.interval_level"}
    assert used["stability.risk_set"]["origin"].startswith("default (defaults.yml)")
    assert "CRC BPC plan" in used["stability.risk_set"]["source"]


def test_a_choice_without_a_default_must_be_filled_in(tmp_path):
    with pytest.raises(SettingsError, match="'time.horizon_years' has no default"):
        settings_with(tmp_path, {"time.horizon_years": "default"})


def test_every_default_is_a_valid_value_with_a_source():
    for choice, d in st.load_defaults().items():
        assert st.SCHEMA_BY_PATH[choice].check(d["value"]) == [] and d["source"].strip()


def test_study_specific_choices_have_no_default():
    for choice in ("cohort.rules", "genotype.groups", "genotype.reference", "genotype.coverage.required",
                   "time.horizon_years", "analyses", "privacy.min_count", "smallest_group.if_short"):
        assert choice not in st.load_defaults()


# ---- No clinical or statistical value in the code ----
FORBIDDEN_NUMBERS = {30, 0.3, 730, 11, 1000, 0.95, 182.625, 89, 20261108, 365.25, 0.05, 0.01}
FORBIDDEN_TEXT = ("KRAS", "NRAS", "BRAF", "V600E", "Stage IV", "Adenocarcinoma", "Squamous", "C181", "Colon",
                  "MSI-H", "MSS", "pMMR", "dMMR", "RAS-mutant", "Wild type")


def code_constants(path: Path):
    """Numbers and strings in the code itself: docstrings and comments left out."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    docstrings = {id(n.body[0].value) for n in ast.walk(tree)
                  if isinstance(n, (ast.Module, ast.FunctionDef, ast.ClassDef, ast.AsyncFunctionDef))
                  and n.body and isinstance(n.body[0], ast.Expr) and isinstance(n.body[0].value, ast.Constant)}
    for n in ast.walk(tree):
        if isinstance(n, ast.Constant) and id(n) not in docstrings:
            yield n.value


def test_no_registered_value_is_written_into_the_code():
    for path in (ROOT / "src" / "delayed_entry").rglob("*.py"):
        for v in code_constants(path):
            if isinstance(v, (int, float)) and not isinstance(v, bool) and path.name != "xlsx.py":   # calendar arithmetic
                assert v not in FORBIDDEN_NUMBERS, f"{path.name}: {v}"
            if isinstance(v, str) and path.name != "check_download.py":   # its SQL checks source units
                assert not any(t in v for t in FORBIDDEN_TEXT), f"{path.name}: {v[:60]}"
