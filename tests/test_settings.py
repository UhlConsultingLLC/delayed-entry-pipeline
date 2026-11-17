"""The settings: schema, loader, the CRC file and the blank template."""
import ast
import re
from pathlib import Path

import pytest
import yaml

from delayed_entry import settings as st
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


def test_the_nsclc_draft_loads_and_every_choice_awaits_review():
    s = load_settings(ROOT / "settings" / "nsclc_bpc.yaml")
    assert s.read("source.genie_bpc.cohort") == "nsclc_v2.0" and s.read("study.registered_plan") == "not registered"
    assert all("awaiting lung oncologist review" in p["made_by"] for p in s.provenance())
    assert [g["name"] for g in s.read("genotype.groups")][-1] == s.read("genotype.reference") == "No listed driver"


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
