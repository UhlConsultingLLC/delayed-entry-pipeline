"""A study's settings: every clinical and statistical choice, with who made it and the source behind it.

    from delayed_entry.settings import load_settings
    s = load_settings("settings/crc_bpc.yaml")
    s.read("smallest_group.patients")      # 30, and the screen records that it read it
    s.provenance()                          # every choice: value, who, source, and where it came from

No clinical or statistical value is written into the code. Each lives in the settings file (YAML), one
entry per choice:

    smallest_group:
      patients:
        value: 30
        made_by: the study's statistician
        source: registered plan, section 5 (D1)

The schema below says which choices exist, what form each takes and, in plain words, what it controls.
Every choice must be written out: a missing choice, a blank value or an unknown key is an error that
names the line. A value may be the word `default`, which takes the value in defaults.yml (shipped with
the package, each with its source); the screen and the locked plan list every default used. Choices
that are specific to a study (who is in the cohort, the groups, the horizon, the privacy minimum) have
no default and must be filled in.

`blank_template()` writes the template in which every choice is blank.
"""
from __future__ import annotations

import copy
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import yaml

from delayed_entry.tables import (BOOLEAN, ENTRY_EVENTS, INTEGER, NUMBER, OUTCOME_FIELDS, RESULT_KINDS, TEXT,
                                  field_spec)

SETTINGS_FORMAT = 1
DEFAULTS_FILE = Path(__file__).with_name("defaults.yml")
DEFAULT_WORD = "default"
ENTRY_KEYS = {"value", "made_by", "source", "note"}


class SettingsError(ValueError):
    """A settings file that can't be used. The message lists every problem, each with its line."""


class _UniqueKeyLoader(yaml.SafeLoader):
    """YAML's safe loader, except that a key written twice in one mapping is an error: plain YAML would
    silently keep the last one, so a second `value:` or a second section would override the first."""

    def construct_mapping(self, node, deep=False):
        seen = set()
        for key_node, _ in node.value:
            key = self.construct_object(key_node, deep=deep)
            if key in seen:
                raise yaml.constructor.ConstructorError(None, None, f"the key '{key}' appears twice",
                                                        key_node.start_mark)
            seen.add(key)
        return super().construct_mapping(node, deep)


# ---------------------------------------------------------------------------
# Validators. Each returns a list of problems ("" for the value itself, or a path inside it).
# ---------------------------------------------------------------------------
Validator = Callable[[Any], list]


def text() -> Validator:
    return lambda v: [] if isinstance(v, str) and v.strip() else ["must be text"]


def number(lo=None, hi=None, integer=False, lo_open=False) -> Validator:
    def check(v):
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            return ["must be a number"]
        if integer and not float(v).is_integer():
            return ["must be a whole number"]
        if lo is not None and (v <= lo if lo_open else v < lo):
            return [f"must be {'more than' if lo_open else 'at least'} {lo}"]
        if hi is not None and v > hi:
            return [f"must be at most {hi}"]
        return []
    return check


def boolean() -> Validator:
    return lambda v: [] if isinstance(v, bool) else ["must be true or false"]


def enum(*options) -> Validator:
    return lambda v: [] if v in options else [f"must be one of: {', '.join(map(str, options))}"]


def list_of(item: Validator, min_len=1, unique=False) -> Validator:
    def check(v):
        if not isinstance(v, list):
            return ["must be a list"]
        if len(v) < min_len:
            return [f"must list at least {min_len}"]
        out = [f"[{i}] {p}" for i, x in enumerate(v) for p in item(x)]
        if unique and len({repr(x) for x in v}) != len(v):
            out.append("has a repeated entry")
        return out
    return check


def mapping_of(key: Validator, value: Validator, min_len=1) -> Validator:
    def check(v):
        if not isinstance(v, dict):
            return ["must be a mapping"]
        if len(v) < min_len:
            return [f"must have at least {min_len} entries"]
        return [f"{k}: {p}" for k, x in v.items() for p in key(k) + value(x)]
    return check


def record(required: dict[str, Validator], optional: dict[str, Validator] | None = None) -> Validator:
    """A mapping with these keys and no others."""
    optional = optional or {}

    def check(v):
        if not isinstance(v, dict):
            return ["must be a mapping"]
        out = [f"missing key '{k}'" for k in required if k not in v]
        out += [f"unknown key '{k}'" for k in v if k not in required and k not in optional]
        for k, f in {**required, **optional}.items():
            if k in v:
                out += [f"{k}: {p}" for p in f(v[k])]
        return out
    return check


def one_of(*records: tuple[str, Validator]) -> Validator:
    """A mapping identified by the key it carries (the first of each pair)."""
    def check(v):
        if not isinstance(v, dict):
            return ["must be a mapping"]
        for marker, f in records:
            if marker in v:
                return f(v)
        return ["must have one of the keys: " + ", ".join(m for m, _ in records)]
    return check


# ---- The parts of complex choices ----
NORMALIZE = ("none", "trim", "letters_digits_upper")
CONDITION_OPS = ("in", "not_in", "starts_with_any", "not_starts_with_any", "ge", "gt", "le", "lt", "known")
STRUCTURAL_CHECKS = {
    "has_entry_test": "the patient has a test of an entry modality",
    "entry_time_known": "the first test's time is known",
    "outcome_recorded": "the patient's outcome is recorded (the screen sees only that it is, never what it is)",
    "entry_before_follow_up_end": "the first test did not come back after follow-up ended, as the source records it",
    "genotype_assignable": "the first test's panel covers what the genotype groups need",
}
PATIENT_FLAGS = ("multiple_cohort_cancers", "other_primary_same_cancer")


def _field_condition(v):
    base = record({"field": text(), "op": enum(*CONDITION_OPS)},
                  {"values": list_of(lambda x: [] if isinstance(x, (str, int, float)) else ["must be text or a number"]),
                   "value": number(), "normalize": enum(*NORMALIZE), "missing": enum("pass", "fail")})(v)
    if base:
        return base
    if v["field"] in OUTCOME_FIELDS:
        return [f"'{v['field']}' is an outcome: it can't decide who is analyzed (the screen never reads outcomes)"]
    spec = field_spec("patients", v["field"])
    kind = spec.kind if spec else None
    values = v.get("values", [])
    if kind in (NUMBER, INTEGER) and any(isinstance(x, bool) or not isinstance(x, (int, float)) for x in values):
        return [f"'{v['field']}' holds numbers, so its values must be numbers"]
    if kind == BOOLEAN and any(not isinstance(x, bool) for x in values):
        return [f"'{v['field']}' holds true or false, so its values must be true or false"]
    if kind != BOOLEAN and any(isinstance(x, bool) for x in values):
        return ["a value was read as true or false: put it in quotes (for example 'Yes')"]
    op = v["op"]
    if op in ("ge", "gt", "le", "lt") and kind in (TEXT, BOOLEAN):
        return [f"'{v['field']}' holds {kind}, which can't be compared with '{op}'"]
    if op in ("in", "not_in", "starts_with_any", "not_starts_with_any") and "values" not in v:
        return [f"op '{op}' needs 'values'"]
    if op in ("ge", "gt", "le", "lt") and "value" not in v:
        return [f"op '{op}' needs 'value'"]
    if op != "known" and "missing" not in v:
        return ["say what a blank means: 'missing: pass' or 'missing: fail'"]
    if op == "known" and "missing" in v:
        return ["op 'known' already fails a blank; leave out 'missing'"]
    return []


condition = one_of(("check", record({"check": enum(*STRUCTURAL_CHECKS)})), ("field", _field_condition))
cohort_rule = record({"id": text(), "label": text(), "all_of": list_of(condition)})

result_match = record({"kind": enum(*RESULT_KINDS), "genes": list_of(text())},
                      {"variant_classes": list_of(text()), "codons": list_of(number(lo=1, integer=True)),
                       "codon_ranges": list_of(list_of(number(lo=1, integer=True), min_len=2)),
                       "variants": list_of(text()), "except_variants": list_of(text())})
genotype_group = one_of(("any_of", record({"name": text(), "any_of": list_of(result_match)})),
                        ("none_of_the_above", record({"name": text(), "none_of_the_above": enum(True)})))
codon_list = lambda v: [] if v == "whole_gene" else list_of(number(lo=1, integer=True))(v)
if_short = enum("describe", "block")

subset_keep = one_of(
    ("within_months", record({"within_months": number(lo=0, lo_open=True)})),
    ("field", _field_condition),
    ("exclude_overlap", record({"exclude_overlap": list_of(text(), min_len=2, unique=True)})),
    ("exclude_from_group", record({"exclude_from_group": text(), "if_result": result_match})),
)
ANALYSIS_KINDS: dict[str, Validator] = {
    "subset": record({"keep": subset_keep, "if_short": if_short}),
    "biomarker": record({
        "results_kind": enum(*RESULT_KINDS),
        "groups": list_of(record({"name": text(), "values": list_of(text())}), min_len=2),
        "unknown_group": text(), "reference": text(),
        "timing": enum("by_entry_test", "any_time"),
        "exclude_if": list_of(enum(*PATIENT_FLAGS), min_len=0),
        "if_short": if_short,
        "fallback": record({"when_panel_sizes_given": enum("mutation_burden", "descriptive_only"),
                            "otherwise": enum("descriptive_only")}),
        "mutation_burden": record({"threshold_per_megabase": number(lo=0, lo_open=True),
                                   "variant_classes": list_of(text())}),
        "not_run_under_fallback": list_of(text(), min_len=0)}),
    "biomarker_variant": record({"of": text()}, {"drop_discordant": boolean(), "timing": enum("by_entry_test", "any_time")}),
    "treatment_clock": record({"first_line": enum("earliest_by_order_on_or_after_diagnosis"), "if_short": if_short}),
    "adjusted_variant": record({"stratify_by": text()}),
    "ph_check": record({}),
    "quasi_independence": record({}),
    "descriptive_delay": record({"by": list_of(text())}),
}
COMMON_ANALYSIS = {"id": text(), "label": text(), "run": boolean(), "kind": enum(*ANALYSIS_KINDS)}


def analysis(v):
    if not isinstance(v, dict):
        return ["must be a mapping"]
    out = [f"missing key '{k}'" for k in COMMON_ANALYSIS if k not in v]
    for k, f in COMMON_ANALYSIS.items():
        if k in v:
            out += [f"{k}: {p}" for p in f(v[k])]
    if out:
        return out
    rest = {k: x for k, x in v.items() if k not in COMMON_ANALYSIS}
    return ANALYSIS_KINDS[v["kind"]](rest)


adjustment_variable = record({"field": text(), "form": enum("linear", "spline", "categorical", "bands")},
                             {"knots_at_percentiles": list_of(number(lo=0, hi=100)), "keep": list_of(text()),
                              "reference": text(), "merge": mapping_of(text(), list_of(text(), unique=True))})


def _adjustment_problems(v: dict) -> list[str]:
    """Checks inside one adjustment variable: the parts its form uses, and its levels."""
    form, out = v.get("form"), []
    name = v.get("field")
    if form != "spline" and "knots_at_percentiles" in v:
        out.append(f"'{name}': knots_at_percentiles is for form spline")
    if form == "spline":
        knots = v.get("knots_at_percentiles")
        if not isinstance(knots, list) or len(knots) < 3:
            out.append(f"'{name}': a spline needs knots_at_percentiles, at least three")
        elif any(b <= a for a, b in zip(knots, knots[1:])):
            out.append(f"'{name}': the knots' percentiles must rise")
    if form != "categorical":
        out += [f"'{name}': {k} is for form categorical" for k in ("keep", "reference", "merge") if k in v]
        return out
    merge = v.get("merge") or {}
    seen: dict[str, str] = {}
    for new, olds in merge.items():
        for old in olds if isinstance(olds, list) else []:
            if old in seen:
                out.append(f"'{name}': the value '{old}' is merged into both '{seen[old]}' and '{new}'")
            seen[old] = new
    levels = list(merge) if merge else None
    keep = v.get("keep")
    if keep is not None and levels is not None:
        out += [f"'{name}': keep names '{k}', which is not one of the merged levels ({', '.join(levels)})"
                for k in keep if k not in levels]
    levels = keep if keep is not None else levels
    ref = v.get("reference")
    if ref is not None and levels is not None and ref not in levels:
        out.append(f"'{name}': the reference '{ref}' is not one of its levels ({', '.join(levels)})")
    return out
year_band = record({"label": text()}, {"from": number(integer=True), "to": number(integer=True)})

IMPUTATION_TERMS = ("event", "cumulative_hazard_corrected", "cumulative_hazard_naive", "entry_time", "exit_time")
imputation_model = record({
    "engine": enum("chained_equations", "mice"),
    "outcome_terms": list_of(enum(*IMPUTATION_TERMS), min_len=0, unique=True),
    "include_group": boolean(),
    "methods": record({"linear": enum("pmm", "norm"), "spline": enum("pmm", "norm"),
                       "categorical": enum("polyreg", "logreg"), "bands": enum("polyreg", "polr")}),
    "pmm_donors": number(lo=1, integer=True),
    "iterations": number(lo=1, integer=True),
    "pooling": enum("rubins_rules"),
    "bootstrap": enum("impute_each_resample_once"),
    "point_estimates": enum("average_across_imputations"),
})


# ---------------------------------------------------------------------------
# The schema: every choice, what form it takes and what it controls.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Choice:
    path: str
    check: Validator
    help: str
    form: str                          # the form, in words, for the template
    when: tuple[str, str] | None = None  # only for this value of another choice
    example: Any = None                # shown in the template for complex choices
    optional: bool = False             # may be left out: whatever needs it then stops and asks for it (no default)


SCHEMA: list[Choice] = [
    Choice("study.name", text(), "A short name for the study.", "text"),
    Choice("study.cancer", text(), "The cancer and the population, in words.", "text"),
    Choice("study.registered_plan", text(), "Where the plan is registered, or 'not registered'.", "text"),

    Choice("source.adapter", enum("genie_bpc", "flat_file", "omop"), "Which adapter reads the data.",
           "genie_bpc, flat_file or omop"),
    Choice("source.genie_bpc.cohort", text(), "Which GENIE BPC cohort and release, for example crc_v2.0.",
           "text", when=("source.adapter", "genie_bpc")),
    Choice("source.genie_bpc.cancer_type_labels", list_of(text()),
           "The source's labels for this cancer: which pathology specimens and which other primaries are of this kind.",
           "a list of labels", when=("source.adapter", "genie_bpc"), example=["Prostate Cancer"]),
    Choice("source.genie_bpc.cancer_site_pattern", text(),
           "A regular expression for this cancer's ICD-O-3 site codes (checked after download).",
           "a regular expression", when=("source.adapter", "genie_bpc")),
    Choice("source.genie_bpc.withheld_survival", enum("censor_at_age_limit", "leave_blank"),
           "When the source withholds survival to hide an age over its limit: censor alive at that age, or leave "
           "the outcome blank (the patient then fails 'outcome_recorded').",
           "censor_at_age_limit or leave_blank", when=("source.adapter", "genie_bpc")),

    Choice("cohort.rules", list_of(cohort_rule),
           "Who is analyzed, in the order the flow diagram counts them. Each rule has an id, a label and the "
           "conditions it needs. Every structural check (" + ", ".join(STRUCTURAL_CHECKS) + ") must appear once.",
           "a list of rules",
           example=[{"id": "adults", "label": "Adults at diagnosis",
                     "all_of": [{"field": "age_at_diagnosis", "op": "ge", "value": 18, "missing": "fail"}]},
                    {"id": "has_test", "label": "Has a sequencing test", "all_of": [{"check": "has_entry_test"}]}]),

    Choice("entry.event", enum(*ENTRY_EVENTS), "The event where the corrected clock starts.", " or ".join(ENTRY_EVENTS)),
    Choice("entry.acceptable_substitutes", list_of(enum(*ENTRY_EVENTS), min_len=0, unique=True),
           "Events that may stand in when the source lacks the entry event; using one is a deviation. "
           "An empty list means none may.", "a list, possibly empty"),
    Choice("entry.modalities", list_of(text(), unique=True),
           "Which kinds of test can be the entry test (the tests table's modality), for example ngs.", "a list"),
    Choice("entry.first_test_order", list_of(enum("earliest_day", "largest_panel", "lowest_test_order"), unique=True),
           "How the first test is chosen: the order of tie-breaks. A test whose time is withheld always ranks last.",
           "a list from: earliest_day, largest_panel, lowest_test_order"),
    Choice("entry.before_diagnosis", enum("enter_at_day_0", "exclude"),
           "What happens when the entry event comes before diagnosis.", "enter_at_day_0 or exclude"),

    Choice("genotype.results_from", enum("entry_test"), "Which test's results define the groups.", "entry_test"),
    Choice("genotype.variant_text", enum("trimmed", "exact"),
           "Whether spaces around a variant's text are ignored when matching it.", "trimmed or exact"),
    Choice("genotype.groups", list_of(genotype_group, min_len=2),
           "The groups, in order of precedence: a patient matching two goes to the first. The last is the "
           "group of everyone else ('none_of_the_above: true').", "a list of groups",
           example=[{"name": "Mutant", "any_of": [{"kind": "mutation", "genes": ["GENE1"], "codons": [12]}]},
                    {"name": "No listed variant", "none_of_the_above": True}]),
    Choice("genotype.reference", text(), "The reference group for every comparison.", "a group's name"),
    Choice("genotype.described_only", list_of(text(), min_len=0, unique=True),
           "Groups that are described, never compared, whatever their size: for example a group gathering several "
           "rare drivers. The reference can't be one. An empty list means every group is compared when it is large "
           "enough.", "a list of groups' names, possibly empty", example=["Rare drivers"]),
    Choice("genotype.coverage.required", mapping_of(text(), codon_list),
           "What a panel must read for a patient's group to be trusted: each gene, with its codons or 'whole_gene'.",
           "gene: list of codons", example={"GENE1": [12, 13]}),
    Choice("genotype.coverage.applies_to", enum("all_groups", "reference_only"),
           "Whether the coverage rule applies to every group or only to the reference (wild type).",
           "all_groups or reference_only"),
    Choice("genotype.coverage.when_only_genes_listed", enum("accept_gene_level", "deviation"),
           "When the panels list genes but not the codons they read, or a group counts a required gene's fusions "
           "and fusion_panels names no panels for it: accept a listed gene as covering its codons (or as detecting "
           "its fusions), or treat it as a deviation.", "accept_gene_level or deviation"),
    Choice("genotype.coverage.fusion_panels", mapping_of(text(), list_of(text(), unique=True), min_len=0),
           "For a required gene whose fusions a group counts: the panels that detect its fusions, from the panels' "
           "documentation. A gene list shows that a panel reads a gene, not that it detects the gene's fusions, so "
           "for a gene named here a panel covers it only if it is on this list. An empty mapping names none.",
           "gene: list of panels; possibly empty", example={"GENE2": ["PANEL-A", "PANEL-B"]}),
    Choice("genotype.reference_drugs", mapping_of(text(), list_of(text(), unique=True), min_len=0),
           "Drugs aimed at a driver, listed by the driver they target, as the source names them. The screen counts, "
           "masked, the reference group's patients who received one at any time: a sign of a driver the first test "
           "didn't find. It only counts; nobody changes group. An empty mapping means no such count.",
           "driver: list of drugs; possibly empty", example={"GENE2 fusion": ["Drug A", "Drug B"]}),

    Choice("time.days_per_year", number(lo=0, lo_open=True), "Days in a year (a month is a twelfth of it).", "a number"),
    Choice("time.horizon_years", number(lo=0, lo_open=True), "The horizon for months alive (RMST).", "years"),
    Choice("summary.survival_at_years", list_of(number(lo=0, lo_open=True), unique=True),
           "The fixed times for survival estimates.", "a list of years"),
    Choice("summary.median_limit_years", number(lo=0, lo_open=True),
           "A median, or a limit of its interval, after this is reported as not reached.", "years"),

    Choice("materiality.rmst_months", number(lo=0, lo_open=True), "The bar for a difference in RMST.", "months"),
    Choice("materiality.median_months", number(lo=0, lo_open=True), "The bar for a difference in median.", "months"),
    Choice("materiality.hr_ratio_change", number(lo=0, lo_open=True),
           "The bar for |HR naive / HR corrected - 1|.", "a share, for example 0.10"),
    Choice("materiality.rmst_gap_change_months", number(lo=0, lo_open=True),
           "The bar for the change in a between-group RMST difference.", "months"),
    Choice("labels.clearly_material", enum("estimate_and_whole_interval"),
           "When a result is labeled clearly material.", "estimate_and_whole_interval"),
    Choice("labels.material", enum("estimate_only"), "When a result is labeled material.", "estimate_only"),

    Choice("smallest_group.patients", number(lo=1, integer=True),
           "The smallest group compared: a smaller one is described, not compared.", "patients"),
    Choice("smallest_group.if_short", mapping_of(text(), if_short),
           "For each genotype group, what happens if it is short: 'describe' it, or 'block' the analysis "
           "(a person then decides, and logs a deviation).", "group: describe or block",
           example={"Mutant": "block", "No listed variant": "block"}),
    Choice("biomarker.unknown_max_share", number(lo=0, hi=1),
           "Above this share of Unknown, a biomarker analysis uses its fallback.", "a share, for example 0.30"),

    Choice("survival_after_entry.window_days", number(lo=0, lo_open=True),
           "The survival-after-sequencing check: deaths within this many days of the entry event.", "days"),
    Choice("survival_after_entry.min_deaths", number(lo=1, integer=True),
           "Fewer deaths than this means the source in effect required survival after the entry event.", "deaths"),
    Choice("survival_after_entry.if_fewer", enum("move_entry_by_window", "warn_only"),
           "What happens then: entry moves later by the window, or a warning only.", "move_entry_by_window or warn_only"),
    Choice("stability.risk_set", number(lo=1, integer=True),
           "The corrected risk set a curve needs before it is trusted (the t* rule).", "patients"),

    Choice("uncertainty.bootstrap_resamples", number(lo=1, integer=True), "Bootstrap resamples.", "a whole number"),
    Choice("uncertainty.seed", number(integer=True), "The random seed for the resamples.", "a whole number"),
    Choice("uncertainty.interval_level", number(lo=0, hi=1, lo_open=True), "The interval's coverage.", "for example 0.95"),
    Choice("uncertainty.when_a_resample_has_no_value", enum("no_interval"),
           "When a resample can't give a value (a median not reached): no interval, and the count is reported.",
           "no_interval"),

    Choice("adjustment.variables", list_of(adjustment_variable),
           "The adjustment variables (patients-table fields, or extras the adapter supplies) and their form. A "
           "categorical variable may list the levels it keeps ('keep': any other value is missing), name its "
           "reference level ('reference': otherwise the first level), and merge the source's values into new "
           "levels ('merge': each new level with the values it takes in; a value in no group is missing, and "
           "'keep' then names new levels).", "a list",
           example=[{"field": "age_at_diagnosis", "form": "spline", "knots_at_percentiles": [10, 50, 90]},
                    {"field": "smoking", "form": "categorical", "reference": "Never",
                     "merge": {"Never": ["Never"], "Former": ["Former, quit over a year ago"],
                               "Current": ["Current", "Former, quit within a year"]}}]),
    Choice("adjustment.year_bands", list_of(year_band, min_len=0),
           "Bands for year of diagnosis, if it is used in bands.", "a list of bands, possibly empty"),
    Choice("adjustment.band_min_patients", number(lo=1, integer=True),
           "A smaller year band is merged before modeling.", "patients"),
    Choice("missing_data.complete_case_max_share", number(lo=0, hi=1),
           "Complete-case analysis if at most this share has a missing covariate; otherwise multiple imputation.",
           "a share"),
    Choice("missing_data.imputations", number(lo=1, integer=True), "Imputations, when imputing.", "a whole number"),
    Choice("missing_data.imputation_model", imputation_model,
           "How multiple imputation is done, when the complete-case rule calls for it (adjustment variables only; the "
           "group and the times are never imputed), by chained equations: the engine (the pipeline's own, in R and "
           "nnet, or the mice package); the outcome terms the imputation model holds besides the adjustment variables "
           "(the event indicator, the Nelson-Aalen cumulative hazard on the corrected or the naive clock, the entry or "
           "exit time); whether it holds the group; the method for each form of variable (pmm: predictive mean "
           "matching; norm: Bayesian linear regression; logreg: logistic, two levels only; polyreg: multinomial; "
           "polr: ordered, in the bands' order); the donors predictive mean matching draws from; and the iterations. "
           "Hazard ratios are pooled by Rubin's rules; for bootstrap intervals, each resample is imputed once; point "
           "estimates are averages across the imputations. The random seed is uncertainty.seed. May be left out if "
           "imputation is never needed: a run that needs it then stops.",
           "a mapping", optional=True,
           example={"engine": "chained_equations",
                    "outcome_terms": ["event", "cumulative_hazard_corrected", "entry_time"], "include_group": True,
                    "methods": {"linear": "pmm", "spline": "pmm", "categorical": "polyreg", "bands": "polyreg"},
                    "pmm_donors": 5, "iterations": 5, "pooling": "rubins_rules",
                    "bootstrap": "impute_each_resample_once", "point_estimates": "average_across_imputations"}),
    Choice("weighting.variables", enum("same_as_adjustment"), "The traits used for weighting.", "same_as_adjustment"),
    Choice("weighting.trim_percentiles", list_of(number(lo=0, hi=100), min_len=2),
           "Weights are trimmed at these percentiles.", "two percentiles"),
    Choice("weighting.balance_target_smd", number(lo=0, lo_open=True),
           "Balance target: absolute standardized mean differences below this.", "a number"),
    Choice("ph_check.p_threshold", number(lo=0, hi=1, lo_open=True),
           "Proportional-hazards check: below this p, hazard ratios are split.", "a p-value"),
    Choice("ph_check.split_months", number(lo=0, lo_open=True), "Where hazard ratios are split.", "months"),
    Choice("assumption_check.estimators", list_of(enum("MB", "IPW1", "IPW2"), unique=True),
           "Estimators of the conditional Kendall's tau for the delayed-entry assumption.", "a list"),

    Choice("analyses", list_of(analysis, min_len=0),
           "The secondary and sensitivity analyses, each switched on or off, with its parameters.", "a list",
           example=[{"id": "SA1", "label": "Short waits only", "run": True, "kind": "subset",
                     "keep": {"within_months": 6}, "if_short": "describe"}]),
    Choice("privacy.min_count", number(lo=1, integer=True),
           "The smallest count shown. It can be raised, never set below the data-use agreement's.", "a whole number"),
]
SCHEMA_BY_PATH = {c.path: c for c in SCHEMA}


def load_defaults(path: Path = DEFAULTS_FILE) -> dict[str, dict]:
    """{choice: {"value": ..., "source": ...}} for the choices that have a documented default."""
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    out = {}
    for choice, entry in raw.items():
        if choice not in SCHEMA_BY_PATH:
            raise SettingsError(f"{path.name}: '{choice}' is not a choice in the schema")
        if not isinstance(entry, dict) or set(entry) != {"value", "source"}:
            raise SettingsError(f"{path.name}: '{choice}' needs exactly 'value' and 'source'")
        problems = SCHEMA_BY_PATH[choice].check(entry["value"])
        if problems:
            raise SettingsError(f"{path.name}: '{choice}': {'; '.join(problems)}")
        out[choice] = entry
    return out


# ---------------------------------------------------------------------------
# Reading the file, keeping line numbers
# ---------------------------------------------------------------------------
def _lines(node, path=(), out=None) -> dict[tuple, int]:
    """The line (1-based) of every key and list item, by its path."""
    out = {} if out is None else out
    if isinstance(node, yaml.MappingNode):
        for k, v in node.value:
            p = path + (k.value,)
            out[p] = k.start_mark.line + 1
            _lines(v, p, out)
    elif isinstance(node, yaml.SequenceNode):
        for i, v in enumerate(node.value):
            out[path + (i,)] = v.start_mark.line + 1
            _lines(v, path + (i,), out)
    return out


@dataclass(frozen=True)
class Resolved:
    choice: str
    value: Any
    made_by: str
    source: str
    origin: str          # "settings file, line N" or "default (defaults.yml)"
    note: str = ""

    @property
    def is_default(self) -> bool:
        return self.origin.startswith("default")


class Settings:
    """A loaded, checked settings file. Read a choice with read(), which records that it was read, so
    a report can list every choice it used and where each came from."""

    def __init__(self, path: Path | None, text_: str, resolved: dict[str, Resolved]):
        self.path = path
        self.text = text_
        # For a file, the text is its bytes decoded as they are (load_settings), so this is the file's own hash.
        self.sha256 = hashlib.sha256(text_.encode("utf-8")).hexdigest()
        self._resolved = resolved
        self._reads: list[str] = []
        # The values in force, defaults included: a changed defaults file changes this, not the file's own hash.
        self.resolved_sha256 = hashlib.sha256(
            json.dumps(self.provenance(), sort_keys=True, default=str).encode("utf-8")).hexdigest()

    def read(self, choice: str) -> Any:
        if choice not in self._resolved:
            raise KeyError(f"no choice '{choice}' in these settings")
        if choice not in self._reads:
            self._reads.append(choice)
        return copy.deepcopy(self._resolved[choice].value)

    def has(self, choice: str) -> bool:
        """Whether the settings hold this choice: always, except an optional choice left out."""
        return choice in self._resolved

    def resolved(self, choice: str) -> Resolved:
        return self._resolved[choice]

    def reads(self) -> list[str]:
        return list(self._reads)

    def forget_reads(self) -> None:
        self._reads.clear()

    def provenance(self, choices: list[str] | None = None) -> list[dict]:
        """Every choice (or the ones named): its value, who made it, the source and where it came from."""
        names = choices if choices is not None else list(self._resolved)
        return [{"choice": c, "value": self._resolved[c].value, "made_by": self._resolved[c].made_by,
                 "source": self._resolved[c].source, "origin": self._resolved[c].origin,
                 **({"note": self._resolved[c].note} if self._resolved[c].note else {})} for c in names]

    def defaults_used(self, choices: list[str] | None = None) -> list[dict]:
        return [p for p in self.provenance(choices) if p["origin"].startswith("default")]

    def analyses(self) -> list[dict]:
        return self.read("analyses")

    def analysis(self, analysis_id: str) -> dict | None:
        return next((a for a in self.read("analyses") if a["id"] == analysis_id), None)


def _get(data, path: str):
    node = data
    for part in path.split("."):
        if not isinstance(node, dict) or part not in node:
            return _MISSING
        node = node[part]
    return node


_MISSING = object()


def _applies(choice: Choice, data) -> bool | None:
    """True or False, or None while the choice it depends on is still blank."""
    if choice.when is None:
        return True
    other = _get(data, choice.when[0])
    if not isinstance(other, dict) or other.get("value") is None:
        return None
    return other.get("value") == choice.when[1]


def _known_prefixes() -> set[tuple]:
    out = set()
    for c in SCHEMA:
        parts = tuple(c.path.split("."))
        for i in range(1, len(parts) + 1):
            out.add(parts[:i])
    return out


def _cross_checks(values: dict[str, Any]) -> list[tuple[str, str]]:
    """Checks between choices: (choice, problem)."""
    out = []
    groups = values.get("genotype.groups")
    if isinstance(groups, list) and all(isinstance(g, dict) and "name" in g for g in groups):
        names = [g["name"] for g in groups]
        if len(set(names)) != len(names):
            out.append(("genotype.groups", "two groups share a name"))
        rest = [i for i, g in enumerate(groups) if g.get("none_of_the_above")]
        if rest != [len(groups) - 1]:
            out.append(("genotype.groups", "exactly one group, the last, must be 'none_of_the_above: true'"))
        ref = values.get("genotype.reference")
        if isinstance(ref, str) and ref not in names:
            out.append(("genotype.reference", f"'{ref}' is not one of the groups ({', '.join(names)})"))
        described = values.get("genotype.described_only")
        described = described if isinstance(described, list) else []
        for g in described:
            if g not in names:
                out.append(("genotype.described_only", f"'{g}' is not one of the groups ({', '.join(names)})"))
        if isinstance(ref, str) and ref in described:
            out.append(("genotype.described_only", f"'{ref}' is the reference, which every comparison needs"))
        compared = [g for g in names if g not in described]
        if described and len(compared) < 2:
            out.append(("genotype.described_only", "at least one group besides the reference must be compared"))
        short = values.get("smallest_group.if_short")
        if isinstance(short, dict) and set(short) != set(compared):
            out.append(("smallest_group.if_short", f"must name each group that is compared exactly once: "
                                                   f"{', '.join(compared)}"
                        + (" (a group described only has no rule here)" if described else "")))
        with_rules = [g["name"] for g in groups if not g.get("none_of_the_above")]
        analyses = values.get("analyses") or []
        for a in analyses if isinstance(analyses, list) else []:
            keep = a.get("keep", {}) if isinstance(a, dict) else {}
            for g in keep.get("exclude_overlap", []) if isinstance(keep, dict) else []:
                if g not in with_rules:
                    out.append(("analyses", f"{a.get('id')}: '{g}' is not a genotype group with its own rules "
                                            f"({', '.join(with_rules)})"))
            g = keep.get("exclude_from_group") if isinstance(keep, dict) else None
            if g is not None and g not in names:
                out.append(("analyses", f"{a.get('id')}: '{g}' is not a genotype group"))
    required, fusion_panels = values.get("genotype.coverage.required"), values.get("genotype.coverage.fusion_panels")
    if isinstance(required, dict) and isinstance(fusion_panels, dict):
        for gene in fusion_panels:
            if gene not in required:
                out.append(("genotype.coverage.fusion_panels",
                            f"'{gene}' is not a required gene ({', '.join(required)}), so its list would do nothing"))
    rules = values.get("cohort.rules")
    if isinstance(rules, list):
        ids = [r.get("id") for r in rules if isinstance(r, dict)]
        if len(set(ids)) != len(ids):
            out.append(("cohort.rules", "two rules share an id"))
        checks = [c.get("check") for r in rules if isinstance(r, dict) for c in r.get("all_of", []) if isinstance(c, dict)]
        for name in STRUCTURAL_CHECKS:
            n = checks.count(name)
            if n != 1:
                out.append(("cohort.rules", f"the structural check '{name}' must appear exactly once (found {n})"))
    analyses = values.get("analyses")
    if isinstance(analyses, list) and all(isinstance(a, dict) for a in analyses):
        ids = [a.get("id") for a in analyses]
        if len(set(ids)) != len(ids):
            out.append(("analyses", "two analyses share an id"))
        kinds = {a.get("id"): a.get("kind") for a in analyses}
        runs = {a.get("id"): a.get("run") for a in analyses}
        parent = {a.get("id"): a.get("of") for a in analyses if a.get("kind") == "biomarker_variant"}
        for a in analyses:
            if a.get("kind") == "biomarker_variant" and kinds.get(a.get("of")) != "biomarker":
                out.append(("analyses", f"{a.get('id')}: 'of' must name a biomarker analysis"))
            elif a.get("kind") == "biomarker_variant" and a.get("run") and runs.get(a.get("of")) is False:
                out.append(("analyses", f"{a.get('id')} runs, but {a.get('of')}, which it varies, is switched off"))
            if a.get("kind") == "biomarker":
                for other in a.get("not_run_under_fallback", []):
                    if parent.get(other) != a.get("id"):
                        out.append(("analyses", f"{a.get('id')}: '{other}' is not one of its own variants "
                                                "(kind biomarker_variant, of: " + str(a.get("id")) + ")"))
                groups = a.get("groups")
                if isinstance(groups, list):
                    names = [g.get("name") for g in groups if isinstance(g, dict)]
                    if a.get("reference") not in names:
                        out.append(("analyses", f"{a.get('id')}: the reference must be one of its groups "
                                                f"({', '.join(str(n) for n in names)})"))
            if a.get("kind") == "adjusted_variant" and a.get("stratify_by") in OUTCOME_FIELDS:
                out.append(("analyses", f"{a.get('id')}: an outcome can't be a stratum"))
    variables = values.get("adjustment.variables") or []
    for v in variables:
        if isinstance(v, dict) and v.get("field") in OUTCOME_FIELDS:
            out.append(("adjustment.variables", f"'{v['field']}' is an outcome, not an adjustment variable"))
        if isinstance(v, dict):
            out += [("adjustment.variables", p) for p in _adjustment_problems(v)]
    fields = [v.get("field") for v in variables if isinstance(v, dict)]
    if len(set(fields)) != len(fields):
        out.append(("adjustment.variables", "a field is listed twice"))
    model = values.get("missing_data.imputation_model")
    if isinstance(model, dict) and isinstance(model.get("methods"), dict) and model["methods"].get("categorical") == "logreg":
        for v in variables:
            if isinstance(v, dict) and v.get("form") == "categorical":
                levels = v.get("keep") or list(v.get("merge") or {}) or None
                if levels is not None and len(levels) != 2:
                    out.append(("missing_data.imputation_model", f"logreg imputes a variable with two levels, and "
                                                                 f"'{v.get('field')}' has {len(levels)}: use polyreg"))
    if any(isinstance(v, dict) and v.get("form") == "bands" for v in variables) and not values.get("adjustment.year_bands"):
        out.append(("adjustment.year_bands", "a variable is used in bands, so the bands must be listed"))
    trim = values.get("weighting.trim_percentiles")
    if isinstance(trim, list) and len(trim) == 2 and all(isinstance(x, (int, float)) for x in trim) and trim[0] >= trim[1]:
        out.append(("weighting.trim_percentiles", "the lower percentile must come first"))
    return out


def loads_settings(text_: str, path: Path | None = None, defaults: dict | None = None) -> Settings:
    """Check settings given as text. Raises SettingsError listing every problem, each with its line."""
    where = path.name if path else "settings"
    try:
        data = yaml.load(text_, Loader=_UniqueKeyLoader)
        lines = _lines(yaml.compose(text_, Loader=yaml.SafeLoader))
    except yaml.YAMLError as e:
        mark = getattr(e, "problem_mark", None)
        problem = getattr(e, "problem", None) or "not valid YAML"
        raise SettingsError(f"{where}, line {mark.line + 1 if mark else '?'}: {problem}") from None
    if not isinstance(data, dict):
        raise SettingsError(f"{where}: must be a mapping of choices")
    defaults = load_defaults() if defaults is None else defaults

    def line(p) -> str:
        return f"line {lines.get(tuple(p.split('.')) if isinstance(p, str) else p, '?')}"

    problems = []
    if data.get("settings_format") != SETTINGS_FORMAT:
        problems.append(f"{where}, {line('settings_format')}: settings_format must be {SETTINGS_FORMAT}")
    known = _known_prefixes() | {("settings_format",)}
    for keys in lines:
        if any(isinstance(k, int) for k in keys) or keys in known:
            continue        # a list belongs to a value; a known key is part of a choice's path
        at = next((i for i in range(1, len(keys)) if ".".join(keys[:i]) in SCHEMA_BY_PATH), None)
        if at is not None and keys[at] in ENTRY_KEYS and (len(keys) == at + 1 or keys[at] == "value"):
            continue        # value, made_by, source or note under a choice, or a key inside a value
        problems.append(f"{where}, line {lines[keys]}: '{'.'.join(map(str, keys))}' is not a known choice "
                        "or entry key (check the spelling)")

    resolved, values = {}, {}
    for c in SCHEMA:
        applies = _applies(c, data)
        entry = _get(data, c.path)
        if applies is False:
            if entry is not _MISSING:
                problems.append(f"{where}, {line(c.path)}: '{c.path}' applies only when {c.when[0]} is {c.when[1]}")
            continue
        if applies is None and entry is _MISSING:
            continue
        if entry is _MISSING:
            if not c.optional:
                problems.append(f"{where}: '{c.path}' is missing. Every choice must be written out ({c.help})")
            continue
        if not isinstance(entry, dict) or "value" not in entry:
            problems.append(f"{where}, {line(c.path)}: '{c.path}' needs 'value', 'made_by' and 'source'")
            continue
        value, made_by, source = entry.get("value"), entry.get("made_by"), entry.get("source")
        if value is None:
            problems.append(f"{where}, {line(c.path + '.value')}: '{c.path}' is not filled in ({c.help})")
            continue
        if not (isinstance(made_by, str) and made_by.strip()):
            problems.append(f"{where}, {line(c.path)}: '{c.path}' needs 'made_by': who made this choice")
        if value == DEFAULT_WORD:
            if c.path not in defaults:
                problems.append(f"{where}, {line(c.path + '.value')}: '{c.path}' has no default; fill in a value")
                continue
            value = defaults[c.path]["value"]
            origin = f"default (defaults.yml), chosen at {line(c.path + '.value')}"
            source = defaults[c.path]["source"] + (f"; {source}" if isinstance(source, str) and source.strip() else "")
        else:
            if not (isinstance(source, str) and source.strip()):
                problems.append(f"{where}, {line(c.path)}: '{c.path}' needs 'source': the source behind the choice")
            origin = f"settings file, {line(c.path + '.value')}"
            bad = c.check(value)
            if bad:
                for p in bad:       # a problem inside a list points at that item's line
                    item = re.match(r"\[(\d+)\]", p)
                    at = (tuple(c.path.split(".")) + ("value", int(item.group(1)))) if item else c.path + ".value"
                    problems.append(f"{where}, {line(at)}: '{c.path}' {p}".rstrip())
                continue
        values[c.path] = value
        resolved[c.path] = Resolved(c.path, value, str(made_by or ""), str(source or ""), origin, str(entry.get("note") or ""))
    for choice, problem in _cross_checks(values):
        problems.append(f"{where}, {line(choice)}: '{choice}': {problem}")
    if problems:
        raise SettingsError(f"{len(problems)} problem(s) in the settings:\n  " + "\n  ".join(problems))
    return Settings(path, text_, resolved)


def load_settings(path: str | Path, defaults: dict | None = None) -> Settings:
    """Load a settings file. It is decoded as it is, line endings included (read_text would turn Windows line
    endings into "\\n"), so the settings' hash is the hash of the file's own bytes: the one the lock's manifest
    and the settings lock compare the file with."""
    path = Path(path)
    return loads_settings(path.read_bytes().decode("utf-8"), path, defaults)


# ---------------------------------------------------------------------------
# The blank template
# ---------------------------------------------------------------------------
def blank_template(defaults: dict | None = None) -> str:
    """The settings template with every choice blank. Each choice says what it controls, its form, an
    example where the form is complex, and its default if it has one (write `value: default` to use it)."""
    return settings_text({}, defaults)


def applies_to(choice: Choice, answers: dict[str, dict]) -> bool | None:
    """Whether a choice applies, given answers ({path: {value, made_by, source}}): True, False, or None while the
    choice it depends on is unanswered."""
    if choice.when is None:
        return True
    other = answers.get(choice.when[0], {}).get("value")
    return None if other is None else other == choice.when[1]


def settings_text(answers: dict[str, dict], defaults: dict | None = None) -> str:
    """The template with these answers filled in ({path: {value, made_by, source, note}}): the blank template when
    there are none. With answers, a choice that doesn't apply is left out, and so is an optional choice left
    unanswered; every other unanswered choice stays blank, for the loader to name."""
    defaults = load_defaults() if defaults is None else defaults
    out = ["# Settings template. Every choice must be filled in: value, made_by (who made the choice) and",
           "# source (the source behind it). A blank, a missing choice or an unknown key stops the tool and",
           "# names the line. Where a choice has a default, writing `value: default` uses it; the locked plan",
           "# lists every default used and its source. Choices marked 'no default' must be filled in.",
           "",
           f"settings_format: {SETTINGS_FORMAT}", ""]
    shown = SCHEMA if not answers else [c for c in SCHEMA if applies_to(c, answers) is not False
                                        and (c.path in answers or not c.optional)]
    open_sections: list[str] = []
    for c in shown:
        parts = c.path.split(".")
        common = 0
        while common < min(len(open_sections), len(parts) - 1) and open_sections[common] == parts[common]:
            common += 1
        if common == 0 and out[-1] != "":
            out.append("")
        for depth in range(common, len(parts) - 1):
            out.append("  " * depth + f"{parts[depth]}:")
        open_sections = parts[:-1]
        indent = "  " * (len(parts) - 1)
        out.append(f"{indent}# {c.help}")
        out.append(f"{indent}# Form: {c.form}." + (f" Only when {c.when[0]} is {c.when[1]}." if c.when else ""))
        if c.path in defaults:
            d = defaults[c.path]
            shown = yaml.safe_dump(d["value"], default_flow_style=True, width=10**6).replace("\n...", "").strip()
            out.append(f"{indent}# Default: {shown} (source: {d['source']})")
        elif c.optional:
            out.append(f"{indent}# No default. Optional: fill it in, or delete this choice if it is never needed.")
        else:
            out.append(f"{indent}# No default: fill this in.")
        if c.example is not None:
            ex = yaml.safe_dump({"value": c.example}, default_flow_style=False, sort_keys=False).rstrip().splitlines()
            out.append(f"{indent}# Example:")
            out += [f"{indent}#   {x}" for x in ex]
        out.append(f"{indent}{parts[-1]}:")
        a = answers.get(c.path)
        if a is None:
            out += [f"{indent}  value: null", f"{indent}  made_by: null", f"{indent}  source: null"]
        else:
            entry = {k: a.get(k) for k in ("value", "made_by", "source")}
            if a.get("note"):
                entry["note"] = a["note"]
            text_ = yaml.safe_dump(entry, sort_keys=False, allow_unicode=True, default_flow_style=False, width=110)
            out += [f"{indent}  {x}" for x in text_.rstrip().splitlines()]
    return "\n".join(out) + "\n"


def write_blank_template(path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(blank_template().encode("utf-8"))
    return path


def unfilled_choices(error: SettingsError) -> list[str]:
    """The choices an error says are not filled in."""
    return re.findall(r"'([\w.]+)' is not filled in", str(error))
