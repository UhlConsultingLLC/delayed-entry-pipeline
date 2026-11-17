"""The cohort, the first test and the genotype groups, built from the outcome-blind tables and the settings.

Generalized from the CRC BPC project's cohort SQL (sql/02_cohort_rules.sql there). Every rule, group,
codon and tie-break comes from the settings; this module only knows how to apply them.

    from delayed_entry.cohort import build_cohort
    c = build_cohort(outcome_blind(tables), settings)
    c.patients        # one row per patient: each rule, where they leave, the first test, the group
    c.flow()          # exact counts at each step: for checking only, never published as is

Like the CRC code, a rule never leaves a patient undecided: each condition says what a blank means.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from delayed_entry.settings import Settings
from delayed_entry.tables import OUTCOME_RECORDED, OutcomeBlindTables

# Codons a protein change spans, in HGVSp short form: "p.G12D" spans 12 to 12, "p.G12_G13del" 12 to 13.
# A change the pattern can't read gives no codon, so it never matches a codon rule.
CODON_FROM = re.compile(r"^p\.[A-Z*](\d+)")
CODON_END = re.compile(r"^p\.[A-Z*]\d+_[A-Z*](\d+)")


class MissingFieldError(KeyError):
    """A rule names a field the tables don't have."""


# ---------------------------------------------------------------------------
# Conditions on patient fields
# ---------------------------------------------------------------------------
def normalize(s: pd.Series, how: str) -> pd.Series:
    if how in (None, "none") or not (pd.api.types.is_string_dtype(s) or s.dtype == object):
        return s
    t = s.astype("string")
    if how == "trim":
        return t.str.strip()
    if how == "letters_digits_upper":
        return t.str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
    raise ValueError(f"unknown normalization {how}")


def field_condition(df: pd.DataFrame, cond: dict) -> pd.Series:
    """True where the condition passes, as a plain boolean Series (never blank)."""
    name = cond["field"]
    if name not in df.columns:
        raise MissingFieldError(f"patients.{name}")
    s = normalize(df[name], cond.get("normalize", "none"))
    op = cond["op"]
    blank = s.isna()
    if op == "known":
        return ~blank
    numeric = pd.api.types.is_numeric_dtype(s)
    if op in ("in", "not_in"):
        # on a field of numbers, a value that isn't a number matches nothing
        values = ([float(v) for v in cond["values"] if isinstance(v, (int, float)) and not isinstance(v, bool)]
                  if numeric else [str(v) for v in cond["values"]])
        hit = s.isin(values)
        res = hit if op == "in" else ~hit
    elif op in ("starts_with_any", "not_starts_with_any"):
        text_ = s.astype("string")
        hit = pd.Series(False, index=s.index)
        for v in cond["values"]:
            hit |= text_.str.startswith(str(v)).fillna(False).astype(bool)
        res = hit if op == "starts_with_any" else ~hit
    else:
        num = pd.to_numeric(s, errors="coerce")
        blank = num.isna()
        v = float(cond["value"])
        res = {"ge": num >= v, "gt": num > v, "le": num <= v, "lt": num < v}[op]
    res = pd.Series(res, index=s.index).astype("boolean").fillna(False).astype(bool)
    res[blank] = cond["missing"] == "pass"
    return res


# ---------------------------------------------------------------------------
# Results: codons and matching
# ---------------------------------------------------------------------------
def variant_text(results: pd.DataFrame, how: str) -> pd.Series:
    v = results["variant"].astype("string") if "variant" in results.columns else pd.Series(pd.NA, index=results.index, dtype="string")
    return v.str.strip() if how == "trimmed" else v


def codon_span(variants: pd.Series) -> tuple[pd.Series, pd.Series]:
    """(first codon, last codon) each change spans; a one-codon change ends where it starts."""
    start = pd.to_numeric(variants.str.extract(CODON_FROM, expand=False), errors="coerce")
    end = pd.to_numeric(variants.str.extract(CODON_END, expand=False), errors="coerce")
    return start, end.fillna(start)


def result_matches(results: pd.DataFrame, match: dict, how: str) -> pd.Series:
    """True for each result row that matches a rule such as {kind, genes, variant_classes, codons, ...}."""
    hit = (results["kind"] == match["kind"]).fillna(False)
    gene = results["gene"] if "gene" in results.columns else pd.Series(pd.NA, index=results.index)
    hit &= gene.isin(match["genes"]).fillna(False)
    v = variant_text(results, how)
    if "variant_classes" in match:
        cls = results["variant_class"] if "variant_class" in results.columns else pd.Series(pd.NA, index=results.index)
        hit &= cls.isin(match["variant_classes"]).fillna(False)
    if "variants" in match:
        hit &= v.isin(match["variants"]).fillna(False)
    if "except_variants" in match:
        hit &= ~v.isin(match["except_variants"]).fillna(False)
    if "codons" in match or "codon_ranges" in match:
        start, end = codon_span(v)
        in_codon = pd.Series(False, index=results.index)
        for c in match.get("codons", []):
            in_codon |= ((start <= c) & (end >= c)).fillna(False)
        for lo, hi in match.get("codon_ranges", []):
            in_codon |= ((start <= hi) & (end >= lo)).fillna(False)
        hit &= in_codon
    return hit.astype(bool)


def drug_matches(regimens: pd.Series, names: list[str]) -> pd.Series:
    """True for each regimen whose drug text names one of the drugs as a whole name, ignoring case: "Drug A" is
    found in "Drug B, Drug A" and in "Drug A Hydrochloride", but not inside a longer word ("Pro-Drug A")."""
    text_ = regimens.astype("string")
    hit = pd.Series(False, index=regimens.index)
    for name in names:
        pattern = rf"(?<![A-Za-z0-9-]){re.escape(name.strip())}(?![A-Za-z0-9-])"
        hit |= text_.str.contains(pattern, case=False, regex=True).fillna(False).astype(bool)
    return hit


# ---------------------------------------------------------------------------
# The first test (the entry event)
# ---------------------------------------------------------------------------
def panel_sizes(panels: pd.DataFrame) -> pd.Series:
    """How many genes each panel reads."""
    if panels is None or len(panels) == 0:
        return pd.Series(dtype="float64")
    return panels.dropna(subset=["panel_id", "gene"]).groupby("panel_id")["gene"].nunique()


def first_tests(blind: OutcomeBlindTables, settings: Settings) -> pd.DataFrame:
    """One row per patient with an entry-modality test: the first test by the settings' tie-breaks, with
    the patient's number of such tests and how many have a withheld time. A withheld time always ranks
    last. If tests are still tied after every tie-break the settings give, the lowest test_id is taken,
    and the row says so (`tied`)."""
    modalities = settings.read("entry.modalities")
    order = settings.read("entry.first_test_order")
    t = blind.tests
    if "modality" not in t.columns:
        raise MissingFieldError("tests.modality")
    t = t[t["modality"].isin(modalities).fillna(False)].copy()
    for col in ("days_to_entry_event", "days_from_birth"):
        if col not in t.columns:
            t[col] = pd.Series(pd.NA, index=t.index, dtype="Float64")     # the source has no such times
    t["_day_missing"] = t["days_to_entry_event"].isna()
    t["_panel_genes"] = t["panel_id"].map(panel_sizes(blind.panels)) if "panel_id" in t.columns else np.nan
    keys, ascending = ["patient_id", "_day_missing"], [True, True]
    tie_cols = []
    for k in order:
        col = {"earliest_day": "days_to_entry_event", "largest_panel": "_panel_genes", "lowest_test_order": "test_order"}[k]
        if col not in t.columns:
            continue
        t[f"_k_{k}"] = pd.to_numeric(t[col], errors="coerce").astype("float64")
        if k == "largest_panel":
            t[f"_k_{k}"] = -t[f"_k_{k}"]
        keys.append(f"_k_{k}")
        ascending.append(True)
        tie_cols.append(f"_k_{k}")
    t = t.sort_values(keys + ["test_id"], ascending=ascending + [True], na_position="last", kind="mergesort")
    grouped = t.groupby("patient_id", sort=False)
    first = grouped.head(1).set_index("patient_id")
    n_tests = grouped.size()
    n_missing = grouped["_day_missing"].sum()
    tie_key = t[["patient_id", "_day_missing"] + tie_cols].astype("string").fillna("<NA>").agg("|".join, axis=1)
    lead = tie_key.groupby(t["patient_id"]).transform("first")
    n_tied = (tie_key == lead).groupby(t["patient_id"]).sum()
    out = pd.DataFrame({
        "test_id": first["test_id"],
        "entry_day": first["days_to_entry_event"].astype("Float64"),
        "entry_from_birth": first["days_from_birth"].astype("Float64"),
        "entry_event": first["entry_event"] if "entry_event" in first.columns else pd.NA,
        "panel_id": first["panel_id"] if "panel_id" in first.columns else pd.NA,
        "after_follow_up_end": first["after_follow_up_end"] if "after_follow_up_end" in first.columns else pd.NA,
        "linked_to_other_cancers": first["linked_to_other_cancers"] if "linked_to_other_cancers" in first.columns else pd.NA,
        "n_tests": n_tests,
        "n_time_missing": n_missing.astype(int),
        "tied": n_tied > 1,
    })
    out.index.name = "patient_id"
    return out


# ---------------------------------------------------------------------------
# Panel coverage and genotype groups
# ---------------------------------------------------------------------------
def panel_coverage(panels: pd.DataFrame, required: dict, when_only_genes: str,
                   fusion_panels: dict | None = None) -> pd.DataFrame:
    """For each panel: whether it covers every required gene (and codon), and whether that rests on a
    gene being listed without the codons it reads (`gene_level`). For a gene in `fusion_panels`, a panel
    covers it only if it is listed there as detecting the gene's fusions."""
    fusion_panels = fusion_panels or {}
    rows = []
    if panels is None or len(panels) == 0:
        return pd.DataFrame(columns=["covered", "gene_level"])
    has_regions = {"codon_from", "codon_to"} <= set(panels.columns)
    for panel, p in panels.groupby("panel_id"):
        covered, gene_level = True, False
        for gene, codons in required.items():
            g = p[p["gene"] == gene]
            if len(g) == 0 or (gene in fusion_panels and panel not in fusion_panels[gene]):
                covered = False
                break
            if codons == "whole_gene":
                continue
            regions = g.dropna(subset=["codon_from", "codon_to"]) if has_regions else g.iloc[0:0]
            if len(regions):
                if not all(((regions["codon_from"] <= c) & (regions["codon_to"] >= c)).any() for c in codons):
                    covered = False
                    break
            else:
                gene_level = True
        rows.append({"panel_id": panel, "covered": covered, "gene_level": covered and gene_level})
    out = pd.DataFrame(rows).set_index("panel_id")
    if when_only_genes not in ("accept_gene_level", "deviation"):
        raise ValueError(when_only_genes)
    return out


def assign_genotype(blind: OutcomeBlindTables, settings: Settings, first: pd.DataFrame,
                    coverage_known: bool = True) -> pd.DataFrame:
    """One row per patient with a first test: whether each group's rules match (before precedence), the
    group by precedence, coverage, and the group after the coverage rule (blank = not assignable).

    coverage_known=False is for a source with no panel coverage at all (its profile says so): groups are
    then assigned without the coverage rule, and every patient is marked unconfirmed (`covered` False)."""
    groups = settings.read("genotype.groups")
    how = settings.read("genotype.variant_text")
    required = settings.read("genotype.coverage.required")
    applies_to = settings.read("genotype.coverage.applies_to")
    when_only_genes = settings.read("genotype.coverage.when_only_genes_listed")
    fusion_panels = settings.read("genotype.coverage.fusion_panels")
    settings.read("genotype.results_from")       # entry_test: only the first test's results count

    res = blind.results
    res = res[res["test_id"].isin(first["test_id"])]
    patient_of_test = first.reset_index().drop_duplicates("test_id").set_index("test_id")["patient_id"]
    out = pd.DataFrame(index=first.index)
    for g in groups:
        if g.get("none_of_the_above"):
            continue
        hit = pd.Series(False, index=res.index)
        for m in g["any_of"]:
            hit |= result_matches(res, m, how)
        matched = set(patient_of_test.reindex(res.loc[hit, "test_id"]).dropna())
        out[f"match:{g['name']}"] = out.index.isin(list(matched))
    rest = next(g["name"] for g in groups if g.get("none_of_the_above"))
    group = pd.Series(rest, index=out.index, dtype="object")
    for g in reversed([g for g in groups if not g.get("none_of_the_above")]):
        group[out[f"match:{g['name']}"]] = g["name"]
    cov = panel_coverage(blind.panels, required, when_only_genes, fusion_panels)
    covered = first["panel_id"].map(cov["covered"]).astype("boolean").fillna(False).astype(bool)
    gene_level = first["panel_id"].map(cov["gene_level"]).astype("boolean").fillna(False).astype(bool)
    out["group_by_rules"] = group
    out["covered"] = covered
    out["gene_level_coverage"] = gene_level & covered
    needs = pd.Series(True, index=out.index) if applies_to == "all_groups" else (group == rest)
    out["genotype"] = group.where(~needs | covered, None) if coverage_known else group
    return out


# ---------------------------------------------------------------------------
# The cohort
# ---------------------------------------------------------------------------
@dataclass
class Cohort:
    patients: pd.DataFrame            # one row per patient in the tables, indexed by patient_id
    rules: list[dict]                 # the settings' rules, in order
    first: pd.DataFrame               # first tests
    genotype: pd.DataFrame            # genotype assignment for patients with a first test
    notes: list[str] = field(default_factory=list)

    @property
    def eligible(self) -> pd.DataFrame:
        return self.patients[self.patients["excluded_at"].isna()]

    def flow(self, start_label: str = "Patients in the standard tables",
             before: list[tuple[str, int]] | None = None) -> pd.DataFrame:
        """Exact counts at each step: excluded and remaining. `before` adds the adapter's own steps
        (for example patients the source holds without a cancer of the study's kind) ahead of the rules.
        Exact counts: for checking only, never published as they stand (privacy.public_flow)."""
        before = before or []
        start = len(self.patients) + sum(n for _, n in before)
        rows = [{"step": "start", "label": start_label, "excluded": 0, "remaining": start}]
        remaining = start
        for label, n in before:
            remaining -= n
            rows.append({"step": "source", "label": label, "excluded": n, "remaining": remaining})
        for r in self.rules:
            n = int((self.patients["excluded_at"] == r["id"]).sum())
            remaining -= n
            rows.append({"step": r["id"], "label": r["label"], "excluded": n, "remaining": remaining})
        return pd.DataFrame(rows)


def build_cohort(blind: OutcomeBlindTables, settings: Settings, coverage_known: bool = True,
                 entry_times_given: bool = True) -> Cohort:
    """Apply the settings' cohort rules in order. Each patient is counted at the first rule they fail.

    entry_times_given=False is for a source with no entry event at all (its profile says so): only the naive
    clock can run, so the entry-time check passes everyone who has a test, and the screen says why."""
    rules = settings.read("cohort.rules")
    before_dx = settings.read("entry.before_diagnosis")
    p = blind.patients.drop_duplicates("patient_id").set_index("patient_id")
    first = first_tests(blind, settings).reindex(p.index)
    has_first = first["test_id"].notna()
    geno = assign_genotype(blind, settings, first[has_first], coverage_known)
    p_geno = geno.reindex(p.index)
    entry_known = first["entry_day"].notna()
    if before_dx == "exclude":
        entry_known &= (first["entry_day"] >= 0).fillna(False)
    if not entry_times_given:
        entry_known = has_first.copy()
    after_end = first["after_follow_up_end"].astype("boolean")
    structural = {
        "has_entry_test": has_first,
        "entry_time_known": entry_known.astype(bool),
        "outcome_recorded": p[OUTCOME_RECORDED].astype("boolean").fillna(False).astype(bool),
        "entry_before_follow_up_end": ~(after_end.fillna(False).astype(bool)),
        "genotype_assignable": p_geno["genotype"].notna(),
    }
    out = pd.DataFrame(index=p.index)
    for r in rules:
        ok = pd.Series(True, index=p.index)
        for cond in r["all_of"]:
            ok &= structural[cond["check"]] if "check" in cond else field_condition(p, cond)
        out[f"r:{r['id']}"] = ok.astype(bool)
    excluded = pd.Series([None] * len(p), index=p.index, dtype="object")
    for r in reversed(rules):
        excluded[~out[f"r:{r['id']}"]] = r["id"]
    out["excluded_at"] = excluded.astype(object).where(excluded.notna(), None)     # None = passes every rule
    out["genotype"] = p_geno["genotype"]
    for c in [c for c in geno.columns if c.startswith("match:")] + ["covered", "gene_level_coverage", "group_by_rules"]:
        out[c] = p_geno[c]
    out["entry_day"] = first["entry_day"]
    out["entry_from_birth"] = first["entry_from_birth"]
    out["first_test_id"] = first["test_id"]
    out["first_panel_id"] = first["panel_id"]
    # flags: patients kept but counted
    out["f:entry_before_diagnosis"] = (first["entry_day"] < 0).fillna(False).astype(bool)
    out["f:unranked_first_test"] = ((first["n_time_missing"] > 0) & first["entry_day"].notna()).fillna(False).astype(bool)
    out["f:first_test_tied"] = first["tied"].astype("boolean").fillna(False).astype(bool)
    out["f:after_follow_up_end_unknown"] = (has_first & after_end.isna()).astype(bool)
    linked = first["linked_to_other_cancers"].astype("boolean")
    out["f:first_test_linked_to_other_cancers"] = linked.fillna(False).astype(bool)
    if "censored_at_age_limit" in p.columns:
        out["f:censored_at_age_limit"] = p["censored_at_age_limit"].astype("boolean").fillna(False).astype(bool)
    return Cohort(patients=out, rules=rules, first=first, genotype=geno)
