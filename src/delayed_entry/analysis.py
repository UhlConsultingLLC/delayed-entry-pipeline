"""The analysis cohort, built after the plan is locked (design plan section 8): the checks that need outcomes,
then who each pass and each analysis keeps, and the adjustment variables coded as the settings say.

    from delayed_entry.analysis import build_analysis
    data = build_analysis(tables, settings, profile, source_steps)
    data.cohort       # one row per eligible patient: group, entry, exit, event, in_analysis
    data.checks       # each check that needed outcomes, with its verdict (texts are masked)
    data.d1           # the smallest-group rule, applied again to the patients each analysis keeps

Generalized from the CRC BPC project's sql/05_analysis_cohort.sql, sql/06_secondary.sql and crcbpc.design (D2
and D1). The cohort rules, groups and first test are phase 1's (delayed_entry.cohort), run on the same tables;
the outcomes are read only here, after the lock. Every rule and threshold comes from the settings.

The checks that need outcomes, in order:
1. time at risk from the times themselves: a patient whose follow-up ends on or before entry is excluded, at the
   cohort rule that holds the structural check `entry_before_follow_up_end` (the same step of the flow diagram);
2. survival after the entry event: if fewer eligible patients than the settings' minimum died within the
   settings' window of their entry event, the source in effect required that survival, and entry moves later by
   the window (or a warning only, as the settings say); a patient left with no time at risk stays eligible but
   leaves the analysis, and is counted;
3. the horizon against follow-up: a warning if, at the horizon, a group's corrected risk set is below the
   stability rule's size;
4. the smallest-group rule, applied again to the patients each analysis keeps.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from delayed_entry.cohort import Cohort, build_cohort, field_condition, result_matches
from delayed_entry.covariates import (Coded, apply_model_bands, code_adjustment, missing_any, missing_data_rule,
                                      model_bands)
from delayed_entry.privacy import show
from delayed_entry.profile import CapabilityProfile
from delayed_entry.screen import (BURDEN_HIGH, BURDEN_LOW, COMPARE, DESCRIBE, NOT_RUN, PRIMARY, biomarker_status,
                                  burden_problems, mutation_burden, panel_sizes_given)
from delayed_entry.settings import Settings
from delayed_entry.tables import StandardTables, outcome_blind

STRUCTURAL_TIME_CHECK = "entry_before_follow_up_end"
SWITCHED_OFF = "switched off"
RUN = "run"
PASS, APPLIED, WARNING = "pass", "applied", "warning"


class RunStop(RuntimeError):
    """The run can't go on as the locked plan says. A person decides what is done instead, and a new plan
    version records it."""


@dataclass(frozen=True)
class Grouping:
    column: str
    levels: tuple[str, ...]          # the reference first
    reference: str

    def comparisons(self) -> list[str]:
        return [comparison(g, self.reference) for g in self.levels[1:]]

    def to_dict(self) -> dict:
        return {"column": self.column, "levels": list(self.levels), "reference": self.reference}


def comparison(group: str, reference: str) -> str:
    return f"{group} vs {reference}"


@dataclass
class PlannedAnalysis:
    """One analysis as the run carries it out."""
    id: str
    kind: str
    label: str
    status: str                                   # run, not run, switched off
    why: str = ""
    grouping: Grouping | None = None
    member_column: str | None = None              # the secondary table's column that says who it keeps
    described: dict = field(default_factory=dict)  # a group or comparison: why it is described, not labeled
    adjusted: bool = False                        # pass 3a too
    stratify_by: str | None = None
    parent: str | None = None

    def to_dict(self) -> dict:
        return {"id": self.id, "kind": self.kind, "label": self.label, "status": self.status, "why": self.why,
                "grouping": self.grouping.to_dict() if self.grouping else None, "member_column": self.member_column,
                "described": dict(self.described), "adjusted": self.adjusted, "stratify_by": self.stratify_by,
                "parent": self.parent}


@dataclass
class OutcomeCheck:
    check: str
    verdict: str          # pass, applied, warning
    detail: str           # masked: safe to show

    def to_dict(self) -> dict:
        return {"check": self.check, "verdict": self.verdict, "detail": self.detail}


@dataclass
class AnalysisData:
    cohort: pd.DataFrame                 # eligible patients, sorted by patient_id
    covariates: pd.DataFrame             # eligible patients: patient_id and one column per adjustment variable
    coded: Coded                         # the same, with the model's terms
    secondary: pd.DataFrame              # patients in the analysis: patient_id and each analysis's member column
    clocks: dict[str, pd.DataFrame]      # treatment-clock cohorts, by analysis id
    analyses: list[PlannedAnalysis]
    d1: pd.DataFrame                     # analysis, comparison, n_reference, n_group, decision (exact counts)
    checks: list[OutcomeCheck]
    covariate_plan: dict
    flow: pd.DataFrame                   # exact counts: never published as they stand
    eligible: pd.DataFrame               # the cohort's own columns (rules, flags) for the eligible patients
    patients: pd.DataFrame               # the patients table, indexed by patient_id (eligible patients)
    d2: dict
    genotype: Grouping
    biomarker_counts: dict               # analysis id: {group: n} over the patients in the analysis, with
                                         # the unknown group and those the analysis sets aside ("not in <id>")
    min_count: int
    blind_cohort: Cohort | None = None

    def analysis(self, analysis_id: str) -> PlannedAnalysis | None:
        return next((a for a in self.analyses if a.id == analysis_id), None)

    def in_analysis(self) -> pd.DataFrame:
        return self.cohort[self.cohort["in_analysis"]]


# ---------------------------------------------------------------------------
# The smallest-group rule
# ---------------------------------------------------------------------------
def d1_rows(analysis: str, counts: dict[str, int], grouping: Grouping, min_group: int, if_short,
            described_only=()) -> list[dict]:
    """The smallest-group rule from counts alone: each comparison with the reference is compared if both groups
    have at least `min_group` patients; otherwise described, or blocked, as `if_short` says (a mapping by group,
    or one word for every group). Groups the settings describe only are never compared."""
    rows = []
    ref = grouping.reference
    for g in grouping.levels[1:]:
        short = [x for x in (g, ref) if counts.get(x, 0) < min_group]
        if g in described_only:
            decision = DESCRIBE
        elif not short:
            decision = COMPARE
        else:
            rule = {x: (if_short[x] if isinstance(if_short, dict) else if_short) for x in short}
            decision = "blocked" if "block" in rule.values() else DESCRIBE
        rows.append({"analysis": analysis, "comparison": comparison(g, ref), "n_reference": int(counts.get(ref, 0)),
                     "n_group": int(counts.get(g, 0)), "decision": decision})
    return rows


def group_counts(frame: pd.DataFrame, column: str, names) -> dict[str, int]:
    vc = frame[column].value_counts()
    return {g: int(vc.get(g, 0)) for g in names}


# ---------------------------------------------------------------------------
# The cohort with outcomes
# ---------------------------------------------------------------------------
def _time_check_rule(c: Cohort) -> str:
    for r in c.rules:
        if any(cond.get("check") == STRUCTURAL_TIME_CHECK for cond in r["all_of"]):
            return r["id"]
    raise RunStop(f"no cohort rule holds the structural check '{STRUCTURAL_TIME_CHECK}'")


def apply_time_check(c: Cohort, no_time: pd.Series) -> int:
    """Fold "no time at risk, judged from the times" into the cohort rule holding the structural time check, and
    recount where each patient leaves. Returns how many eligible patients it removes beyond the source's record."""
    rule = _time_check_rule(c)
    before = c.patients["excluded_at"].isna()
    col = f"r:{rule}"
    c.patients[col] = c.patients[col] & ~no_time.reindex(c.patients.index).fillna(False).astype(bool)
    excluded = pd.Series([None] * len(c.patients), index=c.patients.index, dtype="object")
    for r in reversed(c.rules):
        excluded[~c.patients[f"r:{r['id']}"]] = r["id"]
    c.patients["excluded_at"] = excluded.astype(object).where(excluded.notna(), None)
    return int((before & c.patients["excluded_at"].notna()).sum())


def _outcomes(tables: StandardTables) -> tuple[pd.Series, pd.Series]:
    p = tables.patients.drop_duplicates("patient_id").set_index("patient_id")
    exit_day = pd.to_numeric(p["days_to_death_or_last_followup"], errors="coerce").astype("float64")
    event = pd.to_numeric(p["vital_status"], errors="coerce").astype("float64")
    return exit_day, event


def survival_after_entry(eligible: pd.DataFrame, settings: Settings) -> dict:
    """The CRC BPC plan's design check D2, with the settings' window and minimum: eligible patients who died
    within the window of their entry event, counted from the event itself (a report before diagnosis counts
    from its own day, not from day 0). Exact count: for checking, and shown masked."""
    window = settings.read("survival_after_entry.window_days")
    minimum = settings.read("survival_after_entry.min_deaths")
    if_fewer = settings.read("survival_after_entry.if_fewer")
    early = int(((eligible["event"] == 1) & (eligible["exit_day"] - eligible["report_day"] < window)).sum())
    required = early < minimum
    return {"early_deaths": early, "window_days": window, "required": required,
            "entry_shift_days": window if (required and if_fewer == "move_entry_by_window") else 0,
            "if_fewer": if_fewer}


# ---------------------------------------------------------------------------
# Who each analysis keeps
# ---------------------------------------------------------------------------
def keep_mask(kept: pd.DataFrame, patients: pd.DataFrame, results: pd.DataFrame, keep: dict,
              settings: Settings) -> pd.Series:
    """A subset analysis's rule (the screen's, applied to the patients in the analysis)."""
    if "within_months" in keep:
        limit = keep["within_months"] * settings.read("time.days_per_year") / 12
        return (kept["report_day"] <= limit).fillna(False).astype(bool)
    if "field" in keep:
        return field_condition(patients.loc[kept.index], keep)
    if "exclude_overlap" in keep:
        both = kept[[f"match:{g}" for g in keep["exclude_overlap"]]].astype(bool).all(axis=1)
        return ~both
    if "exclude_from_group" in keep:
        how = settings.read("genotype.variant_text")
        res = results[results["test_id"].isin(kept["first_test_id"].dropna())]
        hit = res.loc[result_matches(res, keep["if_result"], how), "test_id"]
        has = kept["first_test_id"].isin(set(hit))
        return ~((kept["genotype"] == keep["exclude_from_group"]) & has)
    raise ValueError(f"unknown subset rule {keep}")


def first_line(treatments: pd.DataFrame | None) -> pd.DataFrame:
    """Each patient's first regimen by the source's order among those starting on or after diagnosis (a regimen
    whose start is withheld can't be placed and is skipped), as the CRC BPC plan's S2 did."""
    if treatments is None or not {"patient_id", "days_to_start", "regimen_order"} <= set(treatments.columns):
        return pd.DataFrame(columns=["patient_id", "days_to_start"])
    ok = treatments[treatments["days_to_start"].notna() & (treatments["days_to_start"] >= 0)]
    first = ok.sort_values(["patient_id", "regimen_order"], kind="mergesort").drop_duplicates("patient_id")
    return first[["patient_id", "days_to_start"]]


def treatment_clock(cohort: pd.DataFrame, treatments: pd.DataFrame | None) -> pd.DataFrame:
    """The patients in the analysis with time counted from their first-line start: entry is the primary entry
    from that start (a report before the start enters at day 0), exit the follow-up from that start. In this
    clock's analysis only with a first line, a known outcome and entry before exit."""
    x = cohort[cohort["in_analysis"]].copy()
    start = first_line(treatments).set_index("patient_id")["days_to_start"].astype("float64")
    s = x["patient_id"].map(start)
    entry = x["entry_day"] - s
    out = pd.DataFrame({"patient_id": x["patient_id"], "genotype": x["genotype"],
                        "start_day": s, "entry_day": entry.where(~(entry < 0), 0.0),
                        "exit_day": x["exit_day"] - s, "event": x["event"]})
    out["in_analysis"] = (s.notna() & out["exit_day"].notna() & out["event"].notna()
                          & (out["entry_day"] < out["exit_day"])).fillna(False).astype(bool)
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Building it all
# ---------------------------------------------------------------------------
def build_analysis(tables: StandardTables, settings: Settings, profile: CapabilityProfile,
                   source_steps: list[tuple[str, int]] | None = None) -> AnalysisData:
    min_count = max(settings.read("privacy.min_count"), profile.privacy_min_count or 0)
    min_group = settings.read("smallest_group.patients")
    checks: list[OutcomeCheck] = []
    if profile.entry_event == "none":
        raise RunStop("the source has no entry event: the corrected clock can't run, and no bias can be estimated")

    blind = outcome_blind(tables)
    c = build_cohort(blind, settings, coverage_known=profile.panel_coverage != "none", entry_times_given=True)
    exit_day, event = _outcomes(tables)
    raw = c.patients["entry_day"].astype("float64")
    before_dx = settings.read("entry.before_diagnosis")
    entry0 = raw.where(~(raw < 0), 0.0) if before_dx == "enter_at_day_0" else raw
    ex = exit_day.reindex(c.patients.index)
    no_time = raw.notna() & ex.notna() & (ex <= entry0)

    # 1. Time at risk from the times themselves
    removed = apply_time_check(c, no_time)
    checks.append(OutcomeCheck(
        "time at risk from the times themselves", PASS if removed == 0 else APPLIED,
        "every eligible patient's follow-up ends after entry, as the source's own record said" if removed == 0 else
        f"{show(removed, min_count)} more eligible patients have follow-up that ends on or before entry: excluded at "
        f"the cohort rule '{_time_check_rule(c)}' and counted there"))
    flow = c.flow(start_label="Patients in the source", before=list(source_steps or []))

    e = c.eligible
    cohort = pd.DataFrame({
        "patient_id": e.index.astype(str),
        "genotype": e["genotype"].astype(object).to_numpy(),
        "report_day": raw.reindex(e.index).to_numpy(),
        "exit_day": exit_day.reindex(e.index).to_numpy(),
        "event": event.reindex(e.index).to_numpy(),
    })
    cohort["entry_day"] = entry0.reindex(e.index).to_numpy()

    # 2. Survival after the entry event
    d2 = survival_after_entry(cohort, settings)
    d2["reported_as"] = show(d2["early_deaths"], min_count)
    window = d2["window_days"]
    if d2["entry_shift_days"]:
        cohort["entry_day"] = cohort["report_day"] + d2["entry_shift_days"]
    cohort["in_analysis"] = (cohort["exit_day"] > cohort["entry_day"]).fillna(False).astype(bool)
    left = int((~cohort["in_analysis"]).sum())
    deaths = settings.read("survival_after_entry.min_deaths")
    if not d2["required"]:
        checks.append(OutcomeCheck("survival after the entry event", PASS,
                                   f"{d2['reported_as']} eligible patients died within {window:g} days of their entry "
                                   f"event (the settings' minimum is {deaths}): entry stays at the event"))
    elif d2["entry_shift_days"]:
        checks.append(OutcomeCheck("survival after the entry event", APPLIED,
                                   f"fewer than {deaths} eligible patients died within {window:g} days of their entry "
                                   f"event, so the source in effect required that survival: entry moves {window:g} days "
                                   f"later, as planned; {show(left, min_count)} eligible patients are left with no time "
                                   "at risk and leave the analysis"))
    else:
        checks.append(OutcomeCheck("survival after the entry event", WARNING,
                                   f"fewer than {deaths} eligible patients died within {window:g} days of their entry "
                                   "event: the source may have required that survival (the settings ask for a warning "
                                   "only)"))
    cohort = cohort.sort_values("patient_id", kind="mergesort").reset_index(drop=True)
    eligible = e.loc[cohort["patient_id"]]
    patients = blind.patients.drop_duplicates("patient_id").set_index("patient_id").loc[cohort["patient_id"]]

    # The genotype grouping: the reference first, then the settings' order.
    names = [g["name"] for g in settings.read("genotype.groups")]
    ref = settings.read("genotype.reference")
    genotype = Grouping("genotype", tuple([ref] + [g for g in names if g != ref]), ref)
    described_only = settings.read("genotype.described_only")
    inside = cohort[cohort["in_analysis"]]

    # 3. The horizon against follow-up
    dpy = settings.read("time.days_per_year")
    horizon = settings.read("time.horizon_years") * dpy
    risk_set = settings.read("stability.risk_set")
    at_h = (inside["entry_day"] < horizon) & (inside["exit_day"] >= horizon)
    thin = [g for g in genotype.levels if int(at_h[inside["genotype"] == g].sum()) < risk_set]
    checks.append(OutcomeCheck(
        "the horizon against follow-up", PASS if not thin else WARNING,
        f"at the horizon ({settings.read('time.horizon_years'):g} years) every group's corrected risk set has at least "
        f"{risk_set} patients (the stability rule's size)" if not thin else
        f"at the horizon ({settings.read('time.horizon_years'):g} years) the corrected risk set has fewer than {risk_set} "
        f"patients (the stability rule's size) in: {', '.join(thin)}; read their RMST with that in mind"))

    # 4. The smallest-group rule, again: the primary analysis
    d1: list[dict] = []
    counts = group_counts(inside, "genotype", genotype.levels)
    d1 += d1_rows(PRIMARY, counts, genotype, min_group, settings.read("smallest_group.if_short"), described_only)

    # The adjustment variables
    coded = code_adjustment(patients, settings)
    rule, share = missing_data_rule(coded, settings)
    if rule != "complete case" and not settings.has("missing_data.imputation_model"):
        raise RunStop(f"more than {settings.read('missing_data.complete_case_max_share'):.0%} of eligible patients lack "
                      "an adjustment variable, so the plan needs multiple imputation, and the settings don't say how "
                      "(missing_data.imputation_model): add it in a new plan version")
    model_ids = pd.Index(inside["patient_id"])
    if rule == "complete case":                  # the model's patients: those with every variable
        counted = model_ids[~missing_any(coded).reindex(model_ids).to_numpy()]
    else:                                        # imputing: everyone in the analysis, counted where the band is known
        counted = model_ids
    try:
        bands = model_bands(coded, settings, counted)
    except ValueError as err:
        raise RunStop(str(err)) from None
    coded = apply_model_bands(coded, bands)
    covariates = coded.frame.reset_index().rename(columns={"index": "patient_id"})
    covariate_plan = {"missing_data_rule": rule, "missing_share": share, "year_bands": bands,
                      "terms": [t.to_dict() for t in coded.terms]}

    # The other analyses
    secondary = pd.DataFrame({"patient_id": inside["patient_id"].to_numpy()})
    inside_e = eligible.loc[inside["patient_id"]].copy()
    inside_e["report_day"] = inside.set_index("patient_id")["report_day"]
    clocks, analyses, biomarker_counts, states = {}, [], {}, {}
    all_a = settings.read("analyses")
    by_id = {a["id"]: a for a in all_a}
    described = {}
    for r in d1:
        if r["decision"] != COMPARE:
            g = r["comparison"].split(" vs ")[0]
            described[g] = "settings" if g in described_only else "smallest group"
    analyses.append(PlannedAnalysis(PRIMARY, "primary", "the primary analysis", RUN, grouping=genotype,
                                    described=described))
    for a in all_a:
        if not a["run"]:
            analyses.append(PlannedAnalysis(a["id"], a["kind"], a["label"], SWITCHED_OFF, "switched off in the settings"))
            continue
        k = a["kind"]
        if k == "subset":
            mask = keep_mask(inside_e, patients, blind.results, a["keep"], settings)
            col = f"keep_{a['id']}"
            secondary[col] = mask.to_numpy()
            kept = inside[mask.to_numpy()]
            d1 += d1_rows(a["id"], group_counts(kept, "genotype", genotype.levels), genotype, min_group,
                          a["if_short"], described_only)
            analyses.append(PlannedAnalysis(a["id"], k, a["label"], RUN, grouping=genotype, member_column=col))
        elif k == "treatment_clock":
            clock = treatment_clock(cohort, blind.treatments)
            clocks[a["id"]] = clock
            d1 += d1_rows(a["id"], group_counts(clock[clock["in_analysis"]], "genotype", genotype.levels), genotype,
                          min_group, a["if_short"], described_only)
            analyses.append(PlannedAnalysis(a["id"], k, a["label"], RUN, grouping=genotype))
        elif k == "biomarker":
            states[a["id"]] = _biomarker(a, blind, c, cohort, inside, secondary, settings, profile, min_group, d1,
                                         analyses, biomarker_counts)
        elif k == "adjusted_variant":
            if a["stratify_by"] not in coded.frame.columns:
                raise RunStop(f"{a['id']} stratifies by {a['stratify_by']}, which isn't an adjustment variable")
            analyses.append(PlannedAnalysis(a["id"], k, a["label"], RUN, grouping=genotype, stratify_by=a["stratify_by"]))
        elif k in ("ph_check", "quasi_independence", "descriptive_delay"):
            analyses.append(PlannedAnalysis(a["id"], k, a["label"], RUN, grouping=genotype))
    for a in all_a:
        if a["run"] and a["kind"] == "biomarker_variant":
            _biomarker_variant(a, by_id[a["of"]], states.get(a["of"], {}), blind, c, inside, secondary, min_group,
                               d1, analyses)
    d1_frame = pd.DataFrame(d1, columns=["analysis", "comparison", "n_reference", "n_group", "decision"])
    blocked = d1_frame[d1_frame["decision"] == "blocked"]
    if len(blocked):
        raise RunStop("the smallest-group rule, applied again, blocks: " + "; ".join(
            f"{r.analysis}: {r.comparison}" for r in blocked.itertuples()) + ". A person decides what is done instead")
    order = {a.id: i for i, a in enumerate(analyses)}
    analyses.sort(key=lambda x: order[x.id])
    return AnalysisData(
        cohort=cohort, covariates=covariates, coded=coded, secondary=secondary, clocks=clocks, analyses=analyses,
        d1=d1_frame, checks=checks, covariate_plan=covariate_plan, flow=flow, eligible=eligible, patients=patients,
        d2=d2, genotype=genotype, biomarker_counts=biomarker_counts, min_count=min_count, blind_cohort=c)


def _biomarker(a: dict, blind, c: Cohort, cohort, inside, secondary, settings, profile, min_group, d1, analyses,
               biomarker_counts) -> dict:
    """A biomarker analysis: each patient's group by their results (the screen's rule), whether it uses its
    fallback (judged on the eligible patients, as the screen did), and the smallest-group rule. The fallback
    compares mutation burden if the source gives panel sizes and the settings say so, else describes the groups."""
    status_e, _, _ = biomarker_status(blind, c, a, a["timing"], False)
    unknown = int((status_e == a["unknown_group"]).sum())
    fallback = len(status_e) > 0 and unknown / len(status_e) > settings.read("biomarker.unknown_max_share")
    names = [g["name"] for g in a["groups"]]
    grouping = Grouping(f"group_{a['id']}", tuple([a["reference"]] + [g for g in names if g != a["reference"]]),
                        a["reference"])
    ids = inside["patient_id"]
    status = ids.map(status_e).astype(object)
    secondary[grouping.column] = status.where(status.notna(), None).to_numpy()
    counts = {g: int((status == g).sum()) for g in grouping.levels}
    biomarker_counts[a["id"]] = {**counts, a["unknown_group"]: int((status == a["unknown_group"]).sum()),
                                 f"not in {a['id']}": int(status.isna().sum())}
    if fallback:
        form = a["fallback"]["when_panel_sizes_given" if panel_sizes_given(blind, profile) else "otherwise"]
        if form == "mutation_burden":
            _burden(a, blind, c, inside, secondary, min_group, d1, analyses, biomarker_counts)
            _not_run_under_fallback(a, analyses)
            return {"fallback": True}
        rows =[{"analysis": a["id"], "comparison": comparison(g, grouping.reference), "n_reference": counts[grouping.reference],
                 "n_group": counts[g], "decision": DESCRIBE} for g in grouping.levels[1:]]
        d1 += rows
        described = {g: "biomarker fallback" for g in grouping.levels}
        described |= {comparison(g, grouping.reference): "biomarker fallback" for g in grouping.levels[1:]}
        analyses.append(PlannedAnalysis(a["id"], "biomarker", a["label"], RUN, "the biomarker's fallback: described only",
                                        grouping=grouping, member_column=grouping.column, described=described,
                                        adjusted=True))
        _not_run_under_fallback(a, analyses)
        return {"fallback": True}
    rows = d1_rows(a["id"], counts, grouping, min_group, a["if_short"])
    d1 += rows
    described = {g: "smallest group" for g in grouping.levels} if any(r["decision"] != COMPARE for r in rows) else {}
    analyses.append(PlannedAnalysis(a["id"], "biomarker", a["label"], RUN, grouping=grouping,
                                    member_column=grouping.column, described=described, adjusted=True))
    return {"fallback": False, "grouping": grouping}


def _not_run_under_fallback(a: dict, analyses: list) -> None:
    for other in a["not_run_under_fallback"]:
        analyses.append(PlannedAnalysis(other, "biomarker_variant", other, NOT_RUN,
                                        f"{a['id']} uses its fallback, and the settings don't run {other} then",
                                        parent=a["id"]))


def burden_label(a: dict) -> str:
    t = a["mutation_burden"]["threshold_per_megabase"]
    return (f"{BURDEN_HIGH} ({t:g} or more mutations per megabase) vs {BURDEN_LOW}, a hypermutation proxy in place "
            f"of {a['label']}")


def _burden(a: dict, blind, c: Cohort, inside, secondary, min_group, d1, analyses, biomarker_counts) -> None:
    """A biomarker's fallback by mutation burden (the screen's rule, screen.mutation_burden): high against low, the
    reference, on the entry test; the unknown group is described. Not run if the burden can't be counted (the
    screen blocked it too)."""
    problems = burden_problems(blind)
    if problems:
        biomarker_counts.pop(a["id"], None)
        analyses.append(PlannedAnalysis(a["id"], "biomarker", a["label"], NOT_RUN,
                                        "its fallback, mutation burden, can't be counted: " + "; ".join(problems)))
        return
    group_e, _ = mutation_burden(blind, c, a)
    grouping = Grouping(f"group_{a['id']}", (BURDEN_LOW, BURDEN_HIGH), BURDEN_LOW)
    status = inside["patient_id"].map(group_e).astype(object)
    secondary[grouping.column] = status.where(status.notna(), None).to_numpy()
    counts = {g: int((status == g).sum()) for g in grouping.levels}
    biomarker_counts[a["id"]] = {**counts, a["unknown_group"]: int((status == a["unknown_group"]).sum()),
                                 f"not in {a['id']}": int(status.isna().sum())}
    rows = d1_rows(a["id"], counts, grouping, min_group, a["if_short"])
    d1 += rows
    described = {g: "smallest group" for g in grouping.levels} if any(r["decision"] != COMPARE for r in rows) else {}
    analyses.append(PlannedAnalysis(a["id"], "biomarker", burden_label(a), RUN,
                                    "the biomarker's fallback: mutation burden on the entry test",
                                    grouping=grouping, member_column=grouping.column, described=described,
                                    adjusted=True))


def _biomarker_variant(a: dict, parent: dict, state: dict, blind, c, inside, secondary, min_group, d1, analyses):
    if any(x.id == a["id"] for x in analyses):          # not run under the parent's fallback
        return
    if not state or state.get("fallback"):
        return
    timing = a.get("timing", parent["timing"])
    status_e, _, _ = biomarker_status(blind, c, parent, timing, a.get("drop_discordant", False))
    base: Grouping = state["grouping"]
    grouping = Grouping(f"group_{a['id']}", base.levels, base.reference)
    status = inside["patient_id"].map(status_e).astype(object)
    secondary[grouping.column] = status.where(status.notna(), None).to_numpy()
    counts = {g: int((status == g).sum()) for g in grouping.levels}
    d1 += d1_rows(a["id"], counts, grouping, min_group, parent["if_short"])
    analyses.append(PlannedAnalysis(a["id"], "biomarker_variant", a["label"], RUN, grouping=grouping,
                                    member_column=grouping.column, adjusted=True, parent=parent["id"]))


# ---------------------------------------------------------------------------
# The plan's decisions, against the run's
# ---------------------------------------------------------------------------
def compare_with_plan(data: AnalysisData, plan_decisions: dict) -> list[OutcomeCheck]:
    """The smallest-group rule applied again, against the locked plan's decisions. A comparison the plan compares
    and the run can't is a pre-planned change (described, by the settings' rule), stated here. An analysis the plan
    doesn't run must not run."""
    out = []
    changed = []
    for a in data.analyses:
        planned = plan_decisions.get(a.id, {})
        if a.status != RUN:
            if planned and not all(v in (NOT_RUN, SWITCHED_OFF) for v in planned.values()):
                raise RunStop(f"{a.id}: the locked plan runs it, but the run would not ({a.why})")
            continue
        if a.kind not in ("primary", "subset", "treatment_clock", "biomarker", "biomarker_variant"):
            continue
        if any(v in (NOT_RUN, SWITCHED_OFF) for v in planned.values()):
            raise RunStop(f"{a.id}: the locked plan doesn't run it, but the run would")
        rows = data.d1[data.d1["analysis"] == a.id]
        for r in rows.itertuples():
            before = planned.get(r.comparison)
            if before is None:
                continue
            if before.startswith(COMPARE) and r.decision != COMPARE:
                changed.append(f"{a.id}: {r.comparison} is now {r.decision}")
            elif not before.startswith(COMPARE) and r.decision == COMPARE:
                changed.append(f"{a.id}: {r.comparison} stays {before}, as locked")
    if changed:
        out.append(OutcomeCheck("the smallest-group rule on the analysis cohort", APPLIED,
                                "applied again to the patients each analysis keeps: " + "; ".join(changed)))
    else:
        out.append(OutcomeCheck("the smallest-group rule on the analysis cohort", PASS,
                                "applied again to the patients each analysis keeps: every decision is the locked plan's"))
    return out


def keep_plan_decisions(data: AnalysisData, plan_decisions: dict) -> None:
    """A comparison the locked plan describes stays described, even if the patients left would now allow it: the
    plan can't be loosened after the outcomes are read."""
    for i, r in data.d1.iterrows():
        before = plan_decisions.get(r["analysis"], {}).get(r["comparison"])
        if before is not None and not before.startswith(COMPARE) and r["decision"] == COMPARE:
            data.d1.at[i, "decision"] = DESCRIBE
            a = data.analysis(r["analysis"])
            if a is not None and a.id == PRIMARY:
                a.described[r["comparison"].split(" vs ")[0]] = "smallest group"


# ---------------------------------------------------------------------------
# S3: the wait from diagnosis to the entry event
# ---------------------------------------------------------------------------
def delays(data: AnalysisData, settings: Settings, by: list[str]) -> pd.DataFrame:
    """The wait from diagnosis to the entry event, in months, among the patients in the analysis: for all of them
    and by each trait in `by` ('genotype', 'year_band' for the model's year bands, or a patients field, coded as
    the adjusted model codes it when it is an adjustment variable). How many, and the quartiles (linear
    interpolation). A blank trait is 'Missing'. Exact: publish through the safety gate."""
    dpm = settings.read("time.days_per_year") / 12
    x = data.in_analysis().set_index("patient_id")
    cov = data.covariates.set_index("patient_id")
    banded = next((t.field for t in data.coded.terms if t.form == "bands"), None)
    rows = []

    def add(label, level_series):
        for level, d in x.groupby(level_series.reindex(x.index).fillna("Missing"), sort=True):
            q = d["report_day"].quantile([0.25, 0.5, 0.75]).to_numpy() / dpm
            rows.append({"by": label, "level": str(level), "n": len(d), "q1": q[0], "median": q[1], "q3": q[2]})

    add("All patients", pd.Series("All patients", index=x.index))
    for trait in by:
        if trait == "genotype":
            add("Genotype group", x["genotype"])
        elif trait == "year_band":
            if banded is None:
                raise RunStop("S3 is asked by year band, but no adjustment variable is in bands")
            add("Year band", cov[banded].astype(object))
        elif trait in cov.columns:
            add(trait.replace("_", " ").capitalize(), cov[trait].astype(object))
        elif trait in data.patients.columns:
            add(trait.replace("_", " ").capitalize(), data.patients[trait].astype(object))
        else:
            raise RunStop(f"S3 is asked by '{trait}', which the tables don't hold")
    return pd.DataFrame(rows, columns=["by", "level", "n", "q1", "median", "q3"])


def flag_counts(data: AnalysisData, settings: Settings) -> dict[str, int]:
    """The patients kept but counted, among the eligible: each flag phase 1's cohort sets, and patients matching
    more than one group's rules. Exact."""
    e = data.eligible
    out = {}
    for col in [c for c in e.columns if c.startswith("f:")]:
        out[col[2:]] = int(e[col].astype(bool).sum())
    groups = [g["name"] for g in settings.read("genotype.groups") if not g.get("none_of_the_above")]
    if len(groups) >= 2:
        out["overlapping_groups"] = int((e[[f"match:{g}" for g in groups]].astype(bool).sum(axis=1) >= 2).sum())
    return out


__all__ = ["build_analysis", "AnalysisData", "PlannedAnalysis", "Grouping", "OutcomeCheck", "RunStop", "d1_rows",
           "compare_with_plan", "keep_plan_decisions", "delays", "flag_counts", "treatment_clock", "first_line",
           "apply_time_check", "survival_after_entry", "comparison"]
