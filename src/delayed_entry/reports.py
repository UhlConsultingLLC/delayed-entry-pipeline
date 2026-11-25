"""The reports (design plan section 10): the results tables and figure, through the safety gate, and the run report.

Generalized from the CRC BPC project's crcbpc.outputs and crcbpc.secondary.public_s3. Each table is built from the
run's exact results (the R side's Parquet files and the analysis cohort), then passes the gate's rules:

    flow.csv           the cohort's steps, merged so no small count is shown or recoverable
    counts.csv         each analysis's group sizes, and the survival-after-entry count
    descriptives.csv   the patients kept but counted, and the distributions of the fields described
    risk_table.csv     patients at risk at each fixed time, by pass and group
    curves.csv         the Kaplan-Meier steps that may be drawn
    estimates.csv      every estimate: medians, RMST, survival at fixed times, HRs, the assumption checks
    bias.csv           every naive-minus-corrected contrast, with its interval, bar and label
    balance.csv        pass 3b's covariate balance
    s3_delays.csv      the wait from diagnosis to the entry event (when the settings run it)
    figure_km.svg      the naive and corrected curves by group, text kept as text
"""
from __future__ import annotations

import html
import math
from statistics import NormalDist

import duckdb
import pandas as pd

from delayed_entry.analysis import RUN, AnalysisData, comparison, delays, flag_counts
from delayed_entry.cohort import field_condition
from delayed_entry.covariates import humanize
from delayed_entry.gate import (NOT_SHOWN, CountSpec, Gate, GateError, Relation, add_note, balance_precision,
                                binary_smd, count_values, descriptive_values, drawn_steps, hide_across,
                                hide_survival_of_hidden_groups, public_curves, publish_counts, publish_descriptives,
                                require_safe, risk_table, show, small_group_rule, withheld_curves)
from delayed_entry.privacy import public_flow, suppress_small_cells
from delayed_entry.settings import Settings
from delayed_entry.tables import FIELDS

ALL = "All patients"
FLOAT_FORMAT = "%.6g"
SMD_DECIMALS_AFTER = 3          # balance after weighting: shown to three decimals (it can't be worked back to counts)
FLAG_LABELS = {
    "censored_at_age_limit": "Survival withheld by the source's age rule, so censored at the age limit",
    "overlapping_groups": "Matched more than one group's rules, so counted in the first by precedence",
    "overlapping_groups_different_genes": "Matched more than one group's rules through calls in different genes, so "
                                          "counted in the first by precedence",
    "overlapping_groups_same_gene": "Matched more than one group's rules through calls in one gene only, so counted "
                                    "in the first by precedence",
    "entry_before_diagnosis": "Entry event came before diagnosis",
    "first_test_linked_to_other_cancers": "First test is linked to more than one cancer",
    "unranked_first_test": "Another test has a withheld time, so it was ranked last",
    "first_test_tied": "First test tied with another after every tie-break, so the lowest test ID was taken",
    "after_follow_up_end_unknown": "First test has no record of whether it came back after follow-up ended",
}


def read_parquet(path) -> pd.DataFrame:
    return duckdb.execute("SELECT * FROM read_parquet(?)", [str(path)]).df()


def _exp(x):
    return None if x is None or pd.isna(x) else math.exp(x)


def is_true(v) -> bool:
    """A yes/no value from R, where a missing one (NA) is not yes."""
    return False if v is None or pd.isna(v) else bool(v)


def tau_verdict(row) -> str:
    """The note on a Kendall's tau row. With no interval (the estimator failed for the group, or gave no standard
    error) the row is neither flagged nor cleared: it says it wasn't assessed, and why when R said."""
    if row.flagged is None or pd.isna(row.flagged):
        why = getattr(row, "tau_note", None)
        return "not assessed: " + (why if isinstance(why, str) and why else "no interval")
    return "flagged" if row.flagged else "not flagged"


class Results:
    """The run's exact results and everything the tables need. Built once; each table method returns a table
    that has passed the gate's rules for it."""

    def __init__(self, data: AnalysisData, settings: Settings, run_dir, gate: Gate):
        self.d, self.s, self.run_dir, self.gate = data, settings, run_dir, gate
        self.dpy = settings.read("time.days_per_year")
        self.dpm = self.dpy / 12
        self.years = list(settings.read("summary.survival_at_years"))
        self.horizon = settings.read("time.horizon_years")
        self.level = settings.read("uncertainty.interval_level")
        self.z = NormalDist().inv_cdf((1 + self.level) / 2)
        self.groups = list(data.genotype.levels)
        self.ph = next((a for a in data.analyses if a.kind == "ph_check" and a.status == RUN), None)
        self.qi = next((a for a in data.analyses if a.kind == "quasi_independence" and a.status == RUN), None)
        self.s3 = next((a for a in data.analyses if a.kind == "descriptive_delay" and a.status == RUN), None)
        self.notes: list[str] = []
        self._counts()

    def r(self, name) -> pd.DataFrame:
        return read_parquet(self.run_dir / name)

    def has(self, name) -> bool:
        return (self.run_dir / name).exists()

    # ---- counts ----
    def count_spec(self) -> CountSpec:
        d = self.d
        primary = {g: int((d.in_analysis()["genotype"] == g).sum()) for g in self.groups}
        biomarkers, others, first = {}, [], {}
        for a in d.analyses:
            if a.status != RUN:
                continue
            if a.kind == "biomarker":
                spec = self.s.analysis(a.id)
                bc = d.biomarker_counts[a.id]
                biomarkers[a.id] = {"groups": {g: bc[g] for g in a.grouping.levels},
                                    "unknown": (spec["unknown_group"], bc[spec["unknown_group"]]),
                                    "not_in": (f"not in {a.id}", bc[f"not in {a.id}"])}
            elif a.kind in ("subset", "treatment_clock", "biomarker_variant"):
                rows = d.d1[d.d1["analysis"] == a.id]
                groups = {a.grouping.reference: int(rows["n_reference"].iloc[0])}
                groups |= {c.split(" vs ")[0]: int(n) for c, n in zip(rows["comparison"], rows["n_group"])}
                groups = {g: groups[g] for g in a.grouping.levels}
                base = a.parent if a.kind == "biomarker_variant" else "primary"
                others.append((a.id, base, groups))
                if a.kind == "biomarker_variant":
                    spec = self.s.analysis(a.id)
                    parent = self.s.analysis(a.parent)
                    if spec.get("timing", parent["timing"]) == "any_time" and parent["timing"] != "any_time":
                        first[a.id] = parent["groups"][0]["name"]
        return CountSpec(len(d.cohort), primary, biomarkers, others, first)

    def _descriptive_tables(self) -> tuple[list[tuple[str, dict]], list[tuple[str, list[str], dict]]]:
        """The distributions described: each field a subset analysis keeps patients by (with the condition's own
        normalization), then each extra field the adapter carries that the model doesn't use. Also, for each such
        subset, the levels it keeps (for the cross-table relations)."""
        p = self.d.patients
        tables, keeps, seen = [], [], set()
        for a in self.s.read("analyses"):
            keep = a.get("keep") if a["run"] and a["kind"] == "subset" else None
            if not keep or "field" not in keep:
                continue
            f = keep["field"]
            values = p[f]
            if keep.get("normalize") == "trim":
                values = values.astype("string").str.strip()
            levels = values.astype(object).where(values.notna() & (values.astype("string") != ""), None)
            counts = levels.fillna("Missing").value_counts().to_dict()
            counts = {str(k): int(v) for k, v in counts.items()}
            kept = [lv for lv in counts if bool(field_condition(pd.DataFrame({f: [None if lv == "Missing" else lv]}),
                                                                {**keep, "normalize": "none"}).iloc[0])]
            if f not in seen:
                tables.append((humanize(f), counts))
                seen.add(f)
            keeps.append((a["id"], [f"{humanize(f)} | {lv}" for lv in kept], counts))
        standard = {x.name for x in FIELDS if x.table == "patients"}
        model = {v["field"] for v in self.s.read("adjustment.variables")}
        for f in p.columns:
            if f in standard or f in model or f in seen or f == "outcome_recorded":
                continue
            values = p[f]
            if not (pd.api.types.is_string_dtype(values) or values.dtype == object):
                continue
            levels = values.astype("string").str.strip()
            levels = levels.mask(levels == "")
            counts = {str(k): int(v) for k, v in levels.fillna("Missing").value_counts().to_dict().items()}
            tables.append((humanize(f), counts))
            seen.add(f)
        return tables, keeps

    def _counts(self):
        g = self.gate
        self.spec = self.count_spec()
        cv, cr, self.count_order = count_values(self.spec)
        flags = flag_counts(self.d, self.s)
        if "overlapping_groups_same_gene" in flags:
            # Each part is masked on its own (a flag's own rule); the total is shown only when neither part would be
            # masked, so a masked part can't be worked out from the total and the other.
            n = len(self.d.cohort)
            parts = [flags["overlapping_groups_different_genes"], flags["overlapping_groups_same_gene"]]
            if any(g.small(x) or g.small(n - x) for x in parts) or g.small(n - flags["overlapping_groups"]):
                del flags["overlapping_groups"]
        flag_rows = [(self.flag_label(k), flags[k]) for k in FLAG_LABELS if k in flags]
        tables, keeps = self._descriptive_tables()
        dv, dr, self.desc_order = descriptive_values(len(self.d.cohort), flag_rows, tables)
        counts_shown, desc_shown = publish_counts(cv, cr, self.spec, g), publish_descriptives(dv, dr, g)
        cross = []
        for a, kept_keys, _ in keeps:
            parts = tuple(f"{a} | {grp}" for grp in self.groups)
            if not all(p_ in cv for p_ in parts) or not kept_keys:
                continue
            if len(kept_keys) == 1:
                total = kept_keys[0]
            else:
                total = f"{a} | kept levels"
                dv[total] = sum(dv[k] for k in kept_keys)
                dr.append(Relation(total, tuple(kept_keys)))
            closed = sum(cv[p_] for p_ in parts) == dv[total]
            cross.append(Relation(total, parts, closed=closed))
        values = cv | dv
        both = hide_across(values, counts_shown | desc_shown, cr + dr, cross, g)
        self.cv, self.dv = cv, dv
        self.counts_shown, self.desc_shown = both & set(cv), both & set(dv)
        self.n = {}
        for k, x in cv.items():
            if " | " in k:
                a, grp = k.split(" | ", 1)
                self.n[(a, grp)] = x
        for a in {a for a, _ in self.n}:
            self.n[(a, ALL)] = sum(self.n.get((a, grp), 0) for grp in self.analysis_groups(a))
        self.hidden_sizes = {tuple(k.split(" | ", 1)) for k in cv if " | " in k and k not in self.counts_shown}
        g.record("counts: small counts and subtraction", "pass",
                 f"{len(cv) - len(self.counts_shown & set(cv))} of {len(cv)} counts hidden or merged")
        g.record("descriptives: small counts and subtraction", "pass",
                 f"{len([k for k in dv if k not in self.desc_shown])} of {len(dv)} counts hidden")
        if cross:
            g.record("counts and descriptives judged together", "pass",
                     f"{len(cross)} relations between the two tables checked")

    def flag_label(self, key: str) -> str:
        label = FLAG_LABELS[key]
        if key == "entry_before_diagnosis":
            label += {"enter_at_day_0": ", so entry is day 0", "exclude": ", so the patient is excluded"}[
                self.s.read("entry.before_diagnosis")]
        return label

    def analysis_groups(self, a: str) -> list[str]:
        x = self.d.analysis(a)
        if x is not None and x.grouping is not None:
            return list(x.grouping.levels)
        return self.groups

    def counts_table(self) -> pd.DataFrame:
        hidden_extra = set()
        for b, d in self.spec.biomarkers.items():
            un, ni = f"{b} | {d['unknown'][0]}", f"{b} | {d['not_in'][0]}"
            merged = f"{b} | {d['unknown'][0]}, or {d['not_in'][0]}"
            hidden_extra |= {k for k in (un, ni, merged) if k not in self.counts_shown}
        rows = [{"table": k.split(" | ")[0] if " | " in k else "flow", "group": k.split(" | ", 1)[-1],
                 "patients": show(self.cv, k, self.counts_shown, self.gate)}
                for k in self.count_order if k not in hidden_extra]
        window = self.d.d2["window_days"]
        rows.append({"table": "D2", "group": f"eligible patients who died within {window:g} days of their entry event",
                     "patients": self.d.d2["reported_as"]})
        return pd.DataFrame(rows)

    def descriptives_table(self) -> pd.DataFrame:
        return pd.DataFrame([{"table": k.split(" | ")[0], "group": k.split(" | ", 1)[1],
                              "patients": show(self.dv, k, self.desc_shown, self.gate)} for k in self.desc_order])

    def flow_table(self) -> pd.DataFrame:
        return public_flow(self.d.flow[["label", "excluded", "remaining"]], min_count=self.gate.min_count)

    # ---- risk table and curves ----
    def passes(self) -> pd.DataFrame:
        return self.r("r_passes.parquet")

    def group_n(self) -> dict[str, int]:
        return {g: self.n[("primary", g)] for g in [ALL] + self.groups}

    def risk(self) -> pd.DataFrame:
        if not hasattr(self, "_risk"):
            p = self.passes()
            self._risk = risk_table(p[p["pass"].isin(["naive", "corrected"])], ALL, self.groups, self.group_n(), self.gate)
            self.gate.record("risk table: small counts and subtraction", "pass",
                             f"{int((self._risk['at_risk'].astype(str).isin([NOT_SHOWN, self.gate.small_text])).sum())} "
                             f"of {len(self._risk)} cells hidden")
        return self._risk

    def curves(self) -> tuple[pd.DataFrame, dict]:
        steps = drawn_steps(self.r("r_curves.parquet"), self.group_n(), max(self.years) * self.dpy, self.gate)
        withheld = withheld_curves(steps, self.risk(), self.years, self.dpy, self.gate)
        for g, why in withheld.items():
            self.notes.append(f"The curves of {g} are not published: their steps would give away {why}.")
        self.gate.record("curves: steps that would give away a hidden count", "pass" if not withheld else "withheld",
                         "; ".join(f"{g}: withheld" for g in withheld) or "every group's curves may be drawn")
        return public_curves(steps, withheld, self.dpm), withheld

    # ---- estimates ----
    def rmst_label(self) -> str:
        return f"{self.horizon:g}-year RMST, months"

    def estimates_table(self) -> pd.DataFrame:
        rows = []
        z = self.z
        rmst = [self.r("r_rmst.parquet").assign(analysis="primary")]
        if self.has("r_secondary_rmst.parquet"):
            rmst.append(self.r("r_secondary_rmst.parquet"))
        rmst = pd.concat(rmst, ignore_index=True)
        rmst_ci = {(x.analysis, x.pass_, x.quantity, x.group): (x.value, x.lower, x.upper, x.resamples_missing)
                   for x in rmst.rename(columns={"pass": "pass_"}).itertuples()}
        ci_label = f"HR (robust {self.level:.0%} CI)"

        def missing_note(k) -> str:
            m = rmst_ci.get(k, (None, None, None, 0))[3]
            return f"interval not computed: {int(m)} resamples gave no value" if m and m > 0 else ""

        def from_passes(t: pd.DataFrame, analysis: str):
            for (p, grp), d in t.groupby(["pass", "group"], sort=False):
                q = dict(zip(d["quantity"], d["value"]))
                if "median" in q:
                    rows.append((analysis, p, "median, months", grp, q["median"], q.get("median_lower"),
                                 q.get("median_upper"), ""))
                    k = (analysis, p, "rmst", grp)
                    _, lo, hi, _ = rmst_ci.get(k, (None, None, None, 0))
                    rows.append((analysis, p, self.rmst_label(), grp, q["rmst"], lo, hi, missing_note(k)))
                for s in d[d["quantity"] == "survival"].itertuples():
                    ok = s.n_risk >= self.gate.min_count
                    rows.append((analysis, p, f"survival at {float(s.years):g} years", grp, s.value if ok else None,
                                 None, None, "" if ok else f"fewer than {self.gate.min_count} at risk"))
                if "log_hr" in q:
                    b, se = q["log_hr"], q["robust_se"]
                    if "log_hr_lower" in q:          # pooled over imputations: Rubin's interval
                        rows.append((analysis, p, ci_label, grp, _exp(b), _exp(q["log_hr_lower"]),
                                     _exp(q["log_hr_upper"]), self.pooled_note))
                    else:
                        rows.append((analysis, p, ci_label, grp, _exp(b), _exp(b - z * se), _exp(b + z * se), ""))
                    k = (analysis, p, "rmst_difference", grp)
                    if k in rmst_ci:
                        v, lo, hi, _ = rmst_ci[k]
                        rows.append((analysis, p, f"{self.horizon:g}-year RMST difference, months", grp, v, lo, hi,
                                     missing_note(k)))

        passes = self.passes()
        from_passes(passes[passes["pass"].isin(["naive", "corrected"])], "primary")
        from_passes(self.r("r_adjusted.parquet"), "3a")
        names = {"median": "median, months", "rmst": self.rmst_label(), "hr": ci_label,
                 "rmst_difference": f"{self.horizon:g}-year RMST difference, months"}
        for w in self.r("r_weighted.parquet").rename(columns={"pass": "pass_"}).itertuples():
            note = w.note or ""
            m = getattr(w, "resamples_missing", None)
            if m is not None and not pd.isna(m) and m > 0:
                note = "; ".join(x for x in (note, f"interval not computed: {int(m)} resamples gave no value") if x)
            if self.imputing:
                note = "; ".join(x for x in (note, self.pooled_note if w.quantity == "hr" else self.averaged_note) if x)
            rows.append(("3b", w.pass_, names[w.quantity], w.group, w.value, w.lower, w.upper, note))
        if self.ph is not None:
            ph = self.r("r_ph.parquet")
            g = ph[ph["term"] == "GLOBAL"].iloc[0]
            months = self.s.read("ph_check.split_months")
            rows.append((self.ph.id, "corrected", "Schoenfeld test p (global)", ALL, g["p"], None, None,
                         f"rejected: HRs split at {months:g} months added" if g["rejected"] else "not rejected"))
            if g["rejected"]:
                for s in self.r("r_ph_split.parquet").itertuples():
                    rows.append((self.ph.id, f"corrected, {s.period}", ci_label, s.contrast, s.hr, s.lower, s.upper,
                                 f"split at {months:g} months"))
        if self.qi is not None:
            tau = self.r("r_tau.parquet")
            several = tau["estimator"].nunique() > 1
            for s in tau.itertuples():
                q = "conditional Kendall's tau" + (f" ({s.estimator})" if several else "")
                rows.append((self.qi.id, "corrected", q, s.group, s.tau, s.lower, s.upper, tau_verdict(s)))
            for s in tau.drop_duplicates("group").itertuples():
                if s.group != ALL:
                    rows.append((self.qi.id, "transformation model", "median, months", s.group, s.median_transformation,
                                 None, None, s.warning if isinstance(s.warning, str) else ""))
        if self.has("r_secondary_estimates.parquet"):
            for a, t in self.r("r_secondary_estimates.parquet").groupby("analysis", sort=False):
                from_passes(t, a)
        out = pd.DataFrame(rows, columns=["analysis", "pass", "quantity", "group", "estimate", "lower", "upper", "note"])
        out = small_group_rule(out, [self.patients_behind(a, grp) for a, grp in zip(out["analysis"], out["group"])],
                               self.gate)
        out = hide_survival_of_hidden_groups(out, self.hidden_sizes)
        hidden = int((out["note"] == "not shown: counts.csv hides this group's size, and this estimate could give it away").sum())
        self.gate.record("survival estimates that would pin a hidden group size", "pass",
                         f"{hidden} estimates left blank" if hidden else "none to hide")
        corrected = out["analysis"].isin(["primary", "3a", "3b"]) & out["pass"].astype(str).str.endswith("corrected")
        return self.flag_qi(out, corrected)

    def patients_behind(self, analysis: str, group: str) -> int:
        a = analysis.split()[0]
        a = a if (a, ALL) in self.n else "primary"
        if " vs " in group:
            x, y = group.split(" vs ", 1)
            return min(self.n.get((a, x), 0), self.n.get((a, y), 0))
        return self.n.get((a, group), self.n[(a, ALL)])

    def qi_flagged(self) -> set[str]:
        if self.qi is None:
            return set()
        t = self.r("r_tau.parquet")
        return set(t.loc[t["flagged"].map(is_true), "group"])

    def qi_note(self) -> str:
        return f"{self.qi.id}: entry and survival may be dependent (Kendall's tau interval excludes 0)"

    def flag_qi(self, t: pd.DataFrame, rows: pd.Series) -> pd.DataFrame:
        flagged = self.qi_flagged()
        if not flagged:
            return t

        def touches(grp):
            if " vs " in grp:
                x, y = grp.split(" vs ", 1)
                return x in flagged or y in flagged
            return grp in flagged
        return add_note(t, rows & t["group"].map(touches), self.qi_note())

    def bias_table(self) -> pd.DataFrame:
        parts = [("primary", "r_bias.parquet"), ("3a", "r_adjusted_bias.parquet"), ("3b", "r_weighted_bias.parquet"),
                 ("conditional on t*", "r_conditional.parquet")]
        frames = []
        for analysis, f in parts:
            t = self.r(f)
            t.insert(0, "analysis", analysis)
            frames.append(t)
        if self.has("r_secondary_bias.parquet"):
            frames.append(self.r("r_secondary_bias.parquet"))
        t = pd.concat(frames, ignore_index=True)
        t["note"] = t["note"].fillna("") if "note" in t else ""
        t = t[["analysis", "quantity", "group", "estimate", "lower", "upper", "resamples_missing", "bar", "label", "note"]]
        t = small_group_rule(t, [self.patients_behind(a, grp) for a, grp in zip(t["analysis"], t["group"])], self.gate)
        if self.imputing:
            t = add_note(t, t["analysis"].map(self.is_adjusted), self.averaged_note + "; each bootstrap resample imputed once")
        return self.flag_qi(t, t["analysis"].isin(["primary", "3a", "3b", "conditional on t*"]))

    # ---- multiple imputation ----
    @property
    def imputing(self) -> bool:
        return self.d.covariate_plan["missing_data_rule"] == "multiple imputation"

    @property
    def pooled_note(self) -> str:
        return f"pooled over {self.s.read('missing_data.imputations')} imputations by Rubin's rules" if self.imputing else ""

    @property
    def averaged_note(self) -> str:
        return f"averaged over {self.s.read('missing_data.imputations')} imputations" if self.imputing else ""

    def is_adjusted(self, analysis: str) -> bool:
        """An analysis whose models hold the adjustment variables: 3a, 3b, a biomarker's 3a, a stratified model."""
        strat = {a.id for a in self.d.analyses if a.kind == "adjusted_variant"}
        return analysis in ("3a", "3b") or analysis.endswith(" 3a") or analysis in strat

    # ---- balance ----
    def balance_term(self, term: str) -> str:
        for t in sorted(self.d.coded.terms, key=lambda x: -len(x.field)):
            f = t.field
            if term == f:
                return humanize(f)
            if term.startswith(f + "_curve"):
                k = len(t.knots_at_percentiles) - 2
                return humanize(f) + (", spline term" if k == 1 else f", spline term {term[len(f) + 6:]}")
            if term.startswith(f + "_"):
                name = "Year band" if t.form == "bands" else humanize(f)
                return f"{name}: {term[len(f) + 1:]}"
        return term

    def published_traits(self) -> list[str]:
        """Adjustment variables whose per-level totals another table publishes (S3's traits)."""
        if self.s3 is None:
            return []
        by = self.s.analysis(self.s3.id)["by"]
        fields = [t.field for t in self.d.coded.terms if t.form in ("categorical", "bands")]
        banded = [t.field for t in self.d.coded.terms if t.form == "bands"]
        out = [f for f in fields if f in by]
        if "year_band" in by:
            out += [f for f in banded if f not in out]
        return out

    def model_table(self, field: str, values: pd.DataFrame | None = None) -> pd.DataFrame:
        """Patients in the model by group and level of one variable (exact). The model's patients are those in the
        analysis with every adjustment variable; with imputation, `values` holds one imputation's completed values."""
        x = self.d.in_analysis()[["patient_id", "genotype"]].merge(
            self.d.covariates if values is None else values, on="patient_id")
        fields = [t.field for t in self.d.coded.terms]
        x = x[x[fields].notna().all(axis=1)]
        levels = next(t.levels for t in self.d.coded.terms if t.field == field)
        tab = pd.crosstab(x["genotype"], x[field].astype(str)).reindex(index=self.groups, columns=list(levels), fill_value=0)
        return tab.loc[:, (tab.sum(axis=0) > 0)]

    def trait_tables(self, published: list[str], smd_before: pd.DataFrame, values: pd.DataFrame | None) -> dict:
        """For each trait whose totals are published: its group-by-level table and each level's SMD before weighting,
        checked against the gate's own arithmetic."""
        traits = {}
        for f in published:
            tab = self.model_table(f, values)
            smds = {}
            for lv in tab.columns:
                term = smd_before.loc[smd_before["term"] == f"{f}_{lv}", "smd_before"]
                if len(term) == 0 and tab.shape[1] == 2:          # a two-level factor shows one term only
                    other = [x for x in tab.columns if x != lv][0]
                    term = smd_before.loc[smd_before["term"] == f"{f}_{other}", "smd_before"]
                if len(term) == 0:
                    raise GateError(f"balance: no term for {f} = {lv}")
                smds[lv] = float(term.iloc[0])
                mine = binary_smd([int(v) for v in tab[lv]], [int(v) for v in tab.sum(axis=1)])
                if abs(mine - smds[lv]) > 1e-6:
                    raise GateError(f"balance: the gate's arithmetic doesn't reproduce the SMD of {f} = {lv}, so it "
                                    "can't judge its precision")
            traits[humanize(f)] = (tab, smds)
        return traits

    def balance_table(self) -> pd.DataFrame:
        b = self.r("r_balance.parquet")
        target = self.s.read("weighting.balance_target_smd")
        published = self.published_traits()
        if not self.imputing:
            traits = self.trait_tables(published, b, None)
            decimals = balance_precision(traits, SMD_DECIMALS_AFTER, self.gate)
        else:
            # The published SMDs are averages across imputations. Judged on each imputation's own table and own SMDs,
            # as if it were published alone, and the fewest decimals any needs are used.
            by_imp, completed = self.r("r_balance_by_imputation.parquet"), self.r("r_imputed.parquet")
            decimals, traits = SMD_DECIMALS_AFTER, {}
            for i, smd in by_imp.groupby("imputation"):
                values = completed[completed["imputation"] == i].drop(columns="imputation")
                traits = self.trait_tables(published, smd, values)
                decimals = min(decimals, balance_precision(traits, decimals, self.gate))
        self.gate.record("balance: precision that can't rebuild a small cell", "pass",
                         f"before weighting to {decimals} decimals" + (f" (traits whose totals are published: "
                                                                      f"{', '.join(traits)})" if traits else "")
                         + ("; judged on each imputation" if self.imputing else ""))
        after = b["smd_after"].astype(float)
        totals = b["term"].map(lambda t: any(t.startswith(f + "_") for f in published))
        if decimals < SMD_DECIMALS_AFTER:
            self.notes.append(f"Balance before weighting is shown to {decimals} decimals: more would let a small "
                              "group-by-trait count be worked out from published totals.")
        return pd.DataFrame({"term": [self.balance_term(t) for t in b["term"]],
                             "largest_smd_before": b["smd_before"].astype(float).round(decimals),
                             "largest_smd_after": after.round(SMD_DECIMALS_AFTER),
                             f"after_below_{target:g}": after < target,
                             "largest_pair_before": b["pair_before"].where(~totals, NOT_SHOWN),
                             "largest_pair_after": b["pair_after"]})

    # ---- S3 ----
    def s3_table(self) -> pd.DataFrame | None:
        if self.s3 is None:
            return None
        t = delays(self.d, self.s, self.s.analysis(self.s3.id)["by"]).round({"q1": 1, "median": 1, "q3": 1})
        t = suppress_small_cells(t, ["n"], ["q1", "median", "q3"], min_count=self.gate.min_count)
        hidden = t[t["n"] == self.gate.small_text].groupby("by").size()
        lone = sorted(hidden[hidden == 1].index)
        if lone:
            raise GateError(f"{self.s3.id}: one hidden count in {lone} could be worked out from the total; merge its "
                            "groups first")
        self.gate.record(f"{self.s3.id}: small counts", "pass", f"{len(hidden)} hidden")
        return t


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
COLORS = {"naive": "#9A5A00", "corrected": "#0C6B77", "ink": "#15202B", "muted": "#56626E"}


def figure_km(curves: pd.DataFrame, risk: pd.DataFrame, groups: list[str], years: list, withheld, gate: Gate) -> str:
    """The naive and corrected curves for each group, side by side, with the published numbers at risk underneath.
    A group in `withheld` keeps its panel and numbers at risk, with a note where its curves would be. Plain SVG:
    every label is text that can be edited."""
    last = max(years) * 12
    PW, PH, L, T, GAP = 260, 200, 80, 50, 50
    W = L + len(groups) * (PW + GAP) + 20
    H = 400
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="Arial, Helvetica, sans-serif" font-size="12">',
           "<title>Overall survival by group: naive and corrected</title>",
           f"<desc>Kaplan-Meier curves from diagnosis, drawn only where at least {gate.min_count} patients are at risk, "
           "and only for groups whose steps can't give away a count the numbers at risk hide.</desc>",
           f'<rect width="{W}" height="{H}" fill="#FFFFFF"/>']
    ticks = [12 * y for y in range(0, int(math.floor(last / 12)) + 1)]
    for i, g in enumerate(groups):
        x0 = L + i * (PW + GAP)

        def px(m):
            return x0 + m / last * PW

        def py(s):
            return T + (1 - s) * PH
        out.append(f'<text x="{x0}" y="{T - 14}" font-weight="bold" fill="{COLORS["ink"]}">{html.escape(g)}</text>')
        out.append(f'<line x1="{x0}" y1="{T + PH}" x2="{x0 + PW}" y2="{T + PH}" stroke="{COLORS["muted"]}"/>')
        out.append(f'<line x1="{x0}" y1="{T}" x2="{x0}" y2="{T + PH}" stroke="{COLORS["muted"]}"/>')
        for m in ticks:
            out.append(f'<text x="{px(m):.1f}" y="{T + PH + 16}" text-anchor="middle" fill="{COLORS["muted"]}">{m:g}</text>')
        for s in (0, 0.5, 1):
            out.append(f'<text x="{x0 - 6}" y="{py(s) + 4:.1f}" text-anchor="end" fill="{COLORS["muted"]}">{s:g}</text>')
        for p in ("naive", "corrected"):
            d = curves[(curves["pass"] == p) & (curves["group"] == g)].sort_values("months")
            if d.empty:
                continue
            path = f"M{px(d['months'].iloc[0]):.1f},{py(d['survival'].iloc[0]):.1f}"
            for m, s in zip(d["months"].iloc[1:], d["survival"].iloc[1:]):
                path += f"H{px(m):.1f}V{py(s):.1f}"
            out.append(f'<path d="{path}" fill="none" stroke="{COLORS[p]}" stroke-width="2"/>')
        if g in withheld:
            lines = ("Curves not drawn: their steps", f"would let a count below {gate.min_count}", "be worked out")
            for k, line in enumerate(lines):
                out.append(f'<text x="{x0 + PW / 2:.1f}" y="{T + PH / 2 - 12 + 16 * k:.1f}" text-anchor="middle" '
                           f'fill="{COLORS["muted"]}">{html.escape(line)}</text>')
        for j, p in enumerate(("naive", "corrected")):
            y = T + PH + 40 + j * 16
            if i == 0:
                out.append(f'<text x="{x0 - 10}" y="{y}" text-anchor="end" fill="{COLORS[p]}">{p}</text>')
            for yr in years:
                cell = risk[(risk["pass"] == p) & (risk["group"] == g) & (risk["years"].astype(float) == float(yr))]["at_risk"]
                txt = "–" if cell.empty or cell.iloc[0] == NOT_SHOWN else str(cell.iloc[0])
                out.append(f'<text x="{px(12 * yr):.1f}" y="{y}" text-anchor="middle" fill="{COLORS[p]}">{html.escape(txt)}</text>')
    months = ", ".join(f"{12 * y:g}" for y in years)
    note = (f'Months from diagnosis. Numbers at risk at {months} months; "<{gate.min_count}" is fewer than {gate.min_count}, '
            f'"–" is hidden so that a smaller count can\'t be worked out.')
    out.append(f'<text x="{L}" y="{H - 40}" fill="{COLORS["muted"]}">{html.escape(note, quote=False)}</text>')
    out.append(f'<text x="{L}" y="{H - 22}"><tspan fill="{COLORS["naive"]}">naive: from diagnosis, ignoring the wait for the entry event</tspan>'
               f'<tspan dx="24" fill="{COLORS["corrected"]}">corrected: at risk only from the entry event</tspan></text>')
    out.append("</svg>")
    return "\n".join(out) + "\n"


def build_tables(res: Results) -> tuple[dict[str, pd.DataFrame], str, dict]:
    """Every results table and the figure, each through the gate's rules. Returns the tables, the figure and the
    curves withheld."""
    curves, withheld = res.curves()
    tables = {
        "flow.csv": res.flow_table(),
        "counts.csv": res.counts_table(),
        "descriptives.csv": res.descriptives_table(),
        "risk_table.csv": res.risk(),
        "curves.csv": curves,
        "estimates.csv": res.estimates_table(),
        "bias.csv": res.bias_table(),
        "balance.csv": res.balance_table(),
    }
    s3 = res.s3_table()
    if s3 is not None:
        tables["s3_delays.csv"] = s3
    res.gate.record("flow: steps merged so no small count is shown or recoverable", "pass", "")
    figure = figure_km(curves, res.risk(), res.groups, res.years, withheld, res.gate)
    return tables, figure, withheld


__all__ = ["Results", "build_tables", "figure_km", "read_parquet", "FLOAT_FORMAT", "require_safe", "comparison"]
