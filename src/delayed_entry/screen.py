"""The readiness screen (design plan section 6): what the data can support, before any outcome is seen.

    from delayed_entry.tables import outcome_blind
    from delayed_entry.screen import screen
    report = screen(outcome_blind(tables), settings, profile)
    print(report.to_markdown())

The screen is given only the outcome-blind tables: it sees whether each patient's outcome is recorded,
never what it is. A test enforces this: scrambling every outcome value leaves the report unchanged.

It reports three lists:
    ready                 analyses that can run as standard
    missing               what is absent, and which analyses that blocks
    deviations required   each with the rule that triggered it, why, and its effect; severity block (that
                          analysis can't run), deviation (it runs differently, stated in the plan) or note
plus warnings on choices that look implausible for the data (the user decides; the warning stays), the
checks that need outcomes and so run only after the plan is locked, and every setting the screen read,
with where it came from (defaults are listed separately).

Every threshold comes from the settings. Counts in the report are aggregate, and any count from 1 to
the privacy minimum is shown as "<minimum".
"""
from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field as dataclass_field
from pathlib import Path

import pandas as pd
import yaml

from delayed_entry.cohort import (Cohort, MissingFieldError, build_cohort, drug_matches, field_condition,
                                  result_matches)
from delayed_entry.covariates import code_adjustment, missing_data_rule
from delayed_entry.hashing import sha256_file
from delayed_entry.privacy import is_small, public_flow, show
from delayed_entry.profile import CapabilityProfile
from delayed_entry.settings import Settings
from delayed_entry.tables import (ENTRY_EVENTS, FIELDS, OPTIONAL_TABLES, OUTCOME_FIELDS, REQUIRED, RESULT_KINDS,
                                  TABLES, OutcomeBlindTables, StandardTables)

READY, MISSING, DEVIATIONS = "ready", "missing", "deviations required"
BLOCK, DEVIATION, NOTE = "block", "deviation", "note"
COMPARE, DESCRIBE, NOT_RUN, BLOCKED = "compare", "descriptive only", "not run", "blocked"
PRIMARY = "primary"
# A biomarker's fallback by mutation burden compares these two groups; the low one is the reference.
BURDEN_HIGH, BURDEN_LOW = "high mutation burden", "low mutation burden"


def _sentence(text: str) -> str:
    """Text ending in a period, for the Markdown report."""
    text = text.strip()
    return text if not text or text.endswith((".", "!", "?")) else text + "."


@dataclass(frozen=True)
class Item:
    list: str              # READY, MISSING or DEVIATIONS
    severity: str          # BLOCK, DEVIATION or NOTE ("" on the ready list)
    analysis: str          # "primary", an analysis id, or "all"
    topic: str             # the rule or check that found it
    finding: str           # what was found
    effect: str            # what it changes
    settings: tuple = ()   # the choices it read

    def to_dict(self) -> dict:
        return {"severity": self.severity, "analysis": self.analysis, "topic": self.topic, "finding": self.finding,
                "effect": self.effect, "settings": list(self.settings)}


@dataclass
class ReadinessReport:
    study: str
    settings_sha256: str
    profile: dict
    items: list[Item]
    warnings: list[Item]
    decisions: dict                  # analysis: {comparison: decision}
    group_sizes: dict                # analysis: {group: shown count}
    data_checks: list[dict]
    missingness: dict
    flow: list[dict] | None
    deferred: list[dict]
    settings_used: list[dict]
    defaults_used: list[dict]
    min_count: int
    settings_resolved_sha256: str = ""
    inputs: list = dataclass_field(default_factory=list)      # [{"file": name, "sha256": ...}]

    def list(self, name: str) -> list[Item]:
        return [i for i in self.items if i.list == name]

    def find(self, analysis: str | None = None, topic: str | None = None, severity: str | None = None) -> list[Item]:
        return [i for i in self.items + self.warnings
                if (analysis is None or i.analysis == analysis) and (topic is None or i.topic == topic)
                and (severity is None or i.severity == severity)]

    def to_dict(self) -> dict:
        return {
            "study": self.study,
            "settings_sha256": self.settings_sha256,
            "settings_resolved_sha256": self.settings_resolved_sha256,
            "inputs": self.inputs,
            "capability_profile": self.profile,
            "ready": [i.to_dict() for i in self.list(READY)],
            "missing": [i.to_dict() for i in self.list(MISSING)],
            "deviations_required": [i.to_dict() for i in self.list(DEVIATIONS)],
            "warnings": [i.to_dict() for i in self.warnings],
            "decisions": self.decisions,
            "group_sizes": self.group_sizes,
            "flow": self.flow,
            "data_checks": self.data_checks,
            "missingness": self.missingness,
            "checks_after_the_lock": self.deferred,
            "settings_used": self.settings_used,
            "defaults_used": self.defaults_used,
            "privacy_min_count": self.min_count,
        }

    def to_yaml(self) -> str:
        return yaml.safe_dump(self.to_dict(), sort_keys=False, allow_unicode=True, width=120)

    def to_markdown(self) -> str:
        out = [f"# Readiness report: {self.study}", "",
               "Made from the outcome-blind tables: no outcome value was read. Counts from 1 to "
               f"{self.min_count - 1} are shown as \"<{self.min_count}\".", ""]
        for name in (READY, MISSING, DEVIATIONS):
            items = self.list(name)
            out += [f"## {name.capitalize()} ({len(items)})", ""]
            if not items:
                out += ["None.", ""]
            for i in items:
                sev = f"**{i.severity}** " if i.severity else ""
                out.append(f"- {sev}[{i.analysis}] {i.topic}: {_sentence(i.finding)}"
                           + (f" Effect: {_sentence(i.effect)}" if i.effect else ""))
            out.append("")
        out += [f"## Warnings ({len(self.warnings)})", ""]
        out += [f"- [{i.analysis}] {i.topic}: {_sentence(i.finding)} {_sentence(i.effect)}".rstrip()
                for i in self.warnings] or ["None."]
        out += ["", "## Decisions", "", "| Analysis | Comparison | Decision |", "|---|---|---|"]
        for a, comps in self.decisions.items():
            for comp, d in comps.items():
                out.append(f"| {a} | {comp} | {d} |")
        out += ["", "## Group sizes", "", "| Analysis | Group | Patients |", "|---|---|---|"]
        for a, groups in self.group_sizes.items():
            for g, n in groups.items():
                out.append(f"| {a} | {g} | {n} |")
        if self.flow:
            out += ["", "## Flow (steps merged so no small count is shown or recoverable)", "",
                    "| Step | Excluded | Remaining |", "|---|---|---|"]
            out += [f"| {r['label']} | {r['excluded']} | {r['remaining']} |" for r in self.flow]
        out += ["", "## Data checks", "", "| Check | Status | Detail |", "|---|---|---|"]
        out += [f"| {c['check']} | {c['status']} | {c['detail']} |" for c in self.data_checks]
        out += ["", "## Missing values", "", "| Field | Missing |", "|---|---|"]
        out += [f"| {k} | {v} |" for k, v in self.missingness.items()]
        out += ["", "## Checks that need outcomes (run after the plan is locked)", ""]
        out += [f"- {d['check']}: {d['rule']}" for d in self.deferred]
        out += ["", "## Settings the screen read", "", "| Choice | Value | Made by | Source | Where from |", "|---|---|---|---|---|"]
        for s in self.settings_used:
            v = yaml.safe_dump(s["value"], default_flow_style=True, width=10**6).replace("\n...", "").strip()
            if len(v) > 120:
                v = v[:117] + "..."
            out.append(f"| {s['choice']} | `{v}` | {s['made_by']} | {s['source']} | {s['origin']} |")
        out += ["", f"## Defaults applied ({len(self.defaults_used)})", ""]
        out += [f"- {d['choice']} = {d['value']} (source: {d['source']})" for d in self.defaults_used] or ["None."]
        return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
def biomarker_patients(blind: OutcomeBlindTables, c: Cohort, a: dict) -> tuple[pd.DataFrame, list[str]]:
    """The eligible patients a biomarker analysis `a` keeps: those a["exclude_if"] names are left out. Returns
    them, and the flags the tables lack."""
    e = c.eligible
    p = blind.patients.drop_duplicates("patient_id").set_index("patient_id").loc[e.index]
    keep = pd.Series(True, index=e.index)
    absent = []
    for flag in a["exclude_if"]:
        if flag not in p.columns:
            absent.append(flag)
            continue
        keep &= ~p[flag].astype("boolean").fillna(False).astype(bool)
    return e[keep], absent


def panel_sizes_given(blind: OutcomeBlindTables, profile: CapabilityProfile) -> bool:
    """Whether a biomarker's fallback takes its 'when panel sizes are given' form: the capability profile says the
    source gives them, and the panel_sizes table holds at least one."""
    ps = blind.panel_sizes
    return (profile.panel_sizes == "given" and ps is not None and "coding_megabases" in ps.columns
            and bool(ps["coding_megabases"].notna().any()))


def burden_problems(blind: OutcomeBlindTables) -> list[str]:
    """Why mutation burden can't be counted from these tables (labels only), or nothing."""
    out = []
    if "variant_class" not in blind.results.columns:
        out.append("results.variant_class is absent, so no mutation can be counted by its class")
    if "panel_id" not in blind.tests.columns:
        out.append("tests.panel_id is absent, so no entry test has a panel size")
    ps = blind.panel_sizes
    if ps is None or not {"panel_id", "coding_megabases"} <= set(ps.columns):
        out.append("the panel_sizes table, with panel_id and coding_megabases, is absent")
        return out
    sized = ps[ps["coding_megabases"].notna()]
    if (sized["coding_megabases"] <= 0).any():
        out.append("a panel's coding size is zero or less")
    if (sized.groupby("panel_id")["coding_megabases"].nunique() > 1).any():
        out.append("a panel is given two different coding sizes")
    return out


def mutation_burden(blind: OutcomeBlindTables, c: Cohort, a: dict) -> tuple[pd.Series, pd.Series]:
    """Each patient's mutation-burden group, for a biomarker analysis `a` that falls back to it (its settings'
    `mutation_burden`), among the patients the analysis keeps (a["exclude_if"] applies, as for its results).

    The burden is the number of the entry test's mutations whose variant class the settings list (each result row
    once, as the source lists it), over the coding size of the entry test's panel in megabases. High: at or above
    the settings' threshold, judged after rounding to 9 decimals so that a burden exactly at the threshold isn't
    lost to floating point. Low: below it. The analysis's unknown group: the entry test's panel has no size.
    Returns (group, burden per megabase), indexed by patient_id."""
    e, _ = biomarker_patients(blind, c, a)
    mb = a["mutation_burden"]
    res = blind.results
    is_mutation = (res["kind"] == "mutation").fillna(False).astype(bool)
    counted = res["variant_class"].astype("string").str.strip().isin(set(mb["variant_classes"]))
    counted = counted.fillna(False).astype(bool) & is_mutation
    per_test = res.loc[counted.to_numpy(), "test_id"].value_counts()
    mutations = e["first_test_id"].map(per_test).fillna(0).astype(float)
    ps = blind.panel_sizes
    sizes = ps[ps["coding_megabases"].notna()].drop_duplicates("panel_id").set_index("panel_id")["coding_megabases"]
    size = e["first_panel_id"].map(sizes).astype(float)
    burden = (mutations / size).round(9)
    group = pd.Series(a["unknown_group"], index=e.index, dtype="object")
    known = burden.notna()
    high = known & (burden >= mb["threshold_per_megabase"])
    group[high] = BURDEN_HIGH
    group[known & ~high] = BURDEN_LOW
    return group, burden


def biomarker_status(blind: OutcomeBlindTables, c: Cohort, a: dict, timing: str,
                     drop_discordant: bool = False) -> tuple[pd.Series, pd.Series, list[str]]:
    """Each eligible patient's group in a biomarker analysis `a` (a settings entry of kind biomarker), from
    their results of a["results_kind"]: the first group (in the settings' order) any result's value
    belongs to, else the unknown group. With timing "by_entry_test", only results dated on or before the
    entry test count; a result with no date isn't known to, so it doesn't. Patients a["exclude_if"] names
    are left out. Returns (status, discordant, flags the tables lack)."""
    e, absent = biomarker_patients(blind, c, a)
    tests = blind.tests
    # When the source dates events from birth, compare on that clock, as the source's own dates are (the CRC BPC
    # rule); otherwise on days from diagnosis. A blank on either side is not known to come first.
    clock, entry_col = (("days_from_birth", "entry_from_birth") if "days_from_birth" in tests.columns
                        else ("days_to_entry_event", "entry_day"))
    res = blind.results[blind.results["kind"] == a["results_kind"]]
    res = res.merge(tests[[c for c in ("test_id", "patient_id", clock) if c in tests.columns]], on="test_id", how="inner")
    res = res[res["patient_id"].isin(e.index)]
    if timing == "by_entry_test":
        if clock not in res.columns:
            res = res.iloc[0:0]
        else:
            entry = res["patient_id"].map(e[entry_col])
            res = res[(res[clock] <= entry).fillna(False).astype(bool)]
    matched = pd.DataFrame(index=e.index)
    for g in a["groups"]:
        ids = set(res.loc[res["value"].isin(g["values"]).fillna(False), "patient_id"])
        matched[g["name"]] = matched.index.isin(list(ids))
    status = pd.Series(a["unknown_group"], index=e.index, dtype="object")
    for g in reversed(a["groups"]):
        status[matched[g["name"]]] = g["name"]
    discordant = matched.sum(axis=1) >= 2
    if drop_discordant:
        status = status[~discordant]
    return status, discordant, absent


class _Screen:
    def __init__(self, blind: OutcomeBlindTables, settings: Settings, profile: CapabilityProfile,
                 source_steps: list[tuple[str, int]]):
        self.b, self.s, self.p = blind, settings, profile
        self.source_steps = source_steps
        self.items: list[Item] = []
        self.warnings: list[Item] = []
        self.decisions: dict = {}
        self.sizes: dict = {}
        self.checks: list[dict] = []
        # Counts are masked with the larger of the settings' minimum and the data-use agreement's: a setting
        # below the agreement is blocked, and never weakens the masking.
        self.settings_min_count = settings.read("privacy.min_count")
        self.min_count = max(self.settings_min_count, profile.privacy_min_count or 0)
        self.min_group = settings.read("smallest_group.patients")
        self.blocked_analyses: set[str] = set()
        self.banded = False             # True once primary group sizes are shown as ranges only
        self.missingness: dict = {}
        self.entry_times = profile.entry_event != "none"

    # ---- helpers ----
    def add(self, list_, severity, analysis, topic, finding, effect="", settings=()):
        item = Item(list_, severity, analysis, topic, finding, effect, tuple(settings))
        if item not in self.items:
            self.items.append(item)
        if severity == BLOCK:
            self.blocked_analyses.add(analysis)

    def warn(self, analysis, topic, finding, effect, settings=()):
        self.warnings.append(Item(DEVIATIONS, NOTE, analysis, topic, finding, effect, tuple(settings)))

    def n(self, x: int) -> str:
        return show(x, self.min_count)

    def band(self, k: int) -> str:
        """A group size as a range: enough to see the smallest-group decision, too coarse to subtract."""
        if k == 0:
            return "0"
        if is_small(k, self.min_count):
            return f"<{self.min_count}"
        if k < self.min_group:
            return f"{self.min_count} to {self.min_group - 1}"
        return f"at least {max(self.min_group, self.min_count)}"

    def safe_sizes(self, counts: dict, total: int | None = None, relative: dict | None = None,
                   primary: bool = False) -> dict:
        """Group sizes as the report may show them. A small group is masked, and then every size in the
        analysis is shown as a range, so the masked one can't be worked out from a total shown elsewhere.
        So too when the groups fall short of a known total (`total`) by a small number, and, once the
        primary sizes are ranges, everywhere. A group whose count differs from a reference count
        (`relative`, for example the primary analysis's) by a small number is shown as a range."""
        hidden_total = total is not None and is_small(total - sum(counts.values()), self.min_count)
        if self.banded or hidden_total or any(is_small(k, self.min_count) for k in counts.values()):
            if primary:
                self.banded = True
            return {g: self.band(k) for g, k in counts.items()}
        shown = {}
        for g, k in counts.items():
            if relative is not None and g in relative and is_small(abs(k - relative[g]), self.min_count):
                shown[g] = f"{self.band(k)} (differs from its reference by fewer than {self.min_count})"
            else:
                shown[g] = str(k)
        return shown

    def check(self, name, ok, detail=""):
        self.checks.append({"check": name, "status": "pass" if ok is True else ("FAIL" if ok is False else ok),
                            "detail": detail})

    # ---- structure ----
    def structure(self):
        for t in TABLES:
            df = self.b.table(t)
            if df is None:
                if t not in OPTIONAL_TABLES:
                    self.add(MISSING, BLOCK, "all", "standard tables", f"the {t} table is absent",
                             "nothing can be screened without it")
                continue
            for f in [f for f in FIELDS if f.table == t and f.key]:
                if f.name not in df.columns:
                    self.add(MISSING, BLOCK, "all", "standard tables", f"{t}.{f.name} is absent",
                             f"the {t} table can't be used")
        p, t, r = self.b.patients, self.b.tests, self.b.results
        if self.keys_missing():
            return
        dup = int(p["patient_id"].duplicated().sum())
        self.check("One row per patient", dup == 0, f"{self.n(dup)} repeated IDs" if dup else "")
        if dup:
            self.add(DEVIATIONS, BLOCK, "all", "duplicates", f"{self.n(dup)} patient IDs appear more than once",
                     "fix the source or the mapping before going on")
        dup = int(t["test_id"].duplicated().sum())
        self.check("One row per test", dup == 0, f"{self.n(dup)} repeated test IDs" if dup else "")
        if dup:
            self.add(DEVIATIONS, BLOCK, "all", "duplicates", f"{self.n(dup)} test IDs appear more than once",
                     "results can't be linked to one test; fix the source or the mapping")
        orphan = int((~t["patient_id"].isin(p["patient_id"])).sum())
        self.check("Every test belongs to a patient in the patients table", orphan == 0 or "note",
                   f"{self.n(orphan)} tests: not used" if orphan else "")
        orphan = int((~r["test_id"].isin(t["test_id"])).sum())
        self.check("Every result belongs to a test in the tests table", orphan == 0 or "note",
                   f"{self.n(orphan)} results: not used" if orphan else "")
        bad = sorted(set(t["entry_event"].dropna()) - set(ENTRY_EVENTS)) if "entry_event" in t else []
        self.check("Entry events are report, sample_collection or order", not bad, f"{len(bad)} other labels" if bad else "")
        bad = sorted(set(r["kind"].dropna()) - set(RESULT_KINDS)) if "kind" in r else []
        self.check("Result kinds are known", not bad, f"{len(bad)} other kinds" if bad else "")
        if "age_at_diagnosis" in p:
            age = pd.to_numeric(p["age_at_diagnosis"], errors="coerce")
            neg = int((age < 0).sum())
            self.check("No negative age", neg == 0, f"{self.n(neg)} patients" if neg else "")
            if self.p.age_top_code is not None:
                over = int((age > self.p.age_top_code).sum())
                self.check(f"No age above the source's top code ({self.p.age_top_code:g})", over == 0 or "note",
                           f"{self.n(over)} patients" if over else "")
        # Missing values of every field present: the count only, masked. No table total and no share, which
        # together with the flow diagram's counts could give away a step the flow merged to hide it.
        for f in FIELDS:
            df = self.b.table(f.table)
            if df is None or f.name not in df.columns or f.outcome or f.key:
                continue
            self.missingness[f"{f.table}.{f.name}"] = self.missing_text(int(df[f.name].isna().sum()), len(df))
        if self.b.outcome_fields_present:
            rec = self.b.patients["outcome_recorded"].astype("boolean").fillna(False)
            self.missingness["patients: outcome (vital status or time) not recorded"] = \
                self.missing_text(int((~rec).sum()), len(rec))

    def missing_text(self, missing: int, rows: int) -> str:
        present = rows - missing
        if present == 0:
            return "all"
        if is_small(present, self.min_count):
            return f"all but fewer than {self.min_count}"
        return self.n(missing)

    def keys_missing(self) -> bool:
        """True if a table or identity field is absent, so nothing past the structure can be screened."""
        return any(i.topic == "standard tables" and i.severity == BLOCK for i in self.items)

    # ---- what the analyses need ----
    def needs(self):
        """Required fields per component in use, against what the tables hold."""
        components = {"survival": PRIMARY, "corrected_clock": PRIMARY, "genotype": PRIMARY}
        if self.s.read("genotype.coverage.required"):
            components["wild_type_confirmation"] = PRIMARY
        for a in self.s.read("analyses"):
            if a["run"] and a["kind"] == "treatment_clock":
                components["treatment_clock"] = a["id"]
            if a["run"] and a["kind"] == "biomarker":
                components["biomarker"] = a["id"]
        for comp, analysis in components.items():
            for f in FIELDS:
                if f.used_by.get(comp) != REQUIRED:
                    continue
                if f.outcome:
                    if f.name not in self.b.outcome_fields_present:
                        self.add(MISSING, BLOCK, "all", "outcome", f"patients.{f.name} is absent",
                                 "no survival analysis can run")
                    continue
                df = self.b.table(f.table)
                if df is None or f.name not in df.columns:
                    sev = BLOCK
                    if comp == "wild_type_confirmation" and self.p.panel_coverage == "none":
                        continue        # reported once, from the profile
                    self.add(MISSING, sev, analysis, f"needed for {comp.replace('_', ' ')}",
                             f"{f.table}.{f.name} is absent", f"{analysis} can't run as planned")
        # fields the settings name
        named = []
        for r in self.s.read("cohort.rules"):
            named += [("cohort", PRIMARY, c["field"]) for c in r["all_of"] if "field" in c]
        for v in self.s.read("adjustment.variables"):
            named.append(("adjustment", PRIMARY, v["field"]))
        for a in self.s.read("analyses"):
            if a["run"] and a["kind"] == "subset" and "field" in a["keep"]:
                named.append(("subset", a["id"], a["keep"]["field"]))
            if a["run"] and a["kind"] == "adjusted_variant":
                named.append(("adjustment", a["id"], a["stratify_by"]))
        for comp, analysis, name in named:
            if name not in self.b.patients.columns:
                where = "cohort rule" if comp == "cohort" else ("adjusted models" if comp == "adjustment" else "subset")
                self.add(MISSING, BLOCK, analysis, f"needed for the {where}", f"patients.{name} is absent",
                         f"{analysis} can't run as planned")
        uses_flag = any(c.get("check") == "entry_before_follow_up_end" for r in self.s.read("cohort.rules") for c in r["all_of"])
        if uses_flag and "after_follow_up_end" not in self.b.tests.columns:
            self.add(MISSING, NOTE, PRIMARY, "time at risk",
                     "tests.after_follow_up_end is absent: the source doesn't say whether a test came back after "
                     "follow-up ended",
                     "the screen keeps every such patient; the run removes those with no time at risk, so the "
                     "screen's group sizes are upper bounds and the smallest-group rule is checked again at the run")

    # ---- the capability profile against the settings ----
    def profile_rules(self):
        want = self.s.read("entry.event")
        subs = self.s.read("entry.acceptable_substitutes")
        have = self.p.entry_event
        if have == want:
            self.add(READY, "", PRIMARY, "entry event", f"the source gives the {want} time", "the corrected clock starts there",
                     ("entry.event",))
        elif have == "none":
            self.add(DEVIATIONS, BLOCK, PRIMARY, "entry event",
                     "THE SOURCE HAS NO ENTRY EVENT AT ALL: no report, sample collection or order time",
                     "naive clock only; no correction and no bias estimate is possible", ("entry.event",))
        elif have in subs:
            self.add(DEVIATIONS, DEVIATION, PRIMARY, "entry event",
                     f"the source gives the {have} time, not the {want} time",
                     f"entry at {have.replace('_', ' ')}, stated in the plan as a deviation",
                     ("entry.event", "entry.acceptable_substitutes"))
        else:
            self.add(DEVIATIONS, BLOCK, PRIMARY, "entry event",
                     f"the source gives the {have} time; the settings want {want} and accept no substitute",
                     "the corrected clock can't start where planned; change the settings (a new plan version) "
                     "or the source", ("entry.event", "entry.acceptable_substitutes"))
        t = self.b.tests
        if "entry_event" in t and "modality" in t:
            entry_tests = t[t["modality"].isin(self.s.read("entry.modalities")).fillna(False)]
            other = sorted(set(entry_tests["entry_event"].dropna()) - {have})
            if other:
                self.add(DEVIATIONS, NOTE, PRIMARY, "entry event",
                         f"some entry tests are timed at another event ({', '.join(other)}) than the profile says",
                         "check the mapping; the times are used as given")

        required = self.s.read("genotype.coverage.required")
        when_genes = self.s.read("genotype.coverage.when_only_genes_listed")
        needs_codons = any(v != "whole_gene" for v in required.values())
        if self.p.panel_coverage == "none":
            self.add(DEVIATIONS, DEVIATION, PRIMARY, "panel coverage",
                     "the source gives no panel coverage", "wild type is unconfirmed and flagged in every result; "
                     "a sensitivity analysis restricted to confirmable patients is added to the plan",
                     ("genotype.coverage.required",))
        elif self.p.panel_coverage == "genes_only" and needs_codons:
            if when_genes == "accept_gene_level":
                self.add(DEVIATIONS, NOTE, PRIMARY, "panel coverage",
                         "the panels list genes, not the codons they read",
                         "coverage is checked at gene level, as the settings say", ("genotype.coverage.when_only_genes_listed",))
            else:
                self.add(DEVIATIONS, DEVIATION, PRIMARY, "panel coverage",
                         "the panels list genes, not the codons they read",
                         "wild type is confirmed at gene level only, stated as a deviation",
                         ("genotype.coverage.when_only_genes_listed",))
        else:
            self.add(READY, "", PRIMARY, "panel coverage", f"coverage is given as {self.p.panel_coverage.replace('_', ' ')}",
                     "wild type can be confirmed", ("genotype.coverage.required",))
        if self.p.panel_coverage != "none":
            self.fusion_detection(required, when_genes)

        if self.p.death_capture in ("partial", "unknown"):
            self.add(DEVIATIONS, DEVIATION, "all", "death capture", f"death capture is {self.p.death_capture}",
                     "a warning goes with every result, and it is stated as a limit")
        else:
            self.add(READY, "", "all", "death capture", "deaths are linked to a death index", "")

        dua = self.p.privacy_min_count
        if dua is not None and self.settings_min_count < dua:
            self.add(DEVIATIONS, BLOCK, "all", "privacy minimum",
                     f"the settings' minimum count ({self.settings_min_count}) is below the data-use agreement's ({dua})",
                     f"not allowed: raise privacy.min_count. This report masks with the agreement's {dua}, and the "
                     "plan can't be locked until the settings are fixed", ("privacy.min_count",))
        for q in self.p.known_quirks:
            self.add(DEVIATIONS, NOTE, "all", "source quirk", q["name"], q["handling"])

    def fusion_detection(self, required: dict, when_genes: str):
        """A panel's files show the genes it reads, never whether it detects their fusions. For each required gene
        whose fusions a group counts, say where its fusion detection is taken from."""
        counted = {gene for g in self.s.read("genotype.groups") for m in g.get("any_of", [])
                   if m["kind"] == "fusion" for gene in m["genes"]}
        genes = [g for g in required if g in counted]
        lists = self.s.read("genotype.coverage.fusion_panels")
        listed = [g for g in genes if g in lists]
        unlisted = [g for g in genes if g not in lists]
        if listed:
            self.add(DEVIATIONS, NOTE, PRIMARY, "fusion detection",
                     f"whether a panel detects {', '.join(listed)} fusions is taken from the settings' list of panels, "
                     "not from the panel files", "a patient whose first panel isn't on the list is not assignable",
                     ("genotype.coverage.fusion_panels",))
            known = set(self.b.panels["panel_id"].dropna()) if "panel_id" in self.b.panels else set()
            absent = sorted({p for g in listed for p in lists[g]} - known)
            if absent:
                self.warn(PRIMARY, "fusion detection",
                          f"panels the settings list as detecting fusions aren't in the panels table: {', '.join(absent)}",
                          "check each name against the panel files: a misspelled name leaves the patients sequenced "
                          "on that panel unassignable", ("genotype.coverage.fusion_panels",))
        if unlisted:
            accept = when_genes == "accept_gene_level"
            self.add(DEVIATIONS, NOTE if accept else DEVIATION, PRIMARY, "fusion detection",
                     f"the groups count {', '.join(unlisted)} fusions, and the settings name no panels that detect them",
                     "a panel that lists the gene is taken to detect its fusions"
                     + (", as the settings say" if accept else ", stated as a deviation"),
                     ("genotype.coverage.fusion_panels", "genotype.coverage.when_only_genes_listed"))

    # ---- the cohort ----
    def build(self) -> Cohort | None:
        if self.keys_missing():
            return None
        try:
            return build_cohort(self.b, self.s, coverage_known=self.p.panel_coverage != "none",
                                entry_times_given=self.entry_times)
        except MissingFieldError as e:
            self.add(MISSING, BLOCK, "all", "cohort rule", f"{e.args[0]} is absent", "the cohort can't be built")
            return None

    def flags(self, c: Cohort):
        e = c.eligible
        labels = {
            "f:entry_before_diagnosis": ("entry before diagnosis", "patients' first test came before diagnosis",
                                         {"enter_at_day_0": "they enter at day 0", "exclude": "they are excluded"}
                                         [self.s.read("entry.before_diagnosis")]),
            "f:unranked_first_test": ("withheld test time", "patients have another test whose time is withheld",
                                      "it was ranked last without knowing whether it came first"),
            "f:first_test_tied": ("first test tied", "patients have tests tied after every tie-break in the settings",
                                  "the lowest test ID was taken; consider another tie-break"),
            "f:after_follow_up_end_unknown": ("time at risk", "patients' first test has no record of whether it came "
                                              "back after follow-up ended", "kept here; the run checks the times"),
            "f:first_test_linked_to_other_cancers": ("test linked to other cancers", "patients' first test is linked to "
                                                     "more than one cancer", "kept and counted"),
            "f:censored_at_age_limit": ("age limit", "patients' survival is withheld by the source's age rule",
                                        "censored at the age limit, as the settings say"),
        }
        for col, (topic, what, effect) in labels.items():
            if col in e:
                k = int(e[col].sum())
                if k:
                    self.add(DEVIATIONS, NOTE, PRIMARY, topic, f"{self.n(k)} eligible {what}", effect)
        groups = [g["name"] for g in self.s.read("genotype.groups") if not g.get("none_of_the_above")]
        if len(groups) >= 2:
            both = e[[f"match:{g}" for g in groups]].astype(bool).sum(axis=1) >= 2
            k = int(both.sum())
            if k:
                self.add(DEVIATIONS, NOTE, PRIMARY, "overlapping groups",
                         f"{self.n(k)} eligible patients match more than one group's rules",
                         "each is counted in the first group by the settings' precedence")

    def missing_data(self, c: Cohort):
        """The complete-case rule, from counts alone (no outcome). When it calls for multiple imputation, the plan says
        so; if the settings don't say how to impute, the adjusted models can't run, and that blocks. The share itself
        isn't shown: with the eligible count, it could give a small count away."""
        p = self.b.patients.drop_duplicates("patient_id").set_index("patient_id").loc[c.eligible.index]
        try:
            coded = code_adjustment(p, self.s)
        except KeyError:
            return                                  # an absent variable is reported missing already
        if missing_data_rule(coded, self.s)[0] == "complete case":
            return
        what = (f"more than {self.s.read('missing_data.complete_case_max_share'):.0%} of eligible patients lack an "
                "adjustment variable")
        if self.s.has("missing_data.imputation_model"):
            self.add(DEVIATIONS, DEVIATION, PRIMARY, "missing data", what,
                     f"the adjusted and weighted models use multiple imputation ({self.s.read('missing_data.imputations')} "
                     "imputations), as the settings' imputation model says",
                     ("missing_data.complete_case_max_share", "missing_data.imputations", "missing_data.imputation_model"))
        else:
            self.add(DEVIATIONS, BLOCK, PRIMARY, "missing data", what,
                     "the plan needs multiple imputation, and the settings don't say how: add "
                     "missing_data.imputation_model before locking", ("missing_data.complete_case_max_share",))

    def share_of(self, k: int, total: int) -> str:
        """A count out of a total, so that neither it nor what is left of the total gives a small count away."""
        if k == 0:
            return "none"
        if is_small(k, self.min_count):
            return f"<{self.min_count}"
        if is_small(total - k, self.min_count):
            return f"all but fewer than {self.min_count}"
        return str(k)

    def reference_drugs(self, c: Cohort):
        """How many of the reference group's eligible patients received a drug the settings list for a driver,
        at any time: a sign of a driver the first test didn't find. Counted only; nobody changes group."""
        drugs = self.s.read("genotype.reference_drugs")
        if not drugs:
            return
        ref = self.s.read("genotype.reference")
        topic = "targeted drugs in the reference"
        tr = self.b.treatments
        if tr is None or not {"patient_id", "regimen"} <= set(tr.columns) or tr["regimen"].isna().all():
            self.add(MISSING, NOTE, PRIMARY, topic, "the treatments table gives no drug names (treatments.regimen)",
                     f"the {ref} group's patients who received a listed drug can't be counted; stated as a limitation",
                     ("genotype.reference_drugs",))
            return
        refs = c.eligible.index[c.eligible["genotype"] == ref]
        tr = tr[tr["patient_id"].isin(refs)]
        counts = {driver: tr.loc[drug_matches(tr["regimen"], names), "patient_id"].nunique()
                  for driver, names in drugs.items()}
        # Drivers can share drugs, so one count can hold another: if any two differ by a small number, every
        # count is shown as a range, so the difference can't be worked out.
        values = list(counts.values())
        close = any(is_small(abs(a - b), self.min_count) for i, a in enumerate(values) for b in values[i + 1:])
        for driver, names in drugs.items():
            k = counts[driver]
            shown = "none" if k == 0 else (self.band(k) if close else self.share_of(k, len(refs)))
            self.add(DEVIATIONS, NOTE, PRIMARY, topic,
                     f"{shown} of the {ref} group's eligible patients received a drug the settings list for {driver} "
                     f"({', '.join(names)}), at any time",
                     "a sign of a driver the first test didn't find, stated as a limitation; they stay in the group",
                     ("genotype.reference_drugs",))

    # ---- the smallest-group rule (D1) ----
    def small_counts(self, analysis: str, counts: dict):
        for g, k in counts.items():
            if is_small(k, self.min_count):
                self.add(DEVIATIONS, NOTE, analysis, "small count",
                         f"the {g} group has fewer than {self.min_count} patients",
                         "flagged now: the safety gate hides it and anything it could be worked out from",
                         ("privacy.min_count",))

    def d1(self, analysis: str, counts: dict, reference: str, if_short: dict | str, settings_read: tuple,
           relative: dict | None = None, total: int | None = None, primary: bool = False,
           extra: dict | None = None, described: list | None = None) -> dict:
        """Decide each comparison with the reference from counts alone. Records sizes and items. `extra`
        holds groups that are shown but not compared (an unknown group); `described`, groups the settings
        say are described only, whatever their size."""
        decisions = {}
        self.small_counts(analysis, {**counts, **(extra or {})})
        self.sizes[analysis] = self.safe_sizes({**counts, **(extra or {})}, total=total, relative=relative,
                                               primary=primary)
        for g in counts:
            if g == reference:
                continue
            comp = f"{g} vs {reference}"
            if g in (described or []):
                decisions[comp] = DESCRIBE
                self.add(DEVIATIONS, NOTE, analysis, "described only", f"the settings describe {g}, never compare it",
                         f"{comp} is described, not compared", ("genotype.described_only",))
                continue
            short = [x for x in (g, reference) if counts[x] < self.min_group]
            if not short:
                decisions[comp] = COMPARE
                continue
            rule = {x: (if_short[x] if isinstance(if_short, dict) else if_short) for x in short}
            if "block" in rule.values():
                decisions[comp] = BLOCKED
                self.add(DEVIATIONS, BLOCK, analysis, "smallest group",
                         f"{', '.join(short)} has fewer than {self.min_group} patients",
                         "this comparison can't run as planned: a person decides what is done instead, and logs it "
                         "as a deviation, before any outcome is seen", settings_read)
            else:
                decisions[comp] = DESCRIBE
                self.add(DEVIATIONS, DEVIATION, analysis, "smallest group",
                         f"{', '.join(short)} has fewer than {self.min_group} patients",
                         f"{comp} is described, not compared", settings_read)
        self.decisions[analysis] = decisions
        return decisions

    def group_counts(self, frame: pd.DataFrame, column: str, names: list[str]) -> dict:
        vc = frame[column].value_counts()
        return {g: int(vc.get(g, 0)) for g in names}

    # ---- analyses ----
    def primary(self, c: Cohort) -> dict:
        names = [g["name"] for g in self.s.read("genotype.groups")]
        ref = self.s.read("genotype.reference")
        counts = self.group_counts(c.eligible, "genotype", names)
        self.d1(PRIMARY, counts, ref, self.s.read("smallest_group.if_short"),
                ("smallest_group.patients", "smallest_group.if_short", "genotype.groups"),
                total=len(c.eligible), primary=True, described=self.s.read("genotype.described_only"))
        return counts

    def keep_mask(self, c: Cohort, keep: dict) -> pd.Series:
        e = c.eligible
        if "within_months" in keep:
            limit = keep["within_months"] * self.s.read("time.days_per_year") / 12
            return (e["entry_day"] <= limit).fillna(False).astype(bool)
        if "field" in keep:
            p = self.b.patients.drop_duplicates("patient_id").set_index("patient_id").loc[e.index]
            return field_condition(p, keep)
        if "exclude_overlap" in keep:
            both = e[[f"match:{g}" for g in keep["exclude_overlap"]]].astype(bool).all(axis=1)
            return ~both
        if "exclude_from_group" in keep:
            how = self.s.read("genotype.variant_text")
            res = self.b.results[self.b.results["test_id"].isin(e["first_test_id"].dropna())]
            hit = res.loc[result_matches(res, keep["if_result"], how), "test_id"]
            has = e["first_test_id"].isin(set(hit))
            return ~((e["genotype"] == keep["exclude_from_group"]) & has)
        raise ValueError(f"unknown subset rule {keep}")

    def subset(self, c: Cohort, a: dict, primary_counts: dict):
        names = [g["name"] for g in self.s.read("genotype.groups")]
        ref = self.s.read("genotype.reference")
        if "within_months" in a["keep"] and not self.entry_times:
            self.add(MISSING, BLOCK, a["id"], "entry event", "the source has no entry times",
                     f"{a['id']} keeps patients by the time of their entry event, so it can't run")
            self.decisions[a["id"]] = {"all comparisons": NOT_RUN}
            return
        try:
            kept = c.eligible[self.keep_mask(c, a["keep"])]
        except MissingFieldError as err:
            self.add(MISSING, BLOCK, a["id"], "subset", f"{err.args[0]} is absent", f"{a['id']} can't run")
            self.decisions[a["id"]] = {"all comparisons": NOT_RUN}
            return
        counts = self.group_counts(kept, "genotype", names)
        self.d1(a["id"], counts, ref, a["if_short"], ("smallest_group.patients", "analyses"), relative=primary_counts,
                described=self.s.read("genotype.described_only"))

    def treatment_clock(self, c: Cohort, a: dict, primary_counts: dict):
        tr = self.b.treatments
        if tr is None or not {"patient_id", "days_to_start", "regimen_order"} <= set(tr.columns):
            self.add(MISSING, BLOCK, a["id"], "treatments", "the treatments table, with patient_id, days_to_start "
                     "and regimen_order, is absent", f"{a['id']} can't run")
            self.decisions[a["id"]] = {"all comparisons": NOT_RUN}
            return
        ok = tr[tr["days_to_start"].notna() & (tr["days_to_start"] >= 0)]
        first = ok.sort_values(["patient_id", "regimen_order"], kind="mergesort").drop_duplicates("patient_id")
        kept = c.eligible[c.eligible.index.isin(first["patient_id"])]
        names = [g["name"] for g in self.s.read("genotype.groups")]
        counts = self.group_counts(kept, "genotype", names)
        self.d1(a["id"], counts, self.s.read("genotype.reference"), a["if_short"],
                ("smallest_group.patients", "analyses"), relative=primary_counts,
                described=self.s.read("genotype.described_only"))

    def biomarker_status(self, c: Cohort, a: dict, timing: str, drop_discordant: bool) -> tuple[pd.Series, pd.Series]:
        status, discordant, absent = biomarker_status(self.b, c, a, timing, drop_discordant)
        for flag in absent:
            self.add(MISSING, DEVIATION, a["id"], "biomarker exclusions", f"patients.{flag} is absent",
                     f"{a['id']} can't leave out the patients this rule names; stated as a deviation")
        return status, discordant

    def biomarker_blocked(self, a: dict) -> str | None:
        """Why a biomarker analysis can't run with these tables, or None."""
        if "value" not in self.b.results.columns:
            return "results.value is absent"
        if a["timing"] == "by_entry_test" and not self.entry_times:
            return "the source has no entry times, so no result can be placed before the entry test"
        return None

    def biomarker(self, c: Cohort, a: dict) -> dict:
        """S1-type analysis. Returns, for its variants, whether it uses its fallback, whether it can run at
        all, and its exact counts (which its variants' shown counts must not give away)."""
        why = self.biomarker_blocked(a)
        if why:
            self.add(MISSING, BLOCK, a["id"], "biomarker", why, f"{a['id']} can't run")
            self.decisions[a["id"]] = {"all comparisons": NOT_RUN}
            return {"blocked": True}
        status, _ = self.biomarker_status(c, a, a["timing"], False)
        unknown = int((status == a["unknown_group"]).sum())
        total = len(status)
        cutoff = self.s.read("biomarker.unknown_max_share")
        fallback = total > 0 and unknown / total > cutoff
        names = [g["name"] for g in a["groups"]]
        counts = {g: int((status == g).sum()) for g in names}
        sizes_given = panel_sizes_given(self.b, self.p)
        share_text = f"more than {cutoff:.0%}" if fallback else f"{cutoff:.0%} or fewer"
        if fallback:
            form = a["fallback"]["when_panel_sizes_given"] if sizes_given else a["fallback"]["otherwise"]
            if form == "mutation_burden":
                self.burden(c, a, share_text)
            else:
                self.add(DEVIATIONS, DEVIATION, a["id"], "biomarker unknown",
                         f"{share_text} of {a['id']}'s patients have {a['unknown_group']} status, and the source "
                         "gives no panel sizes", f"{a['id']} is descriptive only (the settings' fallback)",
                         ("biomarker.unknown_max_share", "analyses"))
                self.decisions[a["id"]] = {f"{g} vs {a['reference']}": DESCRIBE for g in names if g != a["reference"]}
                # Its groups and Unknown add up to the eligible patients less those it leaves out, so the sizes are
                # ranges if that difference is small.
                all_groups = {**counts, a["unknown_group"]: unknown}
                self.small_counts(a["id"], all_groups)
                self.sizes[a["id"]] = self.safe_sizes(all_groups, total=len(c.eligible))
                # the smallest-group rule is still reported, as the registered design check did
                short = [g for g in names if counts[g] < self.min_group]
                if short:
                    self.add(DEVIATIONS, DEVIATION, a["id"], "smallest group",
                             f"{', '.join(short)} has fewer than {self.min_group} patients",
                             f"{a['id']}'s groups are described, not compared", ("smallest_group.patients",))
            for other in a["not_run_under_fallback"]:
                self.add(DEVIATIONS, BLOCK, other, "biomarker unknown",
                         f"{a['id']} uses its fallback", f"{other} is not run (the settings name it as not run then)",
                         ("analyses",))
                self.decisions[other] = {"all comparisons": NOT_RUN}
        else:
            self.d1(a["id"], counts, a["reference"], a["if_short"], ("smallest_group.patients", "analyses"),
                    total=len(c.eligible), extra={a["unknown_group"]: unknown})
            self.add(DEVIATIONS, NOTE, a["id"], "biomarker unknown",
                     f"{share_text} of {a['id']}'s patients have {a['unknown_group']} status",
                     f"the {a['unknown_group']} group is described, not compared", ("biomarker.unknown_max_share",))
        return {"fallback": fallback, "counts": {**counts, a["unknown_group"]: unknown}}

    def burden(self, c: Cohort, a: dict, share_text: str):
        """A biomarker's fallback by mutation burden: high against low (the reference), counted on the entry test
        (`mutation_burden`), with the smallest-group rule as for any comparison. Blocked if it can't be counted."""
        mb = a["mutation_burden"]
        problems = burden_problems(self.b)
        if problems:
            self.add(MISSING, BLOCK, a["id"], "mutation burden", "; ".join(problems),
                     f"{a['id']} falls back to mutation burden and can't count it, so it can't run: a person decides "
                     "what is done instead", ("analyses",))
            self.decisions[a["id"]] = {"all comparisons": NOT_RUN}
            return
        self.add(DEVIATIONS, DEVIATION, a["id"], "biomarker unknown",
                 f"{share_text} of {a['id']}'s patients have {a['unknown_group']} status",
                 f"{a['id']} compares {BURDEN_HIGH} ({mb['threshold_per_megabase']:g} or more mutations per megabase "
                 f"on the entry test) with {BURDEN_LOW}, labeled a hypermutation proxy (the settings' fallback)",
                 ("biomarker.unknown_max_share", "analyses"))
        group, _ = mutation_burden(self.b, c, a)
        counts = {g: int((group == g).sum()) for g in (BURDEN_LOW, BURDEN_HIGH)}
        unknown = int((group == a["unknown_group"]).sum())
        self.d1(a["id"], counts, BURDEN_LOW, a["if_short"], ("smallest_group.patients", "analyses"),
                total=len(c.eligible), extra={a["unknown_group"]: unknown})
        self.add(DEVIATIONS, NOTE, a["id"], "mutation burden",
                 f"a patient whose entry test's panel has no coding size has no mutation burden: {a['unknown_group']}",
                 f"the {a['unknown_group']} group is described, not compared", ("analyses",))

    def biomarker_variant(self, c: Cohort, a: dict, parent: dict, parent_state: dict):
        if a["id"] in self.decisions:
            return
        if parent_state.get("blocked"):
            self.add(MISSING, BLOCK, a["id"], "biomarker", f"{parent['id']}, which it varies, can't run",
                     f"{a['id']} can't run either")
            self.decisions[a["id"]] = {"all comparisons": NOT_RUN}
            return
        if parent_state.get("fallback"):
            return
        timing = a.get("timing", parent["timing"])
        status, _ = self.biomarker_status(c, parent, timing, a.get("drop_discordant", False))
        names = [g["name"] for g in parent["groups"]]
        counts = {g: int((status == g).sum()) for g in names}
        unknown = int((status == parent["unknown_group"]).sum())
        self.d1(a["id"], counts, parent["reference"], parent["if_short"], ("smallest_group.patients", "analyses"),
                relative=parent_state.get("counts"), extra={parent["unknown_group"]: unknown})

    def analyses(self, c: Cohort, primary_counts: dict):
        all_a = self.s.read("analyses")
        by_id = {a["id"]: a for a in all_a}
        state = {}
        for a in all_a:
            if not a["run"]:
                self.add(READY, "", a["id"], "switched off", f"{a['label']} is switched off in the settings", "not run")
                self.decisions[a["id"]] = {"all comparisons": "switched off"}
                continue
            k = a["kind"]
            if k == "subset":
                self.subset(c, a, primary_counts)
            elif k == "treatment_clock":
                self.treatment_clock(c, a, primary_counts)
            elif k == "biomarker":
                state[a["id"]] = self.biomarker(c, a)
        for a in all_a:
            if a["run"] and a["kind"] == "biomarker_variant":
                self.biomarker_variant(c, a, by_id[a["of"]], state.get(a["of"], {}))

    # ---- guardrails: warnings, never overrides ----
    def guardrails(self, c: Cohort | None):
        dpy = self.s.read("time.days_per_year")
        horizon_days = self.s.read("time.horizon_years") * dpy
        if self.min_group < self.min_count:
            self.warn("all", "smallest group", f"the smallest group compared ({self.min_group}) is below the privacy "
                      f"minimum count ({self.min_count})", "a group's own counts could be compared but not shown",
                      ("smallest_group.patients", "privacy.min_count"))
        window = self.s.read("survival_after_entry.window_days")
        if window >= horizon_days:
            self.warn(PRIMARY, "survival-after-entry check", f"its window ({window:g} days) reaches the horizon "
                      f"({horizon_days:g} days)", "if entry moves by the window, no time at risk is left inside the horizon",
                      ("survival_after_entry.window_days", "time.horizon_years"))
        for a in self.s.read("analyses"):
            if a["run"] and a["kind"] == "subset" and "within_months" in a["keep"]:
                if a["keep"]["within_months"] * dpy / 12 >= horizon_days:
                    self.warn(a["id"], "short-wait cutoff", "its cutoff reaches the horizon", "it keeps almost everyone",
                              ("analyses", "time.horizon_years"))
        if self.p.age_top_code is not None:
            for r in self.s.read("cohort.rules"):
                for cond in r["all_of"]:
                    if cond.get("field") == "age_at_diagnosis" and cond.get("op") in ("ge", "gt") \
                            and cond["value"] > self.p.age_top_code:
                        self.warn(PRIMARY, "cohort rule", f"rule {r['id']} needs an age above the source's top code "
                                  f"({self.p.age_top_code:g})", "no patient can pass it", ("cohort.rules",))
        limit = self.s.read("summary.median_limit_years")
        if limit < self.s.read("time.horizon_years"):
            self.warn(PRIMARY, "median limit", "medians are cut off before the RMST horizon",
                      "medians may be 'not reached' where RMST is still estimated", ("summary.median_limit_years",))
        if c is not None and len(c.eligible):
            share = c.eligible["f:entry_before_diagnosis"].mean()
            if share > 0.5:
                self.warn(PRIMARY, "entry before diagnosis", "most eligible patients' first test came before diagnosis",
                          "check that the time zero and the entry event are mapped the right way round", ("entry.event",))

    def deferred(self) -> list[dict]:
        window = self.s.read("survival_after_entry.window_days")
        deaths = self.s.read("survival_after_entry.min_deaths")
        if_fewer = self.s.read("survival_after_entry.if_fewer")
        return [
            {"check": "survival after the entry event",
             "rule": f"if fewer than {deaths} eligible patients died within {window:g} days of their entry event, the "
                     "source in effect required that survival; then "
                     + ("entry moves later by the window, a pre-planned change that is logged, and the smallest-group "
                        "rule is applied again to the patients left" if if_fewer == "move_entry_by_window"
                        else "a warning is given")},
            {"check": "time at risk from the times themselves",
             "rule": "a patient whose follow-up ends on or before entry has no time at risk and is excluded and counted; "
                     "the screen used only the source's own record of a test that came back after follow-up ended"},
            {"check": "the horizon against follow-up",
             "rule": f"a warning if follow-up rarely reaches the horizon ({self.s.read('time.horizon_years'):g} years)"},
            {"check": "the smallest-group rule on the analysis cohort",
             "rule": f"applied again after any change above, with the same minimum ({self.min_group})"},
        ]

    def ready(self):
        parents = {x["id"]: x["of"] for x in self.s.read("analyses") if x["kind"] == "biomarker_variant"}
        for a in [PRIMARY] + [x["id"] for x in self.s.read("analyses") if x["run"]]:
            mine = [i for i in self.items if i.analysis in (a, "all") and i.list in (MISSING, DEVIATIONS)
                    and i.severity in (BLOCK, DEVIATION)]
            # a variant can't be readier than the analysis it varies
            mine += [i for i in self.items if a in parents and i.analysis == parents[a] and i.severity == BLOCK]
            if not mine:
                label = "the primary analysis" if a == PRIMARY else a
                self.add(READY, "", a, "runs as standard", f"{label} can run as planned", "")


def screen(blind: OutcomeBlindTables, settings: Settings, profile: CapabilityProfile,
           source_steps: list[tuple[str, int]] | None = None,
           inputs: list[str | Path] | None = None) -> ReadinessReport:
    """The readiness report. `source_steps` are the adapter's own exclusions, as (label, count), shown at
    the top of the flow (for example patients the source holds without a cancer of the study's kind).
    `inputs`, the files the tables were read from, are recorded by name and SHA-256, so the plan can be
    locked only with the inputs it was screened from."""
    if isinstance(blind, StandardTables) or not isinstance(blind, OutcomeBlindTables):
        raise TypeError("the screen takes outcome_blind(tables), never the tables with their outcomes")
    if any(c in blind.patients.columns for c in OUTCOME_FIELDS):
        raise TypeError("outcome columns reached the screen")
    settings.forget_reads()
    sc = _Screen(blind, settings, profile, list(source_steps or []))
    sc.structure()
    sc.needs()
    sc.profile_rules()
    c = sc.build()
    flow_rows = None
    if c is not None:
        sc.flags(c)
        sc.missing_data(c)
        sc.reference_drugs(c)
        counts = sc.primary(c)
        sc.analyses(c, counts)
        flow = c.flow(before=sc.source_steps)
        try:
            flow_rows = public_flow(flow[["label", "excluded", "remaining"]], min_count=sc.min_count).to_dict("records")
        except ValueError:
            sc.add(DEVIATIONS, NOTE, PRIMARY, "flow diagram", "no merging of steps makes the flow diagram safe to show",
                   "it is not shown")
    sc.guardrails(c)
    deferred = sc.deferred()
    sc.ready()
    study = settings.read("study.name")
    used = settings.reads()
    hashed = [{"file": Path(p).name, "sha256": sha256_file(Path(p))} for p in sorted({Path(x).resolve() for x in inputs})] \
        if inputs else []
    return ReadinessReport(
        study=study, settings_sha256=settings.sha256, profile=profile.to_dict(), items=sc.items,
        warnings=sc.warnings, decisions=sc.decisions, group_sizes=sc.sizes, data_checks=sc.checks,
        missingness=sc.missingness, flow=flow_rows, deferred=deferred,
        settings_used=settings.provenance(used), defaults_used=settings.defaults_used(used), min_count=sc.min_count,
        settings_resolved_sha256=settings.resolved_sha256, inputs=hashed)


__all__ = ["screen", "ReadinessReport", "Item", "READY", "MISSING", "DEVIATIONS", "BLOCK", "DEVIATION", "NOTE",
           "COMPARE", "DESCRIBE", "NOT_RUN", "BLOCKED", "PRIMARY"]
