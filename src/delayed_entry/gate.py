"""The safety gate (design plan section 9): nothing reaches the output folder unless it passes.

Generalized from the CRC BPC project's crcbpc.outputs, with the minimum count from the settings and the data-use
agreement (the larger of the two). The privacy protections:
1. a count from 1 to the minimum less one is hidden ("<minimum"), and so is any statistic from fewer patients;
2. no hidden count can be worked out by subtracting published counts, within one table or across tables
   (disclosed() checks every relation together, and the gate stops if hiding more can't make it safe);
3. a curve's steps give away its numbers at risk, so a group's curves are withheld if they would give away a count
   the risk table hides (withheld_curves);
4. a group whose size is hidden shows no survival estimate at a fixed time, which could give the size away
   (hide_survival_of_hidden_groups);
5. balance statistics are published at a precision that can't, with published totals, rebuild a table holding a
   small count (balance_precision).
The security checks (check_outputs): only tables, figures and reports; nothing row-level (an ID from the data, a
row-level column, a date that isn't the run's own); no count column holding a small count; no curve drawn past a
hidden number at risk; no survival estimate for a hidden group size.
"""
from __future__ import annotations

import itertools
import math
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from delayed_entry.privacy import is_small as _is_small

NOT_SHOWN = "not shown"                # a count at or above the minimum, hidden so a smaller one can't be worked out
ROW_LEVEL_COLUMNS = {"patient_id", "test_id", "record_id", "sample_id", "tumor_sample_barcode", "cpt_genie_sample_id"}
COUNT_COLUMNS = {"patients", "excluded", "remaining", "at_risk", "n"}
PUBLISHABLE_SUFFIXES = (".csv", ".svg", ".md", ".html", ".json")
DATE = re.compile(r"\b(?:19|20)\d\d-[01]\d-[0-3]\d\b|\b[01]?\d/[0-3]?\d/(?:19|20)\d\d\b")
JOINT_TABLES_LIMIT = 10 ** 6          # beyond this many candidate tables, a trait is checked term by term only


class GateError(ValueError):
    """An output can't be published safely as it stands. A person merges or hides more, and logs it."""


@dataclass
class Gate:
    """The minimum count and a log of every check and what it hid."""
    min_count: int

    def __post_init__(self):
        self.log: list[dict] = []

    def small(self, n) -> bool:
        return _is_small(int(n), self.min_count)

    @property
    def small_text(self) -> str:
        return f"<{self.min_count}"

    def record(self, check: str, verdict: str, detail: str = "") -> None:
        self.log.append({"check": check, "verdict": verdict, "detail": detail})


# ---------------------------------------------------------------------------
# 2. Counts found by subtraction
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Relation:
    """total = the sum of parts. Closed: exactly. Open: plus a remainder that is never published (for example the
    patients not yet entered: naive at risk = corrected at risk + not yet entered)."""
    total: str
    parts: tuple[str, ...]
    closed: bool = True


def disclosed(values: dict[str, int], shown: set[str], relations: list[Relation], gate: Gate) -> list[str]:
    """Every small count the shown counts give away: each hidden count found by subtraction, in as many steps as
    it takes, then what each known total leaves over its known parts. Empty means the counts are safe."""
    known = {k: values[k] for k in shown}
    changed = True
    while changed:
        changed = False
        for r in relations:
            unknown = [k for k in (r.total, *r.parts) if k not in known]
            if r.closed and len(unknown) == 1:
                u = unknown[0]
                known[u] = (sum(known[p] for p in r.parts) if u == r.total
                            else known[r.total] - sum(known[p] for p in r.parts if p != u))
                changed = True
    found = [f"{k}" for k, v in known.items() if k not in shown and gate.small(v)]
    for r in relations:
        if r.total not in known:
            continue
        hidden = [p for p in r.parts if p not in known]
        rest = known[r.total] - sum(known[p] for p in r.parts if p in known)
        if (hidden or not r.closed) and gate.small(rest):
            found.append(f"{r.total} less its parts")
    return found


def require_safe(values, shown, relations, what: str, gate: Gate) -> None:
    left = disclosed(values, shown, relations, gate)
    if left:
        raise GateError(f"{what}: counts below {gate.min_count} could be worked out ({'; '.join(left)}); "
                        "a person must merge or hide more, and log it")


def show(values: dict[str, int], key: str, shown: set[str], gate: Gate):
    """A count as published: the number, "<minimum", or "not shown"."""
    if key in shown:
        return values[key]
    return gate.small_text if gate.small(values[key]) else NOT_SHOWN


# ---------------------------------------------------------------------------
# counts.csv: each analysis's group sizes
# ---------------------------------------------------------------------------
@dataclass
class CountSpec:
    """What count_values needs from the run: exact counts, and which analyses rest on which."""
    eligible: int
    primary: dict[str, int]                       # group: patients in the analysis
    biomarkers: dict[str, dict]                   # id: {"groups": {g: n}, "unknown": (name, n), "not_in": (name, n)}
    others: list[tuple[str, str, dict[str, int]]]  # (analysis, base analysis, {group: n})
    first_groups: dict[str, str]                  # variant id with any-time results: its parent's first group


def count_values(spec: CountSpec) -> tuple[dict[str, int], list[Relation], list[str]]:
    v = {"eligible": spec.eligible, "primary | in the analysis": sum(spec.primary.values())}
    order = ["eligible", "primary | in the analysis"]
    for g, n in spec.primary.items():
        v[f"primary | {g}"] = n
        order.append(f"primary | {g}")
    rel = [Relation("primary | in the analysis", tuple(f"primary | {g}" for g in spec.primary))]
    for b, d in spec.biomarkers.items():
        keys = []
        for g, n in d["groups"].items():
            v[f"{b} | {g}"] = n
            keys.append(f"{b} | {g}")
        (un, n_un), (ni, n_ni) = d["unknown"], d["not_in"]
        v[f"{b} | {un}"], v[f"{b} | {ni}"] = n_un, n_ni
        merged = f"{b} | {un}, or {ni}"
        v[merged] = n_un + n_ni
        order += keys + [f"{b} | {un}", f"{b} | {ni}", merged]
        rel += [Relation("primary | in the analysis", tuple(keys + [f"{b} | {un}", f"{b} | {ni}"])),
                Relation(merged, (f"{b} | {un}", f"{b} | {ni}"))]
    rel.append(Relation("eligible", ("primary | in the analysis",), closed=False))
    for a, base, groups in spec.others:
        for g, n in groups.items():
            v[f"{a} | {g}"] = n
            order.append(f"{a} | {g}")
        if a in spec.first_groups:                  # every result: the parent's first group can only grow
            g = spec.first_groups[a]
            rel.append(Relation(f"{a} | {g}", (f"{base} | {g}",), closed=False))
        else:                                       # the analysis keeps some of its base's patients
            rel += [Relation(f"{base} | {g}", (f"{a} | {g}",), closed=False) for g in groups]
    return v, rel, order


def publish_counts(values, relations, spec: CountSpec, gate: Gate) -> set[str]:
    """Which counts are shown: none from 1 to the minimum less one; a biomarker's unknown and set-aside patients
    merged if either is small, and if the two together are small too, its largest group hidden as well (else the
    groups would leave them over from the analysis's total); an analysis's group hidden if it differs from its
    base's by a small number; then every relation checked together."""
    shown = {k for k, n in values.items() if not gate.small(n)}
    bases = {"primary", *spec.biomarkers}
    for b, d in spec.biomarkers.items():
        un, ni = f"{b} | {d['unknown'][0]}", f"{b} | {d['not_in'][0]}"
        if gate.small(values[un]) or gate.small(values[ni]):
            shown -= {un, ni}
            if gate.small(values[f"{b} | {d['unknown'][0]}, or {d['not_in'][0]}"]):
                largest = max(d["groups"], key=lambda g: d["groups"][g])
                shown.discard(f"{b} | {largest}")
        else:
            shown.discard(f"{b} | {d['unknown'][0]}, or {d['not_in'][0]}")
    for r in relations:
        if not r.closed and r.total != "eligible" and gate.small(values[r.total] - values[r.parts[0]]):
            shown.discard(r.parts[0] if r.total.split(" | ")[0] in bases else r.total)
    base_of = {a: base for a, base, _ in spec.others}
    for k, x in values.items():
        a, sep, g = k.partition(" | ")
        if sep and a in base_of:
            base = f"{base_of[a]} | {g}"
            if base in values and gate.small(abs(x - values[base])):
                shown.discard(k)
    require_safe(values, shown, relations, "counts", gate)
    return shown


# ---------------------------------------------------------------------------
# descriptives.csv: the flags and the distributions
# ---------------------------------------------------------------------------
def descriptive_values(eligible: int, flags: list[tuple[str, int]], tables: list[tuple[str, dict[str, int]]]
                       ) -> tuple[dict[str, int], list[Relation], list[str]]:
    v = {"eligible": eligible}
    rel, order = [], []
    for label, n in flags:
        k = f"Kept but counted | {label}"
        v[k] = n
        order.append(k)
        rel.append(Relation("eligible", (k,), closed=False))
    for table, levels in tables:
        keys = []
        for level, n in sorted(levels.items()):
            k = f"{table} | {level}"
            v[k] = n
            keys.append(k)
        order += keys
        rel.append(Relation("eligible", tuple(keys)))
    return v, rel, order


def publish_descriptives(values, relations, gate: Gate) -> set[str]:
    """Which descriptive counts are shown. A small count is hidden, and in a distribution the next smallest is
    hidden too, as often as it takes, so the hidden ones can't be worked out from the total. A flag is hidden if
    the patients without it are few. Then every relation is checked together."""
    shown = {k for k, n in values.items() if not gate.small(n)} | {"eligible"}
    for r in relations:
        if not r.closed:
            if gate.small(values[r.total] - values[r.parts[0]]):
                shown.discard(r.parts[0])
            continue
        while True:
            hidden = [p for p in r.parts if p not in shown]
            if not hidden or (len(hidden) > 1 and not gate.small(sum(values[p] for p in hidden))):
                break
            more = sorted((p for p in r.parts if p in shown and values[p] > 0), key=lambda p: values[p])
            if not more:
                break
            shown.discard(more[0])
    require_safe(values, shown, relations, "descriptives", gate)
    return shown


def hide_across(values, shown, relations, cross: list[Relation], gate: Gate) -> set[str]:
    """Hide the parts of the tables' cross relations, the smallest first, until together the shown counts give
    nothing away. Raises if hiding every part isn't enough."""
    shown = set(shown)
    for r in cross:
        while disclosed(values, shown, relations + cross, gate):
            more = sorted((p for p in r.parts if p in shown and values[p] > 0), key=lambda p: values[p])
            if not more:
                break
            shown.discard(more[0])
    require_safe(values, shown, relations + cross, "counts and descriptives together", gate)
    return shown


# ---------------------------------------------------------------------------
# The risk table
# ---------------------------------------------------------------------------
def risk_key(p: str, g: str, y) -> str:
    return f"{p} | {g} | {y:g} y"


def group_key(g: str) -> str:
    """A group's size in the risk table's arithmetic: published in counts.csv, not in the risk table."""
    return f"patients | {g}"


def risk_relations(all_label: str, groups: list[str], years: list) -> list[Relation]:
    rel = []
    for y in years:
        rel += [Relation(risk_key(p, all_label, y), tuple(risk_key(p, g, y) for g in groups)) for p in ("naive", "corrected")]
        rel += [Relation(risk_key("naive", g, y), (risk_key("corrected", g, y),), closed=False) for g in [all_label] + groups]
    for g in [all_label] + groups:
        marks = [group_key(g)] + [risk_key("naive", g, y) for y in years]
        rel += [Relation(a, (b,), closed=False) for i, a in enumerate(marks) for b in marks[i + 1:]]
        rel += [Relation(group_key(g), (risk_key("corrected", g, y),), closed=False) for y in years]
    return rel


def publish_risk(values, all_label: str, groups: list[str], years: list, gate: Gate) -> set[str]:
    """Which numbers at risk are shown: not small; not a naive count that is a small number more than the
    corrected one; not one a small number fewer than the group's size or than a shown earlier naive count; not a
    total if any of its groups is hidden. Then every relation is checked together."""
    shown = {k for k, n in values.items() if not gate.small(n)}
    for y in years:
        for g in [all_label] + groups:
            if gate.small(values[risk_key("naive", g, y)] - values[risk_key("corrected", g, y)]):
                shown.discard(risk_key("naive", g, y))
    for g in [all_label] + groups:
        marks = [k for k in [group_key(g)] + [risk_key("naive", g, y) for y in years] if k in values]
        for i, b in enumerate(marks):
            if b in shown and any(a in shown and gate.small(values[a] - values[b]) for a in marks[:i]):
                shown.discard(b)
        if group_key(g) in shown:
            for y in years:
                if gate.small(values[group_key(g)] - values[risk_key("corrected", g, y)]):
                    shown.discard(risk_key("corrected", g, y))
    for y in years:
        for p in ("naive", "corrected"):
            if any(risk_key(p, g, y) not in shown for g in groups):
                shown.discard(risk_key(p, all_label, y))
    require_safe(values, shown, risk_relations(all_label, groups, years), "risk table", gate)
    return shown


def risk_table(passes: pd.DataFrame, all_label: str, groups: list[str], group_n: dict[str, int] | None,
               gate: Gate) -> pd.DataFrame:
    """Patients at risk at each fixed time, by pass and group. `group_n`, each group's size as counts.csv
    publishes it, lets the rules check the two tables together."""
    r = passes[passes["quantity"] == "survival"]
    values = {risk_key(x.pass_, x.group, float(x.years)): int(x.n_risk)
              for x in r.rename(columns={"pass": "pass_"}).itertuples()}
    years = sorted({float(y) for y in r["years"]})
    if group_n:
        values |= {group_key(g): int(group_n[g]) for g in [all_label] + groups}
    shown = publish_risk(values, all_label, groups, years, gate)
    return pd.DataFrame([{"pass": p, "group": g, "years": _num(y), "at_risk": show(values, risk_key(p, g, y), shown, gate)}
                         for p in ("naive", "corrected") for g in [all_label] + groups for y in years])


def _num(y: float):
    return int(y) if float(y).is_integer() else y


# ---------------------------------------------------------------------------
# 3. The curves
# ---------------------------------------------------------------------------
STEP_COLUMNS = ["pass", "group", "day", "surv", "n_risk", "n_event"]


def drawn_steps(curves: pd.DataFrame, group_n: dict[str, int], last_day: float, gate: Gate) -> pd.DataFrame:
    """The KM step functions as they may be drawn: only where at least the minimum is at risk, to the last summary
    time, and only the steps (the times survival drops), plus where each drawn stretch ends. The naive curve starts
    at 1 on day 0, with the whole group at risk. The numbers at risk stay only to judge the steps."""
    out = []
    for (p, g), d in curves.groupby(["pass", "group"], sort=False):
        d = d[["day", "surv", "n_risk", "n_event"]].sort_values("day")
        if p == "naive":
            d = pd.concat([pd.DataFrame([{"day": 0.0, "surv": 1.0, "n_risk": group_n[g], "n_event": 0}]), d])
        d = d[(d["n_risk"] >= gate.min_count) & (d["day"] <= last_day)]
        if d.empty:
            continue
        keep = (d["n_event"] > 0) | (d["day"] == d["day"].max()) | ((p == "naive") & (d["day"] == 0))
        out.append(d[keep].assign(**{"pass": p, "group": g}))
    return pd.concat(out, ignore_index=True)[STEP_COLUMNS] if out else pd.DataFrame(columns=STEP_COLUMNS)


def years_text(y) -> str:
    return f"{y:g} year" if y == 1 else f"{y:g} years"


def withheld_curves(steps: pd.DataFrame, risk: pd.DataFrame, years: list, days_per_year: float, gate: Gate) -> dict[str, str]:
    """The groups whose curves are not published, and why. Each step gives away its number at risk (survival
    falls by the share of those at risk who died). So, taking every drawn step's number at risk as known, a
    group's curves (both passes) are withheld if they would give away a number at risk the risk table hides (the
    curve is drawn past that time), or a small count of patients alive whose entry event hasn't come yet (naive
    minus corrected at risk at a step drawn in both passes)."""
    why = {}
    for g, d in steps.groupby("group", sort=False):
        reasons = []
        for p in ("naive", "corrected"):
            drawn_to = d.loc[d["pass"] == p, "day"].max()
            for y in years:
                cell = risk.loc[(risk["pass"] == p) & (risk["group"] == g) & (risk["years"].astype(float) == float(y)), "at_risk"]
                if len(cell) and str(cell.iloc[0]) in (NOT_SHOWN, gate.small_text) and drawn_to > y * days_per_year:
                    reasons.append(f"the {p} number at risk at {years_text(y)}, which the risk table hides")
        ev = d[d["n_event"] > 0]
        naive = ev[ev["pass"] == "naive"].groupby("day")["n_risk"].first()
        corrected = ev[ev["pass"] == "corrected"].groupby("day")["n_risk"].first()
        if any(gate.small(naive[t] - corrected[t]) for t in naive.index.intersection(corrected.index)):
            reasons.append(f"a count below {gate.min_count} of patients alive whose entry event hadn't come")
        if reasons:
            why[g] = "; ".join(reasons)
    return why


def public_curves(steps: pd.DataFrame, withheld, days_per_month: float) -> pd.DataFrame:
    """The drawn steps as published: months and survival only, without the groups in `withheld`."""
    s = steps[~steps["group"].isin(list(withheld))]
    return pd.DataFrame({"pass": s["pass"], "group": s["group"], "months": (s["day"] / days_per_month).round(2),
                         "survival": s["surv"].astype(float).round(4)}).reset_index(drop=True)


# ---------------------------------------------------------------------------
# 1 and 4. Estimates from few patients, and survival estimates that pin a hidden size
# ---------------------------------------------------------------------------
def small_group_rule(t: pd.DataFrame, behind: list[int], gate: Gate) -> pd.DataFrame:
    """Blank any estimate from fewer than the minimum."""
    t = t.copy()
    few = [n < gate.min_count for n in behind]
    t.loc[few, ["estimate", "lower", "upper"]] = None
    t.loc[few, "note"] = f"fewer than {gate.min_count} patients"
    return t


SIZE_NOTE = "not shown: counts.csv hides this group's size, and this estimate could give it away"


def hide_survival_of_hidden_groups(t: pd.DataFrame, hidden) -> pd.DataFrame:
    """Blank the survival estimates at fixed times of each (analysis, group) in `hidden`, the groups whose size
    counts.csv doesn't show: when nobody was censored before it, survival is survivors over the group's size, so
    at six significant digits it can pin the size. Medians, RMST and contrasts stay: none is a fraction of the group."""
    t = t.copy()
    in_hidden = pd.Series([(a, g) in hidden for a, g in zip(t["analysis"], t["group"])], index=t.index)
    rows = t["quantity"].str.startswith("survival at") & in_hidden & t["estimate"].notna()
    t.loc[rows, ["estimate", "lower", "upper"]] = None
    t.loc[rows, "note"] = SIZE_NOTE
    return t


def add_note(t: pd.DataFrame, rows: pd.Series, note: str) -> pd.DataFrame:
    t = t.copy()
    old = t.loc[rows, "note"].fillna("")
    t.loc[rows, "note"] = [f"{o}; {note}" if o else note for o in old]
    return t


# ---------------------------------------------------------------------------
# 5. Balance statistics
# ---------------------------------------------------------------------------
def binary_smd(counts: list[int], sizes: list[int]) -> float:
    """The largest absolute standardized mean difference across pairs of groups for a yes/no term, as cobalt
    computes it for several groups: the spread of the groups' shares over the square root of the mean of their
    p(1 - p)."""
    p = [a / n for a, n in zip(counts, sizes)]
    sd = math.sqrt(sum(x * (1 - x) for x in p) / len(p))
    return 0.0 if sd == 0 else (max(p) - min(p)) / sd


def _columns(total: int, sizes: list[int], target: float, decimals: int) -> list[tuple[int, ...]]:
    """Every column (a count per group, adding up to `total`) whose SMD rounds to `target` at `decimals`."""
    out = []

    def walk(prefix, left, i):
        if i == len(sizes) - 1:
            if 0 <= left <= sizes[i]:
                col = prefix + (left,)
                if round(binary_smd(list(col), sizes), decimals) == target:
                    out.append(col)
            return
        for a in range(max(0, left - sum(sizes[i + 1:])), min(sizes[i], left) + 1):
            walk(prefix + (a,), left - a, i + 1)

    walk((), total, 0)
    return out


def pinned_small_cells(table: pd.DataFrame, smds: dict[str, float], decimals: int, gate: Gate) -> list[str]:
    """For one trait whose per-level totals are published (rows: groups, published sizes; columns: levels), the
    small cells that the SMDs at `decimals`, with the totals, pin: every table that fits gives them the same value.
    Whole tables are enumerated when few enough fit; otherwise each level is judged alone."""
    sizes = [int(x) for x in table.sum(axis=1)]
    cand = {}
    for lvl in table.columns:
        col = tuple(int(x) for x in table[lvl])
        target = round(smds[lvl], decimals)
        cand[lvl] = _columns(sum(col), sizes, target, decimals) or [col]
    levels = list(table.columns)
    fits: list[list[tuple[int, ...]]]
    product = math.prod(len(cand[lv]) for lv in levels[:-1])
    if product <= JOINT_TABLES_LIMIT:
        fits = []
        last = set(cand[levels[-1]])
        for combo in itertools.product(*(cand[lv] for lv in levels[:-1])):
            rest = tuple(n - sum(c[i] for c in combo) for i, n in enumerate(sizes))
            if rest in last:
                fits.append(list(combo) + [rest])
        per_level = {lv: {f[j] for f in fits} for j, lv in enumerate(levels)}
    else:
        per_level = {lv: set(cand[lv]) for lv in levels}
    pinned = []
    for lv in levels:
        cols = per_level[lv]
        for i, g in enumerate(table.index):
            values = {c[i] for c in cols}
            if len(values) == 1:
                v = next(iter(values))
                if gate.small(v) or gate.small(sizes[i] - v):
                    pinned.append(f"{g} x {lv}")
    return pinned


def balance_precision(traits: dict[str, tuple[pd.DataFrame, dict[str, float]]], most: int, gate: Gate) -> int:
    """The most decimals, up to `most`, at which no trait's before-weighting SMDs pin a small cell. `traits`: for
    each trait whose per-level totals another table publishes, its group-by-level table of patients in the model
    and each level's SMD before weighting. Raises if even whole numbers would."""
    for d in range(most, -1, -1):
        problems = [f"{t}: {c}" for t, (table, smds) in traits.items() for c in pinned_small_cells(table, smds, d, gate)]
        if not problems:
            return d
    raise GateError("balance before weighting: even whole numbers would let a small cell be worked out")


# ---------------------------------------------------------------------------
# The security checks: read back everything the gate would publish
# ---------------------------------------------------------------------------
def curves_past_hidden_counts(curves: pd.DataFrame, risk: pd.DataFrame, gate: Gate) -> list[str]:
    problems = []
    for (p, g), d in curves.groupby(["pass", "group"], sort=False):
        for y in sorted({float(x) for x in risk["years"]}):
            cell = risk.loc[(risk["pass"] == p) & (risk["group"] == g) & (risk["years"].astype(float) == y), "at_risk"]
            if len(cell) and str(cell.iloc[0]) in (NOT_SHOWN, gate.small_text) and d["months"].astype(float).max() > 12 * y:
                problems.append(f"curves.csv: the {p} curve for {g} is drawn past {years_text(y)}, where the risk "
                                "table hides the number at risk")
    return problems


def survival_of_hidden_sizes(counts: pd.DataFrame, est: pd.DataFrame, gate: Gate) -> list[str]:
    hidden = {(r.table, r.group) for r in counts.itertuples() if r.patients in (NOT_SHOWN, gate.small_text)}
    shown = est[est["quantity"].str.startswith("survival at") & (est["estimate"] != "")]
    return [f"estimates.csv: {r.analysis}, {r.group}, {r.quantity} is shown, but counts.csv hides the group's size"
            for r in shown.itertuples() if (r.analysis, r.group) in hidden]


def check_outputs(folder: Path, gate: Gate, ids: set[str], allowed_dates: set[str]) -> list[str]:
    """Every problem in a folder about to be published: a file that isn't a table, figure or report; an ID from
    the data (`ids`: each as a whole word); a row-level column; a date not in `allowed_dates` (the run's own and
    the plan's); a count column holding a small count; a curve drawn past a hidden number at risk; a survival
    estimate for a hidden group size. Empty means it may be published."""
    problems = []
    folder = Path(folder)

    def read(name):
        return pd.read_csv(folder / name, dtype=str, keep_default_na=False)
    if (folder / "curves.csv").exists() and (folder / "risk_table.csv").exists():
        problems += curves_past_hidden_counts(pd.read_csv(folder / "curves.csv"), read("risk_table.csv"), gate)
    if (folder / "counts.csv").exists() and (folder / "estimates.csv").exists():
        problems += survival_of_hidden_sizes(read("counts.csv"), read("estimates.csv"), gate)
    id_pattern = None
    if ids:
        id_pattern = re.compile(r"(?<![\w-])(?:" + "|".join(re.escape(i) for i in sorted(ids, key=len, reverse=True))
                                + r")(?![\w-])")
    for f in sorted(folder.rglob("*")):
        if f.is_dir():
            continue
        name = f.relative_to(folder).as_posix()
        if f.suffix.lower() not in PUBLISHABLE_SUFFIXES:
            problems.append(f"{name}: not a table, figure or report")
            continue
        text = f.read_text(encoding="utf-8")
        if id_pattern is not None and id_pattern.search(text):
            problems.append(f"{name}: holds an ID from the data")
        for m in DATE.findall(text):
            if m not in allowed_dates:
                problems.append(f"{name}: holds a date that isn't the run's or the plan's")
                break
        if f.suffix.lower() == ".csv":
            t = pd.read_csv(f, dtype=str, keep_default_na=False)
            for c in t.columns:
                if c.lower() in ROW_LEVEL_COLUMNS:
                    problems.append(f"{name}: row-level column {c}")
                if c.lower() in COUNT_COLUMNS:
                    small = [x for x in t[c] if x.strip().isdigit() and gate.small(int(x))]
                    if small:
                        problems.append(f"{name}: {c} shows a count below {gate.min_count}")
    return problems


def plain_ids(ids) -> tuple[set[str], int]:
    """The IDs a text scan can look for: those with a character that isn't a digit (a plain number can't be told
    from a result). Also how many are plain numbers."""
    keep = {str(i) for i in ids if i is not None and str(i).strip() and not str(i).strip().isdigit()}
    return keep, sum(1 for i in ids if i is not None and str(i).strip().isdigit())


__all__ = ["Gate", "GateError", "Relation", "disclosed", "require_safe", "show", "CountSpec", "count_values",
           "publish_counts", "descriptive_values", "publish_descriptives", "hide_across", "risk_table",
           "publish_risk", "drawn_steps", "withheld_curves", "public_curves", "small_group_rule",
           "hide_survival_of_hidden_groups", "add_note", "binary_smd", "pinned_small_cells", "balance_precision",
           "check_outputs", "plain_ids", "NOT_SHOWN", "SIZE_NOTE"]
