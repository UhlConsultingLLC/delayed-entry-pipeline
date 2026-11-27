"""Phase 2 acceptance (CRC BPC only): does the run reproduce the CRC BPC study's results (not yet published)?

    python scripts/acceptance_phase2.py <CRC repository> <the run's output folder>

It compares every table of the CRC study's results (its outputs/ folder, not yet published) with the run's, row by
row, and checks
the numbers its README's Results section quotes against the run's tables. The README's numbers are read from the
README itself when the script runs, so none is written here.

It prints only labels and a verdict for each: "match", "match, worded differently" (the general pipeline's wording
for the same thing, listed below), "not produced" or "pipeline only" (with the reason), or "DIFFERS". It never
prints a count or an estimate. The last line is "RESULT: PASS" when nothing DIFFERS.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

# ---- The general pipeline's wording for the same thing ----
REFERENCE = ("vs wild type", "vs Wild type")
FLAGS = {
    "Survival time withheld, so censored at age 89":
        "Survival withheld by the source's age rule, so censored at the age limit",
    "RAS variant and BRAF V600E together, counted as RAS-mutant":
        "Matched more than one group's rules, so counted in the first by precedence",
    "First report came back before diagnosis, so entry is day 0": "Entry event came before diagnosis, so entry is day 0",
    "First report is a test linked to more than one cancer": "First test is linked to more than one cancer",
    "Another report has a withheld time, so it was ranked last": "Another test has a withheld time, so it was ranked last",
}
TABLES = {"Race and ethnicity": "Race ethnicity"}
LEVELS = {("Histology", "Blank"): "Missing"}
NOT_PRODUCED = {
    "Year of diagnosis worked out from": "the standard tables don't record how a source worked out the year of diagnosis",
    "Specimen type of the first NGS sample": "the standard tables don't carry the specimen type",
}
PIPELINE_ONLY = {
    "First test tied with another after every tie-break, so the lowest test ID was taken":
        "a flag the general cohort code counts (phase 1)",
    "First test has no record of whether it came back after follow-up ended": "a flag the general cohort code counts (phase 1)",
}
LABELS = {
    "positive control: described, not labeled": "all patients: described, not labeled",
    "not compared (D1)": "not compared (smallest group)",
    "described, not labeled (D1)": "described, not labeled (smallest group)",
    "described, not labeled (section 4b)": "described, not labeled (biomarker fallback)",
    "described: corrected risk set never reaches 20 (section 6)": "described: corrected risk set never reaches 20",
}
NOTES = {"SA6: entry and survival may be dependent (Kendall's tau CI excludes 0)":
         "SA6: entry and survival may be dependent (Kendall's tau interval excludes 0)"}
MISSING_INTERVAL = re.compile(r"interval not computed: \d+ resamples gave no value")
D2 = ("eligible patients who died within 730 days of their first report",
      "eligible patients who died within 730 days of their entry event")
FIRST_FLOW_LABEL = ("Patients in the release", "Patients in the source")
WORDING = [
    "A comparison names its reference as the settings do: 'vs Wild type', not 'vs wild type'.",
    "The flags kept but counted, the D2 row, the first flow step and the race table are named in general terms.",
    "Labels name the general rule ('smallest group', 'biomarker fallback'), not the CRC plan's section numbers.",
    "'All patients' is described, not labeled; the pipeline doesn't know it is the CRC study's positive control.",
    "The quasi-independence note says 'interval' for the CRC study's 'CI'.",
    "The figure's title and legend are in general terms; its curves and numbers at risk are compared.",
    "Where a pass 3b interval is blank, the note also says why (resamples with no value), as the CRC study's "
    "passes 1 and 2 rows already did.",
]


def read(folder: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(folder / name, dtype=str, keep_default_na=False)


def ref(text: str) -> str:
    return text.replace(*REFERENCE)


def same_number(a: str, b: str) -> bool:
    if a == b:
        return True
    try:
        x, y = float(a), float(b)
    except ValueError:
        return False
    return abs(x - y) <= 1e-9 * max(1.0, abs(x), abs(y))


class Tally:
    def __init__(self):
        self.differs = 0

    def say(self, what: str, verdict: str):
        if verdict == "DIFFERS":
            self.differs += 1
        print(f"  {what}: {verdict}")


def compare_rows(t: Tally, name: str, crc: pd.DataFrame, run: pd.DataFrame, key: list[str], values: list[str],
                 mapping=None, label_map=None, note_map=None):
    """Each CRC row against the run's row with the same key (after mapping the CRC key to the pipeline's wording)."""
    mapping = mapping or (lambda k: k)
    worded = False
    index = {tuple(r[k] for k in key): r for _, r in run.iterrows()}
    seen = set()
    bad, added = [], []
    for _, r in crc.iterrows():
        k0 = tuple(r[k] for k in key)
        k = mapping(k0)
        if k is None:
            continue
        if k != k0:
            worded = True
        mine = index.get(k)
        if mine is None:
            bad.append(" | ".join(k0) + ": missing from the run")
            continue
        seen.add(k)
        for v in values:
            a, b = r[v], mine[v]
            if v == "label" and label_map:
                worded |= a in label_map
                a = label_map.get(a, a)
            if v == "note" and note_map:
                for old, new in note_map.items():
                    if old in a:
                        a, worded = a.replace(old, new), True
            interval_blank = r.get("lower") == mine.get("lower") == "" and r.get("upper") == mine.get("upper") == ""
            if v == "note" and a != b and interval_blank and \
                    MISSING_INTERVAL.fullmatch(b.replace(a, "", 1).lstrip("; ")):
                added.append(" | ".join(k0))
                continue
            if not same_number(a, b):
                bad.append(" | ".join(k0) + f": {v}")
    extra = [" | ".join(k) for k in index if k not in seen]
    for k in added:
        t.say(f"{name}: {k}: note", "match, with the reason for its missing interval added")
    if bad:
        for b in bad:
            t.say(f"{name}: {b}", "DIFFERS")
    else:
        t.say(f"{name}: every row", "match, worded differently" if worded else "match")
    return extra


def compare_outputs(t: Tally, crc: Path, run: Path):
    print("Tables (outputs/):")
    # flow
    a, b = read(crc, "flow.csv"), read(run, "flow.csv")
    ok = len(a) == len(b) and all(same_number(x, y) for c in ("excluded", "remaining") for x, y in zip(a[c], b[c]))
    labels_ok = list(a["label"][1:]) == list(b["label"][1:]) and (a["label"][0], b["label"][0]) == FIRST_FLOW_LABEL
    t.say("flow.csv: every step, its exclusions and what remains", ("match, worded differently" if labels_ok else "DIFFERS")
          if ok else "DIFFERS")
    # counts
    def counts_key(k):
        table, group = k
        if (table, group) == ("D2", D2[0]):
            return ("D2", D2[1])
        return (table, group)
    compare_rows(t, "counts.csv", read(crc, "counts.csv"), read(run, "counts.csv"), ["table", "group"], ["patients"],
                 counts_key)
    # descriptives
    def desc_key(k):
        table, group = k
        if table in NOT_PRODUCED:
            return None
        if table == "Kept but counted":
            return (table, FLAGS.get(group, group))
        return (TABLES.get(table, table), LEVELS.get((table, group), group))
    a = read(crc, "descriptives.csv")
    for table in sorted(set(a["table"]) & set(NOT_PRODUCED)):
        t.say(f"descriptives.csv: {table}", f"not produced ({NOT_PRODUCED[table]})")
    extra = compare_rows(t, "descriptives.csv", a, read(run, "descriptives.csv"), ["table", "group"], ["patients"], desc_key)
    for e in extra:
        group = e.split(" | ", 1)[1]
        t.say(f"descriptives.csv: {e}", f"pipeline only ({PIPELINE_ONLY[group]})" if group in PIPELINE_ONLY else "DIFFERS")
    compare_rows(t, "risk_table.csv", read(crc, "risk_table.csv"), read(run, "risk_table.csv"), ["pass", "group", "years"],
                 ["at_risk"])
    a, b = read(crc, "curves.csv"), read(run, "curves.csv")
    for g in sorted(set(a["group"]) | set(b["group"])):
        x = a[a["group"] == g].reset_index(drop=True)
        y = b[b["group"] == g].reset_index(drop=True)
        same = len(x) == len(y) and all(same_number(p, q) for c in ("months", "survival") for p, q in zip(x[c], y[c])) \
            and list(x["pass"]) == list(y["pass"])
        t.say(f"curves.csv: {g}", "match" if same else "DIFFERS")
    for g in ("All patients", "BRAF V600E"):
        if g not in set(a["group"]) and g not in set(b["group"]):
            t.say(f"curves.csv: {g}", "match (withheld in both)")
    key = lambda k: tuple(ref(x) for x in k)  # noqa: E731
    compare_rows(t, "estimates.csv", read(crc, "estimates.csv"), read(run, "estimates.csv"),
                 ["analysis", "pass", "quantity", "group"], ["estimate", "lower", "upper", "note"], key, note_map=NOTES)
    compare_rows(t, "bias.csv", read(crc, "bias.csv"), read(run, "bias.csv"), ["analysis", "quantity", "group"],
                 ["estimate", "lower", "upper", "resamples_missing", "bar", "label", "note"], key, LABELS, NOTES)
    a, b = read(crc, "balance.csv"), read(run, "balance.csv")
    for c in ("largest_pair_before", "largest_pair_after"):
        a[c], b[c] = a[c].str.lower(), b[c].str.lower()
    compare_rows(t, "balance.csv", a, b, ["term"], ["largest_smd_before", "largest_smd_after", "after_below_0.1",
                                                     "largest_pair_before", "largest_pair_after"])
    compare_rows(t, "s3_delays.csv", read(crc, "s3_delays.csv"), read(run, "s3_delays.csv"), ["by", "level"],
                 ["n", "q1", "median", "q3"])
    paths = lambda f: re.findall(r'<path d="([^"]+)"', f.read_text(encoding="utf-8"))  # noqa: E731
    texts = lambda f: re.findall(r'text-anchor="middle" fill="#[0-9A-F]+">([^<]*)</text>', f.read_text(encoding="utf-8"))  # noqa: E731
    same = paths(crc / "figure_km.svg") == paths(run / "figure_km.svg") and \
        texts(crc / "figure_km.svg") == texts(run / "figure_km.svg")
    t.say("figure_km.svg: curves, axes and numbers at risk", "match, worded differently" if same else "DIFFERS")


# ---- The README's Results section ----
def readme_claims(readme: str, run: Path) -> list[tuple[str, list[str], list[str]]]:
    """(claim, the README's numbers, the run's numbers formatted the same way) for each number the README quotes."""
    text = readme.split("## Results", 1)[1].split("## How to reproduce", 1)[0].replace("\u2212", "-")
    text = re.sub(r"(\d),(\d{3})\b", r"\1\2", text)            # a thousands separator: 1,234 is 1234
    est, bias, counts = read(run, "estimates.csv"), read(run, "bias.csv"), read(run, "counts.csv")
    flow, s3, bal = read(run, "flow.csv"), read(run, "s3_delays.csv"), read(run, "balance.csv")

    def e(a, p, q, g, col="estimate"):
        r = est[(est["analysis"] == a) & (est["pass"] == p) & (est["quantity"] == q) & (est["group"] == g)]
        return r[col].iloc[0] if len(r) else ""

    def b(a, q, g, col="estimate"):
        r = bias[(bias["analysis"] == a) & (bias["quantity"] == q) & (bias["group"] == g)]
        return r[col].iloc[0] if len(r) else ""

    def c(table, g):
        r = counts[(counts["table"] == table) & (counts["group"] == g)]
        return r["patients"].iloc[0] if len(r) else ""

    def f(x, d=1):
        return "" if x == "" else f"{float(x):.{d}f}"

    def interval(get, d=1):
        return [f(get("estimate"), d), f(get("lower"), d), f(get("upper"), d)]

    claims = []

    def claim(name, pattern, mine):
        m = re.search(pattern, text, flags=re.S)
        claims.append((name, list(m.groups()) if m else ["(not found in the README)"], mine))

    H, RMST, RMSTD, ALL = "HR (robust 95% CI)", "3-year RMST, months", "3-year RMST difference, months", "All patients"
    WT, RAS, BRAF = "Wild type", "RAS-mutant", "BRAF V600E"
    claim("cohort: the release, stage IV, left out at the time check, analyzed",
          r"(\d+) patients are in the release\. (\d+) had stage IV.*?(\d+) of those were left out.*?which leaves (\d+)",
          [flow["remaining"].iloc[0], flow.loc[flow["label"] == "Stage IV at diagnosis", "remaining"].iloc[0],
           flow.loc[flow["label"] == "First report back before follow-up ended", "excluded"].iloc[0], c("flow", "eligible")])
    claim("cohort: the three groups", r"- (\d+) wild type;\s*- (\d+) RAS-mutant;\s*- (\d+) BRAF V600E\.",
          [c("primary", WT), c("primary", RAS), c("primary", BRAF)])
    s = s3[s3["by"] == "All patients"].iloc[0]
    claim("the wait to the first report: median and quartiles",
          r"first report was ([\d.]+) months \(quartiles ([\d.]+) to ([\d.]+)\)", [s["median"], s["q1"], s["q3"]])
    claim("the wait for patients diagnosed in 2012 or earlier", r"2012 or earlier, it was ([\d.]+) months",
          [s3.loc[(s3["by"] == "Year band") & (s3["level"] == "2012 or earlier"), "median"].iloc[0]])
    claim("main finding: the naive median", r"half the patients were alive at ([\d.]+) months",
          [f(e("primary", "naive", "median, months", ALL))])
    claim("main finding: the corrected median", r"half had died by ([\d.]+) months",
          [f(e("primary", "corrected", "median, months", ALL))])
    claim("main finding: the positive control and its interval",
          r"The ([\d.]+)-month difference \(95% CI ([\d.]+) to ([\d.]+)\)",
          interval(lambda col: b("primary", "median_diff", ALL, col)))
    for g in (WT, RAS, BRAF):
        claim(f"RMST table: {g}", rf"\| {g} \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) \(([\d.]+) to ([\d.]+)\) \| ([a-z ]+) \|",
              [f(e("primary", "naive", RMST, g)), f(e("primary", "corrected", RMST, g))]
              + interval(lambda col: b("primary", "rmst_diff", g, col)) + [b("primary", "rmst_diff", g, "label")])
    for g in (RAS, BRAF):
        comp = f"{g} vs {WT}"
        claim(f"comparisons table: {g}, gap in 3-year RMST",
              rf"\| {g} vs wild type: gap in 3-year RMST \| ([-\d.]+) months \| ([-\d.]+) months \| understated by "
              r"([-\d.]+) \(([-\d.]+) to ([-\d.]+)\) \| ([a-z ]+) \|",
              [f(e("primary", "naive", RMSTD, comp)), f(e("primary", "corrected", RMSTD, comp))]
              + interval(lambda col: b("primary", "rmst_contrast_change", comp, col))
              + [b("primary", "rmst_contrast_change", comp, "label")])
        claim(f"comparisons table: {g}, hazard ratio",
              rf"\| {g} vs wild type: hazard ratio \| ([\d.]+) \| ([\d.]+) \| ratio ([\d.]+) \(([\d.]+) to ([\d.]+)\) "
              r"\| ([a-z ]+) \|",
              [f(e("primary", "naive", H, comp), 2), f(e("primary", "corrected", H, comp), 2)]
              + interval(lambda col: b("primary", "hr_ratio", comp, col), 2) + [b("primary", "hr_ratio", comp, "label")])
    claim("pass 3b: RAS-mutant's understatement", r"pass 3b\), the understatement is ([\d.]+) months \(([\d.]+) to ([\d.]+)\)",
          interval(lambda col: b("3b", "rmst_contrast_change", f"{RAS} vs {WT}", col)))
    claim("pass 3a: the hazard-ratio ratios",
          r"the ratios are ([\d.]+) \(([\d.]+) to ([\d.]+)\) for RAS-mutant and ([\d.]+) \(([\d.]+) to ([\d.]+)\) for BRAF",
          interval(lambda col: b("3a", "hr_ratio", f"{RAS} vs {WT}", col), 2)
          + interval(lambda col: b("3a", "hr_ratio", f"{BRAF} vs {WT}", col), 2))
    claim("pass 3a: the primary adjusted hazard ratios",
          r"pass 3a\) are ([\d.]+) \(([\d.]+) to ([\d.]+)\) for RAS-mutant and ([\d.]+) \(([\d.]+) to ([\d.]+)\) for BRAF",
          [f(e("3a", "3a corrected", H, f"{RAS} vs {WT}", col), 2) for col in ("estimate", "lower", "upper")]
          + [f(e("3a", "3a corrected", H, f"{BRAF} vs {WT}", col), 2) for col in ("estimate", "lower", "upper")])
    claim("SA1: its patients", r"diagnosis: (\d+) wild type, (\d+) RAS-mutant, (\d+) BRAF V600E",
          [c("SA1", WT), c("SA1", RAS), c("SA1", BRAF)])
    claim("SA1: the largest overstatement", r"RAS-mutant's 3-year RMST, at ([\d.]+) months",
          [f(b("SA1", "rmst_diff", RAS), 2)])
    claim("SA1: the corrected median", r"median survival from first report was ([\d.]+) months",
          [f(e("SA1", "corrected", "median, months", ALL))])
    claim("S2: all patients' RMST overstatement, beside the primary's", r"overstated by ([\d.]+) months against ([\d.]+)",
          [f(b("S2", "rmst_diff", ALL)), f(b("primary", "rmst_diff", ALL))])
    gaps = [float(b(a, "rmst_contrast_change", f"{RAS} vs {WT}")) for a in ("SA3", "SA4", "SA7")]
    claim("SA3, SA4, SA7: RAS-mutant's gap understated", r"still understated, by ([\d.]+) to ([\d.]+) months",
          [f"{min(gaps):.1f}", f"{max(gaps):.1f}"])
    claim("SA5: the Schoenfeld test", r"Schoenfeld test gives p = ([\d.]+)",
          [f(e("SA5", "corrected", "Schoenfeld test p (global)", ALL), 2)])
    claim("stability rule: RAS-mutant", r"stays clearly material, at ([\d.]+) months of 3-year RMST \(([\d.]+) to ([\d.]+)\)",
          interval(lambda col: b("conditional on t*", "rmst_diff", RAS, col)))
    claim("stability rule: BRAF V600E's RMST", r"shrinks to ([\d.]+) months \(([\d.]+) to ([\d.]+)\)",
          interval(lambda col: b("conditional on t*", "rmst_diff", BRAF, col)))
    claim("stability rule: BRAF V600E's median", r"rises to ([\d.]+) months",
          [f(b("conditional on t*", "median_diff", BRAF))])
    claim("SA6: tau for all patients", r"tau is ([-\d.]+) \(95% CI ([-\d.]+) to ([-\d.]+)\)",
          [f(e("SA6", "corrected", "conditional Kendall's tau", ALL, col), 3) for col in ("estimate", "lower", "upper")])
    claim("BRAF V600E: its size", r"BRAF V600E is small \((\d+) patients\)", [c("primary", BRAF)])
    claim("BRAF V600E: its hazard-ratio bias, then adjusted", r"hazard ratio \(([\d.]+)\).*?drops to ([\d.]+) after adjustment",
          [f(b("primary", "hr_ratio", f"{BRAF} vs {WT}"), 2), f(b("3a", "hr_ratio", f"{BRAF} vs {WT}"), 2)])
    un, mss, msi = int(c("S1", "Unknown")), int(c("S1", "MSS/pMMR")), int(c("S1", "MSI-H/dMMR"))
    claim("S1: the share with no result", r"(\d+)% of S1 patients had no MSI", [f"{100 * un / (un + mss + msi):.0f}"])
    claim("S1: MSI-H/dMMR patients", r"Only (\d+) patients are MSI-H/dMMR", [c("S1", "MSI-H/dMMR")])
    after = bal["largest_smd_after"].astype(float)
    claim("balance: terms at or above 0.1 after weighting", r"After weighting, (\d+) of (\d+) covariate terms",
          [str(int((after >= 0.1).sum())), str(len(bal))])
    claim("balance: terms above 0.2", r"and (two) by more than 0.2", ["two" if int((after > 0.2).sum()) == 2 else "other"])
    return claims


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("crc_repository", type=Path)
    ap.add_argument("run_outputs", type=Path)
    a = ap.parse_args(argv)
    t = Tally()
    compare_outputs(t, a.crc_repository / "outputs", a.run_outputs)
    print("README Results:")
    for name, theirs, mine in readme_claims((a.crc_repository / "README.md").read_text(encoding="utf-8"), a.run_outputs):
        t.say(f"README: {name}", "match" if theirs == mine else "DIFFERS")
    print("Wording that differs on purpose:")
    for w in WORDING:
        print(f"  - {w}")
    print("RESULT: PASS" if t.differs == 0 else "RESULT: FAIL (see each DIFFERS above)")
    return 0 if t.differs == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
