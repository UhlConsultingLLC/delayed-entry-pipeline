"""The run report (design plan section 10): what was done and why, in Markdown and HTML.

It holds the inputs and their checksums, the attribution the source's terms require (when its capability profile
gives one), the capability profile, the readiness screen's three lists from the
locked plan, the plan's hash, each check with its verdict (the lock, the checks that needed outcomes, the R steps,
the safety gate), and a plain summary of the results with their labels. It is built from the published tables
only, after the gate, so it can't hold a number the tables don't.
"""
from __future__ import annotations

import html

import pandas as pd

ALL = "All patients"


class Doc:
    """A small document: headings, paragraphs, lists and tables, rendered as Markdown or HTML."""

    def __init__(self, title: str):
        self.title = title
        self.parts: list[tuple] = []

    def h(self, text: str):
        self.parts.append(("h", text))

    def p(self, text: str):
        self.parts.append(("p", text))

    def ul(self, items: list[str]):
        if items:
            self.parts.append(("ul", list(items)))

    def table(self, header: list[str], rows: list[list]):
        self.parts.append(("table", header, [[("" if v is None else str(v)) for v in r] for r in rows]))

    def markdown(self) -> str:
        out = [f"# {self.title}", ""]
        for part in self.parts:
            if part[0] == "h":
                out += [f"## {part[1]}", ""]
            elif part[0] == "p":
                out += [part[1], ""]
            elif part[0] == "ul":
                out += [f"- {x}" for x in part[1]] + [""]
            else:
                _, header, rows = part
                cell = lambda v: str(v).replace("|", "\\|").replace("\n", " ")  # noqa: E731
                out += ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
                out += ["| " + " | ".join(cell(v) for v in r) + " |" for r in rows] + [""]
        return "\n".join(out).rstrip() + "\n"

    def html(self) -> str:
        e = html.escape
        out = ["<!doctype html>", '<html lang="en">', "<head>", '<meta charset="utf-8">',
               '<meta name="viewport" content="width=device-width, initial-scale=1">', f"<title>{e(self.title)}</title>",
               "<style>",
               ":root { --ink: #15202B; --muted: #56626E; --line: #D2D9E0; --bg: #FFFFFF; }",
               "@media (prefers-color-scheme: dark) { :root { --ink: #E6EAEE; --muted: #A9B3BD; --line: #3A4550; --bg: #12181E; } }",
               "body { font-family: Arial, Helvetica, sans-serif; color: var(--ink); background: var(--bg); "
               "max-width: 1100px; margin: 2rem auto; padding: 0 16px; line-height: 1.45; }",
               "h1 { font-size: 1.6rem; } h2 { font-size: 1.2rem; margin-top: 2rem; border-bottom: 1px solid var(--line); }",
               "table { border-collapse: collapse; display: block; overflow-x: auto; margin: 0.5rem 0 1rem; }",
               "th, td { border: 1px solid var(--line); padding: 4px 8px; text-align: left; vertical-align: top; font-size: 0.9rem; }",
               "th { color: var(--muted); }", "</style>", "</head>", "<body>", f"<h1>{e(self.title)}</h1>"]
        for part in self.parts:
            if part[0] == "h":
                out.append(f"<h2>{e(part[1])}</h2>")
            elif part[0] == "p":
                out.append(f"<p>{e(part[1])}</p>")
            elif part[0] == "ul":
                out.append("<ul>" + "".join(f"<li>{e(x)}</li>" for x in part[1]) + "</ul>")
            else:
                _, header, rows = part
                out.append("<table><thead><tr>" + "".join(f"<th>{e(h)}</th>" for h in header) + "</tr></thead><tbody>")
                out += ["<tr>" + "".join(f"<td>{e(v)}</td>" for v in r) + "</tr>" for r in rows]
                out.append("</tbody></table>")
        out += ["</body>", "</html>"]
        return "\n".join(out) + "\n"


def _fmt(x, digits=1) -> str:
    if x is None or (isinstance(x, float) and pd.isna(x)) or x == "":
        return ""
    try:
        return f"{float(x):.{digits}f}"
    except (TypeError, ValueError):
        return str(x)


def _interval(row, digits=1) -> str:
    lo, hi = _fmt(row.get("lower"), digits), _fmt(row.get("upper"), digits)
    return f"{lo} to {hi}" if lo and hi else "no interval"


def results_summary(est: pd.DataFrame, bias: pd.DataFrame, horizon: float, level: float) -> tuple[list[str], list[list]]:
    """A plain summary of the primary analysis, from the published tables: each group's survival on both clocks
    and how much the naive clock overstates it, and how the comparisons change."""
    lines, rows = [], []
    e = est[est["analysis"] == "primary"]

    def val(p, q, g):
        x = e[(e["pass"] == p) & (e["quantity"] == q) & (e["group"] == g)]
        return x.iloc[0].to_dict() if len(x) else {}

    rmst_q = f"{horizon:g}-year RMST, months"
    b = bias[bias["analysis"] == "primary"]
    for g in e["group"].drop_duplicates():
        if " vs " in g:
            continue
        rn, rc = val("naive", rmst_q, g), val("corrected", rmst_q, g)
        d = b[(b["quantity"] == "rmst_diff") & (b["group"] == g)]
        if rn and rc and len(d):
            d = d.iloc[0].to_dict()
            rows.append([g, _fmt(rn.get("estimate")), _fmt(rc.get("estimate")),
                         f"{_fmt(d.get('estimate'))} ({_interval(d)})", d.get("label", "")])
    for _, r in b[b["quantity"].isin(["hr_ratio", "rmst_contrast_change"])].iterrows():
        what = "the ratio of the naive to the corrected hazard ratio" if r["quantity"] == "hr_ratio" else \
            f"the change in the {horizon:g}-year RMST difference, in months"
        lines.append(f"{r['group']}: {what} is {_fmt(r['estimate'], 2)} ({level:.0%} interval {_interval(r, 2)}); "
                     f"{r['label']}.")
    return lines, rows


def build_report(info: dict, tables: dict[str, pd.DataFrame], gate_log: list[dict], notes: list[str]) -> Doc:
    """The run report from the run's record (`info`), the published tables, the gate's log and the notes."""
    doc = Doc(f"Run report: {info['study']}")
    doc.p("What was run, on what, and why each check came out as it did. Every number below is from the published "
          f"tables, after the safety gate: counts from 1 to {info['min_count'] - 1} are shown as \"<{info['min_count']}\".")
    doc.h("The run")
    doc.table(["Item", "Value"], [
        ["Study", info["study"]], ["Registered plan", info.get("registered_plan", "")],
        ["Settings file", f"{info['settings_file']} (SHA-256 {info['settings_sha256']})"],
        ["Locked plan", f"version {info['plan_version']}, locked {info['locked_at']} by {info['locked_by']}"],
        ["Plan files' SHA-256", "; ".join(f"{k}: {v}" for k, v in info["plan_files"].items())],
        ["Tool", info["tool"]], ["Run at (UTC)", info["run_at"]],
        ["R", "; ".join(f"{k} {v}" for k, v in info.get("r_versions", {}).items())],
        ["Bootstrap", f"{info['resamples']} resamples, seed {info['seed']}"],
    ])
    doc.h("Inputs")
    doc.p("Each input file, with its SHA-256 as the lock recorded it. The run stopped unless every one still matched.")
    doc.table(["File", "SHA-256"], [[i["file"], i["sha256"]] for i in info["inputs"]])
    attribution = info["profile"].get("attribution")
    if attribution:
        doc.h("Data source and attribution")
        doc.p("The source's terms require work that uses or references these data to carry this attribution.")
        blank = "NOT FILLED IN: add it to the capability profile's attribution before publishing"
        doc.ul([f"Citation: {attribution.get('citation') or blank}",
                f"Dataset version: {attribution.get('version') or blank}",
                f"Acknowledgment: {attribution.get('acknowledgment') or blank}"])
    doc.h("Capability profile of the source")
    doc.table(["Item", "Value"], [[k, v] for k, v in info["profile"].items() if k not in ("sources", "attribution")])
    doc.h("The readiness screen, from the locked plan")
    for name, key in (("Ready", "ready"), ("Missing", "missing"), ("Deviations required", "deviations_required"),
                      ("Warnings", "warnings")):
        items = info["screen"].get(key, [])
        doc.p(f"{name} ({len(items)})" + ("" if items else ": none."))
        doc.ul([f"{('[' + i['severity'] + '] ') if i['severity'] else ''}{i['analysis']}: {i['topic']}: {i['finding']}"
                + (f" Effect: {i['effect']}" if i["effect"] else "") for i in items])
    doc.h("Checks, with their verdicts")
    doc.table(["Stage", "Check", "Verdict", "Detail"], [[c["stage"], c["check"], c["verdict"], c["detail"]]
                                                        for c in info["checks"]])
    doc.h("The safety gate")
    doc.p("Nothing reached the output folder before every check below passed. The tables were built in the run's own "
          "folder, read back, and only then copied here.")
    doc.table(["Check", "Verdict", "Detail"], [[c["check"], c["verdict"], c["detail"]] for c in gate_log])
    doc.h("What each analysis did")
    doc.table(["Analysis", "What it is", "Status", "Why"], [[a["id"], a["label"], a["status"], a["why"]]
                                                            for a in info["analyses"]])
    if notes:
        doc.h("Notes")
        doc.ul(notes)
    doc.h("Results in brief")
    lines, rows = results_summary(tables["estimates.csv"], tables["bias.csv"], info["horizon"], info["interval_level"])
    doc.p(f"Months alive over the first {info['horizon']:g} years (RMST), counted from diagnosis (naive) and from the "
          "entry event with delayed entry (corrected), and how much the naive clock overstates them. The full results, "
          "with every analysis, are in estimates.csv and bias.csv.")
    doc.table(["Group", "Naive", "Corrected", f"Overstated by ({info['interval_level']:.0%} interval)", "Label"], rows)
    doc.ul(lines)
    doc.h("Files")
    doc.ul([f"{name}: {what}" for name, what in info["files"].items()])
    return doc


__all__ = ["Doc", "build_report", "results_summary"]
