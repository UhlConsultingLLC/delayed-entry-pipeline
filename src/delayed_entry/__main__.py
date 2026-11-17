"""Command line: python -m delayed_entry <command> ...

    template <out.yaml>                        write the blank settings template
    check-settings <settings>                  load a settings file and list every problem, with its line
    lock-settings <settings> <folder> --locked-by NAME
                                               lock the settings alone, before any data are downloaded
    draft-mapping <out.yaml> TABLE=FILE ...    draft a flat-file mapping (a person must confirm it)
    omop-template <out.yaml>                   write a blank OMOP CDM mapping (a person fills it in and confirms it)
    screen <settings> <work> [source]          read the source and write the readiness report
    lock <settings> <work> --locked-by NAME [source] [--reason WHY]
                                               screen, then lock the plan (a person has reviewed it)
    verify-lock <manifest.json>                say whether anything changed since the lock
    run <settings> <work> <out> [source] --r PROGRAM [--r-library DIR] [--manifest FILE] [--reuse-r]
                                               the run, the safety gate and the reports, after the lock

source: --clinical DIR [--genomic DIR] for GENIE BPC, or --mapping FILE for flat files or an OMOP CDM database
(the settings' source.adapter says which). The work folder must be on local disk, not in a cloud-synced folder. The run's
output folder must be outside the work folder. R: --r names Rterm.exe (run as `Rterm -f`) or Rscript; it
defaults to the environment variable DELAYED_ENTRY_R, then to Rscript on the PATH. --r-library puts a folder of
R packages first on R's library path (default: the environment variable DELAYED_ENTRY_R_LIBRARY).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from delayed_entry.lock import lock_plan, lock_settings, verify_lock
from delayed_entry.privacy import check_not_cloud
from delayed_entry.screen import screen
from delayed_entry.settings import SettingsError, load_settings, write_blank_template
from delayed_entry.tables import outcome_blind

REPORT_MD, REPORT_YAML = "readiness_report.md", "readiness_report.yaml"


def read_source(settings, a):
    from delayed_entry.run import read_source as read_
    return read_(settings, a.clinical, a.genomic, a.mapping)


def run_screen(a):
    settings = load_settings(a.settings)
    work = check_not_cloud(Path(a.work).resolve())
    work.mkdir(parents=True, exist_ok=True)
    result = read_source(settings, a)
    report = screen(outcome_blind(result.tables), settings, result.profile, result.source_steps, result.input_files)
    (work / REPORT_MD).write_bytes(report.to_markdown().encode("utf-8"))
    (work / REPORT_YAML).write_bytes(report.to_yaml().encode("utf-8"))
    return settings, result, report, work


def summary(report) -> str:
    lines = [f"ready {len(report.list('ready'))}, missing {len(report.list('missing'))}, "
             f"deviations required {len(report.list('deviations required'))}, warnings {len(report.warnings)}"]
    for analysis, comps in report.decisions.items():
        for comp, d in comps.items():
            lines.append(f"  {analysis}: {comp}: {d}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m delayed_entry", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="command", required=True)
    t = sub.add_parser("template")
    t.add_argument("out")
    c = sub.add_parser("check-settings")
    c.add_argument("settings")
    ls = sub.add_parser("lock-settings")
    ls.add_argument("settings")
    ls.add_argument("folder")
    ls.add_argument("--locked-by", required=True)
    d = sub.add_parser("draft-mapping")
    d.add_argument("out")
    d.add_argument("files", nargs="+", help="TABLE=FILE or TABLE=FILE#SHEET")
    o = sub.add_parser("omop-template")
    o.add_argument("out")
    for name in ("screen", "lock"):
        p = sub.add_parser(name)
        p.add_argument("settings")
        p.add_argument("work")
        p.add_argument("--clinical")
        p.add_argument("--genomic")
        p.add_argument("--mapping")
        if name == "lock":
            p.add_argument("--locked-by", required=True)
            p.add_argument("--reason")
    v = sub.add_parser("verify-lock")
    v.add_argument("manifest")
    r = sub.add_parser("run")
    r.add_argument("settings")
    r.add_argument("work")
    r.add_argument("out")
    r.add_argument("--clinical")
    r.add_argument("--genomic")
    r.add_argument("--mapping")
    r.add_argument("--manifest", help="the locked plan's manifest.json (default: the latest plan_v<n> in the work folder)")
    r.add_argument("--r", default=os.environ.get("DELAYED_ENTRY_R", "Rscript"), help="Rterm.exe or Rscript")
    r.add_argument("--r-library", default=os.environ.get("DELAYED_ENTRY_R_LIBRARY"))
    r.add_argument("--reuse-r", action="store_true",
                   help="keep the R results in the run folder when their inputs are byte for byte the same")
    a = ap.parse_args(argv)

    if a.command == "template":
        print(write_blank_template(a.out))
    elif a.command == "check-settings":
        try:
            s = load_settings(a.settings)
        except SettingsError as e:
            print(e)
            sys.exit(1)
        defaults = s.defaults_used()
        print(f"settings OK: {len(s.provenance())} choices, {len(defaults)} taken from the defaults file")
        for d_ in defaults:
            print(f"  default: {d_['choice']} = {d_['value']} ({d_['source']})")
    elif a.command == "lock-settings":
        load_settings(a.settings)                   # never lock a file that doesn't load
        print("settings locked:", lock_settings(a.settings, a.folder, a.locked_by))
    elif a.command == "draft-mapping":
        from delayed_entry.adapters import flatfile
        files = {}
        for item in a.files:
            table, _, spec = item.partition("=")
            path, _, sheet = spec.partition("#")
            files[table] = {"path": path, "sheet": sheet} if sheet else path
        print(flatfile.write_draft(files, a.out))
    elif a.command == "omop-template":
        from delayed_entry.adapters import omop
        print(omop.write_template(a.out))
    elif a.command == "screen":
        _, _, report, work = run_screen(a)
        print(summary(report))
        print("report:", work / REPORT_MD)
    elif a.command == "lock":
        settings, result, report, work = run_screen(a)
        res = lock_plan(report, settings, result.input_files, work, a.locked_by, a.reason, result.read_notes)
        print(summary(report))
        print(f"locked: version {res.version}, manifest {res.manifest}")
    elif a.command == "verify-lock":
        problems = verify_lock(a.manifest)
        print("lock matches" if not problems else "CHANGED: " + "; ".join(problems))
        if problems:
            sys.exit(1)
    elif a.command == "run":
        from delayed_entry.run import ROptions, run
        summary_ = run(a.settings, a.work, a.out, clinical=a.clinical, genomic=a.genomic, mapping=a.mapping,
                       manifest=a.manifest, r=ROptions(a.r, a.r_library, a.reuse_r))
        for c in summary_["checks"]:
            print(f"  {c['stage']}: {c['check']}: {c['verdict']}")
        for c in summary_["gate"]:
            print(f"  gate: {c['check']}: {c['verdict']}")
        print("published:", ", ".join(summary_["files"]))


if __name__ == "__main__":
    main()
