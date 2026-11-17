"""The run (design plan stages 4 to 6): the analysis, the safety gate and the reports, after the plan is locked.

    python -m delayed_entry run <settings> <work folder> <output folder> --clinical DIR [--genomic DIR] --r PATH

In order, it:
1. checks the lock first (delayed_entry.lock.require_lock) and stops if the plan, the settings, the defaults or
   any input has changed;
2. screens again, without outcomes, and stops unless the screen still gives the locked plan's decisions;
3. reads the outcomes and runs the checks that need them (delayed_entry.analysis): time at risk from the times,
   survival after the entry event, the horizon against follow-up, and the smallest-group rule again;
4. writes the run's row-level files to a folder inside the work folder, and runs the R part there: passes 1 and 2
   on both clocks, the stability rule, pass 3a, pass 3b with its balance, the assumption checks, and the secondary
   and sensitivity analyses the locked plan runs (R for the Cox models' robust standard errors and the bootstrap);
5. builds every results table, the figure and the run report through the safety gate (delayed_entry.gate) in that
   folder, reads them all back, and only then copies them to the output folder.

The work folder holds row-level files: it must be local, and the output folder must be outside it.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd
import yaml

from delayed_entry import __version__
from delayed_entry.analysis import RUN, AnalysisData, RunStop, build_analysis, compare_with_plan, keep_plan_decisions
from delayed_entry.gate import DATE, Gate, GateError, check_outputs, plain_ids
from delayed_entry.hashing import sha256_file
from delayed_entry.lock import MANIFEST, PLAN_MD, PLAN_YAML, LockError, require_lock, verify_lock
from delayed_entry.privacy import check_not_cloud
from delayed_entry.reports import FLOAT_FORMAT, Results, build_tables
from delayed_entry.run_report import build_report
from delayed_entry.screen import BLOCK, PRIMARY, screen
from delayed_entry.settings import Settings, load_settings
from delayed_entry.tables import outcome_blind

R_DIR = Path(__file__).with_name("R")
R_SCRIPT = R_DIR / "run_analysis.R"
R_STEPS = ("primary", "adjusted", "weighted", "assumptions", "secondary")
RUN_SPEC = "run_spec.json"
STAGING = "staging"
OUTPUTS_MANIFEST = "outputs_manifest.json"
PACKAGE_ROOT = Path(__file__).resolve().parents[2]
ALL = "All patients"
FILES = {
    "flow.csv": "the cohort's steps, merged so that no small count is shown or can be worked out",
    "counts.csv": "each analysis's group sizes, and the survival-after-entry count",
    "descriptives.csv": "the patients kept but counted, and the distributions described",
    "risk_table.csv": "patients at risk at each fixed time, by pass and group",
    "curves.csv": "the Kaplan-Meier steps that may be drawn, in months",
    "estimates.csv": "every estimate: medians, RMST, survival at fixed times, hazard ratios, the assumption checks",
    "bias.csv": "every naive-minus-corrected contrast, with its interval, bar and label",
    "balance.csv": "pass 3b's covariate balance, before and after weighting",
    "s3_delays.csv": "the wait from diagnosis to the entry event",
    "figure_km.svg": "the naive and corrected curves by group, with the numbers at risk",
    "run_report.md": "this report",
    "run_report.html": "this report, as a web page",
    OUTPUTS_MANIFEST: "the SHA-256 of every input and every output",
}


@dataclass
class ROptions:
    """How to run R: the program (Rterm.exe, run as `Rterm -f`, or Rscript), and a library folder to put first on
    R's library path (optional)."""
    program: str
    library: str | None = None
    reuse: bool = False          # keep the R results already in the run folder if its inputs are byte for byte the same


def r_command(program: str, step: str, run_dir: Path) -> list[str]:
    name = Path(program).name.lower()
    args = [str(R_DIR), str(run_dir), step]
    if name.startswith("rterm"):
        return [program, "--no-echo", "--no-restore", "-f", str(R_SCRIPT), "--args", *args]
    return [program, str(R_SCRIPT), *args]


def run_r_step(opts: ROptions, step: str, run_dir: Path, log=print) -> dict:
    env = os.environ.copy()
    if opts.library:
        env["DELAYED_ENTRY_R_LIBRARY"] = opts.library
    log(f"R step {step}: running")
    proc = subprocess.run(r_command(opts.program, step, run_dir), cwd=run_dir, env=env, capture_output=True,
                          text=True, encoding="utf-8", errors="replace")
    (run_dir / f"r_{step}.log").write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")
    if proc.returncode != 0:
        raise RunStop(f"R step {step} failed (exit code {proc.returncode}); its log is r_{step}.log in the run folder")
    return json.loads((run_dir / f"r_{step}.json").read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# The spec and the files R reads
# ---------------------------------------------------------------------------
def labels(settings: Settings) -> dict:
    return {
        "clearly_material": "clearly material",
        "material": "material",
        "not_material": "not material",
        "not_compared": "not compared (smallest group)",
        "median_not_reached": f"median not reached within {settings.read('summary.median_limit_years'):g} years: "
                              "RMST is the summary",
        "unstable": f"described: corrected risk set never reaches {settings.read('stability.risk_set')}",
        "described": "described, not labeled (%s)",
        "all_patients": "all patients: described, not labeled",
        "no_death_before_t_star": "described: no death before t*",
        "imperfectly_balanced": "imperfectly balanced",
    }


def build_spec(data: AnalysisData, settings: Settings) -> dict:
    """Every value the R part uses, from the locked settings and the run's decisions."""
    terms = data.coded.terms
    primary = data.analysis(PRIMARY)
    secondary = []
    for a in data.analyses:
        if a.status != RUN or a.kind not in ("subset", "treatment_clock", "biomarker", "biomarker_variant",
                                              "adjusted_variant"):
            continue
        secondary.append({"id": a.id, "kind": a.kind, "member_column": a.member_column,
                          "grouping": a.grouping.to_dict() if a.grouping else None,
                          "described": dict(a.described), "stratify_by": a.stratify_by})
    return {
        "days_per_year": settings.read("time.days_per_year"),
        "survival_at_years": list(settings.read("summary.survival_at_years")),
        "horizon_years": settings.read("time.horizon_years"),
        "median_limit_years": settings.read("summary.median_limit_years"),
        "all_patients": ALL,
        "resamples": settings.read("uncertainty.bootstrap_resamples"),
        "seed": settings.read("uncertainty.seed"),
        "interval_level": settings.read("uncertainty.interval_level"),
        "bars": {"median_diff": settings.read("materiality.median_months"),
                 "rmst_diff": settings.read("materiality.rmst_months"),
                 "hr_ratio": settings.read("materiality.hr_ratio_change"),
                 "rmst_contrast_change": settings.read("materiality.rmst_gap_change_months")},
        "stable_risk_set": settings.read("stability.risk_set"),
        "genotype": data.genotype.to_dict(),
        "d1": data.d1.to_dict("records"),
        "primary_described": dict(primary.described),
        "labels": labels(settings),
        "adjustment": {"fields": [t.field for t in terms], "forms": [t.form for t in terms],
                       "levels": {t.field: list(t.levels) for t in terms if t.form in ("categorical", "bands")},
                       "knots": {t.field: list(t.knots_at_percentiles) for t in terms if t.form == "spline"}},
        "weighting": {"trim_percentiles": list(settings.read("weighting.trim_percentiles")),
                      "balance_target_smd": settings.read("weighting.balance_target_smd")},
        "ph": {"p_threshold": settings.read("ph_check.p_threshold"), "split_months": settings.read("ph_check.split_months")},
        "tau_estimators": list(settings.read("assumption_check.estimators")),
        "missing_data": {"rule": data.covariate_plan["missing_data_rule"],
                         "imputations": settings.read("missing_data.imputations"),
                         "model": (settings.read("missing_data.imputation_model")
                                   if settings.has("missing_data.imputation_model") else None)},
        "run_ph_check": any(a.kind == "ph_check" and a.status == RUN for a in data.analyses),
        "run_quasi_independence": any(a.kind == "quasi_independence" and a.status == RUN for a in data.analyses),
        "secondary": secondary,
    }


def plain(df: pd.DataFrame) -> pd.DataFrame:
    """Columns as plain types for Parquet: numbers as float64 (NaN for blank), yes/no as bool, text as object."""
    out = pd.DataFrame(index=df.index)
    for c in df.columns:
        s = df[c]
        if pd.api.types.is_bool_dtype(s):
            out[c] = s.astype("boolean").fillna(False).astype(bool) if s.isna().any() else s.astype(bool)
        elif pd.api.types.is_numeric_dtype(s):
            out[c] = pd.to_numeric(s, errors="coerce").astype("float64")
        else:
            out[c] = s.astype(object).where(s.notna(), None)
    return out.reset_index(drop=True)


def write_parquet(df: pd.DataFrame, path: Path) -> None:
    con = duckdb.connect()
    try:
        con.register("frame", plain(df))
        con.execute(f"COPY frame TO '{path.as_posix()}' (FORMAT parquet)")
    finally:
        con.close()


def write_run_files(data: AnalysisData, spec: dict, run_dir: Path) -> list[Path]:
    cohort = data.cohort[["patient_id", "genotype", "report_day", "entry_day", "exit_day", "event", "in_analysis"]]
    files = {"analysis_cohort.parquet": cohort, "covariates.parquet": data.covariates.sort_values("patient_id"),
             "secondary.parquet": data.secondary}
    for a, clock in data.clocks.items():
        files[f"clock_{a}.parquet"] = clock.sort_values("patient_id", kind="mergesort")
    written = []
    for name, df in files.items():
        write_parquet(df, run_dir / name)
        written.append(run_dir / name)
    (run_dir / RUN_SPEC).write_bytes((json.dumps(spec, indent=2) + "\n").encode("utf-8"))
    written.append(run_dir / RUN_SPEC)
    return written


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------
def latest_manifest(work: Path) -> Path:
    versions = sorted((int(p.name.split("_v")[1]), p) for p in work.glob("plan_v*")
                      if p.is_dir() and p.name.split("_v")[1].isdigit())
    if not versions:
        raise LockError(f"no locked plan in {work}: lock the plan first (python -m delayed_entry lock)")
    return versions[-1][1] / MANIFEST


def read_source(settings: Settings, clinical=None, genomic=None, mapping=None):
    if settings.read("source.adapter") == "genie_bpc":
        from delayed_entry.adapters import genie_bpc
        if clinical is None:
            raise SystemExit("GENIE BPC needs --clinical (and --genomic if the files are in two folders)")
        return genie_bpc.read(clinical, genomic, settings)
    if mapping is None:
        raise SystemExit("flat files and OMOP need --mapping")
    if settings.read("source.adapter") == "omop":
        from delayed_entry.adapters import omop
        return omop.read(mapping, settings)
    from delayed_entry.adapters import flatfile
    return flatfile.read(mapping, settings)


def _check_folders(work: Path, out: Path) -> None:
    if out == work or work in out.parents:
        raise ValueError("the output folder must not be inside the work folder, which holds row-level files")
    if out == PACKAGE_ROOT or PACKAGE_ROOT in out.parents or work == PACKAGE_ROOT or PACKAGE_ROOT in work.parents:
        raise ValueError("results and work files must stay outside the pipeline's own folder, which is published")
    if out.exists() and any(out.iterdir()) and not (out / OUTPUTS_MANIFEST).exists():
        raise ValueError("the output folder holds other files: give an empty folder, or one a run wrote before")


def _same_bytes(a: Path, b: Path) -> bool:
    return a.exists() and b.exists() and sha256_file(a) == sha256_file(b)


def run(settings_path, work, out, *, clinical=None, genomic=None, mapping=None, manifest=None,
        r: ROptions, log=print, now: datetime | None = None) -> dict:
    """The whole run. Returns a summary (labels and verdicts only). Raises LockError, RunStop or GateError when it
    stops; nothing reaches the output folder then."""
    started = (now or datetime.now(timezone.utc)).astimezone(timezone.utc)
    stamp = started.isoformat(timespec="seconds")
    settings = load_settings(settings_path)
    work = check_not_cloud(Path(work).resolve())
    out = check_not_cloud(Path(out).resolve())
    _check_folders(work, out)
    manifest = Path(manifest) if manifest else latest_manifest(work)
    checks: list[dict] = []

    # 1. The lock, before anything else
    problems = verify_lock(manifest, settings)
    if problems:
        raise LockError("the locked plan doesn't match: " + "; ".join(problems))
    result = read_source(settings, clinical, genomic, mapping)
    require_lock(manifest, settings, result.input_files)
    m = json.loads(manifest.read_text(encoding="utf-8"))
    checks.append({"stage": "lock", "check": "the locked plan, settings, defaults and inputs are unchanged",
                   "verdict": "pass", "detail": f"plan version {m['plan_version']}: every SHA-256 matches"})
    log("lock: matches")
    if result.profile.attribution is not None:
        gaps = result.profile.attribution_gaps()
        checks.append({"stage": "source", "check": "the attribution the source's terms require",
                       "verdict": "warning" if gaps else "pass",
                       "detail": (f"blank in the capability profile: {', '.join(gaps)}; fill it in before publishing"
                                  if gaps else "printed in the run report")})
        log(f"attribution: {'blank: ' + ', '.join(gaps) if gaps else 'complete'}")
    plan = yaml.safe_load((manifest.parent / PLAN_YAML).read_text(encoding="utf-8"))
    screened = plan["readiness_screen"]

    # 2. The screen again, without outcomes
    report = screen(outcome_blind(result.tables), settings, result.profile, result.source_steps, result.input_files)
    if report.decisions != screened["decisions"] or report.group_sizes != screened["group_sizes"]:
        raise RunStop("the screen no longer gives the locked plan's decisions: a new plan version is needed")
    blocks = [i for i in screened["deviations_required"] + screened["missing"]
              if i["severity"] == BLOCK and i["analysis"] in (PRIMARY, "all")]
    if blocks:
        raise RunStop("the locked plan has a block on the primary analysis: " + "; ".join(i["topic"] for i in blocks))
    checks.append({"stage": "screen", "check": "the screen, run again without outcomes, gives the locked plan's "
                   "decisions and group sizes", "verdict": "pass", "detail": ""})
    log("screen: the locked plan's decisions")

    # 3. The checks that need outcomes, and the analysis cohort
    data = build_analysis(result.tables, settings, result.profile, result.source_steps)
    outcome_checks = data.checks + compare_with_plan(data, screened["decisions"])
    keep_plan_decisions(data, screened["decisions"])
    checks += [{"stage": "outcomes", **c.to_dict()} for c in outcome_checks]
    plan_cov = data.covariate_plan
    bands = sorted({v for v in plan_cov["year_bands"].values()}, key=list(plan_cov["year_bands"].values()).index)
    imputing = plan_cov["missing_data_rule"] == "multiple imputation"
    checks.append({"stage": "outcomes", "check": "the adjustment variables", "verdict": "applied" if imputing else "pass",
                   "detail": f"{plan_cov['missing_data_rule']} (the settings' rule)" +
                             (f", with {settings.read('missing_data.imputations')} imputations as the settings' "
                              "imputation model says" if imputing else "") +
                             (f"; year bands in the model: {', '.join(bands)}" if bands else "")})
    for c in outcome_checks:
        log(f"{c.check}: {c.verdict}")

    # 4. The run folder and the R part
    run_dir = work / f"run_plan_v{m['plan_version']}"
    run_dir.mkdir(exist_ok=True)
    spec = build_spec(data, settings)
    fresh = run_dir / "next"
    if fresh.exists():
        shutil.rmtree(fresh)
    fresh.mkdir()
    new_files = write_run_files(data, spec, fresh)
    unchanged = all(_same_bytes(f, run_dir / f.name) for f in new_files)
    r_info = {}
    if r.reuse and unchanged and all((run_dir / f"r_{s}.json").exists() for s in R_STEPS):
        shutil.rmtree(fresh)
        for s in R_STEPS:
            r_info[s] = json.loads((run_dir / f"r_{s}.json").read_text(encoding="utf-8"))
        checks.append({"stage": "R", "check": "the R part", "verdict": "pass",
                       "detail": "its inputs are byte for byte those of the results already in the run folder, which were kept"})
        log("R: inputs unchanged, results kept")
    else:
        for f in run_dir.glob("r_*"):
            f.unlink()
        for f in new_files:
            shutil.move(str(f), run_dir / f.name)
        shutil.rmtree(fresh)
        for s in R_STEPS:
            r_info[s] = run_r_step(r, s, run_dir, log)
            checks.append({"stage": "R", "check": f"R step: {s}", "verdict": "pass",
                           "detail": f"{r_info[s]['seconds']:.0f} seconds" + _warnings_text(r_info[s])})

    # 5. The gate and the reports, built in the run folder, read back, then published
    gate = Gate(data.min_count)
    res = Results(data, settings, run_dir, gate)
    tables, figure, withheld = build_tables(res)
    staging = run_dir / STAGING
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir()
    for name, t in tables.items():
        t.to_csv(staging / name, index=False, lineterminator="\n", float_format=FLOAT_FORMAT)
    (staging / "figure_km.svg").write_bytes(figure.encode("utf-8"))
    published = {name: pd.read_csv(staging / name, dtype=str, keep_default_na=False) for name in tables}
    info = _report_info(settings, result, m, plan, data, r_info, checks, stamp, gate)
    info["files"] = {k: v for k, v in FILES.items() if (staging / k).exists() or k.startswith("run_report")
                     or k == OUTPUTS_MANIFEST}
    doc = build_report(info, published, gate.log + [{"check": "security scan of every file to be published (IDs, "
                                                     "dates, row-level columns, small counts, file types)",
                                                     "verdict": "pass", "detail": "read back before copying"}],
                       res.notes + _limits(data, r_info))
    (staging / "run_report.md").write_bytes(doc.markdown().encode("utf-8"))
    (staging / "run_report.html").write_bytes(doc.html().encode("utf-8"))
    outputs = {f.name: sha256_file(f) for f in sorted(staging.iterdir())}
    (staging / OUTPUTS_MANIFEST).write_bytes((json.dumps(
        {"run_at_utc": stamp, "plan_version": m["plan_version"], "tool": f"delayed-entry-pipeline {__version__}",
         "inputs": [{"file": Path(i["path"]).name, "sha256": i["sha256"]} for i in m["inputs"]],
         "settings": {"file": Path(m["settings"]["path"]).name, "sha256": m["settings"]["sha256"]},
         "plan_files": m["plan_files"], "outputs": outputs}, indent=2) + "\n").encode("utf-8"))
    ids, numeric = plain_ids(set(result.tables.patients["patient_id"].astype(str))
                             | set(result.tables.tests["test_id"].astype(str)))
    allowed = _allowed_dates(settings, manifest, stamp)
    problems = check_outputs(staging, gate, ids, allowed)
    if problems:
        raise GateError("not safe to publish, so nothing was copied to the output folder: " + "; ".join(problems))
    out.mkdir(parents=True, exist_ok=True)
    for f in staging.iterdir():
        shutil.copy2(f, out / f.name)
    log(f"safety gate: passed; {len(list(staging.iterdir()))} files published to {out}")
    return {"output_folder": str(out), "files": sorted(f.name for f in staging.iterdir()),
            "checks": [{k: c[k] for k in ("stage", "check", "verdict")} for c in checks],
            "gate": [{k: c[k] for k in ("check", "verdict")} for c in gate.log], "withheld_curves": sorted(withheld),
            "ids_not_scannable": numeric}


def _warnings_text(info: dict) -> str:
    w = info.get("warnings") or {}
    n = sum(int(v) for part in w.values() for v in (part.values() if isinstance(part, dict) else [part]))
    return f"; {n} model warnings, listed in its log" if n else ""


def _allowed_dates(settings: Settings, manifest: Path, stamp: str) -> set[str]:
    texts = [settings.text, (manifest.parent / PLAN_MD).read_text(encoding="utf-8"), manifest.read_text(encoding="utf-8"),
             stamp]
    return {d for t in texts for d in DATE.findall(t)}


def _limits(data: AnalysisData, r_info: dict) -> list[str]:
    out = []
    if data.covariate_plan["missing_data_rule"] == "multiple imputation":
        m = r_info.get("adjusted", {}).get("model", {}).get("imputations")
        out.append(f"More eligible patients lacked an adjustment variable than the settings' complete-case share, so "
                   f"every adjusted and weighted model used multiple imputation ({m} imputations, as the settings' "
                   "imputation model says). Hazard ratios are pooled by Rubin's rules; each bootstrap resample was "
                   "imputed once; other point estimates, and pass 3b's balance, are averages across the imputations.")
    if r_info.get("weighted", {}).get("balanced") is False:
        out.append("Pass 3b's weighting left at least one term at or above the balance target, so its results are "
                   "labeled imperfectly balanced and pass 3a stays the main adjusted result.")
    if r_info.get("assumptions", {}).get("ph_rejected"):
        out.append("The proportional-hazards check rejected one hazard ratio per comparison: split hazard ratios are "
                   "added, and RMST differences are the main comparison summary.")
    return out


def _report_info(settings, result, m, plan, data, r_info, checks, stamp, gate) -> dict:
    versions = {}
    for s in r_info.values():
        if isinstance(s.get("packages"), dict):
            versions |= s["packages"]
    return {
        "study": settings.read("study.name"), "registered_plan": settings.read("study.registered_plan"),
        "settings_file": Path(m["settings"]["path"]).name, "settings_sha256": m["settings"]["sha256"],
        "plan_version": m["plan_version"], "locked_at": m["locked_at_utc"], "locked_by": m["locked_by"],
        "plan_files": m["plan_files"], "tool": f"delayed-entry-pipeline {__version__}", "run_at": stamp,
        "r_versions": versions, "resamples": settings.read("uncertainty.bootstrap_resamples"),
        "seed": settings.read("uncertainty.seed"),
        "inputs": [{"file": Path(i["path"]).name, "sha256": i["sha256"]} for i in m["inputs"]],
        "profile": result.profile.to_dict(), "screen": plan["readiness_screen"], "checks": checks,
        "analyses": [a.to_dict() for a in data.analyses], "min_count": gate.min_count,
        "horizon": settings.read("time.horizon_years"), "interval_level": settings.read("uncertainty.interval_level"),
    }


__all__ = ["run", "ROptions", "build_spec", "write_run_files", "r_command", "latest_manifest", "labels", "R_DIR"]
