"""Locking the plan (design plan section 7).

    from delayed_entry.lock import lock_plan, require_lock
    report = screen(outcome_blind(tables), settings, profile, inputs=files)
    result = lock_plan(report, settings, files, out_dir, locked_by="name of the reviewer")
    require_lock(result.manifest, settings, files)     # what the run calls first: stops if anything changed

lock_plan writes the filled-in plan twice, in Markdown for people and in YAML for the run. It combines the
settings (every choice, who made it and its source), the adapter's capability profile and the readiness
screen's three lists. A person reviews it first: the lock records who. The manifest holds a SHA-256 hash
of the plan files, the settings file, the defaults file and every input file, with a timestamp.

The lock refuses a plan it can't stand behind: no input files, a report screened from other inputs or
other settings, or a privacy minimum below the data-use agreement's. Nothing is written until every check
has passed, and the version folder appears whole or not at all.

A locked plan changes only by a new version: locking again into the same folder needs a reason, and adds a
line (who, when, why) to the folder's deviations log. Registering the plan elsewhere (for example on OSF)
is a manual step this module prepares for but doesn't do.

lock_settings locks the settings file alone, before any data are downloaded.
"""
from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import yaml

from delayed_entry import __version__
from delayed_entry.hashing import sha256_file
from delayed_entry.screen import BLOCK, ReadinessReport
from delayed_entry.settings import DEFAULTS_FILE, Settings

MANIFEST = "manifest.json"
PLAN_MD, PLAN_YAML, SETTINGS_COPY = "plan.md", "plan.yaml", "settings.yaml"
DEVIATIONS_LOG = "deviations_log.md"
SETTINGS_LOCK = "settings_lock.json"

__all__ = ["lock_plan", "verify_lock", "require_lock", "lock_settings", "settings_lock_matches", "sha256_file",
           "LockError", "LockResult"]


class LockError(RuntimeError):
    pass


def _now(now: datetime | None) -> str:
    return (now or datetime.now(timezone.utc)).astimezone(timezone.utc).isoformat(timespec="seconds")


def _write(path: Path, text: str) -> None:
    path.write_bytes(text.encode("utf-8"))


@dataclass(frozen=True)
class LockResult:
    folder: Path
    version: int
    manifest: Path


def _versions(out_dir: Path) -> list[int]:
    return sorted(int(p.name.split("_v")[1]) for p in out_dir.glob("plan_v*") if p.is_dir() and p.name.split("_v")[1].isdigit())


def plan_yaml(report: ReadinessReport, settings: Settings, version: int, locked_by: str, reason: str | None,
              read_notes: list[str] | None) -> str:
    doc = {
        "plan_version": version,
        "study": report.study,
        "reviewed_and_locked_by": locked_by,
        "reason_for_this_version": reason,
        "settings": settings.provenance(),
        "defaults_used": settings.defaults_used(),
        "capability_profile": report.profile,
        "reading_the_source": list(read_notes or []),
        "readiness_screen": report.to_dict(),
    }
    return yaml.safe_dump(doc, sort_keys=False, allow_unicode=True, width=120)


def plan_markdown(report: ReadinessReport, settings: Settings, version: int, locked_by: str, reason: str | None,
                  read_notes: list[str] | None) -> str:
    out = [f"# Analysis plan, version {version}: {report.study}", "",
           f"Reviewed and locked by: {locked_by}.",
           *([f"Reason for this version: {reason}."] if reason else []),
           "", "Every choice below is locked with this plan. It changes only by a new plan version, which records "
           "who changed it, when and why, and adds a line to the deviations log.", "",
           "## Settings", "", "| Choice | Value | Made by | Source | Where from |", "|---|---|---|---|---|"]
    for s in settings.provenance():
        v = yaml.safe_dump(s["value"], default_flow_style=True, width=10**6).replace("\n...", "").strip()
        out.append(f"| {s['choice']} | `{v}` | {s['made_by']} | {s['source']} | {s['origin']} |")
    defaults = settings.defaults_used()
    out += ["", f"Defaults used: {len(defaults)}."] + [f"- {d['choice']}: {d['source']}" for d in defaults]
    out += ["", "## Capability profile of the source", ""]
    for k, v in report.profile.items():
        out.append(f"- {k}: {v}")
    if read_notes:
        out += ["", "## Reading the source", ""] + [f"- {n}" for n in read_notes]
    out += ["", "## Input files", ""] + [f"- {i['file']}: {i['sha256']}" for i in report.inputs]
    out += ["", report.to_markdown().replace("# Readiness report", "## Readiness report", 1).replace("\n## ", "\n### ")]
    return "\n".join(out) + "\n"


def _input_entries(inputs: list[str | Path]) -> list[dict]:
    entries = []
    for p in sorted({Path(x).resolve() for x in inputs}):
        if not p.is_file():
            raise LockError(f"input file not found: {p.name}")
        entries.append({"path": str(p), "bytes": p.stat().st_size, "sha256": sha256_file(p)})
    return entries


def lock_plan(report: ReadinessReport, settings: Settings, inputs: list[str | Path], out_dir: str | Path,
              locked_by: str, reason: str | None = None, read_notes: list[str] | None = None,
              now: datetime | None = None) -> LockResult:
    """Write plan_v<n>/ (plan.md, plan.yaml, a copy of the settings, manifest.json) into out_dir."""
    if not (isinstance(locked_by, str) and locked_by.strip()):
        raise LockError("a person reviews the plan before it is locked: give locked_by")
    if report.settings_sha256 != settings.sha256 or \
            (report.settings_resolved_sha256 and report.settings_resolved_sha256 != settings.resolved_sha256):
        raise LockError("the report was made with other settings: screen again with these settings")
    if any(i.topic == "privacy minimum" and i.severity == BLOCK for i in report.items):
        raise LockError("the settings' privacy minimum is below the data-use agreement's: fix the settings first")
    if not inputs:
        raise LockError("give the input files: a plan locked without them can't show the data didn't change")
    entries = _input_entries(inputs)
    if not report.inputs:
        raise LockError("the report doesn't record the files it was screened from: screen again with inputs=")
    if sorted(e["sha256"] for e in entries) != sorted(i["sha256"] for i in report.inputs):
        raise LockError("the report was screened from other input files than the ones given")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    versions = _versions(out_dir)
    version = (versions[-1] + 1) if versions else 1
    if versions and not (isinstance(reason, str) and reason.strip()):
        raise LockError(f"a plan is already locked here (version {versions[-1]}): a new version needs a reason")
    stamp = _now(now)
    folder = out_dir / f"plan_v{version}"
    partial = out_dir / f".plan_v{version}.partial"
    if partial.exists():
        shutil.rmtree(partial)
    partial.mkdir()
    try:
        _write(partial / PLAN_MD, plan_markdown(report, settings, version, locked_by, reason, read_notes))
        _write(partial / PLAN_YAML, plan_yaml(report, settings, version, locked_by, reason, read_notes))
        _write(partial / SETTINGS_COPY, settings.text)
        manifest = {
            "plan_version": version,
            "locked_at_utc": stamp,
            "locked_by": locked_by,
            "reason": reason,
            "tool": f"delayed-entry-pipeline {__version__}",
            "plan_files": {name: sha256_file(partial / name) for name in (PLAN_MD, PLAN_YAML, SETTINGS_COPY)},
            "settings": {"path": str(settings.path.resolve()) if settings.path else None, "sha256": settings.sha256,
                         "resolved_sha256": settings.resolved_sha256},
            "defaults": {"path": str(DEFAULTS_FILE.resolve()), "sha256": sha256_file(DEFAULTS_FILE)},
            "inputs": entries,
        }
        _write(partial / MANIFEST, json.dumps(manifest, indent=2) + "\n")
        partial.rename(folder)
    except BaseException:
        shutil.rmtree(partial, ignore_errors=True)
        raise
    if version > 1:
        log = out_dir / DEVIATIONS_LOG
        head = "" if log.exists() else "# Deviations log\n\n| Version | When (UTC) | Who | Why |\n|---|---|---|---|\n"
        with open(log, "ab") as f:
            f.write((head + f"| {version} | {stamp} | {locked_by} | {reason} |\n").encode("utf-8"))
    return LockResult(folder, version, folder / MANIFEST)


def verify_lock(manifest_path: str | Path, settings: Settings | None = None,
                inputs: list[str | Path] | None = None) -> list[str]:
    """What has changed since the lock: [] if nothing. The plan files, the settings file, the defaults file and
    every input are hashed again where the manifest recorded them. With `settings` or `inputs`, also checks
    that they are the ones locked."""
    manifest_path = Path(manifest_path)
    m = json.loads(manifest_path.read_text(encoding="utf-8"))
    folder = manifest_path.parent
    problems = []
    for name, digest in m["plan_files"].items():
        p = folder / name
        if not p.exists():
            problems.append(f"{name} is missing")
        elif sha256_file(p) != digest:
            problems.append(f"{name} has changed since the lock")
    for label, rec in (("the settings file", m["settings"]), ("the defaults file", m.get("defaults"))):
        if rec and rec.get("path"):
            p = Path(rec["path"])
            if not p.exists():
                problems.append(f"{label} is missing")
            elif sha256_file(p) != rec["sha256"]:
                problems.append(f"{label} has changed since the lock")
    for e in m["inputs"]:
        p = Path(e["path"])
        if not p.exists():
            problems.append(f"input {p.name} is missing")
        elif p.stat().st_size != e["bytes"] or sha256_file(p) != e["sha256"]:
            problems.append(f"input {p.name} has changed since the lock")
    if settings is not None and (settings.sha256 != m["settings"]["sha256"]
                                 or settings.resolved_sha256 != m["settings"].get("resolved_sha256")):
        problems.append("the settings are not the ones locked")
    if inputs is not None:
        now = {str(Path(x).resolve()) for x in inputs}
        locked = {e["path"] for e in m["inputs"]}
        if now != locked:
            problems.append(f"the inputs are not the ones locked ({len(now - locked)} added, {len(locked - now)} left out)")
    return problems


def require_lock(manifest_path: str | Path, settings: Settings, inputs: list[str | Path]) -> None:
    """Stop unless the lock still matches these settings and these inputs: what the run calls first."""
    problems = verify_lock(manifest_path, settings, inputs)
    if problems:
        raise LockError("the locked plan doesn't match: " + "; ".join(problems))


def lock_settings(settings_path: str | Path, out_dir: str | Path, locked_by: str, now: datetime | None = None) -> Path:
    """Lock a settings file alone, before any data are downloaded: its hash, who and when."""
    if not (isinstance(locked_by, str) and locked_by.strip()):
        raise LockError("give locked_by: the person who reviewed the settings")
    settings_path = Path(settings_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / SETTINGS_LOCK
    if target.exists():
        raise LockError(f"{SETTINGS_LOCK} already exists: a locked choice changes only by a new plan version")
    record = {"settings": settings_path.name, "sha256": sha256_file(settings_path), "locked_at_utc": _now(now),
              "locked_by": locked_by}
    _write(target, json.dumps(record, indent=2) + "\n")
    return target


def settings_lock_matches(lock_file: str | Path, settings_path: str | Path) -> bool:
    record = json.loads(Path(lock_file).read_text(encoding="utf-8"))
    return record["sha256"] == sha256_file(Path(settings_path))
