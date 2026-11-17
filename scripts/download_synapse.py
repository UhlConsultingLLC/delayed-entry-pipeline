"""Download one part of a GENIE BPC release from Synapse, and write a checksum manifest.

    python scripts/download_synapse.py <release> <part> <data folder> [--settings FILE --settings-lock FILE]

    release      a release in src/delayed_entry/adapters/genie_bpc/synapse_releases.yml, for example nsclc_v2.0
    part         documentation, clinical or genomic
    data folder  the release's local data folder; the files go into <data folder>/<part>

The documentation can be downloaded at any time. The clinical and genomic parts hold patient-level data, so
this script downloads them only when given the study's settings and the lock made before the download
(python -m delayed_entry lock-settings): the settings must load, must name this release in
source.genie_bpc.cohort, and must be unchanged since the lock.

The data folder must be on local disk, not in a cloud-synced folder, and outside this package.

Your Synapse personal access token is read by the Synapse client itself, from the SYNAPSE_AUTH_TOKEN
environment variable or from ~/.synapseConfig. It is never passed on the command line, printed or stored
here. Needs the Synapse client (pip install "synapseclient>=4.0", or this package's [synapse] extra).

Based on the CRC BPC project's scripts/download_synapse_folder.py.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import os
import sys
from pathlib import Path

import yaml

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE / "src"))

from delayed_entry.lock import settings_lock_matches  # noqa: E402
from delayed_entry.privacy import check_not_cloud  # noqa: E402
from delayed_entry.settings import SettingsError, load_settings  # noqa: E402

RELEASES = PACKAGE / "src" / "delayed_entry" / "adapters" / "genie_bpc" / "synapse_releases.yml"


class Refused(SystemExit):
    """The download may not go ahead."""


def releases(path: Path = RELEASES) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def destination(data_folder: str | Path, part: str) -> Path:
    """Where the files go: <data folder>/<part>, on local disk and outside this package."""
    folder = Path(data_folder).resolve()
    try:
        check_not_cloud(folder)
    except ValueError as e:
        raise Refused(f"Stopped: {e}") from None
    if folder == PACKAGE or PACKAGE in folder.parents:
        raise Refused("Stopped: the data folder must be outside this package, which holds code only.")
    return folder / part


def check_allowed(release: str, part: str, settings: str | Path | None, lock: str | Path | None,
                  known: dict | None = None) -> str:
    """The Synapse folder to download, or Refused saying why not."""
    known = releases() if known is None else known
    if release not in known:
        raise Refused(f"Stopped: unknown release '{release}' (known: {', '.join(known)}).")
    parts = known[release]["parts"]
    if part not in parts:
        raise Refused(f"Stopped: unknown part '{part}' (parts: {', '.join(parts)}).")
    spec = parts[part]
    if not spec["patient_data"]:
        return spec["folder"]
    if settings is None or lock is None:
        raise Refused(f"Stopped: the {part} files hold patient-level data. Give --settings and --settings-lock: "
                      "the study's settings are locked before its data are downloaded.")
    try:
        s = load_settings(settings)
    except SettingsError:
        raise Refused("Stopped: the settings file doesn't load (run: python -m delayed_entry check-settings).") from None
    if not s.has("source.genie_bpc.cohort") or s.read("source.genie_bpc.cohort") != release:
        raise Refused(f"Stopped: these settings are not for {release} (their source.genie_bpc.cohort differs).")
    if not Path(lock).is_file() or not settings_lock_matches(lock, settings):
        raise Refused("Stopped: the settings have no matching lock (run: python -m delayed_entry lock-settings), "
                      "or have changed since they were locked.")
    return spec["folder"]


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(syn_id: str, dest: Path, label: str) -> bool:
    """Download the folder and write its manifest. True when files arrived and every checksum matches Synapse."""
    import synapseclient
    import synapseutils

    dest.mkdir(parents=True, exist_ok=True)
    syn = synapseclient.Synapse()
    syn.login()                                   # uses SYNAPSE_AUTH_TOKEN or ~/.synapseConfig
    files = synapseutils.syncFromSynapse(syn, syn_id, path=str(dest))
    rows = []
    for f in files:
        p = Path(f.path)
        rows.append({"synapse_id": f.id, "version": f.versionNumber, "name": p.name, "bytes": p.stat().st_size,
                     "md5_local": md5(p), "md5_synapse": getattr(f, "md5", "")})
    stamp = datetime.date.today().isoformat()
    manifest = dest / f"manifest_{label}_{stamp}.csv"
    with open(manifest, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["synapse_id"])
        w.writeheader()
        w.writerows(rows)
    mism = [r["name"] for r in rows if r["md5_synapse"] and r["md5_local"] != r["md5_synapse"]]
    print(f"{len(rows)} files -> {dest}\nmanifest: {manifest}\nchecksum mismatches: {mism or 'none'}")
    if not rows or mism:
        print("Download failed: " + ("no files arrived." if not rows else "download those files again."))
        return False
    print("Now record the download: date, release, Synapse ID, manifest name, and the date you read the terms.")
    return True


def main(argv: list[str] | None = None) -> bool:
    ap = argparse.ArgumentParser(description="Download one part of a GENIE BPC release from Synapse.")
    ap.add_argument("release")
    ap.add_argument("part")
    ap.add_argument("data_folder")
    ap.add_argument("--settings", default=None)
    ap.add_argument("--settings-lock", default=None)
    a = ap.parse_args(argv)
    syn_id = check_allowed(a.release, a.part, a.settings, a.settings_lock)
    dest = destination(a.data_folder, a.part)
    print(f"{a.release} {a.part}: Synapse {syn_id} -> {dest}")
    return download(syn_id, dest, f"{a.release}_{a.part}")


if __name__ == "__main__":
    ok = main()
    # The download is finished and checked at this point. On Windows the Synapse client can then crash while Python
    # shuts down (an access violation), after a complete and correct download, so leave without that shutdown: the
    # exit code reports the download itself.
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0 if ok else 1)
