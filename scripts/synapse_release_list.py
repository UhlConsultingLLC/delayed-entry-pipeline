"""List a GENIE BPC release's files with their sizes and checksums, from Synapse's metadata: no file is downloaded.

    python scripts/synapse_release_list.py <release>

    release   a release in src/delayed_entry/adapters/genie_bpc/synapse_releases.yml, for example nsclc_v2.0

It writes src/delayed_entry/adapters/genie_bpc/cohorts/<release>/release_files.yml: for every file in the
release's clinical and genomic folders, its folder, name, Synapse ID, version, size in bytes and MD5
checksum, and whether the study reads it (`used`). The check after downloading compares the downloaded files
with this list. File names and sizes are not patient data.

Reading sizes and checksums needs an approved Synapse account. Your Synapse personal access token is read by
the Synapse client itself, from the SYNAPSE_AUTH_TOKEN environment variable or ~/.synapseConfig; it is never
passed on the command line, printed or stored here. Needs the Synapse client (this package's [synapse] extra).
"""
from __future__ import annotations

import datetime
import fnmatch
import os
import sys
from pathlib import Path

import yaml

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE / "src"))

from delayed_entry.adapters.genie_bpc import cohort_spec, raw  # noqa: E402

RELEASES = PACKAGE / "src" / "delayed_entry" / "adapters" / "genie_bpc" / "synapse_releases.yml"
PARTS = ("clinical", "genomic")


def is_used(name: str, cohort: str) -> bool:
    """Whether the study reads this file: a raw file of the adapter's spec, a gene-panel file, or a file the
    cohort lists under also_used."""
    base = name.rsplit("/", 1)[-1]
    spec_files = {f for f, _, _ in raw.raw_files(cohort).values()}
    return base in spec_files or fnmatch.fnmatch(base, raw.PANEL_FILES) or base in (cohort_spec(cohort).get("also_used") or [])


def release_file_entry(part: str, name: str, syn_id: str, version: int, size: int, md5: str, cohort: str) -> dict:
    return {"folder": part, "name": name, "id": syn_id, "version": int(version), "bytes": int(size), "md5": md5,
            "used": is_used(name, cohort)}


def write_release_files(entries: list[dict], cohort: str, path: Path | None = None) -> Path:
    path = path or raw.cohort_dir(cohort) / "release_files.yml"
    head = (f"# The files of {cohort} with their sizes and MD5 checksums, written by scripts/synapse_release_list.py\n"
            f"# from Synapse's metadata on {datetime.date.today().isoformat()}, without downloading any file.\n"
            "# used: true marks the files the study reads; a difference in one of them stops the run.\n")
    body = yaml.safe_dump({"files": entries}, sort_keys=False, default_flow_style=None, width=200)
    path.write_bytes((head + body).encode("utf-8"))
    return path


def size_and_md5(entity) -> tuple[int | None, str | None]:
    """The file's size and MD5 from its metadata, wherever this Synapse client version keeps them."""
    handle = getattr(entity, "_file_handle", None) or {}
    get = handle.get if hasattr(handle, "get") else (lambda k: None)
    return (get("contentSize") or getattr(entity, "contentSize", None),
            get("contentMd5") or getattr(entity, "md5", None))


def list_release(release: str) -> list[dict]:
    import synapseclient

    folders = yaml.safe_load(RELEASES.read_text(encoding="utf-8"))[release]["parts"]
    syn = synapseclient.Synapse()
    syn.login()                                   # uses SYNAPSE_AUTH_TOKEN or ~/.synapseConfig

    def walk(parent: str, prefix: str = ""):
        for child in syn.getChildren(parent, includeTypes=["file", "folder"]):
            if child["type"].endswith("Folder"):
                yield from walk(child["id"], prefix + child["name"] + "/")
            else:
                yield prefix + child["name"], child["id"]

    entries = []
    for part in PARTS:
        for name, syn_id in walk(folders[part]["folder"]):
            e = syn.get(syn_id, downloadFile=False)
            size, md5 = size_and_md5(e)
            if size is None or not md5:
                raise SystemExit(f"Stopped: Synapse gave no size or checksum for {name}: is access to this release approved?")
            entries.append(release_file_entry(part, name, syn_id, e.versionNumber, size, md5, release))
    return entries


def main(argv: list[str] | None = None) -> bool:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        sys.exit(__doc__)
    release = argv[0]
    entries = list_release(release)
    path = write_release_files(entries, release)
    used = sum(e["used"] for e in entries)
    print(f"{len(entries)} files listed ({used} read by the study) -> {path}")
    return bool(entries)


if __name__ == "__main__":
    ok = main()
    # As in the download script: on Windows the Synapse client can crash while Python shuts down, after the work
    # is done, so leave without that shutdown.
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0 if ok else 1)
