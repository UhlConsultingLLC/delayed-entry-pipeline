"""Small counts and local folders. Copied from the CRC BPC project (crcbpc.privacy, crcbpc.cohort.public_flow
and crcbpc.paths.check_not_cloud), with the minimum count passed in from the settings instead of fixed.

A count from 1 to min_count - 1 is small: it is shown as "<min_count". Zero is shown.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

CLOUD_MARKERS = ("onedrive", "dropbox", "google drive", "googledrive", "icloud", "box sync")


def check_not_cloud(p: Path) -> Path:
    """Return p, or stop if it looks cloud-synced. Every folder that receives row-level data is checked."""
    low = str(p).lower().replace("\\", "/")
    for m in CLOUD_MARKERS:
        if m in low:
            raise ValueError(f"{p} looks cloud-synced ('{m}'). Row-level data must stay on local disk only.")
    return p


def is_small(n: int, min_count: int) -> bool:
    """A count the privacy rule hides: 1 to min_count - 1. Zero is shown."""
    return 0 < n < min_count


def show(n: int, min_count: int) -> str:
    """A count as a report prints it."""
    return f"<{min_count}" if is_small(int(n), min_count) else str(int(n))


def suppress_small_cells(table: pd.DataFrame, count_cols, dependent_cols=(), *, min_count: int,
                         label: str | None = None) -> pd.DataFrame:
    """Return a copy with small counts masked.

    count_cols: columns holding counts of people. A value from 1 to min_count - 1 is replaced by `label`.
    dependent_cols: statistics computed from those people (rates, medians, intervals). They are blanked on
        any row where a count was suppressed, so the hidden count can't be back-calculated."""
    label = label or f"<{min_count}"
    out = table.copy()
    small = pd.Series(False, index=out.index)
    for c in list(count_cols):
        v = pd.to_numeric(out[c], errors="coerce")
        hit = v.between(1, min_count - 1)
        small |= hit
        out[c] = out[c].astype(object)
        out.loc[hit, c] = label
    for c in dependent_cols:
        out[c] = out[c].astype(object)
        out.loc[small, c] = ""
    return out


def public_flow(flow_table: pd.DataFrame, *, min_count: int) -> pd.DataFrame:
    """The flow table with steps merged so that no small count is shown or can be worked out.

    Hiding a small count isn't enough in a flow diagram: if the counts remaining before and after a step
    are shown, the hidden count is their difference. So a step that excludes a small number is merged
    with the next step until the merged count is 0 or at least min_count. If the last steps are still
    too small, they merge into the step before. `flow_table` has columns label, excluded, remaining, and
    its first row is the start."""
    start, steps = flow_table.iloc[0], flow_table.iloc[1:]
    blocks: list[dict] = []
    for row in steps.itertuples(index=False):
        if blocks and is_small(blocks[-1]["excluded"], min_count):
            blocks[-1]["labels"].append(row.label)
            blocks[-1]["excluded"] += int(row.excluded)
        else:
            blocks.append({"labels": [row.label], "excluded": int(row.excluded)})
    while len(blocks) > 1 and is_small(blocks[-1]["excluded"], min_count):
        last = blocks.pop()
        blocks[-1]["labels"] += last["labels"]
        blocks[-1]["excluded"] += last["excluded"]

    rows, remaining = [{"label": start["label"], "excluded": 0, "remaining": int(start["remaining"])}], int(start["remaining"])
    for b in blocks:
        remaining -= b["excluded"]
        rows.append({"label": "; ".join(b["labels"]), "excluded": b["excluded"], "remaining": remaining})
    out = pd.DataFrame(rows)
    unsafe = [n for n in pd.concat([out["excluded"], out["remaining"]]) if is_small(int(n), min_count)]
    if unsafe:
        raise ValueError(f"can't publish this flow safely: {len(unsafe)} small counts would be shown or recoverable")
    return out
