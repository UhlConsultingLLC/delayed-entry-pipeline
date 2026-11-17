"""Panel sizes from a GENIE genomic_information file.

GENIE's genomic_information.txt lists, for each panel (SEQ_ASSAY_ID), the regions its center targets:
Chromosome, Start_Position and End_Position (from the center's BED file, on GRCh37, so a region is End_Position -
Start_Position bases long), Hugo_Symbol, Feature_Type ("exon", "intron" or "intergenic") and includeInPanel (the
regions GENIE uses to define the panel's gene list). GENIE labels a submitted region "exon" when it overlaps an
exon, and keeps the region whole; a region that overlaps several exons can be listed once for each. So a panel's
regions are merged where they overlap before they are measured. Sources: the GENIE data guide (its template,
templates/data_guide_template.Rnw) and its processing code (genie_registry/bed.py), in Sage-Bionetworks/Genie.

GENIE BPC releases have no such file; the main GENIE releases do, for the same panels. Which regions count is the
settings' rule (source.genie_bpc.panel_sizes), never the code's.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PATTERN = "genomic_information*"
COLUMNS = ("SEQ_ASSAY_ID", "Chromosome", "Start_Position", "End_Position", "Feature_Type", "includeInPanel")
MARKED = {"true", "1", "yes", "t"}
BASES_PER_MEGABASE = 1_000_000


def find(genomic_dir: Path) -> list[Path]:
    return sorted(Path(genomic_dir).glob(PATTERN))


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, encoding="utf-8")


def problems(df: pd.DataFrame) -> list[str]:
    missing = [c for c in COLUMNS if c not in df.columns]
    return [f"no column {', '.join(missing)}"] if missing else []


def merged_length(starts, ends) -> int:
    """The bases a set of half-open regions covers, each base counted once."""
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(zip(starts, ends)):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    return total + (cur_e - cur_s if cur_e is not None else 0)


def panel_sizes(df: pd.DataFrame, rule: dict) -> tuple[pd.DataFrame, int]:
    """Each panel's size in megabases by the rule ({feature_types, include_in_panel}), and how many rows of those
    feature types couldn't be measured (a position that isn't a whole number, or an end not after its start)."""
    rows = df[df["Feature_Type"].str.strip().str.lower().isin([f.lower() for f in rule["feature_types"]])]
    if rule["include_in_panel"] == "only_marked":
        rows = rows[rows["includeInPanel"].str.strip().str.lower().isin(MARKED)]
    start = pd.to_numeric(rows["Start_Position"], errors="coerce")
    end = pd.to_numeric(rows["End_Position"], errors="coerce")
    ok = start.notna() & end.notna() & (end > start)
    rows = rows.assign(_s=start, _e=end)[ok]
    sizes = [(panel, sum(merged_length(g["_s"].astype(int), g["_e"].astype(int)) for _, g in per.groupby("Chromosome"))
              / BASES_PER_MEGABASE) for panel, per in rows.groupby("SEQ_ASSAY_ID")]
    out = pd.DataFrame(sizes, columns=["panel_id", "coding_megabases"])
    return out.astype({"panel_id": "string", "coding_megabases": "Float64"}), int((~ok).sum())
