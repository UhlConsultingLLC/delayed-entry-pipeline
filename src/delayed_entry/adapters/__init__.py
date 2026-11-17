"""Adapters: each reads one kind of source into the standard tables, with its capability profile.

    genie_bpc   GENIE BPC releases (built from the CRC BPC project's cohort code)
    flatfile    CSV/TSV, Excel and Parquet files, with a mapping a person has confirmed
    omop        an OMOP CDM database (a folder of table files or a DuckDB file), with a mapping a person has
                confirmed
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from delayed_entry.profile import CapabilityProfile
from delayed_entry.tables import StandardTables


@dataclass
class AdapterResult:
    tables: StandardTables
    profile: CapabilityProfile
    input_files: list[Path]                                   # every file read, for the lock's manifest
    source_steps: list[tuple[str, int]] = field(default_factory=list)   # exclusions made while reading
    read_notes: list[str] = field(default_factory=list)       # what reading found (counts are masked)
