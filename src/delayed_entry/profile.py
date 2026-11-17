"""The capability profile: what a data source can supply (design plan section 5).

Every adapter declares one. The GENIE BPC adapter ships its own (adapters/genie_bpc/profile.yml); for
flat files the person who confirms the mapping declares it in the mapping file. The readiness screen
reads it against what each analysis needs.
"""
from __future__ import annotations

from dataclasses import dataclass, field

ENTRY_EVENT_OPTIONS = ("report", "sample_collection", "order", "none")
COVERAGE_OPTIONS = ("regions", "genes_only", "none")
PANEL_SIZE_OPTIONS = ("given", "not_given")
DEATH_CAPTURE_OPTIONS = ("linked", "partial", "unknown")
TIME_FORMAT_OPTIONS = ("real_dates", "shifted_dates", "intervals")
ITEMS = {
    "entry_event": ENTRY_EVENT_OPTIONS,
    "panel_coverage": COVERAGE_OPTIONS,
    "panel_sizes": PANEL_SIZE_OPTIONS,
    "death_capture": DEATH_CAPTURE_OPTIONS,
    "time_format": TIME_FORMAT_OPTIONS,
}


class ProfileError(ValueError):
    pass


@dataclass(frozen=True)
class CapabilityProfile:
    adapter: str
    entry_event: str
    panel_coverage: str
    panel_sizes: str
    death_capture: str
    time_format: str
    privacy_min_count: int | None            # the data-use agreement's minimum count; None if it sets none
    histology_vocabulary: str                 # "ICD-O-3" or a description of the source's own categories
    age_top_code: float | None = None         # ages above this are withheld or top-coded
    privacy_other_limits: tuple[str, ...] = ()
    known_quirks: tuple[dict, ...] = ()        # each {"name": ..., "handling": ...}
    sources: dict = field(default_factory=dict)   # item: the source behind its value
    # What the source's terms require of work that uses or references its data: {citation, version, acknowledgment}.
    # None if the terms require nothing. The run report prints it; a run warns while a part of it is blank.
    attribution: dict | None = None

    def to_dict(self) -> dict:
        out = {"adapter": self.adapter, "entry_event": self.entry_event, "panel_coverage": self.panel_coverage,
               "panel_sizes": self.panel_sizes, "death_capture": self.death_capture, "time_format": self.time_format,
               "privacy_min_count": self.privacy_min_count, "histology_vocabulary": self.histology_vocabulary,
               "age_top_code": self.age_top_code, "privacy_other_limits": list(self.privacy_other_limits),
               "known_quirks": [dict(q) for q in self.known_quirks], "sources": dict(self.sources)}
        if self.attribution is not None:
            out["attribution"] = dict(self.attribution)
        return out

    def attribution_gaps(self) -> list[str]:
        """The parts of the attribution still blank (none if the source requires no attribution)."""
        if self.attribution is None:
            return []
        return [k for k in ATTRIBUTION_KEYS if not (isinstance(self.attribution.get(k), str) and self.attribution[k].strip())]


ATTRIBUTION_KEYS = ("citation", "version", "acknowledgment")


def _attribution_problems(a) -> list[str]:
    if a is None:
        return []
    if not isinstance(a, dict) or set(a) != set(ATTRIBUTION_KEYS):
        return ["attribution must be null, or {citation: ..., version: ..., acknowledgment: ...} (a blank part is null)"]
    if not (isinstance(a["citation"], str) and a["citation"].strip()):
        return ["attribution.citation must be given"]
    return [f"attribution.{k} must be text or null" for k in ("version", "acknowledgment")
            if a[k] is not None and not isinstance(a[k], str)]


def profile_from_dict(d: dict, adapter: str | None = None) -> CapabilityProfile:
    """Check a profile written as a mapping. Every item must be given; none is guessed."""
    if not isinstance(d, dict):
        raise ProfileError("the capability profile must be a mapping")
    problems = []
    for item, options in ITEMS.items():
        if d.get(item) not in options:
            problems.append(f"{item} must be one of: {', '.join(options)}")
    if "privacy_min_count" not in d:
        problems.append("privacy_min_count must be given (a whole number, or null if the agreement sets none)")
    elif d["privacy_min_count"] is not None and not (isinstance(d["privacy_min_count"], int) and d["privacy_min_count"] >= 1):
        problems.append("privacy_min_count must be a whole number of at least 1, or null")
    if not (isinstance(d.get("histology_vocabulary"), str) and d["histology_vocabulary"].strip()):
        problems.append("histology_vocabulary must say what the histology field holds")
    top = d.get("age_top_code")
    if top is not None and (isinstance(top, bool) or not isinstance(top, (int, float))):
        problems.append("age_top_code must be a number or null")
    quirks = d.get("known_quirks", [])
    if not isinstance(quirks, list) or any(not isinstance(q, dict) or set(q) != {"name", "handling"} for q in quirks):
        problems.append("known_quirks must be a list of {name, handling}")
    sources = d.get("sources", {})
    if not isinstance(sources, dict):
        problems.append("sources must map each item to its source")
    problems += _attribution_problems(d.get("attribution"))
    unknown = set(d) - {*ITEMS, "adapter", "privacy_min_count", "histology_vocabulary", "age_top_code",
                        "privacy_other_limits", "known_quirks", "sources", "attribution"}
    if unknown:
        problems.append(f"unknown items: {', '.join(sorted(unknown))}")
    if problems:
        raise ProfileError("capability profile: " + "; ".join(problems))
    return CapabilityProfile(
        adapter=adapter or d.get("adapter") or "unnamed", entry_event=d["entry_event"],
        panel_coverage=d["panel_coverage"], panel_sizes=d["panel_sizes"], death_capture=d["death_capture"],
        time_format=d["time_format"], privacy_min_count=d["privacy_min_count"],
        histology_vocabulary=d["histology_vocabulary"], age_top_code=top,
        privacy_other_limits=tuple(d.get("privacy_other_limits", []) or []),
        known_quirks=tuple(dict(q) for q in quirks), sources=dict(sources),
        attribution=dict(d["attribution"]) if d.get("attribution") is not None else None)
