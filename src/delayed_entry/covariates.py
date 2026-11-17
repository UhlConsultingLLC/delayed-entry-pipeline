"""The adjustment variables, coded as the settings say, for the adjusted and weighted models (passes 3a and 3b).

Generalized from the CRC BPC project's sql/04_covariates.sql and crcbpc.covariates. Every level, band, knot
and threshold comes from the settings; this module only applies them.

    from delayed_entry.covariates import code_adjustment, missing_data_rule, merge_bands
    coded = code_adjustment(patients, settings)      # one column per variable, levels as the settings say
    missing_data_rule(coded, settings)               # "complete case" or "multiple imputation"

A variable's form decides its coding:
    linear, spline   the number itself (R adds the spline's curved terms, from knots set on the model's patients)
    categorical      text, trimmed; `merge` gathers the source's values into new levels, `keep` lists the levels
                     kept, and any other value is missing; `reference` is the first level (else the first listed)
    bands            the year of diagnosis, in the settings' year bands; a year in no band is missing

Small year bands are merged before modeling (merge_bands), from counts alone: no outcome is read.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from delayed_entry.settings import Settings

NUMERIC_FORMS = ("linear", "spline")


@dataclass(frozen=True)
class Term:
    """One adjustment variable as the models use it."""
    field: str                       # the patients-table column, and the model's column name
    form: str                        # linear, spline, categorical or bands
    levels: tuple[str, ...] = ()     # for categorical and bands: the reference first
    knots_at_percentiles: tuple[float, ...] = ()

    def to_dict(self) -> dict:
        return {"field": self.field, "form": self.form, "levels": list(self.levels),
                "knots_at_percentiles": list(self.knots_at_percentiles)}


def humanize(name: str) -> str:
    """A field's name as a reader would write it: 'tumor_location' -> 'Tumor location'."""
    text = name.replace("_", " ").strip()
    return text[:1].upper() + text[1:]


# ---------------------------------------------------------------------------
# Year bands
# ---------------------------------------------------------------------------
def band_of(year: pd.Series, bands: list[dict]) -> pd.Series:
    """Each year's band label: the first band whose `from` and `to` (both included) hold it. Blank if none does."""
    y = pd.to_numeric(year, errors="coerce")
    out = pd.Series(pd.NA, index=year.index, dtype="string")
    for b in reversed(bands):
        hit = y.notna()
        if "from" in b:
            hit &= y >= b["from"]
        if "to" in b:
            hit &= y <= b["to"]
        out = out.mask(hit.fillna(False).astype(bool), b["label"])
    return out


def merged_label(block: list[dict]) -> str:
    """The name of merged bands: '2015 or earlier' when the first is open below, '2016 or later' when the last is
    open above, otherwise '2013-2018'. One band keeps its own label."""
    if len(block) == 1:
        return block[0]["label"]
    first, last = block[0], block[-1]
    if "from" not in first and "to" in last:
        return f"{last['to']} or earlier"
    if "to" not in last and "from" in first:
        return f"{first['from']} or later"
    return f"{first.get('from', '')}-{last.get('to', '')}"


def merge_bands(counts: dict[str, int], bands: list[dict], min_n: int) -> dict[str, str]:
    """Which model band each year band goes into (the CRC BPC plan's rule, section 4c, for any number of bands).

    First, an end band with fewer than `min_n` patients merges into its neighbor (the earliest end first). Then,
    while a band (with whatever joined it) still has fewer than `min_n`, it merges with the smaller of its
    neighbors, the earlier one on a tie. Stops if only one band would be left: the year could then not be
    modeled, and a person decides. Counts only: no outcome is read."""
    blocks = [[b] for b in bands]

    def size(i):
        return sum(counts.get(b["label"], 0) for b in blocks[i])

    def join(i, j):                  # merge block j into block i, keeping calendar order
        lo, hi = sorted((i, j))
        blocks[lo] = blocks[lo] + blocks[hi]
        del blocks[hi]

    if len(blocks) > 1 and size(0) < min_n:
        join(0, 1)
    if len(blocks) > 1 and size(len(blocks) - 1) < min_n:
        join(len(blocks) - 2, len(blocks) - 1)
    while len(blocks) > 1:
        small = [i for i in range(len(blocks)) if size(i) < min_n]
        if not small:
            break
        i = small[0]
        neighbors = [j for j in (i - 1, i + 1) if 0 <= j < len(blocks)]
        target = min(neighbors, key=size)           # min() keeps the first of a tie: the earlier
        join(i, target)
    if len(blocks) < 2:
        raise ValueError("every year band is small: the year can't be modeled as planned, and a person decides")
    return {b["label"]: merged_label(block) for block in blocks for b in block}


# ---------------------------------------------------------------------------
# Coding
# ---------------------------------------------------------------------------
@dataclass
class Coded:
    frame: pd.DataFrame                     # indexed by patient_id: one column per variable
    terms: list[Term]
    notes: list[str] = field(default_factory=list)


def _categorical(s: pd.Series, v: dict) -> tuple[pd.Series, list[str]]:
    text = s.astype("string").str.strip()
    text = text.mask(text == "")
    merge = v.get("merge")
    if merge:
        mapping = {old: new for new, olds in merge.items() for old in olds}
        coded = text.map(lambda x: mapping.get(x, pd.NA) if not pd.isna(x) else pd.NA).astype("string")
        natural = list(merge)
    else:
        coded = text
        natural = list(v["keep"]) if v.get("keep") else sorted(text.dropna().unique())
    if v.get("keep"):
        coded = coded.where(coded.isin(v["keep"]))
        natural = list(v["keep"])
    ref = v.get("reference") or (natural[0] if natural else None)
    levels = ([ref] if ref is not None else []) + [x for x in natural if x != ref]
    return coded, levels


def code_adjustment(patients: pd.DataFrame, settings: Settings) -> Coded:
    """Every adjustment variable for the patients given (indexed by patient_id), coded as the settings say.
    A year band here is the settings' band; merge_bands decides the model's bands later, from counts."""
    variables = settings.read("adjustment.variables")
    bands = settings.read("adjustment.year_bands")
    out, terms = pd.DataFrame(index=patients.index), []
    for v in variables:
        name, form = v["field"], v["form"]
        if name not in patients.columns:
            raise KeyError(f"patients.{name} is absent: the adjusted models need it")
        s = patients[name]
        if form in NUMERIC_FORMS:
            out[name] = pd.to_numeric(s, errors="coerce").astype("Float64")
            terms.append(Term(name, form, knots_at_percentiles=tuple(v.get("knots_at_percentiles", ()))))
        elif form == "categorical":
            out[name], levels = _categorical(s, v)
            terms.append(Term(name, form, tuple(levels)))
        else:
            out[name] = band_of(s, bands)
            terms.append(Term(name, form, tuple(b["label"] for b in bands)))
    return Coded(out, terms)


def missing_any(coded: Coded) -> pd.Series:
    """True for each patient missing any adjustment variable."""
    return coded.frame.isna().any(axis=1)


def missing_data_rule(coded: Coded, settings: Settings, patients: pd.Index | None = None) -> tuple[str, float]:
    """The CRC BPC plan's rule (section 9), with the settings' share: "complete case" if at most that share of the
    patients (the eligible ones) lacks any adjustment variable, otherwise "multiple imputation". Also the share."""
    limit = settings.read("missing_data.complete_case_max_share")
    m = missing_any(coded) if patients is None else missing_any(coded).reindex(patients)
    share = float(m.mean()) if len(m) else 0.0
    return ("complete case" if share <= limit else "multiple imputation"), share


def model_bands(coded: Coded, settings: Settings, in_model: pd.Index) -> dict[str, str]:
    """For a variable in bands: which model band each year band goes into, from the counts of the patients in the
    model (`in_model`: in the analysis, with every variable under the complete-case rule). {} if none is banded."""
    banded = [t for t in coded.terms if t.form == "bands"]
    if not banded:
        return {}
    counts = coded.frame.loc[in_model, banded[0].field].value_counts().to_dict()
    return merge_bands({str(k): int(n) for k, n in counts.items()}, settings.read("adjustment.year_bands"),
                       settings.read("adjustment.band_min_patients"))


def apply_model_bands(coded: Coded, mapping: dict[str, str]) -> Coded:
    """The coded variables with year bands replaced by the model's (merged) bands, and the levels to match."""
    if not mapping:
        return coded
    frame, terms = coded.frame.copy(), []
    for t in coded.terms:
        if t.form == "bands":
            frame[t.field] = frame[t.field].map(lambda x: mapping.get(x, pd.NA) if not pd.isna(x) else pd.NA).astype("string")
            levels = tuple(dict.fromkeys(mapping[b] for b in t.levels))
            terms.append(Term(t.field, t.form, levels))
        else:
            terms.append(t)
    return Coded(frame, terms, coded.notes)


__all__ = ["Term", "Coded", "code_adjustment", "missing_data_rule", "missing_any", "merge_bands", "model_bands",
           "apply_model_bands", "band_of", "merged_label", "humanize"]
