"""The generic flat-file adapter: CSV/TSV, Excel (.xlsx) and Parquet files, read through a mapping file.

    from delayed_entry.adapters import flatfile
    flatfile.write_draft({"patients": "patients.csv", "tests": "tests.xlsx", ...}, "mapping.yaml")
    # a person reads mapping.yaml, fixes every field, fills in the capability profile, and sets
    # status: confirmed, confirmed_by and confirmed_on
    result = flatfile.read("mapping.yaml", settings)

The tool drafts the mapping by matching column names (flatfile_synonyms.yml). A draft is never used: read()
refuses a mapping whose status isn't "confirmed", that names no person, or that leaves any field undecided.
Each standard field is mapped to one of:
    {column: NAME}                    a column, as it is
    {column: NAME, transform: ...}    a column, changed: {map: {from: to}}, {multiply: n} (for example months
                                      to days) or {days_since_diagnosis: {format: "%Y-%m-%d"}} for a date
    {constant: VALUE}                 the same value on every row (for example modality: ngs)
    {absent: true}                    the source doesn't have it
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import duckdb
import pandas as pd
import yaml

from delayed_entry.adapters import AdapterResult
from delayed_entry.adapters.xlsx import read_xlsx
from delayed_entry.privacy import show
from delayed_entry.profile import ITEMS as PROFILE_ITEMS, profile_from_dict
from delayed_entry.settings import Settings
from delayed_entry.tables import OPTIONAL_TABLES, TABLES, StandardTables, conform, table_fields

MAPPING_FORMAT = 1
DRAFT, CONFIRMED = "draft", "confirmed"
SYNONYMS = Path(__file__).with_name("flatfile_synonyms.yml")
SUFFIXES = {".csv": "csv", ".tsv": "tsv", ".txt": "tsv", ".xlsx": "xlsx", ".parquet": "parquet"}
ENTRY_KEYS = {"column", "constant", "absent", "transform"}


class MappingError(ValueError):
    pass


class UnconfirmedMappingError(MappingError):
    """A mapping no person has confirmed. It is never used."""


# ---------------------------------------------------------------------------
# Reading the files
# ---------------------------------------------------------------------------
def file_kind(path: Path) -> str:
    kind = SUFFIXES.get(path.suffix.lower())
    if kind is None:
        raise MappingError(f"{path.name}: can't read this kind of file (CSV, TSV, .xlsx or Parquet)")
    return kind


def read_source(path: str | Path, sheet: str | None = None, missing_values=("",)) -> pd.DataFrame:
    """A source file as a table. CSV and TSV values are read as text; the mapping converts them. Values
    listed in missing_values become blank."""
    path = Path(path)
    kind = file_kind(path)
    if kind in ("csv", "tsv"):
        df = pd.read_csv(path, sep="," if kind == "csv" else "\t", dtype=str, keep_default_na=False, encoding="utf-8")
    elif kind == "xlsx":
        df = read_xlsx(path, sheet)
    else:
        con = duckdb.connect()
        try:
            df = con.execute("SELECT * FROM read_parquet(?)", [str(path)]).df()
        finally:
            con.close()
    missing = set(missing_values)
    return df.apply(lambda s: s.map(lambda v: pd.NA if isinstance(v, str) and v in missing else v)) if len(df) else df


def _norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


# ---------------------------------------------------------------------------
# The draft
# ---------------------------------------------------------------------------
def _file_entry(spec) -> tuple[str, str | None]:
    if isinstance(spec, dict):
        return spec["path"], spec.get("sheet")
    return str(spec), None


def draft_mapping(files: dict[str, str | dict], base_dir: str | Path | None = None) -> str:
    """The text of a draft mapping for these files ({standard table: path or {path, sheet}}), with each
    standard field matched to a column by name where one fits, and the capability profile left blank."""
    synonyms = yaml.safe_load(SYNONYMS.read_text(encoding="utf-8"))
    base = Path(base_dir) if base_dir else Path(".")
    q = json.dumps
    out = ["# Draft mapping, made by matching column names. It is NOT used until a person has checked every",
           "# field and the profile, and set status: confirmed, confirmed_by and confirmed_on.",
           "# Each field: {column: NAME}, optionally with transform: {map: {...}} | {multiply: n} |",
           "# {days_since_diagnosis: {format: ...}}; or {constant: VALUE}; or {absent: true}.",
           f"mapping_format: {MAPPING_FORMAT}", f"status: {DRAFT}", "confirmed_by: null", "confirmed_on: null",
           'missing_values: ["", "NA"]',
           "# For dates: the patients file's diagnosis date, {column: NAME, format: \"%Y-%m-%d\"}, or null.",
           "diagnosis_date: null", "files:"]
    for t in TABLES:
        if t in files:
            p, sheet = _file_entry(files[t])
            out.append(f"  {t}: {{path: {q(p)}" + (f", sheet: {q(sheet)}" if sheet else "") + "}")
        else:
            out.append(f"  {t}: null" + ("    # optional table" if t in OPTIONAL_TABLES else "    # required: give its file"))
    out += ["# The capability profile: what this source can supply. Every item must be filled in.", "profile:"]
    for item, options in PROFILE_ITEMS.items():
        out.append(f"  {item}: null    # {' / '.join(options)}")
    out += ["  privacy_min_count: null    # the data-use agreement's minimum count, or null if it sets none",
            "  histology_vocabulary: null    # what the histology field holds, for example ICD-O-3",
            "  age_top_code: null    # ages above this are withheld or top-coded; null if none",
            "  known_quirks: []    # each {name: ..., handling: ...}",
            "  sources: {}    # each item: the source behind its value",
            "  attribution: null    # or {citation: ..., version: ..., acknowledgment: ...}, if the source's terms "
            "require one",
            "tables:"]
    for t in TABLES:
        if t not in files:
            continue
        p, sheet = _file_entry(files[t])
        cols = list(read_source(base / p, sheet).columns)
        by_norm = {_norm(c): c for c in cols}
        out.append(f"  {t}:")
        for f in table_fields(t):
            names = [f.name] + synonyms.get(t, {}).get(f.name, [])
            hit = next((by_norm[_norm(n)] for n in names if _norm(n) in by_norm), None)
            if hit is not None:
                out.append(f"    {f.name}: {{column: {q(hit)}}}    # drafted: matched the column name; check it")
            else:
                out.append(f"    {f.name}: {{column: null}}    # not found: give a column, a constant, or absent: true")
    return "\n".join(out) + "\n"


def write_draft(files: dict[str, str | dict], out_path: str | Path) -> Path:
    out_path = Path(out_path)
    out_path.write_bytes(draft_mapping(files, out_path.parent).encode("utf-8"))
    return out_path


# ---------------------------------------------------------------------------
# The confirmed mapping
# ---------------------------------------------------------------------------
def load_mapping(path: str | Path) -> dict:
    """A confirmed, complete mapping. Raises UnconfirmedMappingError for a draft, MappingError otherwise."""
    path = Path(path)
    m = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(m, dict) or m.get("mapping_format") != MAPPING_FORMAT:
        raise MappingError(f"{path.name}: not a mapping file (mapping_format {MAPPING_FORMAT})")
    if m.get("status") != CONFIRMED or not (isinstance(m.get("confirmed_by"), str) and m["confirmed_by"].strip()) \
            or not m.get("confirmed_on"):
        raise UnconfirmedMappingError(f"{path.name} has not been confirmed: a person must check it and set status: "
                                      "confirmed, confirmed_by and confirmed_on. A drafted mapping is never used.")
    problems = []
    files = m.get("files") or {}
    tables = m.get("tables") or {}
    if not isinstance(files, dict) or not isinstance(tables, dict):
        raise MappingError(f"{path.name}: files and tables must be mappings")
    for section, names in (("files", files), ("tables", tables)):
        for t in set(names) - set(TABLES):
            problems.append(f"{section}: '{t}' is not a standard table ({', '.join(TABLES)})")
    for t in TABLES:
        if t not in OPTIONAL_TABLES and not files.get(t):
            problems.append(f"files: the {t} table needs a file")
    for t, spec in files.items():
        if spec and not (isinstance(spec, str) or (isinstance(spec, dict) and isinstance(spec.get("path"), str))):
            problems.append(f"files.{t}: give a path, or {{path: ..., sheet: ...}}")
    dd = m.get("diagnosis_date")
    if dd is not None and not (isinstance(dd, dict) and isinstance(dd.get("column"), str)
                               and isinstance(dd.get("format", ""), str) and set(dd) <= {"column", "format"}):
        problems.append("diagnosis_date: give {column: NAME, format: \"%Y-%m-%d\"}")
    if dd is not None and not (isinstance((tables.get("patients") or {}).get("patient_id"), dict)
                               and tables["patients"]["patient_id"].get("column")):
        problems.append("diagnosis_date: the patients table's patient_id must be a column")
    for t, spec in files.items():
        if not spec or t not in TABLES:
            continue
        fields = tables.get(t) or {}
        for f in table_fields(t):
            e = fields.get(f.name)
            if not isinstance(e, dict) or not (set(e) <= ENTRY_KEYS):
                problems.append(f"tables.{t}.{f.name}: must be {{column: ...}}, {{constant: ...}} or {{absent: true}}")
                continue
            decided = (e.get("column") is not None) + ("constant" in e) + (e.get("absent") is True)
            if decided != 1:
                problems.append(f"tables.{t}.{f.name}: not decided (give exactly one of a column, a constant, "
                                "or absent: true)")
            if f.key and e.get("absent") is True:
                problems.append(f"tables.{t}.{f.name}: the table can't be used without it")
            tr = e.get("transform")
            where = f"tables.{t}.{f.name}"
            if tr is not None and not (isinstance(tr, dict) and len(tr) == 1
                                       and next(iter(tr)) in ("map", "multiply", "days_since_diagnosis")):
                problems.append(f"{where}: transform must be one of map, multiply, days_since_diagnosis")
                continue
            if tr and "map" in tr:
                mp = tr["map"]
                if not isinstance(mp, dict) or not mp:
                    problems.append(f"{where}: map must be a mapping {{source value: standard value}}")
                elif any(isinstance(k, bool) or not isinstance(k, (str, int, float)) for k in mp):
                    problems.append(f"{where}: a map key was read as true or false: put it in quotes (for example 'Yes')")
            if tr and "multiply" in tr and (isinstance(tr["multiply"], bool) or not isinstance(tr["multiply"], (int, float))):
                problems.append(f"{where}: multiply needs a number")
            if tr and "days_since_diagnosis" in tr:
                arg = tr["days_since_diagnosis"]
                if arg is not None and not (isinstance(arg, dict) and set(arg) <= {"format"}):
                    problems.append(f"{where}: days_since_diagnosis takes {{format: ...}}")
                if not m.get("diagnosis_date"):
                    problems.append(f"{where}: days_since_diagnosis needs diagnosis_date")
                if t != "patients" and not (isinstance(fields.get("patient_id"), dict) and fields["patient_id"].get("column")):
                    problems.append(f"{where}: days_since_diagnosis needs this table's patient_id as a column")
        unknown = set(fields) - {f.name for f in table_fields(t)}
        if unknown:
            problems.append(f"tables.{t}: unknown fields {sorted(unknown)}")
    try:
        profile_from_dict(m.get("profile"), adapter="flat_file")
    except ValueError as e:
        problems.append(str(e))
    if problems:
        raise MappingError(f"{path.name}: " + "; ".join(problems))
    return m


def _apply(series: pd.Series, transform: dict | None, diag: pd.Series | None) -> tuple[pd.Series, int]:
    """The column after its transform, and how many filled values the transform couldn't handle."""
    if not transform:
        return series, 0
    (kind, arg), = transform.items()
    filled = series.notna()
    if kind == "map":
        out = series.map(lambda v: arg.get(v, arg.get(str(v), pd.NA)) if pd.notna(v) else pd.NA)
    elif kind == "multiply":
        out = pd.to_numeric(series, errors="coerce") * float(arg)
    else:
        fmt = (arg or {}).get("format")
        when = pd.to_datetime(series, format=fmt, errors="coerce")
        out = (when - diag).dt.days.astype("Float64") if diag is not None else pd.Series(pd.NA, index=series.index)
    return out, int((filled & pd.Series(out).isna()).sum())


def read(mapping_path: str | Path, settings: Settings) -> AdapterResult:
    """Read the files a confirmed mapping names into the standard tables. The settings give the privacy
    minimum used to mask the counts in what reading reports."""
    mapping_path = Path(mapping_path)
    m = load_mapping(mapping_path)
    base = mapping_path.parent
    prof = profile_from_dict(m["profile"], adapter="flat_file")
    min_count = max(settings.read("privacy.min_count"), prof.privacy_min_count or 0)

    def count(n: int) -> str:
        return show(n, min_count)
    missing = tuple(m.get("missing_values") or ("",))
    raw: dict[str, pd.DataFrame] = {}
    used: list[Path] = [mapping_path]
    for t, spec in (m.get("files") or {}).items():
        if spec:
            p, sheet = _file_entry(spec)
            path = (base / p) if not Path(p).is_absolute() else Path(p)
            raw[t] = read_source(path, sheet, missing)
            used.append(path)
    notes = []
    diag_by_patient = None
    dd = m.get("diagnosis_date")
    if dd:
        pid_col = m["tables"]["patients"]["patient_id"]["column"]
        pt = raw["patients"]
        for col in (pid_col, dd["column"]):
            if col not in pt.columns:
                raise MappingError(f"no column '{col}' in the patients file")
        dup = int(pt[pid_col].astype(str).duplicated().sum())
        if dup:
            raise MappingError(f"{count(dup)} patient IDs appear more than once in the patients file, so their "
                               "diagnosis dates can't be told apart")
        diag_by_patient = pd.Series(pd.to_datetime(pt[dd["column"]], format=dd.get("format"), errors="coerce").to_numpy(),
                                    index=pt[pid_col].astype(str).to_numpy())
    out: dict[str, pd.DataFrame | None] = {t: None for t in TABLES}
    for t, df in raw.items():
        spec = m["tables"][t]
        built = {}
        diag = None
        if diag_by_patient is not None and "patient_id" in spec and spec["patient_id"].get("column"):
            diag = df[spec["patient_id"]["column"]].astype(str).map(diag_by_patient)
        for f in table_fields(t):
            e = spec[f.name]
            if e.get("absent") is True:
                continue
            if "constant" in e:
                built[f.name] = pd.Series([e["constant"]] * len(df), index=df.index, dtype="object")
                continue
            if e["column"] not in df.columns:
                raise MappingError(f"tables.{t}.{f.name}: no column '{e['column']}' in the {t} file")
            values, failed = _apply(df[e["column"]], e.get("transform"), diag)
            if failed:
                notes.append(f"{t}.{f.name}: {count(failed)} values the transform couldn't convert were read as blank")
            built[f.name] = values
        out[t] = pd.DataFrame(built, index=df.index).reset_index(drop=True)
    tables, failures = conform(StandardTables(**out))
    notes += [f"{k}: {count(n)} values didn't convert to the field's type and were read as blank"
              for k, n in failures.items()]
    return AdapterResult(tables=tables, profile=prof, input_files=used, source_steps=[], read_notes=notes)
