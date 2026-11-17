"""The generic flat-file adapter: CSV/TSV, Excel (.xlsx), Parquet, SAS (.sas7bdat, .xpt), Stata (.dta), MAF and
annotated VCF (.vcf, .vcf.gz) files, read through a mapping file.

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
                                      to days), {days_since_diagnosis: {format: "%Y-%m-%d"}} for a date, or
                                      {hgvsp_short: true}: a protein change in HGVSp short form (p.Gly12Asp,
                                      or VEP's ENSP...:p.Gly12Asp, read as p.G12D)
    {constant: VALUE}                 the same value on every row (for example modality: ngs)
    {absent: true}                    the source doesn't have it

A date is read as it is when the file stores it as a date (Parquet, Stata, and SAS's .sas7bdat for a column with a
SAS date format); from text with `format`; or, where the file holds a number, with `days_from: "1960-01-01"` (or
`seconds_from` for a date and time), which is how SAS and Stata count. A SAS transport file (.xpt) always holds
dates as such numbers.

A file is a path, or {path: ..., sheet: ...} for an Excel sheet. SAS and Stata files take two more options:
    encoding: TEXT        how the file's text is encoded. A .sas7bdat file records its own (the default reads it);
                          a .xpt file doesn't, and is read as ISO-8859-1 unless this says otherwise
    value_labels: false   a Stata file's coded values as codes. By default each coded value is read as its label,
                          the text Stata itself shows, so that a mapping's `map` names readable values
SAS stores every number as a decimal, so a SAS or Stata column whose numbers are all whole is read as whole
numbers: an ID such as 1001 is then the text "1001", not "1001.0".

A MAF is read as a tab-separated table, after the comment lines (#version ...) before its header. A VCF is read
by delayed_entry.adapters.vcf into one row per variant, sample and annotation, and needs one more option, the
two reading rules for its annotation (no default):
    annotation: {pick: canonical | first | all, filters: pass | all}   (and field: CSQ or ANN if it has both)
A path may be a pattern, such as vcf/*.vcf, to read several files of one kind into one table; each row's
source_file column then names its file, without the extension.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import duckdb
import pandas as pd
import yaml

from delayed_entry.adapters import AdapterResult
from delayed_entry.adapters import vcf as vcf_reader
from delayed_entry.adapters.hgvs import short_protein_change
from delayed_entry.adapters.xlsx import read_xlsx
from delayed_entry.privacy import show
from delayed_entry.profile import ITEMS as PROFILE_ITEMS, profile_from_dict
from delayed_entry.settings import Settings
from delayed_entry.tables import OPTIONAL_TABLES, TABLES, StandardTables, conform, table_fields

MAPPING_FORMAT = 1
DRAFT, CONFIRMED = "draft", "confirmed"
SYNONYMS = Path(__file__).with_name("flatfile_synonyms.yml")
SUFFIXES = {".csv": "csv", ".tsv": "tsv", ".txt": "tsv", ".xlsx": "xlsx", ".parquet": "parquet",
            ".sas7bdat": "sas7bdat", ".xpt": "xport", ".dta": "stata", ".maf": "maf", ".vcf": "vcf"}
ENTRY_KEYS = {"column", "constant", "absent", "transform"}
FILE_OPTIONS = {"sheet", "encoding", "value_labels", "annotation"}
DATE_FORMS = ("format", "days_from", "seconds_from")
TRANSFORMS = ("map", "multiply", "days_since_diagnosis", "hgvsp_short")
PATTERN = re.compile(r"[*?\[]")


class MappingError(ValueError):
    pass


class UnconfirmedMappingError(MappingError):
    """A mapping no person has confirmed. It is never used."""


# ---------------------------------------------------------------------------
# Reading the files
# ---------------------------------------------------------------------------
def file_kind(path: Path) -> str:
    kind = "vcf" if path.name.lower().endswith(".vcf.gz") else SUFFIXES.get(path.suffix.lower())
    if kind is None:
        raise MappingError(f"{path.name}: can't read this kind of file (CSV, TSV, .xlsx, Parquet, SAS .sas7bdat or "
                           ".xpt, Stata .dta, MAF, or VCF)")
    return kind


def _comment_lines(path: Path) -> int:
    """How many lines starting with # come before a MAF's header."""
    n = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            n += 1
    return n


def _typed(df: pd.DataFrame) -> pd.DataFrame:
    """A SAS or Stata table made plain: labeled (categorical) columns as their values, and decimal columns whose
    numbers are all whole as whole numbers."""
    out = {}
    for c in df.columns:
        s = df[c]
        if isinstance(s.dtype, pd.CategoricalDtype):
            s = s.astype(object).where(s.notna(), None)
        elif pd.api.types.is_float_dtype(s) and s.notna().any() and bool((s.dropna() % 1 == 0).all()):
            s = s.astype("Int64")
        out[c] = s
    return pd.DataFrame(out, index=df.index)


def read_source(path: str | Path, sheet: str | None = None, missing_values=("",), encoding: str | None = None,
                value_labels: bool = True, annotation: dict | None = None) -> pd.DataFrame:
    """A source file as a table. CSV, TSV, MAF and VCF values are read as text; the mapping converts them. Values
    listed in missing_values become blank. encoding and value_labels apply to SAS and Stata files, annotation to a
    VCF (module notes)."""
    path = Path(path)
    kind = file_kind(path)
    if kind in ("csv", "tsv"):
        df = pd.read_csv(path, sep="," if kind == "csv" else "\t", dtype=str, keep_default_na=False, encoding="utf-8")
    elif kind == "maf":
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, encoding="utf-8",
                         skiprows=_comment_lines(path))
    elif kind == "vcf":
        a = annotation or {}
        if not a.get("pick") or not a.get("filters"):
            raise MappingError(f"{path.name}: a VCF needs its annotation's reading rules: "
                               "annotation: {pick: canonical | first | all, filters: pass | all}")
        try:
            df = vcf_reader.read_vcf(path, a["pick"], a["filters"], a.get("field"))
        except vcf_reader.VcfError as e:
            raise MappingError(str(e)) from e
    elif kind == "xlsx":
        df = read_xlsx(path, sheet)
    elif kind in ("sas7bdat", "xport"):
        df = _typed(pd.read_sas(path, format=kind, encoding=encoding or ("infer" if kind == "sas7bdat" else "ISO-8859-1")))
    elif kind == "stata":
        try:
            df = _typed(pd.read_stata(path, convert_categoricals=value_labels, convert_dates=True))
        except ValueError as e:
            raise MappingError(f"{path.name}: {e}. Give the file value_labels: false to read its codes instead") from e
    else:
        con = duckdb.connect()
        try:
            df = con.execute("SELECT * FROM read_parquet(?)", [str(path)]).df()
        finally:
            con.close()
    missing = set(missing_values)
    if len(df):
        df = df.copy()
        for c in df.columns:       # only text can be a missing-value word; other columns keep their types
            if df[c].dtype == object or pd.api.types.is_string_dtype(df[c]):
                df[c] = df[c].map(lambda v: pd.NA if isinstance(v, str) and v in missing else v)
    return df


def to_dates(series: pd.Series, spec: dict) -> pd.Series:
    """A column as dates: dates as they are; text read with spec's `format`; or a count of days (`days_from`) or of
    seconds (`seconds_from`) after an origin date. A value that doesn't convert is blank."""
    if spec.get("days_from") or spec.get("seconds_from"):
        unit, origin = ("D", spec["days_from"]) if spec.get("days_from") else ("s", spec["seconds_from"])
        return pd.Timestamp(str(origin)) + pd.to_timedelta(pd.to_numeric(series, errors="coerce"), unit=unit)
    if pd.api.types.is_datetime64_any_dtype(series) or \
            pd.api.types.infer_dtype(series, skipna=True) in ("datetime", "datetime64", "date"):
        return pd.to_datetime(series, errors="coerce")
    return pd.to_datetime(series, format=spec.get("format"), errors="coerce")


def _date_form_problems(where: str, spec: dict) -> list[str]:
    forms = [k for k in DATE_FORMS if spec.get(k) is not None]
    if len(forms) > 1:
        return [f"{where}: give one of format, days_from or seconds_from"]
    for k in ("days_from", "seconds_from"):
        if spec.get(k) is not None:
            try:
                pd.Timestamp(str(spec[k]))
            except ValueError:
                return [f"{where}: {k} must be a date, for example \"1960-01-01\""]
    return []


def _norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


# ---------------------------------------------------------------------------
# The draft
# ---------------------------------------------------------------------------
def _file_entry(spec) -> tuple[str, dict]:
    """A file's path and its options (sheet, encoding, value_labels, annotation)."""
    if isinstance(spec, dict):
        return spec["path"], {k: v for k, v in spec.items() if k in FILE_OPTIONS}
    return str(spec), {}


def file_paths(p: str, base: Path) -> list[Path]:
    """The file a path names, or every file a pattern (with *, ? or [) matches, in name order."""
    path = Path(p) if Path(p).is_absolute() else base / p
    if not PATTERN.search(p):
        return [path]
    found = sorted(f for f in Path(path.anchor).glob(str(path.relative_to(path.anchor))) if f.is_file())
    if not found:
        raise MappingError(f"no file matches {p}")
    return found


def _source_name(path: Path) -> str:
    name = path.name
    return vcf_reader.stem(path) if file_kind(path) == "vcf" else name[: -len(path.suffix)] if path.suffix else name


def _read_entry(path: Path, opts: dict, missing_values=("",), drafting: bool = False) -> pd.DataFrame:
    annotation = opts.get("annotation")
    if drafting and annotation is None:                # a draft only needs the columns
        annotation = {"pick": "all", "filters": "all"}
    return read_source(path, opts.get("sheet"), missing_values, encoding=opts.get("encoding"),
                       value_labels=opts.get("value_labels", True), annotation=annotation)


def read_files(p: str, opts: dict, base: Path, missing_values=("",), drafting: bool = False) -> tuple[pd.DataFrame, list[Path]]:
    """The table a file entry names, and the files read. A pattern's files are read into one table, with each row's
    file named in source_file."""
    paths = file_paths(p, base)
    if not PATTERN.search(p):
        return _read_entry(paths[0], opts, missing_values, drafting), paths
    parts = [_read_entry(f, opts, missing_values, drafting).assign(source_file=_source_name(f)) for f in paths]
    return pd.concat(parts, ignore_index=True), paths


def draft_mapping(files: dict[str, str | dict], base_dir: str | Path | None = None) -> str:
    """The text of a draft mapping for these files ({standard table: path or {path, sheet}}), with each
    standard field matched to a column by name where one fits, and the capability profile left blank."""
    synonyms = yaml.safe_load(SYNONYMS.read_text(encoding="utf-8"))
    base = Path(base_dir) if base_dir else Path(".")
    q = json.dumps
    out = ["# Draft mapping, made by matching column names. It is NOT used until a person has checked every",
           "# field and the profile, and set status: confirmed, confirmed_by and confirmed_on.",
           "# Each field: {column: NAME}, optionally with transform: {map: {...}} | {multiply: n} |",
           "# {days_since_diagnosis: {format: ...}} | {hgvsp_short: true}; or {constant: VALUE}; or {absent: true}.",
           "# A date held as a number (SAS .xpt files, and SAS or Stata dates read as numbers) takes",
           "# {days_from: \"1960-01-01\"} or {seconds_from: ...} in place of format.",
           f"mapping_format: {MAPPING_FORMAT}", f"status: {DRAFT}", "confirmed_by: null", "confirmed_on: null",
           'missing_values: ["", "NA"]',
           "# For dates: the patients file's diagnosis date, {column: NAME, format: \"%Y-%m-%d\"} (or days_from or",
           "# seconds_from), or null.",
           "diagnosis_date: null", "files:"]
    for t in TABLES:
        if t in files:
            p, opts = _file_entry(files[t])
            vcf_rules = file_kind(Path(p)) == "vcf" and "annotation" not in opts
            if vcf_rules:                 # the reading rules a person must choose: there is no default
                opts = {**opts, "annotation": {"pick": None, "filters": None}}
            out.append(f"  {t}: {{path: {q(p)}" + "".join(f", {k}: {q(v)}" for k, v in opts.items()) + "}"
                       + ("    # pick: canonical / first / all; filters: pass / all" if vcf_rules else ""))
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
        p, opts = _file_entry(files[t])
        cols = list(read_files(p, opts, base, drafting=True)[0].columns)
        by_norm = {_norm(c): c for c in cols}
        out.append(f"  {t}:")
        from_vcf = t == "results" and file_kind(Path(p)) == "vcf"
        for f in table_fields(t):
            names = synonyms["vcf_results"].get(f.name, []) if from_vcf else \
                [f.name] + synonyms.get(t, {}).get(f.name, [])
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
def _annotation_problems(where: str, p: str, a) -> list[str]:
    """A VCF must give its annotation's reading rules; nothing else may give any."""
    try:
        is_vcf = file_kind(Path(p)) == "vcf"
    except MappingError as e:
        return [f"{where}: {e}"]
    if not is_vcf:
        return [f"{where}: annotation applies only to a VCF"] if a is not None else []
    rules = f"annotation: {{pick: {' | '.join(vcf_reader.PICKS)}, filters: {' | '.join(vcf_reader.FILTERS)}}}"
    if not isinstance(a, dict) or not set(a) <= {"pick", "filters", "field"}:
        return [f"{where}: a VCF needs its annotation's reading rules, {rules} (and field: CSQ or ANN if it has both)"]
    out = []
    if a.get("pick") not in vcf_reader.PICKS:
        out.append(f"{where}: annotation.pick must be one of {', '.join(vcf_reader.PICKS)} (no default: it is a "
                   "reading rule the mapping's author chooses)")
    if a.get("filters") not in vcf_reader.FILTERS:
        out.append(f"{where}: annotation.filters must be one of {', '.join(vcf_reader.FILTERS)}")
    if a.get("field") is not None and a["field"] not in vcf_reader.ANNOTATION_FIELDS:
        out.append(f"{where}: annotation.field must be one of {', '.join(vcf_reader.ANNOTATION_FIELDS)}")
    return out


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
        if spec and not (isinstance(spec, str) or (isinstance(spec, dict) and isinstance(spec.get("path"), str)
                                                   and set(spec) <= {"path"} | FILE_OPTIONS)):
            problems.append(f"files.{t}: give a path, or {{path: ..., and the file's options: sheet, encoding, "
                            "value_labels, annotation}}")
        elif spec:
            p, opts = _file_entry(spec)
            if "value_labels" in opts and not isinstance(opts["value_labels"], bool):
                problems.append(f"files.{t}: value_labels must be true or false")
            if "encoding" in opts and not (isinstance(opts["encoding"], str) and opts["encoding"].strip()):
                problems.append(f"files.{t}: encoding must name an encoding, for example latin-1")
            problems += _annotation_problems(f"files.{t}", p, opts.get("annotation"))
    dd = m.get("diagnosis_date")
    if dd is not None and not (isinstance(dd, dict) and isinstance(dd.get("column"), str)
                               and isinstance(dd.get("format", ""), str) and set(dd) <= {"column", *DATE_FORMS}):
        problems.append("diagnosis_date: give {column: NAME, format: \"%Y-%m-%d\"} (or days_from or seconds_from "
                        "in place of format)")
    elif dd is not None:
        problems += _date_form_problems("diagnosis_date", dd)
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
            if tr is not None and not (isinstance(tr, dict) and len(tr) == 1 and next(iter(tr)) in TRANSFORMS):
                problems.append(f"{where}: transform must be one of {', '.join(TRANSFORMS)}")
                continue
            if tr and "hgvsp_short" in tr and tr["hgvsp_short"] is not True:
                problems.append(f"{where}: hgvsp_short takes true")
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
                if arg is not None and not (isinstance(arg, dict) and set(arg) <= set(DATE_FORMS)):
                    problems.append(f"{where}: days_since_diagnosis takes {{format: ...}}, {{days_from: ...}} or "
                                    "{seconds_from: ...}")
                elif arg is not None:
                    problems += _date_form_problems(f"{where}: days_since_diagnosis", arg)
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
    elif kind == "hgvsp_short":
        out = series.map(short_protein_change)
    else:
        when = to_dates(series, arg or {})
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
            p, opts = _file_entry(spec)
            raw[t], paths = read_files(p, opts, base, missing)
            used += paths
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
        diag_by_patient = pd.Series(to_dates(pt[dd["column"]], dd).to_numpy(), index=pt[pid_col].astype(str).to_numpy())
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
