"""The OMOP CDM adapter: an OMOP Common Data Model (version 5.4) database into the standard tables, read
through a mapping a person has confirmed.

    from delayed_entry.adapters import omop
    omop.write_template("omop_mapping.yaml")   # a person fills it in, checks it and confirms it
    result = omop.read("omop_mapping.yaml", settings)

OMOP codes everything with concepts. Which concepts mean the study's cancer, a sequencing test or an MSI
result is the study's choice, so the mapping names them as concept sets: nothing is assumed. The mapping also
says where each patient field comes from, which genomic extension table holds the variants, and the source's
capability profile. A mapping is never used until a person has confirmed it.

The database is either a folder holding one file per CDM table (person.csv, death.parquet, ...; CSV, TSV or
Parquet) or a DuckDB file with the tables in one schema.

What the adapter works out itself, from the dates:
- the index diagnosis: each person's earliest condition in the diagnosis concept set (the lowest
  condition_occurrence_id on a tie). Persons with none are not patients: the flow diagram's first step;
- every time, as days from that diagnosis;
- the outcome: death from the death table; otherwise alive at the end of the person's last observation period;
- a test's after_follow_up_end: its date is on or after the end of follow-up;
- regimens: a person's drug exposures in the drug concept set that start on the same day form one regimen,
  numbered by start date, since OMOP keeps no regimen number (the CRC plan, section 4d).

A patient field read from the measurement or observation table takes the person's row dated closest to the
diagnosis, the earlier on a tie, then the lowest row ID.
"""
from __future__ import annotations

import re
from pathlib import Path

import duckdb
import pandas as pd
import yaml

from delayed_entry.adapters import AdapterResult
from delayed_entry.adapters.flatfile import read_source as read_flat_file
from delayed_entry.adapters.hgvs import AMINO_ACIDS, short_protein_change  # noqa: F401 (AMINO_ACIDS: the tests' fixture)
from delayed_entry.privacy import show
from delayed_entry.profile import ITEMS as PROFILE_ITEMS, profile_from_dict
from delayed_entry.settings import Settings
from delayed_entry.tables import (ENTRY_EVENTS, RESULT_KINDS, TABLES, StandardTables, conform, empty_table,
                                  table_fields)

MAPPING_FORMAT = 1
KIND = "omop_cdm"
DRAFT, CONFIRMED = "draft", "confirmed"
CDM_SUFFIXES = (".csv", ".tsv", ".parquet")
# The patient fields the adapter always works out itself, never from the mapping.
DERIVED = ("patient_id", "vital_status", "days_to_death_or_last_followup")
MAPPED = tuple(f.name for f in table_fields("patients") if f.name not in DERIVED)
COMPUTABLE = {"age_at_diagnosis": "the year of diagnosis less the year of birth (OHDSI's convention)",
              "year_of_diagnosis": "the calendar year of the diagnosis date"}
SOURCE_KEYS = ("person", "care_site", "index_condition", "measurement", "observation", "computed", "absent")
EVENT_TABLES = {"measurement": ("measurement_concept_id", "measurement_date", "measurement_id"),
                "observation": ("observation_concept_id", "observation_date", "observation_id")}
DATE_COLUMNS = {
    "condition_occurrence": ("condition_start_date",),
    "death": ("death_date",),
    "observation_period": ("observation_period_start_date", "observation_period_end_date"),
    "procedure_occurrence": ("procedure_date",),
    "measurement": ("measurement_date",),
    "observation": ("observation_date",),
    "drug_exposure": ("drug_exposure_start_date",),
}
VARIANT_KEYS = {"table": True, "test_link": True, "gene": True, "protein_change": True, "partner_gene": False,
                "fusion_description": False, "variant_class": False, "variant_class_labels": False,
                "copy_number": False}
READING_KINDS = tuple(k for k in RESULT_KINDS if k not in ("mutation", "fusion"))
VALUE_COLUMNS = ("value_as_concept_id", "value_as_string", "value_as_number")
DRUG_NAMES = ("drug_source_value", "concept_name")


class MappingError(ValueError):
    pass


class UnconfirmedMappingError(MappingError):
    """A mapping no person has confirmed. It is never used."""


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def _key(v):
    """A code or ID as text, so that 8507, 8507.0 and "8507" are the same key."""
    if v is None or (not isinstance(v, (list, dict)) and pd.isna(v)):
        return pd.NA
    if isinstance(v, bool):
        return str(v).lower()
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v).strip()


def _keys(s: pd.Series) -> pd.Series:
    return s.map(_key).astype("string")


def _days(when: pd.Series, dx: pd.Series) -> pd.Series:
    return (when - dx).dt.days.astype("Float64")


# ---------------------------------------------------------------------------
# The template
# ---------------------------------------------------------------------------
def template_text() -> str:
    """A blank mapping, with every choice to make listed and explained."""
    out = [
        "# An OMOP CDM mapping. It is NOT used until a person has filled in every choice, checked it against the",
        "# database, and set status: confirmed, confirmed_by and confirmed_on.",
        f"mapping_format: {MAPPING_FORMAT}", f"kind: {KIND}", f"status: {DRAFT}", "confirmed_by: null",
        "confirmed_on: null",
        "# Where the CDM tables are: a folder with one file per table (person.csv, ...; CSV, TSV or Parquet), with",
        "# the format of the dates in its CSV or TSV files; or a DuckDB file and the schema holding the tables.",
        "cdm: {folder: null, date_format: \"%Y-%m-%d\"}    # or {duckdb: FILE, schema: cdm}",
        "# Concept sets, by name: each {concepts: [concept IDs], include_descendants: true or false}.",
        "# include_descendants needs the concept_ancestor table.",
        "concept_sets: {}",
        "patient_id: null    # person_id or person_source_value",
        "diagnosis: null    # the concept set of the study's cancer: the earliest such condition is the diagnosis",
        "# Where each patient field comes from: one of",
        "#   {person: COLUMN}, {care_site: COLUMN} (through person.care_site_id), {index_condition: COLUMN},",
        "#   {measurement: CONCEPT_SET, value: COLUMN} or {observation: CONCEPT_SET, value: COLUMN} (the row dated",
        "#   closest to diagnosis), {computed: true} (age_at_diagnosis, year_of_diagnosis only), {absent: true};",
        "# optionally with extract: REGEX (its one group is kept) and labels: {source value: standard value}.",
        "# labels must list every value found; list one as null to read it as blank.",
        "patients:",
    ]
    for name in MAPPED:
        hint = "    # or {computed: true}" if name in COMPUTABLE else ""
        out.append(f"  {name}: null{hint}")
    out += [
        "extras: {}    # more patient fields the settings name: {name: a source as above}",
        "# The sequencing tests: each procedure_occurrence row in the concept set is a test.",
        "sequencing:",
        "  procedures: null    # concept set",
        "  modality: null    # as the settings' entry.modalities name it, for example ngs",
        f"  entry_event: null    # what the procedure date is: {' / '.join(ENTRY_EVENTS)}",
        "  panel_id: null    # {procedure: COLUMN} with extract or labels if needed, or {absent: true}",
        "  # The genomic extension table with one row per variant, linked to the test; null if there is none.",
        "  variants:",
        "    table: null",
        "    test_link: null    # its column holding the procedure_occurrence_id",
        "    gene: null    # its column holding the gene symbol (HGNC)",
        "    protein_change: null    # its column holding the protein change (HGVSp, either form)",
        "    partner_gene: null    # a fusion's second gene, or null",
        "    fusion_description: null    # or null: a fusion is then described as GENE1-GENE2",
        "    variant_class: null    # or null",
        "    variant_class_labels: null    # {source value: class}, for example Missense_Mutation; or null",
        "    copy_number: null    # a column whose filled value makes the row a copy-number result; or null",
        "# Readings such as MSI or MMR status, from the measurement table: each row is a test with one result.",
        f"# Each: {{measurement: CONCEPT_SET, kind: {' / '.join(READING_KINDS)}, gene: null or a gene,",
        f"#        value: {' / '.join(VALUE_COLUMNS)}, labels: {{...}} or null, modality: ..., entry_event: ...}}",
        "biomarker_readings: []",
        "panels_file: null    # a file with panel_id, gene (and codon_from, codon_to), or null",
        "panel_sizes_file: null    # a file with panel_id, coding_megabases, or null",
        f"treatments: null    # or {{drugs: CONCEPT_SET, name: {' or '.join(DRUG_NAMES)}}}",
        "# The capability profile: what this source can supply. Every item must be filled in.",
        "profile:",
        "  adapter: omop",
    ]
    for item, options in PROFILE_ITEMS.items():
        out.append(f"  {item}: null    # {' / '.join(options)}")
    out += ["  privacy_min_count: null    # the data-use agreement's minimum count, or null if it sets none",
            "  histology_vocabulary: null    # what the histology field holds, for example ICD-O-3",
            "  age_top_code: null    # ages above this are withheld or top-coded; null if none",
            "  known_quirks: []    # each {name: ..., handling: ...}",
            "  sources: {}    # each item: the source behind its value",
            "  attribution: null    # or {citation: ..., version: ..., acknowledgment: ...}, if the source's terms "
            "require one"]
    return "\n".join(out) + "\n"


def write_template(out_path: str | Path) -> Path:
    out_path = Path(out_path)
    out_path.write_bytes(template_text().encode("utf-8"))
    return out_path


# ---------------------------------------------------------------------------
# The confirmed mapping
# ---------------------------------------------------------------------------
def _source_problems(where: str, spec, sets: dict, field: str | None) -> list[str]:
    if not isinstance(spec, dict):
        return [f"{where}: must be a mapping, for example {{person: gender_concept_id}} or {{absent: true}}"]
    sources = [k for k in SOURCE_KEYS if k in spec]
    if len(sources) != 1:
        return [f"{where}: give exactly one of {', '.join(SOURCE_KEYS)}"]
    src = sources[0]
    allowed = {src, "extract", "labels"} | ({"value"} if src in EVENT_TABLES else set())
    out = [f"{where}: unknown keys {sorted(set(spec) - allowed)}"] if set(spec) - allowed else []
    if src == "absent" and spec["absent"] is not True:
        out.append(f"{where}: absent must be true")
    if src == "computed":
        if spec["computed"] is not True:
            out.append(f"{where}: computed must be true")
        if field not in COMPUTABLE:
            out.append(f"{where}: only {', '.join(COMPUTABLE)} can be computed")
    if src in ("person", "care_site", "index_condition") and not isinstance(spec[src], str):
        out.append(f"{where}: {src} needs a column name")
    if src in EVENT_TABLES:
        if spec[src] not in sets:
            out.append(f"{where}: '{spec[src]}' is not a concept set")
        if not isinstance(spec.get("value"), str):
            out.append(f"{where}: give the {src} table's column to read as value")
    if src in ("absent", "computed") and ({"extract", "labels"} & set(spec)):
        out.append(f"{where}: extract and labels apply only to a value read from a table")
    if "extract" in spec:
        try:
            if re.compile(str(spec["extract"])).groups != 1:
                out.append(f"{where}: extract needs exactly one group in parentheses")
        except re.error as e:
            out.append(f"{where}: extract is not a regular expression ({e})")
    out += _labels_problems(where, spec.get("labels"))
    return out


def _labels_problems(where: str, labels) -> list[str]:
    if labels is None:
        return []
    if not isinstance(labels, dict) or not labels:
        return [f"{where}: labels must be a mapping {{source value: standard value}}"]
    if any(isinstance(k, bool) for k in labels):
        return [f"{where}: a labels key was read as true or false: put it in quotes (for example 'Yes')"]
    return []


def _file_spec_ok(spec) -> bool:
    return spec is None or isinstance(spec, str) or (isinstance(spec, dict) and isinstance(spec.get("path"), str)
                                                     and set(spec) <= {"path", "sheet"})


def load_mapping(path: str | Path) -> dict:
    """A confirmed, complete OMOP mapping. Raises UnconfirmedMappingError for a draft, MappingError otherwise."""
    path = Path(path)
    m = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(m, dict) or m.get("mapping_format") != MAPPING_FORMAT or m.get("kind") != KIND:
        raise MappingError(f"{path.name}: not an OMOP mapping file (mapping_format {MAPPING_FORMAT}, kind {KIND})")
    if m.get("status") != CONFIRMED or not (isinstance(m.get("confirmed_by"), str) and m["confirmed_by"].strip()) \
            or not m.get("confirmed_on"):
        raise UnconfirmedMappingError(f"{path.name} has not been confirmed: a person must check it and set status: "
                                      "confirmed, confirmed_by and confirmed_on. A drafted mapping is never used.")
    problems = []
    known = {"mapping_format", "kind", "status", "confirmed_by", "confirmed_on", "cdm", "concept_sets", "patient_id",
             "diagnosis", "patients", "extras", "sequencing", "biomarker_readings", "panels_file", "panel_sizes_file",
             "treatments", "profile"}
    if set(m) - known:
        problems.append(f"unknown keys {sorted(set(m) - known)}")
    cdm = m.get("cdm")
    if not isinstance(cdm, dict) or (("folder" in cdm) == ("duckdb" in cdm)):
        problems.append("cdm: give {folder: FOLDER, date_format: ...} or {duckdb: FILE, schema: ...}")
    elif "folder" in cdm and not (isinstance(cdm["folder"], str) and isinstance(cdm.get("date_format"), str)
                                  and set(cdm) <= {"folder", "date_format"}):
        problems.append("cdm: a folder needs its path and date_format (the format of dates in CSV and TSV files)")
    elif "duckdb" in cdm and not (isinstance(cdm["duckdb"], str) and isinstance(cdm.get("schema"), str)
                                  and set(cdm) <= {"duckdb", "schema"}):
        problems.append("cdm: a DuckDB file needs its path and the schema holding the tables")
    sets = m.get("concept_sets")
    if not isinstance(sets, dict):
        problems.append("concept_sets: must be a mapping of names to {concepts: [...], include_descendants: ...}")
        sets = {}
    for name, cs in sets.items():
        if not (isinstance(cs, dict) and set(cs) == {"concepts", "include_descendants"}
                and isinstance(cs["concepts"], list) and cs["concepts"]
                and all(isinstance(c, int) and not isinstance(c, bool) for c in cs["concepts"])
                and isinstance(cs["include_descendants"], bool)):
            problems.append(f"concept_sets.{name}: give {{concepts: [concept IDs, whole numbers], "
                            "include_descendants: true or false}")
    if m.get("patient_id") not in ("person_id", "person_source_value"):
        problems.append("patient_id: person_id or person_source_value")
    if m.get("diagnosis") not in sets:
        problems.append("diagnosis: must name a concept set")
    patients = m.get("patients")
    if not isinstance(patients, dict):
        problems.append("patients: must be a mapping")
        patients = {}
    for name in MAPPED:
        if name not in patients or patients[name] is None:
            problems.append(f"patients.{name}: not decided (give a source, or absent: true)")
        else:
            problems += _source_problems(f"patients.{name}", patients[name], sets, name)
    for name in set(patients) - set(MAPPED):
        problems.append(f"patients.{name}: not a patient field the mapping gives"
                        + (" (the adapter works it out)" if name in DERIVED else ""))
    extras = m.get("extras") or {}
    if not isinstance(extras, dict):
        problems.append("extras: must be a mapping")
        extras = {}
    standard = {f.name for t in TABLES for f in table_fields(t)}
    for name, spec in extras.items():
        if name in standard:
            problems.append(f"extras.{name}: a standard field's name")
        problems += _source_problems(f"extras.{name}", spec, sets, None)
    seq = m.get("sequencing")
    if not isinstance(seq, dict):
        problems.append("sequencing: must be a mapping")
        seq = {}
    if set(seq) - {"procedures", "modality", "entry_event", "panel_id", "variants"}:
        problems.append(f"sequencing: unknown keys {sorted(set(seq) - {'procedures', 'modality', 'entry_event', 'panel_id', 'variants'})}")
    if seq.get("procedures") not in sets:
        problems.append("sequencing.procedures: must name a concept set")
    if not (isinstance(seq.get("modality"), str) and seq["modality"].strip()):
        problems.append("sequencing.modality: give the tests' modality")
    if seq.get("entry_event") not in ENTRY_EVENTS:
        problems.append(f"sequencing.entry_event: one of {', '.join(ENTRY_EVENTS)}")
    panel = seq.get("panel_id")
    if not isinstance(panel, dict) or not (panel.get("absent") is True or isinstance(panel.get("procedure"), str)):
        problems.append("sequencing.panel_id: give {procedure: COLUMN} or {absent: true}")
    else:
        bad = set(panel) - ({"absent"} if panel.get("absent") is True else {"procedure", "extract", "labels"})
        if bad:
            problems.append(f"sequencing.panel_id: unknown keys {sorted(bad)}")
        if "extract" in panel:
            try:
                if re.compile(str(panel["extract"])).groups != 1:
                    problems.append("sequencing.panel_id: extract needs exactly one group in parentheses")
            except re.error as e:
                problems.append(f"sequencing.panel_id: extract is not a regular expression ({e})")
        problems += _labels_problems("sequencing.panel_id", panel.get("labels"))
    variants = seq.get("variants")
    if variants is not None:
        if not isinstance(variants, dict):
            problems.append("sequencing.variants: a mapping, or null if there is no variant table")
        else:
            for k, required in VARIANT_KEYS.items():
                v = variants.get(k)
                if k == "variant_class_labels":
                    problems += _labels_problems("sequencing.variants.variant_class_labels", v)
                elif (required and not isinstance(v, str)) or (not required and v is not None and not isinstance(v, str)):
                    problems.append(f"sequencing.variants.{k}: " + ("a column name" if k != "table" else "the table's name")
                                    + ("" if required else ", or null"))
            if set(variants) - set(VARIANT_KEYS):
                problems.append(f"sequencing.variants: unknown keys {sorted(set(variants) - set(VARIANT_KEYS))}")
            if variants.get("variant_class_labels") and not variants.get("variant_class"):
                problems.append("sequencing.variants.variant_class_labels: needs variant_class")
    readings = m.get("biomarker_readings")
    if not isinstance(readings, list):
        problems.append("biomarker_readings: a list (empty if there are none)")
        readings = []
    for i, r in enumerate(readings):
        where = f"biomarker_readings[{i}]"
        if not isinstance(r, dict):
            problems.append(f"{where}: must be a mapping")
            continue
        need = {"measurement", "kind", "gene", "value", "labels", "modality", "entry_event"}
        if set(r) != need:
            problems.append(f"{where}: give exactly {', '.join(sorted(need))}")
        if r.get("measurement") not in sets:
            problems.append(f"{where}.measurement: must name a concept set")
        if r.get("kind") not in READING_KINDS:
            problems.append(f"{where}.kind: one of {', '.join(READING_KINDS)}")
        if r.get("gene") is not None and not isinstance(r.get("gene"), str):
            problems.append(f"{where}.gene: a gene symbol, or null")
        if r.get("value") not in VALUE_COLUMNS:
            problems.append(f"{where}.value: one of {', '.join(VALUE_COLUMNS)}")
        if not (isinstance(r.get("modality"), str) and r["modality"].strip()):
            problems.append(f"{where}.modality: give the tests' modality")
        if r.get("entry_event") not in ENTRY_EVENTS:
            problems.append(f"{where}.entry_event: one of {', '.join(ENTRY_EVENTS)}")
        problems += _labels_problems(where, r.get("labels"))
    for key in ("panels_file", "panel_sizes_file"):
        if not _file_spec_ok(m.get(key)):
            problems.append(f"{key}: a path, {{path: ..., sheet: ...}}, or null")
    tr = m.get("treatments")
    if tr is not None and not (isinstance(tr, dict) and set(tr) == {"drugs", "name"} and tr["drugs"] in sets
                               and tr["name"] in DRUG_NAMES):
        problems.append(f"treatments: null, or {{drugs: CONCEPT_SET, name: {' or '.join(DRUG_NAMES)}}}")
    try:
        prof = profile_from_dict(m.get("profile"), adapter="omop")
    except ValueError as e:
        problems.append(str(e))
    else:
        if prof.time_format == "intervals":
            problems.append("profile.time_format: an OMOP database holds dates (real_dates or shifted_dates)")
        year = patients.get("year_of_diagnosis")
        if prof.time_format == "shifted_dates" and isinstance(year, dict) and year.get("computed"):
            problems.append("patients.year_of_diagnosis: shifted dates don't give the calendar year, so it can't be computed")
        if prof.entry_event != seq.get("entry_event"):
            problems.append("profile.entry_event: must say what sequencing.entry_event says")
        if (m.get("panels_file") is None) != (prof.panel_coverage == "none"):
            problems.append("profile.panel_coverage: none exactly when there is no panels_file")
        if (m.get("panel_sizes_file") is None) != (prof.panel_sizes == "not_given"):
            problems.append("profile.panel_sizes: not_given exactly when there is no panel_sizes_file")
    if problems:
        raise MappingError(f"{path.name}: " + "; ".join(problems))
    return m


# ---------------------------------------------------------------------------
# The CDM tables
# ---------------------------------------------------------------------------
class _Cdm:
    """The CDM's tables, each read once, with column names in lower case and dates as dates."""

    def __init__(self, spec: dict, base: Path, count):
        self.files: list[Path] = []
        self.notes: list[str] = []
        self._count = count
        self._cache: dict[str, pd.DataFrame | None] = {}
        if "folder" in spec:
            folder = Path(spec["folder"]) if Path(spec["folder"]).is_absolute() else base / spec["folder"]
            if not folder.is_dir():
                raise MappingError(f"cdm.folder: no folder {folder}")
            self.date_format = spec["date_format"]
            self.by_name: dict[str, Path] = {}
            for p in sorted(folder.iterdir()):
                if p.suffix.lower() in CDM_SUFFIXES:
                    if p.stem.lower() in self.by_name:
                        raise MappingError(f"cdm.folder: two files hold the table {p.stem.lower()}")
                    self.by_name[p.stem.lower()] = p
            self.con = None
        else:
            db = Path(spec["duckdb"]) if Path(spec["duckdb"]).is_absolute() else base / spec["duckdb"]
            if not db.is_file():
                raise MappingError(f"cdm.duckdb: no file {db}")
            self.con = duckdb.connect(str(db), read_only=True)
            self.catalog = self.con.execute("SELECT current_database()").fetchone()[0]
            self.schema = spec["schema"]
            self.files.append(db)
            self.date_format = None

    def close(self):
        if self.con is not None:
            self.con.close()

    def has(self, name: str) -> bool:
        if self.con is None:
            return name in self.by_name
        return bool(self.con.execute("SELECT count(*) FROM information_schema.tables WHERE table_catalog = ? "
                                     "AND table_schema = ? AND lower(table_name) = ?",
                                     [self.catalog, self.schema, name]).fetchone()[0])

    def table(self, name: str, columns: tuple[str, ...] = (), required: bool = True) -> pd.DataFrame | None:
        name = name.lower()
        if name not in self._cache:
            self._cache[name] = self._read(name) if self.has(name) else None
        df = self._cache[name]
        if df is None:
            if required:
                raise MappingError(f"the CDM has no {name} table")
            return None
        missing = [c for c in columns if c not in df.columns]
        if missing:
            raise MappingError(f"the {name} table has no column {', '.join(missing)}")
        return df

    def _read(self, name: str) -> pd.DataFrame:
        if self.con is None:
            path = self.by_name[name]
            self.files.append(path)
            if path.suffix.lower() == ".parquet":
                con = duckdb.connect()
                try:
                    df = con.execute("SELECT * FROM read_parquet(?)", [str(path)]).df()
                finally:
                    con.close()
            else:
                df = pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t", dtype=str,
                                 keep_default_na=False, encoding="utf-8")
                df = df.apply(lambda s: s.mask(s.str.strip() == "")) if len(df) else df
        else:
            df = self.con.execute(f'SELECT * FROM "{self.catalog}"."{self.schema}"."{name}"').df()
        df.columns = [str(c).lower() for c in df.columns]
        for col in DATE_COLUMNS.get(name, ()):
            if col in df.columns:
                df[col] = self._dates(name, col, df[col])
        return df

    def _dates(self, table: str, col: str, s: pd.Series) -> pd.Series:
        if pd.api.types.is_datetime64_any_dtype(s):
            return s.dt.normalize()
        filled = s.notna()
        if self.date_format is None:
            out = pd.to_datetime(s, errors="coerce")
        else:
            out = pd.to_datetime(s, format=self.date_format, errors="coerce")
        bad = int((filled & out.isna()).sum())
        if bad:
            self.notes.append(f"{table}.{col}: {self._count(bad)} dates didn't match the date format and were read "
                              "as blank")
        return out.dt.normalize()


def _concepts(cdm: _Cdm, sets: dict) -> dict[str, set[str]]:
    out = {}
    for name, cs in sets.items():
        ids = {str(c) for c in cs["concepts"]}
        if cs["include_descendants"]:
            anc = cdm.table("concept_ancestor", ("ancestor_concept_id", "descendant_concept_id"))
            hit = _keys(anc["ancestor_concept_id"]).isin(ids)
            ids |= set(_keys(anc.loc[hit.to_numpy(), "descendant_concept_id"]).dropna())
        out[name] = ids
    return out


# ---------------------------------------------------------------------------
# Patient fields
# ---------------------------------------------------------------------------
class _Patients:
    """Each patient's person row, index diagnosis and follow-up, indexed by person_id (as text)."""

    def __init__(self, cdm: _Cdm, m: dict, concepts: dict[str, set[str]], count):
        person = cdm.table("person", ("person_id", "year_of_birth", m["patient_id"]))
        cond = cdm.table("condition_occurrence", ("condition_occurrence_id", "person_id", "condition_concept_id",
                                                  "condition_start_date"))
        person = person.assign(_pid=_keys(person["person_id"]))
        if person["_pid"].duplicated().any():
            raise MappingError("person: a person_id appears more than once")
        cond = cond.assign(_pid=_keys(cond["person_id"]), _cid=pd.to_numeric(cond["condition_occurrence_id"], errors="coerce"))
        cond = cond[_keys(cond["condition_concept_id"]).isin(concepts[m["diagnosis"]]).to_numpy()
                    & cond["condition_start_date"].notna().to_numpy()]
        index = cond.sort_values(["_pid", "condition_start_date", "_cid"], kind="mergesort").drop_duplicates("_pid")
        self.no_diagnosis = int((~person["_pid"].isin(index["_pid"])).sum())
        self.not_in_person = int((~index["_pid"].isin(person["_pid"])).sum())
        self.person = person[person["_pid"].isin(index["_pid"])].set_index("_pid", drop=False)
        self.index = index.set_index("_pid").reindex(self.person.index)
        self.dx = self.index["condition_start_date"]
        ids = _keys(self.person[m["patient_id"]])
        if ids.isna().any() or ids.duplicated().any():
            raise MappingError(f"person.{m['patient_id']}: blank or repeated for a patient, so it can't be the patient ID")
        self.patient_id = ids
        op = cdm.table("observation_period", ("person_id", "observation_period_end_date"))
        last_seen = op.assign(_pid=_keys(op["person_id"])).groupby("_pid")["observation_period_end_date"].max()
        death = cdm.table("death", ("person_id", "death_date"), required=False)
        if death is not None:
            death = death.assign(_pid=_keys(death["person_id"]))
            died_on = death.groupby("_pid")["death_date"].min()
            died = self.person.index.isin(death["_pid"])
        else:
            died_on, died = pd.Series(dtype="datetime64[ns]"), [False] * len(self.person)
        died = pd.Series(died, index=self.person.index)
        seen = self.person.index.to_series().map(last_seen)
        self.end = self.person.index.to_series().map(died_on).where(died, seen)
        vital = pd.Series(pd.NA, index=self.person.index, dtype="Int64")
        vital[died] = 1
        vital[~died & seen.notna()] = 0
        self.vital = vital
        self.no_follow_up = int((~died & seen.isna()).sum())

    def days(self, when: pd.Series) -> pd.Series:
        return _days(when, self.dx.reindex(when.index))


def _labeled(values: pd.Series, spec: dict, where: str) -> tuple[pd.Series, int]:
    """Values after extract and labels, and how many filled values extract didn't match."""
    unmatched = 0
    if "extract" in spec:
        filled = values.notna()
        values = values.astype("string").str.extract(str(spec["extract"]), expand=False)
        unmatched = int((filled & values.isna()).sum())
    labels = spec.get("labels")
    if labels:
        table = {_key(k): v for k, v in labels.items()}
        keys = values.map(_key)
        unlisted = sorted({k for k in keys.dropna() if k not in table})
        if unlisted:
            shown = ", ".join(unlisted[:10]) + (", ..." if len(unlisted) > 10 else "")
            raise MappingError(f"{where}: values the labels don't list: {shown}. List each one (null reads it as blank).")
        values = keys.map(lambda k: pd.NA if pd.isna(k) else table[k])
        values = values.map(lambda v: pd.NA if v is None else v)
    return values, unmatched


def _event_value(cdm: _Cdm, table: str, ids: set[str], column: str, pts: _Patients) -> pd.Series:
    concept_col, date_col, id_col = EVENT_TABLES[table]
    df = cdm.table(table, ("person_id", concept_col, date_col, id_col, column))
    rows = df[_keys(df[concept_col]).isin(ids).to_numpy()].assign(_pid=lambda d: _keys(d["person_id"]))
    rows = rows[rows["_pid"].isin(pts.person.index)]
    gap = (rows[date_col] - rows["_pid"].map(pts.dx)).dt.days.abs()
    rows = rows.assign(_gap=gap.fillna(float("inf")), _id=pd.to_numeric(rows[id_col], errors="coerce"))
    first = rows.sort_values(["_pid", "_gap", date_col, "_id"], kind="mergesort").drop_duplicates("_pid")
    return pts.person.index.to_series().map(first.set_index("_pid")[column])


def _patient_field(cdm: _Cdm, name: str, spec: dict, pts: _Patients, concepts, where: str) -> tuple[pd.Series, int]:
    idx = pts.person.index
    if spec.get("absent") is True:
        return None, 0
    if spec.get("computed") is True:
        if name == "year_of_diagnosis":
            return pts.dx.dt.year.astype("Int64"), 0
        return (pts.dx.dt.year - pd.to_numeric(pts.person["year_of_birth"], errors="coerce")).astype("Float64"), 0
    if "person" in spec:
        col = spec["person"].lower()
        cdm.table("person", (col,))
        values = pts.person[col]
    elif "care_site" in spec:
        col = spec["care_site"].lower()
        cdm.table("person", ("care_site_id",))
        cs = cdm.table("care_site", ("care_site_id", col))
        by_id = cs.assign(_k=_keys(cs["care_site_id"])).drop_duplicates("_k").set_index("_k")[col]
        values = _keys(pts.person["care_site_id"]).map(by_id)
    elif "index_condition" in spec:
        col = spec["index_condition"].lower()
        cdm.table("condition_occurrence", (col,))
        values = pts.index[col]
    else:
        table = "measurement" if "measurement" in spec else "observation"
        values = _event_value(cdm, table, concepts[spec[table]], spec["value"].lower(), pts)
    values = pd.Series(values.to_numpy(), index=idx, dtype="object").map(lambda v: pd.NA if v is None else v)
    return _labeled(values, spec, where)


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------
def read(mapping_path: str | Path, settings: Settings) -> AdapterResult:
    """Read the OMOP database a confirmed mapping names into the standard tables. The settings give the
    privacy minimum used to mask the counts in what reading reports."""
    mapping_path = Path(mapping_path)
    m = load_mapping(mapping_path)
    base = mapping_path.parent
    prof = profile_from_dict(m["profile"], adapter="omop")
    min_count = max(settings.read("privacy.min_count"), prof.privacy_min_count or 0)

    def count(n: int) -> str:
        return show(n, min_count)

    cdm = _Cdm(m["cdm"], base, count)
    try:
        tables, notes, steps = _build(cdm, m, base, count)
    finally:
        cdm.close()
    files = [mapping_path] + cdm.files + [p for p in (_file_path(m.get(k), base) for k in ("panels_file", "panel_sizes_file"))
                                          if p is not None]
    return AdapterResult(tables=tables, profile=prof, input_files=files, source_steps=steps,
                         read_notes=cdm.notes + notes)


def _file_path(spec, base: Path) -> Path | None:
    if spec is None:
        return None
    p = Path(spec if isinstance(spec, str) else spec["path"])
    return p if p.is_absolute() else base / p


def _side_file(spec, base: Path, columns: tuple[str, ...], optional: tuple[str, ...], what: str) -> pd.DataFrame:
    df = read_flat_file(_file_path(spec, base), None if isinstance(spec, str) else spec.get("sheet"), ("", "NA"))
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise MappingError(f"{what}: no column {', '.join(missing)}")
    for c in optional:
        if c not in df.columns:
            df[c] = pd.NA
    return df[list(columns + optional)].reset_index(drop=True)


def _build(cdm: _Cdm, m: dict, base: Path, count) -> tuple[StandardTables, list[str], list[tuple[str, int]]]:
    notes: list[str] = []
    concepts = _concepts(cdm, m["concept_sets"])
    pts = _Patients(cdm, m, concepts, count)
    notes.append(f"{count(pts.no_diagnosis)} persons have no condition in the concept set '{m['diagnosis']}': left out "
                 "at the flow diagram's first step")
    if pts.not_in_person:
        notes.append(f"{count(pts.not_in_person)} persons with a diagnosis are not in the person table: left out")
    if pts.no_follow_up:
        notes.append(f"{count(pts.no_follow_up)} patients have neither a death nor an observation period: their outcome "
                     "is blank")
    pid = pts.patient_id
    patients = pd.DataFrame({"patient_id": pid.to_numpy(), "vital_status": pts.vital.to_numpy(),
                             "days_to_death_or_last_followup": pts.days(pts.end).to_numpy()})
    for name, spec in list((n, m["patients"][n]) for n in MAPPED) + list((m.get("extras") or {}).items()):
        where = f"patients.{name}" if name in MAPPED else f"extras.{name}"
        values, unmatched = _patient_field(cdm, name, spec, pts, concepts, where)
        if values is None:
            continue
        if unmatched:
            notes.append(f"{where}: {count(unmatched)} values didn't match extract and were read as blank")
        patients[name] = values.to_numpy()

    # ---- sequencing tests ----
    seq = m["sequencing"]
    proc = cdm.table("procedure_occurrence", ("procedure_occurrence_id", "person_id", "procedure_concept_id",
                                              "procedure_date"))
    proc = proc.assign(_pid=_keys(proc["person_id"]), _tid=_keys(proc["procedure_occurrence_id"]))
    proc = proc[_keys(proc["procedure_concept_id"]).isin(concepts[seq["procedures"]]).to_numpy()
                & proc["_pid"].isin(pts.person.index).to_numpy()]
    if proc["_tid"].duplicated().any():
        raise MappingError("procedure_occurrence: a procedure_occurrence_id appears more than once")
    proc = proc.set_index("_tid", drop=False)
    ngs = _tests(proc["_pid"], "procedure:" + proc["_tid"], proc["procedure_date"], seq["modality"],
                 seq["entry_event"], pts)
    panel = seq["panel_id"]
    if panel.get("absent") is not True:
        cdm.table("procedure_occurrence", (panel["procedure"].lower(),))
        values = pd.Series(proc[panel["procedure"].lower()].to_numpy(), index=proc.index, dtype="object")
        values, unmatched = _labeled(values.map(lambda v: pd.NA if v is None else v), panel, "sequencing.panel_id")
        if unmatched:
            notes.append(f"sequencing.panel_id: {count(unmatched)} values didn't match extract and were read as blank")
        ngs["panel_id"] = values.to_numpy()
    tests = [ngs]
    results = [empty_table("results")]
    variants = seq.get("variants")
    if variants:
        v_rows, unlinked = _variant_results(cdm, variants, proc)
        results.append(v_rows)
        notes.append(f"{count(unlinked)} variant rows are not on a patient's sequencing test: not used")

    # ---- readings, each a test with one result ----
    seen: set[str] = set()
    for i, r in enumerate(m["biomarker_readings"]):
        where = f"biomarker_readings[{i}]"
        meas = cdm.table("measurement", ("measurement_id", "person_id", "measurement_concept_id", "measurement_date",
                                         r["value"]))
        rows = meas.assign(_pid=_keys(meas["person_id"]), _mid=_keys(meas["measurement_id"]))
        rows = rows[_keys(rows["measurement_concept_id"]).isin(concepts[r["measurement"]]).to_numpy()
                    & rows["_pid"].isin(pts.person.index).to_numpy()]
        twice = seen & set(rows["_mid"])
        if twice or rows["_mid"].duplicated().any():
            raise MappingError(f"{where}: a measurement row is read twice (two readings' concept sets overlap, or a "
                               "measurement_id repeats)")
        seen |= set(rows["_mid"])
        test_ids = "measurement:" + rows["_mid"]
        tests.append(_tests(rows["_pid"], test_ids, rows["measurement_date"], r["modality"], r["entry_event"], pts))
        values = pd.Series(rows[r["value"]].to_numpy(), dtype="object").map(lambda v: pd.NA if v is None else v)
        values, _ = _labeled(values, {"labels": r["labels"]} if r["labels"] else {}, where)
        results.append(pd.DataFrame({"test_id": test_ids.to_numpy(), "kind": r["kind"], "gene": r["gene"],
                                     "variant": pd.NA, "variant_class": pd.NA, "value": values.to_numpy()}))
        undated = int(rows["measurement_date"].isna().sum())
        if undated:
            notes.append(f"{where}: {count(undated)} readings have no date: not known to come before the entry test")
    tests_df = pd.concat(tests, ignore_index=True)
    results_df = pd.concat([r for r in results if len(r)] or [empty_table("results")], ignore_index=True)
    undated_tests = int(ngs["days_to_entry_event"].isna().sum())
    if undated_tests:
        notes.append(f"{count(undated_tests)} sequencing tests have no date")

    # ---- panels, sizes and treatments ----
    panels = (_side_file(m["panels_file"], base, ("panel_id", "gene"), ("codon_from", "codon_to"), "panels_file")
              if m.get("panels_file") else empty_table("panels"))
    sizes = (_side_file(m["panel_sizes_file"], base, ("panel_id", "coding_megabases"), (), "panel_sizes_file")
             if m.get("panel_sizes_file") else None)
    treatments = None
    if m.get("treatments"):
        treatments, undated = _treatments(cdm, m["treatments"], concepts, pts)
        if undated:
            notes.append(f"{count(undated)} drug exposures in the drug concept set have no start date: left out")
    tables, failures = conform(StandardTables(patients=patients, tests=tests_df, results=results_df, panels=panels,
                                              treatments=treatments, panel_sizes=sizes))
    notes += [f"{k}: {count(n)} values didn't convert to the field's type and were read as blank"
              for k, n in failures.items()]
    steps = [(f"Has a condition in the concept set '{m['diagnosis']}'", pts.no_diagnosis)]
    return tables, notes, steps


def _tests(person: pd.Series, test_id: pd.Series, when: pd.Series, modality: str, entry_event: str,
           pts: _Patients) -> pd.DataFrame:
    dx = person.map(pts.dx)
    end = person.map(pts.end)
    after = pd.Series((when >= end).to_numpy(), dtype="boolean")
    after[(when.isna() | end.isna()).to_numpy()] = pd.NA
    return pd.DataFrame({"patient_id": person.map(pts.patient_id).to_numpy(), "test_id": test_id.to_numpy(),
                         "modality": modality, "entry_event": entry_event,
                         "days_to_entry_event": _days(when, dx).to_numpy(), "after_follow_up_end": after.to_numpy()})


def _variant_results(cdm: _Cdm, v: dict, proc: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    cols = tuple(v[k] for k in ("test_link", "gene", "protein_change", "partner_gene", "fusion_description",
                                "variant_class", "copy_number") if v.get(k))
    df = cdm.table(v["table"], tuple(c.lower() for c in cols))
    col = {k: (v[k].lower() if v.get(k) else None) for k in VARIANT_KEYS if k not in ("table", "variant_class_labels")}
    link = _keys(df[col["test_link"]])
    on_test = link.isin(proc.index)
    df, link = df[on_test.to_numpy()], link[on_test]
    test_id = ("procedure:" + link).to_numpy()
    gene = df[col["gene"]].astype("string").str.strip()
    partner = df[col["partner_gene"]].astype("string").str.strip() if col["partner_gene"] else pd.Series(pd.NA, index=df.index, dtype="string")
    cn = df[col["copy_number"]] if col["copy_number"] else pd.Series(pd.NA, index=df.index)
    is_cn = cn.notna().to_numpy()
    is_fusion = partner.notna().to_numpy() & ~is_cn
    is_mutation = ~is_cn & ~is_fusion
    vclass = pd.Series(pd.NA, index=df.index, dtype="object")
    if col["variant_class"]:
        vclass, _ = _labeled(pd.Series(df[col["variant_class"]].to_numpy(), index=df.index, dtype="object"),
                              {"labels": v["variant_class_labels"]} if v.get("variant_class_labels") else {},
                              "sequencing.variants.variant_class")
    out = [pd.DataFrame({"test_id": test_id[is_mutation], "kind": "mutation", "gene": gene[is_mutation].to_numpy(),
                         "variant": df.loc[is_mutation, col["protein_change"]].map(short_protein_change).to_numpy(),
                         "variant_class": vclass[is_mutation].to_numpy(), "value": pd.NA})]
    if is_fusion.any():
        desc = (df.loc[is_fusion, col["fusion_description"]].astype("string").to_numpy() if col["fusion_description"]
                else (gene[is_fusion] + "-" + partner[is_fusion]).to_numpy())
        for genes in (gene[is_fusion].to_numpy(), partner[is_fusion].to_numpy()):      # one row per partner gene
            out.append(pd.DataFrame({"test_id": test_id[is_fusion], "kind": "fusion", "gene": genes, "variant": desc,
                                     "variant_class": pd.NA, "value": pd.NA}))
    if is_cn.any():
        out.append(pd.DataFrame({"test_id": test_id[is_cn], "kind": "copy_number", "gene": gene[is_cn].to_numpy(),
                                 "variant": pd.NA, "variant_class": pd.NA,
                                 "value": cn[is_cn].map(_key).to_numpy()}))
    return pd.concat(out, ignore_index=True), int((~on_test).sum())


def _treatments(cdm: _Cdm, spec: dict, concepts, pts: _Patients) -> tuple[pd.DataFrame, int]:
    need = ("person_id", "drug_concept_id", "drug_exposure_start_date") + (("drug_source_value",)
                                                                          if spec["name"] == "drug_source_value" else ())
    de = cdm.table("drug_exposure", need)
    rows = de.assign(_pid=_keys(de["person_id"]))
    rows = rows[_keys(rows["drug_concept_id"]).isin(concepts[spec["drugs"]]).to_numpy()
                & rows["_pid"].isin(pts.person.index).to_numpy()]
    if spec["name"] == "concept_name":
        concept = cdm.table("concept", ("concept_id", "concept_name"))
        names = concept.assign(_k=_keys(concept["concept_id"])).drop_duplicates("_k").set_index("_k")["concept_name"]
        rows = rows.assign(_name=_keys(rows["drug_concept_id"]).map(names))
    else:
        rows = rows.assign(_name=rows["drug_source_value"])
    undated = int(rows["drug_exposure_start_date"].isna().sum())
    rows = rows[rows["drug_exposure_start_date"].notna()]
    rows = rows.assign(_name=rows["_name"].astype("string").str.strip())
    grouped = (rows.groupby(["_pid", "drug_exposure_start_date"])["_name"]
               .agg(lambda s: ", ".join(sorted(set(s.dropna()))) or pd.NA).reset_index())
    grouped = grouped.sort_values(["_pid", "drug_exposure_start_date"], kind="mergesort")
    grouped["regimen_order"] = grouped.groupby("_pid").cumcount() + 1
    dx = grouped["_pid"].map(pts.dx)
    out = pd.DataFrame({"patient_id": grouped["_pid"].map(pts.patient_id).to_numpy(),
                        "regimen": grouped["_name"].to_numpy(), "regimen_order": grouped["regimen_order"].to_numpy(),
                        "days_to_start": _days(grouped["drug_exposure_start_date"], dx).to_numpy()})
    return out, undated
