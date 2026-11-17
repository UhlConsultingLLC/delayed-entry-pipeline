"""Builders for the tests: made-up data only, never real counts."""
from __future__ import annotations

import dataclasses
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

import duckdb
import numpy as np
import pandas as pd
import yaml

from delayed_entry.adapters import genie_bpc
from delayed_entry.profile import CapabilityProfile
from delayed_entry.settings import Settings, load_settings
from delayed_entry.tables import StandardTables, conform

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "genie_bpc_crc_synthetic"
CRC_SETTINGS = ROOT / "settings" / "crc_bpc.yaml"


def crc_settings() -> Settings:
    return load_settings(CRC_SETTINGS)


def settings_with(tmp_path: Path, changes: dict, name: str = "changed.yaml") -> Settings:
    """The CRC settings file with some values changed, written to a new file and loaded from it, as a user
    would: no code changes."""
    data = yaml.safe_load(CRC_SETTINGS.read_text(encoding="utf-8"))
    for path, value in changes.items():
        node = data
        for part in path.split("."):
            node = node[part]
        node["value"] = value
    p = tmp_path / name
    p.write_text(yaml.safe_dump(data, sort_keys=False, allow_unicode=True), encoding="utf-8")
    return load_settings(p)


def bpc_profile(**changes) -> CapabilityProfile:
    return dataclasses.replace(genie_bpc.profile(), **changes)


def scrambled(t: StandardTables, seed: int) -> StandardTables:
    """Every recorded outcome value replaced by a random one; blanks stay blank."""
    rng = np.random.default_rng(seed)
    p = t.patients.copy()
    known = p["vital_status"].notna()
    p.loc[known, "vital_status"] = rng.integers(0, 2, int(known.sum()))
    known = p["days_to_death_or_last_followup"].notna()
    p.loc[known, "days_to_death_or_last_followup"] = rng.uniform(1, 5000, int(known.sum()))
    return dataclasses.replace(t, patients=p)


# ---- Standard tables by hand ----
RAS_BRAF = ["KRAS", "NRAS", "BRAF"]


def make_tables(groups: dict[str, int], *, extra_patients: list[dict] | None = None, msi: dict | None = None,
                treated: dict[str, int] | None = None, report_day: float = 100.0) -> StandardTables:
    """Standard tables with `groups` patients per CRC genotype group: each stage IV, age 60, one NGS report on
    `report_day` on a panel reading KRAS, NRAS and BRAF, and one call that puts them in their group.
    msi: {group: (n MSI-H, n MSS, n unknown)} gives that many patients of the group a pathology result dated
    before the report. treated: {group: n} gives n patients of the group a first-line regimen."""
    call = {"RAS-mutant": ("KRAS", "p.G12D", "Missense_Mutation"), "BRAF V600E": ("BRAF", "p.V600E", "Missense_Mutation"),
            "Wild type": ("TP53", "p.R175H", "Missense_Mutation")}
    patients, tests, results, treatments = [], [], [], []
    for g, n in groups.items():
        hi, ms, unk = (msi or {}).get(g, (0, 0, n))
        readings = ["MSI-H"] * hi + ["MSS"] * ms + [None] * (n - hi - ms)
        for i in range(n):
            pid = f"{g[:3]}-{i:04d}"
            patients.append({"patient_id": pid, "vital_status": 1, "days_to_death_or_last_followup": 900.0,
                             "age_at_diagnosis": 60.0, "sex": "Female", "stage_at_diagnosis": "Stage IV",
                             "histology": "Adenocarcinoma", "primary_site_code": "C18.7", "year_of_diagnosis": 2015,
                             "center": "MSK", "censored_at_age_limit": False, "multiple_cohort_cancers": False,
                             "other_primary_same_cancer": False, "tumor_location": "Left colon"})
            tests.append({"patient_id": pid, "test_id": f"{pid}-T01", "modality": "ngs", "entry_event": "report",
                          "days_to_entry_event": report_day, "panel_id": "P8", "test_order": 1,
                          "after_follow_up_end": False, "linked_to_other_cancers": False})
            gene, variant, cls = call[g]
            results.append({"test_id": f"{pid}-T01", "kind": "mutation", "gene": gene, "variant": variant,
                            "variant_class": cls, "value": None})
            if readings[i] is not None:
                tests.append({"patient_id": pid, "test_id": f"{pid}-P01", "modality": "pathology", "entry_event": "report",
                              "days_to_entry_event": report_day - 10, "panel_id": None, "test_order": None,
                              "after_follow_up_end": None, "linked_to_other_cancers": None})
                results.append({"test_id": f"{pid}-P01", "kind": "msi_mmr", "gene": None, "variant": None,
                                "variant_class": None, "value": readings[i]})
            if i < (treated or {}).get(g, 0):
                treatments.append({"patient_id": pid, "regimen": None, "regimen_order": 1, "line": None,
                                   "days_to_start": 30.0})
    patients += extra_patients or []
    panels = [{"panel_id": "P8", "gene": g, "codon_from": None, "codon_to": None}
              for g in RAS_BRAF + ["TP53", "APC", "GENE6", "GENE7", "GENE8"]]
    t = StandardTables(pd.DataFrame(patients), pd.DataFrame(tests), pd.DataFrame(results), pd.DataFrame(panels),
                       pd.DataFrame(treatments, columns=["patient_id", "regimen", "regimen_order", "line", "days_to_start"]),
                       None)
    return conform(t)[0]


# ---- Raw GENIE BPC tables by hand (ported from the CRC BPC project's tests) ----
DOB_DX = 60 * 365.25        # a hand-made report's days from birth assume diagnosis at age 60
HAND_MADE_PANELS = {        # made-up panels: name, genes read
    "P8": RAS_BRAF + [f"GENE{i}" for i in range(4, 9)],      # 8 genes
    "P10": RAS_BRAF + [f"GENE{i}" for i in range(4, 11)],    # 10 genes
    "NO_NRAS": ["KRAS", "BRAF"],
    "NO_BRAF": ["KRAS", "NRAS"],
}
PATIENT_DEFAULTS = {"institution": "MSK", "birth_year": 1955, "naaccr_sex_code": "Male",
                    "race_ethnicity": "Non-Hispanic White"}


def cancer(record_id, ca_seq=0, index="Yes", stage="Stage IV", histology="Adenocarcinoma", age=60.0, site="C18.7",
           os_days=1000.0, os_status=1, crc_type="Left colon"):
    """One cancer. Unless given, it passes every rule, and the patient died 1,000 days after diagnosis."""
    dob_dx = None if age is None else age * 365.25
    return (record_id, ca_seq, index, stage, histology, age, site, os_days, os_status, dob_dx, crc_type)


def report(record_id, cpt_number, day, panel="P8", ca_seq=0, n_ca_seq=1, post_death=0, post_last_alive=0,
           sample_type="Primary tumor", seq_year="2015"):
    """One sequencing report. Its sample ID is the patient's ID plus the test number."""
    dob_report = None if day is None else DOB_DX + day
    months = None if day is None else day / (365.25 / 12)
    return (record_id, ca_seq, cpt_number, n_ca_seq, day, panel, f"{record_id}-T{cpt_number:02d}",
            post_death, post_last_alive, dob_report, months, sample_type, seq_year)


def path_report(record_id, n, specimens=("Colon Cancer",), msi=(), mmr=()):
    """One pathology report. msi and mmr: up to three (label, days from birth) each."""
    spec = list(specimens) + [None] * (3 - len(specimens))
    msi = list(msi) + [(None, None)] * (3 - len(msi))
    mmr = list(mmr) + [(None, None)] * (3 - len(mmr))
    return (record_id, 1, n, *spec, *(x for pair in msi for x in pair), *(x for pair in mmr for x in pair))


def hand_made_db(cases, redacted=frozenset(), calls=(), patients=None, pathology=(), others=(), regimens=()):
    """A tiny database of raw BPC tables holding only the columns the adapter reads.

    cases: {record_id: (cancer rows, report rows)}. redacted: record_ids whose redacted flag is "Yes".
    calls: (sample ID, gene, variant class, protein change). pathology: path_report rows.
    others: (record_id, ca_seq, redcap_ca_index, ca_type). regimens: (record_id, ca_seq, number, start day)."""
    c = duckdb.connect()
    c.execute("""CREATE TABLE raw_patient (record_id VARCHAR, redacted VARCHAR, institution VARCHAR, birth_year INTEGER,
                 naaccr_sex_code VARCHAR, race_ethnicity VARCHAR)""")
    c.executemany("INSERT INTO raw_patient VALUES (?, ?, ?, ?, ?, ?)",
                  [(rid, "Yes" if rid in redacted else "No", *({**PATIENT_DEFAULTS, **(patients or {}).get(rid, {})}).values())
                   for rid in cases])
    c.execute("""CREATE TABLE raw_cancer (record_id VARCHAR, ca_seq INTEGER, redcap_ca_index VARCHAR, stage_dx VARCHAR,
                 ca_hist_adeno_squamous VARCHAR, age_dx DOUBLE, ca_d_site VARCHAR, tt_os_dx_days DOUBLE,
                 os_dx_status INTEGER, dob_ca_dx_days DOUBLE, crc_type VARCHAR)""")
    c.executemany("INSERT INTO raw_cancer VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                  [r for cancers, _ in cases.values() for r in cancers])
    c.execute("""CREATE TABLE raw_report (record_id VARCHAR, ca_seq INTEGER, cpt_number INTEGER, cpt_n_ca_seq INTEGER,
                 dx_cpt_rep_days DOUBLE, cpt_seq_assay_id VARCHAR, cpt_genie_sample_id VARCHAR,
                 cpt_report_post_death INTEGER, cpt_report_post_last_alive INTEGER, dob_cpt_report_days DOUBLE,
                 dx_cpt_rep_mos DOUBLE, sample_type VARCHAR, cpt_seq_date VARCHAR)""")
    c.executemany("INSERT INTO raw_report VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                  [r for _, reports in cases.values() for r in reports])
    c.execute("CREATE TABLE raw_panel_gene (panel VARCHAR, gene VARCHAR)")
    c.executemany("INSERT INTO raw_panel_gene VALUES (?, ?)",
                  [(panel, gene) for panel, genes in HAND_MADE_PANELS.items() for gene in genes])
    c.execute("""CREATE TABLE raw_mutation (Tumor_Sample_Barcode VARCHAR, Hugo_Symbol VARCHAR,
                 Variant_Classification VARCHAR, HGVSp_Short VARCHAR)""")
    if calls:
        c.executemany("INSERT INTO raw_mutation VALUES (?, ?, ?, ?)", list(calls))
    c.execute("CREATE TABLE raw_cancer_other (record_id VARCHAR, ca_seq INTEGER, redcap_ca_index VARCHAR, ca_type VARCHAR)")
    if others:
        c.executemany("INSERT INTO raw_cancer_other VALUES (?, ?, ?, ?)", list(others))
    c.execute("""CREATE TABLE raw_pathology (record_id VARCHAR, path_proc_number INTEGER, path_rep_number INTEGER,
                 path_ca_type1 VARCHAR, path_ca_type2 VARCHAR, path_ca_type3 VARCHAR,
                 msi_result VARCHAR, msi_prepaint DOUBLE, msi_result_2 VARCHAR, msi_prepaint_2 DOUBLE,
                 msi_result_3 VARCHAR, msi_prepaint_3 DOUBLE,
                 mmr_result VARCHAR, mmr_prepaint DOUBLE, mmr_result_2 VARCHAR, mmr_prepaint_2 DOUBLE,
                 mmr_result_3 VARCHAR, mmr_prepaint_3 DOUBLE)""")
    if pathology:
        c.executemany(f"INSERT INTO raw_pathology VALUES ({', '.join(['?'] * 18)})", list(pathology))
    c.execute("""CREATE TABLE raw_regimen (record_id VARCHAR, ca_seq INTEGER, regimen_number_within_cancer INTEGER,
                 redcap_ca_index VARCHAR, dx_reg_start_int DOUBLE, os_g_status INTEGER, tt_os_g_days DOUBLE)""")
    if regimens:
        c.executemany("INSERT INTO raw_regimen VALUES (?, ?, ?, 'Yes', ?, 0, 100)", list(regimens))
    return c


def hand_made_tables(settings: Settings, **kw) -> StandardTables:
    con = hand_made_db(**kw)
    return genie_bpc.tables_from_connection(con, settings)[0]


# ---- A minimal .xlsx writer, for the reader's tests ----
def write_xlsx(path: Path, sheets: dict[str, list[list]], date_cells: set[tuple[str, int, int]] = frozenset()) -> Path:
    """Write sheets of rows. Text is stored as shared strings, numbers as numbers, booleans as booleans, and
    the cells named in date_cells (sheet, row, column) get a date style (their value is an Excel day number)."""
    shared: list[str] = []
    def sidx(s):
        if s not in shared:
            shared.append(s)
        return shared.index(s)
    def col(i):
        out = ""
        i += 1
        while i:
            i, r = divmod(i - 1, 26)
            out = chr(65 + r) + out
        return out
    parts = {}
    for n, (name, rows) in enumerate(sheets.items(), start=1):
        xml_rows = []
        for ri, row in enumerate(rows):
            cells = []
            for ci, v in enumerate(row):
                ref = f"{col(ci)}{ri + 1}"
                style = ' s="1"' if (name, ri, ci) in date_cells else ""
                if v is None:
                    continue
                if isinstance(v, bool):
                    cells.append(f'<c r="{ref}" t="b"><v>{int(v)}</v></c>')
                elif isinstance(v, (int, float)):
                    cells.append(f'<c r="{ref}"{style}><v>{v}</v></c>')
                else:
                    cells.append(f'<c r="{ref}" t="s"><v>{sidx(str(v))}</v></c>')
            xml_rows.append(f'<row r="{ri + 1}">{"".join(cells)}</row>')
        parts[f"xl/worksheets/sheet{n}.xml"] = (
            '<?xml version="1.0" encoding="UTF-8"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            f'<sheetData>{"".join(xml_rows)}</sheetData></worksheet>')
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    parts["xl/workbook.xml"] = (f'<?xml version="1.0" encoding="UTF-8"?><workbook {ns}><sheets>'
                                + "".join(f'<sheet name="{escape(name)}" sheetId="{n}" r:id="rId{n}"/>'
                                          for n, name in enumerate(sheets, start=1)) + "</sheets></workbook>")
    parts["xl/_rels/workbook.xml.rels"] = (
        '<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        + "".join(f'<Relationship Id="rId{n}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
                  f'Target="worksheets/sheet{n}.xml"/>' for n in range(1, len(sheets) + 1)) + "</Relationships>")
    parts["xl/sharedStrings.xml"] = (
        '<?xml version="1.0" encoding="UTF-8"?><sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        + "".join(f"<si><t>{escape(s)}</t></si>" for s in shared) + "</sst>")
    parts["xl/styles.xml"] = (
        '<?xml version="1.0" encoding="UTF-8"?><styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<cellXfs count="2"><xf numFmtId="0"/><xf numFmtId="14"/></cellXfs></styleSheet>')
    parts["[Content_Types].xml"] = '<?xml version="1.0" encoding="UTF-8"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"/>'
    with zipfile.ZipFile(path, "w") as z:
        for name, xml in parts.items():
            z.writestr(name, xml)
    return path
