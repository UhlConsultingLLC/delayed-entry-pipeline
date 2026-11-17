"""Load the GENIE BPC files into DuckDB tables, exactly as released, and check the links between them.

Copied from the CRC BPC project (crcbpc.db), with the list of files taken from the cohort's own spec
(cohorts/<cohort>/raw_files.yml) instead of one fixed list.

Each file becomes one raw_ table: no rows dropped, no values changed. Every column the adapter uses
gets a fixed type, so DuckDB never has to guess it. Other columns are read with DuckDB's own guesses.

    from delayed_entry.adapters.genie_bpc import raw
    con, files = raw.load("crc_v2.0", clinical_dir, genomic_dir)
    print(raw.run_checks(con))
"""
from __future__ import annotations

import re
from pathlib import Path

import duckdb
import yaml

HERE = Path(__file__).parent
COHORTS = HERE / "cohorts"
CHECKS_SQL = HERE / "link_checks.sql"
PANEL_FILES = "data_gene_panel_*.txt"
# A cohort whose spec reads no fusion file gets this empty table, so the steps and checks that read fusions run
# the same way and find none. It does nothing when the fusion file was read.
EMPTY_FUSION = "CREATE TABLE IF NOT EXISTS raw_fusion (Hugo_Symbol VARCHAR, Tumor_Sample_Barcode VARCHAR, Fusion VARCHAR)"


def cohort_dir(cohort: str) -> Path:
    d = COHORTS / cohort
    if not d.is_dir():
        known = ", ".join(sorted(p.name for p in COHORTS.iterdir() if p.is_dir()))
        raise ValueError(f"no GENIE BPC cohort '{cohort}' (known: {known})")
    return d


def raw_files(cohort: str) -> dict[str, tuple[str, str, dict[str, str]]]:
    """{table: (file name, which folder, {column: type})}."""
    spec = yaml.safe_load((cohort_dir(cohort) / "raw_files.yml").read_text(encoding="utf-8"))
    return {table: (s["file"], s["folder"], s["types"]) for table, s in spec.items()}


def sql_text(s: str) -> str:
    """A string written as SQL text: wrapped in single quotes, with any quote inside doubled."""
    return "'" + s.replace("'", "''") + "'"


def read_table_sql(table: str, path: Path, types: dict[str, str]) -> str:
    """The SQL that makes one raw table from one file. A blank or the text NA is a missing value, as in
    genieBPC's reader. sample_size = -1 makes DuckDB look at every row before guessing a column's type,
    so a column that turns to text deep in a file can't stop the load."""
    delim = "\t" if path.suffix == ".txt" else ","
    type_list = ", ".join(f"{sql_text(col)}: {sql_text(t)}" for col, t in types.items())
    return (f"CREATE OR REPLACE TABLE {table} AS\n"
            f"SELECT * FROM read_csv({sql_text(path.as_posix())}, header = true, delim = {sql_text(delim)},\n"
            f"                       nullstr = ['', 'NA'], sample_size = -1,\n"
            f"                       types = {{{type_list}}})")


def read_panel_file(path: Path) -> tuple[str, list[str]]:
    """A gene-panel file ("stable_id: X", "gene_list: A B C ...", tab-separated genes) as (panel name, genes)."""
    fields = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        key, _, value = line.partition(":")
        fields[key.strip()] = value
    genes = [g.strip() for g in fields["gene_list"].split("\t") if g.strip()]
    return fields["stable_id"].strip(), genes


def load(cohort: str, clinical_dir: str | Path, genomic_dir: str | Path | None = None,
         database: str | Path = ":memory:") -> tuple[duckdb.DuckDBPyConnection, list[Path]]:
    """Make one raw_ table per file, plus raw_panel_gene (one row per panel and gene). Returns the
    connection and every file read. For made-up data both kinds of file may sit in one folder. A cohort
    that reads no fusion file gets an empty raw_fusion table."""
    folders = {"clinical": Path(clinical_dir), "genomic": Path(genomic_dir or clinical_dir)}
    con = duckdb.connect(str(database))
    read = []
    for table, (name, folder, types) in raw_files(cohort).items():
        path = folders[folder] / name
        if not path.exists():
            raise FileNotFoundError(f"{table}: no file {name} in the {folder} folder")
        con.execute(read_table_sql(table, path, types))
        read.append(path)
    con.execute(EMPTY_FUSION)
    panels = sorted(folders["genomic"].glob(PANEL_FILES))
    if not panels:
        raise FileNotFoundError(f"no gene-panel files ({PANEL_FILES}) in the genomic folder")
    con.execute("CREATE OR REPLACE TABLE raw_panel_gene (panel VARCHAR, gene VARCHAR)")
    for path in panels:
        panel, genes = read_panel_file(path)
        con.executemany("INSERT INTO raw_panel_gene VALUES (?, ?)", [(panel, g) for g in genes])
        read.append(path)
    return con, read


def read_checks(sql_file: Path = CHECKS_SQL) -> dict[str, str]:
    """Split a checks file into {check name: query}. Each check starts with a line '-- check: name'."""
    parts = re.split(r"^-- check: *(\S+) *$", sql_file.read_text(encoding="utf-8"), flags=re.M)
    return {name: query.strip().rstrip(";") for name, query in zip(parts[1::2], parts[2::2])}


def scalar(con: duckdb.DuckDBPyConnection, query: str, params=None):
    """The single value a query returns. fetchall() reads the result to the end, which DuckDB requires
    before the connection can do anything else."""
    return con.execute(query, params or []).fetchall()[0][0]


def run_checks(con: duckdb.DuckDBPyConnection, sql_file: Path = CHECKS_SQL) -> dict[str, int]:
    """Run every link check and count the rows each one returns. 0 means the check passes."""
    return {name: scalar(con, f"SELECT count(*) FROM ({query})") for name, query in read_checks(sql_file).items()}
