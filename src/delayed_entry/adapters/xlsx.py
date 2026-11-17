"""Read one sheet of an Excel workbook (.xlsx) with the standard library only.

An .xlsx file is a zip of XML parts. This reads the cached values of one sheet: text (shared or inline),
numbers, true/false, and dates (numbers with a date format, given as YYYY-MM-DD or YYYY-MM-DD HH:MM:SS).
Formulas are read as their last saved value. The first row is the header.
"""
from __future__ import annotations

import re
import zipfile
from datetime import datetime, timedelta
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd

M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
DATE_FORMAT_IDS = set(range(14, 23)) | {45, 46, 47}


def _col_index(ref: str) -> int:
    letters = re.match(r"[A-Z]+", ref).group(0)
    n = 0
    for ch in letters:
        n = n * 26 + ord(ch) - 64
    return n - 1


def _text(el) -> str:
    return "".join(t.text or "" for t in el.iter(f"{M}t"))


def _is_date_format(code: str) -> bool:
    code = re.sub(r'"[^"]*"|\[[^\]]*\]|\\.', "", code)     # drop quoted text, colors and escapes
    return bool(re.search(r"[dmyhs]", code, re.I))


def sheet_names(path: str | Path) -> list[str]:
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
    return [s.get("name") for s in wb.find(f"{M}sheets")]


def read_xlsx(path: str | Path, sheet: str | None = None) -> pd.DataFrame:
    """The sheet (the first if none is named) as a table of text, numbers, booleans and blanks."""
    with zipfile.ZipFile(path) as z:
        names = set(z.namelist())
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        pr = wb.find(f"{M}workbookPr")
        epoch = datetime(1904, 1, 1) if pr is not None and pr.get("date1904") in ("1", "true") else datetime(1899, 12, 30)
        sheets = [(s.get("name"), s.get(f"{R}id")) for s in wb.find(f"{M}sheets")]
        if sheet is None:
            name, rid = sheets[0]
        else:
            match = [x for x in sheets if x[0] == sheet]
            if not match:
                raise ValueError(f"no sheet '{sheet}' (sheets: {', '.join(n for n, _ in sheets)})")
            name, rid = match[0]
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = next(r.get("Target") for r in rels if r.get("Id") == rid)
        part = target.lstrip("/") if target.startswith("/") else "xl/" + target
        shared = []
        if "xl/sharedStrings.xml" in names:
            shared = [_text(si) for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(f"{M}si")]
        date_styles = set()
        if "xl/styles.xml" in names:
            st = ET.fromstring(z.read("xl/styles.xml"))
            custom = {int(f.get("numFmtId")): f.get("formatCode", "") for f in st.iter(f"{M}numFmt")}
            xfs = st.find(f"{M}cellXfs")
            for i, xf in enumerate(xfs if xfs is not None else []):
                fid = int(xf.get("numFmtId", "0"))
                if fid in DATE_FORMAT_IDS or (fid in custom and _is_date_format(custom[fid])):
                    date_styles.add(i)
        root = ET.fromstring(z.read(part))
    rows: dict[int, dict[int, object]] = {}
    for r_i, row in enumerate(root.iter(f"{M}row")):
        ri = int(row.get("r", r_i + 1)) - 1
        for c in row.findall(f"{M}c"):
            ci = _col_index(c.get("r")) if c.get("r") else len(rows.get(ri, {}))
            kind, v = c.get("t", "n"), c.find(f"{M}v")
            if kind == "s":
                value = shared[int(v.text)] if v is not None else None
            elif kind == "inlineStr":
                is_ = c.find(f"{M}is")
                value = _text(is_) if is_ is not None else None
            elif kind == "b":
                value = None if v is None else v.text == "1"
            elif kind in ("str", "e"):
                value = None if v is None else v.text
            else:
                value = None if v is None or v.text is None else float(v.text)
                if value is not None and int(c.get("s", "0")) in date_styles:
                    d = epoch + timedelta(days=value)
                    value = d.strftime("%Y-%m-%d") if d.time() == datetime.min.time() else d.strftime("%Y-%m-%d %H:%M:%S")
                elif value is not None and value.is_integer():
                    value = int(value)
            rows.setdefault(ri, {})[ci] = value
    if not rows:
        return pd.DataFrame()
    first = min(rows)
    header = rows[first]
    width = max(header) + 1
    cols = [str(header.get(i)) if header.get(i) is not None else f"column_{i + 1}" for i in range(width)]
    body = [[rows[r].get(i) for i in range(width)] for r in sorted(rows) if r != first]
    return pd.DataFrame(body, columns=cols, dtype="object")
