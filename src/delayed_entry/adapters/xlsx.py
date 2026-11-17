"""Read and write Excel workbooks (.xlsx) with the standard library only.

An .xlsx file is a zip of XML parts. read_xlsx reads the cached values of one sheet: text (shared or inline),
numbers, true/false, and dates (numbers with a date format, given as YYYY-MM-DD or YYYY-MM-DD HH:MM:SS).
Formulas are read as their last saved value. The first row is the header. write_xlsx writes sheets of text,
numbers and true/false, with the parts Excel needs to open the file.
"""
from __future__ import annotations

import re
import zipfile
from datetime import datetime, timedelta
from pathlib import Path
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape

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


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------
_CT = "application/vnd.openxmlformats-officedocument.spreadsheetml"
_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
DATE, WRAP, HEADER = 1, 2, 3          # cell styles: a date; wrapped text, top-aligned; bold wrapped text


def _col_name(i: int) -> str:
    out, i = "", i + 1
    while i:
        i, r = divmod(i - 1, 26)
        out = chr(65 + r) + out
    return out


def write_xlsx(path: str | Path, sheets: dict[str, list[list]], date_cells=frozenset(),
               widths: dict[str, list[float]] | None = None, wrapped: bool = False) -> Path:
    """Write sheets of rows ({sheet name: [row, ...]}). Text is stored as shared strings, numbers as numbers and
    true/false as booleans; a blank (None) is left empty. The cells named in date_cells ((sheet, row, column),
    counted from 0) get a date style; their value is an Excel day number. widths gives each sheet's column widths
    in characters. wrapped: every cell's text wraps, top-aligned, with the first row bold and kept in view."""
    shared: list[str] = []
    index: dict[str, int] = {}

    def sidx(s: str) -> int:
        if s not in index:
            index[s] = len(shared)
            shared.append(s)
        return index[s]
    parts = {}
    for n, (name, rows) in enumerate(sheets.items(), start=1):
        xml_rows = []
        for ri, row in enumerate(rows):
            cells = []
            for ci, v in enumerate(row):
                if v is None:
                    continue
                ref = f"{_col_name(ci)}{ri + 1}"
                style = DATE if (name, ri, ci) in date_cells else (HEADER if ri == 0 else WRAP) if wrapped else 0
                s = f' s="{style}"' if style else ""
                if isinstance(v, bool):
                    cells.append(f'<c r="{ref}" t="b"{s}><v>{int(v)}</v></c>')
                elif isinstance(v, (int, float)):
                    cells.append(f'<c r="{ref}"{s}><v>{v}</v></c>')
                else:
                    cells.append(f'<c r="{ref}" t="s"{s}><v>{sidx(str(v))}</v></c>')
            xml_rows.append(f'<row r="{ri + 1}">{"".join(cells)}</row>')
        view = ('<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" '
                'state="frozen"/></sheetView></sheetViews>') if wrapped else ""
        cols = "".join(f'<col min="{i + 1}" max="{i + 1}" width="{w}" customWidth="1"/>'
                       for i, w in enumerate((widths or {}).get(name, [])))
        parts[f"xl/worksheets/sheet{n}.xml"] = (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            f'{view}{"<cols>" + cols + "</cols>" if cols else ""}<sheetData>{"".join(xml_rows)}</sheetData></worksheet>')
    ns = f'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="{_REL}"'
    parts["xl/workbook.xml"] = (
        f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook {ns}><sheets>'
        + "".join(f'<sheet name="{escape(name)}" sheetId="{n}" r:id="rId{n}"/>' for n, name in enumerate(sheets, start=1))
        + "</sheets></workbook>")
    k = len(sheets)
    parts["xl/_rels/workbook.xml.rels"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        + "".join(f'<Relationship Id="rId{n}" Type="{_REL}/worksheet" Target="worksheets/sheet{n}.xml"/>'
                  for n in range(1, k + 1))
        + f'<Relationship Id="rId{k + 1}" Type="{_REL}/styles" Target="styles.xml"/>'
        + f'<Relationship Id="rId{k + 2}" Type="{_REL}/sharedStrings" Target="sharedStrings.xml"/>'
        + "</Relationships>")
    parts["xl/sharedStrings.xml"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        f'<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" count="{len(shared)}" '
        f'uniqueCount="{len(shared)}">' + "".join(f'<si><t xml:space="preserve">{escape(s)}</t></si>' for s in shared)
        + "</sst>")
    parts["xl/styles.xml"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
        '<font><b/><sz val="11"/><name val="Calibri"/></font></fonts>'
        '<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>'
        '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="4"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
        '<xf numFmtId="14" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1">'
        '<alignment wrapText="1" vertical="top"/></xf>'
        '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1" applyAlignment="1">'
        '<alignment wrapText="1" vertical="top"/></xf></cellXfs>'
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>')
    parts["_rels/.rels"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f'<Relationship Id="rId1" Type="{_REL}/officeDocument" Target="xl/workbook.xml"/></Relationships>')
    parts["[Content_Types].xml"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        f'<Override PartName="/xl/workbook.xml" ContentType="{_CT}.sheet.main+xml"/>'
        + "".join(f'<Override PartName="/xl/worksheets/sheet{n}.xml" ContentType="{_CT}.worksheet+xml"/>'
                  for n in range(1, k + 1))
        + f'<Override PartName="/xl/styles.xml" ContentType="{_CT}.styles+xml"/>'
        + f'<Override PartName="/xl/sharedStrings.xml" ContentType="{_CT}.sharedStrings+xml"/></Types>')
    path = Path(path)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for name in ("[Content_Types].xml", "_rels/.rels", "xl/workbook.xml", "xl/_rels/workbook.xml.rels",
                     "xl/styles.xml", "xl/sharedStrings.xml") + tuple(f"xl/worksheets/sheet{n}.xml" for n in range(1, k + 1)):
            z.writestr(name, parts[name])
    return path
