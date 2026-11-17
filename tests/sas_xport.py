"""A small writer of SAS transport files (XPORT version 5, as SAS's technical note TS-140 lays it out), so that the
tests can make made-up .xpt files: pandas reads them but can't write them.

    write_xport(path, {"PT_ID": ["a", "b"], "AGE": [61, None]}, name="PATIENTS")

A column whose values are all numbers (or blank) is numeric, stored as 8-byte IBM floating point, with a blank as
SAS's missing value; any other column is text, padded with spaces. Each record must be longer than 80 bytes, so that
the reader can't mistake a record's trailing spaces for the file's padding.
"""
from __future__ import annotations

import math
import struct
from pathlib import Path

STAMP = "01JAN26:00:00:00"
LIBRARY = "HEADER RECORD*******LIBRARY HEADER RECORD!!!!!!!000000000000000000000000000000  "
MEMBER = "HEADER RECORD*******MEMBER  HEADER RECORD!!!!!!!000000000000000001600000000140  "
DESCRIPTOR = "HEADER RECORD*******DSCRPTR HEADER RECORD!!!!!!!000000000000000000000000000000  "
OBSERVATIONS = "HEADER RECORD*******OBS     HEADER RECORD!!!!!!!000000000000000000000000000000  "


def _card(text: str) -> bytes:
    data = text.encode("ascii")
    assert len(data) <= 80, text
    return data.ljust(80, b" ")


def ibm_double(x) -> bytes:
    """A number as an 8-byte IBM floating-point value; blank (None or NaN) as SAS's missing value."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return b"\x2e" + b"\x00" * 7
    x = float(x)
    if x == 0:
        return b"\x00" * 8
    sign = 0x80 if x < 0 else 0
    x, exponent = abs(x), 0
    while x >= 1.0:
        x, exponent = x / 16.0, exponent + 1
    while x < 0.0625:
        x, exponent = x * 16.0, exponent - 1
    fraction = int(round(x * (1 << 56)))
    if fraction >= 1 << 56:
        fraction, exponent = fraction >> 4, exponent + 1
    return bytes([sign | (exponent + 64)]) + fraction.to_bytes(7, "big")


def _numeric(values) -> bool:
    return all(v is None or (isinstance(v, (int, float)) and not isinstance(v, bool)) for v in values)


def write_xport(path: Path, columns: dict[str, list], name: str = "DATA", encoding: str = "latin-1") -> Path:
    names = list(columns)
    n = len(next(iter(columns.values())))
    assert all(len(v) == n for v in columns.values()) and all(len(c) <= 8 for c in names)
    fields, pos = [], 0
    for i, c in enumerate(names):
        values = columns[c]
        if _numeric(values):
            fields.append((c, 1, 8, pos))
        else:
            width = max([len(str(v).encode(encoding)) for v in values if v is not None] + [1])
            fields.append((c, 2, width, pos))
        pos += fields[-1][2]
    assert pos > 80, "each record must be longer than 80 bytes"
    out = [_card(LIBRARY),
           _card("SAS     SAS     SASLIB  9.4     X64_10PR" + " " * 24 + STAMP),
           _card(STAMP),
           _card(MEMBER), _card(DESCRIPTOR),
           _card("SAS     " + name.upper().ljust(8)[:8] + "SASDATA 9.4     X64_10PR" + " " * 24 + STAMP),
           _card(STAMP + " " * 16 + "made-up data".ljust(40) + "DATA    "),
           _card(f"HEADER RECORD*******NAMESTR HEADER RECORD!!!!!!!000000{len(fields):04d}00000000000000000000  ")]
    namestrs = b"".join(
        struct.pack(">hhhh8s40s8shhh2s8shhl52s", kind, 0, width, i + 1, c.encode("ascii").ljust(8),
                    b"".ljust(40), b"".ljust(8), 0, 0, 0, b"  ", b"".ljust(8), 0, 0, start, b"\x00" * 52)
        for i, (c, kind, width, start) in enumerate(fields))
    out.append(namestrs + b" " * (-len(namestrs) % 80))
    out.append(_card(OBSERVATIONS))
    rows = b""
    for r in range(n):
        for c, kind, width, _ in fields:
            v = columns[c][r]
            rows += ibm_double(v) if kind == 1 else ("" if v is None else str(v)).encode(encoding).ljust(width, b" ")
    out.append(rows + b" " * (-len(rows) % 80))
    path = Path(path)
    path.write_bytes(b"".join(out))
    return path
