"""Protein changes in HGVSp short form, as the standard tables hold them (p.G12D), from whatever form a source
gives: the long three-letter form, with or without a reference sequence or parentheses, or as VEP writes them.

    short_protein_change("NP_004976.2:p.(Gly12Asp)")  ->  "p.G12D"
"""
from __future__ import annotations

import re

import pandas as pd

AMINO_ACIDS = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
               "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
               "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Sec": "U", "Pyl": "O", "Ter": "*", "Xaa": "X"}
_THREE_LETTER = re.compile("|".join(AMINO_ACIDS))


def short_protein_change(v):
    """A protein change in HGVSp short form: p.Gly12Asp -> p.G12D, NP_004976.2:p.(Arg1450Ter) -> p.R1450*, and VEP's
    p.Leu10%3D (an escaped =) -> p.L10=. A change already in short form, or anything that isn't a protein change
    (no p.), is left as it is; a blank stays blank."""
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return pd.NA
    s = str(v).strip().split(":")[-1].replace("%3D", "=")
    if not s:
        return pd.NA
    if s.startswith("p.(") and s.endswith(")"):
        s = "p." + s[3:-1]
    return _THREE_LETTER.sub(lambda m: AMINO_ACIDS[m.group(0)], s) if s.startswith("p.") else s
