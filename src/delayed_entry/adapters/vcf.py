"""A VCF's annotated variants as a table the flat-file mapping can use.

A VCF holds no protein change unless an annotation tool has added one, so the reader takes the annotation the file
carries in its INFO column: VEP's CSQ or SnpEff's ANN, whose subfields the file's header names. The table has one
row per variant record, per sample that carries it, per annotation kept:

    SAMPLE                     the sample's column name; for a file with no sample columns, the file's name
                               without its extension
    CHROM, POS, ID, REF, ALT, QUAL, FILTER   as the record gives them (ALT: the allele the sample carries)
    GT                         the sample's genotype (blank when the file has no sample columns)
    each annotation subfield   under the header's name for it: for VEP, Allele, Consequence, SYMBOL, HGVSp,
                               CANONICAL...; for SnpEff, Allele, Annotation, Gene_Name, HGVS.p...

A sample carries an allele when its genotype (GT) names it; a file with no GT counts every allele as carried. An
annotation belongs to the allele it names: the ALT itself, or VEP's form of it with a shared first base removed
("-" for a deletion of everything after it). A record with one ALT keeps its annotations even when no allele
matches. Two reading rules come from the mapping, never from code:
    pick      canonical (VEP's CANONICAL=YES), first (the allele's first annotation listed) or all
    filters   pass (records whose FILTER is PASS or .) or all
Every value is text. A carried variant with no annotation kept still gets a row, with its annotation blank.
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd

ANNOTATION_FIELDS = ("CSQ", "ANN")
PICKS = ("canonical", "first", "all")
FILTERS = ("pass", "all")
FIXED = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
BASE = ["SAMPLE", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "GT"]


class VcfError(ValueError):
    pass


def stem(path: Path) -> str:
    """The file's name without .vcf or .vcf.gz."""
    name = Path(path).name
    for ending in (".vcf.gz", ".vcf"):
        if name.lower().endswith(ending):
            return name[: -len(ending)]
    return Path(path).stem


def _open(path: Path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def subfields(meta_line: str) -> list[str]:
    """The subfield names an ##INFO line for CSQ or ANN describes: after "Format:" (VEP) or the colon (SnpEff)."""
    m = re.search(r'Description="([^"]*)"', meta_line)
    text = m.group(1) if m else ""
    text = text.split("Format:", 1)[1] if "Format:" in text else (text.split(":", 1)[1] if ":" in text else text)
    return [f.strip() for f in text.strip().strip("'").strip().split("|")]


def _annotations(info: str, field: str, names: list[str]) -> list[dict]:
    for item in info.split(";"):
        if item.startswith(field + "="):
            out = []
            for entry in item[len(field) + 1:].split(","):
                values = entry.split("|")
                out.append({n: (values[i] if i < len(values) else "") for i, n in enumerate(names)})
            return out
    return []


def _allele_forms(ref: str, alt: str) -> set[str]:
    forms = {alt}
    if ref and alt and ref[0] == alt[0] and (len(ref) > 1 or len(alt) > 1):
        forms.add(alt[1:] or "-")
    return forms


def _choose(anns: list[dict], ref: str, alt: str, n_alts: int, pick: str) -> list[dict]:
    forms = _allele_forms(ref, alt)
    mine = [a for a in anns if a.get("Allele", "") in forms]
    if not mine and n_alts == 1:
        mine = anns
    if pick == "canonical":
        return [a for a in mine if a.get("CANONICAL", "") == "YES"]
    return mine[:1] if pick == "first" else mine


def _carried(gt: str) -> set[int]:
    return {int(a) for a in re.split(r"[/|]", gt) if a.isdigit() and int(a) > 0}


def read_vcf(path: str | Path, pick: str, filters: str, field: str | None = None) -> pd.DataFrame:
    """The table above. field: CSQ or ANN, needed only when the header describes both."""
    path = Path(path)
    if pick not in PICKS or filters not in FILTERS:
        raise VcfError(f"{path.name}: pick must be one of {', '.join(PICKS)}, filters one of {', '.join(FILTERS)}")
    described, header, records = {}, None, []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("##"):
                m = re.match(r"##INFO=<ID=(CSQ|ANN),", line)
                if m:
                    described[m.group(1)] = subfields(line)
            elif line.startswith("#"):
                header = line[1:].split("\t")
            elif line.strip():
                if header is None:
                    raise VcfError(f"{path.name}: a record comes before the #CHROM header line")
                records.append(line.split("\t"))
    if header is None or header[:8] != FIXED:
        raise VcfError(f"{path.name}: no VCF header line (#CHROM POS ID REF ALT QUAL FILTER INFO)")
    if field is None:
        if not described:
            raise VcfError(f"{path.name}: the header describes no CSQ (VEP) or ANN (SnpEff) annotation, so the file "
                           "holds no protein changes: annotate it first")
        if len(described) > 1:
            raise VcfError(f"{path.name}: the header describes both CSQ and ANN: say which to read (annotation.field)")
        field = next(iter(described))
    if field not in described:
        raise VcfError(f"{path.name}: the header describes no {field} annotation")
    names = described[field]
    if pick == "canonical" and "CANONICAL" not in names:
        raise VcfError(f"{path.name}: pick canonical needs a CANONICAL subfield (VEP's --canonical), and this {field} "
                       "has none: pick first or all")
    samples = header[9:]
    rows = []
    for r in records:
        r = r + [""] * (len(header) - len(r))
        chrom, pos, vid, ref, alt, qual, flt, info = r[:8]
        if filters == "pass" and flt not in ("PASS", "."):
            continue
        alts = alt.split(",")
        anns = _annotations(info, field, names)
        if samples:
            fmt = r[8].split(":")
            carriers = []
            for name, value in zip(samples, r[9:]):
                parts = value.split(":")
                gt = parts[fmt.index("GT")] if "GT" in fmt and fmt.index("GT") < len(parts) else ""
                alleles = _carried(gt) if "GT" in fmt else set(range(1, len(alts) + 1))
                carriers += [(name, gt, k) for k in sorted(alleles) if k <= len(alts)]
        else:
            carriers = [(stem(path), "", k) for k in range(1, len(alts) + 1)]
        for name, gt, k in carriers:
            a = alts[k - 1]
            base = dict(zip(BASE, (name, chrom, pos, vid, ref, a, qual, flt, gt)))
            chosen = _choose(anns, ref, a, len(alts), pick)
            rows += [base | ann for ann in chosen] or [base | {n: "" for n in names}]
    columns = BASE + [n for n in names if n not in BASE]
    return pd.DataFrame(rows, columns=columns).astype(object)
