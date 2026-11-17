"""Phase 1 acceptance on the real CRC BPC data. Runs on your machine only, against the data already there.

    python scripts/acceptance_phase1.py <clinical folder> <genomic folder> [<work folder>]

Reads the release with the GENIE BPC adapter and the CRC settings, screens it, and compares the screen's
decisions and reasons with the registered project's design gate (crc_bpc_phase1_expected.yml). It also
checks, in memory, that scrambling every outcome value leaves the readiness report unchanged.

It prints only labels and pass or fail: no count, no patient. With a work folder (local, not
cloud-synced, outside this package), it also writes the readiness report there.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))

from delayed_entry.adapters import genie_bpc  # noqa: E402
from delayed_entry.privacy import check_not_cloud  # noqa: E402
from delayed_entry.screen import screen  # noqa: E402
from delayed_entry.settings import load_settings  # noqa: E402
from delayed_entry.tables import StandardTables, outcome_blind  # noqa: E402

SETTINGS = HERE.parent / "settings" / "crc_bpc.yaml"
EXPECTED = HERE / "crc_bpc_phase1_expected.yml"


def scrambled(t: StandardTables, seed: int) -> StandardTables:
    """Every recorded outcome value replaced by a random one; blanks stay blank."""
    rng = np.random.default_rng(seed)
    p = t.patients.copy()
    known = p["vital_status"].notna()
    p.loc[known, "vital_status"] = rng.integers(0, 2, int(known.sum()))
    known = p["days_to_death_or_last_followup"].notna()
    p.loc[known, "days_to_death_or_last_followup"] = rng.uniform(1, 4000, int(known.sum()))
    return StandardTables(p, t.tests, t.results, t.panels, t.treatments, t.panel_sizes)


def main(argv: list[str]) -> int:
    clinical, genomic = Path(argv[0]), Path(argv[1])
    work = Path(argv[2]) if len(argv) > 2 else None
    settings = load_settings(SETTINGS)
    expected = yaml.safe_load(EXPECTED.read_text(encoding="utf-8"))
    result = genie_bpc.read(clinical, genomic, settings)
    report = screen(outcome_blind(result.tables), settings, result.profile, result.source_steps)
    ok_all = True

    print("Decisions (expected from the registered design gate):")
    for analysis, comps in expected["decisions"].items():
        for comp, want in comps.items():
            got = report.decisions.get(analysis, {}).get(comp)
            ok = got == want
            ok_all &= ok
            print(f"  {'PASS' if ok else 'FAIL'}  {analysis}: {comp}: {got}")
    print("Reasons:")
    for r in expected["reasons"]:
        hits = [i for i in report.find(analysis=r["analysis"], topic=r["topic"]) if r["contains"] in i.finding]
        ok = bool(hits)
        ok_all &= ok
        print(f"  {'PASS' if ok else 'FAIL'}  {r['analysis']}: {r['topic']}: \"{r['contains']}\"")
    print("No outcome read:")
    same = all(screen(outcome_blind(scrambled(result.tables, seed)), settings, result.profile,
                      result.source_steps).to_dict() == report.to_dict() for seed in (1, 2))
    ok_all &= same
    print(f"  {'PASS' if same else 'FAIL'}  scrambling every outcome value leaves the readiness report unchanged")
    blocks = [i for i in report.items if i.severity == "block" and i.analysis not in ("SA8", "SA9")]
    ok_all &= not blocks
    print(f"  {'PASS' if not blocks else 'FAIL'}  no block outside SA8 and SA9")
    if work is not None:
        work = check_not_cloud(work.resolve())
        if (HERE.parent.resolve() in work.parents) or work == HERE.parent.resolve():
            raise SystemExit("the work folder must be outside this package")
        work.mkdir(parents=True, exist_ok=True)
        (work / "readiness_report.md").write_bytes(report.to_markdown().encode("utf-8"))
        (work / "readiness_report.yaml").write_bytes(report.to_yaml().encode("utf-8"))
        print("Report written to the work folder.")
    print("RESULT:", "PASS" if ok_all else "FAIL")
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
