"""The run: the lock first, folders kept apart, the R part, and the safety gate before anything is published.

The tests that run R need it configured: DELAYED_ENTRY_R names Rterm.exe or Rscript, and DELAYED_ENTRY_R_LIBRARY
(optional) a folder of R packages with survival, WeightIt, cobalt, tranSurv, duckdb and jsonlite. Without it they
are skipped. Made-up data only: the GENIE BPC fixture, with few resamples."""
import dataclasses
import json
import os
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import pytest
import yaml

from delayed_entry import run as run_module
from delayed_entry.adapters import genie_bpc
from delayed_entry.gate import Gate, GateError, check_outputs
from delayed_entry.lock import LockError, lock_plan
from delayed_entry.run import OUTPUTS_MANIFEST, R_DIR, ROptions, r_command, run
from delayed_entry.screen import BURDEN_HIGH, BURDEN_LOW, screen
from delayed_entry.settings import load_settings
from delayed_entry.tables import conform, outcome_blind
from helpers import CRC_SETTINGS, FIXTURE, ROOT

R = os.environ.get("DELAYED_ENTRY_R")
needs_r = pytest.mark.skipif(not R, reason="R is not configured (set DELAYED_ENTRY_R)")


def locked(folder: Path, resamples: int = 20, imputing: bool = False, unknown_max_share: float | None = None) -> Path:
    """The CRC settings with few resamples, screened on the fixture and locked in `folder`. Returns the settings.
    imputing: the 'Other' tumor location counts as missing, the complete-case share is lowered below its share of
    the made-up patients, and an imputation model with 3 imputations is given, so the run imputes.
    unknown_max_share: the biomarker's fallback share (the fixture's S1 is 26% Unknown)."""
    folder.mkdir(parents=True, exist_ok=True)
    data = yaml.safe_load(CRC_SETTINGS.read_text(encoding="utf-8"))
    data["uncertainty"]["bootstrap_resamples"]["value"] = resamples
    if unknown_max_share is not None:
        data["biomarker"]["unknown_max_share"]["value"] = unknown_max_share
    if imputing:
        from delayed_entry.settings import SCHEMA_BY_PATH
        for v in data["adjustment"]["variables"]["value"]:
            if v["field"] == "tumor_location":
                v["keep"] = ["Left colon", "Right colon", "Rectal"]
        data["missing_data"]["complete_case_max_share"]["value"] = 0.01
        data["missing_data"]["imputations"]["value"] = 3
        data["missing_data"]["imputation_model"] = {"value": SCHEMA_BY_PATH["missing_data.imputation_model"].example,
                                                    "made_by": "a test", "source": "made up"}
    path = folder / "settings.yaml"
    path.write_bytes(yaml.safe_dump(data, sort_keys=False).encode("utf-8"))
    s = load_settings(path)
    res = genie_bpc.read(FIXTURE, None, s)
    lock_plan(screen(outcome_blind(res.tables), s, res.profile, res.source_steps, res.input_files), s, res.input_files,
              folder, "A. Reviewer", read_notes=res.read_notes)
    return path


# ---- No R needed: what the run refuses before anything is computed ----
def test_the_run_stops_if_the_settings_changed_since_the_lock(tmp_path):
    settings = locked(tmp_path / "work")
    settings.write_bytes(settings.read_bytes().replace(b"value: 20\n", b"value: 21\n", 1))
    with pytest.raises(LockError, match="changed since the lock"):
        run(settings, tmp_path / "work", tmp_path / "out", clinical=FIXTURE, r=ROptions("no-such-program"))
    assert not (tmp_path / "out").exists()


def test_the_run_needs_a_locked_plan(tmp_path):
    (tmp_path / "work").mkdir()
    with pytest.raises(LockError, match="no locked plan"):
        run(CRC_SETTINGS, tmp_path / "work", tmp_path / "out", clinical=FIXTURE, r=ROptions("no-such-program"))


@pytest.mark.parametrize("where", ["inside the work folder", "inside the package"])
def test_outputs_stay_out_of_the_work_folder_and_the_package(tmp_path, where):
    out = tmp_path / "work" / "out" if where == "inside the work folder" else ROOT / "outputs_test"
    with pytest.raises(ValueError, match="must not be inside|outside the pipeline's own folder"):
        run(CRC_SETTINGS, tmp_path / "work", out, clinical=FIXTURE, r=ROptions("no-such-program"))
    assert not (ROOT / "outputs_test").exists()


def test_an_output_folder_holding_other_files_is_refused(tmp_path):
    (tmp_path / "out").mkdir()
    (tmp_path / "out" / "notes.txt").write_text("mine", encoding="utf-8")
    with pytest.raises(ValueError, match="holds other files"):
        run(CRC_SETTINGS, tmp_path / "work", tmp_path / "out", clinical=FIXTURE, r=ROptions("no-such-program"))


def test_r_is_run_as_rterm_dash_f_or_with_rscript(tmp_path):
    term = r_command(r"C:\R\bin\x64\Rterm.exe", "primary", tmp_path)
    assert term[1:5] == ["--no-echo", "--no-restore", "-f", str(R_DIR / "run_analysis.R")] and term[-1] == "primary"
    assert r_command("Rscript", "weighted", tmp_path)[:2] == ["Rscript", str(R_DIR / "run_analysis.R")]


def test_the_r_code_holds_no_study_value():
    """Every value comes from run_spec.json: the R files name no group, gene, threshold or label of a study."""
    text = "\n".join(p.read_text(encoding="utf-8") for p in R_DIR.glob("*.R"))
    for word in ("KRAS", "BRAF", "Wild type", "RAS-mutant", "MSI", "MSK", "1000", "20261108", "730", "0.95", "365.25"):
        assert word not in text, word


# ---- With R ----
@pytest.fixture(scope="module")
def whole_run(tmp_path_factory):
    if not R:
        pytest.skip("R is not configured (set DELAYED_ENTRY_R)")
    base = tmp_path_factory.mktemp("phase2")
    settings = locked(base / "work")
    summary = run(settings, base / "work", base / "out", clinical=FIXTURE,
                  r=ROptions(R, os.environ.get("DELAYED_ENTRY_R_LIBRARY")), log=lambda *_: None)
    return base, settings, summary


@needs_r
def test_r_unit_checks(tmp_path):
    script = ROOT / "tests" / "R" / "test_r_units.R"
    if Path(R).name.lower().startswith("rterm"):
        cmd = [R, "--no-echo", "--no-restore", "-f", str(script), "--args", str(R_DIR), str(tmp_path)]
    else:
        cmd = [R, str(script), str(R_DIR), str(tmp_path)]
    env = {**os.environ, **({"DELAYED_ENTRY_R_LIBRARY": os.environ["DELAYED_ENTRY_R_LIBRARY"]}
                            if os.environ.get("DELAYED_ENTRY_R_LIBRARY") else {})}
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert proc.returncode == 0 and "R unit checks:" in proc.stdout, proc.stderr[-2000:]


@needs_r
def test_a_whole_run_publishes_every_file_and_they_pass_the_scan(whole_run):
    base, _, summary = whole_run
    out = base / "out"
    expected = {"flow.csv", "counts.csv", "descriptives.csv", "risk_table.csv", "curves.csv", "estimates.csv", "bias.csv",
                "balance.csv", "s3_delays.csv", "figure_km.svg", "run_report.md", "run_report.html", OUTPUTS_MANIFEST}
    assert set(summary["files"]) == expected == {f.name for f in out.iterdir()}
    res = genie_bpc.read(FIXTURE, None, load_settings(CRC_SETTINGS))
    ids = {i for i in set(res.tables.patients["patient_id"]) | set(res.tables.tests["test_id"]) if not i.isdigit()}
    assert [p for p in check_outputs(out, Gate(11), ids, set()) if "date" not in p] == []
    assert not [p for p in check_outputs(out, Gate(11), set(), set()) if "date" in p and not p.startswith("run_report")]
    assert all(v["verdict"] in ("pass", "warning", "applied") for v in summary["checks"])


@needs_r
def test_the_outputs_manifest_holds_each_files_checksum(whole_run):
    from delayed_entry.hashing import sha256_file
    out = whole_run[0] / "out"
    m = json.loads((out / OUTPUTS_MANIFEST).read_text(encoding="utf-8"))
    assert all(sha256_file(out / name) == digest for name, digest in m["outputs"].items())
    assert all("path" not in i for i in m["inputs"])                # names and checksums, never local paths


@needs_r
def test_the_tables_hold_every_analysis_the_plan_runs(whole_run):
    out = whole_run[0] / "out"
    bias = pd.read_csv(out / "bias.csv", keep_default_na=False)
    assert {"primary", "3a", "3b", "conditional on t*", "S1", "S1 3a", "S2", "SA1", "SA2", "SA3", "SA4", "SA7"} <= set(bias["analysis"])
    est = pd.read_csv(out / "estimates.csv", keep_default_na=False)
    assert {"SA5", "SA6"} <= set(est["analysis"])
    labels = set(bias["label"])
    assert labels <= {"clearly material", "material", "not material", "all patients: described, not labeled",
                      "not compared (smallest group)", "described, not labeled (smallest group)",
                      "described, not labeled (biomarker fallback)", "described: no death before t*",
                      "median not reached within 5 years: RMST is the summary",
                      "described: corrected risk set never reaches 20"}


@needs_r
def test_the_figure_is_well_formed_svg_with_text_kept_as_text(whole_run):
    svg = (whole_run[0] / "out" / "figure_km.svg").read_text(encoding="utf-8")
    root = ET.fromstring(svg)
    texts = "".join(t.text or "" for t in root.iter("{http://www.w3.org/2000/svg}text"))
    assert all(g in texts for g in ("Wild type", "RAS-mutant", "BRAF V600E"))


@needs_r
def test_the_run_report_says_what_was_done_and_why(whole_run):
    md = (whole_run[0] / "out" / "run_report.md").read_text(encoding="utf-8")
    for heading in ("## Inputs", "## Capability profile of the source", "## The readiness screen, from the locked plan",
                    "## Checks, with their verdicts", "## The safety gate", "## Results in brief"):
        assert heading in md
    assert "survival after the entry event" in md and "the locked plan, settings, defaults and inputs are unchanged" in md
    assert "<!doctype html>" in (whole_run[0] / "out" / "run_report.html").read_text(encoding="utf-8")


@needs_r
def test_a_genie_bpc_run_prints_the_citation_and_warns_while_the_acknowledgment_is_blank(whole_run):
    _, _, summary = whole_run
    md = (whole_run[0] / "out" / "run_report.md").read_text(encoding="utf-8")
    assert "## Data source and attribution" in md and "Cancer Discov. 2017;7(8):818-831" in md
    assert "Dataset version: GENIE BPC crc_v2.0 (Synapse syn30991637)" in md
    check = next(c for c in summary["checks"] if c["check"] == "the attribution the source's terms require")
    if genie_bpc.profile().attribution["acknowledgment"]:
        assert check["verdict"] == "pass" and "NOT FILLED IN" not in md
    else:
        assert check["verdict"] == "warning" and "Acknowledgment: NOT FILLED IN" in md


@needs_r
def test_a_second_run_with_the_same_inputs_keeps_the_r_results(whole_run, tmp_path):
    base, settings, _ = whole_run
    summary = run(settings, base / "work", tmp_path / "again", clinical=FIXTURE,
                  r=ROptions(R, os.environ.get("DELAYED_ENTRY_R_LIBRARY"), reuse=True), log=lambda *_: None)
    assert [c for c in summary["checks"] if c["stage"] == "R"] == [{"stage": "R", "check": "the R part", "verdict": "pass"}]
    a = (base / "out" / "bias.csv").read_text(encoding="utf-8")
    assert (tmp_path / "again" / "bias.csv").read_text(encoding="utf-8") == a


@needs_r
def test_nothing_is_published_when_the_scan_fails(whole_run, tmp_path, monkeypatch):
    base, settings, _ = whole_run
    real = run_module.build_tables
    res = genie_bpc.read(FIXTURE, None, load_settings(CRC_SETTINGS))
    an_id = sorted(i for i in res.tables.patients["patient_id"] if not i.isdigit())[0]

    def leaky(results):
        tables, figure, withheld = real(results)
        tables["flow.csv"].loc[0, "label"] = an_id                  # a patient's ID slips into a table
        return tables, figure, withheld
    monkeypatch.setattr(run_module, "build_tables", leaky)
    with pytest.raises(GateError, match="holds an ID from the data"):
        run(settings, base / "work", tmp_path / "leak", clinical=FIXTURE,
            r=ROptions(R, os.environ.get("DELAYED_ENTRY_R_LIBRARY"), reuse=True), log=lambda *_: None)
    assert not (tmp_path / "leak").exists()


# ---- With R and multiple imputation ----
@pytest.fixture(scope="module")
def imputed_run(tmp_path_factory):
    if not R:
        pytest.skip("R is not configured (set DELAYED_ENTRY_R)")
    base = tmp_path_factory.mktemp("imputed")
    settings = locked(base / "work", resamples=10, imputing=True)
    summary = run(settings, base / "work", base / "out", clinical=FIXTURE,
                  r=ROptions(R, os.environ.get("DELAYED_ENTRY_R_LIBRARY")), log=lambda *_: None)
    return base, summary


@needs_r
def test_an_imputed_run_says_so_and_passes_the_gate(imputed_run):
    base, summary = imputed_run
    check = next(c for c in summary["checks"] if c["check"] == "the adjustment variables")
    assert check["verdict"] == "applied"
    md = (base / "out" / "run_report.md").read_text(encoding="utf-8")
    assert "multiple imputation (3 imputations" in md and "Rubin's rules" in md
    assert any(g["check"].startswith("balance") for g in summary["gate"])


@needs_r
def test_imputed_hazard_ratios_are_pooled_and_their_contrasts_averaged(imputed_run):
    out = imputed_run[0] / "out"
    est = pd.read_csv(out / "estimates.csv", keep_default_na=False)
    hr = est[(est["analysis"] == "3a") & est["quantity"].str.startswith("HR")]
    assert len(hr) and (hr["note"] == "pooled over 3 imputations by Rubin's rules").all()
    assert (hr["lower"].astype(float) < hr["estimate"].astype(float)).all()
    assert (hr["estimate"].astype(float) < hr["upper"].astype(float)).all()
    w = est[(est["analysis"] == "3b") & (est["quantity"] == "median, months")]
    assert w["note"].str.contains("averaged over 3 imputations").all()
    bias = pd.read_csv(out / "bias.csv", keep_default_na=False)
    adjusted = bias[bias["analysis"].isin(["3a", "3b", "SA2"])]
    assert adjusted["note"].str.contains("each bootstrap resample imputed once").all()
    primary = bias[bias["analysis"] == "primary"]
    assert not primary["note"].str.contains("imputation").any()     # passes 1 and 2 never use the covariates


@needs_r
def test_an_imputed_run_is_reproducible(imputed_run, tmp_path):
    base, _ = imputed_run
    run(base / "work" / "settings.yaml", base / "work", tmp_path / "again", clinical=FIXTURE,
        r=ROptions(R, os.environ.get("DELAYED_ENTRY_R_LIBRARY")), log=lambda *_: None)
    for name in ("bias.csv", "estimates.csv", "balance.csv"):
        assert (tmp_path / "again" / name).read_bytes() == (base / "out" / name).read_bytes(), name


# ---- A biomarker's fallback by mutation burden ----
def with_panel_sizes(real_read):
    """The GENIE BPC adapter's read, plus a made-up coding size of 0.15 megabases for every panel, and a profile
    saying the source gives them. On the fixture, 3 or more counted mutations reach 20, and every patient in S1 has a
    burden: no one is Unknown, and the few S1 sets aside are a small count the groups would give away."""
    def read(clinical, genomic, settings):
        res = real_read(clinical, genomic, settings)
        t = res.tables
        panels = t.panels["panel_id"].dropna().unique()
        sizes = pd.DataFrame({"panel_id": panels, "coding_megabases": 0.15})
        tables = conform(dataclasses.replace(t, panel_sizes=sizes))[0]
        return dataclasses.replace(res, tables=tables, profile=dataclasses.replace(res.profile, panel_sizes="given"))
    return read


@pytest.fixture(scope="module")
def burden_run(tmp_path_factory):
    if not R:
        pytest.skip("R is not configured (set DELAYED_ENTRY_R)")
    base = tmp_path_factory.mktemp("burden")
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(genie_bpc, "read", with_panel_sizes(genie_bpc.read))
        settings = locked(base / "work", resamples=10, unknown_max_share=0.25)
        summary = run(settings, base / "work", base / "out", clinical=FIXTURE,
                      r=ROptions(R, os.environ.get("DELAYED_ENTRY_R_LIBRARY")), log=lambda *_: None)
    return base, summary


@needs_r
def test_a_burden_fallback_run_compares_high_with_low_and_passes_the_gate(burden_run):
    base, summary = burden_run
    out = base / "out"
    assert all(v["verdict"] in ("pass", "warning", "applied") for v in summary["checks"])
    bias = pd.read_csv(out / "bias.csv", keep_default_na=False)
    s1 = bias[bias["analysis"].isin(["S1", "S1 3a"])]
    assert {"S1", "S1 3a"} <= set(s1["analysis"])
    contrast = f"{BURDEN_HIGH} vs {BURDEN_LOW}"
    assert set(s1["group"]) == {"All patients", BURDEN_HIGH, BURDEN_LOW, contrast}
    labels = s1.loc[s1["group"] == contrast, "label"]
    assert len(labels) and labels.isin(["clearly material", "material", "not material"]).all()   # compared
    assert not {"SA8", "SA9"} & set(bias["analysis"])
    est = pd.read_csv(out / "estimates.csv", keep_default_na=False)
    assert {BURDEN_HIGH, BURDEN_LOW} <= set(est.loc[est["analysis"] == "S1", "group"])
    md = (out / "run_report.md").read_text(encoding="utf-8")
    assert "hypermutation proxy" in md and "mutation burden on the entry test" in md
    # The gate: the larger group's size is hidden, so the few set aside can't be found by subtraction
    counts = pd.read_csv(out / "counts.csv", keep_default_na=False).set_index(["table", "group"])["patients"]
    assert counts[("S1", BURDEN_LOW)] == "not shown" and counts[("S1", BURDEN_HIGH)] != "not shown"
