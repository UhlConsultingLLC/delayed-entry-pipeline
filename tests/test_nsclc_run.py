"""Whole runs with the NSCLC settings on the made-up NSCLC-shaped files (test_nsclc_spec's fixture): the lock, the R
part, the reports and the safety gate.

Two runs. In the fixture as it is, the made-up ALK fusion group is three patients: too few for SA6's transformation
model, which must be reported for that group rather than stop the R part; the safety gate then stops the run, since
counts that small could be worked out (a person must merge or hide). With fifteen more made-up EML4-ALK patients the
run publishes, a Kendall's tau with no standard error is reported as not assessed rather than stopping the reports,
and the described-only group's comparison says the settings describe it (question 13), not that it is short.

The whole runs need R configured, as test_run's do (DELAYED_ENTRY_R, and DELAYED_ENTRY_R_LIBRARY if the packages
are elsewhere); without it they are skipped. Made-up data only, with few resamples."""
import os
import shutil

import duckdb
import pandas as pd
import pytest

from delayed_entry.adapters import genie_bpc
from delayed_entry.cohort import build_cohort
from delayed_entry.gate import Gate, GateError, check_outputs
from delayed_entry.lock import lock_plan
from delayed_entry.reports import is_true, read_parquet, tau_verdict
from delayed_entry.run import OUTPUTS_MANIFEST, ROptions, run
from delayed_entry.screen import screen
from delayed_entry.settings import load_settings
from delayed_entry.tables import outcome_blind
from test_nsclc_spec import FUSION_COLUMNS, dfci_rows, nsclc_files, nsclc_settings_file, read, write  # noqa: F401

R = os.environ.get("DELAYED_ENTRY_R")
MORE_ALK = 15


def locked_run(files, base):
    """The NSCLC settings with few resamples, screened on `files` and locked in base/work; then the run into base/out.
    Returns the settings, what the adapter read, and the run's summary (or the error that stopped it)."""
    work = base / "work"
    work.mkdir(parents=True)
    settings = nsclc_settings_file(work, {"uncertainty.bootstrap_resamples": 20, "missing_data.imputations": 3})
    s = load_settings(settings)
    res = genie_bpc.read(files, None, s)
    lock_plan(screen(outcome_blind(res.tables), s, res.profile, res.source_steps, res.input_files), s,
              res.input_files, work, "A. Reviewer", read_notes=res.read_notes)
    try:
        summary = run(settings, work, base / "out", clinical=files,
                      r=ROptions(R, os.environ.get("DELAYED_ENTRY_R_LIBRARY")), log=lambda *_: None)
    except GateError as err:
        summary = err
    return s, res, summary


def with_more_alk(files, dst, tmp):
    """A copy of the made-up files where the last MORE_ALK made-up reference patients gain an EML4-ALK fusion
    (DFCI-style rows for both genes, on every sample), so the made-up ALK fusion group is large enough to publish."""
    tmp.mkdir(parents=True, exist_ok=True)
    s = load_settings(nsclc_settings_file(tmp, {}))
    c = build_cohort(outcome_blind(genie_bpc.read(files, None, s).tables), s)
    chosen = sorted(c.eligible.index[c.eligible["genotype"] == "No recognized driver"])[-MORE_ALK:]
    shutil.copytree(files, dst)
    reports = read(dst / "cancer_panel_test_level_dataset.csv")
    rows = [r for pid in chosen for r in dfci_rows(reports, pid, ("EML4", "ALK"), ["EML4", "ALK"])]
    fusions = read(dst / "data_fusions.txt", "\t")
    write(pd.concat([fusions, pd.DataFrame(rows, columns=FUSION_COLUMNS)], ignore_index=True),
          dst / "data_fusions.txt", "\t")
    return dst


@pytest.fixture(scope="module")
def nsclc_run(nsclc_files, tmp_path_factory):  # noqa: F811
    if not R:
        pytest.skip("R is not configured (set DELAYED_ENTRY_R)")
    base = tmp_path_factory.mktemp("nsclc_run")
    files = with_more_alk(nsclc_files, base / "files", base / "tmp")
    s, res, summary = locked_run(files, base)
    if isinstance(summary, Exception):
        raise summary
    return base, s, res, summary


@pytest.fixture(scope="module")
def three_alk_run(nsclc_files, tmp_path_factory):  # noqa: F811
    if not R:
        pytest.skip("R is not configured (set DELAYED_ENTRY_R)")
    base = tmp_path_factory.mktemp("three_alk_run")
    return (base, *locked_run(nsclc_files, base))


def sparse_smoking(files, dst, tmp):
    """A copy of the made-up files where every made-up ALK fusion patient's smoking is "Never used" but two, whose
    smoking is "Unknown": a sparse genotype-by-category cell for the multinomial imputation of smoking."""
    tmp.mkdir(parents=True, exist_ok=True)
    s = load_settings(nsclc_settings_file(tmp, {}))
    c = build_cohort(outcome_blind(genie_bpc.read(files, None, s).tables), s)
    alk = sorted(c.eligible.index[c.eligible["genotype"] == "ALK fusion"])
    shutil.copytree(files, dst)
    ca = read(dst / "cancer_level_dataset_index.csv")
    ca.loc[ca["record_id"].isin(alk), "ca_lung_cigarette"] = "Never used"
    ca.loc[ca["record_id"].isin(alk[:2]), "ca_lung_cigarette"] = "Unknown"
    write(ca, dst / "cancer_level_dataset_index.csv")
    return dst, alk


@pytest.fixture(scope="module")
def sparse_run(nsclc_files, tmp_path_factory):  # noqa: F811
    if not R:
        pytest.skip("R is not configured (set DELAYED_ENTRY_R)")
    base = tmp_path_factory.mktemp("sparse_run")
    files, alk = sparse_smoking(with_more_alk(nsclc_files, base / "more_alk", base / "tmp1"), base / "files",
                                base / "tmp2")
    _, _, summary = locked_run(files, base)
    return base, alk, summary


def tau_table(base):
    return read_parquet(next((base / "work").glob("run_*/r_tau.parquet")))


# ---- No R needed ----
def test_a_tau_row_with_no_interval_is_reported_as_not_assessed(tmp_path):
    """A missing flag (R's NA, read back as a missing value) is neither flagged nor cleared, and never stops the
    report. Made-up rows written as R writes them."""
    p = (tmp_path / "r_tau.parquet").as_posix()
    duckdb.execute(f"""COPY (SELECT * FROM (VALUES ('MB', 'A', 0.1, TRUE, ''),
                                ('IPW2', 'A', 0.2, NULL, 'no standard error, so no interval'),
                                ('IPW1', 'B', NULL, NULL, 'the estimator stopped: made up'),
                                ('MB', 'B', -0.1, FALSE, '')) t(estimator, "group", tau, flagged, tau_note))
                       TO '{p}' (FORMAT parquet)""")
    t = read_parquet(p)
    assert [tau_verdict(r) for r in t.itertuples()] == [
        "flagged", "not assessed: no standard error, so no interval",
        "not assessed: the estimator stopped: made up", "not flagged"]
    assert [is_true(v) for v in t["flagged"]] == [True, False, False, False]
    older = t.drop(columns="tau_note")                       # a run folder written before tau_note existed
    assert tau_verdict(next(older.iloc[[1]].itertuples())) == "not assessed: no interval"


# ---- With R ----
def test_a_three_patient_group_gets_through_the_r_part_and_the_gate_then_stops_the_run(three_alk_run):
    """SA6's transformation model can't be fitted to the made-up three-patient ALK fusion group: R records why, and
    the run goes on to the safety gate, which stops it because counts that small could be worked out."""
    base, _, _, summary = three_alk_run
    assert isinstance(summary, GateError) and "counts below 11 could be worked out" in str(summary)
    assert "ALK fusion" in str(summary)
    tau = tau_table(base)
    alk = tau[tau["group"] == "ALK fusion"].iloc[0]
    assert pd.isna(alk["median_transformation"]) and alk["warning"].startswith("transformation model not fitted:")
    assert not (base / "out").exists() or not any((base / "out").iterdir())


def test_an_nsclc_shaped_run_publishes_every_file_and_passes_the_gate(nsclc_run):
    base, _, res, summary = nsclc_run
    out = base / "out"
    assert {"flow.csv", "counts.csv", "estimates.csv", "bias.csv", "figure_km.svg", "run_report.md",
            OUTPUTS_MANIFEST} <= set(summary["files"]) == {f.name for f in out.iterdir()}
    ids = {i for i in set(res.tables.patients["patient_id"]) | set(res.tables.tests["test_id"]) if not i.isdigit()}
    assert [p for p in check_outputs(out, Gate(11), ids, set()) if "date" not in p] == []
    assert all(v["verdict"] in ("pass", "warning", "applied") for v in summary["checks"])


def test_the_tables_hold_the_nsclc_analyses_the_plan_runs(nsclc_run):
    base, s, _, _ = nsclc_run
    est = pd.read_csv(base / "out" / "estimates.csv", keep_default_na=False)
    bias = pd.read_csv(base / "out" / "bias.csv", keep_default_na=False)
    medians = est[(est["analysis"] == "SA6") & (est["quantity"] == "median, months")]
    assert {g["name"] for g in s.read("genotype.groups")} == set(medians["group"])
    assert {"primary", "SA13"} <= set(bias["analysis"])


def test_the_described_only_groups_comparison_says_the_settings_describe_it(nsclc_run):
    """Question 13, option B: Other recognized driver is described, never compared. Its comparison rows say the
    settings are why, in every analysis, never "not compared (smallest group)"; in the primary its own estimates are
    described too, and in the subsets they keep their labels (the settings' note)."""
    base, s, _, _ = nsclc_run
    assert s.read("genotype.described_only") == ["Other recognized driver"]
    bias = pd.read_csv(base / "out" / "bias.csv", keep_default_na=False)
    comp = bias[bias["group"] == "Other recognized driver vs No recognized driver"]
    assert {"primary", "3a", "3b", "SA13"} <= set(comp["analysis"])
    assert (comp["label"] == "described, not labeled (settings)").all()
    own = bias[bias["group"] == "Other recognized driver"]
    assert (own.loc[own["analysis"] == "primary", "label"] == "described, not labeled (settings)").all()
    assert not own.loc[own["analysis"] == "SA13", "label"].str.startswith("described, not labeled").any()
    alk = bias[(bias["group"] == "ALK fusion vs No recognized driver") & (bias["analysis"] == "primary")]
    assert len(alk) and (alk["label"] == "not compared (smallest group)").all()    # 18 made-up patients: short


def test_the_run_splits_the_overlap_count_by_gene(nsclc_run):
    """6a (h): the NSCLC groups share EGFR, so descriptives.csv gives the overlaps through calls in different genes
    and in one gene, each masked on its own; the total only when neither part is masked."""
    base, _, _, _ = nsclc_run
    d = pd.read_csv(base / "out" / "descriptives.csv", keep_default_na=False)
    rows = dict(zip(d.loc[d["table"] == "Kept but counted", "group"], d.loc[d["table"] == "Kept but counted", "patients"]))
    cross = rows["Matched more than one group's rules through calls in different genes, so counted in the first by precedence"]
    one = rows.get("Matched more than one group's rules through calls in one gene only, so counted in the first by precedence")
    total = rows.get("Matched more than one group's rules, so counted in the first by precedence")
    assert one is not None
    if "<11" in (cross, one):
        assert total is None
    else:
        assert int(total) == int(cross) + int(one)


def test_the_nsclc_shaped_run_imputes(nsclc_run):
    """6a (m): about a sixth of the made-up patients have smoking Unknown, so the run imputes (question 11's model)."""
    base, _, _, summary = nsclc_run
    check = next(c for c in summary["checks"] if c["check"] == "the adjustment variables")
    assert check["verdict"] == "applied"
    assert "multiple imputation (3 imputations" in (base / "out" / "run_report.md").read_text(encoding="utf-8")
    completed = read_parquet(next((base / "work").glob("run_*/r_imputed.parquet")))
    assert sorted(completed["imputation"].unique()) == [1, 2, 3]
    assert set(completed["smoking_at_diagnosis"].dropna()) <= {"Never", "Former", "Current"}
    assert completed["smoking_at_diagnosis"].notna().all()


def test_a_sparse_genotype_by_smoking_cell_still_imputes_and_publishes(sparse_run):
    """6a (m): every made-up ALK fusion patient with a known smoking value is "Never used". The multinomial draw's
    guard against perfect prediction (mice's augmenting records; checked in tests/R/test_r_units.R) keeps the other
    categories possible for the cell's two Unknown patients; the run finishes and passes the gate."""
    base, alk, summary = sparse_run
    assert not isinstance(summary, Exception), summary
    assert all(v["verdict"] in ("pass", "warning", "applied") for v in summary["checks"])
    completed = read_parquet(next((base / "work").glob("run_*/r_imputed.parquet")))
    unknown = completed[completed["patient_id"].isin(alk[:2])]
    assert len(unknown) == 2 * 3 and unknown["smoking_at_diagnosis"].notna().all()
    assert set(unknown["smoking_at_diagnosis"]) <= {"Never", "Former", "Current"}


def test_a_tau_with_no_standard_error_is_published_as_not_assessed(nsclc_run):
    """The made-up run has at least one tau with no standard error (the fault this guards against stopped the
    reports): its published row says it wasn't assessed, unless the small-group rule blanks the group's rows."""
    base, _, _, _ = nsclc_run
    tau = tau_table(base)
    missing = tau[tau["flagged"].isna()]
    assert len(missing) and (missing["tau_note"] == "no standard error, so no interval").all()
    est = pd.read_csv(base / "out" / "estimates.csv", keep_default_na=False)
    tau_rows = est[(est["analysis"] == "SA6") & est["quantity"].str.startswith("conditional Kendall's tau")]
    for r in missing.itertuples():
        row = tau_rows[(tau_rows["group"] == r.group) & tau_rows["quantity"].str.endswith(f"({r.estimator})")]
        assert row["note"].iloc[0] in ("not assessed: no standard error, so no interval", "fewer than 11 patients")
    assert (tau_rows["note"] == "not assessed: no standard error, so no interval").any()
