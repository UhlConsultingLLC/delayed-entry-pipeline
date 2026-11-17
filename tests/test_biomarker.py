"""MSI/MMR status for a biomarker analysis (S1 in the CRC settings). Ported from the CRC BPC project's
tests/test_msi.py: the same hand-made patients and fixture counts, now through the GENIE BPC adapter
(which reads the pathology results) and the settings' biomarker groups (which say what each reading means)."""
import dataclasses

import pandas as pd
import pytest
import yaml

from delayed_entry.adapters import genie_bpc
from delayed_entry.cohort import build_cohort
from delayed_entry.screen import (BURDEN_HIGH, BURDEN_LOW, DESCRIBE, NOT_RUN, biomarker_status, burden_problems,
                                  mutation_burden, panel_sizes_given, screen)
from delayed_entry.tables import conform, outcome_blind
from helpers import (DOB_DX, FIXTURE, bpc_profile, cancer, crc_settings, hand_made_db, make_tables, path_report, report,
                     settings_with)

S = crc_settings()
S1 = S.analysis("S1")
LABELS = yaml.safe_load(genie_bpc.LABELS.read_text(encoding="utf-8"))
MSS, MSI_L, MSI_H, UNCLEAR_MSI = (next(k for k, v in LABELS["MSI"].items() if v == r)
                                  for r in ("MSS", "MSI-L", "MSI-H", "indeterminate"))
PMMR, DMMR, UNCLEAR = (next(k for k, v in LABELS["MMR"].items() if v == r) for r in ("pMMR", "dMMR", "indeterminate"))
HIGH, STABLE, UNKNOWN = "MSI-H/dMMR", "MSS/pMMR", "Unknown"


def path(specimens=("Colon Cancer",), msi=(), mmr=()):
    """A pathology report: results are (label, day from diagnosis), None for a blank date."""
    def dob(results):
        return [(label, None if day is None else DOB_DX + day) for label, day in results]
    return list(specimens), dob(msi), dob(mmr)


MSI_CASES = {  # record_id: (pathology reports, other cancer's type, BPC cancers, status by the first report, at any time)
    "no_pathology":            ([], None, 1, UNKNOWN, UNKNOWN),
    "mmr_proficient":          ([path(mmr=[(PMMR, 20)])], None, 1, STABLE, STABLE),
    "mmr_deficient":           ([path(mmr=[(DMMR, 20)])], None, 1, HIGH, HIGH),
    "msi_high":                ([path(msi=[(MSI_H, 20)])], None, 1, HIGH, HIGH),
    "msi_low":                 ([path(msi=[(MSI_L, 20)])], None, 1, STABLE, STABLE),       # the settings count MSI-L as stable
    "mss":                     ([path(msi=[(MSS, 20)])], None, 1, STABLE, STABLE),
    "indeterminate_only":      ([path(mmr=[(UNCLEAR, 20)])], None, 1, UNKNOWN, UNKNOWN),
    "on_report_day":           ([path(mmr=[(PMMR, 100)])], None, 1, STABLE, STABLE),      # "on or before"
    "day_after_report":        ([path(mmr=[(PMMR, 101)])], None, 1, UNKNOWN, STABLE),
    "date_blank":              ([path(mmr=[(PMMR, None)])], None, 1, UNKNOWN, STABLE),    # not known to be before
    "discordant":              ([path(msi=[(MSI_H, 20)], mmr=[(PMMR, 20)])], None, 1, HIGH, HIGH),
    "retest_after_report":     ([path(mmr=[(UNCLEAR, 20), (DMMR, 300)])], None, 1, UNKNOWN, HIGH),  # slot 2's own date
    "retest_before_report":    ([path(mmr=[(UNCLEAR, 20), (DMMR, 60)])], None, 1, HIGH, HIGH),
    "third_result":            ([path(mmr=[(UNCLEAR, 10), (UNCLEAR, 30), (DMMR, 60)])], None, 1, HIGH, HIGH),
    "cancer_in_specimen_2":    ([path([None, "Colon Cancer"], mmr=[(PMMR, 20)])], None, 1, STABLE, STABLE),
    "uterus_specimen":         ([path(["Uterus Cancer NOS"], mmr=[(DMMR, 20)])], None, 1, UNKNOWN, UNKNOWN),
    "no_cancer_in_specimen":   ([path([None], mmr=[(PMMR, 20)])], None, 1, UNKNOWN, UNKNOWN),
    "specimen_space":          ([path(["Colon Cancer "], mmr=[(PMMR, 20)])], None, 1, STABLE, STABLE),
    "rectal_specimen":         ([path(["Rectum and Rectosigmoid Cancer"], mmr=[(PMMR, 20)])], None, 1, STABLE, STABLE),
    "colon_rectum_specimen":   ([path(["Colon/Rectum Cancer"], mmr=[(PMMR, 20)])], None, 1, STABLE, STABLE),
    "label_space":             ([path(msi=[(MSS + " ", 20)])], None, 1, STABLE, STABLE),
    "later_report_deficient":  ([path(mmr=[(PMMR, 20)]), path(mmr=[(DMMR, 400)])], None, 1, STABLE, HIGH),
    "other_cancer_prostate":   ([path(mmr=[(PMMR, 20)])], "Prostate Cancer", 1, STABLE, STABLE),
    "other_cancer_colon":      ([path(mmr=[(PMMR, 20)])], "Colon Cancer", 1, None, None),          # left out
    "two_bpc_cancers":         ([path(mmr=[(PMMR, 20)])], None, 2, None, None),                    # left out
    "both_exclusions":         ([path(mmr=[(PMMR, 20)])], "Rectum and Rectosigmoid Cancer", 2, None, None),
}
DISCORDANT = {"discordant"}
DISCORDANT_ANY_TIME = {"discordant", "later_report_deficient"}


@pytest.fixture(scope="module")
def msi_cohort():
    cases = {rid: ([cancer(rid, s) for s in range(n_bpc)], [report(rid, 1, 100.0)])
             for rid, (_, _, n_bpc, _, _) in MSI_CASES.items()}
    cases["not_eligible"] = ([cancer("not_eligible", stage="Stage II")], [report("not_eligible", 1, 100.0)])
    reports = [(rid, n, r) for rid, (rs, *_) in MSI_CASES.items() for n, r in enumerate(rs, start=1)]
    reports.append(("not_eligible", 1, path(mmr=[(PMMR, 20)])))
    pathology = [path_report(rid, n, specs, msi, mmr) for rid, n, (specs, msi, mmr) in reports]
    others = [(rid, 5, "No", other) for rid, (_, other, *_) in MSI_CASES.items() if other]
    con = hand_made_db(cases=cases, pathology=pathology, others=others)
    t = genie_bpc.tables_from_connection(con, S)[0]
    blind = outcome_blind(t)
    return blind, build_cohort(blind, S)


@pytest.mark.parametrize("record_id", list(MSI_CASES))
def test_status_on_hand_made_patients(msi_cohort, record_id):
    blind, c = msi_cohort
    *_, primary, any_time = MSI_CASES[record_id]
    by_first, _, _ = biomarker_status(blind, c, S1, "by_entry_test")
    every, _, _ = biomarker_status(blind, c, S1, "any_time")
    assert by_first.get(record_id) == primary and every.get(record_id) == any_time


def test_discordant_flags_on_hand_made_patients(msi_cohort):
    blind, c = msi_cohort
    for timing, expected in (("by_entry_test", DISCORDANT), ("any_time", DISCORDANT_ANY_TIME)):
        _, discordant, _ = biomarker_status(blind, c, S1, timing)
        assert set(discordant.index[discordant]) == expected


def test_dropping_discordant_patients(msi_cohort):
    blind, c = msi_cohort
    status, _, _ = biomarker_status(blind, c, S1, "by_entry_test", drop_discordant=True)
    assert "discordant" not in status.index and "mss" in status.index


def test_only_eligible_patients_get_a_status(msi_cohort):
    blind, c = msi_cohort
    status, _, _ = biomarker_status(blind, c, S1, "by_entry_test")
    assert set(status.index) == {rid for rid, case in MSI_CASES.items() if case[3] is not None}


def test_what_a_reading_means_comes_from_the_settings(tmp_path, msi_cohort):
    """Counting MSI-L as high instead of stable is a settings change, not a code change."""
    blind, c = msi_cohort
    analyses = S.read("analyses")
    s1 = next(a for a in analyses if a["id"] == "S1")
    s1["groups"] = [{"name": HIGH, "values": ["MSI-H", "dMMR", "MSI-L"]}, {"name": STABLE, "values": ["MSS", "pMMR"]}]
    changed = settings_with(tmp_path, {"analyses": analyses})
    status, _, _ = biomarker_status(blind, c, changed.analysis("S1"), "by_entry_test")
    assert status["msi_low"] == HIGH


def test_before_the_first_report_is_judged_on_the_days_from_birth():
    """BPC's days from birth and from diagnosis can differ by a day. The CRC rule compares the days from birth:
    a result dated on the same day from birth as the report counts, whatever the days from diagnosis say."""
    cases = {"clocks_differ": ([cancer("clocks_differ")], [report("clocks_differ", 1, 100.0)]),
             "report_birth_day_blank": ([cancer("report_birth_day_blank")], [report("report_birth_day_blank", 1, 100.0)])}
    pathology = [path_report("clocks_differ", 1, ["Colon Cancer"], [(MSI_H, DOB_DX + 101)]),
                 path_report("report_birth_day_blank", 1, ["Colon Cancer"], [(MSI_H, DOB_DX + 20)])]
    con = hand_made_db(cases=cases, pathology=pathology)
    con.execute("UPDATE raw_report SET dob_cpt_report_days = dob_cpt_report_days + 1 WHERE record_id = 'clocks_differ'")
    con.execute("UPDATE raw_report SET dob_cpt_report_days = NULL WHERE record_id = 'report_birth_day_blank'")
    blind = outcome_blind(genie_bpc.tables_from_connection(con, S)[0])
    status, _, _ = biomarker_status(blind, build_cohort(blind, S), S1, "by_entry_test")
    assert status["clocks_differ"] == HIGH                 # day 101 from diagnosis, but the report's own birth day
    assert status["report_birth_day_blank"] == UNKNOWN     # a blank on either side is not known to come first


def test_an_unknown_label_stops_the_adapter():
    con = hand_made_db(cases={"odd": ([cancer("odd")], [report("odd", 1, 100.0)])},
                       pathology=[path_report("odd", 1, ["Colon Cancer"], [("MSI-H: HIGH (see comment)", DOB_DX + 20)])])
    with pytest.raises(ValueError, match="see comment"):
        genie_bpc.tables_from_connection(con, S)


# ---- On the fixture ----
@pytest.fixture(scope="module")
def fixture_s1(fixture_result):
    blind = outcome_blind(fixture_result.tables)
    return blind, build_cohort(blind, S)


def test_status_matches_the_truth(fixture_s1):
    """A result from either test is right, except where the two tests disagree (discordant); those are called
    MSI-H/dMMR whatever the truth."""
    blind, c = fixture_s1
    status, discordant, _ = biomarker_status(blind, c, S1, "by_entry_test")
    truth = pd.read_csv(FIXTURE / "truth.csv").set_index("record_id")["true_msi"]
    known = status[status != UNKNOWN]
    agree = known == truth.reindex(known.index)
    assert agree[~discordant.reindex(known.index)].all()
    assert (~agree[discordant.reindex(known.index)]).any()


def test_s1_counts_on_the_fixture(fixture_s1):
    """The counts the CRC project's tests found on this made-up fixture."""
    blind, c = fixture_s1
    by_first, discordant, _ = biomarker_status(blind, c, S1, "by_entry_test")
    every, _, _ = biomarker_status(blind, c, S1, "any_time")
    assert len(by_first) == 652                                       # 658 eligible, 6 left out
    assert [int((by_first == g).sum()) for g in (HIGH, STABLE, UNKNOWN)] == [40, 441, 171]
    assert [int((every == g).sum()) for g in (HIGH, STABLE, UNKNOWN)] == [41, 487, 124]
    assert int(discordant.sum()) == 2


def test_the_fallback_switches_on_above_the_settings_share(tmp_path, fixture_result):
    """171 of 652 Unknown is 26%: under the registered 30%, over a 25% set in the settings file."""
    blind = outcome_blind(fixture_result.tables)
    at_30 = screen(blind, S, fixture_result.profile)
    assert at_30.decisions["SA8"] == {"MSI-H/dMMR vs MSS/pMMR": "compare"}
    at_25 = screen(blind, settings_with(tmp_path, {"biomarker.unknown_max_share": 0.25}), fixture_result.profile)
    assert at_25.decisions["S1"] == {"MSI-H/dMMR vs MSS/pMMR": DESCRIBE}
    assert at_25.decisions["SA8"] == at_25.decisions["SA9"] == {"all comparisons": NOT_RUN}
    assert any("no panel sizes" in i.finding for i in at_25.find(analysis="S1", topic="biomarker unknown"))


def sized(t, size, unsized=()):
    """The tables with every panel given `size` megabases, except those named in `unsized`, which have none."""
    panels = t.panels["panel_id"].dropna().unique()
    sizes = pd.DataFrame({"panel_id": panels, "coding_megabases": [None if p in unsized else size for p in panels]})
    return dataclasses.replace(t, panel_sizes=conform(dataclasses.replace(t, panel_sizes=sizes))[0].panel_sizes)


def test_with_panel_sizes_the_fallback_compares_mutation_burden(tmp_path, fixture_result):
    """On the made-up fixture, entry tests carry 0 to 6 counted mutations: over 0.15 megabases, 3 or more reach 20."""
    blind = outcome_blind(sized(fixture_result.tables, 0.15))
    s = settings_with(tmp_path, {"biomarker.unknown_max_share": 0.25})
    rep = screen(blind, s, bpc_profile(panel_sizes="given"))
    assert rep.decisions["S1"] == {f"{BURDEN_HIGH} vs {BURDEN_LOW}": "compare"}
    assert rep.decisions["SA8"] == rep.decisions["SA9"] == {"all comparisons": NOT_RUN}
    assert set(rep.group_sizes["S1"]) == {BURDEN_LOW, BURDEN_HIGH, UNKNOWN}
    assert any("hypermutation proxy" in i.effect for i in rep.find(analysis="S1", topic="biomarker unknown"))
    c = build_cohort(blind, s)
    group, burden = mutation_burden(blind, c, S1)
    assert ((burden >= 20) == (group == BURDEN_HIGH)).all() and (group != UNKNOWN).all()
    # The profile must say the source gives sizes; otherwise the fallback is descriptive only
    assert screen(blind, s, bpc_profile()).decisions["S1"] == {"MSI-H/dMMR vs MSS/pMMR": DESCRIBE}


def test_a_short_burden_group_is_described(tmp_path, fixture_result):
    """Over 1 megabase no entry test reaches 20: the high group is empty, so the comparison is described."""
    blind = outcome_blind(sized(fixture_result.tables, 1.0))
    rep = screen(blind, settings_with(tmp_path, {"biomarker.unknown_max_share": 0.25}), bpc_profile(panel_sizes="given"))
    assert rep.decisions["S1"] == {f"{BURDEN_HIGH} vs {BURDEN_LOW}": DESCRIBE}


def test_a_burden_that_cant_be_counted_blocks_its_analysis(tmp_path, fixture_result):
    t = sized(fixture_result.tables, 0.15)
    t = dataclasses.replace(t, results=t.results.drop(columns="variant_class"))
    rep = screen(outcome_blind(t), settings_with(tmp_path, {"biomarker.unknown_max_share": 0.25}),
                 bpc_profile(panel_sizes="given"))
    assert rep.decisions["S1"] == {"all comparisons": NOT_RUN}
    block = rep.find(analysis="S1", topic="mutation burden")
    assert block and block[0].severity == "block" and "variant_class" in block[0].finding


# ---- Mutation burden on made-up tables (helpers.make_tables: one call per patient on panel P8) ----
BURDEN_GROUPS = {"RAS-mutant": 40, "BRAF V600E": 35, "Wild type": 40}


def burden_case(mutations=None, sizes=None, panel_of=None, later_test=None, exclude=()):
    """mutations: {patient: [(variant class, n)]} adds n result rows to the patient's entry test (its first, T01).
    sizes: {panel: megabases} (default P8 at 0.5). panel_of: {patient: panel} moves an entry test to another panel
    that reads the same genes. later_test: {patient: n} gives the patient a second report, 200 days later, with n
    missense rows. exclude: patients with more than one cancer of the study's kind (S1 leaves them out)."""
    t = make_tables(BURDEN_GROUPS)
    rows, tests = [], [t.tests]
    for pid, items in (mutations or {}).items():
        rows += [{"test_id": f"{pid}-T01", "kind": "mutation", "gene": "GENE4", "variant": f"p.X{cls[:3]}{i}",
                  "variant_class": cls, "value": None} for cls, n in items for i in range(n)]
    for pid, n in (later_test or {}).items():
        tests.append(pd.DataFrame([{"patient_id": pid, "test_id": f"{pid}-T02", "modality": "ngs", "entry_event": "report",
                                    "days_to_entry_event": 300.0, "panel_id": "P8", "test_order": 2,
                                    "after_follow_up_end": False, "linked_to_other_cancers": False}]))
        rows += [{"test_id": f"{pid}-T02", "kind": "mutation", "gene": "GENE5", "variant": f"p.Y{i}",
                  "variant_class": "Missense_Mutation", "value": None} for i in range(n)]
    test_frame = pd.concat(tests, ignore_index=True)
    for pid, panel in (panel_of or {}).items():
        test_frame.loc[test_frame["test_id"] == f"{pid}-T01", "panel_id"] = panel
    patients = t.patients.copy()
    patients.loc[patients["patient_id"].isin(exclude), "multiple_cohort_cancers"] = True
    sizes = {"P8": 0.5} if sizes is None else sizes
    t = dataclasses.replace(
        t, patients=patients, tests=test_frame, results=pd.concat([t.results, pd.DataFrame(rows)], ignore_index=True),
        panels=pd.concat([t.panels, t.panels.assign(panel_id="P10")], ignore_index=True),
        panel_sizes=pd.DataFrame({"panel_id": list(sizes), "coding_megabases": list(sizes.values())}))
    blind = outcome_blind(conform(t)[0])
    return blind, build_cohort(blind, S)


def test_burden_counts_the_listed_classes_on_the_entry_test_over_its_panel_size():
    """S1 counts nine classes at 20 per megabase. P8 is 0.5 megabases, so 10 counted rows reach 20."""
    blind, c = burden_case({"RAS-0000": [("Missense_Mutation", 9)],                       # the call + 9 = 10
                            "RAS-0001": [("Missense_Mutation", 8)],                       # 9: 18 per megabase
                            "RAS-0002": [("Missense_Mutation", 8), ("Silent", 5)],        # silent isn't counted
                            "RAS-0003": [("Frame_Shift_Del", 4), ("Splice_Site", 5)]},    # 10, by other classes
                           later_test={"RAS-0004": 30})                                   # a later test doesn't count
    group, burden = mutation_burden(blind, c, S1)
    assert burden["RAS-0000"] == 20 and group["RAS-0000"] == BURDEN_HIGH                  # at the threshold: high
    assert burden["RAS-0001"] == 18 and group["RAS-0001"] == BURDEN_LOW
    assert burden["RAS-0002"] == 18 and group["RAS-0002"] == BURDEN_LOW
    assert burden["RAS-0003"] == 20 and group["RAS-0003"] == BURDEN_HIGH
    assert burden["RAS-0004"] == 2 and group["RAS-0004"] == BURDEN_LOW
    assert burden["Wil-0007"] == 2                                                         # one call over 0.5


def test_a_burden_exactly_at_the_threshold_is_high_despite_floating_point():
    """7 mutations over 0.56 megabases is 12.5 exactly, but 7 / 0.56 is 12.499999999999998 in floating point."""
    assert 7 / 0.56 < 12.5
    blind, c = burden_case({"BRA-0000": [("Missense_Mutation", 6)]}, sizes={"P8": 0.56})
    a = {**S1, "mutation_burden": {**S1["mutation_burden"], "threshold_per_megabase": 12.5}}
    group, burden = mutation_burden(blind, c, a)
    assert burden["BRA-0000"] == 12.5 and group["BRA-0000"] == BURDEN_HIGH


def test_an_entry_test_on_a_panel_with_no_size_has_no_burden():
    blind, c = burden_case(panel_of={"Wil-0001": "P10", "Wil-0002": "P10"}, sizes={"P8": 0.5, "P10": None})
    group, burden = mutation_burden(blind, c, S1)
    assert group["Wil-0001"] == group["Wil-0002"] == UNKNOWN and pd.isna(burden["Wil-0001"])
    assert (group.drop(["Wil-0001", "Wil-0002"]) != UNKNOWN).all()


def test_burden_leaves_out_the_patients_the_biomarker_leaves_out():
    blind, c = burden_case(exclude=["RAS-0005"])
    group, _ = mutation_burden(blind, c, S1)
    assert "RAS-0005" not in group.index and "RAS-0006" in group.index


@pytest.mark.parametrize("change, problem", [
    (lambda t: dataclasses.replace(t, results=t.results.drop(columns="variant_class")), "variant_class is absent"),
    (lambda t: dataclasses.replace(t, panel_sizes=None), "panel_sizes table"),
    (lambda t: dataclasses.replace(t, panel_sizes=pd.DataFrame({"panel_id": ["P8", "P8"], "coding_megabases": [0.5, 0.6]})),
     "two different coding sizes"),
    (lambda t: dataclasses.replace(t, panel_sizes=pd.DataFrame({"panel_id": ["P8"], "coding_megabases": [0.0]})),
     "zero or less"),
])
def test_burden_problems_are_named(change, problem):
    blind, _ = burden_case()
    assert burden_problems(blind) == []
    assert any(problem in p for p in burden_problems(change(blind)))


def test_sizes_count_only_when_the_profile_says_the_source_gives_them():
    blind, _ = burden_case()
    assert panel_sizes_given(blind, bpc_profile(panel_sizes="given"))
    assert not panel_sizes_given(blind, bpc_profile())
    assert not panel_sizes_given(dataclasses.replace(blind, panel_sizes=None), bpc_profile(panel_sizes="given"))


def test_a_centers_bed_file_stops_the_adapter(tmp_path):
    import shutil
    data = shutil.copytree(FIXTURE, tmp_path / "data")
    (data / "DFCI-1.bed").write_text("1\t100\t200\n", encoding="utf-8")
    with pytest.raises(NotImplementedError, match="BED file"):
        genie_bpc.read(data, None, S)
