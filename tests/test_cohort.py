"""The cohort rules, the first test and the genotype groups, applied from the CRC settings to tables the
GENIE BPC adapter built. Ported from the CRC BPC project's tests/test_cohort.py: the same hand-made
patients and the same made-up fixture, with the same expected results, except where noted (two cases
whose exclusion needs outcome times, which the screen never reads)."""
import pandas as pd
import pytest

from delayed_entry.adapters.genie_bpc import raw
from delayed_entry.cohort import build_cohort, codon_span, field_condition
from delayed_entry.privacy import is_small, public_flow
from delayed_entry.tables import outcome_blind
from helpers import FIXTURE, cancer, crc_settings, hand_made_tables, report

S = crc_settings()


def cohort_of(**kw):
    t = hand_made_tables(S, **kw)
    return build_cohort(outcome_blind(t), S), t


def excluded_at(c, rid):
    return c.patients.loc[rid, "excluded_at"] if rid in c.patients.index else "not in the tables"


# ---- Hand-made patients: one per situation each rule must handle ----
CASES = {  # record_id: (its cancer rows, the step it should be excluded at; None = eligible)
    "ok": ([cancer("ok")], None),
    "two_cancers": ([cancer("two_cancers", 1, stage="Stage II"), cancer("two_cancers", 2)], "stage_iv"),
    "not_index": ([cancer("not_index", index="No")], "not in the tables"),     # left out by the adapter
    "no_cancer": ([], "not in the tables"),
    "stage_space": ([cancer("stage_space", stage="Stage IV ")], None),
    "stage_missing": ([cancer("stage_missing", stage=None)], "stage_iv"),
    "squamous_space": ([cancer("squamous_space", histology="Squamous cell ")], "histology"),
    "sarcoma": ([cancer("sarcoma", histology="Sarcoma")], "histology"),
    "histology_missing": ([cancer("histology_missing", histology=None)], None),
    "age_17": ([cancer("age_17", age=17.0)], "age_18"),
    "age_missing": ([cancer("age_missing", age=None)], "entry_known"),  # passes the age rule; the entry rule takes it
    "appendix_dot": ([cancer("appendix_dot", site="C18.1")], "not_appendix"),
    "appendix_no_dot": ([cancer("appendix_no_dot", site="C181")], "not_appendix"),
    "appendix_text": ([cancer("appendix_text", site="c18.1 Appendix")], "not_appendix"),
    "cecum": ([cancer("cecum", site="C18.0")], None),
    "site_missing": ([cancer("site_missing", site=None)], None),
    "two_failures": ([cancer("two_failures", stage="Stage II", histology="Squamous cell")], "stage_iv"),
}


@pytest.fixture(scope="module")
def hand_made():
    c, _ = cohort_of(cases={rid: (cancers, [report(rid, i, 30.0, ca_seq=row[1]) for i, row in enumerate(cancers, start=1)])
                            for rid, (cancers, _) in CASES.items()})
    return c


@pytest.mark.parametrize("record_id", list(CASES))
def test_each_rule_on_hand_made_patients(hand_made, record_id):
    assert excluded_at(hand_made, record_id) == CASES[record_id][1]


# ---- The first test ----
REPORT_CASES = {  # record_id: (its cancers' ca_seq, its reports, the cpt_number chosen as first, the step excluded at)
    "one_report":            ([0], [(1, 50)], 1, None),
    "earlier_wins":          ([0], [(1, 400), (2, 50)], 2, None),               # the date decides, not the number
    "same_day_bigger_panel": ([0], [(1, 50, "P8"), (2, 50, "P10")], 2, None),
    "date_beats_panel":      ([0], [(1, 50, "P8"), (2, 60, "P10")], 1, None),
    "same_day_same_panel":   ([0], [(3, 50, "P10"), (2, 50, "P10")], 2, None),   # listed out of order on purpose
    "before_diagnosis":      ([0], [(1, 100), (2, -20)], 2, None),               # enters at day 0, flagged
    "time_withheld_and_known": ([0], [(1, None, "P10"), (2, 400)], 2, None),     # flagged: unranked
    "time_withheld_only":    ([0], [(1, None)], 1, "entry_known"),               # chosen, but its time is unknown
    "panel_without_file":    ([0], [(1, 50, "NO-FILE"), (2, 50, "P8")], 2, None),
    "other_cancer_earlier":  ([0, 1], [(1, 10, "P8", 1), (2, 300, "P8", 0)], 2, None),
    "no_report":             ([0], [], None, "has_report"),
    "other_cancer_only":     ([0, 1], [(1, 10, "P8", 1)], None, "has_report"),
    "two_cancer_test_two_rows": ([0, 1], [(1, 50, "P8", 0, 2), (1, 50, "P8", 1, 2)], 1, None),   # flagged
    "two_cancer_test_one_row":  ([0], [(1, 50, "P8", 0, 2)], 1, None),                          # flagged
}
FLAGGED = {"f:first_test_linked_to_other_cancers": {"two_cancer_test_two_rows", "two_cancer_test_one_row"},
           "f:unranked_first_test": {"time_withheld_and_known"},
           "f:entry_before_diagnosis": {"before_diagnosis"}}


@pytest.fixture(scope="module")
def report_cohort():
    return cohort_of(cases={rid: ([cancer(rid, s) for s in seqs], [report(rid, *r) for r in reports])
                            for rid, (seqs, reports, _, _) in REPORT_CASES.items()})[0]


@pytest.mark.parametrize("record_id", list(REPORT_CASES))
def test_first_test_on_hand_made_patients(report_cohort, record_id):
    _, _, first, step = REPORT_CASES[record_id]
    got = report_cohort.patients.loc[record_id, "first_test_id"]
    assert (None if pd.isna(got) else got) == (None if first is None else f"{record_id}-T{first:02d}")
    assert excluded_at(report_cohort, record_id) == step


def test_flags_on_hand_made_patients(report_cohort):
    for name, expected in FLAGGED.items():
        assert set(report_cohort.patients.index[report_cohort.patients[name]]) == expected, name


def test_a_tie_the_settings_leave_open_is_flagged():
    t = hand_made_tables(S, cases={"tie": ([cancer("tie")], [report("tie", 1, 50.0), report("tie", 1, 50.0, ca_seq=0)])})
    t.tests.loc[1, "test_id"] = "tie-T01b"                    # two tests, same day, panel and number
    c = build_cohort(outcome_blind(t), S)
    assert c.patients.loc["tie", "f:first_test_tied"] and c.patients.loc["tie", "first_test_id"] == "tie-T01"


# ---- Entry and exit: the adapter's outcome fields, and the screen's time-at-risk check ----
# One stage IV cancer and one report each. Unless given: age 60, died on day 1,000, report on day 100, BPC's after-death
# flags 0, and redacted "No". None means withheld or missing.
FOLLOW_CASES = {  # record_id: (cancer fields, report fields, redacted?, step excluded at, (vital status, days) if recorded)
    "plain":                    ({}, {}, False, None, (1, 1000)),
    "report_before_dx":         ({}, {"day": -20}, False, None, (1, 1000)),
    "report_on_dx_day":         ({}, {"day": 0}, False, None, (1, 1000)),
    "report_time_withheld":     ({}, {"day": None}, False, "entry_known", (1, 1000)),
    "age_withheld":             ({"age": None, "os_days": None, "os_status": None}, {}, True, "entry_known", None),
    "censored_at_89":           ({"age": 85, "os_days": None, "os_status": None}, {}, True, None, (0, 4 * 365.25)),
    "age_89_report_too_late":   ({"age": 89, "os_days": None, "os_status": None}, {"day": 282}, True,
                                 "entry_before_exit", (0, 0)),
    "survival_missing_not_redacted": ({"os_days": None, "os_status": None}, {}, False, "exit_known", None),
    "status_missing":           ({"os_status": None}, {}, False, "exit_known", None),
    "redacted_survival_known":  ({"age": 70}, {}, True, None, (1, 1000)),
    "report_after_death":       ({}, {"day": 1200, "post_death": 1, "post_last_alive": 1}, False, "entry_before_exit", (1, 1000)),
    "report_on_last_alive_day": ({"os_status": 0}, {"day": 1000, "post_last_alive": 1}, False, "entry_before_exit", (0, 1000)),
    "flag_alone_excludes":      ({}, {"day": 990, "post_death": 1}, False, "entry_before_exit", (1, 1000)),
    # Differs from the CRC project, by design: there the times decided these two, and the screen never reads
    # outcome times. The screen keeps them (the first is flagged as unknown); the run's time check removes both.
    "flags_missing_times_decide": ({}, {"day": 1200, "post_death": None, "post_last_alive": None}, False, None, (1, 1000)),
    "died_on_day_0":            ({"os_days": 0}, {"day": -5}, False, None, (1, 0)),
}


@pytest.fixture(scope="module")
def follow():
    cases = {rid: ([cancer(rid, **c)], [report(rid, 1, **{"day": 100.0, **r})]) for rid, (c, r, _, _, _) in FOLLOW_CASES.items()}
    return cohort_of(cases=cases, redacted={rid for rid, case in FOLLOW_CASES.items() if case[2]})


@pytest.mark.parametrize("record_id", list(FOLLOW_CASES))
def test_entry_and_exit_on_hand_made_patients(follow, record_id):
    c, t = follow
    *_, step, outcome = FOLLOW_CASES[record_id]
    assert excluded_at(c, record_id) == step
    p = t.patients.set_index("patient_id").loc[record_id]
    got = None if pd.isna(p["vital_status"]) else (int(p["vital_status"]), float(p["days_to_death_or_last_followup"]))
    assert got == (None if outcome is None else pytest.approx(outcome))


def test_entry_flags_on_hand_made_patients(follow):
    c, _ = follow
    p = c.patients
    assert set(p.index[p["f:censored_at_age_limit"]]) == {"censored_at_89", "age_89_report_too_late"}
    assert set(p.index[p["f:entry_before_diagnosis"]]) == {"report_before_dx", "died_on_day_0"}
    assert set(p.index[p["f:after_follow_up_end_unknown"]]) == {"flags_missing_times_decide"}


def test_left_blank_instead_of_censored_when_the_settings_say_so(tmp_path):
    from helpers import settings_with
    s = settings_with(tmp_path, {"source.genie_bpc.withheld_survival": "leave_blank"})
    t = hand_made_tables(s, cases={"c89": ([cancer("c89", age=85, os_days=None, os_status=None)], [report("c89", 1, 100.0)])},
                         redacted={"c89"})
    c = build_cohort(outcome_blind(t), s)
    assert t.patients["vital_status"].isna().all() and c.patients.loc["c89", "excluded_at"] == "exit_known"


# ---- The genotype groups ----
ONE_P8 = [(1, 100.0, "P8")]
MISSENSE, DEL, INS = "Missense_Mutation", "In_Frame_Del", "In_Frame_Ins"
GENOTYPE_CASES = {  # record_id: (reports, calls, genotype expected; None = not assignable, so excluded)
    "no_calls":             (ONE_P8, [], "Wild type"),
    "kras_g12d":            (ONE_P8, [(1, "KRAS", MISSENSE, "p.G12D")], "RAS-mutant"),
    "kras_g13d":            (ONE_P8, [(1, "KRAS", MISSENSE, "p.G13D")], "RAS-mutant"),
    "kras_a59t":            (ONE_P8, [(1, "KRAS", MISSENSE, "p.A59T")], "RAS-mutant"),
    "nras_q61k":            (ONE_P8, [(1, "NRAS", MISSENSE, "p.Q61K")], "RAS-mutant"),
    "kras_k117n":           (ONE_P8, [(1, "KRAS", MISSENSE, "p.K117N")], "RAS-mutant"),
    "kras_a146t":           (ONE_P8, [(1, "KRAS", MISSENSE, "p.A146T")], "RAS-mutant"),
    "kras_codon_19":        (ONE_P8, [(1, "KRAS", MISSENSE, "p.L19F")], "Wild type"),       # not a listed codon
    "kras_codon_135":       (ONE_P8, [(1, "KRAS", MISSENSE, "p.R135T")], "Wild type"),      # "135" contains "13"
    "kras_silent_codon_12": (ONE_P8, [(1, "KRAS", "Silent", "p.G12=")], "Wild type"),       # class not counted
    "deletion_spans_12":    (ONE_P8, [(1, "KRAS", DEL, "p.A11_G12del")], "RAS-mutant"),     # starts at 11, ends at 12
    "insertion_12_13":      (ONE_P8, [(1, "KRAS", INS, "p.G12_G13insA")], "RAS-mutant"),
    "deletion_next_to_61":  (ONE_P8, [(1, "KRAS", DEL, "p.E62_E63del")], "Wild type"),      # 62 to 63: misses 61
    "braf_v600e":           (ONE_P8, [(1, "BRAF", MISSENSE, "p.V600E")], "BRAF V600E"),
    "braf_v600e_space":     (ONE_P8, [(1, "BRAF", MISSENSE, "p.V600E ")], "BRAF V600E"),    # trailing space
    "braf_other":           (ONE_P8, [(1, "BRAF", MISSENSE, "p.K601E")], "Wild type"),      # stays wild type
    "co_mutant":            (ONE_P8, [(1, "KRAS", MISSENSE, "p.G12D"), (1, "BRAF", MISSENSE, "p.V600E")], "RAS-mutant"),
    "no_nras_on_panel":     ([(1, 100.0, "NO_NRAS")], [], None),
    "no_nras_kras_found":   ([(1, 100.0, "NO_NRAS")], [(1, "KRAS", MISSENSE, "p.G12D")], None),  # the panel decides
    "no_braf_on_panel":     ([(1, 100.0, "NO_BRAF")], [], None),
    "panel_without_file":   ([(1, 100.0, "NO-FILE")], [], None),
    "kras_on_later_sample": ([(1, 100.0, "P8"), (2, 400.0, "P8")], [(2, "KRAS", MISSENSE, "p.G12D")], "Wild type"),
}


@pytest.fixture(scope="module")
def genotype_cohort():
    cases = {rid: ([cancer(rid)], [report(rid, n, day, panel) for n, day, panel in reports])
             for rid, (reports, _, _) in GENOTYPE_CASES.items()}
    calls = [(f"{rid}-T{n:02d}", *call) for rid, (_, rid_calls, _) in GENOTYPE_CASES.items() for n, *call in rid_calls]
    return cohort_of(cases=cases, calls=calls)[0]


@pytest.mark.parametrize("record_id", list(GENOTYPE_CASES))
def test_genotype_on_hand_made_patients(genotype_cohort, record_id):
    expected = GENOTYPE_CASES[record_id][2]
    p = genotype_cohort.patients.loc[record_id]
    assert (None if pd.isna(p["genotype"]) else p["genotype"]) == expected
    assert p["excluded_at"] == (None if expected else "genotype_assignable")


def test_co_mutants_match_both_groups_rules(genotype_cohort):
    p = genotype_cohort.patients
    both = p["match:RAS-mutant"] & p["match:BRAF V600E"]
    assert set(p.index[both]) == {"co_mutant"}


def test_codon_spans():
    start, end = codon_span(pd.Series(["p.G12D", "p.G12_G13del", "p.E62_E63del", "p.R1450*", "c.35G>A", None],
                                      dtype="string"))
    assert start.tolist()[:4] == [12, 12, 62, 1450] and end.tolist()[:4] == [12, 13, 63, 1450]
    assert pd.isna(start.iloc[4]) and pd.isna(start.iloc[5])


def test_reference_only_coverage_assigns_a_variant_found_on_a_partial_panel(tmp_path):
    from helpers import settings_with
    s = settings_with(tmp_path, {"genotype.coverage.applies_to": "reference_only"})
    cases = {rid: ([cancer(rid)], [report(rid, 1, 100.0, "NO_NRAS")]) for rid in ("kras", "nothing")}
    t = hand_made_tables(s, cases=cases, calls=[("kras-T01", "KRAS", MISSENSE, "p.G12D")])
    p = build_cohort(outcome_blind(t), s).patients
    assert p.loc["kras", "genotype"] == "RAS-mutant" and pd.isna(p.loc["nothing", "genotype"])


# ---- Conditions ----
@pytest.mark.parametrize("cond, expected", [
    ({"field": "x", "op": "in", "values": ["A"], "normalize": "trim", "missing": "fail"}, [True, True, False, False]),
    ({"field": "x", "op": "not_in", "values": ["A"], "normalize": "trim", "missing": "pass"}, [False, False, True, True]),
    ({"field": "x", "op": "in", "values": ["A"], "missing": "fail"}, [True, False, False, False]),
    ({"field": "x", "op": "known"}, [True, True, True, False]),
])
def test_field_conditions_say_what_a_blank_means(cond, expected):
    df = pd.DataFrame({"x": pd.array(["A", "A ", "B", None], dtype="string")})
    assert field_condition(df, cond).tolist() == expected


# ---- The rules on the made-up fixture ----
@pytest.fixture(scope="module")
def fixture_cohort(fixture_result):
    return build_cohort(outcome_blind(fixture_result.tables), S)


TRUTH = pd.read_csv(FIXTURE / "truth.csv")


def test_one_first_test_per_patient(fixture_cohort):
    assert fixture_cohort.first["test_id"].notna().sum() == 1500


def test_first_test_day_matches_the_truth(fixture_cohort):
    """The generator knows each patient's real first report day. Where another report's time is withheld, the
    patient is flagged, and the report chosen is never earlier than the truth."""
    p = fixture_cohort.patients.join(TRUTH.set_index("record_id"))
    p = p[p["entry_day"].notna()]
    flagged = p["f:unranked_first_test"]
    assert (p.loc[~flagged, "entry_day"] == p.loc[~flagged, "true_first_report_day"]).all()
    assert (p.loc[flagged, "entry_day"] >= p.loc[flagged, "true_first_report_day"]).all()
    assert (p.loc[flagged, "entry_day"] > p.loc[flagged, "true_first_report_day"]).any()


def test_genotype_matches_the_truth(fixture_cohort):
    p = fixture_cohort.patients.join(TRUTH.set_index("record_id"))
    p = p[p["genotype"].notna()]
    assert len(p) and (p["genotype"] == p["true_group"]).all()


def test_the_eligible_count_and_flags_on_the_fixture(fixture_cohort):
    """The same counts the CRC project's tests found on this made-up fixture."""
    e = fixture_cohort.eligible
    assert len(e) == 658
    assert (e["match:RAS-mutant"] & e["match:BRAF V600E"]).sum() == 4
    assert e["f:censored_at_age_limit"].sum() == 3 and e["f:entry_before_diagnosis"].sum() == 6
    assert e["f:unranked_first_test"].sum() == 1 and e["f:first_test_linked_to_other_cancers"].sum() == 0


def test_samples_with_no_calls_are_kept_as_wild_type(fixture_result, fixture_cohort):
    """The inner-join trap: a first sample with no calls at all must stay in, as wild type."""
    called = set(fixture_result.tables.results["test_id"])
    e = fixture_cohort.eligible
    none = e[~e["first_test_id"].isin(called)]
    assert len(none) == 14 and set(none["genotype"]) == {"Wild type"}


def pandas_first_reports():
    ca = pd.read_csv(FIXTURE / "cancer_level_dataset_index.csv")
    ca = ca[ca["redcap_ca_index"].str.strip() == "Yes"].sort_values("ca_seq").drop_duplicates("record_id")
    rep = pd.read_csv(FIXTURE / "cancer_panel_test_level_dataset.csv").merge(ca[["record_id", "ca_seq"]])
    sizes = {name: len(set(genes)) for name, genes in map(raw.read_panel_file, FIXTURE.glob(raw.PANEL_FILES))}
    rep["n_genes"] = rep["cpt_seq_assay_id"].map(sizes)
    return (rep.sort_values(["record_id", "dx_cpt_rep_days", "n_genes", "cpt_number"],
                            ascending=[True, True, False, True], na_position="last").drop_duplicates("record_id"))


def test_first_test_matches_a_separate_choice_in_pandas(fixture_cohort):
    first = pandas_first_reports()
    mine = fixture_cohort.patients["first_test_id"]
    assert dict(zip(first["record_id"], first["cpt_genie_sample_id"])) == {k: v for k, v in mine.items() if pd.notna(v)}


def test_flow_adds_up(fixture_cohort):
    f = fixture_cohort.flow()
    assert f.loc[0, "remaining"] == 1500
    assert f["excluded"].sum() + f["remaining"].iloc[-1] == 1500
    assert (f["remaining"].diff().dropna() <= 0).all()
    assert f["remaining"].iloc[-1] == len(fixture_cohort.eligible)


def test_the_flow_has_one_step_per_rule_in_the_settings_order(fixture_cohort):
    f = fixture_cohort.flow(before=[("Has a BPC Project cancer", 0)])
    assert f["step"].tolist() == ["start", "source"] + [r["id"] for r in S.read("cohort.rules")]


# ---- Publishing the flow safely (privacy.public_flow, with the minimum count passed in) ----
def toy_flow(start, excluded):
    rows = [{"label": "start", "excluded": 0, "remaining": start}]
    left = start
    for i, n in enumerate(excluded, start=1):
        left -= n
        rows.append({"label": f"step {i}", "excluded": n, "remaining": left})
    return pd.DataFrame(rows)


def nothing_small_shown_or_recoverable(out, m=11):
    shown = list(out["excluded"]) + list(out["remaining"]) + list(-out["remaining"].diff().dropna())
    return not any(is_small(int(n), m) for n in shown)


def test_small_steps_merge_with_the_next_step():
    out = public_flow(toy_flow(1500, [0, 740, 5, 0, 10]), min_count=11)
    assert out["excluded"].tolist() == [0, 0, 740, 15]
    assert out["label"].iloc[-1] == "step 3; step 4; step 5"
    assert nothing_small_shown_or_recoverable(out)


def test_a_small_last_step_merges_into_the_one_before():
    out = public_flow(toy_flow(500, [100, 3]), min_count=11)
    assert out["excluded"].tolist() == [0, 103]


def test_zero_steps_stay_separate():
    assert public_flow(toy_flow(500, [0, 50, 0]), min_count=11)["excluded"].tolist() == [0, 0, 50, 0]


def test_a_flow_that_cannot_be_made_safe_is_refused():
    with pytest.raises(ValueError):
        public_flow(toy_flow(20, [3, 4]), min_count=11)
    with pytest.raises(ValueError):
        public_flow(toy_flow(100, [95]), min_count=11)


def test_the_minimum_count_comes_from_the_caller():
    assert public_flow(toy_flow(500, [100, 15]), min_count=11)["excluded"].tolist() == [0, 100, 15]
    assert public_flow(toy_flow(500, [100, 15]), min_count=20)["excluded"].tolist() == [0, 115]


def test_the_fixture_flow_is_safe_to_publish(fixture_cohort):
    f = fixture_cohort.flow()[["label", "excluded", "remaining"]]
    assert nothing_small_shown_or_recoverable(public_flow(f, min_count=11))
