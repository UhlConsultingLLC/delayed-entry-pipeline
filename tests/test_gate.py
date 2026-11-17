"""The safety gate: the five privacy protections and the security checks. Made-up numbers only.

Ported from the CRC BPC project's tests of crcbpc.outputs, with the minimum count passed in."""
import random

import pandas as pd
import pytest

from delayed_entry.gate import (NOT_SHOWN, SIZE_NOTE, CountSpec, Gate, GateError, Relation, balance_precision,
                                binary_smd, check_outputs, count_values, disclosed, drawn_steps, group_key,
                                hide_across, hide_survival_of_hidden_groups, pinned_small_cells, plain_ids,
                                public_curves, publish_counts, publish_descriptives, publish_risk, risk_key,
                                small_group_rule, withheld_curves)

G = Gate(11)
NS = NOT_SHOWN
ALL = "All patients"
GROUPS = ["Ref", "Group A", "Group B"]
YEARS = [1.0, 2.0, 3.0, 5.0]


# ---- 2. Counts found by subtraction ----
def test_a_lone_hidden_count_is_found_by_subtraction():
    values = {"T": 30, "A": 24, "B": 6}
    assert disclosed(values, {"T", "A"}, [Relation("T", ("A", "B"))], G) == ["B"]
    assert disclosed(values, {"A"}, [Relation("T", ("A", "B"))], G) == []


def test_two_hidden_counts_are_safe_only_if_their_sum_is_not_small():
    rel = [Relation("T", ("A", "B", "C"))]
    assert disclosed({"T": 40, "A": 27, "B": 6, "C": 7}, {"T", "A"}, rel, G) == []
    assert disclosed({"T": 34, "A": 27, "B": 3, "C": 4}, {"T", "A"}, rel, G) == ["T less its parts"]


def test_an_open_relation_gives_away_its_remainder():
    rel = [Relation("naive", ("corrected",), closed=False)]
    assert disclosed({"naive": 25, "corrected": 22}, {"naive", "corrected"}, rel, G) == ["naive less its parts"]
    assert disclosed({"naive": 104, "corrected": 104}, {"naive", "corrected"}, rel, G) == []


def test_subtraction_can_take_several_steps():
    rel = [Relation("B", ("C", "D")), Relation("T", ("A", "B"))]
    assert disclosed({"T": 100, "A": 60, "B": 40, "C": 35, "D": 5}, {"T", "A", "C"}, rel, G) == ["D"]


def test_the_minimum_comes_from_the_caller():
    values = {"T": 40, "A": 25, "B": 15}
    assert disclosed(values, {"T", "A"}, [Relation("T", ("A", "B"))], G) == []
    assert disclosed(values, {"T", "A"}, [Relation("T", ("A", "B"))], Gate(20)) == ["B"]


# ---- counts.csv ----
def spec(**changes):
    base = dict(eligible=500, primary={"Ref": 200, "Group A": 260, "Group B": 40},
                biomarkers={"B1": {"groups": {"Neg": 300, "Pos": 40}, "unknown": ("Unknown", 140),
                                   "not_in": ("not in B1", 20)}},
                others=[("SUB", "primary", {"Ref": 170, "Group A": 230, "Group B": 33})], first_groups={})
    base.update(changes)
    return CountSpec(**base)


def test_a_subset_count_close_to_its_base_is_hidden():
    s = spec()
    values, rel, _ = count_values(s)
    shown = publish_counts(values, rel, s, G)
    assert "SUB | Group B" not in shown                       # 7 fewer than the primary's
    assert "SUB | Ref" in shown and "SUB | Group A" in shown


def test_the_unknown_and_set_aside_groups_merge_when_either_is_small():
    s = spec(biomarkers={"B1": {"groups": {"Neg": 300, "Pos": 40}, "unknown": ("Unknown", 155), "not_in": ("not in B1", 5)}})
    values, rel, _ = count_values(s)
    shown = publish_counts(values, rel, s, G)
    assert "B1 | Unknown" not in shown and "B1 | not in B1" not in shown and "B1 | Unknown, or not in B1" in shown


def test_when_the_unknown_and_set_aside_together_are_small_the_largest_group_is_hidden_too():
    """With no unknown patients and 5 set aside, the groups alone leave 5 over from the analysis's total."""
    s = spec(biomarkers={"B1": {"groups": {"Neg": 300, "Pos": 195}, "unknown": ("Unknown", 0), "not_in": ("not in B1", 5)}})
    values, rel, _ = count_values(s)
    shown = {k for k, n in values.items() if not G.small(n)} - {"B1 | Unknown"}
    assert disclosed(values, shown, rel, G)                   # what the gate would give away without the rule
    shown = publish_counts(values, rel, s, G)
    assert "B1 | Neg" not in shown and "B1 | Pos" in shown and "B1 | Unknown, or not in B1" not in shown
    assert not disclosed(values, shown, rel, G)


def test_a_count_hidden_in_one_table_cant_be_worked_out_with_another():
    """A subset keeps the patients with one value of a field. The field's distribution shows that value's count,
    and the subset's groups add up to it, so subtraction would find the hidden group and with it a small count."""
    cv = {"primary | Ref": 200, "primary | Group A": 260, "primary | Group B": 40,
          "SUB | Ref": 170, "SUB | Group A": 230, "SUB | Group B": 33}
    cr = [Relation(f"primary | {g}", (f"SUB | {g}",), closed=False) for g in GROUPS]
    dv = {"Field | value": 433}
    shown = (set(cv) - {"SUB | Group B"}) | {"Field | value"}
    cross = [Relation("Field | value", tuple(f"SUB | {g}" for g in GROUPS))]
    values = cv | dv
    assert disclosed(values, shown, cr + cross, G)
    safe = hide_across(values, shown, cr, cross, G)
    assert safe == shown - {"SUB | Ref"} and not disclosed(values, safe, cr + cross, G)


def test_descriptives_hide_the_next_smallest_count_when_one_is_small():
    values = {"eligible": 200, "Field | x": 180, "Field | Missing": 4, "Field | y": 16, "Kept but counted | flag": 195}
    rel = [Relation("eligible", ("Field | x", "Field | Missing", "Field | y")),
           Relation("eligible", ("Kept but counted | flag",), closed=False)]
    shown = publish_descriptives(values, rel, G)
    assert "Field | Missing" not in shown and "Field | y" not in shown and "Field | x" in shown
    assert "Kept but counted | flag" not in shown               # 5 without the flag


# ---- The risk table ----
def made_up_risk_values(sizes=True):
    values = {}
    for g, n, naive, corrected in (("Ref", 300, [260, 220, 170, 60], [150, 150, 120, 40]),
                                   ("Group A", 330, [290, 190, 110, 30], [170, 140, 80, 20]),
                                   ("Group B", 40, [31, 15, 5, 0], [18, 12, 4, 0])):
        if sizes:
            values[group_key(g)] = n
        for y, a, c in zip(YEARS, naive, corrected):
            values[risk_key("naive", g, y)], values[risk_key("corrected", g, y)] = a, c
    if sizes:
        values[group_key(ALL)] = sum(values[group_key(g)] for g in GROUPS)
    for p in ("naive", "corrected"):
        for y in YEARS:
            values[risk_key(p, ALL, y)] = sum(values[risk_key(p, g, y)] for g in GROUPS)
    return values


def test_a_number_at_risk_a_small_number_below_the_group_size_is_hidden():
    b1, a1 = risk_key("naive", "Group B", 1.0), risk_key("naive", ALL, 1.0)
    assert b1 in publish_risk(made_up_risk_values(sizes=False), ALL, GROUPS, YEARS, G)
    shown = publish_risk(made_up_risk_values(), ALL, GROUPS, YEARS, G)
    assert b1 not in shown and a1 not in shown and risk_key("naive", "Ref", 1.0) in shown


def test_the_risk_rules_leave_nothing_to_subtract_in_random_tables():
    rng = random.Random(7)
    drop = lambda: rng.choice([0, rng.randint(1, 10), rng.randint(11, 60)])  # noqa: E731
    for _ in range(300):
        values = {}
        for g in GROUPS:
            n = naive = rng.randint(11, 300)
            values[group_key(g)] = n
            for y in YEARS:
                naive = max(0, naive - drop())
                values[risk_key("naive", g, y)] = naive
                values[risk_key("corrected", g, y)] = max(0, naive - drop())
        values[group_key(ALL)] = sum(values[group_key(g)] for g in GROUPS)
        for p in ("naive", "corrected"):
            for y in YEARS:
                values[risk_key(p, ALL, y)] = sum(values[risk_key(p, g, y)] for g in GROUPS)
        shown = publish_risk(values, ALL, GROUPS, YEARS, G)
        assert not any(G.small(values[k]) for k in shown)


# ---- 3. The curves ----
DPY = 365.25


def made_up_steps(naive_n, corrected_n, group="Group B", days=(100.0, 400.0, 600.0)):
    return pd.DataFrame([{"pass": p, "group": group, "day": d, "surv": 0.9 - 0.1 * i, "n_risk": n, "n_event": 1}
                         for p, ns in (("naive", naive_n), ("corrected", corrected_n))
                         for i, (d, n) in enumerate(zip(days, ns))])


def made_up_risk(hidden=(), group="Group B"):
    return pd.DataFrame([{"pass": p, "group": group, "years": y, "at_risk": NS if (p, y) in hidden else 50}
                         for p in ("naive", "corrected") for y in [1, 2, 3, 4, 5]])


def test_curves_are_drawn_only_where_enough_are_at_risk():
    d = pd.DataFrame({"pass": "corrected", "group": "Ref", "day": [10, 20, 30, 40, 50, 60, 70, 2000],
                      "surv": [0.9, 0.9, 0.8, 0.7, 0.7, 0.6, 0.5, 0.1], "n_risk": [5, 12, 11, 15, 13, 12, 10, 40],
                      "n_event": [1, 0, 1, 1, 0, 1, 1, 1]})
    steps = drawn_steps(d, {"Ref": 100}, 5 * DPY, G)
    assert list(steps["day"]) == [30, 40, 60]
    c = public_curves(steps, {}, DPY / 12)
    assert list(c["survival"]) == [0.8, 0.7, 0.6]


def test_curves_drawn_past_a_hidden_number_at_risk_are_withheld():
    steps = drawn_steps(made_up_steps([45, 35, 25], [30, 18, 12]), {"Group B": 50}, 5 * DPY, G)
    years = [1, 2, 3, 4, 5]
    assert withheld_curves(steps, made_up_risk(), years, DPY, G) == {}
    why = withheld_curves(steps, made_up_risk(hidden={("naive", 1)}), years, DPY, G)
    assert list(why) == ["Group B"] and "naive number at risk at 1 year," in why["Group B"]
    assert withheld_curves(steps, made_up_risk(hidden={("naive", 2)}), years, DPY, G) == {}   # drawn only to day 600


@pytest.mark.parametrize("naive_n, withheld", [([45, 35, 25], False), ([45, 35, 20], True), ([30, 18, 12], False)])
def test_curves_that_give_away_a_small_count_awaiting_entry_are_withheld(naive_n, withheld):
    steps = drawn_steps(made_up_steps(naive_n, [30, 18, 12]), {"Group B": 50}, 5 * DPY, G)
    assert ("Group B" in withheld_curves(steps, made_up_risk(), [1, 2, 3, 4, 5], DPY, G)) == withheld


# ---- 1 and 4. Estimates from few patients, and survival estimates that pin a hidden size ----
def test_a_group_of_few_patients_gets_no_estimate():
    t = pd.DataFrame({"analysis": ["primary"] * 2, "group": ["Group B", "Group B vs Ref"], "estimate": [13.6, 2.8],
                      "lower": [10.0, 1.8], "upper": [16.0, 4.3], "note": ["", ""]})
    out = small_group_rule(t, [8, 8], G)
    assert out["estimate"].isna().all() and (out["note"] == "fewer than 11 patients").all()


def test_a_survival_estimate_can_pin_a_hidden_group_size():
    """With nobody censored before it, survival is survivors / size: a size known to within 1 to 10 of 38 and an
    estimate of 0.774194 leave one size that fits."""
    assert [m for m in range(28, 38) if abs(round(0.774194 * m) / m - 0.774194) < 1e-6] == [31]


def test_survival_estimates_of_a_hidden_size_are_left_out():
    t = pd.DataFrame({"analysis": ["S2"] * 3, "group": ["Group B", "Group B", "Group A"],
                      "quantity": ["survival at 1 years", "median, months", "survival at 1 years"],
                      "estimate": [0.774194, 20.1, 0.8], "lower": None, "upper": None, "note": ["", "", ""]})
    out = hide_survival_of_hidden_groups(t, {("S2", "Group B")})
    assert pd.isna(out["estimate"][0]) and out["note"][0] == SIZE_NOTE
    assert out["estimate"][1] == 20.1 and out["estimate"][2] == 0.8


# ---- 5. Balance statistics ----
def test_the_smd_is_the_spread_of_shares_over_the_pooled_sd():
    assert binary_smd([10, 10], [100, 100]) == 0.0
    p = [0.2, 0.5]
    assert binary_smd([20, 50], [100, 100]) == pytest.approx((0.5 - 0.2) / ((0.2 * 0.8 + 0.5 * 0.5) / 2) ** 0.5)


def test_precise_smds_can_pin_a_small_cell_and_rounder_ones_cant():
    """Made-up: a trait with two levels, published per-level totals and group sizes. At many decimals the SMD fits
    only the true table, which holds a small cell; at fewer it fits many tables."""
    table = pd.DataFrame({"x": [40, 30, 3], "y": [160, 170, 37]}, index=GROUPS)
    smds = {lv: binary_smd(list(table[lv]), list(table.sum(axis=1))) for lv in table.columns}
    assert pinned_small_cells(table, smds, 6, G)
    d = balance_precision({"trait": (table, smds)}, 6, G)
    assert d < 6 and not pinned_small_cells(table, smds, d, G)


def test_a_trait_with_no_small_cell_keeps_every_decimal():
    table = pd.DataFrame({"x": [40, 30, 20], "y": [160, 170, 180]}, index=GROUPS)
    smds = {lv: binary_smd(list(table[lv]), list(table.sum(axis=1))) for lv in table.columns}
    assert balance_precision({"trait": (table, smds)}, 3, G) == 3


# ---- The security checks ----
def write(folder, name, text):
    (folder / name).write_text(text, encoding="utf-8")


@pytest.mark.parametrize("name, text, problem", [
    ("a.csv", "patient_id,months\nX,1\n", "row-level column"),
    ("b.csv", "group,months\nPT-0042,1\n", "an ID from the data"),
    ("c.csv", "group,patients\nRef,5\n", "patients shows a count below 11"),
    ("d.csv", "label,excluded,remaining\nStart,0,658\nStep,7,651\n", "excluded shows"),
    ("e.parquet", "", "not a table, figure or report"),
    ("f.md", "Diagnosed 2015-03-04.\n", "a date that isn't the run's"),
])
def test_the_security_scan_catches_what_must_not_be_published(tmp_path, name, text, problem):
    write(tmp_path, name, text)
    problems = check_outputs(tmp_path, G, {"PT-0042"}, set())
    assert any(problem in p for p in problems)


def test_the_security_scan_passes_safe_files_and_the_runs_own_dates(tmp_path):
    write(tmp_path, "ok.csv", "group,patients,n\nRef,296,0\nGroup B,<11,11\n")
    write(tmp_path, "report.md", "Run at 2026-01-02T03:04:05+00:00. PT-00421 is not an ID here.\n")
    assert check_outputs(tmp_path, G, {"PT-0042"}, {"2026-01-02"}) == []


def test_ids_that_are_plain_numbers_cant_be_scanned_for():
    keep, numeric = plain_ids(["PT-1", "17", "  ", None, "A2"])
    assert keep == {"PT-1", "A2"} and numeric == 1


def test_a_curve_drawn_past_a_hidden_number_at_risk_is_caught(tmp_path):
    write(tmp_path, "risk_table.csv", "pass,group,years,at_risk\nnaive,Group B,1,not shown\n")
    curve = "pass,group,months,survival\nnaive,Group B,0,1\nnaive,Group B,{},0.7\n"
    write(tmp_path, "curves.csv", curve.format(13.2))
    assert any("drawn past 1 year," in p for p in check_outputs(tmp_path, G, set(), set()))
    write(tmp_path, "curves.csv", curve.format(11.5))
    assert check_outputs(tmp_path, G, set(), set()) == []


def test_a_survival_estimate_for_a_hidden_size_is_caught(tmp_path):
    write(tmp_path, "counts.csv", "table,group,patients\nS2,Group B,not shown\nS2,Group A,300\n")
    est = "analysis,pass,quantity,group,estimate,lower,upper,note\nS2,naive,survival at 1 years,Group B,{},,,\n"
    write(tmp_path, "estimates.csv", est.format("0.774194"))
    assert any("hides the group's size" in p for p in check_outputs(tmp_path, G, set(), set()))
    write(tmp_path, "estimates.csv", est.format(""))
    assert check_outputs(tmp_path, G, set(), set()) == []


def test_a_gate_that_cant_be_made_safe_stops():
    with pytest.raises(GateError):
        from delayed_entry.gate import require_safe
        require_safe({"T": 30, "A": 24, "B": 6}, {"T", "A"}, [Relation("T", ("A", "B"))], "made up", G)
