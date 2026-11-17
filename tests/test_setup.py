"""The guided setup: the questionnaire, driven by scripted answers, and the spreadsheet form, filled in and read
back. The CRC settings serve as the answers: going through either form must give the same settings."""
import json
import re

import pandas as pd
import pytest
import yaml

from delayed_entry.__main__ import main
from delayed_entry.adapters.xlsx import read_xlsx, write_xlsx
from delayed_entry.settings import SCHEMA, SettingsError, load_settings, unfilled_choices
from delayed_entry.setup import (COLUMNS, GUIDE, SHEET, answers_from_file, form_to_settings, parse_value,
                                 questionnaire, read_form, write_form)
from helpers import CRC_SETTINGS

CRC = answers_from_file(CRC_SETTINGS)


def same_as_crc(path) -> None:
    """The settings file at path holds every CRC choice's value, made_by and source."""
    a, b = load_settings(CRC_SETTINGS), load_settings(path)
    for p in a.provenance():
        assert b.read(p["choice"]) == a.read(p["choice"]), p["choice"]
    got = {p["choice"]: (p["made_by"], p["source"]) for p in b.provenance()}
    assert got == {p["choice"]: (p["made_by"], p["source"]) for p in a.provenance()}


class Script:
    """A person at the questionnaire: answers each choice it is shown with the CRC value, unless `answers` says
    otherwise ({path: [answer, ...]}, given in turn); quits at the value prompt numbered quit_at."""

    def __init__(self, answers=None, quit_at=None):
        self.answers = {k: list(v) for k, v in (answers or {}).items()}
        self.quit_at, self.prompts, self.shown, self.said, self.path = quit_at, 0, [], [], None

    def say(self, text):
        self.said.append(text)
        m = re.search(r"\[\d+/\d+\] (\S+)", text)
        if m:
            self.path = m.group(1)
            self.shown.append(self.path)

    def ask(self, prompt):
        crc = CRC.get(self.path)
        if prompt.startswith("  Value"):
            self.prompts += 1
            if self.prompts == self.quit_at:
                return "q"
            if self.answers.get(self.path):
                return self.answers[self.path].pop(0)
            if crc is None:
                return "-"                                          # an optional choice the CRC file leaves out
            v = crc["value"]
            return v if isinstance(v, str) else json.dumps(v)
        if prompt.startswith("  Made by"):
            return crc["made_by"] if crc else "A. Reviewer"
        return (crc or {}).get("source") or "made up"


# ---- The questionnaire ----
def test_the_questionnaire_writes_the_crc_settings(tmp_path):
    s = Script()
    questionnaire(tmp_path / "s.yaml", s.ask, s.say)
    same_as_crc(tmp_path / "s.yaml")
    assert "The settings check finds no problem." in s.said[-1]


def test_it_can_stop_and_go_on_where_it_stopped(tmp_path):
    first = Script(quit_at=4)
    questionnaire(tmp_path / "s.yaml", first.ask, first.say)
    assert len(answers_from_file(tmp_path / "s.yaml")) == 3 and "Run setup on it again" in first.said[-1]
    second = Script()
    questionnaire(tmp_path / "s.yaml", second.ask, second.say)
    assert second.shown[0] == first.shown[3]                      # it asks only what is still blank
    same_as_crc(tmp_path / "s.yaml")


def test_an_invalid_value_is_asked_again_with_the_problem(tmp_path):
    s = Script({"time.horizon_years": ["-1", "3"]})
    questionnaire(tmp_path / "s.yaml", s.ask, s.say)
    assert any(t.startswith("  Not valid:") for t in s.said)
    assert load_settings(tmp_path / "s.yaml").read("time.horizon_years") == 3


def test_enter_takes_the_default_with_its_source(tmp_path):
    s = Script({"stability.risk_set": [""]})
    questionnaire(tmp_path / "s.yaml", s.ask, s.say)
    used = {d["choice"] for d in load_settings(tmp_path / "s.yaml").defaults_used()}
    assert used == {"stability.risk_set"}


def test_a_long_value_can_come_from_a_file(tmp_path):
    rules = tmp_path / "rules.yaml"
    rules.write_text(yaml.safe_dump(CRC["cohort.rules"]["value"], sort_keys=False), encoding="utf-8")
    s = Script({"cohort.rules": [f"@{rules}"]})
    questionnaire(tmp_path / "s.yaml", s.ask, s.say)
    same_as_crc(tmp_path / "s.yaml")


def test_a_choice_that_does_not_apply_is_not_asked(tmp_path):
    s = Script({"source.adapter": ["omop"]})
    questionnaire(tmp_path / "s.yaml", s.ask, s.say)
    assert not [p for p in s.shown if p.startswith("source.genie_bpc.")]
    assert "genie_bpc" not in yaml.safe_load((tmp_path / "s.yaml").read_text(encoding="utf-8"))["source"]


def test_an_answer_is_kept_as_text_where_only_the_text_is_valid():
    name = next(c for c in SCHEMA if c.path == "study.name")
    assert parse_value(name, "Yes") == "Yes" and parse_value(name, "2.0") == "2.0"
    horizon = next(c for c in SCHEMA if c.path == "time.horizon_years")
    assert parse_value(horizon, "3") == 3 and parse_value(horizon, "default") == "default"


# ---- The spreadsheet form ----
def test_the_form_has_a_row_for_every_choice_and_a_guide(tmp_path):
    path = write_form(tmp_path / "form.xlsx")
    df = read_xlsx(path, SHEET)
    assert list(df.columns) == COLUMNS and df["Choice"].tolist() == [c.path for c in SCHEMA]
    row = df.set_index("Choice").loc["stability.risk_set"]
    assert row["Default"] is not None and row["Default's source"] and pd.isna(row["Value"])
    assert row["What it controls"] == next(c.help for c in SCHEMA if c.path == "stability.risk_set")
    assert "setup-read" in " ".join(read_xlsx(path, GUIDE).iloc[:, 0].astype(str))


def test_the_crc_settings_come_back_through_the_form(tmp_path):
    problems = form_to_settings(write_form(tmp_path / "form.xlsx", CRC_SETTINGS), tmp_path / "s.yaml")
    assert problems == []
    same_as_crc(tmp_path / "s.yaml")


def edited(tmp_path, change) -> object:
    """The CRC form, with cells changed as a person would ({choice: {column: text}})."""
    path = write_form(tmp_path / "form.xlsx", CRC_SETTINGS)
    df = read_xlsx(path, SHEET).set_index("Choice", drop=False)
    for choice, cells in change.items():
        for col, text in cells.items():
            df.loc[choice, col] = text
    rows = [COLUMNS] + [[None if pd.isna(v) else v for v in r] for r in df[COLUMNS].values.tolist()]
    return write_xlsx(tmp_path / "edited.xlsx", {SHEET: rows})


def test_cells_a_person_edits_are_read_back(tmp_path):
    form = edited(tmp_path, {"time.horizon_years": {"Value": "4"}, "study.name": {"Value": "Yes"},
                             "stability.risk_set": {"Value": "default", "Source": None}})
    assert form_to_settings(form, tmp_path / "s.yaml") == []
    s = load_settings(tmp_path / "s.yaml")
    assert s.read("time.horizon_years") == 4 and s.read("study.name") == "Yes"
    assert [d["choice"] for d in s.defaults_used()] == ["stability.risk_set"]


def test_a_blank_form_names_every_choice_still_to_fill_in(tmp_path):
    problems = form_to_settings(write_form(tmp_path / "form.xlsx"), tmp_path / "s.yaml")
    with pytest.raises(SettingsError) as e:
        load_settings(tmp_path / "s.yaml")
    assert set(unfilled_choices(e.value)) == {c.path for c in SCHEMA} and len(problems) == len(SCHEMA)


def test_a_value_the_settings_refuse_is_named(tmp_path):
    form = edited(tmp_path, {"time.horizon_years": {"Value": "soon"}})
    problems = form_to_settings(form, tmp_path / "s.yaml")
    assert any("time.horizon_years" in p for p in problems)


@pytest.mark.parametrize("change, problem", [
    ({"study.name": {"Choice": "study.title"}}, "'study.title' is not a choice"),
    ({"study.cancer": {"Choice": "study.name"}}, "'study.name' appears twice"),
])
def test_a_form_whose_rows_are_wrong_is_refused(tmp_path, change, problem):
    with pytest.raises(SettingsError, match=problem):
        read_form(edited(tmp_path, change))


# ---- The command line ----
def test_the_commands_write_and_read_the_form(tmp_path, capsys):
    main(["setup-sheet", str(tmp_path / "form.xlsx"), "--from", str(CRC_SETTINGS)])
    main(["setup-read", str(tmp_path / "form.xlsx"), str(tmp_path / "s.yaml")])
    assert "finds no problem" in capsys.readouterr().out
    same_as_crc(tmp_path / "s.yaml")
    main(["setup-sheet", str(tmp_path / "blank.xlsx")])
    with pytest.raises(SystemExit):
        main(["setup-read", str(tmp_path / "blank.xlsx"), str(tmp_path / "b.yaml")])
