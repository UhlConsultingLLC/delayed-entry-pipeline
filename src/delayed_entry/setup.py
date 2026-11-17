"""The guided setup (design plan section 4): a study's settings without editing the settings file by hand, in two
forms that write the same file.

    python -m delayed_entry setup settings/my_study.yaml                  a questionnaire, one choice at a time
    python -m delayed_entry setup-sheet form.xlsx [--from SETTINGS]       a spreadsheet form, one row per choice
    python -m delayed_entry setup-read form.xlsx settings/my_study.yaml   a filled form back to a settings file

Both go through every choice in the settings schema's order, in its own words: what it controls, its form, an
example where the form is complex, and its default with the default's source. Both write the commented settings
file that `check-settings` reads, and say what that check finds. A value is only ever judged by the settings' own
rules.

A value is written as YAML: plain words, a number, [a, list], {a: mapping}, or several lines of YAML (in a sheet's
cell). Where the YAML reading of an answer isn't a valid value but the answer as plain text is (a study named "Yes",
a cohort named 2.0), the text is kept. A choice with a default takes the word default, and its source is then the
default's own. In the questionnaire, a long value (cohort rules, genotype groups, analyses) can be given as @FILE,
a YAML file holding just the value. made_by and source are asked for every choice; Enter repeats the last made_by.

The questionnaire checks each value as it is given, and saves the file after each answer, so it can be stopped
(q) and resumed: run it again on the same file and it asks only what is still blank. A choice that doesn't apply
(a GENIE BPC choice when the source is OMOP) isn't asked, and an optional choice can be passed over with -. The
checks that compare choices with each other run when the file is complete.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import yaml

from delayed_entry.adapters.xlsx import read_xlsx, write_xlsx
from delayed_entry.settings import (DEFAULT_WORD, SCHEMA, SCHEMA_BY_PATH, Choice, SettingsError, applies_to,
                                    load_defaults, loads_settings, settings_text)

SHEET, GUIDE = "Settings", "How to fill it in"
COLUMNS = ["Choice", "What it controls", "Form", "Only when", "Default", "Default's source", "Example",
           "Value", "Made by", "Source", "Note"]
WIDTHS = [32, 48, 32, 24, 22, 34, 40, 40, 20, 34, 24]
QUIT, PASS = "q", "-"
GUIDE_LINES = [
    "Fill in three cells on each row of the Settings sheet: Value, Made by (who made the choice) and Source (the "
    "source behind it). Note is optional. Leave the other columns as they are.",
    "Value: plain words or a number for a simple choice. For a list or a mapping, write YAML, as the Example "
    "column does: [a, b] or one item per line starting with - , and name: value lines for a mapping.",
    "Default: where a choice has one, writing the word default in Value uses it, with the default's own source. "
    "A choice with no default must be filled in.",
    "Only when: a choice that applies only when another choice has a given value. Leave it blank otherwise.",
    "Optional choices (Default says none (optional)) may be left blank.",
    "Then: python -m delayed_entry setup-read <this file> <settings file>. It writes the settings file and lists "
    "every problem the settings check finds, each with the choice it concerns. Fix them here and read it again, or "
    "in the settings file.",
]


def show_value(v) -> str:
    """A value as an answer would give it: text as it is, anything else as YAML."""
    if isinstance(v, str):
        return v
    text = yaml.safe_dump(v, sort_keys=False, allow_unicode=True, default_flow_style=None, width=10**6).strip()
    return text[:-3].strip() if text.endswith("...") else text


def parse_value(choice: Choice, text: str):
    """An answer as a value: its YAML reading, or the answer as plain text where only the text is valid."""
    text = text.strip()
    if text == DEFAULT_WORD:
        return DEFAULT_WORD
    try:
        parsed = yaml.safe_load(text)
    except yaml.YAMLError:
        return text
    if parsed is None or not choice.check(parsed):
        return parsed
    return text if not choice.check(text) else parsed


def answers_from_file(path: str | Path) -> dict[str, dict]:
    """The choices a settings file has filled in, as {path: {value, made_by, source, note}}."""
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    out = {}
    for c in SCHEMA:
        node = data
        for part in c.path.split("."):
            node = node.get(part) if isinstance(node, dict) else None
        if isinstance(node, dict) and node.get("value") is not None:
            out[c.path] = {k: node[k] for k in ("value", "made_by", "source", "note") if node.get(k) is not None}
    return out


def check(text: str, path: Path | None = None) -> list[str]:
    """What the settings check finds in a settings file's text: nothing if it loads."""
    try:
        loads_settings(text, path)
    except SettingsError as e:
        return str(e).splitlines()[1:]
    return []


def save(out: Path, answers: dict[str, dict]) -> str:
    text = settings_text(answers)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(text.encode("utf-8"))
    return text


# ---------------------------------------------------------------------------
# The questionnaire
# ---------------------------------------------------------------------------
def describe(c: Choice, number: int, defaults: dict) -> str:
    lines = ["", f"[{number}/{len(SCHEMA)}] {c.path}", f"  {c.help}", f"  Form: {c.form}."]
    if c.path in defaults:
        d = defaults[c.path]
        lines.append(f"  Default: {show_value(d['value'])} (source: {d['source']}). Enter takes it.")
    else:
        lines.append("  No default." + (" Optional: - passes it over." if c.optional else ""))
    if c.example is not None:
        lines.append("  Example:")
        lines += [f"    {x}" for x in show_value(c.example).splitlines()]
    return "\n".join(lines)


def questionnaire(out: str | Path, ask: Callable[[str], str] = input, say: Callable[[str], None] = print,
                  defaults: dict | None = None) -> dict[str, dict]:
    """Ask each choice still blank in `out` (a new or partly filled settings file), saving after each answer.
    ask(prompt) returns the person's answer; say(text) shows text. Returns the answers."""
    out = Path(out)
    defaults = load_defaults() if defaults is None else defaults
    answers = answers_from_file(out) if out.exists() else {}
    last_by = next((a["made_by"] for a in reversed(list(answers.values())) if a.get("made_by")), None)
    if answers:
        say(f"{out.name}: {len(answers)} choices already filled in; asking the rest.")
    for number, c in enumerate(SCHEMA, start=1):
        if c.path in answers or applies_to(c, answers) is not True:
            continue
        say(describe(c, number, defaults))
        value = _ask_value(c, ask, say, defaults)
        if value is QUIT:
            save(out, answers)
            say(f"Saved to {out}. Run setup on it again to go on.")
            return answers
        if value is PASS:
            continue
        by = ask("  Made by" + (f" (Enter: {last_by})" if last_by else "") + ": ").strip() or last_by
        while not by:
            by = ask("  Made by (who made this choice): ").strip()
        source = ask("  Source" + (" (Enter: the default's own)" if value == DEFAULT_WORD else "") + ": ").strip()
        while not source and value != DEFAULT_WORD:
            source = ask("  Source (the source behind the choice): ").strip()
        answers[c.path] = {"value": value, "made_by": by, "source": source or None}
        last_by = by
        save(out, answers)
    problems = check(save(out, answers), out)
    say(f"\nSaved to {out}. " + ("The settings check finds no problem." if not problems else
                                 "The settings check finds:\n" + "\n".join(problems)))
    return answers


def _ask_value(c: Choice, ask, say, defaults):
    while True:
        raw = ask("  Value: ").strip()
        if raw == QUIT:
            return QUIT
        if raw == PASS and c.optional:
            return PASS
        if not raw:
            if c.path in defaults:
                return DEFAULT_WORD
            say("  A value is needed (or q to stop here and go on later).")
            continue
        if raw.startswith("@"):
            try:
                value = yaml.safe_load(Path(raw[1:].strip()).read_text(encoding="utf-8"))
            except (OSError, yaml.YAMLError) as e:
                say(f"  Can't read {raw[1:].strip()}: {e}")
                continue
        else:
            value = parse_value(c, raw)
        if value == DEFAULT_WORD and c.path not in defaults:
            say("  This choice has no default: give a value.")
            continue
        problems = [] if value == DEFAULT_WORD else c.check(value)
        if not problems:
            return value
        say("  Not valid: " + "; ".join(problems))


# ---------------------------------------------------------------------------
# The spreadsheet form
# ---------------------------------------------------------------------------
def write_form(path: str | Path, from_file: str | Path | None = None, defaults: dict | None = None) -> Path:
    """The form: one row per choice, with its description, form, default and example, and the Value, Made by,
    Source and Note to fill in (already filled from `from_file`, a settings file, if one is given)."""
    defaults = load_defaults() if defaults is None else defaults
    answers = answers_from_file(from_file) if from_file else {}
    rows = [COLUMNS]
    for c in SCHEMA:
        d, a = defaults.get(c.path), answers.get(c.path, {})
        rows.append([c.path, c.help, c.form, f"{c.when[0]} is {c.when[1]}" if c.when else None,
                     show_value(d["value"]) if d else ("none (optional)" if c.optional else "none"),
                     d["source"] if d else None,
                     show_value(c.example) if c.example is not None else None,
                     show_value(a["value"]) if "value" in a else None,
                     a.get("made_by"), a.get("source"), a.get("note")])
    return write_xlsx(path, {SHEET: rows, GUIDE: [["How to fill in the settings form"]] + [[x] for x in GUIDE_LINES]},
                      widths={SHEET: WIDTHS, GUIDE: [120]}, wrapped=True)


def _cell(v):
    if v is None or (isinstance(v, float) and v != v):
        return None
    if isinstance(v, str):
        return v.strip() or None
    return v


def read_form(path: str | Path) -> dict[str, dict]:
    """A filled form's answers, as {path: {value, made_by, source, note}}. Rows whose Value is blank are left out."""
    df = read_xlsx(path, SHEET)
    missing = [c for c in ("Choice", "Value", "Made by", "Source") if c not in df.columns]
    if missing:
        raise SettingsError(f"{Path(path).name}: the {SHEET} sheet has no column {', '.join(missing)}")
    answers, problems = {}, []
    for i, r in df.iterrows():
        name = _cell(r["Choice"])
        if name is None:
            continue
        c = SCHEMA_BY_PATH.get(str(name))
        if c is None:
            problems.append(f"row {i + 2}: '{name}' is not a choice")
            continue
        if c.path in answers:
            problems.append(f"row {i + 2}: '{c.path}' appears twice")
            continue
        v = _cell(r["Value"])
        if v is None:
            continue
        value = parse_value(c, v) if isinstance(v, str) else v
        note = _cell(r["Note"]) if "Note" in df.columns else None
        answers[c.path] = {"value": value, "made_by": _cell(r["Made by"]), "source": _cell(r["Source"])}
        if note is not None:
            answers[c.path]["note"] = note
    if problems:
        raise SettingsError(f"{Path(path).name}: " + "; ".join(problems))
    return answers


def form_to_settings(form: str | Path, out: str | Path) -> list[str]:
    """Write the settings file a filled form gives, and return what the settings check finds in it."""
    out = Path(out)
    return check(save(out, read_form(form)), out)
