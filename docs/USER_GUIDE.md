# Delayed-entry pipeline: setting up, testing and running it

This guide grows with the pipeline. Each step is marked **available** or **not built yet**, with the phase
that will build it. The design is in [DESIGN_PLAN.md](../DESIGN_PLAN.md); [README.md](../README.md) gives
the overview.

| Stage | What it does | Status |
|---|---|---|
| 1. Read | a source into the standard tables | available (phase 1) |
| 2. Screen | what the data support, before any outcome is seen | available (phase 1) |
| 3. Lock | the reviewed plan, with a hash of everything it rests on | available (phase 1) |
| 4. Run | both clocks, the adjusted and weighted models, the checks | available (phase 2) |
| 5. Safety gate | privacy and security checks on every output | available (phase 2) |
| 6. Report | results and the run report | available (phase 2) |

Commands below are for Windows PowerShell. On macOS or Linux, use `/` in paths and
`export PYTHONPATH=src` instead of `$env:PYTHONPATH = "src"`.

---

## Part 1. Setting up

### 1.1 What you need

- **Python 3.10 or later**, with pandas, NumPy, DuckDB and PyYAML; pytest to run the tests.
- **R** (for the run, step 8), for the Cox models with robust standard errors and the bootstrap, with the
  packages survival, WeightIt, cobalt, tranSurv, duckdb and jsonlite. The CRC BPC study used R 4.6.2 with
  survival 3.8-12, WeightIt 2.1.0, cobalt 5.0.0 and tranSurv 1.2.4; the same versions reproduce its results
  exactly.
- **Data you are allowed to use**, under your own data-use agreement. The pipeline doesn't fetch data.
- **A local data folder** that is not inside a cloud-synced folder (OneDrive, Dropbox, Google Drive,
  iCloud). The pipeline refuses one that is. Keep it outside the pipeline's own folder too.
- **A work folder** for reports and locked plans, under the same rules.

### 1.2 Get the code

The pipeline is a local folder for now. It will be published as a public repository once phase 2's
acceptance test passes. Nothing in the folder holds patient data, and nothing should be added that does.

### 1.3 Make a Python environment

On a new machine, from the pipeline folder (this downloads the packages):

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install -e ".[test]"
```

`-e` installs the pipeline from its folder, so edits take effect without reinstalling.

Without installing, put the source folder on Python's path instead, in each new terminal:

```powershell
$env:PYTHONPATH = "src"
```

The tests don't need either: `pyproject.toml` tells pytest where the source is.

### 1.4 Check the setup

```powershell
python -m delayed_entry --help
```

It should list the commands: `template`, `check-settings`, `lock-settings`, `draft-mapping`, `screen`,
`lock`, `verify-lock` and `run`.

---

## Part 2. Testing

### 2.1 Run the tests

From the pipeline folder:

```powershell
python -m pytest -q -p no:cacheprovider
```

- `-q` prints a dot per passing test and a one-line summary.
- `-p no:cacheprovider` stops pytest leaving a cache folder behind.
- It takes about a minute. The last line should read `N passed`, with no failures.

To run one area only, give its file, for example `python -m pytest -q -p no:cacheprovider tests/test_screen.py`.

The tests that run R (a whole run on the made-up fixture, and the R unit checks in `tests/R/`) are skipped
unless R is named. To include them:

```powershell
$env:DELAYED_ENTRY_R = "<path to Rterm.exe or Rscript>"
$env:DELAYED_ENTRY_R_LIBRARY = "<folder of R packages>"     # only if they aren't in R's own library
python -m pytest -q -p no:cacheprovider
```

The last line then reads `N passed`, with nothing skipped.

### 2.2 What the tests use and check

Made-up data only: a made-up GENIE BPC release (`tests/fixtures/genie_bpc_crc_synthetic`), the same
patients written as a made-up OMOP database (`tests/omop_fixture.py`), and small tables built by hand. Never real data, never real counts. They check, among other things, that:

- scrambling every outcome value leaves the readiness report unchanged (the screen never reads outcomes);
- a threshold changed in the settings file changes the screen's result, with no code change;
- no clinical or statistical value is written into the code;
- no count below the privacy minimum can appear in a report, or be worked out from one, and the safety gate
  catches each of the five ways a small count can leak, and an ID, a date or a row-level column in an output;
- the run stops if the lock doesn't match, publishes nothing when the gate finds a problem, and takes every
  value R uses from the settings;
- the cohort, first test, genotype groups and MSI/MMR status match the CRC BPC project's results on the
  same made-up data.

A failing test prints its name and what differed. Don't go on to real data until every test passes.

---

## Part 3. Running a study

The order matters: the settings are written and locked before the data are downloaded, and the plan is
locked before any outcome is analyzed.

### Step 1. Write the study's settings (available)

Every clinical and statistical choice lives in one settings file, with who made it and the source behind
it. Start from the blank template, or copy a ready-made file such as `settings/crc_bpc.yaml`:

```powershell
python -m delayed_entry template settings\my_study.yaml
```

Fill in every choice: `value`, `made_by` and `source`. Where the template shows a default, writing
`value: default` uses it, and the plan will list it with its source. Choices specific to a study (cohort
rules, genotype groups, coverage, horizon, analyses, privacy minimum) have no default.

Then check the file. It lists every problem with its line number, or says it is OK:

```powershell
python -m delayed_entry check-settings settings\my_study.yaml
```

The privacy minimum (`privacy.min_count`) can be raised but never set below the data-use agreement's.
For a new cancer, have a clinician of that specialty review the file before it is locked.

### Step 2. Lock the settings before downloading (available)

```powershell
python -m delayed_entry lock-settings settings\my_study.yaml <work folder>\settings_lock --locked-by "<reviewer>"
```

This records the file's SHA-256 hash, who locked it and when. It refuses a file that doesn't load, and
won't overwrite an existing lock: a locked choice changes only by a new plan version.

### Step 3. Get the data

Download under your own account and data-use agreement, into the local data folder. Never put
credentials in any file in the pipeline folder.

**GENIE BPC (available).** Access is granted on Synapse: one access requirement covers every public BPC
release, so approval for one cohort covers the others. The download script needs the Synapse client
(`pip install -e ".[synapse]"`) and reads your Synapse personal access token itself, from the
`SYNAPSE_AUTH_TOKEN` environment variable or `~/.synapseConfig`; the token is never typed into a command.

The documentation (data guide, variable synopsis) can be downloaded at any time, and is needed to write
the settings:

```powershell
python scripts\download_synapse.py nsclc_v2.0 documentation <data folder>
```

The clinical and genomic files hold patient-level data. The script downloads them only with the study's
settings and the lock from step 2, and only if the settings name this release (`source.genie_bpc.cohort`)
and haven't changed since the lock:

```powershell
python scripts\download_synapse.py nsclc_v2.0 clinical <data folder> --settings settings\my_study.yaml --settings-lock <work folder>\settings_lock\settings_lock.json
python scripts\download_synapse.py nsclc_v2.0 genomic <data folder> --settings settings\my_study.yaml --settings-lock <work folder>\settings_lock\settings_lock.json
```

Each part goes into `<data folder>\<part>`, with a checksum manifest comparing every file with Synapse;
it ends `checksum mismatches: none` when the download is complete. Record each download (date, release,
Synapse ID, manifest name, and the date you read the terms) in the study's provenance log, for example
`docs/provenance/<cohort>.md`. The releases and their Synapse folders are
listed in `src/delayed_entry/adapters/genie_bpc/synapse_releases.yml`: `crc_v2.0` and `nsclc_v2.0` so far.

**Other sources.** Download as the data holder directs; the pipeline reads the files through a mapping
(step 4).

### Step 4. Check the files

**GENIE BPC (available).** The check compares each downloaded file with the release list: its size and
checksum on Synapse. For a release whose list isn't in the package yet (NSCLC v2.0-public so far), make it
once, from Synapse's metadata, without downloading anything (it needs an approved account, like the
download):

```powershell
python scripts\synapse_release_list.py nsclc_v2.0
```

It writes `src\delayed_entry\adapters\genie_bpc\cohorts\nsclc_v2.0\release_files.yml`. Then, right after
downloading:

```powershell
python -m delayed_entry.adapters.genie_bpc.check_download settings\my_study.yaml <clinical folder> <genomic folder> <work folder> --settings-lock <work folder>\settings_lock\settings_lock.json
```

It checks the folders, each file against the release list (size and checksum), that every file reads and
every column the adapter uses is there with values of the right kind, the coded values, the units, the
links between files, and the labels the settings search for. It then builds the standard tables once.
Each check ends `pass`, `FAIL`, `review` or `skipped`. Any `FAIL` stops the run; read each `review` and note
it in the run log. The report is written to the work folder; counts below the minimum show as "<minimum".
Add `--rehearsal` to try it on made-up data.

**Flat files: CSV/TSV, Excel (.xlsx) or Parquet (available).** Draft a mapping from the files' column names:

```powershell
python -m delayed_entry draft-mapping <work folder>\mapping.yaml patients=patients.csv tests=tests.xlsx#Tests results=results.parquet panels=panels.csv
```

The draft is never used as it stands. A person must:
1. check every drafted field, and give each one a column, a constant or `absent: true`;
2. add any transform: `map` for coded values (quote words like `'Yes'`), `multiply` for units (months to
   days), `days_since_diagnosis` for dates (with the patients file's `diagnosis_date`);
3. fill in the capability profile: entry event, panel coverage, panel sizes, death capture, time format,
   the agreement's minimum count, the histology vocabulary, the age top code and known quirks;
4. set `status: confirmed`, `confirmed_by` and `confirmed_on`.

Set the settings' `source.adapter` to `flat_file`.

**An OMOP CDM database (available).** The adapter reads OMOP version 5.4 tables. Either export them to one
folder, one file per table named after it (`person.csv`, `death.parquet`, and so on; CSV, TSV or Parquet), or
give a DuckDB file holding them in one schema. It reads `person`, `condition_occurrence`,
`observation_period` and `procedure_occurrence`, plus `death`, `measurement`, `observation`, `care_site`,
`drug_exposure`, `concept` and `concept_ancestor` where the mapping uses them, and the genomic extension
table it names. Start from the blank mapping:

```powershell
python -m delayed_entry omop-template <work folder>\omop_mapping.yaml
```

A person fills in every choice. The template explains each one.
1. **Where the tables are:** `cdm: {folder: ..., date_format: "%Y-%m-%d"}` or
   `cdm: {duckdb: ..., schema: ...}`.
2. **The concept sets,** by name: the study's cancer, the sequencing procedures, each kind of reading (for
   example MSI and MMR status), the cancer drugs, and any concept a patient field reads. Each set lists
   concept IDs, and can take in their descendants through `concept_ancestor`.
3. **The patient ID** (`person_id` or `person_source_value`) and the diagnosis concept set. The earliest
   such condition is each person's diagnosis, and persons with none are left out at the flow diagram's
   first step.
4. **Each patient field's source:**
   - `{person: COLUMN}`, `{care_site: COLUMN}` or `{index_condition: COLUMN}`;
   - `{measurement: SET, value: COLUMN}` or `{observation: SET, value: COLUMN}`, which reads the row dated
     closest to diagnosis;
   - `{computed: true}`, for the age (the year of diagnosis less the year of birth) and the year of diagnosis;
   - or `{absent: true}`.

   Add `extract` (a regular expression with one group, for example the morphology code in
   `condition_source_value`) or `labels` (each source value and what it reads as) where needed. `labels`
   must list every value found: the read stops and names any it doesn't list. List a value as `null` to read
   it as blank. Extra fields the settings name, such as a tumor location, go under `extras`.
5. **The sequencing tests:**
   - the concept set of their procedures;
   - their modality, as the settings' `entry.modalities` name it;
   - what the procedure date is (`report`, `sample_collection` or `order`);
   - where the panel's name is;
   - the variant table and its columns: the link to the procedure, gene, protein change, and if it has
     them, the fusion partner, the variant class and the copy number.
6. **The readings** from the measurement table (MSI or MMR status, a protein), each with its labels, its
   modality and its entry event. Each reading becomes a test with one result.
7. **Optional files:** a panels file (`panel_id`, `gene`, and optionally `codon_from`, `codon_to`) and a
   panel sizes file (`panel_id`, `coding_megabases`). OMOP holds neither. Without a panels file, wild type
   can't be confirmed from panel coverage.
8. **The cancer drugs' concept set,** and whether drugs are named by `drug_source_value` or by the concept's
   name. Same-day exposures form one regimen, numbered by date.
9. **The capability profile.** It must agree with the rest:
   - `time_format` is `real_dates` or `shifted_dates`, and with shifted dates the year of diagnosis can't be
     computed;
   - `entry_event` is the sequencing tests';
   - the panel coverage is `none` exactly when there is no panels file, and likewise for the panel sizes.
10. Set `status: confirmed`, `confirmed_by` and `confirmed_on`.

Set the settings' `source.adapter` to `omop`, and give the mapping with `--mapping` at steps 5, 7 and 8. Where
`person_source_value` isn't a plain number, prefer it as the patient ID. The safety gate's scan can't look for
plain-number IDs, so it leaves them out and says so.

### Step 5. Screen (available)

```powershell
python -m delayed_entry screen settings\my_study.yaml <work folder> --clinical <clinical folder> --genomic <genomic folder>
```

For flat files or OMOP, use `--mapping <work folder>\mapping.yaml` instead of the two folders. It prints a short
summary and writes `readiness_report.md` and `readiness_report.yaml` to the work folder. It reads no
outcome: it sees only whether each patient's outcome is recorded.

### Step 6. Review the readiness report (a person's step)

Read it in full. In particular:
- **Ready, Missing, Deviations required:** each item says what was found and what it changes. A `block`
  means that analysis can't run as planned; a `deviation` means it runs differently, stated in the plan;
  a `note` is for information.
- **Decisions:** each comparison is compared, described only, not run, or blocked. A blocked comparison
  needs a person to decide what is done instead, logged as a deviation, before any outcome is seen.
- **Warnings:** choices that look implausible for the data. They never change a result; you decide.
- **Checks that need outcomes:** what the run will check after the lock (for example deaths soon after
  the entry event).
- **Settings the screen read** and **Defaults applied:** where each value came from.

If something must change, change the settings (or the mapping), check them, and screen again.

### Step 7. Lock the plan (available)

When the report has been reviewed:

```powershell
python -m delayed_entry lock settings\my_study.yaml <work folder> --clinical <clinical folder> --genomic <genomic folder> --locked-by "<reviewer>"
```

This screens again and writes `plan_v1\` into the work folder: the plan in Markdown and YAML, a copy of the
settings, and `manifest.json` with the SHA-256 hashes of the plan, the settings file, the defaults file and
every input file, and the time. It refuses if the privacy minimum is below the agreement's. A later change
needs a new version: add `--reason "<why>"`, and a line goes into `deviations_log.md`.

Registering the plan (for example on OSF) is a manual step; the Markdown plan is written for it.

To check that nothing has changed since:

```powershell
python -m delayed_entry verify-lock <work folder>\plan_v1\manifest.json
```

### Step 8. Run the analysis (available)

When the plan is locked (step 7):

```powershell
python -m delayed_entry run settings\my_study.yaml <work folder> <output folder> --clinical <clinical folder> --genomic <genomic folder> --r <R program> --r-library <R packages folder>
```

- For flat files or OMOP, use `--mapping <work folder>\mapping.yaml` instead of the two folders.
- `--r` names the R program: `Rterm.exe` (run as `Rterm.exe -f <script>`) or `Rscript`. Without it, the run
  uses the environment variable `DELAYED_ENTRY_R`, then `Rscript` on the PATH. On Windows machines where
  `Rscript.exe` is blocked by an application-control policy, give the full path to `bin\x64\Rterm.exe`.
- `--r-library` puts a folder of R packages first on R's library path, for packages installed outside R's own
  library (default: the environment variable `DELAYED_ENTRY_R_LIBRARY`).
- `--manifest <work folder>\plan_v<n>\manifest.json` picks a plan version; without it, the latest is used.
- The output folder must be empty (or hold an earlier run's outputs), outside the work folder and outside the
  pipeline's folder.

What it does, and prints (labels and verdicts only, never a count):

1. **The lock.** It checks the plan, the settings, the defaults and every input against the lock, and stops
   with `the locked plan doesn't match` if anything changed. It screens again, without outcomes, and stops
   unless the screen gives the locked plan's decisions.
2. **The checks that need outcomes**, each with its verdict (`pass`, `applied` or `warning`):
   - time at risk from the times themselves: a patient whose follow-up ends on or before entry is excluded,
     at the cohort rule holding the check `entry_before_follow_up_end`;
   - survival after the entry event: with too few deaths soon after it, entry moves later by the window
     (`applied`), or a warning only, as the settings say;
   - the horizon against follow-up: a warning if a group has too few patients at risk at the horizon;
   - the smallest-group rule, again: a comparison that falls short is described, as the settings say; one the
     settings mark `block` stops the run for a person to decide, logged as a deviation;
   - the adjustment variables: the missing-data rule, and the year bands after small ones are merged. When
     more patients lack an adjustment variable than the settings' share allows, every adjusted and weighted model
     uses multiple imputation, as `missing_data.imputation_model` says (`applied`); if that choice is left out,
     the run stops. The screen says so before the lock (step 5).
3. **R**, in five steps (`primary`, `adjusted`, `weighted`, `assumptions`, `secondary`). With 1,000 resamples
   it takes a few minutes. Each step's log, results and timing stay in `<work folder>\run_plan_v<n>\`, with the
   row-level files R reads. If a step fails, the run stops and names its log.

Add `--reuse-r` to keep R's results from an earlier run when everything R reads is byte for byte the same (to
rebuild the reports without waiting for R).

### Step 9. The safety gate (available)

It runs as part of step 8; it can't be skipped. Every table, the figure and the run report are built in
`<work folder>\run_plan_v<n>\staging\`, then read back. The run prints each check and its verdict:

- counts below the minimum are hidden, and so are estimates from fewer patients;
- counts that subtraction could recover, within one table or across tables, are hidden too (counts and
  descriptives are judged together);
- a group's curves are withheld when their steps would give away a number at risk the risk table hides
  (`withheld`, with the reason in the run report);
- survival estimates at fixed times are left out for a group whose size is hidden;
- balance before weighting is shown at the most decimals at which no small group-by-trait count can be
  rebuilt from published totals;
- the security scan: no ID from the data, no row-level column, no date but the run's and the plan's, no small
  count, and only tables, figures and reports.

Only when every check passes are the files copied to the output folder. If one fails, the run stops with
`not safe to publish`, names the problem, and the output folder is left untouched; the staged files stay in
the work folder for a person to look at.

### Step 10. The reports (available)

In the output folder:

| File | What it holds |
|---|---|
| `run_report.html`, `run_report.md` | what was done and why: the inputs and their checksums, the attribution the source's terms require (when the profile gives one), the capability profile, the screen's three lists from the locked plan, the plan's hashes, each check with its verdict, what each analysis did, notes (for example curves withheld) and the results in brief |
| `flow.csv` | the cohort's steps, merged so that no small count is shown or can be worked out |
| `counts.csv` | each analysis's group sizes, and the survival-after-entry count |
| `descriptives.csv` | the patients kept but counted, and the distributions described |
| `risk_table.csv` | patients at risk at each fixed time, by clock and group |
| `curves.csv`, `figure_km.svg` | the Kaplan-Meier steps that may be drawn, and the figure (text kept as text) |
| `estimates.csv` | every estimate: medians, RMST, survival at fixed times, hazard ratios, the assumption checks |
| `bias.csv` | every naive-minus-corrected contrast, with its interval, bar and label |
| `balance.csv` | pass 3b's covariate balance, before and after weighting |
| `s3_delays.csv` | the wait from diagnosis to the entry event (when the settings run it) |
| `outputs_manifest.json` | the SHA-256 of every input and every output |

Read the run report first: its checks say what the run changed (for example entry moved later) and why, and
its notes say what was left out to protect small counts. The labels follow the settings: "clearly material"
when the estimate and its whole interval pass the bar, "material" when the estimate does, "not material"
otherwise; a comparison the smallest-group rule doesn't fit is "not compared", and a group it describes is
"described, not labeled".

---

## Part 4. The acceptance checks (CRC BPC only)

**Phase 1.** On a machine holding the GENIE BPC CRC v2.0-public release, this screens it with the CRC settings and
compares the decisions and their reasons with the registered study's design gate. It prints only labels
and pass or fail, never a count:

```powershell
python scripts\acceptance_phase1.py <clinical folder> <genomic folder> [<work folder>]
```

The last line should read `RESULT: PASS`.

### The phase 2 acceptance check

After a run of the CRC plan (step 8), this compares the run's outputs with the CRC BPC study's published
`outputs/` and the numbers its README's Results section quotes. It prints only labels and match or differ:

```powershell
python scripts\acceptance_phase2.py <CRC BPC repository> <the run's output folder>
```

The last line should read `RESULT: PASS`. Lines saying "worded differently", "not produced" or "pipeline only"
name differences of wording or of what the general pipeline can describe; the script lists each.

---

## Part 5. When something goes wrong

| Message | What it means | What to do |
|---|---|---|
| `N problem(s) in the settings`, each with a line | a choice is missing, blank, misspelled or of the wrong kind | fix each line named, run `check-settings` again |
| `the key '...' appears twice` | a key is written twice in the settings | keep one |
| `put it in quotes` | YAML read a word such as Yes or No as true or false | write `'Yes'` |
| `... has not been confirmed` | the flat-file or OMOP mapping is still a draft | complete step 4 and confirm it |
| `values the labels don't list` (OMOP) | a field's `labels` miss a value the database holds | list it, as what it reads as or as `null` for blank |
| `the CDM has no ... table` / `has no column` (OMOP) | the mapping uses a table or column the database lacks | export it, or change the mapping |
| `the attribution the source's terms require: warning` (a run check) | a part of the profile's `attribution` is blank, for GENIE BPC the acknowledgment until it is pasted in | paste it from the source's terms into the profile (for GENIE BPC, `src/delayed_entry/adapters/genie_bpc/profile.yml`), then run again before publishing |
| `looks cloud-synced` | a data or work folder is in OneDrive, Dropbox or similar | move it to local disk |
| `reading it isn't written yet` (GENIE BPC) | the release has a panel coverage or size file | not supported yet: raise it |
| a `block` in the readiness report | that analysis can't run as planned | read its effect; change the settings or the data, or log a deviation |
| `the locked plan doesn't match` / `CHANGED:` | a plan file, the settings, the defaults or an input changed since the lock | restore it, or lock a new version with a reason |
| `the screen no longer gives the locked plan's decisions` | the run's screen differs from the locked one (other code or data) | find out why; a new plan version if it stands |
| `the smallest-group rule, applied again, blocks` | after the checks that need outcomes, a group the settings mark `block` is short | a person decides what is done instead, and locks a new version with the reason |
| `needs multiple imputation, and the settings don't say how` | more patients lack an adjustment variable than the settings' share allows, and `missing_data.imputation_model` is left out | add it (a new plan version if the plan is locked) |
| `R step ... failed`, and its log says `An Application Control policy has blocked this file` | Windows blocked an R package's DLL (for example `lme4.dll`, which the mice engine needs) | use `engine: chained_equations` in the imputation model, which needs only R, nnet and MASS |
| `falls back to mutation burden and can't count it` (a screen block) | the fallback needs `results.variant_class` and one coding size above zero for each panel | fix the tables, or change the analysis's fallback in a new plan version |
| `every year band is small` | the year bands can't be merged into two that are big enough | change the bands in a new plan version |
| `R step ... failed` | R stopped; its log is in `<work folder>\run_plan_v<n>\` | read the log (it stays local) |
| `not safe to publish` | the safety gate found a problem; nothing was copied | read the problem it names; merge or hide more in a new version |
| `must not be inside the work folder` / `holds other files` | the output folder is in the wrong place, or not empty | give an empty folder outside the work folder |

---

## Part 6. On this machine

Notes for one machine (its folders, its R library, where its data and its locked plans live) go in
`docs/local/`, which git leaves out, so they are never published. The machine this pipeline was built
on keeps its notes there.

---

## Changes to this guide

| Date | Change |
|---|---|
| 2026-11-13 | First version, at the end of phase 1: setup, tests, settings, settings lock, download check, flat-file mapping, screen, review, plan lock. |
| 2026-11-13 | Step 3: the GENIE BPC download script, with the CRC and NSCLC v2.0-public releases; patient-level parts only after the settings lock. |
| 2026-11-13 | Step 4: the release list for NSCLC v2.0-public, made with `scripts/synapse_release_list.py`. |
| 2026-11-16 | Part 6: the CRC plan's lock on this machine, and the two files it fingerprints. |
| 2026-11-16 | Phase 2: steps 8-10 (the run, the safety gate, the reports) with their commands; the R requirements; the R tests; the phase 2 acceptance check; new messages in Part 5; the run on this machine in Part 6. Settings: an adjustment variable's levels can be merged under new names (`merge`). |
| 2026-11-16 | Multiple imputation: the optional choice `missing_data.imputation_model` (engine, imputation model, methods, donors, iterations); the screen's block when it is needed and missing; step 8's check; two messages in Part 5. |
| 2026-11-16 | A biomarker's fallback by mutation burden: high against low on the entry test, decided by the smallest-group rule at the screen and the run; the gate hides a biomarker's largest group when its unknown and set-aside patients together are small; one message in Part 5. |
| 2026-11-16 | The lock: the settings file is hashed as it is, line endings included, so a file saved with Windows line endings now verifies (before, its lock always reported the settings file as changed). A lock of a file with Unix line endings, such as the CRC plan's, is unaffected. |
| 2026-11-17 | The settings check now refuses a biomarker analysis whose reference isn't one of its groups (before, that check never ran). |
| 2026-11-17 | The OMOP CDM adapter (step 4): a folder of table files or a DuckDB file, read through a confirmed mapping (`omop-template` writes a blank one); `source.adapter: omop`; `--mapping` at steps 5, 7 and 8; the made-up OMOP database in the tests (2.2); two messages in Part 5. |
| 2026-11-17 | The attribution a source's terms require: an optional capability-profile item (`attribution`), printed in the run report (step 10), with a run check that warns while a part is blank; GENIE BPC's profile carries the consortium citation and the release read; one message in Part 5. |
