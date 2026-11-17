# delayed-entry-pipeline

Survival by genotype with delayed entry, for any cancer and any data holder. The design is in
[DESIGN_PLAN.md](DESIGN_PLAN.md). This folder holds **phase 1**: reading a source into standard tables,
screening it against a study's settings without reading any outcome, and locking the filled-in plan; and
**phase 2**: the run, the safety gate and the reports.

It generalizes the CRC BPC delayed-entry study (registered at osf.io/eyq6b). That project is finished and
unchanged: the code reused here was copied from it, never imported.

## What phase 1 does

1. **Read** a source into the standard tables (`delayed_entry.tables`): patients, tests, results, panels,
   and optionally treatments and panel sizes. All times are days from diagnosis. Each field is marked,
   for each analysis component, as required, optional or needed for one check only.
   - `adapters.genie_bpc`: GENIE BPC releases, built from the CRC BPC project's cohort code. The cohort's
     own files (raw file list, release list, extra columns) live in `adapters/genie_bpc/cohorts/<cohort>/`.
     Mutations are always read; fusions and regimen drug names are read where the cohort's raw file list
     names them (NSCLC does, CRC doesn't). `adapters.genie_bpc.check_download` is the check to run right
     after downloading a release.
   - `adapters.flatfile`: CSV/TSV, Excel (.xlsx), Parquet, SAS (.sas7bdat and .xpt transport files), Stata
     (.dta), MAF and annotated VCF (.vcf, .vcf.gz), through a mapping file. The tool drafts the mapping from
     column names; a person must confirm it, and a draft is never used.
     - **MAF:** read as a tab-separated table, after its `#version` comment lines. The drafted mapping matches
       its standard columns (Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short, Variant_Classification).
     - **VCF:** a VCF holds protein changes only once VEP (its CSQ field) or SnpEff (ANN) has annotated it, so
       the reader (`adapters.vcf`) reads that annotation. It makes one row per variant record, per sample that
       carries it (by its genotype), per annotation kept, with the sample's name and each annotation subfield
       as columns. Two reading rules come from the mapping, with no default:
       - which annotations count: `canonical` (VEP's CANONICAL=YES), `first` (the allele's first listed) or
         `all`;
       - which records: `pass` (FILTER PASS or .) or `all`.

       An annotation belongs to the allele it names. VEP's form of an allele, with a shared first base
       removed and "-" for a deletion, counts too.
     - **Long protein changes:** the `hgvsp_short` transform puts any protein change in short form, for
       example VEP's `ENSP...:p.Gly12Asp` as `p.G12D`.
     - **Several files:** a path may be a pattern such as `vcf/*.vcf`, to read one file per sample into one
       table; its `source_file` column names each row's file.
     - **Dates:** a date stored as a date is read as it is (Parquet, Stata, and .sas7bdat columns with a SAS
       date format). Text dates take a `format`. A date held as a number, as a .xpt file always holds it,
       takes `days_from: "1960-01-01"` (or `seconds_from` for a date and time), the origin SAS and Stata
       count from.
     - **Stata value labels:** a coded value is read as its label, the text Stata shows, so a mapping's `map`
       names readable values. `value_labels: false` on the file reads the codes instead, which is also the way
       out when pandas refuses labels that aren't unique.
     - **Text encoding:** a .sas7bdat file's text is decoded with the encoding it records. A .xpt file records
       none and is read as ISO-8859-1, unless the file's `encoding` says otherwise.
     - **Whole numbers:** SAS stores every number as a decimal, so a column whose numbers are all whole is
       read as whole numbers. A numeric ID is then "1001", and still joins across files.
   - `adapters.omop`: an OMOP Common Data Model (version 5.4) database, either a folder with one file per
     CDM table (CSV, TSV or Parquet) or a DuckDB file, through a mapping a person confirms
     (`python -m delayed_entry omop-template` writes a blank one). OMOP codes everything with concepts, so the
     mapping names the concept sets: the study's cancer, the sequencing tests, readings such as MSI, the
     cancer drugs. It also says where each patient field comes from: a person, care site or index-condition
     column, or the measurement or observation row dated closest to diagnosis, with a regular expression to
     pick part of a value and labels to rename values (labels must list every value found). Variants come
     from the genomic extension table the mapping names, with one row per variant linked to the sequencing
     test. Fusions (two genes) and copy numbers can come from it too, and protein changes are put in HGVSp
     short form. The adapter works out the rest from the dates:
     - the diagnosis is the earliest condition in its concept set;
     - every time is in days from it;
     - the outcome is death, or alive at the end of the last observation period;
     - a test's `after_follow_up_end` comes from the same dates;
     - same-day drug exposures form one regimen, numbered by date, as the CRC plan's section 4d did.
     Panels and panel sizes, which OMOP doesn't hold, can come from files.
   - Each adapter declares a **capability profile** (`delayed_entry.profile`): entry event, panel
     coverage, panel sizes, death capture, time format, the data-use agreement's minimum count, and known
     quirks. It may also give the **attribution** the source's terms require: `{citation, version,
     acknowledgment}`. The run report prints it, and each run warns while a part of it is blank. GENIE BPC's
     profile gives the consortium citation, the adapter fills in the release the settings name with its
     Synapse ID, and the acknowledgment is pasted in from the GENIE terms of use
     (`adapters/genie_bpc/profile.yml`).
2. **Screen** (`delayed_entry.screen`): three lists (ready, missing, deviations required), each item with
   a severity (block, deviation or note), plus warnings on implausible choices (never overrides), the
   checks that need outcomes, and every setting read with where it came from.
3. **Lock** (`delayed_entry.lock`): the plan in Markdown and YAML, with a SHA-256 manifest of the plan,
   the settings and every input file, and a timestamp. A locked plan changes only by a new version with a
   reason, which adds a line to the deviations log. The settings can also be locked alone, before any
   download.

## What phase 2 does

`python -m delayed_entry run` (`delayed_entry.run`) carries out the locked plan, in this order:

1. **The lock first.** It checks the plan files, the settings, the defaults file and every input against the
   lock's SHA-256 hashes, and stops if anything changed. Then it screens again, without outcomes, and stops
   unless the screen still gives the locked plan's decisions and group sizes, or if the plan has a block on
   the primary analysis.
2. **The checks that need outcomes** (`delayed_entry.analysis`), listed in the readiness report's last section:
   - time at risk from the times themselves: a patient whose follow-up ends on or before entry is excluded
     at the cohort rule that holds the structural check `entry_before_follow_up_end` (the same flow step);
   - survival after the entry event (the CRC study's D2), with the settings' `survival_after_entry` window,
     minimum and action: entry moves later by the window, or a warning only;
   - the horizon against follow-up: a warning if a group's corrected risk set at the horizon is below the
     stability rule's size;
   - the smallest-group rule, applied again to the patients each analysis keeps. A comparison the plan
     compares and the run can't is described (or, if the settings say `block`, the run stops for a person to
     decide). A comparison the plan describes stays described.
3. **The adjustment variables** (`delayed_entry.covariates`), coded as the settings say: a spline's knots,
   a categorical variable's kept levels and reference, **levels merged under new names** (`merge`, below),
   and year bands, merged when small by the CRC plan's rule made general. The missing-data rule chooses
   complete cases, or **multiple imputation** (below).
4. **The passes, in R** (`src/delayed_entry/R/`, ported from the CRC study's R code): passes 1 and 2 on both
   clocks (naive from diagnosis, corrected from the entry event), with KM survival at the fixed times,
   medians, RMST to the horizon and Cox models with robust standard errors; the paired bootstrap; the
   stability rule (t*) and its conditional contrasts; pass 3a (adjusted Cox models); pass 3b (weighted
   curves, with the balance table); the proportional-hazards and quasi-independence checks; and every
   secondary and sensitivity analysis the plan runs (subsets, the treatment clock, the biomarker and its
   variants, the stratified model). Every value R uses comes from `run_spec.json`, which the run writes from
   the locked settings.
5. **The safety gate** (`delayed_entry.gate`), generalized from the CRC study's `outputs.py`, with the larger
   of the settings' and the data-use agreement's minimum count:
   1. counts below the minimum hidden, and estimates from fewer patients;
   2. counts that subtraction could recover, within one table and across tables (each table's relations,
      then counts and descriptives together). A biomarker's unknown and set-aside patients are merged if
      either is small; if the two together are small too, its largest group's size is not shown either, since
      the groups would leave them over from the analysis's total;
   3. curves withheld when their steps would give away a number at risk the risk table hides;
   4. survival estimates at fixed times left out for a group whose size is hidden;
   5. balance before weighting published at the most decimals at which no small group-by-trait cell can be
      rebuilt from published totals (the gate checks every table that fits).

   Then the security checks: every file to be published is read back and scanned for the data's own IDs,
   row-level columns, dates that aren't the run's or the plan's, small counts, and files that aren't tables,
   figures or reports. The tables are built in the run's own folder (inside the work folder) and copied to
   the output folder only after every check has passed.
6. **The reports** (`delayed_entry.reports`, `delayed_entry.run_report`): `flow.csv`, `counts.csv`,
   `descriptives.csv`, `risk_table.csv`, `curves.csv`, `estimates.csv`, `bias.csv`, `balance.csv`,
   `s3_delays.csv` and `figure_km.svg` (text kept as text), the run report in Markdown and HTML (the inputs
   and their checksums, the capability profile, the screen's three lists, the plan's hashes, each check with
   its verdict, and a plain summary of the results with their labels), and `outputs_manifest.json` with the
   SHA-256 of every input and output. No local path and no row-level value is written to the output folder.

The run's row-level files (`analysis_cohort.parquet` and the others R reads, and R's results) stay in
`<work folder>/run_plan_v<n>/`.

## Settings: every choice is the user's

No clinical or statistical value is written into the code (a test checks this). Every choice lives in a
settings file, with who made it and the source behind it:

- `settings/template_blank.yaml`: every choice blank. Loading it lists every choice still to fill in,
  each with its line. A missing choice, a blank value or an unknown key is an error naming the line.
- `settings/crc_bpc.yaml`: the registered CRC BPC values, as a worked example.
- `settings/nsclc_bpc.yaml`: a DRAFT for NSCLC (GENIE BPC NSCLC v2.0-public), awaiting a lung
  oncologist's review. The review's questions and answers will be published once the settings are locked.
- A value may be the word `default`, which takes the documented default in
  `src/delayed_entry/defaults.yml`. The screen and the locked plan list every default used, with its
  source. Choices specific to a study (cohort rules, groups, coverage, horizon, analyses, privacy minimum)
  have no default.
- **Multiple imputation** follows the CRC plan's section 9, made general, with every choice in the settings
  (`missing_data.imputations` and `missing_data.imputation_model`). Only adjustment variables are imputed, by
  chained equations; the imputation model holds the settings' outcome terms (for example the event indicator,
  the Nelson-Aalen cumulative hazard with delayed entry, and the entry time) and the group. Hazard ratios are
  pooled by Rubin's rules; each bootstrap resample is imputed once; other point estimates, and pass 3b's balance,
  are averages across the imputations. The engine is the pipeline's own (`chained_equations`, in R with nnet and
  MASS: predictive mean matching, Bayesian linear, logistic, multinomial or ordered draws, as mice does them) or
  the mice package (`mice`). `missing_data.imputation_model` is the one optional choice: a study that never needs
  imputation may leave it out. If the complete-case rule calls for imputation and it is missing, the screen blocks
  the plan, before any outcome is read, and the run stops.
- **A biomarker's fallback by mutation burden** follows the CRC plan's section 4b, made general. When more
  of a biomarker analysis's patients are Unknown than `biomarker.unknown_max_share` allows, and its `fallback`
  says `mutation_burden` when panel sizes are given, the analysis compares high with low mutation burden instead,
  labeled a hypermutation proxy, and the variants it names in `not_run_under_fallback` are not run. The burden
  is the number of the entry test's mutations whose variant class `mutation_burden.variant_classes` lists (each
  result row once), over the coding size of the entry test's panel (`panel_sizes.coding_megabases`). High is at
  or above `threshold_per_megabase`, judged after rounding to 9 decimals, so that a burden exactly at the
  threshold is high despite floating point; low, below it, is the reference. A patient whose entry test's panel
  has no size is Unknown: described, not compared. The smallest-group rule decides the comparison at the screen
  and again at the run, as for any comparison. Panel sizes count only when the capability profile says the
  source gives them (`panel_sizes: given`); otherwise the fallback is descriptive only. If the burden can't be
  counted (no `results.variant_class`, a panel given two sizes, or a size of zero or less), the screen blocks the
  analysis and it is not run.
- An adjustment variable of form `categorical` may merge the source's values into new levels: `merge` maps
  each new level to the values it takes in (a value in no group is missing), `keep` then names new levels,
  and `reference` names the reference level. For example, six smoking categories into never, former and
  current. The settings check refuses a value merged twice, or a kept or reference level that isn't one.

## The screen never reads an outcome

The screen is given `outcome_blind(tables)`: vital status and the time to death or last follow-up are
removed, and only whether each patient's outcome is recorded is kept. A test scrambles every outcome value
and checks the report is unchanged, on the made-up data and, locally, on the real CRC data.

Two checks the design plan lists for the screen need outcome values, so they run after the lock instead,
as the plan's locked rules, and the report lists them under "checks that need outcomes":

- deaths soon after the entry event (the CRC study's design check D2), and
- a test that comes back after follow-up ended, judged from the times themselves.

For the second, the screen uses what the source itself records: `tests.after_follow_up_end`. GENIE BPC
gives it (its report-after-death and report-after-last-contact flags, and the age-89 rule, which needs
only the age). With it, the screen reproduces the registered cohort patient for patient. Without it, the
screen says its group sizes are upper bounds and the smallest-group rule is applied again at the run.

## What the readiness report shows

Aggregate counts only, masked with the larger of the settings' minimum and the data-use agreement's (a
settings minimum below the agreement's is a block, and such a plan can't be locked). So that no hidden
count can be worked out by subtraction:

- if any group in an analysis is masked, or its groups fall short of a total shown elsewhere by a small
  number, that analysis's sizes are shown as ranges ("at least 30", "11 to 29", "<11"); once the primary
  sizes are ranges, every analysis's are;
- a sensitivity analysis's group, or a biomarker variant's, whose count differs from its reference by a
  small number is shown as a range;
- missing values are shown as a masked count, never with the table's total or a share.

Results of a biomarker test count "by the entry test" on the source's own clock: days from birth when the
source dates events that way (GENIE BPC does, and the CRC study compared on it), otherwise days from
diagnosis. The report also records each input file's SHA-256, and the lock accepts only those inputs.

## Running it

The step-by-step guide to setting up, testing and running the pipeline is
[docs/USER_GUIDE.md](docs/USER_GUIDE.md). In short:

The tests use made-up data only. They need Python 3.10+ with pandas, numpy, duckdb, PyYAML and pytest
(`pip install -e ".[test]"`). From this folder:

```bash
python -m pytest -q -p no:cacheprovider
```

The tests that run R (the whole run on the made-up fixture, and the R unit checks) are skipped unless the
environment variable `DELAYED_ENTRY_R` names Rterm.exe or Rscript (and `DELAYED_ENTRY_R_LIBRARY`, if the R
packages live in a folder of their own). The run needs R with survival, WeightIt, cobalt, tranSurv, duckdb
and jsonlite.

The command line (run from this folder, with `src` on the path or the package installed):

```bash
python -m delayed_entry template settings/my_study.yaml
python -m delayed_entry check-settings settings/my_study.yaml
python -m delayed_entry lock-settings settings/my_study.yaml <work folder>/settings_lock --locked-by "<reviewer>"
python -m delayed_entry draft-mapping mapping.yaml patients=patients.csv tests=tests.xlsx#Tests results=results.parquet panels=panels.csv
python -m delayed_entry omop-template omop_mapping.yaml
python -m delayed_entry screen settings/my_study.yaml <work folder> --mapping omop_mapping.yaml
python -m delayed_entry screen settings/crc_bpc.yaml <work folder> --clinical <clinical folder> --genomic <genomic folder>
python -m delayed_entry lock settings/crc_bpc.yaml <work folder> --clinical <clinical folder> --genomic <genomic folder> --locked-by "<reviewer>"
python -m delayed_entry verify-lock <work folder>/plan_v1/manifest.json
python -m delayed_entry run settings/crc_bpc.yaml <work folder> <output folder> --clinical <clinical folder> --genomic <genomic folder> --r <Rterm.exe or Rscript> [--r-library <R packages folder>]
```

The work folder must be on local disk, outside any cloud-synced folder and outside this package. The
output folder must be outside the work folder and outside this package. Reports hold aggregate counts
only, with any count below the privacy minimum shown as "<minimum".

**Phase 1 acceptance** (local only, prints labels and pass or fail, never a count):

```bash
python scripts/acceptance_phase1.py <clinical folder> <genomic folder> [<work folder>]
```

It screens the real GENIE BPC CRC v2.0-public data with the CRC settings and compares the decisions and
their reasons with the registered study's design gate (`scripts/crc_bpc_phase1_expected.yml`).

**Phase 2 acceptance** (local only, prints labels and match or differ, never a count or an estimate):

```bash
python scripts/acceptance_phase2.py <CRC BPC repository> <the run's output folder>
```

It compares every table of the CRC study's published `outputs/` with the run's, row by row, and checks each
number its README's Results section quotes (read from the README when it runs, so none is stored here).
Passed on 2026-11-16, on the CRC plan locked that day: every count, estimate, interval, label, privacy
decision, curve and the figure's curves and numbers at risk match. What differs is wording only (the general
pipeline's names for the same rows and labels; the script lists each), two descriptive tables the standard
tables can't make (how the year of diagnosis was worked out, and the specimen type), two flags only the
general cohort code counts, and one note that says why a blank pass 3b interval is blank.

## Layout

```
settings/                 the CRC settings and the blank template
scripts/                  the phase 1 and phase 2 acceptance checks (local data only), the Synapse download
                          and release list, and the Markdown-to-Word converter for documents sent outside
src/delayed_entry/
  tables.py               the standard tables and the outcome-blind view
  settings.py             the settings schema and loader; defaults.yml
  profile.py              capability profiles
  cohort.py               cohort rules, the first test, genotype groups, panel coverage
  screen.py               the readiness screen
  lock.py                 the plan lock
  privacy.py              small counts, the public flow diagram, local folders
  analysis.py             the run's analysis cohort: the checks that need outcomes, who each analysis keeps
  covariates.py           the adjustment variables: merged and kept levels, year bands, the missing-data rule
  run.py                  the run: the lock, the outcome checks, R, the gate, publishing
  R/                      the passes, the bootstrap and the models in R (from the CRC study's R code)
  gate.py                 the safety gate: the five privacy protections and the security scan
  reports.py              the results tables and the figure
  run_report.py           the run report
  adapters/               genie_bpc (with check_download), flatfile, xlsx, vcf, hgvs, omop
tests/                    made-up data only; fixtures/ holds the CRC project's made-up BPC files; R/ the R checks;
                          omop_fixture.py writes those files' tables as a made-up OMOP database
```

## Known limits

- Panel coverage files (codon regions) and panel sizes are not read from GENIE BPC yet: the adapter stops
  if it finds one, as the CRC project did. The flat-file adapter reads codon regions and panel sizes.
- No panel file shows whether a panel detects a gene's fusions, so that comes from the settings
  (`genotype.coverage.fusion_panels`), with the panels' documentation as its source.
- Excel files are read with the standard library (cached cell values, first row as header).
- A VCF's structural variants and copy-number records are read like any other record, but only an annotated
  protein change places a patient in a mutation group. A combined consequence such as
  `missense_variant&splice_region_variant` must be listed whole in a `map`. Fusions and copy numbers come from a
  table the mapping names, not from a VCF.
- SPSS (.sav) files aren't read: that needs the pyreadstat package, which isn't installed. The SAS .sas7bdat
  reading is pandas' own, and its test files can't be made without such a package, so only .xpt and Stata
  files are tested here, on made-up files a small writer in the tests makes (`tests/sas_xport.py`).
- GENIE BPC gives no panel sizes (above), so with its releases a biomarker's fallback is descriptive only.
  Mutation burden needs a source that gives them, such as a flat file, or an OMOP database with a panel sizes
  file.
- OMOP has no standard home for genomic variants, so the OMOP adapter reads them from the extension table the
  mapping names, with its columns named there too. It has been tested on made-up data only: the made-up GENIE
  BPC fixture written as an OMOP database, which it reads back field for field and screens to the same
  decisions as the GENIE BPC adapter.
- OMOP keeps no test order, so tests on the same day are told apart by their ID. It keeps no regimen number or
  line either, so `line` is blank. A condition, test or drug can't be tied to one of a person's cancers rather
  than another, so `multiple_cohort_cancers` and `other_primary_same_cancer` come from the mapping or are
  absent.
- OMOP's person_id is a plain number, and the safety gate's scan can't tell a plain number from a count, so it
  leaves such IDs out and says so. Where person_source_value isn't a plain number, using it as the patient ID
  lets the scan look for it.
- On Windows with Smart App Control on, the mice engine may not run: Smart App Control can block `lme4.dll`,
  which mice loads through mitml and jomo. The pipeline's own engine (`chained_equations`) needs only R, nnet
  and MASS. Its imputations follow mice's algorithms but are not mice's number for number.
- The descriptive tables are the flags kept but counted, each field a subset analysis keeps patients by, and
  each extra field the adapter carries that the model doesn't use. A source-specific description (for
  example how the year of diagnosis was worked out) needs that field in the standard tables first.
- R's point estimates are not re-derived in Python (the CRC study did that for passes 1 and 2); the phase 2
  acceptance check against the CRC study's published results stands in for it.

## Data and attribution

This repository holds no GENIE data. The tests run on made-up files shaped like a GENIE BPC release. To run
the pipeline on GENIE BPC, request access on Synapse and download the data yourself under the AACR Project
GENIE terms of use; the data must not be redistributed. The files here that describe the releases (Synapse
IDs, file names, sizes, checksums and column names) are descriptions, not data.

Work that uses or references the GENIE data must cite the consortium and the release, and include the
acknowledgment the terms of use give:
- **Citation:** AACR Project GENIE Consortium. AACR Project GENIE: Powering Precision Medicine Through an
  International Consortium. Cancer Discov. 2017;7(8):818-831. Releases referenced here: GENIE BPC CRC
  v2.0-public (the acceptance checks) and GENIE BPC NSCLC v2.0-public (the draft NSCLC settings).
- **Acknowledgment:** The authors would like to acknowledge the American Association
  for Cancer Research and its financial and material support in the development of the
  AACR Project GENIE registry, as well as members of the consortium for their commitment
  to data sharing. Interpretations are the responsibility of study authors.

A run on GENIE BPC prints this attribution in its run report, with the release it read, from the GENIE BPC
capability profile (`src/delayed_entry/adapters/genie_bpc/profile.yml`).

## Security

How to report a security problem, and what never to post in an issue or pull request (patient data, small
counts, credentials), is in [SECURITY.md](SECURITY.md). The tests run on every push, on Windows
(`.github/workflows/tests.yml`): with Python 3.10 and 3.14, and once more with R and the packages in
`.github/r-packages`, so the whole runs, the R unit checks and the mice engine run there too. A guard among them
fails if a committed file looks like data, holds a local path or looks like a credential.

## License

MIT: see [LICENSE](LICENSE).
