# Design plan: a general pipeline for survival by genotype with delayed entry

Status: DRAFT for review, 2026-11-13. Nothing is built until this plan is approved. The five open questions were settled on 2026-11-13 (section 14).

## 0. What it is, and what it is not

**What it is.** A tool that takes survival data with genomic results, for any cancer and from any data holder, and:
1. reads it, whatever the file format, into one standard set of tables;
2. screens it, before any outcome is looked at, and reports what is present, what is missing, and which deviations from the standard analysis the data force;
3. locks a filled-in analysis plan;
4. runs the two-clock analysis (counted from diagnosis vs counted from first report, with delayed entry);
5. runs the privacy and security checks, and publishes nothing that fails them;
6. reports the results and a run report saying exactly what was done and why.

**What it is not.**
- Not a change to the CRC BPC project. That analysis is registered (osf.io/eyq6b) and finished; its repository stays as registered, with its deviations log. The new tool lives in its own folder and reuses the project's code by copying it into a new package with its tests, so the registered repository never changes.
- Not automatic clinical judgment. Which patients, which genotype groups, which time horizon and which thresholds count as material are set by a person who knows the cancer, and locked before the data are seen.
- Not a way to get data. Users bring data they are allowed to use, under their own data-use agreement.

## 1. Ground rules carried over from the CRC BPC project

- **Every choice is the user's.** No clinical or statistical choice is written into the code. Each lives in the settings file, with a documented default and a plain explanation, and is locked with the plan. The CRC settings are one ready-made example, not the tool's rules (section 4).
- **Plan before outcomes.** The screen uses only structure, completeness and counts, never survival. The plan is locked (with a hash and a timestamp) before any outcome is estimated, and the run refuses to start without a lock that matches the current settings and inputs.
- **Privacy by design**, with the five protections the CRC project needed, generalized and configurable (section 9).
- **Data stay local.** No network calls; inputs read-only; no data folder inside a cloud-synced folder (the project's `check_not_cloud`); outputs are aggregate only.
- **Tests on made-up data only**, never real counts. Every adapter, rule and check has tests.
- **Every published number traceable** to a file and a rule.

## 2. How it fits together

Six stages, and two layers of settings that make it general.

```
 data files ──► 1 READ ──► standard tables ──► 2 SCREEN ──► readiness report
   (any format)    ▲                               ▲              │
                   │                               │              ▼
           Layer B: data-holder           Layer A: cancer    3 LOCK ──► locked plan (hash)
           adapter + capability profile   settings                         │
                                                                           ▼
                        6 REPORT ◄── 5 SAFETY GATE ◄── 4 RUN ◄─────────────┘
                   results + run report   (privacy, security)   (both clocks, checks)
```

- **Layer A, cancer settings** (section 4): who is in the cohort, how genotype groups are defined, which adjustment variables, the time horizon and the materiality thresholds, which secondary and sensitivity analyses.
- **Layer B, data-holder adapters** (section 5): how to read a source into the standard tables, plus a profile of what the source can supply.

## 3. The standard tables

Every adapter produces these. All times are kept internally as **days from diagnosis**, whatever the source gives: real dates, shifted dates, or intervals (as GENIE BPC gives).

| Table | Fields (required in bold) | Used for |
|---|---|---|
| patients | **patient_id**, **vital status**, **days to death or last follow-up**, age at diagnosis, sex, stage at diagnosis, histology code (ICD-O-3), primary site code, year of diagnosis, center | cohort rules, outcomes, adjustment |
| tests | **patient_id**, **test_id**, **days to the entry event** and **which event it is** (report, sample collection, or order), panel_id | the corrected clock's entry time |
| results | **test_id**, gene, kind (mutation, fusion, copy number, protein test, MSI/MMR), variant or value | genotype groups |
| panels | panel_id, gene, regions covered (codons or exons) | confirming "wild type" (a gene the panel never tested isn't wild type) |
| treatments (optional) | patient_id, regimen, line, days to start | the treatment-clock analysis |

Each field is marked, per analysis component, as required, optional, or "needed for this one check only", so the screen can say exactly which analyses a missing field blocks.

## 4. Layer A: cancer settings

One settings file per study (YAML), validated against a schema, locked with the plan. What it holds:

- **Cohort:** stage at diagnosis, histologies to include or exclude (ICD-O-3), primary sites to include or exclude, minimum age, and the order the rules apply in (for the flow diagram).
- **Genotype groups:** rules that map results to groups, with precedence (for example "RAS variant and BRAF V600E together count as RAS-mutant"), the reference group, and the regions a panel must cover for a patient to count as wild type.
- **Adjustment variables** for the adjusted model, and the traits used for weighting.
- **Time horizon** for months alive (3 years for metastatic colorectal cancer; shorter for pancreatic, longer for breast or prostate) and the **materiality thresholds** (months for RMST and median, the ratio change for hazard ratios).
- **Secondary and sensitivity analyses**, each switched on or off with its parameters: the short-wait check (and its cutoff), the treatment clock, histology restrictions, a biomarker described only when too much is unknown, and so on.

Shipped with the tool:
- **CRC settings** that reproduce the registered CRC BPC analysis exactly (the acceptance test, section 12).
- **An NSCLC (lung) settings file** as the second example, with groups such as EGFR (exon 19 deletions, L858R), ALK fusions and KRAS G12C (since question 5, one KRAS group, with G12C as a sensitivity analysis). You draft each choice from published sources, each cited in the file, and decide it, as for the CRC study. There is no clinical review (question E, reversed 2026-11-20), stated as a limitation.

Every settings file records, for each choice, who made it and the source behind it.

### What users decide

Every item below is a setting. The CRC column gives the registered CRC BPC plan's value, which the CRC settings file uses. A new study starts from a blank template, where every choice must be filled in, or from a ready-made settings file such as the CRC one. A default is never applied silently: the locked plan lists every value used and where it came from.

| Choice | CRC value | What it controls |
|---|---|---|
| Cohort rules and their order | stage IV at diagnosis; no squamous, small cell or sarcoma histologies; no appendix; age 18+ | who is analyzed, and the flow diagram |
| Genotype groups, precedence, reference group | wild type, RAS-mutant, BRAF V600E; both variants count as RAS-mutant; wild type is the reference | the groups compared |
| Panel coverage needed for wild type | the listed KRAS, NRAS and BRAF codons | when "wild type" can be trusted |
| Entry event | first sequencing report | where the corrected clock starts |
| Time horizon for months alive | 3 years | the RMST window |
| Summary times | survival at 1 to 5 years; medians past 5 years counted as not reached | the fixed-time results |
| Materiality thresholds | 2 months (RMST), 3 months (median), a 10% change in a hazard ratio, 2 months (change in an RMST gap) | which results are material |
| Label rules | clearly material = estimate and whole interval past the threshold; material = estimate past it | how results are labeled |
| Smallest group compared | 30 patients | below it, a group is described, not compared |
| Biomarker unknown cutoff | over 30% unknown | the biomarker is described only |
| Survival-after-sequencing check | deaths within 730 days of the report; fewer than 11 moves entry | whether entry moves later |
| Stability rule | 20 patients in the corrected risk set | the uncertain-start restart (t*) |
| Uncertainty | 1,000 bootstrap resamples, fixed seed, 95% intervals | the intervals |
| Adjustment variables | age (with a spline), sex, center, tumor location, year of diagnosis | the adjusted model |
| Weighting | the same traits; weights trimmed at the 1st and 99th percentiles; balance target 0.1 | the weighted curves |
| Proportional-hazards check | bar p < 0.01; hazard ratios split at 12 months if it fails | when one hazard ratio isn't enough |
| Assumption check | the MB estimator (registered); the censoring-weighted IPW1 and IPW2 were added post hoc | the delayed-entry assumption |
| Sensitivity analyses | short waits (report within 6 months), treatment clock, other choices; each on or off | which checks run |
| Privacy minimum count | 11 | the smallest count shown |

### The guardrails users can't switch off

- **No outcome is read before the plan is locked.** The screen and the settings never see survival.
- **The safety gate always runs.** The minimum count can be raised but never set below the data-use agreement's.
- **A locked choice changes only by a new plan version.** That records who changed it, when and why, and adds a line to the study's deviations log.
- **Warnings, not overrides.** The screen flags a choice that looks implausible for the data, such as a horizon longer than follow-up or a smallest group under 11. The user decides, and the warning stays in the run report.

### How users set them

- **Phase 2:** a commented settings template (YAML), in two forms: a blank one where every choice must be filled in, and the CRC file as a worked example. A schema check rejects a malformed file and names the line.
- **Phase 4:** a guided setup for people who'd rather not edit a settings file. It asks each choice in plain words, shows the default and why, and writes the file. It comes as a short command-line questionnaire and as a spreadsheet form a clinician can fill in, which the tool reads back and checks.

## 5. Layer B: data-holder adapters and capability profiles

An adapter is a reader, a mapping (source columns to standard fields, with their transforms), and a **capability profile** saying what the source can supply:

| Profile item | Options | Why it matters |
|---|---|---|
| Entry event | report date / sample collection date / order date / none | defines the corrected clock; anything but the report date is a deviation; none means no correction is possible |
| Panel coverage | regions / genes only / none | whether "wild type" can be confirmed |
| Death capture | linked to a national death index / partial / unknown | incomplete death records bias survival itself |
| Time format | real dates / shifted dates / intervals | how times are converted to days |
| Privacy rules | the data-use agreement's minimum count (5, 11, 20...) and other limits | the safety gate's settings |
| Known quirks | e.g. withheld values, ages top-coded at 89, several reports per patient | handled in the adapter, reported in the screen |

Planned adapters, in order:
1. **Generic flat files:** CSV/TSV, Excel and Parquet first; SAS, Stata and SPSS later. A mapping file, which the tool drafts by matching column names and which a person confirms. Guessed mappings are never used unconfirmed.
2. **GENIE BPC, every cohort**, built from the CRC project's cohort code. Per the genieBPC package's release list (version 2.2.0), public releases exist for NSCLC (v2.0), CRC (v2.0), breast (v1.0), pancreatic (v1.0) and prostate (v1.0); bladder, kidney and newer releases are consortium-only.
3. **OMOP CDM** (the common data model many health systems use), built from the CRC project's OMOP module, with genomic results where the source's OMOP extension has them.

Later: genomic file formats (MAF, VCF), REDCap exports. Commercial holders' formats (for example Flatiron, Tempus, Foundation Medicine) only with a license; until then, tested on made-up data only.

## 6. Stage 2: the readiness screen

Runs without reading any outcome. A test enforces this: scrambling the outcome columns must leave the readiness report unchanged.

What it checks:
- required fields present, of the right type and units;
- dates or times in order; duplicates; impossible values (negative times, ages over the top-code);
- a first report that comes back after follow-up ends (no time at risk on the corrected clock);
- missingness of every field, and group sizes (counts only);
- panel coverage for each patient's genotype assignment;
- the capability profile against what each analysis needs.

What it reports, in three lists:
- **Ready:** analyses that can run as standard.
- **Missing:** what is absent, and which analyses that blocks.
- **Deviations required:** each with the rule that triggered it, why, and its effect. Severity: *block* (that analysis can't run), *deviation* (it runs differently, stated in the plan), or *note*.

Examples, generalized from the CRC plan's decision rules. The numbers shown are the CRC settings' values; the screen reads each one from the study's settings file, never from code.

| Condition found | Deviation |
|---|---|
| No report dates, but sample collection dates | Entry at sample collection, stated as a deviation |
| No entry event at all | Naive clock only; no bias estimate is possible; loud warning |
| A genotype group under 30 patients | Described, not compared (as design check D1) |
| More than 30% unknown for a biomarker | That biomarker described only |
| No panel coverage | Wild type unconfirmed; flagged; sensitivity analysis added |
| Death capture partial or unknown | Warning in every result; stated as a limit |
| Few deaths soon after the first report | The registry may have required survival after sequencing; entry moved (as design check D2) |
| A report dated before diagnosis | Entry at day 0 |
| Any count of 1 to the minimum | Flagged now for the safety gate |

## 7. Stage 3: locking the plan

The tool writes the filled-in plan in Markdown, for people, and YAML, for the run. It combines the cancer settings, the adapter's profile and the screen's deviations. A person reviews it. Locking records a SHA-256 hash of the plan, the settings and every input file, with a timestamp, in a manifest. The run checks the lock first and stops if anything changed. Registering the plan externally (for example on OSF) is a manual step the tool prepares for but doesn't do.

## 8. Stage 4: the run

Reused from the CRC project, with every number checked against it (section 12):
- **Passes 1 and 2:** Kaplan–Meier on both clocks, with the plan's median rule, months alive to the horizon (RMST) and Cox models.
- **Adjusted and weighted models**, with a balance report.
- **Paired bootstrap:** 1,000 resamples by default, seeded, with no interval when any resample has no value.
- **The stability rule (t\*)** for the uncertain start of the corrected curve.
- **The delayed-entry assumption check**, with more than one estimator, reporting when they disagree (the CRC project's SA6 lesson).
- **Labels** by the locked thresholds, and the sensitivity analyses switched on in the settings.

**One language limit shapes open question B.** In the CRC project, two things came from R alone:
- **Robust standard errors for the Cox models.** The Python survival library (lifelines) can't compute them with delayed entry.
- **The bootstrap resamples**, which were drawn in R.

So a Python-only tool can reproduce the point estimates exactly, but its intervals would be drawn fresh, and its Cox intervals would need another method. Keeping R for those parts reproduces everything exactly.

**Settled (question B): Python + R.**
- Python reads, screens, locks, gates and reports.
- R's survival package fits the Cox models with robust standard errors and draws the bootstrap resamples, as in the CRC project.
- Users need R installed. On this machine it runs through `Rterm.exe`, because `Rscript.exe` is blocked by an application-control policy.

## 9. Stage 5: the safety gate

Nothing reaches the output folder unless it passes. Generalized from the CRC project's `outputs.py`, with the minimum count from the data-use agreement:
1. counts below the minimum hidden;
2. counts that subtraction could recover, within one table and across tables;
3. curve steps that would reveal a hidden number at risk (curves withheld or cut);
4. survival estimates that would pin a hidden group size;
5. precision limits on balance statistics that, with published totals, would rebuild small tables.

Security checks: inputs opened read-only; outputs scanned for row-level patterns (IDs, dates, single-patient rows); data folders not in cloud-synced locations; no network access; checksums of inputs and outputs in the manifest; every check and its verdict in the run log.

## 10. Stage 6: the reports

- **Results:** tables (CSV) and figures (SVG, text kept as text), aggregate only.
- **Run report** (HTML and Markdown):
  - the inputs and their checksums;
  - the capability profile;
  - the screen's three lists;
  - the locked plan's hash;
  - each check with its verdict;
  - a plain summary of the results with their labels.

## 11. Where it lives

A new package folder, `Portfolio_Projects/delayed-entry-pipeline/`, built locally. It is published as its own public repository under the UhlConsultingLLC GitHub account once phase 2's acceptance test passes (question A). Only code and made-up test data are ever published, never patient data. The CRC BPC repository isn't restructured. The new package copies what it reuses (the passes, the privacy functions, the cohort and OMOP code, `reanalyze.py`) with their tests. The CRC project's own code stays as registered.

## 12. Testing

- Made-up data for every adapter, rule, deviation and check, built with the CRC project's synthetic-data approach.
- **The screen never reads outcomes:** scrambled outcomes give the same readiness report.
- **The safety gate:** a test for each of the five protections, built from the cases the CRC project found.
- **Acceptance test:** with the CRC settings and the GENIE BPC adapter, the pipeline reproduces the CRC project's published `outputs/`. It runs on your machine only, against the data already there, and prints only pass or fail.
  - Point estimates, counts, labels and the privacy decisions must match exactly.
  - Bootstrap intervals and robust Cox intervals match exactly only if the R parts are kept (section 8, question B). Otherwise they are checked for closeness, within what resampling noise allows.

## 13. Phases, each with what "done" means

1. **Standard tables, reading and screening.** The generic flat-file and GENIE BPC adapters, the readiness screen and the plan lock. The screen's thresholds come from a settings file (the CRC one), never from code, so phase 2's settings layer builds on it.
   *Done when:* the CRC BPC data, screened, give the deviations the registered project found (S1 descriptive only; BRAF V600E described, not compared, where under 30), with no outcome read.
2. **Settings, run, gate and report.** The cancer settings layer with the CRC file, the run (reusing the passes), the safety gate and the run report.
   *Done when:* the acceptance test passes, as section 12 defines it for the language chosen.
3. **A second cancer and a second holder.** NSCLC settings on GENIE BPC NSCLC v2.0-public, and the OMOP adapter on made-up OMOP data.
   - Access: you request NSCLC access on Synapse while phases 1 and 2 are built (question D).
   - The NSCLC settings are drafted, decided by the study author and locked before any NSCLC file is downloaded, the same order the CRC project used.

   *Done when:* both run end to end and their screens report correctly.
4. **More formats, a command line and the guided setup.** SAS, Stata and SPSS; MAF and VCF; a command-line interface; the guided setup (questionnaire and spreadsheet form, section 4); documentation.

## 14. Settled questions (2026-11-13)

- **A. Name and home:** "delayed-entry-pipeline". Built locally in `Portfolio_Projects/`; published as a public repository under the UhlConsultingLLC GitHub account once phase 2's acceptance test passes.
- **B. Language:** Python + R. R fits the Cox models with robust standard errors and draws the bootstrap resamples, so the CRC numbers reproduce exactly (section 8).
- **C. Second cancer:** lung, GENIE BPC NSCLC v2.0-public. It has well-defined driver groups with matching drugs, and its ALK fusions test the non-mutation result types.
- **D. Data access:** you request NSCLC access on Synapse now, in parallel with phases 1 and 2. Downloading waits until the NSCLC settings are locked. Requesting access is yours to do; it can't be done for you.
- **E. Clinical review:** you draft the lung settings from published sources, each choice cited. A lung oncologist reviews them before locking, and the settings file records who made each choice and its source.
  **Reversed 2026-11-20:** there is no clinical review. The study author decides each choice from the cited sources, as the CRC study's author did; the settings file records the study author as having made each choice, with its source; and the lack of clinical review is stated as a limitation. An ALK sensitivity analysis (SA12) is added to the NSCLC plan before any download: the reference group kept only where the entry test's panel can detect ALK fusions, since in 2014-2017 ALK fusions were often found by FISH or IHC rather than DNA panels, and a DNA panel detects them only if it tiles ALK's introns. (2026-11-23: question 6 was decided as option A, so only patients whose first panel detects ALK fusions are assignable; SA12 would then repeat the primary analysis and is switched off for NSCLC, while the rule stays for other studies.)
