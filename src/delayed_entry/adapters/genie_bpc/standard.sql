-- From the raw GENIE BPC tables (raw.load) to the standard tables. Run by
-- delayed_entry.adapters.genie_bpc.read(), one step at a time, with these parameters:
--   $censor          BOOLEAN    censor a withheld survival time at the age limit (settings: withheld_survival)
--   $age_limit       DOUBLE     the source's age limit (the capability profile's age_top_code)
--   $days_per_year   DOUBLE     settings: time.days_per_year
--   $labels          VARCHAR[]  this cancer's labels in the source (settings: cancer_type_labels)
-- and three pieces filled in before it runs: {first_order}, the settings' order for the first test;
-- {extras}, the cohort's own patient columns (cohorts/<cohort>/cohort.yml); and {regimen}, the regimen's drug
-- names (r.regimen_drugs when the cohort's raw_files.yml reads that column, otherwise blank).
--
-- Adapted from the CRC BPC project's sql/02_cohort_rules.sql, 03_msi_mmr.sql, 04_covariates.sql and
-- 06_secondary.sql. No cohort rule is applied here: who is in the cohort is the settings' business.
-- Values are kept as released (no trimming); the settings' rules say how to compare them.

-- step: index_cancer
-- The index cancer: the patient's first BPC Project cancer (lowest ca_seq).
CREATE OR REPLACE TABLE index_cancer AS
SELECT *
FROM raw_cancer
WHERE trim(redcap_ca_index) = 'Yes'
QUALIFY row_number() OVER (PARTITION BY record_id ORDER BY ca_seq) = 1;

-- step: panel_size
-- How many genes each panel reads.
CREATE OR REPLACE TABLE panel_size AS
SELECT panel, count(DISTINCT gene) AS n_genes
FROM raw_panel_gene
GROUP BY panel;

-- step: follow_up
-- The source's age rule. BPC withholds a redacted patient's survival time when it would reveal an age
-- over its limit. If the age at diagnosis is known, the patient was alive at the limit: with $censor,
-- they are censored there. Only the age and whether the time is withheld are used.
CREATE OR REPLACE TABLE follow_up AS
SELECT
    c.record_id,
    coalesce($censor AND trim(p.redacted) = 'Yes' AND c.tt_os_dx_days IS NULL AND c.age_dx IS NOT NULL, FALSE)
        AS censored_at_age_limit,
    ($age_limit - c.age_dx) * $days_per_year AS age_limit_day
FROM index_cancer AS c
JOIN raw_patient AS p ON p.record_id = c.record_id;

-- step: std_tests_ngs
-- One row per sequencing report for the index cancer. The join needs record_id AND ca_seq: a report for
-- the patient's other cancer doesn't count, and a test linked to two cancers is listed once, here.
-- after_follow_up_end: BPC's own flags (report on or after death, or after the last known alive date),
-- and for a patient censored at the age limit, the report against that limit (or a limit reached by
-- diagnosis). Blank when the source doesn't say.
CREATE OR REPLACE TABLE std_tests_ngs AS
SELECT
    r.record_id AS patient_id,
    r.cpt_genie_sample_id AS test_id,
    'ngs' AS modality,
    'report' AS entry_event,
    r.dx_cpt_rep_days AS days_to_entry_event,
    r.dob_cpt_report_days AS days_from_birth,
    r.cpt_seq_assay_id AS panel_id,
    r.cpt_number AS test_order,
    CASE
        WHEN r.cpt_report_post_death = 1 OR r.cpt_report_post_last_alive = 1 THEN TRUE
        WHEN u.censored_at_age_limit AND (r.dx_cpt_rep_days >= u.age_limit_day OR u.age_limit_day <= 0) THEN TRUE
        WHEN u.censored_at_age_limit AND r.dx_cpt_rep_days IS NOT NULL THEN FALSE
        WHEN r.cpt_report_post_death = 0 AND r.cpt_report_post_last_alive = 0 THEN FALSE
    END AS after_follow_up_end,
    r.cpt_n_ca_seq > 1 AS linked_to_other_cancers,
    r.cpt_seq_date
FROM index_cancer AS c
JOIN raw_report AS r ON r.record_id = c.record_id AND r.ca_seq = c.ca_seq
JOIN follow_up AS u ON u.record_id = c.record_id;

-- step: std_patients
-- One row per patient with a BPC Project cancer. The year of diagnosis is approximated as in the CRC BPC
-- plan, section 4c: from the birth year, assuming a mid-year birthday; if that is withheld, from the year
-- of the first sequencing report, assuming a mid-year report. Blank if both are withheld.
CREATE OR REPLACE TABLE std_patients AS
WITH first_ngs AS (
    SELECT t.*, g.n_genes
    FROM std_tests_ngs AS t
    LEFT JOIN panel_size AS g ON g.panel = t.panel_id
    QUALIFY row_number() OVER (PARTITION BY t.patient_id ORDER BY {first_order}) = 1
),
bpc_cancers AS (
    SELECT record_id, count(*) AS n
    FROM raw_cancer
    WHERE trim(redcap_ca_index) = 'Yes'
    GROUP BY record_id
),
other_same AS (     -- another primary of this cancer type: its pathology reports could be mistaken for the index cancer's
    SELECT DISTINCT record_id
    FROM raw_cancer_other
    WHERE list_contains($labels, trim(ca_type))
)
SELECT
    c.record_id AS patient_id,
    CASE WHEN u.censored_at_age_limit THEN 0 ELSE c.os_dx_status END AS vital_status,
    CASE WHEN u.censored_at_age_limit THEN u.age_limit_day ELSE c.tt_os_dx_days END AS days_to_death_or_last_followup,
    c.age_dx AS age_at_diagnosis,
    p.naaccr_sex_code AS sex,
    c.stage_dx AS stage_at_diagnosis,
    c.ca_hist_adeno_squamous AS histology,
    c.ca_d_site AS primary_site_code,
    CASE
        WHEN p.birth_year IS NOT NULL
            THEN CAST(floor(p.birth_year + 0.5 + c.dob_ca_dx_days / $days_per_year) AS INTEGER)
        ELSE CAST(floor(TRY_CAST(trim(f.cpt_seq_date) AS INTEGER) + 0.5 - f.days_to_entry_event / $days_per_year) AS INTEGER)
    END AS year_of_diagnosis,
    p.institution AS center,
    u.censored_at_age_limit,
    coalesce(b.n, 0) > 1 AS multiple_cohort_cancers,
    o.record_id IS NOT NULL AS other_primary_same_cancer{extras}
FROM index_cancer AS c
JOIN raw_patient AS p ON p.record_id = c.record_id
JOIN follow_up AS u ON u.record_id = c.record_id
LEFT JOIN first_ngs AS f ON f.patient_id = c.record_id
LEFT JOIN bpc_cancers AS b ON b.record_id = c.record_id
LEFT JOIN other_same AS o ON o.record_id = c.record_id;

-- step: pathology_result
-- One row per MSI or MMR result on a pathology report with at least one specimen of this cancer. UNPIVOT
-- turns the specimen columns into rows, so every specimen is checked. A report holds up to three MSI and
-- three MMR results, each with its OWN date: msi_result_2 goes with msi_prepaint_2, never with msi_prepaint.
CREATE OR REPLACE TABLE pathology_result AS
WITH path AS (
    SELECT rowid AS path_row, * FROM raw_pathology
),
of_this_cancer AS (
    SELECT DISTINCT path_row
    FROM (UNPIVOT path ON COLUMNS('^path_ca_type\d+$') INTO NAME specimen VALUE ca_type)
    WHERE list_contains($labels, trim(ca_type))
),
slot AS (
    SELECT path_row, record_id, 'MSI' AS test, 1 AS slot, msi_result AS label, msi_prepaint AS dob_result_days FROM path
    UNION ALL SELECT path_row, record_id, 'MSI', 2, msi_result_2, msi_prepaint_2 FROM path
    UNION ALL SELECT path_row, record_id, 'MSI', 3, msi_result_3, msi_prepaint_3 FROM path
    UNION ALL SELECT path_row, record_id, 'MMR', 1, mmr_result, mmr_prepaint FROM path
    UNION ALL SELECT path_row, record_id, 'MMR', 2, mmr_result_2, mmr_prepaint_2 FROM path
    UNION ALL SELECT path_row, record_id, 'MMR', 3, mmr_result_3, mmr_prepaint_3 FROM path
)
SELECT s.path_row, s.record_id, s.test, s.slot, trim(s.label) AS label, s.dob_result_days
FROM slot AS s
JOIN of_this_cancer AS k ON k.path_row = s.path_row
WHERE nullif(trim(s.label), '') IS NOT NULL;

-- step: std_tests_pathology
-- Each result is its own test, dated in days from diagnosis of the index cancer, and from birth as the source
-- dates it. "By the first report" compares the days from birth, as the CRC BPC plan does (section 4b): BPC's
-- days from birth and from diagnosis can differ by a day. A withheld date stays blank.
CREATE OR REPLACE TABLE std_tests_pathology AS
SELECT
    r.record_id AS patient_id,
    'pathology:' || r.path_row || ':' || r.test || ':' || r.slot AS test_id,
    'pathology' AS modality,
    'report' AS entry_event,
    r.dob_result_days - c.dob_ca_dx_days AS days_to_entry_event,
    r.dob_result_days AS days_from_birth
FROM pathology_result AS r
JOIN index_cancer AS c ON c.record_id = r.record_id;

-- step: std_results
-- Mutation calls and fusions on the index cancer's sequenced samples, and the MSI/MMR readings (bpc_label holds
-- the source's labels and the standard readings they map to). A fusion has a row for each fused gene the panel
-- reads (DFCI writes each with its own gene first: "A-B" under A, "B-A" under B), so a rule on a gene the panel
-- reads finds the fusion whichever partner the gene is, and a fusion whose partner is off the panel has one row;
-- its description is kept as the variant.
CREATE OR REPLACE TABLE std_results AS
SELECT m.Tumor_Sample_Barcode AS test_id, 'mutation' AS kind, m.Hugo_Symbol AS gene, m.HGVSp_Short AS variant,
       m.Variant_Classification AS variant_class, CAST(NULL AS VARCHAR) AS value
FROM raw_mutation AS m
WHERE m.Tumor_Sample_Barcode IN (SELECT test_id FROM std_tests_ngs)
UNION ALL
SELECT f.Tumor_Sample_Barcode, 'fusion', f.Hugo_Symbol, f.Fusion, NULL, NULL
FROM raw_fusion AS f
WHERE f.Tumor_Sample_Barcode IN (SELECT test_id FROM std_tests_ngs)
UNION ALL
SELECT 'pathology:' || r.path_row || ':' || r.test || ':' || r.slot, 'msi_mmr', NULL, NULL, NULL, l.reading
FROM pathology_result AS r
JOIN index_cancer AS c ON c.record_id = r.record_id
JOIN bpc_label AS l ON l.test = r.test AND l.label = r.label;

-- step: std_treatments
-- The index cancer's regimens, with their drug names as released when the cohort reads them.
CREATE OR REPLACE TABLE std_treatments AS
SELECT
    r.record_id AS patient_id,
    {regimen} AS regimen,
    r.regimen_number_within_cancer AS regimen_order,
    CAST(NULL AS INTEGER) AS line,
    r.dx_reg_start_int AS days_to_start
FROM index_cancer AS c
JOIN raw_regimen AS r ON r.record_id = c.record_id AND r.ca_seq = c.ca_seq;
