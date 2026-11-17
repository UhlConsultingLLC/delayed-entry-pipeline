-- Link checks on the raw GENIE BPC tables (made by delayed_entry.adapters.genie_bpc.raw.load). Copied from
-- the CRC BPC project (sql/01_link_checks.sql). Each check returns the rows that break one link; zero rows
-- means the check passes. raw.run_checks() runs them all and counts the rows.
--
-- Most checks are anti-joins: LEFT JOIN, then keep the rows that found no partner
-- (the right-hand key IS NULL). Those are the rows the link leaves stranded.

-- check: patient_id_unique
-- One row per patient in the patient file.
SELECT record_id, count(*) AS n_rows
FROM raw_patient
GROUP BY record_id
HAVING count(*) > 1;

-- check: patient_has_cancer
-- Every patient has at least one BPC Project cancer.
SELECT p.record_id
FROM raw_patient AS p
LEFT JOIN raw_cancer AS c ON c.record_id = p.record_id
WHERE c.record_id IS NULL;

-- check: cancer_has_patient
-- Every cancer belongs to a patient in the patient file.
SELECT c.record_id, c.ca_seq
FROM raw_cancer AS c
LEFT JOIN raw_patient AS p ON p.record_id = c.record_id
WHERE p.record_id IS NULL;

-- check: cancer_key_unique
-- A patient's cancers have different numbers, so record_id + ca_seq names one cancer.
SELECT record_id, ca_seq, count(*) AS n_rows
FROM raw_cancer
GROUP BY record_id, ca_seq
HAVING count(*) > 1;

-- check: report_has_cancer
-- Every report points at a cancer: the same record_id AND the same ca_seq.
SELECT r.record_id, r.ca_seq, r.cpt_genie_sample_id
FROM raw_report AS r
LEFT JOIN raw_cancer AS c ON c.record_id = r.record_id AND c.ca_seq = r.ca_seq
WHERE c.record_id IS NULL;

-- check: report_sample_unique
-- One row per sequenced sample and cancer. A test linked to two cancers may be listed
-- once for each (the plan keeps and counts those), but never twice for the same cancer.
SELECT cpt_genie_sample_id, ca_seq, count(*) AS n_rows
FROM raw_report
GROUP BY cpt_genie_sample_id, ca_seq
HAVING count(*) > 1;

-- check: report_number_unique
-- A cancer's reports have different test numbers. The first-report rule breaks its last
-- tie on cpt_number, so this is what makes its choice the same every time it runs.
SELECT record_id, ca_seq, cpt_number, count(*) AS n_rows
FROM raw_report
GROUP BY record_id, ca_seq, cpt_number
HAVING count(*) > 1;

-- check: report_has_panel
-- Every report names a panel that has a gene-panel file.
SELECT r.cpt_genie_sample_id, r.cpt_seq_assay_id
FROM raw_report AS r
LEFT JOIN (SELECT DISTINCT panel FROM raw_panel_gene) AS g ON g.panel = r.cpt_seq_assay_id
WHERE g.panel IS NULL;

-- check: mutation_has_report
-- Every mutation call belongs to a sample in the report file.
SELECT m.Tumor_Sample_Barcode, count(*) AS n_calls
FROM raw_mutation AS m
LEFT JOIN raw_report AS r ON r.cpt_genie_sample_id = m.Tumor_Sample_Barcode
WHERE r.cpt_genie_sample_id IS NULL
GROUP BY m.Tumor_Sample_Barcode;

-- check: call_gene_on_panel
-- Every call is in a gene its sample's panel reads.
SELECT m.Tumor_Sample_Barcode, m.Hugo_Symbol, r.cpt_seq_assay_id
FROM raw_mutation AS m
JOIN raw_report AS r ON r.cpt_genie_sample_id = m.Tumor_Sample_Barcode
LEFT JOIN raw_panel_gene AS g ON g.panel = r.cpt_seq_assay_id AND g.gene = m.Hugo_Symbol
WHERE g.gene IS NULL;

-- check: fusion_has_report
-- Every fusion row belongs to a sample in the report file (a cohort that reads no fusion file has none).
SELECT f.Tumor_Sample_Barcode, count(*) AS n_rows
FROM raw_fusion AS f
LEFT JOIN raw_report AS r ON r.cpt_genie_sample_id = f.Tumor_Sample_Barcode
WHERE r.cpt_genie_sample_id IS NULL
GROUP BY f.Tumor_Sample_Barcode;

-- check: other_cancer_has_patient
-- Every non-index cancer belongs to a patient in the patient file.
SELECT o.record_id, o.ca_seq
FROM raw_cancer_other AS o
LEFT JOIN raw_patient AS p ON p.record_id = o.record_id
WHERE p.record_id IS NULL;

-- check: other_cancer_not_index
-- The non-index file holds no BPC Project cancer, and none of its cancers shares a number with one. Any index
-- flag other than "No" fails, except a blank one, which is the next check's case.
SELECT o.record_id, o.ca_seq
FROM raw_cancer_other AS o
LEFT JOIN raw_cancer AS c ON c.record_id = o.record_id AND c.ca_seq = o.ca_seq
WHERE (nullif(trim(o.redcap_ca_index), '') IS NOT NULL AND trim(o.redcap_ca_index) <> 'No')
   OR c.record_id IS NOT NULL;

-- check: other_cancer_flag_blank
-- A non-index cancer whose index flag is blank. The release has some (seen 2026-11-11, with a blank cancer type
-- too). Only the index file defines BPC Project cancers, so such a row can't add one. The one rule that reads this
-- file, S1's other colorectal primary, is covered by check_download's check on colorectal site codes.
SELECT o.record_id, o.ca_seq
FROM raw_cancer_other AS o
WHERE nullif(trim(o.redcap_ca_index), '') IS NULL;

-- check: pathology_has_patient
-- Every pathology report belongs to a patient in the patient file.
SELECT a.record_id, a.path_proc_number, a.path_rep_number
FROM raw_pathology AS a
LEFT JOIN raw_patient AS p ON p.record_id = a.record_id
WHERE p.record_id IS NULL;

-- check: pathology_key_unique
-- record_id + procedure number + report number names one pathology report.
SELECT record_id, path_proc_number, path_rep_number, count(*) AS n_rows
FROM raw_pathology
GROUP BY record_id, path_proc_number, path_rep_number
HAVING count(*) > 1;

-- check: regimen_number_unique
-- A cancer's regimens have different numbers, so S2 can take the first one (plan section 7).
SELECT record_id, ca_seq, regimen_number_within_cancer, count(*) AS n_rows
FROM raw_regimen
GROUP BY record_id, ca_seq, regimen_number_within_cancer
HAVING count(*) > 1;

-- check: report_seq_year_is_a_year
-- The year of sequencing is stored as text. Where it is given, it must be a four-digit year,
-- or the year of diagnosis worked out from it (sql/04_covariates.sql) would quietly go blank.
SELECT record_id, cpt_number, cpt_seq_date
FROM raw_report
WHERE nullif(trim(cpt_seq_date), '') IS NOT NULL
  AND NOT regexp_full_match(trim(cpt_seq_date), '[0-9]{4}');
