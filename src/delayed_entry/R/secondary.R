# The secondary and sensitivity analyses that rerun the passes (design plan section 8). Ported from the CRC BPC
# project's R/secondary.R. Which patients each one keeps, and which groups it compares, now come from the Python
# side: a member column per analysis in the secondary table (a TRUE/FALSE column for a subset, the patient's group
# for a biomarker analysis), the treatment-clock cohorts, and the smallest-group rule's decisions in the spec.
#
#   subset             passes 1 and 2 on the patients it keeps
#   treatment_clock    passes 1 and 2, with time counted from the first-line start
#   biomarker          passes 1, 2 and 3a, comparing the biomarker's groups (and its variants)
#   adjusted_variant   pass 3a with one variable as strata instead of HRs

# The analysis cohort's patients in the analysis, with the secondary table's columns.
with_secondary <- function(cohort, secondary) {
  x <- merge(cohort[cohort$in_analysis, ], secondary, by = "patient_id")
  if (nrow(x) != sum(cohort$in_analysis)) stop("the secondary table doesn't match the analysis cohort")
  x
}

# One analysis's patients: a TRUE/FALSE member column keeps the TRUE ones; a group column keeps the patients in
# one of the grouping's levels.
subset_for <- function(x, column, grouping = NULL) {
  keep <- if (is.logical(x[[column]])) x[[column]] else x[[column]] %in% grouping$levels
  x[!is.na(keep) & keep, ]
}

# Stops unless R's counts for one analysis are the ones the Python side decided the smallest-group rule from:
# both must have chosen the same patients. Returns the decisions.
check_counts <- function(analysis, x, grouping, d1) {
  n <- table(factor(x[[grouping$column]][x$in_analysis], levels = grouping$levels))
  ours <- data.frame(comparison = contrast_name(grouping$levels[-1], grouping),
                     n_reference = as.integer(n[[grouping$levels[1]]]), n_group = as.integer(n[grouping$levels[-1]]))
  mine <- d1[d1$analysis == analysis, ]
  theirs <- data.frame(comparison = mine$comparison, n_reference = as.integer(mine$n_reference),
                       n_group = as.integer(mine$n_group))
  if (!isTRUE(all.equal(ours, theirs, check.attributes = FALSE))) {
    stop(analysis, ": R's patients differ from the ones the smallest-group rule counted")
  }
  mine[, c("comparison", "decision")]
}

# Passes 1 and 2 on one analysis's patients: every estimate, the contrasts with intervals and labels, and every
# RMST with its interval.
run_passes <- function(x, grouping, d1, described = character()) {
  est <- pass_estimates(x, grouping)
  boot <- bootstrap_bias(x, grouping = grouping)
  list(estimates = est,
       bias = bias_table(bias_estimates(est, grouping), boot, d1, grouping, described),
       rmst = rmst_table(est, boot, grouping))
}

# Pass 3a on one analysis's patients (the knots are set from these patients). With imputation, the analysis's own
# patients are imputed (R/impute.R).
run_adjusted <- function(x, covariates, grouping, d1, stratify_by = NULL, described = character()) {
  if (imputing()) {
    r <- adjusted_mi(x, covariates, grouping, stratify_by)
    return(list(estimates = r$estimates, bias = bias_table(r$bias, r$boot, d1, grouping, described)))
  }
  frame <- model_frame(x, covariates, grouping = grouping)
  est <- adjusted_estimates(frame, grouping, stratify_by)
  boot <- bootstrap_adjusted(frame, grouping = grouping, stratify_by = stratify_by)
  list(estimates = est, bias = bias_table(adjusted_bias(est), boot, d1, grouping, described))
}

# Every secondary analysis the spec runs. Returns a list of results by name ("S1", "S1 3a", ...) and the warnings.
run_secondary <- function(con, run_dir, cohort, covariates) {
  x <- with_secondary(cohort, read_table(con, run_dir, "secondary.parquet"))
  d1 <- as.data.frame(SPEC$d1)
  results <- list()
  warned <- list()
  keep <- function(name, expr) {
    r <- with_warnings(expr)
    results[[name]] <<- r$value
    if (length(r$warnings)) warned[[name]] <<- r$warnings
  }
  for (a in SPEC$secondary) {
    grouping <- if (is.null(a$grouping)) GENOTYPE else grouping_from(a$grouping)
    described <- named_text(a$described)
    if (a$kind %in% c("biomarker", "biomarker_variant")) {
      p <- subset_for(x, a$member_column, grouping)
      decisions <- check_counts(a$id, p, grouping, d1)
      keep(a$id, run_passes(p, grouping, decisions, described))
      keep(paste(a$id, "3a"), run_adjusted(p, covariates, grouping, decisions, described = described))
    } else if (a$kind == "subset") {
      p <- subset_for(x, a$member_column)
      keep(a$id, run_passes(p, grouping, check_counts(a$id, p, grouping, d1)))
    } else if (a$kind == "treatment_clock") {
      p <- read_table(con, run_dir, paste0("clock_", a$id, ".parquet"))
      keep(a$id, run_passes(p, grouping, check_counts(a$id, p, grouping, d1)))
    } else if (a$kind == "adjusted_variant") {
      keep(a$id, run_adjusted(cohort, covariates, grouping, d1[d1$analysis == "primary", ], stratify_by = a$stratify_by))
    }
  }
  list(results = results, warnings = warned)
}

# One part of every analysis's results, with the analysis's name in front.
tag <- function(results, part) {
  named <- names(results)[!vapply(results, function(r) is.null(r[[part]]), TRUE)]
  if (!length(named)) return(NULL)
  do.call(rbind, lapply(named, function(a) data.frame(analysis = a, results[[a]][[part]])))
}
