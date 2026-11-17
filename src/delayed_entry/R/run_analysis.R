# One step of the run's R part. The Python side (delayed_entry.run) writes the run folder first: run_spec.json
# (every value, from the locked settings), analysis_cohort.parquet, covariates.parquet, secondary.parquet and any
# clock_<id>.parquet; then it calls each step, in order, as
#
#   Rterm.exe --no-echo --no-restore -f run_analysis.R --args <this folder> <run folder> <step>
#
# (or Rscript run_analysis.R <this folder> <run folder> <step>). Steps, ported from the CRC BPC project's scripts:
#   primary       passes 1 and 2, the contrasts, the stability rule and the curves      (scripts/r_passes.R)
#   adjusted      pass 3a                                                                (scripts/r_adjusted.R)
#   weighted      pass 3b and its balance                                                (scripts/r_weighted.R)
#   assumptions   proportional hazards and quasi-independence                            (scripts/r_assumptions.R)
#   secondary     the secondary and sensitivity analyses that rerun the passes           (scripts/r_secondary.R)
# Each writes its results as r_*.parquet in the run folder, and r_<step>.json: what it did, how long it took, its
# warnings, and the versions of R and the packages. The run folder holds row-level files: it must be local.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: run_analysis.R <R folder> <run folder> <step>")
lib <- Sys.getenv("DELAYED_ENTRY_R_LIBRARY")
if (nzchar(lib)) .libPaths(c(strsplit(lib, .Platform$path.sep, fixed = TRUE)[[1]], .libPaths()))
r_dir <- args[1]
step <- args[3]
for (f in c("common.R", "passes.R", "bias.R", "stability.R", "adjusted.R", "weighted.R", "impute.R",
            "assumptions.R", "secondary.R")) source(file.path(r_dir, f))
run_dir <- check_not_cloud(normalizePath(args[2], winslash = "/", mustWork = TRUE))
load_spec(run_dir)
started <- Sys.time()
con <- connect()
cohort <- read_table(con, run_dir, "analysis_cohort.parquet")
d1 <- as.data.frame(SPEC$d1)
d1_primary <- d1[d1$analysis == "primary", ]
described <- named_text(SPEC$primary_described)
warned <- list()
done <- list(step = step)

if (step == "primary") {
  est <- pass_estimates(cohort)
  bias <- bias_estimates(est)
  boot <- bootstrap_bias(cohort)
  stab <- stability(cohort)
  table <- bias_table(bias, boot, d1_primary, described = described, unstable = unstable_groups(stab))
  cond <- conditional_estimates(cohort, stab)
  cond_bias <- conditional_bias(cond)
  cond_table <- conditional_table(cond_bias, bootstrap_conditional(cohort, stab), stab, described)
  write_table(con, rbind(est, bias, stability_rows(stab), cond, cond_bias), run_dir, "r_passes.parquet")
  write_table(con, table, run_dir, "r_bias.parquet")
  write_table(con, boot, run_dir, "r_bootstrap.parquet")
  write_table(con, stab, run_dir, "r_stability.parquet")
  write_table(con, cond_table, run_dir, "r_conditional.parquet")
  write_table(con, rmst_table(est, boot), run_dir, "r_rmst.parquet")
  write_table(con, rbind(km_curves(cohort, "naive"), km_curves(cohort, "corrected")), run_dir, "r_curves.parquet")
  done$packages <- as.list(versions(c("survival", "duckdb")))
} else if (step == "adjusted" && imputing()) {
  covariates <- read_table(con, run_dir, "covariates.parquet")
  r <- with_warnings(adjusted_mi(cohort, covariates))
  if (length(r$warnings)) warned$adjusted <- r$warnings
  write_table(con, r$value$estimates, run_dir, "r_adjusted.parquet")
  write_table(con, bias_table(r$value$bias, r$value$boot, d1_primary), run_dir, "r_adjusted_bias.parquet")
  done$model <- list(patients = r$value$patients, imputations = imputations())
  done$packages <- as.list(versions(c("survival", "mice", "duckdb")))
} else if (step == "adjusted") {
  covariates <- read_table(con, run_dir, "covariates.parquet")
  frame <- model_frame(cohort, covariates)
  est <- adjusted_estimates(frame)
  boot <- bootstrap_adjusted(frame)
  write_table(con, est, run_dir, "r_adjusted.parquet")
  write_table(con, bias_table(adjusted_bias(est), boot, d1_primary), run_dir, "r_adjusted_bias.parquet")
  done$model <- list(patients = nrow(frame), left_out = attr(frame, "left_out"), knots = attr(frame, "knots"))
  done$packages <- as.list(versions(c("survival", "duckdb")))
} else if (step == "weighted") {
  covariates <- read_table(con, run_dir, "covariates.parquet")
  unstable <- unstable_groups(stability(cohort))
  if (imputing()) {
    r <- with_warnings(weighted_mi(cohort, covariates))
    mi_fit <- r$value
    r$value <- mi_fit$boot
    est <- mi_fit$estimates
    bal <- mi_fit$balance
    n_model <- mi_fit$patients
    write_table(con, mi_fit$by_imputation, run_dir, "r_balance_by_imputation.parquet")
    write_table(con, mi_fit$completed, run_dir, "r_imputed.parquet")
  } else {
    frame <- model_frame(cohort, covariates)
    w <- iptw(frame)
    bal <- balance(w)
    est <- weighted_estimates(frame, w$weights)
    r <- with_warnings(bootstrap_weighted(frame))
    n_model <- nrow(frame)
  }
  if (length(r$warnings)) warned$weighted <- r$warnings
  tables <- weighted_tables(est, r$value, d1_primary, bal, described, unstable)
  hrs <- tables$hrs
  hrs$resamples_missing <- NA_real_
  hrs$note <- tables$hrs$note
  write_table(con, bal, run_dir, "r_balance.parquet")
  write_table(con, rbind(tables$curves, hrs[, names(tables$curves)]), run_dir, "r_weighted.parquet")
  write_table(con, tables$bias, run_dir, "r_weighted_bias.parquet")
  done$balanced <- balanced(bal)
  done$weights <- list(patients = n_model, imputations = if (imputing()) imputations() else 0L)
  done$packages <- as.list(versions(c("survival", "WeightIt", "cobalt", if (imputing()) "mice", "duckdb")))
} else if (step == "assumptions") {
  if (isTRUE(SPEC$run_ph_check)) {
    ph <- ph_test(cohort)
    rejected <- ph_rejected(ph)
    write_table(con, data.frame(ph, rejected = rejected), run_dir, "r_ph.parquet")
    if (rejected) write_table(con, split_hrs(cohort), run_dir, "r_ph_split.parquet")
    done$ph_rejected <- rejected
  }
  if (isTRUE(SPEC$run_quasi_independence)) {
    tau <- quasi_independence(cohort)
    tr <- transformation_medians(cohort)
    sa6 <- merge(tau, tr, by = "group", all.x = TRUE)
    sa6 <- sa6[order(match(sa6$estimator, SPEC$tau_estimators), match(sa6$group, c(ALL_PATIENTS, GENOTYPE$levels))), ]
    write_table(con, sa6, run_dir, "r_tau.parquet")
  }
  done$packages <- as.list(versions(c("survival", "tranSurv", "duckdb")))
} else if (step == "secondary") {
  covariates <- read_table(con, run_dir, "covariates.parquet")
  out <- run_secondary(con, run_dir, cohort, covariates)
  warned <- out$warnings
  for (part in c("estimates", "bias", "rmst")) {
    t <- tag(out$results, part)
    if (!is.null(t)) write_table(con, t, run_dir, sprintf("r_secondary_%s.parquet", part))
  }
  done$analyses <- names(out$results)
  done$packages <- as.list(versions(c("survival", if (imputing()) "mice", "duckdb")))
} else {
  stop("unknown step: ", step)
}
DBI::dbDisconnect(con, shutdown = TRUE)
done$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
done$resamples <- BOOT_RESAMPLES
done$seed <- BOOT_SEED
done$warnings <- lapply(warned, function(w) {
  u <- unique(w)
  stats::setNames(as.list(vapply(u, function(m) sum(w == m), numeric(1))), gsub("\\s+", " ", u))
})
jsonlite::write_json(done, file.path(run_dir, sprintf("r_%s.json", step)), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(sprintf("R step %s: done in %.0f seconds\n", step, done$seconds))
