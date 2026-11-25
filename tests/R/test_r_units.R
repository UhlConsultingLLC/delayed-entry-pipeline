# Unit checks of the run's R code on made-up numbers. Run by tests/test_run.py when R is configured:
#   Rterm.exe --no-echo --no-restore -f tests/R/test_r_units.R --args <the package's R folder> <a temporary folder>
args <- commandArgs(trailingOnly = TRUE)
lib <- Sys.getenv("DELAYED_ENTRY_R_LIBRARY")
if (nzchar(lib)) .libPaths(c(strsplit(lib, .Platform$path.sep, fixed = TRUE)[[1]], .libPaths()))
for (f in c("common.R", "passes.R", "bias.R", "stability.R", "adjusted.R", "weighted.R", "impute.R", "assumptions.R",
            "secondary.R")) {
  source(file.path(args[1], f))
}
spec <- list(days_per_year = 365.25, survival_at_years = c(1, 2), horizon_years = 1, median_limit_years = 2,
             all_patients = "All patients", resamples = 5, seed = 7, interval_level = 0.9,
             bars = list(median_diff = 3, rmst_diff = 2, hr_ratio = 0.1, rmst_contrast_change = 2), stable_risk_set = 3,
             genotype = list(column = "genotype", levels = c("Ref", "A"), reference = "Ref"),
             labels = list(clearly_material = "clearly material", material = "material", not_material = "not material"),
             adjustment = list(fields = c("age", "trait"), forms = c("linear", "categorical"),
                               levels = list(trait = c("a", "b", "c")), knots = setNames(list(), character())),
             missing_data = list(rule = "multiple imputation", imputations = 3,
                                 model = list(engine = "chained_equations", outcome_terms = list("event"),
                                              include_group = TRUE, pmm_donors = 3, iterations = 4,
                                              methods = list(linear = "pmm", spline = "pmm", categorical = "polyreg",
                                                             bands = "polyreg"))),
             secondary = list())
jsonlite::write_json(spec, file.path(args[2], "run_spec.json"), auto_unbox = TRUE, digits = NA)
load_spec(args[2])
ok <- 0
check <- function(what, cond) {
  if (!isTRUE(cond)) stop("FAILED: ", what)
  ok <<- ok + 1
}

# The values come from the spec
check("the interval's percentiles come from the level", isTRUE(all.equal(PROBS, c(0.05, 0.95))))
check("the horizon in days", RMST_DAYS == 365.25)

# Labels by the bars
check("clearly material: estimate and whole interval past the bar", label_bias(3, 2.5, 4, 2) == "clearly material")
check("material: the interval crosses the bar", label_bias(3, 1, 4, 2) == "material")
check("not material: the estimate falls short", label_bias(1.5, 1, 4, 2) == "not material")
check("a ratio is judged against 1", label_bias(1.15, 0.93, 1.39, 0.1, center = 1) == "material")
check("at the bar despite floating point", label_bias(0.9, 0.85, 0.9, 0.1, center = 1) == "clearly material")

# The spline: with 3 knots, the CRC BPC project's one curved term exactly
x <- c(40, 55, 61, 70, 82)
k <- c(45, 60, 75)
crc <- function(x, k) {
  pos3 <- function(v) pmax(v, 0)^3
  (pos3(x - k[1]) - pos3(x - k[2]) * (k[3] - k[1]) / (k[3] - k[2]) + pos3(x - k[3]) * (k[2] - k[1]) / (k[3] - k[2])) /
    (k[3] - k[1])^2
}
check("rcs with 3 knots is the CRC term", identical(as.vector(rcs_terms(x, k)), crc(x, k)))
check("rcs with 4 knots has 2 curved terms", ncol(rcs_terms(x, c(45, 55, 65, 75))) == 2)

# The risk set and t*
d <- data.frame(entry_day = c(0, 10, 20, 30, 40), exit_day = c(100, 100, 100, 15, 100), event = c(1, 0, 1, 1, 0))
check("at risk in (entry, exit]", identical(at_risk(d, "corrected", c(5, 25)), c(1, 3)))
check("t* is the first entry after which the risk set reaches the size", t_star(d, 3) == 20)
check("no t* when the risk set never reaches the size", is.na(t_star(d, 6)))

# Percentile intervals: none when a resample has no value
boot <- data.frame(quantity = "q", group = "g", value = c(1:19, NA))
check("no interval when a resample has no value", is.na(percentile_ci(boot)$lower) && percentile_ci(boot)$resamples_missing == 1)
# Rubin's rules: mean 2, within 1, between 1, total 1 + (4/3) * 1, Rubin's df (m - 1)(1 + 1 / (4/3))^2 = 6.125
r <- rubin(c(1, 2, 3), c(1, 1, 1))
check("Rubin: the pooled estimate is the mean", r[["estimate"]] == 2)
check("Rubin: the total variance", isTRUE(all.equal(r[["se"]]^2, 1 + 4 / 3)))
check("Rubin: the degrees of freedom", isTRUE(all.equal(r[["df"]], 6.125)))
check("Rubin: the interval on t", isTRUE(all.equal(r[["upper"]], 2 + stats::qt(0.95, 6.125) * sqrt(7 / 3))))
check("Rubin: imputations that agree give a normal interval", is.infinite(rubin(c(2, 2), c(1, 1))[["df"]]))

# The cumulative hazard with delayed entry: one death at day 50 among the 2 who had entered by then
d <- data.frame(entry_day = c(0, 10, 60), exit_day = c(50, 100, 100), event = c(1, 0, 0))
check("Nelson-Aalen with delayed entry", isTRUE(all.equal(cumulative_hazard(d, "corrected"), c(0.5, 0.5, 0.5))))
check("Nelson-Aalen from diagnosis", isTRUE(all.equal(cumulative_hazard(d, "naive"), c(1 / 3, 1 / 3, 1 / 3))))
# The pipeline's own imputation engine, on made-up patients: age follows the group, the trait is random
set.seed(1)
n <- 120
x <- data.frame(patient_id = sprintf("P%03d", 1:n), group = factor(rep(c("Ref", "A"), each = n / 2), levels = c("Ref", "A")),
                entry_day = 0, exit_day = stats::runif(n, 50, 500), event = stats::rbinom(n, 1, 0.5))
x$age <- ifelse(x$group == "A", 70, 50) + stats::rnorm(n)
x$trait <- factor(sample(c("a", "b", "c"), n, replace = TRUE), levels = c("a", "b", "c"))
gaps <- c(3, 17, 64, 90, 111)
x$age[gaps] <- NA
x$trait[c(5, 66)] <- NA
done <- impute(x, 3, seed = 11)
check("three completed copies", length(done) == 3)
check("every gap filled", all(vapply(done, function(y) !anyNA(y$age) && !anyNA(y$trait), logical(1))))
check("the observed values are untouched", all(vapply(done, function(y) identical(y$age[-gaps], x$age[-gaps]), logical(1))))
check("pmm imputes observed values only", all(done[[1]]$age[gaps] %in% x$age[-gaps]))
check("pmm follows a strong predictor", all(abs(done[[1]]$age[gaps] - ifelse(x$group[gaps] == "A", 70, 50)) < 5))
check("categories stay within their levels", identical(levels(done[[1]]$trait), c("a", "b", "c")))
check("the same seed gives the same imputations", identical(impute(x, 3, seed = 11), done))
check("another seed gives others", !identical(impute(x, 3, seed = 12)[[1]]$age, done[[1]]$age))

# The guard against perfect prediction: in group A every observed trait is "a". The augmenting records keep the other
# categories possible there (on these numbers about 1 in 100; a plain fit leaves a few in a million)
set.seed(5)
cell <- data.frame(group = factor(rep(c("Ref", "A"), each = 30), levels = c("Ref", "A")), event = stats::rbinom(60, 1, 0.5))
cell_trait <- factor(c(sample(c("a", "b", "c"), 30, replace = TRUE), rep("a", 30)), levels = c("a", "b", "c"))
cell_trait[c(31, 32)] <- NA
cell_miss <- is.na(cell_trait)
cell_x <- design(cell)[, -1, drop = FALSE]
pr <- polyreg_probs(droplevels(cell_trait[!cell_miss]), cell_x, cell_miss)
check("a row of probabilities for each missing case, summing to 1", nrow(pr) == 2 && isTRUE(all.equal(rowSums(pr), c(1, 1))))
check("a perfectly predicted cell keeps the other categories possible", all(pr[, 2:3] > 0.001) && all(pr[, 1] > 0.9))
one_gap <- cell_miss
one_gap[32] <- FALSE
one_trait <- cell_trait
one_trait[32] <- "a"
pr1 <- polyreg_probs(droplevels(one_trait[!one_gap]), cell_x, one_gap)
check("and with one missing case (where mice leaves the records out)", nrow(pr1) == 1 && all(pr1[, 2:3] > 0.001))
cell_two <- factor(ifelse(cell_trait == "a", "a", "b"), levels = c("a", "b"))
pr2 <- polyreg_probs(droplevels(cell_two[!cell_miss]), cell_x, cell_miss)
check("and with two levels", ncol(pr2) == 2 && all(pr2[, 2] > 0.001) && isTRUE(all.equal(rowSums(pr2), c(1, 1))))
set.seed(9)
drawn <- unlist(lapply(1:400, function(i) as.character(draw_logreg(cell_two, design(cell), cell_miss))))
check("the logistic draw is guarded too: the cell's other level is drawn, but rarely",
      any(drawn == "b") && mean(drawn == "b") < 0.1)
pr0 <- polyreg_probs(droplevels(cell_trait[!cell_miss]), cell_x[, 0, drop = FALSE], cell_miss)
check("with no predictors, the observed shares",
      isTRUE(all.equal(pr0[1, ], as.vector(table(cell_trait[!cell_miss])) / sum(!cell_miss))))
check("multinomial draws on a sparse cell stay within the levels",
      all(vapply(1:20, function(i) all(draw_polyreg(cell_trait, cell, cell_miss) %in% c("a", "b", "c")), logical(1))))

# SA6 on a group too small for its estimators: no value and a note saying why, never a stop
tiny <- data.frame(entry_day = c(10, 50, 200), exit_day = c(400, 300, 900), event = c(0, 1, 0))
tr <- transformation_row(tiny, "A")
check("a group too small for the transformation model gets no values", is.na(tr$median_transformation) && is.na(tr$a))
check("and a warning that says so", startsWith(tr$warning, "transformation model not fitted:"))
kr <- kendall_row(tiny, "no such estimator", "A")
check("an estimator that stops gives no tau and no flag", is.na(kr$tau) && is.na(kr$flagged))
check("and says why", startsWith(kr$tau_note, "the estimator stopped:"))
kr <- kendall_row(tiny[0, ], "MB", "A")
check("no standard error: neither flagged nor cleared", is.na(kr$flagged) && kr$tau_note == "no standard error, so no interval")
g <- data.frame(entry_day = c(0, 5, 12, 20, 33, 41, 58, 64, 77, 90, 102, 130),
                exit_day = c(300, 95, 410, 260, 120, 515, 200, 380, 160, 610, 240, 450),
                event = c(1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0))
kr <- kendall_row(g, "MB", "A")
check("a group with an interval is flagged or cleared, with no note", identical(kr$flagged, FALSE) && kr$tau_note == "")
check("and the transformation model is fitted", !is.na(transformation_row(g, "A")$median_transformation))

# A group the settings describe only: its comparison says so, not "smallest group"; its own estimates keep labels
LABELS <- c(LABELS, list(not_compared = "not compared (smallest group)", described = "described, not labeled (%s)",
                         all_patients = "all patients", median_not_reached = "not reached", unstable = "unstable"))
set.seed(3)
b <- data.frame(quantity = "rmst_contrast_change", group = c("A", "A vs Ref"), value = c(3, 3))
bt <- data.frame(quantity = "rmst_contrast_change", group = rep(c("A", "A vs Ref"), each = 20), value = 3 + rnorm(40, 0, 0.1))
d1_short <- data.frame(comparison = "A vs Ref", decision = "describe")
SPEC$settings_described <- list(`A vs Ref` = "settings")
lab <- bias_table(b, bt, d1_short)
check("a comparison the settings describe says so", lab$label[lab$group == "A vs Ref"] == "described, not labeled (settings)")
check("and the group's own estimate keeps its label", lab$label[lab$group == "A"] == "clearly material")
SPEC$settings_described <- NULL
lab <- bias_table(b, bt, d1_short)
check("otherwise a comparison not compared is the smallest group's", lab$label[lab$group == "A vs Ref"] == "not compared (smallest group)")

# The mice engine, where mice loads (on Windows, Smart App Control can block a DLL it needs)
if (isTRUE(tryCatch(requireNamespace("mice", quietly = TRUE), error = function(e) FALSE))) {
  SPEC$missing_data$model$engine <- "mice"
  by_mice <- impute(x, 2, seed = 11)
  check("mice fills every gap", all(vapply(by_mice, function(y) !anyNA(y$age) && !anyNA(y$trait), logical(1))))
  check("mice's pmm imputes observed values only", all(by_mice[[1]]$age[gaps] %in% x$age[-gaps]))
  SPEC$missing_data$model$engine <- "chained_equations"
} else {
  cat("mice doesn't load here: its engine's checks skipped\n")
}
cat("R unit checks:", ok, "passed\n")
