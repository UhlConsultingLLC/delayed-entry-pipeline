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
