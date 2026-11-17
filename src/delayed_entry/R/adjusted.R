# Pass 3a (design plan section 8): Cox models adjusted for the settings' adjustment variables, fitted naive and
# corrected, so the adjustment is the same on both sides. The adjusted bias in each comparison is
# HR(naive) / HR(corrected) within 3a, with a bootstrap interval and a label.
#
# Ported from the CRC BPC project's R/adjusted.R. Its fixed covariates (age, sex, center, tumor location, year
# band), their levels and the spline's knots now come from the spec's `adjustment`: each variable's form, levels
# (the reference first) and knots. The Python side codes the values (levels merged, kept or missing; year bands
# merged when small), so a value outside the levels is already missing here.

# Harrell's restricted cubic spline with k knots: k - 2 curved terms, each a cubic between the knots and a
# straight line beyond the outer two, divided by (last knot - first knot)^2 as rms does. With 3 knots this is the
# CRC BPC project's one curved term exactly.
rcs_terms <- function(x, knots) {
  k <- sort(knots)
  n <- length(k)
  pos3 <- function(v) pmax(v, 0)^3
  m <- vapply(seq_len(n - 2), function(j) {
    (pos3(x - k[j]) - pos3(x - k[n - 1]) * (k[n] - k[j]) / (k[n] - k[n - 1]) +
       pos3(x - k[n]) * (k[n - 1] - k[j]) / (k[n] - k[n - 1])) / (k[n] - k[1])^2
  }, numeric(length(x)))
  matrix(m, nrow = length(x))
}

# The knots: the spec's percentiles of the variable among the patients in the model (type 7, as rms).
spline_knots <- function(x, percentiles) unname(stats::quantile(x, percentiles / 100, type = 7))

adjustment <- function() SPEC$adjustment

# A name as it may appear in a formula.
formula_name <- function(f) if (make.names(f) == f) f else paste0("`", f, "`")

# The model's terms, in the settings' order: a spline is the variable and its curved terms; a categorical or
# banded variable is a factor. `stratify_by` names a variable that gets its own baseline hazard instead of HRs.
model_terms <- function(stratify_by = NULL) {
  adj <- adjustment()
  out <- character()
  for (i in seq_along(adj$fields)) {
    f <- adj$fields[i]
    form <- adj$forms[i]
    if (form == "spline") {
      out <- c(out, formula_name(f), formula_name(paste0(f, "_curve", seq_len(length(adj$knots[[f]]) - 2))))
    } else if (!is.null(stratify_by) && identical(f, stratify_by)) {
      out <- c(out, sprintf("strata(%s)", formula_name(f)))
    } else {
      out <- c(out, formula_name(f))
    }
  }
  out
}

# The patients pass 3a models, with their variables coded as the spec says, and their group (as in
# analysis_frame) in the column `group`. Patients with no time at risk are left out, and so is anyone missing a
# variable (the complete-case rule, which the Python side checked applies). The knots are set from the patients in
# the model unless given.
model_frame <- function(cohort, covariates, knots = NULL, grouping = GENOTYPE) {
  adj <- adjustment()
  fields <- adj$fields
  x <- model_data(cohort, covariates, grouping)
  keep <- stats::complete.cases(x[, fields, drop = FALSE])
  x <- x[keep, ]
  x$group <- droplevels(x$group)
  for (i in seq_along(fields)) {
    if (adj$forms[i] %in% c("categorical", "bands")) x[[fields[i]]] <- droplevels(x[[fields[i]]])
  }
  splines <- fields[adj$forms == "spline"]
  if (is.null(knots)) {
    knots <- lapply(stats::setNames(splines, splines), function(f) spline_knots(x[[f]], adj$knots[[f]]))
  }
  x <- add_splines(x, knots)
  attr(x, "left_out") <- sum(!keep)
  x
}

# Everyone in the analysis with their adjustment variables, before any is imputed or anyone left out: the categorical
# and banded variables as factors with the spec's levels (the reference first).
model_data <- function(cohort, covariates, grouping = GENOTYPE) {
  adj <- adjustment()
  fields <- adj$fields
  x <- merge(analysis_frame(cohort, grouping), covariates[, c("patient_id", fields)], by = "patient_id")
  for (i in seq_along(fields)) {
    f <- fields[i]
    if (adj$forms[i] %in% c("categorical", "bands")) x[[f]] <- factor(x[[f]], levels = adj$levels[[f]])
  }
  x
}

# Each spline variable's curved terms, from the given knots.
add_splines <- function(x, knots) {
  for (f in names(knots)) {
    curves <- rcs_terms(x[[f]], knots[[f]])
    for (j in seq_len(ncol(curves))) x[[paste0(f, "_curve", j)]] <- curves[, j]
  }
  attr(x, "knots") <- knots
  x
}

# The adjusted model for one pass: the group plus the adjustment variables.
adjusted_formula <- function(pass, stratify_by = NULL) {
  stats::update(pass_formula(pass), stats::as.formula(paste(". ~ . +", paste(model_terms(stratify_by), collapse = " + "))))
}

adjusted_cox <- function(frame, pass, robust = TRUE, stratify_by = NULL) {
  if (robust) {
    survival::coxph(adjusted_formula(pass, stratify_by), data = frame, ties = "efron", cluster = patient_id)
  } else {
    survival::coxph(adjusted_formula(pass, stratify_by), data = frame, ties = "efron")
  }
}

# Both passes' group HRs: log HR, model-based SE and robust SE for each comparison with the reference.
adjusted_estimates <- function(frame, grouping = GENOTYPE, stratify_by = NULL) {
  rows <- list()
  for (p in PASSES) {
    fit <- adjusted_cox(frame, p, stratify_by = stratify_by)
    g <- grep("^group", names(stats::coef(fit)), value = TRUE)
    contrast <- contrast_name(sub("^group", "", g), grouping)
    est <- list(log_hr = stats::coef(fit)[g], model_se = sqrt(diag(fit$naive.var))[match(g, names(stats::coef(fit)))],
                robust_se = sqrt(diag(fit$var))[match(g, names(stats::coef(fit)))])
    for (q in names(est)) {
      rows[[length(rows) + 1]] <- data.frame(pass = paste("3a", p), quantity = q, group = contrast, years = NA_real_,
                                             value = unname(est[[q]]), n_risk = NA_real_)
    }
  }
  do.call(rbind, rows)
}

# The adjusted bias in each comparison: HR(naive) / HR(corrected), both adjusted.
adjusted_bias <- function(est) {
  contrasts <- unique(est$group[est$quantity == "log_hr"])
  val <- function(p, g) est$value[est$pass == p & est$quantity == "log_hr" & est$group == g]
  data.frame(pass = "3a bias", quantity = "hr_ratio", group = contrasts, years = NA_real_,
             value = vapply(contrasts, function(g) exp(val("3a naive", g) - val("3a corrected", g)), 1),
             n_risk = NA_real_)
}

# The patient-level bootstrap for the adjusted HR ratios. The coding (the knots and the bands) is fixed from the
# patients in the model; each resample refits both models.
bootstrap_adjusted <- function(frame, resamples = BOOT_RESAMPLES, seed = BOOT_SEED, grouping = GENOTYPE,
                               stratify_by = NULL) {
  others <- levels(frame$group)[-1]
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  out <- lapply(seq_len(resamples), function(b) {
    xb <- frame[sample.int(nrow(frame), nrow(frame), replace = TRUE), ]
    lhr <- vapply(PASSES, function(p) {
      cf <- stats::coef(adjusted_cox(xb, p, robust = FALSE, stratify_by = stratify_by))
      cf[paste0("group", others)]
    }, numeric(length(others)))
    lhr <- matrix(lhr, ncol = length(PASSES), dimnames = list(NULL, PASSES))
    data.frame(resample = b, quantity = "hr_ratio", group = contrast_name(others, grouping),
               value = unname(exp(lhr[, "naive"] - lhr[, "corrected"])))
  })
  do.call(rbind, out)
}
