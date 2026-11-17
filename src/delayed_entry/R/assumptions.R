# Checks of the assumptions behind the corrected estimates (design plan section 8). Ported from the CRC BPC
# project's R/assumptions.R; the p-value bar, where hazard ratios are split, and the estimators now come from the
# spec.
#
# Proportional hazards: the global Schoenfeld test on the corrected Cox model (cox.zph). Below the spec's bar,
#   HRs split at the spec's time are added.
# Quasi-independence: the conditional Kendall's tau between entry and survival (tranSurv::cKendall), with its
#   interval, for the whole cohort and each group, by each estimator the spec names; a group whose interval
#   excludes 0 is flagged, and estimators that disagree are reported. Each group's survival is also re-estimated
#   with the transformation-model estimator for dependent truncation (tranSurv::trSurvfit): medians only.
MEDIAN_TOL <- sqrt(.Machine$double.eps)   # the tolerance survival uses for a curve "exactly" at 0.5

# ---- Proportional hazards ----
ph_test <- function(cohort, grouping = GENOTYPE) {
  x <- analysis_frame(cohort, grouping)
  fit <- survival::coxph(pass_formula("corrected"), data = x, ties = "efron", cluster = patient_id, model = TRUE)
  z <- survival::cox.zph(fit)$table
  data.frame(term = rownames(z), chisq = z[, "chisq"], df = z[, "df"], p = z[, "p"], row.names = NULL)
}

ph_rejected <- function(ph) ph$p[ph$term == "GLOBAL"] < SPEC$ph$p_threshold

# HRs against the reference in two periods, before and after the split, from the corrected pass. Each patient's
# time at risk is split at the cut, and each group gets its own coefficient in each period. Robust SEs.
split_hrs <- function(cohort, cut = SPEC$ph$split_months * DAYS_PER_MONTH, grouping = GENOTYPE) {
  x <- analysis_frame(cohort, grouping)
  split_formula <- stats::as.formula("Surv(entry_day, exit_day, event) ~ .", env = asNamespace("survival"))
  s <- survival::survSplit(split_formula, data = x, cut = cut, start = "tstart", end = "tstop", episode = "period")
  terms <- character()
  for (g in levels(x$group)[-1]) for (k in 1:2) {
    name <- sprintf("g%d_p%d", match(g, levels(x$group)), k)
    s[[name]] <- as.numeric(s$group == g & s$period == k)
    terms <- c(terms, name)
  }
  formula <- stats::as.formula(paste("survival::Surv(tstart, tstop, event) ~", paste(terms, collapse = " + ")))
  fit <- survival::coxph(formula, data = s, ties = "efron", cluster = patient_id)
  b <- stats::coef(fit)
  se <- sqrt(diag(fit$var))
  g <- levels(x$group)[as.integer(sub("^g(\\d+)_p\\d$", "\\1", names(b)))]
  months <- format(SPEC$ph$split_months)
  data.frame(contrast = contrast_name(g, grouping),
             period = ifelse(grepl("_p1$", names(b)), paste("first", months, "months"), paste("after", months, "months")),
             hr = unname(exp(b)), lower = unname(exp(b - Z * se)), upper = unname(exp(b + Z * se)))
}

# ---- Quasi-independence ----
quasi_independence <- function(cohort, grouping = GENOTYPE) {
  x <- analysis_frame(cohort, grouping)
  groups <- c(list(x), lapply(levels(x$group), function(g) x[x$group == g, ]))
  names(groups) <- c(ALL_PATIENTS, levels(x$group))
  rows <- list()
  for (est in SPEC$tau_estimators) for (g in names(groups)) {
    d <- groups[[g]]
    k <- tranSurv::cKendall(trun = d$entry_day, obs = d$exit_day, delta = d$event, method = est)
    rows[[length(rows) + 1]] <- data.frame(estimator = est, group = g, tau = k$PE, se = k$SE,
                                           lower = k$PE - Z * k$SE, upper = k$PE + Z * k$SE, p = k$p.value,
                                           flagged = (k$PE - Z * k$SE > 0) | (k$PE + Z * k$SE < 0))
  }
  do.call(rbind, rows)
}

# The median read off any step curve: the earliest time at or below 0.5, or the midpoint of a stretch at exactly
# 0.5 (within survival's tolerance).
median_from_curve <- function(times, surv, tol = MEDIAN_TOL) {
  i <- which(surv <= 0.5 + tol)[1]
  if (is.na(i)) return(NA_real_)
  if (abs(surv[i] - 0.5) <= tol) {
    later <- which(surv < 0.5 - tol & seq_along(surv) > i)[1]
    return(if (is.na(later)) times[i] else (times[i] + times[later]) / 2)
  }
  times[i]
}

# Each group's median from the transformation-model estimator, beside the corrected KM median (in months), and
# the fitted transformation parameter a. tranSurv's warnings are kept.
transformation_medians <- function(cohort, grouping = GENOTYPE) {
  x <- analysis_frame(cohort, grouping)
  do.call(rbind, lapply(levels(x$group), function(g) {
    d <- x[x$group == g, ]
    r <- with_warnings(tranSurv::trSurvfit(trun = d$entry_day, obs = d$exit_day, delta = d$event))
    s <- r$value$surv
    data.frame(group = g, a = r$value$byTau$par[1],
               median_transformation = within_cap(median_from_curve(s$Time, s$trSurv)) / DAYS_PER_MONTH,
               median_km = within_cap(median_from_curve(s$Time, s$kmSurv)) / DAYS_PER_MONTH,
               warning = paste(unique(r$warnings), collapse = "; "))
  }))
}
