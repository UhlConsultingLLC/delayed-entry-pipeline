# Passes 1 and 2 (design plan section 8), with R's survival package. Ported from the CRC BPC project's
# R/passes.R; its fixed values (the groups, the years, the horizon, the median limit) now come from the spec.
#
# Pass 1, naive: Kaplan-Meier (KM) curves and Cox models from diagnosis, as if every patient had been followed
# from day 0, Surv(exit_day, event). Pass 2, corrected: the same patients in counting-process form,
# Surv(entry_day, exit_day, event): a patient is in the risk set only after entry, (entry, exit].
#
# The estimates come back in one long table (pass, quantity, group, years, value, n_risk). Each pass gives KM
# survival at the fixed times, the median and its interval, the RMST to the horizon, and the Cox model. Medians
# and RMST are in months.

# A comparison's label: "Group vs Reference".
contrast_name <- function(level, grouping = GENOTYPE) paste(level, "vs", grouping$reference)

# The survival outcome for each pass, with `rhs` on the right of the ~.
pass_formula <- function(pass, rhs = "group") {
  lhs <- switch(pass,
                naive = "Surv(exit_day, event)",                # at risk from diagnosis (day 0)
                corrected = "Surv(entry_day, exit_day, event)", # at risk only after entry
                stop("unknown pass: ", pass))
  stats::as.formula(paste(lhs, "~", rhs), env = asNamespace("survival"))
}

# The patients the passes analyze: those with time at risk. Their group goes in the column `group`, with the
# reference first. A group with no patients is dropped. A patient in none of the groups stops the analysis.
analysis_frame <- function(cohort, grouping = GENOTYPE) {
  x <- cohort[cohort$in_analysis, ]
  x$group <- droplevels(factor(x[[grouping$column]], levels = grouping$levels))
  if (anyNA(x$group)) stop(sum(is.na(x$group)), " patients aren't in any of the groups in ", grouping$column)
  x
}

# The number at risk at each of `days`, counted directly: from day 0 in the naive pass, and in (entry, exit]
# in the corrected one, so the groups add up to all patients.
at_risk <- function(d, pass, days) {
  vapply(days, function(t) {
    if (pass == "naive") sum(d$exit_day >= t) else sum(d$entry_day < t & d$exit_day >= t)
  }, numeric(1))
}

# KM survival at each of `years`, for all patients and for each group, with the number at risk.
km_survival <- function(cohort, pass, years = YEARS, grouping = GENOTYPE) {
  x <- analysis_frame(cohort, grouping)
  one <- function(d, group) {
    fit <- survival::survfit(pass_formula(pass, "1"), data = d)
    s <- summary(fit, times = years * DAYS_PER_YEAR, extend = TRUE)
    data.frame(pass = pass, quantity = "survival", group = group, years = years,
               value = s$surv, n_risk = at_risk(d, pass, years * DAYS_PER_YEAR))
  }
  groups <- lapply(levels(x$group), function(g) one(x[x$group == g, ], g))
  do.call(rbind, c(list(one(x, ALL_PATIENTS)), groups))
}

# The area under a KM curve from day 0 to `tau` days: the restricted mean survival time.
rmst_days <- function(fit, tau = RMST_DAYS) {
  keep <- fit$time < tau
  edges <- c(0, fit$time[keep], tau)
  heights <- c(1, fit$surv[keep])
  sum(diff(edges) * heights)
}

# A time beyond the median limit, or never reached, is "not reached" (NA).
within_cap <- function(days) ifelse(is.na(days) | days > MEDIAN_CAP_DAYS, NA_real_, days)

# The median, its interval, and the RMST, all in months, for all patients and each group. The median is
# survival's own: the earliest time the curve is at or below 0.5, or the midpoint of a stretch where it equals
# 0.5 exactly. Its interval uses log-log limits for the curve, at the spec's level.
km_summaries <- function(cohort, pass, grouping = GENOTYPE) {
  x <- analysis_frame(cohort, grouping)
  one <- function(d, group) {
    fit <- survival::survfit(pass_formula(pass, "1"), data = d, conf.type = "log-log", conf.int = SPEC$interval_level)
    q <- stats::quantile(fit, probs = 0.5, conf.int = TRUE)
    v <- c(median = within_cap(unname(q$quantile)), median_lower = within_cap(unname(q$lower)),
           median_upper = within_cap(unname(q$upper)), rmst = rmst_days(fit)) / DAYS_PER_MONTH
    data.frame(pass = pass, quantity = names(v), group = group, years = NA_real_, value = unname(v), n_risk = NA_real_)
  }
  groups <- lapply(levels(x$group), function(g) one(x[x$group == g, ], g))
  do.call(rbind, c(list(one(x, ALL_PATIENTS)), groups))
}

# Every step of each KM curve, for the figure: the day, the survival just after it, and how many were at risk
# and died then. Row-level in effect: the safety gate decides what may be shown.
km_curves <- function(cohort, pass, grouping = GENOTYPE) {
  x <- analysis_frame(cohort, grouping)
  one <- function(d, group) {
    fit <- survival::survfit(pass_formula(pass, "1"), data = d, conf.type = "none")
    data.frame(pass = pass, group = group, day = fit$time, surv = fit$surv, n_risk = fit$n.risk, n_event = fit$n.event)
  }
  groups <- lapply(levels(x$group), function(g) one(x[x$group == g, ], g))
  do.call(rbind, c(list(one(x, ALL_PATIENTS)), groups))
}

# The Cox model: each group against the reference, Efron's method for ties, robust (sandwich) standard errors
# (clustered on the patient), and the model-based ones kept too.
cox_hr <- function(cohort, pass, grouping = GENOTYPE) {
  x <- analysis_frame(cohort, grouping)
  fit <- survival::coxph(pass_formula(pass), data = x, ties = "efron", cluster = patient_id)
  b <- stats::coef(fit)
  est <- list(log_hr = b,
              model_se = sqrt(diag(fit$naive.var)),
              robust_se = sqrt(diag(fit$var)))
  contrast <- contrast_name(sub("^group", "", names(b)), grouping)
  do.call(rbind, lapply(names(est), function(q)
    data.frame(pass = pass, quantity = q, group = contrast, years = NA_real_, value = unname(est[[q]]), n_risk = NA_real_)))
}

# Both passes, every estimate, in one long table.
pass_estimates <- function(cohort, grouping = GENOTYPE) {
  do.call(rbind, lapply(PASSES, function(p) rbind(km_survival(cohort, p, grouping = grouping),
                                                  km_summaries(cohort, p, grouping), cox_hr(cohort, p, grouping))))
}
