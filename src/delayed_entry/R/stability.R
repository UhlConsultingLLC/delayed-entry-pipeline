# The stability rule for left-truncated curves (design plan section 8). Ported from the CRC BPC project's
# R/stability.R; the risk set's size now comes from the spec.
#
#   t*(g): the earliest time at which group g's corrected risk set reaches the spec's size.
#   If any death in group g comes at or before t*(g), the conditional contrast becomes co-primary for that group:
#   naive and corrected survival, both conditional on surviving to t*(g).
#
# A curve conditional on surviving to t is the KM curve of the patients still followed after t, with nobody at
# risk before t: in the naive pass everyone starts at t; in the corrected pass each starts at the later of their
# entry and t. Its median and RMST, measured from diagnosis like every other summary, come from R/passes.R.

# The corrected risk set at time t: patients who have entered (entry < t) and are still followed (exit >= t).
risk_set <- function(d, t) sum(d$entry_day < t & d$exit_day >= t)

# t*: the earliest entry day after which the corrected risk set reaches `size`. NA if it never does.
t_star <- function(d, size = STABLE_RISK_SET) {
  for (e in sort(unique(d$entry_day))) {
    if (sum(d$entry_day <= e & d$exit_day > e) >= size) return(e)
  }
  NA_real_
}

# For each group: its t* (in days), the risk set at its first death, how many deaths come at or before t*, and
# whether the conditional contrast is co-primary. Exact counts: for checking.
stability <- function(cohort, grouping = GENOTYPE) {
  x <- analysis_frame(cohort, grouping)
  do.call(rbind, lapply(levels(x$group), function(g) {
    d <- x[x$group == g, ]
    ts <- t_star(d)
    deaths <- sort(d$exit_day[d$event == 1])
    early <- if (is.na(ts)) length(deaths) else sum(deaths <= ts)
    data.frame(group = g, patients = nrow(d), t_star_days = ts,
               risk_set_at_first_death = if (length(deaths)) risk_set(d, deaths[1]) else NA_integer_,
               deaths_by_t_star = early, co_primary = !is.na(ts) & early > 0)
  }))
}

# The groups whose corrected risk set never reaches the size: described, not labeled.
unstable_groups <- function(stab) stab$group[is.na(stab$t_star_days)]

# The cohort for curves conditional on surviving to `ts`.
condition_on <- function(d, ts, pass) {
  x <- d[d$in_analysis & d$exit_day > ts, ]
  x$entry_day <- if (pass == "naive") rep(ts, nrow(x)) else pmax(x$entry_day, ts)
  x
}

# The conditional median and RMST, in months, for each group with a t*, both passes.
conditional_estimates <- function(cohort, stab, grouping = GENOTYPE) {
  x <- analysis_frame(cohort, grouping)
  rows <- list()
  for (i in which(!is.na(stab$t_star_days))) {
    g <- stab$group[i]
    for (p in PASSES) {
      d <- condition_on(x[x$group == g, ], stab$t_star_days[i], p)
      fit <- survival::survfit(survival::Surv(entry_day, exit_day, event) ~ 1, data = d, conf.type = "none")
      m <- within_cap(unname(stats::quantile(fit, probs = 0.5, conf.int = FALSE)))
      rows[[length(rows) + 1]] <- data.frame(pass = paste("conditional", p), quantity = c("median", "rmst"),
                                             group = g, years = NA_real_,
                                             value = c(m, rmst_days(fit)) / DAYS_PER_MONTH, n_risk = NA_real_)
    }
  }
  if (!length(rows)) return(data.frame(pass = character(), quantity = character(), group = character(),
                                       years = numeric(), value = numeric(), n_risk = numeric()))
  do.call(rbind, rows)
}

# The conditional contrasts: naive minus corrected conditional median and RMST, for each group.
conditional_bias <- function(est) {
  val <- function(p, q, g) {
    v <- est$value[est$pass == paste("conditional", p) & est$quantity == q & est$group == g]
    if (length(v) == 1) v else NA_real_
  }
  groups <- unique(est$group)
  if (!length(groups)) return(data.frame(pass = character(), quantity = character(), group = character(),
                                         years = numeric(), value = numeric(), n_risk = numeric()))
  do.call(rbind, lapply(groups, function(g) data.frame(
    pass = "conditional bias", quantity = c("median_diff", "rmst_diff"), group = g, years = NA_real_,
    value = c(val("naive", "median", g) - val("corrected", "median", g),
              val("naive", "rmst", g) - val("corrected", "rmst", g)), n_risk = NA_real_)))
}

# t* as rows of the long estimates table.
stability_rows <- function(stab) {
  rbind(data.frame(pass = "stability", quantity = "t_star_days", group = stab$group, years = NA_real_,
                   value = stab$t_star_days, n_risk = NA_real_),
        data.frame(pass = "stability", quantity = "deaths_by_t_star", group = stab$group, years = NA_real_,
                   value = stab$deaths_by_t_star, n_risk = NA_real_))
}

# The bootstrap for the conditional contrasts. t* and which groups are co-primary are set once, from the analysis
# cohort, and kept in every resample.
bootstrap_conditional <- function(cohort, stab, resamples = BOOT_RESAMPLES, seed = BOOT_SEED, grouping = GENOTYPE) {
  x <- analysis_frame(cohort, grouping)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  out <- lapply(seq_len(resamples), function(b) {
    xb <- x[sample.int(nrow(x), nrow(x), replace = TRUE), ]
    est <- conditional_bias(conditional_estimates(xb, stab, grouping))
    data.frame(resample = rep(b, nrow(est)), quantity = est$quantity, group = est$group, value = est$value)
  })
  do.call(rbind, out)
}

# The conditional contrasts with intervals and labels. A group whose conditional contrast isn't co-primary is
# described only; a co-primary one is labeled against the same bars. A described group stays described.
conditional_table <- function(bias, boot, stab, described = character(), grouping = GENOTYPE) {
  t <- bias_table(bias, boot, data.frame(comparison = character(), decision = character()), grouping,
                  described = described)
  no_death <- t$group %in% stab$group[!stab$co_primary] & !(t$group %in% names(described))
  t$label[no_death] <- LABELS$no_death_before_t_star
  t$co_primary <- t$group %in% stab$group[stab$co_primary] & !(t$group %in% names(described))
  t
}
