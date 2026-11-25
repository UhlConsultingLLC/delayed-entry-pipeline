# The main estimates (design plan section 8): how far the naive pass is from the corrected one, with percentile
# intervals from a patient-level bootstrap, and a label for each. Ported from the CRC BPC project's R/bias.R; the
# bars, the resamples, the seed, the interval level and the labels' words now come from the spec.
#
#   In a group's survival:  naive minus corrected median and RMST.
#   In a comparison:        HR(naive) / HR(corrected), and the change in the between-group difference in RMST.
#
# Labels: "clearly material" if the estimate and its whole interval are at or beyond the bar on the same side;
# "material" if the estimate is, but the interval crosses the bar; "not material" if the estimate falls short.

# The naive-minus-corrected contrasts, from one table of pass estimates. One row per contrast.
bias_estimates <- function(est, grouping = GENOTYPE) {
  val <- function(p, q, g) {
    v <- est$value[est$pass == p & est$quantity == q & est$group == g]
    if (length(v) == 1) v else NA_real_
  }
  rows <- list()
  add <- function(q, g, v) rows[[length(rows) + 1]] <<- data.frame(pass = "bias", quantity = q, group = g,
                                                                    years = NA_real_, value = v, n_risk = NA_real_)
  reference <- grouping$levels[1]
  for (g in c(ALL_PATIENTS, grouping$levels)) {
    add("median_diff", g, val("naive", "median", g) - val("corrected", "median", g))
    add("rmst_diff", g, val("naive", "rmst", g) - val("corrected", "rmst", g))
  }
  for (g in grouping$levels[-1]) {
    contrast <- contrast_name(g, grouping)
    add("hr_ratio", contrast, exp(val("naive", "log_hr", contrast) - val("corrected", "log_hr", contrast)))
    gap <- function(p) val(p, "rmst", g) - val(p, "rmst", reference)
    add("rmst_contrast_change", contrast, gap("naive") - gap("corrected"))
  }
  do.call(rbind, rows)
}

# Each group's RMST and each comparison's RMST difference (the group minus the reference), for each pass,
# keyed "naive rmst", "corrected rmst_difference" and so on.
rmst_rows <- function(est, grouping = GENOTYPE) {
  val <- function(p, g) {
    v <- est$value[est$pass == p & est$quantity == "rmst" & est$group == g]
    if (length(v) == 1) v else NA_real_
  }
  reference <- grouping$levels[1]
  rows <- list()
  for (p in PASSES) {
    for (g in c(ALL_PATIENTS, grouping$levels)) {
      rows[[length(rows) + 1]] <- data.frame(quantity = paste(p, "rmst"), group = g, value = val(p, g))
    }
    for (g in grouping$levels[-1]) {
      rows[[length(rows) + 1]] <- data.frame(quantity = paste(p, "rmst_difference"), group = contrast_name(g, grouping),
                                             value = val(p, g) - val(p, reference))
    }
  }
  do.call(rbind, rows)
}

# The RMSTs and RMST differences of passes 1 and 2 with their bootstrap intervals, from the same resamples as the
# contrasts. In months.
rmst_table <- function(est, boot, grouping = GENOTYPE) {
  t <- merge(rmst_rows(est, grouping), percentile_ci(boot), by = c("quantity", "group"), all.x = TRUE)
  t$pass <- sub(" .*$", "", t$quantity)
  t$quantity <- sub("^\\S+ ", "", t$quantity)
  t <- t[order(match(t$pass, PASSES), t$quantity == "rmst_difference",
               match(t$group, c(ALL_PATIENTS, grouping$levels, contrast_name(grouping$levels[-1], grouping)))), ]
  rownames(t) <- NULL
  t[, c("pass", "quantity", "group", "value", "lower", "upper", "resamples_missing")]
}

# Point estimates only, as fast as possible, for one bootstrap resample: medians and RMST without intervals,
# and Cox models without robust SEs (the coefficients are the same either way).
quick_estimates <- function(x, grouping = GENOTYPE) {
  rows <- list()
  for (p in PASSES) {
    groups <- c(list(x), lapply(levels(x$group), function(g) x[x$group == g, ]))
    names(groups) <- c(ALL_PATIENTS, levels(x$group))
    for (g in names(groups)) {
      fit <- survival::survfit(pass_formula(p, "1"), data = groups[[g]], conf.type = "none")
      m <- within_cap(unname(stats::quantile(fit, probs = 0.5, conf.int = FALSE))) / DAYS_PER_MONTH
      rows[[length(rows) + 1]] <- data.frame(pass = p, quantity = c("median", "rmst"), group = g,
                                             value = c(m, rmst_days(fit) / DAYS_PER_MONTH))
    }
    b <- stats::coef(survival::coxph(pass_formula(p), data = x, ties = "efron"))
    rows[[length(rows) + 1]] <- data.frame(pass = p, quantity = "log_hr",
                                           group = contrast_name(sub("^group", "", names(b)), grouping), value = unname(b))
  }
  do.call(rbind, rows)
}

# The patient-level bootstrap: draw as many patients as there are, with replacement, and recompute every
# contrast, every RMST and every RMST difference. One row per resample and estimate.
bootstrap_bias <- function(cohort, resamples = BOOT_RESAMPLES, seed = BOOT_SEED, grouping = GENOTYPE) {
  x <- analysis_frame(cohort, grouping)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  out <- lapply(seq_len(resamples), function(b) {
    xb <- x[sample.int(nrow(x), nrow(x), replace = TRUE), ]
    xb$group <- droplevels(xb$group)
    quick <- quick_estimates(xb, grouping)
    est <- rbind(bias_estimates(quick, grouping)[, c("quantity", "group", "value")], rmst_rows(quick, grouping))
    data.frame(resample = b, quantity = est$quantity, group = est$group, value = est$value)
  })
  do.call(rbind, out)
}

# Percentile intervals at the spec's level. If any resample couldn't give a value (a median not reached, say),
# the interval isn't estimated, and the number of such resamples is reported.
percentile_ci <- function(boot) {
  keys <- unique(boot[, c("quantity", "group")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    v <- boot$value[boot$quantity == keys$quantity[i] & boot$group == keys$group[i]]
    missing <- sum(is.na(v))
    lim <- if (missing == 0) stats::quantile(v, PROBS, names = FALSE) else c(NA_real_, NA_real_)
    data.frame(quantity = keys$quantity[i], group = keys$group[i], lower = lim[1], upper = lim[2],
               resamples = length(v), resamples_missing = missing)
  }))
}

# The label for one estimate and its interval. `center` is 1 for a ratio, 0 for a difference.
label_bias <- function(estimate, lower, upper, bar, center = 0) {
  d <- estimate - center
  if (is.na(d)) return(NA_character_)
  if (abs(d) < bar - BAR_TOLERANCE) return(LABELS$not_material)
  beyond <- function(v) !is.na(v) && if (d > 0) v - center >= bar - BAR_TOLERANCE else v - center <= -bar + BAR_TOLERANCE
  if (beyond(lower) && beyond(upper)) LABELS$clearly_material else LABELS$material
}

# Every contrast with its interval, bar and label. `d1` holds the smallest-group rule's decision for each
# comparison (the Python side decides it from counts): a comparison not compared gets no label, and neither does
# a group or comparison in `described` (named: why it is described) or a group in `unstable` (its corrected risk
# set never reaches the stability size). A group whose median isn't reached is summarized by RMST instead.
# "All patients" is described, not labeled. The comparison of a group the settings describe only
# (SPEC$settings_described, by comparison name) is described for that reason in every analysis, not "not compared
# (smallest group)"; the group's own estimates keep whatever `described` says.
bias_table <- function(bias, boot, d1, grouping = GENOTYPE, described = character(), unstable = character()) {
  by_settings <- named_text(SPEC$settings_described)
  described <- c(described, by_settings[setdiff(names(by_settings), names(described))])
  t <- merge(bias[, c("quantity", "group", "value")], percentile_ci(boot), by = c("quantity", "group"), all.x = TRUE)
  names(t)[names(t) == "value"] <- "estimate"
  t$bar <- unname(BARS[t$quantity])
  t$label <- mapply(label_bias, t$estimate, t$lower, t$upper, t$bar,
                    ifelse(t$quantity == "hr_ratio", 1, 0))
  t$label[t$group %in% d1$comparison[d1$decision != "compare"]] <- LABELS$not_compared
  t$label[t$quantity == "median_diff" & is.na(t$estimate) & t$group != ALL_PATIENTS] <- LABELS$median_not_reached
  t$label[t$group %in% unstable] <- LABELS$unstable
  for (g in names(described)) t$label[t$group == g] <- sprintf(LABELS$described, described[[g]])
  t$label[t$group == ALL_PATIENTS] <- LABELS$all_patients
  order <- c("median_diff", "rmst_diff", "hr_ratio", "rmst_contrast_change")
  t[order(match(t$quantity, order), match(t$group, c(ALL_PATIENTS, grouping$levels,
                                                      contrast_name(grouping$levels[-1], grouping)))), ]
}

# The groups described in the primary analysis: those whose comparison the smallest-group rule didn't fit, or the
# settings describe only (named: why).
described_from <- function(x) named_text(x)
