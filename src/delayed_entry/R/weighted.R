# Pass 3b (design plan section 8): survival curves standardized by inverse probability of treatment weighting
# (IPTW), weighted naive and weighted corrected. Ported from the CRC BPC project's R/weighted.R; the traits, the
# trimming percentiles and the balance target now come from the spec.
#
# The "treatment" is the genotype group. A multinomial propensity model, fitted with WeightIt on pass 3a's terms,
# gives each patient's chance of being in their own group. Weighting by its inverse makes every group look like
# the whole analyzed cohort in those traits (estimand "ATE": standardized to the cohort, not a causal effect). The
# weights are stabilized, then trimmed at the spec's percentiles. Balance is checked with cobalt and reported
# whatever it shows; the model is never re-specified. If any SMD stays at the target or above, 3a remains the main
# adjusted result and the 3b results are labeled imperfectly balanced.

ps_formula <- function() stats::as.formula(paste("group ~", paste(model_terms(), collapse = " + ")))

# Stabilized ATE weights, trimmed at the spec's percentiles: both ends (lower = 100 - upper), or the upper only
# (lower = 0). WeightIt takes each as one patient's own weight.
iptw <- function(frame) {
  w <- WeightIt::weightit(ps_formula(), data = frame, method = "glm", estimand = "ATE", stabilize = TRUE)
  trim <- SPEC$weighting$trim_percentiles
  if (isTRUE(all.equal(trim[1], 100 - trim[2]))) {
    suppressMessages(WeightIt::trim(w, at = trim[2] / 100, lower = TRUE))
  } else if (trim[1] == 0) {
    suppressMessages(WeightIt::trim(w, at = trim[2] / 100))
  } else {
    stop("weights can be trimmed at both ends alike, or at the upper end only")
  }
}

# Balance on every term of the propensity model, before and after weighting: for each term, the largest absolute
# SMD across the pairs of groups, and which pair that is. Binary terms are standardized too.
balance <- function(w) {
  b <- cobalt::bal.tab(w, stats = "m", binary = "std", un = TRUE)$Balance.Across.Pairs
  out <- data.frame(term = rownames(b), smd_before = abs(b$Max.Diff.Un), smd_after = abs(b$Max.Diff.Adj),
                    row.names = NULL)
  pairs <- cobalt::bal.tab(w, stats = "m", binary = "std", un = TRUE, which.treat = .all)$Pair.Balance
  largest <- function(column) {
    m <- vapply(pairs, function(p) abs(p$Balance[out$term, column]), numeric(nrow(out)))
    m <- matrix(m, nrow = nrow(out))
    pair_label(names(pairs))[max.col(m, ties.method = "first")]
  }
  out$pair_before <- largest("Diff.Un")
  out$pair_after <- largest("Diff.Adj")
  out
}

# cobalt names a pair "A vs. B"; the tables say "A vs B".
pair_label <- function(x) sub(" vs\\. ", " vs ", x)

balanced <- function(bal) all(bal$smd_after < SPEC$weighting$balance_target_smd)

# The weighted median and RMST for each group and pass, in months; no intervals (they come from the bootstrap).
weighted_curves <- function(frame, weights, pass) {
  x <- frame
  x$w <- weights
  do.call(rbind, lapply(levels(x$group), function(g) {
    d <- x[x$group == g, ]
    fit <- survival::survfit(pass_formula(pass, "1"), data = d, weights = w, conf.type = "none")
    m <- within_cap(unname(stats::quantile(fit, probs = 0.5, conf.int = FALSE)))
    data.frame(pass = paste("3b", pass), quantity = c("median", "rmst"), group = g,
               value = c(m, rmst_days(fit)) / DAYS_PER_MONTH)
  }))
}

# The weighted Cox model for one pass: the group only, weighted, with robust SEs.
weighted_cox <- function(frame, weights, pass, robust = TRUE, grouping = GENOTYPE) {
  x <- frame
  x$w <- weights
  fit <- if (robust) {
    survival::coxph(pass_formula(pass), data = x, weights = w, ties = "efron", cluster = patient_id)
  } else {
    survival::coxph(pass_formula(pass), data = x, weights = w, ties = "efron")
  }
  b <- stats::coef(fit)
  contrast <- contrast_name(sub("^group", "", names(b)), grouping)
  out <- data.frame(pass = paste("3b", pass), quantity = "log_hr", group = contrast, value = unname(b))
  if (robust) out <- rbind(out, data.frame(pass = paste("3b", pass), quantity = "robust_se", group = contrast,
                                           value = unname(sqrt(diag(fit$var)))))
  out
}

# Each comparison's weighted RMST difference (the group minus the reference), from one pass's curves.
rmst_differences <- function(curves, grouping = GENOTYPE) {
  r <- curves[curves$quantity == "rmst", ]
  reference <- r$value[r$group == grouping$levels[1]]
  others <- r[r$group != grouping$levels[1], ]
  if (length(reference) != 1 || nrow(others) == 0) return(curves[0, ])
  data.frame(pass = others$pass, quantity = "rmst_difference", group = contrast_name(others$group, grouping),
             value = others$value - reference)
}

weighted_estimates <- function(frame, weights, robust = TRUE, grouping = GENOTYPE) {
  do.call(rbind, lapply(PASSES, function(p) {
    curves <- weighted_curves(frame, weights, p)
    rbind(curves, rmst_differences(curves, grouping), weighted_cox(frame, weights, p, robust, grouping))
  }))
}

# The adjusted bias within 3b: weighted naive minus weighted corrected median and RMST for each group, and the
# change in the between-group RMST difference.
weighted_bias <- function(est, grouping = GENOTYPE) {
  val <- function(p, q, g) {
    v <- est$value[est$pass == paste("3b", p) & est$quantity == q & est$group == g]
    if (length(v) == 1) v else NA_real_
  }
  rows <- list()
  for (g in grouping$levels) {
    rows[[length(rows) + 1]] <- data.frame(quantity = c("median_diff", "rmst_diff"), group = g,
      value = c(val("naive", "median", g) - val("corrected", "median", g),
                val("naive", "rmst", g) - val("corrected", "rmst", g)))
  }
  for (g in grouping$levels[-1]) {
    gap <- function(p) val(p, "rmst", g) - val(p, "rmst", grouping$levels[1])
    rows[[length(rows) + 1]] <- data.frame(quantity = "rmst_contrast_change", group = contrast_name(g, grouping),
                                           value = gap("naive") - gap("corrected"))
  }
  out <- do.call(rbind, rows)
  out$pass <- "3b bias"
  out[, c("pass", "quantity", "group", "value")]
}

# The patient-level bootstrap for pass 3b: the weights are refitted in each resample, the coding is kept.
bootstrap_weighted <- function(frame, resamples = BOOT_RESAMPLES, seed = BOOT_SEED, grouping = GENOTYPE) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  out <- lapply(seq_len(resamples), function(b) {
    xb <- frame[sample.int(nrow(frame), nrow(frame), replace = TRUE), ]
    xb$group <- droplevels(xb$group)
    w <- iptw(xb)$weights
    est <- weighted_estimates(xb, w, robust = FALSE, grouping)
    rbind(data.frame(resample = b, est), data.frame(resample = b, weighted_bias(est, grouping)))
  })
  do.call(rbind, out)
}

# The results tables: each weighted median and RMST with its bootstrap interval, each HR with its robust one, each
# contrast with its interval, bar and label. If balance missed the target, every 3b row says so.
weighted_tables <- function(est, boot, d1, bal, described = character(), unstable = character(), grouping = GENOTYPE) {
  note <- if (balanced(bal)) "" else LABELS$imperfectly_balanced
  curves <- est[est$quantity %in% c("median", "rmst", "rmst_difference"), ]
  keyed <- boot[boot$pass %in% paste("3b", PASSES), ]
  keyed$quantity <- paste(keyed$pass, keyed$quantity)
  ci <- percentile_ci(keyed[, c("resample", "quantity", "group", "value")])
  curves$key <- paste(curves$pass, curves$quantity)
  curves <- merge(curves, ci, by.x = c("key", "group"), by.y = c("quantity", "group"), all.x = TRUE)
  curves$key <- NULL
  hr <- est[est$quantity == "log_hr", ]
  se_rows <- est[est$quantity == "robust_se", ]
  se <- se_rows$value[match(paste(hr$pass, hr$group), paste(se_rows$pass, se_rows$group))]
  hrs <- data.frame(pass = hr$pass, quantity = "hr", group = hr$group, value = exp(hr$value),
                    lower = exp(hr$value - Z * se), upper = exp(hr$value + Z * se))
  # With imputation, the interval is Rubin's (its limits come with the pooled estimates).
  limit <- function(q) {
    rows <- est[est$quantity == q, ]
    rows$value[match(paste(hr$pass, hr$group), paste(rows$pass, rows$group))]
  }
  if (any(est$quantity == "log_hr_lower")) {
    hrs$lower <- exp(limit("log_hr_lower"))
    hrs$upper <- exp(limit("log_hr_upper"))
  }
  contrasts <- boot[boot$pass == "3b bias", c("resample", "quantity", "group", "value")]
  bias <- bias_table(weighted_bias(est, grouping), contrasts, d1, grouping, described = described, unstable = unstable)
  curves$note <- note; hrs$note <- note; bias$note <- note
  list(curves = curves[, c("pass", "quantity", "group", "value", "lower", "upper", "resamples_missing", "note")],
       hrs = hrs, bias = bias)
}
