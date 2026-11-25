# Multiple imputation of the adjustment variables, for pass 3a, pass 3b and every adjusted model, when the
# complete-case rule calls for it (the CRC BPC plan's section 9, made general). Every choice comes from the spec's
# missing_data: the number of imputations, and the imputation model (its outcome terms, whether it holds the group,
# mice's method for each form of variable, and its iterations). The seed is the run's.
#
#   Only adjustment variables are imputed; the group and the times never are.
#   Hazard ratios are pooled by Rubin's rules (rubin, pool_adjusted).
#   For bootstrap intervals, each resample is imputed once, and percentile intervals are taken across resamples.
#   Other point estimates are averages across the imputations.
#   A spline's knots are set once, from the observed values, and its curved terms are computed after imputing.

imputing <- function() identical(SPEC$missing_data$rule, "multiple imputation")
imputations <- function() as.integer(SPEC$missing_data$imputations)

# The Nelson-Aalen cumulative hazard at each patient's exit, on one clock (the corrected one has delayed entry).
cumulative_hazard <- function(x, pass) {
  fit <- survival::survfit(pass_formula(pass, "1"), data = x, conf.type = "none")
  stats::stepfun(fit$time, c(0, fit$cumhaz))(x$exit_day)
}

# The data mice imputes: the adjustment variables, the group if the spec says so, and the spec's outcome terms.
imputation_data <- function(x) {
  adj <- adjustment()
  model <- SPEC$missing_data$model
  d <- x[, adj$fields, drop = FALSE]
  for (f in names(d)) if (is.factor(d[[f]])) d[[f]] <- droplevels(d[[f]])
  if (isTRUE(model$include_group)) d$mi_group <- droplevels(x$group)
  for (term in unlist(model$outcome_terms)) {
    d[[paste0("mi_", term)]] <- switch(term,
      event = x$event,
      cumulative_hazard_corrected = cumulative_hazard(x, "corrected"),
      cumulative_hazard_naive = cumulative_hazard(x, "naive"),
      entry_time = x$entry_day,
      exit_time = x$exit_day,
      stop("unknown outcome term: ", term))
  }
  d
}

# mice's method for each incomplete adjustment variable, by its form; "" for the rest (used as predictors only).
imputation_methods <- function(d) {
  adj <- adjustment()
  out <- stats::setNames(rep("", ncol(d)), names(d))
  for (i in seq_along(adj$fields)) {
    f <- adj$fields[i]
    if (anyNA(d[[f]])) out[[f]] <- SPEC$missing_data$model$methods[[adj$forms[i]]]
  }
  out
}

# `m` completed copies of the model's patients (model_data), each with its missing adjustment variables imputed, by
# the spec's engine. A seed of NA leaves the random stream as it is: inside the bootstrap, so the resamples stay
# reproducible.
impute <- function(x, m, seed = NA) {
  d <- imputation_data(x)
  method <- imputation_methods(d)
  if (all(method == "")) return(rep(list(x), m))
  model <- SPEC$missing_data$model
  done <- if (identical(model$engine, "mice")) {
    imp <- mice::mice(d, m = m, method = method, maxit = as.integer(model$iterations), seed = seed,
                      printFlag = FALSE, donors = as.integer(model$pmm_donors))
    lapply(seq_len(m), function(i) mice::complete(imp, i))
  } else {
    if (!is.na(seed)) set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    lapply(seq_len(m), function(i) impute_chain(d, method))
  }
  lapply(done, function(dd) {
    y <- x
    for (f in adjustment()$fields) {
      y[[f]] <- if (is.factor(x[[f]])) factor(as.character(dd[[f]]), levels = levels(x[[f]])) else dd[[f]]
    }
    y
  })
}

# ---- The pipeline's own engine: chained equations, as mice does them ----
# One chain: each incomplete variable starts with random draws from its observed values; then, for the spec's
# iterations, each in turn is drawn again from a model of it on every other column, with their current values.
impute_chain <- function(d, method) {
  todo <- names(method)[method != ""]
  missing <- lapply(stats::setNames(todo, todo), function(v) is.na(d[[v]]))
  for (v in todo) {
    seen <- d[[v]][!missing[[v]]]
    d[[v]][missing[[v]]] <- seen[sample.int(length(seen), sum(missing[[v]]), replace = TRUE)]
  }
  for (iteration in seq_len(as.integer(SPEC$missing_data$model$iterations))) {
    for (v in todo) d[[v]][missing[[v]]] <- draw_values(d, v, missing[[v]], method[[v]])
  }
  d
}

# The other columns that can predict: a factor needs two levels present, a number some spread.
predictors <- function(d, v) {
  p <- d[setdiff(names(d), v)]
  keep <- vapply(p, function(z) if (is.factor(z)) nlevels(droplevels(z)) > 1 else isTRUE(stats::var(z) > 0),
                 logical(1))
  p <- p[keep]
  for (f in names(p)) if (is.factor(p[[f]])) p[[f]] <- droplevels(p[[f]])
  p
}

design <- function(p) stats::model.matrix(~ ., data = p)

# A draw from the posterior of a linear regression (mice's .norm.draw): the least-squares fit with a tiny ridge for
# numerical stability, sigma drawn from its scaled chi-squared, and beta from its normal.
RIDGE <- 1e-5   # numerical: keeps X'X invertible, as mice does
norm_draw <- function(y, X) {
  xtx <- crossprod(X)
  v <- solve(xtx + diag(RIDGE * diag(xtx), nrow = ncol(X)))
  coef <- v %*% crossprod(X, y)
  df <- max(length(y) - ncol(X), 1)
  sigma <- sqrt(sum((y - X %*% coef)^2) / stats::rchisq(1, df))
  beta <- coef + t(chol((v + t(v)) / 2)) %*% stats::rnorm(ncol(X)) * sigma
  list(coef = coef, beta = beta, sigma = sigma)
}

draw_values <- function(d, v, miss, method) {
  y <- d[[v]]
  p <- predictors(d, v)
  switch(method,
         pmm = draw_pmm(y, design(p), miss),
         norm = draw_norm(y, design(p), miss),
         logreg = draw_logreg(y, design(p), miss),
         polyreg = draw_polyreg(y, p, miss),
         polr = draw_polr(y, p, miss),
         stop("unknown imputation method: ", method))
}

# Predictive mean matching (mice's type 1 matching): the observed cases' predictions from the fitted coefficients,
# the missing cases' from the drawn ones; each missing value is the observed value of one of the spec's donors,
# drawn at random among the observed cases whose predictions are closest.
draw_pmm <- function(y, X, miss) {
  obs <- !miss
  p <- norm_draw(y[obs], X[obs, , drop = FALSE])
  fitted <- X[obs, , drop = FALSE] %*% p$coef
  wanted <- X[miss, , drop = FALSE] %*% p$beta
  donors <- min(as.integer(SPEC$missing_data$model$pmm_donors), sum(obs))
  seen <- y[obs]
  vapply(as.vector(wanted), function(h) {
    nearest <- order(abs(fitted - h))[seq_len(donors)]
    seen[nearest[sample.int(donors, 1)]]
  }, numeric(1))
}

# Bayesian linear regression: the drawn coefficients' prediction plus noise from the drawn sigma.
draw_norm <- function(y, X, miss) {
  p <- norm_draw(y[!miss], X[!miss, , drop = FALSE])
  as.vector(X[miss, , drop = FALSE] %*% p$beta) + stats::rnorm(sum(miss)) * p$sigma
}

# A category drawn for each missing case from its probabilities (rows of `probs`, columns in `levels`' order).
draw_category <- function(probs, levels, all_levels) {
  probs <- matrix(probs, ncol = length(levels))
  cum <- t(apply(probs, 1, cumsum))
  if (nrow(probs) == 1) cum <- matrix(cum, nrow = 1)
  pick <- pmin(1 + rowSums(stats::runif(nrow(probs)) > cum), length(levels))
  factor(levels[pick], levels = all_levels)
}

# The guard against perfect prediction (mice's augment, after White, Daniel and Royston, Comput Stat Data Anal
# 2010): when a predictor separates a category (a group in which every observed value is one category, say), a
# logistic or multinomial fit drives that category's probability to 1 there. Pseudo-records keep the fit finite:
# for each predictor and each observed category, two records of that category with the predictor half a standard
# deviation above and below its mean (kept within its observed range) and every other predictor at its mean.
# Together they weigh as much as p + 1 patients, for p predictors. X is the design without its intercept, and the
# means, standard deviations and ranges are over all its rows, as in mice. mice leaves them out when only one value
# is missing; here they are always added.
augment_records <- function(X, lv) {
  p <- ncol(X)
  k <- length(lv)
  rows <- 2 * p * k
  if (p == 0) return(list(X = X[0, , drop = FALSE], y = factor(character(0), levels = lv), w = numeric(0)))
  at <- matrix(colMeans(X), nrow = rows, ncol = p, byrow = TRUE)
  half <- apply(X, 2, stats::sd) / 2
  for (j in seq_len(p)) {
    r <- (j - 1) * 2 * k + seq_len(2 * k)
    at[r, j] <- at[r, j] + rep(c(1, -1), k) * half[j]
  }
  at <- pmin(pmax(at, matrix(apply(X, 2, min), nrow = rows, ncol = p, byrow = TRUE)),
             matrix(apply(X, 2, max), nrow = rows, ncol = p, byrow = TRUE))
  colnames(at) <- colnames(X)
  list(X = at, y = factor(rep(rep(lv, each = 2), p), levels = lv), w = rep((p + 1) / rows, rows))
}

# Logistic regression for a two-level variable (mice's logreg), fitted with the guard's records: coefficients drawn
# from their approximate normal posterior, then each missing case drawn from its probability.
draw_logreg <- function(y, X, miss) {
  seen <- droplevels(y[!miss])
  lv <- levels(seen)
  if (length(lv) == 1) return(factor(rep(lv, sum(miss)), levels = levels(y)))
  if (length(lv) > 2) stop("logreg is for a variable with two levels; use polyreg")
  aug <- augment_records(X[, -1, drop = FALSE], lv)
  Xo <- rbind(X[!miss, , drop = FALSE], cbind(1, aug$X))
  yo <- as.integer(c(as.character(seen), as.character(aug$y)) == lv[2])
  fit <- suppressWarnings(stats::glm.fit(Xo, yo, weights = c(rep(1, length(seen)), aug$w),
                                         family = stats::binomial()))
  b <- fit$coefficients
  b[is.na(b)] <- 0
  info <- crossprod(Xo, Xo * fit$weights)
  v <- solve(info + diag(RIDGE * pmax(diag(info), 1), nrow = ncol(Xo)))
  beta <- b + t(chol((v + t(v)) / 2)) %*% stats::rnorm(length(b))
  prob <- stats::plogis(as.vector(X[miss, , drop = FALSE] %*% beta))
  draw_category(cbind(1 - prob, prob), lv, levels(y))
}

# Multinomial regression (mice's polyreg, with nnet), fitted with the guard's records: each missing case's
# probabilities of the observed categories (a row each, in the order of `seen`'s levels). X is the design without
# its intercept, over every row; with no predictors, the probabilities are the observed shares.
polyreg_probs <- function(seen, X, miss) {
  lv <- levels(seen)
  if (ncol(X) == 0) {
    return(matrix(as.vector(table(seen)) / length(seen), nrow = sum(miss), ncol = length(lv), byrow = TRUE))
  }
  colnames(X) <- paste0("x", seq_len(ncol(X)))           # plain names for the formula
  aug <- augment_records(X, lv)
  train <- data.frame(mi_y = factor(c(as.character(seen), as.character(aug$y)), levels = lv),
                      rbind(X[!miss, , drop = FALSE], aug$X))
  w <- c(rep(1, length(seen)), aug$w)
  fit <- nnet::multinom(mi_y ~ ., data = train, weights = w, trace = FALSE,
                        MaxNWts = (ncol(X) + 2) * (length(lv) + 1) * 2)
  probs <- matrix(stats::predict(fit, newdata = as.data.frame(X[miss, , drop = FALSE]), type = "probs"),
                  nrow = sum(miss))
  if (length(lv) == 2) probs <- cbind(1 - probs, probs)
  probs
}

# Each missing case drawn from its multinomial probabilities.
draw_polyreg <- function(y, p, miss) {
  seen <- droplevels(y[!miss])
  lv <- levels(seen)
  if (length(lv) == 1) return(factor(rep(lv, sum(miss)), levels = levels(y)))
  draw_category(polyreg_probs(seen, design(p)[, -1, drop = FALSE], miss), lv, levels(y))
}

# Ordered logistic regression (mice's polr, with MASS), the levels in their listed order; multinomial if it fails,
# as mice does.
draw_polr <- function(y, p, miss) {
  seen <- droplevels(y[!miss])
  lv <- levels(seen)
  if (length(lv) < 3) return(draw_polyreg(y, p, miss))
  train <- data.frame(mi_y = factor(seen, levels = lv, ordered = TRUE), p[!miss, , drop = FALSE])
  fit <- tryCatch(suppressWarnings(MASS::polr(mi_y ~ ., data = train, Hess = TRUE)), error = function(e) NULL)
  if (is.null(fit)) return(draw_polyreg(y, p, miss))
  probs <- stats::predict(fit, newdata = p[miss, , drop = FALSE], type = "probs")
  draw_category(probs, lv, levels(y))
}

# Each spline variable's knots, from its observed values among the model's patients.
observed_knots <- function(x) {
  adj <- adjustment()
  splines <- adj$fields[adj$forms == "spline"]
  lapply(stats::setNames(splines, splines), function(f) spline_knots(x[[f]][!is.na(x[[f]])], adj$knots[[f]]))
}

# A completed copy made ready for the models: unused levels dropped, the spline terms from the fixed knots.
finish_frame <- function(y, knots) {
  adj <- adjustment()
  y$group <- droplevels(y$group)
  for (i in seq_along(adj$fields)) {
    if (adj$forms[i] %in% c("categorical", "bands")) y[[adj$fields[i]]] <- droplevels(y[[adj$fields[i]]])
  }
  add_splines(y, knots)
}

# The model's patients, imputed `m` times and made ready: a list of frames.
imputed_frames <- function(x, knots, m = imputations(), seed = BOOT_SEED) {
  lapply(impute(x, m, seed), finish_frame, knots = knots)
}

# Rubin's rules for one estimate: the mean across imputations, the total variance (the mean within-imputation
# variance plus (1 + 1/m) times the between-imputation variance), and the interval at the spec's level on Rubin's
# degrees of freedom (a normal interval when the imputations agree exactly).
rubin <- function(q, u) {
  m <- length(q)
  qbar <- mean(q)
  ubar <- mean(u)
  b <- if (m > 1) stats::var(q) else 0
  total <- ubar + (1 + 1 / m) * b
  df <- if (b > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * b))^2 else Inf
  half <- stats::qt(PROBS[2], df) * sqrt(total)
  c(estimate = qbar, se = sqrt(total), lower = qbar - half, upper = qbar + half, df = df)
}

# One value from each imputation's table of estimates (NA where an imputation has none).
values_across <- function(per, pass, quantity, group) {
  vapply(per, function(e) {
    v <- e$value[e$pass == pass & e$quantity == quantity & e$group == group]
    if (length(v) == 1) v else NA_real_
  }, numeric(1))
}

# Pool the log hazard ratios of `per` (one table of estimates per imputation) by Rubin's rules: the pooled log HR,
# its model-based and robust standard errors, the interval's limits on the log scale, and Rubin's degrees of freedom.
pool_hrs <- function(per) {
  keys <- unique(per[[1]][per[[1]]$quantity == "log_hr", c("pass", "group")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    p <- keys$pass[i]
    g <- keys$group[i]
    lhr <- values_across(per, p, "log_hr", g)
    rob <- rubin(lhr, values_across(per, p, "robust_se", g)^2)
    out <- c(log_hr = rob[["estimate"]], robust_se = rob[["se"]], log_hr_lower = rob[["lower"]],
             log_hr_upper = rob[["upper"]], rubin_df = rob[["df"]])
    model <- values_across(per, p, "model_se", g)
    if (!all(is.na(model))) out <- c(out, model_se = rubin(lhr, model^2)[["se"]])
    data.frame(pass = p, quantity = names(out), group = g, years = NA_real_, value = unname(out), n_risk = NA_real_)
  }))
}

# The average across imputations of every other estimate (medians, RMST, RMST differences, contrasts). An estimate
# missing in any imputation (a median not reached, say) is missing.
average_across <- function(per) {
  first <- per[[1]]
  first$value <- vapply(seq_len(nrow(first)), function(i) {
    mean(values_across(per, first$pass[i], first$quantity[i], first$group[i]))
  }, numeric(1))
  first
}

# ---- Pass 3a, and every adjusted model, with imputation ----
# Returns the pooled estimates, the HR ratios (averaged across imputations) and the bootstrap, in which each
# resample is imputed once.
adjusted_mi <- function(cohort, covariates, grouping = GENOTYPE, stratify_by = NULL) {
  x <- model_data(cohort, covariates, grouping)
  knots <- observed_knots(x)
  per <- lapply(imputed_frames(x, knots), adjusted_estimates, grouping = grouping, stratify_by = stratify_by)
  bias <- average_across(lapply(per, adjusted_bias))
  others <- levels(droplevels(x$group))[-1]
  set.seed(BOOT_SEED, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  boot <- do.call(rbind, lapply(seq_len(BOOT_RESAMPLES), function(b) {
    xb <- x[sample.int(nrow(x), nrow(x), replace = TRUE), ]
    yb <- finish_frame(impute(xb, 1)[[1]], knots)
    lhr <- vapply(PASSES, function(p) {
      cf <- stats::coef(adjusted_cox(yb, p, robust = FALSE, stratify_by = stratify_by))
      cf[paste0("group", others)]
    }, numeric(length(others)))
    lhr <- matrix(lhr, ncol = length(PASSES), dimnames = list(NULL, PASSES))
    data.frame(resample = b, quantity = "hr_ratio", group = contrast_name(others, grouping),
               value = unname(exp(lhr[, "naive"] - lhr[, "corrected"])))
  }))
  list(estimates = pool_hrs(per), bias = bias, boot = boot, patients = nrow(x))
}

# ---- Pass 3b with imputation ----
# Each term's absolute SMD for each pair of groups, before and after weighting (terms x pairs).
pair_smds <- function(w) {
  terms <- rownames(cobalt::bal.tab(w, stats = "m", binary = "std", un = TRUE)$Balance.Across.Pairs)
  pairs <- cobalt::bal.tab(w, stats = "m", binary = "std", un = TRUE, which.treat = .all)$Pair.Balance
  get <- function(column) matrix(vapply(pairs, function(p) abs(p$Balance[terms, column]), numeric(length(terms))),
                                 nrow = length(terms), dimnames = list(terms, pair_label(names(pairs))))
  list(before = get("Diff.Un"), after = get("Diff.Adj"))
}

# The balance table from (averaged) pair SMDs: each term's largest SMD across the pairs, and which pair.
balance_from_pairs <- function(before, after) {
  largest <- function(m) colnames(m)[max.col(m, ties.method = "first")]
  data.frame(term = rownames(before), smd_before = apply(before, 1, max), smd_after = apply(after, 1, max),
             pair_before = largest(before), pair_after = largest(after), row.names = NULL)
}

# Pass 3b with imputation: in each completed copy, the weights, the balance and the weighted estimates; then the
# HRs pooled by Rubin's rules, every other estimate and the balance averaged across imputations, and the bootstrap,
# in which each resample is imputed once and its weights refitted. Also each imputation's balance before weighting
# and its completed values, which the safety gate reads to judge how precisely balance may be published.
weighted_mi <- function(cohort, covariates, grouping = GENOTYPE) {
  x <- model_data(cohort, covariates, grouping)
  knots <- observed_knots(x)
  frames <- imputed_frames(x, knots)
  fits <- lapply(frames, function(f) {
    w <- iptw(f)
    list(est = weighted_estimates(f, w$weights, grouping = grouping), pairs = pair_smds(w))
  })
  est_per <- lapply(fits, `[[`, "est")
  curves <- average_across(lapply(est_per, function(e) e[e$quantity %in% c("median", "rmst", "rmst_difference"), ]))
  hrs <- pool_hrs(lapply(est_per, function(e) data.frame(e[e$quantity %in% c("log_hr", "robust_se"), ], years = NA_real_,
                                                           n_risk = NA_real_)))
  est <- rbind(curves, hrs[, c("pass", "quantity", "group", "value")])
  mean_of <- function(part) Reduce(`+`, lapply(fits, function(f) f$pairs[[part]])) / length(fits)
  bal <- balance_from_pairs(mean_of("before"), mean_of("after"))
  by_imputation <- do.call(rbind, lapply(seq_along(fits), function(i) {
    data.frame(imputation = i, term = rownames(fits[[i]]$pairs$before), smd_before = apply(fits[[i]]$pairs$before, 1, max))
  }))
  completed <- do.call(rbind, lapply(seq_along(frames), function(i) {
    data.frame(imputation = i, frames[[i]][, c("patient_id", adjustment()$fields)])
  }))
  set.seed(BOOT_SEED, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  boot <- do.call(rbind, lapply(seq_len(BOOT_RESAMPLES), function(b) {
    xb <- x[sample.int(nrow(x), nrow(x), replace = TRUE), ]
    yb <- finish_frame(impute(xb, 1)[[1]], knots)
    w <- iptw(yb)$weights
    e <- weighted_estimates(yb, w, robust = FALSE, grouping)
    rbind(data.frame(resample = b, e), data.frame(resample = b, weighted_bias(e, grouping)))
  }))
  list(estimates = est, boot = boot, balance = bal, by_imputation = by_imputation, completed = completed,
       patients = nrow(x))
}
