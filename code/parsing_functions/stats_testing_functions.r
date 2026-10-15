# Shared statistical calculations for co-accessibility analyses.
# Source this file explicitly before testing or calling plot_coaccess_pair().
# data.table is required for donor regression and grouped summaries; logistf
# is required only for separated donor fits. No analysis runs when sourced.

# BH families are caller-defined. Failed fits may count in n without being
# assigned a p-value; outputs retain NA in the original positions.
adjust_bh <- function(p, n = length(p)) {
  valid <- is.finite(p)
  if (length(n) != 1L || !is.finite(n) || n < sum(valid) || n != floor(n))
    stop("n must be an integer at least as large as the number of valid tests.")
  if (any(p[valid] < 0 | p[valid] > 1)) stop("p-values must be in [0, 1].")
  out <- rep(NA_real_, length(p))
  if (any(valid)) out[valid] <- stats::p.adjust(p[valid], method = "BH", n = n)
  out
}

# Expected joint accessibility for each supplied 2x2 table; no pseudocount.
# Sum these expectations AFTER calculating them when conditioning on donors.
expected_coaccess <- function(closed, first, second, both) {
  n <- closed + first + second + both
  ifelse(n > 0, (both + first) * (both + second) / n, NA_real_)
}

# Monte Carlo correction, including ties. Smaller model p-values use "less";
# positive correlation tests use "greater". Invalid null draws are excluded.
empirical_pvalue <- function(observed, null, alternative = c("greater", "less")) {
  alternative <- match.arg(alternative)
  if (length(observed) != 1L || !is.finite(observed)) return(NA_real_)
  null <- null[is.finite(null)]
  if (!length(null)) return(NA_real_)
  extreme <- if (alternative == "greater") null >= observed else null <= observed
  (1 + sum(extreme)) / (1 + length(null))
}

# Special-case exact one-sided rank-separation probability, not a general
# two-group test. Returns NA unless every value in higher exceeds every lower.
complete_separation_pvalue <- function(higher, lower) {
  if (!length(higher) || !length(lower) || any(!is.finite(c(higher, lower))))
    return(NA_real_)
  if (min(higher) > max(lower))
    1 / choose(length(higher) + length(lower), length(higher)) else NA_real_
}

# Vectorized two-sided Fisher p-values for supplied counts (no implicit pseudocount).
fisher_p_two_sided <- function(a, b, c, d) {
  m <- a + b; n <- c + d; k <- a + c
  out <- numeric(length(a))
  for (idx in split(seq_along(a), paste(m, n, k, sep = "/"))) {
    i <- idx[1L]
    support <- max(0L, k[i] - n[i]):min(k[i], m[i])
    dens <- dhyper(support, m[i], n[i], k[i])
    
    # Sorts the probabilities from smallest to largest and takes a running sum. A table’s two-sided p-value is the total probability of all tables no more likely than it
    o <- order(dens)
    ds <- dens[o]
    grp <- cumsum(c(TRUE, diff(ds) > ds[-length(ds)] * 1e-7))
    pv <- numeric(length(ds))
    pv[o] <- as.numeric(tapply(cumsum(ds), grp, max))[grp]
    out[idx] <- pmin(1, pv[match(a[idx], support)])
  }
  out
}

# Canonical table interface. Columns: neither, CRE1-only, CRE2-only, both.
# Conditional MLE is distinct from a pseudocount cross-product odds ratio.
fisher_test_counts <- function(counts, pseudocount = 1L, estimate = TRUE,
                               include_raw = TRUE) {
  counts <- as.matrix(counts)
  if (ncol(counts) != 4L || any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts))) stop("counts must have four nonnegative integer columns.")
  if (length(pseudocount) != 1L || !is.finite(pseudocount) ||
      pseudocount < 0 || pseudocount != floor(pseudocount))
    stop("pseudocount must be a nonnegative integer.")
  x <- counts + pseudocount
  out <- data.frame(pval = fisher_p_two_sided(x[, 1], x[, 3], x[, 2], x[, 4]))
  if (estimate) out$fisher_estimate <- vapply(seq_len(nrow(x)), function(i) {
    tab <- matrix(x[i, c(1, 3, 2, 4)], nrow = 2, byrow = TRUE)
    unname(stats::fisher.test(tab)$estimate)
  }, numeric(1))
  if (include_raw) out$pval_raw <- fisher_p_two_sided(
    counts[, 1], counts[, 3], counts[, 2], counts[, 4])
  out
}

# Adapter from label_reads() output to the canonical count-table test.
coaccess_2x2 <- function(labels) {
  levels <- c("neither", "CRE1 only", "CRE2 only", "both accessible")
  n <- table(factor(labels$config[!is.na(labels$config)], levels = levels))
  fit <- fisher_test_counts(matrix(as.numeric(n), nrow = 1), include_raw = FALSE)
  data.frame(both = n[[4]], cre1_only = n[[2]], cre2_only = n[[3]],
    neither = n[[1]], n_shared = sum(n), fisher_or = fit$fisher_estimate,
    fisher_p = fit$pval)
}

# Pooled and sample-centered rank correlations under within-sample shuffling.
permute_pooled_spearman <- function(fibers, n_permutations = 10000L,
                                    seed = 20261015L) {
  fibers <- as.data.frame(fibers)
  required <- c("sample_name", "prop1", "prop2")
  if (!all(required %in% names(fibers))) {
    stop("fibers must contain sample_name, prop1, and prop2.")
  }
  if (!is.numeric(fibers$prop1) || !is.numeric(fibers$prop2)) {
    stop("prop1 and prop2 must be numeric.")
  }
  if (length(n_permutations) != 1L || !is.finite(n_permutations) ||
      n_permutations < 1 || n_permutations != floor(n_permutations) ||
      n_permutations > .Machine$integer.max) {
    stop("n_permutations must be a positive integer.")
  }
  if (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed != floor(seed) || seed > .Machine$integer.max) {
    stop("seed must be an integer between zero and .Machine$integer.max.")
  }
  n_permutations <- as.integer(n_permutations)
  seed <- as.integer(seed)
  finite_pair <- is.finite(fibers$prop1) & is.finite(fibers$prop2)
  d <- as.data.frame(fibers[finite_pair, required, drop = FALSE])
  d$sample_name <- as.character(d$sample_name)
  if (anyNA(d$sample_name) || any(!nzchar(d$sample_name))) {
    stop("Every finite fiber pair must retain a nonmissing sample_name.")
  }
  groups <- split(seq_len(nrow(d)), d$sample_name, drop = TRUE)
  varies1 <- vapply(groups, function(i) length(unique(d$prop1[i])) > 1L,
                   logical(1))
  varies2 <- vapply(groups, function(i) length(unique(d$prop2[i])) > 1L,
                   logical(1))
  sample_sizes <- lengths(groups)
  stats <- data.frame(
    n_input_fibers = nrow(fibers),
    n_paired_finite_fibers = nrow(d),
    n_excluded_nonfinite_pairs = sum(!finite_pair),
    n_samples = length(groups),
    n_singleton_samples = sum(sample_sizes == 1L),
    n_samples_varying_prop1 = sum(varies1),
    n_samples_varying_prop2 = sum(varies2),
    n_samples_varying_both = sum(varies1 & varies2),
    n_permutations_requested = n_permutations,
    n_permutations_valid = 0L,
    n_adjusted_permutations_valid = 0L,
    seed = seed,
    spearman_rho = NA_real_,
    p_positive = NA_real_,
    within_sample_adjusted_rho = NA_real_,
    adjusted_p_positive = NA_real_,
    null_mean = NA_real_,
    null_sd = NA_real_,
    adjusted_null_mean = NA_real_,
    adjusted_null_sd = NA_real_,
    status = "ok",
    adjusted_status = "ok",
    stringsAsFactors = FALSE
  )
  null <- data.frame(permutation = seq_len(n_permutations),
                     spearman_rho = rep(NA_real_, n_permutations),
                     within_sample_adjusted_rho = rep(NA_real_, n_permutations))
  if (nrow(d) < 2L || length(unique(d$prop1)) < 2L ||
      length(unique(d$prop2)) < 2L) {
    stats$status <- if (nrow(d) < 2L) "fewer_than_two_finite_fibers" else
      if (length(unique(d$prop1)) < 2L) "constant_CRE1_proportion" else
        "constant_CRE2_proportion"
    stats$adjusted_status <- stats$status
    return(list(statistics = stats, null = null))
  }

  r1 <- rank(d$prop1, ties.method = "average")
  r2 <- rank(d$prop2, ties.method = "average")
  c1 <- r1 - mean(r1)
  c2 <- r2 - mean(r2)
  a1 <- r1
  a2 <- r2
  for (i in groups) {
    a1[i] <- r1[i] - mean(r1[i])
    a2[i] <- r2[i] - mean(r2[i])
  }
  pooled_denominator <- sqrt(sum(c1^2) * sum(c2^2))
  adjusted_denominator <- sqrt(sum(a1^2) * sum(a2^2))
  # Both observed and permuted statistics use precisely this dot-product form.
  # Clipping protects only against machine rounding outside [-1, 1].
  rank_correlation <- function(x, y, denominator) {
    pmax(-1, pmin(1, sum(x * y) / denominator))
  }
  stats$spearman_rho <- rank_correlation(c1, c2, pooled_denominator)
  adjusted_defined <- is.finite(adjusted_denominator) && adjusted_denominator > 0
  if (adjusted_defined) {
    stats$within_sample_adjusted_rho <- rank_correlation(a1, a2,
                                                       adjusted_denominator)
  } else {
    stats$adjusted_status <- "no_within_sample_variation_in_one_or_both_CREs"
  }
  shuffle_groups <- groups[varies2]
  if (!length(shuffle_groups)) {
    stats$status <- "no_within_sample_CRE2_variation; null_is_fixed"
  }

  # Leave the caller's random-number state unchanged while making each pair
  # reproducible from its exported seed.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  for (b in seq_len(n_permutations)) {
    permutation_index <- seq_len(nrow(d))
    for (i in shuffle_groups) {
      # sample.int avoids sample(x)'s special behavior for length-one integers.
      permutation_index[i] <- i[sample.int(length(i))]
    }
    null$spearman_rho[b] <- rank_correlation(c1, c2[permutation_index],
                                            pooled_denominator)
    if (adjusted_defined) {
      # Sample means of ranks are invariant to these within-sample shuffles.
      null$within_sample_adjusted_rho[b] <- rank_correlation(
        a1, a2[permutation_index], adjusted_denominator)
    }
  }
  stats$n_permutations_valid <- sum(is.finite(null$spearman_rho))
  stats$p_positive <- empirical_pvalue(stats$spearman_rho, null$spearman_rho)
  stats$null_mean <- mean(null$spearman_rho)
  stats$null_sd <- stats::sd(null$spearman_rho)
  if (adjusted_defined) {
    stats$n_adjusted_permutations_valid <- sum(is.finite(
      null$within_sample_adjusted_rho))
    stats$adjusted_p_positive <- empirical_pvalue(
      stats$within_sample_adjusted_rho, null$within_sample_adjusted_rho)
    stats$adjusted_null_mean <- mean(null$within_sample_adjusted_rho)
    stats$adjusted_null_sd <- stats::sd(null$within_sample_adjusted_rho)
  }
  list(statistics = stats, null = null)
}

# Asymptotic positive-tail Spearman log-p-values for diagnostics.
spearman_nominal_logp <- function(rho, n) {
  out <- rep(NA_real_, length(rho))
  valid <- is.finite(rho) & is.finite(n) & n > 2 & abs(rho) <= 1
  r <- rho[valid]
  t_stat <- r * sqrt((n[valid] - 2) / pmax(0, 1 - r^2))
  out[valid] <- pt(t_stat, df = n[valid] - 2, lower.tail = FALSE, log.p = TRUE)
  out
}

# Grouped descriptive summaries; these do not perform tests.
aggregate_donors <- function(d) {
  cells <- c("co_closed", "CRE1_access", "CRE2_access", "co_access")
  out <- d[, lapply(.SD, sum), by = .(CRE_pair, chrom, donor), .SDcols = cells]
  setorder(out, chrom, CRE_pair, donor)
  out
}

mean_props <- function(d, by_columns) {
  d[, .(n_fibers = .N, mean_prop1 = mean(prop1), mean_prop2 = mean(prop2)), by = by_columns]
}

complete_summary <- function(grid, summary, by_columns) {
  z <- merge(grid, summary, by = by_columns, all.x = TRUE, sort = FALSE)
  z[is.na(n_fibers), n_fibers := 0L]
  z
}

# Donor-adjusted binomial LRT / separated-data Rao test with Firth effects.
capture_fit <- function(expr) {
  warnings <- character()
  error <- NULL
  value <- tryCatch(withCallingHandlers(suppressMessages(expr), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) { error <<- conditionMessage(e); NULL })
  list(value = value, error = error, warnings = unique(warnings))
}

fit_pair <- function(d, min_reads = 2L) {
  cells <- c("co_closed", "CRE1_access", "CRE2_access", "co_access")
  n <- d$co_closed + d$CRE1_access + d$CRE2_access + d$co_access
  x_open <- d$CRE1_access + d$co_access
  y_open <- d$CRE2_access + d$co_access
  informative <- n >= min_reads & x_open > 0 & x_open < n & y_open > 0 & y_open < n
  z <- d[informative]
  ni <- n[informative]
  xi <- x_open[informative]
  yi <- y_open[informative]
  sums <- vapply(cells, function(column) sum(z[[column]]), numeric(1))
  expected <- sum(expected_coaccess(z$co_closed, z$CRE1_access, z$CRE2_access, z$co_access))
  result <- data.table(CRE_pair = d$CRE_pair[1],
    n_shared_reads = sum(ni), n_shared_reads_before_filter = sum(n),
    n_donors_total = nrow(d), n_donors_with_reads = sum(n > 0),
    n_donors_informative = nrow(z),
    co_closed = sums[["co_closed"]], CRE1_access = sums[["CRE1_access"]],
    CRE2_access = sums[["CRE2_access"]], co_access = sums[["co_access"]],
    expected_co_access_adj = if (nrow(z)) expected else NA_real_,
    observed_expected_ratio_adj = if (expected > 0) sums[["co_access"]] / expected else NA_real_,
    beta = NA_real_, OR_adj = NA_real_, OR_adj_ci_low = NA_real_, OR_adj_ci_high = NA_real_,
    ci_method = NA_character_, or_method = NA_character_, test = NA_character_,
    lrt_stat = NA_real_, score_stat = NA_real_, test_df = NA_integer_,
    p_LRT = NA_real_, p_Rao = NA_real_, p_adj_one_sided = NA_real_,
    fit_status = if (sum(n) == 0) "no_coverage" else "no_informative_donors",
    separation = "none", fit_warning = "", fit_error = "", ci_warning = "")
  if (!nrow(z)) return(result)
  separation_direction <- if (all(z$CRE1_access == 0 | z$CRE2_access == 0)) "positive" else
    if (all(z$co_access == 0 | z$co_closed == 0)) "negative" else "none"
  separated <- separation_direction != "none"
  result[, `:=`(separation = separation_direction,
                 or_method = if (separated) "firth" else "mle",
                 test = if (separated) "Rao_score" else "LRT")]
  fail <- function(reason) {
    result[, `:=`(fit_status = "fit_failed", fit_error = reason)]
    result
  }
  binomial_rows <- data.frame(donor = factor(rep(z$donor, each = 2L)),
    x = rep(c(0, 1), nrow(z)),
    y1 = as.vector(rbind(z$CRE2_access, z$co_access)),
    y0 = as.vector(rbind(z$co_closed, z$CRE1_access)))
  null_formula <- if (nrow(z) == 1L) cbind(y1, y0) ~ 1 else cbind(y1, y0) ~ factor(donor)
  full_formula <- if (nrow(z) == 1L) cbind(y1, y0) ~ x else cbind(y1, y0) ~ factor(donor) + x
  null_fit <- capture_fit(glm(null_formula, data = binomial_rows, family = binomial(),
                              control = glm.control(maxit = 100L)))
  full_fit <- capture_fit(glm(full_formula, data = binomial_rows, family = binomial(),
                              control = glm.control(maxit = 100L)))
  result[, fit_warning := paste(unique(c(null_fit$warnings, full_fit$warnings)), collapse = " | ")]
  if (!is.null(null_fit$error) || !is.null(full_fit$error))
    return(fail(paste(c(null_fit$error, full_fit$error), collapse = " | ")))
  null <- null_fit$value
  full <- full_fit$value
  # Rao uses the converged NULL residuals/weights and FULL design matrix;
  # it does not require a finite full-model MLE. No separated GLM effect is reported.
  if (!isTRUE(null$converged) || (!separated && !isTRUE(full$converged)))
    return(fail("GLM did not converge"))
  if (full$rank - null$rank != 1L) return(fail("Full/null model rank difference is not one"))
  comparison <- capture_fit(anova(null, full, test = if (separated) "Rao" else "LRT"))
  if (!is.null(comparison$error)) return(fail(comparison$error))
  tab <- comparison$value
  p_column <- grep("^Pr\\(", names(tab), value = TRUE)
  statistic <- if (separated) tab$Rao[2] else tab$Deviance[2]
  if (length(p_column) != 1L || !is.finite(tab[[p_column]][2]) ||
      !is.finite(statistic) || statistic < 0 || tab$Df[2] != 1L)
    return(fail("Invalid likelihood-ratio or Rao test"))
  p_two_sided <- tab[[p_column]][2]
  if (separated) {
    weighted <- data.frame(donor = factor(rep(z$donor, each = 4L)),
      x = rep(c(0, 0, 1, 1), nrow(z)), y = rep(c(0, 1, 0, 1), nrow(z)),
      frequency = as.vector(rbind(z$co_closed, z$CRE2_access, z$CRE1_access, z$co_access)))
    weighted <- weighted[weighted$frequency > 0, , drop = FALSE]
    firth_formula <- if (nrow(z) == 1L) y ~ x else y ~ factor(donor) + x
    x_index <- match("x", colnames(model.matrix(firth_formula, weighted)))
    firth <- capture_fit(logistf::logistf(firth_formula, data = weighted, weights = frequency,
      firth = TRUE, pl = TRUE, plconf = x_index, alpha = 0.05,
      control = logistf::logistf.control(maxit = 100L),
      plcontrol = logistf::logistpl.control(maxit = 100L)))
    result[, fit_warning := paste(c(fit_warning, firth$warnings), collapse = " | ")]
    if (!is.null(firth$error)) return(fail(paste("Firth:", firth$error)))
    f <- firth$value
    beta_firth <- unname(f$coefficients["x"])
    limits <- c(f$ci.lower["x"], f$ci.upper["x"])
    if (f$iter[["full"]] >= 100L || any(f$pl.iter[x_index, 1:2] >= 100L) ||
        !is.finite(beta_firth) || anyNA(limits) ||
        limits[1] > beta_firth || limits[2] < beta_firth ||
        f$method.ci[x_index] != "Profile Likelihood")
      return(fail("Firth fit or penalized profile interval did not converge"))
    result[, `:=`(beta = beta_firth, OR_adj = exp(beta_firth),
      OR_adj_ci_low = exp(limits[1]), OR_adj_ci_high = exp(limits[2]),
      ci_method = "firth_penalized_profile", score_stat = statistic, test_df = 1L,
      p_Rao = p_two_sided,
      p_adj_one_sided = if (separation_direction == "positive") p_two_sided / 2 else 1 - p_two_sided / 2,
      fit_status = paste0("separation_", separation_direction))]
    return(result)
  }
  beta_hat <- unname(coef(full)["x"])
  if (!is.finite(beta_hat)) return(fail("Coefficient x is not estimable"))
  profile <- capture_fit(stats::confint(full, parm = "x", level = 0.95, test = "LRT"))
  interval <- as.numeric(profile$value)
  profile_ok <- is.null(profile$error) && !length(profile$warnings) &&
    length(interval) == 2L && !anyNA(interval) && interval[1] <= beta_hat && interval[2] >= beta_hat
  ci_kind <- "profile_likelihood"
  if (!profile_ok) {
    ci_kind <- "wald_fallback"
    result[, ci_warning := paste(c(profile$error, profile$warnings,
      if (is.null(profile$error) && !length(profile$warnings)) "Invalid profile interval"), collapse = " | ")]
    wald <- capture_fit(sqrt(vcov(full)["x", "x"]))
    if (!is.null(wald$error) || length(wald$value) != 1L || !is.finite(wald$value)) {
      return(fail("Profile and Wald confidence intervals failed"))
    }
    interval <- beta_hat + c(-1, 1) * qnorm(0.975) * wald$value
  }
  result[, `:=`(beta = beta_hat, OR_adj = exp(beta_hat),
    OR_adj_ci_low = exp(interval[1]), OR_adj_ci_high = exp(interval[2]),
    ci_method = ci_kind, lrt_stat = statistic, test_df = 1L, p_LRT = p_two_sided,
    p_adj_one_sided = if (beta_hat > 0) p_two_sided / 2 else 1 - p_two_sided / 2,
    fit_status = if (nzchar(fit_warning)) "ok_with_fit_warnings" else "ok")]
  result
}
