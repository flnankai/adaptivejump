.aj_m4k <- function(k) {
  16 * k * (2 * k^2 - k - 1) / 3
}

.aggregate_disjoint <- function(x, k) {
  n <- length(x)
  usable <- n - n %% k
  if (usable < k) {
    stop("The path is too short for the requested block size.", call. = FALSE)
  }
  colSums(matrix(x[seq_len(usable)], nrow = k))
}

.lm_window_rule <- function(n_returns, exponent = 0.6) {
  as.integer(max(3L, min(n_returns - 1L, ceiling(n_returns^exponent))))
}

.lm_gumbel_normalization <- function(n_gumbel, c_eff) {
  n_gumbel <- max(2L, as.integer(n_gumbel))
  cn <- sqrt(2 * log(n_gumbel)) / c_eff -
    (log(pi) + log(log(n_gumbel))) / (2 * c_eff * sqrt(2 * log(n_gumbel)))
  sn <- 1 / (c_eff * sqrt(2 * log(n_gumbel)))
  list(cn = cn, sn = sn)
}

.lm_components <- function(log_returns, K) {
  n <- length(log_returns)
  if (K < 3L || K >= n) {
    stop("`K` must satisfy 3 <= K < number of returns.", call. = FALSE)
  }
  bipower <- abs(log_returns[-1L]) * abs(log_returns[-n])
  window <- K - 2L
  csum_bp <- c(0, cumsum(bipower))
  roll <- csum_bp[(window + 1L):length(csum_bp)] -
    csum_bp[seq_len(length(csum_bp) - window)]
  sigma2 <- roll[seq_len(n - K + 1L)] / window

  csum_ret <- c(0, cumsum(log_returns))
  start <- seq_len(n - K + 1L)
  end <- (K - 1L):(n - 1L)
  local_mean <- (csum_ret[end + 1L] - csum_ret[start]) / (K - 1L)
  list(sigma2 = sigma2, local_mean = local_mean)
}

.lm_statistic <- function(log_returns, K, n_gumbel = length(log_returns)) {
  comp <- .lm_components(log_returns, K)
  l_values <- (log_returns[K:length(log_returns)] - comp$local_mean) /
    sqrt(pmax(comp$sigma2, .Machine$double.eps))
  c_eff <- sqrt(2 / pi) * sqrt((K - 1) / K)
  normalization <- .lm_gumbel_normalization(n_gumbel, c_eff)
  max_at <- which.max(abs(l_values))
  list(
    statistic = (abs(l_values[max_at]) - normalization$cn) / normalization$sn,
    max_standardized_return = abs(l_values[max_at]),
    max_index = as.integer(K - 1L + max_at),
    standardized_returns = l_values,
    sigma2 = comp$sigma2,
    local_mean = comp$local_mean,
    cn = normalization$cn,
    sn = normalization$sn
  )
}

.lm_bootstrap_statistics <- function(log_returns, K, n_gumbel, bootstrap_rep, seed = NULL) {
  bootstrap_rep <- as.integer(bootstrap_rep)
  if (bootstrap_rep < 20L) {
    stop("`bootstrap_rep` must be at least 20.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(as.integer(abs(seed) %% (.Machine$integer.max - 1L)))
  }
  comp <- .lm_components(log_returns, K)
  sigma_tail <- sqrt(pmax(comp$sigma2, .Machine$double.eps))
  sigma_full <- c(rep(sigma_tail[1L], K - 1L), sigma_tail)

  # All bootstrap paths are generated together; no batching changes the draw.
  z <- matrix(stats::rnorm(bootstrap_rep * length(log_returns)), nrow = bootstrap_rep)
  boot_returns <- sweep(z, 2L, sigma_full, `*`)
  vapply(
    seq_len(bootstrap_rep),
    function(i) .lm_statistic(boot_returns[i, ], K, n_gumbel)$statistic,
    numeric(1L)
  )
}

#' Ait-Sahalia--Jacod power-variation jump test
#'
#' Implements the (p=4) frictionless power-variation ratio statistic used in
#' the paper.  The alternative is a discontinuous path, so the p-value is from
#' the lower tail of the standardized ratio.
#'
#' @param prices A numeric price path, or a log-price path when
#'   `log_prices = FALSE`.
#' @param delta Sampling interval on the normalized time scale.  By default it
#'   is `1 / (length(prices) - 1)`.
#' @param k Integer aggregation block size, at least two.
#' @param log_prices Should raw prices be log transformed first?
#' @param truncation_constant Constant multiplying `delta^truncation_exponent`
#'   in the truncated quarticity estimators.
#' @param truncation_exponent Exponent in the truncation threshold.
#' @param alpha Optional decision level.
#'
#' @return An object of class `jump_test` with statistic `S_n`, standardized
#'   statistic `z_statistic`, and lower-tail p-value.
#' @export
aj_test <- function(prices,
                    delta = NULL,
                    k = 2L,
                    log_prices = TRUE,
                    truncation_constant = 2,
                    truncation_exponent = 0.47,
                    alpha = NULL) {
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 2L) {
    stop("`k` must be an integer of at least two.", call. = FALSE)
  }
  x <- .as_numeric_path(prices, log_prices = log_prices)
  returns <- diff(x)
  usable <- length(returns) - length(returns) %% k
  returns <- returns[seq_len(usable)]
  delta <- .resolve_delta(length(returns), delta)

  b4 <- sum(abs(returns)^4)
  if (!is.finite(b4) || b4 <= 0) {
    stop("The path has no nonzero fourth power variation.", call. = FALSE)
  }
  coarse_returns <- .aggregate_disjoint(returns, k)
  s_statistic <- sum(abs(coarse_returns)^4) / b4

  threshold <- truncation_constant * delta^truncation_exponent
  keep <- abs(returns) <= threshold
  m4 <- .normal_abs_moment(4)
  m8 <- .normal_abs_moment(8)
  a4 <- delta^(-1) * sum(abs(returns)^4 * keep) / m4
  a8 <- delta^(-3) * sum(abs(returns)^8 * keep) / m8
  variance_hat <- delta * .aj_m4k(k) * a8 / pmax(a4^2, .Machine$double.eps)
  z_statistic <- (s_statistic - k) / sqrt(pmax(variance_hat, .Machine$double.eps))
  p_value <- .clip_p(stats::pnorm(z_statistic))

  .new_jump_test(
    "AJ", s_statistic, p_value, alpha = alpha,
    z_statistic = z_statistic, k = k, delta = delta,
    variance_hat = variance_hat, truncation_threshold = threshold,
    n_returns = length(returns)
  )
}

#' Lee--Mykland max jump test
#'
#' Computes the demeaned, locally standardized Lee--Mykland statistic.  Its
#' calibration may use the asymptotic Gumbel distribution or the pathwise
#' parametric bootstrap used in the paper's finite-sample simulations.
#'
#' @param prices A numeric price path, or a log-price path when
#'   `log_prices = FALSE`.
#' @param K Local volatility and demeaning window.  Defaults to
#'   `ceiling(n^window_exponent)`.
#' @param log_prices Should raw prices be log transformed first?
#' @param window_exponent Exponent used by the default `K` rule.
#' @param calibration Either `"bootstrap"` or `"gumbel"`.
#' @param bootstrap_rep Number of pathwise bootstrap replications.
#' @param seed Optional random seed for the bootstrap.
#' @param alpha Optional decision level.
#'
#' @return An object of class `jump_test`.
#' @export
lm_test <- function(prices,
                    K = NULL,
                    log_prices = TRUE,
                    window_exponent = 0.6,
                    calibration = c("bootstrap", "gumbel"),
                    bootstrap_rep = 199L,
                    seed = NULL,
                    alpha = NULL) {
  calibration <- match.arg(calibration)
  x <- .as_numeric_path(prices, log_prices = log_prices)
  returns <- diff(x)
  n_returns <- length(returns)
  if (is.null(K)) {
    K <- .lm_window_rule(n_returns, window_exponent)
  }
  K <- as.integer(K)
  stat <- .lm_statistic(returns, K = K, n_gumbel = n_returns)

  if (calibration == "gumbel") {
    p_value <- .clip_p(1 - exp(-exp(-stat$statistic)))
    critical_values <- stats::qlogis(c(0.90, 0.95), location = 0, scale = 1)
    boot_n <- 0L
  } else {
    boot_stats <- .lm_bootstrap_statistics(
      returns, K = K, n_gumbel = n_returns,
      bootstrap_rep = bootstrap_rep, seed = seed
    )
    p_value <- .clip_p((1 + sum(boot_stats >= stat$statistic)) / (length(boot_stats) + 1))
    critical_values <- as.numeric(stats::quantile(
      boot_stats, probs = c(0.90, 0.95), names = FALSE, type = 8
    ))
    boot_n <- length(boot_stats)
  }

  .new_jump_test(
    "LM", stat$statistic, p_value, alpha = alpha,
    max_standardized_return = stat$max_standardized_return,
    max_index = stat$max_index, K = K, n_returns = n_returns,
    calibration = calibration, critical_10 = critical_values[1L],
    critical_05 = critical_values[2L], bootstrap_rep = boot_n,
    gumbel_center = stat$cn, gumbel_scale = stat$sn
  )
}

#' Adaptive frictionless jump test
#'
#' Applies the AJ and LM tests to one path and combines their p-values using
#' the equal-weight Cauchy combination rule.
#'
#' @inheritParams aj_test
#' @param lm_K Local window for `lm_test()`.
#' @param lm_calibration Calibration for `lm_test()`.
#' @param bootstrap_rep Number of LM bootstrap replications when applicable.
#' @param seed Optional LM bootstrap seed.
#'
#' @return A list with components `aj`, `lm`, and `cc`.
#' @export
adaptive_jump_test <- function(prices,
                               delta = NULL,
                               k = 2L,
                               lm_K = NULL,
                               log_prices = TRUE,
                               lm_calibration = c("bootstrap", "gumbel"),
                               bootstrap_rep = 199L,
                               seed = NULL,
                               alpha = NULL) {
  lm_calibration <- match.arg(lm_calibration)
  aj <- aj_test(
    prices, delta = delta, k = k, log_prices = log_prices, alpha = alpha
  )
  lm <- lm_test(
    prices, K = lm_K, log_prices = log_prices,
    calibration = lm_calibration, bootstrap_rep = bootstrap_rep,
    seed = seed, alpha = alpha
  )
  cc <- cauchy_combine(c(aj$p_value, lm$p_value), alpha = alpha)
  list(aj = aj, lm = lm, cc = cc)
}
