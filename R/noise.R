.build_rho <- function(p) {
  if (length(p) != 1L || p %% 2L != 0L || p < 4L) {
    stop("`p` must be an even integer of at least four.", call. = FALSE)
  }
  half_p <- as.integer(p / 2L)
  rho <- numeric(half_p + 1L)
  rho[1L] <- 1
  moments <- vapply(0:p, .normal_abs_moment, numeric(1L))
  for (j in seq_len(half_p)) {
    value <- 0
    for (ell in 0:(j - 1L)) {
      value <- value + 2^ell * moments[2L * j - 2L * ell + 1L] *
        choose(p - 2L * ell, p - 2L * j) * rho[ell + 1L]
    }
    rho[j + 1L] <- -value / 2^j
  }
  rho
}

.g_weight <- function(x) pmax(0.5 - abs(x - 0.5), 0)
.h_weight <- function(x) .g_weight(2 * x)

.weight_integral <- function(weight, power) {
  stats::integrate(
    function(u) abs(weight(u))^power, lower = 0, upper = 1,
    subdivisions = 200L, rel.tol = 1e-10
  )$value
}

.make_ajj_setup <- function(kn, p = 4L) {
  kn <- as.integer(kn)
  if (kn < 3L) stop("`kn` must be at least three.", call. = FALSE)
  g_grid <- .g_weight((0:kn) / kn)
  h_grid <- .h_weight((0:kn) / kn)
  g2 <- .weight_integral(.g_weight, 2)
  h2 <- .weight_integral(.h_weight, 2)
  gp <- .weight_integral(.g_weight, p)
  hp <- .weight_integral(.h_weight, p)
  gamma <- g2 / h2
  gamma_prime <- gp / hp
  list(
    kn = kn,
    p = as.integer(p),
    rho = .build_rho(p),
    g_bar_rev = rev(g_grid[2:kn]),
    h_bar_rev = rev(h_grid[2:kn]),
    g_dsq_rev = rev(diff(g_grid)^2),
    h_dsq_rev = rev(diff(h_grid)^2),
    gamma = gamma,
    gamma_prime = gamma_prime,
    gamma_doubleprime = gamma^(p / 2) / gamma_prime,
    gamma_p2 = gamma^(p / 2)
  )
}

.preaverage_components <- function(dy, setup) {
  emb <- embed(dy, setup$kn)
  list(
    ybar_g = as.vector(emb[, -1L, drop = FALSE] %*% setup$g_bar_rev),
    yhat_g = as.vector((emb^2) %*% setup$g_dsq_rev),
    ybar_h = as.vector(emb[, -1L, drop = FALSE] %*% setup$h_bar_rev),
    yhat_h = as.vector((emb^2) %*% setup$h_dsq_rev)
  )
}

.preaverage_psi <- function(ybar, yhat, rho, p) {
  value <- numeric(length(ybar))
  for (ell in 0:(p / 2L)) {
    value <- value + rho[ell + 1L] * abs(ybar)^(p - 2L * ell) * yhat^ell
  }
  value
}

.compute_ajj <- function(y, delta, rn, setup) {
  dy <- diff(y)
  n <- length(dy)
  if (setup$kn + 5L >= n) {
    stop("The path is too short for the selected pre-averaging window.", call. = FALSE)
  }
  rn <- as.integer(rn)
  if (rn <= setup$kn) {
    stop("`rn` must exceed `kn`.", call. = FALSE)
  }
  comp <- .preaverage_components(dy, setup)
  psi_g <- .preaverage_psi(comp$ybar_g, comp$yhat_g, setup$rho, setup$p)
  psi_h <- .preaverage_psi(comp$ybar_h, comp$yhat_h, setup$rho, setup$p)
  vg <- sum(psi_g)
  vh <- sum(psi_h)
  if (!is.finite(vh) || vh <= 0) {
    stop("The pre-averaged denominator is nonpositive.", call. = FALSE)
  }
  raw_ratio <- vg / (setup$gamma_prime * vh)
  psi_g_scaled <- delta^(setup$p / 4 - 1) * psi_g
  psi_h_scaled <- delta^(setup$p / 4 - 1) * psi_h
  denominator <- setup$gamma_prime * sum(psi_h_scaled)
  if (!is.finite(denominator) || denominator == 0) {
    stop("The pre-averaged linearization denominator is zero.", call. = FALSE)
  }
  q <- delta^(-1 / 4) * (psi_g_scaled - setup$gamma_p2 * psi_h_scaled) / denominator
  n_blocks <- floor(length(q) / rn)
  if (n_blocks < 2L) {
    stop("The selected `rn` leaves fewer than two self-normalization blocks.", call. = FALSE)
  }
  gamma_blocks <- vapply(seq_len(n_blocks), function(b) {
    left <- (b - 1L) * rn + 1L
    right <- min(b * rn - setup$kn, length(q))
    if (right < left) 0 else sum(q[left:right])
  }, numeric(1L))
  denominator_hat <- sqrt(sum(gamma_blocks^2))
  if (!is.finite(denominator_hat) || denominator_hat <= 0) {
    stop("The AJJ self-normalizer is nonpositive.", call. = FALSE)
  }
  statistic <- sum(gamma_blocks) / denominator_hat
  list(
    statistic = statistic,
    p_value = .clip_p(2 * stats::pnorm(-abs(statistic))),
    raw_ratio = raw_ratio,
    gamma_doubleprime = setup$gamma_doubleprime,
    n_blocks = n_blocks,
    self_normalizer = denominator_hat
  )
}

#' Noise-robust pre-averaged AJJ test
#'
#' Implements the pre-averaged ratio statistic used for the noisy sum
#' component.  The default tuning rates follow the paper; exact simulation
#' settings can be supplied through `kn` and `rn`.
#'
#' @param prices A numeric price path, or a log-price path when
#'   `log_prices = FALSE`.
#' @param delta Sampling interval on the normalized time scale.
#' @param kn Pre-averaging window.  Defaults to `floor(sqrt(n_returns))`.
#' @param rn Block length of the self-normalizer.  Defaults to
#'   `floor(n_returns^0.85)` subject to feasibility.
#' @param log_prices Should raw prices be log transformed first?
#' @param p Even power, currently implemented for the paper's default `p = 4`.
#' @param alpha Optional decision level.
#'
#' @return An object of class `jump_test`.
#' @export
ajj_test <- function(prices,
                     delta = NULL,
                     kn = NULL,
                     rn = NULL,
                     log_prices = TRUE,
                     p = 4L,
                     alpha = NULL) {
  y <- .as_numeric_path(prices, log_prices = log_prices)
  n_returns <- length(y) - 1L
  delta <- .resolve_delta(n_returns, delta)
  if (is.null(kn)) kn <- floor(sqrt(n_returns))
  kn <- as.integer(kn)
  if (is.null(rn)) rn <- max(kn + 1L, floor(n_returns^0.85))
  rn <- as.integer(min(rn, floor((n_returns - kn + 1L) / 2L)))
  setup <- .make_ajj_setup(kn, p = p)
  result <- .compute_ajj(y, delta = delta, rn = rn, setup = setup)
  .new_jump_test(
    "AJJ", result$statistic, result$p_value, alpha = alpha,
    raw_ratio = result$raw_ratio, gamma_doubleprime = result$gamma_doubleprime,
    kn = kn, rn = rn, n_blocks = result$n_blocks, delta = delta,
    self_normalizer = result$self_normalizer
  )
}

.rolling_mean_left <- function(x, k) {
  if (k < 1L || k > length(x)) {
    stop("Invalid rolling-mean window.", call. = FALSE)
  }
  cs <- c(0, cumsum(x))
  (cs[(k + 1L):(length(x) + 1L)] - cs[seq_len(length(x) - k + 1L)]) / k
}

.subsample_prices <- function(y, lag) {
  n <- length(y) - 1L
  n_sub <- as.integer(floor(n / lag))
  y[1L + (0L:n_sub) * lag]
}

.estimate_noise_sd_path <- function(y, lag = 1L) {
  lag <- as.integer(lag)
  n <- length(y) - 1L
  if (lag < 1L || n <= lag) return(NA_real_)
  d <- y[(1L + lag):(n + 1L)] - y[seq_len(n + 1L - lag)]
  sqrt(mean(d^2) / 2)
}

#' Estimate additive noise standard deviation
#'
#' Uses the Lee--Mykland (2012) lag-difference estimator
#' \eqn{\widehat q = \{(2(n-k))^{-1}\sum(Y_{i+k}-Y_i)^2\}^{1/2}}.
#'
#' @param prices A numeric price path, or a log-price path when
#'   `log_prices = FALSE`.
#' @param lag Difference lag, fixed at one for the i.i.d. noise setting in the
#'   paper.
#' @param log_prices Should raw prices be log transformed first?
#'
#' @return Estimated noise standard deviation.
#' @export
estimate_noise_sd <- function(prices, lag = 1L, log_prices = TRUE) {
  y <- .as_numeric_path(prices, log_prices = log_prices)
  .estimate_noise_sd_path(y, lag = lag)
}

.lm2012_optimal_c <- function(q_percent) {
  q_grid <- c(0.01, 0.03, 0.05, 0.07, 0.09, 0.10, 0.20, 0.30,
              0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)
  c_grid <- c(1 / 19, 1 / 19, 1 / 19, 1 / 18, 1 / 18, 1 / 18,
              1 / 16, 1 / 16, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9,
              1 / 8, 1 / 8)
  q_percent <- min(max(q_percent, q_grid[1L]), q_grid[length(q_grid)])
  stats::approx(q_grid, c_grid, xout = q_percent, method = "linear", rule = 2)$y
}

.select_lm_block <- function(n, lag, q_hat) {
  n_sub <- as.integer(max(1L, floor(n / lag)))
  q_percent <- 100 * q_hat
  c_value <- .lm2012_optimal_c(q_percent)
  block <- as.integer(floor(c_value * sqrt(n_sub)))
  block <- max(2L, min(block, max(2L, floor(n_sub / 4L))))
  list(block = block, c_value = c_value, q_percent = q_percent, n_sub = n_sub)
}

.local_average_differences <- function(y_sub, block) {
  n_sub <- length(y_sub) - 1L
  pb <- .rolling_mean_left(y_sub, block)
  if (length(pb) <= block) {
    return(list(L = numeric(), grid = integer(), n_sub = n_sub))
  }
  grid <- seq.int(1L, length(pb) - block, by = 2L * block)
  list(L = pb[grid + block] - pb[grid], grid = grid - 1L, n_sub = n_sub)
}

.bound_tsrsv <- function(n_sub, window_n, K) {
  window_n <- as.integer(min(max(20L, window_n), max(21L, n_sub - 1L)))
  K <- as.integer(min(max(2L, K), max(2L, floor(window_n / 4L))))
  list(window_n = window_n, K = K)
}

.select_tsrsv_rate <- function(n_sub, delta_sub, window_n = NULL, K = NULL,
                               window_scale = 1 / 4, K_scale = 1) {
  horizon <- n_sub * delta_sub
  if (is.null(window_n)) {
    window_n <- floor(window_scale * (horizon / delta_sub)^(5 / 6))
  }
  if (is.null(K)) {
    K <- floor(K_scale * (horizon / delta_sub)^(5 / 9))
  }
  .bound_tsrsv(n_sub, window_n, K)
}

.tsrsv_iv <- function(y_window, K) {
  m <- length(y_window) - 1L
  if (m <= K || K < 1L) return(NA_real_)
  rv_fast <- sum(diff(y_window)^2)
  rv_slow <- sum((y_window[(K + 1L):(m + 1L)] - y_window[seq_len(m + 1L - K)])^2) / K
  nbar <- (m - K + 1L) / K
  rv_slow - (nbar / m) * rv_fast
}

.spot_sigma2_vec <- function(y_sub, delta_sub, centers, window_n, K) {
  n_sub <- length(y_sub) - 1L
  if (!length(centers) || n_sub <= window_n || window_n <= K) {
    return(rep(NA_real_, length(centers)))
  }
  left <- pmax(0L, pmin(as.integer(centers) - floor(window_n / 2L), n_sub - window_n))
  right <- left + window_n
  r1_sq <- diff(y_sub)^2
  cs_fast <- c(0, cumsum(r1_sq))
  rv_fast <- cs_fast[right + 1L] - cs_fast[left + 1L]
  dK_sq <- (y_sub[(K + 1L):(n_sub + 1L)] - y_sub[seq_len(n_sub - K + 1L)])^2
  cs_slow <- c(0, cumsum(dK_sq))
  rv_slow <- (cs_slow[right - K + 2L] - cs_slow[left + 1L]) / K
  nbar <- (window_n - K + 1L) / K
  sigma2 <- (rv_slow - (nbar / window_n) * rv_fast) / (window_n * delta_sub)
  sigma2[!is.finite(sigma2) | sigma2 <= 0] <- NA_real_
  sigma2
}

.fill_positive <- function(x, fallback = NA_real_) {
  valid <- is.finite(x) & x > 0
  if (!any(valid)) return(rep(fallback, length(x)))
  if (sum(valid) == 1L) return(rep(x[valid][1L], length(x)))
  stats::approx(which(valid), x[valid], xout = seq_along(x), method = "linear", rule = 2)$y
}

.integrated_quarticity_sparse <- function(y_sub, delta_sub, seconds_per_day = 23400) {
  n_sub <- length(y_sub) - 1L
  if (n_sub < 2L) return(NA_real_)
  step <- as.integer(max(1L, round(300 / (delta_sub * seconds_per_day))))
  if (n_sub <= step) return(NA_real_)
  idx <- seq.int(1L, n_sub + 1L - step, by = step)
  r_sparse <- y_sub[idx + step] - y_sub[idx]
  sum(r_sparse^4) / (3 * step * delta_sub)
}

.prelim_lambda2 <- function(y_sub, delta_sub, window_n, K) {
  sigma2 <- .spot_sigma2_vec(y_sub, delta_sub, 0:(length(y_sub) - 2L), window_n, K)
  sigma2 <- .fill_positive(sigma2, fallback = NA_real_)
  if (length(sigma2) < 2L || any(!is.finite(sigma2))) return(NA_real_)
  sum(diff(sigma2)^2)
}

.select_tsrsv <- function(n_sub, y_sub, delta_sub, q_hat,
                          window_n = NULL, K = NULL,
                          window_scale = 1 / 4, K_scale = 1,
                          tuning = c("rate", "plugin_eq56"),
                          seconds_per_day = 23400) {
  tuning <- match.arg(tuning)
  rate <- .select_tsrsv_rate(n_sub, delta_sub, window_n, K, window_scale, K_scale)
  rate$tuning <- if (!is.null(window_n) || !is.null(K)) "fixed_override" else "rate"
  rate$K_star_hat <- NA_real_
  rate$h_star_hat <- NA_real_
  rate$iq_hat <- NA_real_
  rate$lambda2_hat <- NA_real_
  if (tuning == "rate" || !is.finite(q_hat) || q_hat <= 0 ||
      !is.null(window_n) || !is.null(K)) return(rate)

  omega2 <- q_hat^2
  iq_hat <- .integrated_quarticity_sparse(y_sub, delta_sub, seconds_per_day)
  lambda2_hat <- .prelim_lambda2(y_sub, delta_sub, rate$window_n, rate$K)
  if (!is.finite(iq_hat) || iq_hat <= 0 || !is.finite(lambda2_hat) || lambda2_hat <= 0) {
    rate$tuning <- "rate_fallback"
    rate$iq_hat <- iq_hat
    rate$lambda2_hat <- lambda2_hat
    return(rate)
  }
  K_star <- (12 * omega2^2 / iq_hat)^(1 / 3)
  h_star <- sqrt((8 * K_star^2 * omega2^2 + (4 / 3) * K_star * iq_hat) /
                   ((1 / 3) * lambda2_hat))
  if (!is.finite(K_star) || K_star <= 0 || !is.finite(h_star) || h_star <= 0) {
    rate$tuning <- "rate_fallback"
    rate$K_star_hat <- K_star
    rate$h_star_hat <- h_star
    rate$iq_hat <- iq_hat
    rate$lambda2_hat <- lambda2_hat
    return(rate)
  }
  horizon <- n_sub * delta_sub
  tuned <- .bound_tsrsv(
    n_sub,
    floor(h_star * (horizon / delta_sub)^(5 / 6)),
    floor(K_star * (horizon / delta_sub)^(2 / 3))
  )
  tuned$tuning <- "plugin_eq56"
  tuned$K_star_hat <- K_star
  tuned$h_star_hat <- h_star
  tuned$iq_hat <- iq_hat
  tuned$lambda2_hat <- lambda2_hat
  tuned
}

#' Estimate a TSRSV spot-volatility path
#'
#' Computes the filtering two-scale realized spot-variance estimator used by
#' the noisy LM standardization.
#'
#' @param prices A numeric price path, or a log-price path when
#'   `log_prices = FALSE`.
#' @param delta Sampling interval on the normalized time scale.
#' @param log_prices Should raw prices be log transformed first?
#' @param window_n Optional local window length.
#' @param K Optional two-scale subsampling lag.
#' @param window_scale Default rate constant for `window_n`.
#' @param K_scale Default rate constant for `K`.
#'
#' @return A list containing the estimated spot variance vector and tuning
#'   values.
#' @export
estimate_spot_volatility <- function(prices,
                                     delta = NULL,
                                     log_prices = TRUE,
                                     window_n = NULL,
                                     K = NULL,
                                     window_scale = 1 / 4,
                                     K_scale = 1) {
  y <- .as_numeric_path(prices, log_prices = log_prices)
  n_sub <- length(y) - 1L
  delta <- .resolve_delta(n_sub, delta)
  tune <- .select_tsrsv_rate(n_sub, delta, window_n, K, window_scale, K_scale)
  raw <- .spot_sigma2_vec(y, delta, 0:(n_sub - 1L), tune$window_n, tune$K)
  valid <- is.finite(raw) & raw > 0
  fallback <- if (any(valid)) mean(raw[valid]) else {
    iv <- .tsrsv_iv(y, tune$K)
    if (is.finite(iv) && iv > 0) iv / (n_sub * delta) else mean(diff(y)^2) / delta
  }
  list(
    sigma2 = .fill_positive(raw, fallback = max(fallback, .Machine$double.eps)),
    sigma2_raw = raw, window_n = tune$window_n, K = tune$K,
    valid_rate = mean(valid), delta = delta
  )
}

.noisy_lm_statistic <- function(y_sub, delta_sub, seconds_per_day = 23400,
                                dependence_lag = 1L,
                                spot_tuning = c("rate", "plugin_eq56"),
                                spot_window = NULL, spot_K = NULL,
                                spot_window_scale = 1 / 4, spot_K_scale = 1) {
  spot_tuning <- match.arg(spot_tuning)
  n_sub <- length(y_sub) - 1L
  q_hat <- .estimate_noise_sd_path(y_sub, lag = dependence_lag)
  if (!is.finite(q_hat) || q_hat <= 0) return(list(statistic = NA_real_))
  block_info <- .select_lm_block(n_sub, dependence_lag, q_hat)
  la <- .local_average_differences(y_sub, block_info$block)
  if (length(la$L) < 2L) return(list(statistic = NA_real_))

  tune <- .select_tsrsv(
    la$n_sub, y_sub, delta_sub, q_hat,
    window_n = spot_window, K = spot_K,
    window_scale = spot_window_scale, K_scale = spot_K_scale,
    tuning = spot_tuning, seconds_per_day = seconds_per_day
  )
  sigma2 <- .spot_sigma2_vec(y_sub, delta_sub, la$grid, tune$window_n, tune$K)
  valid <- is.finite(sigma2) & sigma2 > 0
  L <- la$L[valid]
  grid <- la$grid[valid]
  sigma2 <- sigma2[valid]
  if (length(L) < 2L) return(list(statistic = NA_real_))

  sigma2_boot_path <- .fill_positive(
    stats::approx(grid + 1L, sigma2, xout = seq_len(la$n_sub), method = "linear", rule = 2)$y,
    fallback = mean(sigma2)
  )
  vhat2 <- 2 * q_hat^2 + (2 / 3) * sigma2 * block_info$block^2 * delta_sub
  if (any(!is.finite(vhat2)) || any(vhat2 <= 0)) return(list(statistic = NA_real_))
  standardized <- sqrt(block_info$block) * L / sqrt(vhat2)
  max_at <- which.max(abs(standardized))
  max_std <- abs(standardized[max_at])
  n_theory <- as.integer(max(2L, floor((la$n_sub - 2L * block_info$block) /
                                        (2L * block_info$block)) + 1L))
  center <- sqrt(2 * log(n_theory)) -
    (log(pi) + log(log(n_theory))) / (2 * sqrt(2 * log(n_theory)))
  scale <- 1 / sqrt(2 * log(n_theory))
  list(
    statistic = (max_std - center) / scale,
    max_standardized_return = max_std,
    max_index = as.integer(grid[max_at] + 1L),
    q_hat = q_hat,
    q_percent = block_info$q_percent,
    dependence_lag = dependence_lag,
    block_size = block_info$block,
    C = block_info$c_value,
    sigma2_hat = mean(sigma2), sigma2_min = min(sigma2), sigma2_max = max(sigma2),
    omega2_hat = q_hat^2, spot_window = tune$window_n, spot_K = tune$K,
    spot_tuning = tune$tuning, K_star_hat = tune$K_star_hat,
    h_star_hat = tune$h_star_hat, iq_hat = tune$iq_hat,
    lambda2_hat = tune$lambda2_hat, n_tests = length(L), n_tests_theory = n_theory,
    vhat2 = mean(vhat2), delta_sub = delta_sub, sigma2_boot_path = sigma2_boot_path,
    gumbel_center = center, gumbel_scale = scale
  )
}

.noisy_lm_bootstrap <- function(y_sub, delta_sub, q_hat, sigma2_path,
                                bootstrap_rep, seconds_per_day,
                                dependence_lag, spot_tuning, spot_window, spot_K,
                                spot_window_scale, spot_K_scale, seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(abs(seed) %% (.Machine$integer.max - 1L)))
  bootstrap_rep <- as.integer(bootstrap_rep)
  if (bootstrap_rep < 20L) stop("`bootstrap_rep` must be at least 20.", call. = FALSE)
  n_sub <- length(y_sub) - 1L
  latent_sd <- sqrt(pmax(sigma2_path * delta_sub, .Machine$double.eps))
  result <- numeric(bootstrap_rep)
  accepted <- 0L
  attempts <- 0L
  max_attempts <- 10L * bootstrap_rep
  while (accepted < bootstrap_rep && attempts < max_attempts) {
    latent <- c(0, cumsum(stats::rnorm(n_sub, sd = latent_sd)))
    observed <- latent + stats::rnorm(n_sub + 1L, sd = q_hat)
    stat <- .noisy_lm_statistic(
      observed, delta_sub, seconds_per_day, dependence_lag, spot_tuning,
      spot_window, spot_K, spot_window_scale, spot_K_scale
    )
    attempts <- attempts + 1L
    if (is.finite(stat$statistic)) {
      accepted <- accepted + 1L
      result[accepted] <- stat$statistic
    }
  }
  result[seq_len(accepted)]
}

#' Noise-robust Lee--Mykland local-average jump test
#'
#' Implements the local-average noisy LM statistic with LM (2012) Table 5
#' interpolation for the local block size and TSRSV spot-volatility
#' standardization.  Every bootstrap path regenerates the latent process and
#' noise and re-estimates both noise variance and spot volatility.
#'
#' @param prices A numeric price path, or a log-price path when
#'   `log_prices = FALSE`.
#' @param delta Sampling interval on the normalized time scale.
#' @param log_prices Should raw prices be log transformed first?
#' @param seconds_per_day Number of seconds in a regular trading day.
#' @param dependence_lag Noise-dependence lag.  The paper's i.i.d. setting uses
#'   one.
#' @param spot_tuning Either `"rate"` or `"plugin_eq56"`.
#' @param spot_window Optional TSRSV local-window override.
#' @param spot_K Optional TSRSV subsampling-lag override.
#' @param spot_window_scale Rate constant for the TSRSV local window.
#' @param spot_K_scale Rate constant for the TSRSV subsampling lag.
#' @param bootstrap_rep Number of fully recursive parametric bootstrap paths.
#' @param critical_value_c Uses the finite-sample multiplier
#'   `1 + critical_value_c * sqrt(delta)` on bootstrap statistics.  Set to zero
#'   for an uncorrected bootstrap.
#' @param seed Optional bootstrap seed.
#' @param alpha Optional decision level.
#'
#' @return An object of class `jump_test`.
#' @export
lm_noise_test <- function(prices,
                          delta = NULL,
                          log_prices = TRUE,
                          seconds_per_day = 23400,
                          dependence_lag = 1L,
                          spot_tuning = c("rate", "plugin_eq56"),
                          spot_window = NULL,
                          spot_K = NULL,
                          spot_window_scale = 1 / 4,
                          spot_K_scale = 1,
                          bootstrap_rep = 199L,
                          critical_value_c = 4,
                          seed = NULL,
                          alpha = NULL) {
  spot_tuning <- match.arg(spot_tuning)
  y <- .as_numeric_path(prices, log_prices = log_prices)
  delta <- .resolve_delta(length(y) - 1L, delta)
  dependence_lag <- as.integer(dependence_lag)
  if (dependence_lag < 1L) stop("`dependence_lag` must be at least one.", call. = FALSE)
  y_sub <- .subsample_prices(y, dependence_lag)
  delta_sub <- delta * dependence_lag
  observed <- .noisy_lm_statistic(
    y_sub, delta_sub, seconds_per_day, dependence_lag = 1L,
    spot_tuning = spot_tuning, spot_window = spot_window, spot_K = spot_K,
    spot_window_scale = spot_window_scale, spot_K_scale = spot_K_scale
  )
  if (!is.finite(observed$statistic)) {
    stop("Noisy LM statistic could not be computed from this path.", call. = FALSE)
  }
  boot <- .noisy_lm_bootstrap(
    y_sub, delta_sub, observed$q_hat, observed$sigma2_boot_path,
    bootstrap_rep, seconds_per_day, dependence_lag = 1L, spot_tuning,
    spot_window, spot_K, spot_window_scale, spot_K_scale, seed
  )
  boot <- boot[is.finite(boot)]
  if (length(boot) < 20L) {
    stop("Fewer than 20 valid recursive bootstrap paths were produced.", call. = FALSE)
  }
  critical_factor <- 1 + critical_value_c * sqrt(delta_sub)
  calibrated_boot <- critical_factor * boot
  p_value <- .clip_p((1 + sum(calibrated_boot >= observed$statistic)) /
                        (length(calibrated_boot) + 1))
  critical <- as.numeric(stats::quantile(
    calibrated_boot, probs = c(0.90, 0.95), names = FALSE, type = 8
  ))
  .new_jump_test(
    "LM", observed$statistic, p_value, alpha = alpha,
    max_standardized_return = observed$max_standardized_return,
    max_index = observed$max_index, q_hat = observed$q_hat,
    q_percent = observed$q_percent, dependence_lag = dependence_lag,
    block_size = observed$block_size, C = observed$C,
    sigma2_hat = observed$sigma2_hat, sigma2_min = observed$sigma2_min,
    sigma2_max = observed$sigma2_max, omega2_hat = observed$omega2_hat,
    spot_window = observed$spot_window, spot_K = observed$spot_K,
    spot_tuning = observed$spot_tuning, n_tests = observed$n_tests,
    n_tests_theory = observed$n_tests_theory, vhat2 = observed$vhat2,
    critical_10 = critical[1L], critical_05 = critical[2L],
    bootstrap_rep = length(boot), critical_value_factor = critical_factor,
    delta = delta, delta_sub = delta_sub
  )
}

#' Adaptive jump test under additive microstructure noise
#'
#' Applies AJJ, noisy LM, and their equal-weight Cauchy combination to one
#' observed price path.
#'
#' @inheritParams lm_noise_test
#' @param kn Pre-averaging window for `ajj_test()`.
#' @param rn AJJ self-normalizer block length.
#'
#' @return A list with components `ajj`, `lm`, and `cc`.
#' @export
adaptive_jump_test_noise <- function(prices,
                                     delta = NULL,
                                     kn = NULL,
                                     rn = NULL,
                                     log_prices = TRUE,
                                     seconds_per_day = 23400,
                                     dependence_lag = 1L,
                                     spot_tuning = c("rate", "plugin_eq56"),
                                     spot_window = NULL,
                                     spot_K = NULL,
                                     spot_window_scale = 1 / 4,
                                     spot_K_scale = 1,
                                     bootstrap_rep = 199L,
                                     critical_value_c = 4,
                                     seed = NULL,
                                     alpha = NULL) {
  spot_tuning <- match.arg(spot_tuning)
  ajj <- ajj_test(
    prices, delta = delta, kn = kn, rn = rn, log_prices = log_prices, alpha = alpha
  )
  lm <- lm_noise_test(
    prices, delta = delta, log_prices = log_prices,
    seconds_per_day = seconds_per_day, dependence_lag = dependence_lag,
    spot_tuning = spot_tuning, spot_window = spot_window, spot_K = spot_K,
    spot_window_scale = spot_window_scale, spot_K_scale = spot_K_scale,
    bootstrap_rep = bootstrap_rep, critical_value_c = critical_value_c,
    seed = seed, alpha = alpha
  )
  cc <- cauchy_combine(c(ajj$p_value, lm$p_value), alpha = alpha)
  list(ajj = ajj, lm = lm, cc = cc)
}
