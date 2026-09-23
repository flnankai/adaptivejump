# This late-loaded definition keeps the Gumbel critical values on the
# standard-Gumbel scale rather than the logistic scale.
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
  if (is.null(K)) K <- .lm_window_rule(n_returns, window_exponent)
  K <- as.integer(K)
  stat <- .lm_statistic(returns, K = K, n_gumbel = n_returns)

  if (calibration == "gumbel") {
    p_value <- .clip_p(1 - exp(-exp(-stat$statistic)))
    critical_values <- -log(-log(c(0.90, 0.95)))
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
