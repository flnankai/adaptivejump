.draw_standardized_noise <- function(n, distribution = c("normal", "t8")) {
  distribution <- match.arg(distribution)
  if (distribution == "normal") return(stats::rnorm(n))
  # A t_8 variate has variance 4/3, so this has unit variance.
  stats::rt(n, df = 8) * sqrt(3 / 4)
}

.draw_marks <- function(n, sd, distribution = c("normal", "binary")) {
  distribution <- match.arg(distribution)
  if (n == 0L) return(numeric())
  if (distribution == "normal") return(stats::rnorm(n, sd = sd))
  sample(c(-sd, sd), size = n, replace = TRUE)
}

#' Simulate a stochastic-volatility jump path
#'
#' Simulates the Heston-style log-price design used throughout the paper:
#' \deqn{dX_t = -V_t dt/2 + sqrt(V_t)dW_t + dJ_t,}
#' with a CIR-type variance process and leverage.  The dense alternative uses
#' a per-step Poisson intensity `theta * sqrt(delta)` and jump size
#' `sqrt(V) * sqrt(delta) * Y`; the sparse alternative uses a compound Poisson
#' process with daily intensity `lambda`.
#'
#' @param model One of `"null"`, `"dense"`, or `"sparse"`.
#' @param days Number of trading days in the path.
#' @param seconds_per_day Length of a trading day in seconds.
#' @param sampling_seconds Sampling interval in seconds.
#' @param x0 Initial log price.
#' @param v0 Initial variance.
#' @param beta_bar Long-run variance.
#' @param kappa Variance mean-reversion coefficient.
#' @param gamma_vol Volatility of variance.
#' @param rho Correlation between price and variance Brownian motions.
#' @param theta Dense-alternative intensity coefficient.
#' @param lambda Sparse-alternative daily jump intensity.
#' @param mark_sd Standard deviation of jump marks.
#' @param mark_distribution Either `"normal"` or symmetric `"binary"`.
#' @param noise_sd Additive observation-noise standard deviation.  Set to zero
#'   for the frictionless path.
#' @param noise_distribution Either `"normal"` or variance-standardized `"t8"`.
#' @param noise_model Either homoskedastic `"iid"` or `"ajl_gaussian"`, which
#'   scales noise with instantaneous volatility.
#' @param seed Optional seed for a single direct simulation.
#'
#' @return A list containing latent and observed log prices, variance, jump
#'   increments, and the normalized sampling interval.
#' @export
simulate_jump_path <- function(model = c("null", "dense", "sparse"),
                               days = 1,
                               seconds_per_day = 23400,
                               sampling_seconds = 1,
                               x0 = log(25),
                               v0 = 0.16,
                               beta_bar = 0.16,
                               kappa = 5,
                               gamma_vol = 0.5,
                               rho = -0.5,
                               theta = 1,
                               lambda = 1,
                               mark_sd = 1,
                               mark_distribution = c("normal", "binary"),
                               noise_sd = 0,
                               noise_distribution = c("normal", "t8"),
                               noise_model = c("iid", "ajl_gaussian"),
                               seed = NULL) {
  model <- match.arg(model)
  mark_distribution <- match.arg(mark_distribution)
  noise_distribution <- match.arg(noise_distribution)
  noise_model <- match.arg(noise_model)
  if (!is.null(seed)) set.seed(as.integer(abs(seed) %% (.Machine$integer.max - 1L)))
  if (days <= 0 || seconds_per_day <= 0 || sampling_seconds <= 0 ||
      seconds_per_day %% sampling_seconds != 0) {
    stop("`days`, `seconds_per_day`, and `sampling_seconds` must define a regular grid.",
         call. = FALSE)
  }
  n <- as.integer(days * seconds_per_day / sampling_seconds)
  if (n < 12L) stop("The simulation grid is too short.", call. = FALSE)
  delta <- days / n
  x <- numeric(n + 1L)
  v <- numeric(n + 1L)
  jump_returns <- numeric(n)
  x[1L] <- x0
  v[1L] <- v0
  z1 <- stats::rnorm(n)
  z2 <- stats::rnorm(n)
  dW <- sqrt(delta) * z1
  dB <- sqrt(delta) * (rho * z1 + sqrt(1 - rho^2) * z2)

  for (i in seq_len(n)) {
    v_now <- max(v[i], 1e-10)
    sigma_now <- sqrt(v_now)
    continuous <- -0.5 * v_now * delta + sigma_now * dW[i]
    jump <- 0
    if (model == "dense") {
      count <- stats::rpois(1L, lambda = theta * sqrt(delta))
      jump <- sigma_now * sqrt(delta) * sum(.draw_marks(count, mark_sd, mark_distribution))
    } else if (model == "sparse") {
      count <- stats::rpois(1L, lambda = lambda * delta)
      jump <- sum(.draw_marks(count, mark_sd, mark_distribution))
    }
    jump_returns[i] <- jump
    x[i + 1L] <- x[i] + continuous + jump
    v[i + 1L] <- max(
      v[i] + kappa * (beta_bar - v_now) * delta + gamma_vol * sigma_now * dB[i],
      1e-10
    )
  }

  if (noise_sd == 0) {
    noise <- rep(0, n + 1L)
  } else if (noise_model == "iid") {
    noise <- noise_sd * .draw_standardized_noise(n + 1L, noise_distribution)
  } else {
    noise <- noise_sd * sqrt(v / beta_bar) *
      .draw_standardized_noise(n + 1L, noise_distribution)
  }
  list(
    latent_log_price = x,
    observed_log_price = x + noise,
    variance = v,
    jump_returns = jump_returns,
    noise = noise,
    n_returns = n,
    delta = delta,
    days = days,
    sampling_seconds = sampling_seconds,
    model = model
  )
}

.simulation_one_rep <- function(i, model, test, simulation_args, test_args) {
  sim <- do.call(simulate_jump_path, c(list(model = model), simulation_args))
  if (test == "frictionless") {
    result <- do.call(
      adaptive_jump_test,
      c(list(prices = exp(sim$latent_log_price), log_prices = TRUE), test_args)
    )
    c(AJ = result$aj$p_value, LM = result$lm$p_value, CC = result$cc$p_value)
  } else {
    result <- do.call(
      adaptive_jump_test_noise,
      c(list(prices = exp(sim$observed_log_price), log_prices = TRUE), test_args)
    )
    c(AJJ = result$ajj$p_value, LM = result$lm$p_value, CC = result$cc$p_value)
  }
}

#' Run a Monte Carlo jump-test experiment
#'
#' Runs a null, dense, or sparse simulation design and returns raw p-values and
#' rejection frequencies for the three procedures.  The function never splits
#' bootstrap replications into batches; each test call receives the full number
#' of requested bootstrap paths.
#'
#' @param model Simulation alternative passed to `simulate_jump_path()`.
#' @param test `"frictionless"` for AJ/LM/CC or `"noise"` for AJJ/noisy-LM/CC.
#' @param n_rep Number of Monte Carlo replications.
#' @param simulation_args Named list forwarded to `simulate_jump_path()`.
#' @param test_args Named list forwarded to the selected adaptive test.
#' @param alpha Vector of rejection levels.
#' @param workers Number of parallel workers.  Use one for serial execution.
#' @param seed Optional master seed.
#'
#' @return A list with a p-value data frame, rejection-rate data frame, and the
#'   settings used.
#' @export
run_jump_simulation <- function(model = c("null", "dense", "sparse"),
                                test = c("frictionless", "noise"),
                                n_rep = 1000L,
                                simulation_args = list(),
                                test_args = list(),
                                alpha = c(0.10, 0.05),
                                workers = 1L,
                                seed = NULL) {
  model <- match.arg(model)
  test <- match.arg(test)
  n_rep <- as.integer(n_rep)
  workers <- as.integer(workers)
  if (n_rep < 1L || workers < 1L) {
    stop("`n_rep` and `workers` must be positive integers.", call. = FALSE)
  }
  worker <- function(i) .simulation_one_rep(i, model, test, simulation_args, test_args)
  if (workers == 1L) {
    if (!is.null(seed)) set.seed(as.integer(abs(seed) %% (.Machine$integer.max - 1L)))
    values <- lapply(seq_len(n_rep), worker)
  } else {
    workers <- min(workers, n_rep)
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(adaptivejump))
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = seed)
    values <- parallel::parLapplyLB(cl, seq_len(n_rep), worker)
  }
  p_values <- as.data.frame(do.call(rbind, values))
  p_values$replication <- seq_len(n_rep)
  p_values <- p_values[c("replication", setdiff(names(p_values), "replication"))]
  methods <- setdiff(names(p_values), "replication")
  rates <- do.call(rbind, lapply(alpha, function(a) {
    data.frame(
      alpha = a,
      method = methods,
      rejection_rate = vapply(methods, function(m) mean(p_values[[m]] <= a), numeric(1L)),
      stringsAsFactors = FALSE
    )
  }))
  list(
    p_values = p_values,
    rejection_rates = rates,
    settings = list(model = model, test = test, n_rep = n_rep,
                    simulation_args = simulation_args, test_args = test_args,
                    alpha = alpha, workers = workers, seed = seed)
  )
}
