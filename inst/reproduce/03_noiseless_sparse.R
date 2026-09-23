#!/usr/bin/env Rscript
# Reproduces the frictionless sparse finite-activity alternatives in the paper.

library(adaptivejump)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1L]] else "noiseless_sparse_results"
n_rep <- if (length(args) >= 2L) as.integer(args[[2L]]) else 1000L
workers <- if (length(args) >= 3L) as.integer(args[[3L]]) else max(1L, parallel::detectCores() - 1L)
bootstrap_rep <- if (length(args) >= 4L) as.integer(args[[4L]]) else 199L
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

lambda_values <- c(1, 1.5, 2, 2.5)
frequencies <- c(1L, 5L, 10L, 15L, 30L)
all_rates <- list()
index <- 0L

for (seconds in frequencies) {
  k_values <- if (seconds %in% c(1L, 5L)) c(2L, 3L, 4L) else 2L
  for (lambda in lambda_values) {
    for (k in k_values) {
      index <- index + 1L
      message(sprintf("%d-second, lambda=%g, k=%d", seconds, lambda, k))
      result <- run_jump_simulation(
        model = "sparse", test = "frictionless", n_rep = n_rep, workers = workers,
        seed = 20261000L + index,
        simulation_args = list(days = 1, sampling_seconds = seconds, lambda = lambda,
                               mark_sd = sqrt(0.05)),
        test_args = list(k = k, lm_calibration = "bootstrap", bootstrap_rep = bootstrap_rep),
        alpha = 0.05
      )
      rates <- result$rejection_rates
      rates$sampling_seconds <- seconds
      rates$lambda <- lambda
      rates$k <- k
      all_rates[[index]] <- rates
    }
  }
}

summary <- do.call(rbind, all_rates)
utils::write.csv(summary, file.path(output_dir, "noiseless_sparse_power_summary.csv"), row.names = FALSE)
print(summary)
