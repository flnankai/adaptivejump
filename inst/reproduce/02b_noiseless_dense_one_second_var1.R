#!/usr/bin/env Rscript
# Reproduces the second one-second dense-local panel with Y ~ N(0, 1).

library(adaptivejump)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1L]] else "noiseless_dense_var1_results"
n_rep <- if (length(args) >= 2L) as.integer(args[[2L]]) else 1000L
workers <- if (length(args) >= 3L) as.integer(args[[3L]]) else max(1L, parallel::detectCores() - 1L)
bootstrap_rep <- if (length(args) >= 4L) as.integer(args[[4L]]) else 199L
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theta_values <- c(100, 200, 300, 400, 500)
all_rates <- list()
index <- 0L
for (theta in theta_values) {
  for (k in c(2L, 3L, 4L)) {
    index <- index + 1L
    result <- run_jump_simulation(
      model = "dense", test = "frictionless", n_rep = n_rep, workers = workers,
      seed = 20260980L + index,
      simulation_args = list(days = 1, sampling_seconds = 1L, theta = theta, mark_sd = 1),
      test_args = list(k = k, lm_calibration = "bootstrap", bootstrap_rep = bootstrap_rep),
      alpha = 0.05
    )
    rates <- result$rejection_rates
    rates$theta <- theta
    rates$k <- k
    rates$mark_variance <- 1
    all_rates[[index]] <- rates
  }
}
summary <- do.call(rbind, all_rates)
utils::write.csv(summary, file.path(output_dir, "noiseless_dense_1s_var1_summary.csv"), row.names = FALSE)
print(summary)
