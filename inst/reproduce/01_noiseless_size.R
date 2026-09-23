#!/usr/bin/env Rscript
# Reproduces the frictionless null-size design in the paper.

library(adaptivejump)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1L]] else "noiseless_size_results"
n_rep <- if (length(args) >= 2L) as.integer(args[[2L]]) else 5000L
workers <- if (length(args) >= 3L) as.integer(args[[3L]]) else max(1L, parallel::detectCores() - 1L)
bootstrap_rep <- if (length(args) >= 4L) as.integer(args[[4L]]) else 199L
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

grid <- do.call(rbind, lapply(c(1L, 5L, 10L, 15L, 30L), function(seconds) {
  k_values <- if (seconds %in% c(1L, 5L)) c(2L, 3L, 4L) else 2L
  data.frame(sampling_seconds = seconds, k = k_values)
}))
all_rates <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  design <- grid[i, ]
  message(sprintf("[%d/%d] %d-second data, k=%d", i, nrow(grid),
                  design$sampling_seconds, design$k))
  result <- run_jump_simulation(
    model = "null", test = "frictionless", n_rep = n_rep, workers = workers,
    seed = 20260901L + i,
    simulation_args = list(days = 1, sampling_seconds = design$sampling_seconds),
    test_args = list(k = design$k, lm_calibration = "bootstrap",
                     bootstrap_rep = bootstrap_rep)
  )
  result$p_values$sampling_seconds <- design$sampling_seconds
  result$p_values$k <- design$k
  utils::write.csv(result$p_values,
                   file.path(output_dir, sprintf("pvalues_%ds_k%d.csv", design$sampling_seconds, design$k)),
                   row.names = FALSE)
  rates <- result$rejection_rates
  rates$sampling_seconds <- design$sampling_seconds
  rates$k <- design$k
  all_rates[[i]] <- rates
}

summary <- do.call(rbind, all_rates)
summary <- summary[order(summary$sampling_seconds, summary$k, summary$alpha, summary$method), ]
utils::write.csv(summary, file.path(output_dir, "noiseless_size_summary.csv"), row.names = FALSE)
print(summary)
