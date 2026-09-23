#!/usr/bin/env Rscript
# Reproduces the noisy null-size scenarios with the recursive LM bootstrap.

library(adaptivejump)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1L]] else "noisy_size_results"
n_rep <- if (length(args) >= 2L) as.integer(args[[2L]]) else 5000L
workers <- if (length(args) >= 3L) as.integer(args[[3L]]) else max(1L, parallel::detectCores() - 1L)
bootstrap_rep <- if (length(args) >= 4L) as.integer(args[[4L]]) else 199L
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- expand.grid(
  sampling_seconds = c(5L, 10L, 15L),
  noise_sd = c(0.005, 0.01),
  noise_distribution = c("normal", "t8"),
  stringsAsFactors = FALSE
)
all_rates <- vector("list", nrow(scenarios))

for (i in seq_len(nrow(scenarios))) {
  scenario <- scenarios[i, ]
  message(sprintf("[%d/%d] %ds, q=%.3f, %s", i, nrow(scenarios),
                  scenario$sampling_seconds, scenario$noise_sd, scenario$noise_distribution))
  result <- run_jump_simulation(
    model = "null", test = "noise", n_rep = n_rep, workers = workers,
    seed = 20261100L + i,
    simulation_args = list(
      days = 5, sampling_seconds = scenario$sampling_seconds,
      noise_sd = scenario$noise_sd, noise_distribution = scenario$noise_distribution
    ),
    test_args = list(
      kn = 100L, rn = 1000L, seconds_per_day = 23400,
      bootstrap_rep = bootstrap_rep, critical_value_c = 4,
      spot_tuning = "rate", spot_window_scale = 1 / 4, spot_K_scale = 1
    )
  )
  result$p_values$sampling_seconds <- scenario$sampling_seconds
  result$p_values$noise_sd <- scenario$noise_sd
  result$p_values$noise_distribution <- scenario$noise_distribution
  utils::write.csv(result$p_values,
                   file.path(output_dir, sprintf("pvalues_%ds_%s_q%s.csv",
                                                  scenario$sampling_seconds,
                                                  scenario$noise_distribution,
                                                  gsub("[.]", "p", scenario$noise_sd))),
                   row.names = FALSE)
  rates <- result$rejection_rates
  rates$sampling_seconds <- scenario$sampling_seconds
  rates$noise_sd <- scenario$noise_sd
  rates$noise_distribution <- scenario$noise_distribution
  all_rates[[i]] <- rates
}

summary <- do.call(rbind, all_rates)
utils::write.csv(summary, file.path(output_dir, "noisy_size_summary.csv"), row.names = FALSE)
print(summary)
