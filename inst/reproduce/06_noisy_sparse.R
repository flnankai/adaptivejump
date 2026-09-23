#!/usr/bin/env Rscript
# Reproduces four-panel noisy sparse finite-activity power curves.

library(adaptivejump)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1L]] else "noisy_sparse_results"
n_rep <- if (length(args) >= 2L) as.integer(args[[2L]]) else 1000L
workers <- if (length(args) >= 3L) as.integer(args[[3L]]) else max(1L, parallel::detectCores() - 1L)
bootstrap_rep <- if (length(args) >= 4L) as.integer(args[[4L]]) else 199L
frequencies <- if (length(args) >= 5L) as.integer(strsplit(args[[5L]], ",")[[1L]]) else c(5L, 10L)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

lambda_values <- seq(0.5, 2.5, by = 0.5)
scenarios <- data.frame(
  noise_distribution = c("normal", "normal", "t8", "t8"),
  noise_sd = c(0.005, 0.01, 0.005, 0.01),
  label = c("Normal, q = 0.005", "Normal, q = 0.01",
            "t8, q = 0.005", "t8, q = 0.01"),
  stringsAsFactors = FALSE
)
all_rates <- list()
index <- 0L

for (seconds in frequencies) {
  for (s in seq_len(nrow(scenarios))) {
    for (lambda in lambda_values) {
      index <- index + 1L
      scenario <- scenarios[s, ]
      message(sprintf("%ds, %s, lambda=%.1f", seconds, scenario$label, lambda))
      result <- run_jump_simulation(
        model = "sparse", test = "noise", n_rep = n_rep, workers = workers,
        seed = 20261300L + index,
        simulation_args = list(
          days = 5, sampling_seconds = seconds, lambda = lambda,
          mark_sd = sqrt(0.05), noise_sd = scenario$noise_sd,
          noise_distribution = scenario$noise_distribution
        ),
        test_args = list(kn = 100L, rn = 1000L, seconds_per_day = 23400,
                         bootstrap_rep = bootstrap_rep, critical_value_c = 4,
                         spot_tuning = "rate", spot_window_scale = 1 / 4,
                         spot_K_scale = 1),
        alpha = 0.05
      )
      rates <- result$rejection_rates
      rates$sampling_seconds <- seconds
      rates$lambda <- lambda
      rates$noise_sd <- scenario$noise_sd
      rates$noise_distribution <- scenario$noise_distribution
      rates$label <- scenario$label
      all_rates[[index]] <- rates
    }
  }
}

summary <- do.call(rbind, all_rates)
utils::write.csv(summary, file.path(output_dir, "noisy_sparse_power_summary.csv"), row.names = FALSE)

for (seconds in frequencies) {
  pdf(file.path(output_dir, sprintf("noisy_sparse_power_%ds.pdf", seconds)), width = 12, height = 9)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (s in seq_len(nrow(scenarios))) {
    panel <- summary[summary$sampling_seconds == seconds &
                       summary$noise_distribution == scenarios$noise_distribution[s] &
                       summary$noise_sd == scenarios$noise_sd[s], ]
    plot(range(lambda_values), c(0, 1), type = "n", xlab = expression(lambda),
         ylab = "Power (5% level)", main = scenarios$label[s])
    abline(h = 0.05, lty = 2, col = "grey70")
    for (method in c("AJJ", "LM", "CC")) {
      values <- panel[panel$method == method, ]
      values <- values[order(values$lambda), ]
      lines(values$lambda, values$rejection_rate,
            col = c(AJJ = "#1b9e77", LM = "#d95f02", CC = "#7570b3")[[method]],
            type = "o", pch = c(AJJ = 16, LM = 17, CC = 15)[[method]], lwd = 2)
    }
    legend("topright", legend = c("AJJ", "LM", "CC"),
           col = c("#1b9e77", "#d95f02", "#7570b3"), pch = c(16, 17, 15),
           lwd = 2, bty = "n")
  }
  dev.off()
}
print(summary)
