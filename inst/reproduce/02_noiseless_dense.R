#!/usr/bin/env Rscript
# Reproduces the frictionless dense-local alternatives in the paper.

library(adaptivejump)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1L]] else "noiseless_dense_results"
n_rep <- if (length(args) >= 2L) as.integer(args[[2L]]) else 1000L
workers <- if (length(args) >= 3L) as.integer(args[[3L]]) else max(1L, parallel::detectCores() - 1L)
bootstrap_rep <- if (length(args) >= 4L) as.integer(args[[4L]]) else 199L
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Mark variances are the values used in the paper's 1--30 second designs.
designs <- data.frame(
  sampling_seconds = c(1L, 5L, 10L, 15L, 30L),
  mark_variance = c(0.5, 12.5, 50, 112.5, 450),
  stringsAsFactors = FALSE
)
theta_values <- c(100, 200, 300, 400, 500)
all_rates <- list()
index <- 0L

for (d in seq_len(nrow(designs))) {
  seconds <- designs$sampling_seconds[d]
  k_values <- if (seconds %in% c(1L, 5L)) c(2L, 3L, 4L) else 2L
  for (theta in theta_values) {
    for (k in k_values) {
      index <- index + 1L
      message(sprintf("%d-second, theta=%g, k=%d", seconds, theta, k))
      result <- run_jump_simulation(
        model = "dense", test = "frictionless", n_rep = n_rep, workers = workers,
        seed = 20260950L + index,
        simulation_args = list(
          days = 1, sampling_seconds = seconds, theta = theta,
          mark_sd = sqrt(designs$mark_variance[d])
        ),
        test_args = list(k = k, lm_calibration = "bootstrap", bootstrap_rep = bootstrap_rep),
        alpha = 0.05
      )
      rates <- result$rejection_rates
      rates$sampling_seconds <- seconds
      rates$mark_variance <- designs$mark_variance[d]
      rates$theta <- theta
      rates$k <- k
      all_rates[[index]] <- rates
    }
  }
}

summary <- do.call(rbind, all_rates)
utils::write.csv(summary, file.path(output_dir, "noiseless_dense_power_summary.csv"), row.names = FALSE)

pdf(file.path(output_dir, "noiseless_dense_power.pdf"), width = 12, height = 9)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (seconds in designs$sampling_seconds) {
  panel <- summary[summary$sampling_seconds == seconds & summary$k == 2L, ]
  plot(range(theta_values), c(0, 1), type = "n", xlab = expression(theta),
       ylab = "Power (5% level)", main = sprintf("%d-second", seconds))
  abline(h = 0.05, lty = 2, col = "grey70")
  for (method in c("AJ", "LM", "CC")) {
    values <- panel[panel$method == method, ]
    values <- values[order(values$theta), ]
    lines(values$theta, values$rejection_rate,
          col = c(AJ = "#1b9e77", LM = "#d95f02", CC = "#7570b3")[[method]],
          type = "o", pch = c(AJ = 16, LM = 17, CC = 15)[[method]], lwd = 2)
  }
}
plot.new()
legend("center", legend = c("AJ", "LM", "CC"),
       col = c("#1b9e77", "#d95f02", "#7570b3"), pch = c(16, 17, 15), lwd = 2, bty = "n")
dev.off()
print(summary)
