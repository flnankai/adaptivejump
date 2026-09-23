#!/usr/bin/env Rscript
# Runs the paper's three-second noise-robust empirical analysis.
# Download the archive manually from:
# https://www.kaggle.com/datasets/brtnsmth/intraday-market-data

library(adaptivejump)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("Usage: 07_empirical_3s.R <archive.zip> [output_dir] [workers] [bootstrap_rep]")
}
zipfile <- args[[1L]]
output_dir <- if (length(args) >= 2L) args[[2L]] else "empirical_3s_results"
workers <- if (length(args) >= 3L) as.integer(args[[3L]]) else 30L
bootstrap_rep <- if (length(args) >= 4L) as.integer(args[[4L]]) else 199L
tickers <- c("SPY", "QQQ", "IWM", "AAPL", "MSFT", "NVDA")
clean_dir <- file.path(output_dir, "clean_data")
result_dir <- file.path(output_dir, "daily_tests")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

expected <- file.path(clean_dir, paste0(tolower(tickers), "_clean_rth_3s.csv"))
if (!all(file.exists(expected))) {
  extraction <- read_tos_archive(zipfile, tickers = tickers, output_dir = clean_dir)
  utils::write.csv(extraction$summary, file.path(output_dir, "cleaning_summary.csv"), row.names = FALSE)
}

all_rejection <- list()
all_pvalues <- list()
for (i in seq_along(tickers)) {
  ticker <- tickers[[i]]
  message(sprintf("[%d/%d] %s", i, length(tickers), ticker))
  path <- file.path(clean_dir, paste0(tolower(ticker), "_clean_rth_3s.csv"))
  data <- utils::read.csv(path, stringsAsFactors = FALSE)
  result <- daily_jump_tests(
    data, method = "noise", workers = workers, clean = FALSE,
    test_args = list(
      kn = 100L, rn = 1000L, seconds_per_day = 23400,
      bootstrap_rep = bootstrap_rep, critical_value_c = 4,
      spot_tuning = "rate", spot_window_scale = 1 / 4, spot_K_scale = 1
    )
  )
  result$p_values$ticker <- ticker
  result$rejection_rates$ticker <- ticker
  utils::write.csv(result$p_values, file.path(result_dir, paste0(tolower(ticker), "_pvalues.csv")), row.names = FALSE)
  utils::write.csv(result$rejection_rates, file.path(result_dir, paste0(tolower(ticker), "_rejection_rates.csv")), row.names = FALSE)
  selection <- bh_select(result$p_values$lm_p[is.finite(result$p_values$lm_p)], fdr = 0.2)
  selected_dates <- result$p_values$date[is.finite(result$p_values$lm_p)][selection$selected]
  utils::write.csv(data.frame(date = selected_dates),
                   file.path(result_dir, paste0(tolower(ticker), "_lm_bh_fdr20.csv")), row.names = FALSE)
  all_pvalues[[i]] <- result$p_values
  all_rejection[[i]] <- result$rejection_rates
}

pvalues <- do.call(rbind, all_pvalues)
rejection <- do.call(rbind, all_rejection)
utils::write.csv(pvalues, file.path(output_dir, "alltime_3s_daily_pvalues.csv"), row.names = FALSE)
utils::write.csv(rejection, file.path(output_dir, "alltime_3s_rejection_rates.csv"), row.names = FALSE)

pdf(file.path(output_dir, "alltime_3s_qhat_hist_panels_3x2.pdf"), width = 12, height = 9)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (ticker in tickers) {
  x <- pvalues$q_hat[pvalues$ticker == ticker & is.finite(pvalues$q_hat)]
  hist(x, breaks = 30, main = ticker, xlab = expression(hat(q)), col = "grey80", border = "white")
}
dev.off()

pdf(file.path(output_dir, "alltime_3s_daily_pvalue_panels_6x1.pdf"), width = 12, height = 15)
par(mfrow = c(6, 1), mar = c(3, 4, 2, 1))
for (ticker in tickers) {
  panel <- pvalues[pvalues$ticker == ticker, ]
  panel <- panel[order(panel$date), ]
  plot(seq_len(nrow(panel)), panel$aj_p, type = "l", ylim = c(0, 1),
       xlab = "Trading day", ylab = "p-value", main = ticker, col = "#1b9e77")
  lines(panel$lm_p, col = "#d95f02")
  lines(panel$cc_p, col = "#7570b3")
  abline(h = 0.05, lty = 2, col = "grey50")
  legend("topright", legend = c("AJJ", "LM", "CC"),
         col = c("#1b9e77", "#d95f02", "#7570b3"), lty = 1, bty = "n")
}
dev.off()

print(rejection[rejection$alpha == 0.05, ])
