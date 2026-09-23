#!/usr/bin/env Rscript
# Runs the paper's 1-, 2-, and 5-minute frictionless empirical analysis.

library(adaptivejump)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("Usage: 08_empirical_coarse.R <clean_3s_directory> [output_dir] [workers] [bootstrap_rep]")
}
clean_dir <- args[[1L]]
output_dir <- if (length(args) >= 2L) args[[2L]] else "empirical_coarse_results"
workers <- if (length(args) >= 3L) as.integer(args[[3L]]) else 30L
bootstrap_rep <- if (length(args) >= 4L) as.integer(args[[4L]]) else 199L
tickers <- c("SPY", "QQQ", "IWM", "AAPL", "MSFT", "NVDA")
frequencies <- c(`1min` = 60L, `2min` = 120L, `5min` = 300L)
data_dir <- file.path(output_dir, "resampled_data")
result_dir <- file.path(output_dir, "daily_tests")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

all_rates <- list()
index <- 0L
for (ticker in tickers) {
  raw_path <- file.path(clean_dir, paste0(tolower(ticker), "_clean_rth_3s.csv"))
  if (!file.exists(raw_path)) stop("Missing clean 3-second file: ", raw_path)
  raw <- utils::read.csv(raw_path, stringsAsFactors = FALSE)
  for (label in names(frequencies)) {
    index <- index + 1L
    message(sprintf("%s, %s", ticker, label))
    coarse <- resample_last_tick(raw, frequencies[[label]])
    utils::write.csv(coarse, file.path(data_dir, sprintf("%s_%s.csv", tolower(ticker), label)), row.names = FALSE)
    result <- daily_jump_tests(
      coarse, method = "frictionless", workers = workers, clean = FALSE,
      test_args = list(k = 2L, lm_calibration = "bootstrap", bootstrap_rep = bootstrap_rep)
    )
    result$p_values$ticker <- ticker
    result$p_values$frequency <- label
    result$rejection_rates$ticker <- ticker
    result$rejection_rates$frequency <- label
    utils::write.csv(result$p_values,
                     file.path(result_dir, sprintf("%s_%s_pvalues.csv", tolower(ticker), label)),
                     row.names = FALSE)
    all_rates[[index]] <- result$rejection_rates
  }
}

rates <- do.call(rbind, all_rates)
utils::write.csv(rates, file.path(output_dir, "coarse_rejection_rates.csv"), row.names = FALSE)

rate05 <- rates[rates$alpha == 0.05, ]
wide <- data.frame(ticker = tickers, stringsAsFactors = FALSE)
for (label in names(frequencies)) {
  for (method in c("AJ", "LM", "CC")) {
    key <- paste(label, method, sep = "_")
    values <- rate05[rate05$frequency == label & rate05$method == method,
                     c("ticker", "rejection_rate")]
    wide[[key]] <- values$rejection_rate[match(tickers, values$ticker)]
  }
}
utils::write.csv(wide, file.path(output_dir, "coarse_rejection_rates_alpha05_wide.csv"), row.names = FALSE)

tex_path <- file.path(output_dir, "coarse_rejection_rates_alpha05.tex")
lines <- c(
  "\\begin{tabular}{lrrrrrrrrr}",
  "\\toprule",
  "Asset & 1min AJ & 1min LM & 1min CC & 2min AJ & 2min LM & 2min CC & 5min AJ & 5min LM & 5min CC \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(wide))) {
  lines <- c(lines, paste(c(wide$ticker[i], sprintf("%.4f", as.numeric(wide[i, -1L]))), collapse = " & "), "\\\\")
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, tex_path)
print(wide)
