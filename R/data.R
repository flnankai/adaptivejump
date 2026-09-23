.parse_timestamps <- function(x) {
  formats <- c(
    "%Y-%m-%d %H:%M:%OS", "%Y-%m-%d %H:%M:%S",
    "%m/%d/%y %H:%M:%OS", "%m/%d/%y %H:%M:%S",
    "%m/%d/%Y %H:%M:%OS", "%m/%d/%Y %H:%M:%S"
  )
  out <- as.POSIXct(rep(NA_real_, length(x)), origin = "1970-01-01", tz = "UTC")
  remaining <- rep(TRUE, length(x))
  for (format in formats) {
    if (!any(remaining)) break
    parsed <- as.POSIXct(x[remaining], format = format, tz = "UTC")
    hit <- !is.na(parsed)
    if (any(hit)) {
      index <- which(remaining)[hit]
      out[index] <- parsed[hit]
      remaining[index] <- FALSE
    }
  }
  out
}

.time_to_seconds <- function(time) {
  as.integer(substr(time, 1L, 2L)) * 3600L +
    as.integer(substr(time, 4L, 5L)) * 60L +
    as.integer(substr(time, 7L, 8L))
}

.is_weekday <- function(date) {
  as.POSIXlt(as.Date(date))$wday %in% 1:5
}

.classify_intraday_day <- function(day, min_obs = 111L) {
  price <- day$price[is.finite(day$price)]
  reasons <- character()
  if (any(!is.finite(day$price))) reasons <- c(reasons, "contains_NA")
  if (length(unique(price)) <= 1L) reasons <- c(reasons, "constant_price")
  if (nrow(day) < min_obs) reasons <- c(reasons, "too_few_obs")
  list(
    keep = !length(reasons),
    reason = paste(reasons, collapse = ";"),
    n_obs = nrow(day),
    n_unique_prices = length(unique(price)),
    start_time = if (nrow(day)) day$time[1L] else NA_character_,
    end_time = if (nrow(day)) day$time[nrow(day)] else NA_character_
  )
}

.append_csv <- function(data, path) {
  if (!nrow(data)) return(invisible(NULL))
  utils::write.table(
    data, path, sep = ",", row.names = FALSE,
    col.names = !file.exists(path), append = file.exists(path), qmethod = "double"
  )
  invisible(path)
}

#' Clean regular-session intraday observations
#'
#' Filters to weekdays and the U.S. regular cash session, sorts and removes
#' duplicate timestamps, and removes days with missing, constant, or too few
#' prices.  It supports either a timestamp column or separate date and time
#' columns.
#'
#' @param data A data frame.
#' @param price_col Price-column name.
#' @param timestamp_col Optional timestamp-column name.
#' @param date_col Date-column name when timestamps are separate.
#' @param time_col Time-column name when timestamps are separate.
#' @param session_start Session start time in `HH:MM:SS` form.
#' @param session_end Session end time in `HH:MM:SS` form.
#' @param min_obs Minimum number of observations retained in a day.
#'
#' @return A list containing clean `data` with columns `date`, `time`, `price`
#'   and an `excluded` data frame.
#' @export
clean_intraday_data <- function(data,
                                price_col = "price",
                                timestamp_col = NULL,
                                date_col = "date",
                                time_col = "time",
                                session_start = "09:30:00",
                                session_end = "16:00:00",
                                min_obs = 111L) {
  if (!is.data.frame(data) || !price_col %in% names(data)) {
    stop("`data` must be a data frame containing `price_col`.", call. = FALSE)
  }
  if (!is.null(timestamp_col)) {
    if (!timestamp_col %in% names(data)) {
      stop("`timestamp_col` is not a column of `data`.", call. = FALSE)
    }
    stamp <- .parse_timestamps(trimws(as.character(data[[timestamp_col]])))
    out <- data.frame(
      date = format(stamp, "%Y-%m-%d"), time = format(stamp, "%H:%M:%S"),
      price = suppressWarnings(as.numeric(data[[price_col]])),
      stringsAsFactors = FALSE
    )
  } else {
    if (!all(c(date_col, time_col) %in% names(data))) {
      stop("Separate date and time columns are required when `timestamp_col` is NULL.",
           call. = FALSE)
    }
    out <- data.frame(
      date = as.character(data[[date_col]]), time = as.character(data[[time_col]]),
      price = suppressWarnings(as.numeric(data[[price_col]])),
      stringsAsFactors = FALSE
    )
  }
  valid_time <- !is.na(out$date) & !is.na(out$time) &
    out$time >= session_start & out$time <= session_end & .is_weekday(out$date)
  out <- out[valid_time, , drop = FALSE]
  out <- out[order(out$date, out$time), , drop = FALSE]
  out <- out[!duplicated(paste(out$date, out$time)), , drop = FALSE]
  by_day <- split(out, out$date)
  accepted <- list()
  excluded <- list()
  ai <- 0L
  ei <- 0L
  for (date in names(by_day)) {
    day <- by_day[[date]]
    info <- .classify_intraday_day(day, min_obs = min_obs)
    if (info$keep) {
      ai <- ai + 1L
      accepted[[ai]] <- day
    } else {
      ei <- ei + 1L
      excluded[[ei]] <- data.frame(
        date = date, reason = info$reason, n_obs = info$n_obs,
        n_unique_prices = info$n_unique_prices, start_time = info$start_time,
        end_time = info$end_time, stringsAsFactors = FALSE
      )
    }
  }
  list(
    data = if (ai) do.call(rbind, accepted) else out[0, , drop = FALSE],
    excluded = if (ei) do.call(rbind, excluded) else data.frame(
      date = character(), reason = character(), n_obs = integer(),
      n_unique_prices = integer(), start_time = character(), end_time = character()
    )
  )
}

#' Resample intraday data using the last tick in each interval
#'
#' @param data A cleaned data frame with `date`, `time`, and `price` columns.
#' @param interval_seconds Positive resampling interval in seconds.
#' @param session_start Start of the regular session.
#'
#' @return A data frame with one last-price observation per nonempty interval.
#' @export
resample_last_tick <- function(data,
                               interval_seconds,
                               session_start = "09:30:00") {
  if (!all(c("date", "time", "price") %in% names(data))) {
    stop("`data` must contain date, time, and price columns.", call. = FALSE)
  }
  interval_seconds <- as.integer(interval_seconds)
  if (is.na(interval_seconds) || interval_seconds < 1L) {
    stop("`interval_seconds` must be a positive integer.", call. = FALSE)
  }
  data <- data[order(data$date, data$time), c("date", "time", "price"), drop = FALSE]
  start_second <- .time_to_seconds(session_start)
  pieces <- lapply(split(data, data$date), function(day) {
    seconds <- .time_to_seconds(day$time)
    bucket <- (seconds - start_second) %/% interval_seconds
    last <- tapply(seq_along(bucket), bucket, tail, n = 1L)
    day[as.integer(last), , drop = FALSE]
  })
  out <- do.call(rbind, pieces)
  out[order(out$date, out$time), , drop = FALSE]
}

.read_tos_entry <- function(zipfile, entry_name, tickers) {
  header <- tryCatch(
    names(utils::read.csv(unz(zipfile, entry_name), nrows = 0L, check.names = FALSE)),
    error = function(e) NULL
  )
  reader <- utils::read.csv
  if (is.null(header) || !all(c("TimeStamp", tickers) %in% header)) {
    header <- tryCatch(
      names(utils::read.delim(unz(zipfile, entry_name), nrows = 0L, check.names = FALSE)),
      error = function(e) NULL
    )
    reader <- utils::read.delim
  }
  if (is.null(header) || !all(c("TimeStamp", tickers) %in% header)) return(NULL)
  col_classes <- rep("NULL", length(header))
  names(col_classes) <- header
  col_classes["TimeStamp"] <- "character"
  col_classes[tickers] <- "numeric"
  raw <- tryCatch(
    reader(unz(zipfile, entry_name), stringsAsFactors = FALSE, check.names = FALSE,
           colClasses = col_classes),
    error = function(e) NULL
  )
  if (is.null(raw)) return(NULL)
  stamp <- .parse_timestamps(trimws(as.character(raw$TimeStamp)))
  keep <- !is.na(stamp) & .is_weekday(format(stamp, "%Y-%m-%d")) &
    format(stamp, "%H:%M:%S") >= "09:30:00" & format(stamp, "%H:%M:%S") <= "16:00:00"
  if (!any(keep)) return(NULL)
  out <- data.frame(date = format(stamp[keep], "%Y-%m-%d"),
                    time = format(stamp[keep], "%H:%M:%S"), stringsAsFactors = FALSE)
  for (ticker in tickers) out[[ticker]] <- raw[[ticker]][keep]
  out
}

#' Extract selected symbols from the public ThinkOrSwim archive
#'
#' The archive is available from Kaggle at
#' <https://www.kaggle.com/datasets/brtnsmth/intraday-market-data>.  This
#' function reads a local archive path and never downloads or bundles data.
#'
#' @param zipfile Local path to the downloaded archive.
#' @param tickers Character vector of columns to extract.
#' @param output_dir Optional directory to receive one cleaned CSV and one
#'   exclusion CSV per ticker.  Supplying it avoids holding the full sample in
#'   memory.
#' @param overwrite Whether existing output files may be replaced.
#' @param min_obs Minimum observations for a valid day.
#'
#' @return A list with output paths and a per-ticker summary.  If `output_dir`
#'   is `NULL`, the list also contains the extracted data in memory.
#' @export
read_tos_archive <- function(zipfile,
                             tickers = c("SPY", "QQQ", "IWM", "AAPL", "MSFT", "NVDA"),
                             output_dir = NULL,
                             overwrite = FALSE,
                             min_obs = 111L) {
  if (!file.exists(zipfile)) stop("Archive not found: ", zipfile, call. = FALSE)
  tickers <- as.character(tickers)
  entries <- utils::unzip(zipfile, list = TRUE)$Name
  entries <- sort(entries[grepl("^Zipped data .*\\.csv$", entries)])
  if (!length(entries)) stop("No weekly CSV entries were found in the archive.", call. = FALSE)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    for (ticker in tickers) {
      prefix <- tolower(ticker)
      for (suffix in c("_clean_rth_3s.csv", "_excluded_days_rth_3s.csv")) {
        path <- file.path(output_dir, paste0(prefix, suffix))
        if (overwrite && file.exists(path)) unlink(path)
      }
    }
  }
  accepted <- setNames(vector("list", length(tickers)), tickers)
  excluded <- setNames(vector("list", length(tickers)), tickers)
  for (entry in entries) {
    week <- .read_tos_entry(zipfile, entry, tickers)
    if (is.null(week) || !nrow(week)) next
    for (ticker in tickers) {
      cleaned <- clean_intraday_data(
        data.frame(date = week$date, time = week$time, price = week[[ticker]]),
        min_obs = min_obs
      )
      if (!is.null(output_dir)) {
        prefix <- tolower(ticker)
        .append_csv(cleaned$data, file.path(output_dir, paste0(prefix, "_clean_rth_3s.csv")))
        .append_csv(cleaned$excluded, file.path(output_dir, paste0(prefix, "_excluded_days_rth_3s.csv")))
      } else {
        accepted[[ticker]] <- c(accepted[[ticker]], list(cleaned$data))
        excluded[[ticker]] <- c(excluded[[ticker]], list(cleaned$excluded))
      }
    }
  }
  if (is.null(output_dir)) {
    accepted <- lapply(accepted, function(x) if (length(x)) do.call(rbind, x) else data.frame())
    excluded <- lapply(excluded, function(x) if (length(x)) do.call(rbind, x) else data.frame())
  }
  summary <- do.call(rbind, lapply(tickers, function(ticker) {
    if (is.null(output_dir)) {
      data <- accepted[[ticker]]
      bad <- excluded[[ticker]]
      data.frame(ticker = ticker, rows_saved = nrow(data),
                 valid_days = length(unique(data$date)), excluded_days = nrow(bad))
    } else {
      data_path <- file.path(output_dir, paste0(tolower(ticker), "_clean_rth_3s.csv"))
      bad_path <- file.path(output_dir, paste0(tolower(ticker), "_excluded_days_rth_3s.csv"))
      data <- if (file.exists(data_path)) utils::read.csv(data_path) else data.frame()
      bad <- if (file.exists(bad_path)) utils::read.csv(bad_path) else data.frame()
      data.frame(ticker = ticker, rows_saved = nrow(data),
                 valid_days = if (nrow(data)) length(unique(data$date)) else 0L,
                 excluded_days = nrow(bad))
    }
  }))
  out <- list(summary = summary, output_dir = output_dir)
  if (is.null(output_dir)) {
    out$data <- accepted
    out$excluded <- excluded
  }
  out
}

.daily_test_one <- function(day, method, test_args) {
  base <- list(prices = day$price, log_prices = TRUE)
  test_fun <- if (method == "noise") adaptive_jump_test_noise else adaptive_jump_test
  result <- tryCatch(do.call(test_fun, c(base, test_args)), error = function(e) e)
  if (inherits(result, "error")) {
    return(data.frame(
      date = day$date[1L], start_time = day$time[1L], end_time = day$time[nrow(day)],
      n_obs = nrow(day), aj_p = NA_real_, lm_p = NA_real_, cc_p = NA_real_,
      q_hat = NA_real_, error = conditionMessage(result), stringsAsFactors = FALSE
    ))
  }
  if (method == "noise") {
    data.frame(
      date = day$date[1L], start_time = day$time[1L], end_time = day$time[nrow(day)],
      n_obs = nrow(day), aj_p = result$ajj$p_value, lm_p = result$lm$p_value,
      cc_p = result$cc$p_value, q_hat = result$lm$q_hat, error = "",
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      date = day$date[1L], start_time = day$time[1L], end_time = day$time[nrow(day)],
      n_obs = nrow(day), aj_p = result$aj$p_value, lm_p = result$lm$p_value,
      cc_p = result$cc$p_value, q_hat = NA_real_, error = "", stringsAsFactors = FALSE
    )
  }
}

#' Apply daily jump tests to an intraday data set
#'
#' @param data A data frame containing `date`, `time`, and `price`.
#' @param method `"noise"` for AJJ/noisy-LM/CC or `"frictionless"` for
#'   AJ/LM/CC.
#' @param test_args Named list passed to the chosen adaptive test.
#' @param alpha Rejection levels to summarize.
#' @param workers Number of parallel day-level workers.
#' @param clean Whether to apply `clean_intraday_data()` first.
#'
#' @return A list with daily p-values, excluded days, and rejection rates.
#' @export
daily_jump_tests <- function(data,
                             method = c("noise", "frictionless"),
                             test_args = list(),
                             alpha = c(0.01, 0.05, 0.10),
                             workers = 1L,
                             clean = TRUE) {
  method <- match.arg(method)
  workers <- as.integer(workers)
  if (workers < 1L) stop("`workers` must be positive.", call. = FALSE)
  cleaned <- if (clean) clean_intraday_data(data) else list(data = data, excluded = data.frame())
  daily <- split(cleaned$data, cleaned$data$date)
  daily <- daily[sort(names(daily))]
  worker <- function(day) .daily_test_one(day, method, test_args)
  if (!length(daily)) {
    return(list(p_values = data.frame(), excluded = cleaned$excluded, rejection_rates = data.frame()))
  }
  if (workers == 1L) {
    values <- lapply(daily, worker)
  } else {
    cl <- parallel::makeCluster(min(workers, length(daily)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(adaptivejump))
    values <- parallel::parLapplyLB(cl, daily, worker)
  }
  p_values <- do.call(rbind, values)
  methods <- c("aj", "lm", "cc")
  rejection <- do.call(rbind, lapply(alpha, function(a) {
    data.frame(
      alpha = a,
      method = toupper(methods),
      rejection_rate = vapply(methods, function(m) mean(p_values[[paste0(m, "_p")]] <= a,
                                                          na.rm = TRUE), numeric(1L)),
      stringsAsFactors = FALSE
    )
  }))
  list(p_values = p_values, excluded = cleaned$excluded, rejection_rates = rejection)
}
