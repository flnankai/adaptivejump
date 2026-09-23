# Internal helpers shared by the frictionless and noisy procedures.

.clip_p <- function(p, eps = 1e-15) {
  pmin(pmax(p, eps), 1 - eps)
}

.normal_abs_moment <- function(r) {
  2^(r / 2) * gamma((r + 1) / 2) / sqrt(pi)
}

.as_numeric_path <- function(x, log_prices = TRUE) {
  x <- as.numeric(x)
  if (length(x) < 12L || any(!is.finite(x))) {
    stop("`prices` must contain at least 12 finite observations.", call. = FALSE)
  }
  if (isTRUE(log_prices)) {
    if (any(x <= 0)) {
      stop("Raw prices must be strictly positive when `log_prices = TRUE`.", call. = FALSE)
    }
    x <- log(x)
  }
  x
}

.resolve_delta <- function(n_returns, delta = NULL) {
  if (is.null(delta)) {
    return(1 / n_returns)
  }
  delta <- as.numeric(delta)
  if (length(delta) != 1L || !is.finite(delta) || delta <= 0) {
    stop("`delta` must be one positive finite number.", call. = FALSE)
  }
  delta
}

.default_workers <- function() {
  n <- parallel::detectCores(logical = TRUE)
  if (is.na(n) || n < 2L) 1L else as.integer(n - 1L)
}

.new_jump_test <- function(method, statistic, p_value, alpha = NULL, ...) {
  out <- list(
    method = method,
    statistic = as.numeric(statistic),
    p_value = as.numeric(p_value),
    alpha = alpha,
    reject = if (is.null(alpha)) NA else is.finite(p_value) && p_value <= alpha,
    ...
  )
  class(out) <- "jump_test"
  out
}

#' Print a jump-test result
#'
#' @param x A `jump_test` object.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.jump_test <- function(x, ...) {
  cat(sprintf("%s jump test\n", x$method))
  cat(sprintf("  statistic: %.6f\n", x$statistic))
  cat(sprintf("  p-value:   %.6g\n", x$p_value))
  if (!is.null(x$alpha)) {
    cat(sprintf("  reject at %.3g: %s\n", x$alpha, if (isTRUE(x$reject)) "yes" else "no"))
  }
  invisible(x)
}

#' Combine p-values with the Cauchy combination rule
#'
#' Computes the weighted Cauchy statistic
#' \eqn{\sum_j w_j\tan\{\pi(1/2-p_j)\}} and its standard-Cauchy upper-tail
#' p-value.  The p-values need not be independent for the usual Cauchy
#' combination approximation.
#'
#' @param p_values Numeric vector of p-values.
#' @param weights Optional nonnegative weights.  Equal weights are used by
#'   default and supplied weights are normalized to sum to one.
#' @param alpha Optional decision level.
#'
#' @return An object of class `jump_test` with method `CC`.
#' @export
cauchy_combine <- function(p_values, weights = NULL, alpha = NULL) {
  p_values <- as.numeric(p_values)
  if (!length(p_values) || any(!is.finite(p_values))) {
    stop("`p_values` must be a nonempty vector of finite values.", call. = FALSE)
  }
  if (is.null(weights)) {
    weights <- rep(1 / length(p_values), length(p_values))
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != length(p_values) || any(!is.finite(weights)) ||
        any(weights < 0) || sum(weights) <= 0) {
      stop("`weights` must be nonnegative, finite, and sum to a positive value.", call. = FALSE)
    }
    weights <- weights / sum(weights)
  }

  statistic <- sum(weights * tan(pi * (0.5 - .clip_p(p_values))))
  p_value <- .clip_p(0.5 - atan(statistic) / pi)
  .new_jump_test("CC", statistic, p_value, alpha = alpha,
                 component_p_values = p_values, weights = weights)
}

#' Benjamini--Hochberg FDR selection
#'
#' @param p_values Numeric vector of p-values.
#' @param fdr Target false discovery rate in `(0, 1)`.
#'
#' @return A list containing the logical selection vector, BH adjusted p-values,
#'   and the selected p-value cutoff.
#' @export
bh_select <- function(p_values, fdr = 0.2) {
  p_values <- as.numeric(p_values)
  if (any(!is.finite(p_values))) {
    stop("`p_values` must be finite.", call. = FALSE)
  }
  if (!is.finite(fdr) || fdr <= 0 || fdr >= 1) {
    stop("`fdr` must lie strictly between 0 and 1.", call. = FALSE)
  }
  n <- length(p_values)
  ord <- order(p_values)
  threshold <- fdr * seq_len(n) / n
  pass <- p_values[ord] <= threshold
  cutoff <- if (any(pass)) max(p_values[ord][pass]) else NA_real_
  list(
    selected = if (is.na(cutoff)) rep(FALSE, n) else p_values <= cutoff,
    adjusted_p = stats::p.adjust(p_values, method = "BH"),
    cutoff = cutoff,
    fdr = fdr
  )
}
