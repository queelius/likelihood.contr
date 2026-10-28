#' Weight a contribution by observed frequencies
#'
#' Wrap an existing contribution to work with counted observations, such as
#' histogram bins. Each row contributes its frequency times its log-likelihood.
#' No expanded data frame is constructed. Analytical derivatives are weighted
#' in the same way; absent derivatives remain available through the model's
#' numerical differentiation fallback.
#'
#' @details Frequencies must be finite, nonnegative integers. This represents
#' repeated independent observations, not survey weights or importance weights.
#' The wrapped contribution must be additive over independent rows. Zero-count
#' rows are skipped, including rows whose likelihood would otherwise be undefined.
#' Multinomial histogram likelihoods differ from this likelihood by a constant
#' independent of the distribution parameters. Rows must describe disjoint bins
#' and counts must be per-bin frequencies, not cumulative histogram counts.
#'
#' The generic fitting layer counts rows, not expanded observations. Do not
#' use its default BIC or row bootstrap for frequency data: use the sum of
#' frequencies for sample size and resample counts under the sampling design.
#' Direct likelihood optimization is illustrated in the timeout benchmarking
#' vignette.
#' @param contribution A contribution object created by [contr_fn()] or
#'   [contr_name()].
#' @param count_col Name of the column containing observed frequencies.
#' @return A `"contr"` object usable in [likelihood_contr()].
#' @examples
#' bins <- data.frame(lo = c(0, 1, 2), hi = c(1, 2, Inf), count = c(12, 5, 3))
#' interval <- contr_name("exp", "interval", "lo", "hi")
#' counted <- contr_freq(interval, "count")
#' counted$loglik(bins, c(rate = 0.7))
#' @export
contr_freq <- function(contribution, count_col) {
  if (!inherits(contribution, "contr")) {
    stop("'contribution' must be a 'contr' object")
  }
  if (!is.character(count_col) || length(count_col) != 1L ||
      is.na(count_col) || !nzchar(count_col)) {
    stop("'count_col' must be a single non-empty column name")
  }
  weighted <- function(fn, zero) {
    if (is.null(fn)) return(NULL)
    function(df, par, ...) {
      if (!is.data.frame(df) || !count_col %in% names(df)) {
        stop("The count column must be present in a data frame")
      }
      counts <- df[[count_col]]
      if (!is.numeric(counts) || any(!is.finite(counts)) ||
          any(counts < 0 | counts != floor(counts))) {
        stop("The count column must contain finite nonnegative integer frequencies")
      }
      result <- zero(par)
      for (i in which(counts > 0)) {
        result <- result + counts[i] * fn(df[i, , drop = FALSE], par, ...)
      }
      result
    }
  }
  contr_fn(
    loglik = weighted(contribution$loglik, function(par) 0),
    score = weighted(contribution$score, function(par) {
      stats::setNames(numeric(length(par)), names(par))
    }),
    hess = weighted(contribution$hess, function(par) {
      matrix(0, length(par), length(par))
    })
  )
}
