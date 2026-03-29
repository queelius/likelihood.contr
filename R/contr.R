#' Create a likelihood contribution from user-supplied functions
#'
#' Wraps user-provided log-likelihood, score, and Hessian functions into
#' a `"contr"` S3 object suitable for use in [likelihood_contr()].
#'
#' The `loglik` function must have signature `function(df, par, ...)` and
#' return a scalar log-likelihood value. The optional `score` and `hess`
#' functions must have the same signature and return a numeric vector
#' (score) or matrix (Hessian), respectively.
#'
#' @param loglik A function `function(df, par, ...)` returning the
#'   log-likelihood contribution (scalar).
#' @param score Optional function `function(df, par, ...)` returning the
#'   score vector. If `NULL`, numerical differentiation via
#'   [numDeriv::grad()] is used.
#' @param hess Optional function `function(df, par, ...)` returning the
#'   Hessian matrix. If `NULL`, numerical differentiation via
#'   [numDeriv::hessian()] is used.
#' @return A `"contr"` S3 object (a list) with elements `$loglik` (function),
#'   `$score` (function or `NULL`), and `$hess` (function or `NULL`). Pass
#'   this object to [likelihood_contr()] to include it in a composed model.
#' @examples
#' # Exponential exact-observation contribution
#' my_contr <- contr_fn(
#'   loglik = function(df, par, ...) {
#'     sum(dexp(df$x, rate = par[1], log = TRUE))
#'   },
#'   score = function(df, par, ...) {
#'     n <- nrow(df)
#'     c(rate = n / par[1] - sum(df$x))
#'   }
#' )
#' @export
contr_fn <- function(loglik, score = NULL, hess = NULL) {
  if (!is.function(loglik)) {
    stop("'loglik' must be a function")
  }
  if (!is.null(score) && !is.function(score)) {
    stop("'score' must be a function or NULL")
  }
  if (!is.null(hess) && !is.function(hess)) {
    stop("'hess' must be a function or NULL")
  }
  structure(
    list(loglik = loglik, score = score, hess = hess),
    class = "contr"
  )
}


#' Create a likelihood contribution from a named R distribution
#'
#' Generates a `"contr"` object for a standard R distribution using the
#' `d<name>` (PDF) and `p<name>` (CDF) functions. Exact observations use
#' the PDF, right-censored use the survival function, left-censored use
#' the CDF, and interval-censored use the CDF difference.
#'
#' @details
#' The log-likelihood form for each contribution type:
#'
#' - `"exact"`: PDF via `d<name>(x, ..., log = TRUE)`
#' - `"right"`: Survival function via `p<name>(x, ..., lower.tail = FALSE, log.p = TRUE)`
#' - `"left"`: CDF via `p<name>(x, ..., log.p = TRUE)`
#' - `"interval"`: `log(F(upper) - F(lower))` via `p<name>`, computed in
#'   log-space for numerical stability
#'
#' @param dist_name The distribution name (e.g., `"norm"`, `"weibull"`, `"exp"`).
#'   Must have both `d<name>` and `p<name>` functions available.
#' @param type One of `"exact"`, `"right"`, `"left"`, `"interval"`.
#' @param ob_col Name of the observation column in the data frame.
#'   For interval censoring, this is the lower bound column.
#' @param ob_col_upper Name of the upper bound column for interval censoring.
#'   Required when `type = "interval"`.
#' @return A `"contr"` S3 object (a list) with a `$loglik` function
#'   derived from the named distribution. The `$score` and `$hess` elements
#'   are `NULL`; numerical differentiation is used automatically when the
#'   contribution is part of a [likelihood_contr()] model.
#' @examples
#' # Exact Weibull contribution
#' exact <- contr_name("weibull", "exact", ob_col = "t")
#'
#' # Right-censored Weibull contribution
#' right <- contr_name("weibull", "right", ob_col = "t")
#'
#' # Interval-censored normal contribution
#' interval <- contr_name("norm", "interval", ob_col = "lo", ob_col_upper = "hi")
#' @export
contr_name <- function(dist_name, type, ob_col = "x", ob_col_upper = NULL) {
  stopifnot(is.character(dist_name), length(dist_name) == 1L)
  stopifnot(is.character(type), length(type) == 1L)
  stopifnot(is.character(ob_col), length(ob_col) == 1L)

  valid_types <- c("exact", "right", "left", "interval")
  if (!type %in% valid_types) {
    stop(sprintf(
      "Invalid type '%s'. Must be one of: %s",
      type, paste(sQuote(valid_types), collapse = ", ")
    ))
  }

  pdf_name <- paste0("d", dist_name)
  cdf_name <- paste0("p", dist_name)

  if (!exists(pdf_name, mode = "function")) {
    stop(sprintf(
      "Distribution '%s' not found: no function '%s'",
      dist_name, pdf_name
    ))
  }
  if (!exists(cdf_name, mode = "function")) {
    stop(sprintf(
      "Distribution '%s' not found: no function '%s'",
      dist_name, cdf_name
    ))
  }

  if (type == "interval" && is.null(ob_col_upper)) {
    stop("'ob_col_upper' is required when type = 'interval'")
  }

  pdf_fn <- get(pdf_name, mode = "function")
  cdf_fn <- get(cdf_name, mode = "function")

  ll <- switch(type,
    exact = function(df, par, ...) {
      x <- df[[ob_col]]
      par_args <- .prepare_args_list(par, pdf_fn)
      sum(do.call(pdf_fn, c(list(x), par_args, list(log = TRUE))))
    },
    right = function(df, par, ...) {
      x <- df[[ob_col]]
      par_args <- .prepare_args_list(par, cdf_fn)
      sum(do.call(cdf_fn, c(list(x), par_args,
                             list(lower.tail = FALSE, log.p = TRUE))))
    },
    left = function(df, par, ...) {
      x <- df[[ob_col]]
      par_args <- .prepare_args_list(par, cdf_fn)
      sum(do.call(cdf_fn, c(list(x), par_args, list(log.p = TRUE))))
    },
    interval = function(df, par, ...) {
      lo <- df[[ob_col]]
      hi <- df[[ob_col_upper]]
      par_args <- .prepare_args_list(par, cdf_fn)
      log_F_hi <- do.call(cdf_fn, c(list(hi), par_args, list(log.p = TRUE)))
      log_F_lo <- do.call(cdf_fn, c(list(lo), par_args, list(log.p = TRUE)))
      log_S_hi <- do.call(cdf_fn, c(list(hi), par_args,
                                     list(lower.tail = FALSE, log.p = TRUE)))
      log_S_lo <- do.call(cdf_fn, c(list(lo), par_args,
                                     list(lower.tail = FALSE, log.p = TRUE)))
      # CDF form is stable when F < 0.5 (left tail);
      # survival form is stable when F >= 0.5 (right tail).
      right_tail <- log_F_hi > -log(2)
      result <- ifelse(
        right_tail,
        log_S_lo + log1p(-exp(log_S_hi - log_S_lo)),
        log_F_hi + log1p(-exp(log_F_lo - log_F_hi))
      )
      sum(result)
    }
  )

  contr_fn(loglik = ll)
}


#' Map a parameter vector to distribution function formals
#'
#' Matches elements of `par` to the formals of `fn` by position or name,
#' excluding standard distribution formals (`x`, `q`, `p`, `n`, `log`,
#' `log.p`, `lower.tail`, `...`).
#'
#' @param par Numeric parameter vector (named or unnamed).
#' @param fn A distribution function whose formals define the parameter names.
#' @return A named list suitable for [do.call()].
#' @keywords internal
.prepare_args_list <- function(par, fn) {
  fn_args <- names(formals(fn))
  fn_args <- fn_args[!fn_args %in% c("x", "q", "p", "n", "log", "log.p",
                                      "lower.tail", "...")]

  if (length(par) > length(fn_args)) {
    stop(sprintf(
      "Too many parameters (%d) for distribution (expects at most %d: %s)",
      length(par), length(fn_args), paste(fn_args, collapse = ", ")
    ))
  }

  if (is.null(names(par))) {
    names(par) <- fn_args[seq_along(par)]
  }

  as.list(par)
}
