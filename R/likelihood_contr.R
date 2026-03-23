#' Compose a likelihood model from heterogeneous observation contributions
#'
#' Constructs a `"likelihood_contr"` model by combining named `"contr"`
#' objects. The model splits a data frame by observation type, evaluates
#' each type's contribution, and sums the results (under the i.i.d.
#' assumption).
#'
#' Observation type dispatch is controlled by `obs_type`:
#' - If a **string**, it names a column in the data frame whose values
#'   are matched against the contribution names.
#' - If a **function**, it is called as `obs_type(df)` and must return a
#'   character vector of the same length as `nrow(df)`.
#'
#' @param obs_type Either a string (column name) or function for
#'   determining observation types. See Details.
#' @param ... Named `"contr"` objects, one per observation type.
#' @param rdata_fn Optional function `function(theta, n, ...)` that
#'   generates random data from the model's data-generating process.
#'   Required by the default [fim.likelihood_model()] method.
#' @param assumptions Character vector of model assumptions. `"iid"` is
#'   always included.
#' @return A `"likelihood_contr"` likelihood model object.
#'
#'
#'
#' @examples
#' # Weibull model with exact and right-censored observations
#' model <- likelihood_contr(
#'   obs_type = "status",
#'   exact = contr_name("weibull", "exact", ob_col = "t"),
#'   right = contr_name("weibull", "right", ob_col = "t"),
#'   assumptions = c("Weibull distribution", "non-informative censoring")
#' )
#'
#' @export
likelihood_contr <- function(obs_type, ..., rdata_fn = NULL,
                             assumptions = character(0)) {
  if (!is.character(obs_type) && !is.function(obs_type)) {
    stop("'obs_type' must be a string (column name) or a function")
  }
  if (is.character(obs_type)) {
    if (length(obs_type) != 1L || is.na(obs_type) || nchar(obs_type) == 0L) {
      stop("'obs_type' string must be a single non-empty column name")
    }
  }

  contrs <- list(...)
  if (length(contrs) == 0L) {
    stop("At least one contribution must be provided")
  }
  if (is.null(names(contrs)) || any(names(contrs) == "")) {
    stop("All contributions must be named")
  }
  for (nm in names(contrs)) {
    if (!inherits(contrs[[nm]], "contr")) {
      stop(sprintf("Contribution '%s' is not a 'contr' object", nm))
    }
  }

  if (!is.null(rdata_fn) && !is.function(rdata_fn)) {
    stop("'rdata_fn' must be a function or NULL")
  }

  assumptions <- unique(c("iid", as.character(assumptions)))

  structure(
    list(
      obs_type = obs_type,
      contrs = contrs,
      rdata_fn = rdata_fn,
      assumptions = assumptions
    ),
    class = c("likelihood_contr", "likelihood_model")
  )
}


#' Log-likelihood for likelihood_contr
#'
#' Returns a closure that splits the data frame by observation type,
#' evaluates each contribution's log-likelihood, and sums them.
#' The data-frame split is cached so repeated calls with the same `df`
#' (e.g., during optimization) skip the split step.
#'
#' @param model A `likelihood_contr` model.
#' @param ... Additional arguments (currently unused).
#' @return A function `function(df, par, ...)` returning the total
#'   log-likelihood (scalar).
#' @examples
#' model <- likelihood_contr(
#'   obs_type = "status",
#'   exact = contr_name("exp", "exact", ob_col = "t"),
#'   right = contr_name("exp", "right", ob_col = "t")
#' )
#' df <- data.frame(t = c(1, 2, 3, 4), status = c("exact", "exact", "right", "right"))
#' ll_fn <- loglik(model)
#' ll_fn(df, par = c(rate = 0.5))
#' @method loglik likelihood_contr
#' @export
loglik.likelihood_contr <- function(model, ...) {
  cache <- .new_cache()
  function(df, par, ...) {
    cache <- .ensure_cache(cache, df, model)
    total <- 0
    for (type in names(cache$dfs)) {
      sub_df <- cache$dfs[[type]]
      if (nrow(sub_df) == 0L) next
      ctr <- model$contrs[[type]]
      total <- total + ctr$loglik(sub_df, par, ...)
    }
    total
  }
}


#' Score for likelihood_contr
#'
#' Returns a closure that splits the data frame by observation type and
#' sums the per-type score vectors. If a contribution does not provide
#' an analytical score, numerical differentiation via [numDeriv::grad()]
#' is used.
#'
#' @param model A `likelihood_contr` model.
#' @param ... Additional arguments (currently unused).
#' @return A function `function(df, par, ...)` returning the total score
#'   vector.
#' @examples
#' model <- likelihood_contr(
#'   obs_type = "status",
#'   exact = contr_name("exp", "exact", ob_col = "t")
#' )
#' df <- data.frame(t = c(0.5, 1.0, 1.5), status = "exact")
#' score_fn <- score(model)
#' score_fn(df, par = c(rate = 2))
#' @method score likelihood_contr
#' @export
score.likelihood_contr <- function(model, ...) {
  cache <- .new_cache()
  function(df, par, ...) {
    cache <- .ensure_cache(cache, df, model)
    total <- rep(0, length(par))
    names(total) <- names(par)
    for (type in names(cache$dfs)) {
      sub_df <- cache$dfs[[type]]
      if (nrow(sub_df) == 0L) next
      ctr <- model$contrs[[type]]
      if (!is.null(ctr$score)) {
        total <- total + ctr$score(sub_df, par, ...)
      } else {
        ll_type <- ctr$loglik
        total <- total + numDeriv::grad(
          func = function(p) ll_type(sub_df, p, ...),
          x = par,
          method.args = list(r = 6)
        )
      }
    }
    total
  }
}


#' Hessian of the log-likelihood for likelihood_contr
#'
#' Returns a closure that splits the data frame by observation type and
#' sums the per-type Hessian matrices. If a contribution does not provide
#' an analytical Hessian, numerical differentiation via
#' [numDeriv::hessian()] is used.
#'
#' @param model A `likelihood_contr` model.
#' @param ... Additional arguments (currently unused).
#' @return A function `function(df, par, ...)` returning the total Hessian
#'   matrix.
#' @examples
#' model <- likelihood_contr(
#'   obs_type = "status",
#'   exact = contr_name("exp", "exact", ob_col = "t")
#' )
#' df <- data.frame(t = c(0.5, 1.0, 1.5), status = "exact")
#' hess_fn <- hess_loglik(model)
#' hess_fn(df, par = c(rate = 2))
#' @method hess_loglik likelihood_contr
#' @export
hess_loglik.likelihood_contr <- function(model, ...) {
  cache <- .new_cache()
  function(df, par, ...) {
    cache <- .ensure_cache(cache, df, model)
    p <- length(par)
    total <- matrix(0, nrow = p, ncol = p)
    for (type in names(cache$dfs)) {
      sub_df <- cache$dfs[[type]]
      if (nrow(sub_df) == 0L) next
      ctr <- model$contrs[[type]]
      if (!is.null(ctr$hess)) {
        total <- total + ctr$hess(sub_df, par, ...)
      } else {
        ll_type <- ctr$loglik
        total <- total + numDeriv::hessian(
          func = function(p) ll_type(sub_df, p, ...),
          x = par,
          method.args = list(r = 6)
        )
      }
    }
    total
  }
}


#' Random data generation for likelihood_contr
#'
#' Returns the user-supplied `rdata_fn` or errors if none was provided.
#'
#' @param model A `likelihood_contr` model.
#' @param ... Additional arguments (currently unused).
#' @return A function `function(theta, n, ...)` returning a data frame.
#' @examples
#' rdata_fn <- function(theta, n, ...) {
#'   data.frame(t = rexp(n, rate = theta[1]), status = "exact")
#' }
#' model <- likelihood_contr(
#'   obs_type = "status",
#'   exact = contr_name("exp", "exact", ob_col = "t"),
#'   rdata_fn = rdata_fn
#' )
#' gen <- rdata(model)
#' gen(theta = c(rate = 2), n = 5)
#' @method rdata likelihood_contr
#' @export
rdata.likelihood_contr <- function(model, ...) {
  if (is.null(model$rdata_fn)) {
    stop("This likelihood_contr model does not have an rdata function")
  }
  model$rdata_fn
}


#' Assumptions for likelihood_contr
#'
#' Returns the character vector of model assumptions stored in the
#' `likelihood_contr` object. The `"iid"` assumption is always included.
#'
#' @param model A `likelihood_contr` model.
#' @param ... Additional arguments (currently unused).
#' @return Character vector of assumptions.
#' @examples
#' model <- likelihood_contr(
#'   obs_type = "status",
#'   exact = contr_name("exp", "exact", ob_col = "t"),
#'   assumptions = c("exponential distribution", "non-informative censoring")
#' )
#' assumptions(model)
#' @method assumptions likelihood_contr
#' @export
assumptions.likelihood_contr <- function(model, ...) {
  model$assumptions
}


#' Print method for likelihood_contr
#'
#' Prints a summary of the model's observation types, dispatch method
#' (column name or function), and assumptions.
#'
#' @param x A `likelihood_contr` model.
#' @param ... Additional arguments (ignored).
#' @return The model object, invisibly.
#' @examples
#' model <- likelihood_contr(
#'   obs_type = "status",
#'   exact = contr_name("weibull", "exact", ob_col = "t"),
#'   right = contr_name("weibull", "right", ob_col = "t"),
#'   assumptions = c("Weibull distribution")
#' )
#' print(model)
#' @export
print.likelihood_contr <- function(x, ...) {
  cat("Likelihood Contribution Model\n")
  cat("-----------------------------\n")

  type_names <- names(x$contrs)
  cat("Observation types:", paste(type_names, collapse = ", "), "\n")

  if (is.function(x$obs_type)) {
    cat("Dispatch method: function\n")
  } else {
    cat("Dispatch method: column '", x$obs_type, "'\n", sep = "")
  }

  cat("Assumptions:\n")
  for (a in x$assumptions) {
    cat(" -", trimws(a), "\n")
  }

  invisible(x)
}
