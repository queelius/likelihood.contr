#' Create a new cache environment for likelihood_contr
#'
#' The cache stores the most recent data frame split to avoid redundant
#' `split()` and observation-type dispatch during repeated optimizer calls.
#'
#' @return An environment with `df_ref` and `dfs` slots.
#' @keywords internal
.new_cache <- function() {
  env <- new.env(parent = emptyenv())
  env$df_ref <- NULL
  env$dfs <- NULL
  env
}


#' Ensure the cache is populated for the given data frame
#'
#' If `df` is a different object from the cached one (checked via
#' [identical()]), splits the data frame by observation type and validates
#' the split. Otherwise returns the existing cache unchanged.
#'
#' @param cache Cache environment from [.new_cache()].
#' @param df Data frame.
#' @param model A `likelihood_contr` object.
#' @return The cache environment (invisibly updated), with `$dfs` populated
#'   as a named list of data-frame subsets keyed by type.
#' @keywords internal
.ensure_cache <- function(cache, df, model) {
  if (!is.data.frame(df)) stop("Observations must be a data frame")
  if (identical(df, cache$df_ref)) {
    return(cache)
  }

  # Determine types
  if (is.function(model$obs_type)) {
    types <- model$obs_type(df)
    if (!is.character(types) || length(types) != nrow(df)) {
      stop("Observation callback must return a character vector with one entry per row")
    }
  } else {
    col <- model$obs_type
    if (!col %in% names(df)) {
      stop(sprintf("Column '%s' not found in data frame", col))
    }
    types <- df[[col]]
  }

  if (any(is.na(types))) {
    stop("NA values in observation types are not allowed")
  }

  # Check all types have matching contributions
  unique_types <- unique(types)
  contr_names <- names(model$contrs)
  unmatched <- setdiff(unique_types, contr_names)
  if (length(unmatched) > 0L) {
    stop(sprintf(
      "No contribution defined for observation type(s): %s",
      paste(sQuote(unmatched), collapse = ", ")
    ))
  }

  cache$df_ref <- df
  cache$dfs <- split(df, types)
  cache
}
