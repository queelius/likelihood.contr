#' @keywords internal
"_PACKAGE"

#' @importFrom likelihood.model loglik score hess_loglik assumptions rdata
NULL

#' Re-exported generics from likelihood.model
#'
#' These generics are re-exported so that users can call them directly
#' after loading \code{likelihood.contr}, without loading
#' \code{likelihood.model} separately.
#'
#' \describe{
#'   \item{\code{loglik(model, ...)}}{Returns a closure
#'     \code{function(df, par, ...)} that evaluates the total log-likelihood
#'     (a scalar).}
#'   \item{\code{score(model, ...)}}{Returns a closure
#'     \code{function(df, par, ...)} that evaluates the score vector
#'     (gradient of the log-likelihood).}
#'   \item{\code{hess_loglik(model, ...)}}{Returns a closure
#'     \code{function(df, par, ...)} that evaluates the Hessian matrix of the
#'     log-likelihood.}
#'   \item{\code{assumptions(model, ...)}}{Returns a character vector of model
#'     assumptions.}
#'   \item{\code{rdata(model, ...)}}{Returns a closure
#'     \code{function(theta, n, ...)} that generates a random data frame from
#'     the model.}
#' }
#'
#' See \code{\link[likelihood.model]{loglik}} for full documentation.
#'
#' @usage loglik(model, ...)
#' score(model, ...)
#' hess_loglik(model, ...)
#' assumptions(model, ...)
#' rdata(model, ...)
#' @param model A likelihood model object (e.g., a \code{"likelihood_contr"}
#'   object).
#' @param ... Additional arguments passed to methods.
#' @return The return type depends on the generic; see Details.
#'   \code{loglik}, \code{score}, \code{hess_loglik}, and \code{rdata} each
#'   return a \strong{closure} (function). \code{assumptions} returns a
#'   character vector.
#'
#' @name loglik
#' @aliases loglik score hess_loglik assumptions rdata
#' @export loglik score hess_loglik assumptions rdata
NULL
