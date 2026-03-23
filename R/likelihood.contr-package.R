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
#'   \item{\code{loglik}}{Returns a log-likelihood closure for the model.}
#'   \item{\code{score}}{Returns a score (gradient) closure for the model.}
#'   \item{\code{hess_loglik}}{Returns a Hessian closure for the model.}
#'   \item{\code{assumptions}}{Returns the character vector of model assumptions.}
#'   \item{\code{rdata}}{Returns a random data generation closure.}
#' }
#'
#' See \code{\link[likelihood.model]{loglik}} for full documentation.
#'
#' @name loglik
#' @aliases loglik score hess_loglik assumptions rdata
#' @export loglik score hess_loglik assumptions rdata
NULL
