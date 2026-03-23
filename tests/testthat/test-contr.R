# ---- contr_fn tests --------------------------------------------------------

test_that("contr_fn creates a contr object from functions", {
  ll <- function(df, par, ...) sum(dexp(df$x, rate = par[1], log = TRUE))
  ctr <- contr_fn(loglik = ll)
  expect_s3_class(ctr, "contr")
  expect_identical(ctr$loglik, ll)
  expect_null(ctr$score)
  expect_null(ctr$hess)
})

test_that("contr_fn stores score and hess when provided", {
  ll <- function(df, par, ...) sum(dexp(df$x, rate = par[1], log = TRUE))
  sc <- function(df, par, ...) c(nrow(df) / par[1] - sum(df$x))
  hs <- function(df, par, ...) matrix(-nrow(df) / par[1]^2, 1, 1)
  ctr <- contr_fn(loglik = ll, score = sc, hess = hs)
  expect_identical(ctr$score, sc)
  expect_identical(ctr$hess, hs)
})

test_that("contr_fn validates loglik is a function", {
  expect_error(contr_fn(loglik = "not a function"), "'loglik' must be a function")
  expect_error(contr_fn(loglik = 42), "'loglik' must be a function")
  expect_error(contr_fn(loglik = NULL), "'loglik' must be a function")
})

test_that("contr_fn validates score is a function or NULL", {
  ll <- function(df, par, ...) 0
  expect_error(contr_fn(loglik = ll, score = "bad"), "'score' must be a function")
})

test_that("contr_fn validates hess is a function or NULL", {
  ll <- function(df, par, ...) 0
  expect_error(contr_fn(loglik = ll, hess = 123), "'hess' must be a function")
})

test_that("contr_fn loglik computes correct values", {
  ll <- function(df, par, ...) sum(dnorm(df$x, mean = par[1], sd = par[2], log = TRUE))
  ctr <- contr_fn(loglik = ll)
  df <- data.frame(x = c(0, 1, -1))
  result <- ctr$loglik(df, c(0, 1))
  expected <- sum(dnorm(c(0, 1, -1), 0, 1, log = TRUE))
  expect_equal(result, expected)
})


# ---- contr_name tests ------------------------------------------------------

test_that("contr_name creates exact contribution", {
  ctr <- contr_name("norm", "exact", ob_col = "x")
  expect_s3_class(ctr, "contr")
  df <- data.frame(x = c(0, 1, -1))
  result <- ctr$loglik(df, c(0, 1))
  expected <- sum(dnorm(c(0, 1, -1), 0, 1, log = TRUE))
  expect_equal(result, expected)
})

test_that("contr_name creates right-censored contribution", {
  ctr <- contr_name("weibull", "right", ob_col = "t")
  expect_s3_class(ctr, "contr")
  df <- data.frame(t = c(1, 2, 3))
  result <- ctr$loglik(df, c(shape = 2, scale = 1))
  expected <- sum(pweibull(c(1, 2, 3), shape = 2, scale = 1,
                           lower.tail = FALSE, log.p = TRUE))
  expect_equal(result, expected)
})

test_that("contr_name creates left-censored contribution", {
  ctr <- contr_name("norm", "left", ob_col = "x")
  df <- data.frame(x = c(0.5, 1.0))
  result <- ctr$loglik(df, c(0, 1))
  expected <- sum(pnorm(c(0.5, 1.0), 0, 1, log.p = TRUE))
  expect_equal(result, expected)
})

test_that("contr_name creates interval-censored contribution", {
  ctr <- contr_name("norm", "interval", ob_col = "lo", ob_col_upper = "hi")
  df <- data.frame(lo = c(0, 1), hi = c(1, 2))
  result <- ctr$loglik(df, c(0, 1))
  expected <- sum(log(pnorm(c(1, 2), 0, 1) - pnorm(c(0, 1), 0, 1)))
  expect_equal(result, expected)
})

test_that("contr_name interval censoring is numerically stable in tails", {
  ctr <- contr_name("norm", "interval", ob_col = "lo", ob_col_upper = "hi")
  # Far right tail: naive log(F_hi - F_lo) underflows to -Inf
  df <- data.frame(lo = 37, hi = 38)
  result <- ctr$loglik(df, c(0, 1))
  expect_true(is.finite(result))
  expect_true(result < 0)
})

test_that("contr_name validates type", {
  expect_error(
    contr_name("norm", "unknown"),
    "Invalid type 'unknown'"
  )
})

test_that("contr_name validates distribution exists", {
  expect_error(
    contr_name("bogus_dist_xyz", "exact"),
    "Distribution 'bogus_dist_xyz' not found"
  )
})

test_that("contr_name errors when CDF function is missing", {
  # Create a distribution with d<name> but no p<name>
  assign("dorphan_test_dist", dnorm, envir = globalenv())
  on.exit(rm("dorphan_test_dist", envir = globalenv()))
  expect_error(
    contr_name("orphan_test_dist", "exact"),
    "no function 'porphan_test_dist'"
  )
})

test_that("contr_name requires ob_col_upper for interval type", {
  expect_error(
    contr_name("norm", "interval", ob_col = "x"),
    "'ob_col_upper' is required"
  )
})

test_that("contr_name does not provide score or hess (numerical fallback)", {
  ctr <- contr_name("norm", "exact", ob_col = "x")
  expect_null(ctr$score)
  expect_null(ctr$hess)
})


# ---- .prepare_args_list tests -----------------------------------------------

test_that(".prepare_args_list maps unnamed parameters by position", {
  args <- likelihood.contr:::.prepare_args_list(c(2, 1), stats::dweibull)
  expect_equal(args, list(shape = 2, scale = 1))
})

test_that(".prepare_args_list preserves named parameters", {
  args <- likelihood.contr:::.prepare_args_list(c(scale = 3, shape = 2), stats::dweibull)
  expect_equal(args, list(scale = 3, shape = 2))
})

test_that(".prepare_args_list errors on too many parameters", {
  expect_error(
    likelihood.contr:::.prepare_args_list(c(1, 2, 3), stats::dexp),
    "Too many parameters"
  )
})

test_that(".prepare_args_list works with single-parameter distribution", {
  args <- likelihood.contr:::.prepare_args_list(c(2), stats::dexp)
  expect_equal(args, list(rate = 2))
})
