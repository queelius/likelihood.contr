# Integration tests using likelihood.model generics (fit, lrt, score,
# hess_loglik, fim) on likelihood_contr models.

# fit, lrt, and fim are not re-exported by likelihood.contr, so load explicitly.
library(likelihood.model)

# ---- fit() on a simple exponential model -----------------------------------

test_that("fit() works with likelihood_contr model (exponential exact)", {
  set.seed(42)
  n <- 200
  true_rate <- 2
  df <- data.frame(t = rexp(n, rate = true_rate), status = "exact")

  model <- likelihood_contr(
    obs_type = "status",
    exact = contr_name("exp", "exact", ob_col = "t")
  )

  solver <- fit(model)
  result <- solver(df, par = c(rate = 1))

  # MLE should be close to the true rate
  expect_equal(coef(result), c(rate = true_rate), tolerance = 0.3)
  expect_true(result$converged)
})


# ---- fit() with mixed exact + right-censored data --------------------------

test_that("fit() works with mixed censoring types", {
  set.seed(123)
  n <- 300
  true_rate <- 1.5
  censor_time <- 1.0

  raw_times <- rexp(n, rate = true_rate)
  obs_times <- pmin(raw_times, censor_time)
  status <- ifelse(raw_times <= censor_time, "exact", "right")

  df <- data.frame(t = obs_times, status = status)

  model <- likelihood_contr(
    obs_type = "status",
    exact = contr_name("exp", "exact", ob_col = "t"),
    right = contr_name("exp", "right", ob_col = "t")
  )

  solver <- fit(model)
  result <- solver(df, par = c(rate = 1))

  expect_equal(coef(result), c(rate = true_rate), tolerance = 0.4)
  expect_true(result$converged)
})


# ---- score() via likelihood.model generic ----------------------------------

test_that("score() returns correct gradient via numerical differentiation", {
  model <- likelihood_contr(
    obs_type = "status",
    exact = contr_name("exp", "exact", ob_col = "t")
  )

  df <- data.frame(t = c(0.5, 1.0, 1.5, 2.0), status = "exact")
  par <- c(rate = 2)

  score_fn <- score(model)
  grad_val <- score_fn(df, par)

  # Analytical score for exponential exact: d/dlambda = n/lambda - sum(t)
  expected <- unname(length(df$t) / par[1] - sum(df$t))
  expect_equal(unname(grad_val), expected, tolerance = 1e-5)
})


test_that("score() uses analytical score when provided in contribution", {
  analytical_score <- function(df, par, ...) {
    c(rate = nrow(df) / par[1] - sum(df$t))
  }
  model <- likelihood_contr(
    obs_type = "status",
    exact = contr_fn(
      loglik = function(df, par, ...) sum(dexp(df$t, rate = par[1], log = TRUE)),
      score  = analytical_score
    )
  )

  df <- data.frame(t = c(0.5, 1.0, 2.0), status = "exact")
  score_fn <- score(model)
  result <- score_fn(df, c(rate = 2))
  expected <- c(rate = 3 / 2 - sum(c(0.5, 1.0, 2.0)))
  expect_equal(result, expected)
})


# ---- hess_loglik() via likelihood.model generic ----------------------------

test_that("hess_loglik() returns correct Hessian via numerical differentiation", {
  model <- likelihood_contr(
    obs_type = "status",
    exact = contr_name("exp", "exact", ob_col = "t")
  )

  df <- data.frame(t = c(0.5, 1.0, 1.5), status = "exact")
  par <- c(rate = 2)

  hess_fn <- hess_loglik(model)
  H <- hess_fn(df, par)

  # Analytical Hessian for exponential exact: -n/lambda^2
  expected <- matrix(-3 / 4, 1, 1)
  expect_equal(H, expected, tolerance = 1e-5)
})


test_that("hess_loglik() uses analytical hess when provided", {
  analytical_hess <- function(df, par, ...) {
    matrix(-nrow(df) / par[1]^2, 1, 1)
  }
  model <- likelihood_contr(
    obs_type = "status",
    exact = contr_fn(
      loglik = function(df, par, ...) sum(dexp(df$t, rate = par[1], log = TRUE)),
      hess   = analytical_hess
    )
  )

  df <- data.frame(t = c(0.5, 1.0, 2.0), status = "exact")
  hess_fn <- hess_loglik(model)
  H <- hess_fn(df, c(rate = 2))
  expected <- matrix(-3 / 4, 1, 1)
  expect_equal(H, expected)
})


# ---- lrt() between nested models ------------------------------------------

test_that("lrt() works with likelihood_contr models", {
  set.seed(99)
  n <- 200
  true_shape <- 2
  true_scale <- 3

  df <- data.frame(t = rweibull(n, shape = true_shape, scale = true_scale),
                   status = "exact")

  # Null: exponential (shape = 1)
  null_model <- likelihood_contr(
    obs_type = "status",
    exact = contr_name("exp", "exact", ob_col = "t")
  )

  # Alt: Weibull (2 parameters)
  alt_model <- likelihood_contr(
    obs_type = "status",
    exact = contr_name("weibull", "exact", ob_col = "t")
  )

  null_solver <- fit(null_model)
  alt_solver <- fit(alt_model)

  # NaN warnings expected as optimizer explores negative parameter space
  null_mle <- suppressWarnings(null_solver(df, par = c(rate = 0.5)))
  alt_mle <- suppressWarnings(alt_solver(df, par = c(shape = 1.5, scale = 2)))

  result <- lrt(
    null = null_model,
    alt  = alt_model,
    data = df,
    null_par = coef(null_mle),
    alt_par  = coef(alt_mle),
    dof = 1
  )

  # Weibull with shape != 1 should be significantly better than exponential
  expect_lt(result$p.value, 0.05)
  expect_gt(result$stat, 0)
})


# ---- fim() with rdata_fn ---------------------------------------------------

test_that("fim() works when rdata_fn is provided", {
  rdata_fn <- function(theta, n, ...) {
    data.frame(
      t = rexp(n, rate = theta[1]),
      status = "exact"
    )
  }

  model <- likelihood_contr(
    obs_type = "status",
    exact = contr_name("exp", "exact", ob_col = "t"),
    rdata_fn = rdata_fn
  )

  fim_fn <- fim(model)
  set.seed(42)
  I <- fim_fn(theta = c(rate = 2), n_obs = 100, n_samples = 500)

  # Analytical FIM for exponential: n / lambda^2 = 100 / 4 = 25
  expect_equal(I[1, 1], 25, tolerance = 5)
})


# ---- Multiparameter: Weibull fit -------------------------------------------

test_that("fit() works with 2-parameter Weibull model", {
  set.seed(77)
  n <- 300
  true_shape <- 2.5
  true_scale <- 4.0

  df <- data.frame(t = rweibull(n, true_shape, true_scale), status = "exact")

  model <- likelihood_contr(
    obs_type = "status",
    exact = contr_name("weibull", "exact", ob_col = "t")
  )

  solver <- fit(model)
  result <- solver(df, par = c(shape = 2, scale = 3))

  expect_equal(unname(coef(result)[1]), true_shape, tolerance = 0.5)
  expect_equal(unname(coef(result)[2]), true_scale, tolerance = 0.5)
})


# ---- is_likelihood_model ---------------------------------------------------

test_that("likelihood_contr satisfies is_likelihood_model", {
  ctr <- contr_fn(loglik = function(df, par, ...) 0)
  model <- likelihood_contr(obs_type = "type", a = ctr)
  expect_true(is_likelihood_model(model))
})


# ---- Mixed analytical + numerical contributions ---------------------------

test_that("score sums analytical and numerical contributions correctly", {
  # Type 'a' has analytical score, type 'b' does not
  model <- likelihood_contr(
    obs_type = "type",
    a = contr_fn(
      loglik = function(df, par, ...) sum(dnorm(df$x, par[1], par[2], log = TRUE)),
      score  = function(df, par, ...) {
        x <- df$x
        n <- length(x)
        mu <- par[1]; sigma <- par[2]
        c(sum(x - mu) / sigma^2,
          -n / sigma + sum((x - mu)^2) / sigma^3)
      }
    ),
    b = contr_fn(
      loglik = function(df, par, ...) sum(dnorm(df$x, par[1], par[2], log = TRUE))
    )
  )

  df <- data.frame(x = c(1, 2, 3, 4), type = c("a", "a", "b", "b"))
  score_fn <- score(model)
  result <- score_fn(df, c(2.5, 1))

  # Should equal the full analytical score for all data
  x <- c(1, 2, 3, 4)
  n <- 4; mu <- 2.5; sigma <- 1
  expected <- c(sum(x - mu) / sigma^2,
                -n / sigma + sum((x - mu)^2) / sigma^3)
  expect_equal(result, expected, tolerance = 1e-4)
})
