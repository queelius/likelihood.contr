# ---- Constructor validation -------------------------------------------------

test_that("likelihood_contr validates obs_type", {
  ctr <- contr_fn(loglik = function(df, par, ...) 0)
  expect_error(likelihood_contr(obs_type = 42, a = ctr), "'obs_type' must be")
  expect_error(likelihood_contr(obs_type = NULL, a = ctr), "'obs_type' must be")
})

test_that("likelihood_contr validates obs_type string is non-empty", {
  ctr <- contr_fn(loglik = function(df, par, ...) 0)
  expect_error(likelihood_contr(obs_type = "", a = ctr), "non-empty")
  expect_error(likelihood_contr(obs_type = NA_character_, a = ctr), "non-empty")
})

test_that("likelihood_contr requires at least one contribution", {
  expect_error(likelihood_contr(obs_type = "type"), "At least one contribution")
})

test_that("likelihood_contr requires all contributions to be named", {
  ctr <- contr_fn(loglik = function(df, par, ...) 0)
  expect_error(likelihood_contr(obs_type = "type", ctr), "must be named")
})

test_that("likelihood_contr validates contributions are contr objects", {
  expect_error(
    likelihood_contr(obs_type = "type", a = "not_contr"),
    "not a 'contr' object"
  )
})

test_that("likelihood_contr always includes iid assumption", {
  ctr <- contr_fn(loglik = function(df, par, ...) 0)
  model <- likelihood_contr(obs_type = "type", a = ctr)
  expect_true("iid" %in% model$assumptions)
})

test_that("likelihood_contr includes user assumptions", {
  ctr <- contr_fn(loglik = function(df, par, ...) 0)
  model <- likelihood_contr(obs_type = "type", a = ctr,
                            assumptions = c("Weibull dist", "non-informative"))
  expect_true("Weibull dist" %in% model$assumptions)
  expect_true("non-informative" %in% model$assumptions)
  expect_true("iid" %in% model$assumptions)
})

test_that("likelihood_contr deduplicates iid in assumptions", {
  ctr <- contr_fn(loglik = function(df, par, ...) 0)
  model <- likelihood_contr(obs_type = "type", a = ctr,
                            assumptions = c("iid", "iid", "other"))
  expect_equal(sum(model$assumptions == "iid"), 1L)
})

test_that("likelihood_contr has correct class", {
  ctr <- contr_fn(loglik = function(df, par, ...) 0)
  model <- likelihood_contr(obs_type = "type", a = ctr)
  expect_s3_class(model, "likelihood_contr")
  expect_s3_class(model, "likelihood_model")
})

test_that("likelihood_contr validates rdata_fn", {
  ctr <- contr_fn(loglik = function(df, par, ...) 0)
  expect_error(
    likelihood_contr(obs_type = "type", a = ctr, rdata_fn = "bad"),
    "'rdata_fn' must be a function"
  )
})


# ---- loglik dispatch: column-based -----------------------------------------

test_that("loglik.likelihood_contr dispatches by column", {
  exact_ll <- function(df, par, ...) {
    sum(dexp(df$t, rate = par[1], log = TRUE))
  }
  right_ll <- function(df, par, ...) {
    sum(pexp(df$t, rate = par[1], lower.tail = FALSE, log.p = TRUE))
  }
  model <- likelihood_contr(
    obs_type = "status",
    exact = contr_fn(loglik = exact_ll),
    right = contr_fn(loglik = right_ll)
  )

  df <- data.frame(
    t      = c(1, 2, 3, 4),
    status = c("exact", "exact", "right", "right")
  )
  ll_fn <- loglik(model)
  result <- ll_fn(df, par = c(0.5))

  # Manual computation
  expected <- sum(dexp(c(1, 2), 0.5, log = TRUE)) +
              sum(pexp(c(3, 4), 0.5, lower.tail = FALSE, log.p = TRUE))
  expect_equal(result, expected)
})


# ---- loglik dispatch: function-based ---------------------------------------

test_that("loglik.likelihood_contr dispatches by function", {
  obs_fn <- function(df) ifelse(df$delta == 1, "exact", "right")
  exact_ll <- function(df, par, ...) {
    sum(dexp(df$t, rate = par[1], log = TRUE))
  }
  right_ll <- function(df, par, ...) {
    sum(pexp(df$t, rate = par[1], lower.tail = FALSE, log.p = TRUE))
  }
  model <- likelihood_contr(
    obs_type = obs_fn,
    exact = contr_fn(loglik = exact_ll),
    right = contr_fn(loglik = right_ll)
  )

  df <- data.frame(t = c(1, 2, 3), delta = c(1, 0, 1))
  ll_fn <- loglik(model)
  result <- ll_fn(df, par = c(0.5))

  expected <- sum(dexp(c(1, 3), 0.5, log = TRUE)) +
              sum(pexp(c(2), 0.5, lower.tail = FALSE, log.p = TRUE))
  expect_equal(result, expected)
})


# ---- loglik caching --------------------------------------------------------

test_that("loglik caches the data frame split", {
  call_count <- 0L
  obs_fn <- function(df) {
    call_count <<- call_count + 1L
    rep("a", nrow(df))
  }
  model <- likelihood_contr(
    obs_type = obs_fn,
    a = contr_fn(loglik = function(df, par, ...) -sum(par^2))
  )
  ll_fn <- loglik(model)
  df <- data.frame(x = 1:3)

  ll_fn(df, c(1))
  ll_fn(df, c(2))  # same df, should not re-split
  # obs_fn should have been called only once (on initial split)
  expect_equal(call_count, 1L)
})


# ---- Edge cases: missing column --------------------------------------------

test_that("loglik errors on missing dispatch column", {
  model <- likelihood_contr(
    obs_type = "status",
    a = contr_fn(loglik = function(df, par, ...) 0)
  )
  ll_fn <- loglik(model)
  expect_error(ll_fn(data.frame(x = 1), c(1)), "Column 'status' not found")
})


# ---- Edge cases: unmatched type --------------------------------------------

test_that("loglik errors when data has unmatched observation type", {
  model <- likelihood_contr(
    obs_type = "type",
    a = contr_fn(loglik = function(df, par, ...) 0)
  )
  ll_fn <- loglik(model)
  df <- data.frame(x = 1, type = "b")
  expect_error(ll_fn(df, c(1)), "No contribution defined.*'b'")
})


# ---- Edge cases: NA in types -----------------------------------------------

test_that("loglik errors on NA in observation types", {
  model <- likelihood_contr(
    obs_type = "type",
    a = contr_fn(loglik = function(df, par, ...) 0)
  )
  ll_fn <- loglik(model)
  df <- data.frame(x = 1:2, type = c("a", NA))
  expect_error(ll_fn(df, c(1)), "NA values")
})


# ---- Edge cases: empty subset ----------------------------------------------

test_that("loglik handles empty subsets gracefully", {
  model <- likelihood_contr(
    obs_type = "type",
    a = contr_fn(loglik = function(df, par, ...) sum(dnorm(df$x, par[1], par[2], log = TRUE))),
    b = contr_fn(loglik = function(df, par, ...) {
      stop("should not be called for empty subset")
    })
  )
  # All data is type 'a', type 'b' is empty.
  # split() with factor only gives levels present, so we won't get an
  # empty 'b' subset from split. This test verifies no error occurs.
  df <- data.frame(x = c(0, 1), type = c("a", "a"))
  ll_fn <- loglik(model)
  result <- ll_fn(df, c(0, 1))
  expected <- sum(dnorm(c(0, 1), 0, 1, log = TRUE))
  expect_equal(result, expected)
})


# ---- assumptions -----------------------------------------------------------

test_that("assumptions.likelihood_contr returns stored assumptions", {
  ctr <- contr_fn(loglik = function(df, par, ...) 0)
  model <- likelihood_contr(obs_type = "type", a = ctr,
                            assumptions = c("Weibull"))
  expect_equal(assumptions(model), c("iid", "Weibull"))
})


# ---- rdata -----------------------------------------------------------------

test_that("rdata.likelihood_contr returns rdata_fn when present", {
  rdf <- function(theta, n, ...) data.frame(x = rnorm(n, theta[1], theta[2]))
  ctr <- contr_fn(loglik = function(df, par, ...) 0)
  model <- likelihood_contr(obs_type = "type", a = ctr, rdata_fn = rdf)
  fn <- rdata(model)
  expect_identical(fn, rdf)
  result <- fn(c(0, 1), 5)
  expect_equal(nrow(result), 5L)
})

test_that("rdata.likelihood_contr errors when no rdata_fn", {
  ctr <- contr_fn(loglik = function(df, par, ...) 0)
  model <- likelihood_contr(obs_type = "type", a = ctr)
  expect_error(rdata(model), "does not have an rdata function")
})


# ---- print method ----------------------------------------------------------

test_that("print.likelihood_contr prints expected output", {
  ctr1 <- contr_fn(loglik = function(df, par, ...) 0)
  ctr2 <- contr_fn(loglik = function(df, par, ...) 0)
  model <- likelihood_contr(
    obs_type = "type",
    exact = ctr1,
    right = ctr2,
    assumptions = c("Weibull dist")
  )
  output <- capture.output(print(model))
  expect_true(any(grepl("Likelihood Contribution Model", output)))
  expect_true(any(grepl("exact, right", output)))
  expect_true(any(grepl("column", output)))
  expect_true(any(grepl("iid", output)))
  expect_true(any(grepl("Weibull dist", output)))
})

test_that("print.likelihood_contr shows function dispatch", {
  ctr <- contr_fn(loglik = function(df, par, ...) 0)
  model <- likelihood_contr(
    obs_type = function(df) rep("a", nrow(df)),
    a = ctr
  )
  output <- capture.output(print(model))
  expect_true(any(grepl("function", output)))
})
