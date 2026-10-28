test_that("frequency contributions equal expanded observations", {
  d <- data.frame(x = c(1, 2, 3), n = c(2, 0, 3))
  expanded <- d[rep(seq_len(nrow(d)), d$n), , drop = FALSE]
  ctr <- contr_fn(
    function(df, par, ...) sum(dexp(df$x, par[1], log = TRUE)),
    function(df, par, ...) c(rate = nrow(df) / par[1] - sum(df$x)),
    function(df, par, ...) matrix(-nrow(df) / par[1]^2, 1)
  )
  freq <- contr_freq(ctr, "n")
  expect_s3_class(freq, "contr")
  expect_equal(freq$loglik(d, 2), ctr$loglik(expanded, 2))
  expect_equal(freq$score(d, 2), ctr$score(expanded, 2))
  expect_equal(freq$hess(d, 2), ctr$hess(expanded, 2))
  expect_equal(freq$loglik(d[FALSE, ], 2), 0)
  expect_equal(freq$score(d[FALSE, ], c(rate = 2)), c(rate = 0))
  expect_equal(freq$hess(d[FALSE, ], 2), matrix(0, 1))
})

test_that("frequency weighting composes with mixed censoring and numeric derivatives", {
  model <- likelihood_contr("kind",
    exact = contr_freq(contr_name("exp", "exact"), "n"),
    right = contr_freq(contr_name("exp", "right"), "n"))
  d <- data.frame(x = c(1, 2, 3), n = c(2, 3, 4), kind = c("exact", "exact", "right"))
  expect_equal(loglik(model)(d, c(rate = 0.5)), 5 * log(0.5) - 0.5 * 20)
  expect_equal(unname(score(model)(d, c(rate = 0.5))), 5 / 0.5 - 20, tolerance = 1e-6)
  expect_equal(hess_loglik(model)(d, c(rate = 0.5)), matrix(-5 / 0.5^2, 1), tolerance = 1e-5)
})

test_that("invalid counts are rejected and zero-frequency rows are not evaluated", {
  ctr <- contr_name("exp", "exact")
  expect_error(contr_freq(list(), "n"), "contr")
  for (col in list(NA_character_, "", c("n", "m"), 1))
    expect_error(contr_freq(ctr, col), "column")
  freq <- contr_freq(ctr, "n")
  expect_error(freq$loglik(data.frame(x = 1), 1), "count")
  for (n in list(-1, 0.5, NA_real_, Inf, "2", TRUE))
    expect_error(freq$loglik(data.frame(x = 1, n = n), 1), "count")
  expect_equal(freq$loglik(data.frame(x = -1, n = 0), 1), 0)
  expect_null(freq$score)
  expect_null(freq$hess)
})

test_that("interval frequencies match a multinomial log likelihood up to a constant", {
  d <- data.frame(lo = c(0, 1, 2), hi = c(1, 2, Inf), n = c(12, 5, 3))
  ctr <- contr_freq(contr_name("exp", "interval", "lo", "hi"), "n")
  probs <- diff(pexp(c(0, 1, 2, Inf), rate = 0.7))
  expect_equal(ctr$loglik(d, c(rate = 0.7)), sum(d$n * log(probs)))
})
