test_that("observation dispatch cannot silently recycle", {
  ctr <- contr_name("exp", "exact")
  d <- data.frame(x = c(1, 2))
  for (callback in list(function(df) "exact", function(df) character(),
                       function(df) c(1, 1))) {
    model <- likelihood_contr(callback, exact = ctr)
    expect_error(loglik(model)(d, 1), "character vector.*row")
  }
  model <- likelihood_contr(function(df) rep("exact", nrow(df)), exact = ctr)
  expect_equal(loglik(model)(d, 1), -3)
  expect_error(loglik(model)(NULL, 1), "data frame")
  d$status <- factor(c("exact", "exact"))
  expect_equal(loglik(likelihood_contr("status", exact = ctr))(d, 1), -3)
})

test_that("contribution identities must be unique and complete", {
  ctr <- contr_name("exp", "exact")
  expect_error(likelihood_contr("status", exact = ctr, exact = ctr), "unique")
  args <- list("status", ctr)
  names(args) <- c("obs_type", "")
  expect_error(do.call(likelihood_contr, args), "named")
})
