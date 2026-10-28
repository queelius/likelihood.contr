# Synthetic benchmark. No provider ranking or empirical performance claim.
library(likelihood.contr)
library(likelihood.model)

local({
  set.seed(20261028)
  deadline <- 2
  true_rate <- 1
  latent <- rexp(1000, rate = true_rate)
  requests <- data.frame(seconds = pmin(latent, deadline),
                         status = ifelse(latent <= deadline, "exact", "right"))
  model <- likelihood_contr("status",
    exact = contr_name("exp", "exact", "seconds"),
    right = contr_name("exp", "right", "seconds"))
  # Closed-form start also supplies an independent correctness check.
  rate <- sum(requests$status == "exact") / sum(requests$seconds)
  fitted <- fit(model)(requests, c(rate = rate))
  stopifnot(fitted$converged, abs(unname(coef(fitted)) - rate) < 1e-5)
  completed <- requests$seconds[requests$status == "exact"]
  report <- data.frame(
    quantity = c("timeout fraction", "completed-only empirical p95 (seconds)",
                 "exponential model p95 (seconds)", "simulation true p95 (seconds)"),
    value = c(mean(requests$status == "right"), unname(quantile(completed, .95)),
              qexp(.95, rate = coef(fitted)), qexp(.95, rate = true_rate)))
  print(report, row.names = FALSE)
  cat("\nThe model p95 extrapolates beyond the deadline; it depends on the exponential assumption.\n")

  # Counted disjoint bins avoid expanding telemetry into individual requests.
  edges <- c(0, .5, 1, 1.5, deadline, Inf)
  counts <- tabulate(findInterval(requests$seconds, edges, left.open = TRUE),
                     nbins = length(edges) - 1L)
  # Timed-out requests were recorded at the deadline; move them to the open tail.
  timed_out <- sum(requests$status == "right")
  counts[length(counts) - 1L] <- counts[length(counts) - 1L] - timed_out
  counts[length(counts)] <- timed_out
  bins <- data.frame(lo = head(edges, -1), hi = tail(edges, -1),
                     count = counts, status = "interval")
  grouped <- likelihood_contr("status",
    interval = contr_freq(contr_name("exp", "interval", "lo", "hi"), "count"))
  grouped_ll <- loglik(grouped)
  # Optimize log(rate) to respect the parameter domain explicitly.
  solution <- optimize(function(log_rate) grouped_ll(bins, c(rate = exp(log_rate))),
                       c(-5, 5), maximum = TRUE)
  print(data.frame(grouped_rate = exp(solution$maximum), request_rate = rate))
})
