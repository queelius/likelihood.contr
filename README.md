# likelihood.contr

Compose likelihood models from heterogeneous observation types.

Many real datasets mix exact observations with censored ones. A sensor
might record a precise failure time or only that a unit survived past
inspection. Each observation type contributes differently to the
log-likelihood, but they share the same parameters. `likelihood.contr`
lets you declare each type's contribution separately, then compose them
into a single model with full inference support: MLE, standard errors,
confidence intervals, likelihood ratio tests, and Fisher information.

## Installation

```r
# Stable release from CRAN
install.packages("likelihood.contr")

# Development builds (with CRAN dependencies)
install.packages("likelihood.contr", repos = c(
  "https://queelius.r-universe.dev", "https://cloud.r-project.org"))

# Or from GitHub
remotes::install_github("queelius/likelihood.contr")
```

## Start with a decision: is this service fast enough?

A timed-out request is evidence that completion took longer than the deadline.
Dropping it can make a slow service look fast. Combine completed requests and
timeouts in one likelihood, and keep measured SLO rates separate from
model-dependent tail extrapolations.

The [timeout benchmarking walkthrough](vignettes/timeout-benchmarking.Rmd)
includes a reproducible service/AI example, an independently checked exponential
MLE, and histogram inference without expanding counts into individual requests.

```r
# Available in version 0.2.0
demo("timeout-benchmark", package = "likelihood.contr")
```

`contr_freq()` weights existing contributions by observed counts, so the same
model machinery handles raw measurements and disjoint telemetry bins. It is
for independent replicated observations, not survey weights. Standard survival
workflows are also well served by `survival` and `fitdistrplus`; use this library
when observation mechanisms need explicit custom composition.

## Quick Example

Fit a Weibull model to a mix of exact and right-censored failure times:

```r
library(likelihood.contr)
library(likelihood.model)

# Define model: Weibull with two observation types
model <- likelihood_contr(
  obs_type = "status",
  exact = contr_name("weibull", "exact", ob_col = "t"),
  right = contr_name("weibull", "right", ob_col = "t")
)

# Simulate data with 50% censoring
set.seed(42)
raw <- rweibull(300, shape = 2, scale = 5)
df <- data.frame(
  t      = pmin(raw, 4),
  status = ifelse(raw <= 4, "exact", "right")
)

# Fit and inspect
result <- fit(model)(df, par = c(shape = 1.5, scale = 4))
summary(result)
#> Coefficients:
#>       Estimate Std. Error   2.5% 97.5%
#> shape   2.0672     0.1611 1.7514 2.383
#> scale   4.9571     0.2309 4.5046 5.410
```

Three lines to define the model. One line to fit it. Full asymptotic
inference for free.

## How It Works

`likelihood.contr` provides two contribution constructors:

- **`contr_name()`** builds contributions from any R distribution
  (`"norm"`, `"exp"`, `"weibull"`, ...) for exact, right-censored,
  left-censored, or interval-censored data.
- **`contr_fn()`** wraps user-supplied log-likelihood, score, and
  Hessian functions for custom observation types.

`likelihood_contr()` composes these into a model that implements the
[likelihood.model](https://github.com/queelius/likelihood.model)
interface: `loglik()`, `score()`, `hess_loglik()`, `fit()`, `lrt()`,
`fim()`, and more.

## Vignettes

- `vignette("likelihood-contr-introduction")` -- Core workflow: build,
  fit, interpret
- `vignette("custom-contributions")` -- Analytical derivatives and
  model comparison
- `vignette("masked-series-system")` -- Competing risks with masked
  failure causes
