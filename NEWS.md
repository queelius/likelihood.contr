# likelihood.contr 0.2.0

* Document usage of re-exported generics, resolving the Rd-contents note on current CRAN R-devel checks.
* Reject duplicate contribution names and malformed observation callbacks before they can silently recycle or misroute observations.
* Add `contr_freq()` for counted independent observations and histogram bins, with analytical derivative preservation and numerical fallback.
* Add a reproducible timeout-aware latency benchmarking vignette and demo, with explicit censoring and extrapolation assumptions.

# likelihood.contr 0.1.1

* Improved documentation: enriched `\value` tags on all exported Rd files
  per CRAN reviewer feedback.
* Added `\arguments` and `\value` to re-exported generics page (`loglik.Rd`).
* Removed single quotes around function names in DESCRIPTION.

# likelihood.contr 0.1.0

* Initial release.
* `contr_fn()` and `contr_name()` for building likelihood contributions.
* `likelihood_contr()` composes contributions into a `likelihood_model`.
* Support for exact, right-censored, left-censored, and interval-censored data.
* Automatic numerical differentiation fallback via `numDeriv`.
