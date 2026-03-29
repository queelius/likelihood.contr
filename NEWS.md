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
