## likelihood.contr 0.2.0

Adds frequency-weighted additive likelihood contributions for counted
observations and histogram bins, with analytical derivative preservation and
numerical fallback. Documentation explicitly distinguishes frequencies from
survey weights and warns that generic row-based nobs/BIC/bootstrap semantics
are inappropriate for frequency data.

Adds a synthetic timeout-aware service/AI latency demonstration. The example
checks the exponential MLE against its closed form and explicitly distinguishes
observed SLO rates from model-dependent extrapolation beyond timeouts.

Also validates observation callback lengths/types and duplicate contribution
names. Adds usage documentation for re-exported generics, addressing the
Rd-contents NOTE on the existing CRAN R-devel Debian check pages.

## Checks

Ubuntu 24.04.3 LTS, R 4.3.3, x86_64 Linux, 2026-10-28:
Final R CMD check --as-cran: 0 errors, 0 warnings, 2 environment notes.

* Unable to verify current time.
* HTML validation skipped because tidy is unavailable.

112 unit/integration assertions pass, with 100% measured line coverage.
Examples, vignette rebuilding and PDF manual pass. The initial remote incoming
checks completed; remote incoming checks were disabled only on the final
repeat. Windows R-devel validation has been requested separately.
