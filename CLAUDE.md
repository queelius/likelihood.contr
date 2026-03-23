# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Package Does

`likelihood.contr` composes likelihood models from heterogeneous observation types (exact, left-censored, right-censored, interval-censored, custom). Each observation type contributes independently to the total log-likelihood under an i.i.d. assumption. The composed model implements the `likelihood.model` interface (`loglik`, `score`, `hess_loglik`, `fim`, `fit`, `lrt`).

## Development Commands

```r
devtools::load_all()          # Load package for development
devtools::test()              # Run all tests
devtools::test(filter = "contr")  # Run single test file
devtools::document()          # Regenerate NAMESPACE + man/
devtools::check()             # Full R CMD check
covr::package_coverage()      # Test coverage
```

## Architecture

Three source files, two S3 classes:

- **`R/contr.R`**: Contribution constructors. `contr_fn()` wraps user functions into `"contr"` objects. `contr_name()` auto-generates contributions from named R distributions (`d<name>`/`p<name>`) for the four censoring types. `.prepare_args_list()` maps flat parameter vectors to distribution formals by introspecting `formals(fn)`.

- **`R/likelihood_contr.R`**: Model composition and `likelihood_model` S3 methods. `likelihood_contr()` constructs the model; `loglik/score/hess_loglik` methods return **closures** `function(df, par, ...)` that split data by observation type and sum per-type contributions. Score/Hessian fall back to `numDeriv` when analytical functions aren't provided.

- **`R/cache.R`**: Split caching via R environments. `identical(df, cache$df_ref)` skips redundant `split()` calls during optimizer iterations.

**Class hierarchy:** `"likelihood_contr"` extends `"likelihood_model"`. All `likelihood.model` generics (`fit`, `lrt`, `fim`, `fisher_mle`) work on these objects.

**Closure-returning pattern:** `loglik(model)` returns a *function*, not a value. This is the convention across the entire rlang ecosystem. The returned closure closes over cached state and the model structure.

## Dependencies

- **`likelihood.model`** (Imports): provides the generics this package implements
- **`numDeriv`** (Imports): numerical score/Hessian fallback when contributions lack analytical derivatives

## Key Design Decisions

- `obs_type` can be a column name (string) or a function, enabling both column-dispatch and computed-dispatch patterns
- `contr_name()` contributions have `score = NULL` and `hess = NULL`, so they always use numerical differentiation. Only `contr_fn()` with explicitly supplied functions gets analytical derivatives
- Parameter vectors are flat numerics mapped to distribution formals by position (unnamed) or name (named)
