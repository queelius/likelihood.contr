# CRAN Submission Comments - likelihood.contr 0.1.1

## Resubmission

This is a resubmission. In response to Benjamin Altmann's review (2026-03-29):

- Removed single quotes around function names in the DESCRIPTION field.
- Added `\arguments` and `\value` tags to the re-exported generics page
  (`loglik.Rd`), with descriptions of return types for each generic.
- Enriched `\value` descriptions on all other exported `.Rd` files to
  explain the structure (class) and meaning of return values.

## Test Environments

- Ubuntu 24.04 (local), R 4.3.3
- GitHub Actions: ubuntu-latest (R release, devel, oldrel-1),
  windows-latest (R release), macos-latest (R release)
- win-builder (R-devel)

## R CMD check Results

0 errors | 0 warnings | 0 notes
