## Submission

This is a new submission. The package implements the methods of Liu (2026),
"Cohort-Anchored HonestDiD" <doi:10.48550/arXiv.2509.01829>.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

  Any words flagged as possibly misspelled in DESCRIPTION are author names
  (Callaway, Liu, Rambachan, Roth, Sant'Anna) or 'HonestDiD', the name of the
  method in the title of the paper that the package implements.

## Test environments

* Local: Windows 11 x64, R 4.4.0, `R CMD check --as-cran`
* win-builder: R-devel and R-release

## Code adapted from 'HonestDiD'

The R files in `inst/upstream/R` contain routines adapted from the 'HonestDiD'
package (MIT license). Its authors, Ashesh Rambachan and Jonathan Roth, are
listed in `Authors@R` with the roles "ctb" and "cph". The original license
and a notice that describes the changes are kept in `inst/upstream`.

'HonestDiD' does not export these routines, and this package modifies them
(positive semidefinite covariances, scalar targets, scale-relative tolerances,
and restoring the user's random-number state), so it carries adapted copies
instead of importing 'HonestDiD'. The copies are loaded into a private
environment on first use. Nothing is attached to the search path, and nothing
is written to the user's file system.
