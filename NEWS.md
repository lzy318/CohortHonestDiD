# CohortHonestDiD 1.1.1

First CRAN release. Results are unchanged from 1.1.0.

* The help pages of `cohort_input()`, `paper_palette()`, and the plot methods
  have runnable examples.
* `DESCRIPTION` lists Ashesh Rambachan and Jonathan Roth as copyright holders
  of the routines adapted from HonestDiD.

# CohortHonestDiD 1.1.0

## New restrictions

* `cohort_sensitivity()` accepts seven more restrictions, `"SDB"`, `"SDM"`,
  `"RMB"`, `"RMM"`, `"SDRM"`, `"SDRMB"`, and `"SDRMM"`. They work with the
  Conditional and C-LF procedures, and with cohort and aggregate inference. C-LF
  is the default. The FLCI stays available under SD only.
* The new argument `direction` sets the sign (`"positive"` or `"negative"`) or
  the monotonicity (`"increasing"` or `"decreasing"`) for every cohort. The sign
  is measured from relative period 0 for every estimator.
* `full.path = TRUE` works under SD, SDB, and SDM. The other six restrictions do
  not take it.
* SDRM, SDRMB, and SDRMM need three pre-treatment coefficients per cohort,
  counting period 0.
* `print()`, `summary()`, and `plot()` name the restriction and its direction.

The paper's simulations do not cover these seven restrictions. The detailed guide
has a section on them.

## Bug fix

* `control = list(flci.solver = "nested")` works with current CVXR. CVXR 1.8 and
  later no longer export `solve()`, so the task used to end with status
  `flci_solver_error`. The default solver `"ECOS"` needs the ECOSolveR package,
  and `flci.control = list(solver = "CLARABEL")` works without it.

## Other changes

* Results that ran under 1.0.0 are unchanged.
* `restriction = "S"` and `restriction = "R"` no longer work as abbreviations,
  because they are now ambiguous. Write `"SD"` or `"RM"`.

# CohortHonestDiD 1.0.0

First public release. It estimates cohort-period effects with the imputation and
Callaway-Sant'Anna estimators, and it applies the SD and RM restrictions under
cohort or aggregate inference with the Conditional, C-LF, and FLCI procedures. It
also provides the CR2 covariance, balanced targets and balanced aggregation, and
`flci_weights()`.
