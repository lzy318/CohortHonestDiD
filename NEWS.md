# CohortHonestDiD 0.2.0

* Renamed the package from HonestDiDbyCohort to CohortHonestDiD, the name used
  in the accompanying paper, Liu (2026), "Cohort-Anchored HonestDiD". Use
  `library(CohortHonestDiD)`; code that referred to `HonestDiDbyCohort::` must
  use `CohortHonestDiD::`.
* New default procedure, following the paper's recommendation:
  `cohort_sensitivity()` now uses C-LF under SD as well as under RM, and under
  SD the default scope is the post-only restriction (`full.path = FALSE`). In
  HonestDiDbyCohort 0.1.1 the SD default was the full-path FLCI, so calls that
  omit `method` under SD now return different confidence sets. Use
  `method = "FLCI"` for the previous default; its `full.path` default remains
  TRUE.
* New option `cohort_fit(covariance = "CR2")`, the bias-reduced CR2 covariance
  of Bell and McCaffrey (2002). It drops the CR1 factor and divides each
  cluster's score by sqrt(1 - p), where p is the cluster's share of the units
  of its cohort (of the unit weights when weights are supplied), the
  never-treated units forming their own group. A cluster that makes up a whole
  cohort keeps a zero score, and a cluster with units in several cohorts has
  each cohort's part of its score adjusted with its share of that cohort. CR2
  fits record the shares and factors in `cluster_shares`. The default remains
  `covariance = "CR1"`.
* New arguments `balanced` and `cohort.paths` of `cohort_sensitivity()`. With
  `balanced = TRUE`, the target uses only the post-treatment relative periods
  observed for every cohort, with cohort weights fixed across periods, and
  aggregate inference forms the aggregated event study only on the relative
  periods observed for every cohort (balanced aggregation). Under cohort
  inference, `cohort.paths = "full"` (the default) keeps each cohort's whole
  path, and `cohort.paths = "window"` applies the restriction only inside the
  balanced window. `settings` records the window, and balanced aggregate
  results contain the aggregated coefficients in `aggregation`. The default,
  `balanced = FALSE`, leaves every result unchanged.
* `print()` and `summary()` of fits and sensitivity results now name the
  covariance: CR1, CR2, or supplied to `cohort_input()`.
* Calls that state `method` explicitly and use the default CR1 covariance give
  the same numerical results as in HonestDiDbyCohort 0.1.1. Dependencies are
  unchanged.
* Added `citation("CohortHonestDiD")`, repository links in DESCRIPTION, a
  README for the GitHub repository, and a source reference on `?minimum_wage`.
* The detailed guide is now built as a second vignette,
  `vignette("detailed-guide", package = "CohortHonestDiD")`.

# HonestDiDbyCohort 0.1.1

* For Conditional initial search radii above 8, automatically increase the
  initial grid count to preserve the standard radius-8 reference density.
  Requested endpoints stay unchanged; integer rounding gives a spacing upper
  bound rather than exact nesting across radius requests. Radius 8 and smaller
  keep their previous computations. Report effective initial counts and
  reference density; wider searches may reach the existing evaluation cap.
* Search progressively wider bilateral domains for Conditional inference,
  keeping mesh spacing fixed within each refinement and requiring two
  successive adjacent-domain agreements plus endpoint checks. The initial
  contrast and sampling scale respect coefficient normalization. Share one
  cumulative point budget across all domains, audits, and refinements.
  `max.expansions` includes required stability expansions; values below 2
  report `domain_stability_limit` with no set. Print, summary, plot captions,
  and detailed diagnostics disclose searched ranges and the lack of a global
  guarantee. Disabling topology audits retains domain and endpoint checks.
  Interior audits also check a combined grid of regular-cell quarter and
  three-quarter phases. Disagreement triggers the existing bounded refinement
  policy; finite mesh checks do not certify every arbitrarily narrow gap.
* Validate numerical and solver controls before inference, including integer
  ranges, scalar values, and unique control names.
* Preserve bias-map and target weights under extreme finite cohort-mass scales.
* Make external covariance and normalization checks robust to large scales;
  normalization row labels no longer affect inference.
* Support the exact zero-covariance FLCI limit without regularization, and keep
  positive-noise FLCI results invariant to changes in outcome units.
  Compare curvature in finite covariance units to prevent construction-bound
  overflow, and preserve exact covariance-null adjustment directions.
  Guard exact-null shortcuts against entry and product underflow; average
  symmetric matrices without overflowing finite large entries or erasing
  equal subnormal entries.
  Scale FLCI curvature and M=0 variance-minimizing weights using only
  contributing covariance blocks, preserving small positive pre variances
  beside large unrelated post variances.
* Report unsupported deterministic Conditional/C-LF cases with an explicit
  failure status and no confidence set, and clarify insufficient prehistory.
  The public failure warning identifies how many tasks have unsupported
  deterministic or numerically degenerate variances and reports that they
  return no confidence set.
* Rescale extreme panel weights when required for finite arithmetic, recording
  the common divisor in `data_info$weight_scale`; preserve exact zero scores
  when a balanced cohort/time cell lies within one cluster.
* Preserve the regular search spacing when profile anchors nearly coincide
  with grid points, so adaptive confidence-set searches can expand as intended.
* Retry eligible unresolved searches with finer grids, up to
  `control$max.refinements` (default 3), sharing one point-evaluation budget and
  the same inference setup. Record all attempts and retain explicit failure
  with no confidence set when the budget or retry allowance is exhausted.
* Resolve private numerical resources from the loaded package directory in
  both installed packages and development source trees.
* Bound SD C-LF searches using certified LP extrema of their own first-stage
  acceptance region, preserving the point test and calibration. Distinguish
  certified empty sets, finite outer support, and unresolved support programs.

# HonestDiDbyCohort 0.1.0

* Initial preliminary release with imputation, CS-NYT, and CS-NT estimation.
* RM Conditional/C-LF and SD Conditional/C-LF/FLCI, with explicit SD scope.
* ATT, dynamic, and cohort-dynamic targets; optional aggregate benchmarks.
* Reproducible parallel inference and the paper's plotting palette.
* Unbalanced-panel approximation with treatment-status carry-forward diagnostics.
