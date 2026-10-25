# CohortHonestDiD

CohortHonestDiD implements cohort-anchored HonestDiD sensitivity analysis for
difference-in-differences designs with staggered adoption. Instead of
restricting event-study coefficients aggregated across cohorts, it places the
sensitivity restrictions of Rambachan and Roth (2023) on each cohort's block
bias: the cohort's untreated trend relative to the units that are untreated
when it first adopts treatment. A known linear map carries block biases into
the biases of the estimated coefficients, so the HonestDiD procedures apply
directly. The package estimates cohort-period effects with the imputation
estimator and the Callaway-Sant'Anna estimators (not-yet-treated or
never-treated controls) and reports confidence sets for average, dynamic, or
cohort-specific dynamic effects. Inference on aggregated coefficients is
available as a benchmark.

The package was previously distributed as HonestDiDbyCohort (versions 0.1.0 and
0.1.1). Version 0.2.0 follows the paper's recommendations: under SD the default
procedure is now C-LF with the post-only restriction (previously the full-path
FLCI), and `cohort_fit(covariance = "CR2")` adds the bias-reduced CR2
covariance; see [NEWS](NEWS.md).

## The paper

Liu, Ziyi (2026). "Cohort-Anchored HonestDiD." Working paper.

The paper develops the method, compares cohort and aggregate inference in
simulations, and applies it to the minimum-wage study of Callaway and
Sant'Anna (2021) and the community-health-center study of Bailey and
Goodman-Bacon (2015). The bundled `minimum_wage` data are the sample of the
first application.

## Installation

CohortHonestDiD is not on CRAN. Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("YOUR-GITHUB-USER/CohortHonestDiD")
library(CohortHonestDiD)
```

To also build the two vignettes (this takes a few minutes and installs the
suggested packages, including knitr and rmarkdown):

```r
remotes::install_github("YOUR-GITHUB-USER/CohortHonestDiD",
                        dependencies = TRUE, build_vignettes = TRUE)
```

The package requires R 4.1.0 or later. The linear-programming dependency
Rglpk needs the GLPK library; on Debian or Ubuntu, install `libglpk-dev`
first. The optional `"nested"` FLCI solver uses CVXR.

## A first analysis

```r
data(minimum_wage)
fit <- cohort_fit(Y ~ D, data = minimum_wage,
                  index = c("id", "year"), estimator = "imp")
print(fit)
plot(fit)

# SD defaults to C-LF under the post-only restriction.
att <- cohort_sensitivity(fit, restriction = "SD", M = c(0, .005, .01, .015, .02),
                           target = "att", seed = 123)
plot(att)

# Overlay robust sets on the conventional event study.
dynamic <- cohort_sensitivity(fit, restriction = "SD", M = .01,
                               target = "dynamic", seed = 123)
plot(dynamic)

# Separate cohort-period effects.
by_cohort <- cohort_sensitivity(fit, restriction = "SD", M = .01,
                                 target = "cohort-dynamic", seed = 123)
plot(by_cohort)
```

The equivalent column-name interface is
`cohort_fit(data = minimum_wage, Y = "Y", D = "D", index = c("id", "year"))`.
Estimator options are `"imp"`, `"cs_nyt"`, and `"cs_nt"`. Supply a clustering
column with `cluster = "state"`; the default is the unit identifier.

## Methods and defaults

| Restriction | Methods | Default | Default `full.path` |
|---|---|---|---|
| SD | FLCI, Conditional, C-LF | C-LF | FALSE (post-only) for C-LF and Conditional; TRUE for FLCI |
| RM | Conditional, C-LF | C-LF | Not applicable |

Following the paper, the default is cohort inference with the C-LF procedure,
under the post-only restriction for SD. The FLCI requires the stronger
full-path restriction, and its coverage deteriorates when some cohorts contain
only a few units; use it with caution.

```r
conditional <- cohort_sensitivity(fit, M = .01, method = "Conditional", seed = 123)
# The full-path FLCI was the SD default before version 0.2.0.
flci <- cohort_sensitivity(fit, M = .01, method = "FLCI")
flci_post <- cohort_sensitivity(fit, M = .01, method = "FLCI", full.path = FALSE)
aggregate <- cohort_sensitivity(fit, M = c(0, .005, .01, .015, .02),
                                 framework = "aggregate", seed = 123)

# RM can be computationally expensive: this design has 320 signed components.
rm <- cohort_sensitivity(fit, restriction = "RM", M = c(0, .5, 1),
                          parallel = TRUE, cores = 2, seed = 123)
```

`framework = "cohort"` is the default. Aggregate inference supports `"att"`
and `"dynamic"`, with the same target weights; it does not support
`"cohort-dynamic"`. SD bounds are in outcome units, while RM bounds are
dimensionless. Relative period **1 is the first treated period**, and 0 is the
last untreated period. Dynamic and cohort-dynamic intervals are pointwise.

`balanced = TRUE` restricts the target, and under aggregate inference the
aggregated event study, to the relative periods observed for every cohort; see
`?cohort_sensitivity`.

## Covariance: CR1 and CR2

`cohort_fit()` estimates the joint covariance of all cohort-period
coefficients from cluster scores of outcomes centered at their cohort-period
cell means. The default, `covariance = "CR1"`, multiplies the sum of the
scores' outer products by the usual finite-sample factor.
`covariance = "CR2"` uses the bias-reduced adjustment of Bell and McCaffrey
(2002) instead: each cluster's score is divided by `sqrt(1 - p)`, where `p` is
the cluster's share of the units of its cohort (of the unit weights when
`weights` is supplied), and the never-treated units form their own group. On a
balanced panel, a cluster that makes up a whole cohort has a zero score under
both options; under CR2 it has a zero score on any panel, while under CR1 it
generally does not when that cohort has missing outcomes. The paper recommends
pairing Conditional or C-LF inference with CR2 when some cohorts contain only a
few units or clusters. `cohort_input()` uses the covariance it is given.

```r
fit_cr2 <- cohort_fit(Y ~ D, data = minimum_wage, index = c("id", "year"),
                      cluster = "state", covariance = "CR2")
print(fit_cr2)  # the header names the covariance
head(fit_cr2$cluster_shares)
clf_cr2 <- cohort_sensitivity(fit_cr2, M = c(0, .01), seed = 123)
```

## Parallel computation and reproducibility

Set `parallel = TRUE`. If `cores` is omitted, the default is the detected
logical core count minus one, with a minimum of one and no more workers than
available tasks. Examples use two workers. Independent target/bound tasks run
in parallel; for a single RM task, component calculations can run in parallel.
An explicit `seed` gives the same results across worker counts and preserves
the caller's RNG state. Avoid requesting more workers than available memory.

## Data handling and diagnostics

The package supports absorbing binary treatment without covariates,
periods interpreted as equally spaced, one-way clustering, and optional positive unit
weights that remain fixed over time. Both treated and never-treated units are
required, with enough pre-treatment support to construct the restrictions.
Irregularly spaced time labels trigger a warning about this interpretation.
Sensitivity analysis needs at least two pre-treatment periods for each included
cohort/path; an otherwise estimable one-pre-period fit has status
`insufficient_prehistory` at the inference step.

Extreme weights are divided by a common positive factor only when needed for
finite arithmetic. `fit$data_info$weight_scale` records that divisor, which is
1 for ordinary weights. Stored working weights and cohort masses then use the
rescaled units; multiply by the divisor to recover the original units.

Missing treatment status is carried forward from the most recent observed
status within unit, including absent unit-time rows. Leading unknown statuses
are not backfilled. Outcomes are never imputed. Thus, a treatment start hidden
inside a missing stretch is dated at the first observed treated period under
this convention. Inspect the treatment-fill and observation-pattern diagnostics.

Unbalanced panels proceed using actual observed histories and the **balanced-panel
bias map as an approximation**, with a warning. The balanced-panel coverage and
Imp/CS-NYT equivalence results do not automatically extend to this approximation.
Unavailable required cohort-period comparisons produce an explanatory error.

`print(fit)` and `summary(fit)` name the covariance (CR1 or CR2) and report
cohort/cluster support and covariance diagnostics.
Few treated clusters can make covariance estimates unreliable despite numerical
success. The covariance uses the paper's common cohort-period-centered score
construction; it is not a generic fixed-effects regression covariance.
Singular covariance alone does not prevent sensitivity analysis. If a target
or required tested moment is deterministic and the Conditional/C-LF procedure
cannot handle it, the package warns, sets `success = FALSE`, and reports
`unsupported_zero_variance` or `unsupported_moment_variance` with `NA` bounds.
FLCI supports its zero-noise limit, whose half-length is the worst-case bias.

## Inspect and customize results

`result$intervals` retains each connected confidence-set component and its
status. Empty, unbounded, disconnected, and failed/unresolved results remain
distinct. A finite search ending at an accepted boundary is not reported as a
complete finite confidence set. Plot methods return editable `ggplot` objects
and use the paper's colors. For dynamic results with several bounds,
`plot(result, M = .01)` selects a computed bound.

For eligible unresolved Conditional/C-LF searches, the package retries a finer
grid up to `control$max.refinements` times (default 3; 0 disables retries).
Every attempt shares the target/bound's `max.evaluations` budget for distinct
confidence-set point evaluations and keeps the same calibration draws and
numerical checks. Calibration and support-optimization work are outside this
counter; it is not an elapsed-time limit. If the budget or refinement limit
is exhausted, the result warns and reports no confidence set. Detailed results
record the attempted grids, statuses and evaluation counts.

SD C-LF uses its own first-stage support bounds to cover the outer search
domain, with the same simulated critical value used by its point tests.
Certified outer bounds do not certify every component or gap inside them:
nonempty finite-grid results retain `topology.certified = FALSE` and
`inversion.exact = FALSE`. Plain Conditional does not use this C-LF gate;
its finite-search diagnostics do not establish a global outer bound.

Conditional inference searches progressively wider bilateral ranges with the
same mesh spacing, and requires two successive agreements between adjacent
ranges plus endpoint checks. `control$max.expansions` (default 6) includes
these stability expansions; values below 2 return no set. Finer-grid retries
restart the same range sequence within the shared evaluation budget.
`print(result)` and `summary(result)` show each target/bound's searched range.
Bounded numerical search does not rule out components outside that range.
Interior audits include regular-cell quarter and three-quarter phases; finite
meshes can still miss arbitrarily narrow interior components or gaps.
`topology.audit = FALSE` skips interior mesh checks but retains the mandatory
range and endpoint checks. Full range histories are available in each task's
`diagnostics$domain.history`.

For Conditional `grid.radius > 8`, `grid.points` now sets the reference density
for radius 8. The initial count grows automatically so a wider initial range
does not coarsen that spacing. Radius 8 and smaller retain the literal initial
count. Integer rounding preserves a spacing upper bound, not exact nesting
between separate radius requests. Wider ranges cost more and can reach the
same cumulative evaluation cap. Diagnostics expose `initial.grid.points`,
`grid.density.adjusted`, `reference.domain.radius` and `reference.grid.spacing`;
settings retain the requested count. This remains finite numerical checking.

The detailed guide documents the numerical controls, `cohort_input()` for
external cohort-level coefficients and their full joint covariance, and
`build_W()`. See `vignette("detailed-guide", package = "CohortHonestDiD")`.

## Documentation

- [Getting started](vignettes/getting-started.Rmd): a short walkthrough of the
  minimum-wage application.
- [Detailed guide](vignettes/detailed-guide.Rmd): assumptions, options,
  executed comparisons, external estimates, and numerical controls.

After installing with `build_vignettes = TRUE`, open them with
`vignette("getting-started", package = "CohortHonestDiD")` and
`vignette("detailed-guide", package = "CohortHonestDiD")`. Help pages:
`?cohort_fit`, `?cohort_sensitivity`, `?cohort_input`, `?build_W`,
`?plot_cohort`, `?paper_palette`, and `?minimum_wage`.

## Citation

If you use CohortHonestDiD, please cite the paper:

> Liu, Ziyi (2026). "Cohort-Anchored HonestDiD." Working paper.

```bibtex
@unpublished{liu2026cohort,
  title  = {Cohort-Anchored {HonestDiD}},
  author = {Liu, Ziyi},
  year   = {2026},
  note   = {Working paper}
}
```

In R, `citation("CohortHonestDiD")` returns this entry and an entry for the
package.

## License and acknowledgments

CohortHonestDiD is released under the MIT license; see [LICENSE](LICENSE) and
the full license text in `LICENSE.md`. The inference routines adapt code from
the MIT-licensed HonestDiD package by Ashesh Rambachan and Jonathan Roth
(Rambachan and Roth 2023, Review of Economic Studies,
doi:10.1093/restud/rdad018); their copyright and permission notices are kept
in [inst/upstream](inst/upstream/NOTICE). The `minimum_wage` data come from the
minimum-wage application of Callaway and Sant'Anna (2021, Journal of
Econometrics 225(2), 200-230).
