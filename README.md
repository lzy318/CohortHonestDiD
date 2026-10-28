# CohortHonestDiD

<!-- badges: start -->
[![Lifecycle: stable](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/mit)
<!-- badges: end -->

**R** package for cohort-anchored HonestDiD sensitivity analysis in
difference-in-differences designs with staggered adoption. It places the
smoothness (SD) and relative-magnitude (RM) restrictions of Rambachan and Roth
(2023) on each adoption cohort's block bias, the cohort's untreated trend
relative to its comparison group at adoption, rather than on the aggregated
event study. It works with the imputation and Callaway-Sant'Anna estimators.
Confidence sets for average, dynamic, and cohort-specific effects come from the
Conditional and C-LF procedures or, under SD, the fixed-length confidence
interval (FLCI); inference on the aggregated event study is available as a
benchmark.

**Installation:**

```r
# install.packages("remotes")
remotes::install_github("lzy318/CohortHonestDiD")
library(CohortHonestDiD)
```

On Debian or Ubuntu, install the GLPK library (`libglpk-dev`) first; the Rglpk
dependency needs it.

**Source Code:** [GitHub](https://github.com/lzy318/CohortHonestDiD)

**User Manual:** [Getting started](vignettes/getting-started.Rmd) and
[Detailed guide](vignettes/detailed-guide.Rmd). To read them in R, install with
`remotes::install_github("lzy318/CohortHonestDiD", dependencies = TRUE, build_vignettes = TRUE)`
and run `browseVignettes("CohortHonestDiD")`.

**Main References:**

Liu, Ziyi (2026). [Cohort-Anchored HonestDiD](https://arxiv.org/abs/2509.01829). arXiv:2509.01829.

Rambachan, Ashesh, and Jonathan Roth (2023). [A More Credible Approach to Parallel Trends](https://doi.org/10.1093/restud/rdad018). *Review of Economic Studies* 90 (5): 2555-2591.

The inference routines adapt code from the
[HonestDiD](https://github.com/asheshrambachan/HonestDiD) package by Ashesh
Rambachan and Jonathan Roth.

**Report bugs:** Open a
[GitHub issue](https://github.com/lzy318/CohortHonestDiD/issues) or email
zyliu2023 \[at\] berkeley.edu, and include a minimal reproducible example: a
short script, a small data set that shows the problem, and the output of
`sessionInfo()`.
