# Architecture — CohortHonestDiD

> Architecture of CohortHonestDiD 0.2.0. The package was formerly distributed as HonestDiDbyCohort (versions 0.1.0
> and 0.1.1). Version 0.2.0 renames it, prepares the repository for GitHub, makes the default inference procedure
> follow the paper's recommendation (C-LF under SD and RM, with the post-only restriction under SD), adds the
> bias-reduced CR2 covariance as an option of `cohort_fit()`, adds balanced targets and balanced aggregation to
> `cohort_sensitivity()` (`balanced = TRUE`, with the cohort-side option `cohort.paths`), and adds the extractor
> `flci_weights()` and the plot type `plot(result, type = "weights")`, which show the linear weights of FLCI results;
> see [NEWS](NEWS.md). Calls that state `method` explicitly, use the default CR1 covariance and keep the default
> `balanced = FALSE` return the same numerical results as version 0.1.1; the FLCI weights add no element to any
> result. Blue nodes in the diagrams mark the code changed in 0.2.0.

## Overview

CohortHonestDiD implements cohort-anchored HonestDiD sensitivity analysis for difference-in-differences designs with
staggered adoption, as proposed in Liu (2026), "Cohort-Anchored HonestDiD" (arXiv:2509.01829). It estimates full
cohort-period coefficients and their joint cluster covariance, then applies the smoothness (SD) and
relative-magnitude (RM) restrictions of Rambachan and Roth (2023) to each cohort's block bias through an
estimator-specific bias map `W`. The package accepts either a panel (`cohort_fit()`) or named external estimates
(`cohort_input()`). Both routes produce the same `cohort_fit` representation; `cohort_sensitivity()` returns every
confidence-set component together with explicit task statuses. Following the paper, `cohort_sensitivity()` uses the
C-LF procedure by default, under the post-only restriction for SD; FLCI and Conditional remain available by name.
The covariance is the cluster-robust CR1 by default; `cohort_fit(covariance = "CR2")` gives the bias-reduced CR2
covariance of Bell and McCaffrey (2002), which the paper recommends with Conditional or C-LF inference when some
cohorts contain only a few units or clusters. With `cohort_sensitivity(balanced = TRUE)` the targets use only the
balanced window, the relative periods observed for every cohort, and aggregate inference forms the aggregated event
study only on that window, with cohort weights fixed across periods (balanced aggregation). Under cohort inference,
`cohort.paths` chooses whether the restriction uses each cohort's full path (`"full"`, the default) or only its part
inside the window (`"window"`). For FLCI results, `flci_weights()` tabulates the linear estimator's weights: the
weights `v` on the pre-treatment coefficients, in the paper's form `v = v^P + D'gamma` (the post-only weights plus a
combination of the pure pre-treatment second differences), or the weights `gamma` on those second differences, with
each weight's contribution to the FLCI center, by element or summed by cohort; `plot(result, type = "weights")` draws
them as a heatmap against the bound. The package is written in R (>= 4.1.0). `Matrix` supports sparse estimation,
`ggplot2` supplies plots, `parallel` and `withr` run reproducible worker tasks, and a bundled private numerical engine
under `inst/upstream/R/` uses the linear-programming (`Rglpk`, `lpSolveAPI`), normal-probability (`mvtnorm`,
`TruncatedNormal`) and matrix packages listed in Imports. The optional nested FLCI solver uses `CVXR` (Suggests).

### Repository layout

| Path | Contents | In the built package |
| --- | --- | --- |
| `DESCRIPTION`, `NAMESPACE` | Package metadata (MIT license); 6 exported functions, 10 registered S3 methods, imports | Yes |
| `R/` (24 files) | Package code; see Module Structure | Yes |
| `inst/upstream/R/` (5 files) | Bundled numerical engine, loaded at run time into a private environment | Yes |
| `inst/upstream/NOTICE`, `LICENSE-HonestDiD`, `LICENSE-HonestDiD-MIT` | Attribution and MIT notices for the routines adapted from HonestDiD | Yes |
| `inst/CITATION` | `citation("CohortHonestDiD")`: the paper and the package | Yes |
| `man/` (8 Rd files) | Help pages `cohort_fit`, `cohort_input`, `cohort_sensitivity`, `flci_weights`, `build_W`, `plot_cohort`, `paper_palette`, `minimum_wage`; the argument reference of the package | Yes |
| `data/minimum_wage.rda` | County panel of the minimum-wage application (15,988 rows, 28 columns) | Yes |
| [`vignettes/getting-started.Rmd`](vignettes/getting-started.Rmd) | Short walkthrough of the minimum-wage application with the defaults: estimation, the default C-LF analysis under SD, the other procedures and scopes (shown, not run; the FLCI paragraph points to `plot(flci, type = "weights")` and the detailed guide), dynamic and cohort-dynamic targets, the aggregate benchmark with a short paragraph on balanced aggregation (no chunk), RM (shown, not run), data requirements and diagnostics, a CR2 fit, and a short note on numerical results with pointers to the detailed guide and help pages for numerical controls and external estimates | Yes |
| [`vignettes/detailed-guide.Rmd`](vignettes/detailed-guide.Rmd) | Assumptions and terms; estimation with the three estimators; weights, clusters, panel histories and the CR2 covariance; targets and representations, with a full-path FLCI comparison, an executed balanced example on the minimum-wage fit (subsection "Balanced targets and balanced aggregation": cohort inference with full and window paths, and aggregate inference) and an executed FLCI weights example (subsection "FLCI weights": the post-only and full-path weights `v`, the weights `gamma` with cell labels, a table of centers and sums of contributions, and the cohort sums); procedures and SD scopes; RM; reading results (with a pointer to `flci_weights()`) and outcome statuses; numerical controls and FLCI solvers; external estimates and `W`; plots (with the weights heatmap). The procedure, RM, result-reading, control and external-input examples use a small artificial panel | Yes |
| `tests/testthat.R`, `tests/testthat/` (16 files) | testthat (edition 3) suite: fits, external input, sensitivity, default procedures, CR2 covariance, balanced targets and balanced aggregation, FLCI weights, RM, refinement, Conditional domains, boundaries, plots, method equivalence, session state, adversarial oracles | Yes |
| [`README.md`](README.md), [`NEWS.md`](NEWS.md), `LICENSE` | Repository front page, change log, CRAN-style license stub | Yes |
| `LICENSE.md` | Full MIT license text | No (build-ignored) |
| `Architecture.md` | This file | No (build-ignored) |
| `data-raw/prepare-minimum-wage.R` | Adds the convenience aliases to the study input to produce `data/minimum_wage.rda`; the input file is not in the repository | No (build-ignored) |
| `.Rbuildignore`, `.gitignore` | Build and Git exclusions | No |

---

## Module Structure

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    subgraph API["Public interface"]
        Input["input.R"]
        External["external_input.R"]
        Sens["sensitivity.R"]
        Methods["methods.R"]
        Plot["plot / palette"]
        Weights["weights.R"]
    end
    subgraph Estimate["Estimation"]
        Imp["estimate_imp.R"]
        CS["estimate_cs_nyt / cs_nt"]
        Cov["covariance.R"]
        Map["mapping.R"]
        Diag["diagnostics.R"]
    end
    subgraph Dispatch["Representation and tasks"]
        Validate["validate.R"]
        Targets["targets.R"]
        Parallel["parallel.R"]
        Utils["utils-engine.R"]
        Search["utils-search.R"]
        Conditional["utils-conditional-search.R"]
    end
    subgraph Inference["Inference adapters"]
        Support["inference_support.R"]
        SD["inference_sd.R"]
        RM["inference_rm.R"]
        FLCI["flci.R"]
    end
    subgraph Private["Bundled private engine"]
        Geometry["restriction_geometry.R"]
        Adapter["inference_adapter.R"]
        Engine["inference_engine.R"]
        Core["flci_core.R"]
        ARP["arp-nuisance.R"]
    end
    Input --> Validate
    Input --> Imp
    Input --> CS
    Input --> Cov
    Input --> Map
    Input --> Diag
    Cov --> Imp
    Map --> Validate
    External --> Validate
    External --> Utils
    Sens --> Validate
    Sens --> Targets
    Sens --> Utils
    Sens --> Parallel
    Parallel --> Support
    Support --> Utils
    Support --> SD
    Support --> RM
    SD --> FLCI
    SD --> Search
    SD --> Conditional
    SD --> Geometry
    SD --> Adapter
    SD --> Engine
    RM --> Search
    RM --> Conditional
    RM --> Geometry
    RM --> Adapter
    RM --> Engine
    FLCI --> Engine
    Conditional --> Engine
    Utils -- loads --> Geometry
    Utils -- loads --> Adapter
    Utils -- loads --> Engine
    Adapter --> Geometry
    Adapter --> Engine
    Engine --> Core
    Engine --> ARP
    Input --> Methods
    Sens --> Methods
    Input --> Plot
    Sens --> Plot
    Sens --> Weights
    FLCI -- stored weights --> Weights
    Plot --> Weights
    Weights --> Targets
    style Input fill:#1e90ff,stroke:#1565c0,color:#fff
    style Sens fill:#1e90ff,stroke:#1565c0,color:#fff
    style Methods fill:#1e90ff,stroke:#1565c0,color:#fff
    style Cov fill:#1e90ff,stroke:#1565c0,color:#fff
    style Parallel fill:#1e90ff,stroke:#1565c0,color:#fff
    style Utils fill:#1e90ff,stroke:#1565c0,color:#fff
    style RM fill:#1e90ff,stroke:#1565c0,color:#fff
    style Engine fill:#1e90ff,stroke:#1565c0,color:#fff
    style Targets fill:#1e90ff,stroke:#1565c0,color:#fff
    style Plot fill:#1e90ff,stroke:#1565c0,color:#fff
    style SD fill:#1e90ff,stroke:#1565c0,color:#fff
    style Geometry fill:#1e90ff,stroke:#1565c0,color:#fff
    style Adapter fill:#1e90ff,stroke:#1565c0,color:#fff
    style Weights fill:#1e90ff,stroke:#1565c0,color:#fff
```

Arrows describe calls or data dependencies; `covariance.R` calls the operator-score helper in `estimate_imp.R`, which
serves all three estimators. `weights.R` reads FLCI results: the weights that `flci.R` stores in each result
("stored weights") and the representation that `targets.R` rebuilds; the weights plot in `plot.R` calls it. Bundled
engine files live under `inst/upstream/R/`; all other modules live under `R/`. The inference adapters reach engine
functions through the private environment that `utils-engine.R` builds ("loads" edges). Blue nodes are the files
whose code changed in 0.2.0:

- `input.R`, `sensitivity.R`, `covariance.R` and `methods.R` carry the new defaults, the covariance option, the
  `balanced` and `cohort.paths` options, and the printed covariance and balanced-window lines;
- `targets.R` computes the balanced window, the balanced targets, the window-path cohort representation and the
  balanced aggregate representation; `plot.R` draws balanced results on the window and, with `type = "weights"`, the
  FLCI weights (`palette.R`, in the same node, is unchanged);
- `weights.R` (new) holds `flci_weights()` and its helpers (see "FLCI weights detail" below);
- `inference_sd.R`, `inference_rm.R`, `restriction_geometry.R` and `inference_adapter.R` pass the window's lower end
  to the restriction operators under `cohort.paths = "window"` (see "Balanced detail" below);
- `parallel.R`, `utils-engine.R`, `inference_rm.R` and `inference_engine.R` are the four files that name the package
  at run time (see "Package-name lookups" below).

Apart from the balanced options, the numerical code of the blue inference and engine files is unchanged. With the
default `balanced = FALSE` the balanced code is not entered: targets, representations and restriction operators are
built by the same operations as before the options existed. The FLCI weights changed no inference code: `flci.R`
and the engine store the same weights as before, and `weights.R` only reads them.

### Module Reference

| Module / File | Layer | Purpose | Key Exports | Changed in 0.2.0 |
| --- | --- | --- | --- | --- |
| `R/input.R` | API | Validate and complete the panel registry; forward-fill D; prepare working weights and observations; validate the covariance option; documents `minimum_wage` | `cohort_fit` | Yes: `covariance` argument (after `...`), `covariance_convention` built from the type, `cluster_shares` on CR2 fits, covariance type named in the two-cluster error; help text, including the `@source` of `minimum_wage` |
| `R/external_input.R` | API | Align named coefficients, covariance, `W` and normalization constraints | `cohort_input` | Help text only (`vcov` is used as supplied) |
| `R/sensitivity.R` | API | Validate options and resolve the method and scope defaults; set up the balanced window; form target/bound tasks; assemble public status and interval rows, the settings and, for balanced aggregate results, the `aggregation` element | `cohort_sensitivity` | Yes: `method = NULL` resolves to C-LF under SD and RM (before 0.2.0 the SD default was the full-path FLCI); formals `balanced = FALSE` and `cohort.paths = NULL` appended after `control`, their validation after all existing checks, the window set-up, and the settings fields and `aggregation` element of balanced results; help text, including the `control` pointer to the detailed guide's "Numerical controls" section, the section "Balanced targets" and a value sentence pointing from `results[[k]]$weights` to `flci_weights()` |
| `R/weights.R` | API | FLCI weights: for each target with weights, rebuild the representation the result used, build its second-difference operators, the post-only weights `v^P`, the estimated second differences and the post-only center in closed form, check them against the stored `gamma` and coefficient weights, and return `v = v^P + D'gamma` or `gamma` with estimates and contributions, by element or summed by cohort | `flci_weights` | Yes (new file): `flci_weights`, `.hc_flci_operators`, `.hc_flci_reconstruct`, `.hc_flci_cohort_sums`, `.hc_flci_weights_method`, `.hc_flci_weights_no_cohorts`, `.hc_flci_weights_fail`, `.hc_flci_weights_empty`; constants `.hc_weights_rank_tolerance` (1e-10) and `.hc_weights_check_tolerance` (1e-8) |
| `R/estimate_imp.R` | Estimation | Sparse outcome operator and untreated unit/time fixed-effect correction; operator application and cluster scores | None | No |
| `R/estimate_cs_nyt.R`, `R/estimate_cs_nt.R` | Estimation | Weighted paired long differences with not-yet-treated or never-treated comparisons | None | No |
| `R/covariance.R` | Estimation | Covariance option validation; cohort-period centering, cluster scores, and the CR1 or CR2 covariance; exact zero scores for single-cluster cells; CR2 shares and factors of (cluster, cohort group) pairs | None | Yes: new `.hc_covariance_type` and `.hc_cr2_adjustment`; `type` argument of `.hc_covariance` |
| `R/mapping.R` | Estimation | Full bias map using scale-stable cohort-mass ratios | `build_W` | No |
| `R/diagnostics.R` | Estimation | Observation, cohort, cluster, normalization and covariance diagnostics; stores the CR1 factor (`NA` under CR2) | None | No |
| `R/validate.R` | Representation | Input assertions, estimator aliases, coordinate registry, normalization and cell checks | None | No |
| `R/targets.R` | Representation | Balanced window; ATT/dynamic/cohort-dynamic weights; the cohort (full or window paths) and aggregate representations | None | Yes: new `.hc_balanced_window` and `.hc_check_window_paths`; `window` argument of `.hc_targets`; `window` and `paths` arguments of `.hc_representation` |
| `R/parallel.R` | Tasks | Fixed task seeds, serial/socket dispatch and worker cleanup | None | Package-name string (2 sites) |
| `R/utils-engine.R` | Representation | Numerical controls, alignment to canonical coordinates, regular grid with optional certified bounds, result/failure constructors, private engine loading from the namespace path | None | Package-name string (1 site) |
| `R/utils-search.R` | Tasks | Shared point-evaluation budget and cache; eligible refinements and attempt diagnostics | None | No |
| `R/utils-conditional-search.R` | Tasks | Conditional domain construction, bilateral expansion, interior and endpoint audits, remembered acceptance checks and range diagnostics | None | No |
| `R/inference_support.R` | Inference | Common task entry; prehistory and deterministic-target checks | None | No |
| `R/inference_sd.R` | Inference | SD restriction and plug-in set; FLCI, Conditional or C-LF; C-LF support bounds prepared once for all attempts | None | Yes: passes the representation's `restriction_window` to `v3_sd_restriction` as `pre_window` (one line) |
| `R/flci.R` | Inference | FLCI rule and interval adapter; stores each result's `weights` (`gamma` and the coefficient weights), which `flci_weights()` reads | None | No |
| `R/inference_rm.R` | Inference | Signed RM components, optional component workers and component union | None | Package-name string (3 sites); `.hc_infer_rm` counts the signed components from the benchmark pairs inside the window and passes `pre_window` to the scenario table and the plug-in restriction |
| `inst/upstream/R/restriction_geometry.R` | Private engine | Coordinate parsing, first/second differences, SD and RM restriction matrices, RM signed scenarios | None | Yes: optional `pre_window = NULL` in `v3_first_difference_operator`, `v3_second_difference_operator`, `v3_sd_restriction`, `v3_rm_scenario_table`, `v3_rm_component_restriction` and `v3_rm_plugin_restriction`; it drops only pure pre-treatment rows that reach before the window |
| `inst/upstream/R/inference_adapter.R` | Private engine | Normalization reduction, plug-in sets, ARP geometry, RM union calculations, reusable support context and budgeted sentinel tests | None | Header comment; three lines in `v3_prepare_rm_union` and `v3_clf_rm_union` pass `restriction_window` to the RM scenario table and plug-in restriction (the file keeps its CRLF line endings and line count) |
| `inst/upstream/R/inference_engine.R` | Private engine | FLCI preparation, Gaussian tests, confidence-set inversion with explicit expansion spacing, endpoint and topology audits; loads the FLCI core and HonestDiD routines | None | Header comment; package-name string (1 site) |
| `inst/upstream/R/flci_core.R` | Private engine | Folded-normal objective, direct/nested optimization and PSD limits | None | Header comment only |
| `inst/upstream/R/arp-nuisance.R` | Private engine | Gaussian moment-inequality (ARP) calculations adapted from HonestDiD | None | No |
| `R/methods.R` | Output | S3 accessors; print/summary with the covariance label and, for balanced results, the balanced-window line; Conditional searched-range tables | `coef`, `vcov`, `print` and `summary` methods | Yes: new `.hc_covariance_label`; one `Covariance:` line in the four print methods; summary objects gain `covariance`; new `.hc_balanced_label` and one balanced-window line in `print.cohort_sensitivity` and `print.summary_cohort_sensitivity` |
| `R/plot.R`, `R/palette.R` | Output | Component-preserving plots and paper colors; bounded-search plot captions; FLCI weights heatmap | `plot` methods; `paper_palette` | `plot.cohort_sensitivity` only: balanced-window caption sentence; conventional points limited to the window for dynamic and cohort-dynamic targets; formals `type = c("sets", "weights")`, `weights`, `by`, `targets = NULL`, `labels = FALSE` appended after `color`, with `type = "sets"` (default) running the earlier body unchanged; new internal `.hc_plot_weights` |
| `R/imports.R`, `NAMESPACE`, `DESCRIPTION` | Package | Import declarations, exports, S3 registrations and metadata | 6 exported functions | `DESCRIPTION` (name, version, title, description, URLs); `NAMESPACE` gains `export(flci_weights)` |
| `data/`, `man/`, `vignettes/`, `README.md`, `NEWS.md`, `inst/CITATION` | Data/docs | Application data, help pages, two vignettes, front page, change log, citation | `minimum_wage` data | Package name in vignettes; rewritten README; new NEWS entry and CITATION; `minimum_wage` help source; new defaults and the CR2 option in the `cohort_fit`, `cohort_sensitivity` and `cohort_input` help pages, README, NEWS and both vignettes; both vignettes rewritten and shortened, and the detailed guide's argument reference removed (the help pages document every argument); the `balanced` and `cohort.paths` options in the `cohort_sensitivity` help page (arguments, section "Balanced targets", example), one README paragraph, a NEWS bullet, a getting-started paragraph, and a detailed-guide subsection `{#balanced}` with an executed example and a `result$aggregation` row in its "Outcomes and statuses" table; the FLCI weights in a new help page `flci_weights`, the `plot_cohort` page (new arguments, details, value) and the `cohort_sensitivity` value section, one README sentence, a NEWS bullet, a getting-started pointer in the FLCI paragraph, and a detailed-guide subsection `{#flci-weights}` with executed minimum-wage chunks plus two cross-references (sections 6 and 9) |
| `tests/testthat.R`, `tests/testthat/` | Tests | Regression and independent public-interface checks | None | Package name in `testthat.R` and `test-session-state.R`; `method = "FLCI"` added to the calls that relied on the earlier SD default; new `test-defaults.R`, `test-covariance-cr2.R`, `test-balanced.R` and `test-flci-weights.R` |

---

## Function Call Graph

### Estimation and dispatch

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    Fit["cohort_fit"] --> Type[".hc_covariance_type"]
    Fit --> Prepare[".hc_prepare_data"]
    Prepare --> Work[".hc_working_data"]
    Work --> Est[".hc_estimate_* (3)"]
    Est --> Apply[".hc_apply_operator"]
    Est --> Cov[".hc_covariance"]
    Type --> Cov
    Fit --> Map["build_W"]
    Fit --> Diag[".hc_fit_diagnostics"]
    Cov --> Diag
    Apply --> Object["cohort_fit object"]
    Cov --> Object
    Map --> Object
    Diag --> Object
    External["cohort_input"] --> Object
    Object --> Sens["cohort_sensitivity"]
    Object --> Methods["print / summary"]
    Sens --> Methods
    Methods --> Label[".hc_covariance_label"]
    Methods --> BLabel[".hc_balanced_label"]
    Sens --> Controls[".hc_engine_controls"]
    Sens --> Window[".hc_balanced_window"]
    Window --> Targets[".hc_targets"]
    Sens --> Targets
    Targets --> Rep[".hc_representation"]
    Sens --> Seeds[".hc_task_seeds"]
    Rep --> Run[".hc_run_tasks"]
    Seeds --> Run
    Run --> Task[".hc_task"]
    Task --> One[".hc_infer_one"]
    style Fit fill:#1e90ff,stroke:#1565c0,color:#fff
    style Type fill:#1e90ff,stroke:#1565c0,color:#fff
    style Prepare fill:#1e90ff,stroke:#1565c0,color:#fff
    style Cov fill:#1e90ff,stroke:#1565c0,color:#fff
    style Sens fill:#1e90ff,stroke:#1565c0,color:#fff
    style Methods fill:#1e90ff,stroke:#1565c0,color:#fff
    style Label fill:#1e90ff,stroke:#1565c0,color:#fff
    style BLabel fill:#1e90ff,stroke:#1565c0,color:#fff
    style Window fill:#1e90ff,stroke:#1565c0,color:#fff
    style Targets fill:#1e90ff,stroke:#1565c0,color:#fff
    style Rep fill:#1e90ff,stroke:#1565c0,color:#fff
    style Run fill:#1e90ff,stroke:#1565c0,color:#fff
```

Arrows follow the order of calls and the objects they produce. `cohort_fit` passes the validated covariance type to
`.hc_covariance` and the covariance list it returns to `.hc_fit_diagnostics`. `.hc_balanced_window` runs only with
`balanced = TRUE`; its window is passed to `.hc_targets` and `.hc_representation` (see "Balanced detail").

| Function | Defined In | Called By | Calls | Changed in 0.2.0 | Purpose |
| --- | --- | --- | --- | --- | --- |
| `cohort_fit` | `input.R` | User | Covariance-type check, preparation, estimation, covariance, mapping and diagnostics | Yes: `covariance = "CR1"` after `...`; `covariance_convention`; `cluster_shares` on CR2 fits | Construct the fitted representation once |
| `.hc_covariance_type` | `covariance.R` | `cohort_fit` | `.hc_assert` | Yes (new) | Accept `"CR1"` or `"CR2"` in any case and return it in upper case; anything else errors with "covariance must be 'CR1' or 'CR2'." |
| `.hc_prepare_data`, `.hc_working_data` | `input.R` | `cohort_fit` | Validation and registry helpers | `.hc_prepare_data`: error message only (names CR1 or CR2) | Forward-fill D, retain usable Y, record any common weight divisor, require two clusters, compress balanced panels to cohort-period cells |
| `.hc_estimate_imp`, `.hc_estimate_cs_nyt`, `.hc_estimate_cs_nt` | Estimator files | `cohort_fit` | Sparse operator helpers | No | Construct the outcome contrast operator |
| `.hc_covariance` | `covariance.R` | `cohort_fit` | `.hc_cr2_adjustment` (CR2 only), `.hc_operator_scores` | Yes: `type` argument | Full joint covariance from centered cluster scores, CR1 or CR2 (see "Covariance detail") |
| `build_W` | `mapping.R` | User / `cohort_fit` | `create_delta_to_Delta_mapping_matrix` | No | Map primitive deviations to coefficient bias |
| `.hc_fit_diagnostics` | `diagnostics.R` | `cohort_fit` | Cohort, cell and eigenvalue summaries | No | Support and covariance diagnostics stored on the fit, computed from whichever `Sigma` was estimated |
| `cohort_input` | `external_input.R` | User | `.hc_check_mapping_support` and validation | No | Construct the same object from external estimates; the supplied `vcov` is used as given |
| `cohort_sensitivity` | `sensitivity.R` | User | Controls, balanced window, targets, representations, seeds and task runner | Yes: default dispatch; `balanced` and `cohort.paths` | Resolve `method` (NULL means C-LF) and, for SD, `full.path` (NULL means TRUE only for FLCI); with `balanced = TRUE`, validate `cohort.paths` (NULL means `"full"` under cohort inference), compute the window and record it; return estimates and all confidence-set components |
| `.hc_balanced_window` | `targets.R` | `cohort_sensitivity` (balanced only) | None | Yes (new) | The balanced window `c(lower = L, upper = S)` from `coordinate_info` (see "Balanced detail") |
| `.hc_check_window_paths` | `targets.R` | `cohort_sensitivity` (`cohort.paths = "window"` only) | None | Yes (new) | Reject a `W` or normalization that links a kept coefficient to one after the window |
| `.hc_targets`, `.hc_representation` | `targets.R` | `cohort_sensitivity` | Mass normalization and matrix transforms | Yes: `window` argument (both) and `paths` argument (`.hc_representation`); with `window = NULL` the same operations as before | Preserve the requested estimand across frameworks; with a window, restrict targets to periods 1 to `S`, form the aggregate path on `L` to `S`, or drop the post-treatment coefficients after `S` under window paths |
| `.hc_task_seeds` | `parallel.R` | `cohort_sensitivity` | None | No | Deterministic per-task seeds, independent of worker count |
| `.hc_run_tasks`, `.hc_task` | `parallel.R` | `cohort_sensitivity` | `.hc_infer_one`, `withr`, `parallel` | `.hc_run_tasks`: package-name string (2 sites) | Reproducible execution with worker cleanup |
| `print.cohort_fit`, `summary.cohort_fit`, `print.summary_cohort_fit`, `print.cohort_sensitivity`, `summary.cohort_sensitivity`, `print.summary_cohort_sensitivity` | `methods.R` | User | `.hc_covariance_label`, `.hc_balanced_label`; searched-range helpers | Yes: covariance line and summary element; balanced-window line | Printed headers name the procedure, scope and covariance, and for balanced results the window |
| `.hc_covariance_label` | `methods.R` | Print and summary methods | None | Yes (new) | `CR1 cluster-robust, G clusters`, `CR2 bias-reduced cluster-robust, G clusters`, `supplied to cohort_input()` for external input, or `not recorded` |
| `.hc_balanced_label` | `methods.R` | `print.cohort_sensitivity`, `print.summary_cohort_sensitivity` (balanced only) | None | Yes (new) | `Balanced window: relative periods L to S; target period 1` (or `target periods 1 to S`), then `; cohort paths: full`, `; cohort paths: window` or `; aggregate path on the window` |

### Covariance detail

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    Fit["cohort_fit"] --> Type[".hc_covariance_type"]
    Fit --> Cov[".hc_covariance"]
    Cov --> Resid["cell-centered residuals"]
    Resid --> Zero["zero single-cluster cells"]
    Cov --> Adj[".hc_cr2_adjustment"]
    Adj --> Factor["pair factor per unit"]
    Adj --> Shares["fit cluster_shares"]
    Zero -- CR2 --> Scale["scale residual rows"]
    Factor --> Scale
    Zero -- CR1 --> Scores[".hc_operator_scores"]
    Scale --> Scores
    Scores --> CR1["CR1: crossprod x factor"]
    Scores --> CR2["CR2: crossprod"]
    CR1 --> Sigma["symmetrized Sigma"]
    CR2 --> Sigma
    Sigma --> Diag[".hc_fit_diagnostics"]
    style Type fill:#1e90ff,stroke:#1565c0,color:#fff
    style Cov fill:#1e90ff,stroke:#1565c0,color:#fff
    style Adj fill:#1e90ff,stroke:#1565c0,color:#fff
    style Scale fill:#1e90ff,stroke:#1565c0,color:#fff
    style CR2 fill:#1e90ff,stroke:#1565c0,color:#fff
```

The CR1 branch is the covariance of version 0.1.1, with the same operations in the same order. The blue steps exist
only for `covariance = "CR2"`; `.hc_cr2_adjustment` runs before the residual matrix is built, and its shares become
the fit's `cluster_shares`.

| Step | Defined In | Called By | Calls | Changed in 0.2.0 | Purpose |
| --- | --- | --- | --- | --- | --- |
| Cell-centered residuals | `covariance.R` | Step of `.hc_covariance` | `rowsum` | No | Residual = outcome minus the weighted mean of its cohort-period cell, over usable rows |
| Residual-by-cluster matrix | `covariance.R` | Step of `.hc_covariance` | `Matrix::sparseMatrix` | No | Balanced panels: one row per cohort-period cell of the compressed working data, holding each cluster's weighted residual sum divided by the cell mass; unbalanced panels: one row per usable observed row, holding its residual |
| Zero single-cluster cells | `covariance.R` | Step of `.hc_covariance` | None | No | Balanced (compressed) path only: a cell supported by a single cluster gets exactly zero entries, without a variance cutoff. The observed-row path for unbalanced panels has no such step |
| `.hc_cr2_adjustment` | `covariance.R` | `.hc_covariance` (CR2) | Base R (`match`, `tabulate`, `rowsum`) | Yes (new) | Pairs (cluster, cohort group) from `prepared$unit_registry`, the never-treated units forming group 0; pair share p = pair weight / group weight (unit counts when unweighted); factor `1/sqrt(1 - p)`, computed as `sqrt(total/rest)`; factor 0 for a pair that holds every unit of its group; returns the `cluster_shares` data frame and one factor per unit |
| Scale residual rows | `covariance.R` | Step of `.hc_covariance` (CR2) | None | Yes (new) | Multiply each row by the factor of its unit's pair. A compressed row is one cell of one group, so every entry of the matrix is scaled by exactly one pair factor |
| `.hc_operator_scores` | `estimate_imp.R` | `.hc_covariance` | Sparse products with the operator and its fixed-effect correction | No | Apply the estimator's linear operator to every cluster column; rows of the result are the cluster scores (`fit$cluster_scores`) |
| CR1 / CR2 sum | `covariance.R` | Step of `.hc_covariance` | `crossprod` | CR2 branch new | CR1: cross-product times `G/(G - 1) * (n - 1)/(n - p)` (G clusters, n usable rows, p nonempty cells); CR2: cross-product with no factor and `cr1 = NA` |
| Symmetrized `Sigma` | `covariance.R` | Step of `.hc_covariance` | None | No | `(Sigma + t(Sigma))/2` with coordinate names; returned with the scores, the factor, the counts and, under CR2, the shares |

### Balanced detail

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    Sens["cohort_sensitivity"] --> Valid["validate the 2 options"]
    Valid --> Bal{"balanced = TRUE?"}
    Bal -- yes --> Win[".hc_balanced_window"]
    Win --> WP{"window paths?"}
    WP -- yes --> Check[".hc_check_window_paths"]
    Bal -- no --> Targets[".hc_targets"]
    WP -- no --> Targets
    Check --> Targets
    Targets --> Rep[".hc_representation"]
    Rep -- cohort, full --> Full["fit representation"]
    Rep -- cohort, window --> Trunc["cells to S + lower end"]
    Rep -- aggregate --> Agg["A on periods L..S"]
    Full --> Infer[".hc_infer_sd / _rm"]
    Trunc --> Infer
    Agg --> Infer
    Infer --> Ops["operators, pre_window"]
    Ops --> Out["settings, aggregation"]
    Out --> BLabel[".hc_balanced_label"]
    Out --> Plot["plot.cohort_sensitivity"]
    style Valid fill:#1e90ff,stroke:#1565c0,color:#fff
    style Win fill:#1e90ff,stroke:#1565c0,color:#fff
    style Check fill:#1e90ff,stroke:#1565c0,color:#fff
    style Targets fill:#1e90ff,stroke:#1565c0,color:#fff
    style Rep fill:#1e90ff,stroke:#1565c0,color:#fff
    style Trunc fill:#1e90ff,stroke:#1565c0,color:#fff
    style Agg fill:#1e90ff,stroke:#1565c0,color:#fff
    style Ops fill:#1e90ff,stroke:#1565c0,color:#fff
    style Out fill:#1e90ff,stroke:#1565c0,color:#fff
    style BLabel fill:#1e90ff,stroke:#1565c0,color:#fff
    style Plot fill:#1e90ff,stroke:#1565c0,color:#fff
```

The task runner and `.hc_infer_one` between the representation and the inference adapters are omitted. With
`balanced = FALSE`, `.hc_targets` and `.hc_representation` receive `window = NULL` and `paths = "full"` and return
what they returned before the options existed; the aggregate branch then uses every observed relative period. Under
`cohort.paths = "full"` the cohort representation is the same as with `balanced = FALSE`; only the target differs.

| Step | Defined In | Called By | Calls | Changed in 0.2.0 | Purpose |
| --- | --- | --- | --- | --- | --- |
| Option validation | `sensitivity.R` | `cohort_sensitivity`, after all existing checks and `.hc_engine_controls` | None | Yes (new) | `balanced` must be one non-missing logical; `cohort.paths` may be supplied only with `balanced = TRUE` and `framework = "cohort"`, and must then be exactly `"full"` or `"window"` (no partial matching); `NULL` means `"full"` under balanced cohort inference |
| `.hc_balanced_window` | `targets.R` | `cohort_sensitivity` | None | Yes (new) | From `coordinate_info` (so also for `cohort_input()` objects and unbalanced panels): `S` = the largest `s >= 1` such that periods `1..s` are observed for every cohort, `L` = the smallest `s <= 0` such that periods `s..0` are; returns the named integer `c(lower = L, upper = S)`. Errors if some cohort lacks relative period 1 (checked first) or relative period 0, which is possible only for `cohort_input()` objects. For a `cohort_fit` with adoption indices `t_1 < ... < t_G` on `T` periods, `L = 2 - t_1` and `S = T - t_G + 1` |
| `.hc_check_window_paths` | `targets.R` | `cohort_sensitivity` (`cohort.paths = "window"`) | None | Yes (new) | Error if a kept row of `W`, or a normalization row, loads on a post-treatment coefficient after `S` (scale-relative tolerance 1e-12). The bias maps that `cohort_fit()` builds never do, so only `cohort_input()` objects can fail |
| `.hc_targets(fit, target, window)` | `targets.R` | `cohort_sensitivity` | Its local `make()` (relative masses, then normalized) | Yes: `window` | Keeps the post cells with `s <= S` before the existing weight arithmetic: the ATT weights each cell of periods 1 to `S` by `m_g / (S * sum m)`, dynamic targets exist for `s = 1..S` with weights `m_g / sum m`, and cohort-dynamic targets for the cells with `s <= S`. Weight vectors still cover every coefficient of `fit`. ATT label `Average treatment effect, relative period 1` (`S = 1`) or `..., relative periods 1 to S` |
| `.hc_representation(fit, ell, framework, window, paths)` | `targets.R` | `cohort_sensitivity`, once per target | `qr.solve`, null-space embedding | Yes: `window`, `paths` | Cohort, `"full"`: the fit's own representation. Cohort, `"window"`: keeps every pre-treatment coefficient and the post-treatment coefficients up to `S` (`betahat`, `Sigma`, `W`, `normalization` and the target subset to them) and adds the element `restriction_window = L`. Aggregate: the existing construction on the period set `L..S` |
| `.hc_infer_sd`, `.hc_infer_rm` | `inference_sd.R`, `inference_rm.R` | `.hc_infer_one` | Restriction operators of the engine | Yes | Read `representation[["restriction_window"]]` (absent, hence `NULL`, except under window paths) and pass it as `pre_window`; `.hc_infer_rm` also counts the signed components from the pre-treatment cells at or after `L` |
| Restriction operators | `restriction_geometry.R` | SD/RM adapters; `v3_prepare_rm_union`, `v3_clf_rm_union` | None | Yes: `pre_window = NULL` | `v3_first_difference_operator` drops pre-treatment first differences whose earlier period is before `L` (RM benchmark candidates); `v3_second_difference_operator` drops pure pre-treatment second differences whose first period is before `L` (full-path SD). Rows that involve a post-treatment coefficient are never dropped. `v3_sd_restriction`, `v3_rm_scenario_table` and `v3_rm_plugin_restriction` pass the argument on; `v3_rm_component_restriction` uses it only when it builds the operator itself (the RM adapter passes it the windowed operator of the scenario table) |
| `v3_prepare_rm_union`, `v3_clf_rm_union` | `inference_adapter.R` | RM C-LF adapter | Scenario table, plug-in restriction | Yes (three lines) | Pass `representation[["restriction_window"]]`; `v3_align_representation` keeps the element |
| Settings and `aggregation` | `sensitivity.R` | `cohort_sensitivity` | None | Yes (new) | Only when `balanced = TRUE`: `settings` gains `balanced`, `cohort.paths` and `window` at the end, and aggregate results gain the element `aggregation` (see Notes). Estimates use the full-length target with the fit's `betahat` and `Sigma` in every case |
| `.hc_balanced_label` | `methods.R` | Print methods | None | Yes (new) | One line after `Covariance:` in `print()`, and the second line of `print(summary())` followed by a blank line |
| `plot.cohort_sensitivity` | `plot.R` | User | `.hc_event_data` (unchanged) | Yes | Caption sentence `Balanced window: relative periods L to S; cohort paths: P.` (cohort) or `Balanced window: relative periods L to S.` (aggregate); dynamic and cohort-dynamic plots draw the conventional points only for `s = L..S`; on the window the aggregated points of a dynamic plot are the balanced aggregated coefficients (`aggregation$betahat` for aggregate results) |

### Inference detail

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    One[".hc_infer_one"] --> Align["engine representation"]
    Align --> SD[".hc_infer_sd"]
    Align --> RM[".hc_infer_rm"]
    SD --> Reduce["restriction / reduction"]
    RM --> Reduce
    Reduce --> Plugin["plug-in identified set"]
    SD --> FLCI[".hc_infer_flci"]
    FLCI --> Prepare["sd_flci_prepare"]
    Prepare --> Solve["direct / nested solver"]
    Solve --> Interval["sd_flci_interval"]
    SD --> ARP["sd_bind_arp_problem"]
    RM --> Union["signed component union"]
    SD --> Budget[".hc_search_budget"]
    RM --> Budget
    Budget -- C-LF --> Retry[".hc_adaptive_search"]
    Budget -- Conditional --> Conditional[".hc_conditional_search"]
    Retry --> Invert
    Retry --> Union
    ARP -- Conditional --> Conditional
    Union -- Conditional --> Conditional
    Conditional --> Fixed["sd_invert_test: fixed domain"]
    Fixed --> Audit
    ARP -- C-LF --> LFSupport["C-LF outer support"]
    LFSupport --> Invert["sd_arp_invert"]
    Invert --> Audit["endpoint / grid audits"]
    Union --> Audit
    Interval --> Result[".hc_engine_result"]
    Audit --> Result
    Plugin --> Result
    One --> Failure[".hc_engine_failure"]
    Result --> Rows["public interval rows"]
    Failure --> Rows
    style SD fill:#1e90ff,stroke:#1565c0,color:#fff
    style RM fill:#1e90ff,stroke:#1565c0,color:#fff
    style Reduce fill:#1e90ff,stroke:#1565c0,color:#fff
```

The blue steps changed in 0.2.0 only to pass the window's lower end (`restriction_window`) to the restriction
operators under `cohort.paths = "window"`; for every other representation the element is absent and the operators
build the same rows as before. The new default changes which branch a call without `method` takes (C-LF instead of
the full-path FLCI under SD), not the branches themselves. Inference uses the covariance only through the fit's
`Sigma`, so CR1 and CR2 fits use the same code. A balanced aggregate representation enters this diagram like any
other aggregate representation.

| Function | Defined In | Called By | Calls | Changed in 0.2.0 | Purpose |
| --- | --- | --- | --- | --- | --- |
| `.hc_infer_one` | `inference_support.R` | Task runner | SD/RM adapters or explicit failure constructor | No | Check prehistory and supported stochastic quantities |
| `.hc_engine_representation`, `.hc_backend` | `utils-engine.R` | Inference adapters | Private alignment; engine loading from the namespace path | `.hc_backend`: package-name string (1 site) | Keep helpers private; support both installed and source-loaded resources |
| `.hc_infer_sd` | `inference_sd.R` | Common task entry | Restriction/reduction, plug-in, FLCI or shared-budget Gaussian inversion; C-LF support LP | Yes: `pre_window` passed to `v3_sd_restriction` | Apply either SD scope; prepare, calibrate and bound support once for all search attempts |
| `.hc_infer_rm` | `inference_rm.R` | Common task entry | Signed-component geometry, prepared support context and shared-budget union inversion | Yes: windowed component count; `pre_window` passed to the scenario table and plug-in restriction | Preserve disconnected RM sets, calibration and worker state across retries |
| `.hc_search_budget` | `utils-search.R` | SD/RM adapters | Cached point-test callback | No | Count distinct test evaluations across searches, probes, audits and retries |
| `.hc_adaptive_search` | `utils-search.R` | SD/RM C-LF adapters | Per-attempt inversion callback | No | Double regular grid intervals for eligible unresolved statuses within the limits |
| `.hc_conditional_domain`, `.hc_conditional_search` | `utils-conditional-search.R` | SD/RM Conditional adapters | Cached point tests, fixed-domain inversion and endpoint audits | No | Search concentric domains with fixed regular spacing and require two successive adjacent-domain agreements; finer retries reuse the same cache and visited-range history |
| `.hc_engine_grid` | `utils-engine.R` | SD/RM adapters | Regular grid/profile anchors, with optional certified bounds | No | Carry regular spacing separately; pad supplied C-LF support by the endpoint probe margin |
| `.hc_infer_flci` | `flci.R` | SD adapter | `sd_flci_prepare`, `sd_flci_interval` | No | Expose the optimized rule and interval diagnostics (including `tie.break`); store `weights$gamma` and `weights$coefficients`, the weights mapped back from the reduced space to the representation's coordinates (read by `flci_weights()`, see "FLCI weights detail") |
| `sd_flci_prepare` | Private `inference_engine.R` | FLCI adapter | Private direct/nested FLCI solvers | No | Handle PSD and exact zero-noise limits without jitter |
| `sd_bind_arp_problem`, `sd_arp_invert` | Private `inference_engine.R` | SD/RM adapters | Gaussian tests, calibration, explicit expansion spacing and shared point tests | No | Compute Conditional/C-LF sets without resimulating calibration on retry |
| `v3_clf_rm_union`, `v3_memoize_rm_union_test`, `v3_rm_verify_support_sentinels` | Private `inference_adapter.R` | RM adapter | Prepared component/support context, cached union tests and separately keyed sentinel tests | `v3_clf_rm_union` (and the `v3_prepare_rm_union` it calls): `restriction_window` passed to the scenario table and plug-in restriction | Reuse RM calibration and account for actual point work across refinements |
| `v3_rm_lf_component_support` | Private `inference_adapter.R` | RM adapter and SD C-LF | Primal/dual-checked support endpoint LPs | No | Bound the necessary LF first-stage support with the realized cutoff |
| `.hc_engine_result`, `.hc_engine_failure` | `utils-engine.R` | Inference adapters | Component/status assembly | No | Keep empty, unbounded, disconnected and unsupported results distinct |

### FLCI weights detail

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    Plot["plot.cohort_sensitivity"] --> Type{"type?"}
    Type -- sets --> Sets["set plot body"]
    Type -- weights --> PW[".hc_plot_weights"]
    PW --> Method[".hc_flci_weights_method"]
    PW --> NoCoh[".hc_flci_weights_no_cohorts"]
    PW --> FW["flci_weights"]
    FW --> Method
    FW --> NoCoh
    FW --> Rep[".hc_representation"]
    FW --> Ops[".hc_flci_operators"]
    FW --> Rec[".hc_flci_reconstruct"]
    Infer[".hc_infer_flci"] --> Stored["stored weights"]
    Stored --> Rec
    Ops --> Rec
    Rec --> Fail[".hc_flci_weights_fail"]
    FW --> Fail
    FW --> Sums[".hc_flci_cohort_sums"]
    FW --> Empty[".hc_flci_weights_empty"]
    PW --> Labels[".hc_cohort_labels"]
    PW --> Palette["paper_palette"]
    PW --> Heat["tile and text layers"]
    style Plot fill:#1e90ff,stroke:#1565c0,color:#fff
    style Type fill:#1e90ff,stroke:#1565c0,color:#fff
    style Sets fill:#1e90ff,stroke:#1565c0,color:#fff
    style PW fill:#1e90ff,stroke:#1565c0,color:#fff
    style Rep fill:#1e90ff,stroke:#1565c0,color:#fff
    style Method fill:#1e90ff,stroke:#1565c0,color:#fff
    style NoCoh fill:#1e90ff,stroke:#1565c0,color:#fff
    style FW fill:#1e90ff,stroke:#1565c0,color:#fff
    style Ops fill:#1e90ff,stroke:#1565c0,color:#fff
    style Rec fill:#1e90ff,stroke:#1565c0,color:#fff
    style Fail fill:#1e90ff,stroke:#1565c0,color:#fff
    style Sums fill:#1e90ff,stroke:#1565c0,color:#fff
    style Empty fill:#1e90ff,stroke:#1565c0,color:#fff
    style Heat fill:#1e90ff,stroke:#1565c0,color:#fff
```

`flci_weights()` is the public entry point; the weights plot calls it after its own argument checks. The extractor
only reads the result: it calls the unchanged `.hc_representation` with the result's own settings, so it sees the
representation the FLCI used, and it accepts the weights that `.hc_infer_flci` stored ("stored weights") only after
reconstructing them. The set plot (`type = "sets"`, the default) runs the earlier body. `.hc_representation` and
the set-plot body are blue for their balanced changes (see "Balanced detail"); the FLCI weights did not change them.
`.hc_infer_flci`, `.hc_cohort_labels` and `paper_palette` are unchanged.

**What is computed.** For one target, the representation gives the coefficients `b` (`betahat`), the bias map `W`,
the target weights `l` and the normalization `N`: under cohort inference those of the fit (or, under
`cohort.paths = "window"`, of the window-path cells, with the window's lower end `L`); under aggregate inference
those of the aggregated path (one internal path `g1`, `W = I`, target `l_A`). `P` and `Q` are its pre- and
post-treatment cells, each ordered by cohort and then by relative period.

- Second differences along each path, whose relative periods are consecutive: from the third cell on, the row
  `d_{g,e}` has `+1` at `(g, e)`, `-2` at `(g, e - 1)` and `+1` at `(g, e - 2)`. Rows with `e <= 0` are the pure
  pre-treatment second differences (under `cohort.paths = "window"` only those with `e - 2 >= L`, the rows of the
  windowed restriction); their columns on `P` form `D_PP` (`J` rows, ordered by cohort and then `e`, which is the
  engine's order of `gamma`; element names `d2_g<g>_s<e>`). Rows with `e >= 1` are post-connected, one per post
  cell, ordered like `Q`; their columns form `D_QP` and the unit lower-triangular `D_QQ`.
- `xi = solve(t(D_QQ), t(W_QQ) %*% l_Q)` and `vP = t(D_QP) %*% xi - t(W_QP) %*% l_Q`: the post-only weights on the
  pre-treatment coefficients (paper `v^P`). They subtract each cohort's linear extrapolation from its last two
  pre-treatment periods and, under CS-NYT, undo the base-period correction of `W`; they do not depend on `M`.
- `zhat = D_PP %*% b_P` (the estimated pure pre-treatment second differences) and the post-only center
  `c0 = sum(vP * b_P) + sum(l_Q * b_Q)`.
- For each task with weights: `gamma` is the stored `weights$gamma` (zeros for the post-only FLCI) and
  `v = vP + t(D_PP) %*% gamma` (paper `v = v^P + D'gamma`).
- Checks before anything is reported: `length(gamma) == J`; if the engine named `gamma` (it does so on its analytic
  `M = 0` branch), the names are `d2_g<g>_s<e>_pos` in the same order; and `Pi %*% v` equals the stored
  `weights$coefficients` on `P` within `1e-8 * max(1, max|w_P|)`, where `Pi` is the orthogonal projector onto the null
  space of the normalization's pre-treatment columns (QR rank tolerance 1e-10). Any failure stops with the
  reconstruction error.
- Contributions are `weight * estimate`. They add up to `center - c0` over the `gamma` rows and to
  `center - estimate` (the conventional estimate) over the `v` rows of a task.

| Function | Defined In | Called By | Calls | Changed in 0.2.0 | Purpose |
| --- | --- | --- | --- | --- | --- |
| `flci_weights(x, weights, by)` | `weights.R` | User; `.hc_plot_weights` | The checks below, `.hc_representation`, `.hc_flci_operators`, `.hc_flci_reconstruct`, `.hc_flci_cohort_sums`, `.hc_flci_weights_empty`, `.hc_flci_weights_fail` | Yes (new, exported) | Check the object, method and framework; loop over targets `j` and bounds `m` (task `k = (j - 1) * length(settings$M) + m`); list tasks without weights in `missing`; rebuild the operators once per target that has weights; one `centers` row and the element rows per task with weights; return the data frame with its attributes (see Notes). It never messages or warns |
| `.hc_flci_weights_method` | `weights.R` | `flci_weights`, `.hc_plot_weights` | None | Yes (new) | Error unless `settings$method` is `"flci"`; the message names Conditional or C-LF |
| `.hc_flci_weights_no_cohorts` | `weights.R` | `flci_weights`, `.hc_plot_weights` | None | Yes (new) | Error for `by = "cohort"` under aggregate inference |
| `.hc_representation` | `targets.R` | `cohort_sensitivity`; `flci_weights` | See "Balanced detail" | Yes for the balanced options (see "Balanced detail"); unchanged by the FLCI weights, which only add a caller | Called with `settings$window` and, for balanced cohort inference, `settings$cohort.paths` (otherwise `"full"`), as `cohort_sensitivity()` called it |
| `.hc_flci_operators` | `weights.R` | `flci_weights`, once per target with weights | `solve`, `qr`, `crossprod` | Yes (new) | Second-difference rows (pure, post-connected, or outside the window), `D_PP`, `xi`, `vP`, `zhat`, `c0`, `Pi`, and the labels of the coefficient and second-difference rows. Errors, which `flci_weights` turns into the reconstruction error, if a post cell has no post-connected row, the triangular solve fails or a value is not finite |
| `.hc_flci_reconstruct` | `weights.R` | `flci_weights`, once per task with weights | `.hc_flci_weights_fail` | Yes (new) | The checks above; returns `v` named by the pre-treatment coordinates |
| `.hc_flci_cohort_sums` | `weights.R` | `flci_weights` (`by = "cohort"`), once per task | None | Yes (new) | Per cohort with at least one row: `elements`, `weight` (sum), `l1` (sum of absolute weights), `share` (`l1` over the task's total, `NA` when that is 0), `contribution` (sum), `tie.break` |
| `.hc_flci_weights_empty`, `.hc_flci_weights_fail` | `weights.R` | `flci_weights`, `.hc_flci_reconstruct` | None | Yes (new) | The 0-row data frame with the documented columns and types; the reconstruction error |
| `.hc_infer_flci` | `flci.R` | SD adapter | Engine FLCI routines | No | Stores `weights$gamma` (the engine's `gamma`; zeros for the post-only FLCI) and `weights$coefficients` (the weights mapped back from the reduced space: the projection of `v` orthogonal to the normalization, with the target weights on the post cells), and `diagnostics$tie.break` |
| `plot.cohort_sensitivity` | `plot.R` | User | `.hc_plot_weights` or the set-plot body | Yes | Match `type` first. `"weights"`: error if `identified.set` or `color` is supplied, then `.hc_plot_weights` (a missing `ylab` means the default label). `"sets"`: error if `weights`, `by`, `targets` or `labels` is supplied, then the earlier body statement for statement |
| `.hc_plot_weights` | `plot.R` | `plot.cohort_sensitivity` | `.hc_flci_weights_method`, `.hc_flci_weights_no_cohorts`, `flci_weights`, `.hc_cohort_labels`, `paper_palette` | Yes (new) | Check `labels`, `targets`, the method, `by` and `cohorts` under aggregate inference; for `weights = "gamma"` on a post-only result print a message and return `NULL` invisibly; check `M`; select tasks by `targets` and by `M` (the set plot's matching rule, relative 1e-10); error if nothing is selected or every selected task failed; keep the selected rows of `flci_weights(x, weights, by)` (message and invisible `NULL` if no `gamma` rows exist), filter `cohorts`, and build the heatmap |
| Tile and text layers | `plot.R` | `.hc_plot_weights` | `ggplot2` | Yes (new) | `geom_tile` (grey borders), optional `geom_text` labels, the diverging fill, discrete axes without expansion, facets, theme and caption (see Notes) |

### Package-name lookups

Socket workers start as empty R sessions, and the private engine is located on disk through the package namespace.
These are the only places where the code names its own package; the rename in 0.2.0 changed exactly these seven
string literals.

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    Sens["cohort_sensitivity"] --> Run[".hc_run_tasks"]
    Run --> Load1["worker: loadNamespace"]
    Run --> Get1["worker: .hc_task by name"]
    RM[".hc_infer_rm"] --> Cluster[".hc_rm_cluster"]
    Cluster --> Load2["worker: loadNamespace"]
    RM --> Job[".hc_rm_worker_job"]
    Job --> Get2["getFromNamespace: job"]
    RM --> WTest[".hc_rm_worker_test"]
    WTest --> Get3["getFromNamespace: batch"]
    Adapters["SD/RM/FLCI adapters"] --> Backend[".hc_backend"]
    Backend --> Path1["namespace path"]
    Path1 --> Source["sys.source 5 engine files"]
    Source --> Deps["engine dependency loader"]
    Deps --> Root[".sd_ie_find_project_root"]
    Root --> Path2["namespace path"]
    style Run fill:#1e90ff,stroke:#1565c0,color:#fff
    style Cluster fill:#1e90ff,stroke:#1565c0,color:#fff
    style Job fill:#1e90ff,stroke:#1565c0,color:#fff
    style WTest fill:#1e90ff,stroke:#1565c0,color:#fff
    style Backend fill:#1e90ff,stroke:#1565c0,color:#fff
    style Root fill:#1e90ff,stroke:#1565c0,color:#fff
```

| Function | Defined In | Called By | Calls | Changed in 0.2.0 | Purpose |
| --- | --- | --- | --- | --- | --- |
| `.hc_run_tasks` | `R/parallel.R` | `cohort_sensitivity` | `makePSOCKcluster`; `clusterCall` with `loadNamespace("CohortHonestDiD")`; `parLapplyLB` with `asNamespace("CohortHonestDiD")` | Yes (2 literals) | Run target/bound tasks on socket workers that load the package from the parent's library paths |
| `.hc_rm_cluster` | `R/inference_rm.R` | `.hc_infer_rm` when `component.cores > 1` | `loadNamespace("CohortHonestDiD")` on each worker | Yes (1) | Start component workers for a single RM task |
| `.hc_rm_worker_job` | `R/inference_rm.R` | `parLapply` in `.hc_infer_rm` | `getFromNamespace(".hc_rm_component_job", "CohortHonestDiD")` | Yes (1) | Prepare one signed RM component on a worker |
| `.hc_rm_worker_test` | `R/inference_rm.R` | RM point test | `getFromNamespace(".hc_rm_test_batch", "CohortHonestDiD")` | Yes (1) | Test one batch of components at a point on a worker |
| `.hc_backend` | `R/utils-engine.R` | All inference adapters | `getNamespaceInfo(asNamespace("CohortHonestDiD"), "path")`; `sys.source` | Yes (1) | Load the five engine files once into a cached private environment, from `upstream/R` (installed) or `inst/upstream/R` (source-loaded) |
| `.sd_ie_find_project_root` | `inst/upstream/R/inference_engine.R` | `sd_inference_load_dependencies` | `getNamespaceInfo(asNamespace("CohortHonestDiD"), "path")` | Yes (1) | Locate the same resource directory for the engine's own dependency environment |

---

## Data Flow

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
    Start{"Panel or estimates?"}
    Start -- Panel --> Panel["Registry and estimation"]
    Start -- Estimates --> Ext["Named input validation"]
    Panel --> CovType{"covariance?"}
    CovType -- CR1 --> CR1["Scores x CR1 factor"]
    CovType -- CR2 --> CR2["Scores x pair factors"]
    CR1 --> Fit["beta, Sigma, W, N"]
    CR2 --> Fit
    Ext --> Fit
    Fit --> Given{"method given?"}
    Given -- No --> Default["C-LF; SD post-only"]
    Given -- Yes --> Explicit["Requested method/scope"]
    Default --> Target["Target and framework"]
    Explicit --> Target
    Target --> Balanced{"balanced = TRUE?"}
    Balanced -- Yes --> Window["Window L..S; targets"]
    Balanced -- No --> Tasks["M grid and fixed seeds"]
    Window --> Tasks
    Tasks --> Support{"History supported?"}
    Support -- Yes --> Method{"Inference method?"}
    Support -- No --> Fail["Warning and NA bounds"]
    Method -- FLCI --> FLCI["Optimized folded normal"]
    Method -- Conditional or C-LF --> ARP["Gaussian test inversion"]
    FLCI --> Result["Components and status"]
    ARP --> Refine["Shared-budget refinement"]
    Refine --> Result
    Fail --> Result
    Result --> Output{"Which output?"}
    Output -- sets, summary --> Plot["Summary and paper plots"]
    Output -- FLCI weights --> Weights["v or gamma per element"]
    Weights --> Heat["Table or heatmap"]
    style CovType fill:#1e90ff,stroke:#1565c0,color:#fff
    style CR2 fill:#1e90ff,stroke:#1565c0,color:#fff
    style Given fill:#1e90ff,stroke:#1565c0,color:#fff
    style Default fill:#1e90ff,stroke:#1565c0,color:#fff
    style Balanced fill:#1e90ff,stroke:#1565c0,color:#fff
    style Window fill:#1e90ff,stroke:#1565c0,color:#fff
    style Output fill:#1e90ff,stroke:#1565c0,color:#fff
    style Weights fill:#1e90ff,stroke:#1565c0,color:#fff
    style Heat fill:#1e90ff,stroke:#1565c0,color:#fff
```

Blue stages changed in 0.2.0: the covariance option, the default procedure, the balanced option and the FLCI weights.
Before 0.2.0 a call without `method` under SD took the FLCI branch with the full-path restriction. The "Window" stage
restricts the targets to periods 1 to `S` and, depending on the framework and `cohort.paths`, forms the aggregate
path on `L..S` or truncates the cohort paths (see "Balanced detail"). The FLCI weights branch exists only for results
computed with `method = "FLCI"`; it reads a finished result and changes nothing in it (see "FLCI weights detail").

| Stage | Purpose | Key Dependencies | Changed in 0.2.0 |
| --- | --- | --- | --- |
| Input | Preserve observed outcomes and forward-fill treatment only; align named external quantities | Input/validation modules | No |
| Covariance | CR1 multiplies the cross-product of the cluster scores by the CR1 factor; CR2 scales each residual by its (cluster, cohort group) factor before the estimator's operator and applies no factor; external input keeps the supplied `vcov` | `covariance.R`, estimator operator | Yes (CR2 option; CR1 unchanged) |
| Representation | Carry the full covariance and deterministic normalization; retain target weights under aggregation | Mapping/targets | Balanced representations only (window paths; aggregate path on the window) |
| Defaults | `method = NULL` resolves to C-LF under SD and RM; under SD `full.path = NULL` resolves to TRUE only for FLCI, so the default SD scope is post-only | `sensitivity.R` | Yes |
| Balanced window | With `balanced = TRUE`: compute `L` and `S` from the coefficient labels; restrict the targets to periods 1 to `S`; aggregate inference on `L..S` with fixed cohort weights; under cohort inference keep the full paths or truncate them to the window | `sensitivity.R`, `targets.R`, restriction operators | Yes (new; default `balanced = FALSE` skips it) |
| Tasks | One task per target and bound `M`, each with a fixed seed; optional socket workers | `parallel.R` | No (workers load the package under its new name) |
| Inference | Apply SD or RM and the chosen procedure to each target/bound task | Adapters/private engine | Only under `cohort.paths = "window"`: the restriction uses the rows inside the window |
| Search accounting | Reuse prepared tests and calibration; cache evaluations across all attempts and audits | `utils-search.R`, SD/RM adapters, private engine/adapter | No |
| C-LF outer domain | Scan the necessary first-stage support using its certified finite bounds plus the endpoint probe margin | SD adapter and the LF support LP helper | No |
| Conditional search domain | Invariant initial contrast/scale, bilateral expansion and mandatory adjacent-domain agreements; no certified global support | `utils-conditional-search.R` and the shared point cache | No |
| Output | Keep each connected component and report failures without successful bounds; printed headers name the procedure, scope and covariance | S3 methods/plots | Covariance line in print/summary; for balanced results, the balanced-window line, settings fields, `aggregation` element and plot caption |
| FLCI weights | For FLCI results: rebuild each target's representation, reconstruct `v = v^P + D'gamma` (or take `gamma`) with estimates and contributions, by element or by cohort; tabulate or draw a heatmap against the bound | `weights.R`, `plot.R`, `targets.R` | Yes (new; reads results only) |

---

## Architectural Patterns

- **One representation for two input routes.** A `cohort_fit` holds named coefficients, covariance, bias mapping,
  normalization, coordinate metadata and cohort masses, whether it comes from `cohort_fit()` or `cohort_input()`.
- **Estimate once, vary assumptions.** Sensitivity tasks reuse the fitted object; they do not re-estimate the panel.
- **Defaults resolved in one place.** `cohort_sensitivity()` keeps the formals `method = NULL` and `full.path = NULL`
  and resolves them at dispatch: `method` becomes C-LF under both restrictions, and under SD `full.path` becomes
  TRUE exactly when the method is FLCI. The default SD analysis is therefore C-LF under the post-only restriction, an
  explicit `method = "FLCI"` still uses the full path, and a call that states `method` follows the same code path as
  in 0.1.1. `settings$method` and `settings$full.path` record the resolved values, and the printed header shows them
  (for example `cohort inference | SD | C-LF | post-only`).
- **The covariance type belongs to the fit.** `covariance` is an argument of `cohort_fit()` only. It follows `...`,
  so it matches only by its full name and no existing positional or abbreviated argument changes meaning. Everything
  downstream (conventional standard errors, representations, inference, fit plots) reads `fit$Sigma`, so a CR2 fit
  needs no other option. `cohort_input()` keeps the covariance it is given, and printing describes it as supplied.
  The printed label is derived from `covariance_convention` and the cluster count in `fit$diagnostics`.
- **CR2 as a row scaling before the linear operator.** Each estimator is linear in the outcomes, and the cluster
  scores are its operator applied to the residual-by-cluster matrix. CR2 multiplies every residual entry by the
  factor of its unit's (cluster, cohort group) pair before that operator. For a cluster nested in one cohort group
  this gives exactly `s_c / sqrt(1 - p_c)`, the rule of the paper; a cluster spanning several groups has each group's
  part of its score scaled by that group's factor. On a balanced panel the scores depend on the residuals only
  through cluster-cell sums, so this is the Bell-McCaffrey CR2 of the weighted cohort-period cell-mean model, with a
  working variance inversely proportional to the weight. The factors depend only on clusters, cohort groups and
  weights, not on the estimator, so the balanced Imp/CS-NYT equivalence also holds under CR2. The CR1 path performs
  exactly the floating-point operations it performed before the option existed.
- **Shares from the fixed unit registry.** Cohort groups are the adoption cohorts plus the never-treated units
  (cohort 0). A pair's share is its units' stored weight divided by its group's weight (unit counts when unweighted),
  on balanced and unbalanced panels alike, and is unchanged by the common rescaling of extreme weights. A pair that
  holds every unit of its group is detected by unit counts, not by a floating tolerance, and gets factor 0 without a
  square root. `1 - p` is computed as `rest / total`; `rest` is summed directly for the one pair per group, if any,
  that holds more than half of the group's weight, so no cancellation occurs.
- **Balanced options as opt-in additions.** `balanced = FALSE` and `cohort.paths = NULL` are appended after `control`,
  validated after every existing check, and with the defaults `window` stays `NULL`, so targets, representations and
  engine calls (`pre_window = NULL`) perform exactly the operations they performed before. The new `settings` fields
  and the `aggregation` element are added only when `balanced = TRUE`, so default results keep their element and
  settings names. The window is read from the coefficient labels, which makes it defined for `cohort_input()` objects
  and unbalanced panels too; in the minimum-wage application (cohorts 2004, 2006 and 2007 in 2001-2007) it runs from
  `L = -2` to `S = 1`.
- **Window paths keep every pre-treatment coefficient.** Under `cohort.paths = "window"` the post-treatment
  coefficients after `S` are dropped, and the restriction uses only rows inside the window: RM benchmarks are the
  first differences between periods `L` and 0, and the full-path SD restriction bounds only the pure pre-treatment
  second differences whose three periods are at or after `L`. Pre-treatment coefficients before `L` stay in the
  representation because the CS-NYT bias map and the imputation normalization use them, but they enter no
  restriction row. Dropping the later post-treatment coefficients is the same as leaving them unrestricted: they
  enter neither the target nor the bias of a kept coefficient, because in the bias maps that `cohort_fit()` builds an
  adjustment term always refers to a cell with a smaller relative period or to a pre-treatment cell
  (`.hc_check_window_paths` enforces this for external input). The engine needs them dropped rather than
  unrestricted, because the FLCI uses one post-connected row per post-treatment coefficient. Consequences: under
  post-only SD both path choices give the same plug-in identified set and the same FLCI, and C-LF sets can differ
  slightly because the least-favorable critical value is simulated in a different dimension; under RM the benchmark
  is the largest change inside the window, so the plug-in set is narrower or equal, and the union has
  `prod_g 2 * (number of window benchmark pairs of g)` signed components (64 instead of 320 in the minimum-wage
  application); under full-path SD the pure pre-treatment rows outside the window are no longer checked, so a
  plug-in set can be nonempty where the full path gives an empty one. On a balanced panel imputation and CS-NYT give
  the same balanced cohort inference under either choice. A truncation of the post-treatment horizon alone would
  reproduce the `"full"` identified sets for every restriction, which is why `"window"` truncates both ends.
- **Balanced aggregation reuses the aggregate construction.** `.hc_representation` keeps its weight arithmetic,
  representability check (`qr.solve` with the 1e-9 check), null-space embedding of the normalization and names, and
  only replaces the period set by `L..S`. Every cohort is observed at every window period, so each row of `A` holds
  the same cohort weights `m_g / sum m`, and the aggregate target `l_A` is unique (the ATT puts `1/S` on each window
  post period; the dynamic target at `s` is the unit vector). The aggregate normalization is the left null space of
  `A E`, with `E` a basis of the null space of the fit's normalization: for imputation it has no row unless the window
  covers every cohort's whole pre-treatment path; for CS-NYT and CS-NT it is the `s = 0` row. The aggregated
  coefficients and covariance are `A betahat` and the symmetrized `A Sigma A'`.
- **FLCI weights as a read-only view.** `flci_weights()` and the weights plot add no element to results and change
  neither `cohort_sensitivity()` nor the engine. The extractor recomputes, for each target, the representation and
  its second-difference operators from the result's own settings, computes the post-only weights and center in
  closed form (one triangular solve and a few products), and combines them with the stored `gamma`. It reports
  weights only after checking them against what the FLCI stored (the length and names of `gamma`, and the
  projection of `v` compared with the stored coefficient weights); a result it cannot reproduce gives an error, not
  a guess.
- **The paper's representative of `v`.** Each cohort's pre-treatment coefficients satisfy one exact normalization
  (CS: the coefficient at `s = 0` is zero; imputation: they sum to zero), so the weights on them are unique only up
  to that direction. The FLCI stores the representative orthogonal to the normalization; `flci_weights()` reports
  the paper's form `v = v^P + D'gamma`. In that form the post-only weights sit only on each cohort's last two
  pre-treatment coefficients (`s = -1` and `s = 0`), plus, under CS-NYT, the cells of the base-period correction
  (imputation and CS-NYT differ only in those cells); they are constant in `M`, so the full-path and post-only
  heatmaps can be compared cell by cell. The estimate, the center, `gamma` and each cohort's sum of contributions
  are the same for every representative; the per-cell weights and the cohort sums `weight`, `l1` and `share` of `v`
  depend on it. Under imputation `D_PP' gamma` already has zero cohort sums, so the full-path and post-only stored
  weights differ by exactly `D_PP' gamma`.
- **Heatmaps follow the paper's figure conventions.** Columns are the computed bounds as a discrete axis, rows the
  relative periods with the earliest at the top, cohorts in adoption order from the top, one panel per target. The
  fill is diverging and symmetric around zero, blue (`#0072B2`) to white to vermilion (`#D55E00`) from
  `paper_palette("representation")`, so exact zeros are white. Captions state what is shown, the SD scope and the
  framework, and add the tie-break, all-zero, failed-result, unbalanced-panel and balanced-window sentences when
  they apply.
- **Type dispatch keeps the set plot.** The five new formals of `plot.cohort_sensitivity` follow `color`, so
  positional calls are unaffected; `type` is matched first, and each type rejects the other type's arguments when
  they are supplied explicitly. With `type = "sets"` the earlier statements run unchanged.
- **Actual histories in estimation.** Sparse operators encode observed comparisons. Unbalanced samples keep a warning
  because their fixed-mass bias map is an approximation.
- **Explicit numerical support.** Singular covariance alone is allowed. Unsupported deterministic Conditional/C-LF
  cases have failure statuses and no public bounds; FLCI uses its defined zero-noise limit.
- **Reproducible task boundaries.** Seeds are assigned before dispatch, so results do not depend on the worker count.
  The public `parallel`/`cores` settings govern both pools: target/bound tasks run on socket workers, and only a
  single RM task passes `cores` on to its component workers, so worker pools are never nested. Clusters are closed
  on exit.
- **Self-reference by name.** Workers receive the parent's library paths and load the package by name, and the
  engine is found through the package namespace path (`upstream/R` when installed, `inst/upstream/R` when
  source-loaded). The seven literals listed under "Package-name lookups" are the only places the code names its own
  package.
- **Private numerical resources.** `.hc_backend()` sources the five engine files into one cached environment whose
  parent is the base environment, so nothing is attached to the user's session; `sd_inference_load_dependencies()`
  keeps a second private environment for the FLCI core and the HonestDiD routines. Upstream notices remain under
  `inst/upstream/`.
- **One search budget per target/bound.** The shared cache counts distinct confidence-set point evaluations across
  original searches, endpoint/topology/sentinel probes and refinements. A separate sentinel context distinguishes
  ungated RM checks from regular support-gated tests. Calibration and support-optimization LP work are outside this
  counter; it does not bound elapsed time.
- **Bounded refinement with preserved checks.** `max.refinements` defaults to 3; 0 disables retries. Eligible
  unresolved statuses (`open_search_boundary`, `topology_audit_failure`, `endpoint_audit_failure`,
  `uncertified_empty_search`, `topology_miss`) double the regular grid intervals, reusing the prepared problems,
  calibration draws and point cache. Input, unsupported-deterministic and solver failures do not trigger retries.
  Tolerances and certification checks stay fixed; a final failure warns and reports no set.
- **Spacing independent of anchors.** Expansion uses the registered regular-grid spacing, including the
  corresponding half spacing for finer audits. A profile anchor nearly equal to a grid point cannot collapse the
  expansion step.
- **Conditional finite domains.** The initial domain is built from the mapped primitive-trend contrast and its
  sampling scale, with a first-difference fallback. Concentric bilateral expansions keep the mesh spacing and require
  two successive adjacent-domain agreements, plus endpoint checks and optional interior mesh audits. The Conditional
  driver uses `.hc_search_budget` and runs its own refinement loop with the same eligible statuses and doubling rule
  as `.hc_adaptive_search`. A finer retry revisits the widest previously evaluated domain and checks all cached
  accepted and rejected points before success. `max.expansions` and `max.refinements` keep separate roles; neither
  resets the cumulative budget.
- **Scheduled and evaluated ranges.** Conditional `search.domain` is the range of coordinates actually evaluated,
  including endpoint probes; it is empty before any evaluation. `requested.domain` and the scheduled history endpoints
  identify blocked requests. `domain.history` records whether and where each domain was evaluated, its regular
  spacing, local status, agreements and cumulative evaluations. Print/summary show compact per-target/M ranges, and
  plot captions identify bounded numerical search. No Conditional result claims global support or exhaustive
  topology.
- **Independent regular-cell audit phases.** The Conditional interior audit also inverts one combined grid at the
  quarter and three-quarter points of every registered regular cell, independently of profile-anchor subintervals.
  The whole component union must agree with the base mesh; disagreement enters the same refinement/no-set path under
  the same cache and budget. Remembered acceptances and rejections constrain the final answer. Finite meshes can still
  miss arbitrarily narrow interior components or gaps.
- **Initial range and reference density.** Only Conditional searches with `grid.radius > 8` (and at least two
  expansions) call `.hc_conditional_density`. It keeps the requested endpoints and raises the initial interval count
  so the spacing is no coarser than the radius-8 reference spacing, including endpoint rounding of the initial
  domain. This is an upper bound on spacing, not exact lattice nesting across calls. Radius 8 and smaller use the
  literal grid count and skip the reference arithmetic. The effective initial count is not a retry; refinement
  retries and domain doublings multiply it under the same cumulative cap. Diagnostics distinguish requested and
  effective counts and the reference radius/spacing; a count limit records the requested domain with no fictitious
  evaluations.
- **C-LF support is a necessary gate of the same test.** SD C-LF uses the selected-moment LP support calculation with
  the realized LF cutoff. It scans both outward-certified endpoints plus the probe margin and disables
  endpoint-triggered tail expansion on this bounded path. Calibration and support preparation happen once, outside
  retries. RM keeps its component support implementation; plain Conditional does not use the LF gate.
- **Certification has separate meanings.** `support.certified` certifies outer support, not exhaustive interior
  topology. Nonempty finite-grid results keep `topology.certified = FALSE` and `inversion.exact = FALSE`. A
  no-feasible-support result is returned as empty only when the global profile dual bound independently certifies LF
  rejection; other unresolved or unbounded support outcomes return a failure and no set.

---

## Notes

- **Covariance objects.** CR1 fits have exactly the elements of version 0.1.1: `covariance_convention` is
  `"cohort_period_centered_CR1"`, `cluster_scores` holds the raw scores, and `diagnostics$cr1` is the CR1 factor.
  CR2 fits have `covariance_convention = "cohort_period_centered_CR2"`, adjusted `cluster_scores` whose symmetrized
  cross-product is `Sigma`, `diagnostics$cr1 = NA`, and a final element `cluster_shares` with columns `cluster`,
  `cohort` (0 = never treated), `units`, `mass` (stored weight units), `share` and `factor`: one row per cluster and
  cohort group, ordered by cluster as in `rownames(cluster_scores)` and then by cohort, with share 1 and factor 0 for
  a cluster that holds a whole group. Summary objects of fits and sensitivity results gain a final `covariance`
  element; the covariance option leaves the structure of sensitivity objects and the plots unchanged.
- **Balanced result objects.** With `balanced = TRUE`, `settings` ends with `balanced = TRUE`, `cohort.paths`
  (`"full"` or `"window"`; `NA_character_` under aggregate inference) and `window` (named integer
  `c(lower = L, upper = S)`). Aggregate results also end with `aggregation`, a list of `aggregation_matrix` (rows
  `delta_pre_g1_s<s>` for `s <= 0` and `delta_post_g1_s<s>` for `s >= 1`, `s = L..S`; columns the fit's coefficient
  names), `betahat`, `Sigma`, `normalization` (possibly with no row), `coordinate_info` and `targets` (a matrix with
  one row per target id and the aggregate coordinates as columns). Target ids are unchanged; only the ATT label
  names the periods. `targets[[j]]$target` keeps one weight per coefficient of the fit, and the conventional
  estimates use it with the fit's `betahat` and `Sigma` under both frameworks and both path choices. Example of the
  printed line: `Balanced window: relative periods -2 to 1; target period 1; cohort paths: full`.
- **Balanced messages.** After all earlier validation, in this order: `balanced must be TRUE or FALSE.`;
  `cohort.paths applies only when balanced = TRUE.`; `cohort.paths applies only to framework = 'cohort'.`;
  `cohort.paths must be 'full' or 'window'.`; `balanced = TRUE requires relative period 1 for every cohort: no
  post-treatment period is observed for all cohorts.`; `balanced = TRUE requires relative period 0 for every cohort:
  no pre-treatment period is observed for all cohorts.`; and `cohort.paths = 'window' requires that the bias map W
  and the normalization link no coefficient inside the balanced window to a coefficient after it.` The last three
  can occur only for `cohort_input()` objects. A cohort with fewer than two pre-treatment periods keeps the per-task
  status `insufficient_prehistory`; when every cohort has at least two, `L <= -1`.
- **Tests of the balanced options.** `tests/testthat/test-balanced.R` (11 blocks) computes its references from the
  data and the fits' public fields: the window and settings on the helper panels, the balanced target weights and
  estimates, the `aggregation` element against a matrix built from the cohort masses, equal post-only plug-in sets
  and FLCI under both path choices, the RM component counts under window paths, imputation against CS-NYT under both
  path choices, a single-cohort fit with and without `balanced`, a `cohort_input()` round trip, every message above,
  the printed lines and the plot captions and points. The earlier test files are unchanged and run with the default
  `balanced = FALSE`.
- **What balanced aggregation does not remove.** It fixes the cohorts averaged at each period but changes the
  estimand (only periods 1 to `S`), and it does not remove the change in the not-yet-treated comparison groups, which
  enters the imputation and CS-NYT biases through `W`; balanced aggregate and cohort inference can therefore still
  differ. Balanced aggregate imputation and CS-NYT results also differ in general, although they coincide when every
  cohort's window ends before a later cohort adopts (as in the minimum-wage application, where `S = 1`): each
  cohort's difference between the two estimators' coefficients is then constant across the window, and SD and RM
  do not respond to a constant shift of the aggregated path. The help page and the detailed guide therefore say
  "in general".
- **FLCI weights objects.** With `by = "element"`, `flci_weights()` returns one row per task with weights and
  element, with columns `target` and `label` (character), `M` (double), `cohort` (integer; `NA` under aggregate
  inference), `relative_period` (integer: the coefficient's period, or the last period `e` of the second
  difference), `element` (character: the coefficient name such as `delta_pre_g4_s-2`, or `d2_g<g>_s<e>`, with
  `g1` under aggregate inference), `weight`, `estimate` (`betahat` or `zhat`), `contribution` (doubles) and
  `tie.break` (logical). Coefficient rows cover every pre-treatment coefficient of the representation, including
  each cohort's normalized cell (CS: `s = 0`, whose coefficient and contribution are zero) and, under
  `cohort.paths = "window"`, the coefficients before `L`, on which CS-NYT weights can load; `gamma` rows are the `J`
  pure pre-treatment second differences, all zero for a post-only FLCI. With `by = "cohort"` (cohort inference only)
  the columns are `target`, `label`, `M`, `cohort`, `elements`, `weight`, `l1`, `share`, `contribution` and
  `tie.break`. Rows follow the targets, then the sorted bounds, then cohort and period. Attributes: `weights`, `by`,
  `framework`, `full.path`, `centers` (`target`, `label`, `M`, `estimate`, `post.only.center`, `center`; one row per
  task with weights) and `missing` (`target`, `label`, `M`, `status`; one row per task without weights). When
  nothing qualifies (every task failed, or `gamma` with `J = 0`) the data frame has 0 rows with the same columns and
  types. `tie.break` is `TRUE` when the task's `diagnostics$tie.break` is `"minimum_euclidean_norm_gamma"`: several
  weight vectors give the same FLCI length (for example at `M = 0` when the covariance of `zhat` is singular), and
  the FLCI took the one with the smallest Euclidean norm of `gamma`. In the minimum-wage application (imputation,
  ATT) the post-only weights are 0.699 and -0.979 at `s = -1` and `s = 0` of cohort 2004 (0.489 and -0.815 for
  2006, 0.476 and -0.951 for 2007) and zero elsewhere; the full-path `gamma` is positive on all eight second
  differences at `M = 0`, nonzero only for 2004's at `M = 0.005` and `0.01`, and zero from `M = 0.015`, where the
  full-path and post-only FLCIs coincide.
- **FLCI weights plot.** The plot data are the selected rows of `flci_weights()` plus `value` (`weight` by element,
  `contribution` by cohort), `column` (the bound labels `format(M, trim = TRUE, drop0trailing = TRUE)`, with 15
  digits if those labels collide), `row` (relative periods or cohort labels, leveled so that the earliest is at the
  top), `cohort_label` (cohort inference), `target_label`, `text` (`formatC(value, format = "fg", digits = 2)`, empty
  for zeros) and `dark` (`|value| > lim / 2`, drawn in white). `lim = max|value|`, or 1 when that is 0. Element plots
  under cohort inference use `facet_grid(cohort_label ~ target_label)` with free row heights and the cohort strips
  on the left; the others have one row of target panels. Defaults: x label `M`; y label `Relative period`,
  `Last period of the second difference` or `Cohort`; legend `Weight v`, `Weight gamma` or `Contribution`. The
  caption starts with `<what>; <scope>; <framework> inference.`, where the scope is `SD, full path` or
  `SD, post-only` (for example `FLCI weights v on the pre-treatment coefficients; SD, full path; cohort
  inference.`), then `White cells are zero weights.` (element) or `White cells are zero. The cohort sums add up to
  the FLCI center minus the conventional estimate.` (`... minus the post-only center.` for `gamma`), and when they
  apply `At M = <bounds>, several weight vectors give the same FLCI length; the one with the smallest Euclidean
  norm of gamma is shown.`, `All weights gamma are zero at the selected bounds: the full-path FLCI equals the
  post-only FLCI.`, `<n> failed result(s) without weights are not shown; inspect summary().`, `Unbalanced-panel
  approximation.` and the balanced-window sentence of the set plot; it is wrapped at 95 characters.
- **FLCI weights messages.** `flci_weights()`, in this order: `x must be a cohort_sensitivity object.`; the
  `match.arg()` errors for `weights` and `by`; `FLCI weights exist only for results with method = 'FLCI'; this
  result used Conditional.` (or `C-LF.`); `by = 'cohort' requires cohort inference: aggregate weights have no
  cohorts.`; and `flci_weights() could not reconstruct the FLCI weights of this result.` for a malformed object,
  operators that cannot be built, or a failed check. `plot(x, type = "weights")`: the `match.arg()` error for
  `type`; `identified.set and color apply only to type = 'sets'.` (with `type = "sets"` instead: `weights, by,
  targets, and labels apply only to type = 'weights'.`); the `match.arg()` errors for `weights` and `by`; `labels
  must be TRUE or FALSE.`; `targets must be NULL or a character vector of target ids.`; the method and
  `by = 'cohort'` errors above; `cohorts cannot be selected under aggregate inference: aggregate weights have no
  cohorts.`; for `weights = "gamma"` on a post-only result the message `The post-only FLCI puts no weight on the
  pure pre-treatment second differences: its weights gamma are zero by construction. Use weights = 'coefficients'
  to plot its weights v.` and an invisible `NULL`; `M must select finite computed bounds.`; `No computed results
  match the requested plot selection.`; `No FLCI weights to plot: every selected result failed; inspect
  summary().`; for `weights = "gamma"` without pure pre-treatment second differences the message `This result has
  no pure pre-treatment second differences, so the FLCI has no weights gamma to plot.` and an invisible `NULL`; and
  the no-match error again when `cohorts` removes every row. `flci_weights()` itself never messages or warns.
- **Tests of the FLCI weights.** `tests/testthat/test-flci-weights.R` (14 blocks; its helpers, prefixed `fw_`, are
  defined in the file) computes its references from public fields of fits and results (`betahat`, `W`,
  `normalization`, `coordinate_info`, `cohort_masses`, `targets`, `results`) without internal functions: columns,
  types, order and attributes; `v` against the closed form and its projection against the stored weights; `gamma`
  against the stored values and their labels; the two sum identities and the post-only center against a separate
  post-only FLCI; cohort sums; the CS-NYT correction cells; aggregate and balanced (window) results; post-only
  results and paths too short for second differences; failed tasks listed as missing; the minimum-norm tie-break
  (a `cohort_input()` object whose covariance is zero on the pre-treatment block); every error and message; the
  unchanged default plot; and the weights plot's tiles, white zeros, captions and selections.
- **Covariance messages.** `covariance must be 'CR1' or 'CR2'.` for any other value; `At least two independent
  clusters are required for CR1 covariance.` (or `CR2`); and, for CR2, an error if a cluster holds all but a
  numerically negligible share of its cohort's weight. `cov = "CR2"` is not a partial match: it lands in `...` and is
  rejected as an unrecognized argument.
- **Whole-cohort clusters.** On a balanced panel, a cluster that consists of one whole cohort has only single-cluster
  cells, so the zeroing in the compressed path makes its score exactly zero under both covariances, and the cohort's
  sampling variation is omitted from the covariance; if the cluster also holds units of other cohort groups, only that
  cohort's part of its score is zero (CR2 factor 0 for that pair). Under CR2 that pair's factor is 0, so the cohort's
  part of the score is zero on any panel. On an unbalanced panel with missing outcomes in that cohort, the CR1 score
  of a cluster that consists of the whole cohort is generally not zero: the estimator is then no longer a function
  of the cell means alone, and the cluster's centered residuals do not cancel in it.
  `tests/testthat/test-covariance-cr2.R` pins both behaviors on an unbalanced panel. The fit warns about cohorts with
  fewer than 10 clusters and about single-cluster cells.
- **Renaming or forking.** A package with a different name must change the seven literals in `R/parallel.R`,
  `R/inference_rm.R`, `R/utils-engine.R` and `inst/upstream/R/inference_engine.R`, as well as `tests/testthat.R` and
  `tests/testthat/test-session-state.R`. A stale name still lets the package load, but socket workers and the engine
  path then fail (or pick up another installed copy). The two-worker test in `test-session-state.R` is skipped when
  `system.file(package = "CohortHonestDiD")` is empty, so a stale name there would skip it silently; a full test run
  should report no skipped tests.
- **Weights.** Extreme panel weights are divided by a common positive factor only when needed for finite arithmetic;
  `data_info$weight_scale` records it, and ordinary weights keep divisor 1. CR2 shares are ratios of stored weights
  and therefore do not depend on this divisor.
- **Diagnostics.** For Conditional/C-LF searches, result diagnostics expose `search.attempts` with `attempt`,
  `grid.points`, `new.evaluations`, `cumulative.evaluations`, `success` and `status`, plus `final.grid.points`,
  `refinements` and cumulative `test.evaluations`. SD C-LF also exposes `support.certified`, `support` and
  `search.domain`; the globally certified empty branch records `empty.certified = TRUE` with zero point evaluations.
  `settings$control` keeps the control list as supplied, without the filled-in defaults.
- **Engine layout.** The adapter lives in `inst/upstream/R/inference_adapter.R`; there is no `R/inference_adapter.R`.
  `arp-nuisance.R` keeps its original HonestDiD author notice; see `inst/upstream/NOTICE` for attribution.
- **Scope of the theory.** The balanced-panel coverage and estimator-equivalence results of the paper do not
  automatically extend to the unbalanced-panel approximation or to unsupported deterministic designs. On unbalanced
  panels the CR2 adjustment, like the bias map, is an approximation: the shares are fixed by unit and do not follow
  the period-by-period composition of the cells, and the estimator is no longer a function of the cell means alone.
  On unbalanced panels the balanced window is read from the coefficient labels and uses the fixed cohort masses, so
  it shares the same approximation.
