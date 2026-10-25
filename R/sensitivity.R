#' Cohort or aggregate sensitivity analysis
#'
#' Run robust inference without re-estimating the cohort-period coefficients.
#' Relative period one is the first treated period; zero is the last untreated
#' period. Each reported confidence set has pointwise, rather than simultaneous,
#' coverage. On unbalanced panels the balanced-design mapping is approximate.
#'
#' @param fit A fitted object returned by [cohort_fit()] or [cohort_input()].
#' @param restriction Either `"SD"` (smoothness) or `"RM"` (relative magnitudes).
#' @param M Nonnegative sensitivity bounds. SD bounds are in outcome units;
#'   RM bounds are dimensionless multiples of each cohort's pre-period benchmark.
#' @param target `"att"`, `"dynamic"`, or `"cohort-dynamic"`.
#' @param framework `"cohort"` (default) or `"aggregate"`.
#' @param method `"FLCI"`, `"Conditional"`, or `"C-LF"`. The default, `NULL`,
#'   is C-LF under both SD and RM, as the paper recommends. FLCI is available
#'   only under SD.
#' @param full.path For SD, restrict pure pre-treatment curvature as well?
#'   The default, `NULL`, is TRUE for FLCI and FALSE (post-only) for
#'   Conditional and C-LF, so the default SD analysis is C-LF under the
#'   post-only restriction. Not applicable to RM.
#' @param conf.level Pointwise confidence level in (0.5, 1), default 0.95.
#' @param parallel Logical; run independent inference tasks in parallel?
#' @param cores Number of workers. Defaults to detected logical cores minus one,
#'   with a minimum of one. Workers are capped by the number of independent
#'   target/bound tasks, or signed components for a single RM task.
#' @param seed Optional integer seed. Fixed seeds give results independent of
#'   worker count and preserve the caller's random-number state.
#' @param control Named list of numerical controls; see the section
#'   "Numerical controls" of \code{vignette("detailed-guide", package = "CohortHonestDiD")}.
#'   `max.refinements` is a nonnegative integer (default 3; 0 disables retries).
#'   Eligible unresolved searches double the regular grid intervals per retry.
#'   `max.evaluations` is the cumulative limit on distinct confidence-set point
#'   evaluations per target/bound, including endpoint probes, topology audits,
#'   and retries. Simulated least-favorable critical values are reused.
#'   For Conditional inference, `max.expansions` (default 6) caps bilateral
#'   domain doublings, including the two required adjacent-domain agreements.
#'   Values 0 or 1 return `domain_stability_limit` and no set. Refinement
#'   restarts the same domains at a finer mesh within the shared budget.
#'   For Conditional `grid.radius > 8`, `grid.points` specifies the reference
#'   density at radius 8: the effective initial count grows enough to keep the
#'   initial spacing no coarser than that reference. At `grid.radius <= 8`,
#'   `grid.points` remains the literal initial count. Larger ranges require
#'   more evaluations and may exhaust the unchanged cumulative budget.
#' @param balanced Logical. If \code{TRUE}, the target and, under aggregate
#'   inference, the aggregated event study use only the relative periods
#'   observed for every cohort, the balanced window. Default \code{FALSE}. See
#'   the section "Balanced targets".
#' @param cohort.paths Only with \code{balanced = TRUE} under cohort inference:
#'   \code{"full"} (used when \code{NULL}) restricts each cohort's whole path,
#'   and \code{"window"} restricts each cohort's path only inside the balanced
#'   window. Supplying it otherwise is an error.
#' @return A `cohort_sensitivity` object containing `intervals` (one row per
#'   connected confidence-set component), `estimates`, `settings`, and detailed
#'   `results`. Empty sets and numerical failures have explicit status rows.
#'   With \code{balanced = TRUE}, see the section "Balanced targets" for the
#'   added settings and the \code{aggregation} element.
#' @details The defaults follow the paper's recommendation: cohort inference
#'   with the C-LF procedure, under the post-only restriction for SD. Before
#'   version 0.2.0 the SD default was the full-path FLCI; request it with
#'   `method = "FLCI"`. The FLCI needs the stronger full-path restriction and
#'   should be used with caution when some cohorts contain only a few units.
#'   Inference uses the covariance stored in `fit`; use
#'   `cohort_fit(covariance = "CR2")` for the bias-reduced CR2 covariance.
#'
#'   Conditional and C-LF procedures return `success = FALSE`, a warning,
#'   and missing confidence-set bounds when the target is deterministic or a
#'   required moment has zero or numerically degenerate variance. The status is
#'   `unsupported_zero_variance` or `unsupported_moment_variance`. Singular
#'   covariance alone does not trigger this limitation. FLCI supports its exact
#'   zero-covariance limit.
#'
#'   Search diagnostics include `search.attempts` (attempt, grid.points,
#'   new.evaluations, cumulative.evaluations, success, status),
#'   `final.grid.points`, `refinements`, and `test.evaluations`. The requested
#'   grid size remains in settings. Refinement preserves endpoint tolerances
#'   and all certification checks; exhaustion returns a warning and no set.
#'   C-LF searches use certified outer bounds from their own first-stage
#'   acceptance inequalities. Under SD, `support.certified`, `support`, and
#'   `search.domain` record these bounds and the searched domain. Finite-grid
#'   agreement does not certify exhaustive interior topology.
#'
#'   Conditional inference searches concentric finite domains centered at a
#'   primitive-trend contrast, with a normalization-invariant sampling scale.
#'   The initial radius adds the plugin radius and `grid.radius` times this
#'   scale. Wider domains retain the regular mesh spacing. Two successive
#'   adjacent-domain agreements are required in addition to endpoint checks.
#'   `topology.audit = FALSE` skips interior mesh checks only; mandatory domain
#'   and endpoint checks remain. The automatic endpoint tolerance still uses
#'   the conventional target standard error, separately from the search scale.
#'   Diagnostics include `search.center`, `search.scale`, `search.scale.source`,
#'   `initial.domain`, `search.domain`, `requested.domain`, `domain.stable`,
#'   `domain.expansions`, `agreement.streak`, and `domain.history`.
#'   `search.domain` records evaluated coordinates; it is empty before any
#'   evaluation. `requested.domain` records a range blocked by an allocation
#'   or evaluation limit. Print and summary methods report searched ranges.
#'   Bounded numerical search does not rule out components outside those
#'   ranges; `support.certified`, `topology.certified`, and `inversion.exact`
#'   remain FALSE for Conditional results. Limits return a warning and no set.
#'   With `topology.audit = TRUE`, an additional shifted finer grid checks the
#'   quarter and three-quarter points of registered regular cells, independently
#'   of inserted anchors. Disagreement uses the same bounded refinement policy.
#'   Finite range and mesh checks can still miss arbitrarily narrow interior
#'   components or gaps, as well as components outside the searched range.
#'   Density diagnostics are `initial.grid.points`, `grid.density.adjusted`,
#'   `reference.domain.radius`, and `reference.grid.spacing`. The reference
#'   radius adds the same plugin radius to eight times the search scale.
#'   Integer ceiling and any necessary upward rounding preserve an initial
#'   spacing upper bound, not exact lattice nesting across separate radius
#'   requests. Settings retain the requested count; histories and
#'   `final.grid.points` record the effective counts. Starting at a larger
#'   count does not consume a refinement. Reference fields are NA when no
#'   wider-range calculation is needed or reached.
#' @section Balanced targets:
#' The balanced window runs from relative period \code{L} to \code{S}. \code{S}
#' is the last post-treatment period and \code{L} the earliest pre-treatment
#' period that every cohort in \code{fit} has. Both are read from the
#' coefficient labels, so the window is also defined for unbalanced panels and
#' for \code{\link{cohort_input}} objects. With \code{balanced = TRUE} the ATT
#' averages the cells of periods 1 to \code{S}, each cohort weighted by its mass
#' in every period, and dynamic and cohort-dynamic targets are reported for
#' periods 1 to \code{S} only. Under aggregate inference the aggregated event
#' study is formed on periods \code{L} to \code{S}, with the same cohort weights
#' in every period (balanced aggregation). This removes the change in the
#' cohorts averaged at each period, but it changes the estimand. It does not
#' remove the change in the not-yet-treated comparison groups, which enters the
#' imputation and CS-NYT biases through the bias map \code{W}, so balanced
#' aggregate and cohort inference can still differ.
#'
#' Under cohort inference \code{cohort.paths} sets the cohort paths that the
#' restriction uses. \code{"full"} keeps every cohort's whole path, as without
#' \code{balanced}. \code{"window"} uses each path only inside the window:
#' coefficients after period \code{S} are dropped, the RM benchmark is a
#' cohort's largest change between periods \code{L} and 0, and the full-path SD
#' restriction bounds only the pre-treatment second differences inside the
#' window. Pre-treatment coefficients before \code{L} are kept, because the bias
#' map and the imputation normalization use them. Under the post-only SD
#' restriction both choices give the same plug-in identified set; their
#' confidence sets can differ slightly, for example through the simulated C-LF
#' critical value. On a balanced panel, imputation and CS-NYT give the same
#' balanced cohort inference under either choice, but in general not the same
#' aggregate inference.
#'
#' With \code{balanced = TRUE}, \code{settings} also contains \code{balanced},
#' \code{cohort.paths} (\code{NA} under aggregate inference) and \code{window}
#' (\code{lower} = \code{L}, \code{upper} = \code{S}). Aggregate results also
#' contain \code{aggregation}: the aggregation matrix
#' (\code{aggregation_matrix}), the aggregated coefficients (\code{betahat}) and
#' their covariance (\code{Sigma}), the aggregate \code{normalization},
#' \code{coordinate_info}, and the aggregate target weights (\code{targets}, one
#' row per target).
#' @export
#' @examples
#' data(minimum_wage)
#' \donttest{
#' fit <- cohort_fit(Y ~ D, minimum_wage, index = c("id", "year"))
#' # Default: C-LF under the post-only SD restriction.
#' result <- cohort_sensitivity(fit, M = c(0, 0.01), seed = 1)
#' plot(result)
#' # Balanced aggregation on the relative periods observed for every cohort.
#' balanced <- cohort_sensitivity(fit, M = 0.01, framework = "aggregate",
#'                                balanced = TRUE, seed = 1)
#' balanced$settings$window
#' }
cohort_sensitivity <- function(fit, restriction = "SD", M,
                               target = "att", framework = "cohort",
                               method = NULL, full.path = NULL, conf.level = 0.95,
                               parallel = FALSE, cores = NULL, seed = NULL,
                               control = list(), balanced = FALSE, cohort.paths = NULL) {
  if (!inherits(fit, "cohort_fit")) stop("fit must be a cohort_fit object.", call. = FALSE)
  for (key in c("restriction", "target", "framework")) {
    value <- get(key)
    .hc_assert(is.character(value) && length(value) == 1L && !is.na(value),
               paste0(key, " must be one character string."))
  }
  restriction <- match.arg(toupper(restriction), c("SD", "RM"))
  target <- match.arg(target, c("att", "dynamic", "cohort-dynamic"))
  framework <- match.arg(framework, c("cohort", "aggregate"))
  if (target == "cohort-dynamic" && framework == "aggregate") {
    stop("target = 'cohort-dynamic' requires framework = 'cohort'.", call. = FALSE)
  }
  if (missing(M) || !is.numeric(M) || !length(M) || any(!is.finite(M)) || any(M < 0)) {
    stop("M must be a nonempty vector of finite nonnegative bounds.", call. = FALSE)
  }
  M <- sort(unique(as.numeric(M)))
  if (is.null(method)) method <- "c_lf"
  .hc_assert(is.character(method) && length(method) == 1L && !is.na(method),
             "method must be one character string.")
  key <- tolower(gsub("[-_ ]", "", method))
  if (length(key) != 1L || !key %in% c("flci", "conditional", "clf")) {
    stop("method must be FLCI, Conditional, or C-LF.", call. = FALSE)
  }
  method <- c(flci = "flci", conditional = "conditional", clf = "c_lf")[[key]]
  if (restriction == "RM" && method == "flci") stop("FLCI is available only under SD.", call. = FALSE)
  if (restriction == "RM" && !is.null(full.path)) stop("full.path applies only to SD.", call. = FALSE)
  if (restriction == "SD" && is.null(full.path)) full.path <- method == "flci"
  if (restriction == "SD" && (!is.logical(full.path) || length(full.path) != 1L || is.na(full.path))) {
    stop("full.path must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L || !is.finite(conf.level) || conf.level <= 0.5 || conf.level >= 1) {
    stop("conf.level must lie strictly between 0.5 and one.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) stop("parallel must be TRUE or FALSE.", call. = FALSE)
  if (is.null(cores)) cores <- .hc_available_cores()
  cores <- .hc_integer_scalar(cores, "cores")
  .hc_engine_controls(control)
  if (!is.logical(balanced) || length(balanced) != 1L || is.na(balanced))
    stop("balanced must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(cohort.paths)) {
    if (!balanced) stop("cohort.paths applies only when balanced = TRUE.", call. = FALSE)
    if (framework == "aggregate") stop("cohort.paths applies only to framework = 'cohort'.", call. = FALSE)
    if (!is.character(cohort.paths) || length(cohort.paths) != 1L || is.na(cohort.paths) ||
        !cohort.paths %in% c("full", "window"))
      stop("cohort.paths must be 'full' or 'window'.", call. = FALSE)
  }
  window <- NULL
  if (balanced) {
    if (framework == "cohort" && is.null(cohort.paths)) cohort.paths <- "full"
    window <- .hc_balanced_window(fit)
    if (identical(cohort.paths, "window")) .hc_check_window_paths(fit, window)
  }
  paths <- if (is.null(cohort.paths)) "full" else cohort.paths
  targets <- .hc_targets(fit, target, window)
  tasks <- list(); task_info <- list(); estimates <- list(); aggregate.targets <- list()
  for (j in seq_along(targets)) {
    one <- targets[[j]]
    rep <- .hc_representation(fit, one$target, framework, window, paths)
    if (balanced && framework == "aggregate") aggregate.targets[[j]] <- rep$target
    theta <- sum(one$target * fit$betahat)
    se <- sqrt(max(0, as.numeric(t(one$target) %*% fit$Sigma %*% one$target)))
    q <- stats::qnorm((1 + conf.level) / 2)
    estimates[[j]] <- data.frame(target = one$id, label = one$label, cohort = one$cohort,
                                 relative_period = one$relative_period, estimate = theta, se = se,
                                 lower = theta - q * se, upper = theta + q * se)
    for (m in M) {
      k <- length(tasks) + 1L
      tasks[[k]] <- list(representation = rep, restriction = restriction, M = m,
                        method = method, full.path = full.path, alpha = 1 - conf.level,
                        seed = 0L, control = control)
      task_info[[k]] <- list(target = one$id, label = one$label, cohort = one$cohort,
                             relative_period = one$relative_period, M = m)
    }
  }
  seeds <- .hc_task_seeds(length(tasks), seed)
  for (k in seq_along(tasks)) tasks[[k]]$seed <- seeds[k]
  # The public switch controls all worker pools; never nest component workers
  # inside the outer target/bound pool.
  for (k in seq_along(tasks)) tasks[[k]]$control$component.cores <- 1L
  if (parallel && length(tasks) == 1L && restriction == "RM") {
    tasks[[1L]]$control$component.cores <- as.integer(cores)
  }
  results <- .hc_run_tasks(tasks, parallel, as.integer(cores))
  used.cores <- if (parallel) min(cores, length(tasks)) else 1L
  if (parallel && length(tasks) == 1L && restriction == "RM") {
    actual <- results[[1L]]$diagnostics$component.cores
    if (length(actual) == 1L && is.finite(actual)) used.cores <- actual
  }
  rows <- lapply(seq_along(results), function(k) {
    r <- results[[k]]; info <- task_info[[k]]
    comp <- r$components
    if (is.null(comp) || !nrow(comp)) comp <- data.frame(lower = NA_real_, upper = NA_real_)
    plugin <- r$plugin
    if (is.null(plugin) || length(plugin) != 2L) plugin <- c(NA_real_, NA_real_)
    data.frame(target = info$target, label = info$label, cohort = info$cohort,
               relative_period = info$relative_period, M = info$M,
               component = seq_len(nrow(comp)), lower = comp$lower, upper = comp$upper,
               center = if (length(r$center) == 1L) r$center else NA_real_,
               plugin.lower = unname(plugin[1L]), plugin.upper = unname(plugin[2L]),
               success = isTRUE(r$success), status = if (length(r$status)) r$status else "unknown",
               message = if (length(r$message)) paste(r$message, collapse = "; ") else "",
               stringsAsFactors = FALSE)
  })
  intervals <- do.call(rbind, rows)
  failures <- sum(!vapply(results, function(x) isTRUE(x$success), logical(1)))
  if (failures) {
    unsupported <- sum(vapply(results, function(x) !isTRUE(x$success) &&
      isTRUE(x$status %in% c("unsupported_zero_variance", "unsupported_moment_variance")),
      logical(1)))
    detail <- if (unsupported) paste0(" Conditional/C-LF inference is unsupported for ",
      unsupported, " task(s) with deterministic or numerically degenerate target/moment variances; ",
      "no confidence set is returned for those tasks.") else ""
    warning(failures, " inference task(s) failed; inspect intervals$status and intervals$message.",
            detail, call. = FALSE)
  }
  out <- structure(list(call = match.call(), intervals = intervals, estimates = do.call(rbind, estimates),
                 fit = fit, targets = targets, results = results,
                 settings = list(restriction = restriction, method = method, full.path = full.path,
                                 target = target, framework = framework, M = M, conf.level = conf.level,
                                 parallel = parallel, cores = used.cores,
                                 seed = seed, task.seeds = seeds, control = control),
                 diagnostics = fit$diagnostics), class = "cohort_sensitivity")
  if (balanced) {
    out$settings$balanced <- TRUE
    out$settings$cohort.paths <- if (framework == "cohort") cohort.paths else NA_character_
    out$settings$window <- window
    if (framework == "aggregate") {
      # A, the aggregated coefficients, and their normalization do not depend
      # on the target; only the aggregate target weights do.
      target.weights <- do.call(rbind, lapply(aggregate.targets, as.numeric))
      dimnames(target.weights) <- list(vapply(targets, function(z) z$id, character(1)),
                                       names(rep$betahat))
      out$aggregation <- list(aggregation_matrix = rep$aggregation_matrix, betahat = rep$betahat,
                              Sigma = rep$Sigma, normalization = rep$normalization,
                              coordinate_info = rep$coordinate_info, targets = target.weights)
    }
  }
  out
}
