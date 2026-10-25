.hc_infer_sd <- function(representation, M, method, full.path, alpha, seed, control) {
  engine <- .hc_backend()
  scope <- if (full.path) "full_path" else "post_connected"
  restriction <- engine$v3_sd_restriction(representation$coordinate_info, M, pre_window = representation[["restriction_window"]])
  built <- engine$v3_reduce_representation(representation, restriction,
                                           allow_singular_covariance = TRUE)
  plugin <- engine$v3_plugin_identified_set(built, scope = scope, allow_empty = TRUE)
  if (method == "flci") return(.hc_infer_flci(
    representation, built, plugin, scope, alpha, seed, control))
  geometry <- engine$v3_make_arp_geometry(built, scope = scope)
  problem <- engine$sd_bind_arp_problem(
    geometry, built$beta_bar, built$Sigma_bar, procedure = method,
    alpha = alpha, hybrid_kappa = alpha / 10, lf_sims = control$lf.sims,
    lf_seed = seed, lf_max_failure_fraction = 0, lf_engine = "diagnostic")
  if (!isTRUE(problem$success)) return(.hc_engine_failure(
    problem$status, problem$message, list(warnings = problem$warnings),
    c(plugin$lower, plugin$upper)))
  anchor <- engine$sd_arp_profile_anchor(problem)
  if (!isTRUE(anchor$success)) return(.hc_engine_failure(
    anchor$status, anchor$message, list(warnings = anchor$warnings),
    c(plugin$lower, plugin$upper)))
  if (method == "conditional") {
    budget <- .hc_search_budget(control$max.evaluations)
    test <- function(theta) budget$evaluate(theta, function()
      engine$sd_arp_truth_test(problem, theta))
    return(.hc_conditional_search(representation, built, plugin, "SD", control,
      anchor$theta, problem, test, budget))
  }
  support <- NULL
  if (method == "c_lf") {
    support <- engine$v3_rm_lf_component_support(problem)
    if (!isTRUE(support$success)) {
      # The same globally profiled dual lower bound used by sd_arp_invert
      # can certify LF rejection everywhere, independently of a finite grid.
      empty.tolerance <- max(1e-10, anchor$dual_tolerance,
                             anchor$dual_objective_residual)
      empty <- identical(support$status, "lf_support_no_feasible") &&
        is.finite(anchor$eta_lower_bound) &&
        anchor$eta_lower_bound > problem$hybrid_list$lf_cv + empty.tolerance
      diagnostics <- list(lf.critical.value = problem$hybrid_list$lf_cv,
        support.certified = empty, support = support, empty.certified = empty,
        topology.certified = empty, inversion.exact = empty,
        test.evaluations = 0L, refinements = 0L, final.grid.points = control$grid.points,
        search.attempts = data.frame(attempt = 0L, grid.points = control$grid.points,
          new.evaluations = 0L, cumulative.evaluations = 0L, success = empty,
          status = if (empty) "empty_certified" else support$status))
      if (empty) return(.hc_engine_result(list(success = TRUE,
        status = "empty_certified", components = data.frame(lower = numeric(), upper = numeric())),
        plugin, diagnostics))
      return(.hc_engine_failure(support$status, support$message, diagnostics,
                                c(plugin$lower, plugin$upper)))
    }
  }
  budget <- .hc_search_budget(control$max.evaluations)
  test <- function(theta) budget$evaluate(theta, function()
    engine$sd_arp_truth_test(problem, theta))
  run_attempt <- function(points) {
  attempt.control <- control
  attempt.control$grid.points <- points
  search <- .hc_engine_grid(built, plugin, anchor$theta, attempt.control,
    bounds = if (is.null(support)) NULL else c(support$lower, support$upper))
  expansions <- if (is.null(support)) control$max.expansions else 0L
  fit <- engine$sd_arp_invert(problem, search$grid,
                            endpoint_tolerance = search$tolerance,
                            max_expansions = expansions,
                            max_test_evaluations = control$max.evaluations,
                            expansion_step = search$spacing,
                            profile_anchor = anchor, test_fn = test)
  diagnostics <- list(lf.critical.value = problem$hybrid_list$lf_cv,
                      endpoint.tolerance = search$tolerance,
                      test.evaluations = fit$total_test_evaluations,
                      endpoint.certified = FALSE, topology.certified = FALSE,
                      resolution.stable = NA, inversion.exact = FALSE)
  if (!is.null(support)) {
    diagnostics$search.domain <- range(search$grid)
    diagnostics$support.certified <- isTRUE(support$certified)
    diagnostics$support <- support
  }
  if (isTRUE(fit$success) && nrow(fit$components)) {
    endpoint <- engine$sd_arp_audit_endpoints(
      problem, fit, probe_distance = max(10 * search$tolerance, 1e-4 * search$se),
      test_fn = test)
    diagnostics$endpoint.certified <- isTRUE(endpoint$certified)
    if (!isTRUE(endpoint$success)) {
      fit$success <- FALSE
      fit$status <- endpoint$status
      fit$message <- endpoint$message
      if (nrow(endpoint$checks) && any(!endpoint$checks$success &
          endpoint$checks$status != "endpoint_probe_too_close")) fit$status <- "endpoint_test_failure"
    }
    if (isTRUE(fit$success) && isTRUE(control$topology.audit)) {
      audit <- engine$sd_arp_topology_audit(
        problem, fit, search$grid, endpoint_tolerance = search$tolerance,
        max_expansions = expansions,
        max_test_evaluations = control$max.evaluations, profile_anchor = anchor,
        expansion_step = search$spacing, test_fn = test)
      diagnostics$resolution.stable <- isTRUE(audit$success)
      if (!isTRUE(audit$success)) {
        fit$success <- FALSE
        fit$status <- audit$status
        fit$message <- audit$message
        failed <- Filter(function(x) !isTRUE(x$success), list(audit$shifted, audit$fine))
        if (length(failed)) fit$status <- failed[[1L]]$status
      }
    }
  }
  list(fit = fit, diagnostics = diagnostics)
  }
  result <- .hc_adaptive_search(control, budget, run_attempt)
  .hc_engine_result(result$fit, plugin, result$diagnostics)
}
