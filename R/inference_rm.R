# Every RM confidence set is a union over signed pre-period benchmarks.
# Conditional inversion must not borrow the C-LF first-stage support gate.
.hc_rm_component_job <- function(job, representation, scenarios, M, alpha,
                                  seed, sims, method) {
  tryCatch({
  engine <- .hc_backend()
  restriction <- engine$v3_rm_component_restriction(
    representation$coordinate_info, M, scenarios$table[job, , drop = FALSE],
    scenarios$operator)
  built <- engine$v3_reduce_polyhedral_representation(representation, restriction)
  geometry <- engine$v3_make_arp_geometry(built, scope = "post_connected")
  if (method == "c_lf") return(engine$sd_lf_critical_value(
    geometry = geometry, Sigma = built$Sigma_bar, hybrid_kappa = alpha / 10,
    sims = as.integer(sims), seed = as.integer(seed + job),
    max_failure_fraction = 0, engine = "diagnostic"))
  problem <- engine$sd_bind_arp_problem(
    geometry, built$beta_bar, built$Sigma_bar,
    procedure = "conditional", alpha = alpha)
  if (!isTRUE(problem$success)) stop(problem$message, call. = FALSE)
  anchor <- engine$sd_arp_profile_anchor(problem)
  if (!isTRUE(anchor$success)) stop(anchor$message, call. = FALSE)
  list(component = job, problem = problem, profile_anchor = anchor)
  }, sd_unsupported_moment_variance = function(e) list(success = FALSE,
    status = "unsupported_moment_variance", message = conditionMessage(e)))
}

.hc_rm_cluster <- function(cores) {
  cluster <- parallel::makePSOCKcluster(as.integer(cores))
  ok <- FALSE
  on.exit(if (!ok) parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, function(paths) {
    .libPaths(paths)
    loadNamespace("CohortHonestDiD")
    NULL
  }, .libPaths())
  ok <- TRUE
  cluster
}

.hc_rm_worker_job <- function(job, ...) {
  utils::getFromNamespace(".hc_rm_component_job", "CohortHonestDiD")(job, ...)
}

.hc_rm_test_batch <- function(components, theta) {
  prepared <- list(success = TRUE, components = components,
                   component_order = seq_along(components))
  .hc_backend()$v3_rm_union_point_test(prepared, theta, use_lf_support = FALSE)
}

.hc_rm_worker_test <- function(components, theta) {
  utils::getFromNamespace(".hc_rm_test_batch", "CohortHonestDiD")(components, theta)
}

.hc_rm_union_result <- function(results, theta) {
  valid <- vapply(results, function(x) isTRUE(x$success) &&
                    is.logical(x$accept) && length(x$accept) == 1L && !is.na(x$accept),
                  logical(1))
  accepted <- vapply(results, function(x) isTRUE(x$accept), logical(1)) & valid
  warnings <- unique(unlist(lapply(results, `[[`, "warnings"), use.names = FALSE))
  if (any(accepted)) {
    failures <- unique(unlist(lapply(results[!valid], `[[`, "message"), use.names = FALSE))
    if (length(failures)) warnings <- unique(c(warnings, paste0("Earlier RM component failure: ", failures)))
    return(list(success = TRUE, accept = TRUE,
                status = if (all(valid)) "ok" else "accepted_after_component_failure",
                message = NULL, theta = theta, warnings = warnings))
  }
  if (!all(valid)) return(list(success = FALSE, accept = NA,
    status = "rm_union_component_failure", theta = theta, warnings = warnings,
    message = paste(unique(unlist(lapply(results[!valid], `[[`, "message"))), collapse = " | ")))
  list(success = TRUE, accept = FALSE, status = "ok", message = NULL,
       theta = theta, warnings = warnings)
}

.hc_infer_rm <- function(representation, M, method, alpha, seed, control) {
  engine <- .hc_backend()
  info <- representation$coordinate_info
  window <- representation[["restriction_window"]]
  pre.rows <- info$type == "pre"
  if (!is.null(window)) pre.rows <- pre.rows & info$relative_period >= window
  pre.count <- table(info$cohort[pre.rows])
  count <- if (M == 0) 1 else prod(2 * (as.numeric(pre.count) - 1))
  if (!is.finite(count) || count > control$max.components)
    stop("RM requires ", format(count), " signed benchmark components; increase max.components deliberately.", call. = FALSE)
  scenarios <- engine$v3_rm_scenario_table(info, M, pre_window = window)
  plugin.restriction <- engine$v3_rm_plugin_restriction(info, representation$betahat, M, pre_window = window)
  built <- engine$v3_reduce_polyhedral_representation(representation, plugin.restriction)
  plugin <- engine$v3_plugin_identified_set(built, scope = "post_connected")
  search <- if (method == "c_lf") .hc_engine_grid(built, plugin, numeric(), control) else NULL
  n <- nrow(scenarios$table)
  cores <- min(control$component.cores, n)
  cluster <- if (cores > 1L) .hc_rm_cluster(cores) else NULL
  if (!is.null(cluster)) on.exit(parallel::stopCluster(cluster), add = TRUE)
  args <- list(representation = representation, scenarios = scenarios, M = M,
               alpha = alpha, seed = seed, sims = control$lf.sims, method = method)
  jobs <- seq_len(n)
  prepared <- if (is.null(cluster)) {
    do.call(lapply, c(list(X = jobs, FUN = .hc_rm_component_job), args))
  } else {
    do.call(parallel::parLapply, c(list(cl = cluster, X = jobs,
                                      fun = .hc_rm_worker_job), args))
  }
  unsupported <- vapply(prepared, function(x)
    identical(x$status, "unsupported_moment_variance"), logical(1))
  if (any(unsupported)) return(.hc_engine_failure("unsupported_moment_variance",
    prepared[[which(unsupported)[1L]]]$message))
  diagnostics <- list(signed.components = n, component.cores = cores,
                       endpoint.tolerance = search$tolerance,
                       topology.certified = FALSE, inversion.exact = FALSE)
  budget <- .hc_search_budget(control$max.evaluations)
  if (method == "c_lf") {
    context <- new.env(parent = emptyenv())
    run_attempt <- function(points) {
    attempt.diagnostics <- diagnostics
    # Seeds and scenario order agree across serial and socket execution.
    fit <- engine$v3_clf_rm_union(
      representation, M, theta_true = sum(built$target_bar * built$beta_bar),
      alpha = alpha, lf_sims = control$lf.sims, lf_seed = seed,
      lf_cv_by_component = prepared,
      topology_initial_radius_se = control$grid.radius,
      topology_mesh_se = 2 * control$grid.radius / (points - 1),
      endpoint_tolerance_se = search$tolerance / search$se,
      max_expansions = control$max.expansions,
      max_test_evaluations = control$max.evaluations,
      topology_audit = control$topology.audit, project_root = "",
      search_context = context, test_budget = budget)
    attempt.diagnostics$lf.critical.values <- fit$lf_cv_by_component
    attempt.diagnostics$resolution.stable <- fit$resolution_stable
    attempt.diagnostics$support.certified <- fit$support_certified
    attempt.diagnostics$component.diagnostics <- fit$component_diagnostics
    if (identical(fit$resolution_stable, FALSE)) {
      fit$success <- FALSE
      fit$status <- "topology_audit_failure"
      fit$message <- "RM inversion disagreed across diagnostic grids."
      if (!isTRUE(fit$resolution_audit$success)) fit$status <- "resolution_audit_failure"
    }
    list(fit = fit, diagnostics = attempt.diagnostics)
    }
    result <- .hc_adaptive_search(control, budget, run_attempt)
    return(.hc_engine_result(result$fit, plugin, result$diagnostics))
  }
  anchors <- vapply(prepared, function(x) x$profile_anchor$theta, numeric(1))
  # Cache across endpoint, topology, and expansion passes. Numerical failures
  # remain failures unless another signed component proves union membership.
  batches <- split(prepared, rep(seq_len(cores), length.out = n))
  test <- function(theta) budget$evaluate(theta, function() {
    if (is.null(cluster)) .hc_rm_test_batch(prepared, theta) else {
      .hc_rm_union_result(parallel::parLapply(cluster, batches,
                                             .hc_rm_worker_test, theta = theta), theta)
    }
  })
  .hc_conditional_search(representation, built, plugin, "RM", control,
    anchors, prepared[[1L]]$problem, test, budget, diagnostics)
}
