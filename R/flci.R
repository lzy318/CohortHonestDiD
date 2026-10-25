.hc_infer_flci <- function(representation, built, plugin, scope, alpha, seed, control) {
  engine <- .hc_backend()
  args <- control$flci.control
  if (control$flci.solver == "direct") args$seed <- seed
  rule <- engine$sd_flci_prepare(built$model, scope = scope, alpha = alpha,
                                solver = control$flci.solver, solve_args = args)
  if (!isTRUE(rule$success)) return(.hc_engine_failure(
    rule$status, rule$message, list(warnings = rule$warnings),
    c(plugin$lower, plugin$upper)))
  fit <- engine$sd_flci_interval(rule, representation$betahat)
  answer <- .hc_engine_result(fit, plugin, list(
    solver.status = rule$solver_status, standard.error = rule$standard_error,
    maximum.bias = rule$maximum_bias, half.length = rule$half_length,
    folded.critical.value = rule$folded_critical_value,
    tie.break = rule$raw_fit$tie_break, covariance.rank = built$model$diagnostics$Sigma_rank,
    endpoint.certified = TRUE, topology.certified = TRUE))
  full <- as.numeric(built$model$E %*% rule$v)
  names(full) <- representation$.original_names[rownames(built$model$E)]
  answer$weights <- list(gamma = rule$gamma, coefficients = full)
  answer
}
