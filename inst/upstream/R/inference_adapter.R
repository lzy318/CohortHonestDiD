# Common inference adapter for the cohort and aggregate representations.
#
# SD cells use full-path FLCI or post-connected C-LF.  RM cells evaluate the
# pointwise logical union of post-connected C-LF tests over signed benchmark
# scenarios on a certified finite LF-support domain.

v3_validate_representation <- function(representation) {
  required <- c(
    "name", "coordinate_info", "W", "normalization", "target",
    "betahat", "Sigma"
  )
  missing <- setdiff(required, names(representation))
  v3_assert(!length(missing), paste0(
    "Representation is missing: ", paste(missing, collapse = ", ")
  ))
  info <- representation$coordinate_info
  v3_validate_coordinate_info(info)
  names_all <- as.character(info$coordinate)
  W <- as.matrix(representation$W)
  Sigma <- as.matrix(representation$Sigma)
  normalization <- as.matrix(representation$normalization)
  v3_assert(identical(dim(W), c(length(names_all), length(names_all))),
            "W has the wrong dimensions.")
  v3_assert(setequal(rownames(W), names_all) && setequal(colnames(W), names_all),
            "W names do not match the representation coordinates.")
  v3_assert(identical(dim(Sigma), c(length(names_all), length(names_all))),
            "Sigma has the wrong dimensions.")
  v3_assert(setequal(rownames(Sigma), names_all) &&
              setequal(colnames(Sigma), names_all),
            "Sigma names do not match the representation coordinates.")
  v3_assert(ncol(normalization) == length(names_all) &&
              setequal(colnames(normalization), names_all),
            "normalization is not aligned to the representation coordinates.")
  target <- representation$target[names_all]
  beta <- representation$betahat[names_all]
  v3_assert(!anyNA(target) && !anyNA(beta) &&
              all(is.finite(c(target, beta, W, Sigma, normalization))),
            "Representation contains missing or nonfinite values.")
  v3_assert(v3_max_abs(target[info$type == "pre"]) <= 1e-12 &&
              v3_max_abs(target[info$type == "post"]) > 0,
            "The target must load only on post-period coordinates.")
  eigenvalues <- eigen(v3_symmetrize(Sigma), symmetric = TRUE,
                       only.values = TRUE)$values
  scale <- max(1, max(abs(eigenvalues)))
  v3_assert(min(eigenvalues) >= -1e-8 * scale,
            "Sigma is not positive semidefinite.")
  invisible(TRUE)
}

v3_align_representation <- function(representation) {
  v3_validate_representation(representation)
  names_all <- as.character(representation$coordinate_info$coordinate)
  representation$W <- representation$W[names_all, names_all, drop = FALSE]
  representation$Sigma <- v3_symmetrize(
    representation$Sigma[names_all, names_all, drop = FALSE]
  )
  representation$normalization <- representation$normalization[
    , names_all, drop = FALSE
  ]
  representation$target <- representation$target[names_all]
  representation$betahat <- representation$betahat[names_all]
  representation
}

v3_reduce_representation <- function(
    representation, restriction, allow_singular_covariance = FALSE) {
  representation <- v3_align_representation(representation)
  names_all <- as.character(representation$coordinate_info$coordinate)
  A_Delta <- restriction$A_Delta[, names_all, drop = FALSE]
  d_vec <- restriction$d_vec[rownames(A_Delta)]
  model <- flci_reduce_model(
    W = representation$W,
    Sigma = representation$Sigma,
    A_Delta = A_Delta,
    d_vec = d_vec,
    target = representation$target,
    normalization = representation$normalization,
    estimator_type = representation$estimator_type %||% "imputation",
    allow_singular_covariance = allow_singular_covariance
  )
  reduced_names <- colnames(model$E)
  dimnames(model$Sigma_bar) <- list(reduced_names, reduced_names)
  dimnames(model$W_bar) <- list(reduced_names, reduced_names)
  beta_bar <- setNames(
    as.numeric(t(model$E) %*% representation$betahat), reduced_names
  )
  target_bar <- setNames(
    as.numeric(t(model$E) %*% representation$target), reduced_names
  )
  beta_support_residual <- v3_max_abs(
    representation$betahat - model$E %*% beta_bar
  )
  covariance_support_residual <- v3_max_abs(
    representation$Sigma - model$E %*% model$Sigma_bar %*% t(model$E)
  )
  v3_assert(beta_support_residual <= 1e-8 *
              max(1, v3_max_abs(representation$betahat)),
            "Coefficient estimate is outside the mechanical estimator support.")
  v3_assert(covariance_support_residual <= 1e-7 *
              max(1, v3_max_abs(representation$Sigma)),
            "Covariance is outside the mechanical estimator support.")
  list(
    representation = representation,
    restriction = restriction,
    model = model,
    reduced_names = reduced_names,
    beta_bar = beta_bar,
    target_bar = target_bar,
    Sigma_bar = model$Sigma_bar,
    beta_support_residual = beta_support_residual,
    covariance_support_residual = covariance_support_residual
  )
}

# Family-generic quotient reducer.  Unlike flci_reduce_model(), this function
# does not inspect SD row names or require a second-difference block.  It is
# therefore the correct reducer for the signed RM polyhedra.
v3_reduce_polyhedral_representation <- function(representation, restriction,
                                                tol = 1e-10) {
  representation <- v3_align_representation(representation)
  info <- representation$coordinate_info
  names_all <- as.character(info$coordinate)
  pre <- info$type == "pre"
  post <- info$type == "post"
  normalization <- representation$normalization[, names_all, drop = FALSE]
  v3_assert(v3_max_abs(normalization[, post, drop = FALSE]) <= 10 * tol,
            "Mechanical normalization must load only on pre coordinates.")
  E_pre <- flci_null_basis(normalization[, pre, drop = FALSE], tol = tol)
  E <- matrix(0, length(names_all), ncol(E_pre) + sum(post))
  if (ncol(E_pre)) E[pre, seq_len(ncol(E_pre))] <- E_pre
  E[post, ncol(E_pre) + seq_len(sum(post))] <- diag(sum(post))
  reduced_names <- c(
    paste0("pre_basis_", seq_len(ncol(E_pre))), names_all[post]
  )
  dimnames(E) <- list(names_all, reduced_names)
  W_bar <- t(E) %*% representation$W %*% E
  Sigma_bar <- v3_symmetrize(t(E) %*% representation$Sigma %*% E)
  dimnames(W_bar) <- dimnames(Sigma_bar) <- list(reduced_names, reduced_names)
  beta_bar <- setNames(
    as.numeric(t(E) %*% representation$betahat), reduced_names
  )
  target_bar <- setNames(
    as.numeric(t(E) %*% representation$target), reduced_names
  )
  v3_assert(qr(W_bar, tol = tol)$rank == ncol(E),
            "The normalized coefficient-bias map is rank deficient.")
  # The RM plug-in benchmark is calculated from observed pre-path
  # coefficients.  This is valid only when those normalized pre coordinates
  # estimate the corresponding latent pre path without mixing in post bias.
  # Both supported representations satisfy this property (cohort imputation
  # and the induced aggregate event-study map); fail closed for any future
  # representation that does not.
  reduced_pre <- seq_len(ncol(E_pre))
  reduced_post <- ncol(E_pre) + seq_len(sum(post))
  pre_identity_residual <- v3_max_abs(
    W_bar[reduced_pre, reduced_pre, drop = FALSE] - diag(ncol(E_pre))
  )
  pre_post_residual <- v3_max_abs(
    W_bar[reduced_pre, reduced_post, drop = FALSE]
  )
  v3_assert(pre_identity_residual <= 10 * tol &&
              pre_post_residual <= 10 * tol,
            paste0(
              "RM plug-in benchmarking requires an identity pre map and no ",
              "post-to-pre mixing (residuals ",
              signif(pre_identity_residual, 6), " and ",
              signif(pre_post_residual, 6), ")."
            ))
  support_residual <- v3_max_abs(
    representation$Sigma - E %*% Sigma_bar %*% t(E)
  )
  support_scale <- max(1, v3_max_abs(representation$Sigma))
  v3_assert(support_residual <= 1e-7 * support_scale,
            paste0("Covariance is outside the estimator support (residual ",
                   signif(support_residual, 8), ")."))
  beta_support_residual <- v3_max_abs(
    representation$betahat - E %*% beta_bar
  )
  v3_assert(beta_support_residual <= 1e-8 *
              max(1, v3_max_abs(representation$betahat)),
            "Coefficient estimate is outside the mechanical estimator support.")
  eigenvalues <- eigen(Sigma_bar, symmetric = TRUE, only.values = TRUE)$values
  v3_assert(min(eigenvalues) >= -tol * max(abs(eigenvalues)),
            "Reduced covariance is not positive semidefinite.")
  model <- list(
    E = E, W_bar = W_bar, Sigma_bar = Sigma_bar,
    estimator_type = "imputation"
  )
  A_Delta <- restriction$A_Delta[, names_all, drop = FALSE]
  d_vec <- restriction$d_vec[rownames(A_Delta)]
  list(
    representation = representation,
    restriction = list(
      A_Delta = A_Delta, d_vec = d_vec,
      operator = restriction$operator %||% NULL
    ),
    model = model, reduced_names = reduced_names,
    beta_bar = beta_bar, target_bar = target_bar,
    Sigma_bar = Sigma_bar, support_residual = support_residual,
    beta_support_residual = beta_support_residual,
    pre_identity_residual = pre_identity_residual,
    pre_post_residual = pre_post_residual
  )
}

v3_make_arp_geometry <- function(built, scope = "post_connected", tol = 1e-10) {
  model <- built$model
  A_delta_E <- built$restriction$A_Delta %*% model$E
  A_final <- t(solve(t(model$W_bar), t(A_delta_E)))
  dimnames(A_final) <- list(
    rownames(built$restriction$A_Delta), built$reduced_names
  )
  residual <- v3_max_abs(A_final %*% model$W_bar - A_delta_E)
  v3_assert(residual <= tol * max(1, v3_max_abs(A_delta_E)),
            paste0("ARP restriction-map residual is ", signif(residual, 8), "."))
  geometry <- sd_arp_geometry(
    A = A_final,
    d = built$restriction$d_vec,
    target = built$target_bar,
    scope = scope,
    tolerance = tol
  )
  geometry$restriction_map_residual <- residual
  geometry
}

v3_glpk_solve_free <- function(objective, matrix, direction, rhs,
                               maximize = FALSE,
                               allow_infeasible = FALSE) {
  v3_assert(requireNamespace("Rglpk", quietly = TRUE),
            "Package 'Rglpk' is required for identified-set calculations.")
  q <- length(objective)
  bounds <- list(
    lower = list(ind = seq_len(q), val = rep(-Inf, q)),
    upper = list(ind = seq_len(q), val = rep(Inf, q))
  )
  fit <- Rglpk::Rglpk_solve_LP(
    obj = objective, mat = matrix, dir = direction, rhs = rhs,
    bounds = bounds, types = rep("C", q), max = maximize,
    control = list(canonicalize_status = FALSE)
  )
  # Raw GLPK LP status: 4 is GLP_NOFEAS and 5 is GLP_OPT.  Keeping the raw
  # code lets population-set diagnostics distinguish a scientifically valid
  # empty identified set from an undefined, merely feasible, or unbounded
  # optimization result.  Inference calls retain the fail-closed default.
  if (isTRUE(allow_infeasible) && length(fit$status) == 1L &&
      fit$status == 4L) {
    return(structure(fit, class = c("v3_lp_infeasible", class(fit))))
  }
  v3_assert(length(fit$status) == 1L && fit$status == 5L &&
              length(fit$optimum) == 1L && is.finite(fit$optimum) &&
              length(fit$solution) == q && all(is.finite(fit$solution)),
            paste0("Identified-set LP failed with status ",
                   paste(fit$status, collapse = ","), "."))
  fit
}

v3_plugin_identified_set <- function(built, scope = "post_connected",
                                     tol = 1e-9,
                                     allow_empty = FALSE) {
  scope <- match.arg(scope, c("full_path", "post_connected"))
  reduced_names <- built$reduced_names
  pre <- grepl("^pre_basis_", reduced_names)
  post <- grepl("^delta_post_", reduced_names)
  v3_assert(all(pre | post) && any(pre) && any(post),
            "Reduced coordinates are not a pre/post partition.")
  A_u <- built$restriction$A_Delta %*% built$model$E
  colnames(A_u) <- reduced_names
  restriction_rows <- if (scope == "full_path") {
    seq_len(nrow(A_u))
  } else {
    which(rowSums(abs(A_u[, post, drop = FALSE]) > tol) > 0L)
  }
  constraint_matrix <- rbind(
    A_u[restriction_rows, , drop = FALSE],
    built$model$W_bar[pre, , drop = FALSE]
  )
  direction <- c(rep("<=", length(restriction_rows)), rep("==", sum(pre)))
  rhs <- c(
    built$restriction$d_vec[restriction_rows], built$beta_bar[pre]
  )
  bias_objective <- as.numeric(t(built$model$W_bar) %*% built$target_bar)
  raw_target <- sum(built$target_bar * built$beta_bar)
  bias_min <- v3_glpk_solve_free(
    bias_objective, constraint_matrix, direction, rhs, maximize = FALSE,
    allow_infeasible = allow_empty
  )
  if (inherits(bias_min, "v3_lp_infeasible")) {
    return(list(
      lower = NA_real_, upper = NA_real_, width = 0,
      center = NA_real_, raw_target = as.numeric(raw_target),
      bias_min = NA_real_, bias_max = NA_real_, scope = scope,
      n_restriction_rows = length(restriction_rows), empty = TRUE,
      status = "empty_infeasible", glpk_status = 4L
    ))
  }
  bias_max <- v3_glpk_solve_free(
    bias_objective, constraint_matrix, direction, rhs, maximize = TRUE,
    allow_infeasible = allow_empty
  )
  v3_assert(!inherits(bias_max, "v3_lp_infeasible"),
            "The identified-set maximum became infeasible after a feasible minimum.")
  lower <- raw_target - bias_max$optimum
  upper <- raw_target - bias_min$optimum
  list(
    lower = as.numeric(lower), upper = as.numeric(upper),
    width = as.numeric(upper - lower), center = as.numeric((upper + lower) / 2),
    raw_target = as.numeric(raw_target),
    bias_min = as.numeric(bias_min$optimum),
    bias_max = as.numeric(bias_max$optimum),
    scope = scope, n_restriction_rows = length(restriction_rows),
    empty = FALSE, status = "ok", glpk_status = 5L
  )
}

v3_topology_grid <- function(built, identified_set, profile_anchor,
                             initial_radius_se = 8, mesh_se = 0.25) {
  target_se <- sqrt(as.numeric(
    t(built$target_bar) %*% built$Sigma_bar %*% built$target_bar
  ))
  v3_assert(is.finite(target_se) && target_se > 0,
            "Target standard error must be positive.")
  center <- profile_anchor$theta
  lower_needed <- min(identified_set$lower, center) - initial_radius_se * target_se
  upper_needed <- max(identified_set$upper, center) + initial_radius_se * target_se
  step <- mesh_se * target_se
  n_left <- max(2L, ceiling((center - lower_needed) / step))
  n_right <- max(2L, ceiling((upper_needed - center) / step))
  center + seq.int(-n_left, n_right) * step
}

v3_component_hull <- function(components) {
  if (is.null(components) || !nrow(components)) {
    return(list(lower = NA_real_, upper = NA_real_, width = 0,
                center = NA_real_, total_length = 0, n_components = 0L))
  }
  lower <- min(components$lower)
  upper <- max(components$upper)
  list(
    lower = lower, upper = upper, width = upper - lower,
    center = (lower + upper) / 2,
    total_length = sum(components$upper - components$lower),
    n_components = nrow(components)
  )
}

v3_merge_components <- function(components, tolerance = 1e-8) {
  if (is.null(components) || !nrow(components)) {
    return(data.frame(lower = numeric(), upper = numeric()))
  }
  components <- components[is.finite(components$lower) &
                             is.finite(components$upper), c("lower", "upper")]
  if (!nrow(components)) return(data.frame(lower = numeric(), upper = numeric()))
  components <- components[order(components$lower, components$upper), , drop = FALSE]
  merged <- list()
  current <- as.numeric(components[1L, ])
  cursor <- 0L
  if (nrow(components) > 1L) {
    for (j in 2:nrow(components)) {
      one <- as.numeric(components[j, ])
      scale <- max(1, abs(current), abs(one))
      if (one[1L] <= current[2L] + tolerance * scale) {
        current[2L] <- max(current[2L], one[2L])
      } else {
        cursor <- cursor + 1L
        merged[[cursor]] <- current
        current <- one
      }
    }
  }
  cursor <- cursor + 1L
  merged[[cursor]] <- current
  ans <- as.data.frame(do.call(rbind, merged))
  names(ans) <- c("lower", "upper")
  ans
}

v3_flci_sd <- function(representation, M_raw, theta_true = 0,
                       alpha = 0.05, solver = "direct", solve_args = list(),
                       project_root, prepared_rule = NULL,
                       allow_singular_covariance = FALSE) {
  dependency_root <- file.path(project_root, "dependencies")
  sd_inference_load_dependencies(dependency_root)
  restriction <- v3_sd_restriction(representation$coordinate_info, M_raw)
  built <- v3_reduce_representation(
    representation, restriction,
    allow_singular_covariance = allow_singular_covariance
  )
  rule <- if (is.null(prepared_rule)) {
    sd_flci_prepare(
      model = built$model, scope = "full_path", alpha = alpha,
      solver = solver, solve_args = solve_args,
      project_root = dependency_root
    )
  } else {
    v3_assert(is.list(prepared_rule) && !is.null(prepared_rule$success) &&
                isTRUE(prepared_rule$success),
              "prepared_rule must be a successful FLCI rule.")
    v3_assert(
      identical(prepared_rule$scope, "full_path") &&
        isTRUE(all.equal(prepared_rule$alpha, alpha, tolerance = 0)) &&
        isTRUE(all.equal(prepared_rule$model, built$model,
                         tolerance = 1e-12, check.attributes = TRUE)),
      "prepared_rule does not match the requested FLCI geometry."
    )
    prepared_rule
  }
  if (!isTRUE(rule$success)) {
    return(list(success = FALSE, status = rule$status, message = rule$message,
                warnings = rule$warnings, plugin = NULL, rule = rule))
  }
  interval <- sd_flci_interval(rule, representation$betahat, theta = theta_true)
  hull <- list(
    lower = interval$lower, upper = interval$upper,
    width = interval$total_length, total_length = interval$total_length,
    center = interval$center, n_components = 1L
  )
  list(
    success = TRUE, status = interval$status, message = NULL,
    warnings = interval$warnings, procedure = "flci", scope = "full_path",
    plugin = NULL, interval = interval, components = interval$components,
    hull = hull, contains_truth = interval$accept,
    target_se = rule$standard_error, maximum_bias = rule$maximum_bias,
    folded_critical_value = rule$folded_critical_value, rule = rule,
    built = built
  )
}

v3_clf_polyhedron <- function(
    built, theta_true = 0, alpha = 0.05, hybrid_kappa = alpha / 10,
    lf_sims = 1000L, lf_seed = 0L, lf_cv = NULL,
    topology_initial_radius_se = 8, topology_mesh_se = 0.25,
    endpoint_tolerance_se = 1e-5, max_expansions = 4L,
    topology_guide = NULL, compute_plugin = TRUE,
    topology_audit = TRUE,
    project_root) {
  dependency_root <- file.path(project_root, "dependencies")
  sd_inference_load_dependencies(dependency_root)
  geometry <- v3_make_arp_geometry(built, scope = "post_connected")
  if (is.null(lf_cv)) {
    lf_cv <- sd_lf_critical_value(
      geometry = geometry, Sigma = built$Sigma_bar,
      hybrid_kappa = hybrid_kappa, sims = as.integer(lf_sims),
      seed = as.integer(lf_seed), max_failure_fraction = 0,
      engine = "diagnostic",
      project_root = dependency_root
    )
  }
  numeric_lf <- is.numeric(lf_cv) && length(lf_cv) == 1L &&
    is.finite(lf_cv)
  object_lf <- inherits(lf_cv, "sd_lf_critical_value") &&
    isTRUE(lf_cv$success) && is.finite(lf_cv$critical_value)
  if (!numeric_lf && !object_lf) {
    return(list(
      success = FALSE,
      status = if (is.list(lf_cv)) lf_cv$status %||%
        "lf_critical_value_unavailable" else "lf_critical_value_unavailable",
      message = if (is.list(lf_cv)) lf_cv$message %||%
        "The supplied LF critical value is invalid." else
        "The supplied LF critical value is invalid.",
      warnings = if (is.list(lf_cv)) lf_cv$warnings %||% character() else
        character(),
      geometry = geometry, lf = lf_cv,
      resolution_audit = NULL, resolution_stable = NA,
      inversion_exact = FALSE, topology_certified = FALSE
    ))
  }
  problem <- sd_bind_arp_problem(
    geometry = geometry, betahat = built$beta_bar, Sigma = built$Sigma_bar,
    procedure = "c_lf", alpha = alpha, hybrid_kappa = hybrid_kappa,
    lf_cv = lf_cv, project_root = dependency_root
  )
  if (!isTRUE(problem$success)) {
    return(list(success = FALSE, status = problem$status,
                message = problem$message, warnings = problem$warnings,
                geometry = geometry, lf = lf_cv, problem = problem,
                resolution_audit = NULL, resolution_stable = NA,
                inversion_exact = FALSE, topology_certified = FALSE))
  }
  # A signed RM component can be incompatible with the realized pre estimates,
  # even though its confidence set is nonempty after accounting for sampling
  # uncertainty.  Consequently, the component-specific plug-in LP is not a
  # prerequisite for inversion.  SD calls compute their own plug-in set; RM
  # calls supply the feasible realization-level RM hull solely as a guide for
  # constructing the topology grid.
  plugin <- if (isTRUE(compute_plugin)) {
    v3_plugin_identified_set(built, scope = "post_connected")
  } else {
    NULL
  }
  guide <- topology_guide %||% plugin
  v3_assert(!is.null(guide) &&
              all(c("lower", "upper") %in% names(guide)) &&
              all(is.finite(c(guide$lower, guide$upper))) &&
              guide$lower <= guide$upper,
            "A finite ordered topology guide is required for C-LF inversion.")
  anchor <- sd_arp_profile_anchor(problem)
  if (!isTRUE(anchor$success)) {
    return(list(success = FALSE, status = anchor$status,
                message = anchor$message, warnings = anchor$warnings,
                plugin = plugin, geometry = geometry, lf = lf_cv,
                problem = problem, profile_anchor = anchor,
                resolution_audit = NULL, resolution_stable = NA,
                inversion_exact = FALSE, topology_certified = FALSE))
  }
  topology_grid <- v3_topology_grid(
    built, guide, anchor,
    initial_radius_se = topology_initial_radius_se,
    mesh_se = topology_mesh_se
  )
  target_se <- sqrt(as.numeric(
    t(built$target_bar) %*% built$Sigma_bar %*% built$target_bar
  ))
  interval <- sd_arp_invert(
    problem = problem, topology_grid = topology_grid,
    endpoint_tolerance = endpoint_tolerance_se * target_se,
    max_expansions = as.integer(max_expansions), truth = theta_true,
    profile_anchor = anchor
  )
  if (!isTRUE(interval$success) || !isTRUE(interval$complete)) {
    return(list(success = FALSE, status = interval$status,
                message = interval$message, warnings = interval$warnings,
                plugin = plugin, geometry = geometry, lf = lf_cv,
                problem = problem, profile_anchor = anchor,
                interval = interval, topology_grid = topology_grid,
                resolution_audit = NULL, resolution_stable = NA,
                inversion_exact = FALSE, topology_certified = FALSE))
  }
  topology_check <- if (isTRUE(topology_audit)) {
    sd_arp_topology_audit(
      problem = problem, base_interval = interval,
      topology_grid = topology_grid,
      endpoint_tolerance = endpoint_tolerance_se * target_se,
      max_expansions = as.integer(max_expansions),
      profile_anchor = anchor
    )
  } else {
    list(
      success = TRUE, status = "not_requested", message = NULL,
      certified = FALSE, resolution_stable = NA, inversion_exact = FALSE
    )
  }
  topology_check$resolution_stable <- if (isTRUE(topology_audit)) {
    isTRUE(topology_check$success)
  } else NA
  topology_check$inversion_exact <- FALSE
  topology_check$certified <- FALSE
  resolution_stable <- if (isTRUE(topology_audit)) {
    isTRUE(topology_check$resolution_stable)
  } else NA
  resolution_warning <- if (identical(resolution_stable, FALSE)) {
    topology_check$message %||%
      "The SD C-LF resolution diagnostic failed or disagreed across grids."
  } else {
    character()
  }
  components <- interval$components[, c("lower", "upper"), drop = FALSE]
  hull <- v3_component_hull(components)
  list(
    success = TRUE,
    status = if (identical(resolution_stable, FALSE)) {
      "ok_with_resolution_warning"
    } else {
      interval$status
    },
    message = NULL,
    warnings = unique(c(interval$warnings,
                        topology_check$warnings %||% character(),
                        resolution_warning)),
    procedure = "c_lf",
    scope = "post_connected", plugin = plugin, interval = interval,
    topology_guide = guide,
    components = components, hull = hull,
    contains_truth = interval$contains_truth,
    target_se = target_se,
    lf_critical_value = if (numeric_lf) as.numeric(lf_cv) else
      lf_cv$critical_value,
    geometry = geometry, lf = lf_cv, problem = problem,
    profile_anchor = anchor, built = built,
    topology_grid = topology_grid,
    topology_audit = topology_check,
    resolution_audit = topology_check,
    resolution_stable = resolution_stable,
    inversion_exact = FALSE,
    topology_certified = FALSE
  )
}

v3_clf_sd <- function(representation, M_raw, theta_true = 0,
                      allow_singular_covariance = FALSE, ...) {
  restriction <- v3_sd_restriction(representation$coordinate_info, M_raw)
  built <- v3_reduce_representation(
    representation, restriction,
    allow_singular_covariance = allow_singular_covariance
  )
  v3_clf_polyhedron(built = built, theta_true = theta_true, ...)
}

v3_clf_rm_componentwise <- function(
    representation, Mbar, theta_true = 0, alpha = 0.05,
    hybrid_kappa = alpha / 10, lf_sims = 1000L, lf_seed = 0L,
    lf_cv_by_component = NULL,
    topology_initial_radius_se = 8, topology_mesh_se = 0.25,
    endpoint_tolerance_se = 1e-5, max_expansions = 4L,
    topology_audit = TRUE,
    project_root) {
  representation <- v3_align_representation(representation)
  scenarios <- v3_rm_scenario_table(representation$coordinate_info, Mbar)
  plugin_restriction <- v3_rm_plugin_restriction(
    representation$coordinate_info, representation$betahat, Mbar
  )
  plugin_built <- v3_reduce_polyhedral_representation(
    representation, plugin_restriction
  )
  plugin <- v3_plugin_identified_set(plugin_built, scope = "post_connected")
  all_components <- list()
  diagnostics <- list()
  component_truth <- logical(nrow(scenarios$table))
  component_lf_values <- rep(NA_real_, nrow(scenarios$table))
  if (!is.null(lf_cv_by_component)) {
    v3_assert(length(lf_cv_by_component) == nrow(scenarios$table),
              "lf_cv_by_component must have one entry per RM scenario.")
  }
  for (j in seq_len(nrow(scenarios$table))) {
    restriction <- v3_rm_component_restriction(
      representation$coordinate_info, Mbar,
      scenarios$table[j, , drop = FALSE], scenarios$operator
    )
    built <- v3_reduce_polyhedral_representation(representation, restriction)
    component_lf <- if (is.null(lf_cv_by_component)) {
      NULL
    } else if (is.list(lf_cv_by_component)) {
      lf_cv_by_component[[j]]
    } else {
      lf_cv_by_component[j]
    }
    fit <- v3_clf_polyhedron(
      built = built, theta_true = theta_true, alpha = alpha,
      hybrid_kappa = hybrid_kappa, lf_sims = lf_sims,
      lf_seed = as.integer(lf_seed + j),
      lf_cv = component_lf,
      topology_initial_radius_se = topology_initial_radius_se,
      topology_mesh_se = topology_mesh_se,
      endpoint_tolerance_se = endpoint_tolerance_se,
      max_expansions = max_expansions, topology_guide = plugin,
      compute_plugin = FALSE, topology_audit = topology_audit,
      project_root = project_root
    )
    diagnostics[[j]] <- data.frame(
      component = j, success = isTRUE(fit$success),
      status = fit$status %||% "missing",
      lf_critical_value = fit$lf_critical_value %||% NA_real_,
      n_components = if (isTRUE(fit$success)) nrow(fit$components) else NA_integer_,
      contains_truth = if (isTRUE(fit$success)) {
        isTRUE(fit$contains_truth)
      } else {
        NA
      },
      topology_certified = if (isTRUE(fit$success)) {
        isTRUE(fit$topology_certified)
      } else {
        FALSE
      },
      stringsAsFactors = FALSE
    )
    component_lf_values[j] <- fit$lf_critical_value %||%
      if (is.numeric(component_lf) && length(component_lf) == 1L) {
        as.numeric(component_lf)
      } else {
        NA_real_
      }
    if (!isTRUE(fit$success)) {
      return(list(
        success = FALSE, status = "rm_component_failure",
        message = paste0("RM signed-benchmark component ", j,
                         " failed: ", fit$status %||% "unknown", "."),
        warnings = fit$warnings %||% character(), plugin = plugin,
        failed_component = j, component_fit = fit,
        component_diagnostics = do.call(rbind, diagnostics),
        lf_cv_by_component = component_lf_values
      ))
    }
    component_truth[j] <- isTRUE(fit$contains_truth)
    if (nrow(fit$components)) {
      one <- fit$components
      one$signed_benchmark_component <- j
      all_components[[length(all_components) + 1L]] <- one
    }
  }
  raw_components <- if (length(all_components)) {
    do.call(rbind, all_components)
  } else {
    data.frame(lower = numeric(), upper = numeric(),
               signed_benchmark_component = integer())
  }
  components <- v3_merge_components(raw_components)
  hull <- v3_component_hull(components)
  contains_truth <- any(component_truth)
  list(
    success = TRUE, status = if (nrow(components)) "ok" else "empty",
    message = NULL, warnings = character(), procedure = "c_lf",
    inversion_engine = "componentwise",
    scope = "post_connected", plugin = plugin, components = components,
    raw_components = raw_components, hull = hull,
    contains_truth = contains_truth,
    target_se = sqrt(as.numeric(
      t(plugin_built$target_bar) %*% plugin_built$Sigma_bar %*%
        plugin_built$target_bar
    )),
    signed_benchmark_components = nrow(scenarios$table),
    signed_benchmark_components_full = scenarios$scenario_count_full,
    component_diagnostics = do.call(rbind, diagnostics),
    lf_cv_by_component = component_lf_values,
    topology_certified = isTRUE(topology_audit) &&
      all(vapply(diagnostics, function(one) {
        isTRUE(one$topology_certified[[1L]])
      }, logical(1L)))
  )
}

# Select the signed benchmark component implied by the realized largest
# absolute pre-period first difference.  This component is tested first.  The
# ordering affects runtime only: the union point test still visits every other
# component before declaring rejection.
v3_rm_preferred_component <- function(scenarios, beta) {
  table <- scenarios$table
  if (nrow(table) <= 1L) return(1L)
  operator <- scenarios$operator
  beta <- beta[colnames(operator$D1)]
  v3_assert(!anyNA(beta),
            "RM preferred-component beta is not aligned to the path.")
  desired <- list()
  for (g in sort(unique(operator$rows$cohort))) {
    use <- operator$rows$cohort == g & operator$rows$component == "pre"
    rows <- operator$rows$row[use]
    values <- as.numeric(operator$D1[rows, , drop = FALSE] %*% beta)
    best <- which.max(abs(values))
    sign <- if (values[best] >= 0) 1L else -1L
    desired[[paste0("g", g)]] <- paste(rows[best], sign, sep = "@@")
  }
  matches <- rep(TRUE, nrow(table))
  for (name in names(desired)) {
    matches <- matches & as.character(table[[name]]) == desired[[name]]
  }
  found <- which(matches)
  v3_assert(length(found) == 1L,
            "The preferred RM signed-benchmark component is not unique.")
  found
}

v3_rm_component_label <- function(scenario_row) {
  fields <- setdiff(names(scenario_row), "component")
  if (!length(fields)) return("Mbar=0")
  paste(paste0(fields, "=", as.character(scenario_row[1L, fields])),
        collapse = ";")
}

v3_rm_glpk_status_class <- function(status) {
  if (length(status) != 1L || is.na(status)) return("malformed")
  switch(
    as.character(as.integer(status)),
    `1` = "undefined", `2` = "feasible", `3` = "infeasible",
    `4` = "no_feasible", `5` = "optimal", `6` = "unbounded",
    "unknown"
  )
}

# The C-LF first stage accepts theta only if some nuisance vector delta obeys
#
#   base_y - a_one * theta - a_minus * delta <= c_LF * sd
#
# on rows_for_arp.  Optimizing theta over this polyhedron therefore gives an
# outer support for the full conditional test, even though the latter can be
# disconnected inside that support.  GLPK's raw status is retained and every
# nominal optimum is independently checked against primal, dual, and objective
# certificates before it can bound the inversion grid.
v3_rm_lf_support_endpoint <- function(
    problem, maximize = FALSE, certificate_tolerance = 1e-8,
    solver_fn = NULL) {
  v3_assert(inherits(problem, "sd_arp_problem") && isTRUE(problem$success),
            "problem must be a successful sd_arp_problem.")
  v3_assert(identical(problem$procedure, "c_lf"),
            "LF support is defined only for a C-LF problem.")
  v3_assert(length(maximize) == 1L && !is.na(maximize),
            "maximize must be one nonmissing logical value.")
  maximize <- isTRUE(maximize)
  v3_assert(length(certificate_tolerance) == 1L &&
              is.finite(certificate_tolerance) && certificate_tolerance > 0,
            "certificate_tolerance must be positive and finite.")
  if (is.null(solver_fn)) {
    v3_assert(requireNamespace("Rglpk", quietly = TRUE),
              "Package 'Rglpk' is required for RM LF-support calculations.")
    solver_fn <- Rglpk::Rglpk_solve_LP
  }
  v3_assert(is.function(solver_fn), "solver_fn must be a function.")

  geometry <- problem$geometry
  rows <- as.integer(geometry$rows_for_arp)
  lf_cv <- as.numeric(problem$hybrid_list$lf_cv)
  y <- as.numeric(problem$base_y[rows])
  design <- cbind(
    theta = as.numeric(geometry$a_one), geometry$a_minus
  )[rows, , drop = FALSE]
  sigma <- problem$sigma_y[rows, rows, drop = FALSE]
  sd <- sqrt(diag(sigma))
  v3_assert(length(rows) > 0L && length(lf_cv) == 1L && is.finite(lf_cv) &&
              nrow(design) == length(rows) && ncol(design) > 0L &&
              all(is.finite(c(y, design, sd))) && all(sd > 0),
            "The RM LF-support program has malformed inputs.")

  # Standardize rowwise before solving.  This is algebraically identical to
  # the displayed LF inequalities and materially improves LP scaling when
  # moment standard errors differ.
  constraint_matrix <- -design / sd
  rhs <- lf_cv - y / sd
  objective <- c(1, rep(0, ncol(design) - 1L))
  bounds <- list(
    lower = list(ind = seq_len(ncol(design)), rep(-Inf, ncol(design))),
    upper = list(ind = seq_len(ncol(design)), rep(Inf, ncol(design)))
  )
  # Rglpk names this element `val`; spell it explicitly after construction so
  # the free-variable contract is transparent in stored diagnostics.
  names(bounds$lower)[2L] <- names(bounds$upper)[2L] <- "val"
  captured_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      solver_fn(
        obj = objective, mat = constraint_matrix,
        dir = rep("<=", length(rhs)), rhs = rhs, bounds = bounds,
        types = rep("C", ncol(design)), max = maximize,
        control = list(canonicalize_status = FALSE)
      ),
      warning = function(w) {
        captured_warnings <<- unique(c(
          captured_warnings, conditionMessage(w)
        ))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  sense <- if (maximize) "maximum" else "minimum"
  if (inherits(fit, "condition")) {
    return(list(
      success = FALSE, certified = FALSE,
      status = "lf_support_solver_error",
      message = paste0("The LF-support ", sense,
                       " LP raised an error: ", conditionMessage(fit)),
      warnings = captured_warnings, sense = sense, maximize = maximize,
      raw_status = NA_integer_, raw_status_class = "solver_error",
      endpoint = NA_real_
    ))
  }
  raw_status <- if (is.list(fit) && length(fit$status) == 1L) {
    suppressWarnings(as.integer(fit$status))
  } else NA_integer_
  raw_class <- v3_rm_glpk_status_class(raw_status)
  if (!identical(raw_class, "optimal")) {
    failure_status <- if (identical(raw_class, "infeasible")) {
      "lf_support_infeasible_nonoptimal"
    } else if (identical(raw_class, "no_feasible")) {
      "lf_support_no_feasible"
    } else if (identical(raw_class, "unbounded")) {
      "lf_support_unbounded"
    } else if (identical(raw_class, "malformed")) {
      "lf_support_status_malformed"
    } else {
      "lf_support_nonoptimal"
    }
    return(list(
      success = FALSE, certified = FALSE, status = failure_status,
      message = paste0("The LF-support ", sense, " LP returned raw GLPK ",
                       "status ", paste(raw_status, collapse = ","), " (",
                       raw_class, ")."),
      warnings = captured_warnings, sense = sense, maximize = maximize,
      raw_status = raw_status, raw_status_class = raw_class,
      endpoint = NA_real_
    ))
  }

  solution <- suppressWarnings(as.numeric(fit$solution))
  dual <- suppressWarnings(as.numeric(fit$auxiliary$dual))
  auxiliary_primal <- suppressWarnings(as.numeric(fit$auxiliary$primal))
  reduced_cost <- suppressWarnings(as.numeric(fit$solution_dual))
  optimum <- suppressWarnings(as.numeric(fit$optimum))
  malformed <- length(solution) != ncol(design) ||
    length(dual) != length(rhs) ||
    length(auxiliary_primal) != length(rhs) ||
    length(reduced_cost) != ncol(design) || length(optimum) != 1L ||
    any(!is.finite(c(solution, dual, auxiliary_primal,
                     reduced_cost, optimum)))
  if (malformed) {
    return(list(
      success = FALSE, certified = FALSE,
      status = "lf_support_solution_malformed",
      message = paste0("The LF-support ", sense,
                       " LP returned a malformed nominal optimum."),
      warnings = captured_warnings, sense = sense, maximize = maximize,
      raw_status = raw_status, raw_status_class = raw_class,
      endpoint = NA_real_
    ))
  }

  activity <- as.numeric(constraint_matrix %*% solution)
  primal_objective <- sum(objective * solution)
  dual_objective <- sum(rhs * dual)
  stationarity <- as.numeric(t(constraint_matrix) %*% dual - objective)
  primal_residual <- max(0, max(activity - rhs))
  auxiliary_primal_residual <- max(abs(auxiliary_primal - activity))
  dual_sign_residual <- if (maximize) {
    max(0, -min(dual))
  } else {
    max(0, max(dual))
  }
  dual_stationarity_residual <- max(abs(stationarity))
  reduced_cost_residual <- max(abs(reduced_cost))
  complementarity_residual <- max(abs(dual * (activity - rhs)))
  primal_objective_residual <- abs(primal_objective - optimum)
  dual_objective_residual <- abs(dual_objective - optimum)
  primal_scaled_residual <- max(
    pmax(0, activity - rhs) /
      pmax(1, abs(activity), abs(rhs))
  )
  auxiliary_primal_scaled_residual <- max(
    abs(auxiliary_primal - activity) /
      pmax(1, abs(auxiliary_primal), abs(activity))
  )
  dual_sign_scaled_residual <- if (maximize) {
    max(pmax(0, -dual) / pmax(1, abs(dual)))
  } else {
    max(pmax(0, dual) / pmax(1, abs(dual)))
  }
  stationarity_activity <- as.numeric(t(constraint_matrix) %*% dual)
  dual_stationarity_scaled_residual <- max(
    abs(stationarity) /
      pmax(1, abs(stationarity_activity), abs(objective))
  )
  reduced_cost_scaled_residual <- max(
    abs(reduced_cost) / pmax(1, abs(reduced_cost), abs(objective))
  )
  complementarity_scaled_residual <- max(
    abs(dual * (activity - rhs)) /
      pmax(1, abs(dual * activity), abs(dual * rhs))
  )
  primal_objective_scaled_residual <- primal_objective_residual /
    max(1, abs(primal_objective), abs(optimum))
  dual_objective_scaled_residual <- dual_objective_residual /
    max(1, abs(dual_objective), abs(optimum))
  primal_scale <- max(1, abs(activity), abs(rhs), abs(auxiliary_primal))
  dual_scale <- max(1, abs(objective), abs(dual), abs(stationarity),
                    abs(reduced_cost))
  objective_scale <- max(1, abs(optimum), abs(primal_objective),
                         abs(dual_objective),
                         abs(dual * (activity - rhs)))
  primal_tolerance <- certificate_tolerance * primal_scale
  dual_tolerance <- certificate_tolerance * dual_scale
  objective_tolerance <- certificate_tolerance * objective_scale
  certified <- primal_scaled_residual <= certificate_tolerance &&
    auxiliary_primal_scaled_residual <= certificate_tolerance &&
    dual_sign_scaled_residual <= certificate_tolerance &&
    dual_stationarity_scaled_residual <= certificate_tolerance &&
    reduced_cost_scaled_residual <= certificate_tolerance &&
    complementarity_scaled_residual <= certificate_tolerance &&
    primal_objective_scaled_residual <= certificate_tolerance &&
    dual_objective_scaled_residual <= certificate_tolerance
  if (!certified) {
    return(list(
      success = FALSE, certified = FALSE,
      status = "lf_support_certificate_failure",
      message = paste0("The LF-support ", sense,
                       " LP failed its primal/dual/objective certificate."),
      warnings = captured_warnings, sense = sense, maximize = maximize,
      raw_status = raw_status, raw_status_class = raw_class,
      endpoint = NA_real_, optimum = optimum, solution = solution,
      dual = dual, primal_residual = primal_residual,
      auxiliary_primal_residual = auxiliary_primal_residual,
      dual_sign_residual = dual_sign_residual,
      dual_stationarity_residual = dual_stationarity_residual,
      reduced_cost_residual = reduced_cost_residual,
      complementarity_residual = complementarity_residual,
      primal_objective_residual = primal_objective_residual,
      dual_objective_residual = dual_objective_residual,
      primal_scaled_residual = primal_scaled_residual,
      auxiliary_primal_scaled_residual = auxiliary_primal_scaled_residual,
      dual_sign_scaled_residual = dual_sign_scaled_residual,
      dual_stationarity_scaled_residual =
        dual_stationarity_scaled_residual,
      reduced_cost_scaled_residual = reduced_cost_scaled_residual,
      complementarity_scaled_residual = complementarity_scaled_residual,
      primal_objective_scaled_residual =
        primal_objective_scaled_residual,
      dual_objective_scaled_residual = dual_objective_scaled_residual,
      primal_tolerance = primal_tolerance, dual_tolerance = dual_tolerance,
      objective_tolerance = objective_tolerance
    ))
  }

  # Widen the numerical primal/dual bracket by the registered certificate
  # tolerance.  This makes the support used by inversion conservative relative
  # to the independently checked floating-point LP optimum.
  objective_values <- c(optimum, primal_objective, dual_objective)
  endpoint <- if (maximize) {
    max(objective_values) + objective_tolerance
  } else {
    min(objective_values) - objective_tolerance
  }
  list(
    success = TRUE, certified = TRUE,
    status = if (length(captured_warnings)) "ok_with_warning" else "ok",
    message = NULL, warnings = captured_warnings,
    sense = sense, maximize = maximize,
    raw_status = raw_status, raw_status_class = raw_class,
    endpoint = endpoint, optimum = optimum,
    primal_objective = primal_objective, dual_objective = dual_objective,
    solution = solution, dual = dual,
    primal_residual = primal_residual,
    auxiliary_primal_residual = auxiliary_primal_residual,
    dual_sign_residual = dual_sign_residual,
    dual_stationarity_residual = dual_stationarity_residual,
    reduced_cost_residual = reduced_cost_residual,
    complementarity_residual = complementarity_residual,
    primal_objective_residual = primal_objective_residual,
    dual_objective_residual = dual_objective_residual,
    primal_scaled_residual = primal_scaled_residual,
    auxiliary_primal_scaled_residual = auxiliary_primal_scaled_residual,
    dual_sign_scaled_residual = dual_sign_scaled_residual,
    dual_stationarity_scaled_residual = dual_stationarity_scaled_residual,
    reduced_cost_scaled_residual = reduced_cost_scaled_residual,
    complementarity_scaled_residual = complementarity_scaled_residual,
    primal_objective_scaled_residual = primal_objective_scaled_residual,
    dual_objective_scaled_residual = dual_objective_scaled_residual,
    primal_tolerance = primal_tolerance, dual_tolerance = dual_tolerance,
    objective_tolerance = objective_tolerance,
    n_constraints = length(rhs), n_variables = ncol(design),
    theta = solution[1L], nuisance = solution[-1L]
  )
}

v3_rm_lf_component_support <- function(problem, certificate_tolerance = 1e-8,
                                       solver_fn = NULL) {
  lower <- v3_rm_lf_support_endpoint(
    problem, maximize = FALSE,
    certificate_tolerance = certificate_tolerance, solver_fn = solver_fn
  )
  if (!isTRUE(lower$success)) {
    return(list(
      success = FALSE, certified = FALSE, status = lower$status,
      message = lower$message, warnings = lower$warnings,
      lower = NA_real_, upper = NA_real_, minimum = lower, maximum = NULL
    ))
  }
  upper <- v3_rm_lf_support_endpoint(
    problem, maximize = TRUE,
    certificate_tolerance = certificate_tolerance, solver_fn = solver_fn
  )
  warnings <- unique(c(lower$warnings, upper$warnings))
  if (!isTRUE(upper$success)) {
    return(list(
      success = FALSE, certified = FALSE, status = upper$status,
      message = upper$message, warnings = warnings,
      lower = lower$endpoint, upper = NA_real_,
      minimum = lower, maximum = upper
    ))
  }
  ordered <- is.finite(lower$endpoint) && is.finite(upper$endpoint) &&
    lower$endpoint <= upper$endpoint
  if (!ordered) {
    return(list(
      success = FALSE, certified = FALSE,
      status = "lf_support_endpoint_order_failure",
      message = "Certified LF-support endpoints are nonfinite or unordered.",
      warnings = warnings, lower = lower$endpoint, upper = upper$endpoint,
      minimum = lower, maximum = upper
    ))
  }
  list(
    success = TRUE, certified = TRUE,
    status = if (length(warnings)) "ok_with_warning" else "ok",
    message = NULL, warnings = warnings,
    lower = lower$endpoint, upper = upper$endpoint,
    width = upper$endpoint - lower$endpoint,
    minimum = lower, maximum = upper
  )
}

# Build every signed-component ARP problem, LF critical value, and global
# profile anchor once.  These objects are shared by all point tests in the
# union inversion and its topology audit.
v3_prepare_rm_union <- function(
    representation, Mbar, alpha = 0.05,
    hybrid_kappa = alpha / 10, lf_sims = 1000L, lf_seed = 0L,
    lf_cv_by_component = NULL, project_root) {
  representation <- v3_align_representation(representation)
  dependency_root <- file.path(project_root, "dependencies")
  sd_inference_load_dependencies(dependency_root)
  scenarios <- v3_rm_scenario_table(representation$coordinate_info, Mbar, pre_window = representation[["restriction_window"]])
  n_scenarios <- nrow(scenarios$table)
  if (!is.null(lf_cv_by_component)) {
    v3_assert(length(lf_cv_by_component) == n_scenarios,
              "lf_cv_by_component must have one entry per RM scenario.")
  }
  preferred <- v3_rm_preferred_component(scenarios, representation$betahat)
  component_order <- c(preferred, setdiff(seq_len(n_scenarios), preferred))
  prepared <- vector("list", n_scenarios)
  diagnostics <- vector("list", n_scenarios)
  lf_values <- rep(NA_real_, n_scenarios)
  all_warnings <- character()

  failure <- function(j, status, message, warnings = character(), extra = list()) {
    out <- c(list(
      success = FALSE, status = status, message = message,
      warnings = unique(c(all_warnings, warnings)), failed_component = j,
      scenarios = scenarios, preferred_component = preferred,
      component_order = component_order,
      component_diagnostics = if (j > 1L) {
        do.call(rbind, diagnostics[seq_len(j - 1L)])
      } else data.frame(),
      lf_cv_by_component = lf_values
    ), extra)
    out
  }

  for (j in seq_len(n_scenarios)) {
    restriction <- v3_rm_component_restriction(
      representation$coordinate_info, Mbar,
      scenarios$table[j, , drop = FALSE], scenarios$operator
    )
    built <- v3_reduce_polyhedral_representation(representation, restriction)
    geometry <- v3_make_arp_geometry(built, scope = "post_connected")
    supplied_lf <- if (is.null(lf_cv_by_component)) {
      NULL
    } else if (is.list(lf_cv_by_component)) {
      lf_cv_by_component[[j]]
    } else {
      lf_cv_by_component[j]
    }
    lf <- supplied_lf
    lf_source <- if (is.null(supplied_lf)) "calibrated" else if (
      inherits(supplied_lf, "sd_lf_critical_value")
    ) "cached_diagnostic_object" else "cached_numeric"
    if (is.null(lf)) {
      lf <- sd_lf_critical_value(
        geometry = geometry, Sigma = built$Sigma_bar,
        hybrid_kappa = hybrid_kappa, sims = as.integer(lf_sims),
        seed = as.integer(lf_seed + j), max_failure_fraction = 0,
        engine = "diagnostic", project_root = dependency_root
      )
    }
    numeric_lf <- is.numeric(lf) && length(lf) == 1L && is.finite(lf)
    object_lf <- inherits(lf, "sd_lf_critical_value") &&
      isTRUE(lf$success) && is.finite(lf$critical_value)
    if (!numeric_lf && !object_lf) {
      return(failure(
        j, if (is.list(lf)) lf$status %||%
          "lf_critical_value_unavailable" else
          "lf_critical_value_unavailable",
        if (is.list(lf)) lf$message %||%
          "The RM component LF critical value is invalid." else
          "The RM component LF critical value is invalid.",
        if (is.list(lf)) lf$warnings %||% character() else character(),
        list(lf = lf, geometry = geometry)
      ))
    }
    lf_value <- if (numeric_lf) as.numeric(lf) else lf$critical_value
    problem <- tryCatch(
      sd_bind_arp_problem(
        geometry = geometry, betahat = built$beta_bar,
        Sigma = built$Sigma_bar, procedure = "c_lf", alpha = alpha,
        hybrid_kappa = hybrid_kappa, lf_cv = lf,
        project_root = dependency_root
      ),
      error = function(e) list(
        success = FALSE, status = "rm_problem_binding_error",
        message = conditionMessage(e), warnings = character()
      )
    )
    if (!isTRUE(problem$success)) {
      return(failure(
        j, problem$status %||% "rm_problem_binding_failure",
        problem$message %||% "The RM component ARP problem could not be bound.",
        problem$warnings %||% character(),
        list(problem = problem, geometry = geometry, lf = lf)
      ))
    }
    anchor <- sd_arp_profile_anchor(problem)
    if (!isTRUE(anchor$success)) {
      return(failure(
        j, anchor$status %||% "arp_profile_anchor_unavailable",
        anchor$message %||% "The RM component profile anchor is unavailable.",
        anchor$warnings %||% character(),
        list(problem = problem, profile_anchor = anchor, lf = lf)
      ))
    }
    certificate_tolerance <- max(
      1e-10, anchor$dual_tolerance %||% 0,
      anchor$dual_objective_residual %||% 0
    )
    empty_certified <- is.finite(anchor$eta_lower_bound) &&
      anchor$eta_lower_bound > lf_value + certificate_tolerance
    lf_support <- v3_rm_lf_component_support(problem)
    support_no_feasible <- !isTRUE(lf_support$success) &&
      identical(lf_support$status, "lf_support_no_feasible")
    if (support_no_feasible && empty_certified) {
      raw_support <- lf_support
      lf_support <- list(
        success = TRUE, certified = TRUE, empty = TRUE,
        status = "empty_certified", message = NULL,
        warnings = raw_support$warnings %||% character(),
        lower = NA_real_, upper = NA_real_, width = 0,
        minimum = raw_support$minimum,
        maximum = raw_support$maximum
      )
    } else if (!isTRUE(lf_support$success)) {
      return(failure(
        j, lf_support$status %||% "lf_support_failure",
        lf_support$message %||%
          "The RM component LF support could not be certified.",
        lf_support$warnings %||% character(),
        list(problem = problem, geometry = geometry, lf = lf,
             profile_anchor = anchor, lf_support = lf_support)
      ))
    } else if (empty_certified) {
      return(failure(
        j, "lf_support_anchor_inconsistency",
        paste(
          "The profile dual certifies empty LF support, but the endpoint LP",
          "returned a certified feasible solution."
        ),
        unique(c(anchor$warnings %||% character(),
                 lf_support$warnings %||% character())),
        list(problem = problem, geometry = geometry, lf = lf,
             profile_anchor = anchor, lf_support = lf_support)
      ))
    } else {
      lf_support$empty <- FALSE
    }
    lf_min <- lf_support$minimum %||% list()
    lf_max <- lf_support$maximum %||% list()
    max_finite <- function(x) {
      x <- as.numeric(x)
      if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
    }
    lf_values[j] <- lf_value
    all_warnings <- unique(c(all_warnings, problem$warnings,
                             anchor$warnings %||% character()))
    prepared[[j]] <- list(
      component = j, scenario = scenarios$table[j, , drop = FALSE],
      restriction = restriction, built = built, geometry = geometry,
      lf = lf, lf_value = lf_value, problem = problem,
      lf_support = lf_support,
      profile_anchor = anchor, empty_certified = empty_certified,
      empty_certificate_tolerance = certificate_tolerance
    )
    diagnostics[[j]] <- data.frame(
      component = j,
      scenario = v3_rm_component_label(scenarios$table[j, , drop = FALSE]),
      preferred = j == preferred,
      evaluation_order = match(j, component_order),
      preparation_success = TRUE,
      preparation_status = problem$status %||% "ok",
      lf_source = lf_source,
      lf_critical_value = lf_value,
      lf_sims = if (object_lf) as.integer(lf$sims %||% NA_integer_) else
        NA_integer_,
      lf_seed = if (object_lf) as.integer(lf$seed %||% NA_integer_) else
        NA_integer_,
      lf_failures = if (object_lf) as.integer(lf$n_failures %||% NA_integer_) else
        NA_integer_,
      n_moments = nrow(geometry$A),
      n_arp_moments = length(geometry$rows_for_arp),
      profile_theta = anchor$theta,
      profile_eta = anchor$eta,
      profile_eta_lower_bound = anchor$eta_lower_bound,
      empty_certified = empty_certified,
      lf_support_certified = isTRUE(lf_support$certified),
      lf_support_empty = isTRUE(lf_support$empty),
      lf_support_lower = lf_support$lower,
      lf_support_upper = lf_support$upper,
      lf_support_width = lf_support$width,
      lf_support_min_raw_status = lf_min$raw_status %||% NA_integer_,
      lf_support_min_status_class =
        lf_min$raw_status_class %||% NA_character_,
      lf_support_min_primal_residual =
        lf_min$primal_residual %||% NA_real_,
      lf_support_min_dual_residual =
        lf_min$dual_stationarity_residual %||% NA_real_,
      lf_support_min_objective_residual =
        max_finite(c(lf_min$primal_objective_residual %||% NA_real_,
                     lf_min$dual_objective_residual %||% NA_real_)),
      lf_support_max_raw_status = lf_max$raw_status %||% NA_integer_,
      lf_support_max_status_class =
        lf_max$raw_status_class %||% NA_character_,
      lf_support_max_primal_residual =
        lf_max$primal_residual %||% NA_real_,
      lf_support_max_dual_residual =
        lf_max$dual_stationarity_residual %||% NA_real_,
      lf_support_max_objective_residual =
        max_finite(c(lf_max$primal_objective_residual %||% NA_real_,
                     lf_max$dual_objective_residual %||% NA_real_)),
      stringsAsFactors = FALSE
    )
  }
  list(
    success = TRUE, status = "ok", message = NULL,
    warnings = all_warnings, representation = representation,
    scenarios = scenarios, components = prepared,
    preferred_component = preferred, component_order = component_order,
    component_diagnostics = do.call(rbind, diagnostics),
    lf_cv_by_component = lf_values,
    lf_support_by_component = lapply(prepared, `[[`, "lf_support")
  )
}

# A full C-LF acceptance is necessarily inside its signed component's
# certified first-stage LF support.  This gate can therefore avoid an ARP
# solve without imposing any shape assumption on acceptance inside support.
v3_rm_component_support_gate <- function(component, theta) {
  support <- component$lf_support %||% NULL
  if (is.null(support) || !isTRUE(support$certified)) {
    return(list(skip = FALSE, reason = "unavailable", guard = NA_real_))
  }
  if (isTRUE(support$empty)) {
    return(list(skip = TRUE, reason = "empty", guard = NA_real_))
  }
  if (!all(is.finite(c(support$lower, support$upper))) ||
      support$lower > support$upper) {
    return(list(skip = FALSE, reason = "malformed", guard = NA_real_))
  }
  endpoint_tolerances <- c(
    support$minimum$objective_tolerance %||% NA_real_,
    support$maximum$objective_tolerance %||% NA_real_
  )
  endpoint_tolerances <- endpoint_tolerances[is.finite(endpoint_tolerances)]
  scale <- max(1, abs(theta), abs(support$lower), abs(support$upper))
  guard <- max(c(1e-10 * scale, endpoint_tolerances), na.rm = TRUE)
  outside <- theta < support$lower - guard ||
    theta > support$upper + guard
  list(
    skip = outside,
    reason = if (outside) "outside" else "inside_guarded_support",
    guard = guard
  )
}

# Exact point test for a logical union.  Solver failures are remembered while
# later components are tried: a later acceptance proves union membership, but
# if every component fails or rejects, any solver failure makes the union test
# fail closed rather than treating the failed component as a rejection.
v3_rm_union_point_test <- function(prepared, theta,
                                   component_test = sd_arp_truth_test,
                                   use_lf_support = TRUE) {
  v3_assert(isTRUE(prepared$success),
            "prepared must be a successful v3_prepare_rm_union() result.")
  v3_assert(is.function(component_test), "component_test must be a function.")
  theta <- as.numeric(theta)
  v3_assert(length(theta) == 1L && is.finite(theta),
            "theta must be one finite scalar.")
  failures <- character()
  failed_ids <- integer()
  warnings <- character()
  evaluated_ids <- integer()
  support_skipped_ids <- integer()
  empty_skipped_ids <- integer()
  for (j in prepared$component_order) {
    gate <- if (isTRUE(use_lf_support)) {
      v3_rm_component_support_gate(prepared$components[[j]], theta)
    } else {
      list(skip = FALSE, reason = "disabled", guard = NA_real_)
    }
    if (isTRUE(gate$skip)) {
      if (identical(gate$reason, "empty")) {
        empty_skipped_ids <- c(empty_skipped_ids, j)
      } else {
        support_skipped_ids <- c(support_skipped_ids, j)
      }
      next
    }
    evaluated_ids <- c(evaluated_ids, j)
    tested <- tryCatch(
      component_test(prepared$components[[j]]$problem, theta),
      error = function(e) list(
        success = FALSE, accept = NA,
        status = "rm_component_test_error", message = conditionMessage(e),
        warnings = character()
      )
    )
    warnings <- unique(c(warnings, tested$warnings %||% character()))
    valid <- isTRUE(tested$success) && is.logical(tested$accept) &&
      length(tested$accept) == 1L && !is.na(tested$accept)
    if (valid && isTRUE(tested$accept)) {
      if (length(failures)) {
        warnings <- unique(c(
          warnings,
          paste0("Earlier RM component failure: ", unique(failures))
        ))
      }
      return(list(
        success = TRUE, status = if (length(failed_ids))
          "accepted_after_component_failure" else "ok",
        message = NULL, warnings = warnings, theta = theta, accept = TRUE,
        evaluated_components = length(evaluated_ids),
        evaluated_component_ids = evaluated_ids,
        support_skipped_components = length(support_skipped_ids),
        support_skipped_component_ids = support_skipped_ids,
        empty_skipped_components = length(empty_skipped_ids),
        empty_skipped_component_ids = empty_skipped_ids,
        lf_support_skipped_components =
          length(support_skipped_ids) + length(empty_skipped_ids),
        accepting_component = j, failed_components = length(failed_ids),
        failed_component_ids = failed_ids
      ))
    }
    if (!valid) {
      failed_ids <- c(failed_ids, j)
      failures <- c(failures, tested$message %||%
                      "A signed RM component point test failed.")
    }
  }
  if (length(failed_ids)) {
    return(list(
      success = FALSE, status = "rm_union_component_failure",
      message = paste(unique(failures), collapse = " | "),
      warnings = warnings, theta = theta, accept = NA,
      evaluated_components = length(evaluated_ids),
      evaluated_component_ids = evaluated_ids,
      support_skipped_components = length(support_skipped_ids),
      support_skipped_component_ids = support_skipped_ids,
      empty_skipped_components = length(empty_skipped_ids),
      empty_skipped_component_ids = empty_skipped_ids,
      lf_support_skipped_components =
        length(support_skipped_ids) + length(empty_skipped_ids),
      accepting_component = NA_integer_,
      failed_components = length(failed_ids),
      failed_component_ids = failed_ids
    ))
  }
  list(
    success = TRUE, status = "ok", message = NULL, warnings = warnings,
    theta = theta, accept = FALSE,
    evaluated_components = length(evaluated_ids),
    evaluated_component_ids = evaluated_ids,
    support_skipped_components = length(support_skipped_ids),
    support_skipped_component_ids = support_skipped_ids,
    empty_skipped_components = length(empty_skipped_ids),
    empty_skipped_component_ids = empty_skipped_ids,
    lf_support_skipped_components =
      length(support_skipped_ids) + length(empty_skipped_ids),
    accepting_component = NA_integer_, failed_components = 0L,
    failed_component_ids = integer()
  )
}

v3_memoize_rm_union_test <- function(prepared,
                                     component_test = sd_arp_truth_test,
                                     use_lf_support = TRUE, test_budget = NULL) {
  state <- new.env(parent = emptyenv())
  state$results <- new.env(hash = TRUE, parent = emptyenv())
  state$requests <- 0L
  state$hits <- 0L
  test_fn <- function(theta) {
    state$requests <- state$requests + 1L
    key <- sprintf("%.17g", theta)
    if (exists(key, envir = state$results, inherits = FALSE)) {
      state$hits <- state$hits + 1L
      return(get(key, envir = state$results, inherits = FALSE))
    }
    run <- function() v3_rm_union_point_test(
      prepared, theta, component_test = component_test,
      use_lf_support = use_lf_support
    )
    result <- if (is.null(test_budget)) run() else test_budget$evaluate(theta, run)
    assign(key, result, envir = state$results)
    result
  }
  list(test_fn = test_fn, state = state)
}

v3_rm_point_cache_diagnostics <- function(memoized) {
  keys <- ls(envir = memoized$state$results, all.names = TRUE)
  values <- lapply(keys, function(key) get(
    key, envir = memoized$state$results, inherits = FALSE
  ))
  points <- if (length(values)) {
    do.call(rbind, lapply(values, function(one) data.frame(
      theta = as.numeric(one$theta), success = isTRUE(one$success),
      accept = if (isTRUE(one$success)) isTRUE(one$accept) else NA,
      evaluated_components = as.integer(one$evaluated_components),
      evaluated_component_ids = paste(one$evaluated_component_ids,
                                      collapse = ";"),
      support_skipped_components =
        as.integer(one$support_skipped_components %||% 0L),
      support_skipped_component_ids = paste(
        one$support_skipped_component_ids %||% integer(), collapse = ";"
      ),
      empty_skipped_components =
        as.integer(one$empty_skipped_components %||% 0L),
      empty_skipped_component_ids = paste(
        one$empty_skipped_component_ids %||% integer(), collapse = ";"
      ),
      lf_support_skipped_components =
        as.integer(one$lf_support_skipped_components %||% 0L),
      accepting_component = as.integer(one$accepting_component %||% NA_integer_),
      failed_components = as.integer(one$failed_components %||% 0L),
      failed_component_ids = paste(one$failed_component_ids %||% integer(),
                                   collapse = ";"),
      status = one$status %||% "missing",
      message = one$message %||% "",
      stringsAsFactors = FALSE
    )))
  } else {
    data.frame(
      theta = numeric(), success = logical(), accept = logical(),
      evaluated_components = integer(), evaluated_component_ids = character(),
      support_skipped_components = integer(),
      support_skipped_component_ids = character(),
      empty_skipped_components = integer(),
      empty_skipped_component_ids = character(),
      lf_support_skipped_components = integer(),
      accepting_component = integer(), failed_components = integer(),
      failed_component_ids = character(), status = character(),
      message = character(), stringsAsFactors = FALSE
    )
  }
  if (nrow(points)) points <- points[order(points$theta), , drop = FALSE]
  summary <- list(
    requests = as.integer(memoized$state$requests),
    cache_hits = as.integer(memoized$state$hits),
    unique_point_tests = nrow(points),
    total_component_evaluations = if (nrow(points))
      sum(points$evaluated_components) else 0L,
    total_outside_support_skips = if (nrow(points))
      sum(points$support_skipped_components) else 0L,
    total_empty_support_skips = if (nrow(points))
      sum(points$empty_skipped_components) else 0L,
    total_lf_support_skips = if (nrow(points))
      sum(points$lf_support_skipped_components) else 0L,
    mean_components_per_point = if (nrow(points))
      mean(points$evaluated_components) else 0,
    max_components_per_point = if (nrow(points))
      max(points$evaluated_components) else 0L
  )
  list(points = points, summary = summary, raw = values)
}

v3_rm_lf_support_hull <- function(prepared) {
  v3_assert(isTRUE(prepared$success) && length(prepared$components),
            "prepared must contain successful RM components.")
  active <- Filter(function(one) !isTRUE(one$lf_support$empty),
                   prepared$components)
  if (!length(active)) {
    return(list(
      success = TRUE, certified = TRUE, empty = TRUE,
      status = "empty_certified", message = NULL,
      lower = NA_real_, upper = NA_real_, width = 0,
      components = data.frame(
        component = integer(), lower = numeric(), upper = numeric()
      )
    ))
  }
  pieces <- do.call(rbind, lapply(active, function(one) {
    support <- one$lf_support
    v3_assert(isTRUE(support$success) && isTRUE(support$certified) &&
                all(is.finite(c(support$lower, support$upper))) &&
                support$lower <= support$upper,
              "Every RM component must have certified finite LF support.")
    data.frame(
      component = as.integer(one$component),
      lower = as.numeric(support$lower),
      upper = as.numeric(support$upper),
      stringsAsFactors = FALSE
    )
  }))
  list(
    success = TRUE, certified = TRUE, empty = FALSE,
    status = "ok", message = NULL,
    lower = min(pieces$lower), upper = max(pieces$upper),
    width = max(pieces$upper) - min(pieces$lower),
    components = pieces
  )
}

v3_rm_union_resolution_grids <- function(
    support, profile_anchors, target_se,
    initial_radius_se = 8, mesh_se = 0.25,
    endpoint_tolerance = 1e-8, max_grid_points = .Machine$integer.max) {
  v3_assert(isTRUE(support$success) && isTRUE(support$certified) &&
              all(is.finite(c(support$lower, support$upper))) &&
              support$lower <= support$upper,
            "The RM union LF-support hull must be finite and certified.")
  profile_anchors <- as.numeric(profile_anchors)
  v3_assert(length(profile_anchors) && all(is.finite(profile_anchors)),
            "Every RM component must have a finite profile anchor.")
  v3_assert(length(target_se) == 1L && is.finite(target_se) && target_se > 0,
            "The RM target standard error must be positive.")
  v3_assert(length(initial_radius_se) == 1L &&
              is.finite(initial_radius_se) && initial_radius_se >= 0,
            "The RM topology radius must be finite and nonnegative.")
  v3_assert(length(mesh_se) == 1L && is.finite(mesh_se) && mesh_se > 0,
            "The RM resolution mesh multiplier must be positive.")
  v3_assert(length(endpoint_tolerance) == 1L &&
              is.finite(endpoint_tolerance) && endpoint_tolerance > 0,
            "The RM endpoint tolerance must be positive.")
  registered_mesh <- mesh_se * target_se
  v3_assert(is.finite(registered_mesh) && registered_mesh > 0,
            "The RM registered resolution mesh must be positive.")
  numeric_scale <- max(1, abs(support$lower), abs(support$upper), target_se)
  sentinel_margin <- max(
    initial_radius_se * target_se, registered_mesh,
    16 * endpoint_tolerance, sqrt(.Machine$double.eps) * numeric_scale
  )
  sentinels <- c(
    lower = support$lower - sentinel_margin,
    upper = support$upper + sentinel_margin
  )
  v3_assert(all(is.finite(sentinels)) && sentinels[1L] < support$lower &&
              sentinels[2L] > support$upper,
            "The RM LF-support sentinels must be finite and strictly outside.")

  regular <- function(step, shifted = FALSE) {
    span <- diff(sentinels)
    requested_count <- if (shifted) floor(span / step) else ceiling(span / step)
    if (!is.finite(requested_count) || requested_count >= .Machine$integer.max ||
        requested_count > max_grid_points) {
      stop(structure(list(message = "The RM refinement grid exceeds the point-evaluation budget or R's integer range.",
        call = NULL), class = c("hc_search_budget_exhausted", "error", "condition")))
    }
    if (shifted) {
      count <- max(1L, as.integer(floor(span / step)))
      points <- sentinels[1L] + (seq_len(count) - 0.5) * step
      points <- points[points > sentinels[1L] & points < sentinels[2L]]
    } else {
      count <- max(1L, as.integer(ceiling(span / step)))
      points <- sentinels[1L] + seq.int(0L, count) * step
      points <- points[points >= sentinels[1L] & points <= sentinels[2L]]
    }
    c(sentinels, points)
  }
  mandatory <- c(
    support$components$lower, support$components$upper, profile_anchors
  )
  mandatory <- mandatory[
    is.finite(mandatory) & mandatory >= sentinels[1L] &
      mandatory <= sentinels[2L]
  ]
  make_grid <- function(step, shifted = FALSE) sort(unique(c(
    regular(step, shifted = shifted), mandatory
  )))
  list(
    success = TRUE, status = "ok", message = NULL,
    registered_mesh = registered_mesh,
    fine_registered_mesh = registered_mesh / 2,
    ultrafine_registered_mesh = registered_mesh / 4,
    sentinel_margin = sentinel_margin, sentinels = sentinels,
    base = make_grid(registered_mesh),
    shifted = make_grid(registered_mesh, shifted = TRUE),
    fine = make_grid(registered_mesh / 2),
    ultrafine = make_grid(registered_mesh / 4),
    support = support
  )
}

v3_rm_lf_point_exclusion_certificate <- function(
    problem, theta, certificate_tolerance = 1e-8) {
  v3_assert(inherits(problem, "sd_arp_problem") && isTRUE(problem$success),
            "problem must be a successful sd_arp_problem.")
  theta <- as.numeric(theta)
  v3_assert(length(theta) == 1L && is.finite(theta),
            "theta must be one finite scalar.")
  lf_cv <- as.numeric(problem$hybrid_list$lf_cv)
  fixed_problem <- problem
  fixed_problem$base_y <- problem$base_y - problem$geometry$a_one * theta
  # Reuse the fully certified free-variable endpoint solver.  Replacing the
  # first design column by sd makes its first variable eta, so its standardized
  # constraints are -eta - (a_minus/sd) delta <= -y(theta)/sd.
  fixed_problem$geometry$a_one <- sqrt(diag(problem$sigma_y))
  fixed_problem$hybrid_list$lf_cv <- 0
  certified_lp <- v3_rm_lf_support_endpoint(
    fixed_problem, maximize = FALSE,
    certificate_tolerance = certificate_tolerance
  )
  if (!isTRUE(certified_lp$success)) {
    return(list(
      success = FALSE, certified = FALSE, excluded = NA,
      status = certified_lp$status,
      message = certified_lp$message, warnings = certified_lp$warnings,
      theta = theta, lf_critical_value = lf_cv,
      eta = NA_real_, eta_lower_bound = NA_real_, lp = certified_lp
    ))
  }
  guard <- max(
    1e-10 * max(1, abs(lf_cv), abs(certified_lp$optimum)),
    certified_lp$objective_tolerance
  )
  excluded <- certified_lp$endpoint > lf_cv + guard
  list(
    success = TRUE, certified = TRUE, excluded = excluded,
    status = if (excluded) "excluded" else "not_excluded",
    message = if (excluded) NULL else paste(
      "The fixed-theta eta lower bound does not exceed the LF critical value",
      "by the registered certificate guard."
    ),
    warnings = certified_lp$warnings, theta = theta,
    lf_critical_value = lf_cv,
    eta = certified_lp$optimum,
    eta_lower_bound = certified_lp$endpoint,
    exclusion_guard = guard, lp = certified_lp
  )
}

v3_rm_verify_support_sentinels <- function(prepared, memoized, sentinels, test_budget = NULL) {
  sentinels <- as.numeric(sentinels)
  v3_assert(length(sentinels) == 2L && all(is.finite(sentinels)) &&
              sentinels[1L] < sentinels[2L],
            "Two ordered finite LF-support sentinels are required.")
  exclusion <- lapply(seq_along(prepared$components), function(j) {
    lapply(sentinels, function(theta) {
      v3_rm_lf_point_exclusion_certificate(
        prepared$components[[j]]$problem, theta
      )
    })
  })
  exclusion_rows <- do.call(rbind, lapply(seq_along(exclusion), function(j) {
    do.call(rbind, lapply(seq_along(sentinels), function(k) {
      one <- exclusion[[j]][[k]]
      data.frame(
        component = j, side = c("lower", "upper")[k],
        theta = sentinels[k], success = isTRUE(one$success),
        certified = isTRUE(one$certified),
        excluded = if (isTRUE(one$success)) isTRUE(one$excluded) else NA,
        status = one$status %||% "missing",
        raw_status = one$lp$raw_status %||% NA_integer_,
        raw_status_class = one$lp$raw_status_class %||% NA_character_,
        eta = one$eta %||% NA_real_,
        eta_lower_bound = one$eta_lower_bound %||% NA_real_,
        lf_critical_value = one$lf_critical_value %||% NA_real_,
        exclusion_guard = one$exclusion_guard %||% NA_real_,
        message = one$message %||% "", stringsAsFactors = FALSE
      )
    }))
  }))
  if (!all(exclusion_rows$success & exclusion_rows$certified)) {
    return(list(
      success = FALSE, status = "lf_support_sentinel_certificate_failure",
      message = paste(
        "A fixed-theta eta LP failed to certify an RM LF-support sentinel",
        "against every signed component."
      ),
      certified = FALSE, exclusion_diagnostics = exclusion_rows,
      exclusion_raw = exclusion
    ))
  }
  if (!all(exclusion_rows$excluded)) {
    return(list(
      success = FALSE, status = "lf_support_sentinel_not_excluded",
      message = paste(
        "A fixed-theta eta lower bound did not exclude an RM LF-support",
        "sentinel from every signed component."
      ),
      certified = FALSE, exclusion_diagnostics = exclusion_rows,
      exclusion_raw = exclusion
    ))
  }

  # Keep this as an unoptimized sanity check: every signed component's ordinary
  # C-LF point test is actually run, independently of the support gate used by
  # production grid scans.
  tested <- lapply(sentinels, function(theta) {
    run <- function() v3_rm_union_point_test(
      prepared, theta, component_test = sd_arp_truth_test,
      use_lf_support = FALSE
    )
    if (is.null(test_budget)) run() else test_budget$evaluate(theta, run, "sentinel")
  })
  success <- vapply(tested, function(one) isTRUE(one$success), logical(1L))
  reject <- vapply(tested, function(one) {
    isTRUE(one$success) && identical(one$accept, FALSE)
  }, logical(1L))
  diagnostics <- data.frame(
    side = c("lower", "upper"), theta = sentinels,
    success = success,
    accept = vapply(tested, function(one) {
      if (isTRUE(one$success)) isTRUE(one$accept) else NA
    }, logical(1L)),
    status = vapply(tested, function(one) one$status %||% "missing",
                    character(1L)),
    evaluated_components = vapply(
      tested, function(one) as.integer(one$evaluated_components %||% 0L),
      integer(1L)
    ),
    message = vapply(tested, function(one) one$message %||% "",
                     character(1L)),
    stringsAsFactors = FALSE
  )
  if (!all(success)) {
    return(list(
      success = FALSE, status = "lf_support_sentinel_test_failure",
      message = "A direct RM union test failed at an LF-support sentinel.",
      certified = FALSE, diagnostics = diagnostics, raw = tested,
      exclusion_diagnostics = exclusion_rows, exclusion_raw = exclusion
    ))
  }
  if (!all(reject)) {
    return(list(
      success = FALSE, status = "lf_support_sentinel_acceptance",
      message = paste(
        "The full RM union test accepted outside its certified first-stage",
        "LF-support hull."
      ),
      certified = FALSE, diagnostics = diagnostics, raw = tested,
      exclusion_diagnostics = exclusion_rows, exclusion_raw = exclusion
    ))
  }
  list(
    success = TRUE, status = "ok", message = NULL, certified = TRUE,
    diagnostics = diagnostics, raw = tested,
    exclusion_diagnostics = exclusion_rows, exclusion_raw = exclusion
  )
}

v3_rm_apply_empty_certificate <- function(interval, prepared) {
  empty_certified <- all(vapply(
    prepared$components, function(one) isTRUE(one$empty_certified), logical(1L)
  ))
  interval$empty_certified <- empty_certified
  if (isTRUE(interval$success) && identical(interval$status, "empty") &&
      !empty_certified) {
    interval$success <- FALSE
    interval$complete <- FALSE
    interval$status <- "uncertified_empty_search"
    interval$message <- paste(
      "The finite union resolution grid contained no accepted point, but at",
      "least one signed component lacks a global C-LF empty-set certificate."
    )
  }
  interval
}

v3_rm_union_invert_once <- function(
    prepared, memoized, topology_grid, theta_true,
    endpoint_tolerance, max_expansions = 0L,
    max_bisection = 100L, max_test_evaluations = 100000L) {
  v3_assert(identical(as.integer(max_expansions), 0L),
            "RM LF-support inversion forbids adaptive tail expansion.")
  interval <- sd_invert_test(
    test_fn = memoized$test_fn, topology_grid = topology_grid,
    endpoint_tolerance = endpoint_tolerance,
    max_expansions = as.integer(max_expansions),
    expansion_points = length(topology_grid) - 1L,
    max_bisection = as.integer(max_bisection),
    max_test_evaluations = as.integer(max_test_evaluations),
    truth = theta_true
  )
  v3_rm_apply_empty_certificate(interval, prepared)
}

v3_rm_compare_union_inversions <- function(base, candidate,
                                           endpoint_tolerance,
                                           agreement_multiplier = 4) {
  same_count <- isTRUE(base$success) && isTRUE(candidate$success) &&
    nrow(base$components) == nrow(candidate$components)
  endpoint_gap <- if (same_count && nrow(candidate$components)) {
    max(abs(
      as.matrix(base$components[, c("lower", "upper"), drop = FALSE]) -
        as.matrix(candidate$components[, c("lower", "upper"), drop = FALSE])
    ))
  } else if (same_count) 0 else Inf
  tolerance <- agreement_multiplier * endpoint_tolerance
  list(
    pass = same_count && is.finite(endpoint_gap) && endpoint_gap <= tolerance,
    same_count = same_count, endpoint_gap = endpoint_gap,
    tolerance = tolerance
  )
}

v3_rm_union_resolution_audit <- function(
    prepared, memoized, authoritative_interval, resolution_grids, theta_true,
    endpoint_tolerance,
    max_bisection = 100L, max_test_evaluations = 100000L,
    agreement_multiplier = 4) {
  run_one <- function(one_grid) v3_rm_union_invert_once(
    prepared = prepared, memoized = memoized, topology_grid = one_grid,
    theta_true = theta_true, endpoint_tolerance = endpoint_tolerance,
    max_expansions = 0L, max_bisection = max_bisection,
    max_test_evaluations = max_test_evaluations
  )
  base <- run_one(resolution_grids$base)
  shifted <- run_one(resolution_grids$shifted)
  fine <- run_one(resolution_grids$fine)
  base_comparison <- v3_rm_compare_union_inversions(
    authoritative_interval, base, endpoint_tolerance, agreement_multiplier
  )
  shifted_comparison <- v3_rm_compare_union_inversions(
    authoritative_interval, shifted, endpoint_tolerance, agreement_multiplier
  )
  fine_comparison <- v3_rm_compare_union_inversions(
    authoritative_interval, fine, endpoint_tolerance, agreement_multiplier
  )
  audit_runs_successful <- isTRUE(authoritative_interval$success) &&
    isTRUE(base$success) && isTRUE(shifted$success) && isTRUE(fine$success)
  stable <- audit_runs_successful && base_comparison$pass &&
    shifted_comparison$pass && fine_comparison$pass
  list(
    # `success` records whether every diagnostic inversion ran successfully.
    # Agreement is retained separately in `resolution_stable`; D043 makes the
    # h/4 result authoritative and treats the coarser-grid comparisons as
    # diagnostics rather than a veto on an otherwise successful inversion.
    success = audit_runs_successful,
    status = if (!audit_runs_successful) {
      "resolution_audit_failure"
    } else if (stable) {
      "ok"
    } else {
      "resolution_disagreement"
    },
    message = if (!audit_runs_successful) {
      paste(
        "At least one RM C-LF diagnostic inversion failed while comparing",
        "the h, shifted-h, h/2, and h/4 grids."
      )
    } else if (!stable) {
      paste(
        "Union-level h, shifted-h, h/2, and h/4 inversions did not agree at",
        "the registered endpoint tolerance; the authoritative h/4 result",
        "was retained and the disagreement was recorded."
      )
    } else NULL,
    resolution_stable = stable, certified = FALSE, inversion_exact = FALSE,
    base = base, shifted = shifted, fine = fine,
    ultrafine = authoritative_interval,
    base_comparison = base_comparison,
    shifted_comparison = shifted_comparison,
    fine_comparison = fine_comparison,
    registered_mesh = resolution_grids$registered_mesh,
    fine_registered_mesh = resolution_grids$fine_registered_mesh,
    ultrafine_registered_mesh = resolution_grids$ultrafine_registered_mesh,
    base_grid = resolution_grids$base,
    shifted_grid = resolution_grids$shifted,
    fine_grid = resolution_grids$fine,
    ultrafine_grid = resolution_grids$ultrafine
  )
}

v3_rm_add_point_counts <- function(component_diagnostics, point_cache) {
  raw <- point_cache$raw
  component_diagnostics$point_evaluations <- 0L
  component_diagnostics$point_acceptances <- 0L
  component_diagnostics$point_failures <- 0L
  component_diagnostics$point_outside_support_skips <- 0L
  component_diagnostics$point_empty_support_skips <- 0L
  for (one in raw) {
    evaluated <- one$evaluated_component_ids %||% integer()
    if (length(evaluated)) {
      component_diagnostics$point_evaluations[evaluated] <-
        component_diagnostics$point_evaluations[evaluated] + 1L
    }
    accepted <- one$accepting_component %||% NA_integer_
    if (length(accepted) == 1L && !is.na(accepted)) {
      component_diagnostics$point_acceptances[accepted] <-
        component_diagnostics$point_acceptances[accepted] + 1L
    }
    failed <- one$failed_component_ids %||% integer()
    if (length(failed)) {
      component_diagnostics$point_failures[failed] <-
        component_diagnostics$point_failures[failed] + 1L
    }
    outside_skipped <- one$support_skipped_component_ids %||% integer()
    if (length(outside_skipped)) {
      component_diagnostics$point_outside_support_skips[outside_skipped] <-
        component_diagnostics$point_outside_support_skips[outside_skipped] + 1L
    }
    empty_skipped <- one$empty_skipped_component_ids %||% integer()
    if (length(empty_skipped)) {
      component_diagnostics$point_empty_support_skips[empty_skipped] <-
        component_diagnostics$point_empty_support_skips[empty_skipped] + 1L
    }
  }
  component_diagnostics
}

v3_clf_rm_union <- function(
    representation, Mbar, theta_true = 0, alpha = 0.05,
    hybrid_kappa = alpha / 10, lf_sims = 1000L, lf_seed = 0L,
    lf_cv_by_component = NULL,
    topology_initial_radius_se = 8, topology_mesh_se = 0.25,
    endpoint_tolerance_se = 1e-5, max_expansions = 4L,
    max_bisection = 100L, max_test_evaluations = 100000L,
    topology_audit = TRUE, project_root,
    search_context = NULL, test_budget = NULL) {
  v3_assert(length(max_expansions) == 1L && is.finite(max_expansions) &&
              max_expansions >= 0,
            "max_expansions must be one nonnegative finite value.")
  requested_max_expansions <- as.integer(max_expansions)
  representation <- v3_align_representation(representation)
  scenarios <- v3_rm_scenario_table(representation$coordinate_info, Mbar, pre_window = representation[["restriction_window"]])
  plugin_restriction <- v3_rm_plugin_restriction(
    representation$coordinate_info, representation$betahat, Mbar, pre_window = representation[["restriction_window"]]
  )
  plugin_built <- v3_reduce_polyhedral_representation(
    representation, plugin_restriction
  )
  plugin <- v3_plugin_identified_set(plugin_built, scope = "post_connected")
  target_se <- sqrt(as.numeric(
    t(plugin_built$target_bar) %*% plugin_built$Sigma_bar %*%
      plugin_built$target_bar
  ))
  prepared <- if (!is.null(search_context) &&
      exists("prepared", search_context, inherits = FALSE)) search_context$prepared else v3_prepare_rm_union(
    representation = representation, Mbar = Mbar, alpha = alpha,
    hybrid_kappa = hybrid_kappa, lf_sims = lf_sims, lf_seed = lf_seed,
    lf_cv_by_component = lf_cv_by_component, project_root = project_root
  )
  if (!is.null(search_context)) search_context$prepared <- prepared
  if (!isTRUE(prepared$success)) {
    return(c(prepared, list(
      procedure = "c_lf", scope = "post_connected", plugin = plugin,
      target_se = target_se, inversion_engine = "union",
      inversion_exact = FALSE, resolution_stable = NA,
      topology_certified = FALSE,
      requested_max_expansions = requested_max_expansions,
      max_expansions_used = 0L
    )))
  }
  endpoint_tolerance <- endpoint_tolerance_se * target_se
  lf_support <- v3_rm_lf_support_hull(prepared)
  memoized <- if (!is.null(search_context) &&
      exists("memoized", search_context, inherits = FALSE)) search_context$memoized else
    v3_memoize_rm_union_test(prepared, test_budget = test_budget)
  if (!is.null(search_context)) search_context$memoized <- memoized
  if (isTRUE(lf_support$empty)) {
    truth_test <- memoized$test_fn(theta_true)
    truth_valid <- isTRUE(truth_test$success) &&
      identical(truth_test$accept, FALSE)
    point_cache <- v3_rm_point_cache_diagnostics(memoized)
    diagnostics <- v3_rm_add_point_counts(
      prepared$component_diagnostics, point_cache
    )
    diagnostics$union_resolution_stable <- TRUE
    diagnostics$union_topology_certified <- FALSE
    if (!truth_valid) {
      diagnostics$union_resolution_stable <- NA
      return(list(
        success = FALSE,
        status = if (!isTRUE(truth_test$success))
          "truth_test_failure" else "empty_support_truth_contradiction",
        message = if (!isTRUE(truth_test$success)) {
          truth_test$message %||% "The direct RM truth-point test failed."
        } else {
          paste(
            "The direct RM truth-point test accepted despite certified-empty",
            "LF support."
          )
        },
        warnings = unique(c(prepared$warnings,
                            truth_test$warnings %||% character())),
        procedure = "c_lf", scope = "post_connected", plugin = plugin,
        target_se = target_se, lf_support = lf_support,
        support_certified = FALSE, inversion_engine = "union",
        inversion_exact = FALSE, resolution_stable = NA,
        topology_certified = FALSE,
        requested_max_expansions = requested_max_expansions,
        max_expansions_used = 0L,
        contains_truth = if (isTRUE(truth_test$success))
          truth_test$accept else NA,
        component_diagnostics = diagnostics,
        lf_cv_by_component = prepared$lf_cv_by_component,
        lf_support_by_component = prepared$lf_support_by_component,
        point_test_diagnostics = point_cache$points,
        evaluated_component_counts = point_cache$summary
      ))
    }
    components <- data.frame(lower = numeric(), upper = numeric())
    raw_components <- components
    raw_components$signed_benchmark_component <- integer()
    hull <- v3_component_hull(components)
    interval <- list(
      success = TRUE, complete = TRUE, status = "empty_certified",
      message = NULL, warnings = character(), components = components,
      lower = NA_real_, upper = NA_real_, total_length = 0,
      truncated_length = 0, n_components = 0L,
      boundary_accepted = FALSE, expansions = 0L,
      endpoint_tolerance = endpoint_tolerance,
      truth_test = truth_test, contains_truth = FALSE,
      empty_certified = TRUE
    )
    resolution_check <- list(
      success = TRUE, status = "not_required_empty_support", message = NULL,
      resolution_stable = TRUE, certified = FALSE, inversion_exact = FALSE
    )
    return(list(
      success = TRUE, status = "empty_certified", message = NULL,
      warnings = prepared$warnings,
      procedure = "c_lf", scope = "post_connected", plugin = plugin,
      components = components, raw_components = raw_components, hull = hull,
      contains_truth = FALSE, target_se = target_se, interval = interval,
      topology_grid = numeric(), topology_audit = resolution_check,
      resolution_audit = resolution_check, resolution_grids = NULL,
      lf_support = lf_support, support_sentinel_tests = NULL,
      support_certified = TRUE, topology_certified = FALSE,
      resolution_stable = TRUE, inversion_exact = FALSE,
      inversion_engine = "union",
      requested_max_expansions = requested_max_expansions,
      max_expansions_used = 0L,
      signed_benchmark_components = nrow(scenarios$table),
      signed_benchmark_components_full = scenarios$scenario_count_full,
      preferred_component = prepared$preferred_component,
      component_evaluation_order = prepared$component_order,
      component_diagnostics = diagnostics,
      lf_cv_by_component = prepared$lf_cv_by_component,
      lf_support_by_component = prepared$lf_support_by_component,
      point_test_diagnostics = point_cache$points,
      evaluated_component_counts = point_cache$summary
    ))
  }
  active_components <- as.integer(lf_support$components$component)
  anchors <- vapply(prepared$components[active_components], function(one) {
    as.numeric(one$profile_anchor$theta)
  }, numeric(1L))
  resolution_grids <- v3_rm_union_resolution_grids(
    support = lf_support, profile_anchors = anchors, target_se = target_se,
    initial_radius_se = topology_initial_radius_se,
    mesh_se = topology_mesh_se,
    endpoint_tolerance = endpoint_tolerance, max_grid_points = max_test_evaluations
  )
  # The finest registered grid is authoritative.  The coarser h, shifted-h,
  # and h/2 inversions are resolution-stability diagnostics only.
  topology_grid <- resolution_grids$ultrafine
  sentinel_check <- v3_rm_verify_support_sentinels(
    prepared, memoized, resolution_grids$sentinels, test_budget = test_budget
  )
  if (!isTRUE(sentinel_check$success)) {
    point_cache <- v3_rm_point_cache_diagnostics(memoized)
    diagnostics <- v3_rm_add_point_counts(
      prepared$component_diagnostics, point_cache
    )
    diagnostics$union_resolution_stable <- NA
    diagnostics$union_topology_certified <- FALSE
    return(list(
      success = FALSE, status = sentinel_check$status,
      message = sentinel_check$message,
      warnings = prepared$warnings,
      procedure = "c_lf", scope = "post_connected", plugin = plugin,
      target_se = target_se, topology_grid = topology_grid,
      resolution_grids = resolution_grids,
      lf_support = lf_support,
      support_sentinel_tests = sentinel_check,
      support_certified = FALSE, inversion_engine = "union",
      inversion_exact = FALSE, resolution_stable = NA,
      topology_certified = FALSE,
      requested_max_expansions = requested_max_expansions,
      max_expansions_used = 0L,
      signed_benchmark_components = nrow(scenarios$table),
      signed_benchmark_components_full = scenarios$scenario_count_full,
      preferred_component = prepared$preferred_component,
      component_evaluation_order = prepared$component_order,
      component_diagnostics = diagnostics,
      lf_cv_by_component = prepared$lf_cv_by_component,
      lf_support_by_component = prepared$lf_support_by_component,
      point_test_diagnostics = point_cache$points,
      evaluated_component_counts = point_cache$summary
    ))
  }
  interval <- v3_rm_union_invert_once(
    prepared = prepared, memoized = memoized,
    topology_grid = topology_grid, theta_true = theta_true,
    endpoint_tolerance = endpoint_tolerance,
    max_expansions = 0L, max_bisection = max_bisection,
    max_test_evaluations = max_test_evaluations
  )
  if (!isTRUE(interval$success) || !isTRUE(interval$complete)) {
    point_cache <- v3_rm_point_cache_diagnostics(memoized)
    diagnostics <- v3_rm_add_point_counts(
      prepared$component_diagnostics, point_cache
    )
    diagnostics$union_resolution_stable <- NA
    diagnostics$union_topology_certified <- FALSE
    return(list(
      success = FALSE, status = interval$status %||% "rm_union_failure",
      message = interval$message %||% "The RM union inversion failed.",
      warnings = unique(c(prepared$warnings,
                          interval$warnings %||% character())),
      procedure = "c_lf", scope = "post_connected", plugin = plugin,
      target_se = target_se, interval = interval,
      topology_grid = topology_grid,
      resolution_grids = resolution_grids,
      lf_support = lf_support,
      support_sentinel_tests = sentinel_check,
      support_certified = isTRUE(sentinel_check$certified),
      inversion_engine = "union", inversion_exact = FALSE,
      resolution_stable = NA, topology_certified = FALSE,
      requested_max_expansions = requested_max_expansions,
      max_expansions_used = 0L,
      signed_benchmark_components = nrow(scenarios$table),
      signed_benchmark_components_full = scenarios$scenario_count_full,
      preferred_component = prepared$preferred_component,
      component_evaluation_order = prepared$component_order,
      component_diagnostics = diagnostics,
      lf_cv_by_component = prepared$lf_cv_by_component,
      lf_support_by_component = prepared$lf_support_by_component,
      point_test_diagnostics = point_cache$points,
      evaluated_component_counts = point_cache$summary
    ))
  }
  resolution_check <- if (isTRUE(topology_audit)) {
    v3_rm_union_resolution_audit(
      prepared = prepared, memoized = memoized,
      authoritative_interval = interval, resolution_grids = resolution_grids,
      theta_true = theta_true, endpoint_tolerance = endpoint_tolerance,
      max_bisection = max_bisection,
      max_test_evaluations = max_test_evaluations
    )
  } else {
    list(
      success = TRUE, status = "not_requested", message = NULL,
      resolution_stable = NA, certified = FALSE, inversion_exact = FALSE
    )
  }
  point_cache <- v3_rm_point_cache_diagnostics(memoized)
  diagnostics <- v3_rm_add_point_counts(
    prepared$component_diagnostics, point_cache
  )
  resolution_stable <- if (isTRUE(topology_audit)) {
    isTRUE(resolution_check$resolution_stable)
  } else NA
  diagnostics$union_resolution_stable <- resolution_stable
  diagnostics$union_topology_certified <- FALSE
  components <- interval$components[, c("lower", "upper"), drop = FALSE]
  hull <- v3_component_hull(components)
  raw_components <- components
  raw_components$signed_benchmark_component <- NA_integer_
  resolution_warning <- if (identical(resolution_stable, FALSE)) {
    resolution_check$message %||%
      "The RM C-LF resolution diagnostic disagreed across grids."
  } else {
    character()
  }
  list(
    success = TRUE,
    status = if (identical(resolution_stable, FALSE)) {
      "ok_with_resolution_warning"
    } else {
      interval$status
    },
    message = NULL,
    warnings = unique(c(prepared$warnings,
                        interval$warnings %||% character(),
                        resolution_check$warnings %||% character(),
                        resolution_warning)),
    procedure = "c_lf", scope = "post_connected", plugin = plugin,
    components = components, raw_components = raw_components, hull = hull,
    contains_truth = interval$contains_truth, target_se = target_se,
    interval = interval, topology_grid = topology_grid,
    topology_audit = resolution_check,
    resolution_audit = resolution_check,
    resolution_grids = resolution_grids,
    lf_support = lf_support,
    support_sentinel_tests = sentinel_check,
    support_certified = isTRUE(sentinel_check$certified),
    topology_certified = FALSE,
    resolution_stable = resolution_stable,
    inversion_exact = FALSE,
    inversion_engine = "union",
    requested_max_expansions = requested_max_expansions,
    max_expansions_used = 0L,
    signed_benchmark_components = nrow(scenarios$table),
    signed_benchmark_components_full = scenarios$scenario_count_full,
    preferred_component = prepared$preferred_component,
    component_evaluation_order = prepared$component_order,
    component_diagnostics = diagnostics,
    lf_cv_by_component = prepared$lf_cv_by_component,
    lf_support_by_component = prepared$lf_support_by_component,
    point_test_diagnostics = point_cache$points,
    evaluated_component_counts = point_cache$summary
  )
}

# Public RM adapter.  The bounded pointwise-union path is the production
# default; the former component-wise inversion remains available for audits.
v3_clf_rm <- function(...,
                      inversion_engine = c("union", "componentwise")) {
  inversion_engine <- match.arg(inversion_engine)
  fit <- if (inversion_engine == "union") {
    v3_clf_rm_union(...)
  } else {
    v3_clf_rm_componentwise(...)
  }
  fit$inversion_engine <- inversion_engine
  fit
}

v3_path_in_sd_class <- function(path, coordinate_info, M_raw,
                                scope = c("full_path", "post_connected"),
                                tolerance = 1e-10) {
  scope <- match.arg(scope)
  operator <- v3_second_difference_operator(coordinate_info)
  path <- path[colnames(operator$D2)]
  curvature <- as.numeric(operator$D2 %*% path)
  rows <- if (scope == "full_path") {
    rep(TRUE, length(curvature))
  } else {
    operator$rows$component == "post_connected"
  }
  max(abs(curvature[rows])) <= M_raw + tolerance
}

v3_path_in_rm_class <- function(path, coordinate_info, Mbar,
                                tolerance = 1e-10) {
  operator <- v3_first_difference_operator(coordinate_info)
  path <- path[colnames(operator$D1)]
  difference <- as.numeric(operator$D1 %*% path)
  for (g in sort(unique(operator$rows$cohort))) {
    pre <- operator$rows$cohort == g & operator$rows$component == "pre"
    post <- operator$rows$cohort == g & operator$rows$component == "post"
    benchmark <- max(abs(difference[pre]))
    if (max(abs(difference[post])) > Mbar * benchmark + tolerance) return(FALSE)
  }
  TRUE
}
