# Safe inference wrappers for staggered-adoption sensitivity analysis.
#
# This file deliberately does not modify the vendored HonestDiD or FLCI
# sources.  It loads those functions into a private environment and adds an
# API with explicit numerical diagnostics.
#
# Public API
# ----------
#   sd_inference_load_dependencies(project_root = NULL)
#
#   sd_flci_prepare(model, scope, alpha = 0.05, fit = NULL,
#                   solver = "nested", solve_args = list())
#   sd_flci_interval(rule, betahat, theta = NULL)
#   sd_flci_truth_test(rule, betahat, theta)
#
#   sd_arp_geometry(A, d, target, scope)
#   sd_lf_critical_value(geometry, Sigma, hybrid_kappa,
#                        sims = 10000L, seed = 0L)
#   sd_bind_arp_problem(geometry, betahat, Sigma, procedure,
#                       alpha = 0.05, hybrid_kappa = alpha / 10,
#                       lf_cv = NULL)
#   sd_arp_profile_anchor(problem)
#   sd_arp_truth_test(problem, theta)
#   sd_arp_invert(problem, topology_grid, endpoint_tolerance = 1e-8)
#   sd_arp_audit_endpoints(problem, interval, probe_distance = 2e-8)
#   sd_arp_topology_audit(problem, base_interval, topology_grid, ...)
#
#   sd_invert_test(test_fn, topology_grid, endpoint_tolerance = 1e-8)
#
# Conventions
# -----------
# * scope is "full_path" or "post_connected".
# * procedure is "conditional" or "c_lf".
# * A, betahat, Sigma and target are in the quotient coordinates used by ARP.
# * FLCI betahat is in the original coordinate system (the rows of model$E).
# * A numerical failure never becomes a non-rejection.  It is returned with
#   success = FALSE and accept/reject = NA.
# * Inversion retains every accepted component and reports total Lebesgue
#   length.  An unresolved accepted search boundary is not reported as a
#   finite confidence-set length.

.sd_ie_state <- new.env(parent = emptyenv())
.sd_ie_state$dependencies <- NULL
.sd_ie_state$project_root <- NULL

.sd_ie_stop <- function(...) {
  stop(..., call. = FALSE)
}

.sd_ie_scalar <- function(x, what, finite = TRUE) {
  ok <- length(x) == 1L && !is.na(x)
  if (finite) ok <- ok && is.finite(x)
  if (!ok) .sd_ie_stop(what, " must be one ", if (finite) "finite " else "", "scalar.")
  invisible(x)
}

.sd_ie_capture <- function(fun) {
  warnings <- character()
  error_message <- NULL
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      error_message <<- conditionMessage(condition)
      NULL
    }
  )
  list(value = value, warnings = unique(warnings), error = error_message)
}

.sd_ie_preserve_rng <- function(fun) {
  withr::with_preserve_seed(fun())
}

.sd_ie_with_seed <- function(seed, fun) {
  if (is.null(seed)) return(fun())
  .sd_ie_scalar(seed, "seed")
  withr::with_seed(as.integer(seed), fun())
}

.sd_ie_find_project_root <- function(project_root = NULL) {
  package_path <- getNamespaceInfo(asNamespace("CohortHonestDiD"), "path")
  roots <- file.path(package_path, c("upstream", "inst/upstream"))
  root <- roots[dir.exists(file.path(roots, "R"))][1L]
  if (is.na(root)) .sd_ie_stop("Package numerical resources are unavailable.")
  root
}

sd_inference_load_dependencies <- function(project_root = NULL) {
  if (is.null(project_root) && !is.null(.sd_ie_state$dependencies)) {
    return(invisible(.sd_ie_state$dependencies))
  }
  root <- .sd_ie_find_project_root(project_root)
  if (!is.null(.sd_ie_state$dependencies) &&
      identical(root, .sd_ie_state$project_root)) {
    return(invisible(.sd_ie_state$dependencies))
  }

  dep <- new.env(parent = baseenv())
  dep$setNames <- stats::setNames
  dep$head <- utils::head
  dep$tail <- utils::tail
  sys.source(
    file.path(root, "R", "flci_core.R"),
    envir = dep
  )
  sys.source(
    file.path(root, "R", "arp-nuisance.R"),
    envir = dep
  )

  # Keep the vendored LP implementation unchanged.  Wrap the private copy so
  # an ill-scaled moment system gets one deterministic, mathematically
  # equivalent row-normalized attempt.  Every returned solution is certified
  # in the original units with the same relative gates used below by
  # sd_arp_truth_test().
  original_test_delta_lp_fn <- get(
    ".test_delta_lp_fn", envir = dep, inherits = FALSE
  )
  assign(
    ".sd_ie_original_test_delta_lp_fn",
    original_test_delta_lp_fn,
    envir = dep
  )

  test_delta_certificate <- function(raw, y_T, X_T, sigma) {
    y_T <- suppressWarnings(as.numeric(y_T))
    X_T <- if (is.null(dim(X_T))) {
      matrix(suppressWarnings(as.numeric(X_T)), ncol = 1L)
    } else {
      suppressWarnings(matrix(
        as.numeric(X_T), nrow = nrow(X_T), ncol = ncol(X_T)
      ))
    }
    sigma <- suppressWarnings(as.matrix(sigma))
    error_flag <- if (is.list(raw) && length(raw$error_flag) == 1L) {
      suppressWarnings(as.integer(raw$error_flag))
    } else {
      NA_integer_
    }
    eta <- if (is.list(raw) && length(raw$eta_star) == 1L) {
      suppressWarnings(as.numeric(raw$eta_star))
    } else {
      NA_real_
    }
    delta <- if (is.list(raw)) {
      suppressWarnings(as.numeric(raw$delta_star))
    } else {
      numeric()
    }
    lambda <- if (is.list(raw)) {
      suppressWarnings(as.numeric(raw$lambda))
    } else {
      numeric()
    }
    dimensions_ok <- is.matrix(X_T) && is.matrix(sigma) &&
      length(y_T) == nrow(X_T) &&
      identical(dim(sigma), c(length(y_T), length(y_T)))
    sd_selected <- if (dimensions_ok) {
      suppressWarnings(sqrt(diag(sigma)))
    } else {
      numeric()
    }
    malformed <- !dimensions_ok || any(!is.finite(c(y_T, X_T, sigma))) ||
      length(sd_selected) != length(y_T) || any(!is.finite(sd_selected)) ||
      is.na(error_flag) || error_flag != 0L || !is.finite(eta) ||
      length(delta) != ncol(X_T) || any(!is.finite(delta)) ||
      length(lambda) != length(y_T) || any(!is.finite(lambda))
    if (malformed) {
      return(list(
        passed = FALSE,
        malformed = TRUE,
        error_flag = error_flag,
        eta = eta,
        primal_residual = NA_real_,
        primal_tolerance = NA_real_,
        dual_residual = NA_real_,
        dual_objective_residual = NA_real_,
        dual_nonnegativity_residual = NA_real_,
        dual_tolerance = NA_real_
      ))
    }

    fitted <- as.numeric(X_T %*% delta)
    primal_slack <- y_T - fitted - eta * sd_selected
    primal_residual <- max(0, max(primal_slack))
    primal_scale <- max(
      1, abs(y_T), abs(fitted), abs(eta * sd_selected)
    )
    primal_tolerance <- 1e-8 * primal_scale

    dual_stationarity <- c(
      sum(lambda * sd_selected) - 1,
      as.numeric(t(X_T) %*% lambda)
    )
    eta_lower_bound <- sum(lambda * y_T)
    dual_residual <- max(abs(dual_stationarity))
    dual_objective_residual <- abs(eta_lower_bound - eta)
    dual_nonnegativity_residual <- max(0, -min(lambda))
    dual_scale <- max(
      1, abs(eta), abs(eta_lower_bound), abs(dual_stationarity), abs(lambda)
    )
    dual_tolerance <- 1e-8 * dual_scale
    passed <- is.finite(primal_residual) &&
      primal_residual <= primal_tolerance &&
      is.finite(dual_residual) && dual_residual <= dual_tolerance &&
      is.finite(dual_objective_residual) &&
      dual_objective_residual <= dual_tolerance &&
      is.finite(dual_nonnegativity_residual) &&
      dual_nonnegativity_residual <= dual_tolerance
    list(
      passed = passed,
      malformed = FALSE,
      error_flag = error_flag,
      eta = eta,
      primal_residual = primal_residual,
      primal_tolerance = primal_tolerance,
      dual_residual = dual_residual,
      dual_objective_residual = dual_objective_residual,
      dual_nonnegativity_residual = dual_nonnegativity_residual,
      dual_tolerance = dual_tolerance
    )
  }
  environment(test_delta_certificate) <- dep
  assign(
    ".sd_ie_test_delta_lp_certificate",
    test_delta_certificate,
    envir = dep
  )

  wrapped_test_delta_lp_fn <- function(
      y_T, X_T, sigma,
      primary_timeout_seconds = 2L,
      fallback_timeout_seconds = 10L) {
    call_original <- function(y_value, X_value, sigma_value) {
      tryCatch(
        .sd_ie_original_test_delta_lp_fn(
          y_T = y_value,
          X_T = X_value,
          sigma = sigma_value,
          primary_timeout_seconds = primary_timeout_seconds,
          fallback_timeout_seconds = fallback_timeout_seconds
        ),
        error = function(condition) structure(
          list(message = conditionMessage(condition)),
          class = c("sd_ie_lp_error", "list")
        )
      )
    }
    solver_flag <- function(value) {
      if (is.list(value) && length(value$error_flag) == 1L) {
        suppressWarnings(as.integer(value$error_flag))
      } else {
        NA_integer_
      }
    }

    unscaled <- call_original(y_T, X_T, sigma)
    unscaled_certificate <- .sd_ie_test_delta_lp_certificate(
      unscaled, y_T = y_T, X_T = X_T, sigma = sigma
    )
    if (isTRUE(unscaled_certificate$passed)) {
      unscaled$row_normalized_fallback_attempted <- FALSE
      unscaled$row_normalized_fallback_used <- FALSE
      unscaled$original_unit_certificate <- unscaled_certificate
      return(unscaled)
    }

    y_numeric <- suppressWarnings(as.numeric(y_T))
    X_matrix <- if (is.null(dim(X_T))) {
      matrix(suppressWarnings(as.numeric(X_T)), ncol = 1L)
    } else {
      suppressWarnings(matrix(
        as.numeric(X_T), nrow = nrow(X_T), ncol = ncol(X_T)
      ))
    }
    sigma_matrix <- suppressWarnings(as.matrix(sigma))
    dimensions_ok <- is.matrix(X_matrix) && is.matrix(sigma_matrix) &&
      length(y_numeric) == nrow(X_matrix) &&
      identical(dim(sigma_matrix), c(length(y_numeric), length(y_numeric)))
    sd_selected <- if (dimensions_ok) {
      suppressWarnings(sqrt(diag(sigma_matrix)))
    } else {
      numeric()
    }
    scaling_ok <- dimensions_ok &&
      all(is.finite(c(y_numeric, X_matrix, sigma_matrix, sd_selected))) &&
      all(sd_selected > 0)

    normalized <- structure(
      list(message = "Row normalization was unavailable."),
      class = c("sd_ie_lp_error", "list")
    )
    normalized_original_units <- normalized
    normalized_certificate <- list(
      passed = FALSE, malformed = TRUE, error_flag = NA_integer_,
      eta = NA_real_, primal_residual = NA_real_,
      primal_tolerance = NA_real_, dual_residual = NA_real_,
      dual_objective_residual = NA_real_,
      dual_nonnegativity_residual = NA_real_, dual_tolerance = NA_real_
    )
    if (scaling_ok) {
      y_normalized <- y_numeric / sd_selected
      X_normalized <- X_matrix / sd_selected
      sigma_normalized <- sigma_matrix / outer(sd_selected, sd_selected)
      sigma_normalized <- (sigma_normalized + t(sigma_normalized)) / 2
      normalized <- call_original(
        y_normalized, X_normalized, sigma_normalized
      )
      normalized_original_units <- normalized
      if (is.list(normalized_original_units) &&
          !inherits(normalized_original_units, "sd_ie_lp_error") &&
          length(normalized_original_units$lambda) == length(sd_selected)) {
        normalized_original_units$lambda_standardized <-
          suppressWarnings(as.numeric(normalized_original_units$lambda))
        normalized_original_units$lambda <-
          normalized_original_units$lambda_standardized / sd_selected
      }
      normalized_certificate <- .sd_ie_test_delta_lp_certificate(
        normalized_original_units,
        y_T = y_numeric, X_T = X_matrix, sigma = sigma_matrix
      )
    }

    if (isTRUE(normalized_certificate$passed)) {
      normalized_original_units$row_normalized_fallback_attempted <- TRUE
      normalized_original_units$row_normalized_fallback_used <- TRUE
      normalized_original_units$unscaled_solver_error_flag <-
        solver_flag(unscaled)
      normalized_original_units$unscaled_original_unit_certificate <-
        unscaled_certificate
      normalized_original_units$original_unit_certificate <-
        normalized_certificate
      return(normalized_original_units)
    }

    list(
      eta_star = NA_real_,
      delta_star = rep(NA_real_, ncol(X_matrix)),
      lambda = rep(NA_real_, length(y_numeric)),
      error_flag = 1001L,
      primary_error_flag = solver_flag(unscaled),
      fallback_used = TRUE,
      row_normalized_fallback_attempted = TRUE,
      row_normalized_fallback_used = FALSE,
      unscaled_solver_error_flag = solver_flag(unscaled),
      row_normalized_solver_error_flag = solver_flag(normalized),
      unscaled_original_unit_certificate = unscaled_certificate,
      row_normalized_original_unit_certificate = normalized_certificate,
      row_normalization_available = scaling_ok
    )
  }
  environment(wrapped_test_delta_lp_fn) <- dep
  assign(
    ".test_delta_lp_fn",
    wrapped_test_delta_lp_fn,
    envir = dep
  )

  required <- c(
    "flci_nested_solve", "flci_direct_solve", "flci_model_v",
    "flci_model_sd", "flci_model_bias", "flci_model_objective",
    ".construct_Gamma", ".test_delta_lp_fn", ".lp_conditional_test_fn",
    ".compute_least_favorable_cv"
  )
  missing <- required[!vapply(required, exists, logical(1L),
                             envir = dep, inherits = FALSE)]
  if (length(missing) > 0L) {
    .sd_ie_stop(
      "The sourced project files lack required functions: ",
      paste(missing, collapse = ", "), "."
    )
  }
  .sd_ie_state$dependencies <- dep
  .sd_ie_state$project_root <- root
  invisible(dep)
}

.sd_ie_dependencies <- function(project_root = NULL) {
  if (is.null(project_root) && !is.null(.sd_ie_state$dependencies)) {
    return(.sd_ie_state$dependencies)
  }
  sd_inference_load_dependencies(project_root)
  .sd_ie_state$dependencies
}

.sd_ie_scope <- function(scope) {
  .sd_ie_scalar(scope, "scope", finite = FALSE)
  key <- gsub("[- ]", "_", tolower(as.character(scope)))
  aliases <- c(full = "full_path", full_path = "full_path",
               post = "post_connected", post_connected = "post_connected")
  if (!(key %in% names(aliases))) {
    .sd_ie_stop("scope must be 'full_path' or 'post_connected'.")
  }
  unname(aliases[[key]])
}

.sd_ie_procedure <- function(procedure) {
  .sd_ie_scalar(procedure, "procedure", finite = FALSE)
  key <- gsub("[- ]", "_", tolower(as.character(procedure)))
  aliases <- c(conditional = "conditional", arp = "conditional",
               c_lf = "c_lf", clf = "c_lf", lf = "c_lf")
  if (!(key %in% names(aliases))) {
    .sd_ie_stop("procedure must be 'conditional' or 'c_lf'.")
  }
  unname(aliases[[key]])
}

.sd_ie_align_vector <- function(x, reference_names, what) {
  if (!is.numeric(x) || length(x) != length(reference_names)) {
    .sd_ie_stop(what, " must be numeric with length ",
                length(reference_names), ".")
  }
  if (!is.null(names(x))) {
    if (anyDuplicated(names(x)) || !setequal(names(x), reference_names)) {
      .sd_ie_stop(what, " names do not match the required coordinates.")
    }
    x <- x[reference_names]
  }
  x <- as.numeric(x)
  names(x) <- reference_names
  if (any(!is.finite(x))) .sd_ie_stop(what, " contains nonfinite values.")
  x
}

.sd_ie_align_covariance <- function(Sigma, reference_names, what = "Sigma",
                                    tolerance = 1e-10) {
  Sigma <- as.matrix(Sigma)
  d <- length(reference_names)
  if (!is.numeric(Sigma) || !identical(dim(Sigma), c(d, d)) ||
      any(!is.finite(Sigma))) {
    .sd_ie_stop(what, " must be a finite numeric ", d, " by ", d, " matrix.")
  }
  if (!is.null(rownames(Sigma)) || !is.null(colnames(Sigma))) {
    if (is.null(rownames(Sigma)) || is.null(colnames(Sigma)) ||
        anyDuplicated(rownames(Sigma)) || anyDuplicated(colnames(Sigma)) ||
        !setequal(rownames(Sigma), reference_names) ||
        !setequal(colnames(Sigma), reference_names)) {
      .sd_ie_stop(what, " row/column names do not match the coordinates.")
    }
    Sigma <- Sigma[reference_names, reference_names, drop = FALSE]
  }
  scale <- max(1, max(abs(Sigma)))
  if (max(abs(Sigma - t(Sigma))) > tolerance * scale) {
    .sd_ie_stop(what, " is materially asymmetric.")
  }
  Sigma <- (Sigma + t(Sigma)) / 2
  dimnames(Sigma) <- list(reference_names, reference_names)
  Sigma
}

.sd_ie_failure <- function(status, message, warnings = character(), ...) {
  c(
    list(
      success = FALSE,
      status = status,
      message = message,
      warnings = unique(warnings)
    ),
    list(...)
  )
}

# -----------------------------------------------------------------------------
# FLCI
# -----------------------------------------------------------------------------

sd_flci_prepare <- function(
    model,
    scope = c("full_path", "post_connected"),
    alpha = 0.05,
    fit = NULL,
    solver = c("nested", "direct"),
    solve_args = list(),
    project_root = NULL) {
  scope <- .sd_ie_scope(scope[[1L]])
  solver <- match.arg(solver)
  .sd_ie_scalar(alpha, "alpha")
  if (alpha <= 0 || alpha >= 0.5) .sd_ie_stop("alpha must lie in (0, 0.5).")
  if (!is.list(model) || is.null(model$E) || is.null(model$R) ||
      is.null(model$Sigma_bar)) {
    .sd_ie_stop("model is not a reduced FLCI model.")
  }
  if (!is.list(solve_args) || is.null(names(solve_args)) && length(solve_args)) {
    .sd_ie_stop("solve_args must be a named list.")
  }
  if (any(c("model", "alpha") %in% names(solve_args))) {
    .sd_ie_stop("solve_args must not replace model or alpha.")
  }
  dep <- .sd_ie_dependencies(project_root)

  captured <- .sd_ie_capture(function() {
    if (scope == "post_connected") {
      return(list(
        gamma = numeric(nrow(model$R)),
        solver_name = "analytic",
        solver_status = "analytic",
        solver_ok = TRUE,
        raw_fit = NULL
      ))
    }

    # When M = 0, every admissible bias is exactly zero. The FLCI problem
    # therefore reduces to minimizing v(gamma)' Sigma v(gamma), a quadratic
    # in gamma. If its curvature is singular, select the Moore--Penrose
    # minimum-Euclidean-norm gamma among the variance-minimizing rules. This
    # tie-break is fixed by the design and covariance, not by the realized
    # coefficient estimates used to center the interval.
    is_zero_bias <- model$B_Q == 0 && all(model$m_R == 0)
    if (is_zero_bias) {
      k <- nrow(model$R)
      if (k == 0L) {
        gamma <- numeric()
        stationarity_residual <- 0
        stationarity_residual_scaled <- 0
        null_rhs_residual <- 0
        null_rhs_residual_scaled <- 0
        curvature_rank <- 0L
        curvature_tolerance <- 0
        minimum_curvature_eigenvalue <- NA_real_
        tie_break <- "not_needed"
      } else {
        D_gamma <- rbind(
          t(model$R),
          matrix(0, nrow = length(model$ell), ncol = k)
        )
        v_zero <- c(model$q_0, model$ell)
        # Common covariance units preserve the minimizer, rank decisions,
        # and relative residuals without overflowing uncanceled products.
        # Keep all RHS columns, but only covariance rows on which gamma acts.
        active <- rowSums(abs(D_gamma)) > 0
        D_active <- D_gamma[active, , drop = FALSE]
        block <- dep$flci_scale_covariance_block(model$Sigma_bar[active, , drop = FALSE])
        covariance_scale <- block$scale
        covariance <- block$covariance
        product <- dep$flci_checked_product
        D_Sigma <- product(t(D_active), covariance)
        D_Sigma_abs <- product(t(abs(D_active)), abs(covariance))
        curvature <- dep$flci_symmetrize(product(D_Sigma, D_gamma))
        rhs <- as.numeric(product(D_Sigma, v_zero))
        curvature_eigen <- eigen(curvature, symmetric = TRUE)
        curvature_scale <- max(abs(curvature_eigen$values), norm(
          product(D_Sigma_abs, abs(D_gamma)), type = "2"))
        curvature_keep <- curvature_eigen$values >
          1e-10 * curvature_scale
        curvature_rank <- sum(curvature_keep)
        curvature_tolerance <- 1e-10 * curvature_scale
        minimum_curvature_eigenvalue <- min(curvature_eigen$values)
        vectors_keep <- curvature_eigen$vectors[, curvature_keep, drop = FALSE]
        gamma <- if (curvature_rank) {
          -as.numeric(
            vectors_keep %*%
              ((t(vectors_keep) %*% rhs) /
                 curvature_eigen$values[curvature_keep])
          )
        } else {
          numeric(k)
        }
        names(gamma) <- names(model$m_R)
        vectors_null <- curvature_eigen$vectors[, !curvature_keep, drop = FALSE]
        null_rhs_residual <- if (ncol(vectors_null)) {
          max(abs(as.numeric(t(vectors_null) %*% rhs)))
        } else 0
        null_rhs_scale <- max(abs(product(D_Sigma_abs, abs(v_zero))))
        relative_residual <- function(residual, scale) {
          if (scale > 0) residual / scale else if (residual == 0) 0 else Inf
        }
        null_rhs_residual_scaled <- relative_residual(null_rhs_residual, null_rhs_scale)
        tie_break <- if (curvature_rank < k) {
          "minimum_euclidean_norm_gamma"
        } else {
          "not_needed"
        }
        stationarity <- as.numeric(
          product(D_Sigma, v_zero + D_gamma %*% gamma)
        )
        stationarity_residual <- if (length(stationarity)) {
          max(abs(stationarity))
        } else 0
        stationarity_scale <- max(abs(product(D_Sigma_abs,
          abs(v_zero) + abs(D_gamma) %*% abs(gamma))))
        stationarity_residual_scaled <-
          relative_residual(stationarity_residual, stationarity_scale)
        stationarity_residual <- stationarity_residual * covariance_scale
        null_rhs_residual <- null_rhs_residual * covariance_scale
        curvature_tolerance <- curvature_tolerance * covariance_scale
        minimum_curvature_eigenvalue <- minimum_curvature_eigenvalue * covariance_scale
      }
      solver_ok <- all(is.finite(gamma)) &&
        is.finite(stationarity_residual_scaled) &&
        stationarity_residual_scaled <= 1e-8 &&
        is.finite(null_rhs_residual_scaled) &&
        null_rhs_residual_scaled <= 1e-8
      return(list(
        gamma = gamma,
        solver_name = "analytic_zero_M",
        solver_status = if (solver_ok) {
          if (tie_break == "minimum_euclidean_norm_gamma") {
            "analytic_zero_M_minimum_norm_tiebreak"
          } else {
            "analytic_zero_M"
          }
        } else {
          "analytic_zero_M_stationarity_fail"
        },
        solver_ok = solver_ok,
        raw_fit = list(
          convergence = 0L,
          optimizer_convergence = 0L,
          selected_source = "analytic_zero_M",
          objective_residual = 0,
          kkt_residual = stationarity_residual,
          kkt_residual_scaled = stationarity_residual_scaled,
          kkt_residual_relative = stationarity_residual_scaled,
          null_rhs_residual = null_rhs_residual,
          null_rhs_residual_scaled = null_rhs_residual_scaled,
          pre_polish_kkt_residual = stationarity_residual,
          pre_polish_kkt_residual_relative = stationarity_residual_scaled,
          curvature_rank = curvature_rank,
          curvature_dimension = k,
          curvature_tolerance = curvature_tolerance,
          minimum_curvature_eigenvalue = minimum_curvature_eigenvalue,
          tie_break = tie_break,
          polish_ok = TRUE,
          polish_status = "analytic_zero_M",
          polish_iterations = 0L,
          polish_attempted_count = 0L,
          polish_success_count = 0L,
          split_complementarity_residual = 0
        )
      ))
    }

    # If covariance annihilates every gamma direction exactly, variance is
    # constant in gamma and its nonnegative L1 bias is minimized at zero.
    # Check the exact identity in finite covariance units, without a cutoff.
    if (nrow(model$R) > 0L && any(model$Sigma_bar != 0)) {
      D_gamma <- rbind(t(model$R),
        matrix(0, nrow = length(model$ell), ncol = nrow(model$R)))
      covariance <- model$Sigma_bar / max(abs(model$Sigma_bar))
      represented <- all(covariance[model$Sigma_bar != 0] != 0)
      if (represented && all(covariance %*% D_gamma == 0)) {
        # Small, nonzero direction coefficients can underflow products even
        # when every normalized covariance entry remains representable.
        represented <- all(vapply(seq_len(ncol(D_gamma)), function(j) {
          index <- D_gamma[, j] != 0
          entries <- covariance[, index, drop = FALSE]
          products <- sweep(entries, 2L, D_gamma[index, j], `*`)
          all(products[entries != 0] != 0)
        }, logical(1)))
      } else {
        represented <- FALSE
      }
      if (represented) {
        return(list(
          gamma = numeric(nrow(model$R)),
          solver_name = "analytic_zero_gamma_variance",
          solver_status = "analytic_zero_gamma_variance",
          solver_ok = TRUE,
          raw_fit = list(tie_break = "minimum_euclidean_norm_gamma")
        ))
      }
    }

    # With identically zero covariance, the objective is exactly the bias
    # B_Q + sum(m_R * abs(gamma)). Its minimum-norm minimizer is gamma = 0.
    if (all(model$Sigma_bar == 0)) {
      return(list(
        gamma = numeric(nrow(model$R)),
        solver_name = "analytic_zero_covariance",
        solver_status = "analytic_zero_covariance",
        solver_ok = TRUE,
        raw_fit = list(tie_break = "minimum_euclidean_norm_gamma")
      ))
    }

    raw <- fit
    if (is.null(raw)) {
      # Express the optimization objective in standard-error units. Scaling
      # every bias bound and covariance by common outcome units leaves gamma
      # unchanged and prevents absolute solver stopping rules from treating
      # small, positive noise or curvature as zero.
      solver_scale <- sqrt(max(abs(model$Sigma_bar)))
      solver_model <- model
      solver_model$Sigma_bar <- model$Sigma_bar / solver_scale / solver_scale
      solver_model$B_Q <- model$B_Q / solver_scale
      solver_model$m_R <- model$m_R / solver_scale
      solver_model$m_Q <- model$m_Q / solver_scale
      call_args <- c(list(model = solver_model, alpha = alpha), solve_args)
      raw <- if (solver == "nested") {
        do.call(dep$flci_nested_solve, call_args)
      } else {
        .sd_ie_preserve_rng(function() {
          do.call(dep$flci_direct_solve, call_args)
        })
      }
      raw$optimization_outcome_scale <- solver_scale
    }
    if (is.list(raw) && !is.null(raw[[solver]]) && is.list(raw[[solver]])) {
      raw <- raw[[solver]]
    }
    if (!is.list(raw) || is.null(raw$gamma)) {
      .sd_ie_stop("The supplied FLCI fit has no gamma solution.")
    }
    gamma <- as.numeric(raw$gamma)
    if (length(gamma) != nrow(model$R) || any(!is.finite(gamma))) {
      .sd_ie_stop("The FLCI solver returned a malformed gamma solution.")
    }
    solver_ok <- TRUE
    solver_status <- "provided_unverified"
    if (!is.null(raw$inner_status_acceptable)) {
      solver_ok <- isTRUE(raw$inner_status_acceptable)
      solver_status <- as.character(raw$inner_status)
    } else if (!is.null(raw$convergence)) {
      convergence <- suppressWarnings(as.integer(raw$convergence))
      optimizer_status_ok <- length(convergence) == 1L &&
        (is.na(convergence) || convergence == 0L)
      optimizer_convergence <- NA_integer_
      if (!is.null(raw$optimizer_convergence)) {
        optimizer_convergence <- suppressWarnings(
          as.integer(raw$optimizer_convergence)
        )
        optimizer_status_ok <- optimizer_status_ok &&
          length(optimizer_convergence) == 1L &&
          !is.na(optimizer_convergence) && optimizer_convergence == 0L
      }
      objective_ok <- length(raw$objective_residual) == 1L &&
        is.finite(raw$objective_residual) && raw$objective_residual <= 1e-8
      # Keep the pilot's audited raw KKT gate.  The historical field named
      # `kkt_residual_scaled` is retained for compatibility, while the new
      # coordinatewise relative residual is recorded as a diagnostic rather
      # than silently changing the acceptance contract.
      kkt_ok <- length(raw$kkt_residual) == 1L &&
        is.finite(raw$kkt_residual) && raw$kkt_residual <= 2e-7
      relative_diagnostic_ok <-
        length(raw$kkt_residual_relative) == 1L &&
        is.finite(raw$kkt_residual_relative)
      polish_ok <- is.null(raw$polish_ok) || isTRUE(raw$polish_ok)
      polish_status <- suppressWarnings(as.character(raw$polish_status))
      polish_certified <- isTRUE(raw$polish_ok) &&
        length(polish_status) == 1L && !is.na(polish_status) &&
        polish_status %in% c("already_stationary", "converged")
      certified_after_optimizer_status <-
        !optimizer_status_ok && polish_certified
      # A raw optimizer termination code is advisory after the unchanged
      # candidate passes the deterministic active-orthant polish and the
      # independently recomputed objective and full nonsmooth KKT gates.
      # Missing or failed certificates still fail closed.
      solver_ok <- (optimizer_status_ok ||
                      certified_after_optimizer_status) &&
        objective_ok && kkt_ok &&
        relative_diagnostic_ok && polish_ok
      recovered_status <- if (length(optimizer_convergence) == 1L &&
          !is.na(optimizer_convergence) && optimizer_convergence != 0L) {
        optimizer_convergence
      } else if (length(convergence) == 1L && !is.na(convergence) &&
          convergence != 0L) {
        convergence
      } else {
        NA_integer_
      }
      solver_status <- if (certified_after_optimizer_status &&
          identical(recovered_status, 52L)) {
        "direct_polished_after_convergence_52"
      } else if (certified_after_optimizer_status) {
        paste0(
          "direct_polished_after_optimizer_status_",
          if (is.na(recovered_status)) "missing" else recovered_status
        )
      } else if (is.na(convergence)) {
        paste0("direct_", .sd_ie_or(raw$selected_source, "initial_point"))
      } else {
        paste0("direct_convergence_", convergence)
      }
      if (!objective_ok) solver_status <- paste0(solver_status, "_objective_fail")
      if (!kkt_ok) solver_status <- paste0(solver_status, "_kkt_fail")
      if (!relative_diagnostic_ok) {
        solver_status <- paste0(solver_status, "_relative_kkt_missing")
      }
      if (!polish_ok) solver_status <- paste0(solver_status, "_polish_fail")
    }
    list(
      gamma = gamma,
      solver_name = if (is.null(fit)) solver else "provided",
      solver_status = solver_status,
      solver_ok = solver_ok,
      raw_fit = raw
    )
  })

  if (!is.null(captured$error)) {
    ans <- .sd_ie_failure(
      "flci_solver_error", captured$error, captured$warnings,
      scope = scope, alpha = alpha, model = model
    )
    class(ans) <- c("sd_flci_rule", "list")
    return(ans)
  }
  solved <- captured$value
  if (!isTRUE(solved$solver_ok)) {
    ans <- .sd_ie_failure(
      "flci_solver_failure",
      paste0("FLCI solver status was '", solved$solver_status, "'."),
      captured$warnings,
      scope = scope, alpha = alpha, model = model,
      raw_fit = solved$raw_fit
    )
    class(ans) <- c("sd_flci_rule", "list")
    return(ans)
  }

  evaluated <- .sd_ie_capture(function() {
    gamma <- solved$gamma
    v <- dep$flci_model_v(model, gamma)
    standard_error <- dep$flci_model_sd(model, gamma)
    maximum_bias <- dep$flci_model_bias(model, gamma)
    half_length <- dep$flci_model_objective(model, gamma, alpha = alpha)
    if (any(!is.finite(c(v, standard_error, maximum_bias, half_length))) ||
        standard_error < 0 || maximum_bias < 0 || half_length < 0) {
      .sd_ie_stop("FLCI rule contains a nonfinite or negative quantity.")
    }
    list(
      gamma = gamma,
      v = as.numeric(v),
      standard_error = standard_error,
      maximum_bias = maximum_bias,
      half_length = half_length,
      folded_critical_value = if (standard_error > 0) {
        half_length / standard_error
      } else {
        NA_real_
      }
    )
  })
  all_warnings <- unique(c(captured$warnings, evaluated$warnings))
  if (!is.null(evaluated$error)) {
    ans <- .sd_ie_failure(
      "flci_rule_error", evaluated$error, all_warnings,
      scope = scope, alpha = alpha, model = model,
      raw_fit = solved$raw_fit
    )
    class(ans) <- c("sd_flci_rule", "list")
    return(ans)
  }

  ans <- c(
    list(
      success = TRUE,
      status = if (length(all_warnings)) "ok_with_warning" else "ok",
      message = NULL,
      warnings = all_warnings,
      scope = scope,
      alpha = alpha,
      model = model,
      solver = solved$solver_name,
      solver_status = solved$solver_status,
      raw_fit = solved$raw_fit
    ),
    evaluated$value
  )
  class(ans) <- c("sd_flci_rule", "list")
  ans
}

# Local null-coalescing helper; avoids an external dependency and does not
# install a commonly named infix operator in the caller's environment.
.sd_ie_or <- function(x, y) if (is.null(x)) y else x

sd_flci_interval <- function(rule, betahat, theta = NULL) {
  if (!inherits(rule, "sd_flci_rule")) {
    .sd_ie_stop("rule must be returned by sd_flci_prepare().")
  }
  if (!isTRUE(rule$success)) {
    return(.sd_ie_failure(
      "flci_rule_unavailable", rule$message, rule$warnings,
      scope = rule$scope, theta = theta,
      accept = NA, reject = NA
    ))
  }
  model <- rule$model
  original_names <- rownames(model$E)
  if (is.null(original_names)) {
    original_names <- paste0("coordinate_", seq_len(nrow(model$E)))
  }
  beta <- .sd_ie_align_vector(betahat, original_names, "betahat")
  y_bar <- as.numeric(t(model$E) %*% beta)
  center <- sum(rule$v * y_bar)
  lower <- center - rule$half_length
  upper <- center + rule$half_length
  if (any(!is.finite(c(center, lower, upper)))) {
    return(.sd_ie_failure(
      "flci_evaluation_failure", "FLCI endpoints are nonfinite.",
      rule$warnings, scope = rule$scope, theta = theta,
      accept = NA, reject = NA
    ))
  }
  if (!is.null(theta)) {
    .sd_ie_scalar(theta, "theta")
    accept <- theta >= lower && theta <= upper
    reject <- !accept
  } else {
    accept <- reject <- NA
  }
  components <- data.frame(
    lower = lower,
    upper = upper,
    lower_search_boundary = FALSE,
    upper_search_boundary = FALSE
  )
  list(
    success = TRUE,
    status = rule$status,
    message = NULL,
    warnings = rule$warnings,
    scope = rule$scope,
    procedure = "flci",
    theta = theta,
    accept = accept,
    reject = reject,
    center = center,
    lower = lower,
    upper = upper,
    half_length = rule$half_length,
    total_length = 2 * rule$half_length,
    n_components = 1L,
    components = components,
    # Analytic FLCI endpoints do not require a search-grid topology audit.
    # Mark that distinction explicitly for the common confidence-set schema.
    endpoint_certified = TRUE,
    boundary_accepted = FALSE,
    expansions = 0L,
    bisection_evaluations = 0L,
    total_test_evaluations = 0L,
    standard_error = rule$standard_error,
    maximum_bias = rule$maximum_bias,
    folded_critical_value = rule$folded_critical_value,
    gamma = rule$gamma,
    v = rule$v
  )
}

sd_flci_truth_test <- function(rule, betahat, theta) {
  sd_flci_interval(rule = rule, betahat = betahat, theta = theta)
}

# -----------------------------------------------------------------------------
# ARP Conditional and C-LF geometry
# -----------------------------------------------------------------------------

# To construct A and d from sdmc_build_model(design, M, Sigma_bar):
#
#   built <- sdmc_build_model(design, M, Sigma_bar)
#   W_bar <- built$model$W_bar
#   dimnames(W_bar) <- list(design$reduced_names, design$reduced_names)
#   A <- t(solve(t(W_bar), t(built$A_Delta %*% design$E)))
#   dimnames(A) <- list(rownames(built$A_Delta), design$reduced_names)
#   d <- built$d_vec
#
# This is A_Delta E W_bar^{-1}.  It maps quotient coefficient estimates to
# the signed SD moments without passing the singular 21-coordinate covariance
# to the ARP code.

sd_arp_geometry <- function(
    A,
    d,
    target,
    scope = c("full_path", "post_connected"),
    project_root = NULL,
    tolerance = 1e-12) {
  scope <- .sd_ie_scope(scope[[1L]])
  .sd_ie_scalar(tolerance, "tolerance")
  if (tolerance <= 0) .sd_ie_stop("tolerance must be positive.")
  dep <- .sd_ie_dependencies(project_root)

  A <- as.matrix(A)
  if (!is.numeric(A) || nrow(A) == 0L || ncol(A) < 2L ||
      any(!is.finite(A))) {
    .sd_ie_stop("A must be a nonempty finite numeric matrix.")
  }
  if (is.null(colnames(A)) || anyDuplicated(colnames(A))) {
    .sd_ie_stop("A must have unique coordinate column names.")
  }
  if (is.null(rownames(A))) rownames(A) <- paste0("moment_", seq_len(nrow(A)))
  if (anyDuplicated(rownames(A))) .sd_ie_stop("A row names must be unique.")
  coordinate_names <- colnames(A)
  d <- .sd_ie_align_vector(d, rownames(A), "d")
  target <- .sd_ie_align_vector(target, coordinate_names, "target")

  pre <- grepl("pre", coordinate_names)
  post <- grepl("post", coordinate_names)
  if (!all(pre | post) || !any(pre) || !any(post)) {
    .sd_ie_stop("A columns must be classifiable as pre or post coordinates.")
  }
  if (!identical(c(which(pre), which(post)), seq_len(ncol(A)))) {
    .sd_ie_stop("A coordinates must be ordered pre first, then post.")
  }
  if (max(abs(target[pre])) > tolerance || max(abs(target[post])) <= tolerance) {
    .sd_ie_stop("target must be zero on pre coordinates and nonzero on post.")
  }

  Gamma <- dep$.construct_Gamma(target[post])
  gamma_condition <- kappa(Gamma, exact = TRUE)
  if (!is.finite(gamma_condition) || gamma_condition > 1 / tolerance) {
    .sd_ie_stop("The target change-of-basis matrix is ill-conditioned.")
  }
  A_gamma_inverse <- A[, post, drop = FALSE] %*% solve(Gamma)
  a_one <- as.numeric(A_gamma_inverse[, 1L])
  a_minus <- A_gamma_inverse[, -1L, drop = FALSE]
  rows_for_arp <- if (scope == "post_connected") {
    which(rowSums(abs(A[, post, drop = FALSE]) > tolerance) > 0L)
  } else {
    seq_len(nrow(A))
  }
  if (length(rows_for_arp) == 0L) {
    .sd_ie_stop("The selected scope contains no ARP moments.")
  }

  ans <- list(
    A = A,
    d = d,
    target = target,
    scope = scope,
    coordinate_names = coordinate_names,
    pre_columns = which(pre),
    post_columns = which(post),
    rows_for_arp = rows_for_arp,
    moment_names = rownames(A),
    selected_moment_names = rownames(A)[rows_for_arp],
    Gamma = Gamma,
    Gamma_condition = gamma_condition,
    a_one = a_one,
    a_minus = a_minus,
    tolerance = tolerance
  )
  class(ans) <- c("sd_arp_geometry", "list")
  ans
}

.sd_ie_sigma_y <- function(geometry, Sigma) {
  Sigma <- .sd_ie_align_covariance(
    Sigma, geometry$coordinate_names, "Sigma", geometry$tolerance
  )
  sigma_y <- geometry$A %*% Sigma %*% t(geometry$A)
  sigma_y <- (sigma_y + t(sigma_y)) / 2
  dimnames(sigma_y) <- list(geometry$moment_names, geometry$moment_names)
  selected_variances <- diag(sigma_y)[geometry$rows_for_arp]
  variance_scale <- max(abs(diag(sigma_y)))
  if (any(!is.finite(selected_variances)))
    .sd_ie_stop("Selected ARP moments have invalid variance.")
  if (any(selected_variances <= geometry$tolerance * variance_scale))
    stop(structure(list(message = paste0(
      "Conditional and C-LF inference is unsupported for selected moments with zero or numerically degenerate variance; ",
      "no confidence set is reported."), call = NULL),
      class = c("sd_unsupported_moment_variance", "error", "condition")))
  list(Sigma = Sigma, sigma_y = sigma_y)
}

# Simulate the C-LF first-stage critical value.  The default diagnostic engine
# calls the project's LP primitive for every Gaussian draw and records every
# nonzero solver flag.  The optional project engine reproduces the original
# helper exactly, but that helper does not reveal its internal LP failure count.
# With seed = NULL, the diagnostic engine consumes the currently registered
# RNG stream (the runner uses independent L'Ecuyer streams).  A numeric seed is
# local and restores the caller's RNG state after the critical value is built.
sd_lf_critical_value <- function(
    geometry,
    Sigma,
    hybrid_kappa,
    sims = 10000L,
    seed = 0L,
    max_failure_fraction = 0,
    engine = c("diagnostic", "project"),
    project_root = NULL,
    keep_eta = FALSE) {
  if (!inherits(geometry, "sd_arp_geometry")) {
    .sd_ie_stop("geometry must be returned by sd_arp_geometry().")
  }
  engine <- match.arg(engine)
  .sd_ie_scalar(hybrid_kappa, "hybrid_kappa")
  if (hybrid_kappa <= 0 || hybrid_kappa >= 1) {
    .sd_ie_stop("hybrid_kappa must lie in (0, 1).")
  }
  .sd_ie_scalar(sims, "sims")
  sims <- as.integer(sims)
  if (sims < 2L) .sd_ie_stop("sims must be at least two.")
  if (!is.null(seed)) .sd_ie_scalar(seed, "seed")
  seed_label <- if (is.null(seed)) {
    "current_registered_stream"
  } else {
    as.character(as.integer(seed))
  }
  .sd_ie_scalar(max_failure_fraction, "max_failure_fraction")
  if (max_failure_fraction < 0 || max_failure_fraction >= 1) {
    .sd_ie_stop("max_failure_fraction must lie in [0, 1).")
  }
  if (!is.logical(keep_eta) || length(keep_eta) != 1L || is.na(keep_eta)) {
    .sd_ie_stop("keep_eta must be TRUE or FALSE.")
  }
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    .sd_ie_stop("Package 'mvtnorm' is required for C-LF critical values.")
  }
  dep <- .sd_ie_dependencies(project_root)
  covariance <- .sd_ie_sigma_y(geometry, Sigma)
  sigma_y <- covariance$sigma_y
  rows <- geometry$rows_for_arp
  sigma_selected <- sigma_y[rows, rows, drop = FALSE]
  X_selected <- geometry$a_minus[rows, , drop = FALSE]

  if (engine == "project") {
    if (is.null(seed)) {
      .sd_ie_stop(
        "seed = NULL requires engine = 'diagnostic': the original project ",
        "helper unconditionally resets the global RNG."
      )
    }
    captured <- .sd_ie_capture(function() {
      .sd_ie_preserve_rng(function() {
        dep$.compute_least_favorable_cv(
          X_T = geometry$a_minus,
          sigma = sigma_y,
          hybrid_kappa = hybrid_kappa,
          sims = sims,
          rowsForARP = rows,
          seed = seed
        )
      })
    })
    if (!is.null(captured$error) || length(captured$value) != 1L ||
        !is.finite(captured$value)) {
      ans <- .sd_ie_failure(
        "lf_critical_value_failure",
        .sd_ie_or(captured$error, "Project LF helper returned a nonfinite value."),
        captured$warnings,
        critical_value = NA_real_, scope = geometry$scope,
        sims = sims, seed = seed, seed_label = seed_label,
        n_failures = NA_integer_,
        failure_fraction = NA_real_, engine = engine,
        rows_for_arp = rows, sigma_y = sigma_y
      )
    } else {
      ans <- list(
        success = TRUE,
        status = if (length(captured$warnings)) "ok_with_warning" else "ok",
        message = paste(
          "The project engine does not expose internal LF LP failure counts."
        ),
        warnings = captured$warnings,
        critical_value = as.numeric(captured$value),
        scope = geometry$scope,
        hybrid_kappa = hybrid_kappa,
        sims = sims,
        seed = seed,
        seed_label = seed_label,
        n_failures = NA_integer_,
        failure_fraction = NA_real_,
        engine = engine,
        rows_for_arp = rows,
        sigma_y = sigma_y,
        eta = NULL
      )
    }
    class(ans) <- c("sd_lf_critical_value", "list")
    return(ans)
  }

  captured <- .sd_ie_capture(function() {
    draws <- if (is.null(seed)) {
      # The Monte Carlo driver registers independent L'Ecuyer streams.  NULL
      # deliberately consumes (and advances) the currently installed stream.
      mvtnorm::rmvnorm(n = sims, sigma = sigma_selected)
    } else {
      # An explicit seed is a reproducible testing convenience and is fully
      # local: the caller's RNG state is restored on exit.
      .sd_ie_with_seed(seed, function() {
        mvtnorm::rmvnorm(n = sims, sigma = sigma_selected)
      })
    }
    eta <- rep(NA_real_, sims)
    errors <- rep(NA_character_, sims)
    solver_flags <- rep(NA_integer_, sims)

    if (ncol(X_selected) == 0L) {
      eta <- apply(
        sweep(draws, 2L, sqrt(diag(sigma_selected)), "/"),
        1L,
        max
      )
      solver_flags[] <- 0L
    } else {
      for (draw_index in seq_len(sims)) {
        one <- tryCatch(
          dep$.test_delta_lp_fn(
            y_T = as.numeric(draws[draw_index, ]),
            X_T = X_selected,
            sigma = sigma_selected
          ),
          error = function(condition) condition
        )
        if (inherits(one, "condition")) {
          errors[draw_index] <- conditionMessage(one)
          solver_flags[draw_index] <- NA_integer_
        } else {
          solver_flags[draw_index] <- as.integer(one$error_flag)
          if (length(solver_flags[draw_index]) != 1L ||
              is.na(solver_flags[draw_index]) || solver_flags[draw_index] != 0L) {
            errors[draw_index] <- paste0(
              "lpSolveAPI flag ", .sd_ie_or(solver_flags[draw_index], "NA")
            )
          } else if (length(one$eta_star) != 1L || !is.finite(one$eta_star)) {
            errors[draw_index] <- "nonfinite eta"
          } else {
            eta[draw_index] <- as.numeric(one$eta_star)
          }
        }
      }
    }
    list(eta = eta, errors = errors, solver_flags = solver_flags)
  })
  if (!is.null(captured$error)) {
    ans <- .sd_ie_failure(
      "lf_critical_value_error", captured$error, captured$warnings,
      critical_value = NA_real_, scope = geometry$scope,
      sims = sims, seed = seed, seed_label = seed_label,
      n_failures = sims,
      failure_fraction = 1, engine = engine,
      rows_for_arp = rows, sigma_y = sigma_y
    )
    class(ans) <- c("sd_lf_critical_value", "list")
    return(ans)
  }

  eta <- captured$value$eta
  good <- is.finite(eta)
  n_failures <- sum(!good)
  failure_fraction <- n_failures / sims
  critical_value <- if (any(good)) {
    as.numeric(stats::quantile(
      eta[good], probs = 1 - hybrid_kappa, names = FALSE, type = 7
    ))
  } else {
    NA_real_
  }
  success <- is.finite(critical_value) &&
    failure_fraction <= max_failure_fraction
  error_summary <- unique(stats::na.omit(captured$value$errors))
  message <- if (success) {
    NULL
  } else {
    paste0(
      "LF critical-value simulation had ", n_failures, " failed LP draw(s) ",
      "out of ", sims, "."
    )
  }
  ans <- list(
    success = success,
    status = if (success) {
      if (length(captured$warnings)) "ok_with_warning" else "ok"
    } else {
      "lf_solver_failure"
    },
    message = message,
    warnings = captured$warnings,
    critical_value = critical_value,
    scope = geometry$scope,
    hybrid_kappa = hybrid_kappa,
    sims = sims,
    seed = seed,
    seed_label = seed_label,
    n_failures = n_failures,
    failure_fraction = failure_fraction,
    failure_messages = head(error_summary, 20L),
    engine = engine,
    rows_for_arp = rows,
    sigma_y = sigma_y,
    eta = if (keep_eta) eta else NULL
  )
  class(ans) <- c("sd_lf_critical_value", "list")
  ans
}

sd_bind_arp_problem <- function(
    geometry,
    betahat,
    Sigma,
    procedure = c("conditional", "c_lf"),
    alpha = 0.05,
    hybrid_kappa = alpha / 10,
    lf_cv = NULL,
    lf_sims = 10000L,
    lf_seed = 0L,
    lf_max_failure_fraction = 0,
    lf_engine = c("diagnostic", "project"),
    project_root = NULL) {
  if (!inherits(geometry, "sd_arp_geometry")) {
    .sd_ie_stop("geometry must be returned by sd_arp_geometry().")
  }
  procedure <- .sd_ie_procedure(procedure[[1L]])
  lf_engine <- match.arg(lf_engine)
  .sd_ie_scalar(alpha, "alpha")
  .sd_ie_scalar(hybrid_kappa, "hybrid_kappa")
  if (alpha <= 0 || alpha >= 0.5 ||
      (procedure == "c_lf" && (hybrid_kappa <= 0 || hybrid_kappa >= alpha))) {
    .sd_ie_stop("Require 0 < hybrid_kappa < alpha < 0.5.")
  }
  beta <- .sd_ie_align_vector(
    betahat, geometry$coordinate_names, "betahat"
  )
  covariance <- .sd_ie_sigma_y(geometry, Sigma)
  sigma_y <- covariance$sigma_y
  base_y <- as.numeric(geometry$A %*% beta - geometry$d)
  names(base_y) <- geometry$moment_names
  warnings <- character()

  critical <- NULL
  if (procedure == "c_lf") {
    if (is.null(lf_cv)) {
      critical <- sd_lf_critical_value(
        geometry = geometry,
        Sigma = covariance$Sigma,
        hybrid_kappa = hybrid_kappa,
        sims = lf_sims,
        seed = lf_seed,
        max_failure_fraction = lf_max_failure_fraction,
        engine = lf_engine,
        project_root = project_root
      )
      if (!isTRUE(critical$success)) {
        ans <- .sd_ie_failure(
          "lf_critical_value_unavailable", critical$message,
          critical$warnings,
          geometry = geometry, procedure = procedure, alpha = alpha,
          hybrid_kappa = hybrid_kappa, lf_critical = critical
        )
        class(ans) <- c("sd_arp_problem", "list")
        return(ans)
      }
      lf_value <- critical$critical_value
      warnings <- critical$warnings
    } else if (inherits(lf_cv, "sd_lf_critical_value")) {
      if (!isTRUE(lf_cv$success)) {
        ans <- .sd_ie_failure(
          "lf_critical_value_unavailable", lf_cv$message, lf_cv$warnings,
          geometry = geometry, procedure = procedure, alpha = alpha,
          hybrid_kappa = hybrid_kappa, lf_critical = lf_cv
        )
        class(ans) <- c("sd_arp_problem", "list")
        return(ans)
      }
      if (!identical(lf_cv$scope, geometry$scope) ||
          !identical(as.integer(lf_cv$rows_for_arp),
                     as.integer(geometry$rows_for_arp)) ||
          abs(lf_cv$hybrid_kappa - hybrid_kappa) > 1e-14) {
        .sd_ie_stop("Precomputed lf_cv does not match scope or hybrid_kappa.")
      }
      scale <- max(1, max(abs(sigma_y)))
      if (!identical(dim(lf_cv$sigma_y), dim(sigma_y)) ||
          max(abs(lf_cv$sigma_y - sigma_y)) > 1e-9 * scale) {
        .sd_ie_stop("Precomputed lf_cv does not match the supplied Sigma.")
      }
      lf_value <- lf_cv$critical_value
      critical <- lf_cv
      warnings <- lf_cv$warnings
    } else {
      .sd_ie_scalar(lf_cv, "lf_cv")
      lf_value <- as.numeric(lf_cv)
      critical <- list(
        success = TRUE, critical_value = lf_value,
        status = "user_supplied_numeric", diagnostics_verified = FALSE
      )
    }
    if (!is.finite(lf_value)) .sd_ie_stop("lf_cv must be finite.")
  } else {
    lf_value <- NULL
  }

  hybrid_flag <- if (procedure == "conditional") "ARP" else "LF"
  hybrid_list <- list(hybrid_kappa = hybrid_kappa)
  if (procedure == "c_lf") hybrid_list$lf_cv <- lf_value
  ans <- list(
    success = TRUE,
    status = if (length(warnings)) "ok_with_warning" else "ok",
    message = NULL,
    warnings = warnings,
    geometry = geometry,
    procedure = procedure,
    alpha = alpha,
    hybrid_kappa = hybrid_kappa,
    hybrid_flag = hybrid_flag,
    hybrid_list = hybrid_list,
    betahat = beta,
    Sigma = covariance$Sigma,
    base_y = base_y,
    sigma_y = sigma_y,
    lf_critical = critical
  )
  class(ans) <- c("sd_arp_problem", "list")
  ans
}

# Find a data-dependent location around which to search for the inverted ARP
# confidence set.  For a fixed theta, the ARP statistic profiles out the
# nuisance vector delta.  Here theta is profiled out as well:
#
#   min_{eta, theta, delta} eta
#     s.t. base_y - a_one theta - a_minus delta <= eta * sd.
#
# This is the same LP primitive used by the ARP test, with a_one appended to
# the nuisance design.  The resulting theta is not the true parameter and
# uses no DGP information.  It is a nuisance-adjusted, sample-dependent anchor
# that remains meaningful when unrestricted linear paths make the raw target
# estimate far from the ARP confidence set.
sd_arp_profile_anchor <- function(problem) {
  if (!inherits(problem, "sd_arp_problem")) {
    .sd_ie_stop("problem must be returned by sd_bind_arp_problem().")
  }
  if (!isTRUE(problem$success)) {
    return(.sd_ie_failure(
      "arp_problem_unavailable", problem$message, problem$warnings,
      theta = NA_real_, eta = NA_real_, eta_lower_bound = NA_real_
    ))
  }

  dep <- .sd_ie_dependencies()
  geometry <- problem$geometry
  rows <- geometry$rows_for_arp
  y_selected <- as.numeric(problem$base_y[rows])
  X_selected <- cbind(
    theta = geometry$a_one,
    geometry$a_minus
  )[rows, , drop = FALSE]
  sigma_selected <- problem$sigma_y[rows, rows, drop = FALSE]
  sd_selected <- sqrt(diag(sigma_selected))

  captured <- .sd_ie_capture(function() {
    dep$.test_delta_lp_fn(
      y_T = y_selected,
      X_T = X_selected,
      sigma = sigma_selected
    )
  })
  all_warnings <- unique(c(problem$warnings, captured$warnings))
  if (!is.null(captured$error)) {
    return(.sd_ie_failure(
      "arp_profile_anchor_error", captured$error, all_warnings,
      theta = NA_real_, eta = NA_real_, eta_lower_bound = NA_real_
    ))
  }

  raw <- captured$value
  error_flag <- if (is.list(raw) && length(raw$error_flag) == 1L) {
    suppressWarnings(as.integer(raw$error_flag))
  } else NA_integer_
  eta <- if (is.list(raw) && length(raw$eta_star) == 1L) {
    suppressWarnings(as.numeric(raw$eta_star))
  } else NA_real_
  coefficients <- if (is.list(raw)) {
    suppressWarnings(as.numeric(raw$delta_star))
  } else numeric()
  lambda <- if (is.list(raw)) {
    suppressWarnings(as.numeric(raw$lambda))
  } else numeric()
  malformed <- is.na(error_flag) || error_flag != 0L || !is.finite(eta) ||
    length(coefficients) != ncol(X_selected) ||
    any(!is.finite(coefficients)) ||
    length(lambda) != length(y_selected) || any(!is.finite(lambda))
  if (malformed) {
    return(.sd_ie_failure(
      "arp_profile_anchor_solver_failure",
      "The joint theta/nuisance profile LP returned a malformed solution.",
      all_warnings,
      theta = NA_real_, eta = eta, eta_lower_bound = NA_real_,
      error_flag = error_flag, raw = raw
    ))
  }

  fitted <- as.numeric(X_selected %*% coefficients)
  primal_slack <- y_selected - fitted - eta * sd_selected
  primal_residual <- max(0, max(primal_slack))
  primal_scale <- max(
    1, abs(y_selected), abs(fitted), abs(eta * sd_selected)
  )
  primal_tolerance <- 1e-8 * primal_scale

  dual_stationarity <- c(
    sum(lambda * sd_selected) - 1,
    as.numeric(t(X_selected) %*% lambda)
  )
  eta_lower_bound <- sum(lambda * y_selected)
  dual_residual <- max(abs(dual_stationarity))
  dual_objective_residual <- abs(eta_lower_bound - eta)
  dual_nonnegativity_residual <- max(0, -min(lambda))
  dual_scale <- max(
    1, abs(eta), abs(eta_lower_bound), abs(dual_stationarity), abs(lambda)
  )
  dual_tolerance <- 1e-8 * dual_scale
  certified <- is.finite(primal_residual) &&
    primal_residual <= primal_tolerance &&
    is.finite(dual_residual) && dual_residual <= dual_tolerance &&
    is.finite(dual_objective_residual) &&
    dual_objective_residual <= dual_tolerance &&
    is.finite(dual_nonnegativity_residual) &&
    dual_nonnegativity_residual <= dual_tolerance
  if (!certified) {
    return(.sd_ie_failure(
      "arp_profile_anchor_certificate_failure",
      "The joint theta/nuisance profile LP failed its primal/dual certificate.",
      all_warnings,
      theta = coefficients[1L], eta = eta,
      eta_lower_bound = eta_lower_bound,
      primal_residual = primal_residual,
      primal_tolerance = primal_tolerance,
      dual_residual = dual_residual,
      dual_objective_residual = dual_objective_residual,
      dual_nonnegativity_residual = dual_nonnegativity_residual,
      dual_tolerance = dual_tolerance,
      error_flag = error_flag, raw = raw
    ))
  }

  list(
    success = TRUE,
    status = if (length(all_warnings)) "ok_with_warning" else "ok",
    message = NULL,
    warnings = all_warnings,
    theta = unname(coefficients[1L]),
    nuisance = unname(coefficients[-1L]),
    eta = eta,
    eta_lower_bound = eta_lower_bound,
    primal_residual = primal_residual,
    primal_tolerance = primal_tolerance,
    dual_residual = dual_residual,
    dual_objective_residual = dual_objective_residual,
    dual_nonnegativity_residual = dual_nonnegativity_residual,
    dual_tolerance = dual_tolerance,
    error_flag = error_flag,
    scope = geometry$scope
  )
}

sd_arp_truth_test <- function(problem, theta) {
  if (!inherits(problem, "sd_arp_problem")) {
    .sd_ie_stop("problem must be returned by sd_bind_arp_problem().")
  }
  .sd_ie_scalar(theta, "theta")
  if (!isTRUE(problem$success)) {
    return(.sd_ie_failure(
      "arp_problem_unavailable", problem$message, problem$warnings,
      theta = theta, accept = NA, reject = NA,
      procedure = .sd_ie_or(problem$procedure, NA_character_),
      scope = .sd_ie_or(problem$geometry$scope, NA_character_)
    ))
  }
  dep <- .sd_ie_dependencies()
  geometry <- problem$geometry
  captured <- .sd_ie_capture(function() {
    dep$.lp_conditional_test_fn(
      theta = theta,
      y_T = problem$base_y - geometry$a_one * theta,
      X_T = geometry$a_minus,
      sigma = problem$sigma_y,
      alpha = problem$alpha,
      hybrid_flag = problem$hybrid_flag,
      hybrid_list = problem$hybrid_list,
      rowsForARP = geometry$rows_for_arp
    )
  })
  all_warnings <- unique(c(problem$warnings, captured$warnings))
  if (!is.null(captured$error)) {
    return(.sd_ie_failure(
      "arp_test_error", captured$error, all_warnings,
      theta = theta, accept = NA, reject = NA,
      procedure = problem$procedure, scope = geometry$scope
    ))
  }
  solver_warning <- any(grepl(
    "LP for eta did not converge properly|solver.*fail|did not converge",
    captured$warnings,
    ignore.case = TRUE
  ))
  raw <- captured$value
  reject <- if (is.list(raw) && length(raw$reject) == 1L) {
    suppressWarnings(as.integer(raw$reject))
  } else {
    NA_integer_
  }
  eta <- if (is.list(raw) && length(raw$eta) == 1L) {
    suppressWarnings(as.numeric(raw$eta))
  } else {
    NA_real_
  }
  invalid <- solver_warning || is.na(reject) || !(reject %in% c(0L, 1L)) ||
    !is.finite(eta)
  if (invalid) {
    return(.sd_ie_failure(
      "arp_solver_failure",
      if (solver_warning) {
        "The ARP LP reported nonconvergence; its fallback non-rejection was discarded."
      } else {
        "The ARP test returned a malformed rejection indicator or eta value."
      },
      all_warnings,
      theta = theta, accept = NA, reject = NA,
      procedure = problem$procedure, scope = geometry$scope,
      eta = eta, raw = raw
    ))
  }

  # The upstream ARP routine does not return the LP status flag on every
  # branch.  It does return the primal solution (and, except for a valid C-LF
  # first-stage rejection, the dual multipliers), so independently certify
  # the LP feasibility and primal/dual objective identities here.  This also
  # prevents a numerically malformed conservative fallback from being counted
  # as a successful non-rejection.
  rows <- geometry$rows_for_arp
  y_selected <- as.numeric(
    (problem$base_y - geometry$a_one * theta)[rows]
  )
  X_selected <- geometry$a_minus[rows, , drop = FALSE]
  sigma_selected <- problem$sigma_y[rows, rows, drop = FALSE]
  sd_selected <- sqrt(diag(sigma_selected))
  delta <- if (is.list(raw)) suppressWarnings(as.numeric(raw$delta)) else NA_real_
  delta_ok <- length(delta) == ncol(X_selected) && all(is.finite(delta))
  if (!delta_ok) {
    return(.sd_ie_failure(
      "arp_primal_malformed",
      "The ARP LP returned a malformed primal nuisance solution.",
      all_warnings,
      theta = theta, accept = NA, reject = NA,
      procedure = problem$procedure, scope = geometry$scope,
      eta = eta, raw = raw
    ))
  }
  fitted <- as.numeric(X_selected %*% delta)
  primal_slack <- y_selected - fitted - eta * sd_selected
  primal_residual <- max(0, max(primal_slack))
  primal_scale <- max(
    1, abs(y_selected), abs(fitted), abs(eta * sd_selected)
  )
  primal_tolerance <- 1e-8 * primal_scale

  lambda_present <- is.list(raw) && !is.null(raw$lambda)
  lambda <- if (lambda_present) {
    suppressWarnings(as.numeric(raw$lambda))
  } else numeric()
  dual_residual <- NA_real_
  dual_objective_residual <- NA_real_
  dual_nonnegativity_residual <- NA_real_
  dual_tolerance <- NA_real_
  conditioning_branch <- "clf_first_stage"
  dual_ok <- TRUE
  if (lambda_present) {
    dual_ok <- length(lambda) == length(y_selected) && all(is.finite(lambda))
    if (dual_ok) {
      dual_stationarity <- c(
        sum(lambda * sd_selected) - 1,
        as.numeric(t(X_selected) %*% lambda)
      )
      dual_residual <- max(abs(dual_stationarity))
      dual_objective_residual <- abs(sum(lambda * y_selected) - eta)
      dual_nonnegativity_residual <- max(0, -min(lambda))
      dual_scale <- max(
        1, abs(eta), abs(sum(lambda * y_selected)),
        abs(dual_stationarity), abs(lambda)
      )
      dual_tolerance <- 1e-8 * dual_scale

      binding <- lambda > 1e-6
      X_binding <- X_selected[binding, , drop = FALSE]
      full_rank <- nrow(X_binding) > 0L &&
        min(dim(X_binding)) > 0L &&
        as.integer(Matrix::rankMatrix(X_binding)) == min(dim(X_binding))
      degenerate <- sum(binding) != ncol(X_selected) + 1L
      conditioning_branch <- if (!full_rank || degenerate) "dual" else "primal"
      dual_ok <- dual_residual <= dual_tolerance &&
        dual_objective_residual <= dual_tolerance &&
        dual_nonnegativity_residual <= dual_tolerance
    }
  } else {
    # The only intended no-lambda branch is the C-LF first-stage rejection.
    dual_ok <- problem$procedure == "c_lf" && reject == 1L &&
      is.finite(problem$hybrid_list$lf_cv) &&
      eta > problem$hybrid_list$lf_cv
  }
  if (!is.finite(primal_residual) || primal_residual > primal_tolerance ||
      !isTRUE(dual_ok)) {
    return(.sd_ie_failure(
      "arp_lp_certificate_failure",
      "The ARP LP failed its independently recomputed primal/dual certificate.",
      all_warnings,
      theta = theta, accept = NA, reject = NA,
      procedure = problem$procedure, scope = geometry$scope,
      eta = eta, raw = raw,
      primal_residual = primal_residual,
      primal_tolerance = primal_tolerance,
      dual_residual = dual_residual,
      dual_objective_residual = dual_objective_residual,
      dual_nonnegativity_residual = dual_nonnegativity_residual,
      dual_tolerance = dual_tolerance,
      conditioning_branch = conditioning_branch
    ))
  }
  list(
    success = TRUE,
    status = if (length(all_warnings)) "ok_with_warning" else "ok",
    message = NULL,
    warnings = all_warnings,
    theta = theta,
    accept = reject == 0L,
    reject = reject == 1L,
    reject_integer = reject,
    eta = eta,
    delta = delta,
    lambda = if (lambda_present) lambda else NULL,
    primal_residual = primal_residual,
    primal_tolerance = primal_tolerance,
    dual_residual = dual_residual,
    dual_objective_residual = dual_objective_residual,
    dual_nonnegativity_residual = dual_nonnegativity_residual,
    dual_tolerance = dual_tolerance,
    conditioning_branch = conditioning_branch,
    procedure = problem$procedure,
    scope = geometry$scope,
    lf_critical_value = .sd_ie_or(problem$hybrid_list$lf_cv, NA_real_)
  )
}

# -----------------------------------------------------------------------------
# Component-preserving adaptive inversion
# -----------------------------------------------------------------------------

.sd_ie_normalize_test_result <- function(result, theta) {
  if (is.logical(result) && length(result) == 1L && !is.na(result)) {
    return(list(
      success = TRUE, status = "ok", message = NULL, warnings = character(),
      theta = theta, accept = result
    ))
  }
  if (!is.list(result)) {
    return(.sd_ie_failure(
      "malformed_test_result",
      "test_fn must return one logical value or a list with success and accept.",
      theta = theta, accept = NA
    ))
  }
  if (!isTRUE(result$success)) {
    return(.sd_ie_failure(
      .sd_ie_or(result$status, "test_failure"),
      .sd_ie_or(result$message, "The test failed without a diagnostic message."),
      .sd_ie_or(result$warnings, character()),
      theta = theta, accept = NA, raw = result
    ))
  }
  if (!is.logical(result$accept) || length(result$accept) != 1L ||
      is.na(result$accept)) {
    return(.sd_ie_failure(
      "malformed_test_result",
      "A successful test result must contain one nonmissing logical accept.",
      .sd_ie_or(result$warnings, character()),
      theta = theta, accept = NA, raw = result
    ))
  }
  result$theta <- theta
  result
}

sd_invert_test <- function(
    test_fn,
    topology_grid,
    endpoint_tolerance = 1e-8,
    max_expansions = 2L,
    expansion_points = NULL,
    max_bisection = 100L,
    max_test_evaluations = 100000L,
    truth = NULL,
    expansion_step = NULL) {
  if (!is.function(test_fn)) .sd_ie_stop("test_fn must be a function.")
  grid <- sort(unique(as.numeric(topology_grid)))
  if (length(grid) < 3L || any(!is.finite(grid)) || any(diff(grid) <= 0)) {
    .sd_ie_stop("topology_grid must contain at least three distinct finite points.")
  }
  .sd_ie_scalar(endpoint_tolerance, "endpoint_tolerance")
  if (endpoint_tolerance <= 0) .sd_ie_stop("endpoint_tolerance must be positive.")
  .sd_ie_scalar(max_expansions, "max_expansions")
  max_expansions <- as.integer(max_expansions)
  if (max_expansions < 0L) .sd_ie_stop("max_expansions must be nonnegative.")
  .sd_ie_scalar(max_bisection, "max_bisection")
  max_bisection <- as.integer(max_bisection)
  if (max_bisection < 1L) .sd_ie_stop("max_bisection must be positive.")
  .sd_ie_scalar(max_test_evaluations, "max_test_evaluations")
  max_test_evaluations <- as.integer(max_test_evaluations)
  if (max_test_evaluations < length(grid)) {
    return(.sd_ie_failure("evaluation_limit",
      "The actual initial topology grid exceeds the point-evaluation budget.",
      complete = FALSE, components = data.frame(), total_test_evaluations = 0L))
  }
  if (is.null(expansion_points)) expansion_points <- length(grid) - 1L
  .sd_ie_scalar(expansion_points, "expansion_points")
  expansion_points <- as.integer(expansion_points)
  if (expansion_points < 1L) .sd_ie_stop("expansion_points must be positive.")
  if (!is.null(truth)) .sd_ie_scalar(truth, "truth")

  base_step <- if (is.null(expansion_step)) min(diff(grid)) else expansion_step
  .sd_ie_scalar(base_step, "expansion_step")
  if (base_step <= 0) .sd_ie_stop("expansion_step must be positive.")
  cache <- new.env(hash = TRUE, parent = emptyenv())
  n_evaluations <- 0L
  evaluate <- function(theta) {
    key <- sprintf("%.17g", theta)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    if (n_evaluations >= max_test_evaluations) {
      return(.sd_ie_failure(
        "evaluation_limit",
        paste(
          "Inversion reached the registered limit of",
          max_test_evaluations, "test evaluations."
        ),
        theta = theta, accept = NA
      ))
    }
    captured <- .sd_ie_capture(function() test_fn(theta))
    if (!is.null(captured$error)) {
      result <- .sd_ie_failure(
        "test_function_error", captured$error, captured$warnings,
        theta = theta, accept = NA
      )
    } else {
      result <- .sd_ie_normalize_test_result(captured$value, theta)
      result$warnings <- unique(c(result$warnings, captured$warnings))
    }
    assign(key, result, envir = cache)
    n_evaluations <<- n_evaluations + 1L
    result
  }
  evaluate_grid <- function(x) lapply(x, evaluate)
  scan_results <- evaluate_grid(grid)

  failed <- which(!vapply(scan_results, function(x) isTRUE(x$success), logical(1L)))
  if (length(failed) > 0L) {
    return(list(
      success = FALSE,
      complete = FALSE,
      status = "test_failure",
      message = "At least one topology test failed; inversion was not completed.",
      warnings = unique(unlist(lapply(scan_results, `[[`, "warnings"))),
      failed_points = grid[failed],
      failed_results = scan_results[failed],
      scan = data.frame(
        grid = grid,
        accept = vapply(scan_results, function(x) {
          if (isTRUE(x$success)) as.integer(x$accept) else NA_integer_
        }, integer(1L))
      ),
      total_test_evaluations = n_evaluations,
      components = data.frame()
    ))
  }

  accept <- vapply(scan_results, `[[`, logical(1L), "accept")
  expansions <- 0L
  while ((accept[1L] || accept[length(accept)]) &&
         expansions < max_expansions) {
    # Double the newly searched span after every unsuccessful boundary check.
    # Keeping the original mesh while expanding geometrically avoids skipping
    # narrow accept/reject transitions in the added tails.
    expansion_scale <- 2^expansions
    points_this_expansion <- as.integer(expansion_points * expansion_scale)
    add_left <- if (accept[1L]) {
      grid[1L] - base_step * rev(seq_len(points_this_expansion))
    } else {
      numeric()
    }
    add_right <- if (accept[length(accept)]) {
      grid[length(grid)] + base_step * seq_len(points_this_expansion)
    } else {
      numeric()
    }
    if (n_evaluations + length(add_left) + length(add_right) >
        max_test_evaluations) {
      return(list(
        success = FALSE,
        complete = FALSE,
        status = "evaluation_limit",
        message = paste(
          "Adaptive topology expansion would exceed the registered limit of",
          max_test_evaluations, "test evaluations."
        ),
        warnings = unique(unlist(lapply(scan_results, `[[`, "warnings"))),
        scan = data.frame(grid = grid, accept = as.integer(accept)),
        total_test_evaluations = n_evaluations,
        components = data.frame()
      ))
    }
    grid <- c(add_left, grid, add_right)
    scan_results <- evaluate_grid(grid)
    failed <- which(!vapply(
      scan_results, function(x) isTRUE(x$success), logical(1L)
    ))
    if (length(failed) > 0L) {
      return(list(
        success = FALSE,
        complete = FALSE,
        status = "test_failure",
        message = "A test failed while expanding the topology search.",
        warnings = unique(unlist(lapply(scan_results, `[[`, "warnings"))),
        failed_points = grid[failed],
        failed_results = scan_results[failed],
        scan = data.frame(
          grid = grid,
          accept = vapply(scan_results, function(x) {
            if (isTRUE(x$success)) as.integer(x$accept) else NA_integer_
          }, integer(1L))
        ),
        total_test_evaluations = n_evaluations,
        components = data.frame()
      ))
    }
    accept <- vapply(scan_results, `[[`, logical(1L), "accept")
    expansions <- expansions + 1L
  }

  truth_test <- if (is.null(truth)) NULL else evaluate(truth)
  accepted <- which(accept)
  warnings <- unique(unlist(lapply(scan_results, `[[`, "warnings")))
  scan <- data.frame(grid = grid, accept = as.integer(accept))
  if (!is.null(truth_test) && !isTRUE(truth_test$success)) {
    return(list(
      success = FALSE,
      complete = FALSE,
      status = "truth_test_failure",
      message = paste(
        "The explicitly requested truth-point test failed; coverage cannot",
        "be recorded for this inversion."
      ),
      warnings = unique(c(warnings, truth_test$warnings)),
      scan = scan,
      components = data.frame(),
      total_test_evaluations = n_evaluations,
      truth_test = truth_test,
      contains_truth = NA
    ))
  }
  if (length(accepted) == 0L) {
    truth_topology_miss <- !is.null(truth_test) &&
      isTRUE(truth_test$success) && isTRUE(truth_test$accept)
    return(list(
      success = !truth_topology_miss,
      complete = !truth_topology_miss,
      status = if (truth_topology_miss) "topology_miss" else "empty",
      message = if (truth_topology_miss) paste(
        "The supplied diagnostic truth point was accepted but the topology",
        "grid found no corresponding confidence-set component."
      ) else NULL,
      warnings = warnings,
      scan = scan,
      lower = NA_real_,
      upper = NA_real_,
      total_length = 0,
      truncated_length = 0,
      n_components = 0L,
      components = data.frame(
        lower = numeric(), upper = numeric(),
        lower_search_boundary = logical(),
        upper_search_boundary = logical()
      ),
      boundary_accepted = FALSE,
      expansions = expansions,
      topology_mesh = max(diff(grid)),
      endpoint_tolerance = endpoint_tolerance,
      total_test_evaluations = n_evaluations,
      truth_test = truth_test,
      contains_truth = if (is.null(truth_test)) NA else truth_test$accept
    ))
  }

  groups <- split(accepted, cumsum(c(TRUE, diff(accepted) > 1L)))
  components <- data.frame(
    lower = rep(NA_real_, length(groups)),
    upper = rep(NA_real_, length(groups)),
    lower_search_boundary = FALSE,
    upper_search_boundary = FALSE
  )
  bisection_evaluations <- 0L
  bisection_failure <- NULL

  bisect_transition <- function(left, right, accept_left, accept_right) {
    if (identical(accept_left, accept_right)) {
      return(.sd_ie_failure(
        "invalid_bisection_bracket", "Bisection bracket has no transition."
      ))
    }
    iterations <- 0L
    while (right - left > endpoint_tolerance && iterations < max_bisection) {
      midpoint <- (left + right) / 2
      middle <- evaluate(midpoint)
      if (!isTRUE(middle$success)) return(middle)
      if (identical(middle$accept, accept_left)) {
        left <- midpoint
      } else {
        right <- midpoint
      }
      iterations <- iterations + 1L
    }
    if (right - left > endpoint_tolerance) {
      return(.sd_ie_failure(
        "bisection_limit",
        "Endpoint bisection reached max_bisection before its tolerance."
      ))
    }
    list(
      success = TRUE,
      boundary = (left + right) / 2,
      iterations = iterations
    )
  }

  for (component_index in seq_along(groups)) {
    lower_index <- min(groups[[component_index]])
    upper_index <- max(groups[[component_index]])
    if (lower_index == 1L) {
      components$lower[component_index] <- grid[1L]
      components$lower_search_boundary[component_index] <- TRUE
    } else {
      edge <- bisect_transition(
        grid[lower_index - 1L], grid[lower_index],
        accept[lower_index - 1L], accept[lower_index]
      )
      if (!isTRUE(edge$success)) {
        bisection_failure <- edge
        break
      }
      components$lower[component_index] <- edge$boundary
      bisection_evaluations <- bisection_evaluations + edge$iterations
    }
    if (upper_index == length(grid)) {
      components$upper[component_index] <- grid[length(grid)]
      components$upper_search_boundary[component_index] <- TRUE
    } else {
      edge <- bisect_transition(
        grid[upper_index], grid[upper_index + 1L],
        accept[upper_index], accept[upper_index + 1L]
      )
      if (!isTRUE(edge$success)) {
        bisection_failure <- edge
        break
      }
      components$upper[component_index] <- edge$boundary
      bisection_evaluations <- bisection_evaluations + edge$iterations
    }
  }
  if (!is.null(bisection_failure)) {
    return(list(
      success = FALSE,
      complete = FALSE,
      status = bisection_failure$status,
      message = bisection_failure$message,
      warnings = unique(c(warnings, bisection_failure$warnings)),
      scan = scan,
      components = components,
      total_test_evaluations = n_evaluations,
      bisection_evaluations = bisection_evaluations,
      truth_test = truth_test
    ))
  }

  boundary_accepted <- any(components$lower_search_boundary) ||
    any(components$upper_search_boundary)
  truncated_length <- sum(components$upper - components$lower)
  truth_topology_miss <- !is.null(truth_test) &&
    isTRUE(truth_test$success) && isTRUE(truth_test$accept) &&
    !any(truth >= components$lower - endpoint_tolerance &
           truth <= components$upper + endpoint_tolerance)
  complete <- !boundary_accepted && !truth_topology_miss
  list(
    success = complete,
    complete = complete,
    status = if (complete) {
      "ok"
    } else if (truth_topology_miss) {
      "topology_miss"
    } else {
      "open_search_boundary"
    },
    message = if (complete) {
      NULL
    } else if (truth_topology_miss) {
      paste(
        "The supplied diagnostic truth point was accepted but lies outside",
        "every component recovered from the topology grid."
      )
    } else {
      paste(
        "An accepted component still reaches a search boundary after",
        expansions, "expansion(s); total length is unresolved."
      )
    },
    warnings = warnings,
    scan = scan,
    lower = min(components$lower),
    upper = max(components$upper),
    total_length = if (complete) truncated_length else NA_real_,
    truncated_length = truncated_length,
    n_components = nrow(components),
    components = components,
    boundary_accepted = boundary_accepted,
    expansions = expansions,
    topology_mesh = max(diff(grid)),
    endpoint_tolerance = endpoint_tolerance,
    bisection_evaluations = bisection_evaluations,
    total_test_evaluations = n_evaluations,
    truth_test = truth_test,
    contains_truth = if (is.null(truth_test) || !isTRUE(truth_test$success)) {
      NA
    } else {
      truth_test$accept
    }
  )
}

sd_arp_invert <- function(
    problem,
    topology_grid,
    endpoint_tolerance = 1e-8,
    max_expansions = 2L,
    expansion_points = NULL,
    max_bisection = 100L,
    max_test_evaluations = 100000L,
    truth = NULL,
    profile_anchor = NULL,
    expansion_step = NULL,
    test_fn = NULL) {
  if (!inherits(problem, "sd_arp_problem")) {
    .sd_ie_stop("problem must be returned by sd_bind_arp_problem().")
  }
  if (is.null(profile_anchor)) profile_anchor <- sd_arp_profile_anchor(problem)
  if (!is.list(profile_anchor) || !isTRUE(profile_anchor$success) ||
      length(profile_anchor$theta) != 1L ||
      !is.finite(profile_anchor$theta) ||
      length(profile_anchor$eta) != 1L || !is.finite(profile_anchor$eta)) {
    anchor_message <- if (is.list(profile_anchor)) {
      .sd_ie_or(
        profile_anchor$message,
        "The data-dependent profile anchor was unavailable."
      )
    } else {
      "The supplied profile anchor is malformed."
    }
    anchor_warnings <- if (is.list(profile_anchor)) {
      .sd_ie_or(profile_anchor$warnings, character())
    } else character()
    return(.sd_ie_failure(
      "arp_profile_anchor_unavailable",
      anchor_message,
      anchor_warnings,
      complete = FALSE,
      components = data.frame(),
      procedure = .sd_ie_or(problem$procedure, NA_character_),
      scope = .sd_ie_or(problem$geometry$scope, NA_character_),
      profile_anchor = profile_anchor
    ))
  }
  result <- sd_invert_test(
    test_fn = if (is.null(test_fn)) function(theta) sd_arp_truth_test(problem, theta) else test_fn,
    topology_grid = topology_grid,
    endpoint_tolerance = endpoint_tolerance,
    max_expansions = max_expansions,
    expansion_points = expansion_points,
    max_bisection = max_bisection,
    max_test_evaluations = max_test_evaluations,
    truth = truth,
    expansion_step = expansion_step
  )

  # A finite all-rejected grid is not, by itself, a certificate that the
  # inverted confidence set is empty: an accepted component can lie beyond a
  # rejected search boundary.  The C-LF first stage supplies an exact empty-set
  # certificate only when the globally minimized eta is itself above the LF
  # critical value.  Otherwise fail closed and let the caller enlarge or
  # recenter the topology search rather than record a false empty set.
  empty_certificate <- FALSE
  if (isTRUE(result$success) && identical(result$status, "empty") &&
      identical(problem$procedure, "c_lf")) {
    lf_cv <- .sd_ie_or(problem$hybrid_list$lf_cv, NA_real_)
    certificate_tolerance <- max(
      1e-10,
      .sd_ie_or(profile_anchor$dual_tolerance, 0),
      .sd_ie_or(profile_anchor$dual_objective_residual, 0)
    )
    empty_certificate <- is.finite(lf_cv) &&
      is.finite(profile_anchor$eta_lower_bound) &&
      profile_anchor$eta_lower_bound > lf_cv + certificate_tolerance
  }
  if (isTRUE(result$success) && identical(result$status, "empty") &&
      !empty_certificate) {
    result$success <- FALSE
    result$complete <- FALSE
    result$status <- "uncertified_empty_search"
    result$message <- paste(
      "The finite topology grid contained no accepted point, but this does",
      "not certify an empty ARP confidence set."
    )
  }
  result$procedure <- .sd_ie_or(problem$procedure, NA_character_)
  result$scope <- .sd_ie_or(problem$geometry$scope, NA_character_)
  result$profile_anchor <- profile_anchor
  result$empty_certified <- empty_certificate
  result
}

sd_arp_audit_endpoints <- function(
    problem,
    interval,
    probe_distance = 2e-8,
    minimum_probe_multiplier = 2,
    test_fn = NULL) {
  if (is.null(test_fn)) test_fn <- function(theta) sd_arp_truth_test(problem, theta)
  if (!inherits(problem, "sd_arp_problem")) {
    .sd_ie_stop("problem must be returned by sd_bind_arp_problem().")
  }
  if (!is.list(interval) || !isTRUE(interval$success)) {
    return(.sd_ie_failure(
      "endpoint_audit_unavailable",
      "Endpoint auditing requires a successful, complete inversion.",
      checks = data.frame()
    ))
  }
  .sd_ie_scalar(probe_distance, "probe_distance")
  .sd_ie_scalar(minimum_probe_multiplier, "minimum_probe_multiplier")
  if (probe_distance <= 0 || minimum_probe_multiplier <= 1) {
    .sd_ie_stop("Endpoint-audit probe parameters must be positive.")
  }
  components <- interval$components
  if (is.null(components) || !nrow(components)) {
    return(list(
      success = TRUE, status = "empty_no_endpoints", message = NULL,
      warnings = character(), certified = TRUE, checks = data.frame()
    ))
  }
  endpoint_tolerance <- .sd_ie_or(interval$endpoint_tolerance, 0)
  rows <- list()
  warnings <- character()
  all_ok <- TRUE
  for (component_index in seq_len(nrow(components))) {
    lower <- components$lower[component_index]
    upper <- components$upper[component_index]
    width <- upper - lower
    left_gap <- if (component_index == 1L) {
      Inf
    } else lower - components$upper[component_index - 1L]
    right_gap <- if (component_index == nrow(components)) {
      Inf
    } else components$lower[component_index + 1L] - upper
    distances <- c(
      lower_inside = min(probe_distance, width / 4),
      lower_outside = min(probe_distance, left_gap / 4),
      upper_inside = min(probe_distance, width / 4),
      upper_outside = min(probe_distance, right_gap / 4)
    )
    points <- c(
      lower_inside = lower + distances[["lower_inside"]],
      lower_outside = lower - distances[["lower_outside"]],
      upper_inside = upper - distances[["upper_inside"]],
      upper_outside = upper + distances[["upper_outside"]]
    )
    expected <- c(TRUE, FALSE, TRUE, FALSE)
    for (probe_index in seq_along(points)) {
      distance_ok <- is.finite(distances[probe_index]) &&
        distances[probe_index] >
          minimum_probe_multiplier * endpoint_tolerance
      tested <- if (distance_ok) {
        .sd_ie_normalize_test_result(
          test_fn(points[probe_index]),
          points[probe_index]
        )
      } else {
        .sd_ie_failure(
          "endpoint_probe_too_close",
          "A component or gap is too narrow for a reliable endpoint probe.",
          theta = points[probe_index], accept = NA
        )
      }
      observed <- if (isTRUE(tested$success)) tested$accept else NA
      passed <- isTRUE(tested$success) && identical(observed, expected[probe_index])
      all_ok <- all_ok && passed
      warnings <- unique(c(warnings, tested$warnings))
      rows[[length(rows) + 1L]] <- data.frame(
        component_id = component_index,
        probe = names(points)[probe_index],
        theta = unname(points[probe_index]),
        distance = unname(distances[probe_index]),
        expected_accept = expected[probe_index],
        observed_accept = observed,
        success = isTRUE(tested$success),
        pass = passed,
        status = .sd_ie_or(tested$status, "missing"),
        primal_residual = .sd_ie_or(tested$primal_residual, NA_real_),
        dual_residual = .sd_ie_or(tested$dual_residual, NA_real_),
        dual_objective_residual =
          .sd_ie_or(tested$dual_objective_residual, NA_real_),
        stringsAsFactors = FALSE
      )
    }
  }
  checks <- do.call(rbind, rows)
  rownames(checks) <- NULL
  list(
    success = all_ok,
    status = if (all_ok) "ok" else "endpoint_audit_failure",
    message = if (all_ok) NULL else
      "At least one inside/outside endpoint probe failed.",
    warnings = warnings,
    certified = all_ok,
    checks = checks
  )
}

sd_arp_topology_audit <- function(
    problem,
    base_interval,
    topology_grid,
    endpoint_tolerance = 1e-8,
    max_expansions = 2L,
    max_bisection = 100L,
    max_test_evaluations = 100000L,
    agreement_multiplier = 4,
    profile_anchor = NULL,
    expansion_step = NULL,
    test_fn = NULL) {
  if (is.null(profile_anchor)) profile_anchor <- sd_arp_profile_anchor(problem)
  if (!is.list(profile_anchor) || !isTRUE(profile_anchor$success) ||
      length(profile_anchor$theta) != 1L || !is.finite(profile_anchor$theta)) {
    anchor_message <- if (is.list(profile_anchor)) {
      .sd_ie_or(
        profile_anchor$message,
        "The topology audit could not construct its profile anchor."
      )
    } else {
      "The supplied topology-audit profile anchor is malformed."
    }
    anchor_warnings <- if (is.list(profile_anchor)) {
      .sd_ie_or(profile_anchor$warnings, character())
    } else character()
    return(.sd_ie_failure(
      "arp_profile_anchor_unavailable",
      anchor_message,
      anchor_warnings,
      shifted = list(success = FALSE, components = data.frame()),
      fine = list(success = FALSE, components = data.frame())
    ))
  }
  grid <- sort(unique(as.numeric(topology_grid)))
  midpoints <- head(grid, -1L) + diff(grid) / 2
  shifted_grid <- sort(unique(c(
    grid[1L], midpoints, grid[length(grid)], profile_anchor$theta
  )))
  fine_grid <- sort(unique(c(grid, midpoints, profile_anchor$theta)))
  run_one <- function(audit_grid) {
    sd_arp_invert(
      problem = problem,
      topology_grid = audit_grid,
      endpoint_tolerance = endpoint_tolerance,
      max_expansions = max_expansions,
      max_bisection = max_bisection,
      max_test_evaluations = max_test_evaluations,
      truth = NULL,
      profile_anchor = profile_anchor,
      expansion_step = if (is.null(expansion_step)) NULL else expansion_step / 2,
      test_fn = test_fn
    )
  }
  shifted <- run_one(shifted_grid)
  fine <- run_one(fine_grid)
  compare_one <- function(candidate) {
    same_count <- isTRUE(base_interval$success) && isTRUE(candidate$success) &&
      nrow(base_interval$components) == nrow(candidate$components)
    endpoint_gap <- if (same_count && nrow(candidate$components)) {
      max(abs(
        as.matrix(base_interval$components[, c("lower", "upper")]) -
          as.matrix(candidate$components[, c("lower", "upper")])
      ))
    } else if (same_count) 0 else Inf
    tolerance <- agreement_multiplier * endpoint_tolerance
    list(
      pass = same_count && is.finite(endpoint_gap) && endpoint_gap <= tolerance,
      same_count = same_count,
      endpoint_gap = endpoint_gap,
      tolerance = tolerance
    )
  }
  shifted_comparison <- compare_one(shifted)
  fine_comparison <- compare_one(fine)
  success <- shifted_comparison$pass && fine_comparison$pass
  list(
    success = success,
    status = if (success) "ok" else "topology_audit_failure",
    message = if (success) NULL else paste(
      "Base, half-mesh-shifted, and half-mesh-fine inversions did not agree."
    ),
    shifted = shifted,
    fine = fine,
    shifted_comparison = shifted_comparison,
    fine_comparison = fine_comparison,
    shifted_grid = shifted_grid,
    fine_grid = fine_grid
  )
}
