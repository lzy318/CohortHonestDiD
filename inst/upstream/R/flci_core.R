# Full-path SD FLCI utilities.
#
# This file is self-contained. The package loads it into a private
# environment; it does not alter or replace any other package function.
#
# IMPORTANT API SCOPE: these routines implement the paper's mechanically
# normalized, quotient-reduced special case.  They require a pre-only
# normalization, the reduced W = [[I, 0], [H, C]] block form, a
# square/invertible boundary-post block Q_plus, and full row rank of the
# quotient pure-pre block R.  They deliberately stop on more general designs;
# the generic multiplier value-function program is not implemented here.

flci_stop <- function(..., call. = FALSE) {
  stop(..., call. = call.)
}

flci_symmetrize <- function(x) {
  y <- t(x)
  answer <- x / 2 + y / 2
  safe_sum <- sign(x) != sign(y) |
    (abs(x) <= .Machine$double.xmax / 2 & abs(y) <= .Machine$double.xmax / 2)
  answer[safe_sum] <- (x[safe_sum] + y[safe_sum]) / 2
  equal <- x == y
  answer[equal] <- x[equal]
  answer
}

flci_max_abs <- function(x) {
  if (length(x) == 0L) return(0)
  max(abs(x))
}

flci_scale_covariance_block <- function(x) {
  scale <- max(abs(x))
  if (scale == 0) scale <- 1
  scaled <- x / scale
  if (any(x != 0 & scaled == 0)) {
    flci_stop("FLCI curvature scaling cannot preserve a contributing nonzero covariance entry.")
  }
  list(covariance = scaled, scale = scale)
}

flci_checked_product <- function(x, y) {
  x <- as.matrix(x); y <- as.matrix(y)
  if (any(!is.finite(x)) || any(!is.finite(y)))
    flci_stop("FLCI curvature arithmetic has nonfinite factors.")
  for (j in seq_len(ncol(x))) {
    left <- x[, j]; right <- y[j, ]
    left <- left[left != 0]; right <- right[right != 0]
    if (length(left) && length(right)) {
      products <- outer(left, right)
      if (any(products == 0) || any(!is.finite(products)))
        flci_stop("FLCI curvature arithmetic cannot preserve a contributing nonzero product.")
    }
  }
  answer <- x %*% y
  if (any(!is.finite(answer)))
    flci_stop("FLCI curvature arithmetic has a nonfinite product sum.")
  answer
}

flci_require_names <- function(x, expected_length, label) {
  if (is.null(x) || length(x) != expected_length ||
      anyNA(x) || any(!nzchar(x))) {
    flci_stop(label, " must have one nonempty name per entry.")
  }
  if (anyDuplicated(x)) {
    flci_stop(label, " must be unique.")
  }
  invisible(x)
}

flci_require_same_name_set <- function(
    actual,
    expected,
    label,
    expected_label = "W coordinate names") {
  missing <- setdiff(expected, actual)
  extra <- setdiff(actual, expected)
  if (length(actual) != length(expected) ||
      length(missing) > 0L || length(extra) > 0L) {
    detail <- character(0)
    if (length(missing) > 0L) {
      detail <- c(detail, paste0(
        "missing: ", paste(missing, collapse = ", ")
      ))
    }
    if (length(extra) > 0L) {
      detail <- c(detail, paste0(
        "extra: ", paste(extra, collapse = ", ")
      ))
    }
    flci_stop(
      label, " names must equal the ", expected_label,
      if (length(detail) > 0L) paste0(" (", paste(detail, collapse = "; "), ")")
      else "",
      "."
    )
  }
  invisible(actual)
}

flci_finite_numeric_matrix <- function(x, label) {
  x <- as.matrix(x)
  if (!is.numeric(x)) flci_stop(label, " must be numeric.")
  if (any(!is.finite(x))) flci_stop(label, " must contain only finite values.")
  storage.mode(x) <- "double"
  x
}

flci_finite_numeric_vector <- function(x, label) {
  if (!is.numeric(x)) flci_stop(label, " must be numeric.")
  x <- as.numeric(x)
  if (any(!is.finite(x))) flci_stop(label, " must contain only finite values.")
  x
}

flci_scaled_residual_tol <- function(x, tol) {
  tol * max(1, flci_max_abs(x))
}

flci_safe_sqrt_variance <- function(
    variance,
    scale,
    rel_tol = 1e-10,
    label = "variance") {
  if (length(variance) != 1L || !is.finite(variance)) {
    flci_stop(label, " must be one finite scalar.")
  }
  if (length(scale) != 1L || !is.finite(scale) || scale < 0) {
    flci_stop(label, " scale must be one finite nonnegative scalar.")
  }
  clip_tol <- rel_tol * max(scale, .Machine$double.xmin)
  if (variance < -clip_tol) {
    flci_stop(
      label, " is materially negative (", format(variance, digits = 16),
      "; clipping tolerance ", format(clip_tol, digits = 6), ")."
    )
  }
  sqrt(if (variance < 0) 0 else variance)
}

flci_validate_solver_opts <- function(solver_opts) {
  if (!is.list(solver_opts)) flci_stop("solver_opts must be a named list.")
  if (length(solver_opts) == 0L) return(invisible(solver_opts))
  option_names <- names(solver_opts)
  if (is.null(option_names) || anyNA(option_names) ||
      any(!nzchar(option_names)) || anyDuplicated(option_names)) {
    flci_stop("solver_opts must be a uniquely named list.")
  }
  reserved <- intersect(option_names, c("object", "a", "b", "solver"))
  if (length(reserved) > 0L) {
    flci_stop(
      "solver_opts cannot override: ", paste(reserved, collapse = ", "), "."
    )
  }
  invisible(solver_opts)
}

flci_cvxr_solve <- function(problem, solver, solver_opts = list()) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    flci_stop("CVXR is required for this exploratory optimization.")
  }
  if (length(solver) != 1L || !is.character(solver) ||
      is.na(solver) || !nzchar(solver)) {
    flci_stop("solver must be one nonempty character string.")
  }
  flci_validate_solver_opts(solver_opts)
  installed <- CVXR::installed_solvers()
  if (!(solver %in% installed)) {
    flci_stop(
      "Requested CVXR solver '", solver,
      "' is not installed. Installed solvers: ",
      paste(installed, collapse = ", "), "."
    )
  }
  do.call(
    CVXR::solve,
    c(list(problem, solver = solver), solver_opts)
  )
}

flci_assert_cvxr_status <- function(status, context) {
  if (length(status) != 1L || is.na(status)) {
    flci_stop(context, " returned a missing or malformed solver status.")
  }
  normalized <- tolower(as.character(status))
  accepted <- c("optimal", "optimal_inaccurate")
  if (!(normalized %in% accepted)) {
    flci_stop(
      context, " failed with CVXR status '", as.character(status),
      "'. Accepted statuses are: ", paste(accepted, collapse = ", "), "."
    )
  }
  normalized
}

flci_folded_cv <- function(r, alpha = 0.05) {
  if (length(r) != 1L || !is.finite(r) || r < 0) {
    flci_stop("r must be one finite nonnegative scalar.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 0.5) {
    flci_stop("alpha must lie strictly between zero and one half.")
  }
  if (r == 0) return(stats::qnorm(1 - alpha / 2))

  score <- function(c_value) {
    stats::pnorm(c_value - r) +
      stats::pnorm(c_value + r) -
      (2 - alpha)
  }
  upper <- r + stats::qnorm(1 - alpha) + 8
  stats::uniroot(score, interval = c(0, upper), tol = 1e-12)$root
}

flci_half_length <- function(B, s, alpha = 0.05, zero_tol = 1e-14) {
  if (length(B) != 1L || !is.finite(B) || B < -zero_tol) {
    flci_stop("B must be one finite nonnegative scalar.")
  }
  if (length(s) != 1L || !is.finite(s) || s < -zero_tol) {
    flci_stop("s must be one finite nonnegative scalar.")
  }
  B <- max(B, 0)
  s <- max(s, 0)
  if (s == 0) return(B)
  s * flci_folded_cv(B / s, alpha = alpha)
}

flci_null_basis <- function(N, tol = 1e-10) {
  N <- as.matrix(N)
  d <- ncol(N)
  if (is.null(d)) flci_stop("N must be a matrix.")
  if (nrow(N) == 0L) return(diag(d))

  qr_N <- qr(t(N), tol = tol, LAPACK = FALSE)
  rank_N <- qr_N$rank
  Q <- qr.Q(qr_N, complete = TRUE)
  if (rank_N == d) return(matrix(numeric(0), nrow = d, ncol = 0L))
  Q[, seq.int(rank_N + 1L, d), drop = FALSE]
}

flci_canonical_names <- function(T_cohort, T_total) {
  pre <- list()
  post <- list()
  for (g in T_cohort) {
    n_pre <- g - 1L
    if (n_pre > 0L) {
      pre[[as.character(g)]] <- paste0(
        "delta_pre_g", g, "_s", -c((n_pre - 1L):0)
      )
    }
    n_post <- T_total - g + 1L
    if (n_post > 0L) {
      post[[as.character(g)]] <- paste0(
        "delta_post_g", g, "_s", seq_len(n_post)
      )
    }
  }
  c(unlist(pre, use.names = FALSE), unlist(post, use.names = FALSE))
}

flci_format_saved_results <- function(results_obj, T_cohort, T_total) {
  all_effects <- numeric(0)

  for (cohort_id in names(results_obj$point_estimates$cohort_pre_estimates)) {
    cohort_df <- results_obj$point_estimates$cohort_pre_estimates[[cohort_id]]
    for (i in seq_len(nrow(cohort_df))) {
      nm <- paste0(
        "delta_pre_g", cohort_id, "_s", cohort_df$relative_period[i]
      )
      all_effects[nm] <- cohort_df$effect[i]
    }
  }
  for (cohort_id in names(results_obj$point_estimates$cohort_post_estimates)) {
    cohort_df <- results_obj$point_estimates$cohort_post_estimates[[cohort_id]]
    for (i in seq_len(nrow(cohort_df))) {
      nm <- paste0(
        "delta_post_g", cohort_id, "_s", cohort_df$relative_period[i]
      )
      all_effects[nm] <- cohort_df$effect[i]
    }
  }

  Sigma <- results_obj$inference$vcov_cohort_all_effects
  rename_one <- function(x) {
    x <- gsub(
      "cohort_pre_([0-9]+)_s_m([0-9]+)",
      "delta_pre_g\\1_s-\\2",
      x
    )
    x <- gsub(
      "cohort_pre_([0-9]+)_s([0-9]+)",
      "delta_pre_g\\1_s\\2",
      x
    )
    gsub(
      "cohort_post_([0-9]+)_s([0-9]+)",
      "delta_post_g\\1_s\\2",
      x
    )
  }
  rownames(Sigma) <- rename_one(rownames(Sigma))
  colnames(Sigma) <- rename_one(colnames(Sigma))

  canonical <- flci_canonical_names(T_cohort, T_total)
  missing_beta <- setdiff(canonical, names(all_effects))
  missing_sigma <- setdiff(canonical, rownames(Sigma))
  if (length(missing_beta) > 0L) {
    flci_stop("Missing saved coefficients: ", paste(missing_beta, collapse = ", "))
  }
  if (length(missing_sigma) > 0L) {
    flci_stop("Missing covariance coordinates: ",
              paste(missing_sigma, collapse = ", "))
  }

  list(
    betahat = all_effects[canonical],
    Sigma = Sigma[canonical, canonical, drop = FALSE],
    canonical_names = canonical
  )
}

flci_normalization_matrix <- function(coordinate_names, estimator_type) {
  pre_names <- coordinate_names[grepl("^delta_pre_", coordinate_names)]
  cohorts <- unique(sub(
    "^delta_pre_(g[0-9]+)_.*$",
    "\\1",
    pre_names
  ))
  rows <- list()

  if (estimator_type == "imputation") {
    for (g in cohorts) {
      use <- grepl(paste0("^delta_pre_", g, "_"), coordinate_names)
      row <- numeric(length(coordinate_names))
      row[use] <- 1 / sum(use)
      rows[[paste0("mean_pre_", g)]] <- row
    }
  } else if (estimator_type %in%
             c("cs-not-yet-treated", "cs-never-treated")) {
    for (g in cohorts) {
      target <- paste0("delta_pre_", g, "_s0")
      if (target %in% coordinate_names) {
        row <- numeric(length(coordinate_names))
        row[match(target, coordinate_names)] <- 1
        rows[[paste0("reference_", g)]] <- row
      }
    }
  } else {
    flci_stop("Unknown estimator_type: ", estimator_type)
  }

  if (length(rows) == 0L) {
    ans <- matrix(numeric(0), nrow = 0L, ncol = length(coordinate_names))
  } else {
    ans <- do.call(rbind, rows)
  }
  colnames(ans) <- coordinate_names
  ans
}

flci_oriented_sd_blocks <- function(A_Delta, d_vec, tol = 1e-12) {
  if (is.null(rownames(A_Delta))) {
    flci_stop("A_Delta must have row names identifying positive rows.")
  }
  keep <- grepl("_pos$", rownames(A_Delta))
  D2 <- A_Delta[keep, , drop = FALSE]
  bounds <- as.numeric(d_vec[rownames(D2)])
  names(bounds) <- rownames(D2)

  pre <- grepl("^delta_pre_", colnames(D2))
  post <- grepl("^delta_post_", colnames(D2))
  pure_pre <- rowSums(abs(D2[, post, drop = FALSE]) > tol) == 0L

  list(
    R = D2[pure_pre, pre, drop = FALSE],
    Q_minus = D2[!pure_pre, pre, drop = FALSE],
    Q_plus = D2[!pure_pre, post, drop = FALSE],
    m_R = bounds[pure_pre],
    m_Q = bounds[!pure_pre],
    oriented_rows = D2,
    pure_pre_rows = rownames(D2)[pure_pre],
    boundary_post_rows = rownames(D2)[!pure_pre]
  )
}

flci_make_target <- function(coordinate_names, N_cohort) {
  # Project estimand convention: every retained post cell in cohort g receives
  # raw weight N_g, followed by one global normalization.  Thus cohort g's
  # aggregate weight is proportional to N_g times its retained post-cell count;
  # this is not a within-cohort-first average.
  flci_require_names(
    coordinate_names, length(coordinate_names), "coordinate_names"
  )
  if (!is.numeric(N_cohort) || any(!is.finite(N_cohort)) ||
      any(N_cohort < 0)) {
    flci_stop("N_cohort must contain finite nonnegative numeric weights.")
  }
  flci_require_names(names(N_cohort), length(N_cohort), "N_cohort names")
  target <- setNames(numeric(length(coordinate_names)), coordinate_names)
  cohort_names <- paste0("g", names(N_cohort))
  for (j in seq_along(cohort_names)) {
    use <- grepl(
      paste0("^delta_post_", cohort_names[j], "_"),
      coordinate_names
    )
    if (!any(use)) {
      flci_stop(
        "N_cohort entry '", names(N_cohort)[j],
        "' has no retained post coordinate."
      )
    }
    target[use] <- N_cohort[j]
  }
  if (sum(target) <= 0) flci_stop("Target has no positive post weight.")
  target / sum(target)
}

flci_reduce_model <- function(
    W,
    Sigma,
    A_Delta,
    d_vec,
    target,
    normalization,
    estimator_type,
    tol = 1e-9,
    allow_singular_covariance = FALSE) {
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0 || tol >= 1) {
    flci_stop("tol must be one finite scalar strictly between zero and one.")
  }
  if (length(estimator_type) != 1L || !is.character(estimator_type) ||
      is.na(estimator_type) || !nzchar(estimator_type)) {
    flci_stop("estimator_type must be one nonempty character string.")
  }
  if (length(allow_singular_covariance) != 1L ||
      is.na(allow_singular_covariance)) {
    flci_stop("allow_singular_covariance must be TRUE or FALSE.")
  }
  allow_singular_covariance <- isTRUE(allow_singular_covariance)

  # Validate and align every named input before numeric coercion.  The W column
  # order defines the internal coordinate order; every other input is reordered
  # explicitly by name.
  if (length(dim(W)) != 2L || nrow(W) != ncol(W) || nrow(W) == 0L) {
    flci_stop("W must be a nonempty square matrix-like object.")
  }
  flci_require_names(colnames(W), ncol(W), "W column names")
  flci_require_names(rownames(W), nrow(W), "W row names")
  coordinate_names <- colnames(W)
  flci_require_same_name_set(rownames(W), coordinate_names, "W row")
  W <- W[coordinate_names, coordinate_names, drop = FALSE]

  if (length(dim(Sigma)) != 2L ||
      nrow(Sigma) != length(coordinate_names) ||
      ncol(Sigma) != length(coordinate_names)) {
    flci_stop("Sigma must be square with one row and column per W coordinate.")
  }
  flci_require_names(rownames(Sigma), nrow(Sigma), "Sigma row names")
  flci_require_names(colnames(Sigma), ncol(Sigma), "Sigma column names")
  flci_require_same_name_set(rownames(Sigma), coordinate_names, "Sigma row")
  flci_require_same_name_set(colnames(Sigma), coordinate_names, "Sigma column")
  Sigma <- Sigma[coordinate_names, coordinate_names, drop = FALSE]

  if (length(dim(A_Delta)) != 2L ||
      ncol(A_Delta) != length(coordinate_names) ||
      nrow(A_Delta) == 0L) {
    flci_stop(
      "A_Delta must be a nonempty matrix-like object with one column ",
      "per W coordinate."
    )
  }
  flci_require_names(rownames(A_Delta), nrow(A_Delta),
                     "A_Delta row names")
  flci_require_names(colnames(A_Delta), ncol(A_Delta),
                     "A_Delta column names")
  flci_require_same_name_set(
    colnames(A_Delta), coordinate_names, "A_Delta column"
  )
  restriction_names <- rownames(A_Delta)
  A_Delta <- A_Delta[restriction_names, coordinate_names, drop = FALSE]

  if (length(dim(normalization)) != 2L ||
      ncol(normalization) != length(coordinate_names)) {
    flci_stop(
      "normalization must be matrix-like with one column per W coordinate."
    )
  }
  flci_require_names(
    colnames(normalization), ncol(normalization),
    "normalization column names"
  )
  flci_require_same_name_set(
    colnames(normalization), coordinate_names, "normalization column"
  )
  if (nrow(normalization) > 0L) {
    flci_require_names(
      rownames(normalization), nrow(normalization),
      "normalization row names"
    )
  }
  normalization <- normalization[, coordinate_names, drop = FALSE]

  flci_require_names(names(d_vec), length(d_vec), "d_vec names")
  flci_require_same_name_set(
    names(d_vec), restriction_names, "d_vec", "A_Delta row names"
  )
  d_vec <- d_vec[restriction_names]

  flci_require_names(names(target), length(target), "target names")
  flci_require_same_name_set(names(target), coordinate_names, "target")
  target <- target[coordinate_names]

  W <- flci_finite_numeric_matrix(W, "W")
  Sigma <- flci_finite_numeric_matrix(Sigma, "Sigma")
  A_Delta <- flci_finite_numeric_matrix(A_Delta, "A_Delta")
  normalization <- flci_finite_numeric_matrix(
    normalization, "normalization"
  )
  d_vec_numeric <- flci_finite_numeric_vector(d_vec, "d_vec")
  names(d_vec_numeric) <- restriction_names
  d_vec <- d_vec_numeric
  target_numeric <- flci_finite_numeric_vector(target, "target")
  names(target_numeric) <- coordinate_names
  target <- target_numeric

  if (any(d_vec < 0)) {
    flci_stop("d_vec restriction bounds must all be nonnegative.")
  }
  sigma_asymmetry <- flci_max_abs(Sigma - t(Sigma))
  if (sigma_asymmetry >
      flci_scaled_residual_tol(Sigma, tol)) {
    flci_stop(
      "Sigma is materially asymmetric after name alignment (residual ",
      format(sigma_asymmetry, digits = 6), ")."
    )
  }
  Sigma <- flci_symmetrize(Sigma)

  pre <- grepl("^delta_pre_", coordinate_names)
  post <- grepl("^delta_post_", coordinate_names)
  if (!all(pre | post)) flci_stop("Unclassified coordinate names.")
  if (!any(pre) || !any(post)) {
    flci_stop("The project-special-case reducer requires pre and post blocks.")
  }
  target_pre_residual <- flci_max_abs(target[pre])
  if (target_pre_residual > tol) {
    flci_stop(
      "The project target must have a zero pre-treatment block (residual ",
      format(target_pre_residual, digits = 6), ")."
    )
  }
  target[pre] <- 0
  if (flci_max_abs(target[post]) <= tol) {
    flci_stop("The project target must load on at least one post coordinate.")
  }

  N_minus <- normalization[, pre, drop = FALSE]
  normalization_post_residual <- if (nrow(normalization) == 0L) 0 else
    flci_max_abs(normalization[, post, drop = FALSE])
  if (normalization_post_residual > tol) {
    flci_stop(
      "This project-special-case reducer requires a pre-only normalization ",
      "(post residual ", format(normalization_post_residual, digits = 6), ")."
    )
  }
  E_minus <- flci_null_basis(N_minus, tol = tol)
  rownames(E_minus) <- coordinate_names[pre]
  E <- matrix(
    0,
    nrow = length(coordinate_names),
    ncol = ncol(E_minus) + sum(post)
  )
  if (ncol(E_minus) > 0L) {
    E[pre, seq_len(ncol(E_minus))] <- E_minus
  }
  E[post, ncol(E_minus) + seq_len(sum(post))] <- diag(sum(post))

  normalization_residual <- flci_max_abs(N_minus %*% E_minus)
  basis_orthogonality <- flci_max_abs(
    t(E_minus) %*% E_minus - diag(ncol(E_minus))
  )
  if (normalization_residual > tol || basis_orthogonality > tol) {
    flci_stop(
      "Mechanical normalization reduction failed its numerical checks ",
      "(normalization residual ", format(normalization_residual, digits = 6),
      ", basis residual ", format(basis_orthogonality, digits = 6), ")."
    )
  }

  W_bar <- t(E) %*% W %*% E
  Sigma_bar <- flci_symmetrize(t(E) %*% Sigma %*% E)
  d_pre_bar <- ncol(E_minus)
  pre_bar <- seq_len(d_pre_bar)
  post_bar <- d_pre_bar + seq_len(sum(post))
  W_pre_identity <- flci_max_abs(
    W_bar[pre_bar, pre_bar, drop = FALSE] - diag(d_pre_bar)
  )
  W_pre_post_zero <- flci_max_abs(
    W_bar[pre_bar, post_bar, drop = FALSE]
  )
  if (W_pre_identity > tol || W_pre_post_zero > tol) {
    flci_stop(
      "Unsupported design: reduced W does not have the required ",
      "[[I, 0], [H, C]] block form (identity residual ",
      format(W_pre_identity, digits = 6), ", zero-block residual ",
      format(W_pre_post_zero, digits = 6), ")."
    )
  }
  H <- W_bar[post_bar, pre_bar, drop = FALSE]
  C <- W_bar[post_bar, post_bar, drop = FALSE]

  eigen_sigma <- eigen(
    Sigma_bar, symmetric = TRUE, only.values = TRUE
  )$values
  sigma_eigen_scale <- max(abs(eigen_sigma))
  sigma_pd_tol <- tol * sigma_eigen_scale
  sigma_rank <- sum(eigen_sigma > sigma_pd_tol)
  sigma_invalid <- if (allow_singular_covariance) {
    !is.finite(sigma_eigen_scale) || min(eigen_sigma) < -sigma_pd_tol
  } else {
    !is.finite(sigma_eigen_scale) || sigma_eigen_scale <= 0 ||
      min(eigen_sigma) <= sigma_pd_tol
  }
  if (sigma_invalid) {
    flci_stop(
      if (allow_singular_covariance) {
        "Unsupported design: reduced covariance must be positive semidefinite "
      } else {
        "Unsupported design: reduced covariance must be positive definite "
      },
      "(minimum eigenvalue ", format(min(eigen_sigma), digits = 6),
      ", threshold ", format(sigma_pd_tol, digits = 6), ")."
    )
  }

  if (!any(grepl("_pos$", restriction_names))) {
    flci_stop("A_Delta has no oriented '_pos' restriction rows.")
  }
  sd_blocks <- flci_oriented_sd_blocks(A_Delta, d_vec, tol = tol)
  R_bar <- sd_blocks$R %*% E_minus
  Q_minus_bar <- sd_blocks$Q_minus %*% E_minus
  Q_plus <- sd_blocks$Q_plus
  Q_bar <- cbind(Q_minus_bar, Q_plus)
  ell <- target[post]

  if (nrow(Q_plus) == 0L || nrow(Q_plus) != ncol(Q_plus)) {
    flci_stop(
      "Unsupported design: Q_plus must be nonempty and square; got ",
      paste(dim(Q_plus), collapse = "x"), "."
    )
  }
  rank_Q_plus <- qr(Q_plus, tol = tol, LAPACK = FALSE)$rank
  if (rank_Q_plus != ncol(Q_plus)) {
    flci_stop(
      "Unsupported design: Q_plus must have full rank; rank ",
      rank_Q_plus, " for dimension ", paste(dim(Q_plus), collapse = "x"), "."
    )
  }
  Q_plus_condition <- kappa(Q_plus, exact = TRUE)
  condition_limit <- 1 / tol
  if (!is.finite(Q_plus_condition) || Q_plus_condition > condition_limit) {
    flci_stop(
      "Unsupported design: Q_plus is too ill-conditioned (condition number ",
      format(Q_plus_condition, digits = 6), "; limit ",
      format(condition_limit, digits = 6), ")."
    )
  }

  rank_R <- if (nrow(R_bar) == 0L) 0L else
    qr(R_bar, tol = tol, LAPACK = FALSE)$rank
  if (rank_R != nrow(R_bar)) {
    flci_stop(
      "Unsupported design: quotient R must have full row rank for this ",
      "project-special-case API (rank ", rank_R, " for ", nrow(R_bar),
      " rows). The generic multiplier value-function program is not ",
      "implemented in flci_core.R."
    )
  }

  lambda_rhs <- as.numeric(t(C) %*% ell)
  lambda_0 <- qr.solve(t(Q_plus), lambda_rhs, tol = tol)
  if (length(lambda_0) != nrow(Q_plus) || any(!is.finite(lambda_0))) {
    flci_stop("Q_plus solve returned a nonfinite or malformed multiplier.")
  }
  lambda_residual <- flci_max_abs(
    t(Q_plus) %*% lambda_0 - lambda_rhs
  )
  lambda_residual_tol <- flci_scaled_residual_tol(lambda_rhs, tol)
  if (lambda_residual > lambda_residual_tol) {
    flci_stop(
      "Unsupported design: Q_plus multiplier solve residual is too large ",
      "(", format(lambda_residual, digits = 6), "; tolerance ",
      format(lambda_residual_tol, digits = 6), ")."
    )
  }
  q_0 <- as.numeric(t(Q_minus_bar) %*% lambda_0 - t(H) %*% ell)
  B_Q <- sum(sd_blocks$m_Q * abs(lambda_0))
  if (any(!is.finite(q_0)) || !is.finite(B_Q) || B_Q < 0) {
    flci_stop("The reduced center or boundary-post bias is nonfinite/invalid.")
  }

  if (nrow(R_bar) == 0L) {
    gamma_variance_eigen <- numeric(0)
  } else {
    D_gamma <- rbind(
      t(R_bar),
      matrix(0, nrow = length(ell), ncol = nrow(R_bar))
    )
    # Compare curvature and its construction bound in common covariance
    # units: the uncanceled bound can overflow even when curvature is zero.
    # Only rows on which gamma acts contribute; unrelated large variances
    # must not erase a small but positive curvature block.
    active <- rowSums(abs(D_gamma)) > 0
    D_active <- D_gamma[active, , drop = FALSE]
    block <- flci_scale_covariance_block(Sigma_bar[active, active, drop = FALSE])
    gamma_covariance_scale <- block$scale
    gamma_covariance <- block$covariance
    G_gamma <- flci_symmetrize(
      flci_checked_product(flci_checked_product(t(D_active), gamma_covariance), D_active)
    )
    gamma_variance_eigen <- eigen(
      G_gamma,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    gamma_eigen_scale <- max(abs(gamma_variance_eigen))
    # Cancellation can make D' Sigma D mathematically zero. Its roundoff
    # scale is the absolute construction, not the canceled result itself.
    gamma_construction_scale <- norm(
      flci_checked_product(flci_checked_product(t(abs(D_active)), abs(gamma_covariance)),
                           abs(D_active)), type = "2"
    )
    gamma_pd_tol <- tol * max(gamma_eigen_scale, gamma_construction_scale)
    gamma_rank <- if (gamma_eigen_scale > 0) {
      sum(gamma_variance_eigen > gamma_pd_tol)
    } else {
      0L
    }
    gamma_condition <- if (min(gamma_variance_eigen) > 0) {
      max(gamma_variance_eigen) / min(gamma_variance_eigen)
    } else {
      Inf
    }
    gamma_invalid <- if (allow_singular_covariance) {
      !is.finite(gamma_eigen_scale) ||
        min(gamma_variance_eigen) < -gamma_pd_tol
    } else {
      !is.finite(gamma_eigen_scale) || gamma_eigen_scale <= 0 ||
        min(gamma_variance_eigen) <= gamma_pd_tol ||
        !is.finite(gamma_condition) ||
        gamma_condition > condition_limit
    }
    if (gamma_invalid) {
      flci_stop(
        if (allow_singular_covariance) {
          "Unsupported design: gamma variance curvature must be positive semidefinite "
        } else {
          paste0(
            "Unsupported design: gamma variance curvature must be positive ",
            "definite and well-conditioned "
          )
        },
        "(minimum eigenvalue ",
        format(min(gamma_variance_eigen), digits = 6),
        ", condition number ", format(gamma_condition, digits = 6), ")."
      )
    }
    # Report native covariance units where representable. An overflowing
    # construction bound is also retained as a finite factor pair.
    gamma_construction_scale_normalized <- gamma_construction_scale
    gamma_variance_eigen <- gamma_variance_eigen * gamma_covariance_scale
    gamma_pd_tol <- gamma_pd_tol * gamma_covariance_scale
    gamma_construction_scale <- gamma_construction_scale * gamma_covariance_scale
  }

  colnames(E) <- c(
    paste0("pre_basis_", seq_len(d_pre_bar)),
    coordinate_names[post]
  )
  rownames(E) <- coordinate_names
  names(q_0) <- colnames(E)[pre_bar]
  names(ell) <- coordinate_names[post]
  names(lambda_0) <- rownames(Q_plus)
  m_R <- as.numeric(sd_blocks$m_R)
  names(m_R) <- sd_blocks$pure_pre_rows
  m_Q <- as.numeric(sd_blocks$m_Q)
  names(m_Q) <- sd_blocks$boundary_post_rows

  model <- list(
    implementation_scope = paste(
      "Project-specific quotient-reduced full-row-rank case;",
      "generic multiplier value functions are unsupported."
    ),
    estimator_type = estimator_type,
    coordinate_names = coordinate_names,
    pre_original = coordinate_names[pre],
    post_original = coordinate_names[post],
    E_minus = E_minus,
    E = E,
    W_bar = W_bar,
    Sigma_bar = Sigma_bar,
    H = H,
    C = C,
    R = R_bar,
    Q_minus = Q_minus_bar,
    Q_plus = Q_plus,
    Q = Q_bar,
    m_R = m_R,
    m_Q = m_Q,
    lambda_0 = as.numeric(lambda_0),
    q_0 = q_0,
    ell = ell,
    B_Q = B_Q,
    normalization = normalization,
    N_minus = N_minus,
    target = target,
    variance_clip_rel_tol = 1e-10,
    diagnostics = list(
      normalization_rank = if (nrow(N_minus) == 0L) 0L else
        qr(N_minus, tol = tol)$rank,
      normalization_residual = normalization_residual,
      normalization_post_residual = normalization_post_residual,
      basis_orthogonality = basis_orthogonality,
      W_pre_identity = W_pre_identity,
      W_pre_post_zero = W_pre_post_zero,
      target_pre_residual = target_pre_residual,
      Sigma_asymmetry = sigma_asymmetry,
      Q_plus_rank = rank_Q_plus,
      Q_plus_dimension = paste(dim(Q_plus), collapse = "x"),
      Q_plus_condition = Q_plus_condition,
      lambda_residual = lambda_residual,
      lambda_residual_tolerance = lambda_residual_tol,
      R_rank = rank_R,
      R_nrow = nrow(R_bar),
      R_quotient_full_row_rank = rank_R == nrow(R_bar),
      gamma_variance_min_eigenvalue =
        if (length(gamma_variance_eigen) == 0L) NA_real_ else
          min(gamma_variance_eigen),
      gamma_variance_condition =
        if (length(gamma_variance_eigen) == 0L) NA_real_ else
          max(gamma_variance_eigen) / min(gamma_variance_eigen),
      Sigma_min_eigenvalue = min(eigen_sigma),
      Sigma_max_eigenvalue = max(eigen_sigma),
      Sigma_condition = if (min(eigen_sigma) > 0) {
        max(eigen_sigma) / min(eigen_sigma)
      } else {
        Inf
      },
      Sigma_rank = sigma_rank,
      Sigma_dimension = length(eigen_sigma),
      Sigma_singular = sigma_rank < length(eigen_sigma),
      singular_covariance_allowed = allow_singular_covariance,
      gamma_variance_rank = if (length(gamma_variance_eigen) == 0L) {
        0L
      } else {
        gamma_rank
      },
      gamma_variance_tolerance = if (length(gamma_variance_eigen) == 0L) 0 else gamma_pd_tol,
      gamma_variance_construction_scale = if (length(gamma_variance_eigen) == 0L) 0 else gamma_construction_scale,
      gamma_variance_construction_scale_normalized = if (length(gamma_variance_eigen) == 0L) 0 else gamma_construction_scale_normalized,
      gamma_variance_covariance_scale = if (length(gamma_variance_eigen) == 0L) 1 else gamma_covariance_scale,
      maintained_assumptions_passed = TRUE
    )
  )
  class(model) <- "exploratory_flci_model"
  model
}

flci_model_v <- function(model, gamma = numeric(nrow(model$R))) {
  if (length(gamma) != nrow(model$R) || any(!is.finite(gamma))) {
    flci_stop("gamma must have the correct length and contain finite values.")
  }
  q <- model$q_0 + as.numeric(t(model$R) %*% gamma)
  c(q, model$ell)
}

flci_model_sd <- function(model, gamma = numeric(nrow(model$R))) {
  v <- flci_model_v(model, gamma)
  Sigma_v <- as.numeric(model$Sigma_bar %*% v)
  variance <- sum(v * Sigma_v)
  variance_scale <- sum(abs(v) * abs(Sigma_v))
  flci_safe_sqrt_variance(
    variance,
    variance_scale,
    rel_tol = model$variance_clip_rel_tol,
    label = "model variance"
  )
}

flci_model_bias <- function(model, gamma = numeric(nrow(model$R))) {
  if (length(gamma) != nrow(model$R) || any(!is.finite(gamma))) {
    flci_stop("gamma must have the correct length and contain finite values.")
  }
  model$B_Q + sum(model$m_R * abs(gamma))
}

flci_model_objective <- function(
    model,
    gamma = numeric(nrow(model$R)),
    alpha = 0.05) {
  if (length(gamma) != nrow(model$R) || any(!is.finite(gamma))) {
    return(Inf)
  }
  B <- flci_model_bias(model, gamma)
  s <- flci_model_sd(model, gamma)
  if (!is.finite(B) || !is.finite(s)) return(Inf)
  flci_half_length(B, s, alpha = alpha)
}

flci_support_primal <- function(
    model,
    gamma,
    solver = "GLPK",
    solver_opts = list(),
    feasibility_tol = 1e-6,
    agreement_tol = 1e-6) {
  if (length(feasibility_tol) != 1L || !is.finite(feasibility_tol) ||
      feasibility_tol <= 0 ||
      length(agreement_tol) != 1L || !is.finite(agreement_tol) ||
      agreement_tol <= 0) {
    flci_stop("Support feasibility/agreement tolerances must be positive.")
  }
  v <- flci_model_v(model, gamma)
  r <- as.numeric(t(model$W_bar) %*% v)
  d <- length(r)
  d_pre <- length(model$q_0)

  A_parts <- list()
  rhs_parts <- list()
  if (nrow(model$R) > 0L) {
    A_R <- cbind(
      model$R,
      matrix(0, nrow(model$R), length(model$ell))
    )
    A_parts <- c(A_parts, list(A_R, -A_R))
    rhs_parts <- c(rhs_parts, list(model$m_R, model$m_R))
  }
  A_parts <- c(A_parts, list(model$Q, -model$Q))
  rhs_parts <- c(rhs_parts, list(model$m_Q, model$m_Q))
  A <- do.call(rbind, A_parts)
  rhs <- unlist(rhs_parts, use.names = FALSE)

  Delta <- CVXR::Variable(d)
  problem <- CVXR::Problem(
    CVXR::Maximize(t(r) %*% Delta),
    list(A %*% Delta <= rhs)
  )
  result <- flci_cvxr_solve(problem, solver, solver_opts)
  normalized_status <- flci_assert_cvxr_status(
    result$status, "Support LP"
  )
  solver_value <- as.numeric(result$value)
  Delta_value <- as.numeric(result$getValue(Delta))
  if (length(solver_value) != 1L || !is.finite(solver_value) ||
      length(Delta_value) != d || any(!is.finite(Delta_value))) {
    flci_stop("Support LP returned nonfinite or malformed primal values.")
  }
  lhs <- as.numeric(A %*% Delta_value)
  feasibility_scale <- max(1, flci_max_abs(lhs), flci_max_abs(rhs))
  primal_violation <- max(0, max(lhs - rhs))
  if (!is.finite(primal_violation) ||
      primal_violation > feasibility_tol * feasibility_scale) {
    flci_stop(
      "Support LP violates its restriction budget (violation ",
      format(primal_violation, digits = 6), "; tolerance ",
      format(feasibility_tol * feasibility_scale, digits = 6), ")."
    )
  }
  objective_recomputed <- sum(r * Delta_value)
  objective_residual <- abs(solver_value - objective_recomputed)
  objective_scale <- max(
    1, abs(solver_value), abs(objective_recomputed)
  )
  if (!is.finite(objective_recomputed) ||
      objective_residual > agreement_tol * objective_scale) {
    flci_stop(
      "Support LP solver objective does not match the recomputed objective ",
      "(residual ", format(objective_residual, digits = 6), ")."
    )
  }
  formula <- flci_model_bias(model, gamma)
  formula_residual <- abs(objective_recomputed - formula)
  formula_scale <- max(1, abs(objective_recomputed), abs(formula))
  if (!is.finite(formula) ||
      formula_residual > agreement_tol * formula_scale) {
    flci_stop(
      "Support LP does not agree with the project-special-case support ",
      "formula (residual ", format(formula_residual, digits = 6), ")."
    )
  }
  list(
    status = result$status,
    normalized_status = normalized_status,
    status_acceptable = TRUE,
    value = solver_value,
    objective_recomputed = objective_recomputed,
    objective_residual = objective_residual,
    formula = formula,
    difference = objective_recomputed - formula,
    formula_residual = formula_residual,
    primal_violation = primal_violation,
    feasibility_tolerance = feasibility_tol * feasibility_scale,
    r = r,
    Delta = Delta_value,
    d_pre = d_pre
  )
}

flci_inner_frontier <- function(
    model,
    bias_budget,
    solver = "ECOS",
    solver_opts = list(),
    feasibility_tol = 1e-6,
    objective_tol = 1e-6) {
  if (length(bias_budget) != 1L || !is.finite(bias_budget) ||
      bias_budget < -1e-12) {
    flci_stop("bias_budget must be one finite nonnegative scalar.")
  }
  if (length(feasibility_tol) != 1L || !is.finite(feasibility_tol) ||
      feasibility_tol <= 0 ||
      length(objective_tol) != 1L || !is.finite(objective_tol) ||
      objective_tol <= 0) {
    flci_stop("Inner feasibility/objective tolerances must be positive.")
  }
  flci_validate_solver_opts(solver_opts)
  bias_budget <- max(0, bias_budget)
  k <- nrow(model$R)
  if (k == 0L) {
    gamma <- numeric(0)
    sd_value <- flci_model_sd(model, gamma)
    return(list(
      status = "analytic",
      normalized_status = "analytic",
      status_acceptable = TRUE,
      gamma = gamma,
      q = model$q_0,
      sd = sd_value,
      variance = sd_value^2,
      solver_objective = sd_value^2,
      objective_recomputed = sd_value^2,
      objective_residual = 0,
      achieved_bias_budget = 0,
      requested_bias_budget = bias_budget,
      budget_violation = 0,
      budget_tolerance = feasibility_tol * max(1, bias_budget)
    ))
  }

  A_gamma <- rbind(
    t(model$R),
    matrix(0, nrow = length(model$ell), ncol = k)
  )
  v_0 <- c(model$q_0, model$ell)
  G <- flci_symmetrize(t(A_gamma) %*% model$Sigma_bar %*% A_gamma)
  g <- as.numeric(t(A_gamma) %*% model$Sigma_bar %*% v_0)
  c_0 <- as.numeric(t(v_0) %*% model$Sigma_bar %*% v_0)

  gamma_var <- CVXR::Variable(k)
  variance_expression <-
    CVXR::quad_form(gamma_var, G) +
    2 * t(g) %*% gamma_var +
    c_0
  budget_expression <- CVXR::sum_entries(
    model$m_R * abs(gamma_var)
  )
  problem <- CVXR::Problem(
    CVXR::Minimize(variance_expression),
    list(budget_expression <= bias_budget)
  )
  result <- flci_cvxr_solve(problem, solver, solver_opts)
  normalized_status <- flci_assert_cvxr_status(
    result$status, "Inner frontier problem"
  )
  gamma <- as.numeric(result$getValue(gamma_var))
  if (length(gamma) != k || any(!is.finite(gamma))) {
    flci_stop("Inner frontier returned a nonfinite or malformed gamma.")
  }
  q <- model$q_0 + as.numeric(t(model$R) %*% gamma)
  if (any(!is.finite(q))) {
    flci_stop("Inner frontier returned a nonfinite reduced pre center.")
  }
  achieved_budget <- sum(model$m_R * abs(gamma))
  budget_tolerance <- feasibility_tol * max(1, bias_budget)
  budget_violation <- max(0, achieved_budget - bias_budget)
  if (!is.finite(achieved_budget) ||
      budget_violation > budget_tolerance) {
    flci_stop(
      "Inner frontier violates its bias budget (violation ",
      format(budget_violation, digits = 6), "; tolerance ",
      format(budget_tolerance, digits = 6), ")."
    )
  }
  v <- c(q, model$ell)
  Sigma_v <- as.numeric(model$Sigma_bar %*% v)
  variance_recomputed <- sum(v * Sigma_v)
  variance_scale <- sum(abs(v) * abs(Sigma_v))
  sd_value <- flci_safe_sqrt_variance(
    variance_recomputed,
    variance_scale,
    rel_tol = model$variance_clip_rel_tol,
    label = "inner-frontier variance"
  )
  solver_objective <- as.numeric(result$value)
  objective_residual <- abs(solver_objective - variance_recomputed)
  objective_scale <- max(
    1, abs(solver_objective), abs(variance_recomputed)
  )
  if (length(solver_objective) != 1L || !is.finite(solver_objective) ||
      !is.finite(variance_recomputed) ||
      objective_residual > objective_tol * objective_scale) {
    flci_stop(
      "Inner-frontier solver objective does not match the independently ",
      "recomputed variance (residual ",
      format(objective_residual, digits = 6), ")."
    )
  }
  list(
    status = result$status,
    normalized_status = normalized_status,
    status_acceptable = TRUE,
    gamma = gamma,
    q = q,
    sd = sd_value,
    variance = sd_value^2,
    solver_objective = solver_objective,
    objective_recomputed = variance_recomputed,
    objective_residual = objective_residual,
    achieved_bias_budget = achieved_budget,
    requested_bias_budget = bias_budget,
    budget_violation = budget_violation,
    budget_tolerance = budget_tolerance,
    G_min_eigenvalue = min(eigen(G, symmetric = TRUE,
                                 only.values = TRUE)$values)
  )
}

flci_nested_solve <- function(
    model,
    alpha = 0.05,
    solver = "ECOS",
    solver_opts = list(),
    outer_tol = 1e-7,
    grid_size = 15L,
    feasibility_tol = 1e-6,
    objective_tol = 1e-6) {
  if (length(outer_tol) != 1L || !is.finite(outer_tol) || outer_tol <= 0 ||
      length(grid_size) != 1L || !is.finite(grid_size) || grid_size < 3) {
    flci_stop("outer_tol must be positive and grid_size must be at least 3.")
  }
  gamma_zero <- numeric(nrow(model$R))
  post_sd <- flci_model_sd(model, gamma_zero)
  chi_post <- flci_half_length(model$B_Q, post_sd, alpha = alpha)
  upper <- max(0, chi_post - model$B_Q)
  cache <- new.env(parent = emptyenv())

  evaluate <- function(b) {
    b <- min(max(as.numeric(b), 0), upper)
    key <- sprintf("%.13g", b)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    inner <- flci_inner_frontier(
      model,
      b,
      solver = solver,
      solver_opts = solver_opts,
      feasibility_tol = feasibility_tol,
      objective_tol = objective_tol
    )
    if (!is.finite(inner$sd)) {
      ans <- list(
        budget = b,
        outer_value = Inf,
        exact_value = Inf,
        exact_value_recomputed = Inf,
        exact_objective_residual = Inf,
        inner = inner
      )
    } else {
      exact_value <- flci_model_objective(
        model,
        inner$gamma,
        alpha = alpha
      )
      exact_value_recomputed <- flci_half_length(
        flci_model_bias(model, inner$gamma),
        inner$sd,
        alpha = alpha
      )
      ans <- list(
        budget = b,
        outer_value = flci_half_length(
          model$B_Q + b,
          inner$sd,
          alpha = alpha
        ),
        exact_value = exact_value,
        exact_value_recomputed = exact_value_recomputed,
        exact_objective_residual = abs(
          exact_value - exact_value_recomputed
        ),
        inner = inner
      )
    }
    assign(key, ans, envir = cache)
    ans
  }

  if (upper <= outer_tol) {
    candidate_b <- 0
    optimize_result <- NULL
  } else {
    optimize_result <- stats::optimize(
      function(b) evaluate(b)$outer_value,
      interval = c(0, upper),
      tol = outer_tol
    )
    candidate_b <- unique(c(
      0,
      upper,
      optimize_result$minimum,
      seq(0, upper, length.out = max(3L, as.integer(grid_size)))
    ))
  }
  evaluated <- lapply(candidate_b, evaluate)

  achieved <- vapply(
    evaluated,
    function(x) x$inner$achieved_bias_budget,
    numeric(1)
  )
  achieved <- achieved[is.finite(achieved)]
  if (length(achieved) > 0L) {
    evaluated <- c(evaluated, lapply(unique(achieved), evaluate))
  }

  exact_values <- vapply(evaluated, `[[`, numeric(1), "exact_value")
  exact_values_recomputed <- vapply(
    evaluated, `[[`, numeric(1), "exact_value_recomputed"
  )
  exact_objective_residuals <- vapply(
    evaluated, `[[`, numeric(1), "exact_objective_residual"
  )
  outer_values <- vapply(evaluated, `[[`, numeric(1), "outer_value")
  best_exact <- evaluated[[which.min(exact_values)]]
  best_outer <- evaluated[[which.min(outer_values)]]

  frontier <- data.frame(
    requested_budget = vapply(evaluated, `[[`, numeric(1), "budget"),
    achieved_budget = vapply(
      evaluated,
      function(x) x$inner$achieved_bias_budget,
      numeric(1)
    ),
    sd = vapply(evaluated, function(x) x$inner$sd, numeric(1)),
    outer_value = outer_values,
    exact_value = exact_values,
    exact_value_recomputed = exact_values_recomputed,
    exact_objective_residual = exact_objective_residuals,
    inner_objective_residual = vapply(
      evaluated,
      function(x) x$inner$objective_residual,
      numeric(1)
    ),
    budget_violation = vapply(
      evaluated,
      function(x) x$inner$budget_violation,
      numeric(1)
    ),
    budget_tolerance = vapply(
      evaluated,
      function(x) x$inner$budget_tolerance,
      numeric(1)
    ),
    status = vapply(
      evaluated,
      function(x) as.character(x$inner$status),
      character(1)
    ),
    status_acceptable = vapply(
      evaluated,
      function(x) isTRUE(x$inner$status_acceptable),
      logical(1)
    ),
    stringsAsFactors = FALSE
  )
  frontier <- frontier[order(frontier$requested_budget), , drop = FALSE]
  rownames(frontier) <- NULL

  list(
    gamma = best_exact$inner$gamma,
    q = best_exact$inner$q,
    sd = best_exact$inner$sd,
    bias = flci_model_bias(model, best_exact$inner$gamma),
    half_length = best_exact$exact_value,
    post_sd = post_sd,
    post_half_length = chi_post,
    upper_bias_budget = upper,
    best_exact_budget = best_exact$budget,
    best_outer_budget = best_outer$budget,
    best_outer_value = best_outer$outer_value,
    inner_status = best_exact$inner$status,
    inner_status_acceptable = best_exact$inner$status_acceptable,
    budget_violation = best_exact$inner$budget_violation,
    budget_tolerance = best_exact$inner$budget_tolerance,
    inner_objective_residual = best_exact$inner$objective_residual,
    exact_objective_residual = best_exact$exact_objective_residual,
    optimize_result = optimize_result,
    frontier = frontier
  )
}

# Smooth objective derivatives on one fixed gamma orthant.  The objective is
# nonsmooth only when a gamma coordinate crosses zero; holding `orthant`
# fixed makes its gradient and Hessian analytic.  Keeping this calculation in
# one helper lets the direct optimizer, KKT certificate, Newton polish, and
# finite-difference tests use exactly the same derivatives.
flci_direct_gamma_parts <- function(
    model,
    gamma,
    alpha = 0.05,
    orthant = NULL,
    hessian = FALSE) {
  k <- nrow(model$R)
  if (length(gamma) != k || any(!is.finite(gamma))) {
    flci_stop("gamma must have the correct length and contain finite values.")
  }
  gamma <- as.numeric(gamma)
  if (is.null(orthant)) orthant <- sign(gamma)
  if (length(orthant) != k || any(!is.finite(orthant)) ||
      any(!(orthant %in% c(-1, 0, 1)))) {
    flci_stop("orthant must contain one value in {-1, 0, 1} per gamma.")
  }
  orthant <- as.numeric(orthant)
  nonzero <- gamma != 0
  if (any(nonzero & orthant != sign(gamma))) {
    flci_stop("orthant must agree with every nonzero gamma coordinate.")
  }
  if (length(hessian) != 1L || is.na(hessian)) {
    flci_stop("hessian must be TRUE or FALSE.")
  }

  D_gamma <- rbind(
    t(model$R),
    matrix(0, nrow = length(model$ell), ncol = k)
  )
  v <- flci_model_v(model, gamma)
  B <- model$B_Q + sum(model$m_R * abs(gamma))
  Sigma_v <- as.numeric(model$Sigma_bar %*% v)
  variance <- sum(v * Sigma_v)
  variance_scale <- sum(abs(v) * abs(Sigma_v))
  s <- flci_safe_sqrt_variance(
    variance,
    variance_scale,
    rel_tol = model$variance_clip_rel_tol,
    label = "direct gamma-objective variance"
  )
  if (s <= 0) {
    flci_stop(
      "Direct gamma derivatives encountered zero variance despite the ",
      "positive-definite reduced covariance and nonzero target."
    )
  }
  r <- B / s
  c_value <- flci_folded_cv(r, alpha = alpha)
  c_prime <- tanh(r * c_value)
  c_second <- (1 - c_prime^2) * (c_value + r * c_prime)
  H_B <- c_prime
  H_s <- c_value - r * c_prime
  ds_dgamma <- as.numeric(t(D_gamma) %*% model$Sigma_bar %*% v) / s
  variance_gradient <- H_s * ds_dgamma
  bias_marginal <- H_B * model$m_R
  bias_gradient <- model$m_R * orthant
  gradient <- variance_gradient + H_B * bias_gradient

  result <- list(
    value = s * c_value,
    bias = B,
    sd = s,
    ratio = r,
    critical_value = c_value,
    critical_value_first = c_prime,
    critical_value_second = c_second,
    H_B = H_B,
    H_s = H_s,
    ds_dgamma = ds_dgamma,
    variance_gradient = variance_gradient,
    bias_marginal = bias_marginal,
    component_scale = abs(variance_gradient) + abs(bias_marginal),
    gradient = gradient,
    orthant = orthant,
    D_gamma = D_gamma,
    v = v
  )
  if (isTRUE(hessian)) {
    G <- flci_symmetrize(t(D_gamma) %*% model$Sigma_bar %*% D_gamma)
    Hessian_s <- G / s - tcrossprod(ds_dgamma) / s
    result$hessian <- flci_symmetrize(
      (c_second / s) * tcrossprod(bias_gradient) -
        (r * c_second / s) *
          (tcrossprod(bias_gradient, ds_dgamma) +
             tcrossprod(ds_dgamma, bias_gradient)) +
        (r^2 * c_second / s) * tcrossprod(ds_dgamma) +
        H_s * Hessian_s
    )
  }
  result
}

# Full nonsmooth KKT certificate in gamma coordinates.  For a_j equal to the
# variance contribution and b_j equal to the nonnegative bias contribution,
# stationarity is a_j +/- b_j = 0 off zero and |a_j| <= b_j at zero.  The
# coordinatewise relative residual is invariant to outcome-unit scaling and
# to a rescaling of an individual gamma coordinate.
flci_direct_kkt <- function(
    model,
    gamma,
    alpha = 0.05,
    active_tol = NULL) {
  parts <- flci_direct_gamma_parts(
    model,
    gamma,
    alpha = alpha,
    orthant = sign(gamma),
    hessian = FALSE
  )
  gamma <- as.numeric(gamma)
  if (is.null(active_tol)) {
    active_tol <- sqrt(.Machine$double.eps) *
      max(1, flci_max_abs(gamma))
  }
  if (length(active_tol) != 1L || !is.finite(active_tol) || active_tol < 0) {
    flci_stop("active_tol must be one finite nonnegative scalar.")
  }
  active <- abs(gamma) > active_tol
  a <- parts$variance_gradient
  b <- parts$bias_marginal
  components <- ifelse(
    active,
    abs(a + b * sign(gamma)),
    pmax(abs(a) - b, 0)
  )
  component_scale <- parts$component_scale
  relative_components <- numeric(length(components))
  positive_scale <- component_scale > 0
  relative_components[positive_scale] <-
    components[positive_scale] / component_scale[positive_scale]
  relative_components[!positive_scale & components > 0] <- Inf

  # Preserve the historical split-gradient scaling as a compatibility field.
  # It is not genuinely unit-free when every split gradient is below one.
  split_gradient <- c(b + a, b - a)
  residual <- flci_max_abs(components)
  list(
    kkt_residual = residual,
    kkt_residual_scaled = residual /
      max(1, flci_max_abs(split_gradient)),
    kkt_residual_relative = flci_max_abs(relative_components),
    components = components,
    relative_components = relative_components,
    component_scale = component_scale,
    active = active,
    active_tol = active_tol,
    variance_gradient = a,
    bias_marginal = b,
    split_gradient = split_gradient,
    objective = parts$value
  )
}

# Fail-closed active-orthant Newton refinement.  The input candidate is never
# replaced unless the polished point is finite, non-increasing in objective,
# and satisfies a high-accuracy full nonsmooth KKT certificate.  Coordinates
# at zero are activated only when their subgradient inequality is violated;
# coordinates reaching an orthant boundary are dropped and re-certified.
flci_direct_polish <- function(
    model,
    gamma,
    alpha = 0.05,
    maxit = 20L,
    relative_tol = 1e-9,
    absolute_tol = 1e-12,
    line_search_maxit = 50L) {
  k <- nrow(model$R)
  if (length(gamma) != k || any(!is.finite(gamma))) {
    flci_stop("gamma must have the correct length and contain finite values.")
  }
  integer_args <- c(maxit, line_search_maxit)
  if (any(!is.finite(integer_args)) || any(integer_args < 1) ||
      any(integer_args != as.integer(integer_args))) {
    flci_stop("Newton and line-search iteration limits must be positive integers.")
  }
  if (length(relative_tol) != 1L || !is.finite(relative_tol) ||
      relative_tol <= 0 || length(absolute_tol) != 1L ||
      !is.finite(absolute_tol) || absolute_tol <= 0) {
    flci_stop("Newton KKT tolerances must be finite and positive.")
  }
  gamma_names <- names(gamma)
  gamma <- as.numeric(gamma)
  names(gamma) <- gamma_names
  input_gamma <- gamma
  input_certificate <- flci_direct_kkt(
    model,
    input_gamma,
    alpha = alpha
  )
  input_value <- input_certificate$objective
  # Canonicalize numerical zeros before constructing the initial orthant.
  # Otherwise a coordinate can be classified inactive by the KKT tolerance
  # while remaining algebraically nonzero, which is inconsistent with a zero
  # orthant entry in the smooth derivative helper.
  initial_active_tolerance <- sqrt(.Machine$double.eps) *
    max(1, flci_max_abs(gamma))
  gamma[abs(gamma) <= initial_active_tolerance] <- 0
  working_initial_gamma <- gamma
  working_initial_certificate <- flci_direct_kkt(
    model,
    gamma,
    alpha = alpha
  )
  working_initial_value <- working_initial_certificate$objective
  objective_tolerance <- 100 * .Machine$double.eps * max(1, abs(input_value))

  converged <- function(certificate) {
    is.finite(certificate$kkt_residual) &&
      is.finite(certificate$kkt_residual_relative) &&
      (certificate$kkt_residual_relative <= relative_tol ||
         certificate$kkt_residual <=
           absolute_tol * max(1, abs(certificate$objective)))
  }
  canonicalization_nonincreasing <-
    working_initial_value <= input_value + objective_tolerance
  if (converged(working_initial_certificate) &&
      canonicalization_nonincreasing) {
    return(list(
      gamma = working_initial_gamma,
      success = TRUE,
      accepted = TRUE,
      status = "already_stationary",
      iterations = 0L,
      activations = 0L,
      drops = 0L,
      objective_initial = input_value,
      objective_working_initial = working_initial_value,
      objective_final = working_initial_value,
      objective_change = working_initial_value - input_value,
      certificate_initial = input_certificate,
      certificate_working_initial = working_initial_certificate,
      certificate = working_initial_certificate
    ))
  }

  active <- working_initial_certificate$active
  orthant <- ifelse(active, sign(gamma), 0)
  activations <- 0L
  drops <- 0L
  status <- "maximum_iterations"
  iterations <- 0L
  for (iteration in seq_len(as.integer(maxit))) {
    iterations <- iteration
    certificate <- flci_direct_kkt(model, gamma, alpha = alpha)
    if (converged(certificate)) {
      status <- "converged"
      break
    }

    inactive <- which(orthant == 0)
    if (length(inactive)) {
      inactive_components <- certificate$components[inactive]
      inactive_scales <- certificate$component_scale[inactive]
      inactive_scores <- ifelse(
        inactive_scales > 0,
        inactive_components / inactive_scales,
        inactive_components / max(1, abs(certificate$objective))
      )
      eligible <- inactive_components >
        absolute_tol * max(1, abs(certificate$objective)) &
        inactive_scores > relative_tol
      if (any(eligible)) {
        candidates <- inactive[eligible]
        selected <- candidates[
          which.max(certificate$components[candidates])
        ]
        direction_sign <- -sign(certificate$variance_gradient[selected])
        if (direction_sign == 0) {
          status <- "inactive_direction_failure"
          break
        }
        orthant[selected] <- direction_sign
        gamma[selected] <- 0
        activations <- activations + 1L
      }
    }

    active_index <- which(orthant != 0)
    if (!length(active_index)) {
      status <- "empty_active_set_not_stationary"
      break
    }
    parts <- flci_direct_gamma_parts(
      model,
      gamma,
      alpha = alpha,
      orthant = orthant,
      hessian = TRUE
    )
    q <- orthant[active_index]
    eta <- q * gamma[active_index]
    gradient_eta <- q * parts$gradient[active_index]
    Hessian_eta <- parts$hessian[active_index, active_index, drop = FALSE]
    Hessian_eta <- flci_symmetrize(Hessian_eta * tcrossprod(q))
    Hessian_scale <- max(.Machine$double.xmin, flci_max_abs(Hessian_eta))
    minimum_eigenvalue <- min(eigen(
      Hessian_eta,
      symmetric = TRUE,
      only.values = TRUE
    )$values)
    eigen_floor <- sqrt(.Machine$double.eps) * Hessian_scale
    diagonal_shift <- max(0, eigen_floor - minimum_eigenvalue)
    Newton_matrix <- Hessian_eta +
      diag(diagonal_shift, nrow(Hessian_eta))
    step_eta <- tryCatch(
      as.numeric(solve(Newton_matrix, -gradient_eta)),
      error = function(condition) rep(NA_real_, length(active_index))
    )
    slope <- if (all(is.finite(step_eta))) {
      sum(gradient_eta * step_eta)
    } else {
      NA_real_
    }
    if (!is.finite(slope) || slope >= 0) {
      step_eta <- -gradient_eta
      slope <- -sum(gradient_eta^2)
    }
    if (!is.finite(slope) || slope >= 0 || !any(step_eta != 0)) {
      status <- "descent_direction_failure"
      break
    }

    maximum_step <- 1
    toward_boundary <- step_eta < 0
    if (any(toward_boundary)) {
      maximum_step <- min(
        maximum_step,
        min(-eta[toward_boundary] / step_eta[toward_boundary])
      )
    }
    if (!is.finite(maximum_step) || maximum_step <= 0) {
      boundary <- active_index[eta <=
        sqrt(.Machine$double.eps) * max(1, flci_max_abs(gamma))]
      if (!length(boundary)) {
        status <- "orthant_boundary_failure"
        break
      }
      gamma[boundary] <- 0
      orthant[boundary] <- 0
      drops <- drops + length(boundary)
      next
    }

    step_length <- maximum_step
    accepted_step <- FALSE
    candidate <- gamma
    for (unused in seq_len(as.integer(line_search_maxit))) {
      candidate_eta <- pmax(eta + step_length * step_eta, 0)
      candidate <- gamma
      candidate[active_index] <- q * candidate_eta
      candidate_value <- flci_model_objective(model, candidate, alpha = alpha)
      if (is.finite(candidate_value) &&
          candidate_value <= parts$value +
            1e-4 * step_length * slope + objective_tolerance) {
        accepted_step <- TRUE
        break
      }
      step_length <- step_length / 2
    }
    if (!accepted_step) {
      status <- "line_search_failure"
      break
    }
    gamma <- candidate
    active_tolerance <- sqrt(.Machine$double.eps) *
      max(1, flci_max_abs(gamma))
    boundary <- active_index[abs(gamma[active_index]) <= active_tolerance]
    if (length(boundary)) {
      gamma[boundary] <- 0
      orthant[boundary] <- 0
      drops <- drops + length(boundary)
    }
  }

  final_certificate <- flci_direct_kkt(model, gamma, alpha = alpha)
  if (converged(final_certificate)) status <- "converged"
  objective_nonincreasing <- is.finite(final_certificate$objective) &&
    final_certificate$objective <= input_value + objective_tolerance
  success <- identical(status, "converged") && objective_nonincreasing
  if (!success) {
    gamma <- input_gamma
    final_certificate <- input_certificate
  }
  list(
    gamma = gamma,
    success = success,
    accepted = success,
    status = status,
    iterations = iterations,
    activations = activations,
    drops = drops,
    objective_initial = input_value,
    objective_working_initial = working_initial_value,
    objective_final = final_certificate$objective,
    objective_change = final_certificate$objective - input_value,
    certificate_initial = input_certificate,
    certificate_working_initial = working_initial_certificate,
    certificate = final_certificate
  )
}

flci_direct_solve <- function(
    model,
    alpha = 0.05,
    seed = 1234,
    n_random_starts = 8L,
    maxit = 20000L) {
  k <- nrow(model$R)
  if (k == 0L) {
    gamma <- numeric(0)
    objective_value <- flci_model_objective(model, gamma, alpha)
    return(list(
      gamma = gamma,
      half_length = objective_value,
      sd = flci_model_sd(model, gamma),
      bias = flci_model_bias(model, gamma),
      convergence = 0L,
      optimizer_convergence = 0L,
      selected_source = "analytic",
      selected_start = NA_integer_,
      objective_recomputed = objective_value,
      objective_residual = 0,
      kkt_residual = 0,
      kkt_residual_scaled = 0,
      kkt_residual_relative = 0,
      pre_polish_kkt_residual = 0,
      pre_polish_kkt_residual_relative = 0,
      polish_ok = TRUE,
      polish_status = "analytic",
      polish_iterations = 0L,
      polish_attempted_count = 0L,
      polish_success_count = 0L,
      split_complementarity_residual = 0,
      canonical_split = numeric(0),
      starts = data.frame()
    ))
  }

  chi_post <- flci_model_objective(model, numeric(k), alpha)
  bias_room <- max(chi_post - model$B_Q, 0)
  scale <- rep(1, k)
  positive <- model$m_R > 0
  scale[positive] <- pmax(
    bias_room / pmax(model$m_R[positive], .Machine$double.eps),
    0.25
  )

  set.seed(seed)
  gamma_starts <- list(rep(0, k))
  if (n_random_starts > 0L) {
    for (j in seq_len(n_random_starts)) {
      gamma_starts[[length(gamma_starts) + 1L]] <-
        stats::runif(k, -1, 1) * scale
    }
  }

  A_gamma <- rbind(
    t(model$R),
    matrix(0, nrow = length(model$ell), ncol = k)
  )
  split_to_gamma <- function(x) x[seq_len(k)] - x[k + seq_len(k)]
  gamma_to_split <- function(gamma) c(pmax(gamma, 0), pmax(-gamma, 0))
  split_bias_weights <- rep(model$m_R, 2L)

  objective_parts <- function(x) {
    if (length(x) != 2L * k || any(!is.finite(x)) ||
        any(x < -sqrt(.Machine$double.eps))) {
      return(list(
        value = .Machine$double.xmax / 100,
        gradient = rep(0, 2L * k)
      ))
    }
    gamma <- split_to_gamma(x)
    v <- flci_model_v(model, gamma)
    B <- model$B_Q + sum(split_bias_weights * x)
    Sigma_v <- as.numeric(model$Sigma_bar %*% v)
    variance <- sum(v * Sigma_v)
    variance_scale <- sum(abs(v) * abs(Sigma_v))
    if (!is.finite(B) || !is.finite(variance) ||
        !is.finite(variance_scale)) {
      return(list(
        value = .Machine$double.xmax / 100,
        gradient = rep(0, 2L * k)
      ))
    }
    s <- flci_safe_sqrt_variance(
      variance,
      variance_scale,
      rel_tol = model$variance_clip_rel_tol,
      label = "direct-objective variance"
    )
    if (s <= 0) {
      flci_stop(
        "Direct optimization encountered zero variance despite the ",
        "positive-definite reduced covariance and nonzero target."
      )
    }
    r <- B / s
    c_value <- flci_folded_cv(r, alpha = alpha)
    c_prime <- tanh(r * c_value)
    H_B <- c_prime
    H_s <- c_value - r * c_prime
    ds_dgamma <- as.numeric(
      t(A_gamma) %*% model$Sigma_bar %*% v
    ) / s
    list(
      value = s * c_value,
      gradient = c(
        H_B * model$m_R + H_s * ds_dgamma,
        H_B * model$m_R - H_s * ds_dgamma
      )
    )
  }

  objective <- function(x) objective_parts(x)$value
  gradient <- function(x) objective_parts(x)$gradient
  split_starts <- lapply(gamma_starts, gamma_to_split)
  fits <- lapply(split_starts, function(start) {
    suppressWarnings(stats::optim(
      par = start,
      fn = objective,
      gr = gradient,
      method = "L-BFGS-B",
      lower = rep(0, 2L * k),
      control = list(
        maxit = maxit,
        factr = 1e5,
        pgtol = 1e-10
      )
    ))
  })
  fit_values <- vapply(fits, function(x) x$value, numeric(1))
  start_values <- vapply(split_starts, objective, numeric(1))
  canonical_fit_x <- lapply(
    fits,
    function(fit) gamma_to_split(split_to_gamma(fit$par))
  )
  canonical_fit_values <- vapply(canonical_fit_x, objective, numeric(1))
  use_fit <- is.finite(canonical_fit_values) &
    canonical_fit_values <= start_values
  candidate_x <- lapply(seq_along(fits), function(index) {
    if (use_fit[index]) canonical_fit_x[[index]] else split_starts[[index]]
  })
  candidate_gamma <- lapply(candidate_x, split_to_gamma)
  pre_polish_values <- vapply(candidate_x, objective, numeric(1))
  pre_polish_values[!is.finite(pre_polish_values)] <- Inf
  if (!any(is.finite(pre_polish_values))) {
    flci_stop("Direct optimization produced no finite candidate objective.")
  }
  minimum_pre_polish <- min(pre_polish_values)
  near_tie_tolerance <- 1e-8 * max(1, abs(minimum_pre_polish))
  near_tie <- is.finite(pre_polish_values) &
    pre_polish_values <= minimum_pre_polish + near_tie_tolerance
  polishes <- vector("list", length(candidate_gamma))
  for (index in which(near_tie)) {
    polishes[[index]] <- tryCatch(
      flci_direct_polish(
        model,
        candidate_gamma[[index]],
        alpha = alpha
      ),
      error = function(condition) {
        list(
          gamma = candidate_gamma[[index]],
          success = FALSE,
          accepted = FALSE,
          status = paste0("error: ", conditionMessage(condition)),
          iterations = NA_integer_,
          objective_final = pre_polish_values[index],
          certificate_initial = NULL,
          certificate = NULL
        )
      }
    )
  }
  polish_success <- vapply(
    polishes,
    function(x) !is.null(x) && is.list(x) && isTRUE(x$success),
    logical(1)
  )
  polished_values <- pre_polish_values
  for (index in which(polish_success)) {
    polished_values[index] <- polishes[[index]]$objective_final
  }
  eligible <- near_tie & polish_success & is.finite(polished_values)
  best_index <- if (any(eligible)) {
    which(eligible)[which.min(polished_values[eligible])]
  } else {
    which.min(pre_polish_values)
  }
  selected_source <- if (use_fit[best_index]) {
    "optimizer"
  } else {
    "initial_point"
  }
  selected_polish <- polishes[[best_index]]
  polish_ok <- is.list(selected_polish) && isTRUE(selected_polish$success)
  best_gamma <- if (polish_ok) {
    selected_polish$gamma
  } else {
    candidate_gamma[[best_index]]
  }
  best_x <- gamma_to_split(best_gamma)
  final_objective <- flci_model_objective(model, best_gamma, alpha)
  objective_residual <- abs(
    final_objective - if (polish_ok) {
      polished_values[best_index]
    } else {
      pre_polish_values[best_index]
    }
  )
  if (!is.finite(final_objective) ||
      objective_residual >
      1e-8 * max(
        1,
        abs(final_objective),
        if (polish_ok) polished_values[best_index]
        else pre_polish_values[best_index]
      )) {
    flci_stop(
      "Direct candidate objective does not match its independently ",
      "recomputed half-length (residual ",
      format(objective_residual, digits = 6), ")."
    )
  }
  final_certificate <- flci_direct_kkt(model, best_gamma, alpha = alpha)
  pre_polish_certificate <- if (
    is.list(selected_polish) &&
      is.list(selected_polish$certificate_initial)
  ) {
    selected_polish$certificate_initial
  } else {
    flci_direct_kkt(
      model,
      candidate_gamma[[best_index]],
      alpha = alpha
    )
  }
  kkt_residual <- final_certificate$kkt_residual
  kkt_residual_scaled <- final_certificate$kkt_residual_scaled
  kkt_residual_relative <- final_certificate$kkt_residual_relative
  split_complementarity_residual <- flci_max_abs(
    pmin(best_x[seq_len(k)], best_x[k + seq_len(k)])
  )

  start_table <- data.frame(
    start = seq_along(fits),
    initial_value = start_values,
    optimizer_value = fit_values,
    canonical_optimizer_value = canonical_fit_values,
    pre_polish_value = pre_polish_values,
    value = ifelse(polish_success, polished_values, pre_polish_values),
    selected_source = ifelse(use_fit, "optimizer", "initial_point"),
    selected = seq_along(fits) == best_index,
    near_tie = near_tie,
    polish_success = polish_success,
    polish_status = vapply(
      polishes,
      function(x) if (!is.null(x) && is.list(x)) {
        as.character(x$status)
      } else {
        "not_attempted"
      },
      character(1)
    ),
    polish_iterations = vapply(
      polishes,
      function(x) if (!is.null(x) && is.list(x)) {
        as.integer(x$iterations)
      } else {
        NA_integer_
      },
      integer(1)
    ),
    convergence = vapply(fits, function(x) x$convergence, integer(1)),
    iterations = vapply(fits, function(x) x$counts[["function"]], numeric(1))
  )
  list(
    gamma = best_gamma,
    half_length = final_objective,
    sd = flci_model_sd(model, best_gamma),
    bias = flci_model_bias(model, best_gamma),
    convergence = if (selected_source == "optimizer") {
      fits[[best_index]]$convergence
    } else {
      NA_integer_
    },
    optimizer_convergence = fits[[best_index]]$convergence,
    selected_source = selected_source,
    selected_start = best_index,
    objective_recomputed = final_objective,
    objective_residual = objective_residual,
    kkt_residual = kkt_residual,
    kkt_residual_scaled = kkt_residual_scaled,
    kkt_residual_relative = kkt_residual_relative,
    pre_polish_kkt_residual = pre_polish_certificate$kkt_residual,
    pre_polish_kkt_residual_relative =
      pre_polish_certificate$kkt_residual_relative,
    polish_ok = polish_ok,
    polish_status = if (!is.null(selected_polish) && is.list(selected_polish)) {
      as.character(selected_polish$status)
    } else {
      "not_attempted"
    },
    polish_iterations = if (!is.null(selected_polish) && is.list(selected_polish)) {
      as.integer(selected_polish$iterations)
    } else {
      NA_integer_
    },
    polish_attempted_count = sum(near_tie),
    polish_success_count = sum(polish_success),
    split_complementarity_residual = split_complementarity_residual,
    canonical_split = best_x,
    near_tie_tolerance = near_tie_tolerance,
    starts = start_table,
    all_fits = fits,
    polishes = polishes
  )
}

flci_compare_solvers <- function(model, alpha = 0.05) {
  nested <- flci_nested_solve(model, alpha = alpha)
  direct <- flci_direct_solve(model, alpha = alpha)
  direct_at_nested <- flci_model_objective(
    model,
    nested$gamma,
    alpha = alpha
  )
  nested_at_direct <- flci_model_objective(
    model,
    direct$gamma,
    alpha = alpha
  )

  list(
    nested = nested,
    direct = direct,
    summary = data.frame(
      post_half_length = nested$post_half_length,
      nested_half_length = nested$half_length,
      direct_half_length = direct$half_length,
      direct_at_nested = direct_at_nested,
      nested_at_direct = nested_at_direct,
      nested_minus_direct = nested$half_length - direct$half_length,
      full_minus_post = nested$half_length - nested$post_half_length,
      nested_sd = nested$sd,
      direct_sd = direct$sd,
      nested_bias = nested$bias,
      direct_bias = direct$bias
    )
  )
}
