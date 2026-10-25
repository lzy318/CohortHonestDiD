.hc_fit_diagnostics <- function(prepared, info, covariance, normalization, beta) {
  panel <- prepared$panel
  units <- prepared$unit_registry
  groups <- c(0L, prepared$cohort_times)
  counts <- do.call(rbind, lapply(groups, function(g) {
    uu <- units[units$cohort == g, , drop = FALSE]
    data.frame(cohort = g, units = nrow(uu), clusters = length(unique(uu$cluster_id)),
      mass = sum(uu$unit_weight), stringsAsFactors = FALSE)
  }))
  cell_counts <- do.call(rbind, lapply(groups, function(g) {
    do.call(rbind, lapply(seq_along(prepared$time_labels), function(tt) {
      z <- panel[panel$usable & panel$cohort == g & panel$time_index == tt, , drop = FALSE]
      data.frame(cohort = g, time_index = tt, observed_units = nrow(z),
        observed_clusters = length(unique(z$cluster_id)), mass = sum(z$unit_weight))
    }))
  }))
  eig <- eigen(covariance$Sigma, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(abs(eig))
  tolerance <- if (scale > 0) scale * 1e-9 else 0
  rank <- sum(eig > tolerance)
  expected_rank <- nrow(info) - nrow(normalization)
  coefficient_norm <- max(abs(normalization %*% beta))
  covariance_norm <- max(abs(normalization %*% covariance$Sigma))
  flags <- character(0)
  if (!prepared$balanced) flags <- c(flags, "unbalanced_bias_mapping_approximation")
  if (any(counts$clusters < 10L)) {
    flags <- c(flags, "few_clusters_in_cohort")
    warning("Some cohorts or the never-treated group contain fewer than 10 clusters. Cluster covariance and sensitivity intervals may be unreliable; inspect diagnostics$cohort_counts and diagnostics$cell_counts.", call. = FALSE)
  }
  if (any(cell_counts$observed_clusters <= 1L & cell_counts$observed_units > 0L)) {
    flags <- c(flags, "single_cluster_cell")
    warning("Some observed cohort-period cells contain only one cluster. Centered covariance can collapse in directions supported by these cells.", call. = FALSE)
  }
  if (rank < expected_rank) flags <- c(flags, "covariance_rank_below_normalized_dimension")
  if (any(diag(covariance$Sigma)[info$type == "post"] <=
          max(diag(covariance$Sigma), 0) * 1e-12)) flags <- c(flags, "near_zero_post_variance")
  list(flags = unique(flags), cohort_counts = counts, cell_counts = cell_counts,
    covariance_rank = rank, normalized_dimension = expected_rank,
    covariance_eigenvalues = eig, covariance_rank_tolerance = tolerance,
    normalization_error = coefficient_norm, covariance_normalization_error = covariance_norm,
    n_clusters = covariance$n_clusters, n_observations = covariance$n_observations,
    n_cells = covariance$n_cells, cr1 = covariance$cr1,
    balanced = prepared$balanced,
    bias_mapping = if (prepared$balanced) "exact_balanced_panel" else "balanced_panel_approximation",
    missingness = prepared$missingness)
}
