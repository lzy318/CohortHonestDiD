.hc_covariance_type <- function(covariance) {
  .hc_assert(is.character(covariance) && length(covariance) == 1L && !is.na(covariance) &&
               toupper(covariance) %in% c("CR1", "CR2"),
             "covariance must be 'CR1' or 'CR2'.")
  toupper(covariance)
}

# CR2 shares and factors for each (cluster, cohort group) pair. Groups are the
# never-treated units (cohort 0) and each adoption cohort; shares use the fixed
# unit weights of the registry. A pair holding every unit of its group gets a
# zero factor, decided by unit counts rather than a floating tolerance.
.hc_cr2_adjustment <- function(registry, clusters) {
  cl <- match(registry$cluster_id, clusters)
  .hc_assert(!anyNA(cl), "Every unit's cluster must have usable outcomes.")
  grp <- as.integer(registry$cohort)
  w <- registry$unit_weight
  groups <- sort(unique(grp))
  gi <- match(grp, groups)
  nG <- length(groups)
  pair <- (cl - 1L) * nG + gi
  pairs <- sort(unique(pair))
  pi <- match(pair, pairs)
  pair_cluster <- (pairs - 1L) %/% nG + 1L
  pair_group <- (pairs - 1L) %% nG + 1L
  n_pair <- tabulate(pi, nbins = length(pairs))
  m_pair <- as.numeric(rowsum(w, pi, reorder = TRUE))
  N_group <- tabulate(gi, nbins = nG)
  M_group <- as.numeric(rowsum(w, gi, reorder = TRUE))
  whole <- n_pair == N_group[pair_group]
  total <- M_group[pair_group]
  # Weight of the group's units outside the cluster. A pair holding more than
  # half of its group's weight (at most one per group) sums the rest directly
  # to avoid cancellation.
  rest <- total - m_pair
  for (k in which(!whole & m_pair > total / 2)) {
    rest[k] <- sum(w[gi == pair_group[k] & cl != pair_cluster[k]])
  }
  undefined <- "CR2 is undefined: a cluster holds all but a numerically negligible share of its cohort's weight."
  .hc_assert(all(rest[!whole] > 0), undefined)
  factor <- numeric(length(pairs))
  factor[!whole] <- sqrt(total[!whole] / rest[!whole])
  .hc_assert(all(is.finite(factor)), undefined)
  share <- m_pair / total
  share[whole] <- 1
  list(shares = data.frame(cluster = clusters[pair_cluster], cohort = groups[pair_group],
                           units = n_pair, mass = m_pair, share = share, factor = factor,
                           stringsAsFactors = FALSE),
       factor_by_unit = factor[pi])
}

.hc_covariance <- function(operator, prepared, work, info, type = "CR1") {
  observed <- prepared$panel[prepared$panel$usable, , drop = FALSE]
  cell <- interaction(observed$cohort, observed$time_index, drop = TRUE, lex.order = TRUE)
  cell_index <- as.integer(cell)
  cell_mass <- as.numeric(rowsum(observed$unit_weight, cell_index, reorder = TRUE))
  cell_mean <- as.numeric(rowsum(observed$unit_weight * observed$Y,
                                cell_index, reorder = TRUE)) / cell_mass
  residual <- observed$Y - cell_mean[cell_index]
  clusters <- sort(unique(observed$cluster_id))
  cluster_index <- match(observed$cluster_id, clusters)
  n_obs <- nrow(observed); n_cells <- nlevels(cell); n_clusters <- length(clusters)
  .hc_assert(n_clusters >= 2L, "At least two clusters with usable outcomes are required.")
  .hc_assert(n_obs > n_cells, paste0("Cohort-period centering leaves no residual degrees of freedom (",
    n_obs, " observations and ", n_cells, " nonempty cells). Covariance cannot be estimated."))
  if (type == "CR2") adjustment <- .hc_cr2_adjustment(prepared$unit_registry, clusters)
  if (work$compressed) {
    rr <- work$observed_to_working
    xx <- observed$unit_weight * residual / work$cell_mass[rr]
    # A centered cohort/time mean supported by a single cluster has exactly
    # zero score. Enforce this identity before summation, without a variance
    # cutoff that could erase genuine small positive sampling variation.
    single_cluster <- vapply(split(cluster_index, rr), function(z)
      length(unique(z)) == 1L, logical(1))
    xx[as.character(rr) %in% names(single_cluster)[single_cluster]] <- 0
  } else {
    rr <- seq_len(n_obs); xx <- residual
  }
  # CR2: each observed row is scaled by the factor of its unit's (cluster,
  # cohort group) pair. A working row is one cell of one group, so every
  # (working row, cluster) entry of R is its raw entry times one pair factor.
  if (type == "CR2") xx <- xx * adjustment$factor_by_unit[observed$unit_index]
  R <- Matrix::sparseMatrix(i = rr, j = cluster_index, x = xx,
    dims = c(nrow(work$data), n_clusters))
  scores <- .hc_operator_scores(operator, R)
  dimnames(scores) <- list(clusters, info$coordinate)
  if (type == "CR1") {
    cr1 <- n_clusters / (n_clusters - 1L) * (n_obs - 1L) / (n_obs - n_cells)
    Sigma <- crossprod(scores) * cr1
  } else {
    cr1 <- NA_real_
    Sigma <- crossprod(scores)
  }
  Sigma <- (Sigma + t(Sigma)) / 2
  dimnames(Sigma) <- list(info$coordinate, info$coordinate)
  out <- list(Sigma = Sigma, scores = scores, cr1 = cr1,
              n_observations = n_obs, n_cells = n_cells, n_clusters = n_clusters)
  if (type == "CR2") out$cluster_shares <- adjustment$shares
  out
}
