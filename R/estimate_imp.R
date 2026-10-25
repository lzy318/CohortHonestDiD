# Sparse linear operators retain the actual outcome histories used for estimation.
.hc_sparse_rows <- function(rows, n_columns) {
  lengths <- vapply(rows, function(z) length(z$i), integer(1))
  Matrix::sparseMatrix(i = rep(seq_along(rows), lengths),
    j = unlist(lapply(rows, `[[`, "i"), use.names = FALSE),
    x = unlist(lapply(rows, `[[`, "x"), use.names = FALSE),
    dims = c(length(rows), n_columns))
}

.hc_mean_contrast <- function(work, positive, negative = integer(0)) {
  .hc_assert(length(positive) > 0L, "A required treated cohort-period mean is not estimable.")
  p <- work$unit_weight[positive]
  if (!length(negative)) return(list(i = positive, x = p / sum(p)))
  n <- work$unit_weight[negative]
  list(i = c(positive, negative), x = c(p / sum(p), -n / sum(n)))
}

.hc_estimate_imp <- function(work, info) {
  rows <- vector("list", nrow(info))
  for (j in seq_len(nrow(info))) {
    g <- info$cohort[j]; tt <- info$time_index[j]
    treated <- which(work$cohort == g & work$time_index == tt)
    if (info$type[j] == "pre") {
      controls <- which((work$cohort == 0L | work$cohort > g) & work$time_index == tt)
      .hc_assert(length(controls) > 0L, paste0("No initial comparison observations for cohort ", g,
        " at time index ", tt, "."))
      rows[[j]] <- .hc_mean_contrast(work, treated, controls)
    } else rows[[j]] <- .hc_mean_contrast(work, treated)
  }
  A <- .hc_sparse_rows(rows, nrow(work))
  center <- Matrix::Diagonal(nrow(info))
  for (g in unique(info$cohort)) {
    rr <- which(info$cohort == g & info$type == "pre")
    center[rr, rr] <- diag(length(rr)) - 1 / length(rr)
  }
  A <- center %*% A
  units <- sort(unique(work$unit_index))
  n_units <- length(units); n_periods <- max(work$time_index)
  ui <- match(work$unit_index, units)
  untreated <- work$D == 0L
  counts <- tabulate(ui[untreated], nbins = n_units)
  .hc_assert(all(counts > 0L), paste0("Imputation requires at least one observed untreated outcome for every unit. Working unit indices without one: ",
    paste(units[counts == 0L], collapse = ", "), "."))
  nonbase <- which(work$time_index > 1L)
  X <- Matrix::sparseMatrix(i = c(seq_len(nrow(work)), nonbase),
    j = c(ui, n_units + work$time_index[nonbase] - 1L), x = 1,
    dims = c(nrow(work), n_units + n_periods - 1L))
  fit_weights <- work$unit_weight * untreated
  weighted_X <- Matrix::Diagonal(x = sqrt(fit_weights)) %*% X
  gram <- Matrix::crossprod(weighted_X)
  factor <- tryCatch(suppressWarnings(Matrix::Cholesky(gram, LDL = FALSE)),
                     error = function(e) NULL)
  .hc_assert(!is.null(factor), paste0("The untreated unit/time fixed-effect design is not identified. ",
    "Every analysis period must be connected to the reference period through units with observed untreated outcomes."))
  H <- A %*% X
  H[info$type == "pre", ] <- 0
  correction <- t(as.matrix(Matrix::solve(factor, Matrix::t(H), system = "A")))
  .hc_assert(all(is.finite(correction)), "The untreated fixed-effect solve produced nonfinite coefficients.")
  list(A = A, X = X, fit_weights = fit_weights, correction = correction)
}

.hc_apply_operator <- function(operator, y) {
  out <- as.numeric(operator$A %*% y)
  if (!is.null(operator$correction)) {
    out <- out - as.numeric(operator$correction %*%
      Matrix::crossprod(operator$X, operator$fit_weights * y))
  }
  out
}

.hc_operator_scores <- function(operator, residual_by_cluster) {
  out <- operator$A %*% residual_by_cluster
  if (!is.null(operator$correction)) {
    rhs <- Matrix::crossprod(operator$X,
      Matrix::Diagonal(x = operator$fit_weights) %*% residual_by_cluster)
    out <- out - operator$correction %*% rhs
  }
  t(as.matrix(out))
}
