# Internal target construction. All weights refer to the original cohort cells.

# The balanced window: the relative periods observed for every cohort.
.hc_balanced_window <- function(fit) {
  info <- fit$coordinate_info
  s <- as.integer(round(info$relative_period))
  g <- info$cohort
  cohorts <- unique(g)
  observed <- function(p) all(vapply(cohorts, function(k) any(g == k & s == p), logical(1)))
  if (!observed(1L)) stop("balanced = TRUE requires relative period 1 for every cohort: ",
                          "no post-treatment period is observed for all cohorts.", call. = FALSE)
  if (!observed(0L)) stop("balanced = TRUE requires relative period 0 for every cohort: ",
                          "no pre-treatment period is observed for all cohorts.", call. = FALSE)
  upper <- 1L
  while (observed(upper + 1L)) upper <- upper + 1L
  lower <- 0L
  while (observed(lower - 1L)) lower <- lower - 1L
  c(lower = lower, upper = upper)
}

# Under cohort.paths = "window" the coefficients after the window are dropped,
# so no kept row of W and no normalization row may load on them.
.hc_check_window_paths <- function(fit, window) {
  info <- fit$coordinate_info
  drop <- info$type == "post" & info$relative_period > window[["upper"]]
  if (!any(drop)) return(invisible(TRUE))
  W <- as.matrix(fit$W)[info$coordinate, info$coordinate, drop = FALSE]
  N <- as.matrix(fit$normalization)[, info$coordinate, drop = FALSE]
  w.link <- max(abs(W[!drop, drop, drop = FALSE]))
  n.link <- if (nrow(N)) max(abs(N[, drop, drop = FALSE])) else 0
  n.scale <- if (nrow(N)) max(abs(N)) else 0
  if (w.link > 1e-12 * max(abs(W)) || n.link > 1e-12 * n.scale)
    stop("cohort.paths = 'window' requires that the bias map W and the normalization link no ",
         "coefficient inside the balanced window to a coefficient after it.", call. = FALSE)
  invisible(TRUE)
}

.hc_targets <- function(fit, target, window = NULL) {
  info <- fit$coordinate_info
  post <- which(info$type == "post")
  if (!is.null(window)) post <- post[info$relative_period[post] <= window[["upper"]]]
  masses <- fit$cohort_masses
  if (is.null(masses)) masses <- fit$cohort_sizes
  if (is.null(names(masses))) names(masses) <- sort(unique(info$cohort))
  mass <- as.numeric(masses[as.character(info$cohort)])
  if (any(!is.finite(mass)) || any(mass <= 0)) stop("Missing or invalid cohort masses.", call. = FALSE)
  make <- function(indices, id, label, period = NA_integer_, cohort = NA_integer_) {
    ell <- setNames(numeric(nrow(info)), info$coordinate)
    relative.mass <- mass[indices] / max(mass[indices])
    ell[indices] <- relative.mass / sum(relative.mass)
    list(id = id, label = label, target = ell, relative_period = period,
         cohort = cohort, cells = info$coordinate[indices])
  }
  if (target == "att") {
    label <- if (is.null(window)) "Average treatment effect" else if (window[["upper"]] == 1L)
      "Average treatment effect, relative period 1" else
      paste0("Average treatment effect, relative periods 1 to ", window[["upper"]])
    return(list(make(post, "att", label)))
  }
  if (target == "dynamic") return(lapply(sort(unique(info$relative_period[post])), function(s) {
    make(post[info$relative_period[post] == s], paste0("s", s),
         paste0("Relative period ", s), period = s)
  }))
  lapply(post, function(j) make(j, paste0("g", info$cohort[j], "_s", info$relative_period[j]),
                              paste0("Cohort ", info$cohort[j], ", period ", info$relative_period[j]),
                              period = info$relative_period[j], cohort = info$cohort[j]))
}

.hc_representation <- function(fit, ell, framework, window = NULL, paths = "full") {
  core <- fit[c("coordinate_info", "betahat", "Sigma", "W", "normalization")]
  core$name <- "cohort"
  core$target <- ell
  if (framework == "cohort") {
    if (is.null(window) || paths == "full") return(core)
    # Window paths: keep every pre-treatment cell (W and the imputation
    # normalization use them) and the post cells up to the window's end.
    info <- fit$coordinate_info
    keep <- info$type == "pre" | info$relative_period <= window[["upper"]]
    nm <- info$coordinate[keep]
    core$coordinate_info <- info[keep, , drop = FALSE]
    core$betahat <- fit$betahat[nm]
    core$Sigma <- fit$Sigma[nm, nm, drop = FALSE]
    core$W <- fit$W[nm, nm, drop = FALSE]
    core$normalization <- fit$normalization[, nm, drop = FALSE]
    core$target <- ell[nm]
    core$restriction_window <- as.integer(window[["lower"]])
    return(core)
  }
  info <- fit$coordinate_info
  periods <- if (is.null(window)) sort(unique(info$relative_period)) else
    seq.int(window[["lower"]], window[["upper"]])
  nm <- paste0("delta_", ifelse(periods <= 0, "pre", "post"), "_g1_s", periods)
  A <- matrix(0, length(periods), nrow(info), dimnames = list(nm, info$coordinate))
  masses <- fit$cohort_masses
  if (is.null(masses)) masses <- fit$cohort_sizes
  if (is.null(names(masses))) names(masses) <- sort(unique(info$cohort))
  for (j in seq_along(periods)) {
    idx <- which(info$relative_period == periods[j])
    w <- as.numeric(masses[as.character(info$cohort[idx])])
    w <- w / max(w)
    A[j, idx] <- w / sum(w)
  }
  ell_A <- as.numeric(qr.solve(t(A), ell, tol = 1e-10))
  names(ell_A) <- nm
  if (max(abs(as.numeric(t(A) %*% ell_A) - ell)) > 1e-9) {
    stop("The aggregate event study cannot represent this target.", call. = FALSE)
  }
  # Carry mechanical normalizations through aggregation; do not impose a new base.
  null_basis <- function(mat, tol = 1e-10) {
    nc <- ncol(mat)
    if (!nrow(mat)) return(diag(nc))
    q <- qr(t(mat), tol = tol, LAPACK = FALSE)
    Q <- qr.Q(q, complete = TRUE)
    if (q$rank >= nc) return(matrix(0, nc, 0L))
    Q[, seq.int(q$rank + 1L, nc), drop = FALSE]
  }
  embedding <- null_basis(fit$normalization)
  norm_A <- t(null_basis(t(A %*% embedding)))
  colnames(norm_A) <- nm
  if (nrow(norm_A)) rownames(norm_A) <- paste0("aggregate_normalization_", seq_len(nrow(norm_A)))
  Sigma <- A %*% fit$Sigma %*% t(A)
  W <- diag(length(nm)); dimnames(W) <- list(nm, nm)
  list(name = "aggregate", betahat = setNames(as.numeric(A %*% fit$betahat), nm),
       Sigma = (Sigma + t(Sigma)) / 2, W = W, normalization = norm_A, target = ell_A,
       coordinate_info = data.frame(coordinate = nm, type = ifelse(periods <= 0, "pre", "post"),
                                    cohort = 1L, relative_period = periods, time = seq_along(periods)),
       aggregation_matrix = A)
}
