#' Construct the Cohort Bias Mapping Matrix
#'
#' Constructs the full mapping from cohort-specific primitive deviations to
#' estimator bias, using the paper's balanced-panel formulas.
#' @param cohort_times Integer adoption indices on the ordered analysis time grid.
#' @param cohort_sizes Positive cohort sizes or fixed cohort masses, in the same
#'   order as \code{cohort_times}. If named, names must be the adoption indices.
#' @param never_treated_size Positive size or fixed mass of the never-treated group.
#' @param n_periods Number of periods in the analysis window.
#' @param estimator One of \code{"imp"}, \code{"cs_nyt"}, or \code{"cs_nt"}.
#' @return A named square matrix \code{W}, with all pre-treatment coordinates
#'   followed by all post-treatment coordinates. Within each block, cohorts and
#'   periods are ordered increasingly. Relative period one is the adoption
#'   period.
#' @details For an unbalanced panel this fixed-mass mapping is an approximation;
#'   it does not generally reproduce the actual observation-specific bias map.
#' @export
#' @examples
#' W <- build_W(c(4, 6), c(20, 30), 50, 7, estimator = "imp")
#' dim(W)
build_W <- function(cohort_times, cohort_sizes, never_treated_size, n_periods,
                    estimator = "imp") {
  create_delta_to_Delta_mapping_matrix(cohort_times, cohort_sizes,
    never_treated_size, n_periods, estimator)$W
}

create_delta_to_Delta_mapping_matrix <- function(T_cohort, N_cohort,
    N_never_treated, T_total, estimator_type = "imputation") {
  estimator <- .hc_estimator(estimator_type)
  n_periods <- .hc_integer_scalar(T_total, "n_periods", 2L)
  .hc_assert(is.numeric(T_cohort) && length(T_cohort) > 0L &&
    all(is.finite(T_cohort) & T_cohort == round(T_cohort) &
          T_cohort >= 2 & T_cohort <= n_periods) && !anyDuplicated(T_cohort),
    "cohort_times must be distinct integer adoption indices from 2 to n_periods.")
  .hc_assert(is.numeric(N_cohort) && length(N_cohort) == length(T_cohort) &&
    all(is.finite(N_cohort) & N_cohort > 0),
    "cohort_sizes must contain one positive finite size or mass per cohort.")
  .hc_assert(is.numeric(N_never_treated) && length(N_never_treated) == 1L &&
    is.finite(N_never_treated) && N_never_treated > 0,
    "never_treated_size must be positive and finite.")
  if (!is.null(names(N_cohort))) {
    .hc_assert(!anyDuplicated(names(N_cohort)) &&
      setequal(names(N_cohort), as.character(T_cohort)),
      "Named cohort_sizes must use exactly the adoption indices in cohort_times.")
    N_cohort <- N_cohort[match(as.character(T_cohort), names(N_cohort))]
  }
  ord <- order(T_cohort)
  cohort_times <- as.integer(T_cohort[ord]); masses <- as.numeric(N_cohort[ord])
  info <- .hc_coordinate_info(cohort_times, seq_len(n_periods))
  W <- diag(nrow(info))
  dimnames(W) <- list(info$coordinate, info$coordinate)
  weight <- vapply(seq_along(masses), function(j) {
    denominator <- c(N_never_treated, masses[seq.int(j, length(masses))])
    scale <- max(denominator)
    (masses[j] / scale) / sum(denominator / scale)
  }, numeric(1))
  blocks <- list()
  if (estimator != "cs_nt") {
    for (r in which(info$type == "post")) {
      g <- info$cohort[r]; tt <- info$time_index[r]
      later <- which(cohort_times > g & cohort_times <= tt)
      for (j in later) {
        cc <- which(info$cohort == cohort_times[j] & info$time_index == tt)
        W[r, cc] <- W[r, cc] + weight[j]
        if (estimator == "cs_nyt") {
          base <- which(info$cohort == cohort_times[j] & info$time_index == g - 1L)
          W[r, base] <- W[r, base] - weight[j]
        }
      }
    }
    if (estimator == "imp") {
      for (tt in seq_len(n_periods)) {
        rr <- which(info$type == "post" & info$time_index == tt)
        if (length(rr) > 1L) {
          block <- W[rr, rr, drop = FALSE]
          dimnames(block) <- rep(list(paste0("g", info$cohort[rr])), 2L)
          blocks[[as.character(tt)]] <- block
        }
      }
    }
  }
  list(M = solve(W), W = W, W_matrices = blocks)
}
