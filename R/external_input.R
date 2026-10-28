#' Supply externally estimated cohort coefficients
#'
#' Construct the common input object from externally estimated quantities.
#' The pre-treatment coefficients must be the paper's fixed-initial-control
#' block-bias estimates with the stated normalization. Ordinary imputation
#' residuals and an aggregated event study are not sufficient substitutes.
#'
#' @param betahat Named numeric vector of cohort-period coefficients, including
#'   pre-treatment and reference-period coefficients.
#' @param vcov Their joint sampling covariance, with matching row/column names.
#'   It is used as supplied: the \code{covariance} option of
#'   \code{\link{cohort_fit}()} does not apply here, and printed results
#'   describe it as supplied to \code{cohort_input()}. To use a CR2 covariance,
#'   compute it and pass it here.
#' @param coefficient_info Data frame with \code{coordinate}, \code{cohort}
#'   (integer adoption index), \code{relative_period} (one at adoption), and
#'   \code{type} (\code{"pre"}/\code{"post"}). Rows identify the supplied
#'   coefficient names.
#' @param W Square bias map with matching coordinate names: delta = W Delta.
#' @param normalization Matrix of deterministic coefficient constraints, with
#'   named columns. For CS this selects the zero reference coefficients; for
#'   imputation it encodes mean-zero pre-treatment coefficients.
#' @param cohort_sizes Named positive sizes or fixed masses, one per cohort.
#' @param estimator Estimator label: \code{"imp"}, \code{"cs_nyt"}, or
#'   \code{"cs_nt"}.
#' @param diagnostics Optional list describing covariance provenance and sample.
#' @return A \code{cohort_fit} object accepted by
#'   \code{\link{cohort_sensitivity}()}.
#' @export
cohort_input <- function(betahat, vcov, coefficient_info, W, normalization,
                         cohort_sizes, estimator = "imp", diagnostics = list()) {
  .hc_assert(is.character(estimator) && length(estimator) == 1L && !is.na(estimator),
             "estimator must be one character string.")
  estimator <- match.arg(estimator, c("imp", "cs_nyt", "cs_nt"))
  info <- as.data.frame(coefficient_info, stringsAsFactors = FALSE)
  need <- c("coordinate", "cohort", "relative_period", "type")
  if (!all(need %in% names(info))) stop("coefficient_info requires coordinate, cohort, relative_period, and type.", call. = FALSE)
  nm <- as.character(info$coordinate)
  if (!length(nm) || anyNA(nm) || anyDuplicated(nm) || any(!nzchar(nm))) stop("Coefficient labels must be unique and nonempty.", call. = FALSE)
  if (!is.numeric(betahat) || length(betahat) != length(nm) || is.null(names(betahat)) ||
      anyDuplicated(names(betahat)) || !setequal(names(betahat), nm)) stop("betahat must be uniquely named and match coefficient_info.", call. = FALSE)
  if (!is.numeric(info$cohort) || !is.numeric(info$relative_period) ||
      any(!is.finite(info$cohort)) || any(!is.finite(info$relative_period)) || any(info$cohort < 2) ||
      any(info$cohort > .Machine$integer.max) || any(abs(info$relative_period) > .Machine$integer.max) ||
      any(info$cohort != floor(info$cohort)) || any(info$relative_period != floor(info$relative_period))) stop("Cohort and relative-period indices must be representable R integers.", call. = FALSE)
  if (anyNA(info$type) || any(info$type != ifelse(info$relative_period <= 0, "pre", "post"))) stop("type must agree with relative_period.", call. = FALSE)
  square <- function(x, what) {
    x <- as.matrix(x)
    if (!is.numeric(x) || !identical(dim(x), c(length(nm), length(nm))) ||
        !setequal(rownames(x), nm) || !setequal(colnames(x), nm) || any(!is.finite(x))) stop(what, " must be a finite named square matrix.", call. = FALSE)
    x[nm, nm, drop = FALSE]
  }
  Sigma <- square(vcov, "vcov"); W <- square(W, "W")
  sigma.scale <- max(abs(Sigma))
  if (max(abs(Sigma - t(Sigma))) > 1e-8 * sigma.scale) stop("vcov must be symmetric.", call. = FALSE)
  Sigma <- Sigma + (t(Sigma) - Sigma) / 2
  if (sigma.scale > 0) {
    ev <- eigen(Sigma / sigma.scale, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -1e-8 * max(abs(ev))) stop("vcov must be positive semidefinite.", call. = FALSE)
  }
  if (qr(W)$rank != nrow(W)) stop("W must be invertible.", call. = FALSE)
  N <- as.matrix(normalization)
  if (!is.numeric(N) || ncol(N) != length(nm) || !setequal(colnames(N), nm) || any(!is.finite(N))) stop("normalization must have finite named columns matching the coefficients.", call. = FALSE)
  N <- N[, nm, drop = FALSE]
  .hc_check_mapping_support(W, N)
  beta <- betahat[nm]
  if (any(!is.finite(beta))) stop("betahat must be finite.", call. = FALSE)
  if (nrow(N)) {
    row.scale <- apply(abs(N), 1L, max)
    scaled.N <- N / ifelse(row.scale > 0, row.scale, 1)
    norm.tolerance <- 1e-7 * rowSums(abs(scaled.N))
    beta.scale <- max(abs(beta))
    if (beta.scale > 0 && any(abs(as.numeric(scaled.N %*% (beta / beta.scale))) > norm.tolerance)) stop("betahat violates the supplied normalization.", call. = FALSE)
    if (sigma.scale > 0 && any(abs(scaled.N %*% (Sigma / sigma.scale)) > norm.tolerance)) stop("vcov violates the supplied normalization.", call. = FALSE)
  }
  g <- sort(unique(info$cohort))
  if (!is.numeric(cohort_sizes) || is.null(names(cohort_sizes)) ||
      length(cohort_sizes) != length(g) || anyDuplicated(names(cohort_sizes)) ||
      !setequal(as.character(g), names(cohort_sizes)) || any(!is.finite(cohort_sizes)) || any(cohort_sizes <= 0)) stop("cohort_sizes must give one uniquely named positive size for each cohort.", call. = FALSE)
  canonical <- paste0("delta_", info$type, "_g", info$cohort, "_s", info$relative_period)
  if (anyDuplicated(canonical)) stop("Each cohort-period cell must occur once.", call. = FALSE)
  info$coordinate <- canonical
  info$time <- info$cohort + info$relative_period - 1L
  names(beta) <- canonical; dimnames(Sigma) <- dimnames(W) <- list(canonical, canonical)
  colnames(N) <- canonical
  structure(list(call = match.call(), estimator = estimator, betahat = beta, Sigma = Sigma,
                 W = W, normalization = N, coordinate_info = info,
                 cohort_sizes = cohort_sizes, cohort_masses = cohort_sizes,
                 total_periods = max(info$time), diagnostics = diagnostics,
                 data_info = list(source = "external", balanced = NA)), class = "cohort_fit")
}
