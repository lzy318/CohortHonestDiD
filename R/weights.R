#' FLCI weights on the pre-treatment coefficients
#'
#' Tabulate the linear weights of fixed-length confidence intervals (FLCIs):
#' the weights \code{v} on the pre-treatment coefficients or the weights
#' \code{gamma} on the pure pre-treatment second differences, for every target
#' and bound of an FLCI result. \code{plot(x, type = "weights")} draws them.
#'
#' @param x A \code{cohort_sensitivity} object computed with
#'   \code{method = "FLCI"}.
#' @param weights \code{"coefficients"} (default) for the weights \code{v} on
#'   the pre-treatment coefficients, or \code{"gamma"} for the weights on the
#'   pure pre-treatment second differences.
#' @param by \code{"element"} (default) for one row per coefficient or second
#'   difference, or \code{"cohort"} for sums over each cohort's rows (cohort
#'   inference only).
#' @return A data frame. With \code{by = "element"}: one row per target, bound
#'   and element, with columns \code{target} and \code{label} (the target's id
#'   and label), \code{M}, \code{cohort} (\code{NA} under aggregate inference),
#'   \code{relative_period} (the coefficient's period, or the last period of
#'   the second difference), \code{element} (the coefficient name, or
#'   \code{d2_g<cohort>_s<period>} for a second difference), \code{weight},
#'   \code{estimate} (the coefficient estimate, or the estimated second
#'   difference), \code{contribution} (\code{weight} times \code{estimate}), and
#'   \code{tie.break} (\code{TRUE} when the FLCI chose its weights by the
#'   minimum-norm rule described below). With \code{by = "cohort"}: one row per
#'   target, bound and cohort, with \code{elements} (the number of rows
#'   summed), \code{weight} (their sum), \code{l1} (the sum of absolute
#'   weights), \code{share} (the cohort's part of the total \code{l1}),
#'   \code{contribution} (the sum of contributions), and \code{tie.break}. Rows
#'   follow the targets, then the sorted bounds, then cohorts and periods in
#'   increasing order. Attributes: \code{weights}, \code{by}, \code{framework},
#'   \code{full.path}; \code{centers}, one row per target and bound with the
#'   conventional \code{estimate}, the post-only center
#'   \code{post.only.center}, and the FLCI \code{center}; and \code{missing},
#'   the targets and bounds whose FLCI failed and that therefore have no
#'   weights, with their \code{status}.
#' @details The FLCI is a linear estimator: its weights on the post-treatment
#'   coefficients are the target weights, and its weights on the pre-treatment
#'   coefficients are \code{v}. Under the post-only SD restriction, \code{v} is
#'   fixed at the post-only weights, which subtract each cohort's linear
#'   extrapolation from its last two pre-treatment periods (with CS-NYT, they
#'   also undo the base-period correction of the adjustment cohorts); the
#'   resulting estimate is the post-only center. Under the full-path
#'   restriction, \code{v} equals the post-only weights plus
#'   \code{t(D) \%*\% gamma}, where the rows of \code{D} form the pure
#'   pre-treatment second differences \code{zhat = D \%*\% betahat}, and the
#'   FLCI center is the post-only center plus \code{sum(gamma * zhat)}. Hence
#'   the contributions of the \code{gamma} rows add up to the FLCI center
#'   minus the post-only center, and those of the \code{v} rows to the FLCI
#'   center minus the conventional estimate. The post-only FLCI has
#'   \code{gamma = 0}.
#'
#'   Each cohort's pre-treatment coefficients satisfy one normalization (CS:
#'   the coefficient at period 0 is zero; imputation: they sum to zero), so
#'   weights that differ along it give the same estimate.
#'   \code{flci_weights()} reports \code{v} in the form above, the paper's;
#'   \code{results[[k]]$weights$coefficients} stores its projection orthogonal
#'   to the normalization. The cohort sums of the contributions are the same
#'   for both.
#'
#'   Under aggregate inference the rows are the aggregated coefficients and
#'   their second differences. With \code{balanced = TRUE} and
#'   \code{cohort.paths = "window"}, the \code{gamma} rows are the second
#'   differences inside the balanced window. When several weight vectors give
#'   the same FLCI length, as can happen at \code{M = 0}, the FLCI uses the one
#'   with the smallest Euclidean norm of \code{gamma}; \code{tie.break} marks
#'   these rows.
#' @seealso \code{\link{cohort_sensitivity}}, \code{\link{plot_cohort}}
#' @export
#' @examples
#' data(minimum_wage)
#' \donttest{
#' fit <- cohort_fit(Y ~ D, minimum_wage, index = c("id", "year"))
#' full <- cohort_sensitivity(fit, M = c(0, 0.005, 0.01), method = "FLCI",
#'                            seed = 1)
#' post <- cohort_sensitivity(fit, M = c(0, 0.005, 0.01), method = "FLCI",
#'                            full.path = FALSE, seed = 1)
#' head(flci_weights(full))
#' flci_weights(full, weights = "gamma", by = "cohort")
#' plot(post, type = "weights")
#' plot(full, type = "weights", weights = "gamma", labels = TRUE)
#' }
flci_weights <- function(x, weights = c("coefficients", "gamma"), by = c("element", "cohort")) {
  if (!inherits(x, "cohort_sensitivity")) stop("x must be a cohort_sensitivity object.", call. = FALSE)
  weights <- match.arg(weights)
  by <- match.arg(by)
  s <- x$settings
  .hc_flci_weights_method(s)
  if (by == "cohort" && identical(s$framework, "aggregate")) .hc_flci_weights_no_cohorts()
  aggregate <- identical(s$framework, "aggregate")
  bounds <- s$M
  n.bounds <- length(bounds)
  if (!n.bounds || length(x$results) != length(x$targets) * n.bounds) .hc_flci_weights_fail()
  window <- s$window
  paths <- if (isTRUE(s$balanced) && identical(s$framework, "cohort")) s$cohort.paths else "full"
  rows <- list(); centers <- list(); missing <- list()
  for (j in seq_along(x$targets)) {
    one <- x$targets[[j]]
    tasks <- (j - 1L) * n.bounds + seq_len(n.bounds)
    has.weights <- vapply(x$results[tasks], function(r) !is.null(r$weights), logical(1))
    for (m in which(!has.weights)) {
      r <- x$results[[tasks[m]]]
      missing[[length(missing) + 1L]] <- data.frame(target = one$id, label = one$label,
        M = bounds[m], status = if (length(r$status)) as.character(r$status)[1L] else "unknown",
        stringsAsFactors = FALSE)
    }
    if (!any(has.weights)) next
    representation <- .hc_representation(x$fit, one$target, s$framework, window, paths)
    ops <- tryCatch(.hc_flci_operators(representation), error = function(e) NULL)
    if (is.null(ops)) .hc_flci_weights_fail()
    for (m in which(has.weights)) {
      r <- x$results[[tasks[m]]]
      v <- .hc_flci_reconstruct(ops, r$weights)
      gamma <- as.numeric(r$weights$gamma)
      tie <- identical(r$diagnostics$tie.break, "minimum_euclidean_norm_gamma")
      centers[[length(centers) + 1L]] <- data.frame(target = one$id, label = one$label,
        M = bounds[m], estimate = as.numeric(x$estimates$estimate[j]),
        post.only.center = ops$c0,
        center = if (length(r$center) == 1L) as.numeric(r$center) else NA_real_,
        stringsAsFactors = FALSE)
      if (weights == "coefficients") {
        cells <- ops$pre
        value <- v
        estimate <- ops$b.pre
      } else {
        cells <- ops$pure
        value <- gamma
        estimate <- ops$zhat
      }
      n <- nrow(cells)
      if (!n) next
      element <- data.frame(target = rep(one$id, n), label = rep(one$label, n),
        M = rep(bounds[m], n),
        cohort = if (aggregate) rep(NA_integer_, n) else cells$cohort,
        relative_period = cells$relative_period, element = cells$element,
        weight = value, estimate = estimate, contribution = value * estimate,
        tie.break = rep(tie, n), stringsAsFactors = FALSE)
      rows[[length(rows) + 1L]] <- if (by == "element") element else .hc_flci_cohort_sums(element)
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else .hc_flci_weights_empty(by)
  rownames(out) <- NULL
  centers <- if (length(centers)) do.call(rbind, centers) else data.frame(target = character(),
    label = character(), M = numeric(), estimate = numeric(), post.only.center = numeric(),
    center = numeric(), stringsAsFactors = FALSE)
  missing <- if (length(missing)) do.call(rbind, missing) else data.frame(target = character(),
    label = character(), M = numeric(), status = character(), stringsAsFactors = FALSE)
  rownames(centers) <- NULL
  rownames(missing) <- NULL
  attr(out, "weights") <- weights
  attr(out, "by") <- by
  attr(out, "framework") <- s$framework
  attr(out, "full.path") <- s$full.path
  attr(out, "centers") <- centers
  attr(out, "missing") <- missing
  out
}

# Tolerances of the weight reconstruction: the rank tolerance of the
# normalization's QR decomposition and the relative reconstruction check.
.hc_weights_rank_tolerance <- 1e-10
.hc_weights_check_tolerance <- 1e-8

#' @noRd
.hc_flci_weights_method <- function(s) {
  if (identical(s$method, "flci")) return(invisible(TRUE))
  used <- if (identical(s$method, "conditional")) "Conditional" else
    if (identical(s$method, "c_lf")) "C-LF" else paste(s$method, collapse = ", ")
  stop("FLCI weights exist only for results with method = 'FLCI'; this result used ",
       used, ".", call. = FALSE)
}

#' @noRd
.hc_flci_weights_no_cohorts <- function() {
  stop("by = 'cohort' requires cohort inference: aggregate weights have no cohorts.", call. = FALSE)
}

#' @noRd
.hc_flci_weights_fail <- function() {
  stop("flci_weights() could not reconstruct the FLCI weights of this result.", call. = FALSE)
}

#' @noRd
.hc_flci_weights_empty <- function(by) {
  if (by == "element") return(data.frame(target = character(), label = character(),
    M = numeric(), cohort = integer(), relative_period = integer(), element = character(),
    weight = numeric(), estimate = numeric(), contribution = numeric(), tie.break = logical(),
    stringsAsFactors = FALSE))
  data.frame(target = character(), label = character(), M = numeric(), cohort = integer(),
             elements = integer(), weight = numeric(), l1 = numeric(), share = numeric(),
             contribution = numeric(), tie.break = logical(), stringsAsFactors = FALSE)
}

# Second-difference operators, post-only weights, estimated pure pre-treatment
# second differences, and the post-only center of one representation.
#' @noRd
.hc_flci_operators <- function(representation) {
  info <- representation$coordinate_info
  coord <- as.character(info$coordinate)
  cohort <- as.integer(info$cohort)
  period <- as.integer(round(info$relative_period))
  ord <- order(cohort, period)
  pre <- ord[info$type[ord] == "pre"]
  post <- ord[info$type[ord] == "post"]
  P <- coord[pre]
  Q <- coord[post]
  L <- representation[["restriction_window"]]
  d2 <- list(); meta <- list()
  for (g in sort(unique(cohort))) {
    path <- ord[cohort[ord] == g]
    if (length(path) < 3L) next
    for (i in seq.int(3L, length(path))) {
      e <- period[path[i]]
      row <- stats::setNames(numeric(length(coord)), coord)
      row[coord[path[(i - 2L):i]]] <- c(1, -2, 1)
      kind <- if (e >= 1L) "post" else if (is.null(L) || e - 2L >= L) "pure" else "dropped"
      d2[[length(d2) + 1L]] <- row
      meta[[length(meta) + 1L]] <- data.frame(element = paste0("d2_g", g, "_s", e), cohort = g,
        relative_period = e, end = coord[path[i]], kind = kind, stringsAsFactors = FALSE)
    }
  }
  D <- matrix(0, length(d2), length(coord), dimnames = list(NULL, coord))
  for (i in seq_along(d2)) D[i, ] <- d2[[i]]
  meta <- if (length(meta)) do.call(rbind, meta) else data.frame(element = character(),
    cohort = integer(), relative_period = integer(), end = character(), kind = character(),
    stringsAsFactors = FALSE)
  pure <- which(meta$kind == "pure")
  connected <- which(meta$kind == "post")
  connected <- connected[match(Q, meta$end[connected])]
  if (anyNA(connected) || length(connected) != length(Q) || !length(Q)) stop("no post-connected rows")
  D_PP <- D[pure, P, drop = FALSE]
  D_QP <- D[connected, P, drop = FALSE]
  D_QQ <- D[connected, Q, drop = FALSE]
  W <- as.matrix(representation$W)
  l <- representation$target
  b <- representation$betahat
  l_Q <- as.numeric(l[Q])
  xi <- solve(t(D_QQ), crossprod(W[Q, Q, drop = FALSE], l_Q))
  vP <- as.numeric(crossprod(D_QP, xi) - crossprod(W[Q, P, drop = FALSE], l_Q))
  names(vP) <- P
  b.pre <- as.numeric(b[P])
  zhat <- as.numeric(D_PP %*% b.pre)
  c0 <- sum(vP * b.pre) + sum(l_Q * as.numeric(b[Q]))
  N <- as.matrix(representation$normalization)
  N_P <- N[, P, drop = FALSE]
  Pi <- diag(length(P))
  if (nrow(N_P)) {
    qrN <- qr(t(N_P), tol = .hc_weights_rank_tolerance)
    Q1 <- qr.Q(qrN)[, seq_len(qrN$rank), drop = FALSE]
    Pi <- Pi - tcrossprod(Q1)
  }
  if (any(!is.finite(c(vP, zhat, c0, Pi)))) stop("nonfinite operators")
  list(pre = data.frame(cohort = cohort[pre], relative_period = period[pre], element = P,
                        stringsAsFactors = FALSE),
       pure = meta[pure, c("cohort", "relative_period", "element"), drop = FALSE],
       P = P, D_PP = D_PP, vP = vP, b.pre = b.pre, zhat = zhat, c0 = c0, Pi = Pi)
}

# The paper's representative v = vP + t(D_PP) gamma of one task's weights,
# checked against the gamma and the projected weights the FLCI stored.
#' @noRd
.hc_flci_reconstruct <- function(ops, stored) {
  gamma <- stored$gamma
  if (length(gamma) != nrow(ops$pure) || any(!is.finite(as.numeric(gamma))))
    .hc_flci_weights_fail()
  if (!is.null(names(gamma)) && !identical(as.character(names(gamma)), paste0(ops$pure$element, "_pos")))
    .hc_flci_weights_fail()
  v <- ops$vP + as.numeric(crossprod(ops$D_PP, as.numeric(gamma)))
  w <- stored$coefficients
  if (is.null(names(w)) || !all(ops$P %in% names(w))) .hc_flci_weights_fail()
  w.pre <- as.numeric(w[ops$P])
  gap <- max(abs(as.numeric(ops$Pi %*% v) - w.pre))
  if (!is.finite(gap) || gap > .hc_weights_check_tolerance * max(1, max(abs(w.pre))))
    .hc_flci_weights_fail()
  stats::setNames(v, ops$P)
}

# Sums of one task's element rows within each cohort.
#' @noRd
.hc_flci_cohort_sums <- function(element) {
  cohorts <- sort(unique(element$cohort))
  pick <- lapply(cohorts, function(g) element$cohort == g)
  l1 <- vapply(pick, function(i) sum(abs(element$weight[i])), numeric(1))
  total <- sum(l1)
  data.frame(target = rep(element$target[1L], length(cohorts)),
             label = rep(element$label[1L], length(cohorts)),
             M = rep(element$M[1L], length(cohorts)), cohort = as.integer(cohorts),
             elements = vapply(pick, function(i) as.integer(sum(i)), integer(1)),
             weight = vapply(pick, function(i) sum(element$weight[i]), numeric(1)),
             l1 = l1, share = if (total > 0) l1 / total else rep(NA_real_, length(cohorts)),
             contribution = vapply(pick, function(i) sum(element$contribution[i]), numeric(1)),
             tie.break = rep(element$tie.break[1L], length(cohorts)), stringsAsFactors = FALSE)
}
