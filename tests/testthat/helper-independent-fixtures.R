# These fixtures and references do not call package model builders.
independent_panel <- function() {
  d <- expand.grid(time = seq_len(7L), id = seq_len(66L))
  d <- d[order(d$id, d$time), ]
  g <- ifelse(d$id <= 20L, 4L, ifelse(d$id <= 42L, 6L, Inf))
  d$d <- as.integer(d$time >= g)
  d$y <- 0.7 + d$id / 40 + 0.06 * d$time +
    ifelse(is.finite(g), 0.006 * g, 0) * (d$time - 2) +
    0.015 * sin(d$id * d$time * 0.37) +
    0.02 * cos(d$id * 0.23 + d$time * 0.61) +
    0.0008 * (d$id %% 5) * d$time^2 +
    ifelse(d$d == 1L, -0.08 - 0.02 * (d$time - g) +
             0.01 * (d$id %% 3), 0)
  d$w <- 1 + (d$id %% 7) / 5
  d$state <- ceiling(d$id / 2)
  rownames(d) <- NULL
  d
}

independent_short_panel <- function() {
  d <- expand.grid(time = seq_len(3L), id = seq_len(24L))
  d <- d[order(d$id, d$time), ]
  d$d <- as.integer(d$id <= 12L & d$time >= 3L)
  d$y <- d$id / 20 + 0.04 * d$time +
    0.03 * sin(d$id * d$time * 0.7) +
    0.015 * cos(d$id * 0.31 + d$time * 0.4) - 0.1 * d$d
  rownames(d) <- NULL
  d
}

independent_two_cohort_short_panel <- function() {
  d <- expand.grid(time = seq_len(4L), id = seq_len(60L))
  d <- d[order(d$id, d$time), ]
  g <- ifelse(d$id <= 20L, 3L, ifelse(d$id <= 40L, 4L, Inf))
  d$d <- as.integer(d$time >= g)
  d$y <- d$id / 25 + 0.04 * d$time +
    0.03 * sin(d$id * d$time * 0.7) +
    0.015 * cos(d$id * 0.31 + d$time * 0.4) +
    0.001 * (d$id %% 5) * d$time^2 -
    0.1 * d$d - ifelse(d$d == 1L, 0.01 * (d$time - g), 0)
  rownames(d) <- NULL
  d
}

independent_cs <- function(d, info, never_only = TRUE, weights = NULL) {
  ids <- unique(d$id)
  adoption <- vapply(ids, function(i) {
    tt <- d$time[d$id == i & d$d == 1L]
    if (length(tt)) min(tt) else Inf
  }, numeric(1))
  names(adoption) <- as.character(ids)
  wt <- if (is.null(weights)) rep(1, nrow(d)) else d[[weights]]
  mean_at <- function(group, t) {
    j <- d$id %in% group & d$time == t
    sum(d$y[j] * wt[j]) / sum(wt[j])
  }
  vapply(seq_len(nrow(info)), function(j) {
    g <- info$cohort[j]
    t <- info$time_index[j]
    treated <- ids[adoption == g]
    controls <- if (never_only) ids[is.infinite(adoption)] else
      ids[adoption > max(t, g - 1L) & adoption != g]
    mean_at(treated, t) - mean_at(treated, g - 1L) -
      mean_at(controls, t) + mean_at(controls, g - 1L)
  }, numeric(1))
}

independent_folded_half_length <- function(bias, se, level = 0.95) {
  if (se == 0) return(bias)
  f <- function(h) stats::pnorm((h - bias) / se) -
    stats::pnorm((-h - bias) / se) - level
  stats::uniroot(f, c(0, bias + 12 * se), tol = 1e-12)$root
}

ordered_endpoints <- function(x) {
  z <- x$intervals
  ordering <- intersect(c("target", "M", "component", "component_id"), names(z))
  if (length(ordering)) z <- z[do.call(order, unname(z[ordering])), , drop = FALSE]
  unname(as.matrix(z[, c("lower", "upper"), drop = FALSE]))
}
