# FLCI weights: references are computed here from the public fields of fits and
# results (betahat, W, normalization, coordinate_info, targets, results).

fw_cell <- function(info, g, e) info$coordinate[info$cohort == g & info$relative_period == e]

# Closed form of the post-only weights vP, the pure pre-treatment second
# differences, their estimates zhat, and the post-only center c0.
fw_reference <- function(info, b, W, ell, window_lower = NULL) {
  info <- info[order(info$cohort, info$relative_period),
               c("coordinate", "type", "cohort", "relative_period")]
  P <- info$coordinate[info$type == "pre"]
  Q <- info$coordinate[info$type == "post"]
  pure <- list(); connected <- list(); pure_cohort <- integer(); pure_period <- integer()
  for (k in seq_len(nrow(info))) {
    g <- info$cohort[k]; e <- info$relative_period[k]
    one <- fw_cell(info, g, e - 1L); two <- fw_cell(info, g, e - 2L)
    if (!length(one) || !length(two)) next
    r <- setNames(numeric(nrow(info)), info$coordinate)
    r[c(two, one, info$coordinate[k])] <- c(1, -2, 1)
    if (e >= 1L) {
      connected[[info$coordinate[k]]] <- r
    } else if (is.null(window_lower) || e - 2L >= window_lower) {
      pure[[paste0("d2_g", g, "_s", e)]] <- r
      pure_cohort <- c(pure_cohort, g); pure_period <- c(pure_period, e)
    }
  }
  D_pure <- if (length(pure)) do.call(rbind, pure) else
    matrix(0, 0, nrow(info), dimnames = list(NULL, info$coordinate))
  D_post <- do.call(rbind, connected)[Q, , drop = FALSE]
  ell_Q <- as.numeric(ell[Q])
  xi <- solve(t(D_post[, Q, drop = FALSE]), t(W[Q, Q, drop = FALSE]) %*% ell_Q)
  vP <- setNames(as.numeric(t(D_post[, P, drop = FALSE]) %*% xi -
                              t(W[Q, P, drop = FALSE]) %*% ell_Q), P)
  b_P <- as.numeric(b[P])
  list(P = P, Q = Q, D_pure = D_pure[, P, drop = FALSE], D_post = D_post, vP = vP,
       pure = rownames(D_pure), pure_cohort = as.integer(pure_cohort),
       pure_period = as.integer(pure_period), b_P = b_P,
       zhat = as.numeric(D_pure[, P, drop = FALSE] %*% b_P),
       c0 = sum(vP * b_P) + sum(ell_Q * as.numeric(b[Q])))
}

fw_projector <- function(N) {
  if (!nrow(N)) return(diag(ncol(N)))
  q <- qr(t(N), tol = 1e-10)
  Q1 <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  diag(ncol(N)) - Q1 %*% t(Q1)
}

fw_aggregation <- function(fit, periods = NULL) {
  info <- fit$coordinate_info
  if (is.null(periods)) periods <- sort(unique(info$relative_period))
  nm <- paste0("delta_", ifelse(periods <= 0, "pre", "post"), "_g1_s", periods)
  A <- matrix(0, length(periods), nrow(info), dimnames = list(nm, info$coordinate))
  for (j in seq_along(periods)) {
    idx <- which(info$relative_period == periods[j])
    m <- as.numeric(fit$cohort_masses[as.character(info$cohort[idx])])
    A[j, idx] <- m / sum(m)
  }
  list(A = A, info = data.frame(coordinate = nm, type = ifelse(periods <= 0, "pre", "post"),
                                cohort = 1L, relative_period = as.integer(periods),
                                stringsAsFactors = FALSE))
}

fw_flat <- function(x) gsub("\\s+", " ", x)

fw_fit <- cohort_fit(y ~ d, data = independent_panel(), index = c("id", "time"))
fw_bounds <- c(0, 0.01, 0.05)
fw_full <- cohort_sensitivity(fw_fit, M = fw_bounds, method = "FLCI", seed = 1)
fw_post <- cohort_sensitivity(fw_fit, M = fw_bounds, method = "FLCI", full.path = FALSE, seed = 1)

test_that("flci_weights returns the documented columns, types, order, and attributes", {
  v <- flci_weights(fw_full)
  expect_s3_class(v, "data.frame")
  expect_identical(names(v), c("target", "label", "M", "cohort", "relative_period", "element",
                               "weight", "estimate", "contribution", "tie.break"))
  expect_identical(unname(vapply(v, function(z) class(z)[1L], character(1))),
                   c("character", "character", "numeric", "integer", "integer", "character",
                     "numeric", "numeric", "numeric", "logical"))
  expect_identical(rownames(v), as.character(seq_len(nrow(v))))
  pre <- fw_fit$coordinate_info$coordinate[fw_fit$coordinate_info$type == "pre"]
  expect_identical(nrow(v), length(fw_bounds) * length(pre))
  expect_identical(v$M, rep(fw_bounds, each = length(pre)))
  expect_identical(v$element[v$M == 0], pre)
  expect_identical(v$cohort[v$M == 0], c(4L, 4L, 4L, 6L, 6L, 6L, 6L, 6L))
  expect_identical(v$relative_period[v$M == 0], c(-2:0, -4:0))
  expect_identical(unique(v$target), "att")
  expect_identical(unique(v$label), fw_full$targets[[1]]$label)
  expect_equal(v$contribution, v$weight * v$estimate, tolerance = 0)
  expect_identical(attr(v, "weights"), "coefficients")
  expect_identical(attr(v, "by"), "element")
  expect_identical(attr(v, "framework"), "cohort")
  expect_true(attr(v, "full.path"))
  centers <- attr(v, "centers")
  expect_identical(names(centers), c("target", "label", "M", "estimate", "post.only.center", "center"))
  expect_identical(centers$M, fw_bounds)
  expect_equal(centers$estimate, rep(fw_full$estimates$estimate, 3L), tolerance = 0)
  expect_equal(centers$center, fw_full$intervals$center, tolerance = 0)
  missing <- attr(v, "missing")
  expect_identical(names(missing), c("target", "label", "M", "status"))
  expect_identical(nrow(missing), 0L)

  g <- flci_weights(fw_full, weights = "gamma", by = "cohort")
  expect_identical(names(g), c("target", "label", "M", "cohort", "elements", "weight", "l1",
                               "share", "contribution", "tie.break"))
  expect_identical(unname(vapply(g, function(z) class(z)[1L], character(1))),
                   c("character", "character", "numeric", "integer", "integer", "numeric",
                     "numeric", "numeric", "numeric", "logical"))
  expect_identical(g$cohort, rep(c(4L, 6L), 3L))
  expect_identical(g$elements, rep(c(1L, 3L), 3L))
  expect_identical(attr(g, "weights"), "gamma")
  expect_identical(attr(g, "by"), "cohort")
})

test_that("weights v follow the closed form and project onto the stored weights", {
  info <- fw_fit$coordinate_info
  ref <- fw_reference(info, fw_fit$betahat, fw_fit$W, fw_full$targets[[1]]$target)
  Pi <- fw_projector(as.matrix(fw_fit$normalization)[, ref$P, drop = FALSE])
  v <- flci_weights(fw_full)
  for (m in seq_along(fw_bounds)) {
    stored <- fw_full$results[[m]]$weights
    expected <- ref$vP + as.numeric(t(ref$D_pure) %*% as.numeric(stored$gamma))
    got <- v$weight[v$M == fw_bounds[m]]
    expect_equal(got, as.numeric(expected), tolerance = 1e-10)
    expect_equal(as.numeric(Pi %*% got), as.numeric(stored$coefficients[ref$P]), tolerance = 1e-10)
    expect_equal(v$estimate[v$M == fw_bounds[m]], ref$b_P, tolerance = 0)
  }
  # The post-only weights do not depend on M and load only on s = -1 and s = 0
  # here, because the imputation bias map W has no post-pre block.
  vp <- flci_weights(fw_post)
  for (m in fw_bounds) expect_equal(vp$weight[vp$M == m], as.numeric(ref$vP), tolerance = 1e-10)
  first <- vp[vp$M == 0, ]
  expect_true(all(first$weight[first$relative_period < -1L] == 0))
  expect_true(all(first$weight[first$relative_period >= -1L] != 0))
  # The post-only estimator's bias depends on the post-connected second
  # differences only: W'(ell + v) lies in the row space of those rows.
  full_v <- setNames(numeric(nrow(info)), info$coordinate)
  full_v[ref$P] <- first$weight
  bias <- t(fw_fit$W[info$coordinate, info$coordinate]) %*%
    (fw_full$targets[[1]]$target[info$coordinate] + full_v)
  resid <- qr.resid(qr(t(ref$D_post[, info$coordinate])), as.numeric(bias))
  expect_lt(max(abs(resid)), 1e-10)
})

test_that("weights gamma equal the stored gamma with second-difference labels", {
  ref <- fw_reference(fw_fit$coordinate_info, fw_fit$betahat, fw_fit$W, fw_full$targets[[1]]$target)
  g <- flci_weights(fw_full, weights = "gamma")
  expect_identical(g$element[g$M == 0], c("d2_g4_s0", "d2_g6_s-2", "d2_g6_s-1", "d2_g6_s0"))
  expect_identical(g$element[g$M == 0], ref$pure)
  expect_identical(g$cohort[g$M == 0], ref$pure_cohort)
  expect_identical(g$relative_period[g$M == 0], ref$pure_period)
  for (m in seq_along(fw_bounds)) {
    rows <- g[g$M == fw_bounds[m], ]
    expect_equal(rows$weight, as.numeric(fw_full$results[[m]]$weights$gamma), tolerance = 0)
    expect_equal(rows$estimate, ref$zhat, tolerance = 1e-12)
    expect_equal(rows$contribution, rows$weight * rows$estimate, tolerance = 0)
  }
  expect_true(any(g$weight[g$M == 0] != 0))
})

test_that("contributions add up to the center differences", {
  v <- flci_weights(fw_full)
  g <- flci_weights(fw_full, weights = "gamma")
  centers <- attr(g, "centers")
  ref <- fw_reference(fw_fit$coordinate_info, fw_fit$betahat, fw_fit$W, fw_full$targets[[1]]$target)
  expect_equal(centers$post.only.center, rep(ref$c0, 3L), tolerance = 1e-12)
  expect_equal(centers$post.only.center, fw_post$intervals$center, tolerance = 1e-12)
  expect_equal(as.numeric(tapply(g$contribution, g$M, sum)),
               centers$center - centers$post.only.center, tolerance = 1e-12)
  expect_equal(as.numeric(tapply(v$contribution, v$M, sum)),
               centers$center - centers$estimate, tolerance = 1e-12)
  vp <- flci_weights(fw_post)
  cp <- attr(vp, "centers")
  expect_equal(as.numeric(tapply(vp$contribution, vp$M, sum)), cp$center - cp$estimate,
               tolerance = 1e-12)
})

test_that("cohort sums add up each cohort's element rows", {
  for (w in c("coefficients", "gamma")) {
    e <- flci_weights(fw_full, weights = w)
    s <- flci_weights(fw_full, weights = w, by = "cohort")
    key <- paste(e$M, e$cohort)
    expect_equal(s$weight, as.numeric(tapply(e$weight, key, sum)[paste(s$M, s$cohort)]), tolerance = 1e-14)
    expect_equal(s$l1, as.numeric(tapply(abs(e$weight), key, sum)[paste(s$M, s$cohort)]), tolerance = 1e-14)
    expect_equal(s$contribution, as.numeric(tapply(e$contribution, key, sum)[paste(s$M, s$cohort)]),
                 tolerance = 1e-14)
    expect_identical(s$elements, as.integer(table(key)[paste(s$M, s$cohort)]))
    for (m in fw_bounds) {
      share <- s$share[s$M == m]
      if (all(s$l1[s$M == m] == 0)) expect_true(all(is.na(share))) else
        expect_equal(sum(share), 1, tolerance = 1e-14)
    }
  }
  # gamma is zero at the largest bound here, so its shares are undefined.
  s <- flci_weights(fw_full, weights = "gamma", by = "cohort")
  expect_true(all(is.na(s$share[s$M == 0.05])))
  expect_equal(sum(s$contribution[s$M == 0]),
               fw_full$intervals$center[1] - fw_post$intervals$center[1], tolerance = 1e-12)
})

test_that("CS-NYT weights include the bias map's post-pre block", {
  fit <- cohort_fit(y ~ d, data = independent_panel(), index = c("id", "time"), estimator = "cs_nyt")
  x <- cohort_sensitivity(fit, M = c(0, 0.01), method = "FLCI", seed = 1)
  ref <- fw_reference(fit$coordinate_info, fit$betahat, fit$W, x$targets[[1]]$target)
  Pi <- fw_projector(as.matrix(fit$normalization)[, ref$P, drop = FALSE])
  v <- flci_weights(x)
  for (m in 1:2) {
    expected <- ref$vP + as.numeric(t(ref$D_pure) %*% as.numeric(x$results[[m]]$weights$gamma))
    expect_equal(v$weight[v$M == x$settings$M[m]], as.numeric(expected), tolerance = 1e-10)
    expect_equal(as.numeric(Pi %*% expected),
                 as.numeric(x$results[[m]]$weights$coefficients[ref$P]), tolerance = 1e-10)
  }
  expect_true(ref$vP[["delta_pre_g6_s-2"]] != 0)
  # Imputation and CS-NYT give the same gamma and centers on this balanced panel.
  g_imp <- flci_weights(fw_full, weights = "gamma")
  g_nyt <- flci_weights(x, weights = "gamma")
  expect_equal(g_nyt$weight, g_imp$weight[g_imp$M <= 0.01], tolerance = 1e-6)
  expect_equal(attr(g_nyt, "centers")$center, fw_full$intervals$center[1:2], tolerance = 1e-8)
})

test_that("aggregate inference reports the aggregated path's weights", {
  x <- cohort_sensitivity(fw_fit, M = c(0, 0.01), method = "FLCI", framework = "aggregate", seed = 1)
  agg <- fw_aggregation(fw_fit)
  b <- setNames(as.numeric(agg$A %*% fw_fit$betahat), rownames(agg$A))
  ell <- setNames(as.numeric(qr.solve(t(agg$A), x$targets[[1]]$target, tol = 1e-10)), rownames(agg$A))
  W <- diag(nrow(agg$A)); dimnames(W) <- list(rownames(agg$A), rownames(agg$A))
  ref <- fw_reference(agg$info, b, W, ell)
  v <- flci_weights(x)
  g <- flci_weights(x, weights = "gamma")
  expect_true(all(is.na(v$cohort)) && is.integer(v$cohort))
  expect_identical(v$element[v$M == 0], paste0("delta_pre_g1_s", -4:0))
  expect_identical(g$element[g$M == 0], paste0("d2_g1_s", -2:0))
  expect_identical(g$relative_period[g$M == 0], -2:0)
  for (m in 1:2) {
    stored <- x$results[[m]]$weights
    expect_equal(g$weight[g$M == x$settings$M[m]], as.numeric(stored$gamma), tolerance = 0)
    expected <- ref$vP + as.numeric(t(ref$D_pure) %*% as.numeric(stored$gamma))
    expect_equal(v$weight[v$M == x$settings$M[m]], as.numeric(expected), tolerance = 1e-10)
  }
  centers <- attr(g, "centers")
  expect_equal(centers$post.only.center, rep(ref$c0, 2L), tolerance = 1e-12)
  expect_equal(as.numeric(tapply(g$contribution, g$M, sum)),
               centers$center - centers$post.only.center, tolerance = 1e-12)
  expect_equal(as.numeric(tapply(v$contribution, v$M, sum)),
               centers$center - centers$estimate, tolerance = 1e-12)
  expect_error(flci_weights(x, by = "cohort"),
               "by = 'cohort' requires cohort inference: aggregate weights have no cohorts.", fixed = TRUE)
})

test_that("balanced results keep the window's second differences", {
  win <- cohort_sensitivity(fw_fit, M = c(0, 0.01), method = "FLCI", balanced = TRUE,
                            cohort.paths = "window", seed = 1)
  L <- win$settings$window[["lower"]]
  S <- win$settings$window[["upper"]]
  expect_identical(c(L, S), c(-2L, 2L))
  info <- fw_fit$coordinate_info
  keep <- info$type == "pre" | info$relative_period <= S
  nm <- info$coordinate[keep]
  ref <- fw_reference(info[keep, ], fw_fit$betahat[nm], fw_fit$W[nm, nm],
                      win$targets[[1]]$target[nm], window_lower = L)
  g <- flci_weights(win, weights = "gamma")
  v <- flci_weights(win)
  expect_identical(g$element[g$M == 0], c("d2_g4_s0", "d2_g6_s0"))
  expect_identical(v$element[v$M == 0], info$coordinate[info$type == "pre"])
  for (m in 1:2) {
    stored <- win$results[[m]]$weights
    expect_equal(g$weight[g$M == win$settings$M[m]], as.numeric(stored$gamma), tolerance = 0)
    expected <- ref$vP + as.numeric(t(ref$D_pure) %*% as.numeric(stored$gamma))
    expect_equal(v$weight[v$M == win$settings$M[m]], as.numeric(expected), tolerance = 1e-10)
  }
  centers <- attr(g, "centers")
  expect_equal(as.numeric(tapply(g$contribution, g$M, sum)),
               centers$center - centers$post.only.center, tolerance = 1e-12)

  full_paths <- cohort_sensitivity(fw_fit, M = 0, method = "FLCI", balanced = TRUE, seed = 1)
  expect_identical(flci_weights(full_paths, weights = "gamma")$element,
                   c("d2_g4_s0", "d2_g6_s-2", "d2_g6_s-1", "d2_g6_s0"))

  ba <- cohort_sensitivity(fw_fit, M = c(0, 0.01), method = "FLCI", framework = "aggregate",
                           balanced = TRUE, seed = 1)
  a <- ba$aggregation
  W <- diag(length(a$betahat)); dimnames(W) <- list(names(a$betahat), names(a$betahat))
  ref <- fw_reference(a$coordinate_info, a$betahat, W, setNames(a$targets[1, ], names(a$betahat)))
  v <- flci_weights(ba)
  g <- flci_weights(ba, weights = "gamma")
  expect_identical(g$element[g$M == 0], "d2_g1_s0")
  expect_identical(v$element[v$M == 0], paste0("delta_pre_g1_s", -2:0))
  expect_identical(nrow(a$normalization), 0L)
  for (m in 1:2) {
    stored <- ba$results[[m]]$weights
    expected <- ref$vP + as.numeric(t(ref$D_pure) %*% as.numeric(stored$gamma))
    expect_equal(v$weight[v$M == ba$settings$M[m]], as.numeric(expected), tolerance = 1e-10)
    expect_equal(v$weight[v$M == ba$settings$M[m]], as.numeric(stored$coefficients[ref$P]),
                 tolerance = 1e-10)
  }
})

test_that("post-only results have zero gamma and short paths have no gamma rows", {
  g <- flci_weights(fw_post, weights = "gamma")
  expect_identical(nrow(g), 12L)
  expect_true(all(g$weight == 0))
  expect_identical(attr(g, "full.path"), FALSE)
  expect_message(r <- plot(fw_post, type = "weights", weights = "gamma"),
    "The post-only FLCI puts no weight on the pure pre-treatment second differences: its weights gamma are zero by construction. Use weights = 'coefficients' to plot its weights v.",
    fixed = TRUE)
  expect_null(r)
  expect_s3_class(plot(fw_post, type = "weights"), "ggplot")

  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  x <- cohort_sensitivity(fit, M = c(0, 0.01), method = "FLCI", seed = 1)
  expect_length(x$results[[1]]$weights$gamma, 0L)
  for (by in c("element", "cohort")) {
    z <- flci_weights(x, weights = "gamma", by = by)
    expect_identical(nrow(z), 0L)
    expect_identical(names(z), names(flci_weights(fw_full, weights = "gamma", by = by)))
    expect_identical(unname(vapply(z, function(y) class(y)[1L], character(1))),
                     unname(vapply(flci_weights(fw_full, weights = "gamma", by = by),
                                   function(y) class(y)[1L], character(1))))
    expect_identical(nrow(attr(z, "centers")), 2L)
  }
  expect_message(r <- plot(x, type = "weights", weights = "gamma"),
    "This result has no pure pre-treatment second differences, so the FLCI has no weights gamma to plot.",
    fixed = TRUE)
  expect_null(r)
  v <- flci_weights(x)
  expect_identical(v$element, rep(c("delta_pre_g3_s-1", "delta_pre_g3_s0"), 2L))
  expect_equal(v$weight, rep(c(1, -2), 2L), tolerance = 1e-12)
})

test_that("results without weights are listed as missing", {
  short <- cohort_fit(y ~ d, data = subset(independent_two_cohort_short_panel(), time >= 2),
                      index = c("id", "time"))
  x <- suppressWarnings(cohort_sensitivity(short, M = c(0, 0.01), method = "FLCI", seed = 1))
  expect_false(any(x$intervals$success))
  z <- flci_weights(x)
  expect_identical(nrow(z), 0L)
  expect_identical(attr(z, "missing")$status, rep("insufficient_prehistory", 2L))
  expect_identical(attr(z, "missing")$M, c(0, 0.01))
  expect_identical(nrow(attr(z, "centers")), 0L)
  expect_error(plot(x, type = "weights"),
               "No FLCI weights to plot: every selected result failed; inspect summary().", fixed = TRUE)

  part <- fw_full
  part$results[[2]]$weights <- NULL
  part$results[[2]]$status <- "inference_error"
  z <- flci_weights(part)
  expect_identical(unique(z$M), c(0, 0.05))
  expect_identical(attr(z, "missing")$M, 0.01)
  expect_identical(attr(z, "missing")$status, "inference_error")
  expect_identical(attr(z, "centers")$M, c(0, 0.05))
  p <- plot(part, type = "weights")
  expect_match(fw_flat(p$labels$caption),
               "1 failed result(s) without weights are not shown; inspect summary().", fixed = TRUE)
  expect_identical(levels(p$data$column), c("0", "0.05"))
  expect_error(plot(part, type = "weights", M = 0.01),
               "No FLCI weights to plot: every selected result failed; inspect summary().", fixed = TRUE)
})

test_that("the minimum-norm tie-break is marked in the table and the caption", {
  info <- fw_fit$coordinate_info
  pre <- info$type == "pre"
  Sigma <- fw_fit$Sigma
  Sigma[pre, ] <- 0
  Sigma[, pre] <- 0
  input <- cohort_input(betahat = coef(fw_fit), vcov = Sigma, coefficient_info = info,
                        W = fw_fit$W, normalization = fw_fit$normalization,
                        cohort_sizes = fw_fit$cohort_masses, estimator = "imp")
  x <- cohort_sensitivity(input, M = c(0, 0.01), method = "FLCI", seed = 1)
  expect_identical(x$results[[1]]$diagnostics$tie.break, "minimum_euclidean_norm_gamma")
  g <- flci_weights(x, weights = "gamma")
  expect_true(all(g$tie.break))
  expect_true(all(g$weight == 0))
  expect_true(all(flci_weights(x, by = "cohort")$tie.break))
  p <- plot(x, type = "weights", weights = "gamma")
  caption <- fw_flat(p$labels$caption)
  expect_match(caption, "At M = 0, 0.01, several weight vectors give the same FLCI length; the one with the smallest Euclidean norm of gamma is shown.", fixed = TRUE)
  expect_match(caption, "All weights gamma are zero at the selected bounds: the full-path FLCI equals the post-only FLCI.", fixed = TRUE)
  expect_false(any(flci_weights(fw_full)$tie.break))
  expect_false(grepl("several weight vectors", fw_flat(plot(fw_full, type = "weights")$labels$caption)))
})

test_that("errors and messages follow the documented texts", {
  expect_error(flci_weights(fw_full$intervals), "x must be a cohort_sensitivity object.", fixed = TRUE)
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  clf <- cohort_sensitivity(fit, M = 0.01, method = "C-LF", seed = 1,
                            control = list(lf.sims = 200L, grid.points = 21L))
  cond <- cohort_sensitivity(fit, M = 0.01, method = "Conditional", seed = 1,
                             control = list(grid.points = 21L))
  msg_clf <- "FLCI weights exist only for results with method = 'FLCI'; this result used C-LF."
  msg_cond <- "FLCI weights exist only for results with method = 'FLCI'; this result used Conditional."
  expect_error(flci_weights(clf), msg_clf, fixed = TRUE)
  expect_error(plot(clf, type = "weights"), msg_clf, fixed = TRUE)
  expect_error(flci_weights(cond), msg_cond, fixed = TRUE)
  expect_error(plot(cond, type = "weights"), msg_cond, fixed = TRUE)
  expect_error(flci_weights(fw_full, weights = "v"), "should be one of")
  expect_error(flci_weights(fw_full, by = "path"), "should be one of")
  expect_error(plot(fw_full, type = "weights", weights = "v"), "should be one of")
  expect_error(plot(fw_full, type = "weights", by = "path"), "should be one of")
  expect_error(plot(fw_full, type = "heatmap"), "should be one of")

  agg <- cohort_sensitivity(fw_fit, M = 0.01, method = "FLCI", framework = "aggregate", seed = 1)
  expect_error(plot(agg, type = "weights", by = "cohort"),
               "by = 'cohort' requires cohort inference: aggregate weights have no cohorts.", fixed = TRUE)
  expect_error(plot(agg, type = "weights", cohorts = 4),
               "cohorts cannot be selected under aggregate inference: aggregate weights have no cohorts.",
               fixed = TRUE)
  none <- "No computed results match the requested plot selection."
  expect_error(plot(fw_full, type = "weights", M = 0.3), none, fixed = TRUE)
  expect_error(plot(fw_full, type = "weights", targets = "s1"), none, fixed = TRUE)
  expect_error(plot(fw_full, type = "weights", cohorts = 5), none, fixed = TRUE)
  expect_error(plot(fw_full, weights = "gamma"),
               "weights, by, targets, and labels apply only to type = 'weights'.", fixed = TRUE)
  expect_error(plot(fw_full, labels = TRUE),
               "weights, by, targets, and labels apply only to type = 'weights'.", fixed = TRUE)
  expect_error(plot(fw_full, targets = "att"),
               "weights, by, targets, and labels apply only to type = 'weights'.", fixed = TRUE)
  expect_error(plot(fw_full, by = "cohort"),
               "weights, by, targets, and labels apply only to type = 'weights'.", fixed = TRUE)
  expect_error(plot(fw_full, type = "weights", identified.set = FALSE),
               "identified.set and color apply only to type = 'sets'.", fixed = TRUE)
  expect_error(plot(fw_full, type = "weights", color = "red"),
               "identified.set and color apply only to type = 'sets'.", fixed = TRUE)
  for (bad in list(NA, c(TRUE, FALSE), "yes", 1))
    expect_error(plot(fw_full, type = "weights", labels = bad), "labels must be TRUE or FALSE.", fixed = TRUE)
  for (bad in list(1, character(), NA_character_, c("att", NA)))
    expect_error(plot(fw_full, type = "weights", targets = bad),
                 "targets must be NULL or a character vector of target ids.", fixed = TRUE)
  for (bad in list(Inf, NA_real_, "0"))
    expect_error(plot(fw_full, type = "weights", M = bad), "M must select finite computed bounds.", fixed = TRUE)

  broken <- fw_full
  broken$results[[1]]$weights$coefficients[["delta_pre_g6_s-1"]] <- 1
  expect_error(flci_weights(broken), "flci_weights() could not reconstruct the FLCI weights of this result.",
               fixed = TRUE)
  broken <- fw_full
  broken$results[[2]]$weights$gamma <- broken$results[[2]]$weights$gamma[-1]
  expect_error(flci_weights(broken, weights = "gamma"),
               "flci_weights() could not reconstruct the FLCI weights of this result.", fixed = TRUE)
  broken <- fw_full
  names(broken$results[[1]]$weights$gamma) <- rev(names(broken$results[[1]]$weights$gamma))
  expect_error(flci_weights(broken), "flci_weights() could not reconstruct the FLCI weights of this result.",
               fixed = TRUE)
  expect_silent(flci_weights(fw_full))
  expect_silent(flci_weights(fw_post, weights = "gamma"))
})

test_that("the default plot is unchanged and the weights plot draws one tile per row", {
  default <- plot(fw_full)
  sets <- plot(fw_full, type = "sets")
  expect_equal(ggplot2::ggplot_build(default)$data, ggplot2::ggplot_build(sets)$data)
  expect_identical(default$labels$caption, sets$labels$caption)
  expect_false(any(vapply(default$layers, function(z) inherits(z$geom, "GeomTile"), logical(1))))

  p <- plot(fw_full, type = "weights")
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$layers[[1]]$geom, "GeomTile"))
  expect_length(p$layers, 1L)
  table <- flci_weights(fw_full)
  expect_identical(nrow(p$data), nrow(table))
  expect_equal(p$data$value, table$weight, tolerance = 0)
  expect_identical(levels(p$data$column), c("0", "0.01", "0.05"))
  expect_identical(levels(p$data$row), as.character(0:-4))
  expect_identical(levels(p$data$cohort_label), c("Cohort 4", "Cohort 6"))
  built <- ggplot2::ggplot_build(p)
  tiles <- built$data[[1]]
  expect_identical(nrow(tiles), nrow(p$data))
  expect_true(all(toupper(tiles$fill[p$data$value == 0]) == "#FFFFFF"))
  expect_true(any(p$data$value == 0))
  visible <- abs(p$data$value) > 0.01 * max(abs(p$data$value))
  expect_true(all(toupper(tiles$fill[visible]) != "#FFFFFF"))
  expect_identical(p$labels$x, "M")
  expect_identical(p$labels$y, "Relative period")
  expect_identical(p$labels$fill, "Weight v")
  caption <- fw_flat(p$labels$caption)
  expect_identical(caption, paste("FLCI weights v on the pre-treatment coefficients; SD, full path;",
                                  "cohort inference. White cells are zero weights."))

  q <- plot(fw_full, type = "weights", weights = "gamma", labels = TRUE, M = c(0, 0.01),
            main = "Gamma", xlab = "Bound", ylab = "Period")
  expect_length(q$layers, 2L)
  expect_true(inherits(q$layers[[2]]$geom, "GeomText"))
  expect_identical(q$data$text[q$data$value == 0], rep("", sum(q$data$value == 0)))
  expect_identical(q$data$text[q$data$value != 0],
                   trimws(formatC(q$data$value[q$data$value != 0], format = "fg", digits = 2)))
  expect_identical(nrow(ggplot2::ggplot_build(q)$data[[1]]), 8L)
  expect_identical(q$labels$title, "Gamma")
  expect_identical(q$labels$x, "Bound")
  expect_identical(q$labels$y, "Period")
  expect_identical(q$labels$fill, "Weight gamma")
  expect_identical(fw_flat(q$labels$caption),
                   paste("FLCI weights gamma on the pure pre-treatment second differences; SD, full path;",
                         "cohort inference. White cells are zero weights."))
  expect_identical(plot(fw_full, type = "weights", weights = "gamma")$labels$y,
                   "Last period of the second difference")

  r <- plot(fw_full, type = "weights", by = "cohort")
  expect_equal(r$data$value, flci_weights(fw_full, by = "cohort")$contribution, tolerance = 0)
  expect_identical(levels(r$data$row), c("Cohort 6", "Cohort 4"))
  expect_identical(r$labels$y, "Cohort")
  expect_identical(r$labels$fill, "Contribution")
  expect_identical(fw_flat(r$labels$caption), paste(
    "Sum over each cohort's pre-treatment coefficients of v times the estimate; SD, full path;",
    "cohort inference. White cells are zero. The cohort sums add up to the FLCI center minus the",
    "conventional estimate."))
  r <- plot(fw_full, type = "weights", weights = "gamma", by = "cohort", M = 0.05)
  expect_identical(fw_flat(r$labels$caption), paste(
    "Sum over each cohort's second differences of gamma times the estimate; SD, full path;",
    "cohort inference. White cells are zero. The cohort sums add up to the FLCI center minus the",
    "post-only center. All weights gamma are zero at the selected bounds: the full-path FLCI",
    "equals the post-only FLCI."))
  expect_identical(fw_flat(plot(fw_post, type = "weights")$labels$caption),
                   paste("FLCI weights v on the pre-treatment coefficients; SD, post-only;",
                         "cohort inference. White cells are zero weights."))

  s <- plot(fw_full, type = "weights", cohorts = 6)
  expect_identical(unique(s$data$cohort), 6L)
  expect_identical(nrow(s$data), 15L)
})

test_that("weights plots select targets and describe the panel and window", {
  dyn <- cohort_sensitivity(fw_fit, M = c(0, 0.01), method = "FLCI", target = "dynamic", seed = 1)
  p <- plot(dyn, type = "weights")
  expect_identical(levels(p$data$target_label), c("Relative period 1", "Relative period 2",
                                                   "Relative period 3", "Relative period 4"))
  q <- plot(dyn, type = "weights", targets = c("s3", "s1"), M = 0.01)
  expect_identical(unique(q$data$target), c("s1", "s3"))
  expect_identical(levels(q$data$target_label), c("Relative period 1", "Relative period 3"))
  expect_identical(unique(q$data$M), 0.01)

  agg <- cohort_sensitivity(fw_fit, M = c(0, 0.01), method = "FLCI", framework = "aggregate", seed = 1)
  a <- plot(agg, type = "weights")
  expect_false("cohort_label" %in% names(a$data))
  expect_match(fw_flat(a$labels$caption), "SD, full path; aggregate inference.", fixed = TRUE)

  win <- cohort_sensitivity(fw_fit, M = 0.01, method = "FLCI", balanced = TRUE,
                            cohort.paths = "window", seed = 1)
  expect_match(fw_flat(plot(win, type = "weights")$labels$caption),
               "Balanced window: relative periods -2 to 2; cohort paths: window.", fixed = TRUE)
  ba <- cohort_sensitivity(fw_fit, M = 0.01, method = "FLCI", framework = "aggregate",
                           balanced = TRUE, seed = 1)
  expect_match(fw_flat(plot(ba, type = "weights")$labels$caption),
               "Balanced window: relative periods -2 to 2.", fixed = TRUE)

  d <- independent_panel()
  d <- d[!(d$id == 3 & d$time == 2), ]
  fit <- suppressWarnings(cohort_fit(y ~ d, data = d, index = c("id", "time")))
  x <- cohort_sensitivity(fit, M = 0.01, method = "FLCI", seed = 1)
  expect_match(fw_flat(plot(x, type = "weights")$labels$caption), "Unbalanced-panel approximation.",
               fixed = TRUE)
  expect_false(grepl("Unbalanced", fw_flat(plot(fw_full, type = "weights")$labels$caption)))
})
