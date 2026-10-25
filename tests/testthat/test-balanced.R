# Balanced targets, balanced aggregation, and cohort.paths. Reference values
# are computed here from the data and from the fits' public fields.

balanced_panel_fit <- function(estimator = "cs_nt", weights = NULL) {
  cohort_fit(y ~ d, data = independent_panel(), index = c("id", "time"),
             estimator = estimator, weights = weights)
}

balanced_panel_mass <- function(d) {
  c(`4` = sum(d$w[d$time == 1 & d$id <= 20]),
    `6` = sum(d$w[d$time == 1 & d$id > 20 & d$id <= 42]))
}

test_that("the balanced window and options are recorded in settings", {
  fit <- balanced_panel_fit()
  info <- fit$coordinate_info
  # Relative periods observed for both cohorts (adoption in periods 4 and 6 of 7).
  common <- Reduce(intersect, split(info$relative_period, info$cohort))
  expect_equal(sort(common), -2:2)
  att <- cohort_sensitivity(fit, M = 0.01, method = "FLCI", balanced = TRUE, seed = 1)
  expect_identical(att$settings$window, c(lower = -2L, upper = 2L))
  expect_true(att$settings$balanced)
  expect_identical(att$settings$cohort.paths, "full")
  expect_identical(tail(names(att$settings), 3L), c("balanced", "cohort.paths", "window"))
  win <- cohort_sensitivity(fit, M = 0.01, method = "FLCI", balanced = TRUE,
                            cohort.paths = "window", seed = 1)
  expect_identical(win$settings$cohort.paths, "window")
  agg <- cohort_sensitivity(fit, M = 0.01, method = "FLCI", framework = "aggregate",
                            balanced = TRUE, seed = 1)
  expect_identical(agg$settings$cohort.paths, NA_character_)
  expect_identical(agg$settings$window, c(lower = -2L, upper = 2L))
  short <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  one <- cohort_sensitivity(short, M = 0.01, method = "FLCI", balanced = TRUE, seed = 1)
  expect_identical(one$settings$window, c(lower = -1L, upper = 1L))
})

test_that("balanced targets weight cohorts by mass on periods 1 to S", {
  d <- independent_panel()
  fit <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nt", weights = "w")
  info <- fit$coordinate_info
  mass <- balanced_panel_mass(d)
  m <- unname(mass[as.character(info$cohort)])
  S <- 2
  inside <- info$type == "post" & info$relative_period <= S
  check <- function(result, weights) {
    expect_length(result$targets, length(weights))
    for (k in seq_along(weights)) {
      l <- result$targets[[k]]$target
      expect_identical(names(l), info$coordinate)
      expect_equal(unname(l), weights[[k]], tolerance = 1e-12)
      expect_equal(result$estimates$estimate[k], sum(weights[[k]] * fit$betahat), tolerance = 1e-12)
      expect_equal(result$estimates$se[k],
                   sqrt(drop(t(weights[[k]]) %*% fit$Sigma %*% weights[[k]])), tolerance = 1e-12)
    }
  }
  run <- function(target) cohort_sensitivity(fit, M = 0.01, target = target, method = "FLCI",
                                             full.path = FALSE, balanced = TRUE, seed = 2)
  att <- run("att")
  check(att, list(ifelse(inside, m / (S * sum(mass)), 0)))
  expect_identical(att$targets[[1]]$id, "att")
  expect_identical(att$estimates$label, "Average treatment effect, relative periods 1 to 2")
  dyn <- run("dynamic")
  check(dyn, lapply(seq_len(S), function(s)
    ifelse(info$type == "post" & info$relative_period == s, m / sum(mass), 0)))
  expect_identical(dyn$estimates$target, c("s1", "s2"))
  expect_equal(dyn$estimates$relative_period, 1:2)
  expect_identical(dyn$estimates$label, c("Relative period 1", "Relative period 2"))
  cel <- run("cohort-dynamic")
  cells <- which(inside)
  check(cel, lapply(cells, function(j) as.numeric(seq_len(nrow(info)) == j)))
  expect_identical(cel$estimates$target,
                   paste0("g", info$cohort[cells], "_s", info$relative_period[cells]))
  for (x in list(att, dyn, cel)) expect_true(all(x$intervals$success))
  expect_true(all(cel$intervals$relative_period <= S))
  short <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  one <- cohort_sensitivity(short, M = 0.01, method = "FLCI", balanced = TRUE, seed = 2)
  expect_identical(one$estimates$label, "Average treatment effect, relative period 1")
})

test_that("balanced aggregation forms the event study on the window with fixed weights", {
  d <- independent_panel()
  fit <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "imp", weights = "w")
  info <- fit$coordinate_info
  mass <- balanced_panel_mass(d)
  periods <- -2:2
  A <- matrix(0, length(periods), nrow(info), dimnames = list(
    paste0("delta_", ifelse(periods <= 0, "pre", "post"), "_g1_s", periods), info$coordinate))
  for (k in seq_along(periods)) {
    j <- which(info$relative_period == periods[k])
    A[k, j] <- mass[as.character(info$cohort[j])] / sum(mass)
  }
  for (target in c("att", "dynamic")) {
    agg <- cohort_sensitivity(fit, M = 0.01, target = target, framework = "aggregate",
                              method = "FLCI", full.path = FALSE, balanced = TRUE, seed = 3)
    coh <- cohort_sensitivity(fit, M = 0.01, target = target, method = "FLCI",
                              full.path = FALSE, balanced = TRUE, seed = 3)
    a <- agg$aggregation
    expect_identical(names(a), c("aggregation_matrix", "betahat", "Sigma", "normalization",
                                 "coordinate_info", "targets"))
    expect_equal(a$aggregation_matrix, A, tolerance = 1e-12)
    expect_equal(a$betahat, setNames(as.numeric(A %*% fit$betahat), rownames(A)), tolerance = 1e-12)
    expect_equal(a$Sigma, A %*% fit$Sigma %*% t(A), tolerance = 1e-12)
    expect_identical(a$coordinate_info$coordinate, rownames(A))
    expect_identical(a$coordinate_info$type, ifelse(periods <= 0, "pre", "post"))
    expect_equal(a$coordinate_info$relative_period, periods)
    expect_identical(colnames(a$normalization), rownames(A))
    if (nrow(a$normalization)) expect_lt(max(abs(a$normalization %*% a$betahat)), 1e-12)
    expect_identical(rownames(a$targets), vapply(agg$targets, function(z) z$id, character(1)))
    expect_identical(colnames(a$targets), rownames(A))
    for (k in seq_along(agg$targets)) {
      expect_equal(as.numeric(t(A) %*% a$targets[k, ]), unname(coh$targets[[k]]$target),
                   tolerance = 1e-10)
    }
    expect_equal(agg$estimates$estimate, coh$estimates$estimate, tolerance = 1e-12)
    expect_equal(agg$estimates$se, coh$estimates$se, tolerance = 1e-12)
    expect_true(all(agg$intervals$success))
  }
  att <- cohort_sensitivity(fit, M = 0.01, framework = "aggregate", method = "FLCI",
                            balanced = TRUE, seed = 3)
  # The balanced ATT averages the window's post-treatment periods equally.
  expect_equal(unname(att$aggregation$targets["att", ]), c(0, 0, 0, 0.5, 0.5), tolerance = 1e-12)
})

test_that("under post-only SD, full and window cohort paths give the same sets", {
  fit <- balanced_panel_fit("imp")
  for (target in c("att", "dynamic")) {
    full <- cohort_sensitivity(fit, M = c(0, 0.01), target = target, method = "FLCI",
                               full.path = FALSE, balanced = TRUE, seed = 4)
    window <- cohort_sensitivity(fit, M = c(0, 0.01), target = target, method = "FLCI",
                                 full.path = FALSE, balanced = TRUE, cohort.paths = "window", seed = 4)
    expect_true(all(full$intervals$success) && all(window$intervals$success))
    expect_equal(window$intervals$plugin.lower, full$intervals$plugin.lower, tolerance = 1e-9)
    expect_equal(window$intervals$plugin.upper, full$intervals$plugin.upper, tolerance = 1e-9)
    expect_equal(ordered_endpoints(window), ordered_endpoints(full), tolerance = 1e-7)
  }
  ctl <- list(lf.sims = 200L, grid.points = 21L)
  full <- cohort_sensitivity(fit, M = 0.01, balanced = TRUE, seed = 4, control = ctl)
  window <- cohort_sensitivity(fit, M = 0.01, balanced = TRUE, cohort.paths = "window",
                               seed = 4, control = ctl)
  expect_true(all(full$intervals$success) && all(window$intervals$success))
  expect_equal(window$intervals$plugin.lower, full$intervals$plugin.lower, tolerance = 1e-9)
  expect_equal(window$intervals$plugin.upper, full$intervals$plugin.upper, tolerance = 1e-9)
})

test_that("window paths use only the window's pre-treatment changes", {
  fit <- balanced_panel_fit("imp")
  # RM components: prod over cohorts of 2 * (pre-period first differences).
  # Full paths: cohort 4 has 3 pre periods, cohort 6 has 5; the window keeps 3 each.
  ctl <- list(lf.sims = 64L, grid.points = 21L, max.evaluations = 5000L)
  full <- cohort_sensitivity(fit, restriction = "RM", M = 0.5, balanced = TRUE, seed = 5, control = ctl)
  window <- cohort_sensitivity(fit, restriction = "RM", M = 0.5, balanced = TRUE,
                               cohort.paths = "window", seed = 5, control = ctl)
  expect_true(all(full$intervals$success) && all(window$intervals$success))
  expect_identical(full$results[[1]]$diagnostics$signed.components, (2L * 2L) * (2L * 4L))
  expect_identical(window$results[[1]]$diagnostics$signed.components, (2L * 2L) * (2L * 2L))
  # The full-path SD restriction drops cohort 6's second differences before -2.
  a <- cohort_sensitivity(fit, M = 0, method = "FLCI", balanced = TRUE, seed = 5)
  b <- cohort_sensitivity(fit, M = 0, method = "FLCI", balanced = TRUE, cohort.paths = "window", seed = 5)
  expect_true(all(a$intervals$success) && all(b$intervals$success))
  expect_gt(max(abs(ordered_endpoints(a) - ordered_endpoints(b))), 1e-6)
})

test_that("on a balanced panel imputation and CS-NYT agree under both cohort paths", {
  imp <- balanced_panel_fit("imp")
  nyt <- balanced_panel_fit("cs_nyt")
  for (paths in c("full", "window")) {
    for (target in c("att", "dynamic")) {
      a <- cohort_sensitivity(imp, M = c(0, 0.01), target = target, method = "FLCI",
                              balanced = TRUE, cohort.paths = paths, seed = 31)
      b <- cohort_sensitivity(nyt, M = c(0, 0.01), target = target, method = "FLCI",
                              balanced = TRUE, cohort.paths = paths, seed = 31)
      expect_true(all(a$intervals$success) && all(b$intervals$success))
      expect_equal(ordered_endpoints(a), ordered_endpoints(b), tolerance = 1e-7)
    }
  }
})

test_that("a single-cohort fit gives the same numbers with and without balanced", {
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  columns <- c("M", "component", "lower", "upper", "center", "plugin.lower", "plugin.upper",
               "success", "status", "message")
  ctl <- list(lf.sims = 200L, grid.points = 21L)
  cases <- list(list(method = "FLCI"), list(method = "C-LF", control = ctl),
                list(method = "FLCI", framework = "aggregate"),
                list(method = "C-LF", control = ctl, target = "dynamic", framework = "aggregate"))
  for (args in cases) {
    run <- function(...) do.call(cohort_sensitivity,
                                 c(list(fit = fit, M = c(0, 0.01), seed = 6), args, list(...)))
    base <- run()
    expect_true(all(base$intervals$success))
    variants <- list(run(balanced = TRUE))
    if (is.null(args$framework)) variants <- c(variants, list(run(balanced = TRUE, cohort.paths = "window")))
    for (x in variants) {
      expect_identical(x$intervals[columns], base$intervals[columns])
      expect_identical(x$estimates[c("estimate", "se", "lower", "upper")],
                       base$estimates[c("estimate", "se", "lower", "upper")])
    }
  }
})

test_that("cohort_input objects support balanced targets", {
  source_fit <- balanced_panel_fit()
  external <- cohort_input(betahat = source_fit$betahat, vcov = source_fit$Sigma,
    coefficient_info = source_fit$coordinate_info, W = source_fit$W,
    normalization = source_fit$normalization, cohort_sizes = source_fit$cohort_sizes,
    estimator = "cs_nt")
  for (paths in c("full", "window")) {
    a <- cohort_sensitivity(source_fit, M = 0.01, method = "FLCI", balanced = TRUE,
                            cohort.paths = paths, seed = 12)
    b <- cohort_sensitivity(external, M = 0.01, method = "FLCI", balanced = TRUE,
                            cohort.paths = paths, seed = 12)
    expect_identical(b$settings$window, a$settings$window)
    expect_equal(ordered_endpoints(b), ordered_endpoints(a), tolerance = 1e-9)
  }
})

test_that("balanced options reject invalid values and combinations", {
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"),
                    estimator = "cs_nt")
  run <- function(x = fit, ...) cohort_sensitivity(x, M = 0.01, method = "FLCI", seed = 7, ...)
  expect_error(cohort_sensitivity(fit, M = -1, balanced = NA), "nonnegative")
  for (bad in list(NA, "yes", c(TRUE, FALSE), 1, logical()))
    expect_error(run(balanced = bad), "balanced must be TRUE or FALSE.", fixed = TRUE)
  expect_error(run(cohort.paths = "full"), "cohort.paths applies only when balanced = TRUE.", fixed = TRUE)
  expect_error(run(balanced = FALSE, cohort.paths = "window"),
               "cohort.paths applies only when balanced = TRUE.", fixed = TRUE)
  expect_error(run(balanced = TRUE, framework = "aggregate", cohort.paths = "full"),
               "cohort.paths applies only to framework = 'cohort'.", fixed = TRUE)
  for (bad in list("win", "Full", NA_character_, c("full", "window"), 1, TRUE))
    expect_error(run(balanced = TRUE, cohort.paths = bad), "cohort.paths must be 'full' or 'window'.",
                 fixed = TRUE)
  source_fit <- balanced_panel_fit()
  subset_input <- function(drop = character(), W = source_fit$W) {
    info <- source_fit$coordinate_info
    keep <- setdiff(info$coordinate, drop)
    N <- source_fit$normalization[, keep, drop = FALSE]
    N <- N[rowSums(abs(N)) > 0, , drop = FALSE]
    cohort_input(betahat = source_fit$betahat[keep], vcov = source_fit$Sigma[keep, keep],
      coefficient_info = info[info$coordinate %in% keep, , drop = FALSE],
      W = W[keep, keep], normalization = N, cohort_sizes = source_fit$cohort_sizes,
      estimator = "cs_nt")
  }
  expect_error(run(subset_input("delta_post_g6_s1"), balanced = TRUE),
    "balanced = TRUE requires relative period 1 for every cohort: no post-treatment period is observed for all cohorts.",
    fixed = TRUE)
  expect_error(run(subset_input("delta_pre_g6_s0"), balanced = TRUE),
    "balanced = TRUE requires relative period 0 for every cohort: no pre-treatment period is observed for all cohorts.",
    fixed = TRUE)
  W <- source_fit$W
  # A coefficient inside the window (period 1) loaded on one after it (period 3).
  W["delta_post_g4_s1", "delta_post_g4_s3"] <- 0.5
  linked <- subset_input(W = W)
  expect_error(run(linked, balanced = TRUE, cohort.paths = "window"),
    "cohort.paths = 'window' requires that the bias map W and the normalization link no coefficient inside the balanced window to a coefficient after it.",
    fixed = TRUE)
  expect_true(all(run(linked, balanced = TRUE)$intervals$success))
})

test_that("print and summary state the balanced window only for balanced results", {
  fit <- balanced_panel_fit()
  coh <- cohort_sensitivity(fit, M = 0.01, method = "FLCI", balanced = TRUE,
                            cohort.paths = "window", seed = 8)
  out <- capture.output(print(coh))
  expect_match(out[3], "^Covariance: ")
  expect_identical(out[4], "Balanced window: relative periods -2 to 2; target periods 1 to 2; cohort paths: window")
  sm <- capture.output(print(summary(coh)))
  expect_match(sm[1], "^Covariance: ")
  expect_identical(sm[2], "Balanced window: relative periods -2 to 2; target periods 1 to 2; cohort paths: window")
  expect_identical(sm[3], "")
  expect_identical(sm[4], "Conventional estimates")
  agg <- cohort_sensitivity(fit, M = 0.01, method = "FLCI", framework = "aggregate",
                            balanced = TRUE, seed = 8)
  expect_identical(capture.output(print(agg))[4],
    "Balanced window: relative periods -2 to 2; target periods 1 to 2; aggregate path on the window")
  short <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  one <- cohort_sensitivity(short, M = 0.01, method = "FLCI", balanced = TRUE, seed = 8)
  expect_identical(capture.output(print(one))[4],
    "Balanced window: relative periods -1 to 1; target period 1; cohort paths: full")
  plain <- cohort_sensitivity(fit, M = 0.01, method = "FLCI", seed = 8)
  explicit <- cohort_sensitivity(fit, M = 0.01, method = "FLCI", seed = 8, balanced = FALSE)
  expect_identical(names(plain), c("call", "intervals", "estimates", "fit", "targets",
                                   "results", "settings", "diagnostics"))
  expect_identical(names(plain$settings), c("restriction", "method", "full.path", "target",
    "framework", "M", "conf.level", "parallel", "cores", "seed", "task.seeds", "control"))
  expect_identical(plain[-1L], explicit[-1L])
  expect_identical(names(coh), names(plain))
  expect_identical(names(agg), c(names(plain), "aggregation"))
  expect_false(any(grepl("Balanced", capture.output(print(plain)))))
  sm <- capture.output(print(summary(plain)))
  expect_match(sm[1], "^Covariance: ")
  expect_identical(sm[2], "")
  expect_identical(sm[3], "Conventional estimates")
})

test_that("balanced plots state the window and draw conventional points on it", {
  fit <- balanced_panel_fit()
  dyn <- cohort_sensitivity(fit, M = 0.01, target = "dynamic", method = "FLCI",
                            full.path = FALSE, balanced = TRUE, seed = 9)
  p <- plot(dyn)
  caption <- gsub("\n", " ", p$labels$caption, fixed = TRUE)
  expect_match(caption, "Balanced window: relative periods -2 to 2; cohort paths: full.", fixed = TRUE)
  expect_equal(p$data$relative_period, -2:2)
  agg <- cohort_sensitivity(fit, M = 0.01, target = "dynamic", framework = "aggregate",
                            method = "FLCI", full.path = FALSE, balanced = TRUE, seed = 9)
  p <- plot(agg)
  caption <- gsub("\n", " ", p$labels$caption, fixed = TRUE)
  expect_match(caption, "Balanced window: relative periods -2 to 2.", fixed = TRUE)
  expect_false(grepl("cohort paths", caption, fixed = TRUE))
  # With a fixed cohort mix, the conventional points are the aggregated coefficients.
  expect_equal(p$data$estimate, unname(agg$aggregation$betahat), tolerance = 1e-12)
  cel <- cohort_sensitivity(fit, M = 0.01, target = "cohort-dynamic", method = "FLCI",
                            full.path = FALSE, balanced = TRUE, seed = 9)
  p <- plot(cel)
  # Two cohorts with five window periods each, for one bound.
  expect_true(all(p$data$relative_period >= -2 & p$data$relative_period <= 2))
  expect_identical(nrow(p$data), 10L)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  att <- cohort_sensitivity(fit, M = c(0, 0.01), method = "FLCI", full.path = FALSE,
                            balanced = TRUE, seed = 9)
  caption <- gsub("\n", " ", plot(att)$labels$caption, fixed = TRUE)
  expect_match(caption, "Balanced window: relative periods -2 to 2; cohort paths: full.", fixed = TRUE)
})
