independent_plot_result <- function(disconnected = FALSE) {
  intervals <- data.frame(target = "att", label = "Average treatment effect",
    cohort = NA_integer_, relative_period = NA_integer_, M = 1,
    component = 1L, lower = -0.4, upper = 0.6, center = 0.1,
    plugin.lower = -0.3, plugin.upper = -0.1, success = TRUE,
    status = "ok", message = "")
  if (disconnected) {
    intervals <- intervals[c(1, 1), ]
    intervals$component <- 1:2
    intervals$lower <- c(-0.4, 0.2)
    intervals$upper <- c(-0.1, 0.6)
    intervals$status <- "disconnected"
  }
  structure(list(intervals = intervals,
    estimates = data.frame(estimate = -0.2),
    fit = list(data_info = list(balanced = TRUE)),
    settings = list(target = "att", restriction = if (disconnected) "RM" else "SD",
      framework = "cohort", method = if (disconnected) "c_lf" else "flci",
      full.path = TRUE, conf.level = 0.95, M = 1)), class = "cohort_sensitivity")
}

test_that("disconnected confidence sets retain a visible gap and no hull-center marker", {
  x <- independent_plot_result(disconnected = TRUE)
  p <- plot(x, identified.set = FALSE)
  expect_s3_class(p, "ggplot")
  idx <- which(vapply(p$layers, function(z) inherits(z$geom, "GeomLinerange"), logical(1)))
  expect_length(idx, 1L)
  b <- ggplot2::ggplot_build(p)$data[[idx]]
  expect_equal(sort(b$ymin), c(-0.4, 0.2), tolerance = 1e-12)
  expect_equal(sort(b$ymax), c(-0.1, 0.6), tolerance = 1e-12)
  expect_false(any(vapply(p$layers, function(z) inherits(z$geom, "GeomPoint"), logical(1))))
})

test_that("FLCI plots use the optimized center rather than the plug-in midpoint", {
  x <- independent_plot_result()
  p <- plot(x)
  expect_s3_class(p, "ggplot")
  idx <- which(vapply(p$layers, function(z) inherits(z$geom, "GeomPoint"), logical(1)))
  expect_length(idx, 1L)
  expect_equal(ggplot2::ggplot_build(p)$data[[idx]]$y, 0.1, tolerance = 1e-12)
})

test_that("plot dispatch supports fitted and dynamic result objects", {
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  expect_s3_class(plot(fit), "ggplot")
  expect_s3_class(plot(fit, type = "cohort-dynamic"), "ggplot")
  ans <- cohort_sensitivity(fit, M = c(0, 0.01), target = "dynamic", method = "FLCI", seed = 15)
  expect_s3_class(plot(ans), "ggplot")
  expect_s3_class(ggplot2::ggplot_build(plot(ans)), "ggplot_built")
})

test_that("all-failed and all-empty dynamic results retain conventional estimates", {
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  ans <- cohort_sensitivity(fit, M = 0.01, target = "dynamic", method = "FLCI", seed = 15)
  for (is_failure in c(TRUE, FALSE)) {
    x <- ans
    x$intervals$lower <- x$intervals$upper <- NA_real_
    x$intervals$success <- !is_failure
    x$intervals$status <- if (is_failure) "inference_error" else "empty"
    p <- plot(x)
    built <- ggplot2::ggplot_build(p)
    expect_s3_class(built, "ggplot_built")
    points <- which(vapply(p$layers, function(z) inherits(z$geom, "GeomPoint"), logical(1)))
    expect_length(points, 1L)
    expect_equal(sort(built$data[[points]]$y), sort(as.numeric(fit$betahat)), tolerance = 1e-12)
    expect_match(p$labels$caption, if (is_failure) "failed or unresolved" else "empty set")
    expect_match(p$labels$caption, "Black: conventional\\s+estimates")
  }
})

test_that("all-failed ATT results retain the requested sensitivity-bound domain", {
  x <- independent_plot_result()
  x$intervals <- x$intervals[c(1, 1), ]
  x$intervals$M <- c(2, 4)
  x$settings$M <- c(2, 4)
  x$intervals$success <- FALSE
  x$intervals$lower <- x$intervals$upper <- NA_real_
  x$intervals$status <- "inference_error"
  p <- plot(x)
  built <- ggplot2::ggplot_build(p)
  rg <- built$layout$panel_params[[1]]$x.range
  expect_lte(rg[1], 2)
  expect_gte(rg[2], 4)
  expect_match(p$labels$caption, "failed or unresolved")
})
