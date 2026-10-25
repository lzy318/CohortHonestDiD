test_that("one-post-period RM is finite and has the analytic plug-in bounds", {
  fit <- cohort_fit(y ~ d, data = independent_short_panel(),
                    index = c("id", "time"), estimator = "cs_nt")
  ctl <- list(lf.sims = 64L, grid.points = 21L, max.evaluations = 5000L)
  ans <- cohort_sensitivity(fit, restriction = "RM", M = c(0, 0.5), seed = 83, control = ctl)
  expect_identical(ans$settings$method, "c_lf")
  expect_true(all(ans$intervals$success))
  expect_true(all(is.finite(ans$intervals$lower)))
  expect_true(all(is.finite(ans$intervals$upper)))
  i <- fit$coordinate_info
  beta_post <- unname(fit$betahat[i$type == "post"])
  benchmark <- abs(unname(fit$betahat[i$relative_period == -1]))
  expect_equal(ans$intervals$plugin.lower, beta_post - c(0, 0.5) * benchmark, tolerance = 1e-9)
  expect_equal(ans$intervals$plugin.upper, beta_post + c(0, 0.5) * benchmark, tolerance = 1e-9)
  expect_true(all(ans$intervals$lower <= beta_post & ans$intervals$upper >= beta_post))
})

test_that("stochastic RM inference is reproducible across two workers", {
  fit <- cohort_fit(y ~ d, data = independent_short_panel(),
                    index = c("id", "time"), estimator = "cs_nt")
  ctl <- list(lf.sims = 64L, grid.points = 21L, max.evaluations = 5000L)
  set.seed(955)
  rng <- .Random.seed
  a <- cohort_sensitivity(fit, restriction = "RM", M = c(0, 0.5), seed = 91, control = ctl)
  b <- cohort_sensitivity(fit, restriction = "RM", M = c(0, 0.5), seed = 91,
                          parallel = TRUE, cores = 2L, control = ctl)
  expect_true(all(a$intervals$success) && all(b$intervals$success))
  expect_identical(.Random.seed, rng)
  expect_identical(a$settings$task.seeds, b$settings$task.seeds)
  expect_equal(ordered_endpoints(a), ordered_endpoints(b), tolerance = 1e-7)
})

test_that("RM retains its dimensionless bound when outcome units change", {
  d <- independent_short_panel()
  fit <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nt")
  d$y <- 100 * d$y
  scaled <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nt")
  ctl <- list(lf.sims = 64L, grid.points = 21L, max.evaluations = 5000L)
  a <- cohort_sensitivity(fit, restriction = "RM", M = 0.5, seed = 86, control = ctl)
  b <- cohort_sensitivity(scaled, restriction = "RM", M = 0.5, seed = 86, control = ctl)
  expect_true(all(a$intervals$success) && all(b$intervals$success))
  expect_equal(ordered_endpoints(b) / 100, ordered_endpoints(a), tolerance = 1e-7)
  expect_equal(b$intervals$plugin.lower / 100, a$intervals$plugin.lower, tolerance = 1e-10)
  expect_equal(b$intervals$plugin.upper / 100, a$intervals$plugin.upper, tolerance = 1e-10)
})
