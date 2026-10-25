test_that("the SD default is C-LF under the post-only restriction", {
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  ctl <- list(lf.sims = 200L, grid.points = 21L)
  default <- cohort_sensitivity(fit, M = c(0, 0.01), seed = 64, control = ctl)
  explicit <- cohort_sensitivity(fit, M = c(0, 0.01), method = "C-LF",
                                 full.path = FALSE, seed = 64, control = ctl)
  expect_identical(default$settings$method, "c_lf")
  expect_false(default$settings$full.path)
  expect_true(all(default$intervals$success))
  expect_identical(default$intervals, explicit$intervals)
  expect_identical(default$estimates, explicit$estimates)
  expect_identical(default$results, explicit$results)
  expect_identical(default$settings, explicit$settings)
  given_null <- cohort_sensitivity(fit, M = c(0, 0.01), method = NULL,
                                   full.path = NULL, seed = 64, control = ctl)
  expect_identical(given_null$intervals, default$intervals)
  expect_identical(given_null$results, default$results)
  expect_identical(given_null$settings, default$settings)
})

test_that("explicit methods keep their own SD scope defaults", {
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  ctl <- list(lf.sims = 200L, grid.points = 21L)
  flci <- cohort_sensitivity(fit, M = 0.01, method = "FLCI", seed = 65)
  conditional <- cohort_sensitivity(fit, M = 0.01, method = "Conditional",
                                    seed = 65, control = ctl)
  clf <- cohort_sensitivity(fit, M = 0.01, method = "C-LF", seed = 65, control = ctl)
  expect_identical(flci$settings$method, "flci")
  expect_true(flci$settings$full.path)
  expect_identical(conditional$settings$method, "conditional")
  expect_false(conditional$settings$full.path)
  expect_identical(clf$settings$method, "c_lf")
  expect_false(clf$settings$full.path)
})

test_that("the RM default is C-LF", {
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  ctl <- list(lf.sims = 64L, grid.points = 21L, max.evaluations = 5000L)
  default <- cohort_sensitivity(fit, restriction = "RM", M = 0.5, seed = 66, control = ctl)
  explicit <- cohort_sensitivity(fit, restriction = "RM", M = 0.5, method = "C-LF",
                                 seed = 66, control = ctl)
  expect_identical(default$settings$method, "c_lf")
  expect_null(default$settings$full.path)
  expect_true(all(default$intervals$success))
  expect_identical(default$intervals, explicit$intervals)
  expect_identical(default$results, explicit$results)
  expect_identical(default$settings, explicit$settings)
})

test_that("a default SD result prints its procedure and scope", {
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  ans <- cohort_sensitivity(fit, M = 0.01, seed = 67,
                            control = list(lf.sims = 200L, grid.points = 21L))
  out <- capture.output(print(ans))
  expect_match(out[1], "| C-LF | post-only", fixed = TRUE)
})
