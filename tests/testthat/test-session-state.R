test_that("fitting and seeded inference do not write into the working directory", {
  place <- withr::local_tempdir()
  withr::local_dir(place)
  before_options <- options()
  before <- list.files(place, all.files = TRUE, no.. = TRUE, recursive = TRUE)
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  ans <- cohort_sensitivity(fit, M = 0.01, method = "FLCI", seed = 41)
  expect_true(all(ans$intervals$success))
  expect_identical(normalizePath(getwd()), normalizePath(place))
  expect_identical(list.files(place, all.files = TRUE, no.. = TRUE, recursive = TRUE), before)
  expect_identical(options()[names(before_options)], before_options)
})

test_that("two-worker inference reproduces sequential seeded results", {
  # PSOCK workers use the installed namespace, as they do in normal package use.
  skip_if_not(nzchar(system.file(package = "CohortHonestDiD")))
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  set.seed(702)
  rng <- .Random.seed
  serial <- cohort_sensitivity(fit, M = c(0, 0.01), method = "FLCI", seed = 97)
  parallel <- cohort_sensitivity(fit, M = c(0, 0.01), method = "FLCI", parallel = TRUE, cores = 2, seed = 97)
  expect_true(all(parallel$intervals$success))
  expect_equal(ordered_endpoints(parallel), ordered_endpoints(serial), tolerance = 1e-7)
  expect_identical(.Random.seed, rng)
  expect_lte(parallel$settings$cores, 2L)
  single <- cohort_sensitivity(fit, M = 0, method = "FLCI", parallel = TRUE, cores = 2, seed = 97)
  expect_equal(single$settings$cores, 1L)
})

test_that("one RM task uses at most two component workers only when requested", {
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  ctl <- list(lf.sims = 64L, grid.points = 21L, max.evaluations = 5000L, component.cores = 2L)
  serial <- cohort_sensitivity(fit, restriction = "RM", M = 0.5,
    parallel = FALSE, cores = 2L, seed = 93, control = ctl)
  workers <- cohort_sensitivity(fit, restriction = "RM", M = 0.5,
    parallel = TRUE, cores = 2L, seed = 93, control = ctl)
  expect_true(all(serial$intervals$success) && all(workers$intervals$success))
  expect_equal(serial$settings$cores, 1L)
  expect_equal(serial$results[[1]]$diagnostics$component.cores, 1L)
  expect_equal(workers$settings$cores, 2L)
  expect_equal(workers$results[[1]]$diagnostics$component.cores, 2L)
  expect_equal(ordered_endpoints(workers), ordered_endpoints(serial), tolerance = 1e-7)
})

test_that("fixed seeds do not depend on the caller RNG kind or worker count", {
  withr::local_seed(110, .rng_kind = "L'Ecuyer-CMRG")
  kind <- RNGkind()
  rng <- .Random.seed
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  ctl <- list(lf.sims = 64L, grid.points = 21L, max.evaluations = 5000L)
  serial <- cohort_sensitivity(fit, restriction = "RM", M = c(0, 0.5),
    seed = 97, control = ctl)
  workers <- cohort_sensitivity(fit, restriction = "RM", M = c(0, 0.5),
    parallel = TRUE, cores = 2L, seed = 97, control = ctl)
  expect_true(all(serial$intervals$success) && all(workers$intervals$success))
  expect_identical(RNGkind(), kind)
  expect_identical(.Random.seed, rng)
  expect_equal(ordered_endpoints(workers), ordered_endpoints(serial), tolerance = 1e-7)
})
