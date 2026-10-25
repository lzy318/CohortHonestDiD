independent_refinement_fixture <- function() {
  d <- expand.grid(time = 1:4, id = 1:72)
  d <- d[order(d$id, d$time), ]
  d$g <- ifelse(d$id <= 24, 0, ifelse(d$id <= 48, 3, 4))
  d$D <- as.integer(d$g > 0 & d$time >= d$g)
  set.seed(271828)
  d$Y <- .01*d$id + .03*d$time + .007*(d$g > 0)*(d$time-1) +
    rnorm(nrow(d), 0, .04) + d$D*(.08 + .01*(d$time-d$g))
  d$D[(d$id == 25 & d$time == 1) | (d$id == 26 & d$time == 2) |
        (d$id == 50 & d$time == 4)] <- NA
  d$D[d$id == 50 & d$time == 3] <- 1
  expect_warning(fit <- cohort_fit(Y ~ D, d, index = c('id', 'time'),
                                  estimator = 'cs_nyt'), 'unbalanced|approximation')
  fit
}

test_that('bounded refinement preserves a genuinely disconnected confidence set', {
  f <- independent_refinement_fixture()
  run <- function(refinements, budget = 100000L, grid = 21L) {
    cohort_sensitivity(f, M = 0, method = 'Conditional', full.path = TRUE,
      cores = 1L, seed = 20260920L,
      control = list(lf.sims = 100L, grid.points = grid,
        max.refinements = refinements, max.evaluations = budget,
        max.expansions = 6L, topology.audit = TRUE))
  }
  coarse_warnings <- character()
  coarse <- withCallingHandlers(run(0L), warning = function(w) {
    coarse_warnings <<- c(coarse_warnings, conditionMessage(w))
    invokeRestart('muffleWarning')
  })
  expect_identical(coarse$results[[1]]$diagnostics$refinements, 0L)
  if (all(coarse$intervals$success)) {
    expect_equal(nrow(coarse$results[[1]]$components), 2L)
    expect_false(any(coarse$results[[1]]$components$lower < .099 &
                     coarse$results[[1]]$components$upper > .099))
  } else {
    expect_true(length(coarse_warnings) > 0L)
    expect_match(coarse$results[[1]]$status, 'topology|endpoint')
    expect_true(all(is.na(coarse$intervals$lower)))
    expect_true(all(is.na(coarse$intervals$upper)))
  }
  refined <- run(3L)
  fine <- run(0L, grid = 101L)
  expect_true(all(fine$intervals$success))
  expect_equal(nrow(fine$results[[1]]$components), 2L)
  expect_true(all(refined$intervals$success))
  expect_equal(nrow(refined$intervals), 2L)
  expect_lt(refined$intervals$upper[1], refined$intervals$lower[2])
  expect_false(any(refined$intervals$lower < .099 &
                   refined$intervals$upper > .099))
  tol <- 8 * max(refined$results[[1]]$diagnostics$endpoint.tolerance,
                 fine$results[[1]]$diagnostics$endpoint.tolerance) +
    1e-7 * fine$estimates$se[1]
  expect_lte(max(abs(refined$intervals$lower - fine$intervals$lower)), tol)
  expect_lte(max(abs(refined$intervals$upper - fine$intervals$upper)), tol)
  if (all(coarse$intervals$success)) {
    coarse_tol <- 8 * max(coarse$results[[1]]$diagnostics$endpoint.tolerance,
                          fine$results[[1]]$diagnostics$endpoint.tolerance) +
      1e-7 * fine$estimates$se[1]
    expect_lte(max(abs(as.matrix(coarse$results[[1]]$components) -
                        as.matrix(fine$results[[1]]$components))), coarse_tol)
  }
  dg <- refined$results[[1]]$diagnostics
  expect_equal(sum(dg$search.attempts$new.evaluations), dg$test.evaluations)
  required <- as.integer(dg$test.evaluations)
  expect_warning(short <- run(3L, required - 1L), 'failed|limit|budget|evaluation')
  expect_true(all(!short$intervals$success))
  expect_true(all(is.na(short$intervals$lower)))
  expect_true(all(is.na(short$intervals$upper)))
  for (budget in c(required, required + 1L)) {
    exact <- run(3L, budget)
    expect_true(all(exact$intervals$success))
    expect_equal(exact$intervals, refined$intervals)
    expect_lte(exact$results[[1]]$diagnostics$test.evaluations, budget)
  }
})

test_that('refinement count validates numeric scalar integer inputs', {
  f <- independent_refinement_fixture()
  for (bad in list(-1, 1.5, NA_real_, Inf, '3', TRUE, c(1, 2), 2147483648)) {
    expect_error(cohort_sensitivity(f, M = .015, cores = 1L,
      control = list(max.refinements = bad)), 'refinements|integer|numeric|finite|scalar')
  }
})
