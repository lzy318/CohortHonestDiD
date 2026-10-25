independent_d4_scalar_fixture <- function() {
  nm <- c('early', 'reference', 'post')
  beta <- setNames(c(.1, 0, .4), nm)
  Sigma <- matrix(0, 3, 3, dimnames = list(nm, nm))
  Sigma[c(1, 3), c(1, 3)] <- matrix(c(.04, .01, .01, .09), 2)
  W <- diag(3)
  dimnames(W) <- list(nm, nm)
  N <- matrix(c(0, 1, 0), 1, dimnames = list('reference', nm))
  info <- data.frame(coordinate = nm, cohort = 3L,
    relative_period = c(-1L, 0L, 1L), type = c('pre', 'pre', 'post'))
  cohort_input(beta, Sigma, info, W, N, c('3' = 1), 'cs_nt')
}

test_that('Conditional cannot claim stability before two wider-domain agreements', {
  fit <- independent_d4_scalar_fixture()
  for (cap in c(0L, 1L)) {
    expect_warning(answer <- cohort_sensitivity(fit, M = 0,
      method = 'Conditional', full.path = FALSE, cores = 1L,
      seed = 20260920L, control = list(grid.points = 21L,
        max.expansions = cap, max.evaluations = 100000L)),
      'failed|domain|range|expansion|stability')
    expect_true(all(!answer$intervals$success))
    expect_true(all(is.na(answer$intervals$lower)))
    expect_true(all(is.na(answer$intervals$upper)))
    expect_match(paste(answer$intervals$message, collapse = ' '),
      'domain|range|expansion|stability', ignore.case = TRUE)
  }
  answer <- cohort_sensitivity(fit, M = 0, method = 'Conditional',
    full.path = FALSE, cores = 1L, seed = 20260920L,
    control = list(grid.points = 21L, max.expansions = 2L,
      max.evaluations = 100000L))
  expect_true(all(answer$intervals$success))
  expected <- .5 + c(-1, 1) * qnorm(.975) * sqrt(.15)
  tolerance <- 8 * answer$results[[1]]$diagnostics$endpoint.tolerance +
    1e-7 * answer$estimates$se[1]
  expect_lte(max(abs(c(answer$intervals$lower, answer$intervals$upper) -
    expected)), tolerance)
  expect_identical(answer$results[[1]]$diagnostics$support.certified, FALSE)
  expect_identical(answer$results[[1]]$diagnostics$topology.certified, FALSE)
  expect_identical(answer$results[[1]]$diagnostics$inversion.exact, FALSE)
})

independent_d4_balanced_fixture <- function(estimator) {
  d <- expand.grid(time = 1:4, id = 1:72)
  d <- d[order(d$id, d$time), ]
  d$g <- ifelse(d$id <= 24, 0, ifelse(d$id <= 48, 3, 4))
  d$D <- as.integer(d$g > 0 & d$time >= d$g)
  set.seed(271828)
  d$Y <- .01*d$id + .03*d$time + .007*(d$g > 0)*(d$time - 1) +
    rnorm(nrow(d), 0, .04) + d$D*(.08 + .01*(d$time - d$g))
  cohort_fit(Y ~ D, d, index = c('id', 'time'), estimator = estimator)
}

test_that('wider Conditional searches retain the known separated component', {
  # Full-precision complete Imp union saved before D4, independent of new code.
  expected <- matrix(c(-.00951900336365546, .139143511615014,
    .121620125095407, .311141338412962), nrow = 2L)
  answers <- lapply(c('imp', 'cs_nyt'), function(estimator) {
    cohort_sensitivity(independent_d4_balanced_fixture(estimator), M = 0,
      target = 'dynamic', method = 'Conditional', framework = 'cohort',
      full.path = TRUE, cores = 1L, seed = 20260920L,
      control = list(lf.sims = 100L, grid.points = 21L,
        max.expansions = 6L, max.refinements = 3L,
        max.evaluations = 100000L, topology.audit = TRUE))
  })
  for (answer in answers) {
    result <- answer$results[[2]]
    expect_true(result$success)
    expect_equal(nrow(result$components), 2L)
    tolerance <- 8 * max(result$diagnostics$endpoint.tolerance,
      1.37887839366637e-7) + 1e-7 * .0137887839366637
    expect_lte(max(abs(as.matrix(result$components) - expected)), tolerance)
    expect_true(any(result$components$lower < .2 & result$components$upper > .2))
    expect_false(any(result$components$lower < .13 & result$components$upper > .13))
    expect_identical(result$diagnostics$support.certified, FALSE)
    expect_identical(result$diagnostics$topology.certified, FALSE)
    expect_identical(result$diagnostics$inversion.exact, FALSE)
  }
})

test_that('a wider initial range retains reference spacing and the known gap', {
  expected <- matrix(c(-.00951900336365546, .139143511615014,
    .121620125095407, .311141338412962), nrow = 2L)
  for (estimator in c('imp', 'cs_nyt')) {
    answer <- cohort_sensitivity(independent_d4_balanced_fixture(estimator),
      M = 0, target = 'dynamic', method = 'Conditional', full.path = TRUE,
      cores = 1L, seed = 20260920L,
      control = list(lf.sims = 100L, grid.points = 21L, grid.radius = 32,
        max.expansions = 6L, max.refinements = 3L,
        max.evaluations = 100000L, topology.audit = TRUE))
    result <- answer$results[[2L]]
    d <- result$diagnostics
    expect_true(result$success)
    expect_equal(nrow(result$components), 2L)
    tolerance <- 8 * max(d$endpoint.tolerance, 1.37887839366637e-7) +
      1e-7 * .0137887839366637
    expect_lte(max(abs(as.matrix(result$components) - expected)), tolerance)
    expect_false(any(result$components$lower < .13 &
      result$components$upper > .13))
    expect_equal(answer$settings$control$grid.points, 21L)
    expect_true(d$grid.density.adjusted)
    expect_gte(d$initial.grid.points, 81)
    expect_equal(d$reference.domain.radius, 8 * d$search.scale, tolerance = 0)
    expect_lte(diff(d$initial.domain) / (d$initial.grid.points - 1),
      d$reference.grid.spacing)
    h <- d$domain.history
    expect_equal(h$regular.points,
      (d$initial.grid.points - 1) * 2^h$refinement * 2^h$domain + 1,
      tolerance = 0)
  }
})

test_that('density uses the plugin radius and respects the original budget', {
  f <- independent_d4_scalar_fixture()
  control <- list(grid.points = 21L, grid.radius = 13.25,
    max.expansions = 6L, max.refinements = 0L,
    max.evaluations = 100000L, topology.audit = TRUE)
  answer <- cohort_sensitivity(f, M = .2, method = 'Conditional',
    full.path = FALSE, cores = 1L, seed = 20260920L, control = control)
  result <- answer$results[[1L]]
  d <- result$diagnostics
  expect_true(result$success)
  expect_identical(d$refinements, 0L)
  reference_radius <- .2 + 8 * sqrt(.15)
  requested_radius <- .2 + 13.25 * sqrt(.15)
  expect_lte(abs(d$reference.domain.radius - reference_radius),
    1e-10 + 1e-8 * reference_radius)
  expect_lte(max(abs(d$initial.domain -
    (.5 + c(-requested_radius, requested_radius)))),
    1e-10 + 1e-8 * requested_radius)
  intervals <- d$initial.grid.points - 1
  minimum_ratio <- 20 * requested_radius / reference_radius
  expect_gte(intervals, minimum_ratio)
  expect_true(intervals - 1 < minimum_ratio ||
    diff(d$initial.domain) / (intervals - 1) > d$reference.grid.spacing)
  expect_lte(diff(d$initial.domain) / intervals, d$reference.grid.spacing)

  control$max.evaluations <- as.integer(d$initial.grid.points - 1)
  expect_warning(limited <- cohort_sensitivity(f, M = .2,
    method = 'Conditional', full.path = FALSE, cores = 1L,
    seed = 20260920L, control = control), 'failed|budget|evaluat|limit')
  expect_false(limited$results[[1L]]$success)
  expect_identical(limited$results[[1L]]$status, 'evaluation_limit')
  expect_equal(limited$results[[1L]]$diagnostics$test.evaluations, 0)
  expect_true(all(is.na(limited$intervals$lower)))
  expect_true(all(is.na(limited$intervals$upper)))

  control$grid.radius <- 1e200
  control$max.evaluations <- 100000L
  expect_warning(huge <- cohort_sensitivity(f, M = .2,
    method = 'Conditional', full.path = FALSE, cores = 1L,
    seed = 20260920L, control = control), 'failed|grid|budget|limit')
  expect_false(huge$results[[1L]]$success)
  expect_equal(huge$results[[1L]]$diagnostics$test.evaluations, 0)
  expect_true(all(is.na(huge$intervals$lower)))
  expect_true(all(is.na(huge$intervals$upper)))
})
