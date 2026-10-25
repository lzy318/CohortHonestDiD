test_that('public bias maps preserve extreme common finite mass scaling', {
  for (est in c('imp', 'cs_nyt', 'cs_nt')) {
    regular <- build_W(c(3, 4), c(1, 1), 1, 4, estimator = est)
    scaled <- build_W(c(3, 4), c(1e308, 1e308), 1e308, 4,
                      estimator = est)
    expect_equal(scaled, regular, tolerance = 1e-10)
  }
})

independent_scalar_input <- function(zero = FALSE) {
  nm <- c('early', 'reference', 'post')
  b <- setNames(c(.1, 0, .4), nm)
  V <- matrix(0, 3, 3, dimnames = list(nm, nm))
  if (!zero) V[c(1, 3), c(1, 3)] <- matrix(c(.04, .01, .01, .09), 2)
  W <- diag(3)
  dimnames(W) <- list(nm, nm)
  N <- matrix(c(0, 1, 0), 1, dimnames = list('reference', nm))
  info <- data.frame(coordinate = nm, cohort = 3L,
                     relative_period = c(-1L, 0L, 1L),
                     type = c('pre', 'pre', 'post'))
  cohort_input(b, V, info, W, N, c('3' = 1), estimator = 'cs_nt')
}

test_that('single-post FLCI matches an independent folded-normal calculation', {
  input <- independent_scalar_input()
  for (scope in c(FALSE, TRUE)) {
    result <- cohort_sensitivity(input, M = c(0, .015, .2),
                                 method = 'FLCI', full.path = scope,
                                 seed = 20260920L, cores = 1L)
    for (i in seq_len(nrow(result$intervals))) {
      bound <- result$intervals$M[i]
      sigma <- sqrt(.15)
      half <- uniroot(function(h) {
        pnorm((h - bound) / sigma) + pnorm((h + bound) / sigma) - 1.95
      }, c(0, bound + 10 * sigma), tol = 1e-13)$root
      expect_true(result$intervals$success[i])
      expect_equal(result$intervals$lower[i], .5 - half, tolerance = 1e-7)
      expect_equal(result$intervals$upper[i], .5 + half, tolerance = 1e-7)
    }
  }
})

test_that('zero-noise FLCI has its exact deterministic limit', {
  input <- independent_scalar_input(zero = TRUE)
  for (scope in c(FALSE, TRUE)) {
    result <- cohort_sensitivity(input, M = c(0, .015, .2),
                                 method = 'FLCI', full.path = scope,
                                 seed = 20260920L, cores = 1L)
    expect_true(all(result$intervals$success))
    expect_equal(result$intervals$lower, .5 - c(0, .015, .2), tolerance = 1e-9)
    expect_equal(result$intervals$upper, .5 + c(0, .015, .2), tolerance = 1e-9)
  }
})

test_that('hand-derived two-cohort mapping has the specified coefficients', {
  for (est in c('imp', 'cs_nyt', 'cs_nt')) {
    observed <- build_W(c(3, 4), c(2, 3), 5, 4, estimator = est)
    expected <- diag(8)
    if (est != 'cs_nt') expected[7, 8] <- 3 / 8
    if (est == 'cs_nyt') expected[7, 4] <- -3 / 8
    expect_equal(unname(observed), expected, tolerance = 1e-10)
    expect_equal(det(observed), 1, tolerance = 1e-10)
  }
})

test_that('finite large PSD covariance retains its scalar FLCI limit', {
  nm <- c('early', 'later', 'reference', 'post')
  covariance <- matrix(0, 4, 4, dimnames = list(nm, nm))
  covariance[1:2, 1:2] <- 2e307 * matrix(c(4, 2, 2, 1), 2)
  covariance[4, 4] <- 2e307
  W <- diag(4)
  dimnames(W) <- list(nm, nm)
  f <- cohort_input(setNames(c(.1, .2, 0, .4), nm), covariance,
    data.frame(coordinate = nm, cohort = 4L,
      relative_period = c(-2L, -1L, 0L, 1L),
      type = c('pre', 'pre', 'pre', 'post')),
    W, matrix(c(0, 0, 1, 0), 1, dimnames = list('reference', nm)),
    c('4' = 1), estimator = 'cs_nt')
  z <- cohort_sensitivity(f, M = c(0, .2), method = 'FLCI', full.path = TRUE,
                          seed = 20260920L, cores = 1L)
  h <- qnorm(.975) * sqrt(4e307)
  expect_true(all(z$intervals$success))
  expect_equal(z$intervals$center, rep(.6, 2), tolerance = 1e-9)
  expect_equal(z$intervals$lower / h, rep((.6 - h) / h, 2), tolerance = 1e-8)
  expect_equal(z$intervals$upper / h, rep((.6 + h) / h, 2), tolerance = 1e-8)
  for (rule in z$results) {
    expect_equal(as.numeric(rule$weights$gamma), 0, tolerance = 1e-10)
    expect_equal(unname(rule$weights$coefficients[c(1, 2, 4)]), c(0, 1, 1), tolerance = 1e-10)
  }
})

test_that('positive pre-curvature noise survives a huge post variance at M zero', {
  nm <- c('early', 'later', 'reference', 'post')
  covariance <- diag(c(1e-308, 4e-308, 0, 1e308))
  dimnames(covariance) <- list(nm, nm)
  W <- diag(4)
  dimnames(W) <- list(nm, nm)
  f <- cohort_input(setNames(c(.1, .2, 0, .4), nm), covariance,
    data.frame(coordinate = nm, cohort = 4L,
      relative_period = c(-2L, -1L, 0L, 1L),
      type = c('pre', 'pre', 'pre', 'post')),
    W, matrix(c(0, 0, 1, 0), 1, dimnames = list('reference', nm)),
    c('4' = 1), estimator = 'cs_nt')
  z <- cohort_sensitivity(f, M = 0, method = 'FLCI', full.path = TRUE,
                          seed = 20260920L, cores = 1L)
  rule <- z$results[[1]]
  # In the public additive rule c + gamma*k, gamma = -c'Sigma k / k'Sigma k.
  expect_true(rule$success)
  expect_equal(as.numeric(rule$weights$gamma), 8 / 17, tolerance = 1e-10)
  expect_equal(unname(rule$weights$coefficients[c(1, 2, 4)]),
               c(8 / 17, 1 / 17, 1), tolerance = 1e-10)
  expect_equal(rule$center, .4 + .1 * 8 / 17 + .2 / 17, tolerance = 1e-9)
})

