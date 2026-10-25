test_that("external estimates preserve named matrices and reject indefinite covariance", {
  source_fit <- cohort_fit(y ~ d, data = independent_short_panel(),
                           index = c("id", "time"), estimator = "cs_nt")
  make_input <- function(V) cohort_input(
    betahat = source_fit$betahat, vcov = V,
    coefficient_info = source_fit$coordinate_info, W = source_fit$W,
    normalization = source_fit$normalization,
    cohort_sizes = source_fit$cohort_sizes, estimator = "cs_nt")
  fit <- make_input(source_fit$Sigma)
  expect_equal(fit$betahat, source_fit$betahat, tolerance = 1e-12)
  expect_equal(fit$Sigma, source_fit$Sigma, tolerance = 1e-12)
  j <- which(source_fit$coordinate_info$type == "post")[1]
  bad <- source_fit$Sigma
  bad[j, j] <- -1
  expect_error(make_input(bad), "positive semidefinite")
  expect_error(make_input(-source_fit$Sigma * 1e-8), "positive semidefinite")
  expect_error(cohort_input(
    betahat = c(source_fit$betahat, source_fit$betahat[1]), vcov = source_fit$Sigma,
    coefficient_info = source_fit$coordinate_info, W = source_fit$W,
    normalization = source_fit$normalization, cohort_sizes = source_fit$cohort_sizes,
    estimator = "cs_nt"), "unique|duplicat|exactly|length")
  bad_W <- source_fit$W
  reference <- which(source_fit$coordinate_info$relative_period == 0)[1]
  bad_W[reference, j] <- 1
  expect_error(cohort_input(
    betahat = source_fit$betahat, vcov = source_fit$Sigma,
    coefficient_info = source_fit$coordinate_info, W = bad_W,
    normalization = source_fit$normalization, cohort_sizes = source_fit$cohort_sizes,
    estimator = "cs_nt"), "normalization")
  a <- cohort_sensitivity(source_fit, M = 0.01, method = "FLCI", seed = 12)
  b <- cohort_sensitivity(fit, M = 0.01, method = "FLCI", seed = 12)
  expect_equal(ordered_endpoints(a), ordered_endpoints(b), tolerance = 1e-9)
})
