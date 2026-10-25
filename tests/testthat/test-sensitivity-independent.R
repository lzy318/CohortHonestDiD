test_that("three-period FLCI equals an independent folded-normal calculation", {
  d <- independent_short_panel()
  fit <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nt")
  info <- fit$coordinate_info
  v <- c(1, -2, 1)[info$time_index]
  center <- sum(v * fit$betahat)
  se <- sqrt(as.numeric(t(v) %*% fit$Sigma %*% v))
  ans <- cohort_sensitivity(fit, M = c(0, 0.01), method = "FLCI", seed = 42)
  expect_s3_class(ans, "cohort_sensitivity")
  expect_identical(ans$settings$method, "flci")
  expect_true(ans$settings$full.path)
  expect_true(all(ans$intervals$success))
  expect_equal(ans$intervals$center, rep(center, 2), tolerance = 1e-9)
  half <- vapply(c(0, 0.01), independent_folded_half_length, numeric(1), se = se)
  expect_equal(ans$intervals$lower, center - half, tolerance = 1e-8)
  expect_equal(ans$intervals$upper, center + half, tolerance = 1e-8)
  expect_equal(ans$intervals$plugin.lower, center - c(0, 0.01), tolerance = 1e-9)
  expect_equal(ans$intervals$plugin.upper, center + c(0, 0.01), tolerance = 1e-9)
})

test_that("balanced imputation and CS-NYT FLCIs agree", {
  d <- independent_panel()
  imp <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "imp")
  cs <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nyt")
  a <- cohort_sensitivity(imp, M = c(0, 0.01), method = "FLCI", seed = 31)
  b <- cohort_sensitivity(cs, M = c(0, 0.01), method = "FLCI", seed = 31)
  expect_true(all(a$intervals$success))
  expect_true(all(b$intervals$success))
  expect_equal(ordered_endpoints(a), ordered_endpoints(b), tolerance = 1e-7)
})

test_that("target estimates use independent treated-cell weighting definitions", {
  d <- independent_panel()
  fit <- cohort_fit(y ~ d, data = d, index = c("id", "time"),
                    estimator = "cs_nt", weights = "w")
  i <- fit$coordinate_info
  post <- which(i$type == "post")
  mass <- c(`4` = sum(d$w[d$time == 1 & d$id <= 20]),
            `6` = sum(d$w[d$time == 1 & d$id > 20 & d$id <= 42]))
  target_value <- function(j) {
    w <- unname(mass[as.character(i$cohort[j])])
    sum(w * fit$betahat[j]) / sum(w)
  }
  att <- cohort_sensitivity(fit, M = 0.01, target = "att", method = "FLCI", seed = 10)
  dyn <- cohort_sensitivity(fit, M = 0.01, target = "dynamic", method = "FLCI", seed = 10)
  cel <- cohort_sensitivity(fit, M = 0.01, target = "cohort-dynamic", method = "FLCI", seed = 10)
  expect_equal(att$estimates$estimate, target_value(post), tolerance = 1e-10)
  for (s in sort(unique(i$relative_period[post]))) {
    j <- post[i$relative_period[post] == s]
    expect_equal(dyn$estimates$estimate[dyn$estimates$relative_period == s],
                 target_value(j), tolerance = 1e-10)
  }
  for (j in post) {
    z <- cel$estimates$cohort == i$cohort[j] &
      cel$estimates$relative_period == i$relative_period[j]
    expect_equal(cel$estimates$estimate[z], unname(fit$betahat[j]), tolerance = 1e-10)
  }
})

test_that("SD is equivariant to a change of outcome units", {
  d <- independent_short_panel()
  fit <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nt")
  d$y <- 100 * d$y
  scaled <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nt")
  a <- cohort_sensitivity(fit, M = c(0, 0.01), method = "FLCI", seed = 28)
  b <- cohort_sensitivity(scaled, M = c(0, 1), method = "FLCI", seed = 28)
  expect_true(all(a$intervals$success) && all(b$intervals$success))
  expect_equal(ordered_endpoints(b) / 100, ordered_endpoints(a), tolerance = 1e-7)
})

test_that("singular covariance at zero bound gives the exact deterministic target", {
  d <- independent_short_panel()
  # All within-group paths have only a unit-specific intercept.
  d$y <- d$id / 20 + 0.04 * d$time - 0.1 * d$d
  fit <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nt")
  ans <- cohort_sensitivity(fit, M = 0, method = "FLCI", seed = 12)
  expect_true(all(ans$intervals$success))
  expect_equal(ans$intervals$lower, -0.1, tolerance = 1e-8)
  expect_equal(ans$intervals$upper, -0.1, tolerance = 1e-8)
})

test_that("inference validates bounds, restrictions, targets, and scope", {
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  expect_error(cohort_sensitivity(fit, M = -0.1), "nonnegative")
  expect_error(cohort_sensitivity(fit, M = NA_real_), "finite")
  expect_error(cohort_sensitivity(fit, M = 0, restriction = "RM", method = "FLCI"), "SD")
  expect_error(cohort_sensitivity(fit, M = 1, restriction = "RM", full.path = FALSE), "SD")
  expect_error(cohort_sensitivity(fit, M = 0, target = "cohort-dynamic", framework = "aggregate"), "cohort")
  expect_error(cohort_sensitivity(fit, M = 0, conf.level = 1), "zero|one")
  expect_error(cohort_sensitivity(fit, M = 0, cores = 0), "cores.*integer")
  expect_error(cohort_sensitivity(fit, M = 0, method = "bogus"), "method")
  expect_error(cohort_sensitivity(fit, M = 0, control = list(unknown.control = 1)), "Unknown|unknown|unsupported")
})

test_that("aggregate and cohort representations preserve the conventional target", {
  fit <- cohort_fit(y ~ d, data = independent_panel(), index = c("id", "time"))
  for (target in c("att", "dynamic")) {
    a <- cohort_sensitivity(fit, M = 0.01, target = target, method = "FLCI", seed = 19)
    b <- cohort_sensitivity(fit, M = 0.01, target = target, framework = "aggregate", method = "FLCI", seed = 19)
    expect_equal(a$estimates$estimate, b$estimates$estimate, tolerance = 1e-10)
    expect_equal(a$estimates$se, b$estimates$se, tolerance = 1e-10)
  }
})

test_that("seeded inference preserves the caller RNG and has stable task seeds", {
  fit <- cohort_fit(y ~ d, data = independent_short_panel(), index = c("id", "time"))
  set.seed(501)
  rng <- .Random.seed
  a <- cohort_sensitivity(fit, M = c(0, 0.01), method = "FLCI", seed = 52)
  expect_identical(.Random.seed, rng)
  b <- cohort_sensitivity(fit, M = c(0, 0.01), method = "FLCI", seed = 52)
  expect_identical(.Random.seed, rng)
  expect_identical(a$settings$task.seeds, b$settings$task.seeds)
  expect_equal(ordered_endpoints(a), ordered_endpoints(b), tolerance = 1e-12)
})
