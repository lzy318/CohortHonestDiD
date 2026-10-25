independent_boundary_scalar <- function() {
  nm <- c("pre", "reference", "post")
  cohort_input(
    betahat = setNames(c(.1, 0, .4), nm),
    vcov = matrix(c(.04, 0, .01, 0, 0, 0, .01, 0, .09), 3,
                  dimnames = list(nm, nm)),
    coefficient_info = data.frame(
      coordinate = nm, cohort = 3L,
      relative_period = c(-1L, 0L, 1L), type = c("pre", "pre", "post")
    ),
    W = structure(diag(3), dimnames = list(nm, nm)),
    normalization = matrix(c(0, 1, 0), 1,
                           dimnames = list("reference", nm)),
    cohort_sizes = c("3" = 24), estimator = "cs_nt"
  )
}

test_that("oversized worker counts are rejected as worker inputs", {
  f <- independent_boundary_scalar()
  expect_error(cohort_sensitivity(f, M = .015, cores = 2147483648),
               "cores|worker", ignore.case = TRUE)
})

test_that("logical values cannot silently become numeric controls", {
  f <- independent_boundary_scalar()
  numeric_controls <- c("lf.sims", "endpoint.tol", "grid.radius",
                        "max.expansions", "max.evaluations",
                        "max.components", "component.cores")
  for (key in numeric_controls) {
    expect_error(
      cohort_sensitivity(f, M = .015, cores = 1,
                         control = setNames(list(TRUE), key)),
      paste(key, "control|numeric|integer", sep = "|"),
      info = key
    )
  }
})

test_that("FLCI controls validate even when the interval is analytic", {
  f <- independent_boundary_scalar()
  malformed <- list(
    list(unknown = 1), structure(list(1), names = ""),
    list(maxit = 1, maxit = 2)
  )
  for (control in malformed) {
    expect_error(cohort_sensitivity(
      f, M = .015, method = "FLCI", cores = 1, control = list(flci.control = control)
    ), "flci|control|named|duplicate|unknown", ignore.case = TRUE)
  }
  for (key in c("maxit", "n_random_starts")) {
    invalid <- list(-1, 1.5, NA_real_, Inf, c(1, 2), "1", TRUE)
    for (value in invalid) {
      expect_error(cohort_sensitivity(
        f, M = .015, method = "FLCI", cores = 1,
        control = list(flci.control = setNames(list(value), key))
      ), paste(key, "control|numeric|integer|finite|scalar", sep = "|"),
      info = paste(key, paste(value, collapse = ",")))
    }
  }
})

test_that("zero optional optimizer iterations retain certified analytic FLCI", {
  f <- independent_boundary_scalar()
  bound <- .015
  sigma <- sqrt(.15)
  h <- uniroot(function(x) {
    pnorm((x - bound) / sigma) + pnorm((x + bound) / sigma) - 1.95
  }, c(0, bound + 10 * sigma), tol = 1e-13)$root
  z <- cohort_sensitivity(f, M = bound, method = "FLCI", cores = 1, seed = 20260920L,
    control = list(flci.control = list(maxit = 0, n_random_starts = 0)))
  expect_true(all(z$intervals$success))
  expect_true(z$results[[1]]$diagnostics$endpoint.certified)
  expect_equal(z$intervals$lower, .5 - h, tolerance = 1e-7)
  expect_equal(z$intervals$upper, .5 + h, tolerance = 1e-7)
})

test_that("ordinary and extreme common weight scales give the same fit", {
  z <- expand.grid(time = 1:4, id = 1:72)
  z <- z[order(z$id, z$time), ]
  z$g <- ifelse(z$id <= 24, 0, ifelse(z$id <= 48, 3, 4))
  z$D <- as.integer(z$g > 0 & z$time >= z$g)
  set.seed(271828)
  z$Y <- .01 * z$id + .03 * z$time +
    .007 * (z$g > 0) * (z$time - 1) + rnorm(nrow(z), sd = .04) +
    z$D * (.08 + .01 * (z$time - z$g))
  for (estimator in c("imp", "cs_nyt", "cs_nt")) {
    ordinary <- cohort_fit(Y ~ D, z, index = c("id", "time"),
                           estimator = estimator)
    scaled <- cohort_fit(Y ~ D, z, index = c("id", "time"),
                         estimator = estimator,
                         weights = rep(1e308, nrow(z)))
    expect_equal(scaled$betahat, ordinary$betahat, tolerance = 1e-10,
                 info = estimator)
    expect_equal(scaled$Sigma, ordinary$Sigma, tolerance = 1e-10,
                 info = estimator)
    expect_equal(scaled$W, ordinary$W, tolerance = 1e-10,
                 info = estimator)
  }
})

test_that("one pre period is diagnosed as insufficient sensitivity history", {
  nm <- c("reference", "post1", "post2")
  f <- cohort_input(
    setNames(c(0, .1, .2), nm),
    structure(diag(c(0, .04, .09)), dimnames = list(nm, nm)),
    data.frame(coordinate = nm, cohort = 2L,
               relative_period = 0:2, type = c("pre", "post", "post")),
    structure(diag(3), dimnames = list(nm, nm)),
    matrix(c(1, 0, 0), 1, dimnames = list("reference", nm)),
    c("2" = 24), estimator = "cs_nt"
  )
  for (restriction in c("SD", "RM")) {
    result <- tryCatch(
      suppressWarnings(cohort_sensitivity(
        f, restriction = restriction, M = .015, cores = 1, seed = 1
      )), error = identity
    )
    if (inherits(result, "error")) {
      expect_match(conditionMessage(result), "pre|history|period|benchmark",
                   ignore.case = TRUE)
    } else {
      expect_true(all(!result$intervals$success))
      expect_true(all(is.na(result$intervals$lower)))
      expect_true(all(is.na(result$intervals$upper)))
      expect_true(all(grepl("pre|history|period|benchmark",
                            result$intervals$message, ignore.case = TRUE)))
    }
  }
})
