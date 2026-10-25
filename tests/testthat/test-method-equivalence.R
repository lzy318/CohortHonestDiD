test_that("balanced Imp and CS-NYT have equivalent conditional and C-LF sets", {
  d <- independent_two_cohort_short_panel()
  imp <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "imp")
  cs <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nyt")
  ctl <- list(lf.sims = 128L, grid.points = 31L, endpoint.tol = 1e-7,
              max.evaluations = 10000L)
  cases <- list(
    list(restriction = "SD", M = 0.02, method = "Conditional", full.path = FALSE),
    list(restriction = "SD", M = 0.02, method = "Conditional", full.path = TRUE),
    list(restriction = "SD", M = 0.02, method = "C-LF", full.path = FALSE),
    list(restriction = "SD", M = 0.02, method = "C-LF", full.path = TRUE),
    list(restriction = "RM", M = 0.5, method = "Conditional"),
    list(restriction = "RM", M = 0.5, method = "C-LF"))
  for (spec in cases) {
    label <- paste(unlist(spec), collapse = "/")
    a <- do.call(cohort_sensitivity, c(list(fit = imp, seed = 117, control = ctl), spec))
    b <- do.call(cohort_sensitivity, c(list(fit = cs, seed = 117, control = ctl), spec))
    expect_true(all(a$intervals$success), info = label)
    expect_true(all(b$intervals$success), info = label)
    expect_equal(a$intervals$status, b$intervals$status, info = label)
    # Inversion stops at 1e-7 outcome units; allow five such steps between runs.
    expect_equal(ordered_endpoints(a), ordered_endpoints(b), tolerance = 5e-7, info = label)
  }
})
