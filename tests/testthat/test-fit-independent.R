test_that("CS-NT and CS-NYT equal independently formed sample contrasts", {
  d <- independent_panel()
  for (estimator in c("cs_nt", "cs_nyt")) {
    fit <- cohort_fit(y ~ d, data = d, index = c("id", "time"),
                      estimator = estimator)
    expect_s3_class(fit, "cohort_fit")
    expected <- independent_cs(d, fit$coordinate_info,
                               never_only = estimator == "cs_nt")
    expect_equal(unname(fit$betahat), expected, tolerance = 1e-10)
    expect_equal(fit$Sigma, t(fit$Sigma), tolerance = 1e-12)
    expect_gte(min(eigen(fit$Sigma, symmetric = TRUE, only.values = TRUE)$values),
               -1e-10)
  }
})

test_that("weighted CS estimates use population masses and weighted cell means", {
  d <- independent_panel()
  fit <- cohort_fit(y ~ d, data = d, index = c("id", "time"),
                    estimator = "cs_nt", weights = "w")
  expected <- independent_cs(d, fit$coordinate_info, weights = "w")
  expect_equal(unname(fit$betahat), expected, tolerance = 1e-10)
  expect_equal(unname(fit$cohort_masses[c("4", "6")]),
               c(sum(d$w[d$time == 1 & d$id <= 20]),
                 sum(d$w[d$time == 1 & d$id > 20 & d$id <= 42])),
               tolerance = 1e-12)
})

test_that("formula and named-column interfaces fit the same observations", {
  d <- independent_panel()
  a <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nt")
  b <- cohort_fit(data = d, index = c("id", "time"), Y = "y", D = "d", estimator = "cs_nt")
  expect_identical(a$betahat, b$betahat)
  expect_identical(a$Sigma, b$Sigma)
})

test_that("CS covariance equals independent cell-centered cluster scores", {
  d <- independent_panel()
  for (cl in c("id", "state")) {
    fit <- cohort_fit(y ~ d, data = d, index = c("id", "time"),
                      estimator = "cs_nt", cluster = cl, weights = "w")
    info <- fit$coordinate_info
    group <- ifelse(d$id <= 20, 4, ifelse(d$id <= 42, 6, 0))
    cell <- interaction(group, d$time, drop = TRUE)
    mu <- ave(seq_len(nrow(d)), cell, FUN = function(j) sum(d$w[j] * d$y[j]) / sum(d$w[j]))
    residual <- d$y - mu
    op <- matrix(0, nrow(info), nrow(d))
    for (j in seq_len(nrow(info))) {
      g <- info$cohort[j]
      t <- info$time_index[j]
      for (term in seq_len(4L)) {
        gg <- c(g, g, 0, 0)[term]
        tt <- c(t, g - 1, t, g - 1)[term]
        sign <- c(1, -1, -1, 1)[term]
        k <- which(group == gg & d$time == tt)
        op[j, k] <- op[j, k] + sign * d$w[k] / sum(d$w[k])
      }
    }
    clusters <- unique(d[[cl]])
    score <- vapply(clusters, function(cc) {
      k <- which(d[[cl]] == cc)
      as.numeric(op[, k, drop = FALSE] %*% residual[k])
    }, numeric(nrow(info)))
    n <- nrow(d)
    p <- nlevels(cell)
    nc <- length(clusters)
    expected <- nc / (nc - 1) * (n - 1) / (n - p) * tcrossprod(score)
    expect_equal(unname(fit$Sigma), expected, tolerance = 1e-10)
  }
})

test_that("balanced imputation and CS-NYT implied paths differ by cohort constants", {
  d <- independent_panel()
  imp <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "imp")
  cs <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nyt")
  expect_identical(names(imp$betahat), names(cs$betahat))
  di <- as.vector(solve(imp$W, imp$betahat))
  dc <- as.vector(solve(cs$W, cs$betahat))
  for (g in unique(imp$coordinate_info$cohort)) {
    z <- di[imp$coordinate_info$cohort == g] - dc[cs$coordinate_info$cohort == g]
    expect_lt(diff(range(z)), 1e-10)
  }
})

test_that("argument checks reject data contradictions", {
  d <- independent_panel()
  expect_error(cohort_fit(y ~ d, data = rbind(d, d[1, ]), index = c("id", "time")),
               "duplicat|unique", ignore.case = TRUE)
  bad <- d
  bad$d[bad$id == 1 & bad$time == 5] <- 0
  expect_error(cohort_fit(y ~ d, data = bad, index = c("id", "time")),
               "absorb|revers|switch|monoton", ignore.case = TRUE)
  bad <- d
  bad$w[1] <- -1
  expect_error(cohort_fit(y ~ d, data = bad, index = c("id", "time"), weights = "w"),
               "weight|positive|nonnegative", ignore.case = TRUE)
  expect_error(cohort_fit(y ~ d + w, data = d, index = c("id", "time")),
               "covariat|formula|one|only", ignore.case = TRUE)
  expect_error(cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "bogus"))
})

test_that("unbalanced support is allowed and treatment LOCF does not impute Y", {
  d <- independent_panel()
  # Missing row before first treatment carries forward zero, not future one.
  d <- d[!(d$id == 1 & d$time == 3), ]
  d$d[d$id == 1 & d$time == 5] <- NA
  d$y[d$id == 50 & d$time == 2] <- NA
  expect_warning(fit <- cohort_fit(y ~ d, data = d, index = c("id", "time")),
                 "unbalanc|approxim|missing", ignore.case = TRUE)
  p <- fit$data_info$panel
  row_before <- p[p$unit_id == 1 & p$time_index == 3, ]
  row_after <- p[p$unit_id == 1 & p$time_index == 5, ]
  expect_equal(row_before$D, 0)
  expect_equal(row_after$D, 1)
  expect_true(is.na(row_before$Y))
  expect_false(row_before$observed_row)
  expect_true(is.na(p$Y[p$unit_id == 50 & p$time_index == 2]))
  expect_equal(sum(!is.na(p$Y)), sum(!is.na(d$y)))
  expect_true(all(is.finite(fit$betahat)))
})

test_that("leading missing status is never filled from future treatment", {
  d <- independent_panel()
  d$d[d$id == 1 & d$time == 1] <- NA
  fit <- suppressWarnings(cohort_fit(y ~ d, data = d, index = c("id", "time")))
  p <- fit$data_info$panel
  z <- p[p$unit_id == 1 & p$time_index == 1, ]
  expect_true(is.na(z$D))
  expect_false(z$usable)
})

test_that("unbalanced imputation agrees with an independent weighted fixed-effects fit", {
  d <- independent_panel()
  d <- d[!(d$id == 1 & d$time == 3), ]
  d$d[d$id == 1 & d$time == 5] <- NA
  d$y[d$id == 50 & d$time == 2] <- NA
  expect_warning(fit <- cohort_fit(y ~ d, data = d, index = c("id", "time"),
                                   weights = "w"), "Unbalanced")
  treatment <- d$d
  for (i in unique(d$id)) {
    j <- which(d$id == i)
    j <- j[order(d$time[j])]
    previous <- NA_real_
    for (k in j) {
      if (!is.na(d$d[k])) previous <- d$d[k]
      else treatment[k] <- previous
    }
  }
  untreated <- !is.na(d$y) & !is.na(treatment) & treatment == 0
  reference <- stats::lm(y ~ factor(id) + factor(time), data = d[untreated, ], weights = w)
  predicted <- stats::predict(reference, newdata = d)
  info <- fit$coordinate_info
  for (j in which(info$type == "post")) {
    members <- if (info$cohort[j] == 4) 1:20 else 21:42
    k <- which(d$id %in% members & d$time == info$time_index[j] & !is.na(d$y))
    expected <- sum(d$w[k] * (d$y[k] - predicted[k])) / sum(d$w[k])
    expect_equal(unname(fit$betahat[j]), unname(expected), tolerance = 1e-10)
  }
})

test_that("a genuinely missing required cohort-time cell is rejected", {
  d <- independent_panel()
  d$y[d$id <= 20 & d$time == 2] <- NA
  expect_warning(expect_error(cohort_fit(y ~ d, data = d, index = c("id", "time")),
               "estimab|cell|support|observ", ignore.case = TRUE), "Unbalanced")
})
