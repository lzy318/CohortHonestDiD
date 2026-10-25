# CR2 references are computed here from the data and from CR1 fits: a
# cluster's share p is its units' weight (count when unweighted) divided by
# that of its cohort group, the never-treated units forming one group.

cr2_units <- function(d, cluster = NULL, weights = NULL) {
  ids <- sort(unique(d$id))
  first <- match(ids, d$id)
  group <- vapply(ids, function(i) {
    tt <- d$time[d$id == i & d$d == 1L]
    if (length(tt)) min(tt) else 0
  }, numeric(1))
  data.frame(id = ids, group = group,
             cluster = if (is.null(cluster)) as.character(ids) else as.character(d[[cluster]][first]),
             weight = if (is.null(weights)) 1 else d[[weights]][first],
             stringsAsFactors = FALSE)
}

# Factor 1/sqrt(1 - p) of the units in `keep` (one cluster, or one cluster's
# part of one group); zero when they make up their whole group.
cr2_factor <- function(u, keep) {
  g <- unique(u$group[keep])
  stopifnot(length(g) == 1L)
  in_group <- u$group == g
  if (all(keep[in_group])) return(0)
  1 / sqrt(1 - sum(u$weight[keep]) / sum(u$weight[in_group]))
}

expect_rel_close <- function(actual, expected, tolerance = 1e-10) {
  actual <- unname(as.matrix(actual)); expected <- unname(as.matrix(expected))
  expect_identical(dim(actual), dim(expected))
  expect_lte(max(abs(actual - expected)), tolerance * max(abs(expected)))
}

expect_cr2_nested <- function(cr1, cr2, u) {
  s <- cr1$cluster_scores
  a <- vapply(rownames(s), function(cc) cr2_factor(u, u$cluster == cc), numeric(1))
  expect_identical(rownames(cr2$cluster_scores), rownames(s))
  expect_identical(colnames(cr2$cluster_scores), colnames(s))
  expect_rel_close(cr2$cluster_scores, s * a)
  expect_rel_close(cr2$Sigma, crossprod(s * a))
  cp <- crossprod(cr2$cluster_scores)
  expect_rel_close(cr2$Sigma, (cp + t(cp)) / 2, tolerance = 1e-14)
  expect_true(isSymmetric(unname(cr2$Sigma), tol = 0))
  invisible(a)
}

uneven_clusters <- function(id) {
  ifelse(id <= 11, "a", ifelse(id <= 20, paste0("s", id),
    ifelse(id <= 24, "b", ifelse(id <= 26, "b2", ifelse(id <= 42, paste0("s", id),
      ifelse(id <= 45, "c", paste0("s", id)))))))
}

fit_pair <- function(d, estimator, cluster = NULL, weights = NULL) {
  base <- list(y ~ d, data = d, index = c("id", "time"), estimator = estimator,
               cluster = cluster, weights = weights)
  list(cr1 = do.call(cohort_fit, base),
       cr2 = do.call(cohort_fit, c(base, list(covariance = "CR2"))))
}

test_that("the covariance option is validated and CR1 is the unchanged default", {
  d <- independent_panel()
  for (bad in list("CR3", c("CR1", "CR2"), NA_character_, 2)) {
    expect_error(cohort_fit(y ~ d, data = d, index = c("id", "time"), covariance = bad),
                 "covariance must be 'CR1' or 'CR2'", fixed = TRUE)
  }
  expect_error(cohort_fit(y ~ d, data = d, index = c("id", "time"), cov = "CR2"),
               "Unrecognized arguments in ...: cov", fixed = TRUE)
  lower <- cohort_fit(y ~ d, data = d, index = c("id", "time"), covariance = "cr2")
  expect_identical(lower$covariance_convention, "cohort_period_centered_CR2")
  default <- cohort_fit(y ~ d, data = d, index = c("id", "time"), cluster = "state",
                        weights = "w")
  cr1 <- cohort_fit(y ~ d, data = d, index = c("id", "time"), cluster = "state",
                    weights = "w", covariance = "CR1")
  for (field in c("betahat", "Sigma", "W", "cluster_scores", "diagnostics",
                  "covariance_convention", "normalization", "coordinate_info")) {
    expect_identical(cr1[[field]], default[[field]], info = field)
  }
  expect_identical(default$covariance_convention, "cohort_period_centered_CR1")
  expect_null(default$cluster_shares)
  expect_null(cr1$cluster_shares)
  expect_identical(setdiff(names(lower), names(default)), "cluster_shares")
  expect_identical(names(lower)[seq_along(names(default))], names(default))
})

test_that("CR2 rescales nested cluster scores by 1/sqrt(1 - share)", {
  d <- independent_panel()
  d$uneven <- uneven_clusters(d$id)
  cases <- list(
    list(estimator = "imp", cluster = NULL, weights = NULL),
    list(estimator = "cs_nt", cluster = NULL, weights = "w"),
    list(estimator = "imp", cluster = "state", weights = NULL),
    list(estimator = "cs_nyt", cluster = "uneven", weights = NULL),
    list(estimator = "imp", cluster = "uneven", weights = "w"))
  for (k in seq_along(cases)) {
    z <- cases[[k]]
    fits <- fit_pair(d, z$estimator, z$cluster, z$weights)
    u <- cr2_units(d, z$cluster, z$weights)
    a <- expect_cr2_nested(fits$cr1, fits$cr2, u)
    expect_true(all(a > 1))
    expect_identical(fits$cr2$betahat, fits$cr1$betahat)
    expect_true(is.na(fits$cr2$diagnostics$cr1))
    expect_true(is.finite(fits$cr1$diagnostics$cr1))
    expect_identical(fits$cr2$diagnostics$n_clusters, fits$cr1$diagnostics$n_clusters)
    expect_identical(fits$cr2$diagnostics$covariance_rank, fits$cr1$diagnostics$covariance_rank)
  }
  # A cluster holding more than half of its cohort.
  u <- cr2_units(d, "uneven")
  expect_equal(cr2_factor(u, u$cluster == "a"), 1 / sqrt(9 / 20), tolerance = 1e-14)
})

test_that("a whole-cohort cluster keeps a zero CR2 score", {
  d <- independent_panel()
  d$whole <- ifelse(d$id <= 20, "whole", paste0("u", d$id))
  fits <- suppressWarnings(fit_pair(d, "imp", "whole", "w"))
  expect_s3_class(fits$cr2, "cohort_fit")
  sh <- fits$cr2$cluster_shares
  row <- sh[sh$cluster == "whole", ]
  expect_identical(nrow(row), 1L)
  expect_identical(row$cohort, 4L)
  expect_identical(row$units, 20L)
  expect_identical(row$share, 1)
  expect_identical(row$factor, 0)
  expect_true(all(fits$cr2$cluster_scores["whole", ] == 0))
  expect_true(all(fits$cr1$cluster_scores["whole", ] == 0))
  expect_cr2_nested(fits$cr1, fits$cr2, cr2_units(d, "whole", "w"))
  expect_true(all(is.finite(fits$cr2$Sigma)))
})

test_that("clusters spanning cohorts scale each cohort's part of their score", {
  d <- independent_panel()
  d$span <- ceiling((d$id + 1) / 2)
  for (weights in list(NULL, "w")) {
    u <- cr2_units(d, "span", weights)
    group <- u$group[match(d$id, u$id)]
    d$pair <- paste(d$span, group, sep = "_")
    by_pair <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "imp",
                          cluster = "pair", weights = weights)
    cr2 <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "imp",
                      cluster = "span", weights = weights, covariance = "CR2")
    sp <- by_pair$cluster_scores
    unit_pair <- paste(u$cluster, u$group, sep = "_")
    a <- vapply(rownames(sp), function(k) cr2_factor(u, unit_pair == k), numeric(1))
    owner <- sub("_.*$", "", rownames(sp))
    reference <- rowsum(sp * a, owner)[rownames(cr2$cluster_scores), , drop = FALSE]
    expect_setequal(unique(owner), rownames(cr2$cluster_scores))
    expect_rel_close(cr2$cluster_scores, reference)
    expect_rel_close(cr2$Sigma, crossprod(reference))
    sh <- cr2$cluster_shares
    expect_identical(sh$cohort[sh$cluster == "11"], c(4L, 6L))
    expect_identical(sh$cohort[sh$cluster == "22"], c(0L, 6L))
    expect_identical(nrow(sh), nrow(sp))
  }
})

test_that("CR2 uses registry shares on unbalanced panels", {
  d <- independent_panel()
  d$uneven <- uneven_clusters(d$id)
  d$whole <- ifelse(d$id <= 20, "whole", paste0("u", d$id))
  drop <- (d$id %in% c(50, 51) & d$time == 2) | (d$id == 5 & d$time == 7) |
    (d$id == 30 & d$time == 1)
  d <- d[!drop, ]
  fits <- suppressWarnings(fit_pair(d, "cs_nt", "uneven", "w"))
  expect_false(fits$cr2$data_info$balanced)
  expect_cr2_nested(fits$cr1, fits$cr2, cr2_units(d, "uneven", "w"))
  whole <- suppressWarnings(fit_pair(d, "imp", "whole", NULL))
  expect_cr2_nested(whole$cr1, whole$cr2, cr2_units(d, "whole", NULL))
  expect_true(all(whole$cr2$cluster_scores["whole", ] == 0))
  # Under CR1 the whole cohort's score is not zero: that cohort has a missing
  # outcome (id 5, time 7). Under CR2 it is exactly zero (documented behavior).
  expect_gt(max(abs(whole$cr1$cluster_scores["whole", ])),
            1e-8 * max(abs(whole$cr1$cluster_scores)))
})

test_that("balanced imputation and CS-NYT FLCIs agree under CR2", {
  d <- independent_panel()
  imp <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "imp",
                    covariance = "CR2")
  cs <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nyt",
                   covariance = "CR2")
  a <- cohort_sensitivity(imp, M = c(0, 0.01), method = "FLCI", seed = 31)
  b <- cohort_sensitivity(cs, M = c(0, 0.01), method = "FLCI", seed = 31)
  expect_true(all(a$intervals$success))
  expect_true(all(b$intervals$success))
  expect_equal(ordered_endpoints(a), ordered_endpoints(b), tolerance = 1e-7)
})

test_that("printed fits and results name the covariance", {
  d <- independent_short_panel()
  cr1 <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nt")
  cr2 <- cohort_fit(y ~ d, data = d, index = c("id", "time"), estimator = "cs_nt",
                    covariance = "CR2")
  line1 <- "Covariance: CR1 cluster-robust, 24 clusters "
  line2 <- "Covariance: CR2 bias-reduced cluster-robust, 24 clusters "
  expect_identical(capture.output(print(cr1))[3], line1)
  expect_identical(capture.output(print(cr2))[3], line2)
  expect_identical(capture.output(print(summary(cr1)))[2], line1)
  expect_identical(capture.output(print(summary(cr2)))[2], line2)
  expect_identical(summary(cr2)$covariance, "CR2 bias-reduced cluster-robust, 24 clusters")
  ans <- cohort_sensitivity(cr2, M = 0.01, method = "FLCI", seed = 3)
  expect_identical(capture.output(print(ans))[3], line2)
  out <- capture.output(print(summary(ans)))
  expect_identical(out[1:3], c(line2, "", "Conventional estimates"))
  external <- cohort_input(betahat = cr2$betahat, vcov = cr2$Sigma,
                           coefficient_info = cr2$coordinate_info, W = cr2$W,
                           normalization = cr2$normalization,
                           cohort_sizes = cr2$cohort_sizes, estimator = "cs_nt")
  line3 <- "Covariance: supplied to cohort_input() "
  expect_identical(capture.output(print(external))[3], line3)
  ext <- cohort_sensitivity(external, M = 0.01, method = "FLCI", seed = 3)
  expect_identical(capture.output(print(ext))[3], line3)
  expect_equal(ordered_endpoints(ext), ordered_endpoints(ans), tolerance = 1e-9)
})

test_that("cluster_shares records each cluster-cohort pair", {
  d <- independent_panel()
  d$span <- ceiling((d$id + 1) / 2)
  fit <- cohort_fit(y ~ d, data = d, index = c("id", "time"), cluster = "span",
                    weights = "w", covariance = "CR2")
  sh <- fit$cluster_shares
  expect_s3_class(sh, "data.frame")
  expect_identical(names(sh), c("cluster", "cohort", "units", "mass", "share", "factor"))
  expect_identical(vapply(sh, typeof, character(1)),
                   c(cluster = "character", cohort = "integer", units = "integer",
                     mass = "double", share = "double", factor = "double"))
  expect_identical(rownames(sh), as.character(seq_len(nrow(sh))))
  u <- cr2_units(d, "span", "w")
  key <- unique(paste(u$cluster, u$group))
  expect_identical(nrow(sh), length(key))
  expect_identical(unique(sh$cluster), rownames(fit$cluster_scores))
  expect_true(all(diff(match(sh$cluster, rownames(fit$cluster_scores))) >= 0))
  for (k in seq_len(nrow(sh))) {
    keep <- u$cluster == sh$cluster[k] & u$group == sh$cohort[k]
    expect_identical(sh$units[k], sum(keep))
    expect_equal(sh$mass[k], sum(u$weight[keep]), tolerance = 1e-14)
    expect_equal(sh$share[k], sum(u$weight[keep]) / sum(u$weight[u$group == sh$cohort[k]]),
                 tolerance = 1e-14)
  }
  expect_true(all(sh$share > 0 & sh$share <= 1))
  expect_true(all(sh$share < 1))
  expect_equal(sh$factor, 1 / sqrt(1 - sh$share), tolerance = 1e-14)
  expect_identical(sort(unique(sh$cohort)), c(0L, 4L, 6L))
})

test_that("CR2 requires at least two clusters", {
  d <- independent_short_panel()
  d$one <- 1
  expect_error(cohort_fit(y ~ d, data = d, index = c("id", "time"), cluster = "one",
                          covariance = "CR2"),
               "At least two independent clusters are required for CR2 covariance.",
               fixed = TRUE)
  expect_error(cohort_fit(y ~ d, data = d, index = c("id", "time"), cluster = "one"),
               "At least two independent clusters are required for CR1 covariance.",
               fixed = TRUE)
})
