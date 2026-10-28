#' Estimate Cohort Coefficients for Staggered Adoption
#'
#' Fits an imputation or Callaway--Sant'Anna difference-in-differences estimator
#' and constructs the joint covariance and block-bias map used by sensitivity
#' analysis. Outcome missingness is allowed: outcomes are never imputed.
#'
#' @param formula A two-variable formula such as \code{outcome ~ treatment}.
#'   Covariate terms and transformed variables are not supported by this interface.
#' @param data A data frame of unit-period observations.
#' @param index Character vector naming the unit and time columns, in that order.
#' @param Y,D Outcome and treatment column names, as an alternative to
#'   \code{formula}.
#' @param estimator One of \code{"imp"}, \code{"cs_nyt"}, or \code{"cs_nt"}.
#' @param cluster Cluster column name or row-aligned vector, constant within unit.
#'   Defaults to the unit identifier.
#' @param weights Positive weights, supplied as a column name or row-aligned
#'   numeric vector. Weights must be constant within unit; the default is one.
#'   Extreme weights are rescaled by a common positive factor when required for
#'   finite arithmetic. \code{data_info$weight_scale} records that divisor;
#'   stored unit weights and cohort masses then use the rescaled units.
#' @param covariance Covariance estimator: \code{"CR1"} (default) or
#'   \code{"CR2"}, case-insensitive. It follows \code{...}, so give it by name.
#'   See Details.
#' @param ... Reserved for future extensions. Unrecognized arguments cause an error.
#' @return A \code{cohort_fit} object containing full named coefficients
#'   \code{betahat}, their joint covariance \code{Sigma}, the bias map \code{W},
#'   normalization constraints, a cohort-period registry, sample information,
#'   and diagnostic tables. \code{covariance_convention} is
#'   \code{"cohort_period_centered_CR1"} or \code{"cohort_period_centered_CR2"}.
#'   \code{cluster_scores} has one row per cluster: under CR1, \code{Sigma} is
#'   the CR1 factor times their cross-product; under CR2 the rows are the
#'   adjusted scores and \code{Sigma} is their cross-product. CR2 fits also
#'   contain \code{cluster_shares}, a data frame with one row per cluster and
#'   cohort group (\code{cluster}, \code{cohort}, \code{units}, \code{mass},
#'   \code{share}, \code{factor}; cohort 0 is the never-treated group), and
#'   \code{diagnostics$cr1} is \code{NA} because no CR1 factor is applied.
#' @details Treatment must be binary and absorbing. Missing treatment values,
#'   including absent unit-time rows, are filled forward from the previous
#'   observed treatment status within unit. Leading unknown treatment statuses
#'   remain unknown; their outcomes are not used. At least one observed untreated
#'   period is required for every treated unit, and a never-treated group is
#'   required. An unbalanced sample uses actual observed histories for estimation;
#'   the balanced-panel bias map based on fixed cohort masses is an approximation.
#'   Estimator equivalence and the paper's coverage guarantees are not established
#'   for that approximation. A warning and diagnostic flag make this explicit.
#'
#'   The covariance centers outcomes within cohort-period cells (weighted cell
#'   means when weights are supplied) and sums the resulting scores within
#'   clusters. With \code{covariance = "CR1"}, the default, the sum of the outer
#'   products of the cluster scores is multiplied by the CR1 factor
#'   G/(G - 1) * (n - 1)/(n - p), with G clusters, n observed usable rows, and
#'   p nonempty cohort-period cells. \code{covariance = "CR2"} is the
#'   bias-reduced covariance of Bell and McCaffrey (2002): there is no such
#'   factor, and each cluster's score is instead divided by sqrt(1 - p_c), where
#'   p_c is the cluster's share of the units of its cohort, measured by the unit
#'   weights when weights are supplied. The never-treated units form their own
#'   group. In the cohort-period cell-mean model this is the CR2 adjustment for
#'   a working model in which each unit's disturbance variance is inversely
#'   proportional to its weight. On a balanced panel, a cluster that makes up a
#'   whole cohort has a zero score under both rules, so that cohort's sampling
#'   variation is omitted. If a cluster contains units of several cohorts, the
#'   part of its score that comes from each cohort is divided by sqrt(1 - p),
#'   with p its share of that cohort; for clusters nested in cohorts this is
#'   the rule above. On unbalanced panels the CR2 adjustment, like the bias
#'   map, is an approximation: the shares are fixed by unit and do not follow
#'   the period-by-period composition of the cells, and the estimator is no
#'   longer a function of the cell means alone. There, a cluster that makes up
#'   a whole cohort still has a zero score under CR2, but if that cohort has
#'   missing outcomes, its CR1 score is generally not zero because the
#'   cluster's centered residuals no longer cancel in the estimator. The paper
#'   recommends CR2 with Conditional or C-LF inference when some cohorts
#'   contain only a few units or clusters. Few treated clusters can make either
#'   covariance unreliable even when the numerical calculation succeeds.
#' @references Bell, R. M. and McCaffrey, D. F. (2002). Bias reduction in
#'   standard errors for linear regression with multi-stage samples. Survey
#'   Methodology, 28(2), 169-181.
#' @export
#' @examples
#' data(minimum_wage)
#' fit <- cohort_fit(Y ~ D, data = minimum_wage,
#'                   index = c("id", "year"), estimator = "cs_nt")
#' head(fit$betahat)
#' fit_cr2 <- cohort_fit(Y ~ D, data = minimum_wage, index = c("id", "year"),
#'                       estimator = "cs_nt", covariance = "CR2")
#' fit_cr2
cohort_fit <- function(formula = NULL, data, index, Y = NULL, D = NULL,
                       estimator = "imp", cluster = NULL, weights = NULL, ...,
                       covariance = "CR1") {
  call <- match.call()
  dots <- list(...)
  .hc_assert(!length(dots), paste0("Unrecognized arguments in ...: ",
                                   paste(names(dots), collapse = ", ")))
  estimator <- .hc_estimator(estimator)
  covariance <- .hc_covariance_type(covariance)
  prepared <- .hc_prepare_data(formula, data, index, Y, D, cluster, weights, covariance)
  info <- .hc_coordinate_info(prepared$cohort_times, prepared$time_labels)
  work <- .hc_working_data(prepared)
  .hc_check_cells(work$data, info)
  operator <- switch(estimator,
    imp = .hc_estimate_imp(work$data, info),
    cs_nyt = .hc_estimate_cs_nyt(work$data, info),
    cs_nt = .hc_estimate_cs_nt(work$data, info))
  beta <- .hc_apply_operator(operator, work$data$Y)
  names(beta) <- info$coordinate
  vc <- .hc_covariance(operator, prepared, work, info, covariance)
  normalization <- .hc_normalization(info, estimator)
  W <- build_W(prepared$cohort_times, prepared$cohort_masses,
               prepared$never_treated_mass, length(prepared$time_labels), estimator)
  W <- W[info$coordinate, info$coordinate, drop = FALSE]
  diagnostics <- .hc_fit_diagnostics(prepared, info, vc, normalization, beta)
  result <- list(call = call, estimator = estimator,
    betahat = beta, Sigma = vc$Sigma, W = W,
    normalization = normalization, coordinate_info = info,
    cohort_sizes = prepared$cohort_sizes, cohort_masses = prepared$cohort_masses,
    never_treated_mass = prepared$never_treated_mass,
    total_periods = length(prepared$time_labels), diagnostics = diagnostics,
    data_info = list(panel = prepared$panel, unit_registry = prepared$unit_registry,
      time_labels = prepared$time_labels, index = index,
      outcome = prepared$outcome_name, treatment = prepared$treatment_name,
      n_units = nrow(prepared$unit_registry), n_observations = nrow(data),
      n_usable = sum(prepared$panel$usable), balanced = prepared$balanced,
      weight_scale = prepared$weight_scale,
      missingness = prepared$missingness),
    cluster_scores = vc$scores,
    covariance_convention = paste0("cohort_period_centered_", covariance),
    specification = list(formula = formula, cluster = cluster, weights = weights,
      treatment_fill = "forward_from_previous_observed_status",
      bias_map = if (prepared$balanced) "balanced_exact" else "balanced_approximation"))
  if (covariance == "CR2") result$cluster_shares <- vc$cluster_shares
  class(result) <- "cohort_fit"
  result
}

.hc_prepare_data <- function(formula, data, index, Y, D, cluster, weights, covariance = "CR1") {
  .hc_assert(is.data.frame(data) && nrow(data) > 0L,
             "data must be a nonempty data frame.")
  .hc_assert(is.character(index) && length(index) == 2L && !anyNA(index) &&
               !anyDuplicated(index) && all(index %in% names(data)),
             "index must name the unit and time columns, in that order.")
  if (!is.null(formula)) {
    .hc_assert(inherits(formula, "formula") && length(formula) == 3L &&
                 is.symbol(formula[[2L]]) && is.symbol(formula[[3L]]),
      "Use a two-variable formula such as Y ~ D; covariates and transformations are not supported here.")
    .hc_assert(is.null(Y) && is.null(D), "Supply formula or Y and D column names, not both.")
    Y <- as.character(formula[[2L]]); D <- as.character(formula[[3L]])
  }
  y <- .hc_column(Y, data, "Y")
  d <- .hc_column(D, data, "D")
  .hc_assert(is.numeric(y) && all(is.finite(y) | is.na(y)),
             "The outcome must be numeric, with finite values or NA.")
  .hc_assert((is.numeric(d) || is.logical(d)) && all(is.na(d) | d %in% c(0, 1)),
             "Treatment must be binary (0/1), with missing values allowed.")
  id <- data[[index[1L]]]; time <- data[[index[2L]]]
  .hc_assert(!anyNA(id) && !anyNA(time), "Unit and time identifiers cannot be missing.")
  .hc_assert(!is.list(id) && !is.list(time), "Unit and time identifiers must be atomic vectors.")
  .hc_assert(!anyDuplicated(data[index]), "Duplicate unit-time rows are not allowed.")
  time_labels <- sort(unique(time))
  n_periods <- length(time_labels)
  .hc_assert(n_periods >= 2L, "At least two distinct periods are required.")
  if (is.numeric(time_labels) || inherits(time_labels, "Date") || inherits(time_labels, "POSIXt")) {
    gaps <- diff(as.numeric(time_labels))
    .hc_assert(all(is.finite(gaps) & gaps > 0), "Time values must be finite and ordered.")
    if (length(gaps) > 1L && max(abs(gaps - gaps[1L])) > 1e-8 * max(1, abs(gaps[1L]))) {
      warning("Time labels are not equally spaced. Relative periods and curvature bounds use the ordered observed time grid, not elapsed calendar time.", call. = FALSE)
    }
  }
  ids <- sort(unique(as.character(id)))
  ui <- match(as.character(id), ids); ti <- match(time, time_labels)
  w <- if (is.null(weights)) rep(1, nrow(data)) else .hc_column(weights, data, "weights", TRUE)
  .hc_assert(is.numeric(w) && all(is.finite(w) & w > 0), "weights must be positive and finite.")
  unit_weight <- as.numeric(.hc_unit_values(w, ui, "weights"))
  weight_scale <- 1
  weighted_outcome <- w[!is.na(y)] * y[!is.na(y)]
  if (!is.finite(sum(w)) || !is.finite(sum(abs(weighted_outcome))) ||
      any(weighted_outcome == 0 & y[!is.na(y)] != 0) ||
      max(unit_weight) < .Machine$double.xmin) {
    weight_scale <- max(unit_weight)
    unit_weight <- unit_weight / weight_scale
    .hc_assert(all(unit_weight > 0),
               "The range of weights is too large to preserve all positive weights after rescaling.")
  }
  cl <- if (is.null(cluster)) as.character(id) else .hc_column(cluster, data, "cluster", TRUE)
  .hc_assert(is.atomic(cl), "Cluster identifiers must be an atomic vector.")
  unit_cluster <- as.character(.hc_unit_values(cl, ui, "cluster"))
  .hc_assert(length(unique(unit_cluster)) >= 2L, paste0("At least two independent clusters are required for ", covariance, " covariance."))
  n_units <- length(ids)
  full <- data.frame(unit_index = rep(seq_len(n_units), each = n_periods),
    unit_id = rep(ids, each = n_periods), time_index = rep(seq_len(n_periods), n_units),
    calendar_time = rep(time_labels, n_units), Y = NA_real_, D = NA_real_,
    D_original = NA_real_, observed_row = FALSE, treatment_filled = FALSE,
    stringsAsFactors = FALSE)
  rows <- (ui - 1L) * n_periods + ti
  full$Y[rows] <- y; full$D[rows] <- as.numeric(d); full$D_original[rows] <- as.numeric(d)
  full$observed_row[rows] <- TRUE
  cohort <- integer(n_units)
  for (i in seq_len(n_units)) {
    rr <- (i - 1L) * n_periods + seq_len(n_periods)
    observed <- full$D_original[rr]
    known <- which(!is.na(observed))
    .hc_assert(length(known) > 0L, paste0("Unit ", ids[i], " has no observed treatment status."))
    .hc_assert(!any(diff(observed[known]) < 0), paste0("Treatment reverses from 1 to 0 for unit ", ids[i], ". Absorbing treatment is required."))
    .hc_assert(observed[known[1L]] == 0, paste0("Unit ", ids[i], " is treated at its first known status. Always-treated or left-censored adoption is not supported."))
    previous <- NA_real_
    for (j in seq_along(rr)) {
      if (!is.na(observed[j])) previous <- observed[j]
      else if (!is.na(previous)) {
        full$D[rr[j]] <- previous; full$treatment_filled[rr[j]] <- TRUE
      }
    }
    first <- which(observed == 1)
    cohort[i] <- if (length(first)) first[1L] else 0L
  }
  .hc_assert(any(cohort == 0L), "At least one never-treated unit is required.")
  .hc_assert(any(cohort > 0L), "At least one treated unit is required.")
  full$cohort <- cohort[full$unit_index]
  full$unit_weight <- unit_weight[full$unit_index]
  full$cluster_id <- unit_cluster[full$unit_index]
  full$usable <- !is.na(full$Y) & !is.na(full$D)
  .hc_assert(any(full$usable), "No outcome observations have known treatment status.")
  usable_count <- tabulate(full$unit_index[full$usable], nbins = n_units)
  .hc_assert(all(usable_count > 0L), paste0("Every unit must have at least one usable outcome. Units without one: ", paste(ids[usable_count == 0L], collapse = ", ")))
  registry <- data.frame(unit_index = seq_len(n_units), unit_id = ids, cohort = cohort,
    unit_weight = unit_weight, cluster_id = unit_cluster,
    adoption_time = time_labels[pmax(cohort, 1L)], stringsAsFactors = FALSE)
  registry$adoption_time[cohort == 0L] <- NA
  cohort_times <- sort(unique(cohort[cohort > 0L]))
  sizes <- setNames(vapply(cohort_times, function(g) sum(cohort == g), integer(1)), as.character(cohort_times))
  masses <- setNames(vapply(cohort_times, function(g) sum(unit_weight[cohort == g]), numeric(1)), as.character(cohort_times))
  balanced <- all(full$observed_row & full$usable)
  missingness <- list(absent_rows = sum(!full$observed_row),
    missing_outcomes_observed = sum(full$observed_row & is.na(full$Y)),
    treatment_filled = sum(full$treatment_filled),
    treatment_filled_observed_rows = sum(full$treatment_filled & full$observed_row),
    treatment_filled_absent_rows = sum(full$treatment_filled & !full$observed_row),
    leading_unknown_treatment = sum(is.na(full$D)),
    outcomes_excluded_unknown_treatment = sum(!is.na(full$Y) & is.na(full$D)),
    outcomes_imputed = 0L)
  if (!balanced) warning("Unbalanced sample: estimation uses observed outcomes and actual available histories. W uses the balanced-panel formula with fixed cohort masses; estimator equivalence and sensitivity-analysis coverage are not guaranteed for this approximation.", call. = FALSE)
  list(panel = full, unit_registry = registry, time_labels = time_labels,
    cohort_times = cohort_times, cohort_sizes = sizes, cohort_masses = masses,
    weight_scale = weight_scale,
    never_treated_mass = sum(unit_weight[cohort == 0L]), balanced = balanced,
    missingness = missingness, outcome_name = Y, treatment_name = D)
}

.hc_working_data <- function(prepared) {
  observed <- prepared$panel[prepared$panel$usable, , drop = FALSE]
  if (!prepared$balanced) return(list(data = observed, compressed = FALSE))
  groups <- c(0L, prepared$cohort_times)
  n_periods <- length(prepared$time_labels)
  key <- (match(observed$cohort, groups) - 1L) * n_periods + observed$time_index
  mass <- as.numeric(rowsum(observed$unit_weight, key, reorder = TRUE))
  means <- as.numeric(rowsum(observed$unit_weight * observed$Y, key, reorder = TRUE)) / mass
  working <- data.frame(unit_index = rep(seq_along(groups), each = n_periods),
    cohort = rep(groups, each = n_periods), time_index = rep(seq_len(n_periods), length(groups)),
    Y = means, unit_weight = mass)
  working$D <- as.integer(working$cohort > 0L & working$time_index >= working$cohort)
  list(data = working, compressed = TRUE, observed_to_working = key,
       cell_mass = mass)
}

#' County Minimum-Wage Panel Used in the Study
#'
#' The study's seven-year balanced county panel, with convenience aliases for
#' the outcome, absorbing treatment, county identifier, and state identifier.
#' Original study variables are retained without alteration.
#' @format A data frame with 15,988 rows and 28 columns:
#' \describe{
#'   \item{id, countyreal, FIPS}{Numeric county identifiers.}
#'   \item{year}{Calendar year, from 2001 through 2007.}
#'   \item{Y, lemp}{Log employment, the study outcome.}
#'   \item{D}{Absorbing treatment indicator: one from the first treatment year.}
#'   \item{state}{Numeric state FIPS identifier, usable for state clustering.}
#'   \item{cohort, first.treat}{First treatment calendar year; zero for never treated.}
#'   \item{state_name, county_name}{State and county names.}
#'   \item{emp0A01_BS}{Employment count underlying the log outcome.}
#'   \item{quarter}{Quarter of the underlying observation.}
#'   \item{censusdiv, region}{Census division and region factors.}
#'   \item{msa}{Indicator for metropolitan statistical area status.}
#'   \item{pop, lpop}{Population measure and its natural logarithm.}
#'   \item{white, black}{White and Black population shares in the study input.}
#'   \item{hs, col}{High-school and college attainment measures.}
#'   \item{medinc, lmedinc}{Median income measure and its natural logarithm.}
#'   \item{pov}{Poverty measure.}
#'   \item{nssi}{Study input's nssi demographic variable.}
#'   \item{treat}{Indicator for ever belonging to a treated cohort; this is not D.}
#' }
#' @source County panel of the minimum-wage application in Callaway and
#'   Sant'Anna (2021), "Difference-in-Differences with Multiple Time Periods",
#'   Journal of Econometrics, 225(2), 200--230, as used in Liu (2026),
#'   "Cohort-Anchored HonestDiD" (arXiv:2509.01829). The script
#'   \file{data-raw/prepare-minimum-wage.R} in the package's source repository
#'   records all convenience aliases.
#' @usage data(minimum_wage)
"minimum_wage"
