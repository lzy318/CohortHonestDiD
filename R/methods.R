#' @export
coef.cohort_fit <- function(object, ...) object$betahat

#' @export
vcov.cohort_fit <- function(object, ...) object$Sigma

.hc_covariance_label <- function(fit) {
  convention <- fit$covariance_convention
  type <- if (identical(convention, "cohort_period_centered_CR1")) "CR1 cluster-robust" else
    if (identical(convention, "cohort_period_centered_CR2")) "CR2 bias-reduced cluster-robust" else NULL
  if (is.null(type)) {
    if (is.list(fit$data_info) && identical(fit$data_info$source, "external"))
      return("supplied to cohort_input()")
    return("not recorded")
  }
  clusters <- if (is.list(fit$diagnostics)) fit$diagnostics$n_clusters else NULL
  if (is.numeric(clusters) && length(clusters) == 1L && is.finite(clusters))
    paste0(type, ", ", as.integer(clusters), " clusters") else type
}

#' @export
print.cohort_fit <- function(x, ...) {
  cat("Cohort difference-in-differences fit\n")
  cat("Estimator:", x$estimator, "\n")
  cat("Covariance:", .hc_covariance_label(x), "\n")
  cat("Cohorts:", length(unique(x$coordinate_info$cohort)),
      " | Coefficients:", length(x$betahat), "\n")
  if (isFALSE(x$data_info$balanced)) cat("Unbalanced panel: balanced-design W approximation.\n")
  tab <- x$diagnostics$cohort_support
  if (is.null(tab)) tab <- x$diagnostics$cohort_table
  if (is.null(tab)) tab <- x$diagnostics$cohort_counts
  if (!is.null(tab)) print(tab, row.names = FALSE)
  invisible(x)
}

#' @export
summary.cohort_fit <- function(object, ...) {
  se <- sqrt(pmax(0, diag(object$Sigma)))
  structure(list(estimator = object$estimator,
                 coefficients = cbind(object$coordinate_info, estimate = as.numeric(object$betahat), se = se),
                 diagnostics = object$diagnostics, data_info = object$data_info,
                 covariance = .hc_covariance_label(object)),
            class = "summary_cohort_fit")
}

#' @export
print.summary_cohort_fit <- function(x, ...) {
  cat("Cohort fit:", x$estimator, "\n")
  if (!is.null(x$covariance)) cat("Covariance:", x$covariance, "\n")
  print(x$coefficients, row.names = FALSE)
  cat("\nDiagnostics:\n"); print(x$diagnostics)
  invisible(x)
}

.hc_balanced_label <- function(s) {
  w <- s$window
  periods <- if (w[["upper"]] == 1L) "target period 1" else paste0("target periods 1 to ", w[["upper"]])
  paths <- if (identical(s$framework, "cohort")) paste0("cohort paths: ", s$cohort.paths) else
    "aggregate path on the window"
  paste0("Balanced window: relative periods ", w[["lower"]], " to ", w[["upper"]], "; ", periods, "; ", paths)
}

#' @export
print.cohort_sensitivity <- function(x, ...) {
  s <- x$settings
  method <- switch(s$method, flci = "FLCI", conditional = "Conditional", c_lf = "C-LF")
  scope <- if (s$restriction == "RM") {
    if (s$framework == "cohort") "cohort-specific pre-period benchmarks" else "aggregate pre-period benchmark"
  } else if (s$full.path) "full path" else "post-only"
  cat(sprintf("%s inference | %s | %s | %s\n", s$framework, s$restriction, method, scope))
  cat(sprintf("Target: %s | %.1f%% pointwise confidence sets\n", s$target, 100 * s$conf.level))
  cat("Covariance:", .hc_covariance_label(x$fit), "\n")
  if (isTRUE(s$balanced)) cat(.hc_balanced_label(s), "\n", sep = "")
  print(x$intervals[, c("target", "M", "component", "lower", "upper", "status")], row.names = FALSE)
  if (s$method == "conditional") .hc_print_search_ranges(.hc_search_ranges(x))
  invisible(x)
}

.hc_search_ranges <- function(object) {
  keys <- unique(object$intervals[, c("target", "M"), drop = FALSE])
  rows <- lapply(seq_along(object$results), function(k) {
    d <- object$results[[k]]$diagnostics
    domain <- d$search.domain
    if (length(domain) != 2L) domain <- c(NA_real_, NA_real_)
    data.frame(target = keys$target[k], M = keys$M[k],
      searched.lower = domain[1L], searched.upper = domain[2L],
      domain.stable = isTRUE(d$domain.stable),
      resolution.stable = if (length(d$resolution.stable)) d$resolution.stable else NA)
  })
  do.call(rbind, rows)
}

.hc_print_search_ranges <- function(ranges) {
  cat("\nConditional search ranges\n")
  print(ranges, row.names = FALSE)
  cat(.hc_conditional_note(), "\n")
}

#' @export
summary.cohort_sensitivity <- function(object, ...) {
  structure(list(settings = object$settings, estimates = object$estimates,
                 intervals = object$intervals, diagnostics = object$diagnostics,
                 search.ranges = if (object$settings$method == "conditional")
                   .hc_search_ranges(object) else NULL,
                 covariance = .hc_covariance_label(object$fit)),
            class = "summary_cohort_sensitivity")
}

#' @export
print.summary_cohort_sensitivity <- function(x, ...) {
  balanced <- isTRUE(x$settings$balanced)
  if (!is.null(x$covariance)) cat("Covariance:", x$covariance, if (balanced) "\n" else "\n\n")
  if (balanced) cat(.hc_balanced_label(x$settings), "\n\n", sep = "")
  cat("Conventional estimates\n"); print(x$estimates, row.names = FALSE)
  cat("\nSensitivity confidence sets\n"); print(x$intervals, row.names = FALSE)
  if (x$settings$method == "conditional") .hc_print_search_ranges(x$search.ranges)
  cat("\nSettings\n"); print(x$settings)
  cat("\nDiagnostics\n"); print(x$diagnostics)
  invisible(x)
}
