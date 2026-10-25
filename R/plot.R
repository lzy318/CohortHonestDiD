#' Plot event studies and sensitivity confidence sets
#'
#' ATT results are plotted against the sensitivity bound. Dynamic and
#' cohort-dynamic results overlay robust sets on conventional event-study
#' coefficients. Multiple bounds are faceted. All intervals are pointwise.
#' Disconnected confidence-set components are drawn separately; unbounded
#' endpoints use arrows. Numerical failures are reported, not drawn as sets.
#'
#' @param x A `cohort_fit` or `cohort_sensitivity` object.
#' @param type For fits, an aggregated (`"dynamic"`) or cohort-specific
#'   (`"cohort-dynamic"`) event study.
#' @param M Optional subset of computed bounds to plot.
#' @param cohorts Optional subset of integer cohort adoption indices.
#' @param identified.set Show the plug-in identified-set band for ATT plots?
#' @param main Optional plot title.
#' @param xlab,ylab Axis labels.
#' @param color Optional sensitivity-set color, overriding the paper palette.
#' @param ... Reserved for additional graphical arguments.
#' @return An editable `ggplot` object.
#' @name plot_cohort
#' @importFrom rlang .data
NULL

.hc_cohort_labels <- function(fit, cohorts) {
  labels <- as.character(cohorts)
  time_labels <- fit$data_info$time_labels
  if (!is.null(time_labels)) {
    valid <- !is.na(cohorts) & cohorts >= 1L & cohorts <= length(time_labels)
    labels[valid] <- as.character(time_labels[cohorts[valid]])
  }
  paste0("Cohort ", labels)
}

.hc_event_data <- function(fit, by_cohort, conf.level = 0.95) {
  info <- fit$coordinate_info
  crit <- stats::qnorm((1 + conf.level) / 2)
  if (by_cohort) {
    out <- info
    out$estimate <- as.numeric(fit$betahat)
    out$se <- sqrt(pmax(0, diag(fit$Sigma)))
  } else {
    periods <- sort(unique(info$relative_period))
    masses <- fit$cohort_masses
    if (is.null(masses)) masses <- fit$cohort_sizes
    if (is.null(names(masses))) names(masses) <- sort(unique(info$cohort))
    out <- do.call(rbind, lapply(periods, function(s) {
      idx <- which(info$relative_period == s)
      w <- setNames(numeric(nrow(info)), info$coordinate)
      w[idx] <- masses[as.character(info$cohort[idx])]
      w <- w / sum(w)
      data.frame(cohort = NA_integer_, relative_period = s,
                 estimate = sum(w * fit$betahat),
                 se = sqrt(max(0, as.numeric(t(w) %*% fit$Sigma %*% w))))
    }))
  }
  out$lower <- out$estimate - crit * out$se
  out$upper <- out$estimate + crit * out$se
  out$cohort_label <- if (by_cohort) .hc_cohort_labels(fit, out$cohort) else "All cohorts"
  out
}

#' @rdname plot_cohort
#' @export
plot.cohort_fit <- function(x, type = c("dynamic", "cohort-dynamic"), cohorts = NULL,
                            main = NULL, xlab = NULL, ylab = "Effect", ...) {
  type <- match.arg(type)
  d <- .hc_event_data(x, type == "cohort-dynamic")
  if (!is.null(cohorts)) d <- d[d$cohort %in% cohorts, , drop = FALSE]
  if (!nrow(d)) stop("No plotted coefficients match the selection.", call. = FALSE)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$relative_period, y = .data$estimate)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0.5, color = "grey65") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), linewidth = 0.4) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::labs(title = main, x = if (is.null(xlab)) "Relative period (1 = first treated)" else xlab,
                  y = ylab, caption = "Conventional estimates with 95% pointwise intervals.") +
    ggplot2::theme_bw(base_size = 11)
  if (type == "cohort-dynamic") p <- p + ggplot2::facet_wrap(~cohort_label)
  p
}

.hc_interval_layers <- function(p, d, xvar, color, offset = 0) {
  if (!nrow(d)) return(p)
  d$xvalue <- d[[xvar]] + offset
  finite_values <- c(d$lower[is.finite(d$lower)], d$upper[is.finite(d$upper)], 0)
  rg <- range(finite_values)
  pad <- max(diff(rg) * 0.15, max(abs(rg)) * 0.05, 0.01)
  d$lower_plot <- ifelse(is.infinite(d$lower) & d$lower < 0, rg[1L] - pad, d$lower)
  d$upper_plot <- ifelse(is.infinite(d$upper) & d$upper > 0, rg[2L] + pad, d$upper)
  finite <- is.finite(d$lower) & is.finite(d$upper)
  if (any(finite)) p <- p + ggplot2::geom_linerange(data = d[finite, , drop = FALSE],
    ggplot2::aes(x = .data$xvalue, ymin = .data$lower_plot, ymax = .data$upper_plot),
    inherit.aes = FALSE, color = color, linewidth = 0.65)
  for (which_end in c("first", "last", "both")) {
    sel <- switch(which_end, first = is.infinite(d$lower) & is.finite(d$upper),
                  last = is.finite(d$lower) & is.infinite(d$upper),
                  both = is.infinite(d$lower) & is.infinite(d$upper))
    if (any(sel)) p <- p + ggplot2::geom_segment(data = d[sel, , drop = FALSE],
      ggplot2::aes(x = .data$xvalue, xend = .data$xvalue, y = .data$lower_plot, yend = .data$upper_plot),
      inherit.aes = FALSE, color = color, linewidth = 0.65,
      arrow = grid::arrow(ends = which_end, length = grid::unit(0.08, "inches")))
  }
  p
}

#' @rdname plot_cohort
#' @export
plot.cohort_sensitivity <- function(x, M = NULL, cohorts = NULL, identified.set = TRUE,
                                    main = NULL, xlab = NULL, ylab = "Effect", color = NULL, ...) {
  s <- x$settings
  d <- x$intervals
  if (!is.null(M)) {
    if (!is.numeric(M) || any(!is.finite(M))) stop("M must select finite computed bounds.", call. = FALSE)
    keep <- vapply(d$M, function(v) any(abs(v - M) <= 1e-10 * pmax(1, abs(v), abs(M))), logical(1))
    d <- d[keep, , drop = FALSE]
  }
  if (!is.null(cohorts)) d <- d[d$cohort %in% cohorts, , drop = FALSE]
  if (!nrow(d)) stop("No computed results match the requested plot selection.", call. = FALSE)
  selected.bounds <- sort(unique(d$M))
  if (is.null(color)) color <- unname(paper_palette("representation")[s$framework])
  failures <- unique(d[!d$success, c("target", "M", "status")])
  empty <- unique(d[d$success & is.na(d$lower) & is.na(d$upper), c("target", "M")])
  scope <- if (s$restriction == "RM") "RM" else paste("SD,", if (s$full.path) "full path" else "post-only")
  method_label <- switch(s$method, flci = "FLCI", conditional = "Conditional", c_lf = "C-LF")
  caption <- sprintf("%g%% pointwise sets; %s; %s; %s inference.", 100 * s$conf.level, method_label, scope, s$framework)
  if (s$method == "conditional") caption <- paste(caption,
    "Bounded numerical search; components outside the searched range are not ruled out. See summary() for searched ranges.")
  if (nrow(failures)) caption <- paste(caption, nrow(failures), "failed or unresolved result(s); inspect summary().")
  if (nrow(empty)) caption <- paste(caption, nrow(empty), "empty set(s).")
  if (isFALSE(x$fit$data_info$balanced)) caption <- paste(caption, "Unbalanced-panel approximation.")
  if (isTRUE(s$balanced)) caption <- paste(caption, paste0("Balanced window: relative periods ",
    s$window[["lower"]], " to ", s$window[["upper"]],
    if (identical(s$framework, "cohort")) paste0("; cohort paths: ", s$cohort.paths) else "", "."))
  d <- d[d$success & !is.na(d$lower) & !is.na(d$upper), , drop = FALSE]
  shape <- if (s$method == "flci") 17 else if (s$method == "conditional") 15 else 16
  if (s$target == "att") {
    p <- ggplot2::ggplot(data.frame(M = selected.bounds), ggplot2::aes(x = .data$M)) +
      ggplot2::geom_blank(ggplot2::aes(y = 0)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
      ggplot2::geom_hline(yintercept = x$estimates$estimate[1L], color = "grey35", linetype = "dotted")
    if (identified.set && nrow(d)) {
      bands <- d[!duplicated(d$M) & is.finite(d$plugin.lower) & is.finite(d$plugin.upper), , drop = FALSE]
      if (nrow(bands)) p <- p + ggplot2::geom_ribbon(data = bands,
        ggplot2::aes(x = .data$M, ymin = .data$plugin.lower, ymax = .data$plugin.upper), fill = color, alpha = 0.2)
    }
    p <- .hc_interval_layers(p, d, "M", color)
    if (nrow(d)) {
      counts <- table(d$M)
      centers <- d[as.numeric(counts[as.character(d$M)]) == 1L & is.finite(d$center), , drop = FALSE]
      if (nrow(centers)) p <- p + ggplot2::geom_point(data = centers,
        ggplot2::aes(x = .data$M, y = .data$center), color = color, shape = shape, size = 2)
    }
    x_label <- if (s$restriction == "RM") expression(bar(M)) else "M"
  } else {
    event <- .hc_event_data(x$fit, s$target == "cohort-dynamic", s$conf.level)
    if (isTRUE(s$balanced)) event <- event[event$relative_period >= s$window[["lower"]] &
                                             event$relative_period <= s$window[["upper"]], , drop = FALSE]
    if (!is.null(cohorts)) event <- event[event$cohort %in% cohorts, , drop = FALSE]
    bounds <- if (is.null(M)) s$M else s$M[vapply(s$M, function(v) any(abs(v - M) < 1e-10 * pmax(1, abs(v), abs(M))), logical(1))]
    event <- do.call(rbind, lapply(bounds, function(m) { z <- event; z$M <- m; z }))
    event$panel <- paste0(if (s$target == "cohort-dynamic") paste0(event$cohort_label, " | ") else "", "M = ", format(event$M, trim = TRUE))
    d$panel <- if (nrow(d)) paste0(if (s$target == "cohort-dynamic") paste0(.hc_cohort_labels(x$fit, d$cohort), " | ") else "", "M = ", format(d$M, trim = TRUE)) else character()
    p <- ggplot2::ggplot(event, ggplot2::aes(x = .data$relative_period, y = .data$estimate)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey55", linetype = "dashed") +
      ggplot2::geom_vline(xintercept = 0.5, color = "grey65") +
      ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), color = "black", linewidth = 0.35) +
      ggplot2::geom_point(size = 1.5, color = "black")
    p <- .hc_interval_layers(p, d, "relative_period", color, offset = 0.12)
    if (nrow(d)) {
      key <- paste(d$target, d$M)
      counts <- table(key)
      centers <- d[as.numeric(counts[key]) == 1L & is.finite(d$center), , drop = FALSE]
      if (nrow(centers)) p <- p + ggplot2::geom_point(data = centers,
        ggplot2::aes(x = .data$relative_period + 0.12, y = .data$center), inherit.aes = FALSE,
        color = color, shape = shape, size = 2)
    }
    p <- p + ggplot2::facet_wrap(~panel)
    caption <- paste(caption, "Black: conventional estimates. Colored: sensitivity sets.")
    x_label <- "Relative period (1 = first treated)"
  }
  p + ggplot2::labs(title = main, x = if (is.null(xlab)) x_label else xlab, y = ylab,
                    caption = paste(strwrap(caption, width = 95), collapse = "\n")) +
    ggplot2::theme_bw(base_size = 11)
}
