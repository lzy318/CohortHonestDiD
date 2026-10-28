#' Plot event studies and sensitivity confidence sets
#'
#' ATT results are plotted against the sensitivity bound. Dynamic and
#' cohort-dynamic results overlay robust sets on conventional event-study
#' coefficients. Multiple bounds are faceted. All intervals are pointwise.
#' Disconnected confidence-set components are drawn separately; unbounded
#' endpoints use arrows. Numerical failures are reported, not drawn as sets.
#'
#' @param x A \code{cohort_fit} or \code{cohort_sensitivity} object.
#' @param type For fits, an aggregated (\code{"dynamic"}) or cohort-specific
#'   (\code{"cohort-dynamic"}) event study. For sensitivity results,
#'   \code{"sets"} (default) plots the confidence sets and \code{"weights"}
#'   the FLCI weights.
#' @param M Optional subset of computed bounds to plot.
#' @param cohorts Optional subset of integer cohort adoption indices. With
#'   \code{type = "weights"}, it selects the rows of these cohorts.
#' @param identified.set Show the plug-in identified-set band for ATT plots?
#'   Only for \code{type = "sets"}.
#' @param main Optional plot title.
#' @param xlab,ylab Axis labels.
#' @param color Optional sensitivity-set color, overriding the paper palette.
#'   Only for \code{type = "sets"}.
#' @param weights With \code{type = "weights"}: \code{"coefficients"}
#'   (default) or \code{"gamma"}, as in \code{\link{flci_weights}}.
#' @param by With \code{type = "weights"}: \code{"element"} (default) or
#'   \code{"cohort"}, as in \code{\link{flci_weights}}.
#' @param targets With \code{type = "weights"}: optional target ids, as in
#'   \code{x$intervals$target}; by default every target, one panel each.
#' @param labels With \code{type = "weights"}: print the values in the cells?
#' @param ... Reserved for additional graphical arguments.
#' @return An editable \code{ggplot} object. With \code{type = "weights"},
#'   \code{NULL} (invisibly) when there is nothing to draw.
#' @details With \code{type = "weights"}, an FLCI result is drawn as a heatmap
#'   of \code{flci_weights(x, weights, by)}: one row per pre-treatment
#'   coefficient or second difference (labelled by relative period, grouped by
#'   cohort, earliest at the top) or per cohort, one column per computed bound,
#'   and one panel per target. Negative values are blue, positive values
#'   vermilion, and zeros white. \code{M}, \code{cohorts} and \code{targets}
#'   select bounds, rows and panels; \code{identified.set} and \code{color}
#'   apply only to \code{type = "sets"}. For a post-only FLCI, \code{gamma} is
#'   zero by construction, and \code{weights = "gamma"} only prints a message.
#'   Conditional and C-LF results have no such weights.
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
                                    main = NULL, xlab = NULL, ylab = "Effect", color = NULL,
                                    type = c("sets", "weights"),
                                    weights = c("coefficients", "gamma"),
                                    by = c("element", "cohort"), targets = NULL,
                                    labels = FALSE, ...) {
  type <- match.arg(type)
  if (type == "weights") {
    if (!missing(identified.set) || !missing(color))
      stop("identified.set and color apply only to type = 'sets'.", call. = FALSE)
    return(.hc_plot_weights(x, M = M, cohorts = cohorts, main = main, xlab = xlab,
                            ylab = if (missing(ylab)) NULL else ylab, weights = weights,
                            by = by, targets = targets, labels = labels))
  }
  if (!missing(weights) || !missing(by) || !missing(targets) || !missing(labels))
    stop("weights, by, targets, and labels apply only to type = 'weights'.", call. = FALSE)
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

# FLCI weights as a heatmap: rows are elements or cohorts, columns the bounds,
# one panel per target.
#' @noRd
.hc_plot_weights <- function(x, M, cohorts, main, xlab, ylab, weights, by, targets, labels) {
  weights <- match.arg(weights, c("coefficients", "gamma"))
  by <- match.arg(by, c("element", "cohort"))
  if (!is.logical(labels) || length(labels) != 1L || is.na(labels))
    stop("labels must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(targets) && (!is.character(targets) || !length(targets) || anyNA(targets)))
    stop("targets must be NULL or a character vector of target ids.", call. = FALSE)
  s <- x$settings
  .hc_flci_weights_method(s)
  aggregate <- identical(s$framework, "aggregate")
  if (by == "cohort" && aggregate) .hc_flci_weights_no_cohorts()
  if (!is.null(cohorts) && aggregate)
    stop("cohorts cannot be selected under aggregate inference: aggregate weights have no cohorts.",
         call. = FALSE)
  if (weights == "gamma" && !isTRUE(s$full.path)) {
    message("The post-only FLCI puts no weight on the pure pre-treatment second differences: ",
            "its weights gamma are zero by construction. Use weights = 'coefficients' to plot ",
            "its weights v.")
    return(invisible(NULL))
  }
  if (!is.null(M) && (!is.numeric(M) || any(!is.finite(M))))
    stop("M must select finite computed bounds.", call. = FALSE)
  ids <- vapply(x$targets, function(z) as.character(z$id), character(1))
  task.target <- rep(ids, each = length(s$M))
  task.M <- rep(s$M, times = length(ids))
  keep <- if (is.null(targets)) rep(TRUE, length(task.target)) else task.target %in% targets
  if (!is.null(M)) keep <- keep & vapply(task.M, function(v)
    any(abs(v - M) <= 1e-10 * pmax(1, abs(v), abs(M))), logical(1))
  selected <- which(keep)
  if (!length(selected)) stop("No computed results match the requested plot selection.", call. = FALSE)
  has.weights <- vapply(selected, function(k) !is.null(x$results[[k]]$weights), logical(1))
  if (!any(has.weights))
    stop("No FLCI weights to plot: every selected result failed; inspect summary().", call. = FALSE)
  d <- flci_weights(x, weights, by)
  in.selection <- logical(nrow(d))
  for (k in selected) in.selection <- in.selection | (d$target == task.target[k] & d$M == task.M[k])
  d <- d[in.selection, , drop = FALSE]
  if (!nrow(d)) {
    if (weights == "gamma") {
      message("This result has no pure pre-treatment second differences, so the FLCI has no ",
              "weights gamma to plot.")
      return(invisible(NULL))
    }
    stop("No computed results match the requested plot selection.", call. = FALSE)
  }
  if (!is.null(cohorts)) {
    d <- d[d$cohort %in% cohorts, , drop = FALSE]
    if (!nrow(d)) stop("No computed results match the requested plot selection.", call. = FALSE)
  }
  rownames(d) <- NULL
  d$value <- if (by == "element") d$weight else d$contribution
  bounds <- sort(unique(d$M))
  columns <- format(bounds, trim = TRUE, drop0trailing = TRUE)
  if (anyDuplicated(columns)) columns <- format(bounds, digits = 15, trim = TRUE, drop0trailing = TRUE)
  d$column <- factor(columns[match(d$M, bounds)], levels = columns)
  if (by == "element") {
    periods <- sort(unique(d$relative_period), decreasing = TRUE)
    d$row <- factor(as.character(d$relative_period), levels = as.character(periods))
  } else {
    d$row <- factor(.hc_cohort_labels(x$fit, d$cohort),
                    levels = .hc_cohort_labels(x$fit, sort(unique(d$cohort), decreasing = TRUE)))
  }
  if (!aggregate) d$cohort_label <- factor(.hc_cohort_labels(x$fit, d$cohort),
                                           levels = .hc_cohort_labels(x$fit, sort(unique(d$cohort))))
  target.labels <- vapply(x$targets, function(z) as.character(z$label), character(1))
  d$target_label <- factor(d$label, levels = unique(target.labels[target.labels %in% d$label]))
  d$text <- ifelse(d$value == 0, "", trimws(formatC(d$value, format = "fg", digits = 2)))
  lim <- max(abs(d$value))
  if (!is.finite(lim) || lim == 0) lim <- 1
  d$dark <- abs(d$value) > lim / 2
  legend <- if (by == "cohort") "Contribution" else if (weights == "coefficients") "Weight v" else
    "Weight gamma"
  y_label <- if (by == "cohort") "Cohort" else if (weights == "coefficients") "Relative period" else
    "Last period of the second difference"
  what <- if (by == "element") {
    if (weights == "coefficients") "FLCI weights v on the pre-treatment coefficients" else
      "FLCI weights gamma on the pure pre-treatment second differences"
  } else if (weights == "coefficients") {
    "Sum over each cohort's pre-treatment coefficients of v times the estimate"
  } else "Sum over each cohort's second differences of gamma times the estimate"
  scope <- if (isTRUE(s$full.path)) "SD, full path" else "SD, post-only"
  caption <- sprintf("%s; %s; %s inference.", what, scope, s$framework)
  caption <- paste(caption, if (by == "element") "White cells are zero weights." else
    paste("White cells are zero. The cohort sums add up to the FLCI center minus",
          if (weights == "coefficients") "the conventional estimate." else "the post-only center."))
  tied <- sort(unique(d$M[d$tie.break]))
  if (length(tied)) caption <- paste(caption, paste0("At M = ",
    paste(format(tied, trim = TRUE, drop0trailing = TRUE), collapse = ", "),
    ", several weight vectors give the same FLCI length; the one with the smallest Euclidean norm",
    " of gamma is shown."))
  if (weights == "gamma" && all(d$value == 0)) caption <- paste(caption,
    "All weights gamma are zero at the selected bounds: the full-path FLCI equals the post-only FLCI.")
  failed <- sum(!has.weights)
  if (failed) caption <- paste(caption, failed,
                               "failed result(s) without weights are not shown; inspect summary().")
  if (isFALSE(x$fit$data_info$balanced)) caption <- paste(caption, "Unbalanced-panel approximation.")
  if (isTRUE(s$balanced)) caption <- paste(caption, paste0("Balanced window: relative periods ",
    s$window[["lower"]], " to ", s$window[["upper"]],
    if (identical(s$framework, "cohort")) paste0("; cohort paths: ", s$cohort.paths) else "", "."))
  colors <- paper_palette("representation")
  p <- ggplot2::ggplot(d) +
    ggplot2::geom_tile(ggplot2::aes(x = .data$column, y = .data$row, fill = .data$value),
                       color = "grey88", linewidth = 0.15)
  if (labels) p <- p +
    ggplot2::geom_text(ggplot2::aes(x = .data$column, y = .data$row, label = .data$text,
                                    color = .data$dark), size = 2.5) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none")
  p <- p +
    ggplot2::scale_fill_gradient2(low = colors[["cohort"]], mid = "white",
                                  high = colors[["aggregate"]], midpoint = 0,
                                  limits = c(-lim, lim), name = legend) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0))
  p <- p + if (!aggregate && by == "element") {
    ggplot2::facet_grid(cohort_label ~ target_label, scales = "free_y", space = "free_y",
                        switch = "y")
  } else ggplot2::facet_grid(. ~ target_label)
  p + ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   strip.placement = "outside",
                   strip.text.y.left = ggplot2::element_text(angle = 0),
                   strip.background = ggplot2::element_blank(),
                   panel.spacing = grid::unit(2, "pt")) +
    ggplot2::labs(title = main, x = if (is.null(xlab)) "M" else xlab,
                  y = if (is.null(ylab)) y_label else ylab, fill = legend,
                  caption = paste(strwrap(caption, width = 95), collapse = "\n"))
}
