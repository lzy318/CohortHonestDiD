# Conditional inversion has no certified outer support. These finite searches
# diagnose stability over wider domains without changing the point test.
.hc_conditional_note <- function() {
  "Bounded numerical search; components outside the searched range are not ruled out."
}

.hc_conditional_diagnostics <- function() list(
  support.certified = FALSE, topology.certified = FALSE, inversion.exact = FALSE,
  search.note = .hc_conditional_note(), search.center = NA_real_,
  search.scale = NA_real_, search.scale.source = NA_character_,
  initial.domain = numeric(), search.domain = numeric(),
  requested.domain = numeric(), domain.stable = FALSE,
  initial.grid.points = NA_real_, grid.density.adjusted = NA,
  reference.domain.radius = NA_real_, reference.grid.spacing = NA_real_,
  domain.expansions = 0L, agreement.streak = 0L, domain.history = data.frame(),
  search.attempts = data.frame(), refinements = 0L, test.evaluations = 0L,
  endpoint.certified = FALSE, resolution.stable = NA)

.hc_conditional_merge_diagnostics <- function(diagnostics) {
  answer <- .hc_conditional_diagnostics()
  answer[names(diagnostics)] <- diagnostics
  answer
}

# Normalize only the active quadratic block, so unrelated outcome scales do
# not erase the search contrast. The roundoff coefficient is the existing
# quadratic-form coefficient, not an absolute outcome-unit floor.
.hc_conditional_sd <- function(v, Sigma) {
  active <- which(v != 0)
  if (!length(active)) return(0)
  v <- v[active]
  Sigma <- Sigma[active, active, drop = FALSE]
  vs <- max(abs(v)); ss <- max(abs(Sigma))
  if (ss == 0) return(0)
  vn <- v / vs; Sn <- Sigma / ss
  if (any(v != 0 & vn == 0) || any(Sigma != 0 & Sn == 0)) return(0)
  variance <- sum(vn * as.numeric(Sn %*% vn))
  bound <- sum(abs(vn) * as.numeric(abs(Sn) %*% abs(vn)))
  answer <- .hc_backend()$flci_safe_sqrt_variance(
    variance, bound, rel_tol = 1e-10, label = "Conditional search variance")
  if (variance <= 1e-10 * bound) return(0)
  answer * sqrt(ss) * vs
}

.hc_conditional_domain <- function(representation, built, plugin, restriction, control) {
  info <- representation$coordinate_info
  K <- nrow(info)
  H <- diag(K)
  differences <- list()
  for (g in sort(unique(info$cohort))) {
    rows <- which(info$cohort == g)
    rows <- rows[order(info$relative_period[rows])]
    pre <- rows[info$type[rows] == "pre"]
    post <- rows[info$type[rows] == "post"]
    last <- pre[info$relative_period[pre] == 0]
    previous <- pre[info$relative_period[pre] == -1]
    if (length(last) != 1L || (restriction == "SD" && length(previous) != 1L))
      stop("Conditional search continuation requires the last pre-period coordinates.", call. = FALSE)
    H[post, ] <- 0
    H[post, last] <- if (restriction == "SD") 1 + info$relative_period[post] else 1
    if (restriction == "SD") H[post, previous] <- -info$relative_period[post]
    for (j in seq_len(length(rows) - 1L)) {
      q <- numeric(K)
      q[rows[c(j, j + 1L)]] <- c(-1, 1)
      differences[[length(differences) + 1L]] <- q
    }
  }
  W.inverse <- solve(representation$W)
  ell <- representation$target
  v <- as.numeric(ell - crossprod(representation$W %*% H %*% W.inverse, ell))
  center <- sum(v * representation$betahat)
  scale <- .hc_conditional_sd(v, representation$Sigma)
  source <- "target contrast"
  if (!is.finite(scale) || scale <= 0) {
    Q <- do.call(rbind, differences) %*% W.inverse
    scales <- apply(Q, 1L, .hc_conditional_sd, Sigma = representation$Sigma)
    scale <- max(scales)
    source <- "first-difference fallback"
  }
  B <- if (restriction == "SD") built$model$B_Q else
    max(abs(c(plugin$lower, plugin$upper) - center))
  radius <- B + control$grid.radius * scale
  domain <- center + c(-radius, radius)
  if (!is.finite(center) || !is.finite(scale) || scale <= 0 ||
      !is.finite(B) || B < 0 || !is.finite(radius) || radius <= 0 ||
      any(!is.finite(domain)) || domain[1L] >= domain[2L])
    stop("Conditional search has no finite positive numerical search scale/domain after the invariant first-difference fallback.", call. = FALSE)
  se <- sqrt(as.numeric(crossprod(built$target_bar, built$Sigma_bar %*% built$target_bar)))
  tolerance <- if (is.null(control$endpoint.tol)) 1e-5 * se else control$endpoint.tol
  if (!is.finite(se) || se <= 0 || !is.finite(tolerance) || tolerance <= 0)
    stop("Conditional endpoint tolerance requires a finite positive conventional target standard error.", call. = FALSE)
  list(center = center, scale = scale, source = source, radius = radius, bias.radius = B,
       domain = domain, tolerance = tolerance, se = se)
}

# Called only for a wider-than-default initial radius. A radius change must
# not make the requested reference density coarser. The ceiling preserves a
# spacing bound; separately requested radii need not share a lattice phase.
.hc_conditional_density <- function(initial, control) {
  n <- as.double(control$grid.points) - 1
  integer.limit <- as.double(.Machine$integer.max) - 1
  reference.radius <- initial$bias.radius + 8 * initial$scale
  # Divide first so a representable reference spacing does not require an
  # overflowing doubled radius or n*radius intermediate.
  reference.spacing <- reference.radius / (n / 2)
  details <- list(initial.grid.points = NA_real_, grid.density.adjusted = NA,
    reference.domain.radius = reference.radius,
    reference.grid.spacing = reference.spacing)
  failed <- function(status, message) list(success = FALSE, status = status,
    message = message, diagnostics = details)
  if (!is.finite(reference.radius) || reference.radius <= 0 ||
      !is.finite(reference.spacing) || reference.spacing <= 0)
    return(failed("search_scale_failure",
      "The wider Conditional search has no finite positive reference grid spacing."))
  ratio <- initial$radius / reference.radius
  if (!is.finite(ratio) || ratio > integer.limit / n)
    return(failed("evaluation_limit",
      "Preserving Conditional reference density requires an initial grid beyond R's supported integer range."))
  intervals <- max(n, ceiling(n * ratio))
  if (!is.finite(intervals) || intervals > integer.limit)
    return(failed("evaluation_limit",
      "Preserving Conditional reference density requires an initial grid beyond R's supported integer range."))
  details$initial.grid.points <- intervals + 1
  details$grid.density.adjusted <- intervals > n
  if (intervals + 1 > control$max.evaluations)
    return(failed("evaluation_limit", paste0("Preserving Conditional reference density requires ",
      format(intervals + 1, scientific = FALSE), " initial regular points, exceeding max.evaluations=" ,
      control$max.evaluations, ".")))
  width <- initial$domain[2L] - initial$domain[1L]
  if (!is.finite(width) || width <= 0)
    return(failed("search_scale_failure",
      "The wider Conditional domain has no finite positive represented grid width."))
  if (width / intervals > reference.spacing) {
    # Endpoint rounding can make the represented width slightly larger than
    # 2*radius. Calculate the necessary increase rather than looping over a
    # potentially large correction, then enforce the bound without a tolerance.
    if (width / integer.limit > reference.spacing) {
      details$initial.grid.points <- NA_real_; details$grid.density.adjusted <- NA
      return(failed("evaluation_limit",
        "The represented Conditional domain requires too many intervals to preserve reference density."))
    }
    required <- ceiling(width / reference.spacing)
    if (!is.finite(required) || required > integer.limit) {
      details$initial.grid.points <- NA_real_; details$grid.density.adjusted <- NA
      return(failed("evaluation_limit",
        "The represented Conditional domain requires too many intervals to preserve reference density."))
    }
    intervals <- max(intervals, required)
    while (width / intervals > reference.spacing) {
      if (intervals >= integer.limit || intervals + 2 > control$max.evaluations) {
        details$initial.grid.points <- if (intervals >= integer.limit) NA_real_ else intervals + 2
        details$grid.density.adjusted <- if (is.finite(details$initial.grid.points)) TRUE else NA
        return(failed("evaluation_limit",
          "An upward grid-rounding adjustment would exceed the Conditional evaluation or integer limit."))
      }
      intervals <- intervals + 1
    }
    details$initial.grid.points <- intervals + 1
    details$grid.density.adjusted <- intervals > n
    if (intervals + 1 > control$max.evaluations)
      return(failed("evaluation_limit", paste0("Preserving the represented Conditional grid spacing requires ",
        format(intervals + 1, scientific = FALSE), " initial regular points, exceeding max.evaluations=",
        control$max.evaluations, ".")))
  }
  list(success = TRUE, points = intervals + 1, diagnostics = details)
}

.hc_conditional_agrees <- function(a, b, tolerance) {
  if (is.null(a) || !isTRUE(a$success) || !isTRUE(b$success) ||
      !identical(a$status, b$status) || !nrow(a$components) ||
      nrow(a$components) != nrow(b$components)) return(FALSE)
  x <- as.matrix(a$components[, c("lower", "upper"), drop = FALSE])
  y <- as.matrix(b$components[, c("lower", "upper"), drop = FALSE])
  finite <- is.finite(x) & is.finite(y)
  if (any(is.na(x)) || any(is.na(y)) ||
      any(is.finite(x) != is.finite(y)) ||
      any(sign(x[!finite]) != sign(y[!finite]))) return(FALSE)
  !any(abs(x[finite] - y[finite]) > 4 * tolerance)
}

.hc_conditional_search <- function(representation, built, plugin, restriction,
                                   control, anchors, problem, test, budget,
                                   diagnostics = list()) {
  engine <- .hc_backend()
  diagnostics <- .hc_conditional_merge_diagnostics(diagnostics)
  initial <- tryCatch(.hc_conditional_domain(representation, built, plugin,
    restriction, control), error = identity)
  if (inherits(initial, "error")) return(.hc_engine_failure(
    "search_scale_failure", conditionMessage(initial), diagnostics,
    c(plugin$lower, plugin$upper)))
  diagnostics$search.center <- initial$center
  diagnostics$search.scale <- initial$scale
  diagnostics$search.scale.source <- initial$source
  diagnostics$initial.domain <- initial$domain
  diagnostics$endpoint.tolerance <- initial$tolerance
  points <- control$grid.points
  if (control$grid.radius <= 8) {
    # Preserve the C4 numerical path, including tiny-radius cases where an
    # unused radius-eight reference could overflow.
    diagnostics$initial.grid.points <- points
    diagnostics$grid.density.adjusted <- FALSE
  } else if (control$max.expansions >= 2L) {
    density <- .hc_conditional_density(initial, control)
    diagnostics[names(density$diagnostics)] <- density$diagnostics
    if (!isTRUE(density$success)) {
      effective <- diagnostics$initial.grid.points
      diagnostics$requested.domain <- initial$domain
      diagnostics$final.grid.points <- effective
      diagnostics$search.attempts <- data.frame(attempt = 0L, grid.points = effective,
        new.evaluations = 0L, cumulative.evaluations = 0L,
        success = FALSE, status = density$status)
      diagnostics$domain.history <- data.frame(refinement = 0L, domain = 0L,
        lower = initial$domain[1L], upper = initial$domain[2L],
        regular.spacing = if (is.finite(effective))
          initial$radius / ((effective - 1) / 2) else NA_real_,
        regular.points = effective, new.evaluations = 0L, cumulative.evaluations = 0L,
        local.status = density$status, evaluated = FALSE,
        evaluated.lower = NA_real_, evaluated.upper = NA_real_, components = 0L,
        agreement = FALSE, agreement.streak = 0L, endpoint.certified = FALSE,
        resolution.stable = NA)
      return(.hc_engine_failure(density$status, density$message, diagnostics,
        c(plugin$lower, plugin$upper)))
    }
    points <- density$points
  }
  history <- list()
  attempts <- list()
  searched <- numeric()
  requested <- numeric()
  refinement <- 0L
  streak <- 0L
  widest.level <- 0L
  endpoint.ok <- FALSE
  resolution <- NA
  domain.visited <- FALSE
  domain.searched <- numeric()
  partial <- data.frame()
  # Track actual evaluated coordinates even if the private inversion catches
  # a budget exception midway through a domain.
  tracked.test <- function(theta) {
    cached <- exists(paste0("point:", sprintf("%.17g", theta)),
                     budget$state$cache, inherits = FALSE)
    if (cached || budget$state$evaluations < control$max.evaluations) {
      searched <<- range(c(searched, theta))
      domain.searched <<- range(c(domain.searched, theta))
      domain.visited <<- TRUE
    }
    test(theta)
  }
  failure <- function(status, message, components = data.frame())
    list(success = FALSE, status = status, message = message, components = components)
  limit <- function(message) failure("evaluation_limit", message)
  invert <- function(grid, spacing) {
    keys <- paste0("point:", sprintf("%.17g", grid))
    needed <- sum(!vapply(keys, exists, logical(1), envir = budget$state$cache,
                          inherits = FALSE))
    if (needed > control$max.evaluations - budget$state$evaluations)
      return(limit("The requested domain grid exceeds the remaining cumulative point-evaluation budget."))
    fit <- engine$sd_invert_test(tracked.test, grid,
      endpoint_tolerance = initial$tolerance, max_expansions = 0L,
      max_test_evaluations = control$max.evaluations, expansion_step = spacing)
    if (!isTRUE(fit$success) && length(fit$failed_results)) {
      reasons <- unique(vapply(fit$failed_results, function(x)
        paste(x$message, collapse = "; "), character(1)))
      reasons <- reasons[!is.na(reasons) & nzchar(reasons)]
      if (length(reasons)) fit$message <- paste(fit$message,
        "Point-test detail:", paste(reasons, collapse = " | "))
    }
    if (isTRUE(fit$success) && !nrow(fit$components)) {
      fit$success <- FALSE
      fit$status <- "uncertified_empty_search"
      fit$message <- "The finite Conditional search found no accepted point; this does not certify an empty set."
    }
    fit
  }
  local.search <- function(domain, intervals) {
    regular <- seq(domain[1L], domain[2L], length.out = as.integer(intervals + 1))
    grid <- sort(unique(c(regular, anchors[is.finite(anchors) &
      anchors >= domain[1L] & anchors <= domain[2L]])))
    spacing <- (domain[2L] - domain[1L]) / intervals
    if (!is.finite(spacing) || spacing <= 0 || any(diff(regular) <= 0))
      return(failure("search_scale_failure", "The registered regular Conditional mesh is not representable."))
    fit <- invert(grid, spacing)
    if (nrow(fit$components)) partial <<- fit$components
    if (!isTRUE(fit$success)) return(fit)
    endpoint <- engine$sd_arp_audit_endpoints(problem, fit,
      probe_distance = max(10 * initial$tolerance, 1e-4 * initial$se),
      test_fn = tracked.test)
    endpoint.ok <<- isTRUE(endpoint$certified)
    if (!isTRUE(endpoint$success)) {
      fit$success <- FALSE; fit$status <- endpoint$status; fit$message <- endpoint$message
      if (nrow(endpoint$checks) && any(!endpoint$checks$success &
          endpoint$checks$status != "endpoint_probe_too_close")) fit$status <- "endpoint_test_failure"
      return(fit)
    }
    if (isTRUE(control$topology.audit)) {
      # Check the allocation bound before constructing either audit grid.
      if (2 * as.double(length(grid)) - 1 > control$max.evaluations ||
          2 * as.double(length(grid)) - 1 > .Machine$integer.max)
        return(limit("The requested mesh audit exceeds the cumulative point-evaluation budget or R's integer range."))
      midpoints <- utils::head(grid, -1L) + diff(grid) / 2
      bounded.anchors <- anchors[is.finite(anchors) & anchors >= domain[1L] & anchors <= domain[2L]]
      shifted <- sort(unique(c(domain, midpoints, bounded.anchors)))
      fine <- sort(unique(c(grid, midpoints)))
      audits <- lapply(list(shifted, fine), invert, spacing = spacing / 2)
      resolution <<- all(vapply(audits, function(x)
        .hc_conditional_agrees(fit, x, initial$tolerance), logical(1)))
      if (isTRUE(resolution)) {
        # Domain doublings repeat the same interior lattice. Audit the
        # shifted finer lattice too, using registered cells independently
        # of the irregular subintervals introduced by profile anchors.
        resolution <<- NA
        phase.count <- 2 * as.double(intervals) + 2
        if (!is.finite(intervals) || intervals < 1 || intervals != floor(intervals) ||
            !is.finite(phase.count) || phase.count > .Machine$integer.max ||
            phase.count > control$max.evaluations)
          return(limit("The requested quarter-phase audit exceeds the cumulative point-evaluation budget or R's integer range."))
        cells <- seq.int(0, intervals - 1)
        quarter <- domain[1L] + (cells + 1 / 4) * spacing
        three.quarter <- domain[1L] + (cells + 3 / 4) * spacing
        phases <- as.vector(rbind(quarter, three.quarter))
        cell.left <- domain[1L] + cells * spacing
        cell.right <- domain[1L] + (cells + 1) * spacing
        if (any(!is.finite(phases)) || any(diff(phases) <= 0) ||
            any(quarter <= cell.left) || any(three.quarter >= cell.right) ||
            phases[1L] <= domain[1L] || utils::tail(phases, 1L) >= domain[2L])
          return(failure("search_scale_failure",
            "The registered Conditional quarter-phase audit is not representable."))
        quarter.grid <- sort(unique(c(domain, phases, bounded.anchors)))
        quarter.audit <- invert(quarter.grid, spacing / 2)
        audits[[length(audits) + 1L]] <- quarter.audit
        resolution <<- .hc_conditional_agrees(fit, quarter.audit, initial$tolerance)
      }
      if (!resolution) {
        failed <- Filter(function(x) !isTRUE(x$success), audits)
        fit$success <- FALSE
        fit$status <- if (length(failed)) failed[[1L]]$status else "topology_audit_failure"
        fit$message <- if (length(failed)) failed[[1L]]$message else
          "Conditional base, shifted, finer, and quarter-phase meshes did not agree."
      }
    }
    fit
  }
  eligible <- c("open_search_boundary", "topology_audit_failure",
                "endpoint_audit_failure", "uncertified_empty_search", "topology_miss")
  repeat {
    before.attempt <- budget$state$evaluations
    previous <- NULL
    streak <- 0L
    fit <- failure("domain_stability_limit",
      "Conditional search requires two successive adjacent-domain agreements; max.expansions must be at least 2.")
    # This upper bound includes the first grid-allocation failure; avoid
    # allocating a sequence of a user-supplied, potentially huge domain cap.
    allocation.level <- max(0, ceiling(log2(control$max.evaluations /
      (as.double(points) - 1))) + 1)
    last.level <- min(control$max.expansions, allocation.level)
    if (control$max.expansions >= 2L) for (k in seq.int(0L, last.level)) {
      radius <- initial$radius * 2^k
      domain <- initial$center + c(-radius, radius)
      intervals <- (as.double(points) - 1) * 2^k
      before.domain <- budget$state$evaluations
      endpoint.ok <- FALSE; resolution <- NA
      domain.visited <- FALSE
      domain.searched <- numeric()
      partial <- data.frame()
      if (!is.finite(radius) || any(!is.finite(domain)) ||
          domain[1L] >= domain[2L]) {
        fit <- failure("search_scale_failure", "A requested wider Conditional domain is not finite and representable.")
      } else if (!is.finite(intervals) || intervals + 1 > .Machine$integer.max ||
                 intervals + 1 > control$max.evaluations) {
        fit <- limit("The requested wider Conditional grid exceeds the cumulative point-evaluation budget or R's integer range.")
      } else {
        fit <- tryCatch(local.search(domain, intervals),
          hc_search_budget_exhausted = function(e) limit(conditionMessage(e)))
      }
      if (domain.visited) widest.level <- max(widest.level, k)
      if (budget$state$exhausted) {
        if (nrow(fit$components)) partial <- fit$components
        fit <- limit(paste0("The cumulative limit of ", control$max.evaluations,
          " confidence-set point evaluations was reached."))
      }
      if (identical(fit$status, "evaluation_limit") && !nrow(fit$components)) {
        fit$components <- if (nrow(partial)) partial else if (!is.null(previous))
          previous$components else data.frame()
      }
      if (identical(fit$status, "evaluation_limit")) requested <- domain
      agreement <- isTRUE(fit$success) && !isTRUE(fit$boundary_accepted) &&
        .hc_conditional_agrees(previous, fit, initial$tolerance)
      streak <- if (agreement) streak + 1L else 0L
      if (isTRUE(fit$success) && streak >= 2L && k >= widest.level) {
        cached <- as.list(budget$state$cache)
        accepted <- vapply(cached, function(x) isTRUE(x$success) &&
          isTRUE(x$accept), logical(1))
        known <- as.numeric(sub("^point:", "", names(cached)[accepted]))
        covered <- vapply(known, function(theta) any(
          theta >= fit$components$lower - initial$tolerance &
          theta <= fit$components$upper + initial$tolerance), logical(1))
        rejected <- vapply(cached, function(x) isTRUE(x$success) &&
          identical(x$accept, FALSE), logical(1))
        known.rejections <- as.numeric(sub("^point:", "", names(cached)[rejected]))
        # A narrower rejection gap seen on an earlier mesh cannot disappear
        # inside a later component. Exclude only the existing endpoint error
        # bands, where the bisection approximation leaves membership unresolved.
        enclosed <- vapply(known.rejections, function(theta) any(
          theta > fit$components$lower + initial$tolerance &
          theta < fit$components$upper - initial$tolerance), logical(1))
        if (any(!covered)) {
          fit$success <- FALSE; fit$status <- "topology_miss"
          fit$message <- "A previously evaluated accepted point lies outside the recovered Conditional components."
          agreement <- FALSE; streak <- 0L
        } else if (any(enclosed)) {
          fit$success <- FALSE; fit$status <- "topology_miss"
          fit$message <- "A previously evaluated rejected point lies strictly inside a recovered Conditional component."
          agreement <- FALSE; streak <- 0L
        }
      }
      history[[length(history) + 1L]] <- data.frame(refinement = refinement,
        domain = k, lower = domain[1L], upper = domain[2L],
        regular.spacing = 2 * radius / intervals, regular.points = intervals + 1,
        new.evaluations = budget$state$evaluations - before.domain,
        cumulative.evaluations = budget$state$evaluations, local.status = fit$status,
        evaluated = domain.visited,
        evaluated.lower = if (length(domain.searched)) domain.searched[1L] else NA_real_,
        evaluated.upper = if (length(domain.searched)) domain.searched[2L] else NA_real_,
        components = nrow(fit$components), agreement = agreement,
        agreement.streak = streak, endpoint.certified = endpoint.ok,
        resolution.stable = resolution)
      # A finer retry must revisit every previously reached domain. It may
      # not forget a remote component observed on an earlier, coarser pass.
      if (isTRUE(fit$success) && streak >= 2L && k >= widest.level) break
      if (isTRUE(fit$success)) {
        previous <- fit
        if (k == control$max.expansions) fit <- failure("domain_stability_limit",
          "The maximum Conditional domain was reached before two successive adjacent-domain agreements.", fit$components)
      } else {
        previous <- NULL
        # Open or empty domains may become informative on the next wider
        # domain. Other eligible local uncertainty restarts at a finer mesh.
        if (!fit$status %in% c("open_search_boundary", "uncertified_empty_search")) break
      }
    }
    attempts[[length(attempts) + 1L]] <- data.frame(attempt = refinement,
      grid.points = points, new.evaluations = budget$state$evaluations - before.attempt,
      cumulative.evaluations = budget$state$evaluations,
      success = isTRUE(fit$success), status = fit$status)
    if (isTRUE(fit$success) || !fit$status %in% eligible ||
        refinement >= control$max.refinements) break
    next.points <- 2 * (as.double(points) - 1) + 1
    if (!is.finite(next.points) || next.points > .Machine$integer.max ||
        next.points > control$max.evaluations) {
      fit <- limit("A finer Conditional grid exceeds the cumulative point-evaluation budget or R's integer range.")
      requested <- initial$domain
      break
    }
    points <- next.points
    refinement <- refinement + 1L
  }
  diagnostics$search.domain <- searched
  diagnostics$requested.domain <- requested
  diagnostics$domain.stable <- isTRUE(fit$success) && streak >= 2L
  diagnostics$domain.expansions <- widest.level
  diagnostics$agreement.streak <- streak
  diagnostics$domain.history <- if (length(history)) do.call(rbind, history) else data.frame()
  diagnostics$search.attempts <- do.call(rbind, attempts)
  diagnostics$refinements <- refinement
  diagnostics$final.grid.points <- points
  diagnostics$test.evaluations <- budget$state$evaluations
  diagnostics$endpoint.certified <- endpoint.ok
  diagnostics$resolution.stable <- resolution
  .hc_engine_result(fit, plugin, diagnostics)
}
