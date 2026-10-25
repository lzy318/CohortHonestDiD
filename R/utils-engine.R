# Private numerical resources. No source file is loaded into the user's session.
.hc_engine_state <- new.env(parent = emptyenv())

# The same normalization is used for primitive biases and coefficients.
# Consequently W must map its null space into itself before quotient reduction.
.hc_check_mapping_support <- function(W, normalization, tolerance = 1e-9) {
  N <- as.matrix(normalization)
  if (!nrow(N)) return(invisible(TRUE))
  row.scale <- apply(abs(N), 1L, max)
  N <- N[row.scale > 0, , drop = FALSE]
  if (!nrow(N)) return(invisible(TRUE))
  N <- N / row.scale[row.scale > 0]
  decomposition <- qr(t(N), tol = tolerance, LAPACK = FALSE)
  if (decomposition$rank == ncol(N)) return(invisible(TRUE))
  Q <- qr.Q(decomposition, complete = TRUE)
  constrained <- Q[, seq_len(decomposition$rank), drop = FALSE]
  free <- Q[, seq.int(decomposition$rank + 1L, ncol(N)), drop = FALSE]
  W.scale <- max(abs(W))
  if (!is.finite(W.scale) || W.scale <= 0)
    stop("W must have a positive finite scale.", call. = FALSE)
  scaled.W <- W / W.scale
  residual <- norm(crossprod(constrained, scaled.W %*% free), type = "F")
  scale <- norm(scaled.W, type = "F")
  if (!is.finite(residual) || residual > tolerance * scale)
    stop("W does not preserve the supplied normalization: N W E is nonzero for a basis E of ker(N).", call. = FALSE)
  invisible(TRUE)
}

.hc_backend <- function() {
  if (exists("engine", .hc_engine_state, inherits = FALSE))
    return(.hc_engine_state$engine)
  package.path <- getNamespaceInfo(asNamespace("CohortHonestDiD"), "path")
  paths <- file.path(package.path, c("upstream/R", "inst/upstream/R"))
  path <- paths[dir.exists(paths)][1L]
  if (is.na(path)) stop("Package numerical resources are unavailable.", call. = FALSE)
  env <- new.env(parent = baseenv())
  env$setNames <- stats::setNames
  env$tail <- utils::tail
  env$head <- utils::head
  env$`%||%` <- function(x, y) if (is.null(x)) y else x
  for (file in c("flci_core.R", "arp-nuisance.R", "inference_engine.R",
                 "restriction_geometry.R", "inference_adapter.R")) {
    sys.source(file.path(path, file), envir = env)
  }
  env$sd_inference_load_dependencies()
  .hc_engine_state$engine <- env
  env
}

.hc_engine_controls <- function(control = list()) {
  defaults <- list(lf.sims = 2000L, grid.points = 101L, endpoint.tol = NULL,
                   grid.radius = 8, max.expansions = 6L, max.refinements = 3L,
                   max.evaluations = 100000L, topology.audit = TRUE,
                   flci.solver = "direct", flci.control = list(),
                   max.components = 100000L, component.cores = 1L)
  if (!is.list(control)) stop("control must be a list.", call. = FALSE)
  if (length(control) && (is.null(names(control)) || anyNA(names(control)) || any(!nzchar(names(control))) ||
                         anyDuplicated(names(control))))
    stop("control must be a uniquely named list.", call. = FALSE)
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) stop("Unknown control: ", paste(unknown, collapse = ", "), call. = FALSE)
  ans <- utils::modifyList(defaults, control, keep.null = TRUE)
  for (key in c("lf.sims", "grid.points", "max.evaluations", "max.components", "component.cores")) {
    ans[[key]] <- .hc_integer_scalar(ans[[key]], key)
  }
  if (ans$grid.points < 5L) stop("grid.points must be at least 5.", call. = FALSE)
  positive_scalar <- function(x) is.numeric(x) && !is.complex(x) &&
    length(x) == 1L && is.finite(x) && x > 0
  if (!positive_scalar(ans$grid.radius))
    stop("grid.radius must be positive.", call. = FALSE)
  if (!is.null(ans$endpoint.tol) &&
      !positive_scalar(ans$endpoint.tol))
    stop("endpoint.tol must be NULL or positive in outcome units.", call. = FALSE)
  ans$max.expansions <- .hc_integer_scalar(ans$max.expansions, "max.expansions", 0L)
  ans$max.refinements <- .hc_integer_scalar(ans$max.refinements, "max.refinements", 0L)
  if (!is.logical(ans$topology.audit) || length(ans$topology.audit) != 1L ||
      is.na(ans$topology.audit)) stop("topology.audit must be TRUE or FALSE.", call. = FALSE)
  .hc_assert(is.character(ans$flci.solver) && length(ans$flci.solver) == 1L &&
               !is.na(ans$flci.solver), "flci.solver must be one character string.")
  ans$flci.solver <- match.arg(ans$flci.solver, c("direct", "nested"))
  fc <- if ("flci.control" %in% names(control)) control$flci.control else ans$flci.control
  if (!is.list(fc)) stop("flci.control must be a list.", call. = FALSE)
  if (length(fc) && (is.null(names(fc)) || anyNA(names(fc)) ||
                    any(!nzchar(names(fc))) || anyDuplicated(names(fc))))
    stop("flci.control must be a uniquely named list.", call. = FALSE)
  allowed <- if (ans$flci.solver == "direct") c("n_random_starts", "maxit") else
    c("solver", "solver_opts", "outer_tol", "grid_size", "feasibility_tol", "objective_tol")
  unknown <- setdiff(names(fc), allowed)
  if (length(unknown)) stop("Unknown flci.control: ", paste(unknown, collapse = ", "), call. = FALSE)
  for (key in intersect(names(fc), c("n_random_starts", "maxit", "grid_size"))) {
    minimum <- switch(key, n_random_starts = 0L, maxit = 0L, grid_size = 3L)
    fc[[key]] <- .hc_integer_scalar(fc[[key]], paste0("flci.control$", key), minimum)
  }
  for (key in intersect(names(fc), c("outer_tol", "feasibility_tol", "objective_tol"))) {
    if (!positive_scalar(fc[[key]])) stop("flci.control$", key, " must be positive.", call. = FALSE)
  }
  if ("solver" %in% names(fc) && (!is.character(fc$solver) || length(fc$solver) != 1L ||
                                  is.na(fc$solver) || !nzchar(fc$solver)))
    stop("flci.control$solver must be one nonempty string.", call. = FALSE)
  if ("solver_opts" %in% names(fc) && (!is.list(fc$solver_opts) ||
      (length(fc$solver_opts) && (is.null(names(fc$solver_opts)) || anyNA(names(fc$solver_opts)) ||
        any(!nzchar(names(fc$solver_opts))) || anyDuplicated(names(fc$solver_opts))))))
    stop("flci.control$solver_opts must be a uniquely named list.", call. = FALSE)
  ans$flci.control <- fc
  ans
}

.hc_engine_representation <- function(representation) {
  info <- representation$coordinate_info
  required <- c("coordinate", "type", "cohort", "relative_period")
  if (!is.data.frame(info) || !all(required %in% names(info)))
    stop("coordinate_info must identify every coefficient's cohort and relative period.", call. = FALSE)
  old <- as.character(info$coordinate)
  if (anyNA(info[, required]) || anyDuplicated(old))
    stop("Coefficient coordinates must be unique and nonmissing.", call. = FALSE)
  canonical <- paste0("delta_", info$type, "_g", info$cohort, "_s", info$relative_period)
  if (anyDuplicated(canonical)) stop("Cohort-period cells are duplicated.", call. = FALSE)
  for (key in c("betahat", "target")) {
    x <- representation[[key]]
    if (is.null(names(x)) || anyDuplicated(names(x)) || !setequal(names(x), old))
      stop(key, " must be named and aligned with coordinate_info.", call. = FALSE)
    representation[[key]] <- stats::setNames(as.numeric(x[old]), canonical)
  }
  for (key in c("Sigma", "W")) {
    x <- as.matrix(representation[[key]])
    if (!setequal(rownames(x), old) || !setequal(colnames(x), old))
      stop(key, " must have matching row and column names.", call. = FALSE)
    x <- x[old, old, drop = FALSE]
    dimnames(x) <- list(canonical, canonical)
    representation[[key]] <- x
  }
  normal <- as.matrix(representation$normalization)
  if (ncol(normal) != length(old) || !setequal(colnames(normal), old))
    stop("normalization must have a column for every coefficient.", call. = FALSE)
  normal <- normal[, old, drop = FALSE]
  colnames(normal) <- canonical
  rownames(normal) <- if (nrow(normal)) paste0("normalization_", seq_len(nrow(normal))) else NULL
  representation$normalization <- normal
  .hc_check_mapping_support(representation$W, normal)
  representation$coordinate_info$coordinate <- canonical
  representation$.original_names <- stats::setNames(old, canonical)
  .hc_backend()$v3_align_representation(representation)
}

.hc_engine_failure <- function(status, message, diagnostics = list(), plugin = c(NA_real_, NA_real_)) {
  list(success = FALSE, status = status, message = message,
       components = data.frame(lower = numeric(), upper = numeric()),
       center = NA_real_, plugin = stats::setNames(plugin, c("lower", "upper")),
       diagnostics = diagnostics)
}

.hc_engine_result <- function(fit, plugin, diagnostics = list()) {
  if (!isTRUE(fit$success)) {
    diagnostics$partial.components <- fit$components
    diagnostics$warnings <- unique(c(diagnostics$warnings, fit$warnings))
    return(.hc_engine_failure(fit$status, fit$message, diagnostics,
                             c(plugin$lower, plugin$upper)))
  }
  components <- fit$components[, c("lower", "upper"), drop = FALSE]
  center <- if (!is.null(fit$center)) fit$center else if (
    all(is.finite(c(plugin$lower, plugin$upper)))) {
    (plugin$lower + plugin$upper) / 2
  } else NA_real_
  diagnostics$warnings <- unique(c(diagnostics$warnings, fit$warnings))
  diagnostics$plugin.empty <- isTRUE(plugin$empty)
  list(success = TRUE,
       status = if (!nrow(components)) "empty" else if (any(is.infinite(as.matrix(components))))
         "unbounded" else if (nrow(components) > 1L) "disconnected" else "ok",
       message = fit$message, components = components, center = center,
       plugin = c(lower = plugin$lower, upper = plugin$upper), diagnostics = diagnostics)
}

.hc_engine_grid <- function(built, plugin, anchors, control, bounds = NULL) {
  se <- sqrt(as.numeric(crossprod(built$target_bar, built$Sigma_bar %*% built$target_bar)))
  if (!is.finite(se) || se <= 0) stop("The target's sampling variance must be positive.", call. = FALSE)
  points <- c(plugin$lower, plugin$upper, anchors)
  points <- points[is.finite(points)]
  if (!length(points)) points <- sum(built$target_bar * built$beta_bar)
  lower <- min(points) - control$grid.radius * se
  upper <- max(points) + control$grid.radius * se
  tolerance <- if (is.null(control$endpoint.tol)) 1e-5 * se else control$endpoint.tol
  if (!is.null(bounds)) {
    margin <- max(10 * tolerance, 1e-4 * se)
    lower <- bounds[1L] - margin
    upper <- bounds[2L] + margin
    anchors <- anchors[anchors >= lower & anchors <= upper]
  }
  grid <- sort(unique(c(seq(lower, upper, length.out = control$grid.points), anchors)))
  grid <- grid[is.finite(grid)]
  list(grid = grid, spacing = (upper - lower) / (control$grid.points - 1L),
       se = se, tolerance = tolerance)
}
