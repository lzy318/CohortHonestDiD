.hc_restrictions <- c("SD", "RM", "SDB", "SDM", "RMB", "RMM", "SDRM", "SDRMB", "SDRMM")

# base: the smoothness (SD) or union (RM, SDRM) rows. shape: the extra rows.
.hc_restriction_spec <- function(restriction, direction = NULL) {
  base <- if (restriction %in% c("SD", "SDB", "SDM")) "SD" else
    if (restriction %in% c("RM", "RMB", "RMM")) "RM" else "SDRM"
  shape <- if (restriction %in% c("SDB", "RMB", "SDRMB")) "sign" else
    if (restriction %in% c("SDM", "RMM", "SDRMM")) "monotone" else "none"
  list(restriction = restriction, base = base, shape = shape, direction = direction)
}

# A wrong or missing direction stops before any task starts.
.hc_check_direction <- function(spec, direction) {
  if (spec$shape == "none") {
    if (!is.null(direction)) stop("direction applies only to SDB, RMB, and SDRMB (positive or negative) ",
                                  "and to SDM, RMM, and SDRMM (increasing or decreasing).", call. = FALSE)
    return(invisible(NULL))
  }
  allowed <- if (spec$shape == "sign") c("positive", "negative") else c("increasing", "decreasing")
  if (!is.character(direction) || length(direction) != 1L || is.na(direction) || !direction %in% allowed)
    stop(spec$restriction, " requires direction = '", allowed[1L], "' or '", allowed[2L], "'.", call. = FALSE)
  invisible(NULL)
}

# Internal entry point shared by the panel and external-estimate interfaces.
.hc_infer_one <- function(representation, restriction, M, method, full.path,
                          alpha, seed = NULL, control = list(), direction = NULL) {
  answer <- .hc_infer_one_impl(representation, restriction, M, method, full.path,
                               alpha, seed, control, direction)
  if (identical(method, "conditional")) answer$diagnostics <-
    .hc_conditional_merge_diagnostics(answer$diagnostics)
  answer
}

.hc_infer_one_impl <- function(representation, restriction, M, method, full.path,
                               alpha, seed = NULL, control = list(), direction = NULL) {
  tryCatch({
    restriction <- match.arg(toupper(restriction), .hc_restrictions)
    spec <- .hc_restriction_spec(restriction, direction)
    method <- match.arg(method, c("conditional", "c_lf", "flci"))
    if (length(M) != 1L || !is.finite(M) || M < 0)
      stop("The sensitivity bound must be finite and nonnegative.", call. = FALSE)
    if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 0.5)
      stop("alpha must lie in (0, 0.5).", call. = FALSE)
    if (spec$base == "SD" &&
        (length(full.path) != 1L || !is.logical(full.path) || is.na(full.path)))
      stop("full.path must be TRUE or FALSE.", call. = FALSE)
    if (restriction != "SD" && method == "flci")
      stop("FLCI is available for SD, not ", restriction, ".", call. = FALSE)
    ctl <- .hc_engine_controls(control)
    # Sign and monotonicity rows split the Conditional confidence set into many
    # pieces, and its meshes rarely agree. Run the mesh audit only on request.
    if (method == "conditional" && spec$shape != "none" && !"topology.audit" %in% names(control))
      ctl$topology.audit <- FALSE
    if (method == "c_lf" && ctl$grid.points > ctl$max.evaluations)
      return(.hc_engine_failure("evaluation_limit",
        "The requested initial grid exceeds the cumulative point-evaluation budget.",
        list(test.evaluations = 0L, refinements = 0L,
          final.grid.points = ctl$grid.points,
          search.attempts = data.frame(attempt = 0L, grid.points = ctl$grid.points,
            new.evaluations = 0L, cumulative.evaluations = 0L,
            success = FALSE, status = "evaluation_limit"))))
    info <- representation$coordinate_info
    pre.count <- vapply(unique(info$cohort), function(g)
      sum(info$cohort == g & info$type == "pre"), integer(1))
    if (any(pre.count < 2L)) return(.hc_engine_failure(
      "insufficient_prehistory",
      "Sensitivity analysis requires at least two pre-treatment periods for every included cohort/path."))
    if (spec$base == "SDRM" && any(pre.count < 3L)) return(.hc_engine_failure(
      "insufficient_prehistory",
      paste0(restriction, " requires at least three pre-treatment periods, counting period 0, ",
             "for every included cohort/path.")))
    rep <- .hc_engine_representation(representation)
    rep$restriction_spec <- spec
    if (method != "flci" && as.numeric(crossprod(rep$target, rep$Sigma %*% rep$target)) == 0)
      return(.hc_engine_failure("unsupported_zero_variance",
        "Conditional and C-LF inference is unsupported for a deterministic target with zero sampling variance; no confidence set is reported."))
    if (!is.null(seed) && (length(seed) != 1L || !is.finite(seed) || seed < 0 || seed != floor(seed)))
      stop("seed must be NULL or one nonnegative integer.", call. = FALSE)
    run <- function() {
      task.seed <- if (is.null(seed)) sample.int(1000000000L, 1L) else
        as.integer(seed %% 1000000000L)
      answer <- if (spec$base == "SD") {
        .hc_infer_sd(rep, M, method, full.path, alpha, task.seed, ctl)
      } else {
        .hc_infer_rm(rep, M, method, alpha, task.seed, ctl)
      }
      answer$diagnostics$restriction <- restriction
      answer$diagnostics$method <- method
      answer$diagnostics$full.path <- if (spec$base == "SD") full.path else NA
      answer$diagnostics$seed <- task.seed
      if (!is.null(direction)) answer$diagnostics$direction <- direction
      answer
    }
    if (is.null(seed)) run() else withr::with_seed(as.integer(seed %% 1000000000L), run())
  }, error = function(e) .hc_engine_failure(
    if (inherits(e, "sd_unsupported_moment_variance")) "unsupported_moment_variance" else "inference_error",
    conditionMessage(e)))
}
