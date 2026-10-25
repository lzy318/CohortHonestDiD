# Internal entry point shared by the panel and external-estimate interfaces.
.hc_infer_one <- function(representation, restriction, M, method, full.path,
                          alpha, seed = NULL, control = list()) {
  answer <- .hc_infer_one_impl(representation, restriction, M, method, full.path,
                               alpha, seed, control)
  if (identical(method, "conditional")) answer$diagnostics <-
    .hc_conditional_merge_diagnostics(answer$diagnostics)
  answer
}

.hc_infer_one_impl <- function(representation, restriction, M, method, full.path,
                               alpha, seed = NULL, control = list()) {
  tryCatch({
    restriction <- match.arg(toupper(restriction), c("SD", "RM"))
    method <- match.arg(method, c("conditional", "c_lf", "flci"))
    if (length(M) != 1L || !is.finite(M) || M < 0)
      stop("The sensitivity bound must be finite and nonnegative.", call. = FALSE)
    if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 0.5)
      stop("alpha must lie in (0, 0.5).", call. = FALSE)
    if (restriction == "SD" &&
        (length(full.path) != 1L || !is.logical(full.path) || is.na(full.path)))
      stop("full.path must be TRUE or FALSE.", call. = FALSE)
    if (restriction == "RM" && method == "flci")
      stop("FLCI is available for SD, not RM.", call. = FALSE)
    ctl <- .hc_engine_controls(control)
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
    rep <- .hc_engine_representation(representation)
    if (method != "flci" && as.numeric(crossprod(rep$target, rep$Sigma %*% rep$target)) == 0)
      return(.hc_engine_failure("unsupported_zero_variance",
        "Conditional and C-LF inference is unsupported for a deterministic target with zero sampling variance; no confidence set is reported."))
    if (!is.null(seed) && (length(seed) != 1L || !is.finite(seed) || seed < 0 || seed != floor(seed)))
      stop("seed must be NULL or one nonnegative integer.", call. = FALSE)
    run <- function() {
      task.seed <- if (is.null(seed)) sample.int(1000000000L, 1L) else
        as.integer(seed %% 1000000000L)
      answer <- if (restriction == "SD") {
        .hc_infer_sd(rep, M, method, full.path, alpha, task.seed, ctl)
      } else {
        .hc_infer_rm(rep, M, method, alpha, task.seed, ctl)
      }
      answer$diagnostics$restriction <- restriction
      answer$diagnostics$method <- method
      answer$diagnostics$full.path <- if (restriction == "SD") full.path else NA
      answer$diagnostics$seed <- task.seed
      answer
    }
    if (is.null(seed)) run() else withr::with_seed(as.integer(seed %% 1000000000L), run())
  }, error = function(e) .hc_engine_failure(
    if (inherits(e, "sd_unsupported_moment_variance")) "unsupported_moment_variance" else "inference_error",
    conditionMessage(e)))
}
