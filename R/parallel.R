.hc_available_cores <- function() {
  n <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (length(n) != 1L || is.na(n) || !is.finite(n)) return(1L)
  max(1L, as.integer(n) - 1L)
}

.hc_task_seeds <- function(n, seed) {
  if (is.null(seed)) return(sample.int(.Machine$integer.max, n, replace = FALSE))
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("seed must be NULL or a nonnegative integer.", call. = FALSE)
  }
  # A local deterministic schedule makes results independent of worker count.
  as.integer((as.double(seed) + seq_len(n) * 104729) %% (.Machine$integer.max - 1)) + 1L
}

.hc_task <- function(task) {
  # Each task gets a private RNG stream, restored also in serial mode.
  withr::with_seed(task$seed, tryCatch(do.call(.hc_infer_one, task), error = function(e) {
    list(success = FALSE, status = "numerical_failure", message = conditionMessage(e),
         components = data.frame(lower = numeric(), upper = numeric()),
         center = NA_real_, plugin = c(NA_real_, NA_real_), diagnostics = list())
  }), .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
  .rng_sample_kind = "Rejection")
}

.hc_run_tasks <- function(tasks, parallel, cores) {
  if (!parallel || length(tasks) < 2L || cores == 1L) return(lapply(tasks, .hc_task))
  workers <- min(cores, length(tasks))
  cl <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  libs <- .libPaths()
  parallel::clusterCall(cl, function(paths) { .libPaths(paths); loadNamespace("CohortHonestDiD"); NULL }, libs)
  parallel::parLapplyLB(cl, tasks, function(task) {
    get(".hc_task", envir = asNamespace("CohortHonestDiD"))(task)
  })
}
