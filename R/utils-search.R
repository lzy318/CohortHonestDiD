# One cumulative budget for confidence-set point evaluations, including audits.
.hc_search_budget <- function(limit) {
  state <- new.env(parent = emptyenv())
  state$cache <- new.env(parent = emptyenv())
  state$evaluations <- 0L
  state$exhausted <- FALSE
  evaluate <- function(theta, test, context = "point") {
    key <- paste0(context, ":", sprintf("%.17g", theta))
    if (exists(key, state$cache, inherits = FALSE)) return(get(key, state$cache))
    if (state$evaluations >= limit) {
      state$exhausted <- TRUE
      stop(structure(list(message = paste0("The cumulative limit of ", limit,
        " confidence-set point evaluations was reached."), call = NULL),
        class = c("hc_search_budget_exhausted", "error", "condition")))
    }
    state$evaluations <- state$evaluations + 1L
    value <- test()
    assign(key, value, state$cache)
    value
  }
  list(evaluate = evaluate, state = state)
}

.hc_adaptive_search <- function(control, budget, run_attempt) {
  points <- control$grid.points
  attempt <- 0L
  history <- list()
  eligible <- c("open_search_boundary", "topology_audit_failure",
                "endpoint_audit_failure", "uncertified_empty_search", "topology_miss")
  repeat {
    before <- budget$state$evaluations
    allocation_ok <- is.finite(points) && points <= .Machine$integer.max &&
      points <= control$max.evaluations
    one <- if (allocation_ok) tryCatch(run_attempt(as.integer(points)),
      hc_search_budget_exhausted = function(e) list(fit = .hc_engine_failure(
        "evaluation_limit", conditionMessage(e)), diagnostics = list())) else
      list(fit = .hc_engine_failure("evaluation_limit",
        "The requested refinement grid exceeds the cumulative point-evaluation budget or R's integer range."),
        diagnostics = list())
    if (budget$state$exhausted) {
      one$fit$success <- FALSE
      one$fit$status <- "evaluation_limit"
      one$fit$message <- paste0("The cumulative limit of ", control$max.evaluations,
        " confidence-set point evaluations was reached.")
    }
    history[[length(history) + 1L]] <- data.frame(attempt = attempt,
      grid.points = points, new.evaluations = budget$state$evaluations - before,
      cumulative.evaluations = budget$state$evaluations,
      success = isTRUE(one$fit$success), status = one$fit$status)
    if (isTRUE(one$fit$success) || !one$fit$status %in% eligible ||
        attempt >= control$max.refinements) break
    next.points <- 2 * (as.double(points) - 1) + 1
    if (!is.finite(next.points) || next.points > .Machine$integer.max ||
        next.points > control$max.evaluations) {
      one$fit$success <- FALSE
      one$fit$status <- "evaluation_limit"
      one$fit$message <- "A finer grid would exceed the cumulative point-evaluation budget or R's integer range."
      break
    }
    points <- next.points
    attempt <- attempt + 1L
  }
  one$diagnostics$search.attempts <- do.call(rbind, history)
  one$diagnostics$final.grid.points <- points
  one$diagnostics$refinements <- attempt
  one$diagnostics$test.evaluations <- budget$state$evaluations
  one
}
