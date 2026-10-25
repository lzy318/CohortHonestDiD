.hc_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

.hc_estimator <- function(estimator) {
  .hc_assert(is.character(estimator) && length(estimator) == 1L &&
               !is.na(estimator), "estimator must be one character string.")
  aliases <- c(imp = "imp", imputation = "imp", cs_nyt = "cs_nyt",
               `cs-not-yet-treated` = "cs_nyt", cs_nt = "cs_nt",
               `cs-never-treated` = "cs_nt")
  value <- unname(aliases[tolower(estimator)])
  .hc_assert(length(value) == 1L && !is.na(value),
             "estimator must be 'imp', 'cs_nyt', or 'cs_nt'.")
  value
}

.hc_integer_scalar <- function(x, name, minimum = 1L) {
  .hc_assert(is.numeric(x) && !is.complex(x) && length(x) == 1L &&
               is.finite(x) && x >= minimum && x <= .Machine$integer.max &&
               x == floor(x),
             paste0(name, " must be an integer between ", minimum,
                    " and ", .Machine$integer.max, "."))
  as.integer(x)
}

.hc_column <- function(value, data, name, allow_vector = FALSE) {
  if (is.character(value) && length(value) == 1L && value %in% names(data)) {
    return(data[[value]])
  }
  if (isTRUE(allow_vector) && length(value) == nrow(data)) return(value)
  stop(name, " must name a column", if (allow_vector) " or be a row-aligned vector" else "",
       ".", call. = FALSE)
}

.hc_unit_values <- function(x, unit_index, name) {
  groups <- split(x, unit_index)
  .hc_assert(!anyNA(x), paste0(name, " cannot be missing on observed rows."))
  .hc_assert(all(vapply(groups, function(z) length(unique(z)) == 1L, logical(1))),
             paste0(name, " must be constant within each unit."))
  x[match(seq_len(max(unit_index)), unit_index)]
}

.hc_coordinate_info <- function(cohort_times, time_labels) {
  make_block <- function(type) {
    do.call(rbind, lapply(cohort_times, function(g) {
      tt <- if (type == "pre") seq_len(g - 1L) else seq.int(g, length(time_labels))
      s <- tt - g + 1L
      data.frame(coordinate = paste0("delta_", type, "_g", g, "_s", s),
                 type = type, cohort = as.integer(g), relative_period = s,
                 time_index = tt, calendar_time = time_labels[tt],
                 event_time = s - 1L, stringsAsFactors = FALSE)
    }))
  }
  out <- rbind(make_block("pre"), make_block("post"))
  rownames(out) <- NULL
  out
}

.hc_normalization <- function(info, estimator) {
  groups <- sort(unique(info$cohort))
  out <- matrix(0, length(groups), nrow(info),
                dimnames = list(paste0("normalization_g", groups), info$coordinate))
  for (j in seq_along(groups)) {
    use <- which(info$cohort == groups[j] & info$type == "pre")
    if (estimator == "imp") out[j, use] <- 1 / length(use)
    else out[j, use[info$relative_period[use] == 0L]] <- 1
  }
  out
}

.hc_check_cells <- function(work, info) {
  for (j in seq_len(nrow(info))) {
    use <- work$cohort == info$cohort[j] & work$time_index == info$time_index[j]
    .hc_assert(any(use), paste0("No observed outcome with known treatment status in cohort ",
      info$cohort[j], " at time ", as.character(info$calendar_time[j]),
      ". This cohort-time cell is not estimable; restrict the analysis window or supply more data."))
  }
}
