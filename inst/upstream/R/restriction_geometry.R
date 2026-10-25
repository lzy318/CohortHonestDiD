# Restriction geometry shared by the cohort and aggregate representations.
# All matrices are named.  The path coordinates must use the canonical
# delta_(pre|post)_g<cohort>_s<event-time> convention.

v3_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

v3_max_abs <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) return(0)
  if (any(!is.finite(x))) return(Inf)
  max(abs(x))
}

v3_symmetrize <- function(x) {
  y <- t(x)
  answer <- x / 2 + y / 2
  safe_sum <- sign(x) != sign(y) |
    (abs(x) <= .Machine$double.xmax / 2 & abs(y) <= .Machine$double.xmax / 2)
  answer[safe_sum] <- (x[safe_sum] + y[safe_sum]) / 2
  equal <- x == y
  answer[equal] <- x[equal]
  answer
}

v3_parse_coordinates <- function(coordinate_names) {
  coordinate_names <- as.character(coordinate_names)
  type <- ifelse(
    grepl("^delta_pre_", coordinate_names), "pre",
    ifelse(grepl("^delta_post_", coordinate_names), "post", NA_character_)
  )
  cohort <- suppressWarnings(as.integer(sub(
    "^delta_(?:pre|post)_g([0-9]+)_s-?[0-9]+$", "\\1", coordinate_names
  )))
  relative_period <- suppressWarnings(as.integer(sub(
    "^.*_s(-?[0-9]+)$", "\\1", coordinate_names
  )))
  ans <- data.frame(
    coordinate = coordinate_names,
    type = type,
    cohort = cohort,
    relative_period = relative_period,
    stringsAsFactors = FALSE
  )
  v3_assert(!anyNA(ans), "Could not parse canonical path-coordinate names.")
  v3_assert(!anyDuplicated(ans$coordinate), "Path-coordinate names are duplicated.")
  ans
}

v3_canonical_coordinate_info <- function(adoption_times, T_total) {
  adoption_times <- sort(as.integer(adoption_times))
  T_total <- as.integer(T_total)
  v3_assert(length(adoption_times) > 0L, "adoption_times cannot be empty.")
  v3_assert(all(adoption_times >= 2L & adoption_times <= T_total),
            "Each adoption time must lie in 2,...,T_total.")
  names <- flci_canonical_names(adoption_times, T_total)
  v3_parse_coordinates(names)
}

v3_validate_coordinate_info <- function(coordinate_info) {
  required <- c("coordinate", "type", "cohort", "relative_period")
  v3_assert(is.data.frame(coordinate_info) &&
              all(required %in% names(coordinate_info)),
            "coordinate_info is missing required columns.")
  v3_assert(!anyDuplicated(coordinate_info$coordinate),
            "coordinate_info has duplicate coordinates.")
  v3_assert(all(coordinate_info$type %in% c("pre", "post")),
            "coordinate_info$type must be pre or post.")
  parsed <- v3_parse_coordinates(coordinate_info$coordinate)
  v3_assert(
    identical(as.character(coordinate_info$type), as.character(parsed$type)) &&
      identical(as.integer(coordinate_info$cohort), as.integer(parsed$cohort)) &&
      identical(as.integer(coordinate_info$relative_period),
                as.integer(parsed$relative_period)),
    "coordinate_info metadata disagree with the canonical coordinate names."
  )
  for (g in sort(unique(coordinate_info$cohort))) {
    s <- sort(coordinate_info$relative_period[coordinate_info$cohort == g])
    v3_assert(length(s) >= 3L && all(diff(s) == 1L),
              paste0("Cohort/path ", g,
                     " must contain at least three consecutive event times."))
    v3_assert(any(s <= 0L) && any(s >= 1L),
              paste0("Cohort/path ", g, " must contain pre and post periods."))
  }
  invisible(TRUE)
}

# pre_window (optional integer): drop pure pre-treatment rows that reach a
# period before it.  Rows that involve a post-treatment period are never dropped.
v3_first_difference_operator <- function(coordinate_info, pre_window = NULL) {
  v3_validate_coordinate_info(coordinate_info)
  names_all <- as.character(coordinate_info$coordinate)
  rows <- list()
  metadata <- list()
  cursor <- 0L
  for (g in sort(unique(coordinate_info$cohort))) {
    one <- coordinate_info[coordinate_info$cohort == g, , drop = FALSE]
    one <- one[order(one$relative_period), , drop = FALSE]
    for (j in seq.int(2L, nrow(one))) {
      cursor <- cursor + 1L
      row <- setNames(numeric(length(names_all)), names_all)
      row[one$coordinate[j]] <- 1
      row[one$coordinate[j - 1L]] <- -1
      row_name <- paste0(
        "diff_g", g, "_s", one$relative_period[j],
        "_vs_s", one$relative_period[j - 1L]
      )
      rows[[row_name]] <- row
      metadata[[cursor]] <- data.frame(
        row = row_name,
        cohort = as.integer(g),
        current_period = one$relative_period[j],
        previous_period = one$relative_period[j - 1L],
        component = if (one$relative_period[j] <= 0L) "pre" else "post",
        stringsAsFactors = FALSE
      )
    }
  }
  D1 <- do.call(rbind, rows)
  meta <- do.call(rbind, metadata)
  if (!is.null(pre_window)) {
    keep <- !(meta$component == "pre" & meta$previous_period < pre_window)
    D1 <- D1[keep, , drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
    rownames(meta) <- NULL
  }
  list(D1 = D1, rows = meta)
}

v3_second_difference_operator <- function(coordinate_info, pre_window = NULL) {
  v3_validate_coordinate_info(coordinate_info)
  names_all <- as.character(coordinate_info$coordinate)
  rows <- list()
  metadata <- list()
  cursor <- 0L
  for (g in sort(unique(coordinate_info$cohort))) {
    one <- coordinate_info[coordinate_info$cohort == g, , drop = FALSE]
    one <- one[order(one$relative_period), , drop = FALSE]
    for (j in seq.int(3L, nrow(one))) {
      cursor <- cursor + 1L
      row <- setNames(numeric(length(names_all)), names_all)
      row[one$coordinate[(j - 2L):j]] <- c(1, -2, 1)
      row_name <- paste0("d2_g", g, "_s", one$relative_period[j])
      rows[[row_name]] <- row
      metadata[[cursor]] <- data.frame(
        row = row_name,
        cohort = as.integer(g),
        end_period = one$relative_period[j],
        component = if (one$relative_period[j] >= 1L) {
          "post_connected"
        } else {
          "pure_pre"
        },
        stringsAsFactors = FALSE
      )
    }
  }
  D2 <- do.call(rbind, rows)
  meta <- do.call(rbind, metadata)
  if (!is.null(pre_window)) {
    keep <- !(meta$component == "pure_pre" & meta$end_period - 2L < pre_window)
    D2 <- D2[keep, , drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
    rownames(meta) <- NULL
  }
  list(D2 = D2, rows = meta)
}

v3_sd_restriction <- function(coordinate_info, M_raw, pre_window = NULL) {
  M_raw <- as.numeric(M_raw)
  v3_assert(length(M_raw) == 1L && is.finite(M_raw) && M_raw >= 0,
            "M_raw must be one finite nonnegative scalar.")
  operator <- v3_second_difference_operator(coordinate_info, pre_window = pre_window)
  positive <- operator$D2
  negative <- -operator$D2
  rownames(positive) <- paste0(rownames(positive), "_pos")
  rownames(negative) <- paste0(rownames(operator$D2), "_neg")
  A_Delta <- rbind(positive, negative)
  d_vec <- setNames(rep(M_raw, nrow(A_Delta)), rownames(A_Delta))
  list(A_Delta = A_Delta, d_vec = d_vec, operator = operator)
}

v3_rm_scenario_table <- function(coordinate_info, Mbar, pre_window = NULL) {
  Mbar <- as.numeric(Mbar)
  v3_assert(length(Mbar) == 1L && is.finite(Mbar) && Mbar >= 0,
            "Mbar must be one finite nonnegative scalar.")
  operator <- v3_first_difference_operator(coordinate_info, pre_window = pre_window)
  choices <- list()
  for (g in sort(unique(operator$rows$cohort))) {
    pre <- operator$rows$cohort == g & operator$rows$component == "pre"
    candidates <- operator$rows$row[pre]
    v3_assert(length(candidates) > 0L,
              paste0("Path ", g, " has no pre-period RM benchmark."))
    choices[[paste0("g", g)]] <- paste(
      rep(candidates, each = 2L), rep(c(-1, 1), times = length(candidates)),
      sep = "@@"
    )
  }
  if (Mbar == 0) {
    table <- data.frame(component = 1L)
    full_count <- prod(vapply(choices, length, integer(1L)))
  } else {
    table <- expand.grid(choices, KEEP.OUT.ATTRS = FALSE,
                         stringsAsFactors = FALSE)
    table$component <- seq_len(nrow(table))
    full_count <- nrow(table)
  }
  list(table = table, operator = operator, scenario_count_full = full_count)
}

v3_rm_component_restriction <- function(
    coordinate_info, Mbar, scenario_row, operator = NULL, pre_window = NULL) {
  if (is.null(operator)) {
    operator <- v3_first_difference_operator(coordinate_info, pre_window = pre_window)
  }
  D1 <- operator$D1
  meta <- operator$rows
  pieces <- list()
  cursor <- 0L
  for (g in sort(unique(meta$cohort))) {
    post <- meta$cohort == g & meta$component == "post"
    Dpost <- D1[post, , drop = FALSE]
    if (Mbar == 0) {
      benchmark_name <- "zero_Mbar"
      benchmark <- setNames(numeric(ncol(D1)), colnames(D1))
    } else {
      encoded <- as.character(scenario_row[[paste0("g", g)]])
      decoded <- strsplit(encoded, "@@", fixed = TRUE)[[1L]]
      v3_assert(length(decoded) == 2L && decoded[1L] %in% rownames(D1),
                "Malformed RM benchmark scenario.")
      benchmark_name <- decoded[1L]
      benchmark <- as.numeric(decoded[2L]) * D1[benchmark_name, ]
    }
    upper <- sweep(Dpost, 2L, Mbar * benchmark, "-")
    lower <- sweep(-Dpost, 2L, Mbar * benchmark, "-")
    rownames(upper) <- paste0(
      rownames(Dpost), "_upper_bm_", make.names(benchmark_name)
    )
    rownames(lower) <- paste0(
      rownames(Dpost), "_lower_bm_", make.names(benchmark_name)
    )
    cursor <- cursor + 1L
    pieces[[cursor]] <- upper
    cursor <- cursor + 1L
    pieces[[cursor]] <- lower
  }
  A_Delta <- do.call(rbind, pieces)
  v3_assert(!anyDuplicated(rownames(A_Delta)),
            "RM component has duplicate restriction names.")
  list(
    A_Delta = A_Delta,
    d_vec = setNames(rep(0, nrow(A_Delta)), rownames(A_Delta)),
    operator = operator
  )
}

v3_rm_plugin_restriction <- function(coordinate_info, beta, Mbar, pre_window = NULL) {
  operator <- v3_first_difference_operator(coordinate_info, pre_window = pre_window)
  beta <- beta[colnames(operator$D1)]
  v3_assert(!anyNA(beta), "RM plug-in beta is not aligned to path coordinates.")
  D1 <- operator$D1
  meta <- operator$rows
  pieces <- list()
  bounds <- numeric()
  for (g in sort(unique(meta$cohort))) {
    pre <- meta$cohort == g & meta$component == "pre"
    post <- meta$cohort == g & meta$component == "post"
    benchmark <- max(abs(as.numeric(D1[pre, , drop = FALSE] %*% beta)))
    Dpost <- D1[post, , drop = FALSE]
    pieces[[paste0("g", g, "_pos")]] <- Dpost
    pieces[[paste0("g", g, "_neg")]] <- -Dpost
    bounds <- c(bounds, rep(Mbar * benchmark, 2L * nrow(Dpost)))
  }
  A_Delta <- do.call(rbind, pieces)
  rownames(A_Delta) <- make.unique(paste0("rm_plugin_", seq_len(nrow(A_Delta))))
  list(
    A_Delta = A_Delta,
    d_vec = setNames(bounds, rownames(A_Delta)),
    operator = operator
  )
}
