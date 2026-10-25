.hc_estimate_cs_nyt <- function(work, info) {
  .hc_estimate_cs(work, info, never_only = FALSE)
}

.hc_estimate_cs <- function(work, info, never_only) {
  rows <- vector("list", nrow(info))
  counts <- data.frame(coordinate = info$coordinate, treated_pairs = 0L,
                       comparison_pairs = 0L, stringsAsFactors = FALSE)
  for (j in seq_len(nrow(info))) {
    g <- info$cohort[j]; tt <- info$time_index[j]; base <- g - 1L
    if (tt == base) {
      rows[[j]] <- list(i = integer(0), x = numeric(0))
      next
    }
    current <- which(work$time_index == tt)
    baseline <- which(work$time_index == base)
    matched <- match(work$unit_index[current], work$unit_index[baseline])
    keep <- !is.na(matched)
    current <- current[keep]; baseline <- baseline[matched[keep]]
    treated <- which(work$cohort[current] == g)
    threshold <- if (tt < g) g else tt
    comparison <- if (never_only) which(work$cohort[current] == 0L) else
      which(work$cohort[current] == 0L | work$cohort[current] > threshold)
    .hc_assert(length(treated) > 0L, paste0("No paired treated outcomes for cohort ", g,
      " between time indices ", base, " and ", tt, "."))
    .hc_assert(length(comparison) > 0L, paste0("No paired comparison outcomes for cohort ", g,
      " between time indices ", base, " and ", tt, "."))
    tr <- work$unit_weight[current[treated]]; co <- work$unit_weight[current[comparison]]
    wt <- tr / sum(tr); wc <- co / sum(co)
    rows[[j]] <- list(i = c(current[treated], baseline[treated],
                            current[comparison], baseline[comparison]),
                      x = c(wt, -wt, -wc, wc))
    counts$treated_pairs[j] <- length(treated)
    counts$comparison_pairs[j] <- length(comparison)
  }
  list(A = .hc_sparse_rows(rows, nrow(work)), pair_counts = counts)
}
