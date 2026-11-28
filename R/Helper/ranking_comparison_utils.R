# ==============================================================================
# ranking_comparison_utils.R
#
# Reusable functions for comparing two COMPLETE rankings over the SAME domain
# of n candidates. Shared by:
#   - RQ1: ranking difference (AHP-CDS vs FC-only / Centrality-only / Equal-weight / RWR)
#   - RQ3: incremental value of AHP weighting (AHP-CDS vs Equal-weight, same criteria)
#
# CONVENTIONS
#   - rank 1 = best/top candidate (lower rank number = higher priority)
#   - Delta_r_i = r_baseline,i - r_AHP,i
#       Delta_r_i > 0  -> candidate PROMOTED by AHP-CDS (moved toward rank 1)
#       Delta_r_i < 0  -> candidate DEMOTED by AHP-CDS
#       Delta_r_i = 0  -> unchanged
#   - Both rankings passed to any function here MUST already be restricted to
#     the same domain (same set/length of ids). These functions do not
#     re-derive ranks or handle disjoint domains.
#
# RANK CONVENTION (shared by RQ1, RQ2, RQ3; 01_ablation_conditions.R builds the same
# ranks with dplyr::row_number() on CDS-sorted input, and RQ5 reads them from its output)
#   - Ranks are ORDINAL permutations of 1..n: no ties, no fractional ranks.
#   - Ties in the underlying score are broken by INPUT ORDER (first come, first
#     ranked). Every script must therefore pass candidates in the same fixed
#     order: descending AHP-CDS (assert_sorted_desc() checks this). The ordinal
#     ranks are built ONLY with rank_scores(); never call
#     rank(..., ties.method = "average") directly.
#   - compute_delta_r() (and therefore summarize_pair()) stops with an error if
#     either rank vector is not an ordinal permutation, so an accidental
#     average-rank input cannot pass silently.
#   - Consequence to disclose in Methods: where a baseline has many tied scores
#     (Centrality-only = degree: 135/143 candidates sit in tie blocks), the
#     input-order tie-break copies the AHP-CDS order inside each block. This
#     makes that baseline look slightly MORE similar to AHP-CDS than a neutral
#     tie treatment would, i.e. it is conservative for the claim that AHP-CDS
#     re-prioritizes relative to the baseline.
#
# RBO NOTE
#   Because both rankings are COMPLETE permutations over the same finite domain
#   (no unseen tail, no disjointness), RBO reduces to the closed-form finite sum
#   (Webber, Moffat & Zobel 2010, formula 32) with no extrapolation needed:
#
#       RBO(S,T,p) = (1-p) * sum_{d=1}^{n} p^(d-1) * (|S_d intersect T_d| / d)
#
#   OFFICIAL backend = "manual": this closed form, computed directly in base R.
#   For a fully-observed, equal-domain comparison like ours, this is EXACT
#   (not an approximation) -- there is no unobserved tail to extrapolate.
#   No external dependencies required.
#
#   backend = "gespeR" is kept only for reference / in case the package is
#   reinstated. As of Bioconductor 3.21, gespeR was flagged "Unresponsive"
#   and removed (its dependency cellHTS2 is unmaintained), so this backend
#   will typically fail with a clear error pointing back to "manual".
# ==============================================================================

# ------------------------------------------------------------------------------
# Rank construction and guards (see RANK CONVENTION in the header)
# ------------------------------------------------------------------------------

#' Ordinal rank of a score vector; rank 1 = best.
#' Ties are broken by input order, so the caller's row order is part of the
#' definition (use descending AHP-CDS order everywhere). A missing score is
#' ranked worst (ties among missing scores are again broken by input order).
#' @param score numeric vector
#' @param decreasing TRUE if a larger score is better (default)
rank_scores <- function(score, decreasing = TRUE) {
  stopifnot(is.numeric(score), length(score) > 0)
  key <- if (decreasing) -score else score
  key[is.na(score)] <- Inf
  as.numeric(rank(key, ties.method = "first"))
}

#' Stop unless x is an ordinal permutation of 1..length(x).
assert_ordinal_rank <- function(x, label = deparse(substitute(x))) {
  ok <- !anyNA(x) && anyDuplicated(x) == 0 && all(sort(x) == seq_along(x))
  if (!ok) {
    stop(sprintf(
      "[%s] is not an ordinal permutation of 1..%d (NA, ties or fractional ranks). Build ranks with rank_scores().",
      label, length(x)), call. = FALSE)
  }
  invisible(x)
}

#' Stop unless two rank vectors (same candidate order) are identical.
#' Used to guard that a script reproduces the ranks already exported by another RQ.
assert_same_rank <- function(x, ref, label = "rank") {
  stopifnot(length(x) == length(ref))
  bad <- which(x != ref)
  if (length(bad) > 0) {
    stop(sprintf("[%s] differs from reference at %d position(s), first at index %d (%s vs %s).",
                 label, length(bad), bad[1], x[bad[1]], ref[bad[1]]), call. = FALSE)
  }
  invisible(x)
}

#' Stop unless x is in non-increasing order (input-order tie-break key = AHP-CDS).
assert_sorted_desc <- function(x, label = deparse(substitute(x))) {
  if (is.unsorted(rev(x))) {
    stop(sprintf("[%s] is not sorted in descending order; input order is the tie-break key.", label),
         call. = FALSE)
  }
  invisible(x)
}

# ------------------------------------------------------------------------------
# RBO backends
# ------------------------------------------------------------------------------
#' Manual closed-form RBO for two COMPLETE, equal-length, same-domain rankings.
#' Implements Webber, Moffat & Zobel (2010), formula (32), directly (no
#' extrapolation, since both lists are fully observed over the same domain).
#' Intended as a transparent cross-check against gespeR::rbo(), not as the
#' official reported value unless explicitly requested.
.rbo_exact_complete <- function(rank_a, rank_b, ids, p) {
  stopifnot(length(rank_a) == length(rank_b), length(rank_a) == length(ids))
  n <- length(ids)
  ord_a <- ids[order(rank_a)]
  ord_b <- ids[order(rank_b)]
  
  overlap_at_d <- numeric(n)
  seen_a <- character(0)
  seen_b <- character(0)
  for (d in seq_len(n)) {
    seen_a <- c(seen_a, ord_a[d])
    seen_b <- c(seen_b, ord_b[d])
    overlap_at_d[d] <- length(intersect(seen_a, seen_b)) / d
  }
  (1 - p) * sum(p^(seq_len(n) - 1) * overlap_at_d)
}

#' Official RBO via gespeR::rbo() (Webber et al., 2010).
#' Requires Bioconductor package 'gespeR': BiocManager::install("gespeR")
.rbo_gespeR <- function(rank_a, rank_b, ids, p) {
  if (!requireNamespace("gespeR", quietly = TRUE)) {
    stop(
      "Package 'gespeR' is required for backend = 'gespeR'.\n",
      "Install via: BiocManager::install('gespeR')\n",
      "Or use backend = 'manual' for the formula-equivalent closed-form implementation."
    )
  }
  n <- length(ids)
  # Higher score = better rank (rank 1 -> highest score). Named vectors keyed by id.
  score_a <- stats::setNames(-rank_a, ids)
  score_b <- stats::setNames(-rank_b, ids)
  gespeR::rbo(s = score_a, t = score_b, p = p, k = n, side = "top", uneven.lengths = FALSE)
}

#' Compute RBO between two complete rankings over the same domain.
#'
#' @param rank_a,rank_b numeric rank vectors (same length, same domain, rank 1 = best)
#' @param ids character vector of candidate ids matching rank_a/rank_b positionally
#' @param p RBO persistence parameter in (0,1)
#' @param backend "gespeR" (official, requires Bioconductor) or "manual" (closed-form cross-check)
compute_rbo <- function(rank_a, rank_b, ids, p, backend = c("gespeR", "manual")) {
  backend <- match.arg(backend)
  if (backend == "gespeR") {
    .rbo_gespeR(rank_a, rank_b, ids, p)
  } else {
    .rbo_exact_complete(rank_a, rank_b, ids, p)
  }
}

# ------------------------------------------------------------------------------
# Rank displacement (Delta_r)
# ------------------------------------------------------------------------------

#' Per-candidate rank displacement table.
#' Delta_r_i = r_baseline,i - r_AHP,i  (positive = promoted by AHP-CDS)
compute_delta_r <- function(rank_ahp, rank_baseline, ids) {
  assert_ordinal_rank(rank_ahp, "rank_ahp")
  assert_ordinal_rank(rank_baseline, "rank_baseline")
  data.frame(
    id = ids,
    rank_AHP = rank_ahp,
    rank_baseline = rank_baseline,
    delta_r = rank_baseline - rank_ahp,
    stringsAsFactors = FALSE
  )
}

#' Percentage of candidates whose |Delta_r| exceeds a threshold.
pct_change_gt <- function(delta_r, threshold) {
  mean(abs(delta_r) > threshold) * 100
}

# ------------------------------------------------------------------------------
# Top-k overlap / switching
# ------------------------------------------------------------------------------

#' Top-k overlap between two rankings: shared / AHP-only / baseline-only sets.
topk_overlap <- function(rank_ahp, rank_baseline, ids, k) {
  top_ahp <- ids[rank_ahp <= k]
  top_baseline <- ids[rank_baseline <= k]
  shared <- intersect(top_ahp, top_baseline)
  ahp_only <- setdiff(top_ahp, top_baseline)
  baseline_only <- setdiff(top_baseline, top_ahp)
  list(
    k = k,
    n_shared = length(shared),
    n_ahp_only = length(ahp_only),
    n_baseline_only = length(baseline_only),
    shared = shared,
    ahp_only = ahp_only,
    baseline_only = baseline_only
  )
}

# ------------------------------------------------------------------------------
# Full pairwise summary (the main entry point used by RQ1 and RQ3 scripts)
# ------------------------------------------------------------------------------

#' Summarize one AHP-vs-baseline comparison: correlations, RBO at multiple
#' depths, rank displacement distribution, threshold-based change rates, and
#' top-k switching — everything RQ1/RQ3 need in one call.
#'
#' @param ids character vector of candidate ids (same domain for both rankings)
#' @param rank_ahp,rank_baseline ordinal rank vectors from rank_scores(), rank 1 = best
#' @param pair_label string label for this comparison, e.g. "AHP vs FC-only"
#' @param p_values named numeric vector of RBO persistence parameters to report,
#'   e.g. c(top10 = 0.876, top15 = 0.917)
#' @param delta_thresholds thresholds (rank positions) for pct_change_gt reporting
#' @param k_values depths for top-k overlap / switching tables
#' @param n_extremes number of largest promotions/demotions to report
#' @param rbo_backend "gespeR" (official) or "manual" (closed-form cross-check)
#' @param verify_rbo if TRUE, also compute RBO with the other backend and warn
#'   if the two differ by more than a small numerical tolerance
summarize_pair <- function(ids, rank_ahp, rank_baseline, pair_label,
                           p_values = c(top10 = 0.876, top15 = 0.917),
                           delta_thresholds = c(5, 10, 15, 20),
                           k_values = c(10, 15),
                           n_extremes = 5,
                           rbo_backend = c("gespeR", "manual"),
                           verify_rbo = FALSE) {
  
  rbo_backend <- match.arg(rbo_backend)
  stopifnot(length(ids) == length(rank_ahp), length(ids) == length(rank_baseline))
  
  dr <- compute_delta_r(rank_ahp, rank_baseline, ids)
  
  spearman_rho <- stats::cor(rank_ahp, rank_baseline, method = "spearman")
  kendall_tau  <- stats::cor(rank_ahp, rank_baseline, method = "kendall")
  
  rbo_vals <- sapply(p_values, function(p) {
    val <- compute_rbo(rank_ahp, rank_baseline, ids, p, backend = rbo_backend)
    if (verify_rbo) {
      other <- if (rbo_backend == "gespeR") "manual" else "gespeR"
      other_val <- tryCatch(
        compute_rbo(rank_ahp, rank_baseline, ids, p, backend = other),
        error = function(e) NA_real_
      )
      if (!is.na(other_val) && abs(val - other_val) > 1e-6) {
        warning(sprintf(
          "[%s] RBO backends disagree at p=%.4f: %s=%.6f vs %s=%.6f (diff=%.2e)",
          pair_label, p, rbo_backend, val, other, other_val, abs(val - other_val)
        ))
      }
    }
    val
  })
  names(rbo_vals) <- paste0("RBO_", names(p_values))
  
  pct_vals <- sapply(delta_thresholds, function(t) pct_change_gt(dr$delta_r, t))
  names(pct_vals) <- paste0("pct_change_gt", delta_thresholds)
  
  topk_list <- lapply(k_values, function(k) topk_overlap(rank_ahp, rank_baseline, ids, k))
  names(topk_list) <- paste0("k", k_values)
  
  topk_summary <- unlist(lapply(topk_list, function(x) {
    stats::setNames(
      c(x$n_shared, x$n_ahp_only, x$n_baseline_only),
      paste0("k", x$k, c("_shared", "_ahp_only", "_baseline_only"))
    )
  }))
  
  dr_sorted <- dr[order(-dr$delta_r), ]
  top_promoted <- utils::head(dr_sorted, n_extremes)
  top_demoted  <- utils::tail(dr_sorted, n_extremes)[rev(seq_len(min(n_extremes, nrow(dr_sorted)))), ]
  
  summary_row <- data.frame(
    pair = pair_label,
    n = length(ids),
    spearman = spearman_rho,
    kendall = kendall_tau,
    as.list(rbo_vals),
    delta_r_median = stats::median(dr$delta_r),
    delta_r_Q1 = stats::quantile(dr$delta_r, 0.25, names = FALSE),
    delta_r_Q3 = stats::quantile(dr$delta_r, 0.75, names = FALSE),
    delta_r_min = min(dr$delta_r),
    delta_r_max = max(dr$delta_r),
    as.list(pct_vals),
    as.list(topk_summary),
    check.names = FALSE,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  list(
    summary = summary_row,
    delta_r_full = dr,
    topk_detail = topk_list,
    top_promoted = top_promoted,
    top_demoted = top_demoted
  )
}