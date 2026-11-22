# ==============================================================================
# RQ1_ranking_difference.R
#
# Question: To what extent does AHP-CDS alter candidate prioritization relative
# to expression-only (FC-only), topology-only (Centrality-only), equal-weight
# integration (Equal-weight), and network propagation (RWR)?
#
#
# Inputs:
#   - ablation_ranks_partial.csv : 143-candidate AHP pool with rank_AHPCDS,
#       rank_FConly, rank_Centralityonly, rank_EqualMCDA
#   - rwr_ranks_full.csv         : RWR ranks over the full 291-protein network
#
# Domain handling:
#   - Analysis domain = the 143 AHP-CDS candidates (fixed reference frame).
#   - RWR is restricted to this domain. 11 candidates absent from the RWR output
#     (isolated STRING nodes, no propagation path) have a missing score and are
#     ranked worst by rank_scores(); ties among them are broken by input order
#     (descending AHP-CDS), the same ordinal convention as every other ranking.
#
# Outputs (written to the working directory):
#   - RQ1_summary_table.csv         one row per comparison, all metrics
#   - RQ1_extremes_table.csv        top 5 largest promotions/demotions per pair
#   - RQ1_topk_switching_detail.csv candidate-level top-10/top-15 switching sets
# ==============================================================================
library(here)
source(here::here("R", "Helper","ranking_comparison_utils.R"))
RBO_BACKEND <- "manual"

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
ablation <- read.csv(here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv"), stringsAsFactors = FALSE)
rwr      <- read.csv(here::here("Ablation_Table", "results", "tables" ,"rwr_ranks_full.csv"), stringsAsFactors = FALSE)

stopifnot(nrow(ablation) == 143, !any(duplicated(ablation$Symbol)))
stopifnot(all(c("rank_Centralityonly_ceil", "rank_Centralityonly_floor") %in% names(ablation)))

# ------------------------------------------------------------------------------
# 2. Restrict RWR to the 143-candidate AHP domain; tie-handle missing candidates
# ------------------------------------------------------------------------------
domain_symbols <- ablation$Symbol

rwr_domain <- merge(
  data.frame(Symbol = domain_symbols, stringsAsFactors = FALSE),
  rwr[, c("Symbol", "RWR_score")],
  by = "Symbol", all.x = TRUE
)

n_missing_rwr <- sum(is.na(rwr_domain$RWR_score))
message(sprintf(
  "RWR: %d/%d AHP candidates have no RWR score (isolated STRING nodes) -> ranked worst.",
  n_missing_rwr, nrow(rwr_domain)
))

# merge() sorts alphabetically; restore the AHP-CDS input order BEFORE ranking,
# because input order is the tie-break key.
rwr_domain <- rwr_domain[match(domain_symbols, rwr_domain$Symbol), ]
rwr_domain$rank_RWR_domain <- rank_scores(rwr_domain$RWR_score)

ablation <- merge(ablation, rwr_domain[, c("Symbol", "rank_RWR_domain")], by = "Symbol")
stopifnot(nrow(ablation) == 143)

# Restore a stable row order (merge() re-sorts alphabetically by Symbol)
ablation <- ablation[order(ablation$rank_AHPCDS), ]

# ------------------------------------------------------------------------------
# 3. Run all four comparisons
# ------------------------------------------------------------------------------
# Centrality-only has many exact degree ties, so it is run at both extremes of tie-breaking:
# ceiling (ties follow AHP-CDS order) and floor (ties against it). Both are integer
# permutations, so summarize_pair() is used unchanged.
pairs <- list(
  "AHP vs FC-only"                   = "rank_FConly",
  "AHP vs Centrality-only (ceiling)" = "rank_Centralityonly_ceil",
  "AHP vs Centrality-only (floor)"   = "rank_Centralityonly_floor",
  "AHP vs Equal-weight"              = "rank_EqualMCDA",
  "AHP vs RWR"                       = "rank_RWR_domain"
)

results <- lapply(names(pairs), function(label) {
  summarize_pair(
    ids           = ablation$Symbol,
    rank_ahp      = ablation$rank_AHPCDS,
    rank_baseline = ablation[[pairs[[label]]]],
    pair_label    = label,
    p_values      = c(top10 = 0.876, top15 = 0.917),   # ~90% RBO weight at d=10, d=15 (Webber et al. 2010, formula 32)
    delta_thresholds = c(5, 10, 15, 20),
    k_values      = c(10, 15),
    n_extremes    = 5,
    rbo_backend   = RBO_BACKEND,
    verify_rbo    = FALSE   # set TRUE only if gespeR becomes available again, to cross-check against manual
  )
})
names(results) <- names(pairs)

# ------------------------------------------------------------------------------
# 4. Export
# ------------------------------------------------------------------------------
summary_table <- do.call(rbind, lapply(results, function(x) x$summary))

# Expected value over random tie-breaking for Centrality-only. Only the rank correlations
# are linear in the tie-break (exact mean of the two extremes). RBO, top-k overlap and
# delta-r counts are NOT, so they are left NA here: report them as the [floor, ceiling]
# range from the two rows above.
ceil_row  <- summary_table[summary_table$pair == "AHP vs Centrality-only (ceiling)", ]
floor_row <- summary_table[summary_table$pair == "AHP vs Centrality-only (floor)", ]
exp_row <- ceil_row
exp_row$pair <- "AHP vs Centrality-only (expected over tie-breaks)"
num_cols <- setdiff(names(exp_row)[vapply(exp_row, is.numeric, logical(1))], "n")
exp_row[, num_cols] <- NA_real_
for (cc in intersect(c("spearman", "kendall"), names(exp_row))) {
  exp_row[[cc]] <- (ceil_row[[cc]] + floor_row[[cc]]) / 2
}
summary_table <- rbind(summary_table, exp_row)
write.csv(summary_table, here::here("Ablation_Table", "results", "tables", "RQ1_summary_table.csv"), row.names = FALSE)

extremes_table <- do.call(rbind, lapply(names(results), function(label) {
  rbind(
    cbind(pair = label, direction = "promoted", results[[label]]$top_promoted),
    cbind(pair = label, direction = "demoted",  results[[label]]$top_demoted)
  )
}))
write.csv(extremes_table, here::here("Ablation_Table", "results", "tables", "RQ1_extremes_table.csv"), row.names = FALSE)

topk_detail_table <- do.call(rbind, lapply(names(results), function(label) {
  do.call(rbind, lapply(results[[label]]$topk_detail, function(tk) {
    data.frame(
      pair = label,
      k = tk$k,
      n_shared = tk$n_shared,
      n_ahp_only = tk$n_ahp_only,
      n_baseline_only = tk$n_baseline_only,
      ahp_only = paste(tk$ahp_only, collapse = "; "),
      baseline_only = paste(tk$baseline_only, collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(topk_detail_table, here::here("Ablation_Table", "results", "tables", "RQ1_topk_switching_detail.csv"), row.names = FALSE)

# Per-candidate displacement vs Centrality-only at the three tie-break conventions
# (delta_r = r_baseline - r_AHP, positive = promoted by AHP-CDS). tie_block_size = number
# of candidates sharing the same degree; 1 means the rank is not affected by tie-breaking.
dr_cen <- data.frame(
  Symbol = ablation$Symbol,
  rank_AHPCDS = ablation$rank_AHPCDS,
  tie_block_size = as.numeric(ave(ablation$score_Centralityonly, ablation$score_Centralityonly, FUN = length)),
  rank_Centralityonly_ceil  = ablation$rank_Centralityonly_ceil,
  rank_Centralityonly       = ablation$rank_Centralityonly,
  rank_Centralityonly_floor = ablation$rank_Centralityonly_floor,
  stringsAsFactors = FALSE
)
dr_cen$delta_r_ceil  <- dr_cen$rank_Centralityonly_ceil  - dr_cen$rank_AHPCDS
dr_cen$delta_r_mid   <- dr_cen$rank_Centralityonly       - dr_cen$rank_AHPCDS
dr_cen$delta_r_floor <- dr_cen$rank_Centralityonly_floor - dr_cen$rank_AHPCDS
write.csv(dr_cen, here::here("Ablation_Table", "results", "tables", "RQ1_centrality_tiebreak_delta_r.csv"), row.names = FALSE)

message("Done.\n  - RQ1_summary_table.csv\n  - RQ1_extremes_table.csv\n  - RQ1_topk_switching_detail.csv\n  - RQ1_centrality_tiebreak_delta_r.csv")