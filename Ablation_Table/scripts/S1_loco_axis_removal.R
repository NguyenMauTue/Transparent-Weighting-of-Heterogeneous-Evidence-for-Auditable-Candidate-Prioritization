# ==============================================================================
# loco_axis_removal.R
#
# Leave-one-criterion-out (LOCO) and leave-one-AXIS-out variants of AHP-CDS,
# summarized with the same machinery as RQ1 (ranking_comparison_utils.R).
#
#   criterion removed : -FC, -FDR, -Deg, -Bet
#   axis removed      : -(FC+FDR)   expression axis
#                       -(Deg+Bet)  network axis
#   In every variant the remaining weights are RENORMALIZED to sum to 1.
#
# Input : CDS_candidates_final_result.csv  (sorted by descending CDS;
#         columns UNIPROT, Symbol, ..., n_FC, n_FDR, n_Bet, n_Deg, CDS)
# Output: LOCO_axis_summary.csv  (one row per variant, summarize_pair() format)
#
# Rank convention follows ranking_comparison_utils.R: ordinal ranks, ties broken
# by input order, Delta_r = r_variant - r_AHP (positive = promoted by AHP-CDS).
# ==============================================================================

source(here::here("R", "Helper", "ranking_comparison_utils.R"))
source(here::here("R", "Helper", "ahp_weights.R"))

in_file  <- "CDS_candidates_final_result.csv"
out_file <- "LOCO_axis_summary.csv"

df <- read.csv(in_file, stringsAsFactors = FALSE)
assert_sorted_desc(df$CDS)                       # input order is the tie-break key

# AHP weights: AHP_WEIGHTS (named FC/FDR/Bet/Deg) comes from ahp_weights.R,
# at full precision. Rename to the CSV's criterion columns (n_FC, n_FDR, ...).
W <- AHP_WEIGHTS[c("FC", "FDR", "Bet", "Deg")]
names(W) <- paste0("n_", names(W))
stopifnot(all(is.finite(W)), abs(sum(W) - 1) < 1e-6)

score_without <- function(drop = character(0)) {
  keep <- setdiff(names(W), drop)
  w    <- W[keep] / sum(W[keep])                 # renormalize remaining weights
  as.numeric(as.matrix(df[, keep, drop = FALSE]) %*% w)
}

ids      <- df$UNIPROT
rank_ahp <- rank_scores(df$CDS)

# Guard: the weights must reproduce the exported AHP-CDS ranking exactly
assert_same_rank(rank_scores(score_without()), rank_ahp, "recomputed full CDS")

variants <- list(
  "-FC"         = "n_FC",
  "-FDR"        = "n_FDR",
  "-Deg"        = "n_Deg",
  "-Bet"        = "n_Bet",
  "-(FC+FDR)"   = c("n_FC", "n_FDR"),
  "-(Deg+Bet)"  = c("n_Deg", "n_Bet")
)

res <- lapply(names(variants), function(v) {
  summarize_pair(
    ids           = ids,
    rank_ahp      = rank_ahp,
    rank_baseline = rank_scores(score_without(variants[[v]])),
    pair_label    = paste0("AHP vs ", v),
    p_values      = c(top10 = 0.876, top15 = 0.917),
    rbo_backend   = "manual"                     # closed form; gespeR no longer needed
  )$summary
})
summary_tbl <- do.call(rbind, res)

# Guard: the four single-criterion rows must match the values already in Table S2
expected_rho <- c("AHP vs -FC" = 0.796, "AHP vs -FDR" = 0.981,
                  "AHP vs -Deg" = 0.947, "AHP vs -Bet" = 0.992)
got <- setNames(summary_tbl$spearman, summary_tbl$pair)[names(expected_rho)]
stopifnot(all(abs(got - expected_rho) < 6e-4))

write.csv(summary_tbl, out_file, row.names = FALSE)
print(t(summary_tbl[, setdiff(names(summary_tbl), "n")]))
