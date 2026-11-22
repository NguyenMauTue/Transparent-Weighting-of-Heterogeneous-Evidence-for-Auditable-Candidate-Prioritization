## ==============================================================================
## RQ2_mechanism_analysis.R
##
## RQ2 -- Which candidates move, and why?
##   Part A: weighted criterion contribution shares  (C_ij = w_j * x_ij, P_ij)
##   Part B: LOCO (Leave-One-Criterion-Out) (Dogan et al. 2026)
##
##
## SCOPE NOTE ON RWR:
##   RWR is intentionally excluded from every comparison in this script.
##   RWR is a random-walk-with-restart score, not a linear-additive
##   combination of {FC, FDR, Degree, Betweenness}. It has no per-criterion
##   weight to leave out, so it cannot participate in LOCO, and it was not
##   used to build C_ij/P_ij in Part A either.
##
## BASELINES FOR RQ2A GROUPING:
##   - FC-only        : rank by x_FC (equivalently |logFC|) descending
##   - Centrality-only : rank by raw degree descending
##   Two independent promoted/stable/demoted groupings are produced, one
##   per baseline, using threshold T (10 and 15; see note at top of Part A).
##
## LOCO SCOPE:
##   LOCO comparisons are internal only: full AHP-CDS vs each
##   leave-one-criterion-out variant.
## ==============================================================================
library(here)
library(readxl)

suppressWarnings(suppressMessages({
  source(here::here("R", "Helper","ranking_comparison_utils.R"))
  source(here::here("R", "Helper","ahp_weights.R"))
}))

set.seed(24082026) 
# ------------------------------------------------------------------------------
# Part 0: Load data, reconstruct x_ij, validate against existing CDS column
# ------------------------------------------------------------------------------

dat <- protein_group <- read.csv(here::here("PXD056161", "results", "tables", "CDS_candidates_final_result.csv"))

stopifnot(nrow(dat) == 143)

minmax <- function(x) (x - min(x)) / (max(x) - min(x))

x_FC  <- minmax(abs(dat$logFC))
x_FDR <- minmax(-log10(dat$adj.P.Val))
x_Deg <- minmax(log1p(dat$degree))
x_Bet <- minmax(log1p(dat$betweenness))

X <- data.frame(
  id    = dat$UNIPROT,
  symbol = dat$Symbol,
  FC    = x_FC,
  FDR   = x_FDR,
  Deg   = x_Deg,
  Bet   = x_Bet,
  stringsAsFactors = FALSE
)

w <- AHP_WEIGHTS[c("FC", "FDR", "Deg", "Bet")]

cds_reconstructed <- as.numeric(as.matrix(X[, c("FC", "FDR", "Deg", "Bet")]) %*% w)
max_abs_err <- max(abs(cds_reconstructed - dat$CDS))
cat(sprintf("[validation] max |CDS_reconstructed - CDS_original| = %.3e\n", max_abs_err))
stopifnot(max_abs_err < 1e-6)  # hard stop if x_ij reconstruction is ever wrong

dat$CDS_check <- cds_reconstructed  # kept for audit trail in output, not used downstream

ids <- dat$UNIPROT
# Rank convention = RQ1 (see ranking_comparison_utils.R): ordinal ranks, ties
# broken by input order, so dat must be in descending CDS order.
assert_sorted_desc(dat$CDS, "dat$CDS")
rank_AHP <- rank_scores(dat$CDS)

# ------------------------------------------------------------------------------
# Part A: Weighted criterion contribution shares  C_ij = w_j * x_ij
# ------------------------------------------------------------------------------

C <- sweep(as.matrix(X[, c("FC", "FDR", "Deg", "Bet")]), 2, w, "*")
colnames(C) <- paste0("C_", colnames(C))

# P_ij = C_ij / CDS_i  (sum_k C_ik = CDS_i exactly, by construction of CDS
# as a linear-additive AHP score -- so shares sum to 1 per candidate)
P <- sweep(C, 1, dat$CDS, "/")
colnames(P) <- paste0("P_", sub("^C_", "", colnames(C)))

contribution_table <- data.frame(
  UNIPROT = dat$UNIPROT,
  Symbol  = dat$Symbol,
  CDS     = dat$CDS,
  rank_AHP = rank_AHP,
  C, P,
  stringsAsFactors = FALSE
)

# --- Baseline rankings for grouping ---
rank_FC_only         <- rank_scores(x_FC)
# Degree has many exact ties. rank_scores() breaks them by input order (descending CDS),
# which is the CEILING of agreement with AHP-CDS. The FLOOR reverses the tie-break, and the
# MIDPOINT (= average rank of each tie block = expected rank under random tie-breaking) is
# what the grouping below uses. rank_Centrality_only keeps its name but is now the midpoint.
rank_Centrality_ceil  <- rank_scores(dat$degree)
rank_Centrality_floor <- rev(rank_scores(rev(dat$degree)))
rank_Centrality_only  <- (rank_Centrality_ceil + rank_Centrality_floor) / 2

# Guard: RQ2 must reproduce the ranks RQ1 already exported.
ablation_ref <- read.csv(here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv"),
                         stringsAsFactors = FALSE)
m_ref <- match(dat$Symbol, ablation_ref$Symbol)
stopifnot(!anyNA(m_ref))
assert_same_rank(rank_AHP,             ablation_ref$rank_AHPCDS[m_ref],          "rank_AHP vs RQ1")
assert_same_rank(rank_FC_only,         ablation_ref$rank_FConly[m_ref],          "rank_FC_only vs RQ1")
assert_same_rank(rank_Centrality_ceil,  ablation_ref$rank_Centralityonly_ceil[m_ref],  "rank_Centrality_ceil vs RQ1")
assert_same_rank(rank_Centrality_floor, ablation_ref$rank_Centralityonly_floor[m_ref], "rank_Centrality_floor vs RQ1")
stopifnot(isTRUE(all.equal(rank_Centrality_only, ablation_ref$rank_Centralityonly[m_ref])))

dr_FC  <- compute_delta_r(rank_AHP, rank_FC_only, ids)
dr_Cen_ceil  <- compute_delta_r(rank_AHP, rank_Centrality_ceil,  ids)
dr_Cen_floor <- compute_delta_r(rank_AHP, rank_Centrality_floor, ids)
# Midpoint (expected) displacement: the midpoint rank is not an integer permutation, so
# delta_r = r_baseline - r_AHP is computed directly. The check confirms this is the same
# sign convention as compute_delta_r().
stopifnot(isTRUE(all.equal(as.numeric(dr_Cen_ceil$delta_r), as.numeric(rank_Centrality_ceil - rank_AHP))))
dr_Cen <- list(delta_r = rank_Centrality_only - rank_AHP)

# --- Grouping thresholds T ---
# T is a rank-displacement threshold, in rank positions, applied to
#   delta_r = r_baseline - r_AHP   (positive = promoted by AHP-CDS):
#   promoted if delta_r >= T,  demoted if delta_r <= -T,  stable otherwise.
# T is NOT the top-k depth k used in RQ1/RQ5.
# T = 10 and 15 are both run and reported regardless of outcome
# (sensitivity/robustness check, not threshold-shopping): if the contribution-
# share pattern holds at both, that's evidence the mechanism isn't threshold-
# dependent; if it doesn't, that's also reported as-is.
T_thresholds <- c(10, 15)

group_by_delta <- function(delta_r, T) {
  ifelse(delta_r >= T, "promoted",
         ifelse(delta_r <= -T, "demoted", "stable"))
}

summarize_profile <- function(tbl, group_col) {
  agg <- aggregate(
    tbl[, c("P_FC", "P_FDR", "P_Deg", "P_Bet")],
    by = list(group = tbl[[group_col]]),
    FUN = mean
  )
  agg$n <- as.numeric(table(tbl[[group_col]])[agg$group])
  agg[, c("group", "n", "P_FC", "P_FDR", "P_Deg", "P_Bet")]
}

profiles_by_T <- list()   # profiles_by_T[[paste0("T", T)]][["FC_only" / "Centrality_only"]]

for (T in T_thresholds) {
  tag <- paste0("T", T)
  
  contribution_table[[paste0("group_vs_FC_only_", tag)]] <-
    group_by_delta(dr_FC$delta_r, T)
  contribution_table[[paste0("group_vs_Centrality_only_", tag)]] <-
    group_by_delta(dr_Cen$delta_r, T)
  contribution_table[[paste0("group_vs_Centrality_only_", tag, "_ceil")]] <-
    group_by_delta(dr_Cen_ceil$delta_r, T)
  contribution_table[[paste0("group_vs_Centrality_only_", tag, "_floor")]] <-
    group_by_delta(dr_Cen_floor$delta_r, T)
  contribution_table[[paste0("group_same_at_both_tiebreaks_Centrality_only_", tag)]] <-
    contribution_table[[paste0("group_vs_Centrality_only_", tag, "_ceil")]] ==
    contribution_table[[paste0("group_vs_Centrality_only_", tag, "_floor")]]
  
  cat("\n[Part A] Group sizes (T = ", T, ")\n", sep = "")
  cat("  vs FC-only:        ")
  print(table(contribution_table[[paste0("group_vs_FC_only_", tag)]]))
  cat("  vs Centrality-only: ")
  print(table(contribution_table[[paste0("group_vs_Centrality_only_", tag)]]))
  
  profile_vs_FC  <- summarize_profile(contribution_table, paste0("group_vs_FC_only_", tag))
  profile_vs_Cen <- summarize_profile(contribution_table, paste0("group_vs_Centrality_only_", tag))
  profile_vs_Cen_ceil  <- summarize_profile(contribution_table, paste0("group_vs_Centrality_only_", tag, "_ceil"))
  profile_vs_Cen_floor <- summarize_profile(contribution_table, paste0("group_vs_Centrality_only_", tag, "_floor"))
  
  cat("\n[Part A] Mean contribution-share profile (T=", T, "), grouped vs FC-only:\n", sep = "")
  print(profile_vs_FC, digits = 3)
  cat("\n[Part A] Mean contribution-share profile (T=", T, "), grouped vs Centrality-only:\n", sep = "")
  print(profile_vs_Cen, digits = 3)
  
  profiles_by_T[[tag]] <- list(FC_only = profile_vs_FC, Centrality_only = profile_vs_Cen,
                               Centrality_only_ceil = profile_vs_Cen_ceil,
                               Centrality_only_floor = profile_vs_Cen_floor)
}

# delta_r columns are threshold-independent (T only affects grouping), so
# store once
contribution_table$delta_r_vs_FC_only         <- dr_FC$delta_r
contribution_table$delta_r_vs_Centrality_only <- dr_Cen$delta_r
contribution_table$delta_r_vs_Centrality_only_ceil  <- dr_Cen_ceil$delta_r
contribution_table$delta_r_vs_Centrality_only_floor <- dr_Cen_floor$delta_r

# --- Cross-threshold consistency check ---
# For each baseline, does the "promoted" group's mean P_Deg (the criterion
# that mechanistically distinguishes AHP-CDS from FC-only, and vice versa
# for Centrality-only) point the same direction at T=10 and T=15?
consistency_check <- function(baseline_name, criterion_col) {
  p10 <- profiles_by_T[["T10"]][[baseline_name]]
  p15 <- profiles_by_T[["T15"]][[baseline_name]]
  data.frame(
    baseline = baseline_name,
    criterion = criterion_col,
    promoted_T10  = p10[p10$group == "promoted", criterion_col],
    demoted_T10   = p10[p10$group == "demoted", criterion_col],
    promoted_T15  = p15[p15$group == "promoted", criterion_col],
    demoted_T15   = p15[p15$group == "demoted", criterion_col],
    direction_consistent = sign(p10[p10$group == "promoted", criterion_col] -
                                  p10[p10$group == "demoted", criterion_col]) ==
      sign(p15[p15$group == "promoted", criterion_col] -
             p15[p15$group == "demoted", criterion_col])
  )
}

consistency_table <- rbind(
  consistency_check("FC_only", "P_Deg"),
  consistency_check("Centrality_only", "P_FC")
)

cat("\n[Part A] Cross-threshold (T=10 vs T=15) consistency check:\n")
print(consistency_table, digits = 3, row.names = FALSE)

# --- Tie-break sensitivity for the Centrality-only grouping ---
# Promoted/demoted labels vs Centrality-only are computed at the midpoint rank (used above)
# and at both extremes of tie-breaking (ceiling: ties follow AHP-CDS order; floor: against
# it). The pattern "promoted-vs-Centrality-only <-> higher P_FC" is robust if P_FC_diff keeps
# its sign (and rough size) at all three.
tiebreak_variants <- c(ceiling = "Centrality_only_ceil", midpoint = "Centrality_only", floor = "Centrality_only_floor")
tiebreak_sensitivity <- do.call(rbind, lapply(T_thresholds, function(T) {
  tag <- paste0("T", T)
  do.call(rbind, lapply(names(tiebreak_variants), function(v) {
    p <- profiles_by_T[[tag]][[tiebreak_variants[[v]]]]
    pick <- function(g, col) { x <- p[p$group == g, col]; if (length(x) == 0) NA_real_ else x }
    data.frame(T = T, tiebreak = v,
               n_promoted = pick("promoted", "n"), n_demoted = pick("demoted", "n"),
               promoted_P_FC = pick("promoted", "P_FC"), demoted_P_FC = pick("demoted", "P_FC"),
               P_FC_diff = pick("promoted", "P_FC") - pick("demoted", "P_FC"),
               stringsAsFactors = FALSE)
  }))
}))
cat("\n[Part A] Tie-break sensitivity (Centrality-only grouping):\n")
print(tiebreak_sensitivity, digits = 3, row.names = FALSE)

# ------------------------------------------------------------------------------
# Part B: LOCO -- Leave-One-Criterion-Out (internal only, vs full AHP-CDS)
# ------------------------------------------------------------------------------

criteria <- c("FC", "FDR", "Deg", "Bet")

loco_results <- list()
loco_summary_rows <- list()

for (crit in criteria) {
  remaining <- setdiff(criteria, crit)
  w_loco <- w[remaining] / sum(w[remaining])  # renormalize
  
  cds_loco <- as.numeric(as.matrix(X[, remaining]) %*% w_loco)
  rank_loco <- rank_scores(cds_loco)
  
  pair_label <- paste0("Full-AHP vs LOCO-without-", crit)
  res <- summarize_pair(
    ids = ids,
    rank_ahp = rank_AHP,
    rank_baseline = rank_loco,
    pair_label = pair_label,
    p_values = c(top10 = 0.876, top15 = 0.917),
    delta_thresholds = c(5, 10, 15, 20),
    k_values = c(10, 15),
    n_extremes = 5,
    rbo_backend = "manual"   # gespeR unavailable in this environment; see
    # ranking_comparison_utils.R header note
  )
  
  loco_results[[crit]] <- res
  loco_summary_rows[[crit]] <- res$summary
}

loco_summary_table <- do.call(rbind, loco_summary_rows)
rownames(loco_summary_table) <- NULL

cat("\n[Part B] LOCO summary (full-AHP vs each leave-one-out variant):\n")
# Note: k10/k15 shared-count columns come out of summarize_pair() named
# "k10.k10_shared" / "k15.k15_shared" etc. (unlist() outer.inner naming,
# from the already-tested utils -- left as-is rather than modified here).
print(loco_summary_table[, c("pair", "n", "spearman", "kendall",
                             "RBO_top10", "RBO_top15",
                             "delta_r_median", "pct_change_gt15",
                             "k10.k10_shared", "k15.k15_shared")], digits = 3)

# Per-candidate LOCO delta_r, wide table: one delta_r column per left-out criterion
loco_delta_wide <- data.frame(UNIPROT = ids, Symbol = dat$Symbol, rank_AHP = rank_AHP)
for (crit in criteria) {
  loco_delta_wide[[paste0("delta_r_without_", crit)]] <- loco_results[[crit]]$delta_r_full$delta_r
}

# "Which candidates depend strongly on each criterion" = large |delta_r| when
# that criterion is removed. Flag candidates whose |delta_r| >= T_threshold
# for at least one criterion, and report which criterion(s) drive them most.
loco_delta_cols <- paste0("delta_r_without_", criteria)
loco_delta_wide$max_abs_delta_r <- apply(abs(loco_delta_wide[, loco_delta_cols]), 1, max)
loco_delta_wide$most_dependent_on <- criteria[apply(abs(loco_delta_wide[, loco_delta_cols]), 1, which.max)]
loco_delta_wide <- loco_delta_wide[order(-loco_delta_wide$max_abs_delta_r), ]

cat("\n[Part B] Top 10 candidates most sensitive to a single criterion's removal:\n")
print(head(loco_delta_wide[, c("Symbol", "rank_AHP", loco_delta_cols,
                               "most_dependent_on")], 10), digits = 3, row.names = FALSE)

# ------------------------------------------------------------------------------
# Exports
# ------------------------------------------------------------------------------
write.csv(contribution_table, here::here("Ablation_Table", "results", "tables", "RQ2A_contribution_shares.csv"), row.names = FALSE)
for (T in T_thresholds) {
  tag <- paste0("T", T)
  write.csv(profiles_by_T[[tag]]$FC_only,
            sprintf(here::here("Ablation_Table", "results", "tables", "RQ2A_profile_vs_FC_only_T%d.csv"), T), row.names = FALSE)
  write.csv(profiles_by_T[[tag]]$Centrality_only,
            sprintf(here::here("Ablation_Table", "results", "tables", "RQ2A_profile_vs_Centrality_only_T%d.csv"), T), row.names = FALSE)
  write.csv(profiles_by_T[[tag]]$Centrality_only_ceil,
            sprintf(here::here("Ablation_Table", "results", "tables", "RQ2A_profile_vs_Centrality_only_ceil_T%d.csv"), T), row.names = FALSE)
  write.csv(profiles_by_T[[tag]]$Centrality_only_floor,
            sprintf(here::here("Ablation_Table", "results", "tables", "RQ2A_profile_vs_Centrality_only_floor_T%d.csv"), T), row.names = FALSE)
}
write.csv(consistency_table, here::here("Ablation_Table", "results", "tables", "RQ2A_cross_threshold_consistency.csv"), row.names = FALSE)
write.csv(tiebreak_sensitivity, here::here("Ablation_Table", "results", "tables", "RQ2A_centrality_tiebreak_sensitivity.csv"), row.names = FALSE)

write.csv(loco_summary_table, here::here("Ablation_Table", "results", "tables", "RQ2B_LOCO_summary_table.csv"), row.names = FALSE)
write.csv(loco_delta_wide, here::here("Ablation_Table", "results", "tables", "RQ2B_LOCO_delta_r_wide.csv"), row.names = FALSE)

for (crit in criteria) {
  write.csv(
    loco_results[[crit]]$delta_r_full,
    sprintf(here::here("Ablation_Table", "results", "tables", "RQ2B_LOCO_without_%s_delta_r_full.csv"), crit),
    row.names = FALSE
  )
}

cat("\nDone. Outputs written to", here::here("Ablation_Table", "results", "tables"))