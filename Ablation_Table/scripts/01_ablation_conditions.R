library(dplyr)

data056 <- read.csv(here::here("PXD056161", "results", "tables", "CDS_candidates_final_result.csv"))
#Same normalization method as AHP-CDS
norm01 <- function(x) (x - min(x)) / (max(x) - min(x))


data056 <- data056 |>
  mutate(
    n_FC  = norm01(abs(logFC)),
    n_FDR = norm01(-log10(pmax(adj.P.Val, 1e-10))),
    n_Bet = norm01(log1p(betweenness)),
    n_Deg = norm01(log1p(degree))
  )
## ---- Condition: FC-only ----
# rank purely by |logFC| magnitude (direction-agnostic prioritization, matches
# how AHP-CDS treats FC via n_FC = norm01(abs(logFC)))
fc_only <- data056 |>
  mutate(score_FConly = n_FC) |>
  arrange(desc(score_FConly)) |>
  mutate(rank_FConly = row_number()) |>
  dplyr::select(UNIPROT, score_FConly, rank_FConly)

## ---- Condition: Centrality-only ----
# no expression signal at all
centrality_only <- {
  # Degree has many exact ties (few distinct integer values), so the score alone does
  # not define one ordinal rank. Three versions are kept:
  #   rank_Centralityonly_ceil  : ties ordered like AHP-CDS (descending CDS)
  #                               = the ranks reported so far; upper bound of agreement
  #   rank_Centralityonly_floor : ties ordered against AHP-CDS (ascending CDS); lower bound
  #   rank_Centralityonly       : mean of the two = average rank of each tie block
  #                               = expected rank under random tie-breaking
  stopifnot("input must be sorted by descending CDS so that the ceiling reproduces earlier ranks" =
              !is.unsorted(-data056$CDS))
  cen <- data056 |> mutate(score_Centralityonly = n_Deg)
  ord_ceil  <- order(-cen$score_Centralityonly, -cen$CDS)
  ord_floor <- order(-cen$score_Centralityonly,  cen$CDS)
  cen$rank_Centralityonly_ceil  <- match(seq_len(nrow(cen)), ord_ceil)
  cen$rank_Centralityonly_floor <- match(seq_len(nrow(cen)), ord_floor)
  cen$rank_Centralityonly <- (cen$rank_Centralityonly_ceil + cen$rank_Centralityonly_floor) / 2
  stopifnot(isTRUE(all.equal(cen$rank_Centralityonly,
                             rank(-cen$score_Centralityonly, ties.method = "average"))))
  cen |>
    dplyr::select(UNIPROT, score_Centralityonly, rank_Centralityonly,
                  rank_Centralityonly_ceil, rank_Centralityonly_floor)
}

## ---- Condition: Equal-weight MCDA ----
# same 4 criteria as AHP-CDS (FC, FDR, Bet, Deg) but weights = 1/4 each
# instead of AHP-derived weights -> isolates the value AHP weighting adds
# on top of "just use all 4 criteria"
equal_mcda <- data056 |>
  mutate(score_EqualMCDA = 0.25*n_FC + 0.25*n_FDR + 0.25*n_Bet + 0.25*n_Deg) |>
  arrange(desc(score_EqualMCDA)) |>
  mutate(rank_EqualMCDA = row_number()) |>
  dplyr::select(UNIPROT, score_EqualMCDA, rank_EqualMCDA)

## ---- Reference: AHP-CDS (full), already computed in CDS_candidates_annotated.csv ----
ahp_cds <- data056 |>
  distinct(UNIPROT, .keep_all = TRUE) |>
  arrange(desc(CDS)) |>
  mutate(rank_AHPCDS = row_number()) |>
  dplyr::select(UNIPROT, Symbol, CDS, rank_AHPCDS)

## ---- Merge all into one comparison table ----
ablation_ranks <- ahp_cds |>
  full_join(fc_only, by = "UNIPROT") |>
  full_join(centrality_only, by = "UNIPROT") |>
  full_join(equal_mcda, by = "UNIPROT") |>
  arrange(rank_AHPCDS)

write.csv(ablation_ranks, here::here("Ablation_Table", "results", "tables", "ablation_ranks_partial.csv"), row.names = FALSE)

cat("N candidates in pool:", nrow(data056), "\n")
cat("N in AHP-CDS table:", nrow(ahp_cds), "\n\n")
cat("=== Top 10 by AHP-CDS (full) ===\n")
print(head(ablation_ranks[,c("Symbol","rank_AHPCDS","rank_FConly","rank_Centralityonly","rank_EqualMCDA")], 10))

cat("\n=== Rank agreement (Spearman) between AHP-CDS and each ablated condition ===\n")
cat("AHP-CDS vs FC-only:         ", cor(ablation_ranks$rank_AHPCDS, ablation_ranks$rank_FConly, method="spearman", use="pairwise.complete.obs"), "\n")
# Spearman is linear in the tie-break, so the expected value over random tie-breaking is
# exactly the mean of the ceiling and floor values. (cor() on the average-rank column
# rank_Centralityonly would give the tie-corrected Spearman, which is NOT the same number.)
rho_cen_ceil  <- cor(ablation_ranks$rank_AHPCDS, ablation_ranks$rank_Centralityonly_ceil,  method="spearman", use="pairwise.complete.obs")
rho_cen_floor <- cor(ablation_ranks$rank_AHPCDS, ablation_ranks$rank_Centralityonly_floor, method="spearman", use="pairwise.complete.obs")
cat("AHP-CDS vs Centrality-only (ceiling: ties follow CDS):   ", rho_cen_ceil, "\n")
cat("AHP-CDS vs Centrality-only (floor: ties against CDS):    ", rho_cen_floor, "\n")
cat("AHP-CDS vs Centrality-only (expected = mean of the two): ", (rho_cen_ceil + rho_cen_floor) / 2, "\n")
cat("AHP-CDS vs Equal-weight MCDA:", cor(ablation_ranks$rank_AHPCDS, ablation_ranks$rank_EqualMCDA, method="spearman", use="pairwise.complete.obs"), "\n")

