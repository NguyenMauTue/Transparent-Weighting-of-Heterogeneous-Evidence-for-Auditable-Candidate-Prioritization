## ============================================================
## R/Helper/cross_dataset_helper.R
##
## SINGLE SOURCE OF TRUTH for:
##   (a) the PXD012162 reference fractional ranking, and
##   (b) the "true" AHP-CDS cross-dataset r (real network, real
##       label) against that reference.
##
## Before this helper, the reference table was re-parsed from
## Module_tables.xlsx independently in 13, 14, 15, 16 (4 copies),
## and the TRUE cross-dataset r for AHP-CDS was computed by THREE
## different code paths that are not provably identical:
##   - 08_cross_dataset_benchmark.R: from ablation_ranks_partial.csv
##     (network stats pre-computed in limma_network_table.csv)
##   - 13_null_distribution_deepdive.R: re-derives frac rank from
##     the same CSV, own cor() call
##   - 15_network_randomization_null.R: REBUILDS the STRING graph
##     from the raw TSV and recomputes degree/betweenness via igraph
##
## Depends on: R/Helper/ahp_weights.R, R/Helper/network_helper.R
## ============================================================

suppressMessages({ library(dplyr); library(openxlsx) })

## ---- (a) PXD012162 reference, fractional-rank by CDS ----
load_pxd012162_reference <- function(module_tables_path) {
  module_tables_pxd012162 <- read.xlsx(module_tables_path, sheet = "All")
  n <- nrow(module_tables_pxd012162)
  module_tables_pxd012162 |>
    arrange(desc(CDS)) |>
    mutate(fracrank_pxd012162 = row_number() / n) |>
    dplyr::select(UNIPROT, fracrank_pxd012162)
}

norm01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

## ---- Score a candidate pool with AHP-CDS given a network centrality table ----
## pool: data.frame with UNIPROT, Symbol, logFC, adj.P.Val (the 61-protein
##       candidate pool, or a permuted-split equivalent)
## centrality: data.frame with Symbol, degree, betweenness (real network,
##       or a rewired/permuted equivalent)
## weights: named vector w["FC"], w["FDR"], w["Bet"], w["Deg"] -- pass
##       AHP_WEIGHTS from ahp_weights.R unless deliberately testing
##       a different weighting scheme.
score_ahp_cds_pool <- function(pool, centrality, weights) {
  d <- pool |> inner_join(centrality, by = "Symbol")
  n <- nrow(d)
  d |>
    mutate(
      n_FC  = norm01(abs(logFC)),
      n_FDR = norm01(-log10(pmax(adj.P.Val, 1e-10))),
      n_Bet = norm01(log1p(betweenness)),
      n_Deg = norm01(log1p(degree)),
      CDS   = weights["FC"]*n_FC + weights["FDR"]*n_FDR +
              weights["Bet"]*n_Bet + weights["Deg"]*n_Deg
    ) |>
    arrange(desc(CDS)) |>
    mutate(frac_rank = row_number() / n)
}

## ---- Cross-dataset Pearson r on fractional ranks vs the 012162 reference ----
compute_cross_dataset_r <- function(scored_pool, reference) {
  m <- inner_join(scored_pool[, c("UNIPROT", "frac_rank")], reference, by = "UNIPROT")
  list(n_overlap = nrow(m), r = cor(m$frac_rank, m$fracrank_pxd012162, method = "pearson"))
}

## ============================================================
## THE canonical number: true AHP-CDS cross-dataset r
## (real network x real label). Compute this ONCE, cache it, and
## have 08/13/15/16 all read from the cache instead of
## recomputing it four different ways.
##
## Uses the TSV-built network (network_helper.R) as the network
## source of truth, since that is the one 15/16's null
## distributions are measured against -- the true value and its
## null distribution must come from the same network-construction
## code path or the empirical p-value in 15/16 is not valid.
## ============================================================
compute_true_ahpcds_cross_dataset_r <- function(pool, string_tsv_path, module_tables_path,
                                                  weights = AHP_WEIGHTS) {
  g_real <- build_real_string_network(string_tsv_path)
  centrality <- network_centrality_table(g_real)
  reference <- load_pxd012162_reference(module_tables_path)

  scored <- score_ahp_cds_pool(pool, centrality, weights)
  result <- compute_cross_dataset_r(scored, reference)

  list(
    r = result$r,
    n_overlap = result$n_overlap,
    n_pool = nrow(scored),
    g_real = g_real,
    scored_pool = scored
  )
}
