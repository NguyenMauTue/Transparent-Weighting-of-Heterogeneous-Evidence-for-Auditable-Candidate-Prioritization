## ============================================================
## R/Helper/network_helper.R
##
## SINGLE SOURCE OF TRUTH for building the "real" PPI graph from
## the STRING TSV export. Previously this exact block (read.delim
## -> rename node1 -> filter combined_score > 0.7 -> graph_from_
## data_frame -> simplify) was duplicated verbatim in
## 15_network_randomization_null.R and 16_double_null_2x2.R.
##
## Consolidating this matters for more than tidiness: 08 and 13
## derive degree/betweenness from a *pre-computed* CSV
## (limma_network_table.csv), while 15/16 rebuild the graph from
## the raw TSV and recompute degree/betweenness with igraph. If
## those two paths are not provably equivalent (same score
## threshold, same simplify() behaviour, same node set), the
## "true cross-dataset r for AHP-CDS" reported by 08/13 and the
## rewire_id==0 row in 15/16's output can silently disagree.
## This helper is the fix for that: everything downstream should
## call build_real_string_network() rather than re-deriving it.
## ============================================================

suppressMessages(library(igraph))

STRING_SCORE_THRESHOLD <- 0.7  # matches the value used throughout PXD056161 scripts

## Build the real (non-permuted) PPI graph from a STRING interactions
## TSV export. Returns an igraph object with degree/betweenness
## computable directly via igraph::degree()/igraph::betweenness().
build_real_string_network <- function(tsv_path, score_threshold = STRING_SCORE_THRESHOLD) {
  edges <- read.delim(tsv_path, check.names = FALSE)
  colnames(edges)[1] <- "node1"
  edges_hi <- edges[edges$combined_score > score_threshold, ]
  g <- graph_from_data_frame(edges_hi[, c("node1", "node2")], directed = FALSE)
  g <- simplify(g)
  g
}

## Convenience: degree + betweenness as a tidy data.frame keyed by Symbol,
## in the same shape both 15 and 16 need for the CDS scoring join.
network_centrality_table <- function(g) {
  deg <- igraph::degree(g)
  bet <- igraph::betweenness(g)
  data.frame(Symbol = names(deg), degree = deg, betweenness = bet[names(deg)],
             row.names = NULL)
}

## Diagnostic: compare this graph's centrality table against an
## already-computed limma_network_table.csv (the source 08/13 trust).
## Call this once, on the real network, to confirm both code paths
## agree before relying on either downstream. Returns a data.frame
## of proteins where degree or betweenness disagree beyond `tol`.
compare_network_sources <- function(g_from_tsv, limma_network_table_path, tol = 1e-6) {
  built <- network_centrality_table(g_from_tsv)
  precomputed <- read.csv(limma_network_table_path) |>
    dplyr::distinct(Symbol, .keep_all = TRUE) |>
    dplyr::select(Symbol, degree, betweenness)

  merged <- merge(built, precomputed, by = "Symbol", suffixes = c("_tsv", "_csv"))
  merged$degree_diff <- abs(merged$degree_tsv - merged$degree_csv)
  merged$betweenness_diff <- abs(merged$betweenness_tsv - merged$betweenness_csv)

  mismatches <- merged[merged$degree_diff > tol | merged$betweenness_diff > tol, ]

  cat(sprintf(
    "Network source comparison: %d proteins in common, %d mismatched (tol=%g)\n",
    nrow(merged), nrow(mismatches), tol
  ))
  if (nrow(mismatches) > 0) {
    cat("First few mismatches:\n")
    print(head(mismatches[order(-mismatches$degree_diff), ], 10))
  }
  invisible(mismatches)
}
