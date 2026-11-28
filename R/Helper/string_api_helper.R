############################################################
# STRING REST API helper — replaces manual website upload
#
# CHANGE from original: added `return_graph` param. Default FALSE
# preserves the exact original behavior (data.frame of Symbol/degree/
# betweenness) for the null-permutation-test use case. Set TRUE to
# also get back the igraph object + raw edge table for plotting.
############################################################

library(httr)
library(jsonlite)

#' Build a STRING PPI network for a set of gene symbols and return
#' a degree/betweenness summary table (or, with return_graph=TRUE,
#' also the graph itself for plotting).
#'
#' @param gene_symbols character vector of gene symbols
#' @param species NCBI taxon id, 9606 = homo sapiens
#' @param required_score STRING combined_score threshold, 0-1000 scale
#'        (700 == the 0.7 "high confidence" threshold used in script 10)
#' @param caller_identity required by STRING's fair-use policy — put a
#'        short project identifier here (e.g. your email or repo name)
#' @param return_graph if TRUE, returns list(summary, graph, edges)
#'        instead of just the summary data.frame — use this for the
#'        network figure panel, leave FALSE everywhere else unchanged.
build_string_network <- function(gene_symbols,
                                 species = 9606,
                                 required_score = 700,
                                 caller_identity = "AHP_CDS_null_permutation_test",
                                 return_graph = FALSE) {
  
  gene_symbols <- unique(gene_symbols[!is.na(gene_symbols) & gene_symbols != ""])
  stopifnot(length(gene_symbols) >= 2)
  
  # --- Step 1: map symbols -> STRING IDs ---
  id_resp <- httr::POST(
    url = "https://string-db.org/api/tsv/get_string_ids",
    body = list(
      identifiers = paste(gene_symbols, collapse = "\r"),
      species = species,
      limit = 1,
      echo_query = 1,
      caller_identity = caller_identity
    ),
    encode = "form"
  )
  httr::stop_for_status(id_resp)
  id_map <- read.delim(text = httr::content(id_resp, as = "text", encoding = "UTF-8"),
                       stringsAsFactors = FALSE)
  
  empty_summary <- data.frame(Symbol = character(0), degree = numeric(0), betweenness = numeric(0))
  
  if (nrow(id_map) == 0) {
    warning("No STRING IDs resolved for this gene set.")
    if (return_graph) return(list(summary = empty_summary, graph = NULL, edges = NULL))
    return(empty_summary)
  }
  
  string_ids <- unique(id_map$stringId)
  
  # --- Step 2: pull the high-confidence network for those IDs ---
  net_resp <- httr::POST(
    url = "https://string-db.org/api/tsv/network",
    body = list(
      identifiers = paste(string_ids, collapse = "\r"),
      species = species,
      required_score = required_score,
      caller_identity = caller_identity
    ),
    encode = "form"
  )
  httr::stop_for_status(net_resp)
  edges <- read.delim(text = httr::content(net_resp, as = "text", encoding = "UTF-8"),
                      stringsAsFactors = FALSE)
  
  if (nrow(edges) == 0) {
    warning("STRING returned zero edges above required_score for this gene set.")
    fallback_summary <- data.frame(Symbol = gene_symbols, degree = 0, betweenness = 0)
    if (return_graph) return(list(summary = fallback_summary, graph = NULL, edges = NULL))
    return(fallback_summary)
  }
  
  # score here is 0-1 scale (STRING API convention), hence /1000 below
  edges <- edges[edges$score >= required_score / 1000, ]
  
  g <- igraph::graph_from_data_frame(
    edges[, c("preferredName_A", "preferredName_B")],
    directed = FALSE
  )
  g <- igraph::simplify(g)  # drop any duplicate edges / self-loops
  
  deg <- igraph::degree(g)
  bet <- igraph::betweenness(g)
  
  summary_df <- data.frame(Symbol = names(deg), degree = deg, betweenness = bet, row.names = NULL)
  
  if (return_graph) {
    return(list(summary = summary_df, graph = g, edges = edges))
  }
  summary_df
}

# ------------------------------------------------------------------
# Quick local smoke test (run this locally first, not in sandbox):
# test_genes <- c("FN1","ITGA2","ITGB1","AGRN","COL1A1","ITGAV")
# test_net <- build_string_network(test_genes)
# print(test_net)
#
# For the network figure panel:
# test_net_full <- build_string_network(test_genes, return_graph = TRUE)
# test_net_full$graph
# ------------------------------------------------------------------