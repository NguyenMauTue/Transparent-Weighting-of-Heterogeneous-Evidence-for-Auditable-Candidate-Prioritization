get_annotation_raw = function(df) {
  id_col <- grep("Protein.*ID|uniprot", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  
  if (is.na(id_col)) {
    stop("Can't find UniProt ID in data frame.")
  }
  
  uniprot_universe <- df[[id_col]] %>%
    as.character() %>%
    strsplit(";") %>%
    sapply(function(x) x[1]) %>%
    unique()
  
  uniprot_universe <- uniprot_universe[!is.na(uniprot_universe) & uniprot_universe != ""]
  
  annotations = getBM(
    attributes = c("uniprotswissprot", "external_gene_name", "description",
                   "go_id", "name_1006", "namespace_1003"),
    filters = "uniprotswissprot",
    values = uniprot_universe,
    mart = mart
  )
  
  return(annotations)
}
