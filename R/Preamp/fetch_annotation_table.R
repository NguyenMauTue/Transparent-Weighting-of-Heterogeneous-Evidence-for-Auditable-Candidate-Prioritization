############################################################
# fetch_annotation_table.R
# 
############################################################
library(biomaRt)
library(dplyr)
source(here::here("R", "Helper", "get_annotation_raw.R"))
# Loading mart
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  host    = "https://jun2026.archive.ensembl.org"  # Ensembl 116
)
#Raw annotation for PXD056
raw = fread("https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/03/PXD056161/proteinGroups.txt")
annotations = get_annotation_raw(raw)
write.csv(annotations, here::here("PXD056161", "data", "annotation_raw.csv"), row.names = FALSE)
