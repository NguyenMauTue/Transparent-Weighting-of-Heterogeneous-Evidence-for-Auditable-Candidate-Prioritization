## ============================================================
## R/Helper/split_helper.R
##
## Shared 3-vs-3 split enumeration and per-split limma DE runner.
## Previously duplicated verbatim in 05_null_permutations.R,
## 14_fconly_null_test.R, 14b_fconly_go_ecm_null.R,
## 14c_equalweight_go_ecm_null.R, 16_double_null_2x2.R.
##
## Relevant to the group-A fix: 16's real_r used to be computed
## from a DE run on the TRUE split that was built via yet another
## copy of this same split logic (implicitly, since TRUE = full
## data with real labels) -- using one shared implementation
## removes any chance that a stray edit to one of the 5 copies
## changes which 3 samples count as "true GroupA" in one script
## but not another.
## ============================================================

suppressMessages(library(limma))

## Enumerate all 10 unique 3-vs-3 partitions of 6 samples (indices 1:6).
## Convention: split 1 is always the TRUE split, i.e. the one where
## GroupA = c(1,2,3) (matching column order Normal1-3, Tumor1-3 in
## imputed_matrix). Returns list(true_split, null_splits, all_splits).
enumerate_3v3_splits <- function() {
  idx <- 1:6
  combos <- combn(6, 3)
  splits <- list(); seen <- c()
  for (i in 1:ncol(combos)) {
    A <- sort(combos[, i]); B <- sort(setdiff(idx, A))
    if (!(1 %in% A)) { tmp <- A; A <- B; B <- tmp }
    key <- paste(A, collapse = ",")
    if (!(key %in% seen)) {
      seen <- c(seen, key)
      splits[[length(splits) + 1]] <- list(A = A, B = B)
    }
  }
  stopifnot(length(splits) == 10)
  list(
    true_split  = splits[[1]],
    null_splits = splits[-1],
    all_splits  = splits
  )
}

## Run limma DE (GroupB - GroupA contrast) for a given 3-vs-3 sample split.
## imputed_matrix: samples in columns, in the fixed order the split indices
## (1:6) refer to.
run_de_for_split <- function(imputed_matrix, A_idx, B_idx) {
  group <- factor(rep(NA, 6), levels = c("GroupA", "GroupB"))
  group[A_idx] <- "GroupA"; group[B_idx] <- "GroupB"
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  fit <- lmFit(imputed_matrix, design)
  cm <- makeContrasts(GroupB_VS_GroupA = GroupB - GroupA, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  res <- topTable(fit2, coef = "GroupB_VS_GroupA", number = Inf, sort.by = "none")
  res$X <- rownames(res)
  res |> dplyr::mutate(UNIPROT = sapply(strsplit(X, ";"), `[`, 1))
}
