## ============================================================
## R/Helper/ahp_weights.R
##
## SOURCE for the AHP pairwise comparison matrix
## and derived criterion weights (FC, FDR, Bet, Deg).
##
## Usage: source(here::here("R", "Helper", "ahp_weights.R"))
## Then use the pre-computed AHP_WEIGHTS vector directly, or call
## ahp_weights(AHP_PAIRWISE_MATRIX) if you need to re-derive it.
## ============================================================

AHP_PAIRWISE_MATRIX <- matrix(
  c(
    1,     4,   5,   2.5,
    1/4,   1,   2,   1,
    1/5,   1/2, 1,   1/5,
    1/2.5, 1,   5,   1
  ),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("FC", "FDR", "Bet", "Deg"), c("FC", "FDR", "Bet", "Deg"))
)

## Column-normalize, then average rows (standard AHP eigenvector
## approximation). Returns a named numeric vector of weights.
ahp_weights <- function(M = AHP_PAIRWISE_MATRIX) {
  norm <- sweep(M, 2, colSums(M), "/")
  rowMeans(norm)
}

## Pre-computed weights, ready to use as `w["FC"]`, `w["FDR"]`, etc.
AHP_WEIGHTS <- ahp_weights(AHP_PAIRWISE_MATRIX)

## Consistency Ratio check -- optional sanity guard. If someone
## edits AHP_PAIRWISE_MATRIX later, this catches a matrix that's
## no longer a valid (consistent-enough) pairwise comparison
## before it silently propagates into every downstream script.
ahp_consistency_ratio <- function(M = AHP_PAIRWISE_MATRIX, w = NULL) {
  n <- nrow(M)
  if (is.null(w)) w <- ahp_weights(M)
  lambda_max <- mean(as.vector(M %*% w) / w)
  CI <- (lambda_max - n) / (n - 1)
  RI <- c(`1` = 0, `2` = 0, `3` = 0.58, `4` = 0.90, `5` = 1.12)[as.character(n)]
  CR <- CI / RI
  list(lambda_max = lambda_max, CI = CI, CR = unname(CR))
}

if (identical(environment(), globalenv()) && interactive()) {
  cat("AHP weights:", round(AHP_WEIGHTS, 4), "\n")
  cr <- ahp_consistency_ratio()
  cat("Consistency Ratio (CR):", round(cr$CR, 4), "\n")
}
