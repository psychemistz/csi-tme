#!/usr/bin/env Rscript
# Compute pairwise Moran's I with grid weights for a subset of genes
# Usage: Rscript run_r_grid.R <vst_file> <n_genes> <output_file>

args <- commandArgs(trailingOnly = TRUE)
vst_file <- args[1]
n_genes <- as.integer(args[2])
output_file <- args[3]

library(sigdiscov)

cat("Loading:", vst_file, "\n")
vst <- read.delim(vst_file, row.names = 1, check.names = FALSE)
cat("  Shape:", nrow(vst), "genes x", ncol(vst), "spots\n")

# Parse spot coordinates
spot_names <- colnames(vst)
coords <- do.call(rbind, lapply(strsplit(spot_names, "x"), as.integer))
coords_df <- data.frame(row = coords[, 1], col = coords[, 2])

# Create grid weights (matching Python create_grid_weights)
cat("Creating grid weights...\n")
wt <- create_weights_visium(coords_df, radius = 300, type = "grid",
                            platform = "visium", include_self = FALSE)
W <- wt$W
weight_sum <- wt$weight_sum
cat("  Weight sum:", weight_sum, "\n")
cat("  NNZ:", Matrix::nnzero(W), "\n")

# Subset genes
n_genes <- min(n_genes, nrow(vst))
expr <- as.matrix(vst[1:n_genes, ])

# Standardize (z-score per gene, population SD matching Python/C++)
expr_z <- t(apply(expr, 1, function(x) {
  n <- length(x)
  m <- mean(x)
  s <- sqrt(sum((x - m)^2) / n)  # population SD (N, not N-1)
  if (s < 1e-10) return(rep(0, n))
  (x - m) / s
}))

# Compute pairwise Moran's I: I = (Z @ W @ Z') / weight_sum
cat("Computing pairwise Moran's I for", n_genes, "genes...\n")
ZW <- expr_z %*% W
I_matrix <- (ZW %*% t(expr_z)) / weight_sum

# Save as binary for exact comparison
cat("Saving to:", output_file, "\n")
writeBin(as.double(I_matrix), output_file)

# Also save weight_sum
ws_file <- paste0(output_file, ".ws")
writeBin(weight_sum, ws_file)

cat("Done.\n")
