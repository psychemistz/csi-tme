#!/usr/bin/env Rscript
# Extract and save SpaCET cell type proportions for all SR samples
# Saves per-spot proportion matrices as CSV for downstream use

library(SpaCET)
library(Matrix)
library(reticulate)

project_root <- "/vf/users/parks34/projects/3csitme"
sr_dir <- file.path(project_root, "data/zenodo_16905935/human_glioma_sr")
output_dir <- file.path(project_root, "reports/spacet_proportions")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Python helper for loading h5ad
py_run_string("
import numpy as np
import scipy.sparse as sp

def extract_h5ad_data(h5ad_path):
    import anndata as ad
    adata = ad.read_h5ad(h5ad_path)
    if sp.issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = np.asarray(adata.X)
    obs_names = adata.obs_names.tolist()
    var_names = adata.var_names.tolist()
    obs_dict = {}
    for col in adata.obs.columns:
        obs_dict[col] = adata.obs[col].tolist()
    spatial_coords = None
    if 'spatial' in adata.obsm:
        spatial_coords = adata.obsm['spatial'].tolist()
    return {
        'X': X.tolist(),
        'obs_names': obs_names,
        'var_names': var_names,
        'obs': obs_dict,
        'spatial': spatial_coords,
        'n_obs': adata.n_obs,
        'n_vars': adata.n_vars
    }
")

load_h5ad_counts <- function(h5ad_path) {
  result <- py$extract_h5ad_data(h5ad_path)
  X_list <- result$X
  counts_r <- do.call(rbind, lapply(X_list, unlist))
  rownames(counts_r) <- result$obs_names
  colnames(counts_r) <- result$var_names
  counts_r <- t(counts_r)
  counts_sparse <- as(counts_r, "dgCMatrix")

  spatial <- NULL
  obs_data <- result$obs
  if ("array_row" %in% names(obs_data) && "array_col" %in% names(obs_data)) {
    spatial <- data.frame(
      x = as.numeric(unlist(obs_data$array_col)),
      y = as.numeric(unlist(obs_data$array_row)),
      row.names = result$obs_names
    )
  }
  return(list(counts = counts_sparse, spatial = spatial))
}

# Get samples
samples <- list.dirs(sr_dir, recursive = FALSE, full.names = FALSE)
samples <- samples[!grepl("^\\.", samples)]

cat(sprintf("Processing %d SR samples\n", length(samples)))

for (sample_id in sort(samples)) {
  output_file <- file.path(output_dir, paste0(sample_id, "_proportions.csv"))

  if (file.exists(output_file)) {
    cat(sprintf("  %s: already exists, skipping\n", sample_id))
    next
  }

  cat(sprintf("\n=== %s ===\n", sample_id))
  gene_h5ad <- file.path(sr_dir, sample_id, "gene.quant.h5ad")

  if (!file.exists(gene_h5ad)) {
    cat("  Gene counts not found, skipping\n")
    next
  }

  tryCatch({
    cat("  Loading gene counts...\n")
    gene_data <- load_h5ad_counts(gene_h5ad)
    cat(sprintf("  %d genes x %d spots\n", nrow(gene_data$counts), ncol(gene_data$counts)))

    if (is.null(gene_data$spatial)) {
      cat("  No spatial coordinates, skipping\n")
      next
    }

    cat("  Running SpaCET deconvolution...\n")
    SpaCET_obj <- SpaCET::create.SpaCET.object(
      counts = gene_data$counts,
      spotCoordinates = gene_data$spatial,
      platform = "Visium"
    )
    SpaCET_obj <- SpaCET::SpaCET.deconvolution(
      SpaCET_obj, cancerType = "GBM", coreNo = 4
    )

    prop_mat <- SpaCET_obj@results$deconvolution$propMat
    if (nrow(prop_mat) < ncol(prop_mat)) {
      prop_mat <- t(prop_mat)
    }

    prop_df <- as.data.frame(prop_mat)
    prop_df$spot <- rownames(prop_mat)
    write.csv(prop_df, output_file, row.names = FALSE)
    cat(sprintf("  Saved: %s (%d spots x %d cell types)\n",
                basename(output_file), nrow(prop_mat), ncol(prop_mat)))

  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
  })
}

cat("\nDone! Proportions saved to: ", output_dir, "\n")
