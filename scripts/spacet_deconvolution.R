#!/usr/bin/env Rscript
# SpaCET deconvolution analysis to verify isoform switching between tumor and macrophage
# Correlate VCAN isoform ratios with cell type proportions

library(SpaCET)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reticulate)

# Set paths
project_root <- "/vf/users/parks34/projects/3csitme"
sr_dir <- file.path(project_root, "data/zenodo_16905935/human_glioma_sr")
output_dir <- file.path(project_root, "reports")

cat("================================================================\n")
cat("SpaCET Deconvolution Analysis\n")
cat("================================================================\n")

# Import Python modules
ad <- import("anndata")
np <- import("numpy")

# Define Python helper function to extract matrix as list
py_run_string("
import numpy as np
import scipy.sparse as sp

def extract_h5ad_data(h5ad_path):
    import anndata as ad
    adata = ad.read_h5ad(h5ad_path)

    # Get matrix as dense numpy array
    if sp.issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = np.asarray(adata.X)

    # Convert to list of lists for R
    X_list = X.tolist()

    # Get names
    obs_names = adata.obs_names.tolist()
    var_names = adata.var_names.tolist()

    # Get spatial coords from obs
    obs_dict = {}
    for col in adata.obs.columns:
        obs_dict[col] = adata.obs[col].tolist()

    # Get obsm spatial if available
    spatial_coords = None
    if 'spatial' in adata.obsm:
        spatial_coords = adata.obsm['spatial'].tolist()

    return {
        'X': X_list,
        'obs_names': obs_names,
        'var_names': var_names,
        'obs': obs_dict,
        'spatial': spatial_coords,
        'n_obs': adata.n_obs,
        'n_vars': adata.n_vars
    }
")

# Function to load h5ad gene counts using Python helper
load_h5ad_counts <- function(h5ad_path) {
  # Call Python function
  result <- py$extract_h5ad_data(h5ad_path)

  n_obs <- result$n_obs
  n_vars <- result$n_vars
  obs_names <- result$obs_names
  var_names <- result$var_names

  cat(sprintf("    AnnData shape: %d obs x %d vars\n", n_obs, n_vars))

  # Convert X list to R matrix
  X_list <- result$X
  counts_r <- do.call(rbind, lapply(X_list, unlist))

  cat(sprintf("    R matrix shape: %d x %d\n", nrow(counts_r), ncol(counts_r)))

  # Set dimension names
  rownames(counts_r) <- obs_names
  colnames(counts_r) <- var_names

  # Transpose to genes x spots for SpaCET
  counts_r <- t(counts_r)

  # Convert to sparse matrix
  counts_sparse <- as(counts_r, "dgCMatrix")

  cat(sprintf("    Final counts: %d genes x %d spots\n", nrow(counts_sparse), ncol(counts_sparse)))

  # Get spatial coordinates
  obs_data <- result$obs
  spatial <- NULL

  if ("array_row" %in% names(obs_data) && "array_col" %in% names(obs_data)) {
    spatial <- data.frame(
      x = as.numeric(unlist(obs_data$array_col)),
      y = as.numeric(unlist(obs_data$array_row)),
      row.names = obs_names
    )
  } else if (!is.null(result$spatial)) {
    sp_coords <- do.call(rbind, lapply(result$spatial, unlist))
    spatial <- data.frame(
      x = sp_coords[, 1],
      y = sp_coords[, 2],
      row.names = obs_names
    )
  }

  return(list(counts = counts_sparse, spatial = spatial))
}

# Define Python helper for isoform extraction
py_run_string("
def extract_isoform_data(h5ad_path, gene_name):
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp
    import re

    adata = ad.read_h5ad(h5ad_path)

    var_names = adata.var_names.tolist()
    obs_names = adata.obs_names.tolist()

    # Find isoforms for the gene
    gene_isos = []

    # Check for gene column in var
    gene_col = None
    for col in ['gene_symbol', 'gene', 'gene_name']:
        if col in adata.var.columns:
            gene_col = col
            break

    if gene_col is not None:
        mask = adata.var[gene_col] == gene_name
        gene_isos = [v for v, m in zip(var_names, mask) if m]
    else:
        # Try to find by prefix
        pattern = f'^{gene_name}[-_]'
        gene_isos = [v for v in var_names if re.match(pattern, v)]

    if len(gene_isos) < 2:
        return None

    # Get indices
    iso_idx = [var_names.index(iso) for iso in gene_isos]

    # Get matrix
    if sp.issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = np.asarray(adata.X)

    # Extract isoform columns
    iso_counts = X[:, iso_idx].tolist()

    return {
        'iso_counts': iso_counts,
        'gene_isos': gene_isos,
        'obs_names': obs_names,
        'n_isoforms': len(gene_isos)
    }
")

# Function to load isoform data and calculate ratios
load_isoform_ratios <- function(h5ad_path, gene_name) {
  # Call Python function
  result <- py$extract_isoform_data(h5ad_path, gene_name)

  if (is.null(result)) {
    return(NULL)
  }

  gene_isos <- result$gene_isos
  obs_names <- result$obs_names
  n_isoforms <- result$n_isoforms

  # Convert to R matrix
  iso_counts <- do.call(rbind, lapply(result$iso_counts, unlist))
  colnames(iso_counts) <- gene_isos
  rownames(iso_counts) <- obs_names

  # Calculate total and ratios
  total <- rowSums(iso_counts)

  # Calculate ratio of first isoform (or dominant isoform)
  ratios <- data.frame(
    spot = obs_names,
    total_counts = total,
    n_isoforms = n_isoforms,
    isoform_names = paste(gene_isos, collapse = ";")
  )

  # Add individual isoform fractions
  for (i in seq_along(gene_isos)) {
    ratios[[paste0("iso", i, "_frac")]] <- iso_counts[, i] / pmax(total, 1)
    ratios[[paste0("iso", i, "_counts")]] <- iso_counts[, i]
  }

  return(ratios)
}

# Get list of SR samples
samples <- list.dirs(sr_dir, recursive = FALSE, full.names = FALSE)
samples <- samples[!grepl("^\\.", samples)]

cat(sprintf("\nFound %d SR samples:\n", length(samples)))
print(samples)

# Process each sample
all_results <- list()

for (sample_id in samples) {
  cat(sprintf("\n\n========== Processing %s ==========\n", sample_id))

  sample_dir <- file.path(sr_dir, sample_id)
  gene_h5ad <- file.path(sample_dir, "gene.quant.h5ad")
  iso_h5ad <- file.path(sample_dir, "iso.quant.h5ad")

  if (!file.exists(gene_h5ad)) {
    cat("  Gene counts not found, skipping\n")
    next
  }

  tryCatch({
    # 1. Load gene counts
    cat("  Loading gene counts...\n")
    gene_data <- load_h5ad_counts(gene_h5ad)

    counts <- gene_data$counts
    spatial <- gene_data$spatial

    cat(sprintf("  Counts: %d genes x %d spots\n", nrow(counts), ncol(counts)))

    if (is.null(spatial)) {
      cat("  No spatial coordinates found, skipping\n")
      next
    }

    # 2. Create SpaCET object
    cat("  Creating SpaCET object...\n")

    SpaCET_obj <- SpaCET::create.SpaCET.object(
      counts = counts,
      spotCoordinates = spatial,
      platform = "Visium"
    )

    # 3. Run deconvolution (GBM cancer type)
    cat("  Running SpaCET deconvolution (GBM)...\n")
    SpaCET_obj <- SpaCET::SpaCET.deconvolution(
      SpaCET_obj,
      cancerType = "GBM",
      coreNo = 4
    )

    # 4. Extract cell type proportions
    cat("  Extracting cell type proportions...\n")
    prop_mat <- SpaCET_obj@results$deconvolution$propMat

    # Check propMat structure - might be cell_types x spots or spots x cell_types
    cat(sprintf("  PropMat dimensions: %d x %d\n", nrow(prop_mat), ncol(prop_mat)))

    # SpaCET propMat has cell_types as rows, spots as columns
    # We need to transpose to get spots x cell_types
    if (nrow(prop_mat) < ncol(prop_mat)) {
      cat("  Transposing propMat (cell_types x spots -> spots x cell_types)\n")
      prop_mat <- t(prop_mat)
    }

    # Keep only major lineage cell types (Level 1)
    # Major lineage: Malignant, CAF, Endothelial, Plasma, B cell, T CD4, T CD8,
    #                NK, cDC, pDC, Macrophage, Mast, Neutrophil, Unidentifiable
    major_lineage <- c("Malignant", "CAF", "Endothelial", "Plasma", "B cell",
                       "T CD4", "T CD8", "NK", "cDC", "pDC", "Macrophage",
                       "Mast", "Neutrophil", "Unidentifiable")

    available_major <- intersect(colnames(prop_mat), major_lineage)
    if (length(available_major) > 0) {
      prop_mat <- prop_mat[, available_major, drop = FALSE]
      cat(sprintf("  Using %d major lineage cell types only\n", length(available_major)))
    }

    cat("  Cell types detected:\n")
    print(colnames(prop_mat))

    # 5. Load VCAN isoform ratios
    cat("  Loading VCAN isoform data...\n")
    vcan_ratios <- load_isoform_ratios(iso_h5ad, "VCAN")

    if (is.null(vcan_ratios)) {
      cat("  VCAN isoforms not found, trying other genes...\n")

      # Try other validated genes
      for (gene in c("APP", "HMGB1", "CLU", "SPP1")) {
        vcan_ratios <- load_isoform_ratios(iso_h5ad, gene)
        if (!is.null(vcan_ratios)) {
          cat(sprintf("  Using %s instead (%d isoforms)\n", gene, vcan_ratios$n_isoforms[1]))
          break
        }
      }
    }

    if (is.null(vcan_ratios)) {
      cat("  No isoform data found, skipping correlation\n")

      # Still save proportions
      all_results[[sample_id]] <- list(
        proportions = prop_mat,
        vcan_ratios = NULL,
        sample_id = sample_id
      )
      next
    }

    cat(sprintf("  Isoforms: %d detected\n", vcan_ratios$n_isoforms[1]))

    # 6. Merge proportions with isoform ratios
    cat("  Correlating cell types with isoform ratios...\n")

    # Match spots
    common_spots <- intersect(rownames(prop_mat), vcan_ratios$spot)
    cat(sprintf("  Common spots: %d\n", length(common_spots)))

    if (length(common_spots) < 50) {
      cat("  Too few common spots, skipping\n")
      all_results[[sample_id]] <- list(
        proportions = prop_mat,
        vcan_ratios = vcan_ratios,
        sample_id = sample_id
      )
      next
    }

    # Create merged data
    prop_df <- as.data.frame(prop_mat[common_spots, ])
    prop_df$spot <- common_spots

    merged <- merge(prop_df, vcan_ratios, by = "spot")

    # Filter to spots with expression
    merged_expressed <- merged[merged$total_counts > 0, ]
    cat(sprintf("  Spots with expression: %d\n", nrow(merged_expressed)))

    if (nrow(merged_expressed) < 30) {
      cat("  Too few expressed spots\n")
      next
    }

    # 7. Calculate correlations
    cat("\n  Correlations with isoform 1 fraction:\n")

    cors <- list()
    cell_types <- setdiff(colnames(prop_mat), c("spot"))

    for (cell_type in cell_types) {
      if (cell_type %in% colnames(merged_expressed) && "iso1_frac" %in% colnames(merged_expressed)) {
        cor_test <- cor.test(merged_expressed[[cell_type]], merged_expressed$iso1_frac,
                            method = "spearman")
        cors[[cell_type]] <- data.frame(
          sample = sample_id,
          cell_type = cell_type,
          correlation = cor_test$estimate,
          pvalue = cor_test$p.value,
          n_spots = nrow(merged_expressed)
        )

        if (abs(cor_test$estimate) > 0.1 || cor_test$p.value < 0.05) {
          cat(sprintf("    %s: r = %.3f, p = %.2e\n",
                     cell_type, cor_test$estimate, cor_test$p.value))
        }
      }
    }

    # Store results
    all_results[[sample_id]] <- list(
      proportions = prop_mat,
      vcan_ratios = vcan_ratios,
      merged = merged_expressed,
      correlations = bind_rows(cors),
      sample_id = sample_id
    )

    # 8. Create visualization
    if (nrow(merged_expressed) > 50) {
      # Find macrophage-like columns
      macro_cols <- grep("macro|Macro|TAM|Myeloid", colnames(prop_mat), value = TRUE, ignore.case = TRUE)
      tumor_cols <- grep("malig|Malig|tumor|Tumor|cancer|Cancer|Neoplastic", colnames(prop_mat), value = TRUE, ignore.case = TRUE)

      if (length(macro_cols) > 0) {
        macro_col <- macro_cols[1]

        p <- ggplot(merged_expressed, aes_string(x = macro_col, y = "iso1_frac")) +
          geom_point(alpha = 0.3, size = 1) +
          geom_smooth(method = "lm", color = "red") +
          labs(title = paste(sample_id, "- Isoform 1 Fraction vs", macro_col),
               x = paste(macro_col, "Proportion"),
               y = "Isoform 1 Fraction") +
          theme_minimal()

        ggsave(file.path(output_dir, paste0("spacet_", sample_id, "_macro.png")),
               p, width = 6, height = 5)
        cat(sprintf("  Saved: spacet_%s_macro.png\n", sample_id))
      }
    }

  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
  })
}

# Combine all correlations
cat("\n\n================================================================\n")
cat("SUMMARY OF RESULTS\n")
cat("================================================================\n")

if (length(all_results) > 0) {
  # Collect all correlations
  all_cors <- bind_rows(lapply(all_results, function(x) {
    if (!is.null(x$correlations)) x$correlations else NULL
  }))

  if (nrow(all_cors) > 0) {
    # Save correlations
    write.csv(all_cors, file.path(output_dir, "spacet_isoform_correlations.csv"), row.names = FALSE)
    cat(sprintf("\nSaved: spacet_isoform_correlations.csv (%d rows)\n", nrow(all_cors)))

    # Summary by cell type
    cor_summary <- all_cors %>%
      group_by(cell_type) %>%
      summarise(
        mean_cor = mean(correlation, na.rm = TRUE),
        sd_cor = sd(correlation, na.rm = TRUE),
        n_samples = n(),
        n_sig = sum(pvalue < 0.05, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      arrange(desc(abs(mean_cor)))

    cat("\nCorrelation of isoform ratio with cell type proportions:\n")
    print(cor_summary, n = 20)

    write.csv(cor_summary, file.path(output_dir, "spacet_correlation_summary.csv"), row.names = FALSE)

    # Focus on macrophage vs tumor
    macro_cors <- all_cors %>% filter(grepl("macro|Macro|TAM|Myeloid", cell_type, ignore.case = TRUE))
    tumor_cors <- all_cors %>% filter(grepl("malig|tumor|cancer|Neoplastic", cell_type, ignore.case = TRUE))

    if (nrow(macro_cors) > 0) {
      cat("\n\nMacrophage correlation summary:\n")
      cat(sprintf("  Mean r = %.3f (SD = %.3f)\n", mean(macro_cors$correlation), sd(macro_cors$correlation)))
      cat(sprintf("  Significant samples: %d/%d\n", sum(macro_cors$pvalue < 0.05), nrow(macro_cors)))
    }

    if (nrow(tumor_cors) > 0) {
      cat("\nTumor/Malignant correlation summary:\n")
      cat(sprintf("  Mean r = %.3f (SD = %.3f)\n", mean(tumor_cors$correlation), sd(tumor_cors$correlation)))
      cat(sprintf("  Significant samples: %d/%d\n", sum(tumor_cors$pvalue < 0.05), nrow(tumor_cors)))
    }
  }
}

cat("\n================================================================\n")
cat("Analysis complete!\n")
cat("================================================================\n")
