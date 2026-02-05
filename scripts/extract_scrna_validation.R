#!/usr/bin/env Rscript
# Extract scRNA-seq validation data for HSIC-IR candidate genes
# Output: CSV files with PSI values and differential splicing results

library(dplyr)
library(tidyr)

# Set paths
project_root <- "/vf/users/parks34/projects/3csitme"
marvel_dir <- file.path(project_root, "data/zenodo_17055113/zenodo-20250904/MARVEL")
output_dir <- file.path(project_root, "reports")

# Candidate genes from ONT HSIC-IR analysis
candidates <- c('B2M', 'CLU', 'CD74', 'MIF', 'LGALS1', 'PTGDS', 'CST3', 'SPARCL1',
                'TIMP1', 'ITM2B', 'LGALS3', 'SPP1', 'AGT', 'C3', 'SPARC', 'APOE',
                'APP', 'VCAN', 'HMGB1', 'IGFBP5', 'CCL2', 'SERPINE2')

cat("========================================\n")
cat("Extracting scRNA-seq Validation Data\n")
cat("========================================\n")

# 1. Load tumor state differential splicing events
cat("\n1. Loading tumor state differential events...\n")
tumor_state_file <- file.path(marvel_dir, "event_list_tumor_state_diff_wilcox.RData")

if (file.exists(tumor_state_file)) {
  load(tumor_state_file)
  cat("   Loaded: event_list_tumor_state_diff_wilcox.RData\n")

  # List objects loaded
  cat("   Objects in environment:\n")
  for (obj_name in ls()) {
    obj <- get(obj_name)
    if (is.data.frame(obj)) {
      cat(sprintf("     - %s: %d rows, %d cols\n", obj_name, nrow(obj), ncol(obj)))
    } else if (is.list(obj)) {
      cat(sprintf("     - %s: list with %d elements\n", obj_name, length(obj)))
    }
  }
}

# 2. Load PSI data
cat("\n2. Loading PSI merge data...\n")
psi_file <- file.path(marvel_dir, "PSI_merge.RData")

if (file.exists(psi_file)) {
  load(psi_file)
  cat("   Loaded: PSI_merge.RData\n")

  # List objects
  for (obj_name in ls()) {
    obj <- get(obj_name)
    if (is.data.frame(obj) || is.matrix(obj)) {
      cat(sprintf("     - %s: %d rows, %d cols\n", obj_name, nrow(obj), ncol(obj)))
      if (ncol(obj) > 0) {
        cat(sprintf("       First cols: %s\n", paste(head(colnames(obj), 5), collapse=", ")))
      }
    }
  }
}

# 3. Load splice feature annotations
cat("\n3. Loading splice feature data...\n")
splice_file <- file.path(marvel_dir, "splice_feature.RData")

if (file.exists(splice_file)) {
  load(splice_file)
  cat("   Loaded: splice_feature.RData\n")

  for (obj_name in ls()) {
    obj <- get(obj_name)
    if (is.data.frame(obj)) {
      cat(sprintf("     - %s: %d rows, %d cols\n", obj_name, nrow(obj), ncol(obj)))
      # Show column names
      if (ncol(obj) > 0 && ncol(obj) <= 20) {
        cat(sprintf("       Columns: %s\n", paste(colnames(obj), collapse=", ")))
      }
    }
  }
}

# 4. Extract candidate gene data
cat("\n4. Searching for candidate genes in loaded data...\n")

# Try to find gene column and filter for candidates
all_results <- list()

for (obj_name in ls()) {
  obj <- get(obj_name)

  if (is.data.frame(obj)) {
    # Look for gene-related columns
    gene_cols <- grep("gene|Gene|GENE", colnames(obj), value=TRUE)

    for (gene_col in gene_cols) {
      # Check if any candidates are present
      found_genes <- intersect(toupper(obj[[gene_col]]), toupper(candidates))

      if (length(found_genes) > 0) {
        cat(sprintf("   Found %d candidates in %s$%s\n", length(found_genes), obj_name, gene_col))

        # Extract rows for candidates
        candidate_data <- obj[toupper(obj[[gene_col]]) %in% toupper(candidates), ]

        if (nrow(candidate_data) > 0) {
          candidate_data$source_object <- obj_name
          candidate_data$gene_column <- gene_col
          all_results[[paste(obj_name, gene_col, sep="_")]] <- candidate_data

          # Print found genes
          unique_genes <- unique(candidate_data[[gene_col]])
          cat(sprintf("     Genes found: %s\n", paste(unique_genes, collapse=", ")))
        }
      }
    }
  }
}

# 5. Combine and save results
cat("\n5. Saving results...\n")

if (length(all_results) > 0) {
  # Try to combine similar data frames
  for (name in names(all_results)) {
    df <- all_results[[name]]
    output_file <- file.path(output_dir, paste0("scrna_", gsub("[^a-zA-Z0-9]", "_", name), ".csv"))
    write.csv(df, output_file, row.names=FALSE)
    cat(sprintf("   Saved: %s (%d rows)\n", basename(output_file), nrow(df)))
  }
}

# 6. Generate summary of differential splicing events by tumor state
cat("\n6. Looking for differential splicing results...\n")

# Check for event_list objects which contain differential analysis results
event_objs <- ls()[grep("event_list", ls())]

for (obj_name in event_objs) {
  obj <- get(obj_name)

  if (is.list(obj) && !is.data.frame(obj)) {
    cat(sprintf("\n   Processing %s (list with %d elements):\n", obj_name, length(obj)))

    for (elem_name in names(obj)) {
      elem <- obj[[elem_name]]

      if (is.data.frame(elem)) {
        cat(sprintf("     - %s: %d rows\n", elem_name, nrow(elem)))

        # Look for gene column
        gene_cols <- grep("gene|Gene", colnames(elem), value=TRUE)

        if (length(gene_cols) > 0) {
          gene_col <- gene_cols[1]
          found <- elem[toupper(elem[[gene_col]]) %in% toupper(candidates), ]

          if (nrow(found) > 0) {
            cat(sprintf("       Found %d events for candidates\n", nrow(found)))
            output_file <- file.path(output_dir,
                                     paste0("scrna_diff_", gsub("[^a-zA-Z0-9]", "_", paste(obj_name, elem_name)), ".csv"))
            write.csv(found, output_file, row.names=FALSE)
            cat(sprintf("       Saved: %s\n", basename(output_file)))
          }
        }
      }
    }
  } else if (is.data.frame(obj)) {
    cat(sprintf("\n   Processing %s (data.frame with %d rows):\n", obj_name, nrow(obj)))

    gene_cols <- grep("gene|Gene", colnames(obj), value=TRUE)

    if (length(gene_cols) > 0) {
      gene_col <- gene_cols[1]
      found <- obj[toupper(obj[[gene_col]]) %in% toupper(candidates), ]

      if (nrow(found) > 0) {
        cat(sprintf("     Found %d events for candidates\n", nrow(found)))
        output_file <- file.path(output_dir,
                                 paste0("scrna_diff_", gsub("[^a-zA-Z0-9]", "_", obj_name), ".csv"))
        write.csv(found, output_file, row.names=FALSE)
        cat(sprintf("     Saved: %s\n", basename(output_file)))
      }
    }
  }
}

cat("\n========================================\n")
cat("Extraction complete!\n")
cat("========================================\n")

# List all generated CSV files
cat("\nGenerated CSV files:\n")
csv_files <- list.files(output_dir, pattern="scrna.*\\.csv$", full.names=FALSE)
for (f in csv_files) {
  cat(sprintf("  - %s\n", f))
}
