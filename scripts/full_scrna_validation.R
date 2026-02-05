#!/usr/bin/env Rscript
# Full scRNA-seq validation: Extract PSI values and run differential analysis
# for all HSIC-IR candidate genes

library(dplyr)
library(tidyr)

# Set paths
project_root <- "/vf/users/parks34/projects/3csitme"
marvel_dir <- file.path(project_root, "data/zenodo_17055113/zenodo-20250904/MARVEL")
seurat_dir <- file.path(project_root, "data/zenodo_17055113/zenodo-20250904/seurat")
output_dir <- file.path(project_root, "reports")

# Candidate genes from ONT HSIC-IR analysis
candidates <- c('B2M', 'CLU', 'CD74', 'MIF', 'LGALS1', 'PTGDS', 'CST3', 'SPARCL1',
                'TIMP1', 'ITM2B', 'LGALS3', 'SPP1', 'AGT', 'C3', 'SPARC', 'APOE',
                'APP', 'VCAN', 'HMGB1', 'IGFBP5', 'CCL2', 'SERPINE2')

cat("================================================================\n")
cat("Full scRNA-seq Validation Analysis\n")
cat("================================================================\n")

# 1. Load cell metadata
cat("\n1. Loading cell metadata...\n")
load(file.path(seurat_dir, "metadata_seu_filt_all.RData"))
cat("   Loaded metadata\n")
cat(sprintf("   Metadata: %d cells, %d columns\n", nrow(metadata), ncol(metadata)))
cat(sprintf("   Columns: %s\n", paste(colnames(metadata), collapse=", ")))

# 2. Load splice features to map events to genes
cat("\n2. Loading splice features...\n")
load(file.path(marvel_dir, "splice_feature.RData"))

# Filter to candidate genes
candidate_events <- splice_feature %>%
  filter(gene_short_name %in% candidates)

cat(sprintf("   Found %d splice events for %d candidate genes\n",
            nrow(candidate_events),
            length(unique(candidate_events$gene_short_name))))

# 3. Load PSI matrix
cat("\n3. Loading PSI matrix...\n")
load(file.path(marvel_dir, "PSI_merge.RData"))

cat(sprintf("   PSI matrix: %d events x %d cells\n", nrow(PSI_merge), ncol(PSI_merge)))

# Get candidate event IDs
candidate_event_ids <- candidate_events$id

# Filter PSI matrix to candidate events
psi_candidates <- PSI_merge[rownames(PSI_merge) %in% candidate_event_ids, ]
cat(sprintf("   Filtered to %d candidate events with PSI data\n", nrow(psi_candidates)))

# 4. Get tumor state annotations
cat("\n4. Getting tumor state annotations...\n")

# Use cell_idents_tumor column which contains tumor cell states
state_col <- "cell_idents_tumor"
if (state_col %in% colnames(metadata)) {
  cat(sprintf("   Using column: %s\n", state_col))
  cat("\n   Cell state distribution:\n")
  print(table(metadata[[state_col]], useNA = "ifany"))

  # Get non-NA states
  tumor_states <- unique(metadata[[state_col]])
  tumor_states <- tumor_states[!is.na(tumor_states)]
  cat(sprintf("\n   Found %d tumor states\n", length(tumor_states)))
} else {
  # Try alternative columns
  for (col in c("Tumor_state", "cell_idents", "seurat_clusters")) {
    if (col %in% colnames(metadata)) {
      state_col <- col
      cat(sprintf("   Using fallback column: %s\n", state_col))
      break
    }
  }
}

# 5. Run differential PSI analysis
cat("\n5. Running differential PSI analysis...\n")

# Function to run Wilcoxon test for one event between two groups
run_wilcox_test <- function(psi_values, group1_cells, group2_cells) {
  g1_vals <- psi_values[names(psi_values) %in% group1_cells]
  g2_vals <- psi_values[names(psi_values) %in% group2_cells]

  # Remove NAs
  g1_vals <- g1_vals[!is.na(g1_vals)]
  g2_vals <- g2_vals[!is.na(g2_vals)]

  if (length(g1_vals) < 10 || length(g2_vals) < 10) {
    return(list(
      n_g1 = length(g1_vals),
      n_g2 = length(g2_vals),
      mean_g1 = ifelse(length(g1_vals) > 0, mean(g1_vals), NA),
      mean_g2 = ifelse(length(g2_vals) > 0, mean(g2_vals), NA),
      delta_psi = NA,
      pvalue = NA
    ))
  }

  test_result <- tryCatch(
    wilcox.test(g1_vals, g2_vals),
    error = function(e) list(p.value = NA)
  )

  return(list(
    n_g1 = length(g1_vals),
    n_g2 = length(g2_vals),
    mean_g1 = mean(g1_vals),
    mean_g2 = mean(g2_vals),
    delta_psi = mean(g2_vals) - mean(g1_vals),
    pvalue = test_result$p.value
  ))
}

# Get cells by tumor state
tumor_states <- unique(metadata[[state_col]])
tumor_states <- tumor_states[!is.na(tumor_states)]

# Create list of cells per state
cells_by_state <- list()
for (state in tumor_states) {
  state_cells <- rownames(metadata)[metadata[[state_col]] == state & !is.na(metadata[[state_col]])]
  # Only keep cells that are in PSI matrix
  state_cells <- state_cells[state_cells %in% colnames(psi_candidates)]
  cells_by_state[[state]] <- state_cells
  cat(sprintf("   State '%s': %d cells with PSI data\n", state, length(state_cells)))
}

# Filter to states with enough cells
cells_by_state <- cells_by_state[sapply(cells_by_state, length) >= 20]
tumor_states <- names(cells_by_state)
cat(sprintf("\n   States with ≥20 cells: %d\n", length(tumor_states)))

# Run pairwise comparisons
all_results <- list()
comparisons <- combn(tumor_states, 2, simplify = FALSE)

cat(sprintf("   Running %d pairwise comparisons for %d events...\n",
            length(comparisons), nrow(psi_candidates)))

pb_interval <- max(1, floor(nrow(psi_candidates) / 10))

for (i in seq_len(nrow(psi_candidates))) {
  event_id <- rownames(psi_candidates)[i]
  psi_vals <- as.numeric(psi_candidates[i, ])
  names(psi_vals) <- colnames(psi_candidates)

  # Get gene name for this event
  gene_name <- candidate_events$gene_short_name[candidate_events$id == event_id][1]
  event_type <- candidate_events$event_type[candidate_events$id == event_id][1]

  for (comp in comparisons) {
    state1 <- comp[1]
    state2 <- comp[2]

    result <- run_wilcox_test(psi_vals, cells_by_state[[state1]], cells_by_state[[state2]])

    all_results[[length(all_results) + 1]] <- data.frame(
      event_id = event_id,
      gene = gene_name,
      event_type = event_type,
      state1 = state1,
      state2 = state2,
      n_cells_state1 = result$n_g1,
      n_cells_state2 = result$n_g2,
      mean_psi_state1 = result$mean_g1,
      mean_psi_state2 = result$mean_g2,
      delta_psi = result$delta_psi,
      pvalue = result$pvalue,
      stringsAsFactors = FALSE
    )
  }

  if (i %% pb_interval == 0) {
    cat(sprintf("   Processed %d/%d events (%.0f%%)\n", i, nrow(psi_candidates), 100*i/nrow(psi_candidates)))
  }
}

# Combine results
results_df <- bind_rows(all_results)
cat(sprintf("\n   Total comparisons: %d\n", nrow(results_df)))

# Remove NA p-values before FDR correction
results_valid <- results_df %>% filter(!is.na(pvalue))
cat(sprintf("   Valid comparisons (with p-values): %d\n", nrow(results_valid)))

# Add FDR correction
results_valid$padj <- p.adjust(results_valid$pvalue, method = "BH")

# Sort by p-value
results_valid <- results_valid %>% arrange(pvalue)

# Save full results
output_file <- file.path(output_dir, "scrna_full_psi_differential.csv")
write.csv(results_valid, output_file, row.names = FALSE)
cat(sprintf("\n   Saved full results: %s (%d rows)\n", basename(output_file), nrow(results_valid)))

# 6. Summarize results
cat("\n================================================================\n")
cat("6. SUMMARY OF RESULTS\n")
cat("================================================================\n")

sig_results <- results_valid %>%
  filter(padj < 0.05, abs(delta_psi) > 10)

cat(sprintf("\nSignificant events (padj < 0.05, |ΔPSI| > 10%%): %d\n", nrow(sig_results)))

if (nrow(sig_results) > 0) {
  cat("\nTop 30 significant events:\n")
  print(head(sig_results %>% select(gene, event_type, state1, state2, delta_psi, pvalue, padj), 30))

  # Summarize by gene
  gene_summary <- sig_results %>%
    group_by(gene) %>%
    summarise(
      n_sig_events = n(),
      n_comparisons = n_distinct(paste(state1, state2)),
      max_delta_psi = max(abs(delta_psi)),
      min_pvalue = min(pvalue),
      min_padj = min(padj),
      .groups = 'drop'
    ) %>%
    arrange(min_pvalue)

  cat("\nGenes with significant differential splicing:\n")
  print(gene_summary)

  # Save gene summary
  gene_summary_file <- file.path(output_dir, "scrna_validation_gene_summary.csv")
  write.csv(gene_summary, gene_summary_file, row.names = FALSE)
  cat(sprintf("\nSaved: %s\n", basename(gene_summary_file)))
}

# 7. All genes - best result per gene
cat("\n================================================================\n")
cat("7. ALL GENES - BEST RESULT PER GENE\n")
cat("================================================================\n")

all_gene_summary <- results_valid %>%
  group_by(gene) %>%
  summarise(
    n_events_tested = n_distinct(event_id),
    n_comparisons = n(),
    best_delta_psi = delta_psi[which.min(pvalue)],
    min_pvalue = min(pvalue),
    min_padj = min(padj),
    n_sig_padj05 = sum(padj < 0.05, na.rm = TRUE),
    n_sig_with_effect = sum(padj < 0.05 & abs(delta_psi) > 10, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    validation_status = case_when(
      n_sig_with_effect >= 3 ~ "STRONG",
      n_sig_with_effect >= 1 ~ "MODERATE",
      n_sig_padj05 >= 1 ~ "WEAK",
      TRUE ~ "NOT_VALIDATED"
    )
  ) %>%
  arrange(min_pvalue)

print(all_gene_summary, n = 30)

# Save all gene summary
all_gene_file <- file.path(output_dir, "scrna_validation_all_genes.csv")
write.csv(all_gene_summary, all_gene_file, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", basename(all_gene_file)))

# Final summary
cat("\n================================================================\n")
cat("VALIDATION SUMMARY\n")
cat("================================================================\n")

validation_counts <- table(all_gene_summary$validation_status)
cat("\nValidation status counts:\n")
print(validation_counts)

strong_genes <- all_gene_summary$gene[all_gene_summary$validation_status == "STRONG"]
moderate_genes <- all_gene_summary$gene[all_gene_summary$validation_status == "MODERATE"]

cat(sprintf("\nSTRONG validation (≥3 sig events with |ΔPSI|>10%%): %s\n",
            ifelse(length(strong_genes) > 0, paste(strong_genes, collapse=", "), "None")))
cat(sprintf("MODERATE validation (1-2 sig events with |ΔPSI|>10%%): %s\n",
            ifelse(length(moderate_genes) > 0, paste(moderate_genes, collapse=", "), "None")))

cat("\n================================================================\n")
cat("Analysis complete!\n")
cat("================================================================\n")
