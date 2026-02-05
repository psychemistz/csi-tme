#!/usr/bin/env Rscript
# Create boxplots for APP and VCAN PSI by tumor state

library(dplyr)
library(ggplot2)
library(tidyr)

# Load PSI data
cat("Loading PSI data...\n")
load("/vf/users/parks34/projects/3csitme/data/zenodo_17055113/zenodo-20250904/MARVEL/PSI_merge.RData")

# Load metadata
load("/vf/users/parks34/projects/3csitme/data/zenodo_17055113/zenodo-20250904/seurat/metadata_seu_filt_all.RData")

# Load differential results to get gene-event mapping
diff_results <- read.csv("/vf/users/parks34/projects/3csitme/reports/scrna_full_psi_differential.csv")

output_dir <- "/vf/users/parks34/projects/3csitme/reports"

cat(sprintf("PSI matrix: %d events x %d cells\n", nrow(PSI_merge), ncol(PSI_merge)))
cat(sprintf("Metadata: %d cells\n", nrow(metadata)))

# Get state column
state_col <- "cell_idents_tumor"
cat("Cell states:\n")
print(table(metadata[[state_col]]))

# Get event names
event_names <- rownames(PSI_merge)
cat(sprintf("\nTotal events: %d\n", length(event_names)))
cat("Example event names:\n")
print(head(event_names, 5))

# Find APP events from differential results
app_events <- unique(diff_results$event_id[diff_results$gene == "APP"])
cat(sprintf("\nAPP events from diff results: %d\n", length(app_events)))
app_events_in_psi <- intersect(app_events, event_names)
cat(sprintf("APP events in PSI matrix: %d\n", length(app_events_in_psi)))

# Find VCAN events from differential results
vcan_events <- unique(diff_results$event_id[diff_results$gene == "VCAN"])
cat(sprintf("VCAN events from diff results: %d\n", length(vcan_events)))
vcan_events_in_psi <- intersect(vcan_events, event_names)
cat(sprintf("VCAN events in PSI matrix: %d\n", length(vcan_events_in_psi)))

# Function to create boxplot
create_psi_boxplot <- function(psi_mat, gene_events, gene_name, meta, state_col, output_file, title) {
  if (length(gene_events) == 0) {
    cat(sprintf("No events found for %s\n", gene_name))
    return(FALSE)
  }

  # Use the first event with most non-NA values
  best_event <- NULL
  max_cells <- 0

  for (event in gene_events) {
    if (event %in% rownames(psi_mat)) {
      event_row <- as.numeric(psi_mat[event, ])
      n_valid <- sum(!is.na(event_row))
      if (n_valid > max_cells) {
        max_cells <- n_valid
        best_event <- event
      }
    }
  }

  if (is.null(best_event)) {
    cat(sprintf("No valid events for %s\n", gene_name))
    return(FALSE)
  }

  cat(sprintf("Using event: %s (%d cells with data)\n", best_event, max_cells))

  # Extract PSI values - convert data.frame row to named numeric vector
  psi_row <- psi_mat[best_event, ]
  psi_values <- as.numeric(psi_row)
  names(psi_values) <- colnames(psi_mat)

  # Match cells with metadata
  common_cells <- intersect(names(psi_values), rownames(meta))
  cat(sprintf("Common cells: %d\n", length(common_cells)))

  if (length(common_cells) < 100) {
    cat("Too few common cells\n")
    return(FALSE)
  }

  # Create data frame
  plot_data <- data.frame(
    cell = common_cells,
    psi = psi_values[common_cells],
    state = meta[common_cells, state_col]
  )

  # Remove NA PSI values
  plot_data <- plot_data[!is.na(plot_data$psi), ]
  cat(sprintf("Cells with PSI data: %d\n", nrow(plot_data)))

  # Filter to tumor states only (containing "Tumor")
  tumor_data <- plot_data[grepl("Tumor", plot_data$state), ]
  cat(sprintf("Tumor cells: %d\n", nrow(tumor_data)))

  if (nrow(tumor_data) < 50) {
    cat("Too few tumor cells, including Macrophage for comparison\n")
    tumor_data <- plot_data[grepl("Tumor|Macrophage", plot_data$state), ]
  }

  # Calculate summary stats per state
  summary_stats <- tumor_data %>%
    group_by(state) %>%
    summarise(
      n = n(),
      mean_psi = mean(psi * 100),
      median_psi = median(psi * 100),
      sd_psi = sd(psi * 100)
    ) %>%
    filter(n >= 10) %>%
    arrange(desc(mean_psi))

  cat("\nPSI by state (%):\n")
  print(as.data.frame(summary_stats))

  # Filter to states with enough cells
  valid_states <- summary_stats$state
  tumor_data <- tumor_data[tumor_data$state %in% valid_states, ]

  # Reorder states by mean PSI
  tumor_data$state <- factor(tumor_data$state, levels = summary_stats$state)

  # Color palette - red for high PSI, blue for low
  n_states <- length(valid_states)
  colors <- colorRampPalette(c("#e74c3c", "#f39c12", "#3498db"))(n_states)

  # Create boxplot
  p <- ggplot(tumor_data, aes(x = state, y = psi * 100, fill = state)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    geom_jitter(width = 0.2, alpha = 0.05, size = 0.3) +
    labs(
      title = title,
      x = "Cell State",
      y = "PSI (%)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(values = colors) +
    ylim(0, 100)

  # Add sample size labels
  n_labels <- tumor_data %>%
    group_by(state) %>%
    summarise(n = n(), y_pos = 3)

  p <- p + geom_text(data = n_labels,
                     aes(x = state, y = y_pos, label = paste0("n=", n)),
                     size = 3, color = "gray30")

  ggsave(output_file, p, width = 9, height = 6, dpi = 150)
  cat(sprintf("Saved: %s\n", output_file))

  return(TRUE)
}

# Create APP boxplot
cat("\n========== Creating APP boxplot ==========\n")
create_psi_boxplot(
  PSI_merge, app_events_in_psi, "APP", metadata, state_col,
  file.path(output_dir, "fig1_app_psi_boxplot.png"),
  "APP Exon 7/8 Inclusion (PSI) by Tumor State"
)

# Create VCAN boxplot
cat("\n========== Creating VCAN boxplot ==========\n")
create_psi_boxplot(
  PSI_merge, vcan_events_in_psi, "VCAN", metadata, state_col,
  file.path(output_dir, "fig2_vcan_psi_boxplot.png"),
  "VCAN GAG-Domain Exon Inclusion (PSI) by Tumor State"
)

cat("\n========== Done ==========\n")
