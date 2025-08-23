# pca script for sample identification
# Usage: Rscript pca_analysis.R <directory_path>
# The script will look for pca_data.csv in the specified directory
# Samples with Type="unknown" will be identified against samples with Type="ref"
# All output files will be saved in the same directory as the input data

# Function to install and load packages (cross-platform compatible)
install_and_load <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    # Set CRAN mirror for automatic installation
    if (length(getOption("repos")) == 0 || getOption("repos")["CRAN"] == "@CRAN@") {
      options(repos = c(CRAN = "https://cloud.r-project.org/"))
    }
    
    message("Installing package: ", package)
    install.packages(package, quiet = TRUE, dependencies = TRUE)
    
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
      stop("Failed to install required package: ", package)
    }
  }
}

# Load required libraries (suppress messages)
suppressPackageStartupMessages({
  install_and_load("ggplot2")
  install_and_load("gridExtra")
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set default path or use provided path (cross-platform compatible)
if (length(args) > 0) {
  # If argument is provided, treat it as a directory path
  input_dir <- normalizePath(args[1], mustWork = FALSE)
  
  # Check if it's a directory or a file
  if (dir.exists(input_dir)) {
    # It's a directory, look for pca_data.csv in it
    data_path <- file.path(input_dir, "pca_data.csv")
    output_dir <- input_dir
  } else if (file.exists(input_dir)) {
    # It's a file path
    data_path <- input_dir
    output_dir <- dirname(data_path)
  } else {
    # Assume it's a directory path that might not exist yet
    data_path <- file.path(input_dir, "pca_data.csv")
    output_dir <- input_dir
  }
} else {
  # No argument provided, use current directory
  data_path <- file.path(getwd(), "pca_data.csv")
  output_dir <- getwd()
}

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Check if file exists (cross-platform)
if (!file.exists(data_path)) {
  stop("Data file not found: ", data_path, "\nCurrent working directory: ", getwd())
}

# Print system information for debugging
message("System information:")
message("Platform: ", R.version$platform)
message("R version: ", R.version.string)
message("Working directory: ", getwd())
message("Data file path: ", data_path)
message("Output directory: ", output_dir)

# Read the data
data_all <- read.csv(data_path, stringsAsFactors = FALSE)

# Separate unknown samples from reference samples based on Type column
unknown_samples <- data_all[data_all$Type == "unknown", ]
reference_samples <- data_all[data_all$Type == "ref", ]

# Debug information
message("Total rows in data: ", nrow(data_all))
message("Unknown samples found: ", nrow(unknown_samples))
message("Reference samples found: ", nrow(reference_samples))
if (nrow(unknown_samples) == 0) {
  stop("No unknown samples found. Check that Type column contains 'unknown' values.")
}
if (nrow(reference_samples) == 0) {
  stop("No reference samples found. Check that Type column contains 'ref' values.")
}

# Prepare data for PCA (remove sample names and type column)
reference_numeric <- reference_samples[, -c(1, 2)]  # Remove Sample and Type columns
unknown_numeric <- unknown_samples[, -c(1, 2)]      # Remove Sample and Type columns

# Clean column names (remove X prefix that R adds to numeric column names)
colnames(reference_numeric) <- gsub("^X", "", colnames(reference_numeric))
colnames(unknown_numeric) <- gsub("^X", "", colnames(unknown_numeric))

# Check for zero variance elements
element_variances <- apply(reference_numeric, 2, var, na.rm = TRUE)
zero_var_elements <- names(element_variances[element_variances == 0 | is.na(element_variances)])
near_zero_var_elements <- names(element_variances[element_variances < 1e-10 & element_variances > 0])

excluded_elements <- c(zero_var_elements, near_zero_var_elements)
if (length(excluded_elements) > 0) {
  message("Excluded elements due to zero/near-zero variance: ", length(excluded_elements))
  reference_numeric <- reference_numeric[, !names(reference_numeric) %in% excluded_elements]
  unknown_numeric <- unknown_numeric[, !names(unknown_numeric) %in% excluded_elements]
}

# Perform PCA on reference samples
pca_result <- prcomp(reference_numeric, scale. = TRUE, center = TRUE)

# Calculate variance explained
variance_explained <- summary(pca_result)$importance[2, ] * 100
cumulative_variance <- summary(pca_result)$importance[3, ] * 100

# Project unknown samples into PCA space
unknown_scaled <- scale(unknown_numeric, 
                       center = pca_result$center, 
                       scale = pca_result$scale)
unknown_pca <- unknown_scaled %*% pca_result$rotation

# Determine number of PCs needed for >70% cumulative variance (moved up for dynamic calculation)
pcs_needed <- which(cumulative_variance >= 70)[1]
if (is.na(pcs_needed)) pcs_needed <- min(3, length(cumulative_variance))

message("Number of PCs needed for >70% variance: ", pcs_needed)
message("Cumulative variance with ", pcs_needed, " PCs: ", round(cumulative_variance[pcs_needed], 1), " %")

# Create dynamic reference scores dataframe based on required PCs
ref_scores <- data.frame(Sample = reference_samples$Sample)
for (pc in 1:pcs_needed) {
  ref_scores[[paste0("PC", pc)]] <- pca_result$x[, pc]
}

# Calculate position for each unknown sample using dynamic number of PCs
unknown_results <- data.frame()

for (i in 1:nrow(unknown_samples)) {
  sample_id <- unknown_samples$Sample[i]
  
  # Get sample position for required number of PCs
  sample_coords <- unknown_pca[i, 1:pcs_needed]
  
  # Calculate distances to all reference samples using dynamic number of PCs
  distances <- numeric(nrow(ref_scores))
  for (j in 1:nrow(ref_scores)) {
    ref_coords <- as.numeric(ref_scores[j, 2:(pcs_needed+1)])  # Skip Sample column
    distances[j] <- sqrt(sum((sample_coords - ref_coords)^2))
  }
  
  # Find closest match
  closest_idx <- which.min(distances)
  closest_sample <- ref_scores$Sample[closest_idx]
  closest_distance <- distances[closest_idx]
  
  # Calculate correlation with closest sample
  unknown_composition <- as.numeric(unknown_numeric[i, ])
  closest_ref_idx <- which(reference_samples$Sample == closest_sample)[1]
  closest_composition <- as.numeric(reference_numeric[closest_ref_idx, ])
  
  correlation_coef <- cor(unknown_composition, closest_composition)
  
  # Create results dataframe with dynamic PC columns
  result_row <- data.frame(
    Unknown_Sample = sample_id,
    Closest_Match = closest_sample,
    Distance = closest_distance,
    Correlation = correlation_coef
  )
  
  # Add PC coordinates dynamically
  for (pc in 1:pcs_needed) {
    result_row[[paste0("PC", pc)]] <- sample_coords[pc]
  }
  
  unknown_results <- rbind(unknown_results, result_row)
}

# No threshold calculations - pure distance analysis

# Export intermediate files for alternative analysis script
intermediate_path <- file.path(output_dir, "pca_intermediate.RData")
save(reference_numeric, unknown_numeric, pca_result, unknown_results, 
     reference_samples, unknown_samples, variance_explained, cumulative_variance,
     excluded_elements, ref_scores, unknown_pca,
     file = intermediate_path)
message("Intermediate data saved to: ", intermediate_path)

# Store results for output files (no console output)

# Create comprehensive PDF report (cross-platform file paths)
pdf_path <- file.path(output_dir, "pca_report.pdf")
suppressWarnings({
pdf(pdf_path, width = 11, height = 8.5)

# Page 1: Summary Statistics and Results Table
# pcs_needed already calculated above for dynamic distance calculation

# Create summary information plot
pc_info <- paste0("PC", 1:pcs_needed, ": ", round(variance_explained[1:pcs_needed], 1), "%", collapse = " | ")
cumulative_info <- round(cumulative_variance[pcs_needed], 1)

# Create summary plot
summary_plot <- ggplot() + 
  theme_void() +
  ggtitle("PCA Sample Identification - Summary Results") +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) +
  annotate("text", x = 0.05, y = 0.8, 
           label = paste0("ANALYSIS SUMMARY:\n",
                         "• Dataset: ", nrow(data_all), " samples (", nrow(unknown_samples), " unknowns + ", nrow(reference_samples), " references) | ", ncol(reference_numeric), " elements analyzed\n",
                         "• Principal Components: ", pc_info, " | Cumulative: ", cumulative_info, "% (≥70% threshold)\n",
                         "• Method: Distance-based identification using ", pcs_needed, " principal components"),
           hjust = 0, vjust = 1, size = 3.8, family = "mono") +
  xlim(0, 1) + ylim(0, 1)

# Create results table data
results_summary <- unknown_results
results_summary$Distance <- round(results_summary$Distance, 3)
results_summary$Correlation <- round(results_summary$Correlation, 4)
results_summary$PC1 <- round(results_summary$PC1, 2)
results_summary$PC2 <- round(results_summary$PC2, 2)
results_summary$PC3 <- round(results_summary$PC3, 2)

# Sort by distance for ranking
results_sorted <- results_summary[order(results_summary$Distance), ]
results_sorted$Rank <- 1:nrow(results_sorted)

# Create simplified table data (only Distance and Correlation)
table_data <- data.frame(
  Sample = results_sorted$Unknown_Sample,
  Match = results_sorted$Closest_Match,
  Distance = results_sorted$Distance,
  Correlation = results_sorted$Correlation,
  Rank = results_sorted$Rank
)

# Create a clean table visualization
table_plot <- ggplot() + 
  theme_void() +
  annotate("text", x = 0.05, y = 0.95, 
           label = "IDENTIFICATION RESULTS",
           hjust = 0, vjust = 1, size = 4.5, fontface = "bold") +
  annotate("text", x = 0.05, y = 0.85, 
           label = sprintf("%-10s %-15s %-12s %-12s %-8s", 
                          "Sample", "Closest Match", "Distance", "Correlation", "Rank"),
           hjust = 0, vjust = 1, size = 3.5, family = "mono", fontface = "bold") +
  annotate("text", x = 0.05, y = 0.80, 
           label = paste(rep("-", 65), collapse = ""),
           hjust = 0, vjust = 1, size = 3, family = "mono") +
  xlim(0, 1) + ylim(0, 1)

# Add table rows with better spacing
for (i in 1:nrow(table_data)) {
  y_pos <- 0.72 - (i-1) * 0.08
  table_plot <- table_plot + 
    annotate("text", x = 0.05, y = y_pos,
             label = sprintf("%-10s %-15s %-12.3f %-12.4f %-8s",
                            table_data$Sample[i],
                            table_data$Match[i],
                            table_data$Distance[i],
                            table_data$Correlation[i],
                            paste0("#", table_data$Rank[i])),
             hjust = 0, vjust = 1, size = 3.2, family = "mono")
}

# Combine plots using grid.arrange
p1 <- grid.arrange(summary_plot, table_plot, ncol = 1, heights = c(1, 2))

# Export comprehensive summary table to separate text file (cross-platform)
summary_path <- file.path(output_dir, "pca_summary_table.txt")
sink(summary_path)

writeLines("pca sample identification - comprehensive summary table")
writeLines(paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
writeLines("Data file: pca_data.csv")
writeLines("")

# Recreate simplified table_text for file output
table_text <- "SAMPLE IDENTIFICATION RESULTS\n"
table_text <- paste0(table_text, paste(rep("=", 65), collapse = ""), "\n")
table_text <- paste0(table_text, sprintf("%-10s %-15s %-12s %-12s %-8s\n", 
                                        "Sample", "Closest_Match", "Distance", "Correlation", "Rank"))
table_text <- paste0(table_text, paste(rep("-", 65), collapse = ""), "\n")

for (i in 1:nrow(results_sorted)) {
  table_text <- paste0(table_text, sprintf("%-10s %-15s %-12.3f %-12.4f %-8s\n",
                                          results_sorted$Unknown_Sample[i],
                                          results_sorted$Closest_Match[i],
                                          results_sorted$Distance[i],
                                          results_sorted$Correlation[i],
                                          paste0("#", i)))
}

table_text <- paste0(table_text, paste(rep("-", 65), collapse = ""), "\n")
table_text <- paste0(table_text, "ANALYSIS SUMMARY:\n")
table_text <- paste0(table_text, sprintf("• Dataset: %d samples (%d unknowns + %d references) | %d elements\n", 
                                        nrow(data_all), nrow(unknown_samples), nrow(reference_samples), ncol(reference_numeric)))
table_text <- paste0(table_text, sprintf("• Principal Components: %s | Cumulative: %.1f%% (≥70%% threshold)\n", 
                                        pc_info, cumulative_info))

writeLines(table_text)

# Create comprehensive distance matrix for all unknown samples to all references
writeLines("")
writeLines("")
writeLines("comprehensive distance matrix:")
writeLines("Unknown samples (rows) vs Reference samples (columns)")
writeLines(paste(rep("=", 100), collapse = ""))

# Get all reference sample names
ref_sample_names <- unique(reference_samples$Sample)
ref_sample_names <- sort(ref_sample_names)

# Create distance matrix
distance_matrix <- matrix(0, nrow = nrow(unknown_samples), ncol = length(ref_sample_names))
rownames(distance_matrix) <- paste("Sample", unknown_samples$Sample)
colnames(distance_matrix) <- ref_sample_names

# Calculate distances for each unknown sample to each reference using dynamic PCs
for (i in 1:nrow(unknown_samples)) {
  sample_id <- unknown_samples$Sample[i]
  
  # Get sample coordinates for all required PCs
  sample_coords <- numeric(pcs_needed)
  for (pc in 1:pcs_needed) {
    sample_coords[pc] <- unknown_results[[paste0("PC", pc)]][unknown_results$Unknown_Sample == sample_id]
  }
  
  for (ref_name in ref_sample_names) {
    ref_indices <- which(reference_samples$Sample == ref_name)
    if (length(ref_indices) > 0) {
      # Calculate mean position if multiple replicates for all required PCs
      ref_coords <- numeric(pcs_needed)
      for (pc in 1:pcs_needed) {
        ref_coords[pc] <- mean(pca_result$x[ref_indices, pc])
      }
      
      # Calculate distance using all required PCs
      distance <- sqrt(sum((sample_coords - ref_coords)^2))
      distance_matrix[i, ref_name] <- round(distance, 3)
    }
  }
}

# Print distance matrix header
header_line <- sprintf("%-10s", "Sample")
for (ref_name in ref_sample_names) {
  header_line <- paste0(header_line, sprintf("%-10s", substr(ref_name, 1, 9)))
}
writeLines(header_line)
writeLines(paste(rep("-", nchar(header_line)), collapse = ""))

# Print distance matrix rows
for (i in 1:nrow(unknown_samples)) {
  sample_id <- unknown_samples$Sample[i]
  row_line <- sprintf("%-10s", paste("Sample", sample_id))
  for (ref_name in ref_sample_names) {
    row_line <- paste0(row_line, sprintf("%-10s", distance_matrix[i, ref_name]))
  }
  writeLines(row_line)
}

writeLines("")
writeLines("")
writeLines("detailed sample breakdown:")
writeLines(paste(rep("=", 75), collapse = ""))

for (i in 1:nrow(results_sorted)) {
  sample_data <- results_sorted[i, ]
  sample_id <- sample_data$Unknown_Sample
  
  writeLines(sprintf("sample %s:", sample_id))
  writeLines(sprintf("  • closest match: %s", sample_data$Closest_Match))
  writeLines(sprintf("  • distance: %.3f (rank #%d)", sample_data$Distance, i))
  writeLines(sprintf("  • correlation: %.4f", sample_data$Correlation))
  writeLines(sprintf("  • pca coordinates: pc1=%.2f, pc2=%.2f, pc3=%.2f", sample_data$PC1, sample_data$PC2, sample_data$PC3))
  
  # Show distances to all references for this sample
  writeLines("  • distances to all references:")
  sample_row_idx <- which(unknown_samples$Sample == sample_id)
  sample_distances <- distance_matrix[sample_row_idx, ]
  sorted_distances <- sort(sample_distances)
  for (j in 1:length(sorted_distances)) {
    ref_name <- names(sorted_distances)[j]
    distance_val <- sorted_distances[j]
    rank_indicator <- if (j == 1) " (closest)" else ""
    writeLines(sprintf("    %2d. %-12s: %.3f%s", j, ref_name, distance_val, rank_indicator))
  }
  writeLines("")
}

writeLines("reference samples in dataset:")
writeLines(paste(rep("-", 50), collapse = ""))
ref_samples <- unique(reference_samples$Sample)
for (i in 1:length(ref_samples)) {
  writeLines(sprintf("%2d. %s", i, ref_samples[i]))
}

writeLines("")
writeLines("match summary:")
writeLines(paste(rep("-", 30), collapse = ""))
match_counts <- table(results_summary$Closest_Match)
for (ref_name in names(match_counts)) {
  samples_matched <- results_summary$Unknown_Sample[results_summary$Closest_Match == ref_name]
  writeLines(sprintf("%s: %d samples (%s)", ref_name, match_counts[ref_name], paste(samples_matched, collapse = ", ")))
}

sink()

# p1 is already printed by grid.arrange

# Page 2: Scree Plot
scree_data <- data.frame(
  PC = 1:min(10, length(variance_explained)),
  Variance = variance_explained[1:min(10, length(variance_explained))]
)

p2 <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 3) +
  labs(title = "scree plot - variance explained by principal components",
       x = "principal component",
       y = "variance explained (%)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))

print(p2)

# Page 3: PCA Scores Plot with Unknown Projection (PC1 vs PC2 only)
scores_data <- data.frame(
  Sample = reference_samples$Sample,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)

# Add unknown samples to plot
unknown_plot_data <- data.frame(
  Sample = as.character(unknown_results$Unknown_Sample),
  PC1 = unknown_results$PC1,
  PC2 = unknown_results$PC2
)

p3 <- ggplot() +
  geom_point(data = scores_data, aes(x = PC1, y = PC2, color = Sample), 
             alpha = 0.7, size = 2) +
  geom_text(data = scores_data, aes(x = PC1, y = PC2, label = Sample),
            vjust = -0.8, hjust = 0.5, size = 2, color = "darkblue", alpha = 0.8) +
  geom_point(data = unknown_plot_data, aes(x = PC1, y = PC2), 
             color = "red", size = 2.5, shape = 17) +
  geom_text(data = unknown_plot_data, aes(x = PC1, y = PC2, label = Sample),
            vjust = -1.2, hjust = 0.5, size = 2.5, color = "black", fontface = "bold") +
  labs(title = "pca scores plot with unknown sample projections",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.position = "right")

print(p3)

# Dynamic PC plots based on cumulative variance >70% (already calculated above)

# Generate all PC combination plots dynamically
pc_combinations <- list()
if (pcs_needed >= 2) {
  for (i in 1:(pcs_needed-1)) {
    for (j in (i+1):pcs_needed) {
      pc_combinations[[length(pc_combinations) + 1]] <- c(i, j)
    }
  }
}

# Skip the first combination (PC1 vs PC2) as it's already plotted
if (length(pc_combinations) > 1) {
  for (k in 2:length(pc_combinations)) {
    pc_pair <- pc_combinations[[k]]
    pc_x <- pc_pair[1]
    pc_y <- pc_pair[2]
    
    # Create scores data for this PC combination
    scores_data_dynamic <- data.frame(
      Sample = reference_samples$Sample,
      PCx = pca_result$x[, pc_x],
      PCy = pca_result$x[, pc_y]
    )
    
    # Unknown samples data for this PC combination
    unknown_plot_data_dynamic <- data.frame(
      Sample = as.character(unknown_results$Unknown_Sample),
      PCx = unknown_pca[, pc_x],
      PCy = unknown_pca[, pc_y]
    )
    
    # Create plot
    p_dynamic <- ggplot() +
      geom_point(data = scores_data_dynamic, aes(x = PCx, y = PCy, color = Sample), 
                 alpha = 0.7, size = 2) +
      geom_text(data = scores_data_dynamic, aes(x = PCx, y = PCy, label = Sample),
                vjust = -0.8, hjust = 0.5, size = 2, color = "darkblue", alpha = 0.8) +
      geom_point(data = unknown_plot_data_dynamic, aes(x = PCx, y = PCy), 
                 color = "red", size = 2.5, shape = 17) +
      geom_text(data = unknown_plot_data_dynamic, aes(x = PCx, y = PCy, label = Sample),
                vjust = -1.2, hjust = 0.5, size = 2.5, color = "black", fontface = "bold") +
      labs(title = paste0("pca scores plot: PC", pc_x, " vs PC", pc_y, " with unknown sample projections"),
           x = paste0("PC", pc_x, " (", round(variance_explained[pc_x], 1), "%)"),
           y = paste0("PC", pc_y, " (", round(variance_explained[pc_y], 1), "%)")) +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
            legend.position = "right")
    
    print(p_dynamic)
  }
}

# Loadings Plot (page number depends on dynamic PC plots)
loadings_data <- data.frame(
  Element = rownames(pca_result$rotation),
  PC1 = pca_result$rotation[, 1],
  PC2 = pca_result$rotation[, 2],
  PC3 = pca_result$rotation[, 3]
)

# Select top contributors
loadings_data$Contribution <- sqrt(loadings_data$PC1^2 + loadings_data$PC2^2)
top_loadings <- loadings_data[order(loadings_data$Contribution, decreasing = TRUE)[1:20], ]

p6 <- ggplot(top_loadings, aes(x = PC1, y = PC2)) +
  geom_point(color = "blue", size = 3) +
  geom_text(aes(label = Element), vjust = -0.5, hjust = 0.5, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "pca loadings plot - top 20 contributing elements",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))

print(p6)

# Distance Matrix Visualization (page number depends on dynamic PC plots)
# Create comprehensive distance data for all unknown-reference combinations
all_distance_data <- data.frame()

# Get all reference sample names
ref_sample_names <- unique(reference_samples$Sample)

# Calculate all distances using dynamic number of PCs
for (i in 1:nrow(unknown_samples)) {
  sample_id <- unknown_samples$Sample[i]
  
  # Get sample coordinates for all required PCs
  sample_coords <- numeric(pcs_needed)
  for (pc in 1:pcs_needed) {
    sample_coords[pc] <- unknown_results[[paste0("PC", pc)]][unknown_results$Unknown_Sample == sample_id]
  }
  
  for (ref_name in ref_sample_names) {
    ref_indices <- which(reference_samples$Sample == ref_name)
    if (length(ref_indices) > 0) {
      # Calculate mean position if multiple replicates for all required PCs
      ref_coords <- numeric(pcs_needed)
      for (pc in 1:pcs_needed) {
        ref_coords[pc] <- mean(pca_result$x[ref_indices, pc])
      }
      
      # Calculate distance using all required PCs
      distance <- sqrt(sum((sample_coords - ref_coords)^2))
      
      all_distance_data <- rbind(all_distance_data, data.frame(
        Unknown_Sample = paste("Sample", sample_id),
        Reference_Sample = ref_name,
        Distance = distance
      ))
    }
  }
}

# Create heatmap-style visualization
p7 <- ggplot(all_distance_data, aes(x = Reference_Sample, y = Unknown_Sample, fill = Distance)) +
  geom_tile(alpha = 0.8) +
  scale_fill_gradient(low = "blue", high = "red", name = "Distance") +
  geom_text(aes(label = round(Distance, 2)), color = "white", size = 2.5) +
  labs(title = "distance matrix - all unknown samples to all reference samples",
       x = "reference sample",
       y = "unknown sample") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10))

print(p7)

# Individual Analysis for all unknown samples
# Calculate element variability across all reference samples
element_variability <- apply(reference_numeric, 2, var, na.rm = TRUE)
top_variable_elements <- names(sort(element_variability, decreasing = TRUE)[1:20])

# Generate plots for all unknown samples
for (i in 1:nrow(unknown_samples)) {
  sample_id <- unknown_samples$Sample[i]
  # Get composition of unknown sample
  unknown_composition <- as.numeric(unknown_numeric[i, ])
    names(unknown_composition) <- names(unknown_numeric)
    
    # Get closest match composition
    closest_sample <- unknown_results$Closest_Match[unknown_results$Unknown_Sample == sample_id]
    closest_ref_idx <- which(reference_samples$Sample == closest_sample)[1]
    closest_composition <- as.numeric(reference_numeric[closest_ref_idx, ])
    names(closest_composition) <- names(reference_numeric)
    
    # Page A: Correlation Plot for this sample
    correlation_data <- data.frame(
      Unknown = unknown_composition,
      Reference = closest_composition,
      Element = names(unknown_composition)
    )
    
    correlation_coef <- cor(unknown_composition, closest_composition)
    
    # Identify outlier elements for labeling
    correlation_data$residual <- abs(correlation_data$Unknown - correlation_data$Reference)
    top_outliers <- correlation_data[order(correlation_data$residual, decreasing = TRUE)[1:5], ]
    
    # Add small constant to avoid log(0) issues and handle negative values
    correlation_data$Reference_log <- log10(abs(correlation_data$Reference) + 1)
    correlation_data$Unknown_log <- log10(abs(correlation_data$Unknown) + 1)
    top_outliers$Reference_log <- log10(abs(top_outliers$Reference) + 1)
    top_outliers$Unknown_log <- log10(abs(top_outliers$Unknown) + 1)
    
    p_correlation <- ggplot(correlation_data, aes(x = Reference_log, y = Unknown_log)) +
      geom_point(alpha = 0.6, size = 1.5, color = "blue") +
      geom_smooth(method = "lm", color = "red", se = TRUE, alpha = 0.3) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray", alpha = 0.7) +
      geom_text(data = top_outliers, aes(x = Reference_log, y = Unknown_log, label = Element), 
                vjust = -0.5, hjust = 0.5, size = 2.5, color = "darkred") +
      scale_x_continuous(name = paste("reference sample:", closest_sample, "concentration (log10)")) +
      scale_y_continuous(name = paste("unknown sample", sample_id, "concentration (log10)")) +
      labs(title = paste("element correlation - sample", sample_id, "vs", closest_sample),
           subtitle = paste("all", length(correlation_data$Element), "elements shown as dots (log scale) | correlation: r =", round(correlation_coef, 4))) +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(size = 10, hjust = 0.5))
    
    # Add labels to all elements in the correlation plot
    p_correlation <- ggplot(correlation_data, aes(x = Reference_log, y = Unknown_log)) +
      geom_point(alpha = 0.6, size = 1.5, color = "blue") +
      geom_smooth(method = "lm", color = "red", se = TRUE, alpha = 0.3) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray", alpha = 0.7) +
      geom_text(aes(label = Element), vjust = -0.5, hjust = 0.5, size = 1.8, color = "black", alpha = 0.8) +
      scale_x_continuous(name = paste("reference sample:", closest_sample, "concentration (log10)")) +
      scale_y_continuous(name = paste("unknown sample", sample_id, "concentration (log10)")) +
      labs(title = paste("element correlation - sample", sample_id, "vs", closest_sample),
           subtitle = paste("all", length(correlation_data$Element), "elements labeled (log scale) | correlation: r =", round(correlation_coef, 4))) +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(size = 10, hjust = 0.5))
    
    print(p_correlation)
    
    # Page B: Distance Ranking for this sample
    # Calculate distances from this sample to all reference samples using dynamic PCs
    distances_to_refs <- numeric(nrow(ref_scores))
    
    # Get unknown sample coordinates for all required PCs
    unknown_coords <- numeric(pcs_needed)
    for (pc in 1:pcs_needed) {
      unknown_coords[pc] <- unknown_results[[paste0("PC", pc)]][unknown_results$Unknown_Sample == sample_id]
    }
    
    # Calculate distance to each reference sample
    for (j in 1:nrow(ref_scores)) {
      ref_coords <- as.numeric(ref_scores[j, 2:(pcs_needed+1)])  # Skip Sample column
      distances_to_refs[j] <- sqrt(sum((unknown_coords - ref_coords)^2))
    }
    
    sample_distances <- data.frame(
      Reference_Sample = ref_scores$Sample,
      Distance = distances_to_refs
    )
    
    # Sort by distance and assign quality
    sample_distances <- sample_distances[order(sample_distances$Distance), ]
    sample_distances$Rank <- 1:nrow(sample_distances)
    

    
    # Create ordered factor for x-axis to ensure proper sorting (remove duplicates)
    ordered_levels <- unique(sample_distances$Reference_Sample[order(sample_distances$Distance)])
    sample_distances$Reference_Sample_Ordered <- factor(sample_distances$Reference_Sample, 
                                                       levels = ordered_levels)
    
    p_distance <- ggplot(sample_distances, aes(x = Reference_Sample_Ordered, y = Distance)) +
      geom_col(fill = "steelblue", alpha = 0.8) +
      labs(title = paste("distance ranking - sample", sample_id, "to all reference samples"),
           subtitle = paste("closest match:", closest_sample, "| distance:", 
                           round(min(sample_distances$Distance), 3)),
           x = "reference sample",
           y = "distance in pca space") +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(size = 10, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "top")
    
    print(p_distance)
    
    # Page C: Element Profile Comparison for this sample
    comparison_data <- data.frame(
      Element = rep(top_variable_elements, 2),
      Concentration = c(unknown_composition[top_variable_elements],
                       closest_composition[top_variable_elements]),
      Sample_Type = rep(c(paste("Sample", sample_id), closest_sample), each = 20)
    )
    
    # Add small constant to avoid log(0) issues and apply log transformation
    comparison_data$Concentration_log <- log10(comparison_data$Concentration + 1)
    
    p_profile <- ggplot(comparison_data, aes(x = Element, y = Concentration_log, fill = Sample_Type)) +
      geom_col(position = "dodge", alpha = 0.8) +
      scale_fill_manual(values = c("red", "green")) +
      labs(title = paste("element profile comparison - sample", sample_id, "vs", closest_sample),
           subtitle = paste("top 20 most variable elements (log scale) | distance:", round(unknown_results$Distance[unknown_results$Unknown_Sample == sample_id], 3),
                           "| correlation:", round(unknown_results$Correlation[unknown_results$Unknown_Sample == sample_id], 4)),
           x = "element (top 20 most variable)",
           y = "concentration (log10)",
           fill = "sample") +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(size = 10, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "top")
    
    print(p_profile)
  }

dev.off()
})

# Write detailed results to text file (cross-platform)
results_path <- file.path(output_dir, "pca_results.txt")
sink(results_path)

writeLines("=== complete pca results ===")
writeLines(paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
writeLines("")

writeLines("1. input file:")
writeLines("data file: pca_data.csv")
writeLines("")

writeLines("2. dataset information:")
writeLines(paste("total samples:", nrow(data_all)))
writeLines(paste("unknown samples (1-12):", nrow(unknown_samples)))
writeLines(paste("reference samples:", nrow(reference_samples)))
writeLines(paste("number of elements:", ncol(reference_numeric)))
writeLines(paste("pc1 variance explained:", round(variance_explained[1], 2), "%"))
writeLines(paste("pc2 variance explained:", round(variance_explained[2], 2), "%"))
writeLines(paste("cumulative pc1+pc2:", round(cumulative_variance[2], 2), "%"))
writeLines("")

writeLines("3. unknown sample identification results:")
for (i in 1:nrow(unknown_results)) {
  result <- unknown_results[i, ]
  
  writeLines(paste("unknown sample", result$Unknown_Sample, ":"))
  writeLines(paste("  closest match:", result$Closest_Match))
  writeLines(paste("  distance in pca space:", round(result$Distance, 4)))
  writeLines(paste("  correlation coefficient:", round(result$Correlation, 4)))
  writeLines("")
}

writeLines("5. complete distance ranking:")
writeLines("sample\tclosest_match\tdistance\tcorrelation\tquality")
for (i in 1:nrow(unknown_results)) {
  result <- unknown_results[i, ]
  writeLines(paste(result$Unknown_Sample, "\t", result$Closest_Match, "\t", 
      round(result$Distance, 4), "\t", round(result$Correlation, 4)))
}

writeLines("")
writeLines("6. principal components variance:")
for (i in 1:min(10, length(variance_explained))) {
  writeLines(paste("pc", i, ": ", round(variance_explained[i], 2), "% (cumulative: ", 
      round(cumulative_variance[i], 2), "%)"))
}

if (length(excluded_elements) > 0) {
  writeLines("")
  writeLines("7. excluded elements:")
  writeLines(paste("elements excluded due to zero/near-zero variance:", length(excluded_elements)))
  writeLines(paste(excluded_elements, collapse = ", "))
} else {
  writeLines("")
  writeLines("7. excluded elements:")
  writeLines("no elements excluded - all elements used in analysis")
}

sink()

# Analysis complete - print summary of generated files
message("\n=== analysis complete ===")
message("files generated in: ", output_dir)
message("1. pca_report.pdf - comprehensive visual analysis")
message("2. pca_summary_table.txt - detailed summary tables")
message("3. pca_results.txt - complete analysis results")

# Verify files were created
generated_files <- c("pca_report.pdf", "pca_summary_table.txt", "pca_results.txt")
for (file in generated_files) {
  file_path <- file.path(output_dir, file)
  if (file.exists(file_path)) {
    file_size <- file.info(file_path)$size
    message("✓ ", file, " - size: ", file_size, " bytes")
  } else {
    message("✗ ", file, " - not created")
  }
}