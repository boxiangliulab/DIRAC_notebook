# Load necessary libraries
library(SingleCellExperiment)
library(scmap)
library(zellkonverter)  # For reading h5ad files
library(dplyr)
library(Matrix)
library(ggplot2)
library(tibble)
library(yaml)

# Set parameters
data_path <- "/home/project/11003054/changxu/Projects/DIRAC/Section-4"
methods <- "SCMAP"
sim_methods <- "layers"
omics <- "RNA"

# Define color palette
colormaps <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#FFFF33", "#F781BF", "#999999",
               "#E5D8BD", "#B3CDE3", "#CCEBC5", "#FED9A6", "#FBB4AE", "#8DD3C7", "#BEBADA", "#80B1D3", "#B3DE69", "#FCCDE5",
               "#BC80BD", "#FFED6F", "#8DA0CB", "#E78AC3", "#E5C494", "#CCCCCC", "#FB9A99", "#E31A1C", "#CAB2D6", "#6A3D9A", 
               "#B15928")

# Create a list to save all results
results_list <- list()

# Outer loop: i = 0 to 7
for (i in 0:7) {
  cat(sprintf("Processing dataset i = %d...\n", i))
  
  # 1. Load the data
  file_path <- sprintf("%s/scMultiSim_data/sim-%s/sim_%s_%d.h5ad", data_path, sim_methods, omics, i)
  sce_all <- readH5AD(file_path)
  
  # 2. Preprocessing
  logcounts(sce_all) <- log2(counts(sce_all) + 1)
  rowData(sce_all)$feature_symbol <- rownames(sce_all)
  assay(sce_all, "counts") <- assay(sce_all, "logcounts")
  sce_all <- sce_all[!duplicated(rownames(sce_all)), ]
  
  # Inner loop: j = 1 to 5
  for (j in 1:5) {
    cat(sprintf("  Current batch = %d\n", j))
    
    # 3. Split into source and target
    source_sce <- sce_all[, sce_all$batch != j]
    target_sce <- sce_all[, sce_all$batch == j]
    
    # 4. Feature selection (only on the source data)
    source_sce <- selectFeatures(source_sce, suppress_plot = TRUE)
    
    # 5. Build index
    source_sce <- indexCluster(source_sce)
    
    # 6. Predict target
    scmapCluster_results <- scmapCluster(
      projection = target_sce,
      index_list = list(source = metadata(source_sce)$scmap_cluster_index)
    )
    
    # 7. Evaluate results
    pred_labels <- as.character(scmapCluster_results$combined_labs)
    true_labels <- as.character(colData(target_sce)$cell.type)
    
    # Exclude unassigned predictions
    valid_idx <- pred_labels != "unassigned"
    
    if (sum(valid_idx) > 0) {
      accuracy <- sum(pred_labels[valid_idx] == true_labels[valid_idx]) / length(valid_idx)
      precision <- caret::precision(data = factor(pred_labels[valid_idx], levels = unique(true_labels)),
                                    reference = factor(true_labels[valid_idx], levels = unique(true_labels)))
      recall <- caret::recall(data = factor(pred_labels[valid_idx], levels = unique(true_labels)),
                              reference = factor(true_labels[valid_idx], levels = unique(true_labels)))
      f1 <- caret::F_meas(data = factor(pred_labels[valid_idx], levels = unique(true_labels)),
                          reference = factor(true_labels[valid_idx], levels = unique(true_labels)))
    } else {
      accuracy <- NA
      precision <- NA
      recall <- NA
      f1 <- NA
    }
    
    metrics_all <- list(
      "Accuracy Score" = accuracy,
      "Precision Score" = precision,
      "Recall Score" = recall,
      "F1 Score" = f1
    )
    
    print(metrics_all)
    
    # 8. Save results
    save_path <- sprintf("%s/Results/merged_%s_%s_%s_%d_%d", data_path, omics, methods, sim_methods, i, j)
    if (!dir.exists(save_path)) {
      dir.create(save_path, recursive = TRUE)
    }
    
    # Save evaluation metrics
    write_yaml(metrics_all, file.path(save_path, "annotate_metrics.yaml"))
    
    # Save predictions
    output_df <- data.frame(
      Ground_Truth = true_labels,
      Predicted = pred_labels,
      stringsAsFactors = FALSE
    )
    write.csv(output_df, file.path(save_path, sprintf("%s_%s_%s.csv", omics, methods, sim_methods)), row.names = FALSE)
    
    # Plot results (optional)
    plot_df <- output_df %>% count(Ground_Truth, Predicted)
    
    p <- ggplot(plot_df, aes(x = Ground_Truth, y = n, fill = Predicted)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      scale_fill_manual(values = colormaps) +
      labs(title = sprintf("i=%d, j=%d: scmap Results", i, j),
           x = "True Labels", y = "Count", fill = "Predicted Labels")
    
    ggsave(file.path(save_path, sprintf("%s_%s_%s.pdf", omics, methods, sim_methods)), plot = p, width = 8, height = 6)
    
    # Collect results
    results_list[[paste0("i", i, "_j", j)]] <- metrics_all
  }
}

# 9. Summarize all results
summary_df <- bind_rows(results_list, .id = "dataset")
write.csv(summary_df, file.path(data_path, "Results", sprintf("sim_metrics_%s_%s_%s.csv", omics, methods, sim_methods)), row.names = FALSE)

cat("All processing completed!\n")
