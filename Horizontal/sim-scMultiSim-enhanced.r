library(scMultiSim)
library(dplyr)
library(reticulate)  # For Python integration to save .h5ad

# Load Python's anndata package (must be installed in your Python env)
anndata <- import("anndata")

# Load GRN parameters
data(GRN_params_100, envir = environment())

# Ligand parameters
lig_params <- data.frame(
  target    = c(101, 102),
  regulator = c(103, 104),
  effect    = c(5.2, 5.9)
)

# Spatial options
spatial_options <- function(cell_count, layout_type) {
  cci_opt <- list(
    params = lig_params,
    max.neighbors = 6,
    start.layer = cell_count,
    grid.size = 128,
    cell.type.interaction = "random"
  )
  list(
    rand.seed = 0,
    GRN = GRN_params_100,
    num.cells = cell_count,
    num.cifs = 60,
    num.genes = 200,
    cif.sigma = 0.6,
    tree = Phyla5(),
    diff.cif.fraction = 0.9,
    speed.up = TRUE,
    cci = c(cci_opt, list(layout = layout_type))
  )
}

# Save function
save_results <- function(results, dir_path, sim_name) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }

  write.csv(results$counts, file = file.path(dir_path, "rna_counts.csv"), row.names = TRUE)
  write.csv(results$atacseq_data, file = file.path(dir_path, "atac_counts.csv"), row.names = TRUE)
  write.csv(results$cci_locs, file = file.path(dir_path, "cci_locs.csv"), row.names = TRUE)
  write.csv(results$cell_meta, file = file.path(dir_path, "cell_meta.csv"), row.names = TRUE)

  pdf(file.path(dir_path, "tsne_plot.pdf"), width = 8, height = 6)
  plot_tsne(log2(results$counts + 1),
           results$cell_meta$pop,
           legend = 'pop', plot.name = paste('True RNA Counts Tsne -', sim_name))
  dev.off()

  pdf(file.path(dir_path, "correlation_heatmap.pdf"), width = 8, height = 6)
  plot_gene_module_cor_heatmap(results)
  dev.off()

  pdf(file.path(dir_path, "spatial_distribution.pdf"), width = 8, height = 6)
  plot_cell_loc(results, show.arrows = FALSE)
  dev.off()

  spatial_coords <- as.matrix(results$cci_locs)
  rownames(spatial_coords) <- rownames(results$cell_meta)

  adata_rna <- anndata$AnnData(
    X = t(results$counts),
    obs = results$cell_meta,
    obsm = list(spatial = spatial_coords)
  )
  adata_rna$write(file.path(dir_path, "sim_RNA.h5ad"))

  if (!is.null(results$atac_counts)) {
    adata_atac <- anndata$AnnData(
      X = t(results$atac_counts),
      obs = results$cell_meta,
      obsm = list(spatial = spatial_coords)
    )
    adata_atac$write(file.path(dir_path, "sim_ATAC.h5ad"))
  }
}

# ==== Only generate one dataset ====

# Set parameters
cell_count <- 15000
layout_type <- "enhanced"
sim_name <- "sim-enhanced"
dir_path <- file.path("/home/project/11003054/changxu/Projects/DIRAC/Section-4/scMultiSim_data", sim_name)

cat("Generating simulation:", sim_name, "\n")

results <- sim_true_counts(spatial_options(cell_count, layout_type))
add_expr_noise(results, alpha_mean = 1e3)
divide_batches(results, nbatch = 5, effect = 3)
save_results(results, dir_path, sim_name)

cat("Simulation completed and saved to:", dir_path, "\n")
