# Title: Clustering + clustree (Harmony embedding)
# Project: AD snRNA-seq multi-model analysis
# Author: Jung Hyun Park
# Purpose: Choose clustering resolution (clustree) and save an unannotated clustered object.
# Notes:
#   - Assumes integration was done with Harmony on SCT PCs (analysis/02_integration.R).
#   - Clustering is performed on the Harmony embedding; SCT remains the default assay.
#   - Uses ndims stored in obj@misc$ndims_harmony when available for consistency.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(fs)
  library(ggplot2)
  library(clustree)
})

# -----------------------
# Parameters (edit here)
# -----------------------
in_rds  <- "outputs/integration/integrated_harmony.rds"
out_dir <- "outputs/clustering"
fig_dir <- file.path(out_dir, "qc_plots")
dir_create(out_dir)
dir_create(fig_dir)

set.seed(4444)

# Resolution sweep for clustree
res_grid <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)

# Fallback ndims if obj@misc$ndims_harmony is not present
ndims_fallback <- 30

# -----------------------
# Load integrated object
# -----------------------
if (!file.exists(in_rds)) stop("Missing input: ", in_rds, "\nRun analysis/02_integration.R first.")
obj <- readRDS(in_rds)

# Use SCT for clustering geometry; RNA is typically used later for DE
DefaultAssay(obj) <- "SCT"

# -----------------------
# Pick ndims (consistent with integration if stored)
# -----------------------
ndims <- ndims_fallback
if (!is.null(obj@misc$ndims_harmony) && is.numeric(obj@misc$ndims_harmony) && length(obj@misc$ndims_harmony) == 1) {
  ndims <- as.integer(obj@misc$ndims_harmony)
}
message("Using ndims = ", ndims, " (obj@misc$ndims_harmony if available; fallback=", ndims_fallback, ")")

# -----------------------
# Quick QC: elbow plot (optional)
# -----------------------
pdf(file.path(fig_dir, "elbow_pca.pdf"), width = 7, height = 5)
print(ElbowPlot(obj, ndims = 50))
dev.off()

# -----------------------
# Neighbors/UMAP on Harmony embedding
# -----------------------
# Assumes Harmony reduction is named "harmony" (from 02_integration.R)
if (!"harmony" %in% Reductions(obj)) {
  stop("Harmony reduction not found on object. Expected reduction name: 'harmony'.\n",
       "Check your integration script output.")
}

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:ndims, verbose = FALSE)

# Always (re)run UMAP here for reproducible clustering QC plots
obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = 1:ndims,
  reduction.name = "umap.harmony",
  verbose = FALSE
)

# QC UMAP plots (public-safe)
if ("Condition" %in% colnames(obj[[]])) {
  pdf(file.path(fig_dir, "umap_harmony_by_condition.pdf"), width = 8, height = 6)
  print(DimPlot(obj, reduction = "umap.harmony", group.by = "Condition", raster = FALSE))
  dev.off()
}
if ("ID" %in% colnames(obj[[]])) {
  pdf(file.path(fig_dir, "umap_harmony_by_id.pdf"), width = 8, height = 6)
  print(DimPlot(obj, reduction = "umap.harmony", group.by = "ID", raster = FALSE))
  dev.off()
}
if ("Batch" %in% colnames(obj[[]])) {
  pdf(file.path(fig_dir, "umap_harmony_by_batch.pdf"), width = 8, height = 6)
  print(DimPlot(obj, reduction = "umap.harmony", group.by = "Batch", raster = FALSE))
  dev.off()
}

# -----------------------
# Resolution sweep + clustree
# -----------------------
obj <- FindClusters(obj, resolution = res_grid, verbose = FALSE)

# Save object before pruning resolution columns (reproducibility)
saveRDS(obj, file.path(out_dir, "cluster_sweep_before_clustree.rds"))

# clustree expects cluster columns like SCT_snn_res.0.2, SCT_snn_res.0.4, ...
p_tree <- clustree(obj, prefix = "SCT_snn_res.")
ggsave(filename = file.path(fig_dir, "clustree.pdf"), plot = p_tree, width = 6, height = 8)

# -----------------------
# Choose a resolution and set identities
# -----------------------
res_selected <- 0.4
res_col <- paste0("SCT_snn_res.", res_selected)
if (!res_col %in% colnames(obj[[]])) stop("Selected resolution column not found: ", res_col)

Idents(obj) <- res_col

pdf(file.path(fig_dir, paste0("umap_harmony_res", res_selected, ".pdf")), width = 8, height = 6)
print(DimPlot(obj, reduction = "umap.harmony", label = TRUE, raster = FALSE))
dev.off()

# -----------------------
# Cluster representation by sample (sanity check)
# -----------------------
# Useful to catch clusters dominated by one sample/library
if ("ID" %in% colnames(obj[[]])) {
  cluster_sample_table <- table(Idents(obj), obj$ID)
  write.csv(cluster_sample_table, file.path(out_dir, paste0("cluster_by_sample_res", res_selected, ".csv")))
}

# -----------------------
# Clean up: keep only the selected resolution column
# -----------------------
res_cols_all  <- grep("^SCT_snn_res\\.", colnames(obj[[]]), value = TRUE)
res_cols_drop <- setdiff(res_cols_all, res_col)
for (cc in res_cols_drop) obj[[cc]] <- NULL

# record chosen resolution + dims for downstream scripts
obj@misc$cluster_resolution_selected <- res_selected
obj@misc$ndims_clustering <- ndims

# Save clustered (unannotated) object
out_rds <- file.path(out_dir, paste0("unannotated_clustered_res", res_selected, ".rds"))
saveRDS(obj, out_rds)

message("Done.")
message("Saved clustered object: ", out_rds)
message("Clustree + QC plots saved to: ", fig_dir)
