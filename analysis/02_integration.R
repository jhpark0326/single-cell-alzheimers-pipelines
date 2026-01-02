# Title: Data integration (snRNA-seq) - Harmony on the SCT assay
# Project: AD snRNA-seq multi-model analysis
# Author: Jung Hyun Park
# Purpose: Integrate multiple snRNA-seq samples after preprocessing (SoupX/DF handled upstream).
# Notes:
#   - Integration/clustering are performed on the SCT assay (variance-stabilized).
#   - Differential expression is typically performed later on the RNA assay.
#   - Harmony corrects TECHNICAL covariates only (e.g., ID, Batch). Do NOT include Condition.
#   - Expects per-sample objects saved by analysis/01_qc_preprocessing.R:
#       outputs/preprocess/seurat_preprocessed_<sample_id>.rds
#   - Writes integrated objects + QC plots to outputs/integration/

suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(fs)
  library(ggplot2)
  library(harmony)
})

# -----------------------
# Parameters (edit here)
# -----------------------
meta_path <- "data/metadata.csv"      # local file (do not commit)
pre_dir   <- "outputs/preprocess"     # from 01_qc_preprocessing.R

out_dir <- "outputs/integration"
fig_dir <- file.path(out_dir, "qc_plots")
dir_create(out_dir)
dir_create(fig_dir)

# Harmony variables (TECHNICAL only).
# Default: ID (sample/library effect).
# Optional: include Batch too if you truly have multi-batch sequencing.
# Tip: if Batch is confounded with Condition, prefer ID-only to avoid removing biology.
harmony_vars <- c("ID")  # change to c("Batch","ID") if needed

# PCA / clustering defaults
npcs  <- 50
cluster_resolution <- 0.2
choose_ndims <- function(seu, target = 90, min_dims = 10, max_dims = 50) {
  pvar <- seu[["pca"]]@stdev^2
  cumv <- cumsum(pvar / sum(pvar) * 100)
  nd <- min(which(cumv >= target))
  max(min_dims, min(nd, max_dims))
}

# Optional “hard doublet” rule (flags only; does NOT remove cells)
use_hard_doublet <- TRUE
hard_doublet_rule <- function(seu) {
  # Example conservative thresholds; adjust per chemistry/prep
  seu$HardDoublet <- "Singlet"
  seu$HardDoublet[seu$nFeature_RNA >= 8000 | seu$nCount_RNA >= 40000] <- "Doublet"

  # Combine DF (from preprocessing) + HardDoublet into a single flag
  seu$Doublet <- "Singlet"
  if ("DF" %in% colnames(seu[[]])) {
    seu$Doublet[seu$DF == "Doublet" | seu$HardDoublet == "Doublet"] <- "Doublet"
  } else {
    seu$Doublet[seu$HardDoublet == "Doublet"] <- "Doublet"
  }
  seu
}

# -----------------------
# Load metadata
# -----------------------
if (!file.exists(meta_path)) {
  stop(
    "Missing metadata file: ", meta_path, "\n",
    "Tip: copy data/metadata_template.csv to data/metadata.csv locally."
  )
}
meta <- read_csv(meta_path, show_col_types = FALSE)

required_cols <- c("sample_id", "id", "condition")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("metadata.csv missing required columns: ", paste(missing_cols, collapse = ", "))
}

# -----------------------
# Load per-sample objects
# -----------------------
obj_list <- vector("list", nrow(meta))
names(obj_list) <- meta$sample_id

for (i in seq_len(nrow(meta))) {
  sample_id <- as.character(meta$sample_id[i])
  rds_path  <- file.path(pre_dir, paste0("seurat_preprocessed_", sample_id, ".rds"))
  if (!file.exists(rds_path)) stop("Missing preprocessed object: ", rds_path)

  seu <- readRDS(rds_path)

  # Integration/clustering on SCT assay (variance-stabilized);
  # DE later is typically performed on RNA assay.
  DefaultAssay(seu) <- "SCT"

  # Ensure public-safe labels exist (from 01 script or fallback from metadata)
  if (!"Condition" %in% colnames(seu[[]])) seu$Condition <- as.character(meta$condition[i])
  if (!"ID" %in% colnames(seu[[]]))        seu$ID        <- as.character(meta$id[i])
  if ("batch" %in% colnames(meta) && !"Batch" %in% colnames(seu[[]])) {
    seu$Batch <- as.character(meta$batch[i])
  }

  # Optional doublet flagging (does not subset)
  if (use_hard_doublet) seu <- hard_doublet_rule(seu)

  obj_list[[sample_id]] <- seu
}

message("Loaded ", length(obj_list), " preprocessed objects.")

# -----------------------
# Merge + PCA (SCT)
# -----------------------
integ_features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)

obj <- merge(
  x = obj_list[[1]],
  y = obj_list[2:length(obj_list)],
  merge.data = TRUE
)

DefaultAssay(obj) <- "SCT"
VariableFeatures(obj) <- integ_features

obj <- RunPCA(obj, assay = "SCT", npcs = npcs, verbose = FALSE)
saveRDS(obj, file.path(out_dir, "merged_preintegration.rds"))

ndims <- choose_ndims(obj, target = 90, min_dims = 10, max_dims = 50)
# elbow plot for manual sanity check
pdf(file.path(fig_dir, "pca_elbow_preintegration.pdf"), width = 7, height = 5)
print(ElbowPlot(obj, ndims = npcs))
dev.off()
obj@misc$ndims_harmony <- ndims
# -----------------------
# Harmony integration (SCT PCs)
# -----------------------
missing_hvars <- setdiff(harmony_vars, colnames(obj[[]]))
if (length(missing_hvars) > 0) {
  stop(
    "Harmony vars missing from object metadata: ", paste(missing_hvars, collapse = ", "), "\n",
    "Add them to metadata.csv (e.g., batch column) or change harmony_vars."
  )
}

obj <- RunHarmony(
  object = obj,
  group.by.vars = harmony_vars,
  reduction = "pca",
  assay.use = "SCT",
  reduction.save = "harmony"
)

# Cluster/UMAP on Harmony embedding
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:ndims, verbose = FALSE)
obj <- FindClusters(obj, resolution = cluster_resolution, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:ndims, verbose = FALSE, reduction.name="umap.harmony")

saveRDS(obj, file.path(out_dir, "integrated_harmony.rds"))

# -----------------------
# QC plots (public-safe)
# -----------------------
pdf(file.path(fig_dir, "umap_harmony_by_condition.pdf"), width = 8, height = 6)
print(DimPlot(obj, reduction = "umap", group.by = "Condition", raster = FALSE, label = FALSE))
dev.off()

if ("Batch" %in% colnames(obj[[]])) {
  pdf(file.path(fig_dir, "umap_harmony_by_batch.pdf"), width = 8, height = 6)
  print(DimPlot(obj, reduction = "umap", group.by = "Batch", raster = FALSE, label = FALSE))
  dev.off()
}

if ("ID" %in% colnames(obj[[]])) {
  pdf(file.path(fig_dir, "umap_harmony_by_id.pdf"), width = 8, height = 6)
  print(DimPlot(obj, reduction = "umap", group.by = "ID", raster = FALSE, label = FALSE))
  dev.off()
}

if ("Doublet" %in% colnames(obj[[]])) {
  pdf(file.path(fig_dir, "umap_harmony_by_doubletflag.pdf"), width = 8, height = 6)
  print(DimPlot(obj, reduction = "umap", group.by = "Doublet", raster = FALSE, label = FALSE))
  dev.off()
}

message("Done.")
message("Merged object saved to: ", file.path(out_dir, "merged_preintegration.rds"))
message("Integrated object saved to: ", file.path(out_dir, "integrated_harmony.rds"))
message("QC plots saved to: ", fig_dir)
