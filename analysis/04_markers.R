# Title: Cluster markers (FindAllMarkers)
# Project: AD snRNA-seq multi-model analysis
# Author: Jung Hyun Park
# Purpose: Compute cluster markers after clustering and save a CSV.
# Notes:
#   - Works on clustered object from analysis/03_clustering.R
#   - Uses SCT assay for marker finding (with PrepSCTFindMarkers when multiple SCT models exist).
#   - For downstream DE between conditions, use RNA assay later (separate script).

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(fs)
  library(readr)
})

# -----------------------
# Parameters (edit here)
# -----------------------
in_rds  <- "outputs/clustering/unannotated_clustered_res0.4.rds"
out_dir <- "outputs/markers"
dir_create(out_dir)

# Which identity to use for markers:
# - Default: active ident already set by 03_clustering.R to SCT_snn_res.<selected>
# - If you stored it in misc, we can set it explicitly.
use_selected_ident_from_misc <- TRUE

# FindAllMarkers parameters
only_pos <- TRUE
min_pct <- 0.25
logfc_threshold <- 0.25

# -----------------------
# Load clustered object
# -----------------------
if (!file.exists(in_rds)) stop("Missing input: ", in_rds, "\nRun analysis/03_clustering.R first.")
obj <- readRDS(in_rds)

DefaultAssay(obj) <- "SCT"

# Set identities consistently (optional)
if (use_selected_ident_from_misc && !is.null(obj@misc$cluster_resolution_selected)) {
  res_selected <- obj@misc$cluster_resolution_selected
  res_col <- paste0("SCT_snn_res.", res_selected)
  if (res_col %in% colnames(obj[[]])) {
    Idents(obj) <- res_col
  }
}

# PrepSCTFindMarkers is needed when object contains multiple SCT models (common after merging/integration)
obj <- PrepSCTFindMarkers(obj)

# -----------------------
# Find markers
# -----------------------
markers <- FindAllMarkers(
  object = obj,
  only.pos = only_pos,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold
)

# Save marker table
res_tag <- if (!is.null(obj@misc$cluster_resolution_selected)) paste0("res", obj@misc$cluster_resolution_selected) else "resNA"
out_csv <- file.path(out_dir, paste0("cluster_markers_", res_tag, ".csv"))
write.csv(markers, out_csv, row.names = FALSE)

# Save a copy of the object (optional; useful if you want markers embedded later)
saveRDS(obj, file.path(out_dir, paste0("clustered_for_markers_", res_tag, ".rds")))

message("Done.")
message("Markers CSV: ", out_csv)

