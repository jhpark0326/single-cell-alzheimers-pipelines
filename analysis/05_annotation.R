# Title: Annotation (SAHA marker-based + manual mapping)
# Project: AD snRNA-seq multi-model analysis
# Author: Jung Hyun Park
# Purpose: Generate automated marker-based annotation with SAHA and apply manual cluster labels in a reproducible way.
# Notes:
#   - SAHA requires built-in datasets: ISOCTX_Markers and ABC_meta. Refer to https://github.com/neurogenetics/SAHA

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(fs)
  library(ggplot2)
})

set.seed(4444)

# ---- Load clustered object ----
obj <- readRDS("outputs/clustering/unannotated_clustered_res0.4.rds")
DefaultAssay(obj) <- "SCT"

# Ensure identities are the selected clustering identities (numeric clusters)
# (03_clustering.R should already set this, but we make it explicit)
if (!is.null(obj@misc$cluster_resolution_selected)) {
  res_col <- paste0("SCT_snn_res.", obj@misc$cluster_resolution_selected)
  if (res_col %in% colnames(obj[[]])) Idents(obj) <- res_col
}

# Stable cluster id (numeric as character: "0","1",...)
obj$cluster_id <- as.character(Idents(obj))

# ---- Output folders ----
out_dir <- "outputs/annotation"
fig_dir <- file.path(out_dir, "qc_plots")
dir_create(out_dir)
dir_create(fig_dir)

# ============================================================
# Part A) SAHA marker-based annotation (automated suggestions)
# ============================================================
# SAHA helps you *suggest* labels; final labels are recorded in config/ mapping tables (Part B).
run_saha <- TRUE

if (run_saha) {
  suppressPackageStartupMessages(library(SAHA))

  markers_csv <- "outputs/markers/cluster_markers_res0.4.csv"
  if (!file.exists(markers_csv)) stop("Missing markers CSV: ", markers_csv, "\nRun analysis/04_markers.R first.")
  demo_Markers <- read.csv(markers_csv)

  # SAHA built-in references (Marker-based only; AvgExp method ignored)
  data("ISOCTX_Markers", package = "SAHA")
  data("ABC_meta",       package = "SAHA")
  meta <- ABC_meta

  ann <- Create_SAHA_object(query = demo_Markers, db = ISOCTX_Markers, data_type = "Markers")
  ann <- Initialize_Markers(ann)
  ann <- Initialize_Self_Similiarity(ann, slot = "Markers")
  ann <- Run_Marker_Based(ann)

  ann <- Create_SelfSimilarity_Viz(ann, slot = "Markers")
  ann <- Create_MarkerBased_Viz(ann, meta = meta, facet = TRUE)

  pdf(file.path(out_dir, "SAHA_marker_based.pdf"), width = 10, height = 7)
  call_SAHA_plots(ann, plot_type = "self-similarity", data_type = "Markers")
  call_SAHA_plots(ann, plot_type = "Marker-based",    data_type = "Markers")
  dev.off()

  saveRDS(ann, file.path(out_dir, "SAHA_object_markers.rds"))
}

# ============================================================
# Part B) Manual annotation (record your final decisions)
# ============================================================
# After reviewing SAHA outputs + marker QC, record final labels in:
#   config/cluster_to_cluster_label.csv   (fine-grained cluster labels, optional)
#   config/cluster_to_cell_type.csv       (broad cell types, recommended)
#
# These mapping tables are the reproducible record of your interpretation.

dir_create("config")

map_cluster <- "config/cluster_to_cluster_label.csv"
map_cell    <- "config/cluster_to_cell_type.csv"

# If mapping files are missing, auto-create starter files (Option A-friendly)
# so a GitHub viewer can see the expected format immediately.
make_templates_if_missing <- TRUE

if (make_templates_if_missing) {
  cl_ids <- sort(unique(as.character(obj$cluster_id)))

  if (!file.exists(map_cluster)) {
    write_csv(
      tibble(cluster_id = cl_ids, cluster_label = ""),
      map_cluster
    )
    message("Created: ", map_cluster, " (fill in cluster_label column after reviewing SAHA/markers)")
  }

  if (!file.exists(map_cell)) {
    write_csv(
      tibble(cluster_id = cl_ids, cell_type = ""),
      map_cell
    )
    message("Created: ", map_cell, " (fill in cell_type column after reviewing SAHA/markers)")
  }
}

# Load mappings (if present) and apply
apply_map <- function(obj, map_path, value_col) {
  if (!file.exists(map_path)) return(obj)

  m <- read_csv(map_path, show_col_types = FALSE)
  if (!all(c("cluster_id", value_col) %in% colnames(m))) {
    stop("Mapping file must have columns: cluster_id,", value_col, " -> ", map_path)
  }

  m$cluster_id <- as.character(m$cluster_id)
  lookup <- setNames(m[[value_col]], m$cluster_id)

  obj[[value_col]] <- unname(lookup[as.character(obj$cluster_id)])
  obj
}

obj <- apply_map(obj, map_cluster, "cluster_label")
obj <- apply_map(obj, map_cell,    "cell_type")

# ============================================================
# QC plots (public-safe)
# ============================================================
if ("umap.harmony" %in% Reductions(obj)) {
  pdf(file.path(fig_dir, "umap_by_cluster_id.pdf"), width = 9, height = 6)
  print(DimPlot(obj, reduction = "umap.harmony", group.by = "cluster_id",
                label = TRUE, repel = TRUE, raster = FALSE) + NoLegend())
  dev.off()

  if ("cluster_label" %in% colnames(obj[[]]) && any(nzchar(obj$cluster_label))) {
    pdf(file.path(fig_dir, "umap_by_cluster_label.pdf"), width = 10, height = 6)
    print(DimPlot(obj, reduction = "umap.harmony", group.by = "cluster_label",
                  label = TRUE, repel = TRUE, raster = FALSE) + NoLegend())
    dev.off()
  }

  if ("cell_type" %in% colnames(obj[[]]) && any(nzchar(obj$cell_type))) {
    pdf(file.path(fig_dir, "umap_by_cell_type.pdf"), width = 9, height = 6)
    print(DimPlot(obj, reduction = "umap.harmony", group.by = "cell_type",
                  label = TRUE, repel = TRUE, raster = FALSE) + NoLegend())
    dev.off()
  }
}

# Basic marker sanity DotPlot (edit as you like)
genes <- c("Slc17a7","Gad1","Gad2","Mbp","Mog","P2ry12","Aqp4","Aldoc","Vcan","Flt1")
genes <- genes[genes %in% rownames(obj)]

if (length(genes) > 0) {
  pdf(file.path(fig_dir, "dotplot_basic_markers_by_cluster.pdf"), width = 10, height = 5)
  print(DotPlot(obj, features = genes, group.by = "cluster_id") + RotatedAxis())
  dev.off()

  if ("cell_type" %in% colnames(obj[[]]) && any(nzchar(obj$cell_type))) {
    pdf(file.path(fig_dir, "dotplot_basic_markers_by_cell_type.pdf"), width = 10, height = 5)
    print(DotPlot(obj, features = genes, group.by = "cell_type") + RotatedAxis())
    dev.off()
  }
}

# Save annotated object
saveRDS(obj, file.path(out_dir, "annotated.rds"))

message("Done.")
message("Saved: ", file.path(out_dir, "annotated.rds"))
if (run_saha) message("SAHA PDF: ", file.path(out_dir, "SAHA_marker_based.pdf"))
message("QC plots: ", fig_dir)
