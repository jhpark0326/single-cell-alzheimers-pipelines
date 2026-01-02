# Title: QC + preprocessing (SoupX + DoubletFinder)
# Project: AD snRNA-seq multi-model analysis
# Author: Jung Hyun Park
# Purpose: Generalized workflow derived from research analyses (no raw data included).
# Outputs: per-sample Seurat objects + QC plots + preprocessing summary (all in outputs/; ignored by git)

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(DoubletFinder)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(fs)
})

# -----------------------
# Parameters (edit here)
# -----------------------
input_dir  <- "data/10x"          # data/10x/<sample_id>/  OR data/10x/<sample_id>/outs
meta_path  <- "data/metadata.csv" # local file (do not commit). copy from data/metadata_template.csv
is_cellranger_run <- FALSE        # TRUE if sample folder contains an "outs/" directory

out_dir <- "outputs/preprocess"
fig_dir <- file.path(out_dir, "qc_plots")
dir_create(out_dir)
dir_create(fig_dir)

# QC thresholds
min_features <- 200
max_percent_mt <- 1

# PCA dims selection (cumulative variance explained)
target_cumvar <- 90
ndims_min <- 10
ndims_max <- 50

cluster_resolution <- 0.2

# Turn modules on/off
do_soupx <- TRUE
do_doubletfinder <- TRUE

# Expected doublet rate (edit if your lab uses a different heuristic)
expected_doublet_rate <- function(n_cells) (n_cells / 1000) * 0.008

# -----------------------
# Tiny helpers (keep minimal)
# -----------------------
choose_ndims <- function(seu, target = 90, min_dims = 10, max_dims = 50) {
  pvar <- seu[["pca"]]@stdev^2
  cumv <- cumsum(pvar / sum(pvar) * 100)
  nd <- min(which(cumv >= target))
  max(min_dims, min(nd, max_dims))
}

standardize_df_cols <- function(seu) {
  # remove old pANN_* columns (prevents accumulation across reruns)
  seu@meta.data <- seu@meta.data[, !grepl("^pANN_", colnames(seu@meta.data)), drop = FALSE]

  pANN_col <- grep("^pANN_", colnames(seu[[]]), value = TRUE)
  DF_col_1 <- grep("^DF\\.classifications_", colnames(seu[[]]), value = TRUE)
  DF_col_2 <- grep("^DF\\.", colnames(seu[[]]), value = TRUE)

  if (length(pANN_col) == 1) {
    seu$pANN <- seu[[pANN_col]][, 1]
    seu[[pANN_col]] <- NULL
  }
  if (length(DF_col_1) == 1) {
    seu$DF <- seu[[DF_col_1]][, 1]
    seu[[DF_col_1]] <- NULL
  } else if (length(DF_col_2) == 1) {
    seu$DF <- seu[[DF_col_2]][, 1]
    seu[[DF_col_2]] <- NULL
  }
  seu
}

# -----------------------
# Load metadata (public-safe)
# -----------------------
if (!file_exists(meta_path)) stop("Missing metadata file: ", meta_path, "\nTip: copy metadata_template.csv to metadata.csv locally.")
meta <- read_csv(meta_path, show_col_types = FALSE)

required_cols <- c("sample_id", "id", "condition", "sex", "age_months")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) stop("metadata.csv missing required columns: ", paste(missing_cols, collapse = ", "))

# optional: track SoupX rho estimates
soupx_summary_path <- file.path(out_dir, "soupx_summary.csv")
if (do_soupx && file_exists(soupx_summary_path)) file_delete(soupx_summary_path)

# store per-sample summary rows
summary_rows <- list()

# -----------------------
# Main loop
# -----------------------
for (i in seq_len(nrow(meta))) {

  sample_id <- as.character(meta$sample_id[i])
  message("\n=== Processing: ", sample_id, " ===")

  sample_base <- file.path(input_dir, sample_id)
  sample_path <- if (is_cellranger_run) file.path(sample_base, "outs") else sample_base
  if (!dir_exists(sample_path)) stop("Sample directory not found: ", sample_path)

  # -----------------------
  # 1) Load + SoupX (optional)
  # -----------------------
  if (do_soupx) {
    soup <- load10X(sample_path)

    # Log contamination estimate (print output to a file)
    log_file <- file.path(out_dir, paste0("soupx_", sample_id, "_log.txt"))
    sink(log_file)
    soup <- autoEstCont(soup)
    print(soup)
    sink()

    # Save a minimal per-sample SoupX summary (QC only)
    nCells <- ncol(soup$toc)
    rho <- if (!is.null(soup$metaData$rho)) {
      if (length(soup$metaData$rho) == 1) soup$metaData$rho else median(soup$metaData$rho, na.rm = TRUE)
    } else if (!is.null(soup$rho)) {
      as.numeric(soup$rho)
    } else NA_real_

    write.table(
      data.frame(sample_id = sample_id, nCells = nCells, rho = rho),
      file = soupx_summary_path, sep = ",",
      row.names = FALSE, col.names = !file_exists(soupx_summary_path), append = TRUE
    )

    counts <- adjustCounts(soup)
  } else {
    counts <- Read10X(sample_path)
  }

  # -----------------------
  # 2) Seurat object + metadata
  # -----------------------
  seu <- CreateSeuratObject(
    counts = counts,
    project = "ad_snrnaseq_multimodel",
    min.cells = 3,
    min.features = min_features
  )

  seu$sample_id  <- sample_id
  seu$ID         <- as.character(meta$id[i])
  seu$Condition  <- as.character(meta$condition[i])
  seu$Sex        <- as.character(meta$sex[i])
  seu$Age_months <- as.character(meta$age_months[i])
  if ("batch" %in% colnames(meta)) seu$Batch <- as.character(meta$batch[i])

  # -----------------------
  # 3) QC metrics + pre-filter plot
  # -----------------------
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

  pdf(file.path(fig_dir, paste0("qc_violin_prefilter_", sample_id, ".pdf")), width = 10, height = 6)
  print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()

  # -----------------------
  # 4) Filter
  # -----------------------
  seu <- subset(seu, subset = nFeature_RNA > min_features & percent.mt < max_percent_mt)

  # -----------------------
  # 5) SCT + PCA + adaptive ndims
  # -----------------------
  seu <- SCTransform(seu, verbose = FALSE)
  seu <- RunPCA(seu, verbose = FALSE)
  ndims <- choose_ndims(seu, target = target_cumvar, min_dims = ndims_min, max_dims = ndims_max)

  # -----------------------
  # 6) Quick clustering/UMAP (broad, for DF homotypic adjustment)
  # -----------------------
  seu <- FindNeighbors(seu, dims = 1:ndims, verbose = FALSE)
  seu <- FindClusters(seu, resolution = cluster_resolution, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:ndims, verbose = FALSE)

  pdf(file.path(fig_dir, paste0("umap_quick_", sample_id, ".pdf")), width = 8, height = 6)
  print(DimPlot(seu, reduction = "umap", group.by = "seurat_clusters"))
  dev.off()

  # -----------------------
  # 7) DoubletFinder (optional)
  # -----------------------
  df_rate <- NA_real_
  if (do_doubletfinder) {
    # pK sweep
    sweep.res.list <- paramSweep(seu, PCs = 1:ndims, sct = TRUE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    write.csv(bcmvn, file.path(out_dir, paste0("doubletfinder_pK_", sample_id, ".csv")), row.names = FALSE)

    pK <- as.numeric(as.character(bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric, na.rm = TRUE), "pK"]))

    # expected doublets + homotypic adjustment
    n_cells <- ncol(seu)
    rate <- expected_doublet_rate(n_cells)
    nExp_poi <- round(rate * n_cells)
    homotypic.prop <- modelHomotypic(seu@meta.data$seurat_clusters)
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

    # run DF
    seu <- doubletFinder(seu, PCs = 1:ndims, pN = 0.25, pK = pK, nExp = nExp_poi.adj, sct = TRUE)
    seu <- standardize_df_cols(seu)

    # QC plots by DF
    pdf(file.path(fig_dir, paste0("doublet_violin_", sample_id, ".pdf")), width = 10, height = 6)
    print(VlnPlot(seu, features = c("nFeature_SCT", "nCount_SCT", "percent.mt"), group.by = "DF"))
    dev.off()

    pdf(file.path(fig_dir, paste0("doublet_umap_", sample_id, ".pdf")), width = 8, height = 6)
    print(DimPlot(seu, reduction = "umap", group.by = "DF"))
    dev.off()

    if ("DF" %in% colnames(seu[[]])) df_rate <- mean(seu$DF == "Doublet")
  }

  # -----------------------
  # 8) Save per-sample object
  # -----------------------
  saveRDS(seu, file.path(out_dir, paste0("seurat_preprocessed_", sample_id, ".rds")))

  # -----------------------
  # 9) Collect summary
  # -----------------------
  summary_rows[[sample_id]] <- tibble(
    sample_id = sample_id,
    ndims = ndims,
    n_cells = ncol(seu),
    nFeature_median = median(seu$nFeature_RNA),
    percent_mt_median = median(seu$percent.mt),
    doubletfinder = do_doubletfinder,
    df_doublet_rate = df_rate
  )
}

# write summary table
summary_df <- bind_rows(summary_rows)
write_csv(summary_df, file.path(out_dir, "preprocess_summary.csv"))

message("\nDone.")
message("Objects: ", out_dir)
message("QC plots: ", fig_dir)
if (do_soupx) message("SoupX summary: ", soupx_summary_path)
