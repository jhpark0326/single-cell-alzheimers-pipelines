# Title: pySCENIC input export from R Seurat
# Project: AD snRNA-seq multi-model analysis
# Author: Jung Hyun Park
# Purpose: Export raw RNA counts + metadata (and optional embeddings) from a Seurat object
#          for downstream pySCENIC (Python/CLI).
# Notes:
#   - Exports genes x cells as a sparse MatrixMarket file by default (recommended for size),
#     and can optionally export a CSV (not recommended for large datasets).
#   - Keeps metadata column names as-is. Ensure rownames(meta) match cell barcodes.
#   - Embeddings are optional (e.g., Harmony UMAP coordinates).

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(fs)
})

# -----------------------
# Parameters (edit here)
# -----------------------
in_rds  <- "outputs/annotation/annotated.rds"   # your Seurat object
out_dir <- "outputs/pyscenic"             # where to write files
dir_create(out_dir)

ASSAY_USE <- "RNA"
COUNT_SLOT <- "counts"

# Output toggles
WRITE_COUNTS_MTX <- TRUE    # recommended
WRITE_COUNTS_CSV <- FALSE   # huge for large datasets; use only if needed
WRITE_EMBEDDINGS <- TRUE

# Which embedding to export (if present)
EMBED_REDUCTION <- "umap.harmony"  # e.g. "umap", "umap.harmony", "pca", etc.
EMBED_FILE_BASENAME <- "umap_coord.csv"

# Output names
mtx_dir <- file.path(out_dir, "raw_counts_mtx")
csv_counts_path <- file.path(out_dir, "raw_counts_mat.csv")
meta_path <- file.path(out_dir, "meta.csv")
embed_path <- file.path(out_dir, EMBED_FILE_BASENAME)

# -----------------------
# Load Seurat object
# -----------------------
if (!file.exists(in_rds)) stop("Missing input: ", in_rds)
obj <- readRDS(in_rds)

DefaultAssay(obj) <- ASSAY_USE
if (!ASSAY_USE %in% Assays(obj)) stop("Assay not found: ", ASSAY_USE)

# -----------------------
# Export counts
# -----------------------
counts <- GetAssayData(obj, assay = ASSAY_USE, slot = COUNT_SLOT)

# Ensure sparse dgCMatrix (for MTX writing)
if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")

if (WRITE_COUNTS_MTX) {
  dir_create(mtx_dir)

  # MatrixMarket expects a matrix; we keep genes x cells orientation
  Matrix::writeMM(counts, file = file.path(mtx_dir, "matrix.mtx"))
  writeLines(rownames(counts), con = file.path(mtx_dir, "genes.tsv"))
  writeLines(colnames(counts), con = file.path(mtx_dir, "barcodes.tsv"))

  message("Wrote MTX counts to: ", mtx_dir)
}

if (WRITE_COUNTS_CSV) {
  # WARNING: can be extremely large; will coerce to dense
  message("Writing counts CSV (may be huge)...")
  write.csv(as.matrix(counts), csv_counts_path)
  message("Wrote CSV counts to: ", csv_counts_path)
}

# -----------------------
# Export metadata
# -----------------------
meta <- obj@meta.data

# Ensure cell barcode alignment
if (!all(rownames(meta) == colnames(counts))) {
  # Try to reorder meta to match counts
  meta <- meta[colnames(counts), , drop = FALSE]
  if (!all(rownames(meta) == colnames(counts))) {
    stop("Metadata rownames do not match count matrix colnames (cell barcodes).")
  }
}

write.csv(meta, meta_path, row.names = TRUE)
message("Wrote metadata to: ", meta_path)

# -----------------------
# Export embeddings (optional)
# -----------------------
if (WRITE_EMBEDDINGS) {
  if (EMBED_REDUCTION %in% names(obj@reductions)) {
    emb <- obj@reductions[[EMBED_REDUCTION]]@cell.embeddings
    emb <- emb[colnames(counts), , drop = FALSE]  # match cell order
    write.csv(emb, embed_path, row.names = TRUE)
    message("Wrote embeddings (", EMBED_REDUCTION, ") to: ", embed_path)
  } else {
    message("Embedding reduction not found: ", EMBED_REDUCTION, " (skipping).")
  }
}

message("Done.")
