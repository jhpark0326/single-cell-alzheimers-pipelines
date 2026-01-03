# Title: Differential expression 
# Project: AD snRNA-seq multi-model analysis
# Author: Jung Hyun Park
# Purpose: Run DE by annotation label (cluster_label or cell_type) across conditions,
#          then make volcano plots for visualization.
# Notes:
#   - DE is performed on the RNA assay; clustering/integration assumed done upstream.
#   - Uses MAST with TECHNICAL covariates only (latent.vars). Do NOT include biological factors (e.g., Sex).

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(EnhancedVolcano)
  library(ComplexUpset)
  library(fs)
})

# -----------------------
# Parameters (edit here)
# -----------------------
in_rds <- "outputs/annotation/annotated.rds"

out_dir <- "outputs/de"
fig_dir <- file.path(out_dir, "figures")
dir_create(out_dir)
dir_create(fig_dir)

# Choose label used for DE stratification (matches your propeller choices)
# options: "cluster_label" or "cell_type"
LABEL_COL <- "cluster_label"

# Baseline + targets for contrasts 
BASELINE_CONDITION <- "Control"
TARGET_CONDITIONS  <- NULL   # NULL => auto-detect all non-baseline conditions

# Thresholds
ALPHA <- 0.05
FC_CUTOFF <- 1.5
ABS_LOG2FC_CUTOFF <- log2(FC_CUTOFF)   # ~0.585

# DE settings
ASSAY_FOR_DE <- "RNA"
DE_TEST <- "MAST"

# TECHNICAL covariates only (used if present)
LATENT_VARS_CANDIDATES <- c("nCount_RNA")

# -----------------------
# Load annotated object
# -----------------------
if (!file.exists(in_rds)) stop("Missing input: ", in_rds)
obj <- readRDS(in_rds)

# -----------------------
# Validate metadata
# -----------------------
meta_cols <- colnames(obj[[]])

req <- c("Condition", "sample_id", LABEL_COL)
missing <- setdiff(req, meta_cols)
if (length(missing) > 0) {
  stop("Object metadata missing required columns: ", paste(missing, collapse = ", "))
}

obj$Condition <- as.character(obj$Condition)
obj$sample_id <- as.character(obj$sample_id)
obj$label_use <- as.character(obj[[LABEL_COL]][, 1])

# Make a stable factor order for conditions
all_conditions <- sort(unique(obj$Condition))
if (!(BASELINE_CONDITION %in% all_conditions)) {
  stop("BASELINE_CONDITION not found in obj$Condition: ", BASELINE_CONDITION)
}

if (is.null(TARGET_CONDITIONS)) {
  TARGET_CONDITIONS <- setdiff(all_conditions, BASELINE_CONDITION)
}
if (length(TARGET_CONDITIONS) == 0) stop("No TARGET_CONDITIONS found (non-baseline).")

# Make contrasts
contrasts <- lapply(TARGET_CONDITIONS, function(x) c(x, BASELINE_CONDITION))

# Label levels (cell types / cluster labels)
labels <- sort(unique(obj$label_use))

# -----------------------
# Helper: run DE for one label and one contrast
# -----------------------
run_label_de <- function(seu, label, cond1, cond2) {
  DefaultAssay(seu) <- ASSAY_FOR_DE

  # Identity: label + condition (general, no genotype)
  seu$DEG_group <- paste0(seu$label_use, "__", seu$Condition)
  Idents(seu) <- "DEG_group"

  ident_1 <- paste0(label, "__", cond1)
  ident_2 <- paste0(label, "__", cond2)

  lvls <- levels(Idents(seu))
  if (!(ident_1 %in% lvls && ident_2 %in% lvls)) return(NULL)

  m <- FindMarkers(
    seu,
    ident.1 = ident_1,
    ident.2 = ident_2,
    test.use = DE_TEST,
    logfc.threshold = 0,
    only.pos = FALSE,
    latent.vars = LATENT_VARS_CANDIDATES
  )

  m %>%
    rownames_to_column("gene") %>%
    mutate(
      label = label,
      contrast = paste0(cond1, "_vs_", cond2),
      cond1 = cond1,
      cond2 = cond2
    )
}

# -----------------------
# Run DE: label × contrast
# -----------------------
res_list <- list()

for (lb in labels) {
  for (cc in contrasts) {
    cond1 <- cc[1]; cond2 <- cc[2]
    message("DE: ", LABEL_COL, "=", lb, "  ", cond1, " vs ", cond2)

    rr <- run_label_de(obj, lb, cond1, cond2)
    if (!is.null(rr) && nrow(rr) > 0) {
      key <- paste(LABEL_COL, lb, cond1, "vs", cond2, sep = "__")
      res_list[[key]] <- rr
    }
  }
}

deg <- bind_rows(res_list)
if (nrow(deg) == 0) stop("No DE results produced. Check labels/conditions present in metadata.")

# Standardize expected columns (Seurat MAST output uses avg_log2FC and p_val_adj)
stopifnot(all(c("avg_log2FC", "p_val_adj") %in% colnames(deg)))

# Handle p_adj == 0
min_nonzero <- min(deg$p_val_adj[deg$p_val_adj > 0], na.rm = TRUE)
deg$p_val_adj[deg$p_val_adj == 0] <- min_nonzero * 1e-3

# Add direction + significance columns (integrated up+down)
deg <- deg %>%
  mutate(
    sig = !is.na(p_val_adj) & p_val_adj <= ALPHA & !is.na(avg_log2FC) & abs(avg_log2FC) >= ABS_LOG2FC_CUTOFF,
    direction = case_when(
      sig & avg_log2FC >=  ABS_LOG2FC_CUTOFF ~ "Up",
      sig & avg_log2FC <= -ABS_LOG2FC_CUTOFF ~ "Down",
      TRUE ~ "NS"
    )
  )

write.csv(deg, file.path(out_dir, paste0("DE_", LABEL_COL, "_all_contrasts.csv")), row.names = FALSE)

# -----------------------
# Volcano plots (true volcano: plot ALL genes, label sig only)
# -----------------------
vol_dir <- file.path(fig_dir, "volcano")
dir_create(vol_dir)

make_volcano <- function(df, title, out_pdf) {
  sig_genes <- df %>% filter(sig) %>% pull(gene) %>% unique()

  # dynamic xlim
  xlim_use <- range(df$avg_log2FC, na.rm = TRUE)
  pad <- 0.25 * diff(xlim_use)
  if (is.finite(pad) && pad > 0) xlim_use <- c(xlim_use[1] - pad, xlim_use[2] + pad)

  pdf(out_pdf, width = 6.5, height = 6)
  print(
    EnhancedVolcano(
      df,
      lab = df$gene,
      selectLab = sig_genes,
      x = "avg_log2FC",
      y = "p_val_adj",
      title = title,
      pCutoff = ALPHA,
      FCcutoff = ABS_LOG2FC_CUTOFF,
      xlim = xlim_use,
      ylab = expression(-log[10]("adj.p")),
      legendPosition = "none",
      pointSize = 1.1,
      labSize = 3.8,
      max.overlaps = 20
    )
  )
  dev.off()
}

for (ct in unique(deg$contrast)) {
  for (lb in unique(deg$label)) {
    df <- deg %>% filter(contrast == ct, label == lb)
    if (nrow(df) == 0) next

    out_pdf <- file.path(vol_dir, paste0("volcano__", LABEL_COL, "__", lb, "__", ct, ".pdf"))
    make_volcano(df, paste0(lb, " — ", ct), out_pdf)
  }
}


message("Done.")
message("DE table: ", file.path(out_dir, paste0("DE_", LABEL_COL, "_all_contrasts.csv")))
message("Volcano PDFs: ", vol_dir)
