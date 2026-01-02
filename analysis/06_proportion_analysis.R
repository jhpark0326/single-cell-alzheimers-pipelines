# Title: Proportion analysis (cluster + cell type)
# Project: <YOUR_PROJECT_NAME>
# Author: <YOUR_NAME>
#
# Purpose:
#   Quantify composition changes across experimental conditions at:
#     (1) cluster level (e.g., seurat_clusters)
#     (2) cell-type level (e.g., cell_type)
#   using propeller (speckle) with limma modeling.
#
# Notes:
#   - No setwd(); run from repo root.
#   - Uses obj$ID as sample identifier and obj$Condition as the group variable.
#   - Does not reference genotype-specific labels; uses Condition levels found in object.
#   - Requires metadata columns: ID, Condition.
#   - Output paths are relative (created if missing).
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(speckle)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
})

# -----------------------------
# Inputs / Outputs
# -----------------------------
OBJ_PATH <- "outputs/annotation/annotated.rds"

OUT_DIR <- "outputs/proportion"
PLOT_DIR <- file.path(OUT_DIR, "plots")
TAB_DIR  <- file.path(OUT_DIR, "tables")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR,  recursive = TRUE, showWarnings = FALSE)

# Toggle: if you have 3+ conditions, also compute all pairwise contrasts
RUN_PAIRWISE_IF_MULTI <- TRUE

# -----------------------------
# Load object
# -----------------------------
obj <- readRDS(OBJ_PATH)

# Basic checks
stopifnot("ID" %in% colnames(obj@meta.data))
stopifnot("Condition" %in% colnames(obj@meta.data))

obj$ID <- as.character(obj$ID)
obj$Condition <- as.character(obj$Condition)

# Drop cells missing essentials
keep <- !is.na(obj$ID) & obj$ID != "" & !is.na(obj$Condition) & obj$Condition != ""
obj <- subset(obj, cells = colnames(obj)[keep])

# Freeze condition levels as they appear (stable ordering in outputs)
cond_levels <- unique(obj$Condition)
obj$Condition <- factor(obj$Condition, levels = cond_levels)

message("Conditions found: ", paste(levels(obj$Condition), collapse = ", "))
message("N samples (unique IDs): ", length(unique(obj$ID)))

# -----------------------------
# Optional QC / sanity plots
# -----------------------------
# Only run if these columns/reductions exist (doesn't error if not)
if ("umap.harmony" %in% names(obj@reductions)) {
  pdf(file.path(PLOT_DIR, "umap_by_condition.pdf"), width = 7, height = 6)
  print(DimPlot(obj, reduction = "umap.harmony", group.by = "Condition", raster = FALSE) +
          ggtitle("UMAP by Condition"))
  dev.off()

  if ("Sex" %in% colnames(obj@meta.data)) {
    pdf(file.path(PLOT_DIR, "umap_by_sex.pdf"), width = 7, height = 6)
    print(DimPlot(obj, reduction = "umap.harmony", group.by = "Sex", raster = FALSE) +
            ggtitle("UMAP by Sex"))
    dev.off()
  }
}

# If you track doublets with a column like DF/Doublet/etc, plot if present
if ("DF" %in% colnames(obj@meta.data) && "umap.harmony" %in% names(obj@reductions)) {
  pdf(file.path(PLOT_DIR, "umap_by_doublet_flag.pdf"), width = 7, height = 6)
  print(DimPlot(obj, reduction = "umap.harmony", group.by = "DF", raster = FALSE) +
          ggtitle("UMAP by Doublet Flag"))
  dev.off()
}

# -----------------------------
# Helper: run propeller and write results
# -----------------------------
run_propeller <- function(cluster_vec, sample_vec, group_vec,
                          out_prefix,
                          transform = "logit",
                          run_pairwise = TRUE) {

  # speckle expects vectors of same length
  stopifnot(length(cluster_vec) == length(sample_vec))
  stopifnot(length(sample_vec) == length(group_vec))

  props <- getTransformedProps(clusters = cluster_vec,
                               sample   = sample_vec,
                               transform = transform)

  # Ensure group factor aligns with props sample columns
  # props$Proportions columns are samples; sample_vec is per-cell.
  # Build a sample-level group vector in the same order as colnames(props$Proportions).
  sample_order <- colnames(props$Proportions)
  sample_to_group <- data.frame(
    ID = sample_vec,
    Group = group_vec,
    stringsAsFactors = FALSE
  ) %>%
    distinct(ID, Group)

  group_by_sample <- sample_to_group$Group[match(sample_order, sample_to_group$ID)]
  if (any(is.na(group_by_sample))) {
    stop("Some samples in props are missing group labels. Check ID/Condition mapping.")
  }

  group_by_sample <- factor(group_by_sample, levels = levels(group_vec))
  props$group <- group_by_sample

  # Build design
  design <- model.matrix(~ 0 + group, data = props)
  colnames(design) <- gsub("^group", "", colnames(design))  # cleaner column names

  # Decide test:
  n_groups <- nlevels(props$group)

  results <- list()

  if (n_groups == 2) {
    # Single contrast: group2 - group1
    lv <- levels(props$group)
    contrast_str <- paste0(lv[2], "-", lv[1])
    contr <- makeContrasts(contrasts = contrast_str, levels = design)

    res <- as.data.frame(
      propeller.ttest(props, design, contrasts = contr,
                      robust = TRUE, trend = FALSE, sort = TRUE)
    )
    results[[paste0(lv[2], "_vs_", lv[1])]] <- res

  } else {
    # ANOVA-style test across all groups
    res <- as.data.frame(
      propeller.anova(props, design,
                      robust = TRUE, trend = FALSE, sort = TRUE)
    )
    results[["ANOVA_all_groups"]] <- res

    # Optional: all pairwise contrasts
    if (run_pairwise) {
      lv <- levels(props$group)
      pair_names <- combn(lv, 2, simplify = FALSE)

      for (p in pair_names) {
        contrast_str <- paste0(p[2], "-", p[1])
        contr <- makeContrasts(contrasts = contrast_str, levels = design)

        res_pw <- as.data.frame(
          propeller.ttest(props, design, contrasts = contr,
                          robust = TRUE, trend = FALSE, sort = TRUE)
        )
        results[[paste0(p[2], "_vs_", p[1])]] <- res_pw
      }
    }
  }

  # Write XLSX (one sheet per result)
  wb <- createWorkbook()
  for (nm in names(results)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, results[[nm]], rowNames = TRUE)
  }
  xlsx_path <- file.path(TAB_DIR, paste0(out_prefix, "_propeller_results.xlsx"))
  saveWorkbook(wb, file = xlsx_path, overwrite = TRUE)

  # Return both props and results for plotting
  return(list(props = props, results = results, design = design))
}

# -----------------------------
# Helper: barplot mean +/- SE from props
# -----------------------------
plot_props_bar <- function(props_obj, out_pdf,
                           x_label = "Group",
                           y_label = "Proportion of nuclei",
                           title = NULL,
                           rotate_x = TRUE) {

  # props_obj$Proportions: rows = clusters, cols = samples
  props_long <- as.data.frame(props_obj$Proportions) %>%
    tibble::rownames_to_column("Label") %>%
    pivot_longer(cols = -Label, names_to = "ID", values_to = "Proportion") %>%
    left_join(
      data.frame(ID = colnames(props_obj$Proportions),
                 Group = props_obj$group,
                 stringsAsFactors = FALSE),
      by = "ID"
    )

  # Summarize mean/se per Label x Group
  summary_df <- props_long %>%
    group_by(Label, Group) %>%
    summarise(
      mean_value = mean(Proportion),
      se_value   = sd(Proportion) / sqrt(n()),
      .groups = "drop"
    )

  summary_df$Group <- factor(summary_df$Group, levels = levels(props_obj$group))

  p <- ggplot(summary_df, aes(x = Label, y = mean_value, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    geom_errorbar(aes(ymin = pmax(mean_value - se_value, 0),
                      ymax = mean_value + se_value),
                  position = position_dodge(width = 0.8),
                  width = 0.25) +
    labs(x = NULL, y = y_label, fill = x_label, title = title) +
    theme_minimal(base_size = 14)

  if (rotate_x) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  pdf(out_pdf, width = 11, height = 7)
  print(p)
  dev.off()
}

# ============================================================
# Part 1: Cluster-level proportion analysis
# ============================================================

# You can change this if your “cluster identity” column differs
# We keep it exactly like your prior scripts: seurat_clusters
stopifnot("seurat_clusters" %in% colnames(obj@meta.data))

cluster_res <- run_propeller(
  cluster_vec = obj$seurat_clusters,
  sample_vec  = obj$ID,
  group_vec   = obj$Condition,
  out_prefix  = "cluster",
  transform   = "logit",
  run_pairwise = RUN_PAIRWISE_IF_MULTI
)

plot_props_bar(
  props_obj = cluster_res$props,
  out_pdf   = file.path(PLOT_DIR, "barplot_cluster_proportions_by_condition.pdf"),
  x_label   = "Condition",
  y_label   = "Proportion of nuclei",
  title     = "Cluster proportions by condition",
  rotate_x  = TRUE
)

# ============================================================
# Part 2: Cell-type-level proportion analysis
# ============================================================

stopifnot("cell_type" %in% colnames(obj@meta.data))

celltype_res <- run_propeller(
  cluster_vec = obj$cell_type,
  sample_vec  = obj$ID,
  group_vec   = obj$Condition,
  out_prefix  = "cell_type",
  transform   = "logit",
  run_pairwise = RUN_PAIRWISE_IF_MULTI
)

plot_props_bar(
  props_obj = celltype_res$props,
  out_pdf   = file.path(PLOT_DIR, "barplot_celltype_proportions_by_condition.pdf"),
  x_label   = "Condition",
  y_label   = "Proportion of nuclei",
  title     = "Cell-type proportions by condition",
  rotate_x  = TRUE
)

message("Done. Outputs written to: ", OUT_DIR)
