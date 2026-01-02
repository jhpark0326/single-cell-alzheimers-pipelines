# Title: Proportion analysis (cluster + cell type)
# Project: AD snRNA-seq multi-model analysis
# Author: Jung Hyun Park
#
# Purpose: Quantify composition changes across experimental conditions at:
#     (1) cluster level 
#     (2) cell-type level 
#   using propeller (speckle) with limma modeling.

suppressPackageStartupMessages({
  library(Seurat)
  library(speckle)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
})

# ==============================
# USER-SET SECTION
# ==============================

INPUT_RDS <- "outputs/annotation/annotated.rds"

OUT_DIR <- "outputs/proportion"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CONTROL_LEVEL <- "Control"  # must match seu$Condition

# Planned contrasts: each condition vs control
RUN_VS_CONTROL <- TRUE
# Optional global test (provided for readers; OFF by default)
RUN_ANOVA_OPTIONAL <- FALSE

# Transformation for propeller
TRANSFORM <- "logit"  # recommended; alternatives: "asin"

# ==============================
# Load object + sanity checks
# ==============================

if (!file.exists(INPUT_RDS)) {
  stop("Cannot find INPUT_RDS: ", INPUT_RDS,
       "\nTip: run this script from repo root, or update INPUT_RDS.")
}

seu <- readRDS(INPUT_RDS)

req_cols <- c("sample_id", "ID", "Condition", "Sex", "Age_months", "cluster_label", "cell_type")
missing <- setdiff(req_cols, colnames(seu@meta.data))
if (length(missing) > 0) {
  stop("Missing required metadata columns: ", paste(missing, collapse = ", "))
}

seu$Condition <- factor(seu$Condition)
if (!(CONTROL_LEVEL %in% levels(seu$Condition))) {
  stop("CONTROL_LEVEL '", CONTROL_LEVEL, "' not found in seu$Condition levels: ",
       paste(levels(seu$Condition), collapse = ", "))
}

# ==============================
# Helper: plotting + excel
# ==============================

plot_mean_se <- function(props, sample_meta, feature_name, out_pdf) {
  df <- as.data.frame(props$Proportions) |>
    tibble::rownames_to_column("feature") |>
    pivot_longer(-feature, names_to = "sample_id", values_to = "proportion") |>
    left_join(sample_meta |> select(sample_id, Condition), by = "sample_id")

  sumdf <- df |>
    group_by(feature, Condition) |>
    summarise(
      mean_value = mean(proportion, na.rm = TRUE),
      se_value = sd(proportion, na.rm = TRUE) / sqrt(sum(!is.na(proportion))),
      .groups = "drop"
    )

  sumdf$feature <- factor(sumdf$feature, levels = unique(sumdf$feature))

  pdf(out_pdf, width = 11, height = 8)
  print(
    ggplot(sumdf, aes(x = feature, y = mean_value, fill = Condition)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
      geom_errorbar(
        aes(ymin = pmax(mean_value - se_value, 0), ymax = mean_value + se_value),
        position = position_dodge(width = 0.8),
        width = 0.25
      ) +
      labs(x = feature_name, y = "Mean proportion (across samples)", fill = "Condition") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
}

write_xlsx <- function(path, sheets) {
  wb <- createWorkbook()
  for (nm in names(sheets)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, sheets[[nm]], rowNames = TRUE)
  }
  saveWorkbook(wb, path, overwrite = TRUE)
}

# ==============================
# Core runner
# ==============================

run_propeller <- function(feature_vec, feature_name, out_prefix) {

  # propeller tests are sample-level, so we need one row per sample_id
  sample_meta <- seu@meta.data |>
    transmute(
      sample_id  = as.character(sample_id),
      Condition  = Condition,
      ID         = as.character(ID),
      Sex        = as.character(Sex),
      Age_months = as.character(Age_months),
      Batch      = if ("Batch" %in% colnames(seu@meta.data)) as.character(Batch) else NA_character_
    ) |>
    group_by(sample_id) |>
    summarise(
      Condition  = first(Condition),
      ID         = first(ID),
      Sex        = first(Sex),
      Age_months = first(Age_months),
      Batch      = first(Batch),
      .groups = "drop"
    )

  # proportions (features x samples)
  props <- getTransformedProps(
    clusters = feature_vec,
    sample   = as.character(seu$sample_id),
    transform = TRANSFORM
  )

  # align sample_meta to propeller’s sample order
  sample_ids <- colnames(props$Proportions)
  sample_meta <- sample_meta |>
    filter(sample_id %in% sample_ids) |>
    mutate(sample_id = factor(sample_id, levels = sample_ids)) |>
    arrange(sample_id) |>
    mutate(sample_id = as.character(sample_id))

  if (!identical(sample_ids, sample_meta$sample_id)) {
    stop("Sample order mismatch between props and sample_meta after alignment.")
  }

  # design: condition-only (simple + readable)
  design <- model.matrix(~ 0 + Condition, data = sample_meta)
  colnames(design) <- levels(sample_meta$Condition)

  results <- list()

  if (RUN_VS_CONTROL) {
    others <- setdiff(colnames(design), CONTROL_LEVEL)
    if (length(others) == 0) stop("No non-control conditions found to compare.")

    contr <- makeContrasts(
      contrasts = paste0(others, " - ", CONTROL_LEVEL),
      levels = design
    )

    res <- propeller.ttest(
      props = props,
      design = design,
      contrasts = contr,
      robust = TRUE,
      trend = FALSE,
      sort = TRUE
    )

    results[["vs_control"]] <- as.data.frame(res)
    write.csv(results[["vs_control"]],
              file.path(OUT_DIR, paste0(out_prefix, "_propeller_vs_control.csv")),
              row.names = TRUE)
  }

  if (RUN_ANOVA_OPTIONAL) {
    an <- propeller.anova(
      clusters = feature_vec,
      sample   = as.character(seu$sample_id),
      group    = seu$Condition
    )
    results[["anova_all_conditions"]] <- as.data.frame(an)
    write.csv(results[["anova_all_conditions"]],
              file.path(OUT_DIR, paste0(out_prefix, "_propeller_anova_all_conditions.csv")),
              row.names = TRUE)
  }
  
  plot_mean_se(
    props, sample_meta,
    feature_name = feature_name,
    out_pdf = file.path(OUT_DIR, paste0(out_prefix, "_barplot_meanSE.pdf"))
  )

  write_xlsx(
    file.path(OUT_DIR, paste0(out_prefix, "_results.xlsx")),
    results
  )

  invisible(results)
}

# ==============================
# Run cluster_label + cell_type
# ==============================

res_cluster_label <- run_propeller(
  feature_vec = as.character(seu$cluster_label),
  feature_name = "Cluster label",
  out_prefix = "cluster_label"
)

res_cell_type <- run_propeller(
  feature_vec = as.character(seu$cell_type),
  feature_name = "Cell type",
  out_prefix = "cell_type"
)

message("Done. Outputs written to: ", OUT_DIR)
