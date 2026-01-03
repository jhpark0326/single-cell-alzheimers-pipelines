# Title: high-dimensional Weighted Gene Co-expression Network Analysis (hdWGCNA) by clusters or cell types
# Project: AD snRNA-seq multi-model analysis
# Author: Jung Hyun Park
# Purpose:
#   1) Build hdWGCNA co-expression modules using metacells
#   2) Compute harmonized module eigengenes (hMEs) and hub genes scores (kMEs)
#   3) Run DME (Differential Module Eigengene) tests across conditions within each annotation group
#   4) Run module enrichment analysis (Enrichr) per module

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(WGCNA)
  library(hdWGCNA)
  library(fs)
  library(ggplot2)
})

theme_set(theme_cowplot())
set.seed(12345)

# -----------------------
# Parameters (edit here)
# -----------------------
in_rds <- "outputs/annotation/annotated.rds"

out_dir <- "outputs/hdwgcna"
fig_dir <- file.path(out_dir, "figures")
net_dir <- file.path(out_dir, "networks")
dir_create(out_dir)
dir_create(fig_dir)
dir_create(net_dir)

# Which annotation to run WGCNA within (metacell identity)
# options should match your annotated object, e.g. "cluster_label" or "cell_type"
GROUP_COL <- "cluster_label"

# Sample identifier column (technical); used for metacells and harmonized MEs
SAMPLE_COL <- "sample_id"

# Condition column (biological); used ONLY for DME contrasts
CONDITION_COL <- "Condition"

# WGCNA experiment name saved into obj@misc
WGCNA_NAME <- "hdWGCNA_by_group"

# Expression prep for WGCNA (metacells)
ASSAY_FOR_WGCNA <- "RNA"
SLOT_FOR_WGCNA  <- "data"  # normalized data

# Gene selection for WGCNA
GENE_SELECT_METHOD <- "fraction"
GENE_FRACTION <- 1e-7  # very permissive; tune upward if needed (e.g., 1e-4)

# Metacell construction
# (k=20~75 typical; reduce if you have small groups)
METACELL_K <- 35
METACELL_MAX_SHARED <- 10

# Soft power / network type
NETWORK_TYPE <- "signed"

# Eigengenes
ME_GROUP_BY <- SAMPLE_COL   # harmonize by sample/library

# Rename modules to generic names (Mod1, Mod2, ...)
MODULE_PREFIX <- "Mod"

# Enrichr settings (module enrichment)
RUN_ENRICHR <- TRUE
ENRICHR_DBS <- c(
  "GO_Biological_Process_2025",
  "KEGG_2019_Mouse",
  "Reactome_Pathways_2024",
  "WikiPathways_2024_Mouse"
)
ENRICHR_MAX_GENES <- Inf     # per module; set Inf to use all genes

# DME settings (module-level differential test)
RUN_DME <- TRUE
DME_TEST <- "wilcox"
DME_PADJ <- 0.05

# Define baseline + targets (same pattern as 07_DE script)
BASELINE_CONDITION <- "Control"
TARGET_CONDITIONS <- NULL  # NULL => auto-detect all non-baseline conditions

# -----------------------
# Load object + validate metadata
# -----------------------
if (!file.exists(in_rds)) stop("Missing input: ", in_rds)
obj <- readRDS(in_rds)

meta_cols <- colnames(obj[[]])
req <- c(GROUP_COL, SAMPLE_COL, CONDITION_COL)
missing <- setdiff(req, meta_cols)
if (length(missing) > 0) {
  stop("Object metadata missing required columns: ", paste(missing, collapse = ", "))
}

obj[[GROUP_COL]]     <- as.character(obj[[GROUP_COL]][, 1])
obj[[SAMPLE_COL]]    <- as.character(obj[[SAMPLE_COL]][, 1])
obj[[CONDITION_COL]] <- as.character(obj[[CONDITION_COL]][, 1])

# Keep a stable “group_use” column for hdWGCNA
obj$group_use <- obj[[GROUP_COL]] # either cluster or cell type
obj$sample_use <- obj[[SAMPLE_COL]]
obj$condition_use <- obj[[CONDITION_COL]]

# Conditions / contrasts
all_conditions <- sort(unique(obj$condition_use))
if (!(BASELINE_CONDITION %in% all_conditions)) {
  stop("BASELINE_CONDITION not found in ", CONDITION_COL, ": ", BASELINE_CONDITION)
}
if (is.null(TARGET_CONDITIONS)) {
  TARGET_CONDITIONS <- setdiff(all_conditions, BASELINE_CONDITION)
}
if (length(TARGET_CONDITIONS) == 0) stop("No non-baseline TARGET_CONDITIONS found.")

contrasts <- lapply(TARGET_CONDITIONS, function(x) c(x, BASELINE_CONDITION))

# -----------------------
# Basic RNA prep 
# -----------------------
DefaultAssay(obj) <- ASSAY_FOR_WGCNA

# Seurat v5: JoinLayers to ensure RNA is a single layer for downstream steps
if ("RNA" %in% names(obj@assays)) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

# Ensure normalized data exists (idempotent-ish)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)

# Use an existing reduction if present; otherwise run PCA/UMAP
if (!("pca" %in% names(obj@reductions))) {
  obj <- RunPCA(obj, verbose = FALSE)
}
if (!("umap" %in% names(obj@reductions))) {
  obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", verbose = FALSE)
}

pdf(file.path(fig_dir, "umap_by_group.pdf"), width = 8, height = 6)
print(DimPlot(obj, group.by = "group_use", label = TRUE, raster = FALSE) + NoLegend())
dev.off()

pdf(file.path(fig_dir, "umap_by_condition.pdf"), width = 8, height = 6)
print(DimPlot(obj, group.by = "condition_use", raster = FALSE))
dev.off()

# -----------------------
# Setup hdWGCNA
# -----------------------
obj <- SetupForWGCNA(
  obj,
  gene_select = GENE_SELECT_METHOD,
  fraction = GENE_FRACTION,
  wgcna_name = WGCNA_NAME
)

# Choose a reduction for metacell KNN
reduction_use <- if ("harmony" %in% names(obj@reductions)) "harmony" else if ("pca" %in% names(obj@reductions)) "pca" else "pca"

# Construct metacells within (group × sample)
obj <- MetacellsByGroups(
  obj,
  group.by = c("group_use", "sample_use"),
  reduction = reduction_use,
  k = METACELL_K,
  max_shared = METACELL_MAX_SHARED,
  ident.group = "group_use"
)

# Normalize metacell expression
obj <- NormalizeMetacells(obj)

# -----------------------
# Co-expression network
# -----------------------
Idents(obj) <- "group_use"

# SetDatExpr: use metacells expression matrix
# group.by=NULL => uses Idents(obj) (group_use)
obj <- SetDatExpr(
  obj,
  group_name = NULL,
  group.by = NULL,
  assay = ASSAY_FOR_WGCNA,
  slot = SLOT_FOR_WGCNA
)

# Soft power selection
obj <- TestSoftPowers(obj, networkType = NETWORK_TYPE)
plot_list <- PlotSoftPowers(obj)

pdf(file.path(fig_dir, "soft_powers.pdf"), width = 10, height = 6)
print(wrap_plots(plot_list, ncol = 2))
dev.off()

power_table <- GetPowerTable(obj)
readr::write_csv(power_table, file.path(out_dir, "soft_power_table.csv"))

# Construct network + TOM
obj <- ConstructNetwork(
  obj,
  tom_outdir = file.path(net_dir, "TOM"),
  tom_name = WGCNA_NAME,
  overwrite_tom = TRUE
)

pdf(file.path(fig_dir, "dendrogram.pdf"), width = 10, height = 6)
print(PlotDendrogram(obj, main = paste0(WGCNA_NAME, " Dendrogram")))
dev.off()

# -----------------------
# Module eigengenes + connectivity
# -----------------------
obj <- ModuleEigengenes(
  obj,
  group.by.vars = ME_GROUP_BY
)

hMEs <- GetMEs(obj)                  # harmonized
MEs  <- GetMEs(obj, harmonized=FALSE)

readr::write_csv(as.data.frame(hMEs), file.path(out_dir, "hMEs_harmonized.csv"))
readr::write_csv(as.data.frame(MEs),  file.path(out_dir, "MEs_raw.csv"))

obj <- ModuleConnectivity(
  obj,
  group.by = "group_use"
)

obj <- ResetModuleNames(
  obj,
  new_name = MODULE_PREFIX
)

modules_tbl <- GetModules(obj) %>% filter(module != "grey")
readr::write_csv(modules_tbl, file.path(out_dir, "module_assignments.csv"))

hub_df <- GetHubGenes(obj, n_hubs = 50)
readr::write_csv(hub_df, file.path(out_dir, "hub_genes_top50.csv"))

pdf(file.path(fig_dir, "kME_ranked.pdf"), width = 12, height = 8)
print(PlotKMEs(obj, ncol = 4, text_size = 1))
dev.off()

# Add hMEs to metadata so you can use Seurat plotting functions
obj@meta.data <- cbind(obj@meta.data, hMEs)

# Basic module feature plots
pdf(file.path(fig_dir, "module_featureplots_hMEs.pdf"), width = 12, height = 10)
plot_hME <- ModuleFeaturePlot(
  obj,
  features = "hMEs",
  order = TRUE,
  reduction = if ("umap.harmony" %in% names(obj@reductions)) "umap.harmony" else "umap"
)
plot_hME <- lapply(plot_hME, function(x) x + theme(plot.title = element_text(size = 7)))
print(wrap_plots(plot_hME, ncol = 4))
dev.off()

# DotPlot of hMEs by group (helps see module activation patterns)
mods <- colnames(hMEs)
mods <- mods[mods != "grey"]

pdf(file.path(fig_dir, "dotplot_hMEs_by_group.pdf"), width = 10, height = 8)
p <- DotPlot(obj, features = mods, group.by = "group_use", dot.scale = 5) +
  coord_flip() +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8))
print(p)
dev.off()

# -----------------------
# Module enrichment (Enrichr via hdWGCNA)
# -----------------------
if (RUN_ENRICHR) {
  suppressPackageStartupMessages(library(enrichR))
  setEnrichrSite("Enrichr")

  obj <- RunEnrichr(
    obj,
    dbs = ENRICHR_DBS,
    max_genes = ENRICHR_MAX_GENES,
    wgcna_name = WGCNA_NAME
  )

  enr_tbl <- GetEnrichrTable(obj, wgcna_name = WGCNA_NAME)
  readr::write_csv(enr_tbl, file.path(out_dir, "enrichr_module_results.csv"))

  # Simple barplots per module (built-in)
  enr_plot_dir <- file.path(fig_dir, "enrichr_barplots")
  dir_create(enr_plot_dir)

  # One PDF per module generated by hdWGCNA's EnrichrBarPlot
  # (This uses the enrichr results already stored in the object.)
  EnrichrBarPlot(
    obj,
    outdir = enr_plot_dir,
    n_terms = 10,
    plot_size = c(8, 10),
    logscale = TRUE,
    wgcna_name = WGCNA_NAME
  )
}

# -----------------------
# DME: Differential Module Eigengenes (Condition contrasts within each group)
# -----------------------
if (RUN_DME) {
  dme_dir <- file.path(out_dir, "dme")
  dme_fig <- file.path(fig_dir, "dme")
  dir_create(dme_dir)
  dir_create(dme_fig)

  # Ensure modules list after renaming
  modules <- GetModules(obj, wgcna_name = WGCNA_NAME)
  mod_levels <- levels(modules$module)
  mod_levels <- mod_levels[mod_levels != "grey"]

  groups <- sort(unique(obj$group_use))

  DMEs <- list()

  for (cc in contrasts) {
    cond1 <- cc[1]; cond2 <- cc[2]
    ct_name <- paste0(cond1, "_vs_", cond2)

    message("DME contrast: ", ct_name)

    for (g in groups) {
      # barcodes for group + condition
      bar1 <- rownames(obj@meta.data %>% filter(group_use == g, condition_use == cond1))
      bar2 <- rownames(obj@meta.data %>% filter(group_use == g, condition_use == cond2))

      # If either side empty, skip
      if (length(bar1) == 0 || length(bar2) == 0) next

      cur <- FindDMEs(
        obj,
        barcodes1 = bar1,
        barcodes2 = bar2,
        test.use = DME_TEST,
        wgcna_name = WGCNA_NAME
      )

      cur$group <- g
      cur$contrast <- ct_name
      rownames(cur) <- NULL
      DMEs[[paste0(ct_name, "__", g)]] <- cur
    }
  }

  DME_tbl <- bind_rows(DMEs)
  if (nrow(DME_tbl) == 0) {
    message("No DME results produced. Check Condition/group combinations.")
  } else {
    readr::write_csv(DME_tbl, file.path(dme_dir, "DME_all_groups_all_contrasts.csv"))

    # Heatmap-style tile plot per contrast (modules × groups), showing avg_log2FC where significant
    for (ct in unique(DME_tbl$contrast)) {
      plot_df <- DME_tbl %>%
        filter(contrast == ct) %>%
        mutate(
          module = factor(module, levels = mod_levels),
          group = factor(group, levels = groups),
          avg_log2FC_plot = ifelse(p_val_adj <= DME_PADJ, avg_log2FC, NA_real_),
          sig_star = dplyr::case_when(
            is.na(p_val_adj) ~ "",
            p_val_adj <= 0.001 ~ "***",
            p_val_adj <= 0.01  ~ "**",
            p_val_adj <= 0.05  ~ "*",
            TRUE ~ ""
          )
        )

      # cap for visibility
      cap <- 5 # adjust as needed
      plot_df$avg_log2FC_plot <- pmax(pmin(plot_df$avg_log2FC_plot, cap), -cap)

      p <- ggplot(plot_df, aes(x = module, y = group, fill = avg_log2FC_plot)) +
        geom_tile(color = "white", linewidth = 0.2) +
        geom_text(aes(label = sig_star), size = 3) +
        coord_equal() +
        labs(
          title = paste0("DME: ", ct),
          x = NULL, y = NULL, fill = "avg_log2FC"
        ) +
        theme_minimal(base_size = 10) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()
        )

      pdf(file.path(dme_fig, paste0("DME_heat__", ct, ".pdf")), width = 12, height = 7)
      print(p)
      dev.off()
    }
  }
}

# -----------------------
# Save final object
# -----------------------
saveRDS(obj, file.path(out_dir, paste0(WGCNA_NAME, "_seurat_obj.rds")))

message("Done.")
message("Saved hdWGCNA object: ", file.path(out_dir, paste0(WGCNA_NAME, "_seurat_obj.rds")))
message("Outputs in: ", out_dir)
