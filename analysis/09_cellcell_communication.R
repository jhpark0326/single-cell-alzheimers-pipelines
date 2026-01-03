# Title: Cell-cell communication (CellChat)
# Project: AD snRNA-seq multi-model analysis
# Author: Jung Hyun Park
# Purpose: Run CellChat per Condition by clusters,
#          then merge datasets for differential interaction analysis.

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(fs)
  library(future)
  library(ComplexHeatmap)
})

options(stringsAsFactors = FALSE)

# -----------------------
# Parameters (edit here)
# -----------------------
in_rds <- "outputs/annotation/annotated.rds"

out_dir <- "outputs/cellchat"
fig_dir <- file.path(out_dir, "figures")
dir_create(out_dir)
dir_create(fig_dir)

# CellChat species DB
SPECIES <- "mouse"  # "mouse" or "human"

# Label used as the "cell group" in CellChat (must exist in obj@meta.data)
# options: "cluster_label" or "cell_type"
LABEL_COL <- "cluster_label"

# Condition column (dataset split)
CONDITION_COL <- "Condition"

# Minimum cells per group to keep communications
MIN_CELLS <- 10

# Communication probability estimator
COMM_TYPE <- "triMean"  # common choice

# Merge comparison: choose 2 conditions by name (must exist). If NULL, uses first + last.
COMPARE_CONDITIONS <- c("Control", "Disease")  # set to NULL to auto

# DEG mapping (optional)
RUN_DEG_MAPPING <- TRUE
POS_DATASET <- NULL   # e.g., "Disease" (if NULL, uses COMPARE_CONDITIONS[2])
DEA_THRESH_PC <- 0.1
DEA_THRESH_FC <- 0.1
DEA_THRESH_P  <- 0.05
LR_LOGFC_CUTOFF <- 0.1

# Memory
options(future.globals.maxSize = 4 * 1024^3)
plan("sequential")

# -----------------------
# Load object + validate metadata
# -----------------------
if (!file.exists(in_rds)) stop("Missing input: ", in_rds)
obj <- readRDS(in_rds)
DefaultAssay(obj) <- "RNA"

meta_cols <- colnames(obj[[]])
req <- c(LABEL_COL, CONDITION_COL)
missing <- setdiff(req, meta_cols)
if (length(missing) > 0) {
  stop("Object metadata missing required columns: ", paste(missing, collapse = ", "))
}

obj[[LABEL_COL]] <- as.character(obj[[LABEL_COL]][, 1])
obj[[CONDITION_COL]] <- as.character(obj[[CONDITION_COL]][, 1])

# Use label as CellChat grouping (Idents)
Idents(obj) <- LABEL_COL

# Conditions to run
conditions <- sort(unique(obj[[CONDITION_COL]]))
if (length(conditions) < 2) stop("Need at least 2 conditions in ", CONDITION_COL)

if (is.null(COMPARE_CONDITIONS)) {
  COMPARE_CONDITIONS <- c(conditions[1], conditions[length(conditions)])
}
if (!all(COMPARE_CONDITIONS %in% conditions)) {
  stop("COMPARE_CONDITIONS not found in object: ",
       paste(setdiff(COMPARE_CONDITIONS, conditions), collapse = ", "))
}

message("Conditions: ", paste(conditions, collapse = ", "))
message("Comparing: ", paste(COMPARE_CONDITIONS, collapse = " vs "))

# -----------------------
# Helper: run CellChat for one condition
# -----------------------
run_cellchat_one <- function(seu, cond_name) {
  message("Running CellChat for condition: ", cond_name)

  seu_sub <- subset(seu, subset = !!as.name(CONDITION_COL) == cond_name)
  DefaultAssay(seu_sub) <- "RNA"
  Idents(seu_sub) <- LABEL_COL

  # Create CellChat object
  cellchat <- createCellChat(object = seu_sub, group.by = "ident", assay = "RNA")

  # Set DB
  if (SPECIES == "mouse") {
    cellchat@DB <- CellChatDB.mouse
  } else if (SPECIES == "human") {
    cellchat@DB <- CellChatDB.human
  } else {
    stop("SPECIES must be 'mouse' or 'human'")
  }

  # Standard pipeline
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = COMM_TYPE)
  cellchat <- filterCommunication(cellchat, min.cells = MIN_CELLS)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  cellchat
}

# -----------------------
# Run per-condition and save
# -----------------------
object.list <- list()

for (cond in conditions) {
  cc <- run_cellchat_one(obj, cond)
  object.list[[cond]] <- cc
  saveRDS(cc, file.path(out_dir, paste0("cellchat__", gsub("[^A-Za-z0-9]+", "_", cond), ".rds")))
}

# -----------------------
# Lift (align groups) across datasets and merge
# -----------------------
group.new <- levels(object.list[[COMPARE_CONDITIONS[1]]]@idents)

# Lift any condition missing some labels to the reference label set
for (nm in names(object.list)) {
  if (!identical(levels(object.list[[nm]]@idents), group.new)) {
    object.list[[nm]] <- liftCellChat(object.list[[nm]], group.new = group.new)
  }
}

cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(object.list, file.path(out_dir, "cellchat_object_list.rds"))
saveRDS(cellchat_merged, file.path(out_dir, "cellchat_merged.rds"))

# -----------------------
# Core comparisons
# -----------------------
# 1) Total number of interactions + strength
p1 <- compareInteractions(cellchat_merged, show.legend = FALSE)
p2 <- compareInteractions(cellchat_merged, show.legend = FALSE, measure = "weight")

pdf(file.path(fig_dir, "compareInteractions__count_and_weight.pdf"), width = 10, height = 5)
print(p1 + p2)
dev.off()

# 2) Differential interaction circle plots (two-condition comparison)
cmp_idx <- match(COMPARE_CONDITIONS, names(object.list))
if (any(is.na(cmp_idx))) stop("Internal: comparison conditions not found in object.list names.")

pdf(file.path(fig_dir, paste0("diffInteraction__", COMPARE_CONDITIONS[1], "_vs_", COMPARE_CONDITIONS[2], ".pdf")),
    width = 10, height = 5)
netVisual_diffInteraction(cellchat_merged, comparison = cmp_idx, weight.scale = TRUE)
netVisual_diffInteraction(cellchat_merged, comparison = cmp_idx, weight.scale = TRUE, measure = "weight")
dev.off()

# 3) Differential heatmaps (count + weight)
pdf(file.path(fig_dir, paste0("diffHeatmap__", COMPARE_CONDITIONS[1], "_vs_", COMPARE_CONDITIONS[2], ".pdf")),
    width = 11, height = 10)
g1 <- netVisual_heatmap(cellchat_merged)
g2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
print(g1 + g2)
dev.off()

# 4) Compare outgoing / incoming signaling role heatmaps (union pathways)
pathway.union <- Reduce(union, lapply(object.list[c(COMPARE_CONDITIONS[1], COMPARE_CONDITIONS[2])], function(x) x@netP$pathways))

ht_out_1 <- netAnalysis_signalingRole_heatmap(object.list[[COMPARE_CONDITIONS[1]]],
                                             pattern = "outgoing", signaling = pathway.union,
                                             title = COMPARE_CONDITIONS[1], width = 10, height = 18)
ht_out_2 <- netAnalysis_signalingRole_heatmap(object.list[[COMPARE_CONDITIONS[2]]],
                                             pattern = "outgoing", signaling = pathway.union,
                                             title = COMPARE_CONDITIONS[2], width = 10, height = 18)

pdf(file.path(fig_dir, paste0("signalingRoleHeatmap__outgoing__", COMPARE_CONDITIONS[1], "_vs_", COMPARE_CONDITIONS[2], ".pdf")),
    width = 12, height = 10)
ComplexHeatmap::draw(ht_out_1 + ht_out_2, ht_gap = grid::unit(0.5, "cm"))
dev.off()

ht_in_1 <- netAnalysis_signalingRole_heatmap(object.list[[COMPARE_CONDITIONS[1]]],
                                            pattern = "incoming", signaling = pathway.union,
                                            title = COMPARE_CONDITIONS[1], width = 10, height = 18)
ht_in_2 <- netAnalysis_signalingRole_heatmap(object.list[[COMPARE_CONDITIONS[2]]],
                                            pattern = "incoming", signaling = pathway.union,
                                            title = COMPARE_CONDITIONS[2], width = 10, height = 18)

pdf(file.path(fig_dir, paste0("signalingRoleHeatmap__incoming__", COMPARE_CONDITIONS[1], "_vs_", COMPARE_CONDITIONS[2], ".pdf")),
    width = 12, height = 10)
ComplexHeatmap::draw(ht_in_1 + ht_in_2, ht_gap = grid::unit(0.5, "cm"))
dev.off()

# 5) Compare information flow (rankNet)
pdf(file.path(fig_dir, paste0("rankNet__flow__", COMPARE_CONDITIONS[1], "_vs_", COMPARE_CONDITIONS[2], ".pdf")),
    width = 10, height = 10)
gg_flow1 <- rankNet(cellchat_merged, mode = "comparison", comparison = cmp_idx,
                    measure = "weight", stacked = TRUE, do.stat = TRUE)
gg_flow2 <- rankNet(cellchat_merged, mode = "comparison", comparison = cmp_idx,
                    measure = "weight", stacked = FALSE, do.stat = TRUE)
print(gg_flow1 + gg_flow2)
dev.off()

# -----------------------
# Optional: DEG mapping to extract up/down LR pairs between two datasets
# -----------------------
if (RUN_DEG_MAPPING) {
  pos.dataset <- if (!is.null(POS_DATASET)) POS_DATASET else COMPARE_CONDITIONS[2]
  if (!(pos.dataset %in% names(object.list))) stop("POS_DATASET not found: ", pos.dataset)

  features.name <- paste0(gsub("[^A-Za-z0-9]+", "_", pos.dataset), ".DEA")

  cellchat_merged <- identifyOverExpressedGenes(
    cellchat_merged,
    idents.use = cellchat_merged@idents$joint,
    group.dataset = "datasets",
    pos.dataset = pos.dataset,
    features.name = features.name,
    only.pos = FALSE,
    thresh.pc = DEA_THRESH_PC,
    thresh.fc = DEA_THRESH_FC,
    thresh.p = DEA_THRESH_P,
    group.DE.combined = FALSE
  )

  net.de <- netMappingDEG(cellchat_merged, features.name = features.name, variable.all = TRUE)
  saveRDS(net.de, file.path(out_dir, paste0("netMappingDEG__", features.name, ".rds")))

  net.up <- subsetCommunication(cellchat_merged, net = net.de, datasets = pos.dataset,
                                ligand.logFC = LR_LOGFC_CUTOFF, receptor.logFC = LR_LOGFC_CUTOFF)
  net.down <- subsetCommunication(cellchat_merged, net = net.de, datasets = COMPARE_CONDITIONS[1],
                                  ligand.logFC = LR_LOGFC_CUTOFF, receptor.logFC = LR_LOGFC_CUTOFF)

  write.csv(net.up, file.path(out_dir, paste0("LR_up__", pos.dataset, ".csv")), row.names = FALSE)
  write.csv(net.down, file.path(out_dir, paste0("LR_down__", pos.dataset, ".csv")), row.names = FALSE)

  message("Saved DEG-mapped LR tables to: ", out_dir)
}

# -----------------------
# Export per-dataset communication tables (net + netP)
# -----------------------
dir_create(file.path(out_dir, "tables"))

df.net <- subsetCommunication(cellchat_merged)
for (nm in names(df.net)) {
  write.csv(df.net[[nm]],
            file.path(out_dir, "tables", paste0("net__", gsub("[^A-Za-z0-9]+", "_", nm), ".csv")),
            row.names = FALSE)
}

df.netP <- subsetCommunication(cellchat_merged, slot.name = "netP")
for (nm in names(df.netP)) {
  write.csv(df.netP[[nm]],
            file.path(out_dir, "tables", paste0("netP__", gsub("[^A-Za-z0-9]+", "_", nm), ".csv")),
            row.names = FALSE)
}

message("Done.")
message("Saved objects: ", out_dir)
message("Saved figures: ", fig_dir)
