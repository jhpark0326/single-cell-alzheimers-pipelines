# Title: pySCENIC GRN + cisTarget + AUCell (Slurm job)
# Project: AD snRNA-seq multi-model analysis
# Author: Jung Hyun Park
# Purpose: Run pySCENIC steps on HPC:
#     (1) GRN inference (grn) -> adj.tsv
#     (2) cisTarget regulon pruning (ctx) -> regulons.csv
#     (3) AUCell regulon activity scoring (aucell) -> pyscenic_output.loom
# Notes:
#   - Expects a pySCENIC input loom created by 10b (genes x cells).
#   - You MUST provide the correct TF list and cisTarget DB paths for your organism.

#!/bin/bash
#####  Constructed by HPC everywhere #####
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=0-24:00:00
#SBATCH --mem=100gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=pyscenic_grn_ctx_auc
#SBATCH -A YOUR_ALLOCATION
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

# -----------------------
# Modules / environment
# -----------------------
module load python/3.11.4
source /PATH/TO/YOUR_VENV/bin/activate

# -----------------------
# Parameters (edit here)
# -----------------------
WORKDIR="/PATH/TO/PROJECT/outputs/pyscenic"
cd "$WORKDIR"

# Inputs
TF_LIST="allTFs_mm.txt"                 # transcription factors list (mouse example)
LOOM_IN="filtered_scenic_input.loom"    # from 10b

# Outputs
ADJ_OUT="adj.tsv"
REG_OUT="regulons.csv"
AUC_OUT="pyscenic_output.loom"

# cisTarget databases (feather files)
DB_DIR="/PATH/TO/cistargetDB"
DB_GLOB="${DB_DIR}/*.feather"
DBS=( $DB_GLOB )

# Motif annotation table
MOTIF_ANNOT="${DB_DIR}/motifs.tbl"

# Workers
NWORKERS="${SLURM_CPUS_PER_TASK:-20}"

echo "Working directory: $WORKDIR"
echo "Workers: $NWORKERS"

# -----------------------
# Step 1: GRN inference
# -----------------------
pyscenic grn "$LOOM_IN" "$TF_LIST" \
  -o "$ADJ_OUT" \
  --num_workers "$NWORKERS"

# -----------------------
# Step 2-3: cisTarget pruning (ctx)
# -----------------------
pyscenic ctx "$ADJ_OUT" "${DBS[@]}" \
  --annotations_fname "$MOTIF_ANNOT" \
  --expression_mtx_fname "$LOOM_IN" \
  --output "$REG_OUT" \
  --mask_dropouts \
  --num_workers "$NWORKERS"

# -----------------------
# Step 4: AUCell scoring
# -----------------------
pyscenic aucell "$LOOM_IN" "$REG_OUT" \
   --output "$AUC_OUT" \
   --num_workers "$NWORKERS"

echo "Done."
echo "Adjacencies: $ADJ_OUT"
echo "Regulons:    $REG_OUT"
echo "AUCell loom: $AUC_OUT"

