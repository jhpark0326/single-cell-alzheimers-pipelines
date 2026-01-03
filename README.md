# Single-cell transcriptomics analysis pipelines for Alzheimer's disease  

## Project
This repository demonstrates an end-to-end single-cell transcriptomics analysis workflow,
with examples drawn from Alzheimer’s disease mouse model datasets generated 
using 10x Genomics single-nucleus RNA-seq (snRNA-seq).

## Author
Jung Hyun Park

## Purpose
This repository provides a modular and reproducible collection of analysis scripts designed to:
- Standardize common single-cell transcriptomics analysis steps
- Keep each analytical stage transparent and reusable

## Notes
- Sensitive, unpublished, or identifying experimental details are intentionally omitted.
- Metadata fields, condition labels, and sample identifiers are generalized
  (e.g., `sample_id`, `Condition`, `Batch`).
- Scripts are written to prioritize interpretability and adaptability.

---

## Repository structure

```text
analysis/
├── 01_qc_preprocessing.R
├── 02_integration.R
├── 03_clustering.R
├── 04_markers.R
├── 05_annotation.R
├── 06_proportion_analysis.R
├── 07_differential_expression.R
├── 08_network_analysis_hdWGCNA.R
├── 09_cellcell_communication.R
└── 10_pyscenic/
    ├── 10a_pyscenic_export_from_r.R
    ├── 10b_pyscenic_run_and_visualize.py
    └── 10c_pyscenic_slurm_grn_ctx_aucell.sh

---

## Pipeline overview

### 01. Quality control and preprocessing (SoupX + DoubletFinder)
- Filtering low-quality cells or nuclei
- Computing and visualizing QC metrics
- Generating cleaned objects for downstream analysis

**Typical outputs**
- Filtered single-cell objects
- QC plots and summary statistics

---

### 02. Integration (Harmony)
- Integration across samples or batches (when applicable)
- Selection of dimensionality for downstream analyses

**Typical outputs**
- Integrated representations
- Reduced dimensional embeddings

---

### 03. Clustering (Clustree)
- Construction of nearest-neighbor graphs
- Clustering and dimensionality reduction (e.g., UMAP)

**Typical outputs**
- Cluster assignments
- Embedding visualizations

---

### 04. Marker identification (FindAllMarkers)
- Identification of cluster-enriched genes
- Visualization of marker expression

**Typical outputs**
- Marker tables
- Feature and dot plots

---

### 05. Annotation (automated + manual mapping)
- Automated assignment of cell-type labels based on Allen Brain Cell Atlas taxonomy using Semi-automated Hand Annotation for Single Cell and Spatial Datasets (https://github.com/neurogenetics/SAHA)
- Validation using known markers

**Typical outputs**
- Annotated metadata
- Annotation figures

---

### 06. Proportion analysis (Speckle propeller)
- Quantification of cell-type or cluster proportions
- Comparison across samples or conditions

**Typical outputs**
- Proportion tables
- Summary plots and statistical results

---

### 07. Differential expression + enrichment analysis (enrichR)
- Differential expression testing between defined conditions
- Emphasis on biologically interpretable contrasts

**Typical outputs**
- Differential expression tables
- Volcano and summary plots
- Enrichment terms

---

### 08. Network analysis (hdWGCNA)
- Construction of gene co-expression networks
- Identification of co-expression modules

**Typical outputs**
- Gene modules
- Module-level summaries and visualizations

---

### 09. Cell–cell communication (CellChat)
- Inference of intercellular signaling interactions
- Pathway-level communication summaries

**Typical outputs**
- Interaction tables
- Communication network visualizations

---

### 10. Regulon inference (pySCENIC)
This module is divided into three components:

- **10a**: Export of expression data and metadata from R into pySCENIC-compatible formats
- **10b**: Import the Seurat object into Python and visualize pySCENIC results 
- **10c**: Slurm batch script to run GRN inference, cisTarget motif enrichment,
  and AUCell scoring on high-performance computing systems

**Typical outputs**
- Regulon definitions
- AUCell activity matrices
- Regulon activity visualizations

---

## Environment

Analyses were developed and tested using:

### R
- R ≥ 4.3
- Key packages:
  - Seurat (v5 ecosystem)
  - tidyverse
  - ggplot2 / patchwork
  - hdWGCNA

### Python (pySCENIC module)
- Python ≥ 3.10
- pyscenic
- scanpy / anndata
- loompy (for `.loom` export when required)

### Computing
- Local execution for most R-based analyses
- Python-based pySCENIC workflows executed on high-performance
  computing systems (Slurm)

Exact package versions are not pinned to keep the repository lightweight and readable.

---

## Reproducibility considerations
- Scripts are modular and deterministic where possible
- Random seeds are set explicitly for stochastic steps
- Hard-coded system paths are avoided
- Metadata fields are standardized across pipeline stages

---

## Contact
For questions, suggestions, or reuse of this workflow, please contact:
Jung Hyun Park
