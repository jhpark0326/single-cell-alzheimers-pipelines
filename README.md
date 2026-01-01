# AD snRNA-seq multi-model analysis

This repository demonstrates an end-to-end single-nucleus RNA-seq (snRNA-seq)
analysis pipeline for Alzheimer's disease mouse models using 10x Genomics data.

The analysis integrates multiple AD-related mouse models to:
- perform quality control and preprocessing,
- integrate datasets across batches and conditions,
- identify cell types and subpopulations,
- assess differential expression and network-level changes,
- explore cell–cell communication alterations in disease.

The code is adapted from real research analyses but presented here in a
generalized and reproducible form without raw sequencing data.

## Repository structure

- `analysis/`: ordered analysis scripts (01–06) from QC → integration → annotation → DEG → network → cell–cell communication
- `workflows/`: helper runners for LIANA/CellChat or HPC-style execution
- `environment/`: session info and package requirements
- `examples/`: example figures and minimal demonstration outputs

## How to run (high level)

1. Start with `analysis/01_qc_preprocessing.R` to create per-sample objects.
2. Run `analysis/02_integration.R` to integrate samples/batches.
3. Proceed through scripts `03–06` for downstream analyses.
