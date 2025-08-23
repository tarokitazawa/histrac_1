# scDam&Tag / scHisTrac Analysis

This directory contains the single-cell analysis workflows for **scDam&Tag** and **scHisTrac-seq**, as described in our paper:

---

## Overview

The workflow for scDam&Tag and scHisTrac-seq follows the Methods section of the paper:

1. **Demultiplexing and alignment with Nanoscope**  
   - We use [Nanoscope](https://github.com/bartosovic-lab/nanoscope) (Bartosovic lab) for preprocessing.  
   - Input: raw FASTQ files (from 10x Genomics Chromium).  
   - Nanoscope performs:
     - **Sample/modality demultiplexing** (separates Dam and H3K27ac reads by sample barcodes and antibody barcodes).  
     - **FASTQ splitting** into per-modality, per-sample files.  
     - **Cell Ranger ATAC** (10x Genomics) for alignment, cell calling, and QC.  
     - **Peak calling with MACS2** on aggregated pseudo-bulk data.  
   - Environment setup and detailed instructions: see the [Nanoscope GitHub page](https://github.com/bartosovic-lab/nanoscope).

2. **Downstream analysis with Seurat / Signac**  
   - For single-cell analysis, we adapt the **Nanoscope “Analysis using peaks” workflow** (https://fansalon.github.io/vignette_single-cell-nanoCT.html), implemented with:
     - [Seurat](https://satijalab.org/seurat/) (v5.3.0)
     - [Signac](https://stuartlab.org/signac/) (v1.14.0)
   - We provide custom R scripts in this repository.  
   - These include:
     - Cell filtering and QC
     - Clustering and UMAP embedding
     - Gene activity scoring
     - Integration of Dam (historical) vs H3K27ac (present) modalities
     - Visualization of identity jumps

---

## Notes

- **Nanoscope environment**: please follow the official instructions on the [Bartosovic lab GitHub page](https://github.com/bartosovic-lab/nanoscope).  
- **Modifications**: we describe pipeline modifications (barcode pattern changes, etc.) in this README or in `configs/`.  

---
## Data
- single-cell E-MTAB-15341
---

## Citation
Kawamura YK, Khalil V, Kitazawa T (2025).
Whole-genome single-cell multimodal history tracing to reveal cell identity transition.
bioRxiv. https://doi.org/10.1101/2025.08.12.669973
