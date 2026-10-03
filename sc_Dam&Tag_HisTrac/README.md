# scDam&Tag / scHisTrac Analysis

This directory contains the single-cell analysis workflows for **scDam&Tag** and **scHisTrac-seq**, as described in our paper.

---

## Overview

The workflow for scDam&Tag and scHisTrac-seq follows the Methods section of the paper:

### 1. Demultiplexing and alignment with Nanoscope

- We use [Nanoscope](https://github.com/bartosovic-lab/nanoscope) (Bartosovic lab) for preprocessing.
- Input: raw FASTQ files from 10x Genomics Chromium.
- Nanoscope performs:
  - **Sample/modality demultiplexing**, separating Dam and H3K27ac reads using sample and antibody barcodes.
  - **FASTQ splitting** into per-modality and per-sample files.
  - **Cell Ranger ATAC** (10x Genomics) for alignment, cell calling, and QC.
  - **Peak calling with MACS2** on aggregated pseudo-bulk data.
- Configuration files and modifications used in this study are provided in [`nanoscope_implementation/`](nanoscope_implementation/).
- For Nanoscope environment setup and general usage, see the [Nanoscope GitHub page](https://github.com/bartosovic-lab/nanoscope).

### 2. Downstream analysis with Seurat / Signac

- For downstream single-cell analysis, we adapt the [Nanoscope “Analysis using peaks” workflow](https://fansalon.github.io/vignette_single-cell-nanoCT.html), implemented with:
  - [Seurat](https://satijalab.org/seurat/) (v5.3.0)
  - [Signac](https://stuartlab.org/signac/) (v1.14.0)
- Core downstream analysis scripts are provided in [`single_cell_analysis/`](single_cell_analysis/).
- These analyses include:
  - Cell filtering and QC
  - Clustering and UMAP embedding
  - Gene activity scoring
  - Comparison of Dam (historical) and H3K27ac (present) modalities
  - Visualization of cell identity transitions
- Custom helper functions adapted from the Nanoscope workflow are provided in [`helper_functions/`](helper_functions/).
- A compact example of the HisTrac-seq analysis workflow is provided in [`single_cell_analysis_demo/`](single_cell_analysis_demo/).

### 3. In vitro scHisTrac-seq analysis

- Scripts for the **in vitro mESC neurodifferentiation experiments** are provided in [`single_cell_analysis_vitro/`](single_cell_analysis_vitro/).
- These include downstream analyses and statistical testing of **Dam-Leo1, Free-Dam, and Dam-Biv** datasets, as well as analyses of perturbation experiments including **retinoic acid (RA)** and **A395-mediated Polycomb perturbation**.
- The analyses cover cell identity transitions, ambiguity/commitment measurements, modality comparisons, gene activity analyses, and statistical evaluation of identity-jump phenotypes.

### 4. In vivo scHisTrac-seq analysis

- Scripts for the **in vivo mouse cortex Dam-Leo1 history-tracing experiments** are provided in [`single_cell_analysis_vivo/`](single_cell_analysis_vivo/).
- These include analysis of postnatal and adult cortical states and retrospective Dam-Leo1 recordings, including cell-type annotation, clustering, developmental identity comparisons, and history-tracing analyses of cell identity transitions.

---

## Notes

- **Nanoscope environment**: please follow the official instructions on the [Bartosovic lab GitHub page](https://github.com/bartosovic-lab/nanoscope).
- **Nanoscope modifications**: study-specific changes to barcode patterns and preprocessing configurations are described in [`nanoscope_implementation/`](nanoscope_implementation/), including [`nanoscope_modifications.md`](nanoscope_implementation/nanoscope_modifications.md).
- **Custom helper functions**: modified helper functions used for downstream analysis are provided in [`helper_functions/functions_scCT2.R`](helper_functions/functions_scCT2.R), with the major changes from the original Nanoscope helper functions summarized in [`helper_functions/README.md`](helper_functions/README.md).

---

## Data

- Single-cell dataset: **E-MTAB-15341**

---

## Typical install time on a "normal" desktop computer

- 1–2 h

---

## Citation

Kawamura YK, Khalil V, Kitazawa T (2025).  
**Whole-genome single-cell multimodal history tracing to reveal cell identity transition.**  
bioRxiv. https://doi.org/10.1101/2025.08.12.669973
