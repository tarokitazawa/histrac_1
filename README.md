# histrac_1

Code for the paper:

**“Whole-genome single-cell multimodal history tracing to reveal cell identity transition”**

**Authors:** Yumiko K. Kawamura, Valentina Khalil, and Taro Kitazawa  
**Preprint:** [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.08.12.669973v1)

---

## Overview

This repository contains the analysis code for the **HisTrac-seq** platform described in the above study.

HisTrac-seq uses enzymatic adenine methylation to record past genomic regulatory states and combines this molecular history with present-state single-cell profiling. The method enables retrospective analysis of cell-state dynamics over extended timescales and was used here to investigate cell identity transitions during neuronal differentiation and mouse brain development.

The repository is organized into two main analysis workflows:

### [`bulk_Dam_analysis/`](bulk_Dam_analysis/)

Analysis of **bulk DamID-seq and Dam&Tag** data.

This directory contains workflows for:

- DamID-seq and Dam&Tag read preprocessing and alignment;
- generation of genomic GATC reference coordinates;
- QuasR-based quantification over genomic bins, genes, and ATAC-defined regulatory regions;
- CPM, GATC-density, and FreeDam normalization; and
- generation of BigWig tracks for genome-browser visualization.

Detailed workflow descriptions and software requirements are provided in the [`bulk_Dam_analysis/README.md`](bulk_Dam_analysis/README.md).

---

### [`sc_Dam&Tag_HisTrac/`](sc_Dam%26Tag_HisTrac/)

Analysis of **single-cell Dam&Tag and HisTrac-seq** data.

This directory contains:

- [`nanoscope_implementation/`](sc_Dam%26Tag_HisTrac/nanoscope_implementation/) — implementation and study-specific modifications of the Nanoscope preprocessing workflow for demultiplexing, alignment, cell calling, and peak calling;
- [`single_cell_analysis/`](sc_Dam%26Tag_HisTrac/single_cell_analysis/) — core Seurat/Signac workflow for QC, replicate merging, dimensionality reduction, clustering, gene activity analysis, and comparison of historical Dam and present-state H3K27ac modalities;
- [`single_cell_analysis_demo/`](sc_Dam%26Tag_HisTrac/single_cell_analysis_demo/) — a compact example workflow for reproducing the main single-cell analysis steps;
- [`single_cell_analysis_vitro/`](sc_Dam%26Tag_HisTrac/single_cell_analysis_vitro/) — downstream analyses of the in vitro neuronal differentiation experiments, including Dam-Leo1, Free-Dam, Dam-Biv, and perturbation experiments; and
- [`single_cell_analysis_vivo/`](sc_Dam%26Tag_HisTrac/single_cell_analysis_vivo/) — downstream analyses of Dam-Leo1 history tracing in the developing and adult mouse cortex.

Detailed preprocessing, software requirements, and analysis instructions are provided in the corresponding README files within these directories.

---

## Analysis overview

```text
Bulk DamID / Dam&Tag
        |
        +--> read preprocessing and alignment
        +--> QuasR quantification
        +--> GATC / FreeDam normalization
        `--> genome-browser tracks


Single-cell Dam&Tag / HisTrac-seq
        |
        +--> Nanoscope preprocessing
        +--> Seurat / Signac object preparation
        +--> clustering and gene activity analysis
        +--> historical vs present-state comparison
        |
        +--> in vitro HisTrac analyses
        `--> in vivo cortical HisTrac analyses
```

---

## Demo

A compact example of the HisTrac-seq single-cell analysis workflow is provided in:

[`sc_Dam&Tag_HisTrac/single_cell_analysis_demo/`](sc_Dam%26Tag_HisTrac/single_cell_analysis_demo/)

The demo illustrates the main steps from reorganized Nanoscope output to merged Seurat objects and downstream analysis.

---

## Software requirements

Software dependencies and version information are documented separately for the bulk and single-cell workflows:

- [`bulk_Dam_analysis/README.md`](bulk_Dam_analysis/README.md)
- [`sc_Dam&Tag_HisTrac/README.md`](sc_Dam%26Tag_HisTrac/README.md)

---

## Data

Data generated in this study are available through ArrayExpress:

- **Bulk sequencing:** E-MTAB-15338
- **Single-cell sequencing:** E-MTAB-15341

---

## Citation

Kawamura YK, Khalil V, Kitazawa T (2025).  
**Whole-genome single-cell multimodal history tracing to reveal cell identity transition.**  
bioRxiv. https://doi.org/10.1101/2025.08.12.669973
