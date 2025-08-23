# histrac_1
Code for the paper:  
**"Whole-genome single-cell multimodal history tracing to reveal cell identity transition"**  
Authors: Yumiko K. Kawamura, Valentina Khalil, and Taro Kitazawa 
Journal/Preprint: bioRxiv (https://www.biorxiv.org/content/10.1101/2025.08.12.669973v1)

---

## Overview
This repository contains the code for the HisTrac-seq platform introduced in the above preprint deposited in bioRxiv.  
HisTrac-seq enzymatically labels genomic DNA adenine to “bookmark” regulatory states and enables temporal multi-omics profiling of mouse brain across two months, revealing large-scale single-cell identity transitions (“identity jumps”).

---

## Requirements

### Preprocessing
- cutadapt (v3.5)
- Bowtie2 (v2.5.1)
- STAR (v2.7)
- samtools (v1.6)
- MACS2 (v2.2.7.1)
- deepTools (v3.5.5)

### R (>=4.4.0) packages
- QuasR (1.44.0)
- edgeR (4.2.2)
- limma (3.60.0)
- Seurat (5.3.0)
- Signac (1.14.0)
- clusterProfiler (4.12.6)
- ggalluvial (0.12.5)
- motifmatchr + TFBSTools
- JASPAR2024 database
- pheatmap (1.0.13)

### Python
- Included in Cell Ranger ATAC (10x Genomics v2.0.0)
- nanoscope pipeline (from Bartosovic lab)

### Visualization
- IGV (v2.16.0)

## Single-cell analysis with Nanoscope
For single-cell Dam&Tag / HisTrac-seq data, we used the
[Nanoscope pipeline](https://github.com/bartosovic-lab/nanoscope) 
(bartosovic-lab).  
