# Bulk DamID / Dam&Tag Analysis

This directory contains the scripts and pipelines used for **bulk DamID-seq** and **bulk Dam&Tag** analysis

---

## Overview

- Preprocessing of bulk DamID-seq and Dam&Tagreads
- Generation of normalized coverage tracks 
- Quantification of m6A signals over genomic features 
- Visualization of results

---

## Requirements

### Software
- cutadapt (v3.5)
- Bowtie2 (v2.5.1)
- samtools (v1.6)
- deepTools (v3.5.5)
- IGV (v2.16.0)

### R packages
- QuasR (>= 1.44.0)

---

## Workflow

### 1a. Preprocessing for DamID-seq
- Trim adapters (NEB Next) using cutadapt. Filter DamID reads to retain DpnI-derived fragments (5′ TC overhang)
- Align to **mm10** reference genome using Bowtie2. Convert and sort BAM files with samtools

### 1b. Preprocessing for Dam&Tag
- Trim adapters (Illumina Nextera) using cutadapt
- Align to **mm10** reference genome using Bowtie2. Convert and sort BAM files with samtools

### 2. Quantification
- Use **QuasR** to count reads in genomic bins or features (e.g., genes, enahancers)

### 3. Normalization
- Generate CPM (counts per million) matrices
- Normalize by GATC density
- Normalize by FreeDam (except for Dam-RNAPII, Dam-Leo1)

### 4. Visualization
- Generate bigwig by QuasR
- Generate FreeDam-normalized bigwig by deepTools
- Visualization with IGV
---
## examples of raw fastq files (paried end)
- DamID-seq fastq files: /path/to/your/project/raw_fastq_filepaths_dam.txt
- Dam&Tag fastq files: /path/to/your/project/raw_fastq_filepaths_damtag.txt
- Read1: /path/to/your/raw_fastq/sample1_R1_001.fastq.gz
- Read2: /path/to/your/raw_fastq/sample1_R2_001.fastq.gz
---
## GATC regions of mm10
- mm10_GATC.bed
---
---

## Data
- Bulk E-MTAB-15336
---

## Citation
Kawamura YK, Khalil V, Kitazawa T (2025).
Whole-genome single-cell multimodal history tracing to reveal cell identity transition.
bioRxiv. https://doi.org/10.1101/2025.08.12.669973
