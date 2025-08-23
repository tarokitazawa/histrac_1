# Conversion of Nanoscope Output for Replicate Merging

For downstream Seurat/Signac analysis, we followed vignette of nanoscope_base (https://fansalon.github.io/vignette_single-cell-nanoCT.html)
Here, we needed to **merge biological replicates** (rep1, rep2).  
However, the default **Nanoscope output structure** nests each modality directly under the experiment folder, which is not optimal for merging for the above pipeline.

---

## Original Nanoscope Output (per experiment)

```
Day7R_Rep1/
├── Day7R.Dam.rep1_CCTATCCT/
│   ├── barcode_metrics/
│   │   ├── all_barcodes.txt
│   │   └── peaks_barcodes.txt
│   │
│   ├── cellranger/
│   │   └── outs/
│   │       ├── fragments.tsv.gz
│   │       ├── fragments.tsv.gz.tbi
│   │       └── peaks.bed
│   │
│   ├── cell_picking/
│   │   └── metadata.csv
│   │
│   └── peaks/
│       └── macs_broad/
│           └── Day7R.Dam.rep1_peaks.broadPeak
│
└── Day7R.K27ac.rep1_CGTCTAAT/
    ├── barcode_metrics/
    │   ├── all_barcodes.txt
    │   └── peaks_barcodes.txt
    │
    ├── cellranger/
    │   └── outs/
    │       ├── fragments.tsv.gz
    │       ├── fragments.tsv.gz.tbi
    │       └── peaks.bed
    │
    ├── cell_picking/
    │   └── metadata.csv
    │
    └── peaks/
        └── macs_broad/
            └── Day7R.K27ac.rep1_peaks.broadPeak

Day7R_Rep2/
├── Day7R.Dam.rep2_CAGGACGT/
│   ├── barcode_metrics/
│   │   ├── all_barcodes.txt
│   │   └── peaks_barcodes.txt
│   │
│   ├── cellranger/
│   │   └── outs/
│   │       ├── fragments.tsv.gz
│   │       ├── fragments.tsv.gz.tbi
│   │       └── peaks.bed
│   │
│   ├── cell_picking/
│   │   └── metadata.csv
│   │
│   └── peaks/
│       └── macs_broad/
│           └── Day7R.Dam.rep2_peaks.broadPeak
│
└── Day7R.K27ac.rep2_TCGACTAG/
    ├── barcode_metrics/
    │   ├── all_barcodes.txt
    │   └── peaks_barcodes.txt
    │
    ├── cellranger/
    │   └── outs/
    │       ├── fragments.tsv.gz
    │       ├── fragments.tsv.gz.tbi
    │       └── peaks.bed
    │
    ├── cell_picking/
    │   └── metadata.csv
    │
    └── peaks/
        └── macs_broad/
            └── Day7R.K27ac.rep2_peaks.broadPeak
```
Each modality folder name encodes the sample ID, replicate ID, and barcode sequence.

---

## Restructured Format for Merged Analysis

To simplify merging across replicates, we reorganized outputs into the following structure:

```
Day7R_MERGED/
├── rep1/
│   ├── Day7R.Dam_bb/
│   │   ├── barcode_metrics/
│   │   │   ├── all_barcodes.txt
│   │   │   └── peaks_barcodes.txt
│   │   │
│   │   ├── cellranger/
│   │   │   └── outs/
│   │   │       ├── fragments.tsv.gz
│   │   │       ├── fragments.tsv.gz.tbi
│   │   │       └── peaks.bed
│   │   │
│   │   ├── cell_picking/
│   │   │   └── metadata.csv
│   │   │
│   │   └── peaks/
│   │       └── macs_broad/
│   │           └── Day7R.Dam_peaks.broadPeak
│   │
│   └── Day7R.K27ac_aa/
│       ├── barcode_metrics/
│       │   ├── all_barcodes.txt
│       │   └── peaks_barcodes.txt
│       │
│       ├── cellranger/
│       │   └── outs/
│       │       ├── fragments.tsv.gz
│       │       ├── fragments.tsv.gz.tbi
│       │       └── peaks.bed
│       │
│       ├── cell_picking/
│       │   └── metadata.csv
│       │
│       └── peaks/
│           └── macs_broad/
│               └── Day7R.K27ac_peaks.broadPeak
│
└── rep2/
    ├── Day7R.Dam_bb/
    │   ├── barcode_metrics/
    │   │   ├── all_barcodes.txt
    │   │   └── peaks_barcodes.txt
    │   │
    │   ├── cellranger/
    │   │   └── outs/
    │   │       ├── fragments.tsv.gz
    │   │       ├── fragments.tsv.gz.tbi
    │   │       └── peaks.bed
    │   │
    │   ├── cell_picking/
    │   │   └── metadata.csv
    │   │
    │   └── peaks/
    │       └── macs_broad/
    │           └── Day7R.Dam_peaks.broadPeak
    │
    └── Day7R.K27ac_aa/
        ├── barcode_metrics/
        │   ├── all_barcodes.txt
        │   └── peaks_barcodes.txt
        │
        ├── cellranger/
        │   └── outs/
        │       ├── fragments.tsv.gz
        │       ├── fragments.tsv.gz.tbi
        │       └── peaks.bed
        │
        ├── cell_picking/
        │   └── metadata.csv
        │
        └── peaks/
            └── macs_broad/
                └── Day7R.K27ac_peaks.broadPeak
```

### Key Changes
- **Replicate-level directories (`rep1`, `rep2`)** were introduced at the top level.  
- **Modality names** were simplified to `Day7R.Dam_bb` and `Day7R.K27ac_aa` (replacing sample/modality barcode sequences).  
- **Peak files** were renamed consistently (`Day7R.Dam_peaks.broadPeak`, `Day7R.K27ac_peaks.broadPeak`).  
- Internal subfolders (`barcode_metrics`, `cellranger/outs`, `cell_picking`, `peaks/macs_broad`) were preserved.  

---

### Rationale

This restructuring provides:
- Consistent modality names across replicates  
- Easy programmatic iteration over `rep1`, `rep2` for merging  
- Compatibility with downstream Seurat/Signac pipelines (expecting identical directory structure across replicates)  

---

### Notes
- This conversion is a **one-time preprocessing step** to align Nanoscope output with our downstream pipeline.  


---

## Alternative Layout: Merging Across Timepoints

In some analyses (e.g. **medoid calculation** and direct cross-timepoint comparison),  
we merge the data accross **timepoint (Day2, Day7, Day7R)** and **replicates (rep1, rep1)**.  
This allows us to easily  compare across developmental stages.

Example:

```
modality_merged/
|
+---Day2rep1
|   +---Dam_bb
|   |   +---barcode_metrics
|   |   |       all_barcodes.txt
|   |   |       peaks_barcodes.txt
|   |   |
|   |   +---cellranger
|   |   |   \---outs
|   |   |           fragments.tsv.gz
|   |   |           fragments.tsv.gz.tbi
|   |   |           peaks.bed
|   |   |
|   |   +---cell_picking
|   |   |       metadata.csv
|   |   |
|   |   \---peaks
|   |       \---macs_broad
|   |               Dam_peaks.broadPeak
|   |
|   \---K27ac_aa
|       +---barcode_metrics
|       |       all_barcodes.txt
|       |       peaks_barcodes.txt
|       |
|       +---cellranger
|       |   \---outs
|       |           fragments.tsv.gz
|       |           fragments.tsv.gz.tbi
|       |           peaks.bed
|       |
|       +---cell_picking
|       |       metadata.csv
|       |
|       \---peaks
|           \---macs_broad
|                   K27ac_peaks.broadPeak
|
+---Day2rep2
|   +---Dam_bb
|   |   +---barcode_metrics
|   |   |       all_barcodes.txt
|   |   |       peaks_barcodes.txt
|   |   |
|   |   +---cellranger
|   |   |   \---outs
|   |   |           fragments.tsv.gz
|   |   |           fragments.tsv.gz.tbi
|   |   |           peaks.bed
|   |   |
|   |   +---cell_picking
|   |   |       metadata.csv
|   |   |
|   |   \---peaks
|   |       \---macs_broad
|   |               Dam_peaks.broadPeak
|   |
|   \---K27ac_aa
|       +---barcode_metrics
|       |       all_barcodes.txt
|       |       peaks_barcodes.txt
|       |
|       +---cellranger
|       |   \---outs
|       |           fragments.tsv.gz
|       |           fragments.tsv.gz.tbi
|       |           peaks.bed
|       |
|       +---cell_picking
|       |       metadata.csv
|       |
|       \---peaks
|           \---macs_broad
|                   K27ac_peaks.broadPeak
|
+---Day7rep1
|   +---Dam_bb
|   |   +---barcode_metrics
|   |   |       all_barcodes.txt
|   |   |       peaks_barcodes.txt
|   |   |
|   |   +---cellranger
|   |   |   \---outs
|   |   |           fragments.tsv.gz
|   |   |           fragments.tsv.gz.tbi
|   |   |           peaks.bed
|   |   |
|   |   +---cell_picking
|   |   |       metadata.csv
|   |   |
|   |   \---peaks
|   |       \---macs_broad
|   |               Dam_peaks.broadPeak
|   |
|   \---K27ac_aa
|       +---barcode_metrics
|       |       all_barcodes.txt
|       |       peaks_barcodes.txt
|       |
|       +---cellranger
|       |   \---outs
|       |           fragments.tsv.gz
|       |           fragments.tsv.gz.tbi
|       |           peaks.bed
|       |
|       +---cell_picking
|       |       metadata.csv
|       |
|       \---peaks
|           \---macs_broad
|                   K27ac_peaks.broadPeak
|
+---Day7rep2
|   +---Dam_bb
|   |   +---barcode_metrics
|   |   |       all_barcodes.txt
|   |   |       peaks_barcodes.txt
|   |   |
|   |   +---cellranger
|   |   |   \---outs
|   |   |           fragments.tsv.gz
|   |   |           fragments.tsv.gz.tbi
|   |   |           peaks.bed
|   |   |
|   |   +---cell_picking
|   |   |       metadata.csv
|   |   |
|   |   \---peaks
|   |       \---macs_broad
|   |               Dam_peaks.broadPeak
|   |
|   \---K27ac_aa
|       +---barcode_metrics
|       |       all_barcodes.txt
|       |       peaks_barcodes.txt
|       |
|       +---cellranger
|       |   \---outs
|       |           fragments.tsv.gz
|       |           fragments.tsv.gz.tbi
|       |           peaks.bed
|       |
|       +---cell_picking
|       |       metadata.csv
|       |
|       \---peaks
|           \---macs_broad
|                   K27ac_peaks.broadPeak
|
+---Day7Rrep1
|   +---Dam_bb
|   |   +---barcode_metrics
|   |   |       all_barcodes.txt
|   |   |       peaks_barcodes.txt
|   |   |
|   |   +---cellranger
|   |   |   \---outs
|   |   |           fragments.tsv.gz
|   |   |           fragments.tsv.gz.tbi
|   |   |           peaks.bed
|   |   |
|   |   +---cell_picking
|   |   |       metadata.csv
|   |   |
|   |   \---peaks
|   |       \---macs_broad
|   |               Dam_peaks.broadPeak
|   |
|   \---K27ac_aa
|       +---barcode_metrics
|       |       all_barcodes.txt
|       |       peaks_barcodes.txt
|       |
|       +---cellranger
|       |   \---outs
|       |           fragments.tsv.gz
|       |           fragments.tsv.gz.tbi
|       |           peaks.bed
|       |
|       +---cell_picking
|       |       metadata.csv
|       |
|       \---peaks
|           \---macs_broad
|                   K27ac_peaks.broadPeak
|
\---Day7Rrep2
    +---Dam_bb
    |   +---barcode_metrics
    |   |       all_barcodes.txt
    |   |       peaks_barcodes.txt
    |   |
    |   +---cellranger
    |   |   \---outs
    |   |           fragments.tsv.gz
    |   |           fragments.tsv.gz.tbi
    |   |           peaks.bed
    |   |
    |   +---cell_picking
    |   |       metadata.csv
    |   |
    |   \---peaks
    |       \---macs_broad
    |               Dam_peaks.broadPeak
    |
    \---K27ac_aa
        +---barcode_metrics
        |       all_barcodes.txt
        |       peaks_barcodes.txt
        |
        +---cellranger
        |   \---outs
        |           fragments.tsv.gz
        |           fragments.tsv.gz.tbi
        |           peaks.bed
        |
        +---cell_picking
        |       metadata.csv
        |
        \---peaks
            \---macs_broad
                    K27ac_peaks.broadPeak
```

### Why this layout?
- Fixes **modality names** (`Dam_bb`, `K27ac_aa`) across all timepoints  
- Facilitates **cross-timepoint merging** for the same modality  
- Ensures consistency for downstream QC, UMAP, and Seurat/Signac integration pipelines  

---

