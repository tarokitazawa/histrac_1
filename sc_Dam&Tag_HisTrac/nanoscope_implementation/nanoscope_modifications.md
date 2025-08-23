# Nanoscope Pipeline Implementation and Modifications

This document describes how we implemented the **Nanoscope pipeline** (Bartosovic lab) for the analysis of **scDam&Tag** and **scHisTrac-seq** data, and highlights the modifications applied in our implementation.

---

## Workflow Summary

We used the Nanoscope pipeline, originally designed for the analysis of **nanoCUT&Tag (nanoCT)**, to preprocess raw single-cell sequencing data.

1. **Demultiplexing**  

2. **Cell Ranger ATAC (v2.0.0)**  
   - The demultiplexed FASTQs were processed with **Cell Ranger ATAC** for:
     - Alignment to the **mm10** reference genome  
     - Basic quality control (QC)  
     - Preparation of **fragments files** (`fragments.tsv.gz`)  

3. **Peak calling with MACS2**  
   - Broad-peak calling was performed on aggregated pseudo-bulk fragments for each modality using **MACS2**.

4. **Metadata collation**  
   - Finally, Nanoscope collated all outputs into a unified metadata table, annotating each cell with:
     - Total fragments  
     - Unique fragment counts per cell  
     - FRiP (Fraction of Reads in Peaks)  
   - This metadata was used for downstream single-cell analysis (Seurat/Signac).

---

## Modification of Demultiplexing Script

In our implementation, we modified the **antibody-barcode “follow” pattern** in the Nanoscope demultiplexing script (`debarcode.py`):

- **Original default**  
  ```python
  default = "GCGTGGAGACGCTGCCGACGA"
- **Replaced with**  
  ```python
  default = "GACGCTGCCGACGA"

### Rationale
Our sample/modality-barcoded **ME-A adapters** were designed based on **NTT-seq**, which deviates from the nanoCUT&Tag design.  
As a result, the hard-coded 20-mer sequence in the original script no longer matched the region immediately downstream of the antibody barcode.  
Replacing it with the correct 14-mer sequence ensured proper demultiplexing of our libraries.  
(See Supplementary Table 3 of the manuscript for details.)
