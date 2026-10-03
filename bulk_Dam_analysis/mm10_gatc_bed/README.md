### [`mm10_gatc_bed/`](mm10_gatc_bed/)

Reference coordinates for **GATC motifs in the mouse mm10 genome**, used for GATC-aware feature selection and GATC-density normalization.

- [`mm10_GATC.bed`](mm10_gatc_bed/mm10_GATC.bed)
  - BED coordinates of all GATC motifs in mm10.
  - Used by the bin- and peak-based quantification workflows and by GATC-density normalization.

- [`generate_mm10_GATC_bed.py`](mm10_gatc_bed/generate_mm10_GATC_bed.py)
  - Scans a reference FASTA sequence and generates the corresponding GATC BED file.
  - For another genome assembly or species, modify the reference FASTA path and generate the corresponding GATC coordinates before analysis.
