# Custom helper functions

`functions_scCT2.R` contains helper functions adapted from the Nanoscope `functions_scCT.R` script for the scDam&Tag/scHisTrac-seq analyses used in this study.

The majority of the original helper functions are retained without substantial modification. The main changes are:

- **Modified QC visualization (`plotCounts`)**
  - The original function is retained as `plotCounts_original`.
  - The modified `plotCounts` function uses a half-violin and dot-based representation of single-cell QC distributions.

- **Extended cross-modality UMAP visualization (`plotConnectModal`)**
  - Support was added for the UMAP coordinate naming convention used in our Seurat objects.
  - Additional versions (`plotConnectModal_2` to `plotConnectModal_5`) provide:
    - adjustable point and connecting-line display parameters;
    - compatibility with both `UMAP_1`/`UMAP_2` and `umap_1`/`umap_2` coordinate names;
    - custom cluster colors;
    - explicit control of the drawing order of connecting lines; and
    - optional randomization of point drawing order for visualization.

Other helper functions, including `listToMeta`, `DepthCorMulMod`, `getUpsetPeaks`, `plotPassed`, `plotPassedCells`, and `commonCellHistonMarks`, are retained from the original Nanoscope helper script.

These functions are included directly in this repository so that the analyses can be reproduced without modifying the local Nanoscope installation.
