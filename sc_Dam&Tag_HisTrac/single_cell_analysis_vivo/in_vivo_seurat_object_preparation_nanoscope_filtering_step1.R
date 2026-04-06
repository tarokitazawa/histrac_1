# In vivo-specific Seurat object preparation from nanoscope outputs
# Cell filtering and QC
# Example: DamLeo1-K27ac in vivo dataset containing
#   - P8 snapshot samples
#   - Week 8 snapshot samples
#   - Week 8 retrospective / HisTrac samples
#
# This script summarizes the in vivo-specific preprocessing points that differ
# from the default two-replicate HisTrac example:
#   1) a larger number of samples (12 total in this example), and
#   2) cell filtering based on nanoscope QC (passedMB), rather than 10x Cell Ranger
#      cell-calling output.
#
# The script proceeds up to sample-wise modality matching and saving the
# separated QC-passed objects.
#
# Modified from the nanoscope workflow:
# https://fansalon.github.io/vignette_single-cell-nanoCT.html

# Required packages for this script itself.
# Additional packages may be needed internally by functions_scCT2.R,
# depending on the nanoscope helper functions used.
library(Signac)
library(Seurat)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(ggpubr)

# directory where the nanoscope repo was cloned
repodir <- "/path/to/your/nanoscope/"
source(paste0(repodir, "scripts/functions_scCT2.R"))

setwd("/path/to/your/in_vivo_project/")

# -----------------------------------------------------------------------------
# Basic settings
# -----------------------------------------------------------------------------
genome <- "mm10"
assay <- "peaks"
min_cell <- 1
min_feat <- 1

# In vivo dataset with 12 samples:
#   4 x P8 snapshot
#   4 x Week 8 snapshot
#   4 x Week 8 retrospective / HisTrac
samples <- c(
  "P8rep1", "P8rep2", "P8rep3", "P8rep4",
  "W8rep1", "W8rep2", "W8rep3", "W8rep4",
  "W8Rrep1", "W8Rrep2", "W8Rrep3", "W8Rrep4"
)

# Assumes the same modality structure is present in all samples
modalities <- dir(samples[1])
modalities

# -----------------------------------------------------------------------------
# Read peak coordinates for each modality across all samples
# -----------------------------------------------------------------------------
input.ls <- list()
for (smpl in samples) {
  cat("Loading peaks for sample", smpl, "\n")
  mod_dirs <- dir(smpl, full.names = FALSE)
  for (mod in mod_dirs) {
    cat("\tReading:", mod, "\n")
    mod2 <- sub("_[^_]+$", "", mod)
    file_path <- file.path(smpl, mod, "peaks", "macs_broad", paste0(mod2, "_peaks.broadPeak"))
    if (!file.exists(file_path)) {
      cat("\tFile not found:", file_path, "\n")
      next
    }
    df <- read.table(file_path)[, 1:3]
    colnames(df) <- c("chr", "start", "end")
    input.ls[[paste0(mod, "_", smpl)]] <- df
  }
}

# Convert to GRanges
input.ls <- lapply(input.ls, makeGRangesFromDataFrame)
input.ls

# -----------------------------------------------------------------------------
# Merge peak coordinates across all samples for each modality
# -----------------------------------------------------------------------------
# Unlike the default two-replicate example, the in vivo dataset contains many
# samples spanning multiple stages and recording conditions. Therefore, peaks are
# merged across all 12 samples for each modality.
combined.peaks.ls <- list()
for (mod in modalities) {
  peak_list_mod <- lapply(samples, function(s) input.ls[[paste0(mod, "_", s)]])
  peak_list_mod <- peak_list_mod[!vapply(peak_list_mod, is.null, logical(1))]
  combined.peaks.ls[[mod]] <- reduce(do.call(c, peak_list_mod))
}
combined.peaks.ls

# Optional overlap summary across samples
upset_plot <- getUpsetPeaks(
  modalities = modalities,
  samples = samples,
  combined_peaks_ls = combined.peaks.ls,
  input_ls = input.ls
)
upset_plot

# -----------------------------------------------------------------------------
# Metadata loading
# -----------------------------------------------------------------------------
metadata.ls <- list()
for (smpl in samples) {
  cat("Loading metadata for sample", smpl, "\n")
  for (mod in modalities) {
    cat("\t", mod, "\n")
    metadata.ls[[paste0(mod, "_", smpl)]] <- read.csv(
      paste0(smpl, "/", mod, "/cell_picking/metadata.csv"),
      stringsAsFactors = FALSE
    )
    rownames(metadata.ls[[paste0(mod, "_", smpl)]]) <- metadata.ls[[paste0(mod, "_", smpl)]]$barcode
  }
}

head(metadata.ls)
summary(metadata.ls)

plotPassed(metadata.ls, xaxis_text = 9, angle_x = 60)
plotPassedCells(metadata.ls, samples, modalities)

# -----------------------------------------------------------------------------
# Cell filtering: in vivo-specific choice
# -----------------------------------------------------------------------------
# Cell filtering can also be performed using 10x Cell Ranger-based cell calling,
# as used in some in vitro analyses.
# In the present in vivo dataset, however, we used nanoscope-based QC filtering,
# retaining cells with passedMB == TRUE.
# In the nanoscope plots, passedMB corresponds to the nanoscope QC decision,
# and peak_ratio_MB corresponds to FRiP.
metadata.ls <- lapply(metadata.ls, function(x) x[x$passedMB, ])

head(metadata.ls)
summary(metadata.ls)

plotPassed(metadata.ls, xaxis_text = 9, angle_x = 60)
plotPassedCells(metadata.ls, samples, modalities)

# -----------------------------------------------------------------------------
# Fragments
# -----------------------------------------------------------------------------
fragment.ls <- list()
for (smpl in samples) {
  cat("Loading fragments for sample", smpl, "\n")
  for (mod in modalities) {
    cat("\t", mod, "\n")
    fragment.ls[[paste0(mod, "_", smpl)]] <- CreateFragmentObject(
      path = paste0(smpl, "/", mod, "/cellranger/outs/fragments.tsv.gz"),
      cells = metadata.ls[[paste0(mod, "_", smpl)]]$barcode
    )
  }
}

fragment.ls

# -----------------------------------------------------------------------------
# Quantify peaks for each experiment
# -----------------------------------------------------------------------------
counts.ls <- list()
for (experim in names(fragment.ls)) {
  cat("Analysing experiment:", experim, "\n")

  modal <- paste0(
    str_split_fixed(experim, "_", 4)[, 1], "_",
    str_split_fixed(experim, "_", 4)[, 2]
  )
  cat("\tPeaks from modality:", modal, "\n")

  counts.ls[[experim]] <- FeatureMatrix(
    fragments = fragment.ls[[experim]],
    features = combined.peaks.ls[[modal]],
    cells = metadata.ls[[experim]]$barcode,
    process_n = 20000
  )
}

summary(counts.ls)

# -----------------------------------------------------------------------------
# Create one Seurat object per experiment
# -----------------------------------------------------------------------------
obj.ls <- list()
for (experim in names(counts.ls)) {
  smpl <- str_split_fixed(experim, "_", 3)[, 3]
  modality <- str_split_fixed(experim, "_", 2)[, 1]

  chrom.assay <- CreateChromatinAssay(
    counts = counts.ls[[experim]],
    fragments = fragment.ls[[experim]],
    genome = genome,
    min.cells = min_cell,
    min.features = min_feat
  )

  obj.ls[[experim]] <- CreateSeuratObject(
    counts = chrom.assay,
    assay = assay,
    meta.data = metadata.ls[[experim]],
    project = smpl
  )

  obj.ls[[experim]]$dataset <- experim
  obj.ls[[experim]]$modality <- modality
  obj.ls[[experim]]$sample <- smpl
}

obj.ls

# -----------------------------------------------------------------------------
# Optional additional QC based on logUMI and FRiP-like metrics
# -----------------------------------------------------------------------------
# In the in vivo workflow, nanoscope-based filtering (passedMB) is the primary
# cell-retention criterion. The additional quantile-based filtering below is optional
# and can be used to remove extreme outliers in logUMI or peak_ratio_MB.
quant_high <- 0.99
quant_low <- 0.01

obj.ls.qc <- list()
for (experiment in names(obj.ls)) {
  logUMI_cutoff_high <- quantile(obj.ls[[experiment]]$logUMI, quant_high)
  logUMI_cutoff_low  <- quantile(obj.ls[[experiment]]$logUMI, quant_low)
  peak_ratio_MB_cutoff_low <- quantile(obj.ls[[experiment]]$peak_ratio_MB, quant_low)

  obj.ls.qc[[experiment]] <- subset(
    obj.ls[[experiment]],
    logUMI > logUMI_cutoff_low &
      logUMI < logUMI_cutoff_high &
      peak_ratio_MB > peak_ratio_MB_cutoff_low
  )

  old_n_cell <- nrow(obj.ls[[experiment]][[]])
  new_n_cell <- nrow(obj.ls.qc[[experiment]][[]])
  discarded <- old_n_cell - new_n_cell
  cat(experiment, "\n")
  cat("\tdiscarded", discarded, "cells (", round(discarded / old_n_cell * 100, 2), "%)\n")
}

obj.ls.qc

plotCounts(obj = obj.ls, quantiles = c(quant_low, quant_high), feature = "logUMI")
plotCounts(obj = obj.ls.qc, quantiles = c(quant_low, quant_high), feature = "logUMI")
plotCounts(obj = obj.ls, quantiles = c(quant_low), feature = "peak_ratio_MB", ylabel = "% UMI in peaks")
plotCounts(obj = obj.ls.qc, quantiles = c(quant_low), feature = "peak_ratio_MB", ylabel = "% UMI in peaks")
rm(obj.ls)

# -----------------------------------------------------------------------------
# Sample-wise modality matching
# -----------------------------------------------------------------------------
# The in vivo dataset contains many samples, so modality matching is performed
# sample by sample. For each sample, we keep only cells that passed QC in both
# modalities of interest (Dam and K27ac).

# Optional per-sample overlap plots
venn_list <- lapply(samples, function(s) {
  commonCellHistonMarks(
    mod1 = obj.ls.qc[[paste0(modalities[1], "_", s)]],
    name_mod1 = strsplit(modalities[1], "_")[[1]][1],
    mod2 = obj.ls.qc[[paste0(modalities[2], "_", s)]],
    name_mod2 = gsub("H3", "", strsplit(modalities[2], "_")[[1]][1]),
    sample = s
  )
})

# Example display (modify ncol as needed)
do.call(ggarrange, c(venn_list, ncol = 4))

# Intersect cells between modalities within each sample
common_cells_by_sample <- setNames(vector("list", length(samples)), samples)
for (s in samples) {
  obj1 <- obj.ls.qc[[paste0(modalities[1], "_", s)]]
  obj2 <- obj.ls.qc[[paste0(modalities[2], "_", s)]]
  common_cells_by_sample[[s]] <- intersect(Cells(obj1), Cells(obj2))
}

# Subset objects to the shared cells within each sample
obj.ls.qc.subset <- list()
for (s in samples) {
  for (mod in modalities) {
    key <- paste0(mod, "_", s)
    obj.ls.qc.subset[[key]] <- subset(
      x = obj.ls.qc[[key]],
      cells = common_cells_by_sample[[s]]
    )
  }
}

# Check the result
obj.ls.qc.subset

plotCounts(obj = obj.ls.qc.subset, quantiles = c(quant_low, quant_high), feature = "logUMI")
plotCounts(obj = obj.ls.qc.subset, quantiles = c(quant_low), feature = "peak_ratio_MB", ylabel = "% UMI in peaks")

# Save the separated, modality-matched, QC-passed objects.
# At this stage, samples are still kept separate. This is useful when later
# analyses require sample-wise control before merging at the stage or condition level.
saveRDS(obj.ls.qc.subset, file = "nanoscope_separated_nanoscopefilter_in_vivo.rds")
