# Dam normalization by GATC frequency

library (edgeR)
library (rtracklayer)
library (GenomicRanges)

# working directory
path <- '/path/to/your/project/'
setwd(path)
getwd()

# Load the raw count table produced by the corresponding QuasR quantification script:
#   Bins:  Dam_10kb_bin_count.rds from QuasR_quantification_bin.R
#   Genes: Dam_gene_count.rds from QuasR_quantification_genes.R
#   Peaks: Dam_enhancer_count.rds from QuasR_quantification_peaks.R
#
# Use the exact genomic feature set and order used for quantification:
#   Bins:  bins_gatc
#   Genes: genes_mm10
#   Peaks: atac_peak_gatc_enhancer
#
# The corresponding feature object must also be loaded or reconstructed;
# it is not stored in the count table itself.
# The example below uses 10-kb bins containing at least one GATC motif.

Dam_count <- readRDS("/path/to/raw_count_table/Dam_10kb_bin_count.rds")
dim(Dam_count)
head(Dam_count)

#cpm (Dam_count[,1] is gene length)
Dam_cpm <- cpm(Dam_count[,-1])
head(Dam_cpm)
dim(Dam_cpm)

# Select genomic features containing at least one GATC motif.
#
# Feature objects are defined in the corresponding quantification scripts:
#   Bins:  bins in QuasR_quantification_bin.R
#   Genes: genes_mm10 in QuasR_quantification_genes.R
#   Peaks: atac_peak in QuasR_quantification_peaks.R
#
# For gene- or peak-based analysis, use the corresponding feature object
# and retain only features containing at least one GATC motif.
# The retained features must match the count-table rows in both identity and order.
#
# Note for genes:
#   The current gene quantification script counts all genes_mm10.
#   If GATC-free genes are removed here, remove the same rows from Dam_cpm.
#
# Note for peaks:
#   The current peak quantification script also excludes promoter-overlapping peaks.
#   Therefore, use atac_peak_gatc_enhancer to match Dam_enhancer_count.rds,
#   rather than all GATC-containing peaks.
#
# The example below uses bins defined in QuasR_quantification_bin.R.
# This object must already be available in the R session.

mm10_gatc <- import("/path/to/your/project/mm10_GATC.bed", format = "bed")
overlaps <- findOverlaps(bins, mm10_gatc)
overlap_indices <- unique(queryHits(overlaps))
bins_gatc <- bins[overlap_indices]

bins_gatc
length(bins)
length(bins_gatc)

# Check that the number of count-table rows and reference features agree.
# This checks dimensions only; feature identity and order must also match.
stopifnot(nrow(Dam_cpm) == length(bins_gatc))

# DamID normalization by GATC counts
overlaps <- findOverlaps(bins_gatc, mm10_gatc)
gatc_counts <- table(queryHits(overlaps))
head(gatc_counts)
summary(gatc_counts)
sum(gatc_counts == 0)
length(gatc_counts)
nrow(Dam_cpm)
length(bins_gatc)
# Divide CPM by (0.256 × GATC count) for normalization.
# Rationale: GATC occurs ~1/256 bp, i.e. ~3.9 sites per 1 kb.
# This scaling approximates RPKM by adjusting for GATC density.
Dam_gatcNorm <- sweep(Dam_cpm, 1, 0.256*gatc_counts, "/")
head(Dam_gatcNorm)
dim(Dam_gatcNorm)

saveRDS(Dam_gatcNorm,"/path/to/your/project/gatcNrom_counts/Dam_gatcNorm.rds")

