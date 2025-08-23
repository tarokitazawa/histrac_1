# Dam normalization by GATC frequency

library (QuasR)
library (Rbowtie)
library (BSgenome)
library (Rsamtools)
library (rtracklayer)
library (GenomicFeatures)
library (parallel)
library (GenomicRanges)
library (edgeR)
library (Rhisat2)
library (BSgenome.Mmusculus.UCSC.mm10)
library (TxDb.Mmusculus.UCSC.mm10.knownGene)

# working directory
path <- '/path/to/your/project/'
setwd(path)
getwd()

# Count table produced by QuasR (e.g., QuasR_quantification_bin.R)
Dam_count <- readRDS("/path/to/raw_count_table/Dam_10kb_bin_count.rds")
dim(Dam_count)
head(Dam_count)

#cpm (Dam_count[,1] is gene length)
Dam_cpm <- cpm(Dam_count[,-1])
head(Dam_cpm)
dim(Dam_cpm)

# Bins containing GATC motif
mm10_gatc <- import("/path/to/your/project/mm10_GATC.bed", format = "bed")
overlaps <- findOverlaps(bins, mm10_gatc)
overlap_indices <- unique(queryHits(overlaps))
bins_gatc <- bins[overlap_indices]
bins_gatc
length(bins)
length(bins_gatc)
all(rownames(Dam_cpm) == rownames(bins_gatc)) # Confirmation that we are using the identical reference

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

