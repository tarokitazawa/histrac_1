# QuasR quantification of bulk Dam for peaks (from BAMs)

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

# qProject prepration
# working directory
# Dam_filepaths.txt is the bam input file format of QuasR (tab-delineated, 2 columns with headers FileName and SampleName)
path <- '/path/to/your/project/'
setwd(path)
getwd()
cl <- makeCluster(8)
proj_Dam <- qAlign("Dam_filepaths.txt", "BSgenome.Mmusculus.UCSC.mm10", 
                   aligner = "Rbowtie", paired = "fr", splicedAlignment = FALSE, 
                   clObj = cl
)

#Grange of ATAC peaks (called by MACS2).
df <- read.table("/path/to/your/ATAC_peak/ATAC_summits.bed", 
                 header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
colnames(df) <- c("chr", "summit_start", "summit_end", "peakID", "score")
head(df)
# Compute the midpoint of each summit
df$mid <- floor((df$summit_start + df$summit_end) / 2)
# Expand the midpoint to ±500 bp
df$extended_start <- df$mid - 500
df$extended_end   <- df$mid + 500
# we can set it to 1 to avoid negative or 0-based coordinates.
df$extended_start[df$extended_start < 1] <- 1
atac_peak <- GRanges(
  seqnames = df$chr,
  ranges   = IRanges(
    start = df$extended_start,
    end   = df$extended_end
  ),
  peakID   = df$peakID,
  score    = df$score
)
atac_peak

# Select peaks containing GATC motif
mm10_gatc <- import("/path/to/your/project/mm10_GATC.bed", format = "bed")
overlaps <- findOverlaps(atac_peak, mm10_gatc)
overlap_indices <- unique(queryHits(overlaps))
atac_peak_gatc <- atac_peak[overlap_indices]
atac_peak_gatc
length(atac_peak)
length(atac_peak_gatc)

# Remove promoters and identify enhancers
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
columns(txdb)
promoters_mm10 <- promoters(txdb, upstream=1000, downstream=1000)
promoters_mm10
overlaps <- findOverlaps(atac_peak_gatc, promoters_mm10)
overlapping_indices <- unique(queryHits(overlaps))
non_overlapping_indices <- setdiff(seq_along(atac_peak_gatc), overlapping_indices)
atac_peak_gatc_enhancer <- atac_peak_gatc[non_overlapping_indices]
atac_peak_gatc_enhancer
length(atac_peak_gatc)
length(atac_peak_gatc_enhancer) 

Dam_count <- qCount(proj_Dam, atac_peak_gatc_enhancer, shift = "halfInsert", clObj = cl)
saveRDS(Dam_count,"/path/to/your/project/raw_counts/Dam_enhancer_count.rds")
