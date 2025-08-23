# QuasR quantification of bulk Dam for bins (from BAMs)

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

#Bin Grange generation (10Kb at least one GATC)
genome <- BSgenome.Mmusculus.UCSC.mm10
main_chroms <- paste0("chr", c(1:19, "X", "Y"))
chrom_lengths <- seqlengths(genome)[main_chroms]
bin_size <- 10000  # this is for 10kb bin
bins <- tileGenome(seqlengths = chrom_lengths,
                   tilewidth = bin_size,
                   cut.last.tile.in.chrom = TRUE)
bins
length(bins)
head(bins)

# Select bins containing GATC motif
mm10_gatc <- import("/path/to/your/project/mm10_GATC.bed", format = "bed")
overlaps <- findOverlaps(bins, mm10_gatc)
overlap_indices <- unique(queryHits(overlaps))
bins_gatc <- bins[overlap_indices]
bins_gatc
length(bins)
length(bins_gatc)

Dam_count <- qCount(proj_Dam, bins_gatc, shift = "halfInsert", clObj = cl)
saveRDS(Dam_count,"/path/to/your/project/raw_counts/Dam_10kb_bin_count.rds")
