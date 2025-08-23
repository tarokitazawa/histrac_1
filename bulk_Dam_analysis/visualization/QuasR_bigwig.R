# Bigwig preparation from BAM files with QuasR
library (QuasR)
library (Rbowtie)
library (BSgenome)
library (Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(Gviz)
library (parallel)
library(GenomicRanges)
library(edgeR)
library(Rhisat2)
library(BSgenome.Mmusculus.UCSC.mm10)

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

# statistics (mapped reads)
getwd()
alignmentStats(proj_Dam)


# 100 bp bin bigwig (bin size can be changed by binsize = )
path <- '/path/to/your/project/bigwig/'
setwd(path)
getwd()
cl <- makeCluster(8)
qExportWig(proj_Dam, binsize = 100L, scaling = 1000000,
           createBigWig = TRUE, clObj = cl
)

