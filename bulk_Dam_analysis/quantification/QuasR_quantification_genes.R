# QuasR quantification of bulk Dam for genes (from BAMs)

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

# Gene Grange generation
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_mm10 <- genes(txdb)
genes_mm10
starts <- start(genes_mm10)
ends <- end(genes_mm10)
strands <- strand(genes_mm10)
extend_amount <- 3000
starts <- starts - extend_amount
ends <- ends + extend_amount
start(genes_mm10) <- starts
end(genes_mm10) <- ends
genes_mm10

Dam_count <- qCount(proj_Dam, genes_mm10, shift = "halfInsert", clObj = cl)
saveRDS(Dam_count,"/path/to/your/project/raw_counts/Dam_gene_count.rds")
