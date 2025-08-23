# Dam-POI normalization by FreeDam

############################################
# R function to create ratio columns
# from specified pairs of samples
# with custom column names
############################################

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

# GATC-Normalized count table produced by QuasR (e.g., Dam_normalization.R)
Dam_gatcNorm <- readRDS("/path/to/your/normalized_count/Dam_gatcNorm.rds")
dim(Dam_gatcNorm)
head(Dam_gatcNorm)

create_ratio_df <- function(data, sample_info, pseudocount) {
  ratio_df <- data.frame(row.names = rownames(data))
  
  # data: a data frame (e.g. Dam_Chrom_norm)
  # sample_info: a list of character vectors, each of length 3
  #   e.g. list(c("sampleA1","sampleA2","newSampleNameA"), 
  #             c("sampleB1","sampleB2","newSampleNameB"))
  #   where sampleA1, sampleA2 are column names in 'data',
  #         newSampleNameA is the new column name
  #
  # pseudocount: a small number to avoid division by zero;

  for (triple in sample_info) {
    sampleA <- triple[1]
    sampleB <- triple[2]
    new_label <- triple[3]
    
    ratio_df[[new_label]] <- (data[[sampleA]] + pseudocount) / (data[[sampleB]] + pseudocount)
  }
  return(ratio_df)
}

# Example of Taf3 normalization
sample_info <- list(
  c("Taf3_Day2_rep1", "FreeDam_Day2_rep1", "Taf3_Day2_FreeDamNorm_rep1"),
  c("Taf3_Day2_rep2", "FreeDam_Day2_rep2", "Taf3_Day2_FreeDamNorm_rep2"),
  c("Taf3_Day7_rep1", "FreeDam_Day7_rep1", "Taf3_Day7_FreeDamNorm_rep1"),
  c("Taf3_Day7_rep2", "FreeDam_Day7_rep2", "Taf3_Day7_FreeDamNorm_rep2"),
  c("Taf3_Day7_R_rep1", "FreeDam_Day7_R_rep1", "Taf3_Day7_R_FreeDamNorm_rep1"),
  c("Taf3_Day7_R_rep2", "FreeDam_Day7_R_rep2", "Taf3_Day7_R_FreeDamNorm_rep2"),
)

Dam_gatcFreeNorm <- create_ratio_df(
  data = Dam_gatcNorm, 
  sample_info = sample_info,
  pseudocount = 1
)

head(Dam_gatcFreeNorm)
dim(Dam_gatcFreeNorm)
colnames(Dam_gatcFreeNorm)

saveRDS(Dam_gatcFreeNorm,"/path/to/your/project/gatcFreeNorm_counts/Dam_gatcFreeNorm.rds")

# We carry out log2 normalization for the following analysis
# We normally do not normalize RNAPII and Leo1 by FreeDam