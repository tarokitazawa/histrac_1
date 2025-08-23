# Seurat object preparation from nanoscope outputs
# Cell filtering and QC
# Merging of rep1 and rep2 (example - Dam-Biv-K27ac Day7-Retrospective (=HisTrac) sample)
# Restructuring of nanoscope output is explained in nanoscope_output_restructure.md
# Modified pipeline of nanoscope https://fansalon.github.io/vignette_single-cell-nanoCT.html

library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(stringr)
library(ggplot2)
library(gghalves)
library(ggpubr)
library(EnsDb.Mmusculus.v79)
library(ComplexUpset)
library(regioneR)
library(scales)
library(ggVennDiagram)
library(harmony)
library(biovizBase)

# directory where the nanoscope repo was cloned
repodir <- "/path/to/your/nanoscope/"
source(paste0(repodir,"scripts/functions_scCT2.R"))

setwd("path/to/your/nanoscope/project/")

# Set the used genome version
genome <- "mm10"
# Genome annotations
genome_ann <- EnsDb.Mmusculus.v79
# Assay (peaks or bins)
assay <- "peaks"
# Samples
samples <- c("rep1", "rep2")
# Modalities
modalities <- dir(samples[1]) 
# Minimum number of cells a feature must be detected in to be included in the Seurat obj
min_cell <- 1
# Minium number of features a cell must have to be included in the Seurat obj
min_feat <- 1

# Read-in peak coordinates for each modality iterating across samples and modalities
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

#Grange conversion
input.ls <- lapply(input.ls,function(x){ makeGRangesFromDataFrame(x) })
input.ls

#Merging
combined.peaks.ls <- list()
for (mod in modalities) {
  combined.peaks.ls[[mod]] <- reduce(x = c(input.ls[[paste0(mod,"_",samples[1])]],
                                           input.ls[[paste0(mod,"_",samples[2])]]))
  
}
combined.peaks.ls

# Upset plot works only if we have replicates
upset_plot <- getUpsetPeaks(modalities = modalities, 
                            samples = samples, 
                           combined_peaks_ls = combined.peaks.ls, 
                           input_ls = input.ls)
upset_plot

##### Metadata

# Load metadata for each experiment iterating across the sample list
metadata.ls <- list()
for (smpl in samples) {
  cat("Loading metadata for sample",smpl,"\n")
  for (mod in modalities) {
    cat("\t",mod,"\n")
    
    # Load metadata
    metadata.ls[[paste0(mod,"_",smpl)]] <- read.csv(paste0(smpl,"/",mod,"/cell_picking/metadata.csv"),
                                                    stringsAsFactors = FALSE)
    rownames(metadata.ls[[paste0(mod,"_",smpl)]]) <- metadata.ls[[paste0(mod,"_",smpl)]]$barcode
    
  }
}
head(metadata.ls)
summary(metadata.ls)

plotPassed(metadata.ls,xaxis_text = 9,angle_x=60)
plotPassedCells(metadata.ls,samples,modalities)

# 10X-based cell filtering
metadata.ls <- lapply(metadata.ls, function(x) {
  # Convert if necessary
  x$is__cell_barcode <- as.numeric(x$is__cell_barcode) 
  # or as.logical, depending on your goal
  # Keep only rows == 1
  x[x$is__cell_barcode == 1, ]
})
head(metadata.ls)
summary(metadata.ls)
plotPassed(metadata.ls,xaxis_text = 9,angle_x=60) 
plotPassedCells(metadata.ls,samples,modalities)

##### Fragments
# Create one fragment object per each experiment iterating across the sample list
fragment.ls <- list()
for (smpl in samples) {
  cat("Loading peaks for sample",smpl,"\n")
  for (mod in modalities) {
    cat("\t",mod,"\n")
    
    # Load fragments
    fragment.ls[[paste0(mod,"_",smpl)]] <- CreateFragmentObject(
      path = paste0(smpl,"/",mod,"/cellranger/outs/fragments.tsv.gz"),
      cells = metadata.ls[[paste0(mod,"_",smpl)]]$barcode)
    
  }
}

fragment.ls
head(fragment.ls[[paste0("Day7R.K27ac_aa_",samples[1])]])
head(fragment.ls[[paste0("Day7R.Dam_bb_",samples[2])]])

#Quantifying peaks in each experiment
# As usual, append the object to a list
counts.ls <- list()
# Iterate among the experiment names in order to retrieve fragments and cell barcodes from the same experiment
for (experim in names(fragment.ls)) {
  
  cat("Analysing experiment: ",experim,"\n")
  
  # modality name
  modal <- paste0(str_split_fixed(experim,"_",4)[,1],"_",str_split_fixed(experim,"_",4)[,2])
  cat("\tPeaks from modality: ",modal,"\n")
  
  # Create FeatureMatrix
  counts.ls[[experim]] <- FeatureMatrix(fragments = fragment.ls[[experim]],           # fragments
                                        features = combined.peaks.ls[[modal]],        # modality common set of peaks
                                        cells = metadata.ls[[experim]]$barcode,       # metadata
                                        process_n = 20000) # the higher, the faster, the more memory is used
  
}
summary(counts.ls)
counts.ls[[paste0(modal,"_",samples[1])]][1:10,1:25]

#Create one Seurat object for each experiment
obj.ls <- list()
# Iterate among the experiment names in order to retrieve fragments and cell barcodes from the same experiment
for (experim in names(counts.ls)) {
  
  # Get sample and modality name from the experiment variable
  smpl <- str_split_fixed(experim,"_",3)[,3]
  modality <- str_split_fixed(experim,"_",2)[,1]
  
  # Create the chromatin assay
  chrom.assay <- CreateChromatinAssay(counts = counts.ls[[experim]],
                                      fragments = fragment.ls[[experim]],
                                      genome = genome,
                                      min.cells = min_cell,
                                      min.features = min_feat)
  
  # Create the object
  obj.ls[[experim]] <- CreateSeuratObject(counts = chrom.assay,
                                          assay = assay,
                                          meta.data=metadata.ls[[experim]],
                                          project = smpl)
  
  # Add some metadata of origin
  obj.ls[[experim]]$dataset <- experim
  obj.ls[[experim]]$modality <- modality
  obj.ls[[experim]]$sample <- smpl
  
}
obj.ls
name_1 <- names(obj.ls)[[1]]
name_2 <- names(obj.ls)[[2]]
name_3 <- names(obj.ls)[[3]]
name_4 <- names(obj.ls)[[4]]
name_1
name_2
name_3
name_4

#Additional QC-------------------------------------------------------------------
quant_high <- 0.99
quant_low <- 0.01

obj.ls.qc <- list()
for (experiment in names(obj.ls)) {
  
  # Calculate 1th and 99th quantiles for logUMI
  logUMI_cutoff_high <- quantile(obj.ls[[experiment]]$logUMI, quant_high)
  logUMI_cutoff_low  <- quantile(obj.ls[[experiment]]$logUMI, quant_low)
  # Calculate 1th quantile for logUMI for % fragments in peaks (here indicated by peak_ratio_MB)
  peak_ratio_MB_cutoff_low  <- quantile(obj.ls[[experiment]]$peak_ratio_MB, quant_low)
  
  # Filter out cells outside these thresholds
  obj.ls.qc[[experiment]] <- subset(obj.ls[[experiment]],
                                    logUMI > logUMI_cutoff_low &
                                      logUMI < logUMI_cutoff_high &
                                      peak_ratio_MB > peak_ratio_MB_cutoff_low)
  
  # print number of discarded cells
  old_n_cell <- nrow(obj.ls[[experiment]][[]])
  new_n_cell <- nrow(obj.ls.qc[[experiment]][[]])
  discarded <- old_n_cell - new_n_cell
  cat(experiment,"\n")
  cat("\tdiscarded",discarded,"cells (",round(discarded/old_n_cell*100,2),"%)\n")
}
obj.ls.qc

plotCounts(obj = obj.ls, quantiles = c(quant_low,quant_high), feature = "logUMI")
plotCounts(obj = obj.ls.qc, quantiles = c(quant_low,quant_high), feature = "logUMI")
plotCounts(obj = obj.ls, quantiles = c(quant_low), feature = "peak_ratio_MB",ylabel = "% UMI in peaks")
plotCounts(obj = obj.ls.qc, quantiles = c(quant_low), feature = "peak_ratio_MB",ylabel = "% UMI in peaks")
rm(obj.ls)


#-----------------------------------------------------------------------------------------------
# Plot number of cells for each modality for each sample
venn1 <- commonCellHistonMarks(mod1 = obj.ls.qc[[paste0(modalities[1],"_",samples[1])]], 
                               name_mod1 = strsplit(modalities[1],"_")[[1]][1],
                               mod2 = obj.ls.qc[[paste0(modalities[2],"_",samples[1])]], 
                               name_mod2 = gsub("H3","",strsplit(modalities[2],"_")[[1]][1]),
                               mod3 = obj.ls.qc[[paste0(modalities[3],"_",samples[1])]], 
                               name_mod3 = gsub("H3","",strsplit(modalities[3],"_")[[1]][1]),
                               sample = samples[1] )
venn2 <- commonCellHistonMarks(mod1 = obj.ls.qc[[paste0(modalities[1],"_",samples[2])]], 
                               name_mod1 = strsplit(modalities[1],"_")[[1]][1],
                               mod2 = obj.ls.qc[[paste0(modalities[2],"_",samples[2])]], 
                               name_mod2 = gsub("H3","",strsplit(modalities[2],"_")[[1]][1]),
                               mod3 = obj.ls.qc[[paste0(modalities[3],"_",samples[2])]], 
                               name_mod3 = gsub("H3","",strsplit(modalities[3],"_")[[1]][1]),
                               sample = samples[2] )
ggarrange(venn1,venn2,ncol=2)

# make subset here-------------------------------------------------------------------
# Extract the common cell names between the two Seurat objects (we analyze cells that passed QC of both Dam and H3K27ac)
common_cells_rep1 <- intersect(
  Cells(obj.ls.qc[[name_1]]), 
  Cells(obj.ls.qc[[name_2]]) 
)
common_cells_rep2 <- intersect(
  Cells(obj.ls.qc[[name_3]]), 
  Cells(obj.ls.qc[[name_4]]) 
)

# Create a new list to store the subsetted objects
obj.ls.qc.subset <- list()
obj.ls.qc.subset[[name_1]] <- subset(
  x = obj.ls.qc[[name_1]],
  cells = common_cells_rep1
)
obj.ls.qc.subset[[name_2]] <- subset(
  x = obj.ls.qc[[name_2]],
  cells = common_cells_rep1
)
obj.ls.qc.subset[[name_3]] <- subset(
  x = obj.ls.qc[[name_3]],
  cells = common_cells_rep2
)
obj.ls.qc.subset[[name_4]] <- subset(
  x = obj.ls.qc[[name_4]],
  cells = common_cells_rep2
)

# Check the result
obj.ls.qc.subset

plotCounts(obj = obj.ls.qc.subset, quantiles = c(quant_low,quant_high), feature = "logUMI")
plotCounts(obj = obj.ls.qc.subset, quantiles = c(quant_low), feature = "peak_ratio_MB",ylabel = "% UMI in peaks")

#Merge seurat opjects (if only one replicate, no merging happens)--------------------------------------------------------------
combined.obj.ls <- list()
for (mod in modalities) {
  objs <- lapply(samples, function(s) obj.ls.qc.subset[[paste0(mod, "_", s)]])
  if (length(objs) > 1) {
    combined.obj.ls[[mod]] <- merge(
      x = objs[[1]],
      y = objs[[2]],
      add.cell.ids = samples[1:2]
    )
  } else {
    single.obj <- objs[[1]]
    single.obj <- RenameCells(single.obj, add.cell.id = samples[1])
    combined.obj.ls[[mod]] <- single.obj
  }
}

combined.obj.ls
combined.obj.ls[["Day7R.K27ac_aa"]][["peaks"]]
combined.obj.ls[["Day7R.Dam_bb"]][["peaks"]]

saveRDS(combined.obj.ls, file = "combined.obj.ls.rds")
