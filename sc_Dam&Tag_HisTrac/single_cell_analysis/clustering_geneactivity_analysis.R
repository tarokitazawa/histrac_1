# Seurat object preparation from nanoscope outputs
# clustering, umap, gene activity
# Merging of rep1 and rep2 (example - DamLeo1-K27ac Day7R (Retrospective (=HisTrac)) sample)
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

# Seurat object was prepared by Seurat_object_preparation.R
combined.obj.ls <-readRDS("combined.obj.ls.rds")

# Normalisation
combined.obj.ls <- lapply(combined.obj.ls,RunTFIDF)
# a warning (Some features contain 0 total counts) might rise. This is due to the cell filtering performed in the QC step. Should be safe to proceed

# Variable features
combined.obj.ls <- lapply(combined.obj.ls,FindTopFeatures)

# Dim reduction
combined.obj.ls <- lapply(combined.obj.ls,RunSVD)

# Inspect the relationship between sequencing depth and latent semantic indexing (LSI) components.
# This helps identify components that are strongly driven by depth.
plot.list <- lapply(combined.obj.ls, DepthCorMulMod)
ggarrange(plot.list[[1]], plot.list[[2]])

# Inspect the variance explained across LSI dimensions.
# The number of dimensions used in downstream analyses is chosen based on the elbow plots,
# together with the depth-correlation plots above.
plot.list_elbow <- lapply(combined.obj.ls, ElbowPlot, reduction = "lsi", ndims = 50)
ggarrange(plot.list_elbow[[1]], plot.list_elbow[[2]])

# In downstream analyses, the LSI dimensions to use should be selected based on these diagnostics.
# In practice, the first depth-associated component is often excluded, and the remaining dimensions
# are chosen according to the elbow structure.

# Run dimension reduction and clustering on each modality individually

DefaultAssay(combined.obj.ls$Day7R.K27ac_aa) <- "peaks"
DefaultAssay(combined.obj.ls$Day7R.Dam_bb) <- "peaks"

# K27ac (adjust parameter, e.g., dims, resolution)
# Since the first LSI component often reflects depth or other batch-related effects, we set dims = 2:15 as intended.
combined.obj.ls$Day7R.K27ac_aa <- RunUMAP(combined.obj.ls$Day7R.K27ac_aa,reduction = 'lsi', dims = 2:15)
combined.obj.ls$Day7R.K27ac_aa <- FindNeighbors(combined.obj.ls$Day7R.K27ac_aa,reduction = 'lsi', dims = 2:15)
combined.obj.ls$Day7R.K27ac_aa <- FindClusters(combined.obj.ls$Day7R.K27ac_aa, verbose = FALSE, algorithm = 3, resolution = 0.2)

# Dam
# Since the first LSI component often reflects depth or other batch-related effects, we set dims = 2:15 as intended.
combined.obj.ls$Day7R.Dam_bb <- RunUMAP(combined.obj.ls$Day7R.Dam_bb,reduction = 'lsi', dims = 2:15)
combined.obj.ls$Day7R.Dam_bb <- FindNeighbors(combined.obj.ls$Day7R.Dam_bb,reduction = 'lsi', dims = 2:15)
combined.obj.ls$Day7R.Dam_bb <- FindClusters(combined.obj.ls$Day7R.Dam_bb, verbose = FALSE, algorithm = 3, resolution = 0.2)

table(Idents(combined.obj.ls$Day7R.K27ac_aa))
table(Idents(combined.obj.ls$Day7R.Dam_bb))

# Plot each modality separately
p1=DimPlot(combined.obj.ls$Day7R.K27ac_aa, label = TRUE) + NoLegend() +
  ggtitle("Day7R.K27ac_aa") +
  theme(plot.title = element_text(size=15,hjust=0.5,face='bold'))

p2=DimPlot(combined.obj.ls$Day7R.Dam_bb, label = TRUE) + NoLegend() +
  ggtitle("Day7R.Dam_bb") +
  theme(plot.title = element_text(size=15,hjust=0.5,face='bold'))

p3=FeaturePlot(
  object = combined.obj.ls$Day7R.K27ac_aa,
  features = 'passed_filters'
) + ggtitle("Read Depth on UMAP (Peaks Assay)")

p4=FeaturePlot(
  object = combined.obj.ls$Day7R.Dam_bb,
  features = 'passed_filters'
) + ggtitle("Read Depth on UMAP (Peaks Assay)")

p5 <- DimPlot(combined.obj.ls$Day7R.K27ac_aa,  group.by = "sample") + 
  ggtitle("Day7R.K27ac_aa (colored by sample)")

p6 <- DimPlot(combined.obj.ls$Day7R.Dam_bb,  group.by = "sample") + 
  ggtitle("Day7R.Dam_bb (colored by sample)")

p1 + p3 + p5 + p2 + p4 + p6


#--- Cell type tnransfer -------------------------------------------------------------------------

# Transfer of cell type annotation
colnames(combined.obj.ls$Day7R.K27ac_aa@meta.data)
colnames(combined.obj.ls$Day7R.Dam_bb@meta.data)

# Add  annotations to the object
combined.obj.ls$Day7R.K27ac_aa <- AddMetaData(combined.obj.ls$Day7R.K27ac_aa,
                                              Idents(combined.obj.ls$Day7R.K27ac_aa),"K27ac_annotation")
combined.obj.ls$Day7R.Dam_bb <- AddMetaData(combined.obj.ls$Day7R.Dam_bb,
                                            Idents(combined.obj.ls$Day7R.Dam_bb),"Dam_annotation")

idents_1 <- combined.obj.ls$Day7R.K27ac_aa$K27ac_annotation
idents_2 <- combined.obj.ls$Day7R.Dam_bb$Dam_annotation

head(idents_1)
head(idents_2)
summary(idents_1)
summary(idents_2)

combined.obj.ls$Day7R.K27ac_aa <- AddMetaData(combined.obj.ls$Day7R.K27ac_aa,
                                              idents_2,"Dam_annotation")
combined.obj.ls$Day7R.Dam_bb <- AddMetaData(combined.obj.ls$Day7R.Dam_bb,
                                            idents_1,"K27ac_annotation")
# Plot

p1 <- DimPlot(
  object = combined.obj.ls$Day7R.K27ac_aa,
  group.by = "K27ac_annotation",
  label = TRUE
) + ggtitle("Clusters based on K27ac Assay")

p2 <- DimPlot(
  object = combined.obj.ls$Day7R.K27ac_aa,
  group.by = "Dam_annotation",
  label = TRUE
) + ggtitle("Clusters based on DamDam Assay")

p3 <- DimPlot(
  object = combined.obj.ls$Day7R.Dam_bb,
  group.by = "K27ac_annotation",
  label = TRUE
) + ggtitle("Clusters based on K27ac Assay")

p4 <- DimPlot(
  object = combined.obj.ls$Day7R.Dam_bb,
  group.by = "Dam_annotation",
  label = TRUE
) + ggtitle("Clusters based on DamDam Assay")

p1 + p2 + p3 + p4 


# Annotation-------------------------------------------------------------------
# Genome annotations
genome_ann <- EnsDb.Mmusculus.v79
# get gene annotations from Ensembl
annotations <- GetGRangesFromEnsDb(ensdb = genome_ann)

# change to UCSC style since the data was mapped to mm10
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(combined.obj.ls[["Day7R.K27ac_aa"]]) <- annotations
Annotation(combined.obj.ls[["Day7R.Dam_bb"]]) <- annotations

combined.obj.ls
combined.obj.ls[["Day7R.K27ac_aa"]][["peaks"]]
combined.obj.ls[["Day7R.Dam_bb"]][["peaks"]]

# Calculate gene activities for K27ac-----------------------------------------
DefaultAssay(combined.obj.ls$Day7R.K27ac_aa) <- "peaks"
gene.activities <- GeneActivity(combined.obj.ls$Day7R.K27ac_aa,
                                extend.upstream = 2000,
                                extend.downstream = 0)

# Add the gene activity matrix to the Seurat object as a new assay
combined.obj.ls$Day7R.K27ac_aa[['RNA.K27ac']] <- CreateAssayObject(counts = gene.activities)

# normalize RNA data
DefaultAssay(combined.obj.ls$Day7R.K27ac_aa) <- "RNA.K27ac"

combined.obj.ls$Day7R.K27ac_aa <- NormalizeData(
  object = combined.obj.ls$Day7R.K27ac_aa,
  assay = 'RNA.K27ac',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined.obj.ls$Day7R.K27ac_aa$nCount_RNA.K27ac)
)
# Calculate gene activities for Leo1-------------------------------------------
DefaultAssay(combined.obj.ls$Day7R.Dam_bb) <- "peaks"
gene.activities <- GeneActivity(combined.obj.ls$Day7R.Dam_bb,
                                extend.upstream = 2000,
                                extend.downstream = 2000)

# Add the gene activity matrix to the Seurat object as a new assay
combined.obj.ls$Day7R.Dam_bb[['RNA.Dam']] <- CreateAssayObject(counts = gene.activities)

# normalize RNA data
DefaultAssay(combined.obj.ls$Day7R.Dam_bb) <- "RNA.Dam"

combined.obj.ls$Day7R.Dam_bb <- NormalizeData(
  object = combined.obj.ls$Day7R.Dam_bb,
  assay = 'RNA.Dam',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined.obj.ls$Day7R.Dam_bb$nCount_RNA.Dam)
)

combined.obj.ls
combined.obj.ls[["Day7R.K27ac_aa"]][["peaks"]]
combined.obj.ls[["Day7R.K27ac_aa"]][["peaks"]]
combined.obj.ls[["Day7R.K27ac_aa"]][["RNA.K27ac"]]
combined.obj.ls[["Day7R.Dam_bb"]][["RNA.Dam"]]

# We can use plotConnectModal (nanoscope) to connect identical cells between two modalities with lines
plotConnectModal(seurat = combined.obj.ls, group = 'K27ac_annotation')
plotConnectModal(seurat = combined.obj.ls, group = 'Dam_annotation')

saveRDS(combined.obj.ls, file = "combined.obj.ls_analysis.rds")
