# In vivo downstream analysis after sample-wise QC and modality matching
# Merging across 12 samples (P8, W8, W8R) for Dam and K27ac modalities
# Cell filtering and sample-wise modality matching are assumed to have been completed beforehand
# Example input: nanoscope_separated_nanoscopefilter_260218_1percentile.rds
# Modified pipeline based on nanoscope: https://fansalon.github.io/vignette_single-cell-nanoCT.html

# Required packages for this script itself.
# Additional packages may be needed internally by functions_scCT2.R.
library(Signac)
library(Seurat)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(ggpubr)
library(EnsDb.Mmusculus.v79)

# directory where the nanoscope repo was cloned
repodir <- "/path/to/your/nanoscope/"
source(paste0(repodir, "scripts/functions_scCT2.R"))

setwd("path/to/your/nanoscope/project/")

# -------------------------------------------------------------------
# Input
# -------------------------------------------------------------------
# This object is assumed to contain sample-separated Seurat objects that have:
#   1) passed nanoscope-based filtering,
#   2) been matched between Dam and K27ac within each sample.
obj.ls.qc.subset <- readRDS("nanoscope_separated_nanoscopefilter_in_vivo.rds")
obj.ls.qc.subset

# -------------------------------------------------------------------
# Merge Seurat objects across all in vivo samples for each modality
# -------------------------------------------------------------------
# In vivo analysis here includes 12 samples in total:
#   - P8 snapshot:   P8rep1-4
#   - W8 snapshot:   W8rep1-4
#   - W8R HisTrac:   W8Rrep1-4
#
# The same set of samples is merged independently for the two modalities.
modality2samples <- list(
  "K27ac_aa" = c("P8rep1", "P8rep2", "P8rep3", "P8rep4",
                 "W8rep1", "W8rep2", "W8rep3", "W8rep4",
                 "W8Rrep1", "W8Rrep2", "W8Rrep3", "W8Rrep4"),
  "Dam_bb"   = c("P8rep1", "P8rep2", "P8rep3", "P8rep4",
                 "W8rep1", "W8rep2", "W8rep3", "W8rep4",
                 "W8Rrep1", "W8Rrep2", "W8Rrep3", "W8Rrep4")
)

combined.obj.ls <- list()
for (mod in names(modality2samples)) {
  these_samples <- modality2samples[[mod]]

  seurat_list <- lapply(these_samples, function(smpl) {
    key <- paste0(mod, "_", smpl)
    if (!key %in% names(obj.ls.qc.subset)) {
      stop("No object found for: ", key)
    }
    obj.ls.qc.subset[[key]]
  })

  if (length(seurat_list) == 1) {
    merged_obj <- seurat_list[[1]]
  } else {
    merged_obj <- merge(
      x = seurat_list[[1]],
      y = seurat_list[-1],
      add.cell.ids = these_samples
    )
  }

  combined.obj.ls[[mod]] <- merged_obj
}
combined.obj.ls

# -------------------------------------------------------------------
# TF-IDF normalization, variable features, and LSI
# -------------------------------------------------------------------
combined.obj.ls <- lapply(combined.obj.ls, RunTFIDF)
combined.obj.ls <- lapply(combined.obj.ls, FindTopFeatures)
combined.obj.ls <- lapply(combined.obj.ls, RunSVD)

combined.obj.ls
combined.obj.ls[["K27ac_aa"]][["peaks"]]
combined.obj.ls[["Dam_bb"]][["peaks"]]

# Inspect the relationship between sequencing depth and LSI components.
# This helps identify dimensions strongly driven by depth.
plot.list <- lapply(combined.obj.ls, DepthCorMulMod)
ggarrange(plot.list[[1]], plot.list[[2]])

# Inspect the elbow structure of the LSI reduction.
# The number of LSI dimensions used downstream is chosen based on the elbow plots,
# together with the depth-correlation plots above.
plot.list_elbow <- lapply(combined.obj.ls, ElbowPlot, reduction = "lsi", ndims = 50)
ggarrange(plot.list_elbow[[1]], plot.list_elbow[[2]])

# In downstream analyses, the LSI dimensions are selected based on these diagnostics.
# In this example, dims = 2:15 are used, because the first component often reflects
# depth or other technical effects.

# -------------------------------------------------------------------
# UMAP and clustering for each modality
# -------------------------------------------------------------------
DefaultAssay(combined.obj.ls$K27ac_aa) <- "peaks"
DefaultAssay(combined.obj.ls$Dam_bb)   <- "peaks"

combined.obj.ls$K27ac_aa <- RunUMAP(combined.obj.ls$K27ac_aa, reduction = "lsi", dims = 2:15)
combined.obj.ls$K27ac_aa <- FindNeighbors(combined.obj.ls$K27ac_aa, reduction = "lsi", dims = 2:15)
combined.obj.ls$K27ac_aa <- FindClusters(combined.obj.ls$K27ac_aa, verbose = FALSE, algorithm = 3, resolution = 2)

combined.obj.ls$Dam_bb <- RunUMAP(combined.obj.ls$Dam_bb, reduction = "lsi", dims = 2:15)
combined.obj.ls$Dam_bb <- FindNeighbors(combined.obj.ls$Dam_bb, reduction = "lsi", dims = 2:15)
combined.obj.ls$Dam_bb <- FindClusters(combined.obj.ls$Dam_bb, verbose = FALSE, algorithm = 3, resolution = 0.5)

table(Idents(combined.obj.ls$K27ac_aa))
table(Idents(combined.obj.ls$Dam_bb))

# Plot each modality separately
p1 <- DimPlot(combined.obj.ls$K27ac_aa, label = TRUE) + NoLegend() +
  ggtitle("K27ac_aa") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))

p2 <- DimPlot(combined.obj.ls$Dam_bb, label = TRUE) + NoLegend() +
  ggtitle("Dam_bb") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))

p3 <- FeaturePlot(combined.obj.ls$K27ac_aa, features = "passed_filters") +
  ggtitle("Read Depth on UMAP (Peaks Assay)")

p4 <- FeaturePlot(combined.obj.ls$Dam_bb, features = "passed_filters") +
  ggtitle("Read Depth on UMAP (Peaks Assay)")

p5 <- DimPlot(combined.obj.ls$K27ac_aa, group.by = "sample") +
  ggtitle("K27ac_aa (colored by sample)")

p6 <- DimPlot(combined.obj.ls$Dam_bb, group.by = "sample") +
  ggtitle("Dam_bb (colored by sample)")

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3)

# -------------------------------------------------------------------
# Transfer cluster labels between modalities
# -------------------------------------------------------------------
# Cluster identities are first defined independently in each modality.
# These labels are then copied to the paired modality so that each cell carries
# both the K27ac-based and Dam-based cluster annotations.
combined.obj.ls$K27ac_aa <- AddMetaData(
  combined.obj.ls$K27ac_aa,
  Idents(combined.obj.ls$K27ac_aa),
  "K27ac_annotation"
)
combined.obj.ls$Dam_bb <- AddMetaData(
  combined.obj.ls$Dam_bb,
  Idents(combined.obj.ls$Dam_bb),
  "Dam_annotation"
)

idents_1 <- combined.obj.ls$K27ac_aa$K27ac_annotation
idents_2 <- combined.obj.ls$Dam_bb$Dam_annotation

combined.obj.ls$K27ac_aa <- AddMetaData(combined.obj.ls$K27ac_aa, idents_2, "Dam_annotation")
combined.obj.ls$Dam_bb   <- AddMetaData(combined.obj.ls$Dam_bb, idents_1, "K27ac_annotation")

p1 <- DimPlot(combined.obj.ls$K27ac_aa, group.by = "K27ac_annotation", label = TRUE) +
  ggtitle("Clusters based on K27ac assay")
p2 <- DimPlot(combined.obj.ls$K27ac_aa, group.by = "Dam_annotation", label = TRUE) +
  ggtitle("Clusters based on Dam assay")
p3 <- DimPlot(combined.obj.ls$Dam_bb, group.by = "K27ac_annotation", label = TRUE) +
  ggtitle("Clusters based on K27ac assay")
p4 <- DimPlot(combined.obj.ls$Dam_bb, group.by = "Dam_annotation", label = TRUE) +
  ggtitle("Clusters based on Dam assay")
p1 + p2 + p3 + p4

# -------------------------------------------------------------------
# Add stage metadata
# -------------------------------------------------------------------
# Because multiple developmental conditions are merged into one object,
# a simplified stage label is derived from orig.ident.
# Example:
#   P8rep1  -> P8
#   W8rep2  -> W8
#   W8Rrep3 -> W8R
combined.obj.ls$K27ac_aa$stage <- sub("rep\\d+$", "", combined.obj.ls$K27ac_aa$orig.ident)
combined.obj.ls$Dam_bb$stage   <- sub("rep\\d+$", "", combined.obj.ls$Dam_bb$orig.ident)

table(combined.obj.ls$K27ac_aa$stage)
table(combined.obj.ls$Dam_bb$stage)

p1 <- DimPlot(combined.obj.ls$K27ac_aa, group.by = "K27ac_annotation", label = TRUE) +
  ggtitle("Clusters based on K27ac assay")
p2 <- DimPlot(combined.obj.ls$Dam_bb, group.by = "Dam_annotation", label = TRUE) +
  ggtitle("Clusters based on Dam assay")
p3 <- DimPlot(combined.obj.ls$K27ac_aa, group.by = "stage") +
  ggtitle("K27ac_aa (colored by stage)")
p4 <- DimPlot(combined.obj.ls$Dam_bb, group.by = "stage") +
  ggtitle("Dam_bb (colored by stage)")
p1 + p2 + p3 + p4

# -------------------------------------------------------------------
# Gene annotation and gene activity calculation
# -------------------------------------------------------------------
# Gene activity is computed from the peaks assay and stored as a separate assay
# for each modality.
genome_ann <- EnsDb.Mmusculus.v79
annotations <- GetGRangesFromEnsDb(ensdb = genome_ann)
seqlevelsStyle(annotations) <- "UCSC"

Annotation(combined.obj.ls[["K27ac_aa"]]) <- annotations
Annotation(combined.obj.ls[["Dam_bb"]])   <- annotations

# K27ac gene activity
DefaultAssay(combined.obj.ls$K27ac_aa) <- "peaks"
gene.activities <- GeneActivity(
  combined.obj.ls$K27ac_aa,
  extend.upstream = 2000,
  extend.downstream = 0
)
combined.obj.ls$K27ac_aa[["RNA.K27ac"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(combined.obj.ls$K27ac_aa) <- "RNA.K27ac"
combined.obj.ls$K27ac_aa <- NormalizeData(
  object = combined.obj.ls$K27ac_aa,
  assay = "RNA.K27ac",
  normalization.method = "LogNormalize",
  scale.factor = median(combined.obj.ls$K27ac_aa$nCount_RNA.K27ac)
)

# Dam gene activity
DefaultAssay(combined.obj.ls$Dam_bb) <- "peaks"
gene.activities <- GeneActivity(
  combined.obj.ls$Dam_bb,
  extend.upstream = 2000,
  extend.downstream = 2000
)
combined.obj.ls$Dam_bb[["RNA.Dam"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(combined.obj.ls$Dam_bb) <- "RNA.Dam"
combined.obj.ls$Dam_bb <- NormalizeData(
  object = combined.obj.ls$Dam_bb,
  assay = "RNA.Dam",
  normalization.method = "LogNormalize",
  scale.factor = median(combined.obj.ls$Dam_bb$nCount_RNA.Dam)
)

combined.obj.ls
combined.obj.ls[["K27ac_aa"]][["RNA.K27ac"]]
combined.obj.ls[["Dam_bb"]][["RNA.Dam"]]

# -------------------------------------------------------------------
# Example marker visualization in the gene activity assays
# -------------------------------------------------------------------
DefaultAssay(combined.obj.ls$K27ac_aa) <- "RNA.K27ac"
FeaturePlot(
  object = combined.obj.ls$K27ac_aa,
  features = c("Rbfox3", "Map2", "Slc17a7", "Slc17a6", "Neurod6", "Gad1", "Gad2", "Slc1a3", "Foxp2"),
  pt.size = 0.1, max.cutoff = "q95", min.cutoff = "q10", ncol = 3, label = FALSE
)
FeaturePlot(
  object = combined.obj.ls$K27ac_aa,
  features = c("Satb2", "Cux2", "Rorb", "Etv1", "Tshz2", "Bcl11b", "Pou3f2", "Tbr1", "Syt6"),
  pt.size = 0.1, max.cutoff = "q95", ncol = 3, label = FALSE
)
FeaturePlot(
  object = combined.obj.ls$K27ac_aa,
  features = c("Oprk1", "Cdh9", "Pou3f1", "Dlx5", "Dlx6", "Sst", "Reln", "Vip", "Pvalb"),
  pt.size = 0.1, max.cutoff = "q95", ncol = 3, label = FALSE
)
FeaturePlot(
  object = combined.obj.ls$K27ac_aa,
  features = c("Mfge8", "Mag", "P2ry12", "C1qc", "Slc1a3", "Pdgfra", "Mobp", "Enpp6", "Olig1"),
  pt.size = 0.1, max.cutoff = "q95", ncol = 3, label = FALSE
)

DefaultAssay(combined.obj.ls$Dam_bb) <- "RNA.Dam"
FeaturePlot(
  object = combined.obj.ls$Dam_bb,
  features = c("Rbfox3", "Map2", "Slc17a7", "Slc17a6", "Neurod6", "Gad1", "Gad2", "Slc1a3", "Foxp2"),
  pt.size = 0.1, max.cutoff = "q95", min.cutoff = "q10", ncol = 3, label = FALSE
)
FeaturePlot(
  object = combined.obj.ls$Dam_bb,
  features = c("Satb2", "Cux2", "Rorb", "Etv1", "Tshz2", "Bcl11b", "Pou3f2", "Tbr1", "Syt6"),
  pt.size = 0.1, max.cutoff = "q95", ncol = 3, label = FALSE
)
FeaturePlot(
  object = combined.obj.ls$Dam_bb,
  features = c("Oprk1", "Cdh9", "Pou3f1", "Dlx5", "Dlx6", "Sst", "Reln", "Vip", "Pvalb"),
  pt.size = 0.1, max.cutoff = "q95", ncol = 3, label = FALSE
)
FeaturePlot(
  object = combined.obj.ls$Dam_bb,
  features = c("Mfge8", "Mag", "P2ry12", "C1qc", "Slc1a3", "Pdgfra", "Mobp", "Enpp6", "Olig1"),
  pt.size = 0.1, max.cutoff = "q95", ncol = 3, label = FALSE
)

# We can use plotConnectModal (nanoscope) to connect identical cells between two modalities with lines.
plotConnectModal(seurat = combined.obj.ls, group = "K27ac_annotation")
plotConnectModal(seurat = combined.obj.ls, group = "Dam_annotation")

saveRDS(combined.obj.ls, file = "nanoscope_final_peaks_nanoscopefilter_in_vivo.rds")
