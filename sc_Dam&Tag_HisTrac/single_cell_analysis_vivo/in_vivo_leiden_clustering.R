# In vivo cell type annotation and module-score analysis
# Input object: nanoscope_final_peaks_nanoscopefilter_in_vivo.rds
# This script assumes that Dam/K27ac merged objects, clustering, stage metadata,
# and gene-activity assays (RNA.Dam and RNA.K27ac) have already been generated.
#
# Cell-type module scores are derived from adult snapshot DamLeo1 Dam&Tag marker tables.
# The marker tables should be downloaded/exported beforehand and placed in a local directory.
# Example directory layout:
#   marker_tables/
#     adult_snapshot_DamLeo1_DamTag/
#       FindMarkers_Cluster0_vsRest_UP.csv
#       FindMarkers_Cluster0_vsRest_DOWN.csv
#       ...

library(Seurat)
library(ggplot2)
library(patchwork)

# -------------------------------------------------------------------
# Input
# -------------------------------------------------------------------
combined.obj.ls <- readRDS("nanoscope_final_peaks_nanoscopefilter_in_vivo.rds")

# -------------------------------------------------------------------
# Marker-table input from adult snapshot DamLeo1 Dam&Tag
# -------------------------------------------------------------------
marker_dir <- "marker_tables/adult_snapshot_DamLeo1_DamTag"

Dam_L2_3_UP        <- read.csv(file.path(marker_dir, "FindMarkers_Cluster0_vsRest_UP.csv"), row.names = 1)
Dam_L2_3_DOWN      <- read.csv(file.path(marker_dir, "FindMarkers_Cluster0_vsRest_DOWN.csv"), row.names = 1)
Dam_L4_5_UP        <- read.csv(file.path(marker_dir, "FindMarkers_Cluster1_vsRest_UP.csv"), row.names = 1)
Dam_L4_5_DOWN      <- read.csv(file.path(marker_dir, "FindMarkers_Cluster1_vsRest_DOWN.csv"), row.names = 1)
Dam_L5_Pou3f1_UP   <- read.csv(file.path(marker_dir, "FindMarkers_Cluster2_vsRest_UP.csv"), row.names = 1)
Dam_L5_Pou3f1_DOWN <- read.csv(file.path(marker_dir, "FindMarkers_Cluster2_vsRest_DOWN.csv"), row.names = 1)
Dam_L6_Foxp2_UP    <- read.csv(file.path(marker_dir, "FindMarkers_Cluster3_vsRest_UP.csv"), row.names = 1)
Dam_L6_Foxp2_DOWN  <- read.csv(file.path(marker_dir, "FindMarkers_Cluster3_vsRest_DOWN.csv"), row.names = 1)
Dam_Sst_UP         <- read.csv(file.path(marker_dir, "FindMarkers_Cluster4_vsRest_UP.csv"), row.names = 1)
Dam_Sst_DOWN       <- read.csv(file.path(marker_dir, "FindMarkers_Cluster4_vsRest_DOWN.csv"), row.names = 1)
Dam_Vip_Reln_UP    <- read.csv(file.path(marker_dir, "FindMarkers_Cluster5_vsRest_UP.csv"), row.names = 1)
Dam_Vip_Reln_DOWN  <- read.csv(file.path(marker_dir, "FindMarkers_Cluster5_vsRest_DOWN.csv"), row.names = 1)
Dam_Pvalb_UP       <- read.csv(file.path(marker_dir, "FindMarkers_Cluster6_vsRest_UP.csv"), row.names = 1)
Dam_Pvalb_DOWN     <- read.csv(file.path(marker_dir, "FindMarkers_Cluster6_vsRest_DOWN.csv"), row.names = 1)
Dam_Astrocyte_UP   <- read.csv(file.path(marker_dir, "FindMarkers_Cluster7_vsRest_UP.csv"), row.names = 1)
Dam_Astrocyte_DOWN <- read.csv(file.path(marker_dir, "FindMarkers_Cluster7_vsRest_DOWN.csv"), row.names = 1)
Dam_L5_Tshz2_UP    <- read.csv(file.path(marker_dir, "FindMarkers_Cluster8_vsRest_UP.csv"), row.names = 1)
Dam_L5_Tshz2_DOWN  <- read.csv(file.path(marker_dir, "FindMarkers_Cluster8_vsRest_DOWN.csv"), row.names = 1)

# -------------------------------------------------------------------
# Module scores from adult snapshot DamLeo1 Dam&Tag markers
# -------------------------------------------------------------------
DefaultAssay(combined.obj.ls$K27ac_aa) <- "RNA.K27ac"
DefaultAssay(combined.obj.ls$Dam_bb)   <- "RNA.Dam"

score_sets <- list(
  Dam_L2_3      = list(up = rownames(Dam_L2_3_UP),        down = rownames(Dam_L2_3_DOWN)),
  Dam_L4_5      = list(up = rownames(Dam_L4_5_UP),        down = rownames(Dam_L4_5_DOWN)),
  Dam_L5_Pou3f1 = list(up = rownames(Dam_L5_Pou3f1_UP),   down = rownames(Dam_L5_Pou3f1_DOWN)),
  Dam_L6_Foxp2  = list(up = rownames(Dam_L6_Foxp2_UP),    down = rownames(Dam_L6_Foxp2_DOWN)),
  Dam_Sst       = list(up = rownames(Dam_Sst_UP),         down = rownames(Dam_Sst_DOWN)),
  Dam_Vip_Reln  = list(up = rownames(Dam_Vip_Reln_UP),    down = rownames(Dam_Vip_Reln_DOWN)),
  Dam_Pvalb     = list(up = rownames(Dam_Pvalb_UP),       down = rownames(Dam_Pvalb_DOWN)),
  Dam_Astrocyte = list(up = rownames(Dam_Astrocyte_UP),   down = rownames(Dam_Astrocyte_DOWN)),
  Dam_L5_Tshz2  = list(up = rownames(Dam_L5_Tshz2_UP),    down = rownames(Dam_L5_Tshz2_DOWN))
)

add_signed_module_scores <- function(obj, assay_name, score_sets) {
  DefaultAssay(obj) <- assay_name
  for (nm in names(score_sets)) {
    obj <- AddModuleScore(obj, features = list(score_sets[[nm]]$up), name = paste0(nm, "_up_GeneActivity"))
    obj <- AddModuleScore(obj, features = list(score_sets[[nm]]$down), name = paste0(nm, "_down_GeneActivity"))
    obj[[paste0(nm, "_GeneActivity")]] <-
      obj[[paste0(nm, "_up_GeneActivity1")]][, 1] - obj[[paste0(nm, "_down_GeneActivity1")]][, 1]
  }
  obj
}

combined.obj.ls$K27ac_aa <- add_signed_module_scores(combined.obj.ls$K27ac_aa, "RNA.K27ac", score_sets)
combined.obj.ls$Dam_bb   <- add_signed_module_scores(combined.obj.ls$Dam_bb,   "RNA.Dam",   score_sets)

# -------------------------------------------------------------------
# Example visualization of module scores
# -------------------------------------------------------------------
FeaturePlot(
  object = combined.obj.ls$K27ac_aa,
  features = c(
    "Dam_L2_3_GeneActivity", "Dam_L4_5_GeneActivity", "Dam_L5_Pou3f1_GeneActivity",
    "Dam_L5_Tshz2_GeneActivity", "Dam_L6_Foxp2_GeneActivity", "Dam_Sst_GeneActivity",
    "Dam_Vip_Reln_GeneActivity", "Dam_Pvalb_GeneActivity", "Dam_Astrocyte_GeneActivity"
  ),
  pt.size = 0.1, max.cutoff = "q95", min.cutoff = "q01", ncol = 3, label = FALSE
)

FeaturePlot(
  object = combined.obj.ls$Dam_bb,
  features = c(
    "Dam_L2_3_GeneActivity", "Dam_L4_5_GeneActivity", "Dam_L5_Pou3f1_GeneActivity",
    "Dam_L5_Tshz2_GeneActivity", "Dam_L6_Foxp2_GeneActivity", "Dam_Sst_GeneActivity",
    "Dam_Vip_Reln_GeneActivity", "Dam_Pvalb_GeneActivity", "Dam_Astrocyte_GeneActivity"
  ),
  pt.size = 0.1, max.cutoff = "q95", min.cutoff = "q01", ncol = 3, label = FALSE
)

# -------------------------------------------------------------------
# Manual cell-type annotation from cluster identities
# -------------------------------------------------------------------
# K27ac modality
obj <- combined.obj.ls$K27ac_aa
cluster_map <- c(
  "3" = "L2-3", "4" = "L2-3", "17" = "L2-3", "31" = "L2-3",
  "8" = "L4-5", "11" = "L4-5", "14" = "L4-5", "22" = "L4-5", "24" = "L4-5", "38" = "L4-5", "39" = "L4-5", "43" = "L4-5",
  "5" = "L6-Foxp2", "20" = "L6-Foxp2", "23" = "L6-Foxp2", "25" = "L6-Foxp2", "41" = "L6-Foxp2",
  "12" = "L5-Pou3f1", "28" = "L5-Pou3f1", "29" = "L5-Pou3f1", "42" = "L5-Pou3f1",
  "35" = "L5-Tshz2",
  "6" = "Sst", "27" = "Sst",
  "7" = "Pvalb", "21" = "Pvalb",
  "19" = "Vip-Reln",
  "16" = "Astrocyte", "26" = "Astrocyte", "33" = "Astrocyte",
  "15" = "P8_W8R_1", "36" = "P8_W8R_1",
  "18" = "P8_W8R_2", "30" = "P8_W8R_2", "34" = "P8_W8R_2",
  "9" = "W8R_1", "13" = "W8R_1",
  "0" = "P8_1", "1" = "P8_1",
  "2" = "P8_2", "10" = "P8_2", "32" = "P8_2", "37" = "P8_2", "40" = "P8_2"
)
vals <- cluster_map[as.character(obj$seurat_clusters)]
names(vals) <- colnames(obj)
obj <- AddMetaData(obj, metadata = vals, col.name = "K27ac_Annotation_cell_type")
levels_K27ac <- c("L2-3", "L4-5", "L5-Pou3f1", "L5-Tshz2", "L6-Foxp2",
                  "Sst", "Vip-Reln", "Pvalb", "Astrocyte",
                  "P8_W8R_1", "P8_W8R_2", "W8R_1", "P8_1", "P8_2")
obj$K27ac_Annotation_cell_type <- factor(as.character(obj$K27ac_Annotation_cell_type), levels = levels_K27ac)
combined.obj.ls$K27ac_aa <- obj

# Dam modality
obj <- combined.obj.ls$Dam_bb
cluster_map <- c(
  "0" = "P8_W8R_1", "1" = "L2-3", "2" = "L6-Foxp2", "3" = "L5-Pou3f1",
  "4" = "L4-5", "5" = "Sst", "6" = "P8_W8R_2", "7" = "Pvalb",
  "8" = "Astrocyte", "9" = "P8_W8R_3", "10" = "L6-Foxp2", "11" = "L4-5",
  "12" = "Vip-Reln", "13" = "L4-5", "14" = "L2-3"
)
vals <- cluster_map[as.character(obj$seurat_clusters)]
names(vals) <- colnames(obj)
obj <- AddMetaData(obj, metadata = vals, col.name = "Dam_Annotation_cell_type")
levels_Dam <- c("L2-3", "L4-5", "L5-Pou3f1", "L5-Tshz2", "L6-Foxp2",
                "Sst", "Vip-Reln", "Pvalb", "Astrocyte",
                "P8_W8R_1", "P8_W8R_2", "P8_W8R_3")
obj$Dam_Annotation_cell_type <- factor(as.character(obj$Dam_Annotation_cell_type), levels = levels_Dam)
combined.obj.ls$Dam_bb <- obj

# -------------------------------------------------------------------
# Transfer cell-type annotations between modalities
# -------------------------------------------------------------------
idents_1 <- combined.obj.ls$K27ac_aa$K27ac_Annotation_cell_type
idents_2 <- combined.obj.ls$Dam_bb$Dam_Annotation_cell_type
combined.obj.ls$K27ac_aa <- AddMetaData(combined.obj.ls$K27ac_aa, idents_2, "Dam_Annotation_cell_type")
combined.obj.ls$Dam_bb   <- AddMetaData(combined.obj.ls$Dam_bb,   idents_1, "K27ac_Annotation_cell_type")

p1 <- DimPlot(combined.obj.ls$K27ac_aa, group.by = "K27ac_Annotation_cell_type", label = FALSE) + ggtitle("Cell types based on K27ac assay")
p2 <- DimPlot(combined.obj.ls$K27ac_aa, group.by = "Dam_Annotation_cell_type", label = FALSE) + ggtitle("Cell types based on Dam assay")
p3 <- DimPlot(combined.obj.ls$Dam_bb,   group.by = "K27ac_Annotation_cell_type", label = FALSE) + ggtitle("Cell types based on K27ac assay")
p4 <- DimPlot(combined.obj.ls$Dam_bb,   group.by = "Dam_Annotation_cell_type", label = FALSE) + ggtitle("Cell types based on Dam assay")
p1 + p2 + p3 + p4

# Example stage-restricted visualization
for (st in c("P8", "W8", "W8R")) {
  k27_sub <- subset(combined.obj.ls$K27ac_aa, subset = stage == st)
  dam_sub <- subset(combined.obj.ls$Dam_bb,   subset = stage == st)

  p1 <- DimPlot(k27_sub, group.by = "K27ac_Annotation_cell_type", label = FALSE) + ggtitle(paste0(st, ": K27ac labels"))
  p2 <- DimPlot(k27_sub, group.by = "Dam_Annotation_cell_type",   label = FALSE) + ggtitle(paste0(st, ": Dam labels"))
  p3 <- DimPlot(dam_sub, group.by = "K27ac_Annotation_cell_type", label = FALSE) + ggtitle(paste0(st, ": K27ac labels"))
  p4 <- DimPlot(dam_sub, group.by = "Dam_Annotation_cell_type",   label = FALSE) + ggtitle(paste0(st, ": Dam labels"))
  print(p1 + p2 + p3 + p4)
}

saveRDS(combined.obj.ls, file = "nanoscope_final_peaks_nanoscopefilter_in_vivo_leiden_cluster.rds")
