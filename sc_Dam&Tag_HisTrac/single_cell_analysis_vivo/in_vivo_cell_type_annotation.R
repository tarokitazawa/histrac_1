# In vivo reference mapping and final object construction
# Input object: nanoscope_final_peaks_nanoscopefilter_in_vivo_leiden_cluster.rds
#
# This script performs the in vivo-specific post-processing steps after initial
# clustering / gene-activity generation:
#   1) remove ambiguous K27ac-defined clusters that are not retained as adult-like references
#   2) use Dam modality to map P8 cells onto an adult reference built from W8 + W8R
#   3) create a joint reference-space LSI representation (joint.ref.lsi)
#   4) remove ambiguously mapped P8 cells
#   5) transfer adult labels / Cell_Type metadata from Dam to K27ac
#   6) re-normalize and re-cluster the filtered K27ac object for visualization
#   7) harmonize visualization UMAP naming across modalities and save the final object
#
# Biological rationale:
# - Early postnatal AAV injection can label populations that are not recovered in
#   the adult snapshot AAV datasets.
# - Therefore, only unambiguous adult-like populations are retained as references.
# - P8 cells are then projected to the adult reference in Dam LSI space, and
#   ambiguous P8 assignments are removed.
# - In the HisTrac condition, the Dam modality reflects earlier molecular states
#   (here, P8), whereas the paired K27ac modality is read out at Week 8 and
#   therefore provides the later-stage ground truth.
# - This makes it possible to assign the historical Dam profiles of HisTrac cells
#   to unambiguous adult identities using their matched Week 8 endpoint state.
# - In other words, although the Dam signal records an earlier P8-like state,
#   the final identity of each HisTrac cell can be resolved by its Week 8
#   K27ac-defined adult outcome.
#
# Note:
# - The supervised K27ac UMAP is used for visualization only.
# - Downstream quantitative analyses should rely on unbiased low-dimensional
#   representations such as LSI / joint.ref.lsi, not on the supervised UMAP.

library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(uwot)
library(dplyr)

# -------------------------------------------------------------------
# Input
# -------------------------------------------------------------------
combined.obj.ls <- readRDS("nanoscope_final_peaks_nanoscopefilter_in_vivo_leiden_cluster.rds")

# -------------------------------------------------------------------
# Remove ambiguous K27ac-defined populations from the adult reference framework
# -------------------------------------------------------------------
# Rationale:
# Adult-like reference identities are defined from unambiguous K27ac clusters only.
# Mixed / ambiguous populations are excluded before building the adult reference.

k27ac_obj <- combined.obj.ls$K27ac_aa
dam_obj   <- combined.obj.ls$Dam_bb

k27ac_obj_sub <- subset(
  x = k27ac_obj,
  subset =
    (K27ac_Annotation_cell_type != "P8_W8R_1") &
    (K27ac_Annotation_cell_type != "P8_W8R_2") &
    (K27ac_Annotation_cell_type != "W8R_1") &
    (K27ac_Annotation_cell_type != "P8_1") &
    (K27ac_Annotation_cell_type != "P8_2")
)
k27ac_obj_sub <- subset(
  x = k27ac_obj_sub,
  subset = (K27ac_annotation != "33")
)

dam_obj_sub <- subset(
  x = dam_obj,
  subset =
    (K27ac_Annotation_cell_type != "P8_W8R_1") &
    (K27ac_Annotation_cell_type != "P8_W8R_2") &
    (K27ac_Annotation_cell_type != "W8R_1") &
    (K27ac_Annotation_cell_type != "P8_1") &
    (K27ac_Annotation_cell_type != "P8_2")
)
dam_obj_sub <- subset(
  x = dam_obj_sub,
  subset = (K27ac_annotation != "33")
)

combined.obj.ls.2 <- list(
  K27ac_aa = k27ac_obj_sub,
  Dam_bb   = dam_obj_sub
)

# -------------------------------------------------------------------
# Build adult reference in Dam modality and map P8 cells
# -------------------------------------------------------------------
# Adult reference = W8 + W8R
# Query = P8
# The mapped representation is stored in:
#   - ref.umap   (query projected into reference UMAP)
#   - ref.lsi    (query projected into reference LSI)
# These are then used to construct a joint reference-space embedding.

ref_stages  <- c("W8", "W8R")
query_stage <- "P8"

ref   <- subset(combined.obj.ls.2$Dam_bb, subset = stage %in% ref_stages)
query <- subset(combined.obj.ls$Dam_bb,   subset = stage == query_stage)

DefaultAssay(ref)   <- "peaks"
DefaultAssay(query) <- "peaks"

# Create stage-suffixed adult labels in the reference
ref$K27ac_Annotation_cell_type_Adult <- paste0(
  as.character(ref$K27ac_Annotation_cell_type),
  "_",
  as.character(ref$stage)
)

ref <- subset(ref, subset = !is.na(K27ac_Annotation_cell_type_Adult))

# Recompute TFIDF / LSI after subsetting
ref   <- RunTFIDF(ref)
ref   <- FindTopFeatures(ref, min.cutoff = "q0")
ref   <- RunSVD(ref)

query <- RunTFIDF(query)
query <- FindTopFeatures(query, min.cutoff = "q0")
query <- RunSVD(query)

# Reference UMAP
dims_use <- 2:14
ref <- RunUMAP(ref, reduction = "lsi", dims = dims_use, return.model = TRUE)

# Anchors and mapping
anchors <- FindTransferAnchors(
  reference = ref,
  query = query,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = dims_use
)

ref_vec <- ref$K27ac_Annotation_cell_type_Adult
names(ref_vec) <- colnames(ref)

query_mapped <- MapQuery(
  anchorset = anchors,
  reference = ref,
  query = query,
  refdata = ref_vec,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = "umap"
)

# Score column (Seurat-version dependent)
md_cols <- colnames(query_mapped@meta.data)
score_col <- if ("predicted.id.score" %in% md_cols) {
  "predicted.id.score"
} else if ("prediction.score.max" %in% md_cols) {
  "prediction.score.max"
} else {
  stop("Score column not found. meta.data columns: ", paste(md_cols, collapse = ", "))
}

# Threshold mapped adult labels
thr <- 0.38
query_mapped$Adult_pred  <- query_mapped$predicted.id
query_mapped$Adult_score <- query_mapped[[score_col]][, 1]

query_mapped$K27ac_Annotation_cell_type_Adult <- ifelse(
  query_mapped$Adult_score >= thr,
  as.character(query_mapped$Adult_pred),
  "Unknown"
)

# Minimal factor ordering: Unknown last
labels_present <- unique(c(
  as.character(ref$K27ac_Annotation_cell_type_Adult),
  as.character(query_mapped$K27ac_Annotation_cell_type_Adult)
))
labels_present <- labels_present[!is.na(labels_present)]
ordered_levels <- c(setdiff(labels_present, "Unknown"), "Unknown")

ref$K27ac_Annotation_cell_type_Adult <- factor(
  as.character(ref$K27ac_Annotation_cell_type_Adult),
  levels = ordered_levels
)
query_mapped$K27ac_Annotation_cell_type_Adult <- factor(
  as.character(query_mapped$K27ac_Annotation_cell_type_Adult),
  levels = ordered_levels
)

# Plot reference and mapped query
if (!("ref.umap" %in% names(query_mapped@reductions))) {
  stop("query_mapped does not contain reduction 'ref.umap'.")
}

emb_refumap <- Embeddings(query_mapped, "ref.umap")
keep_cells  <- rownames(emb_refumap)[complete.cases(emb_refumap)]
drop_cells  <- setdiff(colnames(query_mapped), keep_cells)

query_mapped$ref_umap_ok <- colnames(query_mapped) %in% keep_cells
query_plot <- subset(query_mapped, cells = keep_cells)

p_ref <- DimPlot(
  ref,
  reduction = "umap",
  group.by = "K27ac_Annotation_cell_type_Adult",
  label = TRUE,
  repel = TRUE
) + NoLegend() +
  ggtitle("Reference (W8 + W8R; adult labels)")

p_query <- DimPlot(
  query_plot,
  reduction = "ref.umap",
  group.by = "K27ac_Annotation_cell_type_Adult",
  label = TRUE,
  repel = TRUE
) + NoLegend() +
  ggtitle(
    paste0(
      "Query (P8) mapped; thr=", thr,
      " | shown=", ncol(query_plot),
      " | dropped(no ref.umap)=", length(drop_cells)
    )
  )

p_ref | p_query

# Split reference into W8 vs W8R for visualization
ref_W8  <- subset(ref, subset = stage == "W8")
ref_W8R <- subset(ref, subset = stage == "W8R")

p_ref_W8 <- DimPlot(
  ref_W8,
  reduction = "umap",
  group.by = "K27ac_Annotation_cell_type_Adult",
  label = TRUE,
  repel = TRUE
) + NoLegend() + ggtitle("Reference: W8")

p_ref_W8R <- DimPlot(
  ref_W8R,
  reduction = "umap",
  group.by = "K27ac_Annotation_cell_type_Adult",
  label = TRUE,
  repel = TRUE
) + NoLegend() + ggtitle("Reference: W8R")

p_query2 <- DimPlot(
  query_plot,
  reduction = "ref.umap",
  group.by = "K27ac_Annotation_cell_type_Adult",
  label = TRUE,
  repel = TRUE
) + NoLegend() +
  ggtitle(
    paste0(
      "Query: P8 mapped (thr=", thr,
      ") | shown=", ncol(query_plot),
      " | dropped=", length(drop_cells)
    )
  )

(p_ref_W8 | p_ref_W8R) / p_query2

# -------------------------------------------------------------------
# Construct merged object in reference space
# -------------------------------------------------------------------
# joint.umap:
#   - reference uses ref@umap
#   - query uses query_mapped@ref.umap
#
# joint.ref.lsi:
#   - reference uses ref@lsi
#   - query uses query_mapped@ref.lsi
#
# This creates a shared reference-space low-dimensional representation suitable
# for downstream quantitative analyses.

stopifnot(
  "K27ac_Annotation_cell_type_Adult" %in% colnames(ref@meta.data),
  "K27ac_Annotation_cell_type_Adult" %in% colnames(query_mapped@meta.data),
  "umap"     %in% names(ref@reductions),
  "lsi"      %in% names(ref@reductions),
  "ref.umap" %in% names(query_mapped@reductions),
  "ref.lsi"  %in% names(query_mapped@reductions)
)

overlap <- intersect(colnames(ref), colnames(query_mapped))
if (length(overlap) > 0) {
  merged_all <- merge(ref, y = query_mapped, add.cell.ids = c("REF", "P8"))
  ref_cells_new <- paste0("REF_", colnames(ref))
  qry_cells_new <- paste0("P8_",  colnames(query_mapped))
} else {
  merged_all <- merge(ref, y = query_mapped)
  ref_cells_new <- colnames(ref)
  qry_cells_new <- colnames(query_mapped)
}

merged_all$Source <- NA_character_
merged_all$Source[ref_cells_new] <- "Reference"
merged_all$Source[qry_cells_new] <- "Query"
merged_all$Source <- factor(merged_all$Source, levels = c("Reference", "Query"))

# joint.umap
emb_joint_umap <- matrix(
  NA_real_,
  nrow = ncol(merged_all),
  ncol = 2,
  dimnames = list(colnames(merged_all), c("UMAP_1", "UMAP_2"))
)

emb_joint_umap[ref_cells_new, ] <- Embeddings(ref, "umap")[, 1:2, drop = FALSE]
emb_joint_umap[qry_cells_new, ] <- Embeddings(query_mapped, "ref.umap")[, 1:2, drop = FALSE]

merged_all[["joint.umap"]] <- CreateDimReducObject(
  embeddings = emb_joint_umap,
  key = "JUMAP_",
  assay = DefaultAssay(ref)
)

# joint.ref.lsi
lsi_ref <- Embeddings(ref, "lsi")
lsi_qry <- Embeddings(query_mapped, "ref.lsi")

k_ref <- ncol(lsi_ref)
k_qry <- ncol(lsi_qry)
k_req <- length(dims_use)
k_use <- min(k_qry, k_req)

if (k_use < k_req) {
  message("NOTE: query_mapped@ref.lsi has only ", k_qry,
          " dims; requested dims_use length is ", k_req,
          ". Using first ", k_use, " dims from dims_use.")
}

dims_use_ref <- dims_use[seq_len(k_use)]

if (max(dims_use_ref) > k_ref) {
  stop("ref@lsi has only ", k_ref, " dims, but dims_use_ref requires up to ", max(dims_use_ref), ".")
}

lsi_ref_use <- lsi_ref[, dims_use_ref, drop = FALSE]
lsi_qry_use <- lsi_qry[, seq_len(k_use), drop = FALSE]

emb_joint_lsi <- matrix(
  NA_real_,
  nrow = ncol(merged_all),
  ncol = k_use,
  dimnames = list(colnames(merged_all), paste0("LSI_", dims_use_ref))
)

emb_joint_lsi[ref_cells_new, ] <- lsi_ref_use
emb_joint_lsi[qry_cells_new, ] <- lsi_qry_use

merged_all[["joint.ref.lsi"]] <- CreateDimReducObject(
  embeddings = emb_joint_lsi,
  key = "JLSI_",
  assay = DefaultAssay(ref)
)
merged_all[["joint.ref.lsi"]]@misc$dims_use_ref <- dims_use_ref

# Create simplified Cell_Type by removing adult-stage suffix
merged_all$Cell_Type <- sub(
  "_(W8R|W8)$",
  "",
  as.character(merged_all$K27ac_Annotation_cell_type_Adult)
)

ct_levels <- unique(merged_all$Cell_Type)
ct_levels <- ct_levels[!is.na(ct_levels)]
ct_levels <- c(setdiff(ct_levels, "Unknown"), "Unknown")
merged_all$Cell_Type <- factor(merged_all$Cell_Type, levels = ct_levels)

# Remove P8 cells mapped as Unknown
cells_keep <- colnames(merged_all)[
  !(merged_all$Source == "Query" &
      merged_all$K27ac_Annotation_cell_type_Adult == "Unknown")
]
merged_all <- subset(merged_all, cells = cells_keep)

p1 <- DimPlot(
  merged_all,
  reduction = "joint.umap",
  group.by = "K27ac_Annotation_cell_type_Adult",
  label = FALSE,
  shuffle = TRUE
)

p2 <- DimPlot(
  merged_all,
  reduction = "joint.umap",
  group.by = "Cell_Type",
  label = FALSE,
  shuffle = TRUE
)

p3 <- DimPlot(
  merged_all,
  reduction = "joint.umap",
  group.by = "stage",
  label = FALSE,
  shuffle = TRUE
)

p4 <- DimPlot(
  merged_all,
  reduction = "joint.umap",
  group.by = "orig.ident",
  label = FALSE,
  shuffle = TRUE
)

p1 + p2 + p3 + p4

# Optional layered stage plot
set.seed(123)
p_base <- DimPlot(merged_all, reduction = "joint.umap", group.by = "stage",
                  label = FALSE, raster = FALSE)

df <- p_base$data
df_w8  <- df[df$stage == "W8", , drop = FALSE]
df_top <- df[df$stage %in% c("P8", "W8R"), , drop = FALSE]
df_top <- df_top[sample(nrow(df_top)), , drop = FALSE]
df_other <- df[!(df$stage %in% c("W8", "P8", "W8R")), , drop = FALSE]

p_layered <- p_base
p_layered$data <- rbind(df_w8, df_top, df_other)

p1 + p2 + p_layered + p4

# -------------------------------------------------------------------
# Transfer adult labels and Cell_Type to K27ac
# -------------------------------------------------------------------
# K27ac keeps only cells present in merged_all and receives the adult labels /
# Cell_Type metadata from the Dam-based reference-mapping result.

k27ac_full <- combined.obj.ls$K27ac_aa
dam_full   <- merged_all

cells_k27 <- colnames(k27ac_full)
cells_dam <- colnames(dam_full)

missing_in_k27ac <- setdiff(cells_dam, cells_k27)
if (length(missing_in_k27ac) > 0) {
  stop(
    "Some cells in merged_all are NOT found in combined.obj.ls$K27ac_aa.\n",
    "n_missing = ", length(missing_in_k27ac), "\n",
    "Examples: ", paste(head(missing_in_k27ac, 20), collapse = ", ")
  )
}

cells_use <- cells_k27[cells_k27 %in% cells_dam]

stopifnot(setequal(cells_use, cells_dam))
stopifnot(!anyDuplicated(cells_use))

K27ac_aa_sub <- subset(k27ac_full, cells = cells_use)
Dam_bb_sub   <- subset(dam_full,   cells = cells_use)

# Drop prediction-related assays from Dam object if present
pred_assays <- grep("^prediction|^predictions", Assays(Dam_bb_sub), value = TRUE)
if (length(pred_assays) > 0) {
  if (DefaultAssay(Dam_bb_sub) %in% pred_assays) {
    DefaultAssay(Dam_bb_sub) <- "peaks"
  }
  for (a in pred_assays) {
    Dam_bb_sub[[a]] <- NULL
  }
}

# Enforce Dam order to match K27ac order
Dam_bb_sub@meta.data <- Dam_bb_sub@meta.data[cells_use, , drop = FALSE]
Dam_bb_sub@active.ident <- Dam_bb_sub@active.ident[cells_use]
names(Dam_bb_sub@active.ident) <- cells_use

for (assay_nm in names(Dam_bb_sub@assays)) {
  a <- Dam_bb_sub@assays[[assay_nm]]
  
  if ("layers" %in% slotNames(a)) {
    for (ly in names(a@layers)) {
      m <- a@layers[[ly]]
      if (is.null(m) || ncol(m) == 0) next
      a@layers[[ly]] <- m[, cells_use, drop = FALSE]
    }
    if ("cells" %in% slotNames(a) && !is.null(a@cells)) {
      if (all(cells_use %in% rownames(a@cells))) {
        a@cells <- a@cells[cells_use, , drop = FALSE]
      }
    }
  } else {
    for (sl in intersect(c("counts", "data", "scale.data"), slotNames(a))) {
      m <- slot(a, sl)
      if (is.null(m) || ncol(m) == 0) next
      slot(a, sl) <- m[, cells_use, drop = FALSE]
    }
  }
  
  Dam_bb_sub@assays[[assay_nm]] <- a
}

if (length(Dam_bb_sub@reductions) > 0) {
  for (red_nm in names(Dam_bb_sub@reductions)) {
    emb <- Dam_bb_sub@reductions[[red_nm]]@cell.embeddings
    if (is.null(emb) || nrow(emb) == 0) next
    
    emb2 <- matrix(
      NA_real_,
      nrow = length(cells_use),
      ncol = ncol(emb),
      dimnames = list(cells_use, colnames(emb))
    )
    common <- intersect(cells_use, rownames(emb))
    emb2[common, ] <- emb[common, , drop = FALSE]
    Dam_bb_sub@reductions[[red_nm]]@cell.embeddings <- emb2
  }
}

K27ac_aa_sub@meta.data <- K27ac_aa_sub@meta.data[cells_use, , drop = FALSE]
K27ac_aa_sub@active.ident <- K27ac_aa_sub@active.ident[cells_use]
names(K27ac_aa_sub@active.ident) <- cells_use

stopifnot(identical(colnames(K27ac_aa_sub), cells_use))
stopifnot(identical(colnames(Dam_bb_sub),   cells_use))
stopifnot(identical(colnames(K27ac_aa_sub), colnames(Dam_bb_sub)))

# Transfer metadata from Dam to K27ac
need_cols <- c("K27ac_Annotation_cell_type_Adult", "Cell_Type")
missing_cols <- setdiff(need_cols, colnames(Dam_bb_sub@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing columns in Dam_bb_sub meta.data: ", paste(missing_cols, collapse = ", "))
}

v_adult <- Dam_bb_sub$K27ac_Annotation_cell_type_Adult
names(v_adult) <- colnames(Dam_bb_sub)
K27ac_aa_sub <- AddMetaData(
  K27ac_aa_sub,
  metadata = v_adult,
  col.name = "K27ac_Annotation_cell_type_Adult"
)

v_ct <- Dam_bb_sub$Cell_Type
names(v_ct) <- colnames(Dam_bb_sub)
K27ac_aa_sub <- AddMetaData(
  K27ac_aa_sub,
  metadata = v_ct,
  col.name = "Cell_Type"
)

combined.obj.ls3 <- list(
  K27ac_aa = K27ac_aa_sub,
  Dam_bb   = Dam_bb_sub
)

# -------------------------------------------------------------------
# Re-normalize and re-cluster K27ac for visualization
# -------------------------------------------------------------------
# The K27ac UMAP below is a supervised / visualization-oriented embedding.
# By contrast, downstream quantitative analyses should continue to use the
# unbiased LSI space, especially the reference-space representation
# stored in joint.ref.lsi.

DefaultAssay(combined.obj.ls3$K27ac_aa) <- "peaks"
obj <- combined.obj.ls3$K27ac_aa
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj)
obj <- RunSVD(obj)

combined.obj.ls3$K27ac_aa <- obj

combined.obj.ls3$K27ac_aa <- RunUMAP(combined.obj.ls3$K27ac_aa, reduction = "lsi", dims = 2:10)
combined.obj.ls3$K27ac_aa <- FindNeighbors(combined.obj.ls3$K27ac_aa, reduction = "lsi", dims = 2:15)
combined.obj.ls3$K27ac_aa <- FindClusters(combined.obj.ls3$K27ac_aa, verbose = FALSE, algorithm = 3, resolution = 2)

# Cell-type supervised UMAP for K27ac visualization
obj <- combined.obj.ls3$K27ac_aa
DefaultAssay(obj) <- "peaks"

dims_use_vis <- 2:8
lsi_mat <- Embeddings(obj, "lsi")[, dims_use_vis, drop = FALSE]

y <- as.character(obj$Cell_Type)
y[y %in% c("Unknown", NA)] <- NA
y <- factor(y)

set.seed(1)
umap_fit <- uwot::umap(
  X = lsi_mat,
  metric = "cosine",
  n_neighbors = 30,
  min_dist = 1,
  y = y,
  target_weight = 0.6,
  target_n_neighbors = 30,
  negative_sample_rate = 3,
  ret_model = TRUE,
  verbose = TRUE
)

emb <- umap_fit$embedding
rownames(emb) <- rownames(lsi_mat)
colnames(emb) <- paste0("CTUMAP_", seq_len(ncol(emb)))

obj[["umap.celltype"]] <- CreateDimReducObject(
  embeddings = emb,
  key = "CTUMAP_",
  assay = DefaultAssay(obj)
)
obj[["umap.celltype"]]@misc$model <- umap_fit$model
combined.obj.ls3$K27ac_aa <- obj

# -------------------------------------------------------------------
# Harmonize UMAP naming across modalities for plotting
# -------------------------------------------------------------------
# Dam modality already uses "joint.umap".
# For K27ac, rename the supervised visualization UMAP to "joint.umap"
# so that the same plotting code can be reused across modalities.
#
# This renaming is for visualization consistency only and does not change the
# LSI / joint.ref.lsi used for downstream quantitative analyses.

if ("umap.celltype" %in% Reductions(combined.obj.ls3$K27ac_aa)) {
  combined.obj.ls3$K27ac_aa[["joint.umap"]] <- combined.obj.ls3$K27ac_aa[["umap.celltype"]]
}
if ("umap" %in% Reductions(combined.obj.ls3$K27ac_aa)) {
  combined.obj.ls3$K27ac_aa[["umap"]] <- NULL
}
if ("umap.celltype" %in% Reductions(combined.obj.ls3$K27ac_aa)) {
  combined.obj.ls3$K27ac_aa[["umap.celltype"]] <- NULL
}

# -------------------------------------------------------------------
# Final visualization and save
# -------------------------------------------------------------------
p1 <- DimPlot(combined.obj.ls3$K27ac_aa, reduction = "joint.umap",
              group.by = "Cell_Type", label = FALSE, shuffle = TRUE)
p2 <- DimPlot(combined.obj.ls3$Dam_bb, reduction = "joint.umap",
              group.by = "Cell_Type", label = FALSE, shuffle = TRUE)
p3 <- DimPlot(combined.obj.ls3$K27ac_aa, reduction = "joint.umap",
              group.by = "stage", label = FALSE, shuffle = TRUE)
p4 <- DimPlot(combined.obj.ls3$Dam_bb, reduction = "joint.umap",
              group.by = "stage", label = FALSE, shuffle = TRUE)
p1 + p2 + p3 + p4

# -------------------------------------------------------------------
# Save final object
# -------------------------------------------------------------------
saveRDS(combined.obj.ls3, file = "nanoscope_final_peaks_nanoscopefilter_in_vivo_Cell_Type.rds")
