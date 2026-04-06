# ============================================================
# W8R Dam: L2-3 vs L4-5 ambiguity and commitment strength
#
# Input object:
#   nanoscope_final_peaks_nanoscopefilter_in_vivo_Cell_Type.rds
#
# Assumptions:
# - object contains:
#     Cell_Type (L2-3, L4-5)
#     joint.ref.lsi (Dam modality)
#     joint.umap
# - analysis is restricted to W8R Dam cells (i.e., Dam modality of P8-Wk8 HisTrac-seq sample)
#
# Note on joint.ref.lsi:
# - Dam medoids are computed in joint.ref.lsi, a shared reference-space LSI in which
#   W8/W8R reference cells retain their original Dam LSI coordinates and P8 cells are
#   projected onto that adult reference.
# - Therefore, joint.ref.lsi uses a re-indexed reference-space representation rather
#   than the original raw LSI numbering.
# - For this reason, Dam analyses based on joint.ref.lsi appropriately start from
#   dimension 1, even though the original Dam LSI workflow typically used dimensions
#   beginning from 2 after excluding the first raw LSI component.
# ============================================================

library(Seurat)
library(FNN)
library(glmnet)

# -------------------------------------------------------------------
# Input
# -------------------------------------------------------------------
combined.obj.ls3 <- readRDS("nanoscope_final_peaks_nanoscopefilter_in_vivo_Cell_Type.rds")
obj <- subset(combined.obj.ls3$Dam_bb, subset = stage == "W8R")

# -------------------------------------------------------------------
# Restrict to L2-3 and L4-5
# -------------------------------------------------------------------
obj <- subset(obj, subset = Cell_Type %in% c("L2-3", "L4-5"))
obj$Cell_Type <- factor(obj$Cell_Type, levels = c("L2-3", "L4-5"))

# -------------------------------------------------------------------
# Use Dam reference-space LSI
# -------------------------------------------------------------------
X <- Embeddings(obj, reduction = "joint.ref.lsi")

# -------------------------------------------------------------------
# 1) kNN ambiguity
# - ambiguity = 1 - purity
# - higher values indicate more mixed local neighborhoods
# -------------------------------------------------------------------
lab <- as.character(obj$Cell_Type)
k <- 30

nn <- get.knn(X, k = k)$nn.index

knn_purity <- sapply(seq_len(nrow(nn)), function(i) {
  mean(lab[nn[i, ]] == lab[i])
})

obj$knn_ambiguity <- 1 - knn_purity

p_ambiguity <- FeaturePlot(
  obj,
  reduction = "joint.umap",
  features = "knn_ambiguity",
  pt.size = 2,
) + ggtitle("kNN ambiguity (= 1 - purity)")

print(p_ambiguity)

# -------------------------------------------------------------------
# 2) Commitment strength
# - fit a classifier in joint.ref.lsi space
# - probability of being L4-5 is converted to logit distance from 0.5
# - absolute logit = commitment strength
# -------------------------------------------------------------------
y <- ifelse(obj$Cell_Type == "L4-5", 1, 0)

set.seed(1)
cv <- cv.glmnet(
  x = as.matrix(X),
  y = y,
  family = "binomial",
  alpha = 1,
  nfolds = 5
)

p_L45 <- as.numeric(
  predict(cv, newx = as.matrix(X), s = "lambda.1se", type = "response")
)

obj$p_L45 <- p_L45
obj$id_logit <- qlogis(p_L45)
obj$commitment_strength <- abs(obj$id_logit)

p_commitment <- FeaturePlot(
  obj,
  reduction = "joint.umap",
  features = "commitment_strength",
  pt.size = 2,
) + ggtitle("Commitment strength")

print(p_commitment)




# ============================================================
# DEG analysis
# Continuation:
# W8R (HisTrac-seq) Dam L2-3 vs L4-5
# - Commitment50
# - pseudobulk edgeR (L2-3 vs L4-5)
# - commitment-stratified gene-set split
# - final gene-level heatmap
# ============================================================

library(Matrix)
library(edgeR)
library(limma)
library(dplyr)
library(pheatmap)

stopifnot("Cell_Type" %in% colnames(obj@meta.data))
stopifnot("p_L45" %in% colnames(obj@meta.data))

# ----------------------------
# 1) Build agreement and commitment groups
# ----------------------------
# For pseudobulk DEG analysis, we further restrict to cells with dam_agree == TRUE,
# i.e. cells whose classifier-based lineage prediction matches their annotated cell type.
# This restriction was applied to reduce contamination from highly ambiguous cells
# near the L2-3/L4-5 boundary and to define a more stringent set of cells for
# lineage-specific differential-expression analysis.

obj$dam_pred <- ifelse(obj$p_L45 >= 0.5, "L4-5", "L2-3")
obj$dam_agree <- obj$dam_pred == as.character(obj$Cell_Type)

obj$dam_support_true <- ifelse(
  obj$Cell_Type == "L4-5",
  obj$p_L45,
  1 - obj$p_L45
)

med_by_type <- tapply(obj$dam_support_true, obj$Cell_Type, median, na.rm = TRUE)

obj$DamCommit50 <- ifelse(
  obj$dam_support_true >= med_by_type[as.character(obj$Cell_Type)],
  "High",
  "Low"
)
obj$DamCommit50 <- factor(obj$DamCommit50, levels = c("Low", "High"))

# ----------------------------
# 2) Pseudobulk L2-3 vs L4-5 (agreement-only)
#    -> Pseudobulk_L2-3_vs_L4-5_edgeR.csv
# ----------------------------
# W8Rrep3 is excluded from pseudobulk analyses, consistent with the final analysis
# used in the manuscript.
assay_use <- "RNA.Dam"
exclude_reps <- c("W8Rrep3")   # set character(0) if not needed

DefaultAssay(obj) <- assay_use

obj_pb0 <- subset(
  obj,
  subset = (dam_agree == TRUE) & (Cell_Type %in% c("L2-3", "L4-5"))
)

if (length(exclude_reps) > 0 && "orig.ident" %in% colnames(obj_pb0@meta.data)) {
  obj_pb0 <- subset(obj_pb0, subset = !(orig.ident %in% exclude_reps))
}

obj_pb0$Cell_Type <- factor(obj_pb0$Cell_Type, levels = c("L2-3", "L4-5"))

rep_col <- if ("sample" %in% colnames(obj_pb0@meta.data)) "sample" else "orig.ident"

cts0 <- tryCatch(
  GetAssayData(obj_pb0, assay = assay_use, layer = "counts"),
  error = function(e) GetAssayData(obj_pb0, assay = assay_use, slot = "counts")
)

meta0 <- obj_pb0@meta.data
meta0$rep <- meta0[[rep_col]]
meta0$ct  <- as.character(meta0$Cell_Type)
meta0$pb_id <- paste(meta0$rep, meta0$ct, sep = "__")
pb_f0 <- factor(meta0$pb_id)

mm0 <- Matrix::sparse.model.matrix(~0 + pb_f0)
colnames(mm0) <- levels(pb_f0)

pb_counts0 <- as.matrix(cts0 %*% mm0)

meta_pb0 <- data.frame(
  pb_id = levels(pb_f0),
  stringsAsFactors = FALSE
)
meta_pb0$rep <- sub("__.*$", "", meta_pb0$pb_id)
meta_pb0$Cell_Type <- sub("^.*__", "", meta_pb0$pb_id)
meta_pb0$n_cells <- as.integer(table(pb_f0)[meta_pb0$pb_id])
rownames(meta_pb0) <- meta_pb0$pb_id

meta_pb0$ct2 <- ifelse(meta_pb0$Cell_Type == "L2-3", "L2_3", "L4_5")
meta_pb0$ct2 <- factor(meta_pb0$ct2, levels = c("L2_3", "L4_5"))
meta_pb0$rep <- factor(meta_pb0$rep)

y0 <- DGEList(counts = pb_counts0, samples = meta_pb0)
design0 <- model.matrix(~ ct2 + rep, data = meta_pb0)

keep0 <- filterByExpr(y0, design = design0)
y0 <- y0[keep0, , keep.lib.sizes = FALSE]

y0 <- calcNormFactors(y0)
y0 <- estimateDisp(y0, design0)
fit0 <- glmQLFit(y0, design0)

coef_name <- "ct2L4_5"
stopifnot(coef_name %in% colnames(design0))

# positive logFC = L2-3 > L4-5
contr0 <- limma::makeContrasts(L2_3_minus_L4_5 = -ct2L4_5, levels = design0)
qlf0 <- glmQLFTest(fit0, contrast = contr0)

res0 <- topTags(qlf0, n = Inf)$table
res0$gene <- rownames(res0)
res0 <- res0[order(res0$FDR), ]

write.csv(
  res0,
  "Pseudobulk_L2-3_vs_L4-5_edgeR.csv",
  row.names = FALSE,
  quote = FALSE
)

L23_edgeR <- res0$gene[res0$FDR < 0.05 & res0$logFC > 0]
L45_edgeR <- res0$gene[res0$FDR < 0.05 & res0$logFC < 0]

# ----------------------------
# 3) Build 4-group pseudobulk (L23/L45 x High/Low)
# ----------------------------
group_levels <- c("L23_High", "L23_Low", "L45_High", "L45_Low")

obj_use <- subset(
  obj,
  subset = (dam_agree == TRUE) &
    (Cell_Type %in% c("L2-3", "L4-5")) &
    (DamCommit50 %in% c("Low", "High"))
)

if (length(exclude_reps) > 0 && "orig.ident" %in% colnames(obj_use@meta.data)) {
  obj_use <- subset(obj_use, subset = !(orig.ident %in% exclude_reps))
}

obj_use$Group4 <- NA_character_
obj_use$Group4[obj_use$Cell_Type == "L2-3" & obj_use$DamCommit50 == "High"] <- "L23_High"
obj_use$Group4[obj_use$Cell_Type == "L2-3" & obj_use$DamCommit50 == "Low"]  <- "L23_Low"
obj_use$Group4[obj_use$Cell_Type == "L4-5" & obj_use$DamCommit50 == "High"] <- "L45_High"
obj_use$Group4[obj_use$Cell_Type == "L4-5" & obj_use$DamCommit50 == "Low"]  <- "L45_Low"
obj_use$Group4 <- factor(obj_use$Group4, levels = group_levels)

obj_use <- subset(obj_use, subset = !is.na(Group4))

rep_col <- if ("sample" %in% colnames(obj_use@meta.data)) "sample" else "orig.ident"

cts <- tryCatch(
  GetAssayData(obj_use, assay = assay_use, layer = "counts"),
  error = function(e) GetAssayData(obj_use, assay = assay_use, slot = "counts")
)

meta <- obj_use@meta.data
meta$rep <- meta[[rep_col]]
meta$group4 <- as.character(obj_use$Group4)
meta$pb_id <- paste(meta$rep, meta$group4, sep = "__")
pb_f <- factor(meta$pb_id)

mm <- Matrix::sparse.model.matrix(~0 + pb_f)
colnames(mm) <- levels(pb_f)

pb_counts <- as.matrix(cts %*% mm)

meta_pb <- data.frame(
  pb_id = colnames(pb_counts),
  stringsAsFactors = FALSE
)
meta_pb$rep <- sub("__.*$", "", meta_pb$pb_id)
meta_pb$group4 <- sub("^.*__", "", meta_pb$pb_id)
meta_pb$group4 <- factor(meta_pb$group4, levels = group_levels)
meta_pb$n_cells <- as.integer(table(pb_f)[meta_pb$pb_id])
rownames(meta_pb) <- meta_pb$pb_id
meta_pb$rep <- factor(meta_pb$rep)

# ----------------------------
# 4) Rep-adjusted means for heatmap
# ----------------------------
y_logcpm <- DGEList(counts = pb_counts)
y_logcpm <- calcNormFactors(y_logcpm)
logCPM <- cpm(y_logcpm, log = TRUE, prior.count = 2)

get_keep_reps <- function(meta_pb, low_group, high_group) {
  reps_low  <- unique(as.character(meta_pb$rep[meta_pb$group4 == low_group]))
  reps_high <- unique(as.character(meta_pb$rep[meta_pb$group4 == high_group]))
  intersect(reps_low, reps_high)
}

keep_L23 <- get_keep_reps(meta_pb, "L23_Low", "L23_High")
keep_L45 <- get_keep_reps(meta_pb, "L45_Low", "L45_High")

calc_rep_adjusted_means <- function(logCPM, meta_pb, groups, keep_reps) {
  cols <- rownames(meta_pb)[meta_pb$group4 %in% groups & meta_pb$rep %in% keep_reps]
  sub_meta <- droplevels(meta_pb[cols, , drop = FALSE])
  sub_log  <- logCPM[, cols, drop = FALSE]
  
  design_g <- model.matrix(~0 + group4, data = sub_meta)
  sub_log_adj <- removeBatchEffect(sub_log, batch = sub_meta$rep, design = design_g)
  
  out <- sapply(groups, function(g) {
    rowMeans(sub_log_adj[, sub_meta$group4 == g, drop = FALSE])
  })
  rownames(out) <- rownames(sub_log_adj)
  out
}

expr23 <- calc_rep_adjusted_means(
  logCPM,
  meta_pb,
  groups = c("L23_High", "L23_Low"),
  keep_reps = keep_L23
)

expr45 <- calc_rep_adjusted_means(
  logCPM,
  meta_pb,
  groups = c("L45_High", "L45_Low"),
  keep_reps = keep_L45
)

expr4_edgeR <- cbind(
  expr23[, c("L23_High", "L23_Low"), drop = FALSE],
  expr45[, c("L45_High", "L45_Low"), drop = FALSE]
)

expr4_edgeR <- expr4_edgeR[, group_levels, drop = FALSE]
expr4_edgeR_z <- t(scale(t(expr4_edgeR)))
expr4_edgeR_z[is.na(expr4_edgeR_z)] <- 0

# ----------------------------
# 5) Within-cell-type edgeR: High vs Low
# ----------------------------
run_commit_edger <- function(pb_counts, meta_pb, high_group, low_group) {
  cols <- rownames(meta_pb)[meta_pb$group4 %in% c(low_group, high_group)]
  sub_counts <- pb_counts[, cols, drop = FALSE]
  sub_meta <- meta_pb[cols, , drop = FALSE]
  
  sub_meta$commit <- ifelse(sub_meta$group4 == high_group, "High", "Low")
  sub_meta$commit <- factor(sub_meta$commit, levels = c("Low", "High"))
  
  tab <- table(sub_meta$rep, sub_meta$commit)
  keep_reps <- rownames(tab)[apply(tab, 1, function(x) all(x > 0))]
  sub_keep <- sub_meta$rep %in% keep_reps
  
  sub_counts <- sub_counts[, sub_keep, drop = FALSE]
  sub_meta   <- droplevels(sub_meta[sub_keep, , drop = FALSE])
  
  y <- DGEList(counts = sub_counts)
  design <- model.matrix(~ commit + rep, data = sub_meta)
  
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  
  qlf <- glmQLFTest(fit, coef = "commitHigh")
  res <- topTags(qlf, n = Inf)$table
  res$gene <- rownames(res)
  res
}

res_L23_commit <- run_commit_edger(pb_counts, meta_pb, "L23_High", "L23_Low")
res_L45_commit <- run_commit_edger(pb_counts, meta_pb, "L45_High", "L45_Low")

lfc_L23 <- res_L23_commit$logFC
names(lfc_L23) <- res_L23_commit$gene

lfc_L45 <- res_L45_commit$logFC
names(lfc_L45) <- res_L45_commit$gene

# ----------------------------
# 6) Split L23_edgeR / L45_edgeR into early / consistent / late
# ----------------------------
lfc_thr <- 0.1

split_early_consistent_late <- function(genes, lfc_vec, thr = 0.1) {
  genes <- unique(genes)
  lfc <- lfc_vec[genes]
  ok  <- !is.na(lfc)
  
  list(
    early      = genes[ok & (lfc < -thr)],
    consistent = genes[ok & (abs(lfc) <= thr)],
    late       = genes[ok & (lfc > thr)],
    missing    = genes[!ok]
  )
}

sp23 <- split_early_consistent_late(L23_edgeR, lfc_L23, thr = lfc_thr)
sp45 <- split_early_consistent_late(L45_edgeR, lfc_L45, thr = lfc_thr)

gene_sets <- list(
  L23_early_edgeR      = sp23$early,
  L23_consistent_edgeR = sp23$consistent,
  L23_late_edgeR       = sp23$late,
  L45_early_edgeR      = sp45$early,
  L45_consistent_edgeR = sp45$consistent,
  L45_late_edgeR       = sp45$late
)

gs_df <- bind_rows(lapply(names(gene_sets), function(nm) {
  data.frame(
    gene = gene_sets[[nm]],
    gene_set = nm,
    stringsAsFactors = FALSE
  )
}))

write.csv(
  gs_df,
  "EdgeR_geneSets_6groups_EarlyConsistentLate.csv",
  row.names = FALSE,
  quote = FALSE
)

# ----------------------------
# 7) Final gene-level heatmap
#    -> Heatmap_Genes_Zscore_4groups_6blocks_EdgeRbased.pdf
# ----------------------------
set_order <- c(
  "L23_early_edgeR",
  "L23_consistent_edgeR",
  "L23_late_edgeR",
  "L45_early_edgeR",
  "L45_consistent_edgeR",
  "L45_late_edgeR"
)

genes_ordered <- unlist(gene_sets[set_order], use.names = FALSE)
genes_ordered <- unique(genes_ordered)
genes_ordered <- intersect(genes_ordered, rownames(expr4_edgeR_z))

stopifnot(length(genes_ordered) > 0)

mat_z <- expr4_edgeR_z[genes_ordered, group_levels, drop = FALSE]

gene2set <- rep(NA_character_, length(genes_ordered))
names(gene2set) <- genes_ordered
for (nm in set_order) {
  g <- intersect(gene_sets[[nm]], genes_ordered)
  gene2set[g] <- nm
}

ann_row <- data.frame(
  GeneSet = factor(gene2set[genes_ordered], levels = set_order)
)
rownames(ann_row) <- genes_ordered

sizes <- sapply(gene_sets[set_order], function(g) length(intersect(g, genes_ordered)))
gaps_row <- cumsum(sizes)
gaps_row <- gaps_row[gaps_row > 0 & gaps_row < nrow(mat_z)]

max_abs <- as.numeric(quantile(abs(mat_z), 0.99, na.rm = TRUE))
bk <- seq(-max_abs, max_abs, length.out = 101)
cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

ann_colors <- list(
  GeneSet = setNames(
    c("#E41A1C", "#FB9A99", "#A50F15", "#1F78B4", "#00A6D6", "#6A3D9A"),
    set_order
  )
)

pdf_h <- max(7, 0.12 * nrow(mat_z) + 2)

ph <- pheatmap::pheatmap(
  mat_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 5,
  fontsize_col = 11,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  gaps_row = gaps_row,
  color = cols,
  breaks = bk,
  main = "Genes (rows) x 4 groups (cols) | gene-wise Z-score | 6 blocks (early/consistent/late)",
  silent = TRUE
)

pdf(
  "Heatmap_Genes_Zscore_4groups_6blocks_EdgeRbased.pdf",
  width = 7.5,
  height = pdf_h,
  useDingbats = FALSE
)
grid::grid.newpage()
grid::grid.draw(ph$gtable)
dev.off()


# ============================================================
# Pseudotime analysis
# Continuation:
# W8R (HisTrac) Dam L2-3 vs L4-5 pseudotime and module-score overlay
#
# Assumptions:
# - obj is the output of the previous script
# - obj already contains:
#     Cell_Type
#     joint.ref.lsi
#     joint.umap
#     RNA.Dam
# - the following files already exist in the working directory:
#     Pseudobulk_L2-3_vs_L4-5_edgeR.csv
#     EdgeR_geneSets_6groups_EarlyConsistentLate.csv
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(SingleCellExperiment)
library(slingshot)

stopifnot("Cell_Type" %in% colnames(obj@meta.data))
stopifnot("joint.ref.lsi" %in% Reductions(obj))
stopifnot("joint.umap" %in% Reductions(obj))
stopifnot("RNA.Dam" %in% names(obj@assays))

# ------------------------------------------------------------
# 1) Read DEG table and build L23 vs L45 delta score
# ------------------------------------------------------------
DefaultAssay(obj) <- "RNA.Dam"

expr_mat <- tryCatch(
  GetAssayData(obj, assay = "RNA.Dam", layer = "data"),
  error = function(e) GetAssayData(obj, assay = "RNA.Dam", slot = "data")
)
features_in_assay <- rownames(expr_mat)

deg_df <- read.csv("Pseudobulk_L2-3_vs_L4-5_edgeR.csv", stringsAsFactors = FALSE)

L23_edgeR <- intersect(
  deg_df$gene[deg_df$FDR < 0.05 & deg_df$logFC > 0],
  features_in_assay
)
L45_edgeR <- intersect(
  deg_df$gene[deg_df$FDR < 0.05 & deg_df$logFC < 0],
  features_in_assay
)

stopifnot(length(L23_edgeR) > 0)
stopifnot(length(L45_edgeR) > 0)

obj <- AddModuleScore(obj, features = list(L23_edgeR), name = "L23_EdgeR_")
obj <- AddModuleScore(obj, features = list(L45_edgeR), name = "L45_EdgeR_")
obj$L23vsL45_EdgeR_delta <- obj$L23_EdgeR_1 - obj$L45_EdgeR_1

# ------------------------------------------------------------
# 2) Read 6 module gene sets and add module scores
# ------------------------------------------------------------
gene_sets_df <- read.csv(
  "EdgeR_geneSets_6groups_EarlyConsistentLate.csv",
  stringsAsFactors = FALSE
)

set_order <- c(
  "L23_early_edgeR",
  "L23_consistent_edgeR",
  "L23_late_edgeR",
  "L45_early_edgeR",
  "L45_consistent_edgeR",
  "L45_late_edgeR"
)

for (nm in set_order) {
  genes_use <- intersect(
    gene_sets_df$gene[gene_sets_df$gene_set == nm],
    features_in_assay
  )
  stopifnot(length(genes_use) > 0)
  obj <- AddModuleScore(obj, features = list(genes_use), name = paste0(nm, "_"))
}

# ------------------------------------------------------------
# 3) Build pseudotime with Slingshot
#    Root/tips are defined by L23vsL45_EdgeR_delta
# ------------------------------------------------------------
obj_pt <- subset(obj, subset = Cell_Type %in% c("L2-3", "L4-5"))
obj_pt$Cell_Type <- factor(obj_pt$Cell_Type, levels = c("L2-3", "L4-5"))

stopifnot("L23vsL45_EdgeR_delta" %in% colnames(obj_pt@meta.data))
delta <- obj_pt$L23vsL45_EdgeR_delta

q_root_abs <- quantile(abs(delta), 0.20, na.rm = TRUE)
q_tip_L23  <- quantile(delta[obj_pt$Cell_Type == "L2-3"], 0.95, na.rm = TRUE)
q_tip_L45  <- quantile(delta[obj_pt$Cell_Type == "L4-5"], 0.05, na.rm = TRUE)

obj_pt$state_delta <- "mid"
obj_pt$state_delta[abs(delta) <= q_root_abs] <- "ambiguous"
obj_pt$state_delta[obj_pt$Cell_Type == "L2-3" & delta >= q_tip_L23] <- "L23_tip"
obj_pt$state_delta[obj_pt$Cell_Type == "L4-5" & delta <= q_tip_L45] <- "L45_tip"

obj_pt$state_delta[obj_pt$state_delta == "mid" & obj_pt$Cell_Type == "L2-3"] <- "L23_mid"
obj_pt$state_delta[obj_pt$state_delta == "mid" & obj_pt$Cell_Type == "L4-5"] <- "L45_mid"

obj_pt$state_delta <- factor(
  obj_pt$state_delta,
  levels = c("ambiguous", "L23_mid", "L23_tip", "L45_mid", "L45_tip")
)

X_lsi <- as.matrix(Embeddings(obj_pt, "joint.ref.lsi"))

sce <- as.SingleCellExperiment(obj_pt)
reducedDims(sce)$LSI <- X_lsi
sce$state_delta <- obj_pt$state_delta

sce <- slingshot(
  sce,
  clusterLabels = "state_delta",
  reducedDim = "LSI",
  start.clus = "ambiguous",
  end.clus = c("L23_tip", "L45_tip")
)

pt <- slingPseudotime(sce)
lin_end <- sapply(slingLineages(sce), function(z) tail(z, 1))
col23 <- which(lin_end == "L23_tip")[1]
col45 <- which(lin_end == "L45_tip")[1]

obj_pt$pt_L23 <- pt[, col23]
obj_pt$pt_L45 <- pt[, col45]
obj_pt$pt_branch <- ifelse(obj_pt$Cell_Type == "L2-3", obj_pt$pt_L23, obj_pt$pt_L45)

# ------------------------------------------------------------
# 4) Rank-based pseudotime for FeaturePlot
#    This is the corrected plotting version
# ------------------------------------------------------------
# Note:
# pt_rank01_warp01 is used for pseudotime visualization on the UMAP.
# Quantitative module-score overlays use pt_fromAmb01, which is centered on
# the ambiguous root population and restricted to the forward branch.

obj_pt$pt_rank <- rank(obj_pt$pt_branch, na.last = "keep", ties.method = "average")
obj_pt$pt_rank01 <- obj_pt$pt_rank / max(obj_pt$pt_rank, na.rm = TRUE)

# ---------------------
# pt_rank01 version
# ---------------------
vmin <- 0
vmid <- 0.25
vmax <- 0.7

mid_width <- 0.25
mid_span  <- c(0.01, 0.80)

stopifnot(vmin < vmid, vmid < vmax)

vlow  <- max(vmin, vmid - mid_width)
vhigh <- min(vmax, vmid + mid_width)

if (vlow >= vmid)  vlow  <- (vmin + vmid) / 2
if (vhigh <= vmid) vhigh <- (vmid + vmax) / 2
if (vlow <= vmin)  vlow  <- vmin + 1e-6
if (vhigh >= vmax) vhigh <- vmax - 1e-6

x_breaks <- c(vmin, vlow, vmid, vhigh, vmax)
y_breaks <- c(0, mid_span[1], 0.5, mid_span[2], 1)

x <- obj_pt$pt_rank01
x_clip <- pmin(pmax(x, vmin), vmax)

obj_pt$pt_rank01_warp01 <- NA_real_
ok <- is.finite(x_clip)
obj_pt$pt_rank01_warp01[ok] <- approx(
  x = x_breaks,
  y = y_breaks,
  xout = x_clip[ok],
  ties = "ordered"
)$y

pal <- viridisLite::plasma(256)

p_pt <- FeaturePlot(
  obj_pt,
  reduction = "joint.umap",
  features = "pt_rank01_warp01",
  min.cutoff = 0,
  max.cutoff = 1,
  pt.size = 2,
  combine = FALSE
)[[1]] +
  ggtitle(sprintf(
    "Pseudotime pt_rank01 (plasma; window %.2f-%.2f; mid %.2f)",
    vmin, vmax, vmid
  )) +
  scale_colour_gradientn(
    colours = pal,
    limits = c(0, 1),
    breaks = y_breaks,
    labels = round(x_breaks, 2),
    oob = scales::squish
  )

print(p_pt)

ggsave(
  "Pseudotime_pt_rank01_warp_plasma.pdf",
  p_pt,
  width = 6,
  height = 5,
  useDingbats = FALSE
)

# ------------------------------------------------------------
# 5) Forward pseudotime from ambiguous state
#    Used for overlay plots
# ------------------------------------------------------------
m0 <- median(obj_pt$pt_branch[obj_pt$state_delta == "ambiguous"], na.rm = TRUE)
obj_pt$pt_centered <- obj_pt$pt_branch - m0
obj_pt$pt_fromAmb <- ifelse(obj_pt$pt_centered < 0, NA_real_, obj_pt$pt_centered)
obj_pt$pt_fromAmb01 <- obj_pt$pt_fromAmb / max(obj_pt$pt_fromAmb, na.rm = TRUE)

# ------------------------------------------------------------
# 6) Prepare module-score table
# ------------------------------------------------------------
modules <- c(
  "L23_early_edgeR_1",
  "L23_consistent_edgeR_1",
  "L23_late_edgeR_1",
  "L45_early_edgeR_1",
  "L45_consistent_edgeR_1",
  "L45_late_edgeR_1"
)

minmax01_robust <- function(x, q = c(0.01, 0.99)) {
  lo <- as.numeric(quantile(x, q[1], na.rm = TRUE))
  hi <- as.numeric(quantile(x, q[2], na.rm = TRUE))
  if (!is.finite(lo) || !is.finite(hi) || hi - lo == 0) {
    return(rep(NA_real_, length(x)))
  }
  y <- (x - lo) / (hi - lo)
  pmin(pmax(y, 0), 1)
}

df_mod <- FetchData(
  obj_pt,
  vars = c("pt_fromAmb01", "Cell_Type", modules)
)

df_mod_long <- df_mod %>%
  pivot_longer(cols = all_of(modules), names_to = "module", values_to = "score") %>%
  mutate(
    family = if_else(str_detect(module, "^L23_"), "L23", "L45"),
    phase = case_when(
      str_detect(module, "early") ~ "early",
      str_detect(module, "consistent") ~ "consistent",
      str_detect(module, "late") ~ "late",
      TRUE ~ module
    ),
    phase = factor(phase, levels = c("early", "consistent", "late")),
    Cell_Type = factor(Cell_Type, levels = c("L2-3", "L4-5"))
  ) %>%
  group_by(module) %>%
  mutate(score01 = minmax01_robust(score)) %>%
  ungroup()

# ------------------------------------------------------------
# 7) Bin pseudotime and compute median trends
# ------------------------------------------------------------
bin_median <- function(dat, pt_col = "pt_fromAmb01", n_bins = 6) {
  dat <- dat %>% filter(is.finite(.data[[pt_col]]))
  
  x <- dat[[pt_col]]
  brks <- quantile(x, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  brks <- unique(brks)
  
  if (length(brks) < 3) {
    stop("Too few unique pseudotime breaks for binning.")
  }
  
  dat %>%
    mutate(pt_bin = cut(.data[[pt_col]], breaks = brks, include.lowest = TRUE)) %>%
    group_by(family, phase, module, Cell_Type, pt_bin) %>%
    summarise(
      pt_med = median(.data[[pt_col]], na.rm = TRUE),
      score_med = median(score01, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    arrange(family, module, Cell_Type, pt_med)
}

df_bin <- bin_median(df_mod_long, pt_col = "pt_fromAmb01", n_bins = 6)

# ------------------------------------------------------------
# 8) Overlay plot function
# ------------------------------------------------------------
smooth_predict_df_range <- function(dat, xlim = c(0, 0.3), span = 1.3, n_grid = 200) {
  dat <- dat %>%
    filter(is.finite(pt_med), is.finite(score_med)) %>%
    arrange(pt_med)
  
  if (nrow(dat) < 3 || length(unique(dat$pt_med)) < 3) {
    return(tibble(pt_grid = numeric(0), yhat = numeric(0)))
  }
  
  xmin <- max(xlim[1], min(dat$pt_med, na.rm = TRUE))
  xmax <- min(xlim[2], max(dat$pt_med, na.rm = TRUE))
  
  if (!is.finite(xmin) || !is.finite(xmax) || xmin >= xmax) {
    return(tibble(pt_grid = numeric(0), yhat = numeric(0)))
  }
  
  grid <- seq(xmin, xmax, length.out = n_grid)
  
  yhat <- tryCatch({
    fit <- loess(
      score_med ~ pt_med,
      data = dat,
      span = span,
      na.action = na.exclude,
      control = loess.control(surface = "direct")
    )
    as.numeric(predict(fit, newdata = data.frame(pt_med = grid)))
  }, error = function(e) {
    as.numeric(approx(dat$pt_med, dat$score_med, xout = grid, ties = mean)$y)
  })
  
  tibble(pt_grid = grid, yhat = yhat)
}

make_overlay_plot <- function(df_bin, family_use = c("L23", "L45"), xlim = c(0, 0.3), span = 1.3) {
  family_use <- match.arg(family_use)
  
  d <- df_bin %>%
    filter(family == family_use, is.finite(pt_med), is.finite(score_med))
  
  d_plot <- d %>% filter(pt_med >= xlim[1], pt_med <= xlim[2])
  
  curves <- d %>%
    group_by(module, phase, Cell_Type) %>%
    group_modify(~ smooth_predict_df_range(.x, xlim = xlim, span = span)) %>%
    ungroup() %>%
    filter(is.finite(pt_grid), is.finite(yhat))
  
  sf <- curves %>%
    group_by(module) %>%
    summarise(
      ymin = min(yhat, na.rm = TRUE),
      ymax = max(yhat, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(ymin), is.finite(ymax), ymax > ymin) %>%
    mutate(yrange = ymax - ymin)
  
  curves2 <- curves %>%
    left_join(sf, by = "module") %>%
    filter(is.finite(yrange)) %>%
    mutate(y01 = (yhat - ymin) / yrange)
  
  points2 <- d_plot %>%
    left_join(sf, by = "module") %>%
    filter(is.finite(yrange)) %>%
    mutate(
      y01 = (score_med - ymin) / yrange,
      y01 = pmin(pmax(y01, 0), 1)
    )
  
  ggplot() +
    geom_line(
      data = curves2,
      aes(
        x = pt_grid,
        y = y01,
        color = phase,
        linetype = Cell_Type,
        group = interaction(module, Cell_Type)
      ),
      linewidth = 1
    ) +
    geom_point(
      data = points2,
      aes(
        x = pt_med,
        y = y01,
        color = phase,
        shape = Cell_Type,
        group = interaction(module, Cell_Type)
      ),
      size = 1,
      alpha = 0.5
    ) +
    theme_classic(base_size = 12) +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      x = "Pseudotime",
      y = "Scaled trend (0-1)",
      color = "Module phase",
      linetype = "Cell type",
      shape = "Cell type"
    ) +
    ggtitle(paste0(family_use, " overlay"))
}

# ------------------------------------------------------------
# 9) Final overlay plots and PDFs
# ------------------------------------------------------------
p_L23 <- make_overlay_plot(df_bin, family_use = "L23", xlim = c(0, 0.3), span = 1.3)
p_L45 <- make_overlay_plot(df_bin, family_use = "L45", xlim = c(0, 0.3), span = 1.3)

print(p_L23)
print(p_L45)

ggsave(
  "L23gene_L23_L45_overlay_xlim.pdf",
  p_L23,
  width = 8,
  height = 6,
  useDingbats = FALSE
)

ggsave(
  "L45gene_L23_L45_overlay_xlim.pdf",
  p_L45,
  width = 8,
  height = 6,
  useDingbats = FALSE
)
