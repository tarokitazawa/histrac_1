# In vivo medoid heatmap examples
# Input object: DamLeo1_vivoHisTrac_rep_merged.rds
#
#   - K27ac: based on lsi
#   - Dam:   based on joint.ref.lsi
# Heatmap values:
#   - Pearson correlation between stage medoids
#
# Clustering:
#   - Euclidean distance between stage medoids
#
# Note on joint.ref.lsi:
# - K27ac medoids are computed in the original K27ac LSI space.
# - Dam medoids are computed in joint.ref.lsi, a shared reference-space LSI in which
#   W8/W8R reference cells retain their original Dam LSI coordinates and P8 cells are
#   projected onto that adult reference.
# - Therefore, joint.ref.lsi uses a re-indexed reference-space representation rather
#   than the original raw LSI numbering.
# - For this reason, Dam analyses based on joint.ref.lsi appropriately start from
#   dimension 1, even though the original Dam LSI workflow typically used dimensions
#   beginning from 2 after excluding the first raw LSI component.

library(Seurat)
library(pheatmap)

# -------------------------------------------------------------------
# Input
# -------------------------------------------------------------------
combined.obj.ls3 <- readRDS("DamLeo1_vivoHisTrac_rep_merged.rds")

# -------------------------------------------------------------------
# Helper: rotate 3-stage clustering to the preferred stage order
# -------------------------------------------------------------------
rotate_hc_best_stage3 <- function(hc, desired = c("P8", "W8R", "W8")) {
  labs <- hc$labels
  desired <- desired[desired %in% labs]
  
  stopifnot(length(labs) == 3)
  
  m1 <- hc$merge[1, ]
  pair_idx  <- abs(m1[m1 < 0])
  pair_labs <- labs[pair_idx]
  third_lab <- setdiff(labs, pair_labs)
  
  cand <- list(
    c(pair_labs, third_lab),
    c(rev(pair_labs), third_lab),
    c(third_lab, pair_labs),
    c(third_lab, rev(pair_labs))
  )
  
  rank_des <- setNames(seq_along(desired), desired)
  inv_cost <- function(o) {
    p <- rank_des[o]
    as.integer(p[1] > p[2]) +
      as.integer(p[1] > p[3]) +
      as.integer(p[2] > p[3])
  }
  
  best <- cand[[which.min(sapply(cand, inv_cost))]]
  hc$order <- match(best, labs)
  hc
}

# -------------------------------------------------------------------
# Helper: compute stage medoids and plot correlation heatmap
# -------------------------------------------------------------------
plot_stage_medoid_correlation <- function(
    obj,
    reduction,
    dims_use,
    cell_type = "L2-3",
    file_pdf,
    plot_title = NULL,
    desired_stage_order = c("P8", "W8R", "W8")) {
  
  stopifnot("Cell_Type" %in% colnames(obj@meta.data))
  stopifnot("stage" %in% colnames(obj@meta.data))
  stopifnot(reduction %in% Reductions(obj))
  
  obj_sub <- subset(obj, subset = Cell_Type == cell_type)
  
  emb_all <- Embeddings(obj_sub, reduction = reduction)
  if (max(dims_use) > ncol(emb_all)) {
    stop("Requested dimensions exceed the number of available dimensions in reduction: ", reduction)
  }
  
  emb <- emb_all[, dims_use, drop = FALSE]
  samp <- as.character(obj_sub$stage)
  
  if (length(unique(samp)) != 3) {
    warning("Expected 3 stages, found: ", paste(unique(samp), collapse = ", "))
  }
  
  # ---------------------------------------------------------------
  # Compute one medoid per stage in the chosen low-dimensional space
  # ---------------------------------------------------------------
  medoids <- do.call(rbind, lapply(unique(samp), function(s) {
    idx <- which(samp == s)
    sub_emb <- emb[idx, , drop = FALSE]
    dmat <- as.matrix(dist(sub_emb, method = "euclidean"))
    medoid_row <- which.min(rowSums(dmat))
    sub_emb[medoid_row, , drop = FALSE]
  }))
  rownames(medoids) <- unique(samp)
  
  # ---------------------------------------------------------------
  # Heatmap values = Pearson correlation between stage medoids
  # ---------------------------------------------------------------
  cor_med <- cor(t(medoids), method = "pearson", use = "pairwise.complete.obs")
  
  # ---------------------------------------------------------------
  # Clustering = Euclidean distance between stage medoids
  # ---------------------------------------------------------------
  dist_med <- dist(medoids, method = "euclidean")
  hc_dist <- hclust(dist_med, method = "ward.D2")
  hc_dist_fix <- rotate_hc_best_stage3(hc_dist, desired = desired_stage_order)
  
  if (is.null(plot_title)) {
    plot_title <- paste0(cell_type, " medoid correlation (", reduction, ")")
  }
  
  pheatmap::pheatmap(
    cor_med,
    cluster_rows = hc_dist_fix,
    cluster_cols = hc_dist_fix,
    main = plot_title,
    filename = file_pdf,
    width = 6,
    height = 6
  )
  
  invisible(list(
    object = obj_sub,
    medoids = medoids,
    cor_mat = cor_med,
    dist_mat = as.matrix(dist_med),
    hc = hc_dist_fix
  ))
}

# -------------------------------------------------------------------
# Example 1: K27ac L4-5 using lsi
# -------------------------------------------------------------------
res_k27_l45 <- plot_stage_medoid_correlation(
  obj = combined.obj.ls3$K27ac_aa,
  reduction = "lsi",
  dims_use = 2:14,
  cell_type = "L4-5",
  file_pdf = "Leo1_K27ac_L45.pdf",
  plot_title = "K27ac L4-5: stage medoid correlation"
)

# -------------------------------------------------------------------
# Example 2: Dam L4-5 using joint.ref.lsi
# -------------------------------------------------------------------
res_dam_l45 <- plot_stage_medoid_correlation(
  obj = combined.obj.ls3$Dam_bb,
  reduction = "joint.ref.lsi",
  dims_use = 1:13,
  cell_type = "L4-5",
  file_pdf = "Leo1_Dam_L45.pdf",
  plot_title = "Dam L4-5: stage medoid correlation"
)

