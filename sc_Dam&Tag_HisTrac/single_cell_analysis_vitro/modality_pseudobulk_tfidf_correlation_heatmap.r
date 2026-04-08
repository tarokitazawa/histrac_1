# ============================================================
# Correlation heatmap script for mean TF-IDF on top LSI peaks
#
# Input example:
#   combined.obj.ls
#
# Analytical premise:
# - combined.obj.ls is a merged Seurat object list containing Day2, Day7, and
#   Day7R samples together.
# - The list contains two modality-specific Seurat objects:
#     Dam_bb
#     K27ac_aa
# - Each modality object already contains a peaks assay and an LSI reduction.
# - This script computes, for each modality separately:
#     1) top LSI peaks
#     2) mean TF-IDF per sample/group
#     3) Pearson correlation between group-level mean TF-IDF profiles
#     4) a clustered correlation heatmap saved as PDF
#
# Notes:
# - This script intentionally keeps only the essential steps used for the final
#   PDF heatmaps.
# - Custom color-scale tuning is intentionally omitted.
# - By default, columns/rows are grouped by 'orig.ident', but this can be
#   changed to another metadata column such as 'sample'.
# ============================================================

library(Signac)
library(Seurat)
library(Matrix)
library(pheatmap)
library(dendextend)

get_modality_object <- function(x, modality = c("Dam", "K27ac")) {
  modality <- match.arg(modality)

  if (inherits(x, "Seurat")) {
    return(x)
  }

  if (!is.list(x)) {
    stop("Input must be either a Seurat object or a list containing modality objects.", call. = FALSE)
  }

  obj_name <- if (modality == "Dam") "Dam_bb" else "K27ac_aa"

  if (!obj_name %in% names(x)) {
    stop("List does not contain ", obj_name, ".", call. = FALSE)
  }

  obj <- x[[obj_name]]
  if (!inherits(obj, "Seurat")) {
    stop(obj_name, " exists but is not a Seurat object.", call. = FALSE)
  }

  obj
}

check_required_elements <- function(obj, assay = "peaks", reduction = "lsi", group_col = "orig.ident") {
  if (!assay %in% names(obj@assays)) {
    stop("Assay '", assay, "' was not found in the Seurat object.", call. = FALSE)
  }
  if (!reduction %in% names(obj@reductions)) {
    stop("Reduction '", reduction, "' was not found in the Seurat object.", call. = FALSE)
  }
  if (!group_col %in% colnames(obj@meta.data)) {
    stop("Metadata column '", group_col, "' was not found in the Seurat object.", call. = FALSE)
  }
}

get_top_lsi_peaks <- function(obj, reduction = "lsi", dims = 2:30, top_n = 1000) {
  loadings <- Loadings(obj[[reduction]])

  if (is.null(loadings) || nrow(loadings) == 0) {
    stop("No feature loadings found in reduction '", reduction, "'.", call. = FALSE)
  }
  if (max(dims) > ncol(loadings)) {
    stop("Requested dims exceed available dimensions in reduction '", reduction, "'.", call. = FALSE)
  }

  importance <- rowSums(loadings[, dims, drop = FALSE]^2)
  importance <- sort(importance, decreasing = TRUE)

  head(names(importance), top_n)
}

compute_mean_tfidf_by_group <- function(
    obj,
    group_col = "orig.ident",
    assay = "peaks",
    slot = "data") {

  tfidf <- GetAssayData(obj, assay = assay, slot = slot)
  meta <- obj@meta.data
  groups <- as.character(meta[[group_col]])
  cells <- rownames(meta)

  valid <- !is.na(groups)
  groups <- groups[valid]
  cells  <- cells[valid]

  group_levels <- unique(groups)

  mean_list <- lapply(group_levels, function(g) {
    g_cells <- cells[groups == g]
    Matrix::rowMeans(tfidf[, g_cells, drop = FALSE])
  })

  mean_mat <- do.call(cbind, mean_list)
  rownames(mean_mat) <- rownames(tfidf)
  colnames(mean_mat) <- group_levels

  mean_mat
}

make_tfidf_correlation_heatmap_pdf <- function(
    x,
    modality = c("Dam", "K27ac"),
    output_pdf,
    group_col = ,
    assay = ,
    reduction = "lsi",
    dims = 2:15,
    top_n = 1000,
    clustering_method = "ward.D",
    reverse_order = FALSE) {

  modality <- match.arg(modality)
  obj <- get_modality_object(x, modality = modality)
  check_required_elements(obj, assay = assay, reduction = reduction, group_col = group_col)

  top_peaks <- get_top_lsi_peaks(
    obj = obj,
    reduction = reduction,
    dims = dims,
    top_n = top_n
  )

  tfidf_mean <- compute_mean_tfidf_by_group(
    obj = obj,
    group_col = group_col,
    assay = assay,
    slot = "data"
  )

  common_peaks <- intersect(rownames(tfidf_mean), top_peaks)
  if (length(common_peaks) < 2) {
    stop("Fewer than two common peaks were found between TF-IDF matrix and top LSI peaks.", call. = FALSE)
  }

  tfidf_top <- tfidf_mean[common_peaks, , drop = FALSE]
  cor_mat <- cor(tfidf_top, method = "pearson")
  hc <- hclust(dist(t(tfidf_top)), method = clustering_method)

  if (reverse_order) {
    hc <- as.hclust(rev(as.dendrogram(hc)))
  }

  plot_title <- paste0(
    modality,
    ": correlation heatmap (mean TF-IDF on top LSI peaks)"
  )

  pdf_width  <- max(6, ncol(cor_mat) * 0.6 + 3)
  pdf_height <- max(6, nrow(cor_mat) * 0.6 + 3)

  pheatmap::pheatmap(
    cor_mat,
    cluster_rows = hc,
    cluster_cols = hc,
    main = plot_title,
    filename = output_pdf,
    width = pdf_width,
    height = pdf_height
  )

  invisible(list(
    object = obj,
    top_peaks = top_peaks,
    tfidf_mean = tfidf_mean,
    tfidf_top = tfidf_top,
    cor_mat = cor_mat,
    hc = hc
  ))
}
