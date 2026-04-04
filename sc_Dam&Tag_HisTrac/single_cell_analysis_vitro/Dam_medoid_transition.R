# ============================================================
# Dam-medoid transition analysis for one stage-specific object
#
# Input example:
#   combined.obj.ls_Day2.no2
#   or directly:
#   combined.obj.ls_Day2.no2$Dam_bb
#
# Analytical premise:
# - A shared low-dimensional space (LSI) was first constructed from a merged
#   object containing all stages (Day2, Day7, and Day7R).
# - Therefore, the LSI used here is not stage-specific; it was learned jointly
#   across Day2, Day7, and Day7R in the merged object.
# - The object analyzed in this script corresponds to a stage-specific subset
#   extracted from that merged object (for example, Day2 only).
# - Accordingly, this script analyzes one stage at a time, but always in the
#   context of the shared LSI originally defined across all stages.
#
# Modality premise:
# - Dam and K27ac modalities were measured for the same cells.
# - Identity assignment was performed independently in each modality, producing:
#     Dam_annotation_new
#     K27ac_annotation_new
# - In the downstream object used here, these two identity calls are already
#   stored together in Dam_bb metadata.
# - Therefore, transition status (homo vs hetero) can be evaluated from Dam_bb
#   alone, without separately using K27ac_aa.
# - The distance calculation in this script uses only Dam medoids, defined in
#   the shared LSI space stored in Dam_bb.
#
# Assumptions:
# - The input object is already filtered to the two target identities.
# - Dam_bb already contains the shared LSI reduction inherited from the
#   merged Day2/Day7/Day7R object.
# - Metadata columns below already exist in Dam_bb:
#     sample
#     Dam_annotation_new
#     K27ac_annotation_new
# - Dam_annotation_new and K27ac_annotation_new contain only 0 / 1
#   (for example, id-C vs id-P).
#
# ============================================================


library(Seurat)
library(dplyr)
library(ggplot2)

check_required_columns <- function(obj, cols) {
  missing_cols <- setdiff(cols, colnames(obj@meta.data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

get_dam_bb <- function(x) {
  if (inherits(x, "Seurat")) {
    return(x)
  }
  
  if (is.list(x) && "Dam_bb" %in% names(x)) {
    if (!inherits(x$Dam_bb, "Seurat")) {
      stop("x$Dam_bb exists but is not a Seurat object.", call. = FALSE)
    }
    return(x$Dam_bb)
  }
  
  stop(
    "Input must be either a Seurat object (Dam_bb) or a list containing Dam_bb.",
    call. = FALSE
  )
}

get_cluster_medoid <- function(embedding, cluster_labels, cluster_id) {
  cluster_cells <- rownames(embedding)[cluster_labels == cluster_id]
  
  if (length(cluster_cells) == 0) {
    stop("No cells found for Dam cluster ", cluster_id, call. = FALSE)
  }
  
  cluster_embedding <- embedding[cluster_cells, , drop = FALSE]
  
  if (nrow(cluster_embedding) == 1) {
    return(as.numeric(cluster_embedding[1, ]))
  }
  
  pam_fit <- cluster::pam(cluster_embedding, k = 1)
  as.numeric(cluster_embedding[pam_fit$id.med, , drop = TRUE])
}

euclidean_to_medoid <- function(embedding, medoid_vec) {
  sqrt(rowSums((sweep(embedding, 2, medoid_vec, FUN = "-"))^2))
}

add_dam_opposite_medoid_distance <- function(
    x,
    reduction = "lsi",
    dims = 2:15,
    dam_col = "Dam_annotation_new",
    k27ac_col = "K27ac_annotation_new") {
  
  dam_obj <- get_dam_bb(x)
  
  check_required_columns(dam_obj, c(dam_col, k27ac_col))
  
  if (!reduction %in% names(dam_obj@reductions)) {
    stop("Reduction '", reduction, "' was not found in Dam_bb.", call. = FALSE)
  }
  
  embedding_full <- Embeddings(dam_obj, reduction = reduction)
  
  if (max(dims) > ncol(embedding_full)) {
    stop(
      "Requested dims exceed available dimensions in reduction '", reduction, "'.",
      call. = FALSE
    )
  }
  
  embedding <- embedding_full[, dims, drop = FALSE]
  dam_labels <- as.character(dam_obj@meta.data[[dam_col]])
  k27ac_labels <- as.character(dam_obj@meta.data[[k27ac_col]])
  
  required_clusters <- c("0", "1")
  if (!all(required_clusters %in% unique(dam_labels))) {
    stop("Dam labels must contain both 0 and 1.", call. = FALSE)
  }
  if (!all(required_clusters %in% unique(k27ac_labels))) {
    stop("K27ac labels must contain both 0 and 1.", call. = FALSE)
  }
  
  # Dam medoids are defined only from Dam_annotation_new.
  dam_medoid_0 <- get_cluster_medoid(embedding, dam_labels, "0")
  dam_medoid_1 <- get_cluster_medoid(embedding, dam_labels, "1")
  
  # Distance of every cell to the two Dam medoids.
  dist_dam_0 <- euclidean_to_medoid(embedding, dam_medoid_0)
  dist_dam_1 <- euclidean_to_medoid(embedding, dam_medoid_1)
  
  # Opposite medoid distance (OMD) based only on Dam identity.
  dam_obj$dist_contra_Dam <- ifelse(dam_labels == "1", dist_dam_0, dist_dam_1)
  
  # Transition status is still defined by comparing Dam and K27ac identity calls.
  dam_obj$transition_type <- factor(
    ifelse(dam_labels == k27ac_labels, "homo", "hetero"),
    levels = c("homo", "hetero")
  )
  dam_obj$transition <- paste0(dam_labels, "_to_", k27ac_labels)
  
  attr(dam_obj, "dam_medoid_0") <- dam_medoid_0
  attr(dam_obj, "dam_medoid_1") <- dam_medoid_1
  
  dam_obj
}

plot_dam_omd_homo_hetero <- function(dam_obj, stage_name = "stage") {
  check_required_columns(dam_obj, c("dist_contra_Dam", "transition_type"))
  
  meta <- dam_obj@meta.data
  
  meta %>%
    mutate(
      transition_type = factor(
        as.character(transition_type),
        levels = c("homo", "hetero")
      )
    ) %>%
    filter(is.finite(dist_contra_Dam), dist_contra_Dam > 0) %>%
    ggplot(aes(x = transition_type, y = log2(dist_contra_Dam), fill = transition_type)) +
    geom_boxplot(width = 0.6) +
    scale_x_discrete(labels = c(homo = "Homo", hetero = "Hetero")) +
    scale_fill_manual(values = c(homo = "darkblue", hetero = "darkorange")) +
    labs(
      title = paste0(stage_name, ": Dam contra-medoid distance"),
      x = NULL,
      y = "log2(distance)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}

run_dam_medoid_transition_analysis <- function(
    x,
    stage_name = "stage",
    outdir = NULL,
    reduction = "lsi",
    dims = 2:15,
    dam_col = "Dam_annotation_new",
    k27ac_col = "K27ac_annotation_new") {
  
  dam_obj <- add_dam_opposite_medoid_distance(
    x = x,
    reduction = reduction,
    dims = dims,
    dam_col = dam_col,
    k27ac_col = k27ac_col
  )
  
  omd_plot <- plot_dam_omd_homo_hetero(dam_obj, stage_name = stage_name)
  
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    
    ggsave(
      filename = file.path(outdir, paste0(stage_name, "_Dam_medoid_distance_homo_vs_hetero.pdf")),
      plot = omd_plot,
      width = 5,
      height = 5
    )
  }
  
  list(
    dam_obj = dam_obj,
    omd_plot = omd_plot,
    dam_medoid_0 = attr(dam_obj, "dam_medoid_0"),
    dam_medoid_1 = attr(dam_obj, "dam_medoid_1")
  )
}

# ============================================================
# Example
# ============================================================
# source("dam_medoid_transition_minimal.R")
#
# Example 1: input is the stage-specific list containing Dam_bb
# res_day2 <- run_dam_medoid_transition_analysis(
#   x = combined.obj.ls_Day2.no2,
#   stage_name = "Day2",
#   outdir = "Day2_dam_medoid_analysis"
# )
# res_day2$omd_plot
# head(res_day2$dam_obj@meta.data[, c("Dam_annotation_new", "K27ac_annotation_new",
#                                     "transition_type", "dist_contra_Dam")])
#
# The same function can be used for Day7 / Day7R.

