# ============================================================
# K27ac-medoid transition analysis for one stage-specific object
#
# Input example:
#   combined.obj.ls_Day2.no2
#   or directly:
#   combined.obj.ls_Day2.no2$K27ac_aa
#
# Analytical premise:
# - LSI was first constructed on an object merged across all stages
#   (Day2, Day7, and Day7R), rather than separately within each stage.
# - Therefore, the embedding used here is not learned from Day2 alone,
#   but from the stage-merged object.
# - The object analyzed in this script is a stage-specific subset extracted
#   from that merged object (for example, Day2 only).
#
# Modality premise:
# - Dam and K27ac modalities were measured for the same cells.
# - Identity assignment was performed independently in each modality, producing:
#     Dam_annotation_new
#     K27ac_annotation_new
# - In the downstream object used here, these two identity calls are already
#   stored together in K27ac_aa metadata.
# - Therefore, transition status (homo vs hetero) can be evaluated from K27ac_aa
#   alone, without separately using Dam_bb.
# - The distance calculation in this script uses only K27ac medoids, defined in
#   the LSI space stored in K27ac_aa.
#
# Assumptions:
# - The input object is already filtered to the two target identities.
# - K27ac_aa already contains the LSI reduction inherited from the merged
#   Day2/Day7/Day7R object.
# - Metadata columns below already exist in K27ac_aa:
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
library(cluster)

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

get_k27ac_aa <- function(x) {
  if (inherits(x, "Seurat")) {
    return(x)
  }

  if (is.list(x) && "K27ac_aa" %in% names(x)) {
    if (!inherits(x$K27ac_aa, "Seurat")) {
      stop("x$K27ac_aa exists but is not a Seurat object.", call. = FALSE)
    }
    return(x$K27ac_aa)
  }

  stop(
    "Input must be either a Seurat object (K27ac_aa) or a list containing K27ac_aa.",
    call. = FALSE
  )
}

get_cluster_medoid <- function(embedding, cluster_labels, cluster_id) {
  cluster_cells <- rownames(embedding)[cluster_labels == cluster_id]

  if (length(cluster_cells) == 0) {
    stop("No cells found for K27ac cluster ", cluster_id, call. = FALSE)
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

add_k27ac_opposite_medoid_distance <- function(
    x,
    reduction = "lsi",
    dims = 2:15,
    dam_col = "Dam_annotation_new",
    k27ac_col = "K27ac_annotation_new") {

  k27ac_obj <- get_k27ac_aa(x)

  check_required_columns(k27ac_obj, c(dam_col, k27ac_col))

  if (!reduction %in% names(k27ac_obj@reductions)) {
    stop("Reduction '", reduction, "' was not found in K27ac_aa.", call. = FALSE)
  }

  embedding_full <- Embeddings(k27ac_obj, reduction = reduction)

  if (max(dims) > ncol(embedding_full)) {
    stop(
      "Requested dims exceed available dimensions in reduction '", reduction, "'.",
      call. = FALSE
    )
  }

  embedding <- embedding_full[, dims, drop = FALSE]
  dam_labels <- as.character(k27ac_obj@meta.data[[dam_col]])
  k27ac_labels <- as.character(k27ac_obj@meta.data[[k27ac_col]])

  required_clusters <- c("0", "1")
  if (!all(required_clusters %in% unique(dam_labels))) {
    stop("Dam labels must contain both 0 and 1.", call. = FALSE)
  }
  if (!all(required_clusters %in% unique(k27ac_labels))) {
    stop("K27ac labels must contain both 0 and 1.", call. = FALSE)
  }

  # K27ac medoids are defined only from K27ac_annotation_new.
  k27ac_medoid_0 <- get_cluster_medoid(embedding, k27ac_labels, "0")
  k27ac_medoid_1 <- get_cluster_medoid(embedding, k27ac_labels, "1")

  # Distance of every cell to the two K27ac medoids.
  dist_k27ac_0 <- euclidean_to_medoid(embedding, k27ac_medoid_0)
  dist_k27ac_1 <- euclidean_to_medoid(embedding, k27ac_medoid_1)

  # Opposite medoid distance (OMD) based only on K27ac identity.
  k27ac_obj$dist_contra_K27ac <- ifelse(k27ac_labels == "1", dist_k27ac_0, dist_k27ac_1)

  # Transition status is still defined by comparing Dam and K27ac identity calls.
  k27ac_obj$transition_type <- factor(
    ifelse(dam_labels == k27ac_labels, "homo", "hetero"),
    levels = c("homo", "hetero")
  )
  k27ac_obj$transition <- paste0(dam_labels, "_to_", k27ac_labels)

  attr(k27ac_obj, "k27ac_medoid_0") <- k27ac_medoid_0
  attr(k27ac_obj, "k27ac_medoid_1") <- k27ac_medoid_1

  k27ac_obj
}

plot_k27ac_omd_homo_hetero <- function(k27ac_obj, stage_name = "stage") {
  check_required_columns(k27ac_obj, c("dist_contra_K27ac", "transition_type"))

  meta <- k27ac_obj@meta.data

  meta %>%
    mutate(
      transition_type = factor(
        as.character(transition_type),
        levels = c("homo", "hetero")
      )
    ) %>%
    filter(is.finite(dist_contra_K27ac), dist_contra_K27ac > 0) %>%
    ggplot(aes(x = transition_type, y = log2(dist_contra_K27ac), fill = transition_type)) +
    geom_boxplot(width = 0.6) +
    scale_x_discrete(labels = c(homo = "Homo", hetero = "Hetero")) +
    scale_fill_manual(values = c(homo = "darkblue", hetero = "darkorange")) +
    labs(
      title = paste0(stage_name, ": K27ac contra-medoid distance"),
      x = NULL,
      y = "log2(distance)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}

run_k27ac_medoid_transition_analysis <- function(
    x,
    stage_name = "stage",
    outdir = NULL,
    reduction = "lsi",
    dims = 2:15,
    dam_col = "Dam_annotation_new",
    k27ac_col = "K27ac_annotation_new") {

  k27ac_obj <- add_k27ac_opposite_medoid_distance(
    x = x,
    reduction = reduction,
    dims = dims,
    dam_col = dam_col,
    k27ac_col = k27ac_col
  )

  omd_plot <- plot_k27ac_omd_homo_hetero(k27ac_obj, stage_name = stage_name)

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    ggsave(
      filename = file.path(outdir, paste0(stage_name, "_K27ac_medoid_distance_homo_vs_hetero.pdf")),
      plot = omd_plot,
      width = 5,
      height = 5
    )
  }

  list(
    k27ac_obj = k27ac_obj,
    omd_plot = omd_plot,
    k27ac_medoid_0 = attr(k27ac_obj, "k27ac_medoid_0"),
    k27ac_medoid_1 = attr(k27ac_obj, "k27ac_medoid_1")
  )
}

# ============================================================
# Example
# ============================================================
#
# Example: input is the stage-specific list containing K27ac_aa
# res_day2 <- run_k27ac_medoid_transition_analysis(
#   x = combined.obj.ls_Day2.no2,
#   stage_name = "Day2",
#   outdir = "Day2_k27ac_medoid_analysis"
# )
# res_day2$omd_plot
# head(res_day2$k27ac_obj@meta.data[, c("Dam_annotation_new", "K27ac_annotation_new",
#                                       "transition_type", "dist_contra_K27ac")])
#
# The same function can be used for Day7 / Day7R.
