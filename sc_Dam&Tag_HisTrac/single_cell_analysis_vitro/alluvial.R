# ============================================================
# Alluvial analysis for one stage-specific object
#
# Input example:
#   combined.obj.ls_Day2.no2
#   or directly:
#   combined.obj.ls_Day2.no2$Dam_bb
#
# Analytical premise:
# - Dam and K27ac modalities were measured for the same cells.
# - Identity assignment was performed independently in each modality, producing:
#     Dam_annotation_new
#     K27ac_annotation_new
# - In the downstream object used here, these two identity calls are already
#   stored together in Dam_bb metadata.
# - Therefore, Dam_bb alone is sufficient to visualize identity transitions
#   between Dam-defined and K27ac-defined states.
# - This script does not separately use K27ac_aa.
#
# Assumptions:
# - The object is already filtered to the two target identities.
# - Metadata columns below already exist in Dam_bb:
#     sample
#     Dam_annotation_new
#     K27ac_annotation_new
# - Dam_annotation_new and K27ac_annotation_new contain only 0 / 1 (id-C vs id-P).
#
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(patchwork)

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

make_dam_alluvial_plot <- function(
    x,
    stage_name = "stage",
    sample_value = NULL,
    sample_col = "sample",
    dam_col = "Dam_annotation_new",
    k27ac_col = "K27ac_annotation_new",
    dam_label = "Dam",
    k27ac_label = "K27ac") {

  dam_obj <- get_dam_bb(x)

  required_cols <- c(dam_col, k27ac_col)
  if (!is.null(sample_value)) {
    required_cols <- c(required_cols, sample_col)
  }
  check_required_columns(dam_obj, required_cols)

  meta <- dam_obj@meta.data

  df <- data.frame(
    cell = rownames(meta),
    Dam_cluster = factor(as.character(meta[[dam_col]]), levels = c("0", "1")),
    K27ac_cluster = factor(as.character(meta[[k27ac_col]]), levels = c("0", "1")),
    stringsAsFactors = FALSE
  )

  if (!is.null(sample_value)) {
    df$sample <- as.character(meta[[sample_col]])
    df <- df[df$sample == sample_value, , drop = FALSE]
  }

  if (nrow(df) == 0) {
    stop("No cells available for plotting.", call. = FALSE)
  }

  df_count <- df %>%
    count(Dam_cluster, K27ac_cluster, name = "n", .drop = FALSE) %>%
    filter(n > 0) %>%
    mutate(
      fraction = 100 * n / sum(n),
      transition_fraction_label = paste0(
        Dam_cluster, " -> ", K27ac_cluster,
        " (", round(fraction, 1), "%)"
      )
    )

  legend_df <- data.frame(
    transition_fraction_label = unique(df_count$transition_fraction_label),
    x = Inf,
    y = Inf
  )
  color_vec <- setNames(
    rep("black", nrow(legend_df)),
    legend_df$transition_fraction_label
  )

  subtitle_text <- if (is.null(sample_value)) {
    NULL
  } else {
    paste0("sample = ", sample_value)
  }

  ggplot(
    data = df_count,
    aes(axis1 = Dam_cluster, axis2 = K27ac_cluster, y = n)
  ) +
    geom_alluvium(aes(fill = Dam_cluster), width = 1 / 12) +
    geom_stratum(width = 1 / 12, fill = "grey85", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    geom_point(
      data = legend_df,
      mapping = aes(x = x, y = y, color = transition_fraction_label),
      alpha = 0,
      show.legend = TRUE,
      size = 5,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = color_vec, name = "Transition (%)") +
    scale_x_discrete(
      limits = c(dam_label, k27ac_label),
      expand = c(0.05, 0.05)
    ) +
    labs(
      title = paste0(stage_name, ": Alluvial plot"),
      subtitle = subtitle_text,
      x = "Annotation",
      y = "Number of cells"
    ) +
    guides(
      fill = guide_legend(title = paste0(dam_label, " cluster")),
      color = guide_legend(
        title = "Transition (%)",
        override.aes = list(shape = 15, size = 5)
      )
    ) +
    theme_minimal()
}

run_dam_alluvial_analysis <- function(
    x,
    stage_name = "stage",
    outdir = NULL,
    sample_col = "sample",
    dam_col = "Dam_annotation_new",
    k27ac_col = "K27ac_annotation_new",
    plot_by_sample = TRUE,
    dam_label = "Dam",
    k27ac_label = "K27ac") {

  dam_obj <- get_dam_bb(x)
  check_required_columns(dam_obj, c(dam_col, k27ac_col))

  overall_plot <- make_dam_alluvial_plot(
    x = dam_obj,
    stage_name = stage_name,
    sample_value = NULL,
    sample_col = sample_col,
    dam_col = dam_col,
    k27ac_col = k27ac_col,
    dam_label = dam_label,
    k27ac_label = k27ac_label
  )

  by_sample_plots <- list()
  by_sample_panel <- NULL

  if (plot_by_sample) {
    if (!sample_col %in% colnames(dam_obj@meta.data)) {
      warning("sample column not found. Per-sample alluvial plots were skipped.")
    } else {
      sample_ids <- unique(as.character(dam_obj@meta.data[[sample_col]]))
      sample_ids <- sample_ids[!is.na(sample_ids)]

      by_sample_plots <- setNames(
        lapply(sample_ids, function(sid) {
          make_dam_alluvial_plot(
            x = dam_obj,
            stage_name = stage_name,
            sample_value = sid,
            sample_col = sample_col,
            dam_col = dam_col,
            k27ac_col = k27ac_col,
            dam_label = dam_label,
            k27ac_label = k27ac_label
          )
        }),
        sample_ids
      )

      if (length(by_sample_plots) > 0) {
        ncol_panel <- min(2, length(by_sample_plots))
        by_sample_panel <- patchwork::wrap_plots(by_sample_plots, ncol = ncol_panel)
      }
    }
  }

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    ggsave(
      filename = file.path(outdir, paste0(stage_name, "_alluvial_overall.pdf")),
      plot = overall_plot,
      width = 7,
      height = 5
    )

    if (!is.null(by_sample_panel)) {
      ncol_panel <- min(2, length(by_sample_plots))
      nrow_panel <- ceiling(length(by_sample_plots) / ncol_panel)

      ggsave(
        filename = file.path(outdir, paste0(stage_name, "_alluvial_by_sample.pdf")),
        plot = by_sample_panel,
        width = 7 * ncol_panel,
        height = 5 * nrow_panel
      )
    }
  }

  list(
    dam_obj = dam_obj,
    overall_plot = overall_plot,
    by_sample_plots = by_sample_plots,
    by_sample_panel = by_sample_panel
  )
}

# ============================================================
# Example
# ============================================================
# source("dam_alluvial_minimal.R")
#
# Example 1: input is the stage-specific list containing Dam_bb
# res_day2 <- run_dam_alluvial_analysis(
#   x = combined.obj.ls_Day2.no2,
#   stage_name = "Day2",
#   outdir = "Day2_alluvial_analysis"
# )
#
# res_day2$overall_plot
# res_day2$by_sample_panel
#
# To disable per-sample plots:
# res_day2 <- run_dam_alluvial_analysis(
#   x = combined.obj.ls_Day2.no2,
#   stage_name = "Day2",
#   plot_by_sample = FALSE
# )
