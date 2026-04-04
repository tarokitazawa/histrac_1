# ============================================================
# RA analysis for Dam_1_0_UPscore_GeneActivity1
#   from a gene list-derived AddModuleScore
#
# Input example:
#   combined.obj.ls_new$Dam_bb
#
# Analytical premise:
# - combined.obj.ls_new is a merged Seurat object containing RA0, RA05, and RA25.
# - The analysis uses combined.obj.ls_new$Dam_bb.
# - Dam_1_0_UPscore_GeneActivity1 is defined from a gene list (not enhancer
#   coordinates) using AddModuleScore on the RNA.Dam assay.
# - In the original analysis context, this score represents a jump-associated
#   gene-activity program.
#
# Goal:
# - Compute Dam_1_0_UPscore_GeneActivity1 from an input gene list.
# - Plot violin + boxplot across RA0 / RA05 / RA25.
# - Test the planned one-sided hypotheses:
#     RA05 > RA0
#     RA25 > RA0
#   using a mixed model with random intercept = sample.
#
# Required metadata in Dam_bb:
#   RA
#   sample   (or orig.ident as fallback)
#   nCount_RNA.Dam  (optional covariate)
# ============================================================

library(Seurat)
library(dplyr)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)

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

read_gene_list <- function(gene_list_path, format = c("csv", "txt")) {
  format <- match.arg(format)

  if (format == "csv") {
    df <- read.csv(gene_list_path, row.names = 1, check.names = FALSE)
    genes <- rownames(df)
  } else {
    genes <- scan(gene_list_path, what = character(), sep = "\n", quiet = TRUE)
  }

  genes <- unique(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (length(genes) == 0) {
    stop("No genes were found in the input gene list.", call. = FALSE)
  }

  genes
}

add_dam_1_0_upscore_from_gene_list <- function(
    dam_obj,
    gene_list_path,
    gene_list_format = c("csv", "txt"),
    assay = "RNA.Dam",
    score_name = "Dam_1_0_UPscore_GeneActivity") {

  if (!inherits(dam_obj, "Seurat")) {
    stop("dam_obj must be a Seurat object.", call. = FALSE)
  }
  if (!assay %in% names(dam_obj@assays)) {
    stop("Assay '", assay, "' was not found in the Seurat object.", call. = FALSE)
  }

  genes <- read_gene_list(gene_list_path, format = match.arg(gene_list_format))

  DefaultAssay(dam_obj) <- assay
  genes_present <- intersect(genes, rownames(dam_obj))

  if (length(genes_present) == 0) {
    stop("None of the genes from the input list were found in the assay.", call. = FALSE)
  }

  dam_obj <- AddModuleScore(
    object = dam_obj,
    features = list(genes_present),
    name = score_name
  )

  list(
    dam_obj = dam_obj,
    genes_input = genes,
    genes_present = genes_present,
    score_column = paste0(score_name, "1")
  )
}

fit_ra_mixed_model <- function(
    dam_obj,
    score_var = "Dam_1_0_UPscore_GeneActivity1",
    group_var = "RA",
    sample_var = "sample",
    umi_var = "nCount_RNA.Dam",
    selected_groups = c("RA0", "RA05", "RA25"),
    use_covariate_log10umi = TRUE) {

  check_required_columns(dam_obj, c(group_var, score_var))

  if (!(sample_var %in% colnames(dam_obj@meta.data))) {
    message("'", sample_var, "' not found. Using 'orig.ident' instead.")
    sample_var <- "orig.ident"
  }

  vars_need <- unique(c(score_var, group_var, sample_var, umi_var))
  vars_need <- vars_need[vars_need %in% colnames(dam_obj@meta.data)]

  df_cell <- FetchData(dam_obj, vars = vars_need)

  col_map <- setNames(colnames(df_cell), colnames(df_cell))
  col_map[score_var] <- "score"
  col_map[group_var] <- "RA"
  col_map[sample_var] <- "sample"
  if (umi_var %in% colnames(df_cell)) col_map[umi_var] <- "nCount"
  colnames(df_cell) <- unname(col_map[colnames(df_cell)])

  if (!("nCount" %in% colnames(df_cell))) df_cell$nCount <- NA_real_

  df_cell <- df_cell %>%
    filter(
      !is.na(score),
      !is.na(RA),
      !is.na(sample),
      RA %in% selected_groups
    ) %>%
    mutate(
      RA = factor(RA, levels = selected_groups),
      sample = factor(sample),
      log10_umi = ifelse(is.na(nCount), NA_real_, log10(nCount + 1))
    )

  if (use_covariate_log10umi && all(!is.na(df_cell$log10_umi))) {
    fit <- lmer(score ~ RA + log10_umi + (1 | sample), data = df_cell, REML = FALSE)
    model_label <- "score ~ RA + log10_umi + (1 | sample)"
  } else {
    fit <- lmer(score ~ RA + (1 | sample), data = df_cell, REML = FALSE)
    model_label <- "score ~ RA + (1 | sample)"
  }

  anova_res <- anova(fit)
  p_global <- anova_res["RA", "Pr(>F)"]

  emm <- emmeans(fit, ~ RA)
  contr_vs_RA0 <- contrast(emm, method = "trt.vs.ctrl", ref = 1)

  two_sided_BH <- as.data.frame(summary(contr_vs_RA0, infer = TRUE, adjust = "BH"))
  tmp <- as.data.frame(summary(contr_vs_RA0, infer = TRUE, adjust = "none"))

  ratio_col <- intersect(c("t.ratio", "z.ratio"), colnames(tmp))[1]
  if (is.na(ratio_col)) stop("No t.ratio/z.ratio column found in contrast output.")

  if (ratio_col == "z.ratio") {
    tmp$p.one.sided <- pnorm(tmp[[ratio_col]], lower.tail = FALSE)
  } else {
    tmp$p.one.sided <- with(tmp, pt(.data[[ratio_col]], df = df, lower.tail = FALSE))
  }

  tmp$p.adj.BH.one.sided <- p.adjust(tmp$p.one.sided, method = "BH")
  tmp$label <- sprintf("BH one-sided p = %.3g", tmp$p.adj.BH.one.sided)

  list(
    df_cell = df_cell,
    fit = fit,
    model_label = model_label,
    anova_res = anova_res,
    p_global = p_global,
    contrasts = tmp,
    contrasts_two_sided_BH = two_sided_BH
  )
}

make_ra_violin_plot <- function(df_cell, contrasts_df, p_global, model_label) {
  parts <- strsplit(as.character(contrasts_df$contrast), " - ")
  treat <- vapply(parts, function(x) x[1], character(1))
  ctrl  <- vapply(parts, function(x) x[2], character(1))

  yr <- range(df_cell$score, na.rm = TRUE)
  dY <- diff(yr)
  if (!is.finite(dY) || dY == 0) dY <- max(abs(yr), na.rm = TRUE) + 1e-6
  y0 <- yr[2] + 0.06 * dY
  ystep <- 0.08 * dY

  pval_df <- data.frame(
    group1 = ctrl,
    group2 = treat,
    y.position = y0 + (seq_along(treat) - 1) * ystep,
    label = contrasts_df$label,
    stringsAsFactors = FALSE
  )

  ggviolin(
    df_cell,
    x = "RA", y = "score",
    trim = FALSE,
    add = "boxplot",
    add.params = list(width = 0.35, outlier.shape = NA)
  ) +
    theme_classic(base_size = 11) +
    ylab("Jump-associated gene-activity score") +
    xlab("RA") +
    labs(
      title = "Dam_1_0_UPscore_GeneActivity1 across RA",
      subtitle = paste0(
        "Mixed model: ", model_label,
        " | Global RA effect p=", signif(p_global, 3),
        " | One-sided: H1 (RA > RA0), BH adjusted"
      )
    ) +
    stat_pvalue_manual(
      pval_df,
      label = "label",
      tip.length = 0.01
    )
}

run_dam_1_0_upscore_ra_analysis <- function(
    dam_obj,
    gene_list_path,
    gene_list_format = c("csv", "txt"),
    out_pdf = "Dam_1_0_UPscore_violin_box_mixed_oneSided_RA0ref.pdf",
    out_csv = "Dam_1_0_UPscore_mixed_oneSided_contrasts_RA0ref.csv",
    assay = "RNA.Dam",
    score_name = "Dam_1_0_UPscore_GeneActivity",
    group_var = "RA",
    sample_var = "sample",
    umi_var = "nCount_RNA.Dam",
    selected_groups = c("RA0", "RA05", "RA25"),
    use_covariate_log10umi = TRUE) {

  score_res <- add_dam_1_0_upscore_from_gene_list(
    dam_obj = dam_obj,
    gene_list_path = gene_list_path,
    gene_list_format = match.arg(gene_list_format),
    assay = assay,
    score_name = score_name
  )

  model_res <- fit_ra_mixed_model(
    dam_obj = score_res$dam_obj,
    score_var = score_res$score_column,
    group_var = group_var,
    sample_var = sample_var,
    umi_var = umi_var,
    selected_groups = selected_groups,
    use_covariate_log10umi = use_covariate_log10umi
  )

  plot_obj <- make_ra_violin_plot(
    df_cell = model_res$df_cell,
    contrasts_df = model_res$contrasts,
    p_global = model_res$p_global,
    model_label = model_res$model_label
  )

  write.csv(model_res$contrasts, out_csv, row.names = FALSE, quote = FALSE)

  pdf(out_pdf, width = 4.2, height = 5.2)
  print(plot_obj)
  dev.off()

  list(
    dam_obj = score_res$dam_obj,
    genes_input = score_res$genes_input,
    genes_present = score_res$genes_present,
    score_column = score_res$score_column,
    plot = plot_obj,
    fit = model_res$fit,
    contrasts = model_res$contrasts,
    contrasts_two_sided_BH = model_res$contrasts_two_sided_BH,
    anova_res = model_res$anova_res,
    out_pdf = out_pdf,
    out_csv = out_csv
  )
}

# ============================================================
# Example
# ============================================================
# source("dam_1_0_upscore_from_gene_list_violin_mixed_RA_minimal.R")
#
# combined.obj.ls_new <- readRDS("path_to_RA_object/RA0_RA05_RA25_merged.rds")
#
# res <- run_dam_1_0_upscore_ra_analysis(
#   dam_obj = combined.obj.ls_new$Dam_bb,
#   gene_list_path = "path_to_gene_list/FindMarkers_Dam_1_0_UP.csv",
#   gene_list_format = "csv",
#   out_pdf = "Dam_1_0_UPscore_violin_box_mixed_oneSided_RA0ref.pdf",
#   out_csv = "Dam_1_0_UPscore_mixed_oneSided_contrasts_RA0ref.csv"
# )
#
# print(res$plot)
# head(res$genes_present)
# res$contrasts
