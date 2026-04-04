# ============================================================
# A395 analysis for Dam_1_0_UPscore_GeneActivity1
#
# Input example:
#   combined.obj.ls_new$Dam_bb
#
# Analytical premise:
# - combined.obj.ls_new is a merged Seurat object containing the A395-series
#   samples together (ctrl and A395).
# - The combined object contains two modalities, Dam_bb and K27ac_aa.
# - This script uses only Dam_bb.
# - Dam_1_0_UPscore_GeneActivity1 is computed from a gene list using
#   AddModuleScore() on the RNA.Dam assay.
# - The gene list corresponds to the Dam 1->0 signature used in the original
#   Day2/Day7/Day7R analysis.
#
# Statistical design:
# - This is a 2-group comparison (ctrl vs A395), unlike the RA analysis which
#   contains three groups and two planned contrasts.
# - The score is modeled with a linear mixed model:
#     score ~ cond + log10_umi + (1 | sample)
#   or, if UMI is unavailable,
#     score ~ cond + (1 | sample)
# - Planned one-sided hypothesis:
#     A395 > ctrl
#
# Output:
# - PDF:
#     Dam_1_0_UPscore_violin_box_mixed_oneSided_ctrlRef.pdf
# - CSV:
#     Dam_1_0_UPscore_mixed_oneSided_contrast_ctrlRef.csv
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)

read_gene_list <- function(gene_list_path, gene_list_format = c("csv", "txt")) {
  gene_list_format <- match.arg(gene_list_format)

  if (gene_list_format == "csv") {
    tab <- read.csv(gene_list_path, row.names = 1, check.names = FALSE)
    genes <- rownames(tab)
  } else {
    genes <- readLines(gene_list_path)
  }

  genes <- unique(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (length(genes) == 0) {
    stop("No genes were read from the provided gene list.", call. = FALSE)
  }

  genes
}

add_dam_1_0_module_score <- function(
    dam_obj,
    gene_list_path,
    gene_list_format = c("csv", "txt"),
    assay = "RNA.Dam",
    score_name = "Dam_1_0_UPscore_GeneActivity") {

  gene_list_format <- match.arg(gene_list_format)

  if (!inherits(dam_obj, "Seurat")) {
    stop("dam_obj must be a Seurat object.", call. = FALSE)
  }
  if (!assay %in% names(dam_obj@assays)) {
    stop("Assay '", assay, "' was not found in dam_obj.", call. = FALSE)
  }

  DefaultAssay(dam_obj) <- assay

  genes <- read_gene_list(
    gene_list_path = gene_list_path,
    gene_list_format = gene_list_format
  )

  genes_present <- intersect(genes, rownames(dam_obj))
  if (length(genes_present) == 0) {
    stop("None of the genes from the gene list were found in the assay.", call. = FALSE)
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
    score_var = paste0(score_name, "1")
  )
}

run_dam_1_0_upscore_a395_analysis <- function(
    dam_obj,
    gene_list_path,
    gene_list_format = c("csv", "txt"),
    assay = "RNA.Dam",
    group_var = "cond",
    sample_var = "sample",
    umi_var = "nCount_RNA.Dam",
    selected_groups = c("ctrl", "A395"),
    use_covariate_log10umi = TRUE,
    out_pdf = "Dam_1_0_UPscore_violin_box_mixed_oneSided_ctrlRef.pdf",
    out_csv = "Dam_1_0_UPscore_mixed_oneSided_contrast_ctrlRef.csv",
    box_width = 0.40) {

  gene_list_format <- match.arg(gene_list_format)

  score_res <- add_dam_1_0_module_score(
    dam_obj = dam_obj,
    gene_list_path = gene_list_path,
    gene_list_format = gene_list_format,
    assay = assay,
    score_name = "Dam_1_0_UPscore_GeneActivity"
  )

  dam_obj <- score_res$dam_obj
  score_var <- score_res$score_var

  if (!(group_var %in% colnames(dam_obj@meta.data))) {
    stop("group_var not found in meta.data: ", group_var, call. = FALSE)
  }

  if (!(sample_var %in% colnames(dam_obj@meta.data))) {
    message("'", sample_var, "' not found in meta.data. Using 'orig.ident' instead.")
    sample_var <- "orig.ident"
  }

  vars_need <- unique(c(score_var, group_var, sample_var, umi_var))
  vars_need <- vars_need[vars_need %in% colnames(dam_obj@meta.data)]

  df_cell <- FetchData(dam_obj, vars = vars_need)

  col_map <- setNames(colnames(df_cell), colnames(df_cell))
  col_map[score_var]  <- "score"
  col_map[group_var]  <- "cond"
  col_map[sample_var] <- "sample"
  if (umi_var %in% colnames(df_cell)) col_map[umi_var] <- "nCount"
  colnames(df_cell) <- unname(col_map[colnames(df_cell)])

  if (!("nCount" %in% colnames(df_cell))) df_cell$nCount <- NA_real_

  df_cell <- df_cell %>%
    filter(
      !is.na(score),
      !is.na(cond),
      !is.na(sample),
      cond %in% selected_groups
    ) %>%
    mutate(
      cond = factor(cond, levels = selected_groups),
      sample = factor(sample),
      log10_umi = ifelse(is.na(nCount), NA_real_, log10(nCount + 1))
    )

  if (use_covariate_log10umi) {
    df_cell <- df_cell %>% filter(!is.na(log10_umi))
  }

  if (nlevels(df_cell$cond) < 2) {
    stop("The filtered data contain fewer than 2 condition levels.", call. = FALSE)
  }
  if (nlevels(df_cell$sample) < 2) {
    stop("The filtered data contain fewer than 2 sample levels.", call. = FALSE)
  }

  if (use_covariate_log10umi) {
    fit <- lmer(score ~ cond + log10_umi + (1 | sample), data = df_cell, REML = FALSE)
    model_label <- "score ~ cond + log10_umi + (1 | sample)"
  } else {
    fit <- lmer(score ~ cond + (1 | sample), data = df_cell, REML = FALSE)
    model_label <- "score ~ cond + (1 | sample)"
  }

  anova_res <- anova(fit)
  p_global <- anova_res["cond", "Pr(>F)"]

  emm <- emmeans(fit, ~ cond)
  contr <- contrast(emm, method = "trt.vs.ctrl", ref = 1)
  tmp <- as.data.frame(summary(contr, infer = TRUE, adjust = "none"))

  ratio_col <- intersect(c("t.ratio", "z.ratio"), colnames(tmp))[1]
  if (is.na(ratio_col)) stop("No t.ratio/z.ratio column found in contrast output.", call. = FALSE)

  if (ratio_col == "z.ratio") {
    tmp$p.one.sided <- pnorm(tmp[[ratio_col]], lower.tail = FALSE)
  } else {
    tmp$p.one.sided <- with(tmp, pt(tmp[[ratio_col]], df = df, lower.tail = FALSE))
  }

  # Only one planned contrast here; BH is numerically identical but retained for consistency.
  tmp$p.adj.BH.one.sided <- p.adjust(tmp$p.one.sided, method = "BH")
  tmp$label <- sprintf("one-sided (A395 > ctrl), BH p = %.3g", tmp$p.adj.BH.one.sided)

  write.csv(tmp, out_csv, row.names = FALSE, quote = FALSE)

  parts <- strsplit(as.character(tmp$contrast), " - ")
  treat <- vapply(parts, function(x) x[1], character(1))
  ctrl  <- vapply(parts, function(x) x[2], character(1))

  yr <- range(df_cell$score, na.rm = TRUE)
  dY <- diff(yr)
  if (!is.finite(dY) || dY == 0) dY <- max(abs(yr), na.rm = TRUE) + 1e-6

  pval_df <- data.frame(
    group1 = ctrl,
    group2 = treat,
    y.position = yr[2] + 0.08 * dY,
    label = tmp$label,
    stringsAsFactors = FALSE
  )

  p <- ggviolin(
    df_cell,
    x = "cond", y = "score",
    trim = FALSE,
    add = "boxplot",
    add.params = list(
      width = box_width,
      outlier.shape = NA
    )
  ) +
    theme_classic(base_size = 11) +
    ylab("Module Score (Dam_1_0_UP)") +
    xlab("cond") +
    labs(
      title = "Dam_1_0_UP module score: ctrl vs A395",
      subtitle = paste0(
        "Mixed model: ", model_label,
        " | Global cond effect p=", signif(p_global, 3),
        " | One-sided: H1 (A395 > ctrl)"
      )
    ) +
    stat_pvalue_manual(
      pval_df,
      label = "label",
      tip.length = 0.01
    )

  pdf(out_pdf, width = 4.2, height = 5.2)
  print(p)
  dev.off()

  cat("\n--- Mixed model ANOVA (global cond effect) ---\n")
  print(anova_res)

  cat("\n--- One-sided contrast (H1: A395 > ctrl) ---\n")
  print(tmp)

  cat("\nSaved:\n  PDF: ", out_pdf, "\n  CSV: ", out_csv, "\n", sep = "")

  invisible(list(
    dam_obj = dam_obj,
    genes_input = score_res$genes_input,
    genes_present = score_res$genes_present,
    score_var = score_var,
    df_cell = df_cell,
    fit = fit,
    anova = anova_res,
    contrasts = tmp,
    plot = p
  ))
}

# ============================================================
# Example
# ============================================================
# source("dam_1_0_upscore_from_gene_list_violin_mixed_A395_minimal.R")
#
# combined.obj.ls_new <- readRDS("path_to_A395_object/ctrl_A395_merged.rds")
#
# res <- run_dam_1_0_upscore_a395_analysis(
#   dam_obj = combined.obj.ls_new$Dam_bb,
#   gene_list_path = "path_to_gene_list/FindMarkers_Dam_1_0_UP.csv",
#   gene_list_format = "csv",
#   out_pdf = "Dam_1_0_UPscore_violin_box_mixed_oneSided_ctrlRef.pdf",
#   out_csv = "Dam_1_0_UPscore_mixed_oneSided_contrast_ctrlRef.csv"
# )
#
# print(res$plot)
# res$contrasts
