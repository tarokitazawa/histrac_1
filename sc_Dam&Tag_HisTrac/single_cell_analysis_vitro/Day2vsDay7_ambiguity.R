# ============================================================
# ambiguity analysis
#
# Analytical premise:
# - combined.obj.ls_new is the merged in vitro object containing Day2, Day7,
#   and Day7R cells.
# - A shared LSI was constructed from this merged object and is used as the
#   common low-dimensional space for downstream analyses.
#
# Modality premise:
# - Dam and K27ac were measured for the same cells.
# - Modality-specific identities were assigned independently, producing:
#     Dam_annotation_new
#     K27ac_annotation_new
# - Dam_annotation_new and K27ac_annotation_new contain only 0 / 1
#   (for example, id-C vs id-P).
# - These identity calls are already stored in the merged Seurat objects,
#   allowing transition and ambiguity analyses to be performed directly from
#   Dam_bb or K27ac_aa, depending on the analysis.
#
# Assumptions:
# - The chosen modality object already contains the shared LSI reduction.
# - Relevant subsets (for example, Day2 only, or Day2 vs Day7) are derived
#   from combined.obj.ls_new.
# - Dam_annotation_new and K27ac_annotation_new contain the identities of
#   interest for the downstream comparison.
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)
library(FNN)
library(glmnet)

# -------------------------------------------------------------------
# Input
# -------------------------------------------------------------------
obj <- combined.obj.ls_new$K27ac_aa

stopifnot(all(c("timepoint", "K27ac_annotation_new") %in% colnames(obj@meta.data)))
stopifnot("lsi" %in% Reductions(obj))

# replicate column
sample_var <- if ("sample" %in% colnames(obj@meta.data)) "sample" else "orig.ident"

# -------------------------------------------------------------------
# Basic setup
# -------------------------------------------------------------------
obj$timepoint <- factor(as.character(obj$timepoint), levels = c("Day2", "Day7"))
obj$anno01 <- factor(as.character(obj$K27ac_annotation_new), levels = c("0", "1"))
obj$grp4 <- factor(
  paste0(obj$timepoint, "_", obj$anno01),
  levels = c("Day2_0", "Day2_1", "Day7_0", "Day7_1")
)

# -------------------------------------------------------------------
# Use K27ac LSI coordinates
# -------------------------------------------------------------------
E <- Embeddings(obj, "lsi")
dims_use <- 2:min(30, ncol(E))   # skip LSI1
X_all <- as.matrix(E[, dims_use, drop = FALSE])

# -------------------------------------------------------------------
# Helper: kNN ambiguity
# -------------------------------------------------------------------
knn_mixing_and_ambiguity <- function(X, lab, k = 30) {
  lab <- as.character(lab)
  n <- nrow(X)
  k_use <- min(k, n - 1)
  
  nn <- FNN::get.knn(X, k = k_use)$nn.index
  
  purity <- vapply(seq_len(n), function(i) {
    mean(lab[nn[i, ]] == lab[i])
  }, numeric(1))
  
  ambiguity <- 1 - purity
  
  entropy <- vapply(seq_len(n), function(i) {
    p <- prop.table(table(lab[nn[i, ]]))
    -sum(p * log(p))
  }, numeric(1))
  
  entropy01 <- entropy / log(2)
  
  list(
    purity = purity,
    ambiguity = ambiguity,
    entropy01 = entropy01
  )
}

# -------------------------------------------------------------------
# Compute the 3 metrics within each timepoint
# -------------------------------------------------------------------
obj$knn_entropy01_01 <- NA_real_
obj$knn_ambiguity_01 <- NA_real_
obj$p_anno1 <- NA_real_
obj$commit_strength <- NA_real_

for (tp in levels(obj$timepoint)) {
  obj_t <- subset(obj, subset = timepoint == tp & anno01 %in% c("0", "1"))
  
  lab <- factor(obj_t$anno01, levels = c("0", "1"))
  Xt <- as.matrix(Embeddings(obj_t, "lsi")[, dims_use, drop = FALSE])
  
  # ambiguity
  knn <- knn_mixing_and_ambiguity(Xt, lab, k = 30)
  
  # commitment
  y <- ifelse(lab == "1", 1, 0)
  set.seed(1)
  cv <- cv.glmnet(Xt, y, family = "binomial", alpha = 1, nfolds = 5)
  
  p1 <- as.numeric(predict(cv, newx = Xt, s = "lambda.1se", type = "response"))
  p1 <- pmin(pmax(p1, 1e-6), 1 - 1e-6)
  commit <- abs(qlogis(p1))
  
  cells <- colnames(obj_t)
  obj@meta.data[cells, "knn_entropy01_01"] <- knn$entropy01
  obj@meta.data[cells, "knn_ambiguity_01"] <- knn$ambiguity
  obj@meta.data[cells, "p_anno1"] <- p1
  obj@meta.data[cells, "commit_strength"] <- commit
}

# -------------------------------------------------------------------
# Prepare dataframe for mixed models
# -------------------------------------------------------------------
metrics <- c(
  entropy    = "knn_entropy01_01",
  ambiguity  = "knn_ambiguity_01",
  commitment = "commit_strength"
)

cov_var <- if ("nCount_peaks" %in% colnames(obj@meta.data)) "nCount_peaks" else NA_character_
use_covariate_log10depth <- !is.na(cov_var)

vars_need <- c("timepoint", sample_var, unname(metrics), cov_var)
vars_need <- vars_need[!is.na(vars_need)]

df <- FetchData(obj, vars = vars_need)

colnames(df)[colnames(df) == sample_var] <- "sample"
if (!is.na(cov_var)) {
  colnames(df)[colnames(df) == cov_var] <- "nCount_peaks"
}

df <- df %>%
  mutate(
    timepoint = factor(as.character(timepoint), levels = c("Day2", "Day7")),
    sample = factor(as.character(sample))
  ) %>%
  filter(!is.na(timepoint), !is.na(sample))

if ("nCount_peaks" %in% colnames(df)) {
  df <- df %>% mutate(log10_depth = log10(nCount_peaks + 1))
} else {
  df$log10_depth <- NA_real_
  use_covariate_log10depth <- FALSE
}

# -------------------------------------------------------------------
# Helper: mixed model Day7 - Day2
# -------------------------------------------------------------------
fit_lmm_day7_vs_day2 <- function(df, y_col,
                                 transform = c("none", "log1p", "logit01"),
                                 use_cov = FALSE, eps = 1e-4) {
  transform <- match.arg(transform)
  
  dat <- df %>%
    transmute(
      timepoint = .data$timepoint,
      sample = .data$sample,
      log10_depth = .data$log10_depth,
      y = .data[[y_col]]
    ) %>%
    filter(!is.na(timepoint), !is.na(sample), !is.na(y))
  
  if (transform == "log1p") {
    dat$y_tr <- log1p(dat$y)
  } else if (transform == "logit01") {
    y2 <- pmin(pmax(dat$y, eps), 1 - eps)
    dat$y_tr <- qlogis(y2)
  } else {
    dat$y_tr <- dat$y
  }
  
  if (use_cov) {
    dat <- dat %>% filter(!is.na(log10_depth))
    fit <- lmer(y_tr ~ timepoint + log10_depth + (1 | sample), data = dat, REML = FALSE)
  } else {
    fit <- lmer(y_tr ~ timepoint + (1 | sample), data = dat, REML = FALSE)
  }
  
  emm <- emmeans(fit, ~ timepoint)
  contr <- contrast(emm, method = list("Day7 - Day2" = c(-1, 1)))
  as.data.frame(summary(contr, infer = TRUE, adjust = "none"))
}

# -------------------------------------------------------------------
# Run the 3 tests
# -------------------------------------------------------------------
eps <- 1e-4

res_entropy <- fit_lmm_day7_vs_day2(
  df, metrics[["entropy"]],
  transform = "logit01",
  use_cov = use_covariate_log10depth,
  eps = eps
)

res_ambig <- fit_lmm_day7_vs_day2(
  df, metrics[["ambiguity"]],
  transform = "logit01",
  use_cov = use_covariate_log10depth,
  eps = eps
)

res_commit <- fit_lmm_day7_vs_day2(
  df, metrics[["commitment"]],
  transform = "log1p",
  use_cov = use_covariate_log10depth,
  eps = eps
)

res_all <- bind_rows(
  cbind(metric = "entropy", res_entropy),
  cbind(metric = "ambiguity", res_ambig),
  cbind(metric = "commitment", res_commit)
) %>%
  mutate(
    p.adj.BH = p.adjust(p.value, method = "BH"),
    label = sprintf("LMM two-sided\nBH p = %.3g", p.adj.BH)
  )

write.csv(res_all, "Day2vsDay7_mixedModel_3metrics.csv", row.names = FALSE, quote = FALSE)

# -------------------------------------------------------------------
# Plot helpers
# -------------------------------------------------------------------
make_pval_df <- function(df, y_col, label_text) {
  yr <- range(df[[y_col]], na.rm = TRUE)
  dY <- diff(yr)
  if (!is.finite(dY) || dY == 0) dY <- max(abs(yr), na.rm = TRUE) + 1e-6
  
  data.frame(
    group1 = "Day2",
    group2 = "Day7",
    y.position = yr[2] + 0.08 * dY,
    label = label_text,
    stringsAsFactors = FALSE
  )
}

plot_box_tp <- function(df, y_col, ylab, p_lab) {
  pval_df <- make_pval_df(df, y_col, p_lab)
  
  ggplot(df, aes(x = timepoint, y = .data[[y_col]])) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic(base_size = 12) +
    labs(x = NULL, y = ylab) +
    stat_pvalue_manual(pval_df, label = "label", tip.length = 0.01)
}

# -------------------------------------------------------------------
# Final figure
# -------------------------------------------------------------------
p1 <- plot_box_tp(
  df,
  metrics["entropy"],
  "kNN entropy (0-1)",
  res_all$label[res_all$metric == "entropy"]
)

p2 <- plot_box_tp(
  df,
  metrics["ambiguity"],
  "kNN ambiguity (1-purity)",
  res_all$label[res_all$metric == "ambiguity"]
)

p3 <- plot_box_tp(
  df,
  metrics["commitment"],
  "Commitment strength",
  res_all$label[res_all$metric == "commitment"]
)

p_all <- p1 | p2 | p3
p_all

ggsave(
  "Day2vsDay7_3metrics_box_LMM.pdf",
  p_all,
  width = 10.5,
  height = 3.8,
  useDingbats = FALSE
)