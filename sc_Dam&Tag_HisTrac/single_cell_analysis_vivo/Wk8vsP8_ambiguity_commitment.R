# ============================================================
# P8 vs W8 comparison for 3 metrics
#
# Assumptions:
# - obj is already restricted to the target cells of interest
#   (e.g. K27ac P8/W8, L2-3 vs L4-5 only)
# - obj@meta.data contains:
#     timepoint            : "P8" / "W8"
#     sample or orig.ident : replicate ID
#     knn_entropy01_01
#     knn_ambiguity_01
#     commit_strength
#     nCount_peaks         : optional depth covariate
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
# Input object
# -------------------------------------------------------------------
K27ac_aa <- subset(combined.obj.ls3$K27ac_aa, subset = stage %in% c("P8", "W8"))
obj <- subset(K27ac_aa, subset = Cell_Type %in% c("L2-3", "L4-5"))

# -------------------------------------------------------------------
# Basic metadata
# -------------------------------------------------------------------
obj$timepoint <- factor(as.character(obj$stage), levels = c("P8", "W8"))
obj$anno01 <- factor(as.character(obj$Cell_Type), levels = c("L2-3", "L4-5"))

# replicate column
sample_var <- if ("sample" %in% colnames(obj@meta.data)) "sample" else "orig.ident"

# -------------------------------------------------------------------
# Compute ambiguity / commitment in K27ac LSI space
# -------------------------------------------------------------------
E <- Embeddings(obj, "lsi")
dmax <- ncol(E)
dims_use <- 2:min(30, dmax)   # drop LSI1
X_all <- as.matrix(E[, dims_use, drop = FALSE])

# helper
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
  
  # 2 classes -> normalize by log(2)
  entropy01 <- entropy / log(2)
  
  list(
    nn = nn,
    purity = purity,
    ambiguity = ambiguity,
    entropy01 = entropy01
  )
}

# compute within each timepoint separately
obj$knn_entropy01_01 <- NA_real_
obj$knn_ambiguity_01 <- NA_real_
obj$p_anno1 <- NA_real_
obj$commit_strength <- NA_real_

for (tp in levels(obj$timepoint)) {
  obj_t <- subset(obj, subset = timepoint == tp & anno01 %in% c("L2-3", "L4-5"))
  
  lab <- factor(obj_t$anno01, levels = c("L2-3", "L4-5"))
  Xt  <- as.matrix(Embeddings(obj_t, "lsi")[, dims_use, drop = FALSE])
  
  # ambiguity
  knn <- knn_mixing_and_ambiguity(Xt, lab, k = 30)
  
  # commitment
  y <- ifelse(lab == "L4-5", 1, 0)
  set.seed(1)
  cv <- cv.glmnet(Xt, y, family = "binomial", alpha = 1, nfolds = 5)
  
  p1 <- as.numeric(predict(cv, newx = Xt, s = "lambda.1se", type = "response"))
  p1 <- pmin(pmax(p1, 1e-6), 1 - 1e-6)
  commit <- abs(qlogis(p1))
  
  cells <- colnames(obj_t)
  obj@meta.data[cells, "knn_entropy01_01"] <- knn$entropy01
  obj@meta.data[cells, "knn_ambiguity_01"] <- knn$ambiguity
  obj@meta.data[cells, "p_anno1"]          <- p1
  obj@meta.data[cells, "commit_strength"]  <- commit
}

# sanity check
table(obj$timepoint)
summary(obj$knn_entropy01_01)
summary(obj$knn_ambiguity_01)
summary(obj$commit_strength)

# -------------------------------------------------------------------
# Settings for mixed models
# -------------------------------------------------------------------
time_var <- "timepoint"

metrics <- c(
  entropy    = "knn_entropy01_01",
  ambiguity  = "knn_ambiguity_01",
  commitment = "commit_strength"
)

cov_var <- if ("nCount_peaks" %in% colnames(obj@meta.data)) "nCount_peaks" else NA_character_
use_covariate_log10depth <- !is.na(cov_var)

eps <- 1e-4

# -------------------------------------------------------------------
# Fetch cell-level dataframe
# -------------------------------------------------------------------
vars_need <- c(time_var, sample_var, unname(metrics), cov_var)
vars_need <- vars_need[!is.na(vars_need)]

df <- FetchData(obj, vars = vars_need)

# rename only the columns that actually exist
colnames(df)[colnames(df) == time_var]   <- "timepoint"
colnames(df)[colnames(df) == sample_var] <- "sample"
if (!is.na(cov_var)) {
  colnames(df)[colnames(df) == cov_var] <- "nCount_peaks"
}

df <- df %>%
  mutate(
    timepoint = factor(as.character(timepoint), levels = c("P8", "W8")),
    sample    = factor(as.character(sample))
  ) %>%
  filter(!is.na(timepoint), !is.na(sample))

if ("nCount_peaks" %in% colnames(df)) {
  df <- df %>% mutate(log10_depth = log10(nCount_peaks + 1))
} else {
  df$log10_depth <- NA_real_
  use_covariate_log10depth <- FALSE
}

# -------------------------------------------------------------------
# Helper: mixed model for W8 - P8
# -------------------------------------------------------------------
fit_lmm_W8_vs_P8 <- function(df, y_col,
                             transform = c("none", "log1p", "logit01"),
                             use_cov = FALSE, eps = 1e-4) {
  transform <- match.arg(transform)
  
  dat <- df %>%
    transmute(
      timepoint   = .data$timepoint,
      sample      = .data$sample,
      log10_depth = .data$log10_depth,
      y           = .data[[y_col]]
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
  contr <- contrast(emm, method = list("W8 - P8" = c(-1, 1)))
  as.data.frame(summary(contr, infer = TRUE, adjust = "none"))
}

# -------------------------------------------------------------------
# Run the 3 tests
# -------------------------------------------------------------------
res_entropy <- fit_lmm_W8_vs_P8(
  df, metrics[["entropy"]],
  transform = "logit01",
  use_cov = use_covariate_log10depth,
  eps = eps
)

res_ambig <- fit_lmm_W8_vs_P8(
  df, metrics[["ambiguity"]],
  transform = "logit01",
  use_cov = use_covariate_log10depth,
  eps = eps
)

res_commit <- fit_lmm_W8_vs_P8(
  df, metrics[["commitment"]],
  transform = "log1p",
  use_cov = use_covariate_log10depth,
  eps = eps
)

res_all <- bind_rows(
  cbind(metric = "entropy",    res_entropy),
  cbind(metric = "ambiguity",  res_ambig),
  cbind(metric = "commitment", res_commit)
) %>%
  mutate(
    p.adj.BH = p.adjust(p.value, method = "BH"),
    label = sprintf("LMM two-sided\nBH p = %.3g", p.adj.BH)
  )

print(res_all)
write.csv(res_all, "P8vsW8_mixedModel_3metrics.csv", row.names = FALSE, quote = FALSE)

# -------------------------------------------------------------------
# Plot helpers
# -------------------------------------------------------------------
make_pval_df <- function(df, y_col, label_text) {
  yr <- range(df[[y_col]], na.rm = TRUE)
  dY <- diff(yr)
  if (!is.finite(dY) || dY == 0) dY <- max(abs(yr), na.rm = TRUE) + 1e-6
  
  data.frame(
    group1 = "P8",
    group2 = "W8",
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
# Build final figure
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