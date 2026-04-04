# ============================================================
# Identity jump GLMM for U0126 treatment
#
# Input example:
#   combined.obj.ls_U0126$Dam_bb
#
# Comparison type:
# - Treatment vs control within HisTrac
#
# Biological question:
# - Does U0126 treatment increase identity-jump frequency relative to control?
#
# Analytical premise:
# - The combined object contains two modalities, Dam_bb and K27ac_aa.
# - Dam_bb metadata already stores both modality-derived identity calls:
#     Dam_annotation_new
#     K27ac_annotation_new
# - Therefore, identity jump can be evaluated directly from Dam_bb alone as:
#     jump = 1(Dam_annotation_new != K27ac_annotation_new)
# - K27ac_aa is not required for the GLMM itself, because Dam_bb already carries
#   both annotations for the same cells.
#
# Dataset-specific note:
# - This is not a Snapshot vs HisTrac comparison.
# - Both control and U0126 groups belong to the HisTrac framework.
#
# Group definition:
# - Control:
#     ctrlrep1, ctrlrep2
# - Treatment:
#     U0126rep1, U0126rep2, U0126rep3
#
# Model:
# - Binomial GLMM
# - Formula:
#     jump ~ group + (1|sample)
# - Random intercept:
#     sample
# - Planned one-sided hypothesis:
#     U0126 > control
#
# Output:
# - Per-replicate jump rates
# - Estimated jump probabilities for control and U0126
# - Odds ratio for U0126 vs control
# ============================================================

library(Seurat)
library(dplyr)
library(lme4)
library(emmeans)

# --------------------
# Input
# --------------------
combined.obj.ls_U0126 <- readRDS("path_to_U0126/U0126_ctrl_U0126.rds")
obj <- combined.obj.ls_U0126$Dam_bb

# --------------------
# Column names
# --------------------
dam_id_var <- "Dam_annotation_new"
k27_id_var <- "K27ac_annotation_new"
rep_var    <- "sample"
id_col     <- "orig.ident"

# --------------------
# Manual group assignment
# --------------------
control_ids <- c("ctrlrep1", "ctrlrep2")
treat_ids   <- c("U0126rep1", "U0126rep2", "U0126rep3")

# --------------------
# Output files
# --------------------
out_csv_rates <- "IdentityJump_rates_U0126.csv"
out_csv_glmm  <- "IdentityJump_GLMM_U0126.csv"

# --------------------
# Checks
# --------------------
meta <- obj@meta.data
stopifnot(all(c(dam_id_var, k27_id_var, rep_var, id_col) %in% colnames(meta)))

# --------------------
# Build dataframe
# --------------------
df <- meta %>%
  transmute(
    cell   = rownames(meta),
    id     = as.character(.data[[id_col]]),
    sample = factor(.data[[rep_var]]),
    dam_id = .data[[dam_id_var]],
    k27_id = .data[[k27_id_var]]
  ) %>%
  filter(!is.na(id), !is.na(sample), !is.na(dam_id), !is.na(k27_id)) %>%
  mutate(
    group = case_when(
      id %in% control_ids ~ "control",
      id %in% treat_ids   ~ "U0126",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  mutate(
    group = factor(group, levels = c("control", "U0126")),
    jump  = as.integer(dam_id != k27_id)
  )

print(table(df$group))
cat("Replicates:", nlevels(df$sample), "\n")

if (nlevels(df$group) < 2) {
  stop("group has fewer than 2 levels. Check sample IDs.")
}
if (nlevels(df$sample) < 2) {
  stop("sample has fewer than 2 levels. Cannot fit (1|sample).")
}

# --------------------
# Per-replicate jump rates
# --------------------
rate_df <- df %>%
  group_by(group, sample) %>%
  summarise(
    n_cells   = n(),
    n_jump    = sum(jump),
    jump_rate = n_jump / n_cells,
    .groups = "drop"
  )

write.csv(rate_df, out_csv_rates, row.names = FALSE, quote = FALSE)

# --------------------
# GLMM
# --------------------
fit <- glmer(
  jump ~ group + (1 | sample),
  data = df,
  family = binomial(),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

prob_df <- as.data.frame(summary(emmeans(fit, ~ group), type = "response"))

emm_link <- emmeans(fit, ~ group)
contr <- contrast(emm_link, method = "trt.vs.ctrl", ref = 1)
s <- as.data.frame(summary(contr, infer = TRUE, adjust = "none"))

z <- s$z.ratio[1]
p_one <- pnorm(z, lower.tail = FALSE)  # one-sided: U0126 > control

logOR <- s$estimate[1]
res_df <- data.frame(
  dataset = "U0126",
  p_jump_control = prob_df$prob[prob_df$group == "control"][1],
  p_jump_U0126    = prob_df$prob[prob_df$group == "U0126"][1],
  logOR_U0126_vs_control = logOR,
  OR_U0126_vs_control    = exp(logOR),
  OR_LCL = exp(s$asymp.LCL[1]),
  OR_UCL = exp(s$asymp.UCL[1]),
  p_two_sided = s$p.value[1],
  p_one_sided = p_one,
  stringsAsFactors = FALSE
)

print(res_df)
write.csv(res_df, out_csv_glmm, row.names = FALSE, quote = FALSE)

cat("\nSaved:\n")
cat("  ", out_csv_glmm, "\n", sep = "")
cat("  ", out_csv_rates, "\n", sep = "")
