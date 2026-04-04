# ============================================================
# Identity jump GLMM for FreeDam
#
# Input example:
#   combined.obj.ls_FreeDam$Dam_bb
#
# Comparison type:
# - Snapshot vs HisTrac
#
# Biological question:
# - Is the identity-jump frequency higher in Day7R HisTrac cells than in
#   Snapshot cells (Day2 + Day7) in the FreeDam dataset?
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
# Group definition:
# - Snapshot:
#     Day2rep1, Day2rep2, Day7rep1
# - HisTrac:
#     Day7Rrep1, Day7Rrep2
#
# Model:
# - Binomial GLMM
# - Formula:
#     jump ~ group + (1|sample)
# - Random intercept:
#     sample
# - Planned one-sided hypothesis:
#     HisTrac > Snapshot
#
# Output:
# - Per-replicate jump rates
# - Estimated jump probabilities for Snapshot and HisTrac
# - Odds ratio for HisTrac vs Snapshot
# Note:
# - In this dataset, Day7rep2 is not included in the Snapshot group.
# ============================================================

library(Seurat)
library(dplyr)
library(lme4)
library(emmeans)

# --------------------
# Input
# --------------------
combined.obj.ls_FreeDam <- readRDS("path_to_FreeDam/Dam_Leo1_D2_D7_D7R.rds")
obj <- combined.obj.ls_FreeDam$Dam_bb

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
his_ids  <- c("Day7Rrep1", "Day7Rrep2")
snap_ids <- c("Day2rep1", "Day2rep2", "Day7rep1")

# --------------------
# Output files
# --------------------
out_csv_rates <- "IdentityJump_rates_FreeDam.csv"
out_csv_glmm  <- "IdentityJump_GLMM_FreeDam.csv"

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
      id %in% snap_ids ~ "Snapshot",
      id %in% his_ids  ~ "HisTrac",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  mutate(
    group = factor(group, levels = c("Snapshot", "HisTrac")),
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

# Predicted probabilities
prob_df <- as.data.frame(summary(emmeans(fit, ~ group), type = "response"))

# Contrast: HisTrac - Snapshot
emm_link <- emmeans(fit, ~ group)
contr <- contrast(emm_link, method = "trt.vs.ctrl", ref = 1)
s <- as.data.frame(summary(contr, infer = TRUE, adjust = "none"))

z <- s$z.ratio[1]
p_one <- pnorm(z, lower.tail = FALSE)  # one-sided: HisTrac > Snapshot

logOR <- s$estimate[1]
res_df <- data.frame(
  dataset = "FreeDam",
  p_jump_Snapshot = prob_df$prob[prob_df$group == "Snapshot"][1],
  p_jump_HisTrac  = prob_df$prob[prob_df$group == "HisTrac"][1],
  logOR_HisTrac_vs_Snapshot = logOR,
  OR_HisTrac_vs_Snapshot    = exp(logOR),
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
