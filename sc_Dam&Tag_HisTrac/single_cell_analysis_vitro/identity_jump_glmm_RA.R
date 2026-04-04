# ============================================================
# Identity jump GLMM for the RA treatment series
#
# Input example:
#   combined.obj.ls_RA$Dam_bb
#
# Comparison type:
# - Treatment vs control within HisTrac
#
# Biological question:
# - Does RA treatment increase identity-jump frequency relative to RA0 control?
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
# - All groups in this dataset belong to the HisTrac framework.
# - The merged object contains three RA conditions:
#     RA0, RA05, RA25
# - This dataset was normalized and merged within the RA series, and is analyzed
#   separately from the Day2/Day7/Day7R datasets.
#
# Group definition:
# - Control:
#     RA0rep1, RA0rep2
# - Treatment groups:
#     RA05rep1, RA05rep2
#     RA25rep1, RA25rep2
#
# Model:
# - Binomial GLMM
# - Formula:
#     jump ~ group + (1|sample)
# - Random intercept:
#     sample
# - Planned contrasts:
#     RA05 vs RA0
#     RA25 vs RA0
# - Planned one-sided hypothesis:
#     treatment > control
#
# Output:
# - Per-replicate jump rates
# - Estimated jump probabilities for each contrast
# - Odds ratios for RA05 vs RA0 and RA25 vs RA0
# ============================================================

library(Seurat)
library(dplyr)
library(lme4)
library(emmeans)

# --------------------
# Input
# --------------------
combined.obj.ls_RA <- readRDS("path_to_RA/RA_RA0_RA05_RA25.rds")
obj <- combined.obj.ls_RA$Dam_bb

# --------------------
# Column names
# --------------------
dam_id_var <- "Dam_annotation_new"
k27_id_var <- "K27ac_annotation_new"
rep_var    <- "sample"
id_col     <- "orig.ident"

# --------------------
# Manual assignment
# --------------------
control_name <- "RA0"
control_ids  <- c("RA0rep1", "RA0rep2")
targets <- list(
  RA05 = c("RA05rep1", "RA05rep2"),
  RA25 = c("RA25rep1", "RA25rep2")
)

# --------------------
# Output files
# --------------------
out_csv_glmm  <- "IdentityJump_GLMM_RA.csv"
out_csv_rates <- "IdentityJump_rates_RA.csv"

# --------------------
# Checks
# --------------------
meta <- obj@meta.data
stopifnot(all(c(dam_id_var, k27_id_var, rep_var, id_col) %in% colnames(meta)))

# --------------------
# Build base dataframe
# --------------------
df0 <- meta %>%
  transmute(
    cell   = rownames(meta),
    id     = as.character(.data[[id_col]]),
    sample = as.character(.data[[rep_var]]),
    dam_id = .data[[dam_id_var]],
    k27_id = .data[[k27_id_var]]
  ) %>%
  filter(!is.na(id), !is.na(sample), !is.na(dam_id), !is.na(k27_id)) %>%
  mutate(
    jump   = as.integer(dam_id != k27_id),
    sample = factor(sample)
  )

fit_one_target <- function(df0, control_ids, target_ids, target_name) {
  df <- df0 %>%
    mutate(
      group = case_when(
        id %in% control_ids ~ "control",
        id %in% target_ids  ~ "target",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group)) %>%
    mutate(group = factor(group, levels = c("control", "target")))

  tab <- table(df$group)
  if (any(tab == 0)) {
    stop("Empty group after filtering for target = ", target_name, ". Check IDs.")
  }
  if (nlevels(df$sample) < 2) {
    stop("sample has fewer than 2 levels after filtering for target = ", target_name)
  }

  rates <- df %>%
    group_by(group, sample) %>%
    summarise(
      n_cells   = n(),
      n_jump    = sum(jump),
      jump_rate = n_jump / n_cells,
      .groups = "drop"
    ) %>%
    mutate(target = target_name)

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
  p_one <- pnorm(z, lower.tail = FALSE)  # one-sided: target > control

  data.frame(
    dataset = "RA",
    target = target_name,
    p_jump_control = prob_df$prob[prob_df$group == "control"][1],
    p_jump_target  = prob_df$prob[prob_df$group == "target"][1],
    logOR_target_vs_control = s$estimate[1],
    OR_target_vs_control    = exp(s$estimate[1]),
    OR_LCL = exp(s$asymp.LCL[1]),
    OR_UCL = exp(s$asymp.UCL[1]),
    p_two_sided = s$p.value[1],
    p_one_sided = p_one,
    stringsAsFactors = FALSE
  ) -> res

  list(res = res, rates = rates, fit = fit)
}

all_res   <- list()
all_rates <- list()
all_fits  <- list()

for (nm in names(targets)) {
  fr <- fit_one_target(df0, control_ids, targets[[nm]], nm)
  all_res[[nm]]   <- fr$res
  all_rates[[nm]] <- fr$rates
  all_fits[[nm]]  <- fr$fit
}

res_df   <- bind_rows(all_res)
rates_df <- bind_rows(all_rates)

# Multiple-testing correction across the two planned contrasts
res_df$p_one_BH   <- p.adjust(res_df$p_one_sided, method = "BH")
res_df$p_one_Holm <- p.adjust(res_df$p_one_sided, method = "holm")
res_df$p_two_BH   <- p.adjust(res_df$p_two_sided, method = "BH")
res_df$p_two_Holm <- p.adjust(res_df$p_two_sided, method = "holm")

print(res_df)
write.csv(res_df, out_csv_glmm, row.names = FALSE, quote = FALSE)
write.csv(rates_df, out_csv_rates, row.names = FALSE, quote = FALSE)

cat("\nSaved:\n")
cat("  ", out_csv_glmm, "\n", sep = "")
cat("  ", out_csv_rates, "\n", sep = "")
