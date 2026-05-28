#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

root <- if (length(args) >= 1) args[[1]] else "Day7R_MERGED"
out  <- if (length(args) >= 2) args[[2]] else "Day7R_MERGED_demo"
n_total <- if (length(args) >= 3) as.integer(args[[3]]) else 500L
seed <- if (length(args) >= 4) as.integer(args[[4]]) else 1L

samples <- c("rep1", "rep2")
mods <- c("Day7R.Dam_bb", "Day7R.K27ac_aa")

set.seed(seed)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

n_per_rep <- ceiling(n_total / length(samples))

read_meta_for_selection <- function(smpl, mod) {
  f <- file.path(root, smpl, mod, "cell_picking", "metadata.csv")
  if (!file.exists(f)) stop("Missing metadata: ", f)

  m <- fread(f)
  if (!"barcode" %in% names(m)) {
    stop("metadata.csv must contain a barcode column: ", f)
  }

  m[, barcode := as.character(barcode)]

  # Prefer Cell Ranger-called cells when available.
  if ("is__cell_barcode" %in% names(m)) {
    suppressWarnings(m[, is__cell_barcode_num := as.numeric(is__cell_barcode)])
    m <- m[is__cell_barcode_num == 1]
    m[, is__cell_barcode_num := NULL]
  }

  # Use existing nanoscope/Cell Ranger QC metrics to avoid extreme low/high cells.
  # These columns are also used later by your current QC code.
  if (all(c("logUMI", "peak_ratio_MB") %in% names(m))) {
    q_logumi_low  <- quantile(m$logUMI, 0.01, na.rm = TRUE)
    q_logumi_high <- quantile(m$logUMI, 0.99, na.rm = TRUE)
    q_frip_low    <- quantile(m$peak_ratio_MB, 0.01, na.rm = TRUE)

    m <- m[
      logUMI > q_logumi_low &
        logUMI < q_logumi_high &
        peak_ratio_MB > q_frip_low
    ]
  }

  m
}

selected_list <- list()

for (smpl in samples) {
  message("Selecting cells for ", smpl)

  metas <- lapply(mods, function(mod) read_meta_for_selection(smpl, mod))
  names(metas) <- mods

  # Important: use only cells present in both Dam and H3K27ac modalities.
  common_barcodes <- Reduce(intersect, lapply(metas, function(x) x$barcode))

  if (length(common_barcodes) == 0) {
    stop("No common barcodes found for ", smpl)
  }

  score_dt <- data.table(barcode = common_barcodes)

  for (mod in mods) {
    m <- metas[[mod]][barcode %in% common_barcodes]
    safe_mod_name <- gsub("[^A-Za-z0-9]+", "_", mod)

    if ("logUMI" %in% names(m)) {
      tmp <- m[, .(barcode, tmp_score = as.numeric(logUMI))]
      setnames(tmp, "tmp_score", paste0("score_", safe_mod_name))
      score_dt <- merge(score_dt, tmp, by = "barcode", all.x = TRUE)
    }
  }

  score_cols <- grep("^score_", names(score_dt), value = TRUE)

  if (length(score_cols) > 0) {
    # Joint quality score: a cell is good only if both modalities are good.
    score_dt[, joint_score := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = score_cols]
    score_dt <- score_dt[order(-joint_score)]

    # Sample from a high-quality pool, rather than always taking the top cells.
    pool_n <- min(nrow(score_dt), max(n_per_rep * 3L, n_per_rep))
    pool <- score_dt[seq_len(pool_n), barcode]
  } else {
    pool <- common_barcodes
  }

  chosen <- if (length(pool) > n_per_rep) {
    sample(pool, n_per_rep)
  } else {
    pool
  }

  selected_list[[smpl]] <- data.table(sample = smpl, barcode = chosen)
}

selected <- rbindlist(selected_list)
fwrite(selected, file.path(out, "selected_cells.tsv"), sep = "\t")

# Write subset metadata to the demo tree.
# Keep original metadata values. Do not overwrite peak_ratio_MB with 100,
# because your current QC code uses strict ">" quantile filtering.
for (smpl in samples) {
  for (mod in mods) {
    src_meta <- file.path(root, smpl, mod, "cell_picking", "metadata.csv")
    dst_dir <- file.path(out, smpl, mod, "cell_picking")
    dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)

    meta <- fread(src_meta)
    meta[, barcode := as.character(barcode)]

    keep <- selected[sample == smpl, barcode]
    meta_sub <- meta[barcode %in% keep]

    if (nrow(meta_sub) == 0) {
      stop("No metadata rows retained for ", smpl, " / ", mod)
    }

    fwrite(meta_sub, file.path(dst_dir, "metadata.csv"))
  }
}

message("Wrote selected cells to: ", file.path(out, "selected_cells.tsv"))
print(selected[, .N, by = sample])