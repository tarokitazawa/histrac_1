# In vivo HisTrac Dam ambiguity analysis across all major cell types
# Input object: nanoscope_final_peaks_nanoscopefilter_in_vivo_Cell_Type.rds
#
# Goal:
# - Identify the most ambiguous cell-type separations in the HisTrac Dam modality
#   (W8R stage only) using two complementary approaches:
#     1) normalized centroid separation (D_norm)
#     2) local kNN-based symmetric mixing (mix_sym)
#
#
# Note on joint.ref.lsi (see in_vivo_cell_type_annotation.R):
# - Dam medoids are computed in joint.ref.lsi, a shared reference-space LSI in which
#   W8/W8R reference cells retain their original Dam LSI coordinates and P8 cells are
#   projected onto that adult reference.
# - Therefore, joint.ref.lsi uses a re-indexed reference-space representation rather
#   than the original raw LSI numbering.
# - For this reason, Dam analyses based on joint.ref.lsi appropriately start from
#   dimension 1, even though the original Dam LSI workflow typically used dimensions
#   beginning from 2 after excluding the first raw LSI component.

library(Seurat)
library(ggplot2)
library(FNN)
library(viridis)
library(dplyr)

# -------------------------------------------------------------------
# Input
# -------------------------------------------------------------------
combined.obj.ls3 <- readRDS("nanoscope_final_peaks_nanoscopefilter_in_vivo_Cell_Type.rds")

# -------------------------------------------------------------------
# Select HisTrac Dam modality (W8R only) and major cell types
# -------------------------------------------------------------------
obj <- subset(combined.obj.ls3$Dam_bb, subset = stage == "W8R")

cts <- c(
  "L2-3", "L4-5", "L5-Pou3f1", "L5-Tshz2", "L6-Foxp2",
  "Sst", "Vip-Reln", "Pvalb", "Astrocyte"
)

obj <- subset(obj, subset = Cell_Type %in% cts)
obj$Cell_Type <- factor(as.character(obj$Cell_Type), levels = cts)

# Dam ambiguity is evaluated in the reference-space Dam LSI
dims_use <- 1:13
X <- Embeddings(obj, "joint.ref.lsi")[, dims_use, drop = FALSE]
lab <- as.character(obj$Cell_Type)

stopifnot(nrow(X) == length(lab))

# -------------------------------------------------------------------
# Method 1: normalized centroid separation
# - Smaller D_norm means more ambiguous separation
# -------------------------------------------------------------------
centroids <- sapply(cts, function(ct) colMeans(X[lab == ct, , drop = FALSE]))

within_rms <- sapply(cts, function(ct) {
  Xi <- X[lab == ct, , drop = FALSE]
  mu <- colMeans(Xi)
  sqrt(mean(rowSums((Xi - matrix(mu, nrow(Xi), length(mu), byrow = TRUE))^2)))
})

D_between <- as.matrix(dist(t(centroids), method = "euclidean"))
D_norm <- D_between / sqrt(outer(within_rms^2, within_rms^2, "+"))

rank_pairs <- function(M, decreasing = FALSE) {
  stopifnot(nrow(M) == ncol(M))
  types <- rownames(M)
  idx <- which(upper.tri(M), arr.ind = TRUE)
  df <- data.frame(
    type1 = types[idx[, 1]],
    type2 = types[idx[, 2]],
    score = M[idx],
    stringsAsFactors = FALSE
  )
  df[order(df$score, decreasing = decreasing), ]
}

sep_rank <- rank_pairs(D_norm, decreasing = FALSE)
write.csv(sep_rank, "ambiguity_rank_Dnorm.csv", row.names = FALSE, quote = FALSE)

# -------------------------------------------------------------------
# Method 2: local kNN mixing
# - Larger mix_sym means more ambiguous separation
# -------------------------------------------------------------------
k <- 30
nn <- FNN::get.knn(X, k = k)$nn.index

from <- rep(lab, each = k)
to   <- lab[as.vector(t(nn))]
edges <- data.frame(from = from, to = to, stringsAsFactors = FALSE)

mix_counts <- xtabs(~ from + to, data = edges)
mix_counts <- as.matrix(mix_counts)

mix_counts2 <- matrix(0, nrow = length(cts), ncol = length(cts), dimnames = list(cts, cts))
mix_counts2[rownames(mix_counts), colnames(mix_counts)] <- mix_counts

mix_prop <- sweep(mix_counts2, 1, rowSums(mix_counts2), "/")
mix_sym <- (mix_prop + t(mix_prop)) / 2
diag(mix_sym) <- NA

mix_rank <- rank_pairs(mix_sym, decreasing = TRUE)
write.csv(mix_rank, "ambiguity_rank_mixSym.csv", row.names = FALSE, quote = FALSE)

# -------------------------------------------------------------------
# Plot helpers
# -------------------------------------------------------------------
make_heatmap_df <- function(M, cts) {
  M2 <- M[cts, cts, drop = FALSE]
  df <- as.data.frame(as.table(M2), stringsAsFactors = FALSE)
  colnames(df) <- c("from", "to", "value")

  df$from <- factor(df$from, levels = rev(cts))
  df$to   <- factor(df$to, levels = cts)
  df
}

plot_heatmap <- function(df, title, direction = 1, trans = "identity") {
  ggplot(df, aes(x = to, y = from, fill = value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    coord_fixed() +
    scale_fill_viridis_c(
      na.value = "grey90",
      direction = direction,
      trans = trans
    ) +
    labs(
      x = "To (neighbor / other class)",
      y = "From (query class)",
      fill = NULL,
      title = title
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(face = "bold")
    )
}

# -------------------------------------------------------------------
# PDF 1: D_norm heatmap
# - smaller values are more ambiguous, so use direction = -1
# -------------------------------------------------------------------
D_heat <- D_norm
diag(D_heat) <- NA
df_D <- make_heatmap_df(D_heat, cts)

p_D <- plot_heatmap(
  df_D,
  title = "HisTrac Dam ambiguity: normalized centroid separation (D_norm)",
  direction = -1,
  trans = "log10"
)

ggsave("Heatmap_Dnorm.pdf", p_D, width = 6.2, height = 5.6, useDingbats = FALSE)

# -------------------------------------------------------------------
# PDF 2: mix_sym heatmap
# - larger values are more ambiguous
# -------------------------------------------------------------------
M_heat <- mix_sym
diag(M_heat) <- NA
df_M <- make_heatmap_df(M_heat, cts)

p_M <- plot_heatmap(
  df_M,
  title = "HisTrac Dam ambiguity: kNN symmetric mixing (mix_sym)",
  direction = 1,
  trans = "identity"
)

ggsave("Heatmap_mixSym.pdf", p_M, width = 6.2, height = 5.6, useDingbats = FALSE)

# -------------------------------------------------------------------
# Console summary of the most ambiguous pairs
# -------------------------------------------------------------------
cat("\nTop ambiguous pairs by D_norm (smaller = more ambiguous):\n")
print(head(sep_rank, 10))

cat("\nTop ambiguous pairs by mix_sym (larger = more ambiguous):\n")
print(head(mix_rank, 10))
