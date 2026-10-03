# convert list of objects to dataframe of metadata
listToMeta <- function(obj) {
  for (i in 1:length(obj)) {
    mdata <- obj[[i]][[]]
    if (i==1) {
      final <- mdata
    } else {
      final <- rbind(final,mdata)
    }
  }
  return(final)
}


# Plot UMI counts per each modality and sample
plotCounts_original <- function(obj,quantiles,feature,ylabel=feature) {
  dataf <- listToMeta(obj)
  # plot
  pp=ggplot(dataf, aes(x="x",y=dataf[,feature],fill=sample)) +
    theme_bw() +
    geom_half_violin(draw_quantiles = .5) +
    geom_half_point(shape=1,aes(color=sample)) +
    facet_grid(sample~modality,scales = "free_x") +
    theme(strip.text.x = element_text(size = 11, colour = "black", angle = 0, face= 'bold')) +
    theme(strip.text.y = element_text(size = 11, colour = "black", face= 'bold')) +
    # x axis
    xlab("Samples") +
    theme(axis.text.x=element_text(size = 0,angle = 0, hjust = .5)) +
    theme(axis.title.x = element_text(size=0)) +
    # y axis
    ylab(ylabel) +
    theme(axis.text.y=element_text(size = 12)) +
    theme(axis.title.y = element_text(size=14)) +
    theme(legend.position = "none") +
    # add quantiles
    stat_summary(fun = "quantile", fun.args = list(probs = quantiles), 
                 geom = "hline", aes(yintercept = ..y..), linetype = "dashed")
    return(pp)
}

#Modified on 27th Oct 2025
plotCounts <- function(obj, quantiles, feature, ylabel = feature) {
  dataf <- listToMeta(obj)
  dataf$x <- "x"
  
  ggplot(dataf, aes(x = x, y = .data[[feature]], fill = sample, color = sample)) +
    theme_bw() +
    gghalves::geom_half_violin(draw_quantiles = .5, side = "l") +
    ggdist::stat_dots(side = "right", scale = 0.6, alpha = 0.7) +
    facet_grid(sample ~ modality, scales = "free_x") +
    theme(
      strip.text.x = element_text(size = 11, colour = "black", face = "bold"),
      strip.text.y = element_text(size = 11, colour = "black", face = "bold"),
      axis.text.x  = element_text(size = 0),
      axis.title.x = element_text(size = 0),
      axis.text.y  = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      legend.position = "none"
    ) +
    xlab("Samples") + ylab(ylabel) +
    stat_summary(
      data = dataf,
      aes(x = x, y = .data[[feature]], yintercept = after_stat(y)),
      fun = function(z) quantile(z, probs = quantiles),
      geom = "hline",
      linetype = "dashed",
      inherit.aes = FALSE
    )
}


# modified version of the Signac DepthCor - it just adds to the plot the modality name
DepthCorMulMod <- function(obj) {
  nm <- unique(obj[[]][,"modality"])
  p1 <- DepthCor(obj) +
    ggtitle(nm) +
    theme(plot.title = element_text(size=14,hjust=0.5,face='bold'))
  return(p1)
}


# Umap with connected modalities 
plotConnectModal <- function(seurat,group) {
  nmodalities <- length(seurat)
  
  # First, get coords of UMAP, modality name, cluster name and cell barcode for each modality
  umap_embeddings <- list()
  i <- 0
  for(mod in names(seurat)) {
    i <- i + 1
    # get modality name without barcode
    mod2 <- strsplit(mod,"_")[[1]][1]
    # get UMAP 1 and 2 coords
    umap_embeddings[[mod2]]               <- as.data.frame(seurat[[mod]]@reductions[['umap']]@cell.embeddings)
    # adjust UMAP1 for plotting pursoses
    umap_embeddings[[mod2]]$umap_1        <- umap_embeddings[[mod2]]$umap_1 + (i-1)*40 
    # add modality name
    umap_embeddings[[mod2]]$modality      <- unique(seurat[[mod]]$modality) 
    # add cluster name
    umap_embeddings[[mod2]]$cluster       <- seurat[[mod]]@meta.data[,group] 
    # add cell barcode
    umap_embeddings[[mod2]]$cell_barcode  <- rownames(umap_embeddings[[mod2]]) 
  }
  
  # convert the list to dataf
  umap.embeddings.merge <- purrr::reduce(umap_embeddings,rbind)
  # get the names of the common cells among the analysed modalities
  common.cells                        <- table(umap.embeddings.merge$cell_barcode)
  common.cells                        <- names(common.cells[common.cells==nmodalities])
  # subset the umap_embeddings, selecting only the info from the common set of cells
  umap.embeddings.merge               <- umap.embeddings.merge[umap.embeddings.merge$cell_barcode %in% common.cells,]
  
  # label modality - get coords
  coords <- aggregate(umap.embeddings.merge$umap_1,by=list(umap.embeddings.merge$modality),max)
  names(coords) <- c("modality","max")
  coords$min <- aggregate(umap.embeddings.merge$umap_1,by=list(umap.embeddings.merge$modality),min)[,"x"]
  coords$mid_point <- (coords$max + coords$min) / 2
  
  plot <- ggplot(data=umap.embeddings.merge,aes(x=umap_1,y=umap_2,col=cluster)) + 
    geom_point(size=0.2) + 
    geom_line(data=umap.embeddings.merge, aes(group=cell_barcode,col=cluster),alpha=0.2,size=0.02) + 
    theme_classic() + NoAxes() +
    guides(color = guide_legend(override.aes = list(size=5),title="")) +
    geom_text(data=coords,aes(label=modality,x=mid_point,y=max(umap.embeddings.merge$umap_2+.1*umap.embeddings.merge$umap_2)),
              colour='black', fontface = "bold",size=4)
  return(plot)
  
}



# Umap with connected modalities adjusted
plotConnectModal_2 <- function(seurat,group) {
  nmodalities <- length(seurat)
  
  # First, get coords of UMAP, modality name, cluster name and cell barcode for each modality
  umap_embeddings <- list()
  i <- 0
  for(mod in names(seurat)) {
    i <- i + 1
    # get modality name without barcode
    mod2 <- strsplit(mod,"_")[[1]][1]
    # get UMAP 1 and 2 coords
    umap_embeddings[[mod2]]               <- as.data.frame(seurat[[mod]]@reductions[['umap']]@cell.embeddings)
    # adjust UMAP1 for plotting pursoses
    umap_embeddings[[mod2]]$umap_1        <- umap_embeddings[[mod2]]$umap_1 + (i-1)*40 
    # add modality name
    umap_embeddings[[mod2]]$modality      <- unique(seurat[[mod]]$modality) 
    # add cluster name
    umap_embeddings[[mod2]]$cluster       <- seurat[[mod]]@meta.data[,group] 
    # add cell barcode
    umap_embeddings[[mod2]]$cell_barcode  <- rownames(umap_embeddings[[mod2]]) 
  }
  
  # convert the list to dataf
  umap.embeddings.merge <- purrr::reduce(umap_embeddings,rbind)
  # get the names of the common cells among the analysed modalities
  common.cells                        <- table(umap.embeddings.merge$cell_barcode)
  common.cells                        <- names(common.cells[common.cells==nmodalities])
  # subset the umap_embeddings, selecting only the info from the common set of cells
  umap.embeddings.merge               <- umap.embeddings.merge[umap.embeddings.merge$cell_barcode %in% common.cells,]
  
  # label modality - get coords
  coords <- aggregate(umap.embeddings.merge$umap_1,by=list(umap.embeddings.merge$modality),max)
  names(coords) <- c("modality","max")
  coords$min <- aggregate(umap.embeddings.merge$umap_1,by=list(umap.embeddings.merge$modality),min)[,"x"]
  coords$mid_point <- (coords$max + coords$min) / 2
  
  plot <- ggplot(data=umap.embeddings.merge,aes(x=umap_1,y=umap_2,col=cluster)) + 
    geom_point(size=0.5) + 
    geom_line(data=umap.embeddings.merge, aes(group=cell_barcode,col=cluster),alpha=0.02,size=0.01) + 
    theme_classic() + NoAxes() +
    guides(color = guide_legend(override.aes = list(size=5),title="")) +
    geom_text(data=coords,aes(label=modality,x=mid_point,y=max(umap.embeddings.merge$umap_2+.1*umap.embeddings.merge$umap_2)),
              colour='black', fontface = "bold",size=4)
  return(plot)
  
}

# color
plotConnectModal_3 <- function(seurat, group, cols = NULL,
                               x_shift = 40,
                               pt.size = 0.5,
                               line.alpha = 0.02,
                               line.size = 0.01,
                               label.size = 4) {
  nmodalities <- length(seurat)
  
  umap_embeddings <- list()
  i <- 0
  
  for (mod in names(seurat)) {
    i <- i + 1
    mod2 <- strsplit(mod, "_")[[1]][1]
    
    emb <- as.data.frame(seurat[[mod]]@reductions[["umap"]]@cell.embeddings)
    
    # If columns are "UMAP_1/UMAP_2", rename them to "umap_1/umap_2" (safety check)
    if (!all(c("umap_1", "umap_2") %in% colnames(emb))) {
      if (all(c("UMAP_1", "UMAP_2") %in% colnames(emb))) {
        colnames(emb)[colnames(emb) == "UMAP_1"] <- "umap_1"
        colnames(emb)[colnames(emb) == "UMAP_2"] <- "umap_2"
      } else {
        stop("UMAP embeddings must have columns umap_1/umap_2 or UMAP_1/UMAP_2")
      }
    }
    
    emb$umap_1 <- emb$umap_1 + (i - 1) * x_shift
    emb$modality <- unique(seurat[[mod]]$modality)
    
    # IMPORTANT: always coerce cluster to character → factor
    # (ensures discrete coloring even if clusters are numeric like 0/1)
    emb$cluster <- as.character(seurat[[mod]]@meta.data[, group])
    
    emb$cell_barcode <- rownames(emb)
    
    umap_embeddings[[mod2]] <- emb
  }
  
  umap.embeddings.merge <- purrr::reduce(umap_embeddings, rbind)
  
  # Keep only cells common to all modalities
  common.cells <- table(umap.embeddings.merge$cell_barcode)
  common.cells <- names(common.cells[common.cells == nmodalities])
  umap.embeddings.merge <- umap.embeddings.merge[
    umap.embeddings.merge$cell_barcode %in% common.cells, ]
  
  # Define cluster levels explicitly (to stabilize color assignment)
  cluster_levels <- sort(unique(umap.embeddings.merge$cluster))
  umap.embeddings.merge$cluster <- factor(
    umap.embeddings.merge$cluster,
    levels = cluster_levels
  )
  
  # Compute label coordinates for each modality
  coords <- aggregate(umap.embeddings.merge$umap_1,
                      by = list(umap.embeddings.merge$modality),
                      max)
  names(coords) <- c("modality", "max")
  coords$min <- aggregate(umap.embeddings.merge$umap_1,
                          by = list(umap.embeddings.merge$modality),
                          min)[, "x"]
  coords$mid_point <- (coords$max + coords$min) / 2
  
  # Y position for modality labels (slightly above the data range)
  yr <- range(umap.embeddings.merge$umap_2, na.rm = TRUE)
  y_lab <- yr[2] + 0.05 * diff(yr)
  coords$y <- y_lab
  
  p <- ggplot(umap.embeddings.merge,
              aes(x = umap_1, y = umap_2, col = cluster)) +
    geom_point(size = pt.size) +
    geom_line(aes(group = cell_barcode, col = cluster),
              alpha = line.alpha, size = line.size) +
    theme_classic() +
    Seurat::NoAxes() +
    guides(color = guide_legend(
      override.aes = list(size = 5),
      title = ""
    )) +
    geom_text(
      data = coords,
      aes(label = modality, x = mid_point, y = y),
      colour = "black",
      fontface = "bold",
      size = label.size
    )
  
  # Manual color specification:
  # If cols is NULL, keep default behavior.
  if (!is.null(cols)) {
    # If cols has no names, assign colors in cluster_levels order
    if (is.null(names(cols))) {
      if (length(cols) != length(cluster_levels)) {
        stop(
          "cols is unnamed, so it must have the same length as the number of cluster levels: ",
          length(cluster_levels)
        )
      }
      names(cols) <- cluster_levels
    } else {
      # If some cluster levels are missing, fill them with default hues (with warning)
      missing_lv <- setdiff(cluster_levels, names(cols))
      if (length(missing_lv) > 0) {
        warning(
          "cols is missing colors for: ",
          paste(missing_lv, collapse = ", "),
          " -> filled with default hue palette."
        )
        extra <- scales::hue_pal()(length(missing_lv))
        names(extra) <- missing_lv
        cols <- c(cols, extra)
      }
      # Extra names in cols are allowed (they will be ignored)
    }
    
    p <- p + scale_color_manual(values = cols, drop = FALSE)
  }
  
  return(p)
}

# color, line order
plotConnectModal_4 <- function(seurat, group, cols = NULL,
                               # NEW: specify drawing order for LINES only
                               line.last  = NULL,   # e.g. "1" to draw on top
                               line.order = NULL,   # e.g. c("0","1") (bottom -> top)
                               x_shift = 40,
                               pt.size = 0.5,
                               line.alpha = 0.02,
                               line.size = 0.01,
                               label.size = 4) {
  
  nmodalities <- length(seurat)
  
  umap_embeddings <- list()
  i <- 0
  
  for (mod in names(seurat)) {
    i <- i + 1
    mod2 <- strsplit(mod, "_")[[1]][1]
    
    emb <- as.data.frame(seurat[[mod]]@reductions[["umap"]]@cell.embeddings)
    
    # If columns are "UMAP_1/UMAP_2", rename to "umap_1/umap_2"
    if (!all(c("umap_1", "umap_2") %in% colnames(emb))) {
      if (all(c("UMAP_1", "UMAP_2") %in% colnames(emb))) {
        colnames(emb)[colnames(emb) == "UMAP_1"] <- "umap_1"
        colnames(emb)[colnames(emb) == "UMAP_2"] <- "umap_2"
      } else {
        stop("UMAP embeddings must have columns umap_1/umap_2 or UMAP_1/UMAP_2")
      }
    }
    
    emb$umap_1 <- emb$umap_1 + (i - 1) * x_shift
    emb$modality <- unique(seurat[[mod]]$modality)
    
    # Use the 'group' column as cluster (always coerce to character)
    emb$cluster <- as.character(seurat[[mod]]@meta.data[, group])
    emb$cell_barcode <- rownames(emb)
    
    umap_embeddings[[mod2]] <- emb
  }
  
  umap.embeddings.merge <- purrr::reduce(umap_embeddings, rbind)
  
  # Keep only cells common to all modalities
  common.cells <- table(umap.embeddings.merge$cell_barcode)
  common.cells <- names(common.cells[common.cells == nmodalities])
  umap.embeddings.merge <- umap.embeddings.merge[
    umap.embeddings.merge$cell_barcode %in% common.cells, ]
  
  # Cluster levels (stabilize color assignment)
  cluster_levels <- sort(unique(umap.embeddings.merge$cluster))
  umap.embeddings.merge$cluster <- factor(umap.embeddings.merge$cluster, levels = cluster_levels)
  
  # ---- Determine line drawing order (do not change point order) ----
  if (!is.null(line.last) && !is.null(line.order)) {
    stop("Use only one of line.last or line.order (not both).")
  }
  
  line_levels <- cluster_levels
  
  if (!is.null(line.order)) {
    line.order <- unique(as.character(line.order))
    bad <- setdiff(line.order, cluster_levels)
    if (length(bad) > 0) stop("line.order has unknown levels: ", paste(bad, collapse = ", "))
    # Levels not specified are drawn first (underneath)
    line_levels <- c(setdiff(cluster_levels, line.order), line.order)
  } else if (!is.null(line.last)) {
    line.last <- unique(as.character(line.last))
    bad <- setdiff(line.last, cluster_levels)
    if (length(bad) > 0) stop("line.last has unknown levels: ", paste(bad, collapse = ", "))
    line_levels <- c(setdiff(cluster_levels, line.last), line.last)
  }
  
  # ---- Modality label positions ----
  coords <- aggregate(umap.embeddings.merge$umap_1,
                      by = list(umap.embeddings.merge$modality),
                      max)
  names(coords) <- c("modality", "max")
  coords$min <- aggregate(umap.embeddings.merge$umap_1,
                          by = list(umap.embeddings.merge$modality),
                          min)[, "x"]
  coords$mid_point <- (coords$max + coords$min) / 2
  
  yr <- range(umap.embeddings.merge$umap_2, na.rm = TRUE)
  coords$y <- yr[2] + 0.05 * diff(yr)
  
  # ---- Plot: points first, then lines (so lines are drawn on top) ----
  p <- ggplot() +
    geom_point(
      data = umap.embeddings.merge,
      aes(x = umap_1, y = umap_2, col = cluster),
      size = pt.size
    )
  
  # Add line layers per cluster in line_levels order (later layers are on top)
  for (lv in line_levels) {
    df_lv <- umap.embeddings.merge[umap.embeddings.merge$cluster == lv, ]
    p <- p + geom_line(
      data = df_lv,
      aes(x = umap_1, y = umap_2, group = cell_barcode, col = cluster),
      alpha = line.alpha, size = line.size,
      show.legend = FALSE  # keep legend from points only
    )
  }
  
  p <- p +
    theme_classic() + Seurat::NoAxes() +
    guides(color = guide_legend(override.aes = list(size = 5), title = "")) +
    geom_text(
      data = coords,
      aes(label = modality, x = mid_point, y = y),
      colour = "black", fontface = "bold", size = label.size
    )
  
  # ---- Manual color specification ----
  if (!is.null(cols)) {
    if (is.null(names(cols))) {
      if (length(cols) != length(cluster_levels)) {
        stop("cols is unnamed, so it must have the same length as the number of cluster levels: ",
             length(cluster_levels))
      }
      names(cols) <- cluster_levels
    } else {
      missing_lv <- setdiff(cluster_levels, names(cols))
      if (length(missing_lv) > 0) {
        warning("cols is missing colors for: ", paste(missing_lv, collapse = ", "),
                "  -> filled with default hue palette.")
        extra <- scales::hue_pal()(length(missing_lv))
        names(extra) <- missing_lv
        cols <- c(cols, extra)
      }
    }
    p <- p + scale_color_manual(values = cols, drop = FALSE)
  }
  
  return(p)
}

# color, line order, randam dot
plotConnectModal_5 <- function(seurat, group, cols = NULL,
                               # Specify drawing order for LINES only
                               line.last  = NULL,   # e.g. "1"
                               line.order = NULL,   # e.g. c("0","1") (bottom -> top)
                               # NEW: randomize drawing order of POINTS
                               randomize.points = FALSE,
                               point.seed = NULL,
                               
                               x_shift = 40,
                               pt.size = 0.5,
                               line.alpha = 0.02,
                               line.size = 0.01,
                               label.size = 4) {
  
  nmodalities <- length(seurat)
  
  umap_embeddings <- list()
  i <- 0
  
  for (mod in names(seurat)) {
    i <- i + 1
    mod2 <- strsplit(mod, "_")[[1]][1]
    
    emb <- as.data.frame(seurat[[mod]]@reductions[["umap"]]@cell.embeddings)
    
    # "UMAP_1/UMAP_2" -> "umap_1/umap_2"
    if (!all(c("umap_1", "umap_2") %in% colnames(emb))) {
      if (all(c("UMAP_1", "UMAP_2") %in% colnames(emb))) {
        colnames(emb)[colnames(emb) == "UMAP_1"] <- "umap_1"
        colnames(emb)[colnames(emb) == "UMAP_2"] <- "umap_2"
      } else {
        stop("UMAP embeddings must have columns umap_1/umap_2 or UMAP_1/UMAP_2")
      }
    }
    
    emb$umap_1 <- emb$umap_1 + (i - 1) * x_shift
    emb$modality <- unique(seurat[[mod]]$modality)
    
    emb$cluster <- as.character(seurat[[mod]]@meta.data[, group])
    emb$cell_barcode <- rownames(emb)
    
    umap_embeddings[[mod2]] <- emb
  }
  
  umap.embeddings.merge <- purrr::reduce(umap_embeddings, rbind)
  
  # Keep only cells common to all modalities
  common.cells <- table(umap.embeddings.merge$cell_barcode)
  common.cells <- names(common.cells[common.cells == nmodalities])
  umap.embeddings.merge <- umap.embeddings.merge[
    umap.embeddings.merge$cell_barcode %in% common.cells, ]
  
  # Levels for stable color mapping
  cluster_levels <- sort(unique(umap.embeddings.merge$cluster))
  umap.embeddings.merge$cluster <- factor(umap.embeddings.merge$cluster, levels = cluster_levels)
  
  # ---- Line drawing order (do not touch point order here) ----
  if (!is.null(line.last) && !is.null(line.order)) {
    stop("Use only one of line.last or line.order (not both).")
  }
  
  line_levels <- cluster_levels
  if (!is.null(line.order)) {
    line.order <- unique(as.character(line.order))
    bad <- setdiff(line.order, cluster_levels)
    if (length(bad) > 0) stop("line.order has unknown levels: ", paste(bad, collapse = ", "))
    line_levels <- c(setdiff(cluster_levels, line.order), line.order)
  } else if (!is.null(line.last)) {
    line.last <- unique(as.character(line.last))
    bad <- setdiff(line.last, cluster_levels)
    if (length(bad) > 0) stop("line.last has unknown levels: ", paste(bad, collapse = ", "))
    line_levels <- c(setdiff(cluster_levels, line.last), line.last)
  }
  
  # ---- Modality label coords ----
  coords <- aggregate(umap.embeddings.merge$umap_1,
                      by = list(umap.embeddings.merge$modality),
                      max)
  names(coords) <- c("modality", "max")
  coords$min <- aggregate(umap.embeddings.merge$umap_1,
                          by = list(umap.embeddings.merge$modality),
                          min)[, "x"]
  coords$mid_point <- (coords$max + coords$min) / 2
  yr <- range(umap.embeddings.merge$umap_2, na.rm = TRUE)
  coords$y <- yr[2] + 0.05 * diff(yr)
  
  # ---- Randomize POINT drawing order only ----
  df_points <- umap.embeddings.merge
  if (randomize.points) {
    if (!is.null(point.seed)) set.seed(point.seed)
    df_points <- df_points[sample.int(nrow(df_points)), ]
  }
  
  # ---- Plot ----
  p <- ggplot() +
    # Points: randomized order (if enabled)
    geom_point(
      data = df_points,
      aes(x = umap_1, y = umap_2, col = cluster),
      size = pt.size
    )
  
  # Lines: add one layer per cluster to control order (last layer is on top)
  for (lv in line_levels) {
    df_lv <- umap.embeddings.merge[umap.embeddings.merge$cluster == lv, ]
    # Sort by cell barcode and x for more stable line rendering
    df_lv <- df_lv[order(df_lv$cell_barcode, df_lv$umap_1), ]
    
    p <- p + geom_line(
      data = df_lv,
      aes(x = umap_1, y = umap_2, group = cell_barcode, col = cluster),
      alpha = line.alpha, size = line.size,
      show.legend = FALSE
    )
  }
  
  p <- p +
    theme_classic() + Seurat::NoAxes() +
    guides(color = guide_legend(override.aes = list(size = 5), title = "")) +
    geom_text(
      data = coords,
      aes(label = modality, x = mid_point, y = y),
      colour = "black", fontface = "bold", size = label.size
    )
  
  # ---- Manual colors ----
  if (!is.null(cols)) {
    if (is.null(names(cols))) {
      if (length(cols) != length(cluster_levels)) {
        stop("cols is unnamed, so it must have the same length as the number of cluster levels: ",
             length(cluster_levels))
      }
      names(cols) <- cluster_levels
    } else {
      missing_lv <- setdiff(cluster_levels, names(cols))
      if (length(missing_lv) > 0) {
        warning("cols is missing colors for: ", paste(missing_lv, collapse = ", "),
                " -> filled with default hue palette.")
        extra <- scales::hue_pal()(length(missing_lv))
        names(extra) <- missing_lv
        cols <- c(cols, extra)
      }
    }
    p <- p + scale_color_manual(values = cols, drop = FALSE)
  }
  
  return(p)
}

# Function to plot upset on peaks of each modality
getUpsetPeaks <- function(modalities, samples, combined_peaks_ls, input_ls) {
  
  list_up <- list()
  for (mod in modalities) {
    i <- 0
    for (smp in samples) {
      i <- i + 1
      
      # overlap
      combined_mod <- as.data.frame(combined_peaks_ls[[mod]])
      overlap <- GenomicRanges::findOverlaps( toGRanges(combined_mod), input_ls[[paste0(mod,"_",smp)]], select = "first")
      # convert NA to 0s and make it binary
      overlap[is.na(overlap)] = 0
      overlap <- ifelse(overlap>0,1,0)
      combined_mod[,smp] <- overlap
      
      # add new values
      if (i==1) {
        final <- combined_mod
      } else {
        final[,smp] <- overlap
      }
      
      # if last sample, append to list
      if (smp==samples[length(samples)]) {
        list_up[[mod]] <- final
      }
      
      
    }
  }
  
  if (length(list_up)==1) {
    pfinal=upset(list_up[[1]][,6:ncol(list_up[[1]])],colnames(list_up[[1]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
                 set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
                 base_annotations=list('Size'=(intersection_size(counts=FALSE))))  + xlab(names(list_up)[1]) + theme(axis.title.x = element_text(size=12))
  } else if (length(list_up)==2) {
    p1=upset(list_up[[1]][,6:ncol(list_up[[1]])],colnames(list_up[[1]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE)))) + xlab(names(list_up)[1])+ theme(axis.title.x = element_text(size=12)) 
    p2=upset(list_up[[2]][,6:ncol(list_up[[1]])],colnames(list_up[[2]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE)))) + xlab(names(list_up)[2]) + theme(axis.title.x = element_text(size=12))
    pfinal=ggarrange(p1,p2,ncol=2)
  } else if (length(list_up)==3) {
    p1=upset(list_up[[1]][,6:ncol(list_up[[1]])],colnames(list_up[[1]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE))))  + xlab(names(list_up)[1]) + theme(axis.title.x = element_text(size=12))
    p2=upset(list_up[[2]][,6:ncol(list_up[[1]])],colnames(list_up[[2]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE)))) + xlab(names(list_up)[2]) + theme(axis.title.x = element_text(size=12))
    p3=upset(list_up[[3]][,6:ncol(list_up[[1]])],colnames(list_up[[3]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE)))) + xlab(names(list_up)[3]) + theme(axis.title.x = element_text(size=12))
    pfinal=ggarrange(p1,p2,p3,ncol=3)
  } else {
    warning("Upset plot function implemented in this vignette is not implemented to work on more than 3 modalities")
  }
  return(pfinal)
}


plotPassed <- function(mdata.list,xaxis_text=NULL,angle_x=NULL) {
  if (is.null(xaxis_text)) { xaxis_text=12 }
  if (is.null(angle_x)) { angle_x=30 }
  i <- 0
  for (exp in names(mdata.list)) {
    i <- i + 1
    df <- mdata.list[[exp]]
    counts <- as.data.frame(table(df$passedMB))
    counts$Perc <- counts$Freq / sum(counts$Freq) * 100
    counts$sample <- exp
    # bind
    if (i==1) {
      metadata.df <- counts
    } else {
      metadata.df <- rbind(metadata.df,counts)
    }
  }
  names(metadata.df)[1] <- "passedMB"
  # now plot
  pl=ggplot(metadata.df,aes(x=sample,y=Freq,fill=passedMB)) +
    theme_bw() +
    geom_bar(stat="identity",color='black',alpha=.7) +
    # x-axis
    theme(axis.text.x=element_text(angle=angle_x, hjust=1, size = xaxis_text)) + 
    theme(axis.title.x = element_text(size= 0)) +
    xlab("") +
    theme(axis.title.x = element_text(size= 14)) +
    # y-axis
    ylab("Number of cells") +
    theme(axis.text.y=element_text(size = 12)) +
    theme(axis.title.y = element_text(size= 14)) 
  return(pl)
}


plotPassedCells <- function(obj,sample_list,mod_list) {
  # convert obj to df
  i <- 0
  for (nm in names(obj)) {
    i <- i + 1
    # get samples and barcode name
    for (s in sample_list) {
      if (grepl(s,nm)) { 
        sample <- s
      }
    }
    for (m in mod_list) {
      if (grepl(m,nm)) { 
        mod <- m
      }
    }
    obj[[nm]]$sample <- sample
    obj[[nm]]$modality <- mod
    if (i==1) {
      df <- obj[[nm]]
    } else {
      df <- rbind(df,obj[[nm]])
    }
  }
  plot <- ggplot(df,aes(x=all_unique_MB,y=peak_ratio_MB,fill=passedMB)) + 
    theme_bw() +
    geom_point(shape=21,size=.4) +
    scale_x_log10(labels=trans_format('log10',math_format(10^.x))) +
    #coord_cartesian(ylim = c(0,1),xlim = c(10,1000000)) +
    facet_grid(modality~sample) +
    theme(strip.text.x = element_text(size = 11, colour = "black", angle = 0, face= 'bold')) +
    theme(strip.text.y = element_text(size = 11, colour = "black", face= 'bold')) + 
    scale_fill_manual(values = c("#F8766D","#00BA38")) +
    # x-axis
    ylab("Fractio reads in peaks") +
    theme(axis.text.x=element_text(angle=0, hjust=.5, size = 12)) + 
    theme(axis.title.x = element_text(size= 14)) +
    # y-axis
    xlab("UMI") +
    theme(axis.text.y=element_text(size = 12)) +
    theme(axis.title.y = element_text(size= 14))+
    guides(fill = guide_legend(override.aes = list(size=5)))
  plot <- ggarrange(plot,legend='bottom')
  return(plot)
}




commonCellHistonMarks <- function(mod1,name_mod1,mod2,name_mod2,mod3=NULL,name_mod3=NULL,sample) {
  
  # for 2 modalities
  x <- list(name_mod1=mod1$barcode,
            name_mod2=mod2$barcode)
  names(x) <- c(name_mod1,name_mod2)
  pp=ggVennDiagram(x) + 
    scale_fill_gradient(low = "white", high = "white") +
    theme(legend.position = "none") +
    ggtitle(sample) +
    theme(plot.title = element_text(size=14,hjust=0.5,face='bold'))
  
  # 3 modalities
  if (!is.null(mod3)) {
    x <- list(name_mod1=mod1$barcode,
              name_mod2=mod2$barcode,
              name_mod3=mod3$barcode)
    names(x) <- c(name_mod1,name_mod2,name_mod3)
    pp=ggVennDiagram(x) + 
      scale_fill_gradient(low = "white", high = "white") +
      theme(legend.position = "none") +
      ggtitle(sample) +
      theme(plot.title = element_text(size=14,hjust=0.5,face='bold'))
  }
  return(pp)
  
}


