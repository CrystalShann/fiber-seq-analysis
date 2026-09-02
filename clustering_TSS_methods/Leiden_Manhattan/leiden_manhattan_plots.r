# leiden_manhattan_plots.r
#
# Plots for the Leiden + Manhattan single-molecule clustering
# (see leiden_manhattan_functions.r). Every function takes the result object
# returned by leiden_manhattan_cluster().
#
#   plot_cluster_met_profiles()  per-cluster m6A methylation proportion at each
#                                m6A site, one panel per cluster - the topic
#                                model's cluster_met_profiles plot, so the
#                                clusterings can be compared panel for panel
#   plot_cluster_composition()   cluster proportions per timepoint and
#                                timepoint proportions per cluster
#   plot_cluster_structure()     one column per read, grouped by timepoint and
#                                coloured by cluster (the hard-assignment
#                                analogue of the topic model's structure plot)
#   plot_cluster_heatmap()       read x feature methylation heatmap, rows split
#                                by cluster, with cluster and timepoint
#                                annotations
#   plot_met_fraction_lines()    m6A fraction per feature, one panel per cluster
#
# The bar panels are drawn on a fixed 0-1 axis so clusters and genes stay
# comparable; the line plot uses the data range, since that is the plot for
# reading profile shape.

suppressMessages({
  requireNamespace("ComplexHeatmap")
})

# 25-colour cluster palette (the topic model's colors_25 / the hamming
# notebook's PHASE_CLUSTER_COLORS), extended by interpolation if Leiden ever
# returns more clusters than that
LEIDEN_CLUSTER_COLORS <- c(
  "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1",
  "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
  "green1", "yellow4", "yellow3", "darkorange4", "brown"
)

# timepoint colours (sequential, as in the spearman notebook)
LEIDEN_TIMEPOINT_COLORS <- c("LPS_0" = "#bdbdbd", "LPS_5" = "#6baed6",
                             "LPS_10" = "#2171b5", "LPS_15" = "#08306b")

cluster_palette <- function(levels) {
  n <- length(levels)
  cols <- if (n <= length(LEIDEN_CLUSTER_COLORS)) {
    LEIDEN_CLUSTER_COLORS[seq_len(n)]
  } else {
    grDevices::colorRampPalette(LEIDEN_CLUSTER_COLORS)(n)
  }
  setNames(cols, levels)
}

# sample_name as an ordered factor, whatever the caller passed in
as_timepoint_factor <- function(x, timepoint_cols = LEIDEN_TIMEPOINT_COLORS) {
  if (is.factor(x)) return(x)
  lv <- intersect(names(timepoint_cols), unique(as.character(x)))
  factor(as.character(x), levels = c(lv, setdiff(unique(as.character(x)), lv)))
}


# ---------------------------------------------------------------------------
# 1. Per-cluster m6A methylation proportion at each m6A site 

# One panel per cluster, x = genomic coordinate, y = methylation proportion.
# ---------------------------------------------------------------------------
plot_cluster_met_profiles <- function(res, met_mat, tss = NULL, main = NULL) {
  # cluster_site_profiles() takes cluster labels and computes mean m6a value at every 
  # position for each cluster
  prof <- cluster_site_profiles(res, met_mat)
  lv   <- levels(res$assignments$cluster)
  pal  <- cluster_palette(lv)

  p_list <- lapply(lv, function(cl) {
    d <- prof[prof$cluster == cl, ]
    gg <- ggplot(d, aes(x = pos, y = met)) +
      geom_col(fill = pal[cl]) +
      ylim(0, 1) +
      ylab("met prop.") + xlab("pos") +
      ggtitle(sprintf("%s (n=%d)", cl, d$n_reads[1])) +
      theme_cowplot(font_size = 10) +
      theme(plot.title = element_text(hjust = 0.5))
    if (!is.null(tss))
      gg <- gg + geom_vline(xintercept = tss, linetype = "dashed", color = "grey40")
    gg
  })
  if (!is.null(main))
    p_list <- c(list(cowplot::ggdraw() +
                       cowplot::draw_label(main, fontface = "bold", size = 12)),
                p_list)
  cowplot::plot_grid(plotlist = p_list, ncol = 1,
                     rel_heights = if (is.null(main)) 1
                                   else c(0.4, rep(1, length(p_list) - 1)))
}


# ---------------------------------------------------------------------------
# 2a. Cluster composition by timepoint: cluster proportions within each
# timepoint (stacked), and each timepoint's reads spread over the clusters.
# ---------------------------------------------------------------------------
plot_cluster_composition <- function(res, main = NULL,
                                     timepoint_cols = LEIDEN_TIMEPOINT_COLORS) {
  df <- res$assignments
  df$sample_name <- as_timepoint_factor(df$sample_name, timepoint_cols)
  pal <- cluster_palette(levels(df$cluster))

  p1 <- ggplot(df, aes(x = sample_name, fill = cluster)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = pal) +
    labs(x = NULL, y = "fraction of reads",
         title = if (is.null(main)) NULL else paste0("Cluster composition per timepoint, ", main)) +
    theme_cowplot(font_size = 10)

  prop <- df %>%
    dplyr::count(cluster, sample_name) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(proportion = n / sum(n)) %>%
    dplyr::ungroup()
  p2 <- ggplot(prop, aes(x = cluster, y = proportion, fill = sample_name)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = timepoint_cols, name = "timepoint") +
    labs(x = "cluster", y = "proportion of that timepoint's reads",
         title = if (is.null(main)) NULL else paste0("Timepoints per cluster, ", main)) +
    theme_cowplot(font_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  cowplot::plot_grid(p1, p2, ncol = 1)
}


# ---------------------------------------------------------------------------
# 2b. Cluster structure: one column per read, grouped by timepoint, coloured by
# cluster (the hard-assignment analogue of the topic model's structure plot).
# ---------------------------------------------------------------------------
plot_cluster_structure <- function(res, main = NULL,
                                   timepoint_cols = LEIDEN_TIMEPOINT_COLORS) {
  df <- res$assignments
  df$sample_name <- as_timepoint_factor(df$sample_name, timepoint_cols)
  df <- df[order(df$sample_name, df$cluster), ]
  df$x <- stats::ave(seq_len(nrow(df)), df$sample_name, FUN = seq_along)

  ggplot(df, aes(x = x, y = 1, fill = cluster)) +
    geom_col(width = 1) +
    facet_grid(~ sample_name, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = cluster_palette(levels(df$cluster))) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = "reads (grouped by timepoint)", y = NULL, title = main) +
    theme_cowplot(font_size = 10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
}


# ---------------------------------------------------------------------------
# 3. Read x feature methylation heatmap: 

# rows = reads split by cluster (each and ordered by read start within cluster,
# columns = features in genomic order, left annotation = cluster + timepoint.
# No colour gradient: a feature with any m6A call is black, none is white,

# ---------------------------------------------------------------------------
plot_cluster_heatmap <- function(res, main = NULL,
                                 timepoint_cols = LEIDEN_TIMEPOINT_COLORS) {
  # takes per read assignment from the clustering result
  df <- res$assignments
  df$sample_name <- as_timepoint_factor(df$sample_name, timepoint_cols)
  # sort by cluster
  # within each cluster, sort by read start position
  o  <- order(df$cluster, df$start)
  df <- df[o, ]
  # Takes the feature matrix used for  clustering 
  # and puts its rows in exactly the same order as df
  mat <- res$feat_mat[df$RID, , drop = FALSE]
  mat <- matrix(ifelse(is.na(mat), NA, ifelse(mat > 0, "m6A", "no m6A")),
                nrow(mat), ncol(mat), dimnames = dimnames(mat))

  # create row annotations
  ha <- ComplexHeatmap::rowAnnotation(
    cluster   = df$cluster,
    timepoint = df$sample_name,
    col = list(cluster   = cluster_palette(levels(df$cluster)),
               timepoint = timepoint_cols[levels(df$sample_name)]))

  # rows are reads, columns are m6a sites
  ComplexHeatmap::Heatmap(
    mat,
    name = "m6A call",
    col  = c("m6A" = "black", "no m6A" = "white"),
    na_col = "grey85",
    show_row_names = FALSE, show_column_names = FALSE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = df$cluster, row_gap = unit(0.6, "mm"),
    row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 8),
    width = unit(11, "cm"), height = unit(14, "cm"),
    use_raster = TRUE,
    column_title = main,
    column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
    left_annotation = ha)
}


# ---------------------------------------------------------------------------
# 4. m6A fraction per feature, one panel per cluster 
# ---------------------------------------------------------------------------
plot_met_fraction_lines <- function(res, tss = NULL, main = NULL, smooth_k = 1) {
  P   <- res$profiles
  pos <- as.numeric(colnames(P))

  smooth_row <- function(v) {
    if (smooth_k <= 1) return(v)
    as.numeric(stats::filter(v, rep(1 / smooth_k, smooth_k), sides = 2))
  }

  df <- do.call(rbind, lapply(rownames(P), function(cl)
    data.frame(cluster = cl, pos = pos, value = smooth_row(P[cl, ]))))
  df$cluster <- factor(df$cluster, levels = rownames(P))
  n <- table(res$assignments$cluster)
  levels(df$cluster) <- sprintf("%s (n=%d)", rownames(P), as.integer(n[rownames(P)]))

  ylab <- if (res$params$window_size == 0) "m6A fraction per site"
          else sprintf("mean m6A per %d-bp window", res$params$window_size)
  sub  <- if (smooth_k > 1) sprintf("rolling mean over %d features", smooth_k) else NULL

  gg <- ggplot(df[!is.na(df$value), ], aes(x = pos, y = value, color = cluster)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = setNames(cluster_palette(rownames(P)), levels(df$cluster)),
                       guide = "none") +
    facet_wrap(~ cluster, ncol = 1, strip.position = "right") +
    labs(x = "genomic position", y = ylab, title = main, subtitle = sub) +
    theme_cowplot(font_size = 10)
  if (!is.null(tss))
    gg <- gg + geom_vline(xintercept = tss, linetype = "dashed", color = "grey40")
  gg
}
