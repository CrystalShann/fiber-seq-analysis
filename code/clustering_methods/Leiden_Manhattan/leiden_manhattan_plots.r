# leiden_manhattan_plots.r
#
# Shared plots and report export for Leiden + Manhattan single-molecule
# clustering, including the Fourier-feature notebook. Plot-specific inputs
# are explicit; Fourier spectrum and phase plots stay in the FFT notebook.
# LCL-only report orchestration and figure exports live in hidden chunks of
# leiden_LCL.Rmd; helpers used by multiple notebooks stay here.
# plot_cluster_heatmap() displays timepoint-annotated clustering features;
# plot_genomic_cluster_heatmap() displays genomic m6A with sample/allele tracks.
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

# Shared report export: retain PNG previews until Pandoc embeds them in HTML.
# A writer owns its output directory and manifest, so notebooks do not need
# global plotting state. Supply draw for objects such as ComplexHeatmap heatmaps.
embed_report_png <- function(path, label, image_dir) {
  stopifnot(file.exists(path), file.info(path)$size > 0)
  dir.create(image_dir, recursive = TRUE, showWarnings = FALSE)
  image_path <- tempfile("plot-", tmpdir = image_dir, fileext = ".png")
  stopifnot(file.copy(path, image_path))
  label_attribute <- as.character(htmltools::htmlEscape(label, attribute = TRUE))
  label <- as.character(htmltools::htmlEscape(label))
  image_path <- as.character(htmltools::htmlEscape(normalizePath(image_path), attribute = TRUE))
  cat('\n\n<figure><img src="', image_path,
      '" alt="', label_attribute, '" style="max-width:100%;height:auto;" />',
      '<figcaption>', label, '</figcaption></figure>\n\n', sep = "")
  invisible(NULL)
}

create_report_plot_writer <- function(image_dir, dpi = 160, bg = "white") {
  dir.create(image_dir, recursive = TRUE, showWarnings = FALSE)
  image_dir <- normalizePath(image_dir)
  plot_files <- character()
  save_plot <- function(filename, plot, width, height, ..., dpi = dpi_default,
                        bg = bg_default, draw = NULL, label = NULL) {
    preview <- tempfile(fileext = ".png")
    on.exit(unlink(preview), add = TRUE)
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    if (is.null(draw)) {
      ggplot2::ggsave(filename, plot = plot, width = width, height = height,
                      ..., dpi = dpi, bg = bg)
      ggplot2::ggsave(preview, plot = plot, width = width, height = height,
                      ..., dpi = dpi, bg = bg)
    } else {
      stopifnot(is.function(draw), tolower(tools::file_ext(filename)) == "pdf")
      grDevices::pdf(filename, width = width, height = height, bg = bg)
      tryCatch(draw(plot), finally = grDevices::dev.off())
      grDevices::png(preview, width = width, height = height, units = "in",
                     res = dpi, type = "cairo", bg = bg)
      tryCatch(draw(plot), finally = grDevices::dev.off())
    }
    if (is.null(label)) label <- tools::file_path_sans_ext(basename(filename))
    embed_report_png(preview, label, image_dir)
    plot_files <<- c(plot_files, filename)
    invisible(filename)
  }
  dpi_default <- dpi
  bg_default <- bg
  list(save = save_plot, files = function() plot_files)
}

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


# ---- Fiber-seq read, allele and composition plots ----
LCL_CLUSTER_COLORS <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1",
  "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
  "green1", "yellow4", "yellow3", "darkorange4", "brown")
cluster_id_palette <- function(levels) {
  idx <- suppressWarnings(as.integer(sub("^cluster", "", levels)))
  if (anyNA(idx) || any(!grepl("^cluster[0-9]+$", levels))) idx <- seq_along(levels)
  cols <- if (max(idx) <= length(LCL_CLUSTER_COLORS)) LCL_CLUSTER_COLORS[idx] else
    grDevices::colorRampPalette(LCL_CLUSTER_COLORS)(max(idx))[idx]
  setNames(cols, levels)
}

LCL_HAPLOTYPE_COLORS <- c(HP1 = "#ADD8E6", HP2 = "#FFF2AE", unphased = "#999999", pooled = "#BBBBBB")


LCL_TRACK_COLORS <- c(m6A = "#800080", "ft_nuc_130-160bp" = "#4d4d4d",
  "FiberHMM_10-30bp" = "#fdae6b", "FiberHMM_40-60bp" = "#f16913",
  "FiberHMM_60-80bp" = "#a63603")
LCL_TRACK_LABELS <- c(m6A = "m6A", "ft_nuc_130-160bp" = "nucleosome footprint (130-160 bp)",
  "FiberHMM_10-30bp" = "TF footprint 10-30 bp", "FiberHMM_40-60bp" = "TF footprint 40-60 bp",
  "FiberHMM_60-80bp" = "TF footprint 60-80 bp")
feature_colors <- function(tracks) {
  stopifnot(all(tracks %in% names(LCL_TRACK_COLORS)))
  LCL_TRACK_COLORS[tracks]
}

plot_smf_reads <- function(result, records, sample_colors, tracks = LCL_FOOTPRINT_TRACKS,
                                include_haplotype = FALSE) {
  tracks <- LCL_FOOTPRINT_TRACKS
  records <- records[records$track %in% tracks, , drop = FALSE]
  prepared <- prepare_read_tracks(result, records, sample_colors)
  reads <- prepared$reads
  features <- prepared$features[prepared$features$track %in% tracks, , drop = FALSE]
  anchor <- prepared$anchor
  width <- max(1, anchor$right - anchor$left)
  feature_colors <- feature_colors(tracks)
  # Draw nucleosomes first and TF footprints on top, as in TNF Fig 1.
  features <- features[order(match(features$track, tracks)), ]
  features$ymin <- features$row - 0.45
  features$ymax <- features$row + 0.45
  features$fill <- unname(feature_colors[features$track])
  cluster_colors <- cluster_id_palette(levels(reads$cluster))
  make_bar <- function(offset, colors) data.frame(cluster = reads$cluster, row = reads$row,
    left = anchor$left - offset * width, right = anchor$left - (offset - 0.025) * width,
    fill = unname(colors))
  bars <- rbind(make_bar(0.075, cluster_colors[as.character(reads$cluster)]),
                 make_bar(0.040, sample_colors[reads$sample_label]))
  allele_colors <- if (identical(result$region$region_type, "top_asfire_het"))
    fiberseq_category_palette(result$allele_display_levels) else LCL_HAPLOTYPE_COLORS
  allele <- if (identical(result$region$region_type, "top_asfire_het")) as.character(reads$allele_display) else reads$haplotype
  if (include_haplotype) bars <- rbind(bars, make_bar(0.110, allele_colors[allele]))
  legend_colors <- c(feature_colors, sample_colors,
    if (include_haplotype) allele_colors)
  legend_labels <- c(unname(LCL_TRACK_LABELS[tracks]), names(sample_colors),
    if (include_haplotype) names(allele_colors))
  counts <- table(reads$cluster)
  cluster_labels <- setNames(paste0(names(counts), " (n=", as.integer(counts), ")"), names(counts))
  ticks <- pretty(c(anchor$left, anchor$right), n = 6)
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = reads, ggplot2::aes(x = left, xend = right, y = row, yend = row),
      color = "grey80", linewidth = 0.25) +
    ggplot2::geom_rect(data = features,
      ggplot2::aes(xmin = left, xmax = right, ymin = ymin, ymax = ymax, fill = fill)) +
    ggplot2::geom_rect(data = bars,
      ggplot2::aes(xmin = left, xmax = right, ymin = row - 0.5, ymax = row + 0.5, fill = fill)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
    ggplot2::scale_fill_identity(guide = "legend", name = NULL,
      breaks = unname(legend_colors), labels = legend_labels) +
    ggplot2::scale_x_continuous(limits = c(anchor$left - (if (include_haplotype) 0.12 * width else 0.085 * width),
      anchor$right + 0.5), breaks = ticks[ticks >= anchor$left & ticks <= anchor$right], expand = c(0, 0)) +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::facet_grid(cluster ~ ., scales = "free_y", space = "free_y",
      labeller = ggplot2::as_labeller(cluster_labels)) +
    ggplot2::labs(x = anchor$x_label, y = "SMF reads",
      title = paste(c(result$region$annotation, result$group_id), collapse = " | "),
      subtitle = "One row per retained molecule; saved m6A-defined Leiden clusters") +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4, byrow = TRUE)) +
    cowplot::theme_cowplot(font_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(angle = 0),
      panel.spacing.y = grid::unit(2, "mm"), legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8), legend.key.size = grid::unit(3, "mm"))
}

# all features overlaid in one panel for every cluster.
plot_signal_profile <- function(profiles, region, clusters,
                                     title = region$annotation) {
  profiles <- profiles[profiles$track %in% names(LCL_TRACK_COLORS), , drop = FALSE]
  anchor <- plot_anchor(region)
  profiles$relative_pos <- anchor$direction * (profiles$pos - anchor$anchor)
  profiles$cluster <- factor(profiles$cluster, levels = clusters)
  profiles$track <- factor(profiles$track, levels = names(LCL_TRACK_COLORS)[names(LCL_TRACK_COLORS) %in% profiles$track])
  profiles <- profiles[order(profiles$cluster, profiles$track, profiles$relative_pos), ]
  colors <- feature_colors(levels(profiles$track))
  counts <- unique(profiles[, c("cluster", "n_reads")])
  labels <- setNames(paste0(counts$cluster, " (n=", counts$n_reads, ")"), counts$cluster)
  plot <- ggplot2::ggplot(profiles,
    ggplot2::aes(relative_pos, fraction, fill = track, color = track, group = track))
  plot <- plot + ggplot2::geom_ribbon(
    data = profiles[profiles$track == "ft_nuc_130-160bp", , drop = FALSE],
    ggplot2::aes(ymin = 0, ymax = fraction), fill = "grey60", color = NA, alpha = .18) +
    ggplot2::geom_line(linewidth = 0.5)
  plot + ggplot2::facet_wrap(~cluster, ncol = 1, labeller = ggplot2::as_labeller(labels)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = colors, labels = LCL_TRACK_LABELS) +
    ggplot2::scale_color_manual(values = colors, labels = LCL_TRACK_LABELS) +
    ggplot2::scale_x_continuous(limits = c(anchor$left - 0.5, anchor$right + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = anchor$x_label, y = "Fraction of cluster reads", title = title,
      subtitle = "m6A calls and footprint occupancy; retained heterozygous reads", fill = "Plot track", color = "Plot track") +
    cowplot::theme_cowplot(font_size = 10) + cowplot::panel_border() +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4), color = ggplot2::guide_legend(ncol = 4)) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(5.5, 16, 5.5, 5.5))
}

sample_palette <- function(sample_names) {
  labels <- sub("_.*$", "", sample_names)
  labels <- unique(labels[order(as.integer(sub("^AL-?([0-9]+).*$", "\\1", labels)))])
  setNames(grDevices::hcl.colors(length(labels), "Dynamic"), labels)
}

plot_genomic_cluster_heatmap <- function(res, region, sample_colors,
                                     include_haplotype = FALSE, variants = NULL,
                                     show_cluster_profiles = FALSE, split_alleles = FALSE,
                                     cluster_label = "Saved m6A-defined clusters") {
  assignments <- res$assignments
  if (split_alleles) stopifnot(include_haplotype, !is.null(assignments$allele_display))
  ordering <- if (split_alleles) order(assignments$allele_display, assignments$cluster, assignments$start, assignments$RID) else
    order(assignments$cluster, assignments$start, assignments$RID)
  assignments <- assignments[ordering, , drop = FALSE]
  sample <- factor(sub("_.*$", "", assignments$sample_name), levels = names(sample_colors))
  annotation <- ComplexHeatmap::rowAnnotation(
    cluster = assignments$cluster, sample = sample,
    col = list(cluster = cluster_id_palette(levels(assignments$cluster)), sample = sample_colors))
  if (include_haplotype && identical(region$region_type, "top_asfire_het")) {
    annotation <- ComplexHeatmap::rowAnnotation(
      cluster = assignments$cluster, sample = sample, allele = assignments$allele_display,
      col = list(cluster = cluster_id_palette(levels(assignments$cluster)), sample = sample_colors,
        allele = fiberseq_category_palette(res$allele_display_levels)))
  } else if (include_haplotype) {
    annotation <- ComplexHeatmap::rowAnnotation(
      cluster = assignments$cluster, sample = sample, haplotype = assignments$haplotype,
      col = list(cluster = cluster_id_palette(levels(assignments$cluster)), sample = sample_colors,
        haplotype = LCL_HAPLOTYPE_COLORS))
  }
  met <- matrix(0L, nrow(assignments), region$width,
                  dimnames = list(assignments$RID, seq.int(region$analysis_start, region$analysis_end)))
  met[, match(colnames(res$site_met_mat), colnames(met))] <- as.matrix(res$site_met_mat[assignments$RID, , drop = FALSE])
  display <- met
  colors <- c("0" = "white", "1" = "black")
  legend <- list(at = c(0, 1), labels = c("no m6A call", "m6A"))
  legend_name <- "m6A"
  top <- NULL
  if (show_cluster_profiles) {
    cluster_colors <- cluster_id_palette(levels(assignments$cluster))
    observed <- match(colnames(res$site_met_mat), colnames(met))
    profile_annotations <- lapply(levels(assignments$cluster), function(cluster) {
      selected <- assignments$cluster == cluster
      ComplexHeatmap::anno_lines(
        stats::approx(observed, colMeans(met[selected, observed, drop = FALSE]),
          xout = seq_len(ncol(met)), rule = 1)$y,
        ylim = c(0, 1), gp = grid::gpar(col = cluster_colors[[cluster]], lwd = 0.7),
        axis_param = list(at = c(0, 0.5, 1), labels = c("0", ".5", "1")),
        height = grid::unit(12, "mm"))
    })
    names(profile_annotations) <- paste0(levels(assignments$cluster), " m6A")
    top <- do.call(ComplexHeatmap::HeatmapAnnotation, c(profile_annotations,
      list(annotation_name_gp = grid::gpar(fontsize = 8), gap = grid::unit(1.5, "mm"))))
  }
  ticks <- unique(round(seq(1, ncol(display), length.out = 5L)))
  bottom_parts <- list(coordinate = ComplexHeatmap::anno_mark(at = ticks,
    labels = format(as.integer(colnames(display)[ticks]), scientific = FALSE, trim = TRUE),
    which = "column", side = "bottom", labels_gp = grid::gpar(fontsize = 8)))
  if (identical(region$region_type, "promoter") && !is.null(region$tss) &&
      !is.na(region$tss) && region$tss >= region$analysis_start && region$tss <= region$analysis_end) {
    bottom_parts$TSS <- ComplexHeatmap::anno_mark(
      at = region$tss - region$analysis_start + 1L,
      labels = paste0("TSS: ", region$chr, ":", region$tss),
      which = "column", side = "bottom", labels_gp = grid::gpar(fontsize = 8, col = "#D55E00"))
  }
  snps <- if (!is.null(region$focal_snp)) data.frame(pos = region$focal_pos,
    label = paste0(region$focal_snp, " ", region$ref, ">", region$alt)) else window_snps(variants, region)
  if (nrow(snps)) {
    bottom_parts$SNP <- ComplexHeatmap::anno_mark(
      at = snps$pos - region$analysis_start + 1L, labels = snps$label,
      which = "column", side = "bottom", labels_gp = grid::gpar(fontsize = 7),
      link_gp = grid::gpar(col = "#984EA3"))
  }
  bottom <- do.call(ComplexHeatmap::HeatmapAnnotation, c(bottom_parts,
    list(annotation_name_gp = grid::gpar(fontsize = 8))))
  ComplexHeatmap::Heatmap(
    display, name = legend_name, col = colors, heatmap_legend_param = legend,
    cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE,
    row_split = if (split_alleles) data.frame(allele = assignments$allele_display,
      cluster = assignments$cluster) else assignments$cluster,
    row_gap = grid::unit(if (split_alleles) 3 else 0.6, "mm"),
    row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 8),
    show_row_names = FALSE, show_column_names = FALSE,
    use_raster = TRUE, raster_quality = 2, raster_resize_mat = FALSE,
    left_annotation = annotation,
    top_annotation = top, bottom_annotation = bottom,
    column_title = paste(c(region$annotation,
      if (include_haplotype) paste0(cluster_label, "; focal allele annotated per sample"),
      if (split_alleles) "Rows split by focal allele, then cluster",
      if (show_cluster_profiles) "Top: m6A call fraction per cluster (all full-span reads)",
      if (include_haplotype && !nrow(snps)) "No phased heterozygous SNP in this window",
      paste0(region$chr, ":", region$analysis_start, "-", region$analysis_end)), collapse = "\n"),
    column_title_gp = grid::gpar(fontsize = 11))
}

plot_category_composition <- function(res, column, colors, main = NULL,
                                      legend_title = column,
                                      y_label = "Fraction of cluster reads",
                                      label_fun = identity) {
  assignments <- res$assignments
  stopifnot(column %in% names(assignments))
  categories <- label_fun(as.character(assignments[[column]]))
  stopifnot(length(categories) == nrow(assignments),
            !anyNA(categories), all(categories %in% names(colors)))
  assignments$category <- factor(categories, levels = names(colors))
  counts <- table(assignments$cluster)
  ggplot2::ggplot(assignments, ggplot2::aes(cluster, fill = category)) +
    ggplot2::geom_bar(position = "fill", width = .75) +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::scale_x_discrete(labels = function(x) paste0(x, "\n(n=", counts[x], ")")) +
    ggplot2::scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, .25),
      expand = ggplot2::expansion(mult = c(0, .02))) +
    ggplot2::labs(x = "Cluster", y = y_label, fill = legend_title, title = main,
      subtitle = "Each bar sums to 100% of retained reads in that cluster") +
    ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "right")
}

plot_sample_composition <- function(res, sample_colors, main = NULL) {
  plot_category_composition(res, "sample_name", sample_colors, main,
    legend_title = "LCL sample", y_label = "Sample composition within cluster",
    label_fun = function(x) sub("_.*$", "", x))
}


# Fiber-seq panels: every composition denominator is the cluster read count.
fiberseq_category_palette <- function(categories) {
  categories <- sort(unique(as.character(categories)))
  categories <- c(setdiff(categories, c("unphased", "Unknown allele")),
    intersect(c("unphased", "Unknown allele"), categories))
  colors <- setNames(grDevices::hcl.colors(length(categories), "Dynamic"), categories)
  focal <- categories[grepl("^rs[0-9]+: [ACGT]$", categories)]
  colors[focal] <- c("#0072B2", "#D55E00", "#CC79A7", "#009E73")[seq_along(focal)]
  unknown <- intersect(c("unphased", "Unknown allele"), categories)
  colors[unknown] <- "#555555"
  colors
}

fiberseq_legend <- function(colors, title, columns = 4L, labels = names(colors)) {
  d <- data.frame(category = factor(names(colors), levels = names(colors)), x = 1, y = 1)
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y, color = category)) +
    ggplot2::geom_point() + ggplot2::scale_color_manual(values = colors, labels = labels, drop = FALSE) +
    ggplot2::guides(color = ggplot2::guide_legend(title = title, title.position = "top", ncol = columns, byrow = TRUE,
      override.aes = list(size = 3))) + ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8))
  cowplot::get_legend(p)
}

fiberseq_legend_height <- function(legend) {
  grid::convertHeight(sum(legend$heights), "in", valueOnly = TRUE) + .2
}

plot_fiberseq_profiles <- function(tables, region, markers, cluster_colors, sample_colors, allele_colors,
                                  title = region$annotation) {
  allele_heading <- if (identical(region$region_type, "top_asfire_het"))
    "Focal SNP allele" else "Allele / local haplotype"
  rows <- lapply(tables$groups, function(cluster) {
    d <- tables$profiles[tables$profiles$cluster == cluster, ]
    n <- tables$counts$n_reads[match(cluster, tables$counts$cluster)]
    nuc <- d[d$track == "ft_nuc_130-160bp", , drop = FALSE]
    met <- d[d$track == "m6A", , drop = FALSE]
    p <- ggplot2::ggplot(d, ggplot2::aes(pos, fraction)) +
      ggplot2::geom_ribbon(data = nuc, ggplot2::aes(ymin = 0, ymax = fraction),
        fill = "#808080", alpha = .18, color = NA) +
      ggplot2::geom_line(data = nuc, color = "#666666", linewidth = .55) +
      ggplot2::geom_line(data = met, color = cluster_colors[[cluster]], linewidth = .45) +
      ggplot2::scale_x_continuous(limits = c(region$analysis_start, region$analysis_end),
        labels = function(x) format(x, scientific = FALSE, trim = TRUE), expand = ggplot2::expansion(mult = 0)) +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
      ggplot2::labs(title = paste0(cluster, " (n = ", n, ")"), x = paste0(region$chr, " (1-based bp)"), y = "Read fraction") +
      ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "none",
        plot.title = ggplot2::element_text(color = cluster_colors[[cluster]], face = "bold"))
    if (nrow(markers)) {
      p <- p +
        ggplot2::geom_vline(data = markers, ggplot2::aes(xintercept = pos),
          inherit.aes = FALSE, color = "#666666", linetype = "dotted", linewidth = .35)
    }
    bar <- function(tab, colors, heading) {
      tab <- tab[tab$cluster == cluster, ]
      tab$category <- factor(tab$category, levels = names(colors))
      ggplot2::ggplot(tab, ggplot2::aes(x = 1, y = fraction, fill = category)) +
        ggplot2::geom_col(width = .6, position = ggplot2::position_stack(reverse = TRUE)) +
        ggplot2::coord_flip() + ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
        ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1), labels = scales::percent,
          expand = ggplot2::expansion(mult = 0)) +
        ggplot2::scale_x_continuous(breaks = NULL) +
        ggplot2::labs(title = heading, subtitle = paste0("n = ", n), x = NULL, y = "Within cluster") +
        ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "none", panel.grid = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(5.5, 16, 5.5, 5.5))
    }
    cowplot::plot_grid(p, bar(tables$sample, sample_colors, "LCL sample"),
      bar(tables$allele, allele_colors, allele_heading), nrow = 1,
      rel_widths = c(2.7, 1, 1.5), align = "h", axis = "tb")
  })
  marker_text <- if (nrow(markers)) paste(paste0(markers$label, ": ", markers$pos),
    collapse = "; ") else "No focal annotation resolved"
  header <- cowplot::ggdraw() + cowplot::draw_label(paste0(title,
    "\nCluster color: m6A calls; grey line / light grey fill: 130-160 bp nucleosome occupancy\n",
    paste(strwrap(marker_text, width = 155), collapse = "\n")), size = 10)
  sample_legend <- fiberseq_legend(sample_colors, "LCL sample", 6, sub("_.*$", "", names(sample_colors)))
  allele_legend <- fiberseq_legend(allele_colors, allele_heading, 3)
  heights <- c(.85, rep(1.35, length(rows)), fiberseq_legend_height(sample_legend), fiberseq_legend_height(allele_legend))
  plot <- cowplot::plot_grid(plotlist = c(list(header), rows, list(sample_legend, allele_legend)),
                            ncol = 1, rel_heights = heights)
  list(plot = plot, height = sum(heights))
}

plot_umap <- function(embedding, column, colors, title, labels = names(colors),
                           show_legend = FALSE) {
  ggplot2::ggplot(embedding, ggplot2::aes(UMAP1, UMAP2, color = .data[[column]])) +
    ggplot2::geom_point(size = .9, alpha = .75) +
    ggplot2::scale_color_manual(values = colors, labels = labels, breaks = names(colors), drop = FALSE) +
    ggplot2::coord_equal() + ggplot2::labs(title = title, color = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = if (show_legend) "right" else "none")
}

plot_knn_graph <- function(result, main = NULL, seed = NULL, vertex_size = 1.5,
                           edge_width = 0.3, edge_alpha = 0.2, color_by = c("cluster", "allele")) {
  color_by <- match.arg(color_by)
  graph <- result$graph
  assignments <- result$assignments
  if (!igraph::is_igraph(graph)) stop("The clustering result must contain its saved KNN graph")
  stopifnot(all(c("RID", "cluster") %in% names(assignments)),
    !igraph::is_directed(graph), igraph::vcount(graph) > 0L,
    !anyNA(assignments$RID), !anyDuplicated(assignments$RID), !anyNA(assignments$cluster))
  read_ids <- igraph::V(graph)$name
  stopifnot(length(read_ids) == igraph::vcount(graph), !anyNA(read_ids),
    !anyDuplicated(read_ids), setequal(read_ids, as.character(assignments$RID)))
  assignments <- assignments[match(read_ids, assignments$RID), , drop = FALSE]
  clusters <- if (is.factor(assignments$cluster)) levels(droplevels(assignments$cluster)) else
    unique(as.character(assignments$cluster))
  colors <- cluster_id_palette(clusters)
  weights <- igraph::edge_attr(graph, "weight")
  if (length(weights)) stopifnot(all(is.finite(weights)), all(weights > 0))
  if (is.null(seed)) seed <- if (is.null(result$params$seed)) 1L else result$params$seed
  stopifnot(length(seed) == 1L, is.finite(seed))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  coordinates <- igraph::layout_with_fr(graph, weights = weights)
  nodes <- data.frame(RID = read_ids, layout_x = coordinates[, 1], layout_y = coordinates[, 2],
    cluster = factor(assignments$cluster, levels = clusters))
  nodes$color_group <- nodes$cluster
  color_levels <- clusters
  legend_title <- "Leiden cluster"
  if (color_by == "allele") {
    summary <- knn_allele_summary(result)
    color_levels <- summary$allele_levels
    nodes$color_group <- factor(assignments$allele_label, levels = color_levels)
    colors <- fiberseq_category_palette(color_levels)
    legend_title <- "Focal SNP allele"
  }
  endpoints <- igraph::as_edgelist(graph, names = FALSE)
  edges <- data.frame(from_x = coordinates[endpoints[, 1], 1],
    from_y = coordinates[endpoints[, 1], 2], to_x = coordinates[endpoints[, 2], 1],
    to_y = coordinates[endpoints[, 2], 2])
  counts <- table(nodes$color_group)
  labels <- paste0(color_levels, " (n=", as.integer(counts[color_levels]), ")")
  if (is.null(main)) main <- paste(c(result$region$annotation,
    paste("Manhattan KNN graph colored by", if (color_by == "allele") "focal SNP allele" else "Leiden cluster")),
    collapse = "\n")
  neighbors <- if (is.null(result$params$k_eff)) result$params$k_neighbors else result$params$k_eff
  settings <- c(paste(nrow(nodes), "reads"), paste(nrow(edges), "edges"),
    if (!is.null(neighbors)) paste0("k = ", neighbors),
    if (!is.null(result$params$resolution)) paste0("resolution = ", result$params$resolution))
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges,
      ggplot2::aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
      color = "grey60", linewidth = edge_width, alpha = edge_alpha) +
    ggplot2::geom_point(data = nodes, ggplot2::aes(layout_x, layout_y, color = color_group),
      size = vertex_size) +
    ggplot2::scale_color_manual(values = colors, breaks = color_levels, labels = labels, drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = main, subtitle = paste(settings, collapse = " | "), color = legend_title) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      plot.margin = ggplot2::margin(12, 12, 12, 12))
}

knn_allele_summary <- function(result) {
  graph <- result$graph
  assignments <- result$assignments
  focal <- result$focal
  stopifnot(igraph::is_igraph(graph), !igraph::is_directed(graph), igraph::is_simple(graph),
    nrow(assignments) > 1L, !is.null(focal), nrow(focal) == 1L,
    all(c("RID", "cluster", "allele_label", "allele_status") %in% names(assignments)),
    !anyNA(assignments$RID), !anyDuplicated(assignments$RID), !anyNA(assignments$cluster),
    !anyNA(assignments$allele_label), all(assignments$allele_status == "phased_focal_genotype"))
  alleles <- sort(paste0(focal$variant_id, ": ", c(focal$ref, focal$alt)))
  stopifnot(length(alleles) == 2L, length(unique(alleles)) == 2L,
    all(grepl("^rs[0-9]+: [ACGT]$", alleles)), all(assignments$allele_label %in% alleles))
  read_ids <- igraph::V(graph)$name
  stopifnot(length(read_ids) == igraph::vcount(graph), !anyNA(read_ids),
    !anyDuplicated(read_ids), setequal(read_ids, as.character(assignments$RID)))
  assignments <- assignments[match(read_ids, assignments$RID), , drop = FALSE]
  clusters <- if (is.factor(assignments$cluster)) levels(droplevels(assignments$cluster)) else
    unique(as.character(assignments$cluster))
  stopifnot(!"All reads" %in% clusters)
  allele <- factor(assignments$allele_label, levels = alleles)
  composition <- as.data.frame(table(group = factor(assignments$cluster, levels = clusters),
    allele = allele), stringsAsFactors = FALSE)
  names(composition)[3] <- "n_reads"
  overall <- data.frame(group = "All reads", allele = alleles, n_reads = as.integer(table(allele)))
  composition <- dplyr::bind_rows(overall, composition)
  composition$group <- factor(composition$group, levels = c("All reads", clusters))
  composition$allele <- factor(composition$allele, levels = alleles)
  composition$total_reads <- ave(composition$n_reads, composition$group, FUN = sum)
  composition$fraction <- composition$n_reads / composition$total_reads
  endpoints <- igraph::as_edgelist(graph, names = FALSE)
  source_nodes <- c(endpoints[, 1], endpoints[, 2])
  target_nodes <- c(endpoints[, 2], endpoints[, 1])
  degree <- tabulate(source_nodes, nbins = length(read_ids))
  same <- tabulate(source_nodes[allele[source_nodes] == allele[target_nodes]], nbins = length(read_ids))
  reads <- data.frame(RID = read_ids, cluster = assignments$cluster, allele = allele,
    degree = degree, same_allele_neighbors = same, opposite_allele_neighbors = degree - same)
  reads$same_fraction <- ifelse(degree > 0L, same / degree, NA_real_)
  reads$opposite_fraction <- ifelse(degree > 0L, (degree - same) / degree, NA_real_)
  connectivity <- dplyr::bind_rows(lapply(alleles, function(label) {
    selected <- reads$allele == label
    connected <- selected & reads$degree > 0L
    data.frame(allele = label, n_reads = sum(selected), n_reads_with_neighbors = sum(connected),
      n_isolated_reads = sum(selected & reads$degree == 0L),
      same_fraction = if (any(connected)) mean(reads$same_fraction[connected]) else NA_real_,
      opposite_fraction = if (any(connected)) mean(reads$opposite_fraction[connected]) else NA_real_,
      expected_same_fraction = if (any(selected)) (sum(selected) - 1L) / (nrow(reads) - 1L) else NA_real_)
  }))
  connectivity$allele <- factor(connectivity$allele, levels = alleles)
  list(allele_levels = alleles, composition = composition, reads = reads, connectivity = connectivity)
}


save_fiberseq_plots <- function(result, footprints, output_dir, sample_colors,
                               embedding = NULL, example_only = FALSE,
                               plot_writer = NULL) {
  output_dir <- lcl_output_path(output_dir)
  result <- fiberseq_display_result(result)
  allele_heading <- if (identical(result$region$region_type, "top_asfire_het"))
    "Focal SNP allele" else "Allele / sample-local haplotype"
  plot_dir <- file.path(output_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- fiberseq_tables(result, footprints)
  cluster_colors <- cluster_id_palette(tables$groups)
  allele_colors <- fiberseq_category_palette(result$allele_display_levels)
  stopifnot(all(result$assignments$sample_name %in% names(sample_colors)))
  save_plot <- function(plot, name, width, height) {
    if (example_only && name != "fiberseq_haplotype_example") return(invisible(NULL))
    if (!is.null(plot_writer)) return(invisible(plot_writer(plot, name, width, height)))
    ggplot2::ggsave(file.path(plot_dir, paste0(name, ".pdf")),
      plot, width = width, height = height, limitsize = FALSE, bg = "white")
  }
  profiles <- plot_fiberseq_profiles(tables, result$region, result$markers,
    cluster_colors, sample_colors, allele_colors,
    result$region$annotation)
  save_plot(profiles$plot, "fiberseq_profiles_composition", 16, profiles$height)
  if (is.null(embedding)) embedding <- fiberseq_umap(result)
  stopifnot(identical(embedding$RID, result$assignments$RID))
  embedding <- data.frame(result$assignments, embedding[, c("UMAP1", "UMAP2")])
  embedding$display_cluster <- factor(result$assignments$cluster, levels = tables$groups)
  cluster_labels <- paste0(tables$groups, " (n=", tables$counts$n_reads, ")")
  pc <- plot_umap(embedding, "display_cluster", cluster_colors, "UMAP: cluster")
  ps <- plot_umap(embedding, "sample_name", sample_colors, "UMAP: LCL sample")
  pa <- plot_umap(embedding, "allele_display", allele_colors, paste0("UMAP: ", allele_heading))
  counts <- tables$counts
  counts$cluster <- factor(counts$cluster, levels = tables$groups)
  counts$label <- sprintf("%s\nn=%d (%.1f%%)", counts$cluster, counts$n_reads, 100 * counts$proportion)
  rgb <- grDevices::col2rgb(cluster_colors) / 255
  text_colors <- setNames(ifelse(colSums(rgb * c(.299, .587, .114)) < .5, "white", "black"), names(cluster_colors))
  pie <- ggplot2::ggplot(counts, ggplot2::aes(x = "", y = proportion, fill = cluster)) +
    ggplot2::geom_col(width = 1, color = "white", position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(proportion >= .03, sprintf("%.1f%%", 100 * proportion), ""), color = cluster),
      position = ggplot2::position_stack(vjust = .5, reverse = TRUE), size = 3) +
    ggplot2::scale_color_manual(values = text_colors, guide = "none") +
    ggplot2::coord_polar(theta = "y") + ggplot2::scale_fill_manual(values = cluster_colors, labels = counts$label) +
    ggplot2::labs(title = paste0("Cluster proportions: ", nrow(embedding), " reads"), fill = NULL) +
    ggplot2::theme_void() + ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size = 8))
  save_plot(pie, "fiberseq_cluster_proportions", 8, 5)
  legends <- list(fiberseq_legend(cluster_colors, "Cluster", 5, cluster_labels),
    fiberseq_legend(sample_colors, "LCL sample", 6, sub("_.*$", "", names(sample_colors))),
    fiberseq_legend(allele_colors, allele_heading, 3))
  lh <- vapply(legends, fiberseq_legend_height, numeric(1))
  heading <- cowplot::ggdraw() + cowplot::draw_label(result$region$annotation, size = 12)
  overview <- cowplot::plot_grid(plotlist = c(list(heading, cowplot::plot_grid(pc, ps, pa, pie, ncol = 2)), legends),
                                ncol = 1, rel_heights = c(.45, 8, lh))
  overview_height <- 8.45 + sum(lh)
  save_plot(overview, "fiberseq_umap_overview", 16, overview_height)
  # One PDF page combines the same UMAP coordinates, cluster pie, aggregate
  # accessibility/nucleosomes, and within-cluster sample/actual-allele bars.
  example <- cowplot::plot_grid(overview, profiles$plot, ncol = 1,
    rel_heights = c(overview_height, profiles$height))
  save_plot(example, "fiberseq_haplotype_example", 16, overview_height + profiles$height)
  invisible(list(tables = tables, embedding = embedding))
}


# ---------------------------------------------------------------------------
# Per-cluster mean m6A call at each m6A site, computed at site resolution
# straight from met_mat 
# ---------------------------------------------------------------------------
cluster_site_profiles <- function(res, met_mat) {
  M  <- as.matrix(met_mat)
  df <- res$assignments
  df <- df[df$RID %in% rownames(M), ]
  do.call(rbind, lapply(levels(df$cluster), function(cl) {
    rids <- df$RID[df$cluster == cl]
    met  <- colMeans(M[rids, , drop = FALSE], na.rm = TRUE)
    met[is.nan(met)] <- NA
    data.frame(cluster = cl, pos = as.numeric(colnames(M)), met = met,
               n_reads = length(rids), row.names = NULL)
  }))
}

extracted_path <- function(sample_table, chromosome, feature) {
  file.path(sample_table$fire_dir, "extracted_results", paste0(feature, "_by_chr"),
            paste0(sample_table$sample_name, ".ft_extracted_", feature, ".", chromosome, ".bed.gz"))
}

# Collapse sample-level variants to exact in-window SNP positions and REF/ALT bases.
window_snps <- function(variants, region) {
  empty <- data.frame(pos = integer(), label = character())
  if (is.null(variants) || !nrow(variants)) return(empty)
  snps <- variants[variants$chr == region$chr & variants$pos >= region$analysis_start &
    variants$pos <= region$analysis_end & nchar(variants$ref) == 1L &
    grepl("^[ACGT](,[ACGT])*$", variants$alt), , drop = FALSE]
  if (!nrow(snps)) return(empty)
  snps$id <- ifelse(is.na(snps$variant_id) | snps$variant_id == ".", "SNP", snps$variant_id)
  snps$label <- paste0(snps$id, " ", snps$pos, " ", snps$ref, ">", snps$alt)
  labels <- tapply(snps$label, snps$pos, function(x) paste(unique(x), collapse = "; "))
  data.frame(pos = as.integer(names(labels)), label = unname(labels))
}


# Plot only the observed focal SNP / or an annotated TSS
focal_markers <- function(result) {
  focal <- result$focal
  if (!is.null(focal) && nrow(focal))
    return(data.frame(pos = focal$pos, end = focal$pos,
      label = paste0(focal$variant_id, " ", focal$ref, ">", focal$alt)))
  region <- result$region
  if (!is.null(region$focal_pos) && !is.na(region$focal_pos))
    return(data.frame(pos = region$focal_pos, end = region$focal_pos, label = region$focal_snp))
  if (!is.null(region$tss) && !is.na(region$tss))
    return(data.frame(pos = region$tss, end = region$tss, label = paste(region$gene, "TSS")))
  data.frame(pos = integer(), end = integer(), label = character())
}

# Keep detailed phasing provenance, but combine unknown focal alleles for display.
fiberseq_display_result <- function(result) {
  result$markers <- focal_markers(result)
  a <- result$assignments
  a$allele_display <- a$allele_label
  anchored <- !is.null(result$region$focal_pos) && !is.na(result$region$focal_pos)
  if (anchored) {
    known <- a$allele_status == "phased_focal_genotype"
    strict <- identical(result$region$region_type, "top_asfire_het")
    if (strict && (!all(known) || !all(a$focal_genotype %in% c("0|1", "1|0"))))
      stop("Top AS-FIRE plots require only resolved phased heterozygous reads")
    a$allele_display[!known] <- "Unknown allele"
    bases <- a$allele_label[known]
    v <- result$variants
    if (!is.null(v) && nrow(v)) {
      v <- v[v$chr == result$region$chr & v$pos == result$region$focal_pos, , drop = FALSE]
      bases <- c(bases, paste0(if (is.null(result$region$focal_snp)) result$region$region_id else result$region$focal_snp, ": ",
        unique(c(v$ref, unlist(strsplit(v$alt, ",", fixed = TRUE))))))
    }
    bases <- sort(unique(bases[grepl("^rs[0-9]+: [ACGT]$", bases)]))
    if (strict) {
      bases <- sort(paste0(result$region$focal_snp, ": ", c(result$region$ref, result$region$alt)))
      stopifnot(length(unique(bases)) == 2L, all(a$allele_label %in% bases))
    }
    result$allele_display_levels <- if (strict) bases else c(bases, "Unknown allele")
  } else {
    result$allele_display_levels <- if (is.null(result$allele_levels))
      sort(unique(a$allele_label)) else result$allele_levels
  }
  a$allele_display <- factor(a$allele_display, levels = result$allele_display_levels)
  stopifnot(!anyNA(a$allele_display))
  result$assignments <- a
  result
}

fiberseq_tables <- function(result, footprints) {
  result <- fiberseq_display_result(result)
  a <- result$assignments
  groups <- if (is.factor(a$cluster)) levels(droplevels(a$cluster)) else unique(a$cluster)
  group <- factor(a$cluster, levels = groups)
  counts <- data.frame(cluster = groups, n_reads = as.integer(table(group)), total_reads = nrow(a))
  counts$proportion <- counts$n_reads / counts$total_reads
  composition <- function(column) {
    tab <- as.data.frame(table(cluster = group, category = a[[column]]), stringsAsFactors = FALSE)
    names(tab)[3] <- "n_reads"
    tab$cluster <- factor(tab$cluster, levels = groups)
    tab$cluster_reads <- counts$n_reads[match(tab$cluster, counts$cluster)]
    tab$fraction <- tab$n_reads / tab$cluster_reads
    stopifnot(all(abs(tapply(tab$fraction, tab$cluster, sum) - 1) < 1e-10))
    tab
  }
  view <- result
  view$assignments$cluster <- group
  nuc <- footprint_profiles(footprints, "ft_nuc_130-160bp", view$assignments, view$region)
  met <- cluster_site_profiles(view, view$site_met_mat)
  profiles <- dplyr::bind_rows(data.frame(cluster = met$cluster, pos = met$pos,
    track = "m6A", fraction = met$met, n_reads = met$n_reads), nuc)
  profiles$cluster <- factor(profiles$cluster, levels = groups)
  profiles$chr <- result$region$chr
  profiles$coordinate_system <- "1-based inclusive"
  stopifnot(all(is.finite(profiles$fraction)), all(profiles$fraction >= 0 & profiles$fraction <= 1))
  list(counts = counts, sample = composition("sample_name"), allele = composition("allele_display"),
       profiles = profiles, groups = groups)
}

fiberseq_umap <- function(result, seed = 1L, n_neighbors = 15L, min_dist = 0.1) {
  stopifnot(requireNamespace("uwot", quietly = TRUE))
  a <- result$assignments
  x <- as.matrix(result$feat_mat[a$RID, , drop = FALSE])
  stopifnot(nrow(x) >= 3L, !anyNA(x), !anyDuplicated(a$RID))
  set.seed(seed)
  xy <- uwot::umap(x, metric = "manhattan", n_neighbors = min(n_neighbors, nrow(x) - 1L),
    min_dist = min_dist, n_threads = 1L, n_sgd_threads = 1L, init = "random", verbose = FALSE)
  data.frame(a, UMAP1 = xy[, 1], UMAP2 = xy[, 2])
}


# ---- Single-molecule footprint processing ----
# Process LCL nucleosome/TF footprints on fixed m6A-defined read clusters.
# BED12 and shared data utilities are defined above in this file.
# Footprint display tracks only; the clustering input remains the m6A matrix.
LCL_FOOTPRINT_TRACKS <- c("ft_nuc_130-160bp", "FiberHMM_10-30bp",
                          "FiberHMM_40-60bp", "FiberHMM_60-80bp")

m6a_intervals <- function(result) {
  sites <- Matrix::summary(as(result$site_met_mat, "dgCMatrix"))
  sites <- sites[sites$x > 0, , drop = FALSE]
  positions <- as.integer(colnames(result$site_met_mat)[sites$j])
  identifiers <- rownames(result$site_met_mat)[sites$i]
  data.frame(RID = identifiers,
    original_RID = result$assignments$original_RID[match(identifiers, result$assignments$RID)],
    chr = rep(result$region$chr, length(identifiers)), start = positions - 1L,
    end = positions, size = rep(1L, length(identifiers)), track = rep("m6A", length(identifiers)))
}

plot_anchor <- function(region) {
  promoter <- region$region_type == "promoter" && !is.na(region$tss)
  focal <- !is.null(region$focal_pos) && !is.na(region$focal_pos)
  anchor <- if (focal) region$focal_pos else if (promoter) region$tss else floor(mean(c(region$analysis_start, region$analysis_end)))
  direction <- if (promoter && region$strand == "-") -1L else 1L
  bounds <- sort(direction * (c(region$analysis_start, region$analysis_end) - anchor))
  list(anchor = anchor, direction = direction, left = bounds[1], right = bounds[2],
    x_label = if (focal) paste0("Position relative to ", region$focal_snp, " (bp)") else if (promoter) "Position relative to canonical TSS (bp)" else "Position relative to region centre (bp)")
}

prepare_read_tracks <- function(result, records, sample_colors) {
  anchor <- plot_anchor(result$region)
  reads <- result$assignments
  reads <- reads[order(reads$cluster, reads$start, reads$RID), , drop = FALSE]
  reads$row <- as.integer(ave(seq_len(nrow(reads)), reads$cluster, FUN = seq_along))
  reads$sample_label <- sub("_.*$", "", reads$sample_name)
  stopifnot(all(reads$sample_label %in% names(sample_colors)))
  if (!"haplotype" %in% names(reads)) reads$haplotype <- "pooled"
  relative_start <- anchor$direction * (reads$start - anchor$anchor)
  relative_end <- anchor$direction * (reads$end - anchor$anchor)
  reads$left <- pmax(pmin(relative_start, relative_end), anchor$left)
  reads$right <- pmin(pmax(relative_start, relative_end), anchor$right)
  matched <- match(records$RID, reads$RID)
  records <- records[!is.na(matched), , drop = FALSE]
  matched <- matched[!is.na(matched)]
  records$row <- reads$row[matched]
  records$cluster <- reads$cluster[matched]
  relative_start <- anchor$direction * (records$start + 1L - anchor$anchor)
  relative_end <- anchor$direction * (records$end - anchor$anchor)
  records$left <- pmax(pmin(relative_start, relative_end) - 0.5, anchor$left - 0.5)
  records$right <- pmin(pmax(relative_start, relative_end) + 0.5, anchor$right + 0.5)
  list(reads = reads, features = records, anchor = anchor, region = result$region)
}


extract_nucleosomes <- function(sample_table, region, assignments, min_size = 130L, max_size = 160L) {
  paths <- extracted_path(sample_table, region$chr, "nuc")
  result <- lapply(seq_len(nrow(sample_table)), function(sample_index) {
    query <- GenomicRanges::GRanges(region$chr, IRanges::IRanges(region$analysis_start, region$analysis_end))
    bed <- read_ft_bed12(paths[sample_index], query, longest_alignment = TRUE)
    if (!nrow(bed)) return(NULL)
    bed$RID <- paste(sample_table$sample_name[sample_index], bed$RID, sep = "::")
    bed <- bed[bed$RID %in% assignments$RID, , drop = FALSE]
    if (!nrow(bed)) return(NULL)
    blocks <- convert_ft_bed12_to_bed6(bed)
    blocks$original_RID <- assignments$original_RID[match(blocks$RID, assignments$RID)]
    blocks$size <- blocks$end - blocks$start
    blocks <- blocks[blocks$size >= min_size & blocks$size <= max_size &
      blocks$start < region$end & blocks$end > region$start, , drop = FALSE]
    blocks[, c("RID", "original_RID", "start", "end", "size"), drop = FALSE]
  })
  result <- dplyr::bind_rows(result)
  if (!ncol(result)) result <- data.frame(RID = character(), original_RID = character(),
    start = integer(), end = integer(), size = integer())
  result$chr <- rep(region$chr, nrow(result))
  result$track <- rep(paste0("ft_nuc_", min_size, "-", max_size, "bp"), nrow(result))
  result
}

cache_footprint_tracks <- function(regions, tracks_dir, cache_dir, workers = 2L, reuse = TRUE) {
  cache_dir <- lcl_output_path(cache_dir)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- unlist(lapply(unique(regions$chr), function(chromosome) {
    files <- list.files(file.path(tracks_dir, chromosome), pattern = "_(10-30|40-60|60-80)bp_fps\\.bed\\.gz$", full.names = TRUE)
    if (!length(files)) stop("No footprint tracks for ", chromosome)
    files
  }), use.names = FALSE)
  jobs <- lapply(paths, function(path) {
    chromosome <- basename(dirname(path))
    list(path = path, regions = regions[regions$chr == chromosome, c("chr", "start", "end")])
  })
  extract_job <- function(job) {
    signature <- list(version = 1L, regions = job$regions, path = job$path)
    cache <- file.path(cache_dir, paste0(basename(job$path), ".rds"))
    if (reuse && file.exists(cache)) {
      previous <- readRDS(cache)$signature
      previous_path <- if (!is.null(previous$path)) previous$path else previous$source$path
      if (identical(previous$version, signature$version) &&
          identical(previous$regions, signature$regions) &&
          identical(previous_path, signature$path)) return(cache)
    }
    conditions <- sprintf("($1 == \"%s\" && $2 < %.0f && $3 > %.0f)",
                            job$regions$chr, job$regions$end, job$regions$start)
    program <- paste0("BEGIN { FS = OFS = \"\\t\" } (", paste(conditions, collapse = " || "), ") { print }")
    temporary <- tempfile(tmpdir = cache_dir, fileext = ".bed")
    on.exit(unlink(temporary), add = TRUE)
    command <- paste("set -o pipefail; gzip -cd", shQuote(job$path), "| awk", shQuote(program), ">", shQuote(temporary))
    message("Extracting regional footprints: ", basename(job$path))
    status <- system2("/bin/bash", c("-c", shQuote(command)))
    if (status != 0L) stop("Footprint extraction failed: ", job$path)
    records <- data.frame(chr = character(), start = integer(), end = integer(), original_RID = character())
    if (file.info(temporary)$size > 0) {
      records <- data.table::fread(temporary, header = FALSE, data.table = FALSE)
      if (ncol(records) != 4L) stop("Expected BED4 footprints: ", job$path)
      names(records) <- c("chr", "start", "end", "original_RID")
    }
    records$size <- records$end - records$start
    records$track <- rep(sub("^combined_chr[^_]+_(.*)_fps\\.bed\\.gz$", "FiberHMM_\\1", basename(job$path)), nrow(records))
    saveRDS(list(records = records, signature = signature), cache)
    cache
  }
  caches <- parallel::mclapply(jobs, extract_job, mc.cores = workers, mc.preschedule = FALSE)
  if (any(vapply(caches, inherits, logical(1), "try-error"))) stop("One or more footprint extractions failed")
  unlist(caches, use.names = FALSE)
}

region_footprints <- function(caches, region, assignments) {
  if (anyDuplicated(assignments$original_RID)) {
    stop("Combined BED4 lacks sample IDs: ambiguous original read names in ", region$region_id)
  }
  caches <- caches[startsWith(basename(caches), paste0("combined_", region$chr, "_"))]
  pieces <- lapply(caches, function(cache) {
    records <- readRDS(cache)$records
    selected <- records$start < region$end & records$end > region$start
    records <- records[selected, , drop = FALSE]
    matched <- match(records$original_RID, assignments$original_RID)
    records$RID <- assignments$RID[matched]
    records[!is.na(matched), , drop = FALSE]
  })
  list(records = dplyr::bind_rows(pieces),
       tracks = vapply(caches, function(cache) {
         sub("^combined_chr[^_]+_(.*)_fps\\.bed\\.gz.rds$", "FiberHMM_\\1", basename(cache))
       }, character(1)))
}

occupancy_matrix <- function(records, assignments, region) {
  occupancy <- matrix(0L, nrow(assignments), region$width,
                       dimnames = list(assignments$RID, seq.int(region$analysis_start, region$analysis_end)))
  if (!nrow(records)) return(occupancy)
  for (record_index in seq_len(nrow(records))) {
    read_index <- match(records$RID[record_index], assignments$RID)
    if (is.na(read_index)) next
    left <- max(records$start[record_index], region$start) - region$start + 1L
    right <- min(records$end[record_index], region$end) - region$start
    if (left <= right) occupancy[read_index, seq.int(left, right)] <- 1L
  }
  occupancy
}

footprint_profiles <- function(records, tracks, assignments, region) {
  dplyr::bind_rows(lapply(tracks, function(track) {
    occupancy <- occupancy_matrix(records[records$track == track, ], assignments, region)
    dplyr::bind_rows(lapply(levels(assignments$cluster), function(cluster) {
      selected <- assignments$cluster == cluster
      data.frame(cluster = cluster, pos = as.integer(colnames(occupancy)), track = track,
                 fraction = colMeans(occupancy[selected, , drop = FALSE]), n_reads = sum(selected))
    }))
  }))
}


# Footprints use per-base occupancy; m6A profiles use the observed-site matrix
signal_profiles <- function(result, footprints) {
  footprints <- footprints[footprints$track %in% LCL_FOOTPRINT_TRACKS, , drop = FALSE]
  records <- dplyr::bind_rows(m6a_intervals(result), footprints)
  fp <- footprint_profiles(footprints, LCL_FOOTPRINT_TRACKS, result$assignments, result$region)
  met <- cluster_site_profiles(result, result$site_met_mat)
  met <- data.frame(cluster = met$cluster, pos = met$pos, track = "m6A",
    fraction = met$met, n_reads = met$n_reads)
  list(records = records, tracks = c("m6A", LCL_FOOTPRINT_TRACKS),
       profiles = dplyr::bind_rows(met, fp))
}
