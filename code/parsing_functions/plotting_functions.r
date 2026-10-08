# Shared plotting palettes 
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

cluster_id_palette <- function(levels) {
  idx <- suppressWarnings(as.integer(sub("^cluster", "", levels)))
  if (anyNA(idx) || any(!grepl("^cluster[0-9]+$", levels))) idx <- seq_along(levels)
  cols <- if (max(idx) <= length(LEIDEN_CLUSTER_COLORS)) LEIDEN_CLUSTER_COLORS[idx] else
    grDevices::colorRampPalette(LEIDEN_CLUSTER_COLORS)(max(idx))[idx]
  setNames(cols, levels)
}



timepoint_palette <- function(x) {
  labels <- unique(as.character(x))
  minutes <- trimws(sub("min$", "", sub("^LPS_", "", labels)))
  colors <- unname(LEIDEN_TIMEPOINT_COLORS[paste0("LPS_", minutes)])
  unknown <- is.na(labels) | is.na(colors)
  if (any(unknown)) {
    shown <- ifelse(is.na(labels[unknown]), "<NA>", paste0('"', labels[unknown], '"'))
    stop("No LEIDEN_TIMEPOINT_COLORS entry for: ", paste(shown, collapse = ", "), call. = FALSE)
  }
  setNames(colors, labels)
}


theme_fiberseq <- function(base = c("classic", "bw", "cowplot"),
                           base_size = 9, base_family = "", overrides = list()) {
  base <- match.arg(base)
  value <- switch(base,
    classic = ggplot2::theme_classic(base_size = base_size, base_family = base_family),
    bw = ggplot2::theme_bw(base_size = base_size, base_family = base_family),
    cowplot = cowplot::theme_cowplot(font_size = base_size, font_family = base_family))
  if (length(overrides)) value <- value + do.call(ggplot2::theme, overrides)
  value
}


save_figure <- function(plot, path, width, height,
                        device = NULL, units = "in", dpi = 300,
                        bg = NULL, limitsize = TRUE,
                        draw = NULL, after_draw = NULL,
                        preview = NULL, on_saved = NULL, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  render <- function(filename, preview_image = FALSE, image_dpi = dpi, image_bg = bg) {
    if (is.null(draw)) {
      ggplot2::ggsave(filename, plot = plot, width = width, height = height,
        device = if (preview_image) "png" else device, units = units,
        dpi = image_dpi, bg = image_bg, limitsize = limitsize, ...)
    } else {
      stopifnot(is.function(draw))
      inches <- switch(units, "in" = 1, "cm" = 1 / 2.54, "mm" = 1 / 25.4,
                       "px" = 1 / image_dpi, stop("Unsupported units: ", units))
      background <- if (is.null(image_bg)) "white" else image_bg
      if (preview_image) {
        grDevices::png(filename, width = width * inches, height = height * inches,
          units = "in", res = image_dpi, type = "cairo", bg = background)
      } else {
        if (tolower(tools::file_ext(filename)) != "pdf")
          stop("A draw callback requires a PDF output")
        pdf_device <- if (is.null(device)) grDevices::pdf else device
        if (!is.function(pdf_device)) stop("For a draw callback, device must be a function")
        pdf_args <- list(filename, width = width * inches, height = height * inches, ...)
  
        if (!is.null(image_bg)) pdf_args$bg <- image_bg
        do.call(pdf_device, pdf_args)
      }
      opened_device <- grDevices::dev.cur()
      tryCatch({
        draw(plot)
        if (!is.null(after_draw)) after_draw(plot)
      }, finally = grDevices::dev.off(opened_device))
    }
  }
  render(path)
  if (!is.null(preview)) {
    temporary <- is.null(preview$path)
    preview_path <- if (temporary) tempfile(fileext = ".png") else preview$path
    if (temporary) on.exit(unlink(preview_path), add = TRUE)
    dir.create(dirname(preview_path), recursive = TRUE, showWarnings = FALSE)
    render(preview_path, preview_image = TRUE,
      image_dpi = if (is.null(preview$dpi)) dpi else preview$dpi,
      image_bg = if (is.null(preview$bg)) bg else preview$bg)
    if (!is.null(preview$embed)) preview$embed(preview_path)
  }
  if (!is.null(on_saved)) on_saved(path)
  invisible(path)
}


plot_stacked_proportion <- function(df, x, fill, colors,
                                    y = NULL, weight = NULL,
                                    position = c("fill", "stack", "dodge"),
                                    reverse = FALSE, horizontal = FALSE,
                                    width = 0.75, border = NULL,
                                    labels = list(), totals = NULL,
                                    scales = list(), legend = list(), theme = NULL) {
  position <- match.arg(position)
  mapping <- ggplot2::aes(x = .data[[x]], fill = .data[[fill]])
  if (!is.null(y)) mapping$y <- ggplot2::aes(y = .data[[y]])$y
  if (!is.null(weight)) mapping$weight <- ggplot2::aes(weight = .data[[weight]])$weight
  placement <- switch(position,
    fill = ggplot2::position_fill(reverse = reverse),
    stack = ggplot2::position_stack(reverse = reverse),
    dodge = ggplot2::position_dodge())
  args <- c(list(width = width, position = placement), border)
  p <- ggplot2::ggplot(df, mapping) +
    do.call(if (is.null(y)) ggplot2::geom_bar else ggplot2::geom_col, args)
  inside <- labels$inside
  labels$inside <- NULL
  if (!is.null(inside)) {
    text <- as.data.frame(df)
    text$.proportion_label <- ifelse(text[[y]] >= inside$min_fraction,
      scales::percent(text[[y]], accuracy = inside$accuracy), "")
    p <- p + ggplot2::geom_text(data = text, ggplot2::aes(label = .proportion_label),
      position = ggplot2::position_stack(vjust = inside$vjust, reverse = reverse),
      color = inside$color, size = inside$size)
  }
  if (!is.null(totals)) {
    text <- as.data.frame(totals$data)
    text$.total_x <- text[[totals$x]]
    text$.total_y <- totals$y
    text$.total_label <- text[[totals$label]]
    p <- p + do.call(ggplot2::geom_text, c(list(data = text,
      mapping = ggplot2::aes(x = .total_x, y = .total_y, label = .total_label),
      inherit.aes = FALSE), totals$style))
  }
  if (horizontal) p <- p + ggplot2::coord_flip()
  p <- p + do.call(ggplot2::scale_fill_manual, c(list(values = colors), legend))
  p <- p + scales + do.call(ggplot2::labs, labels)
  if (!is.null(theme)) p <- p + theme
  p
}


plot_interval_track <- function(df, window, fill_col = NULL, highlight = NULL,
                                start_col = "start", end_col = "end", row_col = NULL,
                                geometry = c("rect", "segment"),
                                colors = NULL, color_scale = NULL,
                                height = 0.7, linewidth = 1,
                                labels = list(), scales = list(), theme = NULL,
                                clip_intervals = FALSE, color = NULL, lineend = "butt") {
  geometry <- match.arg(geometry)
  d <- as.data.frame(df)
  d$.interval_start <- d[[start_col]]
  d$.interval_end <- d[[end_col]]
  d$.interval_row <- if (is.null(row_col)) 1 else d[[row_col]]
  if (clip_intervals) {
    d$.interval_start <- pmax(d$.interval_start, window[1])
    d$.interval_end <- pmin(d$.interval_end, window[2])
  }
  if (geometry == "rect") {
    mapping <- ggplot2::aes(xmin = .interval_start, xmax = .interval_end,
      ymin = .interval_row - height / 2, ymax = .interval_row + height / 2)
    if (!is.null(fill_col)) mapping$fill <- ggplot2::aes(fill = .data[[fill_col]])$fill
    layer <- do.call(ggplot2::geom_rect, c(list(mapping = mapping),
      if (!is.null(color)) list(colour = color)))
  } else {
    mapping <- ggplot2::aes(x = .interval_start, xend = .interval_end,
      y = .interval_row, yend = .interval_row)
    if (!is.null(fill_col)) mapping$colour <- ggplot2::aes(colour = .data[[fill_col]])$colour
    layer <- do.call(ggplot2::geom_segment, c(list(mapping = mapping,
      linewidth = linewidth, lineend = lineend), if (!is.null(color)) list(colour = color)))
  }
  p <- ggplot2::ggplot(d) + layer
  if (!is.null(color_scale)) {
    p <- p + color_scale
  } else if (!is.null(colors)) {
    p <- p + if (geometry == "rect") ggplot2::scale_fill_manual(values = colors) else
      ggplot2::scale_colour_manual(values = colors)
  }
  if (!is.null(highlight)) {
    h <- highlight
    band <- ggplot2::annotate("rect", xmin = h$start, xmax = h$end,
      ymin = -Inf, ymax = Inf, fill = h$fill, alpha = h$alpha)
    if (identical(h$position, "under")) p$layers <- c(list(band), p$layers) else p <- p + band
  }
  coordinates <- scales$coord
  scales$coord <- NULL
  if (is.null(coordinates)) coordinates <- ggplot2::coord_cartesian(xlim = window)
  p <- p + scales + coordinates + do.call(ggplot2::labs, labels)
  if (!is.null(theme)) p <- p + theme
  p
}


plot_group_profile <- function(profile, window = NULL, group_col = NULL,
                               style = c("line", "area", "column", "ribbon"),
                               x_col = "pos", y_col = "fraction",
                               facet_col = NULL, series_col = group_col,
                               colors = NULL, fill_colors = NULL, smooth_k = 1L,
                               layers = NULL, highlight = NULL, markers = NULL,
                               region = NULL, coordinates = c("genomic", "relative"),
                               labels = list(), scales = list(), legend = list(),
                               theme = NULL, mapping = NULL) {
  style <- match.arg(style)
  coordinates <- match.arg(coordinates)
  d <- as.data.frame(profile)
  if (coordinates == "relative") {
    stopifnot(!is.null(region))
    d[[x_col]] <- plot_positions(d[[x_col]], region)
  }
  if (smooth_k > 1L) {
    # Keep input feature order: reversing promoter coordinates must not shift
    # the original rolling-mean bins or alter the treatment of edge NAs.
    groups <- unique(c(facet_col, series_col))
    indices <- if (!length(groups)) list(seq_len(nrow(d))) else
      split(seq_len(nrow(d)), interaction(d[groups], drop = TRUE))
    for (i in indices) d[[y_col]][i] <- as.numeric(stats::filter(
      d[[y_col]][i], rep(1 / smooth_k, smooth_k), sides = 2))
    d <- d[!is.na(d[[y_col]]), , drop = FALSE]
  }
  if (is.null(mapping)) {
    mapping <- ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])
    if (!is.null(group_col)) mapping$colour <- ggplot2::aes(colour = .data[[group_col]])$colour
    if (!is.null(group_col) && style != "line")
      mapping$fill <- ggplot2::aes(fill = .data[[group_col]])$fill
    if (!is.null(series_col)) mapping$group <- ggplot2::aes(group = .data[[series_col]])$group
  }
  p <- ggplot2::ggplot(d, mapping)
  if (is.null(layers)) {
    layers <- switch(style,
      line = list(list(geom = "line")),
      column = list(list(geom = "col")),
      area = list(list(geom = "area"), list(geom = "line")),
      ribbon = list(list(geom = "ribbon", mapping = ggplot2::aes(ymin = 0,
        ymax = .data[[y_col]])), list(geom = "line")))
  }
  if (!is.null(highlight)) {
    band <- do.call(ggplot2::annotate, c(list(geom = "rect", ymin = -Inf, ymax = Inf),
      highlight[setdiff(names(highlight), "position")]))
    if (identical(highlight$position, "under")) p <- p + band
  }
  for (spec in layers) {
    geom <- match.arg(spec$geom, c("line", "area", "col", "ribbon", "vline"))
    spec$geom <- NULL
    p <- p + do.call(getExportedValue("ggplot2", paste0("geom_", geom)), spec)
  }
  if (!is.null(highlight) && !identical(highlight$position, "under")) p <- p + band
  if (!is.null(markers)) p <- p + do.call(ggplot2::geom_vline, markers)
  if (!is.null(colors)) p <- p + do.call(ggplot2::scale_colour_manual,
    c(list(values = colors), legend))
  if (!is.null(fill_colors)) p <- p + do.call(ggplot2::scale_fill_manual,
    c(list(values = fill_colors), legend))
  if (!is.null(facet_col)) p <- p + ggplot2::facet_wrap(
    stats::reformulate(facet_col), ncol = 1)
  p <- p + scales
  if (!is.null(window)) p <- p + ggplot2::coord_cartesian(xlim = window)
  if (length(labels)) p <- p + do.call(ggplot2::labs, labels)
  if (!is.null(theme)) p <- p + theme
  p
}

# Shared read-by-position heatmap. Layouts preserve the original matrix encoding,
# row ordering, annotation styles and ComplexHeatmap defaults. Callers supply
# FIRE window labels/ticks so this definition has no notebook-specific globals.
# By default return the Heatmap object; draw=TRUE also applies FIRE boundaries.
plot_read_heatmap <- function(res, region = res$region, sample_colors = NULL,
                              include_haplotype = FALSE, variants = NULL,
                              show_cluster_profiles = FALSE, split_alleles = FALSE,
                              cluster_label = "Saved m6A-defined clusters",
                              sample_label_column = NULL,
                              layout = c("genomic", "features", "fire"),
                              main = NULL, timepoint_cols = LEIDEN_TIMEPOINT_COLORS,
                              palette = cluster_palette, rows = NULL, window = NULL,
                              cluster_colors = NULL, ticks = NULL,
                              use_raster = NULL, heatmap_options = list(),
                              draw = FALSE, draw_options = list()) {
  layout <- match.arg(layout)
  boundary_positions <- NULL
  if (layout == "features") {
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
    if (plot_anchor(region)$direction == -1L)
      mat <- mat[, rev(seq_len(ncol(mat))), drop = FALSE]
    mat <- matrix(ifelse(is.na(mat), NA, ifelse(mat > 0, "m6A", "no m6A")),
                  nrow(mat), ncol(mat), dimnames = dimnames(mat))
  
    # create row annotations
    ha <- ComplexHeatmap::rowAnnotation(
      cluster   = df$cluster,
      timepoint = df$sample_name,
      col = list(cluster   = palette(levels(df$cluster)),
                 timepoint = timepoint_cols[levels(df$sample_name)]))
  
    # rows are reads, columns are m6a sites
    heatmap_args <- list(
      matrix = mat,
      name = "m6A call",
      col  = c("m6A" = "black", "no m6A" = "white"),
      na_col = "grey85",
      show_row_names = FALSE, show_column_names = FALSE,
      cluster_rows = FALSE, cluster_columns = FALSE,
      row_split = df$cluster, row_gap = grid::unit(0.6, "mm"),
      row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 8),
      width = grid::unit(11, "cm"), height = grid::unit(14, "cm"),
      use_raster = TRUE,
      column_title = main,
      column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
      left_annotation = ha)
  } else if (layout == "fire") {
    lr <- res
    ex <- window
    cl_cols <- cluster_colors
    time_cols <- timepoint_cols
    pos  <- seq.int(ex$window_start, ex$window_end)
    site <- as.matrix(lr$site_met_mat[rows$RID, , drop = FALSE])
    m <- matrix(0L, nrow(rows), length(pos), dimnames = list(rows$RID, pos))
    m[, match(colnames(site), pos)] <- site
    n_cl  <- table(rows$cluster)
    in_region <- ifelse(pos > ex$start & pos <= ex$end, "tested FIRE region", "outside")
  
    heatmap_args <- list(
      matrix = m, name = "m6A", col = c("0" = "white", "1" = "black"),
      heatmap_legend_param = list(at = c(0, 1), labels = c("no m6A call", "m6A call"),
                                  title = "m6A"),
      cluster_rows = FALSE, cluster_columns = FALSE,
      row_split = rows$cluster,
      row_title = sprintf("%s\nn = %d", names(n_cl), as.integer(n_cl)),
      row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 8),
      row_gap = grid::unit(1, "mm"), border = TRUE,
      show_row_names = FALSE, show_column_names = FALSE,
      use_raster = TRUE, raster_quality = 4,
      left_annotation = ComplexHeatmap::rowAnnotation(
        cluster = rows$cluster, time = rows$time,
        col = list(cluster = cl_cols, time = time_cols),
        annotation_name_gp = grid::gpar(fontsize = 8),
        annotation_legend_param = list(cluster = list(title = "Leiden cluster"),
                                       time = list(title = "Time"))),
      bottom_annotation = ComplexHeatmap::HeatmapAnnotation(
        region = ComplexHeatmap::anno_simple(in_region, col = c("tested FIRE region" = "black", outside = "white"),
                             height = grid::unit(1.5, "mm")),
        position = ComplexHeatmap::anno_mark(at = match(ticks, pos), labels = scales::comma(ticks),
                             side = "bottom", labels_rot = 0,
                             labels_gp = grid::gpar(fontsize = 8)),
        show_annotation_name = c(region = TRUE, position = FALSE),
        annotation_name_gp = grid::gpar(fontsize = 8), annotation_label = c("tested region", "")),
      column_title = paste0(main, "\nLeiden + Manhattan clusters, reads pooled 0-15 min (",
                            ex$chrom, ", ", nrow(rows), " molecules)"),
      column_title_gp = grid::gpar(fontsize = 9))
  
    boundary_positions <- match(c(ex$start, ex$end), pos) / length(pos)
    boundary_slices <- seq_along(n_cl)
  } else {
    assignments <- res$assignments
    if (split_alleles) stopifnot(include_haplotype, !is.null(assignments$allele_display))
    ordering <- if (split_alleles) order(assignments$allele_display, assignments$cluster, assignments$start, assignments$RID) else
      order(assignments$cluster, assignments$start, assignments$RID)
    assignments <- assignments[ordering, , drop = FALSE]
    sample_labels <- if (is.null(sample_label_column)) sub("_.*$", "", assignments$sample_name) else
      assignments[[sample_label_column]]
    if (!is.null(sample_label_column))
      stopifnot(length(sample_labels) == nrow(assignments), all(sample_labels %in% names(sample_colors)))
    sample <- factor(sample_labels, levels = names(sample_colors))
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
    anchor <- plot_anchor(region)
    if (anchor$direction == -1L) met <- met[, rev(seq_len(ncol(met))), drop = FALSE]
    display <- met
    colors <- c("0" = "white", "1" = "black")
    legend <- list(at = c(0, 1), labels = c("no m6A call", "m6A"))
    legend_name <- "m6A"
    top <- NULL
    if (show_cluster_profiles) {
      cluster_colors <- cluster_id_palette(levels(assignments$cluster))
      profile_annotations <- lapply(levels(assignments$cluster), function(cluster) {
        selected <- assignments$cluster == cluster
        ComplexHeatmap::anno_lines(
          colMeans(met[selected, , drop = FALSE]),
          ylim = c(0, 1), gp = grid::gpar(col = cluster_colors[[cluster]], lwd = 0.7),
          axis_param = list(at = c(0, 0.5, 1), labels = c("0", ".5", "1")),
          height = grid::unit(12, "mm"))
      })
      names(profile_annotations) <- paste0(levels(assignments$cluster), " m6A")
      top <- do.call(ComplexHeatmap::HeatmapAnnotation, c(profile_annotations,
        list(annotation_name_gp = grid::gpar(fontsize = 8), gap = grid::unit(1.5, "mm"))))
    }
    ticks <- unique(round(seq(1, ncol(display), length.out = 5L)))
    coordinates <- as.integer(colnames(display))
    tick_positions <- if (anchor$promoter) plot_positions(coordinates[ticks], region) else coordinates[ticks]
    bottom_parts <- list(coordinate = ComplexHeatmap::anno_mark(at = ticks,
      labels = format(tick_positions, scientific = FALSE, trim = TRUE),
      which = "column", side = "bottom", labels_gp = grid::gpar(fontsize = 8)))
    if (identical(region$region_type, "promoter") && !is.null(region$tss) &&
        !is.na(region$tss) && region$tss >= region$analysis_start && region$tss <= region$analysis_end) {
      bottom_parts$TSS <- ComplexHeatmap::anno_mark(
        at = match(region$tss, coordinates),
        labels = if (anchor$promoter) "TSS: 0" else paste0("TSS: ", region$chr, ":", region$tss),
        which = "column", side = "bottom", labels_gp = grid::gpar(fontsize = 8, col = "#D55E00"))
    }
    snps <- if (!is.null(region$focal_snp)) data.frame(pos = region$focal_pos,
      label = paste0(region$focal_snp, " ", region$ref, ">", region$alt)) else window_snps(variants, region)
    if (nrow(snps)) {
      bottom_parts$SNP <- ComplexHeatmap::anno_mark(
        at = match(snps$pos, coordinates),
        labels = if (anchor$promoter) paste0(snps$label, " (", plot_positions(snps$pos, region), " bp from TSS)") else snps$label,
        which = "column", side = "bottom", labels_gp = grid::gpar(fontsize = 7),
        link_gp = grid::gpar(col = "#984EA3"))
    }
    bottom <- do.call(ComplexHeatmap::HeatmapAnnotation, c(bottom_parts,
      list(annotation_name_gp = grid::gpar(fontsize = 8))))
    heatmap_args <- list(
      matrix = display, name = legend_name, col = colors, heatmap_legend_param = legend,
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
        if (anchor$promoter) anchor$x_label,
        if (include_haplotype && !nrow(snps)) "No phased heterozygous SNP in this window",
        paste0(region$chr, ":", region$analysis_start, "-", region$analysis_end)), collapse = "\n"),
      column_title_gp = grid::gpar(fontsize = 11))
  }
  if (!is.null(use_raster)) heatmap_args$use_raster <- use_raster
  for (name in names(heatmap_options)) heatmap_args[name] <- heatmap_options[name]
  heatmap <- do.call(ComplexHeatmap::Heatmap, heatmap_args)
  if (!draw) return(heatmap)
  if (layout == "fire") draw_options <- utils::modifyList(list(
    merge_legend = TRUE, heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"), draw_options)
  do.call(ComplexHeatmap::draw, c(list(object = heatmap), draw_options))
  if (!is.null(boundary_positions)) {
    for (slice in boundary_slices) ComplexHeatmap::decorate_heatmap_body("m6A", slice = slice, {
      grid::grid.segments(x0 = boundary_positions, x1 = boundary_positions, y0 = 0, y1 = 1,
        gp = grid::gpar(lty = 2, lwd = 0.8, col = "#D55E00"))
    })
  }
  invisible(NULL)
}
