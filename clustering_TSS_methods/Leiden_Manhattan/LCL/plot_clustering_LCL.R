# All LCL plotting and figure exports. Processing lives in the other two scripts.
# Cluster and footprint colors match the TNF non-LCL examples.
LCL_CLUSTER_COLORS <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1",
  "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
  "green1", "yellow4", "yellow3", "darkorange4", "brown")
lcl_cluster_palette <- function(levels) {
  idx <- suppressWarnings(as.integer(sub("^cluster", "", levels)))
  if (anyNA(idx)) idx <- seq_along(levels)
  cols <- if (max(idx) <= length(LCL_CLUSTER_COLORS)) LCL_CLUSTER_COLORS[idx] else
    grDevices::colorRampPalette(LCL_CLUSTER_COLORS)(max(idx))[idx]
  setNames(cols, levels)
}

LCL_HAPLOTYPE_COLORS <- c(HP1 = "#ADD8E6", HP2 = "#FFF2AE", unphased = "#999999", pooled = "#BBBBBB")

# TNF footprint colors and labels; m6A uses purple outside heatmaps.
LCL_TRACK_COLORS <- c(m6A = "#800080", "ft_nuc_130-160bp" = "#4d4d4d",
  "FiberHMM_10-30bp" = "#fdae6b", "FiberHMM_40-60bp" = "#f16913",
  "FiberHMM_60-80bp" = "#a63603")
LCL_TRACK_LABELS <- c(m6A = "m6A", "ft_nuc_130-160bp" = "nucleosome footprint (130-160 bp)",
  "FiberHMM_10-30bp" = "TF footprint 10-30 bp", "FiberHMM_40-60bp" = "TF footprint 40-60 bp",
  "FiberHMM_60-80bp" = "TF footprint 60-80 bp")
lcl_feature_colors <- function(tracks) {
  stopifnot(all(tracks %in% names(LCL_TRACK_COLORS)))
  LCL_TRACK_COLORS[tracks]
}

plot_lcl_smf_reads <- function(result, records, sample_colors, tracks = LCL_FOOTPRINT_TRACKS,
                                include_haplotype = FALSE) {
  tracks <- LCL_FOOTPRINT_TRACKS
  records <- records[records$track %in% tracks, , drop = FALSE]
  prepared <- prepare_lcl_read_tracks(result, records, sample_colors)
  reads <- prepared$reads
  features <- prepared$features[prepared$features$track %in% tracks, , drop = FALSE]
  anchor <- prepared$anchor
  width <- max(1, anchor$right - anchor$left)
  feature_colors <- lcl_feature_colors(tracks)
  # Draw nucleosomes first and TF footprints on top, as in TNF Fig 1.
  features <- features[order(match(features$track, tracks)), ]
  features$ymin <- features$row - 0.45
  features$ymax <- features$row + 0.45
  features$fill <- unname(feature_colors[features$track])
  cluster_colors <- lcl_cluster_palette(levels(reads$cluster))
  make_bar <- function(offset, colors) data.frame(cluster = reads$cluster, row = reads$row,
    left = anchor$left - offset * width, right = anchor$left - (offset - 0.025) * width,
    fill = unname(colors))
  bars <- rbind(make_bar(0.075, cluster_colors[as.character(reads$cluster)]),
                 make_bar(0.040, sample_colors[reads$sample_label]))
  if (include_haplotype) bars <- rbind(bars, make_bar(0.110, LCL_HAPLOTYPE_COLORS[reads$haplotype]))
  legend_colors <- c(feature_colors, sample_colors,
    if (include_haplotype) LCL_HAPLOTYPE_COLORS[unique(reads$haplotype)])
  legend_labels <- c(unname(LCL_TRACK_LABELS[tracks]), names(sample_colors),
    if (include_haplotype) unique(reads$haplotype))
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
      subtitle = "One row per molecule; pooled m6A-defined Leiden clusters") +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4, byrow = TRUE)) +
    cowplot::theme_cowplot(font_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(angle = 0),
      panel.spacing.y = grid::unit(2, "mm"), legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 7), legend.key.size = grid::unit(3, "mm"))
}

# TNF Fig 2 layout: all features overlaid in one panel for every cluster.
plot_lcl_signal_profile <- function(profiles, region, clusters, type = c("line", "mixed"),
                                     title = region$annotation) {
  type <- match.arg(type)
  profiles <- profiles[profiles$track %in% names(LCL_TRACK_COLORS), , drop = FALSE]
  anchor <- lcl_plot_anchor(region)
  profiles$relative_pos <- anchor$direction * (profiles$pos - anchor$anchor)
  profiles$cluster <- factor(profiles$cluster, levels = clusters)
  profiles$track <- factor(profiles$track, levels = names(LCL_TRACK_COLORS)[names(LCL_TRACK_COLORS) %in% profiles$track])
  profiles <- profiles[order(profiles$cluster, profiles$track, profiles$relative_pos), ]
  colors <- lcl_feature_colors(levels(profiles$track))
  counts <- unique(profiles[, c("cluster", "n_reads")])
  labels <- setNames(paste0(counts$cluster, " (n=", counts$n_reads, ")"), counts$cluster)
  plot <- ggplot2::ggplot(profiles,
    ggplot2::aes(relative_pos, fraction, fill = track, color = track, group = track))
  if (type == "mixed") {
    plot <- plot + ggplot2::geom_col(data = profiles[profiles$track != "m6A", , drop = FALSE],
      width = 1, linewidth = 0, position = "identity", alpha = 0.5) +
      ggplot2::geom_line(data = profiles[profiles$track == "m6A", , drop = FALSE], linewidth = 0.55)
  } else {
    plot <- plot + ggplot2::geom_line(linewidth = 0.5)
  }
  plot + ggplot2::facet_wrap(~cluster, ncol = 1, labeller = ggplot2::as_labeller(labels)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = colors, labels = LCL_TRACK_LABELS) +
    ggplot2::scale_color_manual(values = colors, labels = LCL_TRACK_LABELS) +
    ggplot2::scale_x_continuous(limits = c(anchor$left - 0.5, anchor$right + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = anchor$x_label, y = "Fraction of cluster reads", title = title,
      subtitle = "m6A calls and footprint occupancy; all samples pooled", fill = "Plot track", color = "Plot track") +
    cowplot::theme_cowplot(font_size = 10) + cowplot::panel_border() +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4), color = ggplot2::guide_legend(ncol = 4)) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(5.5, 16, 5.5, 5.5))
}

lcl_sample_palette <- function(sample_names) {
  labels <- sub("_.*$", "", sample_names)
  labels <- unique(labels[order(as.integer(sub("^AL-?([0-9]+).*$", "\\1", labels)))])
  setNames(grDevices::hcl.colors(length(labels), "Dynamic"), labels)
}

plot_lcl_cluster_heatmap <- function(res, region, sample_colors,
                                     include_haplotype = FALSE, variants = NULL,
                                     show_cluster_profiles = FALSE) {
  assignments <- res$assignments
  assignments <- assignments[order(assignments$cluster, assignments$start, assignments$RID), ]
  sample <- factor(sub("_.*$", "", assignments$sample_name), levels = names(sample_colors))
  annotation <- ComplexHeatmap::rowAnnotation(
    cluster = assignments$cluster, sample = sample,
    col = list(cluster = lcl_cluster_palette(levels(assignments$cluster)), sample = sample_colors))
  if (include_haplotype) {
    annotation <- ComplexHeatmap::rowAnnotation(
      cluster = assignments$cluster, sample = sample, haplotype = assignments$haplotype,
      col = list(cluster = lcl_cluster_palette(levels(assignments$cluster)), sample = sample_colors,
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
    cluster_colors <- lcl_cluster_palette(levels(assignments$cluster))
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
  snps <- lcl_window_snps(variants, region)
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
    row_split = assignments$cluster, row_gap = grid::unit(0.6, "mm"),
    row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 8),
    show_row_names = FALSE, show_column_names = FALSE,
    use_raster = TRUE, raster_quality = 2, raster_resize_mat = FALSE,
    left_annotation = annotation,
    top_annotation = top, bottom_annotation = bottom,
    column_title = paste(c(region$annotation,
      if (include_haplotype) "Pooled m6A-defined clusters; haplotype annotated per sample",
      if (show_cluster_profiles) "Top: m6A call fraction per cluster (all full-span reads)",
      if (include_haplotype && !nrow(snps)) "No phased heterozygous SNP in this window",
      paste0(region$chr, ":", region$start, "-", region$end)), collapse = "\n"),
    column_title_gp = grid::gpar(fontsize = 11))
}

plot_lcl_sample_composition <- function(res, sample_colors, main = NULL) {
  assignments <- res$assignments
  assignments$sample <- factor(sub("_.*$", "", assignments$sample_name), levels = names(sample_colors))
  ggplot2::ggplot(assignments, ggplot2::aes(sample, fill = cluster)) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::scale_fill_manual(values = lcl_cluster_palette(levels(assignments$cluster))) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(x = "LCL sample", y = "Fraction of sample reads", title = main) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
}


# One bar per cluster: calls divided by reads x union of observed m6A sites.
plot_lcl_methylation_by_cluster <- function(summary, region) {
  summary$cluster <- factor(summary$cluster, levels = summary$cluster)
  ggplot2::ggplot(summary, ggplot2::aes(cluster, methylation_proportion, fill = cluster)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = lcl_cluster_palette(levels(summary$cluster)), guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.03))) +
    ggplot2::labs(x = "m6A-defined cluster", y = "Mean m6A call proportion", title = region$annotation,
      subtitle = "Mean across observed m6A sites; all cluster reads included") +
    cowplot::theme_cowplot(font_size = 10) + cowplot::panel_border()
}

save_lcl_plots <- function(result, footprints, sample_colors, output_dir, include_haplotype = FALSE) {
  table_dir <- file.path(output_dir, "tables")
  plot_dir <- file.path(output_dir, "plots")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  signals <- lcl_signal_profiles(result, footprints)
  profiles <- signals$profiles
  records <- signals$records
  matched <- match(records$RID, result$assignments$RID)
  for (column in intersect(c("cluster", "sample_name", "haplotype", "phase_set", "group_id"),
                           names(result$assignments))) records[[column]] <- result$assignments[[column]][matched]
  data.table::fwrite(records, file.path(table_dir, "single_molecule_features.tsv.gz"), sep = "\t")
  data.table::fwrite(profiles, file.path(table_dir, "aggregate_signal_profiles.tsv.gz"), sep = "\t")
  met <- cluster_site_profiles(result, result$site_met_mat)
  data.table::fwrite(met, file.path(table_dir, "aggregate_accessibility.tsv"), sep = "\t")
  summary <- lcl_methylation_by_cluster(result)
  data.table::fwrite(summary, file.path(table_dir, "methylation_proportion_by_cluster.tsv"), sep = "\t")
  ggplot2::ggsave(file.path(plot_dir, "methylation_proportion_by_cluster.pdf"),
    plot_lcl_methylation_by_cluster(summary, result$region), width = 8, height = 4.5)
  heatmap_file <- if (include_haplotype) "heatmap_m6a_footprints.pdf" else "heatmap_m6a.pdf"
  heatmap_height <- if (include_haplotype) 11 + 0.55 * result$n_clusters else 11
  heatmap <- plot_lcl_cluster_heatmap(result, result$region, sample_colors,
    include_haplotype = include_haplotype, variants = result$variants,
    show_cluster_profiles = include_haplotype)
  grDevices::pdf(file.path(plot_dir, heatmap_file), width = 14, height = heatmap_height)
  tryCatch(ComplexHeatmap::draw(heatmap, newpage = FALSE), finally = grDevices::dev.off())
  ggplot2::ggsave(file.path(plot_dir, "fig1_read_footprints.pdf"),
    plot_lcl_smf_reads(result, records, sample_colors, include_haplotype = include_haplotype),
    width = 12, height = max(7, 0.03 * nrow(result$assignments) + 3.5 + 0.3 * result$n_clusters), limitsize = FALSE)
  if (include_haplotype) {
    ggplot2::ggsave(file.path(plot_dir, "fig2_occupancy_by_cluster.pdf"),
      plot_lcl_signal_profile(profiles, result$region, levels(result$assignments$cluster), "line"),
      width = 10, height = 1.4 * result$n_clusters + 2.5, limitsize = FALSE)
  }
  aggregate <- plot_lcl_signal_profile(profiles, result$region,
    levels(result$assignments$cluster), if (include_haplotype) "mixed" else "line")
  ggplot2::ggsave(file.path(plot_dir, "aggregate_profile.pdf"), aggregate,
    width = 10, height = 1.4 * result$n_clusters + 2.5, limitsize = FALSE)
  if (!include_haplotype) {
    ggplot2::ggsave(file.path(plot_dir, "cluster_sample_composition.pdf"),
      plot_lcl_sample_composition(result, sample_colors, result$region$annotation), width = 12, height = 5)
  }
  invisible(signals)
}

# One row per region; Pandoc embeds the PNG files when writing self-contained HTML.
# RStudio/knitr may have a shorter PATH than an interactive shell.
lcl_ghostscript <- function() {
  candidates <- unique(c(Sys.getenv("R_GSCMD", unset = ""), unname(Sys.which("gs")),
    "/usr/bin/gs", "/usr/local/bin/gs"))
  candidates <- candidates[nzchar(candidates)]
  available <- candidates[file.exists(candidates) & file.access(candidates, mode = 1L) == 0L]
  if (!length(available)) stop(
    "Ghostscript was not found. Set R_GSCMD to the full path of the Ghostscript executable before knitting.")
  available[[1L]]
}

lcl_report_gallery <- function(result_paths, include_haplotype = FALSE) {
  gs <- lcl_ghostscript()
  files <- c(if (include_haplotype) "heatmap_m6a_footprints.pdf" else "heatmap_m6a.pdf",
    "fig1_read_footprints.pdf", if (include_haplotype) "fig2_occupancy_by_cluster.pdf", "aggregate_profile.pdf",
    "methylation_proportion_by_cluster.pdf",
    if (!include_haplotype) "cluster_sample_composition.pdf")
  labels <- c("m6A heatmap", "Read footprints", if (include_haplotype) "Occupancy by cluster", "Aggregate profile",
    "Methylation by cluster", if (!include_haplotype) "Sample composition")
  rows <- lapply(names(result_paths), function(rid) {
    plot_dir <- file.path(dirname(dirname(result_paths[[rid]])), "plots")
    cells <- lapply(seq_along(files), function(i) {
      pdf <- file.path(plot_dir, files[i])
      # Keep rendering files in the R session's temporary directory, not the outputs.
      png <- tempfile(pattern = "lcl-report-", fileext = ".png")
      stopifnot(file.exists(pdf))

      status <- system2(gs, c("-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m",
        "-r120", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4", "-dFirstPage=1", "-dLastPage=1",
        shQuote(paste0("-sOutputFile=", png)), shQuote(pdf)))
      stopifnot(status == 0L, file.exists(png))
      gateway <- paste0("https://gate.rcc.uchicago.edu/spott/file_show?path=",
        utils::URLencode(normalizePath(pdf), reserved = TRUE))
      htmltools::tags$td(
        htmltools::tags$button(type = "button", class = "lcl-preview", title = "Enlarge preview",
          onclick = "var d=document.getElementById('lcl-zoom'); d.querySelector('img').src=this.querySelector('img').src; d.showModal();",
          htmltools::tags$img(src = normalizePath(png),
            alt = paste(rid, labels[i]), loading = "lazy")),
        htmltools::tags$div(
          htmltools::tags$a(href = gateway, target = "_blank", rel = "noopener", "Open PDF")))
    })
    do.call(htmltools::tags$tr, c(list(htmltools::tags$th(scope = "row", rid)), cells))
  })
  header <- do.call(htmltools::tags$tr, lapply(c("Region", labels), htmltools::tags$th))
  html <- htmltools::tags$div(class = "lcl-gallery", htmltools::tags$table(
    htmltools::tags$thead(header), do.call(htmltools::tags$tbody, rows)))
  knitr::asis_output(paste0("\n```{=html}\n", as.character(html), "\n```\n"))
}
