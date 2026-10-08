source("/project/spott/cshan/fiber-seq/code/parsing_functions/plotting_functions.r", local = TRUE)

# Plot helpers extracted from leiden_LCL.Rmd; no notebook evaluation at runtime.
# Shared primitives are loaded by common.R::load_shared().


plot_knn_allele_connectivity <- function(summary, main = NULL) {
  connectivity <- summary$connectivity
  connected <- connectivity[connectivity$n_reads_with_neighbors > 0L, , drop = FALSE]
  if (!nrow(connected)) stop("No reads have neighbors in the saved KNN graph")
  bars <- dplyr::bind_rows(
    data.frame(allele = connected$allele, relationship = "Same allele", fraction = connected$same_fraction),
    data.frame(allele = connected$allele, relationship = "Opposite allele", fraction = connected$opposite_fraction))
  bars$relationship <- factor(bars$relationship, levels = c("Same allele", "Opposite allele"))
  counts <- setNames(connectivity$n_reads_with_neighbors, as.character(connectivity$allele))
  ggplot2::ggplot(bars, ggplot2::aes(allele, fraction, fill = relationship)) +
    ggplot2::geom_col(width = 0.65, position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(fraction >= 0.05, scales::percent(fraction, accuracy = 1), "")),
      position = ggplot2::position_stack(vjust = 0.5, reverse = TRUE), color = "white", size = 3) +
    ggplot2::geom_point(data = connected, ggplot2::aes(x = allele, y = expected_same_fraction),
      inherit.aes = FALSE, shape = 23, size = 3, fill = "white", color = "black") +
    ggplot2::scale_fill_manual(values = c("Same allele" = "#0072B2", "Opposite allele" = "#D55E00"), drop = FALSE) +
    ggplot2::scale_x_discrete(labels = function(labels) paste0(labels, "\n(n=", counts[labels], ")")) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(title = main, x = "Source read's focal SNP allele", y = "Mean fraction of neighbors", fill = "Neighbor relationship",
      subtitle = "Unweighted neighbor fractions, averaged equally across reads within each allele",
      caption = paste0("Uses all undirected graph neighbors, including connections between Leiden clusters.\n",
        "Diamonds: region-wide random-mixing baseline excluding the source read, not a significance test.\n",
        sum(connectivity$n_isolated_reads), " isolated reads excluded from connectivity means; retained in composition bars.")) +
    ggplot2::theme_bw(base_size = 11)
}

plot_panel <- function(plot, width = 10, height = 7, draw = NULL, status = "written") {
  list(plot = plot, width = width, height = height, draw = draw, status = status)
}

plot_message_panel <- function(message, title = NULL) {
  p <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = message, size = 4) +
    ggplot2::labs(title = title) + ggplot2::theme_void()
  plot_panel(p, width = 10, height = 4, status = "not_applicable")
}

plot_table <- function(context, data, filename) {
  dir.create(context$table_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(data, file.path(context$table_dir, filename), sep = "\t", quote = FALSE,
                    row.names = FALSE, fileEncoding = "UTF-8")
}

workflow_plot_context <- function(result, assembled, cfg, ds, table_dir) {
  context <- new.env(parent = emptyenv())
  context$result <- result
  context$cfg <- cfg
  context$ds <- ds
  context$table_dir <- table_dir
  context$footprints <- assembled$footprints
  context$track_spec <- assembled$track_spec
  context$track_spec$tracks <- setdiff(context$track_spec$tracks, "m6A")
  if (is.null(context$footprints) || is.null(context$track_spec))
    stop("Assembly must contain footprints and an explicit track_spec")
  a <- result$assignments
  samples <- assembled$sample_table
  if (is.null(samples)) stop("Assembly must contain the full dataset sample_table")
  if (!"sample_label" %in% names(a)) a$sample_label <- samples$sample_label[match(a$sample_name, samples$sample_name)]
  stopifnot(!anyNA(a$sample_label), !anyNA(a$cluster), !anyDuplicated(a$RID))
  context$annotations <- unlist(ds$annotate_by, use.names = FALSE)
  unknown <- setdiff(context$annotations, c("sample", "timepoint", "allele"))
  if (length(unknown)) stop("Unsupported plot annotations: ", paste(unknown, collapse = ", "))
  context$include_allele <- "allele" %in% context$annotations
  if ("timepoint" %in% context$annotations) {
    stopifnot("timepoint" %in% names(samples))
    colors <- unname(timepoint_palette(samples$timepoint)[as.character(samples$timepoint)])
    context$sample_colors <- setNames(colors, samples$sample_label)
    context$timepoint_colors <- setNames(colors, samples$sample_name)
    context$sample_heading <- "Timepoint"
  } else {
    context$sample_colors <- suppressWarnings(sample_palette(samples$sample_name, samples$sample_label))
    context$timepoint_colors <- setNames(unname(context$sample_colors[samples$sample_label]), samples$sample_name)
    context$sample_heading <- "LCL sample"
  }
  context$sample_colors <- context$sample_colors[names(context$sample_colors) %in% a$sample_label]
  context$timepoint_colors <- context$timepoint_colors[names(context$timepoint_colors) %in% a$sample_name]
  result$assignments <- a
  if (context$include_allele) result <- fiberseq_display_result(result)
  context$result <- result
  context$cluster_colors <- cluster_id_palette(levels(result$assignments$cluster))
  context$allele_colors <- if (context$include_allele)
    fiberseq_category_palette(result$allele_display_levels) else character()
  context
}

context_signals <- function(context) {
  if (!exists("signals", envir = context, inherits = FALSE)) {
    context$signals <- signal_profiles(context$result, context$footprints, context$track_spec$tracks)
    plot_table(context, context$signals$profiles, "signal_profiles.tsv")
  }
  context$signals
}

context_allele_summary <- function(context) {
  if (!context$include_allele) stop("This panel requires annotate_by: allele")
  if (!exists("allele_summary", envir = context, inherits = FALSE)) {
    context$allele_summary <- knn_allele_summary(context$result)
    plot_table(context, context$allele_summary$composition, "knn_allele_composition.tsv")
    plot_table(context, context$allele_summary$connectivity, "knn_allele_connectivity.tsv")
    plot_table(context, context$allele_summary$reads, "knn_read_allele_connectivity.tsv")
  }
  context$allele_summary
}

context_embedding <- function(context) {
  if (!exists("embedding", envir = context, inherits = FALSE)) {
    if (nrow(context$result$assignments) < 3L) return(NULL)
    u <- context$cfg$umap
    context$embedding <- fiberseq_umap(context$result, seed = u$seed,
                                      n_neighbors = u$n_neighbors, min_dist = u$min_dist)
    plot_table(context, context$embedding, "umap_coordinates.tsv")
  }
  context$embedding
}

context_fiberseq_profiles <- function(context) {
  if (!context$include_allele) stop("fiberseq_profiles_composition requires allele annotations")
  if (!exists("fiberseq_profiles", envir = context, inherits = FALSE)) {
    tables <- fiberseq_tables(context$result, context$footprints,
                              nucleosome_track = context$track_spec$nucleosome_track,
                              sample_column = "sample_label")
    plot_table(context, tables$profiles, "fiberseq_profiles.tsv")
    plot_table(context, tables$sample, "fiberseq_sample_composition.tsv")
    plot_table(context, tables$allele, "fiberseq_allele_composition.tsv")
    plot_table(context, tables$counts, "fiberseq_cluster_counts.tsv")
    nuc_label <- unname(context$track_spec$labels[context$track_spec$nucleosome_track])
    context$fiberseq_profiles <- plot_fiberseq_profiles(tables, context$result$region,
      context$result$markers, context$cluster_colors, context$sample_colors, context$allele_colors,
      nucleosome_track = context$track_spec$nucleosome_track, nucleosome_label = nuc_label,
      nucleosome_color = unname(context$track_spec$colors[context$track_spec$nucleosome_track]),
      nucleosome_fill = unname(context$track_spec$colors[context$track_spec$nucleosome_track]),
      sample_labels = names(context$sample_colors))
  }
  context$fiberseq_profiles
}

context_umap_panel <- function(context) {
  embedding <- context_embedding(context)
  if (is.null(embedding)) return(plot_message_panel("UMAP requires at least three retained reads.",
                                                  context$result$region$annotation))
  pc <- plot_umap(embedding, "cluster", context$cluster_colors, "UMAP: cluster")
  ps <- plot_umap(embedding, "sample_label", context$sample_colors,
                  paste0("UMAP: ", context$sample_heading))
  pie <- plot_cluster_pie(context$result)
  if (!context$include_allele) {
    # Promoter notebook layout: sample and cluster UMAPs alongside the pie.
    return(plot_panel(cowplot::plot_grid(
      ps + ggplot2::theme(legend.position = "right"),
      pc + ggplot2::theme(legend.position = "right"), pie,
      ncol = 3, rel_widths = c(1.3, 1, 1)), width = 20, height = 7))
  }
  pa <- plot_umap(embedding, "allele_display", context$allele_colors, "UMAP: Focal SNP allele")
  counts <- table(context$result$assignments$cluster)
  labels <- paste0(names(counts), " (n=", as.integer(counts), ")")
  legends <- list(fiberseq_legend(context$cluster_colors, "Cluster", 5, labels),
    fiberseq_legend(context$sample_colors, context$sample_heading, 6),
    fiberseq_legend(context$allele_colors, "Focal SNP allele", 3))
  heights <- vapply(legends, fiberseq_legend_height, numeric(1))
  heading <- cowplot::ggdraw() + cowplot::draw_label(context$result$region$annotation, size = 12)
  p <- cowplot::plot_grid(plotlist = c(list(heading, cowplot::plot_grid(pc, ps, pa, pie, ncol = 2)), legends),
                          ncol = 1, rel_heights = c(.45, 8, heights))
  plot_panel(p, width = 16, height = 8.45 + sum(heights))
}

context_haplotype_panel <- function(context) {
  overview <- context_umap_panel(context)
  profiles <- context_fiberseq_profiles(context)
  plot_panel(cowplot::plot_grid(overview$plot, profiles$plot, ncol = 1,
    rel_heights = c(overview$height, profiles$height)), width = 16,
    height = overview$height + profiles$height)
}

context_genomic_heatmap <- function(context) {
  r <- context$result
  heatmap <- plot_read_heatmap(r, r$region, context$sample_colors,
    include_haplotype = context$include_allele, variants = r$variants,
    show_cluster_profiles = TRUE, sample_label_column = "sample_label")
  draw <- function(x) {
    ComplexHeatmap::draw(x, newpage = FALSE)
    if (plot_anchor(r$region)$promoter) {
      # Exact promoter notebook decoration, including reverse-strand columns.
      tss_x <- (match(as.character(r$region$tss), colnames(x@matrix)) - 0.5) / r$region$width
      for (slice in seq_along(levels(r$assignments$cluster))) {
        ComplexHeatmap::decorate_heatmap_body("m6A", {
          grid::grid.lines(x = c(tss_x, tss_x), y = c(0, 1),
            gp = grid::gpar(col = "#D55E00", lty = "dashed", lwd = 1.2))
        }, slice = slice)
      }
      for (cluster in levels(r$assignments$cluster)) {
        ComplexHeatmap::decorate_annotation(paste0(cluster, " m6A"), {
          grid::grid.lines(x = c(tss_x, tss_x), y = c(0, 1),
            gp = grid::gpar(col = "#D55E00", lty = "dashed", lwd = 1.2))
        })
      }
    }
  }
  plot_panel(heatmap, width = 14, height = 11 + .55 * r$n_clusters, draw = draw)
}

context_read_footprints <- function(context) {
  r <- context$result
  p <- plot_smf_reads(r, context_signals(context)$records, context$sample_colors,
    tracks = context$track_spec$tracks, include_haplotype = context$include_allele,
    track_colors = context$track_spec$colors, track_labels = context$track_spec$labels,
    sample_label_column = "sample_label")
  if (plot_anchor(r$region)$promoter) p <- p + ggplot2::labs(caption = paste0(
    "Dashed line: canonical TSS ", r$region$chr, ":", r$region$tss, " (0 bp)"))
  plot_panel(p, width = 12, height = max(7, .03 * nrow(r$assignments) + 3.5 + .3 * r$n_clusters))
}

context_occupancy <- function(context) {
  r <- context$result
  p <- local({
    .profile_args <- list(track_colors = context$track_spec$colors,
      track_labels = context$track_spec$labels,
      nucleosome_track = context$track_spec$nucleosome_track,
      profiles = context_signals(context)$profiles,
      region = r$region,
      clusters = levels(r$assignments$cluster))
    profiles <- .profile_args$profiles
    region <- .profile_args$region
    clusters <- .profile_args$clusters
    title <- region$annotation
    track_colors <- .profile_args$track_colors
    track_labels <- .profile_args$track_labels
    nucleosome_track <- .profile_args$nucleosome_track
    profiles <- signal_profile_inputs(profiles, region, clusters, track_colors)
    anchor <- plot_anchor(region)
    colors <- feature_colors(levels(profiles$track), track_colors)
    counts <- unique(profiles[, c("cluster", "n_reads")])
    labels <- setNames(paste0(counts$cluster, " (n=", counts$n_reads, ")"), counts$cluster)
    plot <- plot_group_profile(profiles, group_col = "track", style = "ribbon",
      x_col = "relative_pos",
      mapping = ggplot2::aes(relative_pos, fraction, fill = track, color = track, group = track),
      layers = list(list(geom = "ribbon",
        data = profiles[profiles$track %in% nucleosome_track, , drop = FALSE],
        mapping = ggplot2::aes(ymin = 0, ymax = fraction), fill = "grey60", color = NA, alpha = .18),
        list(geom = "line", linewidth = 0.5)))
    plot + ggplot2::facet_wrap(~cluster, ncol = 1, labeller = ggplot2::as_labeller(labels)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
      ggplot2::scale_fill_manual(values = colors, labels = track_labels) +
      ggplot2::scale_color_manual(values = colors, labels = track_labels) +
      ggplot2::scale_x_continuous(limits = c(anchor$left - 0.5, anchor$right + 0.5), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
        expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::labs(x = anchor$x_label, y = "Fraction of cluster reads", title = title,
        subtitle = "m6A calls and footprint occupancy; retained heterozygous reads", fill = "Plot track", color = "Plot track") +
      cowplot::theme_cowplot(font_size = 10) + cowplot::panel_border() +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4), color = ggplot2::guide_legend(ncol = 4)) +
      ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8),
        plot.margin = ggplot2::margin(5.5, 16, 5.5, 5.5))
  }) +
    ggplot2::labs(subtitle = if (context$include_allele)
      "m6A calls and footprint occupancy; retained heterozygous reads" else
      "m6A calls and footprint occupancy; all full-span reads")
  if (plot_anchor(r$region)$promoter) p <- p + ggplot2::labs(caption = paste0(
    "Dashed line: canonical TSS ", r$region$chr, ":", r$region$tss, " (0 bp)"))
  plot_panel(p, width = 10, height = 1.4 * r$n_clusters + 2.5)
}

context_pairwise_delta <- function(context) {
  smooth <- context$ds$plot_options$delta_smooth_bp
  if (is.null(smooth)) smooth <- 25L
  delta <- m6a_delta_table(context$result, smooth_bp = smooth)
  if (is.null(delta)) {
    # Preserve a valid declared PDF when Leiden returns only one cluster.
    return(plot_message_panel("Pairwise m6A differences require at least two clusters.",
                              context$result$region$annotation))
  }
  plot_table(context, delta, "pairwise_m6a_delta.tsv")
  p <- plot_m6a_delta(delta, context$result$region, smooth_bp = smooth)
  plot_panel(p$plot, width = 10, height = p$height)
}

# Macrophage custom windows use their centre as the anchor. Only the display
# coordinates change; the source matrices and saved clustering remain untouched.
context_profile_result <- function(context) {
  r <- context$result
  anchor <- plot_anchor(r$region)
  if (!anchor$promoter && (is.null(r$region$focal_pos) || is.na(r$region$focal_pos))) {
    r$region$analysis_start <- r$region$analysis_start - anchor$anchor
    r$region$analysis_end <- r$region$analysis_end - anchor$anchor
    colnames(r$site_met_mat) <- as.character(as.numeric(colnames(r$site_met_mat)) - anchor$anchor)
    colnames(r$profiles) <- as.character(as.numeric(colnames(r$profiles)) - anchor$anchor)
  }
  r
}

