# Figures of the pooled Pol II clustering, built with the shared Leiden-Manhattan
# helpers (plot_read_heatmap(), plot_group_profile(), plot_stacked_proportion(),
# stack_rows(), plot_molecules(), cluster_id_palette(), sample_palette(),
# save_figure()). Load them first with pol2_load_helpers(). Clusters are drawn in
# the summary's display order (PI order).
# plot_molecules() calls ggplot2 functions unqualified.
suppressPackageStartupMessages(library(ggplot2))

POL2_X_LABEL <- "Position relative to FANTOM5 dominant TSS (bp); upstream < 0"

# TSS-relative stand-in region for plot_read_heatmap(): positions are already
# strand-oriented, so tss = 0 on "+" leaves them unchanged.
pol2_plot_region <- function() data.frame(chr = "TSS", analysis_start = -1000L, analysis_end = 999L,
  width = 2000L, region_type = "promoter", strand = "+", tss = 0L, annotation = "")

pol2_sample_colors <- function(cfg, present) {
  colors <- sample_palette(fread(cfg$samples)$sample_name)
  colors[names(colors) %in% sub("_.*$", "", present)]
}

# Rows in heatmap order: cluster, then whole-window m6A proportion (passed as
# `start`, the within-cluster sort key of plot_read_heatmap()).
pol2_heatmap_result <- function(x, annotation, max_reads, seed) {
  m <- annotation$molecules
  stopifnot(identical(m$observation_id, rownames(x$sites)))
  set.seed(seed)
  draw <- if (nrow(m) > max_reads) sort(sample(seq_len(nrow(m)), max_reads)) else seq_len(nrow(m))
  a <- data.frame(RID = m$observation_id[draw], cluster = droplevels(m$cluster[draw]),
                  start = m$m6a_proportion[draw], sample_name = m$sample_name[draw])
  a <- a[order(a$cluster, a$start, a$RID), , drop = FALSE]
  list(assignments = a, site_met_mat = x$sites)
}

# Draw a heatmap and mark the TSS (column 1001 of -1000..999) on every slice.
pol2_draw_heatmap <- function(heatmap, clusters, profiles = FALSE) {
  ComplexHeatmap::draw(heatmap, newpage = FALSE)
  tss_x <- (1001 - 0.5) / 2000
  tss_line <- function() grid::grid.lines(x = c(tss_x, tss_x), y = c(0, 1),
    gp = grid::gpar(col = "#D55E00", lty = "dashed", lwd = 1.2))
  for (slice in seq_along(clusters)) ComplexHeatmap::decorate_heatmap_body("m6A", tss_line(), slice = slice)
  if (profiles) for (cl in clusters) ComplexHeatmap::decorate_annotation(paste0(cl, " m6A"), tss_line())
  invisible(NULL)
}

# Leiden-Manhattan read heatmap: per-bp m6A calls, cluster and sample bars,
# per-cluster call fraction lines on top.
pol2_molecule_heatmap <- function(res, sample_colors) {
  plot_read_heatmap(res, pol2_plot_region(), sample_colors, show_cluster_profiles = TRUE,
    heatmap_options = list(column_title = paste0(
      "Pol II clustering: ", nrow(res$assignments), " footprint-positive molecules, 31 LCL samples pooled\n",
      "Top: m6A call fraction per cluster\n", POL2_X_LABEL), column_title_gp = grid::gpar(fontsize = 11)))
}

# Leiden-Manhattan cluster accessibility: 50-bp bin profiles, one panel per cluster.
pol2_plot_profiles <- function(annotation, pi_labels = FALSE) {
  s <- annotation$summary
  d <- as.data.frame(annotation$profiles)
  d$cluster <- factor(d$cluster, levels = s$cluster)
  colors <- cluster_id_palette(s$cluster)
  labels <- setNames(paste0(s$cluster, " (n=", s$n, ")"), s$cluster)
  p <- plot_group_profile(d, group_col = "cluster", style = "ribbon", x_col = "relative_bp", y_col = "m6a_proportion",
      layers = list(list(geom = "ribbon", mapping = ggplot2::aes(ymin = 0, ymax = m6a_proportion), alpha = 0.35, colour = NA),
                    list(geom = "line", linewidth = 0.5)),
      colors = colors, fill_colors = colors) +
    ggplot2::facet_wrap(~cluster, ncol = 1, labeller = ggplot2::as_labeller(labels)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    ggplot2::scale_x_continuous(limits = c(-1000, 1000), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = POL2_X_LABEL, y = "m6A proportion at observed sites (50 bp bins)",
                  title = "Cluster-specific accessibility") +
    cowplot::theme_cowplot(font_size = 10) + cowplot::panel_border() +
    ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(5.5, 16, 5.5, 5.5))
  if (pi_labels) {
    text <- data.frame(cluster = factor(s$cluster, levels = s$cluster),
      label = ifelse(is.na(s$median_pi), "PI: NA", sprintf("PI: %.1f (%d TSSs)", s$median_pi, s$n_tss_with_pi)))
    p <- p + ggplot2::geom_text(data = text, ggplot2::aes(x = 980, y = 0.9, label = label),
                                inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.8)
  }
  p
}

# Pol II ChIP fold-enrichment at the TSS (-50..+299), log2 of cluster mean over
# the all-molecule mean, with BH stars.
pol2_plot_signal_enrichment <- function(summary) {
  d <- as.data.frame(summary)
  d <- d[is.finite(d$log2_ratio), , drop = FALSE]
  # coord_flip draws the first level at the bottom: reverse to match the heatmap.
  d$cluster <- factor(d$cluster, levels = rev(summary$cluster))
  reach <- max(abs(d$log2_ratio), 0.1)
  plot_stacked_proportion(d, x = "cluster", fill = "cluster", colors = cluster_id_palette(summary$cluster),
    y = "log2_ratio", position = "dodge", horizontal = TRUE,
    labels = list(x = NULL, y = "log2(cluster mean / all molecules)", title = "Pol II ChIP at TSS (-50..+299)",
                  caption = "Wilcoxon, cluster vs rest; BH: *** < 0.001, ** < 0.01, * < 0.05"),
    totals = list(data = d, x = "cluster", y = 1.2 * reach, label = "signif", style = list(size = 3)),
    scales = list(ggplot2::geom_hline(yintercept = 0, linewidth = 0.3),
                  ggplot2::scale_y_continuous(limits = c(-1.35 * reach, 1.35 * reach))),
    legend = list(guide = "none"), theme = theme_fiberseq("classic", base_size = 10))
}

pol2_plot_footprints <- function(annotation) {
  s <- annotation$summary
  fp <- as.data.frame(annotation$footprints)
  fp$cluster <- factor(fp$cluster, levels = s$cluster)
  labels <- setNames(paste0(s$cluster, " (n=", s$n, ")"), s$cluster)
  ggplot2::ggplot(fp, ggplot2::aes(midpoint_relative, fill = cluster)) +
    ggplot2::geom_histogram(binwidth = 10, boundary = 0) +
    ggplot2::facet_wrap(~cluster, ncol = 1, scales = "free_y", labeller = ggplot2::as_labeller(labels)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = cluster_id_palette(s$cluster)) +
    # Footprints overlapping a window edge can have midpoints just outside it.
    ggplot2::coord_cartesian(xlim = c(-1000, 1000)) +
    ggplot2::labs(x = "30-60 bp footprint midpoint relative to TSS (bp)", y = "Footprints",
                  title = "Candidate footprint positions") +
    cowplot::theme_cowplot(font_size = 10) + cowplot::panel_border() +
    ggplot2::theme(legend.position = "none")
}

# 30-60 bp footprints of every heatmap molecule, one row per molecule in the
# heatmap's row order, drawn with plot_molecules(). midpoint_relative is already
# strand-oriented, so each footprint spans midpoint +/- length / 2.
pol2_plot_footprint_molecules <- function(res, annotation) {
  rows <- stack_rows(data.frame(RID = res$assignments$RID, cluster = res$assignments$cluster), "cluster")
  fp <- annotation$footprints[observation_id %in% rows$RID]
  fp_class <- "TF footprint 30-60 bp"
  tracks <- data.frame(RID = fp$observation_id,
                       start = fp$midpoint_relative - fp$footprint_length / 2,
                       end = fp$midpoint_relative + fp$footprint_length / 2,
                       class = factor(fp_class, levels = fp_class))
  n <- table(rows$cluster)
  # start = end = 0 puts plot_molecules()' dashed region lines on the TSS.
  ex <- data.frame(chrom = "TSS", start = 0, end = 0, window_start = -1000, window_end = 1000)
  plot_molecules(rows, tracks, ex, facet_col = "cluster",
    facet_labels = setNames(paste0(names(n), "\n(n=", n, ")"), names(n)),
    bar_cols = "cluster", bar_palette = cluster_id_palette(levels(rows$cluster)),
    track_colors = MOLECULE_TRACK_COLORS[fp_class],
    title = paste0("30-60 bp FiberHMM footprints of ", nrow(rows), " clustered molecules"),
    subtitle = "Rows in molecule_heatmap.pdf order; clusters by median PI, lowest first") +
    ggplot2::labs(x = POL2_X_LABEL)
}

# Fig. 1c wedge beside the profiles: clusters run from lowest median PI (top)
# to highest (bottom), so the wedge widens downward.
pol2_plot_pi_wedge <- function() {
  ggplot2::ggplot(data.frame(x = c(0, 0, 1), y = c(1, 0, 0)), ggplot2::aes(x, y)) +
    ggplot2::geom_polygon(fill = "grey70") +
    ggplot2::annotate("text", x = 0.3, y = 0.4, label = "Pausing index", angle = -90, size = 3.2) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::theme_void()
}

# Fig. 1c layout: heatmap with cluster bar, per-cluster profiles with PI, a PI
# wedge, and Pol II ChIP signal at the TSS, clusters in the same (PI) order.
pol2_fig1c <- function(res, annotation, sample_colors) {
  a <- res$assignments
  clusters <- levels(a$cluster)
  heatmap <- plot_read_heatmap(res, pol2_plot_region(), sample_colors,
    heatmap_options = list(
      left_annotation = ComplexHeatmap::rowAnnotation(cluster = a$cluster,
        col = list(cluster = cluster_id_palette(clusters)), show_legend = FALSE),
      column_title = paste0("Single-molecule m6A (", nrow(a), " molecules)"),
      column_title_gp = grid::gpar(fontsize = 10)))
  height <- max(7, 0.9 * length(clusters) + 2)
  grob <- grid::grid.grabExpr(pol2_draw_heatmap(heatmap, clusters), width = 7, height = height)
  profiles <- cowplot::plot_grid(pol2_plot_profiles(annotation, pi_labels = TRUE), pol2_plot_pi_wedge(),
                                 nrow = 1, align = "h", axis = "tb", rel_widths = c(1, 0.06))
  cowplot::plot_grid(grob, profiles, pol2_plot_signal_enrichment(annotation$summary),
                     nrow = 1, rel_widths = c(1.3, 1.06, 0.7))
}

pol2_save_figures <- function(x, annotation, cfg) {
  res <- pol2_heatmap_result(x, annotation, cfg$heatmap_max_reads, cfg$seed)
  colors <- pol2_sample_colors(cfg, res$assignments$sample_name)
  clusters <- levels(res$assignments$cluster)
  k <- length(clusters)
  plots <- file.path(cfg$output, "plots")
  save_figure(pol2_molecule_heatmap(res, colors), file.path(plots, "molecule_heatmap.pdf"),
              width = 14, height = 11 + 0.55 * k, draw = function(h) pol2_draw_heatmap(h, clusters, profiles = TRUE))
  save_figure(pol2_plot_profiles(annotation), file.path(plots, "cluster_profiles.pdf"),
              width = 10, height = 1.4 * k + 2.5)
  save_figure(pol2_plot_footprints(annotation), file.path(plots, "footprint_positions.pdf"),
              width = 10, height = 1.4 * k + 2.5)
  save_figure(pol2_plot_footprint_molecules(res, annotation), file.path(plots, "footprint_molecules.pdf"),
              width = 9, height = 11 + 0.55 * k)
  save_figure(pol2_fig1c(res, annotation, colors), file.path(plots, "fig1c_pol2_clusters.pdf"),
              width = 16, height = max(7, 0.9 * k + 2))
}
