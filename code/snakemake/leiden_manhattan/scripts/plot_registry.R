source("/project/spott/cshan/fiber-seq/code/parsing_functions/plotting_functions.r", local = TRUE)

# The single panel-name -> PDF filename/function registry for every dataset.
# A renderer returns plot_panel(); expensive shared summaries are cached in context.
panel_entry <- function(filename, render, enabled = NULL) {
  if (is.character(filename)) {
    name <- filename
    filename <- function(region_id) name
  }
  list(filename = filename, render = render, enabled = enabled)
}

panel_registry <- list(
  cluster_met_profiles = panel_entry(function(id) paste0("cluster_met_profiles_", id, ".pdf"), function(c) {
    r <- context_profile_result(c)
    p <- local({
      .profile_args <- list(main = r$region$region_id,
      palette = cluster_id_palette,
      x_label = plot_anchor(c$result$region)$x_label,
      res = r,
      met_mat = r$site_met_mat)
      res <- .profile_args$res
      met_mat <- .profile_args$met_mat
      tss <- NULL
      main <- .profile_args$main
      region <- res$region
      palette <- .profile_args$palette
      x_label <- .profile_args$x_label
      # cluster_site_profiles() takes cluster labels and computes mean m6a value at every 
      # position for each cluster
      prof <- cluster_site_profiles(res, met_mat)
      anchor <- plot_anchor(region)
      if (anchor$promoter) {
        prof$pos <- plot_positions(prof$pos, region)
        tss <- 0
      }
      lv   <- levels(res$assignments$cluster)
      pal  <- palette(lv)
    
      p_list <- lapply(lv, function(cl) {
        d <- prof[prof$cluster == cl, ]
        gg <- plot_group_profile(d, style = "column", y_col = "met",
          mapping = aes(x = pos, y = met),
          layers = list(list(geom = "col", fill = pal[cl]))) +
          ylim(0, 1) +
          ylab("met prop.") + xlab(if (!is.null(x_label)) x_label else if (anchor$promoter) anchor$x_label else "pos") +
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
    })
    plot_panel(p, width = 6, height = 1.4 * r$n_clusters + 1.4)
  }),
  cluster_composition = panel_entry(function(id) paste0("cluster_composition_", id, ".pdf"), function(c) {
    inputs <- cluster_composition_inputs(c$result, c$timepoint_colors)
    p1 <- plot_stacked_proportion(inputs$assignments, "sample_name", "cluster",
      cluster_id_palette(levels(inputs$assignments$cluster)), position = "fill", width = .9,
      labels = list(x = NULL, y = "fraction of reads",
        title = paste0("Cluster composition per timepoint, ", c$result$region$region_id)),
      theme = cowplot::theme_cowplot(font_size = 10))
    p2 <- plot_stacked_proportion(inputs$proportions, "cluster", "sample_name", c$timepoint_colors,
      y = "proportion", position = "dodge", width = .7, legend = list(name = "timepoint"),
      labels = list(x = "cluster", y = "proportion of that timepoint's reads",
        title = paste0("Timepoints per cluster, ", c$result$region$region_id)),
      theme = cowplot::theme_cowplot(font_size = 10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
    plot_panel(cowplot::plot_grid(p1, p2, ncol = 1), width = 7, height = 8)
  }),
  heatmap = panel_entry(function(id) paste0("heatmap_", id, ".pdf"), function(c) {
    plot_panel(plot_read_heatmap(c$result, main = c$result$region$region_id,
      timepoint_cols = c$timepoint_colors, palette = cluster_id_palette, layout = "features"), width = 8, height = 9,
      draw = function(x) ComplexHeatmap::draw(x, newpage = FALSE))
  }),
  met_fraction_lines = panel_entry(function(id) paste0("met_fraction_lines_", id, ".pdf"), function(c) {
    r <- context_profile_result(c)
    smooth <- c$ds$plot_options$line_smooth_k
    if (is.null(smooth)) smooth <- 1L
    p <- local({
      .profile_args <- list(main = r$region$region_id,
      smooth_k = smooth,
      palette = cluster_id_palette,
      x_label = plot_anchor(c$result$region)$x_label,
      res = r)
      res <- .profile_args$res
      tss <- NULL
      main <- .profile_args$main
      smooth_k <- .profile_args$smooth_k
      region <- res$region
      palette <- .profile_args$palette
      x_label <- .profile_args$x_label
      P   <- res$profiles
      pos <- as.numeric(colnames(P))
      anchor <- plot_anchor(region)
      if (anchor$promoter) {
        # Column names are genomic window midpoints; keep the original bins.
        pos <- plot_positions(pos, region)
        tss <- 0
      }
    
      df <- met_fraction_inputs(res, smooth_k, region)
    
      ylab <- if (res$params$window_size == 0) "m6A fraction per site"
              else sprintf("mean m6A per %d-bp window", res$params$window_size)
      sub  <- if (smooth_k > 1) sprintf("rolling mean over %d features", smooth_k) else NULL
    
      gg <- plot_group_profile(df, group_col = "cluster", style = "line", y_col = "value",
        mapping = aes(x = pos, y = value, color = cluster),
        layers = list(list(geom = "line", linewidth = 0.5))) +
        scale_color_manual(values = setNames(palette(rownames(P)), levels(df$cluster)),
                           guide = "none") +
        facet_wrap(~ cluster, ncol = 1, strip.position = "right") +
        labs(x = if (!is.null(x_label)) x_label else if (anchor$promoter) anchor$x_label else "genomic position",
             y = ylab, title = main, subtitle = sub) +
        theme_cowplot(font_size = 10)
      if (!is.null(tss))
        gg <- gg + geom_vline(xintercept = tss, linetype = "dashed", color = "grey40")
      gg
    })
    plot_panel(p, width = 8, height = 1.4 * r$n_clusters + 1)
  }),
  fiberseq_profiles_composition = panel_entry("fiberseq_profiles_composition.pdf", function(c) {
    p <- context_fiberseq_profiles(c)
    plot_panel(p$plot, width = 16, height = p$height)
  }),
  cluster_proportions = panel_entry("fiberseq_cluster_proportions.pdf", function(c) {
    plot_panel(plot_cluster_pie(c$result), width = 8, height = 5)
  }),
  umap_overview = panel_entry("fiberseq_umap_overview.pdf", context_umap_panel),
  haplotype_example = panel_entry("fiberseq_haplotype_example.pdf", context_haplotype_panel),
  knn_graph = panel_entry("knn_graph.pdf", function(c) {
    plot_panel(plot_knn_graph(c$result), width = 12, height = 9)
  }),
  knn_graph_by_allele = panel_entry("knn_graph_by_allele.pdf", function(c) {
    context_allele_summary(c)
    plot_panel(plot_knn_graph(c$result, color_by = "allele"), width = 12, height = 9)
  }),
  knn_allele_composition = panel_entry("knn_allele_composition.pdf", function(c) {
    summary <- context_allele_summary(c)
    p <- plot_stacked_proportion(summary$composition, "group", "allele",
      fiberseq_category_palette(summary$allele_levels), y = "fraction",
      position = "stack", reverse = TRUE, width = .75, legend = list(drop = FALSE),
      scales = list(ggplot2::scale_x_discrete(labels = function(x) paste0(x, "\n(n=",
        tapply(summary$composition$n_reads, summary$composition$group, sum)[x], ")")),
        ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1),
          expand = ggplot2::expansion(mult = c(0, .02)))),
      labels = list(title = paste(c$result$region$region_id, "Focal SNP allele proportions", sep = "\n"), x = "Leiden cluster", y = "Proportion of retained reads",
        fill = "Focal SNP allele", subtitle = "All reads: region-wide baseline; each bar sums to 100%",
        inside = list(min_fraction = .05, accuracy = 1, vjust = .5, color = "white", size = 3)),
      theme = ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
    plot_panel(p, width = max(10, .7 * nlevels(summary$composition$group) + 3), height = 6)
  }),
  knn_allele_connectivity = panel_entry("knn_allele_connectivity.pdf", function(c) {
    p <- plot_knn_allele_connectivity(context_allele_summary(c),
      paste(c$result$region$region_id, "Allele mixing in the KNN graph", sep = "\n"))
    plot_panel(p, width = 11, height = 7)
  }),
  methylation_proportion_by_cluster = panel_entry("methylation_proportion_by_cluster.pdf", function(c) {
    summary <- methylation_by_cluster(c$result)
    plot_table(c, summary, "methylation_proportion_by_cluster.tsv")
    plot_panel(plot_methylation_by_cluster(summary, c$result$region), width = 8, height = 4.5)
  }),
  heatmap_m6a_footprints = panel_entry("heatmap_m6a_footprints.pdf", context_genomic_heatmap),
  read_footprints = panel_entry("fig1_read_footprints.pdf", context_read_footprints),
  occupancy_by_cluster = panel_entry("fig2_occupancy_by_cluster.pdf", context_occupancy),
  cluster_sample_composition = panel_entry("cluster_sample_composition.pdf", function(c) {
    composition <- category_composition_inputs(c$result, "sample_label", c$sample_colors)
    p <- plot_stacked_proportion(composition, "cluster", "category", c$sample_colors,
      position = "fill", width = .75, legend = list(drop = FALSE),
      scales = list(ggplot2::scale_x_discrete(labels = function(x) paste0(x, "\n(n=", table(composition$cluster)[x], ")")),
        ggplot2::scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, .25),
          expand = ggplot2::expansion(mult = c(0, .02)))),
      labels = list(x = "Cluster", y = "Sample composition within cluster", fill = "LCL sample",
        title = c$result$region$annotation, subtitle = "Each bar sums to 100% of retained reads in that cluster"),
      theme = ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "right"))
    plot_panel(p, width = 12, height = 5)
  }),
  pairwise_m6a_delta = panel_entry("pairwise_m6a_delta.pdf", context_pairwise_delta,
    enabled = function(c) {
      limit <- c$ds$plot_options$delta_top_n
      rank <- c$result$region$rank
      is.null(limit) || is.null(rank) || (!is.na(rank) && rank <= limit)
    })
)
