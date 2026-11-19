# Visualization only. Jitter changes display positions, never the analyzed cohort.
plot_diversity <- function(d, output_dir, groups) {
  labels <- c(rao_total = "Total Rao diversity per base",
    rao_between = "Between-cluster Rao diversity per base",
    between_cluster_share = "Between-cluster share",
    q_over_g = "Weighted between-cluster distance per base (B/G)",
    n_clusters = "Number of Leiden clusters")
  d$gene_group <- factor(d$gene_group, levels = groups)
  colors <- setNames(c('#1B9E77', '#D95F02', '#7570B3'), groups)
  for (metric in names(labels)) {
    p <- ggplot2::ggplot(d, ggplot2::aes(x = gene_group, y = .data[[metric]], fill = gene_group)) +
      ggplot2::geom_violin(trim = FALSE, alpha = .35, na.rm = TRUE) +
      ggplot2::geom_boxplot(width = .16, outlier.shape = NA, alpha = .75, na.rm = TRUE) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = .10, height = 0, seed = 2026),
                          size = 1.5, alpha = .65, na.rm = TRUE) +
      ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
      ggplot2::scale_x_discrete(drop = FALSE, labels = c('Housekeeping', 'Cell-type-specific', 'Non-expressed')) +
      ggplot2::labs(x = NULL, y = labels[[metric]],
        caption = 'One point per gene; all eligible fibers; original Leiden assignments') +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(legend.position = 'none', plot.caption = ggplot2::element_text(size = 9))
    ggplot2::ggsave(file.path(output_dir, paste0(metric, '.pdf')), p,
      device = 'pdf', width = 7, height = 5, useDingbats = FALSE)
  }
}

plot_rao_go <- function(results, path) {
  # Show up to eight smallest-p terms per category; display FDR even if nonsignificant.
  chosen <- do.call(rbind, lapply(c('low','intermediate','high'), function(category) {
    z <- results[results$rao_category == category & results$n_hits > 0, ]
    head(z[order(z$p_adjusted_BH, z$p_value, -z$fold_enrichment), ], 8L)
  }))
  chosen$rao_category <- factor(chosen$rao_category, levels = c('low','intermediate','high'))
  chosen$term_key <- paste(chosen$rao_category, chosen$go_id, sep = ':')
  chosen$term_key <- factor(chosen$term_key, levels = rev(chosen$term_key))
  labels <- setNames(vapply(chosen$description, function(x)
    paste(strwrap(x, width = 52), collapse = '\n'), character(1)), as.character(chosen$term_key))
  p <- ggplot2::ggplot(chosen, ggplot2::aes(fold_enrichment, term_key)) +
    ggplot2::geom_point(ggplot2::aes(size = n_hits, color = -log10(p_adjusted_BH))) +
    ggplot2::facet_grid(rao_category ~ ., scales = 'free_y', space = 'free_y') +
    ggplot2::scale_y_discrete(labels = labels) +
    ggplot2::scale_color_viridis_c(name = '-log10(BH FDR)') +
    ggplot2::labs(x = 'Fold enrichment relative to all 149 scored genes', y = NULL,
      size = 'Gene count', title = 'GO Biological Process by total Rao diversity',
      subtitle = 'Top eight terms per category; terms shown are not necessarily significant',
      caption = 'Exact overrepresentation tests; BH across all category-term pairs; HLA genes included') +
    ggplot2::theme_bw(base_size = 10)
  ggplot2::ggsave(path, p, device = 'pdf', width = 12, height = 14, useDingbats = FALSE)
}
