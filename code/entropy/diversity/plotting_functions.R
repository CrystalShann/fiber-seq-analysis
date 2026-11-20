# Visualization only. Jitter changes display positions, never the analyzed cohort.
plot_diversity <- function(d, output_dir, groups) {
  labels <- c(rao_total = "Total Rao diversity per base",
    rao_between = "Between-cluster Rao diversity per base",
    between_cluster_share = "Between-cluster share",
    q_over_g = "Weighted between-cluster distance per base (B/G)")
  files <- c(rao_total = "total_rao_diversity",
    rao_between = "between_cluster_rao_diversity",
    between_cluster_share = "between_cluster_share_of_total_diversity",
    q_over_g = "weighted_between_cluster_distance")
  d$gene_group <- factor(d$gene_group, levels = groups)
  colors <- setNames(c('#1B9E77', '#D95F02', '#7570B3'), groups)
  group_labels <- setNames(c('Housekeeping', 'Cell-type-specific', 'Not expressed'), groups)
  # Same style as the nucleosome entropy 01_gene_group_entropy.pdf; black diamond = median
  for (metric in names(labels)) {
    z <- d[is.finite(d[[metric]]), ]
    s <- do.call(rbind, lapply(split(z[[metric]], z$gene_group, drop = TRUE), function(v)
      data.frame(n_genes = length(v), median_value = median(v))))
    s$gene_group <- factor(rownames(s), levels = groups)
    label_y <- max(z[[metric]]) + .05 * diff(range(z[[metric]]))
    p <- ggplot2::ggplot(z, ggplot2::aes(x = gene_group, y = .data[[metric]], colour = gene_group)) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = gene_group), width = .46, outlier.shape = NA,
                            alpha = .15, show.legend = FALSE) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = .13, seed = 12),
                          size = 1.8, alpha = .65, show.legend = FALSE) +
      ggplot2::geom_point(data = s, ggplot2::aes(gene_group, median_value), inherit.aes = FALSE,
                          shape = 18, size = 4) +
      ggplot2::geom_text(data = s, ggplot2::aes(gene_group, label_y,
                         label = sprintf('Median %.4g', median_value)), inherit.aes = FALSE, size = 3.6) +
      ggplot2::scale_colour_manual(values = colors) + ggplot2::scale_fill_manual(values = colors) +
      ggplot2::scale_x_discrete(labels = setNames(sprintf('%s\nn = %d genes',
        group_labels[as.character(s$gene_group)], s$n_genes), as.character(s$gene_group))) +
      ggplot2::labs(x = NULL, y = labels[[metric]], title = paste(labels[[metric]], 'by gene group')) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = 'bold'))
    ggplot2::ggsave(file.path(output_dir, paste0(files[[metric]], '.pdf')), p,
      device = 'pdf', width = 8.5, height = 6.5, useDingbats = FALSE)
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
    ggplot2::labs(x = sprintf('Fold enrichment relative to all %d scored genes', chosen$n_background[1]), y = NULL,
      size = 'Gene count', title = 'GO Biological Process by total Rao diversity',
      subtitle = 'Top eight terms per category; terms shown are not necessarily significant',
      caption = 'Exact overrepresentation tests; BH across all category-term pairs; HLA genes included') +
    ggplot2::theme_bw(base_size = 10)
  ggplot2::ggsave(path, p, device = 'pdf', width = 12, height = 14, useDingbats = FALSE)
}
