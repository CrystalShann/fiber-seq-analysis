# Co-accessible-pair display using the same plot_smf_reads() renderer, colours,
# fonts and SNP-relative coordinates as auto_correlation/06_autocorrelation.Rmd.
# Source leiden_manhattan_plots.r first. Input reads are already filtered;
# cluster is a factor of SNP-oriented haplotype groups, not inferred clusters.
plot_phased_coaccess_tracks <- function(result, records, pair) {
  a <- result$assignments
  stopifnot(nrow(a) > 0L, !anyDuplicated(a$RID),
    all(a$haplotype %in% c("HP1", "HP2")),
    all(a$allele_status == "phased_focal_genotype"),
    all(a$start <= result$region$analysis_start),
    all(a$end >= result$region$analysis_end))
  region <- result$region
  anchor <- plot_anchor(region)
  sample_colors <- sample_palette(a$sample_name)
  bounds <- c(anchor$left - .12 * (anchor$right - anchor$left), anchor$right + .5)
  cre <- data.frame(
    left = c(pair$cre1_start, pair$cre2_start) + 1L - anchor$anchor,
    right = c(pair$cre1_end, pair$cre2_end) - anchor$anchor,
    label = paste0(c("CRE1 (", "CRE2 ("), sub(".*\\.", "", c(pair$CRE1, pair$CRE2)), ")"),
    color = c("#FF0000", "#FFA500"))
  annotation <- ggplot2::ggplot(cre) +
    ggplot2::geom_rect(ggplot2::aes(xmin = left, xmax = right, ymin = .3, ymax = .6,
                                  fill = color)) +
    ggplot2::geom_text(ggplot2::aes(x = (left + right) / 2, y = .8, label = label), size = 3) +
    ggplot2::geom_vline(xintercept = 0, color = "#984EA3", linetype = "dashed", linewidth = .4) +
    ggplot2::annotate("point", x = 0, y = .1, shape = 18, size = 3, color = "#984EA3") +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(limits = bounds, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = paste0(region$annotation, " | full-window, phased molecules"),
      subtitle = sprintf("%s:%s %s>%s | dashed line at SNP (0 bp) | %d molecules",
        region$chr, format(region$focal_pos, big.mark = ",", scientific = FALSE),
        region$ref, region$alt, nrow(a)), x = NULL, y = NULL) +
    cowplot::theme_cowplot(font_size = 10) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank())
  draw <- function(tracks, ylabel) {
    plot_smf_reads(result, records, sample_colors, tracks = tracks, include_haplotype = TRUE) +
      ggplot2::geom_vline(xintercept = 0, color = "#984EA3", linetype = "dashed", linewidth = .4) +
      ggplot2::geom_vline(xintercept = c(cre$left, cre$right), color = "grey60",
                          linetype = "dotted", linewidth = .25) +
      ggplot2::labs(title = NULL, subtitle = NULL, y = ylabel,
        caption = "Left bars: HP1/HP2 (sample-local), SNP group, sample. Dotted lines: cCRE edges.")
  }
  m6a <- draw("m6A", "m6A: one row per molecule")
  footprints <- draw(LCL_FOOTPRINT_TRACKS, "Footprints: same molecule order")
  read_height <- max(3.5, .035 * nrow(a) + 1.8)
  combined <- cowplot::plot_grid(annotation, m6a, footprints, ncol = 1,
    align = "v", axis = "lr", rel_heights = c(1.1, read_height, read_height))
  list(panels = list(combined = combined, m6a = m6a, footprints = footprints),
       height = 1.1 + 2 * read_height)
}
