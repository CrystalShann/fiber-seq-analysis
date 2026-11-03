# Figures of the pooled Pol II clustering, built with the shared Leiden-Manhattan
# helpers (plot_read_heatmap(), plot_group_profile(), plot_stacked_proportion(),
# stack_rows(), plot_molecules(), cluster_id_palette(), sample_palette(),
# save_figure()). Load them first with pol2_load_helpers(). Clusters are drawn in
# numeric cluster order, except the enrichment composite (saved PI order).
# plot_molecules() calls ggplot2 functions unqualified.
suppressPackageStartupMessages(library(ggplot2))


###########################
# Order cluster labels
############################
# 1. Input: Cluster labels.
# 2. Output: Unique labels in numeric order.
#
# 3. Steps:
#    - Remove duplicate labels
#    - Sort by the number in each label.
pol2_numeric_clusters <- function(clusters) {
  clusters <- unique(as.character(clusters))
  clusters[order(as.integer(sub("^cluster", "", clusters)))]
}

# Change display factors on copies only; saved assignments and PI ranks stay intact.
###########################
# Order annotation tables
############################
# 1. Input: Cluster annotation tables.
# 2. Output: Copies with numeric cluster ordering.
#
# 3. Steps:
#    - Copy the tables
#    - Apply consistent cluster levels
#    - Sort the cluster summary.
pol2_numeric_annotation <- function(annotation) {
  out <- lapply(annotation, data.table::copy)
  clusters <- pol2_numeric_clusters(out$summary$cluster)
  for (name in names(out)) if ("cluster" %in% names(out[[name]])) {
    out[[name]][, cluster := factor(cluster, levels = clusters)]
  }
  data.table::setorder(out$summary, cluster)
  out
}

POL2_X_LABEL <- "Position relative to FANTOM5 dominant TSS (bp); upstream < 0"

# TSS-relative stand-in region for plot_read_heatmap(): positions are already
# strand-oriented, so tss = 0 on "+" leaves them unchanged.
###########################
# Describe the promoter axis
############################
# 1. Input: No arguments.
# 2. Output: A promoter description for shared plotting tools.
#
# 3. Steps:
#    - Use the strand-oriented promoter window
#    - Place the TSS at zero.
pol2_plot_region <- function() data.frame(chr = "TSS", analysis_start = -1000L, analysis_end = 999L,
  width = 2000L, region_type = "promoter", strand = "+", tss = 0L, annotation = "")

###########################
# Choose sample colors
############################
# 1. Input: Workflow settings and sample labels in the plot.
# 2. Output: Colors for the represented samples.
#
# 3. Steps:
#    - Read the sample list
#    - Build the shared palette
#    - Keep the needed entries.
pol2_sample_colors <- function(cfg, present) {
  colors <- sample_palette(fread(cfg$samples)$sample_name)
  colors[names(colors) %in% sub("_.*$", "", present)]
}

# Rows in heatmap order: cluster, then whole-window m6A proportion (passed as
# `start`, the within-cluster sort key of plot_read_heatmap()).
###########################
# Prepare heatmap rows
############################
# 1. Input: Saved clustering, annotations, row limit and random seed.
# 2. Output: Ordered molecule annotations and the base-level call matrix.
#
# 3. Steps:
#    - Check molecule identities
#    - Sample rows if needed
#    - Sort by cluster and accessibility.
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
###########################
# Draw the TSS heatmap
############################
# 1. Input: A heatmap, cluster labels and a profile-display choice.
# 2. Output: A drawn heatmap with TSS markers.
#
# 3. Steps:
#    - Draw the heatmap
#    - Mark the TSS in each cluster section
#    - Mark profile panels when requested.
pol2_draw_heatmap <- function(heatmap, clusters, profiles = FALSE) {
  ComplexHeatmap::draw(heatmap, newpage = FALSE)
  tss_x <- (1001 - 0.5) / 2000
  ###########################
  # Mark the TSS
  ############################
  # 1. Input: The TSS position from the surrounding heatmap function.
  # 2. Output: A vertical line on the active panel.
  #
  # 3. Steps:
  #    - Draw a dashed line through the TSS position.
  tss_line <- function() grid::grid.lines(x = c(tss_x, tss_x), y = c(0, 1),
    gp = grid::gpar(col = "#D55E00", lty = "dashed", lwd = 1.2))
  for (slice in seq_along(clusters)) ComplexHeatmap::decorate_heatmap_body("m6A", tss_line(), slice = slice)
  if (profiles) for (cl in clusters) ComplexHeatmap::decorate_annotation(paste0(cl, " m6A"), tss_line())
  invisible(NULL)
}

# Leiden-Manhattan read heatmap: per-bp m6A calls, cluster and sample bars,
# per-cluster call fraction lines on top.
###########################
# Build the accessibility heatmap
############################
# 1. Input: Prepared heatmap data and sample colors.
# 2. Output: A molecule accessibility heatmap.
#
# 3. Steps:
#    - Use the shared heatmap tool
#    - Add sample and cluster bars
#    - Show cluster profiles above the molecules.
pol2_molecule_heatmap <- function(res, sample_colors) {
  plot_read_heatmap(res, pol2_plot_region(), sample_colors, show_cluster_profiles = TRUE,
    heatmap_options = list(column_title = paste0(
      "Pol II clustering: ", nrow(res$assignments), " footprint-positive molecules, 31 LCL samples pooled\n",
      "Top: m6A call fraction per cluster\n", POL2_X_LABEL), column_title_gp = grid::gpar(fontsize = 11)))
}

# Internal enrichment-composite panel: 50-bp accessibility profiles with PI labels.
###########################
# Show cluster accessibility
############################
# 1. Input: Cluster annotations and profiles.
# 2. Output: Stacked accessibility panels for the enrichment figure.
#
# 3. Steps:
#    - Use the saved display order
#    - Plot each cluster's profile
#    - Add molecule counts and pausing-index labels.
pol2_enrichment_profiles <- function(annotation) {
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
  {
    text <- data.frame(cluster = factor(s$cluster, levels = s$cluster),
      label = ifelse(is.na(s$median_pi), "PI: NA", sprintf("PI: %.1f (%d TSSs)", s$median_pi, s$n_tss_with_pi)))
    p <- p + ggplot2::geom_text(data = text, ggplot2::aes(x = 980, y = 0.9, label = label),
                                inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.8)
  }
  p
}

# Pol II ChIP fold-enrichment at the TSS (-50..+299), log2 of cluster mean over
# the all-molecule mean, with BH stars.
###########################
# Show Pol II enrichment
############################
# 1. Input: The cluster summary.
# 2. Output: A horizontal enrichment plot.
#
# 3. Steps:
#    - Compare cluster signal with the overall mean
#    - Match the cluster display order
#    - Add significance labels.
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

###########################
# Show footprint positions
############################
# 1. Input: Cluster annotations and candidate footprints.
# 2. Output: Footprint-position histograms by cluster.
#
# 3. Steps:
#    - Group strand-relative midpoints into bins
#    - Draw one panel per cluster
#    - Mark the TSS.
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

# Nucleosome and TF footprints (tables/molecule_footprints.tsv.gz, `tracks`
# stage) of every heatmap molecule, one row per molecule in the heatmap's row
# order, drawn with plot_molecules() in the MOLECULE_TRACK_COLORS of each class.
###########################
# Show footprints on molecules
############################
# 1. Input: Heatmap row information and extracted footprints.
# 2. Output: Footprint tracks for individual molecules.
#
# 3. Steps:
#    - Match footprints to displayed molecules
#    - Draw intervals in the heatmap row order
#    - Color each footprint class.
pol2_plot_footprint_molecules <- function(res, footprints) {
  rows <- stack_rows(data.frame(RID = res$assignments$RID, cluster = res$assignments$cluster), "cluster")
  fp <- footprints[observation_id %in% rows$RID]
  classes <- intersect(names(MOLECULE_TRACK_COLORS), fp$class)
  stopifnot(length(classes) > 0L, all(fp$class %in% classes))
  tracks <- data.frame(RID = fp$observation_id, start = fp$rel_start, end = fp$rel_end,
                       class = factor(fp$class, levels = classes))
  # nucleosomes first, TF footprints drawn on top
  tracks <- tracks[order(tracks$class), , drop = FALSE]
  n <- table(rows$cluster)
  # start = end = 0 puts plot_molecules()' dashed region lines on the TSS.
  ex <- data.frame(chrom = "TSS", start = 0, end = 0, window_start = -1000, window_end = 1000)
  plot_molecules(rows, tracks, ex, facet_col = "cluster",
    facet_labels = setNames(paste0(names(n), "\n(n=", n, ")"), names(n)),
    bar_cols = "cluster", bar_palette = cluster_id_palette(levels(rows$cluster)),
    track_colors = MOLECULE_TRACK_COLORS[classes],
    title = paste0("Nucleosome and FiberHMM TF footprints of ", nrow(rows), " clustered molecules"),
    subtitle = paste("Rows in molecule_heatmap.pdf order; clusters in numeric order;",
                     "nucleosomes: fibertools nuc calls, TF footprints: FiberHMM v2 recalled_tf")) +
    ggplot2::labs(x = POL2_X_LABEL)
}

# Fig. 1c wedge beside the profiles: clusters run from lowest median PI (top)
# to highest (bottom), so the wedge widens downward.
###########################
# Draw the pausing-index guide
############################
# 1. Input: No arguments.
# 2. Output: A labeled triangular guide.
#
# 3. Steps:
#    - Draw a triangle
#    - Add the pausing-index label.
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
###########################
# Combine cluster views
############################
# 1. Input: Heatmap data, annotations and sample colors.
# 2. Output: A combined accessibility and enrichment figure.
#
# 3. Steps:
#    - Restore the saved pausing-index order
#    - Draw the molecule heatmap and cluster profiles
#    - Add the pausing-index guide and enrichment panel.
pol2_clusters_and_enrichment <- function(res, annotation, sample_colors) {
  annotation <- lapply(annotation, data.table::copy)
  data.table::setorder(annotation$summary, display_order)
  cluster_order <- as.character(annotation$summary$cluster)
  res$assignments$cluster <- factor(res$assignments$cluster, levels = cluster_order)
  res$assignments <- res$assignments[order(res$assignments$cluster, res$assignments$start, res$assignments$RID), , drop = FALSE]
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
  profiles <- cowplot::plot_grid(pol2_enrichment_profiles(annotation), pol2_plot_pi_wedge(),
                                 nrow = 1, align = "h", axis = "tb", rel_widths = c(1, 0.06))
  cowplot::plot_grid(grob, profiles, pol2_plot_signal_enrichment(annotation$summary),
                     nrow = 1, rel_widths = c(1.3, 1.06, 0.7))
}

# Every eligible molecule at one gene's FANTOM5 TSS (its chromosome feature
# cache), as a plot_locus_molecules() result. Clustered molecules keep their
# Leiden cluster; every other molecule gets pol2_knn_clusters()
# (pol2_molecule_tables.R), the rule of the `project` stage.
# Per-bp calls move from the -1000..999 axis to genomic columns.
###########################
# Prepare one gene's molecules
############################
# 1. Input: Workflow settings, saved clusters, annotations and a gene.
# 2. Output: Molecules and genomic calls for one promoter.
#
# 3. Steps:
#    - Find the gene's analysed TSS
#    - Keep saved labels and assign remaining molecules
#    - Convert relative positions to genomic coordinates.
pol2_locus_result <- function(cfg, x, annotation, gene) {
  site <- pol2_read_prepared(cfg)$tss[gene_name == gene]
  if (nrow(site) != 1L) stop("Expected one analysed TSS for ", gene, ", found ", nrow(site))
  f <- readRDS(file.path(cfg$output, "tables", "features", paste0(site$chrom, ".rds")))
  meta <- f$molecules[tss_id == site$tss_id]
  stopifnot(nrow(meta) > 0L, !anyDuplicated(meta$RID))
  m <- annotation$molecules
  leiden <- m$cluster[match(meta$observation_id, m$observation_id)]
  cluster <- as.character(leiden)
  todo <- which(is.na(cluster))
  if (length(todo)) cluster[todo] <- pol2_knn_clusters(f$features[meta$observation_id[todo], , drop = FALSE],
    x$features, x$molecules$cluster, cfg$k_neighbors)
  # numeric copy of the pattern matrix for m6a_intervals()
  sites <- f$sites[meta$observation_id, , drop = FALSE] * 1
  pos <- site$tss + ifelse(site$strand == "+", 1L, -1L) * as.integer(colnames(sites))
  sites <- sites[, order(pos), drop = FALSE]
  dimnames(sites) <- list(meta$RID, sort(pos))
  region <- data.frame(region_id = paste0(gene, "_TSS"), region_type = "promoter", chr = site$chrom,
    start = site$window_start, end = site$window_end, analysis_start = site$window_start + 1L,
    analysis_end = site$window_end, width = site$window_end - site$window_start,
    strand = site$strand, tss = site$tss, annotation = sprintf("%s | FANTOM5 dominant TSS %s:%s (%s), hg38",
      gene, site$chrom, format(site$tss, big.mark = ","), site$strand))
  assignments <- data.frame(RID = meta$RID, original_RID = meta$original_RID, sample_name = meta$sample_name,
    start = meta$read_start, end = meta$read_end, clustered = !is.na(leiden),
    cluster = droplevels(factor(cluster, levels = levels(m$cluster))))
  list(region = region, assignments = assignments, site_met_mat = sites, site = site)
}

# GTEx eVariants of the gene inside the TSS window (eqtl = signif pairs table);
# the most significant one is labelled.
###########################
# Find nearby expression variants
############################
# 1. Input: An optional expression-variant table and promoter data.
# 2. Output: Variants in the promoter window, or no result.
#
# 3. Steps:
#    - Match the gene
#    - Keep variants inside the window
#    - Label the most significant variant.
pol2_locus_eqtl <- function(eqtl, res) {
  if (is.null(eqtl)) return(NULL)
  site <- res$site
  e <- eqtl[sub("\\..*$", "", gene_id) == sub("\\..*$", "", site$gene_id)]
  if (!nrow(e)) return(NULL)
  e[, c("chrom", "pos", "ref", "alt") := tstrsplit(variant_id, "_", fixed = TRUE, keep = 1:4)]
  e[, pos := as.integer(pos)]
  e <- e[chrom == site$chrom & pos >= res$region$analysis_start & pos <= res$region$analysis_end][order(pval_nominal)]
  if (!nrow(e)) return(NULL)
  e[, label := NA_character_]
  e[1L, label := sprintf("lead %s:%s %s>%s, p = %.1e", chrom, format(pos, big.mark = ","), ref, alt, pval_nominal)]
  as.data.frame(e[, .(pos, label, variant_id, pval_nominal)])
}

# plot_locus_molecules() for one gene, without CTCF tracks: genes, its in-window
# GTEx eVariants, cluster proportions at the TSS and in all clustered
# molecules, per-cluster m6A profiles and every molecule.
###########################
# Build one gene's figure
############################
# 1. Input: Workflow settings, clustering, annotations, gene models and optional variants.
# 2. Output: A locus figure with its supporting data.
#
# 3. Steps:
#    - Prepare molecules and variants
#    - Summarize cluster and sample counts
#    - Combine genes, profiles and molecule tracks.
pol2_locus_figure <- function(cfg, x, annotation, gene, models, eqtl = NULL) {
  res <- pol2_locus_result(cfg, x, annotation, gene)
  a <- res$assignments
  marks <- pol2_locus_eqtl(eqtl, res)
  n <- table(a$cluster)
  samples <- tapply(a$sample_name, a$cluster, function(s) length(unique(s)))
  ###########################
  # Write count labels
  ############################
  # 1. Input: A count and a word.
  # 2. Output: A readable singular or plural label.
  #
  # 3. Steps:
  #    - Place the count before the word
  #    - Add a plural ending when needed.
  plural <- function(k, word) paste0(k, " ", word, ifelse(k == 1, "", "s"))
  labels <- setNames(sprintf("%s\n%s (%.0f%%)\n%s", names(n), plural(as.integer(n), "molecule"),
    100 * as.integer(n) / nrow(a), plural(as.integer(samples[names(n)]), "sample")), names(n))
  s <- annotation$summary
  fig <- plot_locus_molecules(res, models, cluster_id_palette(s$cluster), marks = marks,
    strip_label = "GTEx LCL\neVariants", group_labels = labels, show_proportions = TRUE,
    reference_counts = setNames(s$n, s$cluster), region_label = paste(gene, "TSS"),
    reference_label = "All clustered molecules",
    subtitle = paste0("Top: canonical genes, TSS window shaded | bottom: TSS window, dashed line at the TSS",
      if (!is.null(marks)) paste0("; diamonds: GTEx v8 LCL eVariants of ", gene) else "",
      " | clusters in numeric order"),
    caption = sprintf(paste("Left bars: cluster, sample. %d of %d molecules were clustered and keep their",
      "Leiden cluster; the others take the majority cluster of their %d nearest clustered molecules",
      "(Manhattan distance, 40 bins of 50 bp)."), sum(a$clustered), nrow(a), cfg$k_neighbors),
    x_label = POL2_X_LABEL)
  c(fig, list(result = res, eqtl = marks))
}

# Locus figures of cfg$locus_genes, named by gene; empty without locus_genes.
###########################
# Build selected gene figures
############################
# 1. Input: Workflow settings, saved clustering and annotations.
# 2. Output: Figures named by gene, or an empty list.
#
# 3. Steps:
#    - Read gene models and optional variants
#    - Build a figure for each requested gene.
pol2_locus_figures <- function(cfg, x, annotation) {
  if (!length(cfg$locus_genes)) return(list())
  models <- load_canonical_gene_models(cfg$gene_models)
  eqtl <- if (is.null(cfg$eqtl_pairs)) NULL else
    fread(cfg$eqtl_pairs, select = c("variant_id", "gene_id", "pval_nominal"))
  setNames(lapply(cfg$locus_genes, function(g) pol2_locus_figure(cfg, x, annotation, g, models, eqtl)),
           cfg$locus_genes)
}

###########################
# Save the main figures
############################
# 1. Input: Saved clustering, annotations and workflow settings.
# 2. Output: Main clustering, footprint, imprinting, feature and locus figures.
#
# 3. Steps:
#    - Prepare shared colors and row ordering
#    - Build each figure
#    - Save figures and extracted feature tables.
pol2_save_figures <- function(x, annotation, cfg) {
  annotation <- pol2_numeric_annotation(annotation)
  pol2_required(file.path(cfg$output, "tables", "tracks", c("molecule_projection.tsv.gz", "molecule_footprints.tsv.gz")))
  res <- pol2_heatmap_result(x, annotation, cfg$heatmap_max_reads, cfg$seed)
  colors <- pol2_sample_colors(cfg, res$assignments$sample_name)
  clusters <- levels(res$assignments$cluster)
  k <- length(clusters)
  plots <- file.path(cfg$output, "plots")
  masks <- pol2_observed_masks(cfg,x)
  profile <- pol2_observed_profile(x$sites,res$assignments,masks,
    x$molecules$tss_id[match(res$assignments$RID,x$molecules$observation_id)])
  save_figure(pol2_observed_heatmap(res, colors, profile), file.path(plots, "molecule_heatmap.pdf"),
              width = 14, height = 11 + 0.55 * k, draw = function(h) pol2_draw_heatmap(h, clusters, profiles = TRUE))
  save_figure(pol2_plot_footprints(annotation), file.path(plots, "footprint_positions.pdf"),
              width = 10, height = 1.4 * k + 2.5)
  tables <- file.path(cfg$output, "tables")
  footprints <- fread(file.path(tables, "tracks", "molecule_footprints.tsv.gz"))
  save_figure(pol2_plot_footprint_molecules(res, footprints), file.path(plots, "footprint_molecules.pdf"),
              width = 9, height = 11 + 0.55 * k)
  save_figure(pol2_clusters_and_enrichment(res, annotation, colors), file.path(plots, "pol2_clusters_and_enrichment.pdf"),
              width = 16, height = max(7, 0.9 * k + 2))
  prep <- pol2_read_prepared(cfg)
  imprinting <- pol2_imprinting(cfg, fread(file.path(tables, "tracks", "molecule_projection.tsv.gz")), prep$tss, clusters)
  save_figure(pol2_plot_imprinting(imprinting, clusters), file.path(plots, "imprinting_by_cluster.pdf"), width = 11, height = 4.5)
  features <- pol2_cluster_features(cfg, x, annotation, prep)
  dir.create(file.path(tables,"cluster_features"),recursive=TRUE,showWarnings=FALSE)
  fwrite(features$summary, file.path(tables, "cluster_features", "cluster_features.tsv"), sep = "\t")
  fwrite(features$cpg, file.path(tables, "cluster_features", "cluster_cpg_by_tss.tsv"), sep = "\t")
  fwrite(features$molecules[, .(observation_id, cluster, tss_id, m6a_proportion, cpg, grocap, dnase)],
    file.path(tables, "cluster_features", "cluster_feature_molecules.tsv.gz"), sep = "\t")
  save_figure(pol2_plot_cluster_features(features, clusters), file.path(plots, "cluster_features.pdf"), width = 12, height = 8)
  loci <- pol2_observed_loci(cfg, x, annotation, masks)
  for (gene in names(loci))
    save_figure(loci[[gene]]$combined, file.path(plots, "loci", paste0(gene, "_TSS_molecules.pdf")),
                width = 12, height = loci[[gene]]$height, limitsize = FALSE)
}

# Descriptive physical-molecule counts; no inferential tests.
###########################
# Summarize imprinted promoters
############################
# 1. Input: Workflow settings, projected molecules, TSSs and cluster order.
# 2. Output: Imprinting summaries and a promoter table held in memory.
#
# 3. Steps:
#    - Find promoters overlapping imprinted regions
#    - Count distinct molecules and cluster proportions
#    - Add expression, pausing and signal information.
pol2_imprinting <- function(cfg, projection, tss, clusters) {
  dmr <- fread(cfg$imprinted_dmrs, header = FALSE, select = 1:4,
               col.names = c("chrom", "start", "end", "dmr_name"))
  hits <- findOverlaps(pol2_gr(tss$chrom, tss$window_start, tss$window_end),
                       pol2_gr(dmr$chrom, dmr$start, dmr$end), ignore.strand = TRUE)
  sites <- unique(tss$tss_id[queryHits(hits)])
  p <- copy(projection)
  p[, imprinted := tss_id %in% sites]
  per_cluster <- p[, .(n = uniqueN(RID[imprinted]), N = uniqueN(RID)), by = cluster]
  per_cluster[, fraction := n / N]
  # A physical molecule observed at both types of TSS belongs only to imprinted.
  imprinted_rids <- unique(p[imprinted == TRUE, RID])
  pairs <- unique(p[, .(RID, cluster)])
  pairs[, set := fifelse(RID %in% imprinted_rids, "Imprinted TSSs", "Other TSSs")]
  composition <- pairs[, .(n = .N), by = .(set, cluster)]
  composition[, fraction := n / sum(n), by = set]
  overall <- pairs[, mean(cluster == "cluster4")]
  loci <- p[imprinted == TRUE, {
    u <- unique(.SD[, .(RID, cluster)])
    counts <- u[, .N, by = cluster][order(-N, cluster)]
    list(molecules = uniqueN(RID), samples = uniqueN(sample_name),
         cluster4_fraction = uniqueN(RID[cluster == "cluster4"]) / uniqueN(RID),
         dominant_cluster = counts$cluster[1L])
  }, by = tss_id]
  names <- data.table(tss_id = tss$tss_id[queryHits(hits)], dmr_name = dmr$dmr_name[subjectHits(hits)])
  names <- names[, .(dmr_name = paste(sort(unique(dmr_name)), collapse = ";")), by = tss_id]
  loci <- merge(merge(loci, names, by = "tss_id"), pol2_attach_pi(tss, cfg$pi_table), by = "tss_id")
  expr <- rbindlist(lapply(cfg$expression, function(f) fread(f, select = c("gene_id", "TPM"))))
  expr <- expr[grepl("^ENSG", gene_id) & !grepl("_PAR_Y$", gene_id),
               .(tpm = mean(TPM)), by = .(key = sub("\\..*$", "", gene_id))]
  loci[, `:=`(tpm = expr$tpm[match(sub("\\..*$", "", gene_id), expr$key)],
               fold_overall_cluster4 = cluster4_fraction / overall,
               overall_cluster4_fraction = overall,
               locus = paste0(chrom, ":", window_start + 1L, "-", window_end))]
  for (d in list(per_cluster, composition)) d[, cluster := factor(cluster, levels = clusters)]
  list(per_cluster = per_cluster, composition = composition,
       loci = loci[, .(tss_id, gene = gene_name, locus, dmr_name, molecules, samples,
         cluster4_fraction, overall_cluster4_fraction, fold_overall_cluster4, dominant_cluster,
         tpm, pi, pol2_tss_signal)])
}

###########################
# Compare imprinting patterns
############################
# 1. Input: Imprinting summaries and cluster order.
# 2. Output: A two-panel imprinting figure.
#
# 3. Steps:
#    - Show the imprinted fraction within each cluster
#    - Compare cluster composition between promoter groups.
pol2_plot_imprinting <- function(imprinting, clusters) {
  colors <- cluster_id_palette(clusters)
  a <- ggplot(imprinting$per_cluster, aes(cluster, fraction, fill = cluster)) +
    geom_col() + geom_text(aes(label = paste(n, N, sep = " / ")), vjust = -0.4, size = 3) +
    scale_fill_manual(values = colors, guide = "none") +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .18))) +
    labs(x = NULL, y = "Molecules at imprinted TSSs", title = "Imprinted fraction per cluster") +
    theme_fiberseq("classic", base_size = 10)
  composition <- copy(imprinting$composition)
  b <- ggplot(composition, aes(set, fraction, fill = cluster)) + geom_col() +
    geom_text(aes(label = sprintf("%.1f%%", 100 * fraction)), colour = "white",
              position = position_stack(vjust = .5), size = 3) +
    scale_fill_manual(values = colors) + scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Cluster composition", fill = NULL) + theme_fiberseq("classic", base_size = 10)
  cowplot::plot_grid(a, b, nrow = 1, labels = c("A", "B"))
}

###########################
# Measure cluster CpG methylation
############################
# 1. Input: Workflow settings, saved clusters and prepared inputs.
# 2. Output: CpG summaries per promoter and cluster.
#
# 3. Steps:
#    - Extract each cluster's molecule alignments
#    - Run CpG counting and check the output
#    - Summarize methylation within promoter windows.
pol2_cluster_cpg <- function(cfg, x, prep) {
  tmp <- tempfile("pol2_cpg_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  manifest <- merge(x$molecules[, .(RID, original_RID, sample_name, tss_id, cluster)],
    prep$tss[, .(tss_id, chrom, window_start, window_end)], by = "tss_id")
  manifest <- merge(manifest, prep$samples[, .(sample_name, merged_bam_file)], by = "sample_name")
  for (sample in names(cfg$merged_bam_overrides))
    manifest[sample_name == sample, merged_bam_file := cfg$merged_bam_overrides[[sample]]]
  pol2_required(unique(manifest$merged_bam_file))
  stopifnot(nrow(manifest) == nrow(x$molecules), !anyDuplicated(manifest$RID))
  path <- file.path(tmp, "manifest.tsv")
  fwrite(manifest, path, sep = "\t")
  helper <- file.path(cfg$project, "code/pausing/clustering_Pol2/04_pol2_cluster_bams.py")
  status <- system2(cfg$python, shQuote(c(helper, "--manifest", path, "--output", tmp)))
  if (status != 0) stop("Cluster BAM extraction failed")
  cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  results <- lapply(unique(as.character(manifest$cluster)), function(cl) {
    prefix <- file.path(tmp, cl)
    status <- system2(cfg$pb_cpg_tools, shQuote(c("--bam", paste0(prefix, ".bam"),
      "--output-prefix", prefix, "--pileup-mode", "count", "--modsites-mode", "reference",
      "--ref", cfg$reference, "--min-coverage", "1", "--threads", cores)))
    if (status != 0) stop("CpG extraction failed for ", cl)
    path <- paste0(prefix, ".combined.bed.gz")
    con <- gzfile(path, "rt")
    header <- readLines(con, n = 50L)
    close(con)
    header <- header[startsWith(header, "#chrom\t")]
    expected <- c("#chrom", "begin", "end", "mod_score", "type", "cov", "mod_count", "unmod_count",
                  "avg_mod_score", "avg_unmod_score")
    if (length(header) != 1L || !identical(strsplit(header, "\t", fixed = TRUE)[[1]], expected))
      stop("Unexpected count-mode CpG BED schema")
    bed <- fread(cmd = paste("zcat", shQuote(path), "| sed '/^#/d'"), header = FALSE)
    if (ncol(bed) != length(expected)) stop("Unexpected CpG BED columns")
    message("Confirmed count-mode CpG columns: ", header)
    setnames(bed, names(bed)[1:8], c("chrom", "start", "end", "mod_score", "type", "cov", "mod", "unmod"))
    stopifnot(all(bed$cov >= 1), all(bed$mod >= 0), all(bed$mod <= bed$cov),
              all(abs(bed$mod + bed$unmod - bed$cov) < 1e-5))
    windows <- unique(manifest[cluster == cl, .(tss_id, chrom, window_start, window_end)])
    h <- findOverlaps(pol2_gr(bed$chrom, bed$start, bed$end),
      pol2_gr(windows$chrom, windows$window_start, windows$window_end), ignore.strand = TRUE)
    z <- data.table(tss_id = windows$tss_id[subjectHits(h)], site = queryHits(h),
      mod = bed$mod[queryHits(h)], cov = bed$cov[queryHits(h)])
    by_tss <- z[, .(mod = sum(mod), cov = sum(cov), cpg = sum(mod) / sum(cov),
                    n_cpg_sites = .N, median_cpg_coverage = median(as.numeric(cov))), by = tss_id]
    by_tss <- merge(windows[, .(tss_id)], by_tss, by = "tss_id", all.x = TRUE)
    u <- unique(z, by = "site")
    list(by_tss = by_tss[, cluster := cl], summary = data.table(cluster = cl,
      pooled_cpg = sum(u$mod) / sum(u$cov), n_cpg_sites = nrow(u), median_cpg_coverage = median(u$cov)))
  })
  list(by_tss = rbindlist(lapply(results, `[[`, "by_tss")),
       summary = rbindlist(lapply(results, `[[`, "summary")))
}

###########################
# Combine cluster measurements
############################
# 1. Input: Workflow settings, clustering, annotations and prepared inputs.
# 2. Output: Molecule measurements and cluster feature summaries.
#
# 3. Steps:
#    - Measure CpG methylation
#    - Read GRO-cap and DNase signal
#    - Join measurements to molecules and summarize each cluster.
pol2_cluster_features <- function(cfg, x, annotation, prep) {
  cpg <- pol2_cluster_cpg(cfg, x, prep)
  tss <- prep$tss[tss_id %in% x$molecules$tss_id]
  windows <- tss[, .(tss_id, strand, chrom, start = window_start, end = window_end)]
  windows[, grocap := 0]
  for (st in c("+", "-")) {
    d <- pol2_signal(windows[strand == st], if (st == "+") cfg$grocap_plus else cfg$grocap_minus)
    value <- abs(d$pol2_mean) * d$pol2_covered_fraction * (d$end - d$start)
    value[d$pol2_covered_fraction == 0] <- 0
    windows[strand == st, grocap := value]
  }
  windows[, dnase := pol2_signal(windows, cfg$dnase)$pol2_mean]
  m <- merge(annotation$molecules, windows[, .(tss_id, grocap, dnase)], by = "tss_id")
  m <- merge(m, cpg$by_tss[, .(cluster, tss_id, cpg)], by = c("cluster", "tss_id"), all.x = TRUE)
  summary <- m[, .(n_molecules = .N, n_tss = uniqueN(tss_id),
    pooled_m6a = sum(m6a_calls) / sum(n_positions), median_m6a = median(m6a_proportion),
    median_grocap = median(grocap), median_dnase = median(dnase),
    mean_grocap = mean(grocap), mean_dnase = mean(dnase)), by = cluster]
  summary <- merge(summary, cpg$summary, by = "cluster")
  summary[, cluster := factor(cluster, levels = annotation$summary$cluster)]
  setorder(summary, cluster)
  list(molecules = m, summary = summary, cpg = cpg$by_tss)
}

###########################
# Compare cluster measurements
############################
# 1. Input: Feature tables and cluster order.
# 2. Output: A four-panel comparison of accessibility, CpG, GRO-cap and DNase.
#
# 3. Steps:
#    - Draw the molecule distributions
#    - Add pooled or mean values
#    - Arrange the panels together.
pol2_plot_cluster_features <- function(features, clusters) {
  m <- copy(features$molecules)
  s <- copy(features$summary)
  m[, cluster := factor(cluster, levels = clusters)]
  s[, cluster := factor(cluster, levels = clusters)]
  specs <- list(c("m6a_proportion", "pooled_m6a", "m6A fraction"),
    c("cpg", "pooled_cpg", "CpG methylation fraction"),
    c("grocap", "mean_grocap", "GRO-cap sense reads: log10(x + 1)"),
    c("dnase", "mean_dnase", "DNase: log10(x + 1)"))
  panels <- lapply(seq_along(specs), function(i) {
    spec <- specs[[i]]
    d <- m[, .(cluster, value = get(spec[1]))]
    point <- s[, .(cluster, value = get(spec[2]))]
    if (i > 2L) { d[, value := log10(value + 1)]; point[, value := log10(value + 1)] }
    ggplot(d, aes(cluster, value, fill = cluster)) + geom_boxplot(outlier.size = .3, na.rm = TRUE) +
      geom_point(data = point, shape = 23, fill = "white", size = 2.5) +
      scale_fill_manual(values = cluster_id_palette(clusters), guide = "none") +
      labs(x = "Cluster (numeric order)", y = spec[3],
        caption = if (i <= 2L) "Diamond: pooled fraction" else "Diamond: molecule-weighted mean") +
      theme_fiberseq("classic", base_size = 10)
  })
  cowplot::plot_grid(plotlist = panels, ncol = 2)
}


# Plot routing only: reuse the original autocorrelation plotting function.
###########################
# Build autocorrelation views
############################
# 1. Input: Workflow settings and whether to save figures.
# 2. Output: Upstream and downstream overlays, panels and heatmaps.
#
# 3. Steps:
#    - Load the original plotting function and saved summaries
#    - Use numeric cluster order and saved colors
#    - Build the views and optionally save them.
pol2_autocorrelation_plots <- function(cfg, save = FALSE) {
  tables <- file.path(cfg$output, "tables", "autocorrelation")
# Source the original named function chunk, without executing the report's analysis chunks.
original <- readLines(file.path(cfg$project,"code/clustering_methods/auto_correlation/06_autocorrelation.Rmd"))
start <- grep("^```\\{r report-cluster-acf-function,",original)
end <- start + which(original[(start+1L):length(original)] == "```")[1L]
eval(parse(text=original[(start+1L):(end-1L)]))
d <- data.table::fread(file.path(tables,"autocorrelation_cluster_summary.tsv.gz"))
clusters <- pol2_numeric_clusters(d$cluster)
acf_plots <- list()
for (window in c("upstream","downstream")) {
  entry <- list(tables=list(cluster_acf_summary=as.data.frame(d[d$window == window])),region=list(region_id=window))
  title <- if (window=="upstream") "Upstream: TSS -500 to -25 bp (strand-relative)" else "Downstream: TSS +25 to +500 bp (strand-relative)"
  for (kind in c("overlay","panels","heatmap")) {
    acf_plots[[paste(window,kind,sep="_")]] <- plot_report_cluster_acf(entry,lower=15,upper=475,title=title,
      groups=clusters,colors=cluster_id_palette(clusters),panels=kind=="panels",
      display=if(kind=="heatmap") "heatmap" else "curves",show_counts=FALSE,frame=FALSE)
    if (save) {
      suffix <- if(window=="upstream") "upstream_-500_to_-25_bp" else "downstream_+25_to_+500_bp"
      prefix <- if(kind=="heatmap") "autocorrelation_heatmap_" else paste0("autocorrelation_curves_",kind,"_")
      plots <- file.path(cfg$output,"plots","autocorrelation");dir.create(plots,recursive=TRUE,showWarnings=FALSE)
      ggplot2::ggsave(file.path(plots,paste0(prefix,suffix,".pdf")),acf_plots[[paste(window,kind,sep="_")]],
        width=10,height=if(kind=="panels") 10 else 5)
    }
  }
}

  acf_plots
}

# Plot routing only: reuse the original midpoint and profile plotting functions.
###########################
# Build nucleosome views
############################
# 1. Input: Workflow settings and whether to save figures.
# 2. Output: Occupancy, center-count, position and TSS-occupancy plots.
#
# 3. Steps:
#    - Load the original plotting tools and saved tables
#    - Keep numeric cluster order
#    - Build the four views and optionally save them.
pol2_nucleosome_plots <- function(cfg, save = FALSE) {
  nuc_env <- new.env()
  source(file.path(cfg$project,"code/nucleosome_positioning/process_nuc_pos.R"),local=nuc_env)
  source(file.path(cfg$project,"code/nucleosome_positioning/plot_nuc_pos.R"),local=nuc_env)
  tables <- file.path(cfg$output,"tables","nuc")
  positions <- data.table::fread(file.path(tables,"nucleosome_molecules.tsv.gz"))
  clusters <- pol2_numeric_clusters(positions$cluster)
  positions[,cluster:=factor(cluster,levels=clusters)]
  meta <- list(region_type="promoter",rel_view_start=-1000L,rel_view_end=999L,
    position_label="Position relative to FANTOM5 TSS (bp; strand-relative)")
  plots <- list()
  for (kind in c("occupancy","center_counts")) {
    d <- data.table::fread(file.path(tables,paste0("nucleosome_",kind,"_profiles.tsv.gz")))
    d[,cluster:=factor(cluster,levels=clusters)]
    meta$title <- if(kind=="occupancy") "Nucleosome occupancy" else "Nucleosome center counts"
    label <- if(kind=="occupancy") "Occupancy at 10-bp grid positions" else "Raw inferred-center count (25-bp bins)"
    plots[[kind]] <- nuc_env$plot_position_summary(positions,positions,meta,include_gap=FALSE,
      group_col="cluster",colors=cluster_id_palette(clusters),profile=d,profile_y="value",profile_label=label)
  }
  meta$title <- "Nearest upstream/downstream nucleosome by midpoint"
  plots$positions <- nuc_env$plot_position_summary(positions,positions,meta,include_gap=FALSE,
    group_col="cluster",colors=cluster_id_palette(clusters),x_label="Cluster")
  summary <- data.table::fread(file.path(tables,"nucleosome_cluster_summary.tsv.gz"))
  summary[,cluster:=factor(cluster,levels=clusters)]
  meta$title <- "TSS nucleosome occupancy"
  plots$tss_occupancy <- nuc_env$plot_position_summary(positions,positions,meta,include_gap=FALSE,
    colors=cluster_id_palette(clusters),profile=summary,profile_y="tss_occupancy")
  if (save) {
    output <- file.path(cfg$output,"plots","nuc")
    dir.create(output,recursive=TRUE,showWarnings=FALSE)
    for (kind in c("occupancy","center_counts"))
      ggplot2::ggsave(file.path(output,paste0("nucleosome_",kind,"_molecule.pdf")),plots[[kind]],width=10,height=10)
    ggplot2::ggsave(file.path(output,"nucleosome_positions.pdf"),plots$positions,width=11,height=5)
    ggplot2::ggsave(file.path(output,"nucleosome_tss_occupancy.pdf"),plots$tss_occupancy,width=9,height=5)
  }
  plots
}

# Bulk GM12878 histone ChIP (cfg$histone_marks: ENCODE GRCh38 fold change over control)
# around each clustered molecule's TSS. Bins are 50 bp over -2000..+1999 in transcriptional
# orientation; each bin is the coverage-weighted mean signal (pol2_signal()). Cluster curves
# average over the cluster's molecules, so molecules at one TSS share its locus signal.
###########################
# Summarize histone signals
############################
# 1. Input: Workflow settings, molecule metadata, flank length and bin size.
# 2. Output: Mean histone profiles by cluster and across all molecules.
#
# 3. Steps:
#    - Create strand-oriented bins around each TSS
#    - Read each histone signal
#    - Average the signal over molecules.
pol2_histone_profiles <- function(cfg, m, flank = 2000L, bin = 50L) {
  tss <- unique(m[, .(tss_id, chrom, tss, strand)])
  bins <- tss[, .(relative_start = seq(-flank, flank - bin, by = bin)), by = .(tss_id, chrom, tss, strand)]
  # Inclusive offsets relative_start..relative_start + bin - 1 as 0-based BED; tss is 1-based.
  bins[, `:=`(start = ifelse(strand == "+", tss - 1L + relative_start, tss - relative_start - bin),
              end = ifelse(strand == "+", tss - 1L + relative_start + bin, tss - relative_start))]
  rbindlist(lapply(names(cfg$histone_marks), function(mark) {
    signal <- pol2_signal(bins, cfg$histone_marks[[mark]])
    v <- merge(m[, .(tss_id, cluster = as.character(cluster))],
               signal[, .(tss_id, relative_start, fe = pol2_mean)], by = "tss_id", allow.cartesian = TRUE)
    rbind(v[, .(mean_fe = mean(fe, na.rm = TRUE), n_molecules = sum(!is.na(fe))), by = .(cluster, relative_start)],
          v[, .(cluster = "all", mean_fe = mean(fe, na.rm = TRUE), n_molecules = sum(!is.na(fe))), by = relative_start])[
      , `:=`(mark = mark, relative_end = relative_start + bin - 1L)]
  }))[, .(mark, cluster, relative_start, relative_end, mean_fe, n_molecules)]
}

###########################
# Show histone profiles
############################
# 1. Input: Histone profiles, cluster labels and a stacked-layout choice.
# 2. Output: An overlay or vertically stacked histone figure.
#
# 3. Steps:
#    - Omit the all-molecule reference
#    - Preserve numeric cluster order and colors
#    - Show one mark per column with shared scales across its cluster rows.
pol2_plot_histone_profiles <- function(profiles, clusters, stacked = FALSE) {
  clusters <- pol2_numeric_clusters(clusters)
  d <- copy(profiles)[cluster != "all"]
  d[, `:=`(mark = factor(mark, levels = unique(mark)),
             cluster = factor(cluster, levels = clusters),
             x = (relative_start + relative_end) / 2)]
  if (stacked) {
    # One mark per column; fixed scales within each column compare cluster rows.
    columns <- lapply(levels(d$mark), function(mark_name) {
      ggplot(d[mark == mark_name], aes(x, mean_fe, colour = cluster)) +
        annotate("rect", xmin = -1000, xmax = 999, ymin = -Inf, ymax = Inf, fill = "grey93") +
        geom_vline(xintercept = 0, colour = "grey60", linetype = "dotted", linewidth = .3) +
        geom_line(linewidth = .5) +
        facet_grid(cluster ~ ., drop = FALSE) +
        scale_colour_manual(values = cluster_id_palette(clusters), guide = "none") +
        labs(x = "Position relative to TSS (bp)", y = "Mean fold enrichment over control", title = mark_name) +
        theme_fiberseq("bw", base_size = 10)
    })
    return(cowplot::plot_grid(plotlist = columns, nrow = 1, align = "h", axis = "tb"))
  }
  ggplot(d, aes(x, mean_fe, colour = cluster)) +
    annotate("rect", xmin = -1000, xmax = 999, ymin = -Inf, ymax = Inf, fill = "grey93") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dotted", linewidth = .3) +
    geom_line(linewidth = .5) +
    facet_wrap(~mark, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = cluster_id_palette(clusters)) +
    labs(x = POL2_X_LABEL, y = "Mean fold enrichment over control", colour = "Cluster",
         title = "GM12878 histone ChIP around the TSS",
         subtitle = "Mean over each cluster's molecules; 50-bp bins; shaded: clustering window") +
    theme_fiberseq("bw", base_size = 10)
}
