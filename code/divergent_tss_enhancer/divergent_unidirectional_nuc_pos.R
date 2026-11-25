# Nucleosome positioning of the divergent vs unidirectional TSS clusters
# (code/divergent_tss_enhancer/cluster_divergent_unidirectional_tss.R).
# Source without running; Rscript .../divergent_unidirectional_nuc_pos.R
# [--pilot] executes nuc_main().
ROOT <- '/project/spott/cshan/fiber-seq'
source(file.path(ROOT, 'code/nucleosome_positioning/process_nuc_pos.R'))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

########
# Settings: inputs, output folder and bootstrap draws
########

# 1. Read the clustering result (reads, clusters, windows) of the TSS clustering;
#    its orientation is reused (minus-strand TSSs flipped when the clustering
#    was oriented, as for any pooled promoters).
# 2. Read the 31-sample LCL table (fire_dir locates the fibertools nuc_by_chr files).
# 3. Write tables/ and plots/ under nucleosome_positioning.
# 4. bootstrap: crossed locus x sample draws for the centre-occupancy interval.
nuc_config <- function() list(
  clustering = file.path(ROOT, 'LCL_project/divergent_tss_enhancer/clustering_results/tables/clustering.rds'),
  samples = '/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv',
  output = file.path(ROOT, 'LCL_project/divergent_tss_enhancer/nucleosome_positioning'),
  bootstrap = 500L, seed = 1L)

########
# One anchor row per clustered window, matching the clustered bases exactly
########

# 1. Anchor = the TSS / pair-centre base (offset 0 in the clustering,
#    1-based centre + 1).
# 2. Strand = the TSS strand when the clustering was oriented ("+" for
#    divergent pairs), otherwise "+" for every window.
# 3. Relative window = the clustered bases relative to the anchor, so the
#    anchors cover exactly the clustered window (build_region_anchors(), one
#    call per distinct relative window).
# 4. Check every anchor window equals its clustered window.
nuc_anchors <- function(regions, oriented) {
  r <- as.data.table(regions)
  r[, strand := if (oriented) tss_strand else '+']
  r[, anchor := center + 1L]
  r[, c('rel_lo', 'rel_hi') := rel_interval(analysis_start, analysis_end, anchor, strand)]
  table <- r[, .(region_id, chr, anchor, strand, region_type = class, gene = region_id,
                 ensg = NA_character_, rel_lo, rel_hi)]
  a <- rbindlist(lapply(split(table, by = c('rel_lo', 'rel_hi')), function(x)
    build_region_anchors(NULL, NULL, NULL, relative_window = c(x$rel_lo[1], x$rel_hi[1]),
      anchor_table = x[, !c('rel_lo', 'rel_hi')])))
  a <- a[match(r$region_id, region_id)]
  stopifnot(all(a$view_start == r$analysis_start), all(a$view_end == r$analysis_end))
  a
}

########
# Reads with their spans relative to the anchor, clipped to the window
########

# 1. Keep each clustered read's window, sample, raw read name, class, cluster
#    and alignment span (build_read_table()).
# 2. Clip the span to the window; every read spans the whole window.
nuc_reads <- function(assignments, anchors) {
  a <- as.data.table(assignments)
  reads <- build_read_table(a[, .(region_id, RID, sample_name, original_RID, class, cluster,
    read_start = start, read_end = end)], anchors)
  reads[anchors, on = 'region_id', `:=`(rel_read_start = pmax(rel_read_start, i.rel_view_start),
    rel_read_end = pmin(rel_read_end, i.rel_view_end))]
  stopifnot(nrow(reads) == nrow(a),
    all(reads$rel_read_start == anchors$rel_view_start[match(reads$region_id, anchors$region_id)]),
    all(reads$rel_read_end == anchors$rel_view_end[match(reads$region_id, anchors$region_id)]))
  reads
}

########
# Per-read nucleosome positions, occupancy profiles and summaries
########

# 1. Per read: midpoints of the nearest nucleosome left (< 0) and right (>= 0)
#    of the anchor, the gap between them, and whether one covers the anchor
#    (read_position_summary()).
# 2. Medians per class x cluster and per class (all clusters), with a 95%
#    crossed locus x sample bootstrap interval for centre occupancy
#    (position_summary_by_group()).
# 3. Per 10-bp position: fraction of reads with a nucleosome covering it, per
#    window x class x cluster (occupancy_profile()), pooled per class x cluster
#    and per class.
nuc_summaries <- function(footprints, reads, anchors, cfg) {
  pos <- read_position_summary(footprints, reads)
  by_cluster <- position_summary_by_group(pos, by = c('class', 'cluster'),
    bootstrap = cfg$bootstrap, seed = cfg$seed, sample_col = 'sample_name')
  by_class <- position_summary_by_group(pos, by = 'class',
    bootstrap = cfg$bootstrap, seed = cfg$seed, sample_col = 'sample_name')
  summary <- rbind(by_class[, cluster := 'all'], by_cluster, use.names = TRUE)
  occ <- occupancy_profile(footprints, reads, anchors, by = c('region_id', 'class', 'cluster'))
  profile <- rbind(occ[, .(cluster = 'all', n_occupied = sum(n_occupied), n_reads = sum(n_reads)), by = .(class, pos)],
    occ[, .(n_occupied = sum(n_occupied), n_reads = sum(n_reads)), by = .(class, cluster = as.character(cluster), pos)])
  profile[, fraction := n_occupied / n_reads]
  levels <- c('all', levels(reads$cluster))
  summary[, cluster := factor(as.character(cluster), levels = levels)]
  profile[, cluster := factor(cluster, levels = levels)]
  list(molecules = pos, summary = summary, occupancy = profile)
}

########
# Draw and save the nucleosome-positioning figures as PDFs
########

# 1. nucleosome_occupancy_by_cluster.pdf: occupancy along the window, one
#    panel for all reads and one per cluster, divergent vs unidirectional.
# 2. nucleosome_positions_by_cluster.pdf: per-read midpoints of the nearest
#    nucleosome left / right of the anchor and the gap between them.
# 3. centre_occupancy_by_cluster.pdf: fraction of reads with a nucleosome on
#    the anchor base, with the bootstrap interval.
nuc_plots <- function(s, reads, output, oriented) {
  class_colors <- fiberseq_category_palette(reads$class)
  x_label <- if (oriented)
    'Position relative to TSS / pair centre (bp); minus-strand TSSs flipped' else
    'Position relative to window centre (bp), genome orientation'
  n <- rbind(reads[, .(cluster = 'all', n = .N), by = class], reads[, .(n = .N), by = .(class, cluster)])
  n <- n[, .(label = paste0(cluster[1], ' (', paste0(class, ' n=', n, collapse = ', '), ')')), by = cluster]
  labels <- as_labeller(setNames(n$label, n$cluster))
  p1 <- plot_group_profile(s$occupancy, group_col = 'class', style = 'line', colors = class_colors,
      layers = list(list(geom = 'vline', xintercept = 0, linetype = 'dashed', colour = 'grey30', linewidth = .4),
        list(geom = 'line', linewidth = .6))) +
    facet_wrap(~cluster, ncol = 1, labeller = labels) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1), expand = expansion(mult = c(0, .02))) +
    labs(x = x_label, y = 'Fraction of reads with a nucleosome', colour = 'TSS class',
      title = 'Nucleosome occupancy (fibertools nuc calls, all sizes)') +
    cowplot::theme_cowplot(font_size = 10) + cowplot::panel_border() + theme(legend.position = 'bottom')
  ggsave(file.path(output, 'nucleosome_occupancy_by_cluster.pdf'), p1, width = 8, height = 1.4 * nrow(n) + 1.5)
  long <- melt(rbind(copy(s$molecules)[, cluster := 'all'], copy(s$molecules)[, cluster := as.character(cluster)]),
    id.vars = c('RID', 'class', 'cluster'), measure.vars = c('minus1_mid', 'plus1_mid', 'nfr_width'),
    variable.name = 'metric', value.name = 'value')[!is.na(value)]
  long[, cluster := factor(cluster, levels = levels(s$summary$cluster))]
  p2 <- ggplot(long, aes(class, value, fill = class)) +
    geom_boxplot(outlier.shape = NA, width = .6, alpha = .8, linewidth = .3) +
    facet_grid(metric ~ cluster, scales = 'free_y', labeller = labeller(metric = c(
      minus1_mid = 'nearest nucleosome < 0\nmidpoint (bp)', plus1_mid = 'nearest nucleosome >= 0\nmidpoint (bp)',
      nfr_width = 'gap between them (bp)'))) +
    scale_fill_manual(values = class_colors, guide = 'none') +
    labs(x = NULL, y = NULL, title = 'Nucleosomes flanking the window centre, per read',
      subtitle = sub('^Position', 'Positions', x_label)) +
    cowplot::theme_cowplot(font_size = 9) + cowplot::panel_border() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(output, 'nucleosome_positions_by_cluster.pdf'), p2, width = 2 + 1.6 * nlevels(long$cluster), height = 7)
  p3 <- ggplot(s$summary, aes(cluster, frac_anchor_occupied, colour = class)) +
    geom_point(size = 2, position = position_dodge(.5)) +
    geom_errorbar(aes(ymin = tss_lower, ymax = tss_upper), width = .2, position = position_dodge(.5)) +
    scale_colour_manual(values = class_colors) + scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = 'Fraction with a nucleosome on the centre base', colour = 'TSS class',
      caption = '95% crossed locus x sample bootstrap intervals') +
    cowplot::theme_cowplot(font_size = 10) + cowplot::panel_border()
  ggsave(file.path(output, 'centre_occupancy_by_cluster.pdf'), p3, width = 7, height = 4.5)
}

########
# Run the nucleosome-positioning analysis: anchors, reads, nucleosomes, save, plot
########

# 1. Load the clustering result and the sample table; n_loci keeps only the
#    first n windows per class (for testing).
# 2. Build the anchors in the clustering's orientation (nuc_anchors()) and
#    the read table (nuc_reads()).
# 3. Read every read's fibertools nucleosomes in its window
#    (extract_region_footprints(format = "bed12_fibertools")) and attach class
#    and cluster.
# 4. Summarise (nuc_summaries()), save the tables in tables/ and the figures
#    in plots/.
nuc_main <- function(cfg = nuc_config(), n_loci = NULL) {
  data.table::setDTthreads(1L)
  result <- readRDS(cfg$clustering)
  assignments <- result$assignments
  regions <- result$regions[result$regions$region_id %in% assignments$region_id, ]
  if (!is.null(n_loci)) {
    regions <- dplyr::bind_rows(lapply(split(regions, regions$class), head, n_loci))
    assignments <- assignments[assignments$region_id %in% regions$region_id, ]
  }
  samples <- read.csv(cfg$samples, stringsAsFactors = FALSE)
  oriented <- isTRUE(result$params$orient_strand)
  anchors <- nuc_anchors(regions, oriented)
  reads <- nuc_reads(assignments, anchors)
  footprints <- extract_region_footprints(reads, anchors, include_tf = FALSE,
    original_id_col = 'original_RID', format = 'bed12_fibertools', sample_table = samples)
  footprints <- merge(footprints, reads[, .(region_id, RID, class, cluster)], by = c('region_id', 'RID'))
  s <- nuc_summaries(footprints, reads, anchors, cfg)
  tables <- file.path(cfg$output, 'tables'); plots <- file.path(cfg$output, 'plots')
  dir.create(tables, recursive = TRUE, showWarnings = FALSE); dir.create(plots, showWarnings = FALSE)
  fwrite(footprints, file.path(tables, 'nucleosome_footprints.tsv.gz'), sep = '\t')
  fwrite(s$molecules, file.path(tables, 'nucleosome_molecules.tsv.gz'), sep = '\t')
  fwrite(s$summary, file.path(tables, 'nucleosome_position_summary.tsv'), sep = '\t')
  fwrite(s$occupancy, file.path(tables, 'nucleosome_occupancy_profiles.tsv'), sep = '\t')
  saveRDS(list(config = cfg, anchors = anchors, n_loci = n_loci, oriented = oriented,
    source = 'fibertools ft extract nuc_by_chr; longest alignment per read; all sizes',
    clustering_md5 = unname(tools::md5sum(cfg$clustering))), file.path(tables, 'nucleosome_provenance.rds'))
  nuc_plots(s, reads, plots, oriented)
  message('COMPLETE: ', cfg$output)
  invisible(s)
}
if (sys.nframe() == 0L) {
  if ('--pilot' %in% commandArgs(trailingOnly = TRUE)) {
    cfg <- nuc_config(); cfg$output <- tempfile('nuc-pilot-'); cfg$bootstrap <- 20L
    nuc_main(cfg, n_loci = 20L); unlink(cfg$output, recursive = TRUE)
    message('PILOT PASSED; temporary outputs removed')
  } else nuc_main()
}
