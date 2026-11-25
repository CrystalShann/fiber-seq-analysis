# Autocorrelation of TSS-relative nucleosome profiles, upstream and downstream,
# for the divergent vs unidirectional TSS clusters
# (code/divergent_tss_enhancer/cluster_divergent_unidirectional_tss.R).
# Source without running; Rscript .../divergent_unidirectional_autocorrelation.R
# [--pilot] executes acf_main().
ROOT <- '/project/spott/cshan/fiber-seq'
ACF_DIR <- file.path(ROOT, 'code/clustering_methods/auto_correlation')
source(file.path(ACF_DIR, '02_compute_autocorrelations.R'))
source(file.path(ACF_DIR, '04_summarize_autocorrelations.R'))
source(file.path(ROOT, 'code/clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r'))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

########
# Settings: input, output folder and the two TSS-relative windows
########

# 1. Read the clustering result: its per-read nucleosome profiles (fibertools
#    nuc calls, 1 = base covered by a nucleosome) are TSS-relative, -500..+499,
#    with minus-strand TSSs already flipped.
# 2. Write tables/ and plots/ under autocorrelation.
# 3. Windows as in the Pol II report: upstream -500..-25 and downstream
#    +25..+499 bp, leaving out the 50 bp around the TSS. For divergent pairs
#    upstream = minus-TSS side and downstream = plus-TSS side of the centre.
acf_config <- function() list(
  clustering = file.path(ROOT, 'LCL_project/divergent_tss_enhancer/clustering_results/tables/clustering.rds'),
  output = file.path(ROOT, 'LCL_project/divergent_tss_enhancer/autocorrelation'),
  windows = list(upstream = c(-500L, -25L), downstream = c(25L, 499L)))

########
# Per-read autocorrelation of the nucleosome profile in one window
########

# 1. Take every read's 0/1 nucleosome profile (rows = clustered reads).
# 2. Autocorrelation of the profile inside the window, every lag
#    (autocorrelations()); reads with no nucleosome base, or nucleosome over
#    the whole window, have none (status "zero_variance"; "clustered" = has
#    an ACF, the shared functions' name).
# 3. Per read: its nucleosome fraction and m6A calls per bp in the window,
#    and its repeat peak (highest positive local ACF peak at 140-250 bp,
#    repeat_peak()).
acf_window <- function(result, assignments, name, bounds) {
  a <- assignments
  nuc <- result$annotations$nucleosomes[a$RID, , drop = FALSE]
  res <- autocorrelations(nuc, positions = -500:499, window = bounds)
  cols <- as.character(seq.int(bounds[1], bounds[2]))
  records <- data.frame(RID = a$RID, sample_name = a$sample_name, region_id = a$region_id,
    class = a$class, cluster = as.integer(sub('^cluster', '', a$cluster)), window = name,
    status = ifelse(res$valid, 'clustered', 'zero_variance'),
    nucleosome_fraction = rowMeans(nuc[, cols, drop = FALSE]),
    m6a_call_fraction = Matrix::rowMeans(result$site_met_mat[a$RID, cols, drop = FALSE]))
  peaks <- t(vapply(seq_len(nrow(res$profiles)), function(i) repeat_peak(res$profiles[i, ]), numeric(2)))
  records$repeat_peak_bp <- peaks[, 1]; records$peak_acf <- peaks[, 2]
  list(profiles = res$profiles, records = records)
}

########
# Mean ACF curve and its repeat peak per class x cluster (and per class)
########

# 1. For each class, average the per-read ACFs of each cluster
#    (summarize_autocorrelations(), reads with an ACF only).
# 2. Do the same with all clusters of the class pooled (cluster "all").
# 3. Return the mean curves (mean, SD, SEM per lag) and the repeat peak of
#    each mean curve.
acf_summaries <- function(w) {
  out <- list()
  for (cl in c('divergent', 'unidirectional')) for (pooled in c(TRUE, FALSE)) {
    keep <- w$records$class == cl
    rec <- w$records[keep, ]; rec$region_id <- cl   # rows are labelled by region_id
    if (pooled) rec$cluster <- 0L                   # 0 = all clusters of the class
    out[[length(out) + 1L]] <- summarize_autocorrelations(w$profiles[keep, , drop = FALSE], rec)
  }
  tidy <- function(x) {
    x <- as.data.table(x); setnames(x, 'region_id', 'class')
    x[, cluster := ifelse(cluster == '0', 'all', paste0('cluster', cluster))]
    x[, window := w$records$window[1]][]
  }
  list(curves = rbindlist(lapply(out, function(x) tidy(x$averages))),
       peaks = rbindlist(lapply(out, function(x) tidy(x$summaries))))
}

########
# Draw and save the autocorrelation figures as PDFs
########

# 1. nucleosome_acf_curves.pdf: mean ACF (+/- SEM) against lag, upstream and
#    downstream, one row for all reads and one per cluster, divergent vs
#    unidirectional; points mark each curve's repeat peak.
# 2. nucleosome_acf_repeat_peaks.pdf: per-read repeat-peak lag, same layout.
acf_plots <- function(curves, peaks, molecules, output) {
  class_colors <- fiberseq_category_palette(molecules$class)
  levels <- c('all', sort(unique(curves$cluster[curves$cluster != 'all'])))
  window_labels <- c(upstream = 'upstream (-500 to -25 bp; divergent: minus-TSS side)',
                     downstream = 'downstream (+25 to +499 bp; divergent: plus-TSS side)')
  curves[, `:=`(cluster = factor(cluster, levels = levels), window = factor(window, levels = names(window_labels)))]
  peaks[, `:=`(cluster = factor(cluster, levels = levels), window = factor(window, levels = names(window_labels)))]
  p1 <- ggplot(curves[lag_bp > 0], aes(lag_bp, mean_acf, colour = class, fill = class)) +
    geom_hline(yintercept = 0, colour = 'grey60', linewidth = .3) +
    geom_ribbon(aes(ymin = mean_acf - sem_acf, ymax = mean_acf + sem_acf), colour = NA, alpha = .25) +
    geom_line(linewidth = .5) +
    geom_point(data = peaks[!is.na(positive_local_peak_140_250_bp)],
      aes(positive_local_peak_140_250_bp, peak_acf), size = 1.6) +
    facet_grid(cluster ~ window, labeller = labeller(window = window_labels)) +
    scale_colour_manual(values = class_colors) + scale_fill_manual(values = class_colors) +
    labs(x = 'Lag (bp)', y = 'Mean autocorrelation of nucleosome profile', colour = 'TSS class', fill = 'TSS class',
      title = 'Nucleosome-profile autocorrelation (fibertools nuc calls; minus-strand TSSs flipped)',
      caption = 'Points: repeat peak of each mean curve (highest positive local peak, 140-250 bp)') +
    cowplot::theme_cowplot(font_size = 10) + cowplot::panel_border() + theme(legend.position = 'bottom')
  ggsave(file.path(output, 'nucleosome_acf_curves.pdf'), p1, width = 11, height = 1.8 * length(levels) + 2)
  m <- as.data.table(molecules)[!is.na(repeat_peak_bp)][, cluster := paste0('cluster', cluster)]
  m <- rbind(copy(m)[, cluster := 'all'], m)
  m[, `:=`(cluster = factor(cluster, levels = levels), window = factor(window, levels = names(window_labels)))]
  n <- m[, .N, by = .(class, cluster, window)]
  p2 <- ggplot(m, aes(class, repeat_peak_bp, fill = class)) +
    geom_boxplot(outlier.shape = NA, width = .6, alpha = .8, linewidth = .3) +
    geom_text(data = n, aes(class, 255, label = paste0('n=', N)), size = 2.6, inherit.aes = FALSE) +
    facet_grid(window ~ cluster, labeller = labeller(window = c(upstream = 'upstream', downstream = 'downstream'))) +
    scale_fill_manual(values = class_colors, guide = 'none') +
    scale_y_continuous(limits = c(140, 258)) +
    labs(x = NULL, y = 'Per-read repeat peak (bp)', title = 'Per-read nucleosome repeat peak',
      subtitle = 'Reads whose ACF has a positive local peak at 140-250 bp') +
    cowplot::theme_cowplot(font_size = 9) + cowplot::panel_border() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(output, 'nucleosome_acf_repeat_peaks.pdf'), p2, width = 2 + 1.6 * length(levels), height = 5.5)
}

########
# Run the autocorrelation analysis: per-read ACFs, summaries, save, plot
########

# 1. Load the clustering result; it must be strand-oriented (promoters are
#    pooled) and carry the per-read nucleosome profiles. n_loci keeps only the
#    first n windows per class (for testing).
# 2. Per window (upstream, downstream): per-read ACFs (acf_window()) and
#    class x cluster summaries (acf_summaries()).
# 3. Save the per-read table, mean curves, repeat peaks and per-read ACF
#    matrices in tables/, the figures in plots/.
acf_main <- function(cfg = acf_config(), n_loci = NULL) {
  result <- readRDS(cfg$clustering)
  stopifnot(isTRUE(result$params$orient_strand), !is.null(result$site_met_mat),
    identical(colnames(result$annotations$nucleosomes), as.character(-500:499)))
  a <- result$assignments
  if (!is.null(n_loci)) {
    keep <- unlist(lapply(split(unique(a[, c('region_id', 'class')]), ~class), function(x) head(x$region_id, n_loci)))
    a <- a[a$region_id %in% keep, ]
  }
  w <- lapply(names(cfg$windows), function(name) acf_window(result, a, name, cfg$windows[[name]]))
  names(w) <- names(cfg$windows)
  s <- lapply(w, acf_summaries)
  curves <- rbindlist(lapply(s, `[[`, 'curves')); peaks <- rbindlist(lapply(s, `[[`, 'peaks'))
  molecules <- rbindlist(lapply(w, `[[`, 'records'))
  tables <- file.path(cfg$output, 'tables'); plots <- file.path(cfg$output, 'plots')
  dir.create(tables, recursive = TRUE, showWarnings = FALSE); dir.create(plots, showWarnings = FALSE)
  fwrite(molecules, file.path(tables, 'nucleosome_acf_molecules.tsv.gz'), sep = '\t')
  fwrite(curves, file.path(tables, 'nucleosome_acf_cluster_curves.tsv.gz'), sep = '\t')
  fwrite(peaks, file.path(tables, 'nucleosome_acf_cluster_peaks.tsv'), sep = '\t')
  for (name in names(w)) saveRDS(list(acf = w[[name]]$profiles, lag_bp = seq_len(ncol(w[[name]]$profiles)) - 1L,
    relative_bp = seq.int(cfg$windows[[name]][1], cfg$windows[[name]][2]), RID = a$RID, window = name),
    file.path(tables, paste0('nucleosome_acf_per_molecule_', name, '.rds')))
  saveRDS(list(config = cfg, n_loci = n_loci,
    method = 'autocorrelations(), summarize_autocorrelations(), repeat_peak() on 1-bp fibertools nucleosome profiles',
    clustering_md5 = unname(tools::md5sum(cfg$clustering))), file.path(tables, 'nucleosome_acf_provenance.rds'))
  acf_plots(curves, peaks, molecules, plots)
  message('COMPLETE: ', cfg$output)
  invisible(list(curves = curves, peaks = peaks, molecules = molecules))
}
if (sys.nframe() == 0L) {
  if ('--pilot' %in% commandArgs(trailingOnly = TRUE)) {
    cfg <- acf_config(); cfg$output <- tempfile('acf-pilot-')
    acf_main(cfg, n_loci = 20L); unlink(cfg$output, recursive = TRUE)
    message('PILOT PASSED; temporary outputs removed')
  } else acf_main()
}
