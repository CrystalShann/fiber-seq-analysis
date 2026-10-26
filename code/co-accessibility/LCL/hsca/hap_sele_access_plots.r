# Plotting functions and supporting helpers for LCL_hsca.Rmd.
# Sourced into the report environment; shared plotting/parsing dependencies
# are loaded by the report before these functions are called.
# Report data, params, and output directories are supplied by that environment.

pair_imprinted = function(dt) dt$CRE1 %in% imprinted_cres | dt$CRE2 %in% imprinted_cres

# phased-pair-plot-function
plot_phased_pair_tracks = function(result, records, pair, group_colors = NULL,
    sample_colors = NULL, cre_colors = NULL, highlight_color = "#FFFF00",
    highlight_alpha = .22, title = NULL, min_read_height = 3.5, row_height = .035,
    group_column = "allele_group") {
  a = result$assignments
  region = result$region
  stopifnot(group_column %in% names(a), length(pair$cre1_start) == 1L,
    length(pair$cre2_start) == 1L, length(region$focal_pos) == 1L,
    region$focal_pos >= region$analysis_start, region$focal_pos <= region$analysis_end)
  groups = a[[group_column]]
  a$cluster = if (is.factor(groups)) droplevels(groups) else factor(groups, levels = unique(groups))
  result$assignments = a # Local renderer field only; no clustering is performed.
  stopifnot(nrow(a) > 0L, !anyDuplicated(a$RID), is.factor(a$cluster), !anyNA(a$cluster),
    all(a$start <= region$analysis_start), all(a$end >= region$analysis_end),
    length(highlight_alpha) == 1L, is.finite(highlight_alpha), highlight_alpha >= 0,
    highlight_alpha <= 1, min_read_height > 0, row_height > 0)
  anchor = plot_anchor(region)
  if (is.null(sample_colors)) sample_colors = sample_palette(a$sample_name)
  if (is.null(group_colors)) group_colors = pair_allele_colors(region)
  stopifnot(all(levels(a$cluster) %in% names(group_colors)),
    all(sub("_.*$", "", a$sample_name) %in% names(sample_colors)))
  allele_colors = group_colors[levels(a$cluster)]
  # Give this private renderer copy the requested group palette; leave the shared renderer intact.
  renderer = plot_smf_reads
  renderer_env = new.env(parent = environment(renderer))
  renderer_env$cluster_id_palette = function(groups) allele_colors[groups]
  environment(renderer) = renderer_env
  bounds = c(anchor$left - .085 * (anchor$right - anchor$left), anchor$right + .5)
  cre_classes = sub(".*\\.", "", c(pair$CRE1, pair$CRE2))
  if (is.null(cre_colors)) {
    cre_colors = unname(CCRE_COLORS[cre_classes])
    cre_colors[is.na(cre_colors)] = "grey40"
  }
  stopifnot(length(cre_colors) == 2L, !anyNA(cre_colors))
  cre = data.frame(
    left = c(pair$cre1_start, pair$cre2_start) + 1L - anchor$anchor,
    right = c(pair$cre1_end, pair$cre2_end) - anchor$anchor,
    label = paste0(c("CRE1 (", "CRE2 ("), cre_classes, ")"), color = cre_colors)
  highlight = function() geom_rect(data = cre,
    aes(xmin = left, xmax = right, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = highlight_color, alpha = highlight_alpha)
  annotation = ggplot(cre) + highlight() +
    geom_rect(aes(xmin = left, xmax = right, ymin = .3, ymax = .6, fill = color)) +
    geom_text(aes(x = (left + right) / 2, y = .8, label = label), size = 3) +
    geom_vline(xintercept = 0, color = "#984EA3", linetype = "dashed", linewidth = .4) +
    annotate("point", x = 0, y = .1, shape = 18, size = 3, color = "#984EA3") +
    scale_fill_identity() + scale_x_continuous(limits = bounds, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = if (is.null(title)) paste0(region$annotation, " | full-window, phased molecules") else title,
      subtitle = sprintf("%s:%s %s>%s | dashed line at SNP (0 bp) | %d molecules",
        region$chr, format(region$focal_pos, big.mark = ",", scientific = FALSE),
        region$ref, region$alt, nrow(a)), x = NULL, y = NULL) +
    cowplot::theme_cowplot(font_size = 10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
  draw = function(tracks, ylabel) {
    # FALSE omits sample-local HP1/HP2; the existing group and sample bars remain.
    p = renderer(result, records, sample_colors, tracks = tracks, include_haplotype = FALSE)
    legend_colors = c(feature_colors(tracks), allele_colors, sample_colors)
    legend_labels = c(unname(LCL_TRACK_LABELS[tracks]), names(allele_colors), names(sample_colors))
    p = suppressMessages(p + scale_fill_identity(guide = "legend", name = NULL,
      breaks = unname(legend_colors), labels = legend_labels))
    # Put yellow behind the read marks so their original colors remain intact.
    p$layers = c(list(highlight()), p$layers)
    p + geom_vline(xintercept = 0, color = "#984EA3", linetype = "dashed", linewidth = .4) +
      geom_vline(xintercept = c(cre$left, cre$right), color = "grey60",
        linetype = "dotted", linewidth = .25) +
      labs(title = NULL, subtitle = NULL, y = ylabel,
        caption = "Left bars: SNP allele, sample. Shaded bands: tested cCREs.")
  }
  m6a = draw("m6A", "m6A: one row per molecule")
  footprints = draw(LCL_FOOTPRINT_TRACKS, "Footprints: same molecule order")
  read_height = max(min_read_height, row_height * nrow(a) + 1.8)
  combined = cowplot::plot_grid(annotation, m6a, footprints, ncol = 1,
    align = "v", axis = "lr", rel_heights = c(1.1, read_height, read_height))
  list(panels = list(combined = combined, m6a = m6a, footprints = footprints),
    height = 1.1 + 2 * read_height)
}

# pair-allele-acf-functions
pair_allele_colors = function(region) {
  setNames(c('#FDBE85', '#9ECAE1'),
    c(paste0('REF-', region$ref), paste0('ALT-', region$alt)))
}

annotate_pair_alleles = function(result, audit) {
  a = result$assignments
  audit = as.data.frame(audit)
  keys = paste(audit$sample_name, audit$original_RID)
  stopifnot(!anyDuplicated(keys))
  matched = match(paste(a$sample_name, a$original_RID), keys)
  status = audit$allele_status[matched]
  genotype = audit$focal_genotype[matched]
  bases = sub('.*: ', '', audit$allele_label[matched])
  known = !is.na(status) & status == 'phased_focal_genotype' &
    genotype %in% c('0|1', '1|0') & bases %in% c(result$region$ref, result$region$alt)
  allele = rep('Unassigned', nrow(a))
  allele[known] = ifelse(bases[known] == result$region$ref,
    paste0('REF-', result$region$ref), paste0('ALT-', result$region$alt))
  a$allele_group = factor(allele, levels = c(names(pair_allele_colors(result$region)), "Unassigned"))
  a$focal_genotype = genotype
  a$allele_status = ifelse(is.na(status), 'no_phasing_record', status)
  a$allele_assigned = known
  result$assignments = a
  result
}

# Filter before computing either distance matrix or graph. The audit retains exclusions.
cluster_phased_pair = function(input, audit, config) {
  input = annotate_pair_alleles(input, audit)
  keep = input$assignments$allele_assigned
  a = droplevels(input$assignments[keep, , drop = FALSE])
  if (nrow(a) < 3L) stop("At least three REF/ALT-assigned full-window reads are required for clustering.")
  mat = input$site_met_mat[a$RID, , drop = FALSE]
  region = input$region
  result = leiden_manhattan_cluster(mat, a, region$analysis_start, region$analysis_end,
    window_size = 0L, k_neighbors = config$leiden_k_neighbors,
    resolution = config$leiden_resolution, seed = config$leiden_seed)
  meta = a[match(result$assignments$RID, a$RID), , drop = FALSE]
  for (key in setdiff(names(meta), names(result$assignments))) result$assignments[[key]] = meta[[key]]
  result$region = region
  result$site_met_mat = mat[result$assignments$RID, , drop = FALSE]
  result$params$read_filter = "full-window; phased heterozygous focal SNP; REF/ALT assigned"
  result$params$n_input_full_window = nrow(input$assignments)
  result$params$n_unassigned_excluded = sum(!keep)
  stopifnot(all(result$assignments$allele_assigned),
    !anyNA(result$assignments$allele_group),
    !any(result$assignments$allele_group == "Unassigned"))
  result
}

pair_read_annotations = function(assignments, region, include_cluster = TRUE) {
  stopifnot(all(as.character(assignments$allele_group) %in% names(pair_allele_colors(region))))
  samples = sample_palette(assignments$sample_name)
  df = data.frame(sample = factor(sub('_.*$', '', assignments$sample_name), levels = names(samples)),
    allele = factor(assignments$allele_group, levels = names(pair_allele_colors(region))))
  colors = list(sample = samples, allele = pair_allele_colors(region))
  if (include_cluster) {
    df = data.frame(cluster = assignments$cluster, df)
    colors = c(list(cluster = cluster_id_palette(levels(assignments$cluster))), colors)
  }
  stopifnot(!anyNA(df))
  ComplexHeatmap::rowAnnotation(df = df, col = colors,
    annotation_name_gp = grid::gpar(fontsize = 8),
    annotation_legend_param = list(allele = list(title = 'Focal SNP allele')))
}

plot_pair_allele_graph = function(result, main = NULL) {
  # Reuse the shared graph layout. Its color_by='allele' branch requires every
  # read to have an rsID-based allele label, so adapt only the point colors here.
  p = plot_knn_graph(result, main = main)
  point = which(vapply(p$layers, function(layer) inherits(layer$geom, 'GeomPoint'), logical(1)))
  stopifnot(length(point) == 1L)
  nodes = p$layers[[point]]$data
  nodes$color_group = result$assignments$allele_group[match(nodes$RID, result$assignments$RID)]
  p$layers[[point]]$data = nodes
  colors = pair_allele_colors(result$region)
  counts = table(factor(nodes$color_group, levels = names(colors)))
  suppressMessages(p + scale_color_manual(values = colors, breaks = names(colors),
    labels = paste0(names(colors), ' (n=', as.integer(counts), ')'), drop = FALSE) +
    labs(color = 'Focal SNP allele'))
}

pair_cluster_alleles = function(result, title) {
  counts = as.data.table(as.data.frame(table(cluster = result$assignments$cluster,
    allele = result$assignments$allele_group)))
  setnames(counts, 'Freq', 'n_reads')
  counts[, total_reads := sum(n_reads), by = cluster]
  counts[, fraction := n_reads / total_reads]
  p = ggplot(counts, aes(cluster, fraction, fill = allele)) + geom_col() +
    scale_fill_manual(values = pair_allele_colors(result$region), drop = FALSE) +
    scale_y_continuous(labels = scales::percent) + cowplot::theme_cowplot(font_size = 10) +
    labs(title = title, x = 'Cluster', y = 'Fraction of molecules in cluster', fill = 'Focal SNP allele')
  list(counts = counts, plot = p)
}

pair_autocorrelation = function(result, config) {
  # Call the existing ACF calculation and PCA/correlation-neighbor Leiden code
  # directly, with REF/ALT-assigned full-window molecules and no allele features.
  work = tempfile('pair_acf_')
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  a = result$assignments
  stopifnot(all(a$allele_assigned), !anyNA(a$allele_group),
    !any(a$allele_group == 'Unassigned'))
  region = result$region
  binary = matrix(0L, nrow(a), region$width, dimnames = list(a$RID, NULL))
  binary[, as.integer(colnames(result$site_met_mat)) - region$analysis_start + 1L] =
    as.matrix(result$site_met_mat[a$RID, , drop = FALSE])
  fwrite(as.data.table(binary), file.path(work, 'binary.tsv'), sep = '\t', col.names = FALSE)
  runner = r"---(
import importlib.util, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
source, work = map(Path, sys.argv[1:3])
def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, source / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
compute = load('pair_acf_compute', '03_compute_autocorrelations.py')
cluster = load('pair_acf_cluster', '04_cluster_autocorrelations.py')
binary = np.loadtxt(work / 'binary.tsv', delimiter='\t', ndmin=2)
profiles, valid = compute.autocorrelations(binary)
labels, status, embedding, info, edges = cluster.cluster_profiles(
    profiles, valid, n_pcs=int(sys.argv[3]), n_neighbors=int(sys.argv[4]),
    resolution=float(sys.argv[5]), seed=int(sys.argv[6]))
np.savetxt(work / 'acf.tsv', profiles, delimiter='\t', fmt='%.17g')
pd.DataFrame(dict(row_index=np.arange(len(labels)), cluster=labels, status=status,
                  valid=valid, umap1=embedding[:,0], umap2=embedding[:,1])).to_csv(
    work / 'assignments.tsv', sep='\t', index=False)
pd.DataFrame(edges, columns=['source_index','target_index','weight']).to_csv(
    work / 'edges.tsv', sep='\t', index=False)
with open(work / 'parameters.json', 'w') as f:
    json.dump(info, f)
)---"
  script = file.path(work, 'run.py')
  writeLines(runner, script)
  args = c(script, file.path(config$code_root, 'clustering_methods/auto_correlation'), work,
    config$acf_n_pcs, config$acf_neighbors, config$acf_resolution, config$acf_seed)
  log = file.path(work, 'run.log')
  status = system2(config$acf_python, vapply(as.character(args), shQuote, character(1)),
    stdout = log, stderr = log, env = c('OMP_NUM_THREADS=1', 'OPENBLAS_NUM_THREADS=1',
      'MKL_NUM_THREADS=1', 'NUMBA_NUM_THREADS=1', 'PYTHONDONTWRITEBYTECODE=1',
      paste0('NUMBA_CACHE_DIR=', file.path(work, 'numba')),
      paste0('MPLCONFIGDIR=', file.path(work, 'matplotlib'))))
  if (status != 0L) stop('Autocorrelation failed:\n', paste(readLines(log), collapse = '\n'))
  meta = fread(file.path(work, 'assignments.tsv'), colClasses = list(character = c('cluster','status')))
  stopifnot(identical(meta$row_index, seq_len(nrow(a)) - 1L))
  profiles = as.matrix(fread(file.path(work, 'acf.tsv'), header = FALSE, na.strings = 'nan'))
  rownames(profiles) = a$RID
  colnames(profiles) = as.character(seq_len(ncol(profiles)) - 1L)
  stopifnot(nrow(profiles) == nrow(a), ncol(profiles) == region$width,
    all(abs(profiles[meta$valid, 1] - 1) < 1e-10))
  a$m6a_cluster = a$cluster
  a$acf_status = meta$status
  a$acf_valid = meta$valid
  a$umap1 = meta$umap1
  a$umap2 = meta$umap2
  selected = which(meta$status == 'clustered')
  if (!length(selected)) stop('No ACF profiles could be clustered; see the input read variability.')
  cluster = rep(NA_character_, nrow(a))
  cluster[selected] = as.character(relabel_clusters_by_size(setNames(meta$cluster[selected], a$RID[selected])))
  a$cluster = factor(cluster, levels = paste0('cluster', seq_len(length(unique(cluster[selected])))))
  edges = fread(file.path(work, 'edges.tsv'))
  edges[, `:=`(from = a$RID[source_index + 1L], to = a$RID[target_index + 1L])]
  graph = igraph::graph_from_data_frame(as.data.frame(edges[, .(from, to, weight)]),
    directed = FALSE, vertices = a$RID[selected])
  settings = jsonlite::fromJSON(file.path(work, 'parameters.json'))
  settings$k_eff = settings$actual_neighbors
  settings$read_filter = result$params$read_filter
  settings$n_input_full_window = result$params$n_input_full_window
  settings$n_unassigned_excluded = result$params$n_unassigned_excluded
  settings$lag_min = 0L
  settings$lag_max = region$width - 1L
  list(region = region, assignments = a[selected, , drop = FALSE], all_assignments = a,
    acf = profiles[selected, , drop = FALSE], all_acf = profiles,
    site_met_mat = result$site_met_mat[a$RID[selected], , drop = FALSE],
    graph = graph, n_clusters = nlevels(droplevels(a$cluster)), params = settings)
}

pair_acf_heatmap = function(result) {
  a = result$assignments
  ordering = order(a$cluster, a$start, a$RID)
  display = result$acf[a$RID[ordering], -1, drop = FALSE]
  limit = max(as.numeric(stats::quantile(abs(display), .99)), .01)
  ticks = unique(round(seq(1, ncol(display), length.out = min(6L, ncol(display)))))
  ComplexHeatmap::Heatmap(display, name = 'ACF',
    col = circlize::colorRamp2(c(-limit, 0, limit), c('#2166ac', 'white', '#b2182b')),
    cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE,
    row_split = a$cluster[ordering], row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 8),
    show_row_names = FALSE, show_column_names = FALSE,
    use_raster = TRUE, raster_quality = 2, raster_resize_mat = FALSE, raster_by_magick = FALSE,
    left_annotation = pair_read_annotations(a[ordering, , drop = FALSE], result$region),
    column_title = paste0(result$region$annotation, ' | autocorrelation + Leiden | ', nrow(a), ' phased molecules\n',
      'All nonzero lags shown; REF/ALT labels are annotations only'),
    column_title_gp = grid::gpar(fontsize = 10),
    bottom_annotation = ComplexHeatmap::HeatmapAnnotation(lag_bp = ComplexHeatmap::anno_mark(
      at = ticks, labels = as.character(ticks), which = 'column', side = 'bottom')))
}

pair_acf_profiles = function(result, max_lag = 500L) {
  a = result$all_assignments
  lag = seq_len(min(max_lag, ncol(result$all_acf) - 1L))
  summarize = function(group, selection, panel) {
    selection = selection & a$acf_valid
    if (!any(selection)) return(NULL)
    data.table(panel = panel, group = group, lag_bp = lag,
      mean_acf = colMeans(result$all_acf[selection, lag + 1L, drop = FALSE]), n_reads = sum(selection))
  }
  values = rbindlist(c(
    lapply(levels(a$allele_group), function(x) summarize(x, a$allele_group == x, 'Focal SNP allele')),
    lapply(levels(a$cluster), function(x) summarize(x, !is.na(a$cluster) & a$cluster == x, 'ACF cluster'))))
  colors = c(pair_allele_colors(result$region), cluster_id_palette(levels(a$cluster)))
  p = ggplot(values, aes(lag_bp, mean_acf, color = group)) + geom_line(linewidth = .6) +
    facet_grid(panel ~ .) + scale_color_manual(values = colors) +
    cowplot::theme_cowplot(font_size = 10) +
    labs(title = paste0(result$region$annotation, ' | mean per-molecule autocorrelation'),
      subtitle = paste0('Display: lags 1-', max(lag), ' bp; clustering uses all lags. Zero-variance reads excluded.'),
      x = 'Lag (bp)', y = 'Mean ACF', color = 'Group')
  list(values = values, plot = p)
}

write_pair_allele_acf_plots = function(result, pair, config, output_dir, output_prefix) {
  path = function(suffix, extension = '.pdf') file.path(output_dir, paste0(output_prefix, '_', suffix, extension))
  allele_graph = plot_pair_allele_graph(result,
    paste0(result$region$annotation, ' | Manhattan graph colored by focal SNP allele'))
  allele_composition = pair_cluster_alleles(result,
    paste0(result$region$annotation, ' | allele composition of m6A Leiden clusters'))
  acf = pair_autocorrelation(result, config)
  acf_heatmap = pair_acf_heatmap(acf)
  acf_m6a_heatmap = pair_m6a_heatmap(acf, pair, cluster_label = 'Autocorrelation + Leiden clusters')
  acf_profiles = pair_acf_profiles(acf, config$acf_plot_max_lag)
  acf_graph = plot_knn_graph(acf, main = paste0(result$region$annotation, ' | ACF PCA/correlation graph; Leiden clusters'))
  acf_allele_graph = plot_pair_allele_graph(acf,
    paste0(result$region$annotation, ' | ACF PCA/correlation graph; focal SNP allele'))
  acf_composition = pair_cluster_alleles(acf,
    paste0(result$region$annotation, ' | allele composition of ACF Leiden clusters'))
  plots = list(manhattan_leiden_graph_by_allele = allele_graph,
    leiden_allele_composition = allele_composition$plot,
    autocorrelation_profiles = acf_profiles$plot,
    autocorrelation_graph = acf_graph,
    autocorrelation_graph_by_allele = acf_allele_graph,
    autocorrelation_allele_composition = acf_composition$plot)
  for (name in names(plots)) ggsave(path(name), plots[[name]], width = 12, height = 8)
  heatmap_height = max(8, .012 * nrow(acf$assignments) + 3)
  save_figure(acf_heatmap, path('autocorrelation_heatmap'), width = 14, height = heatmap_height,
    draw = function(x) ComplexHeatmap::draw(x, heatmap_legend_side = 'right', annotation_legend_side = 'right'))
  save_figure(acf_m6a_heatmap, path('m6a_heatmap_autocorrelation'), width = 14,
    height = heatmap_height, draw = draw_pair_m6a_heatmap)
  suffixes = c(names(plots), 'autocorrelation_heatmap', 'm6a_heatmap_autocorrelation')
  list(result = acf, heatmap = acf_heatmap, m6a_heatmap = acf_m6a_heatmap,
    profile_plot = acf_profiles$plot, allele_graph = allele_graph,
    allele_composition = allele_composition, plots = plots,
    files = setNames(vapply(suffixes, path, character(1)), suffixes))
}

# pair-m6a-heatmap-functions
# Reusable wrapper with a top track aggregated over the displayed phased molecules.
pair_m6a_heatmap = function(result, pair, split_by_cluster = TRUE,
    cluster_label = "Manhattan + Leiden clusters") {
  region = result$region
  a = result$assignments
  positions = seq.int(region$analysis_start, region$analysis_end)
  stopifnot(identical(rownames(result$site_met_mat), a$RID),
    all(a$start <= region$analysis_start), all(a$end >= region$analysis_end))
  aggregate = numeric(length(positions))
  aggregate[match(as.integer(colnames(result$site_met_mat)), positions)] =
    as.numeric(Matrix::colMeans(result$site_met_mat))
  if (!split_by_cluster) {
    result$assignments$cluster = factor(rep("All phased molecules", nrow(a)))
  }
  cre_track = rep("Outside tested CREs", length(positions))
  cre_track[positions > pair$cre1_start & positions <= pair$cre1_end] = "CRE1"
  cre_track[positions > pair$cre2_start & positions <= pair$cre2_end] = "CRE2"
  top = ComplexHeatmap::HeatmapAnnotation(
    Phased_reads_m6A = ComplexHeatmap::anno_lines(aggregate, ylim = c(0, 1),
      gp = grid::gpar(col = "#333333", lwd = .8),
      axis_param = list(at = c(0, .5, 1), labels = c("0", ".5", "1")),
      height = grid::unit(20, "mm")),
    Tested_CREs = ComplexHeatmap::anno_simple(cre_track,
      col = c(CRE1 = "#FFFF00", CRE2 = "#FFFF00", "Outside tested CREs" = "white"),
      height = grid::unit(3, "mm")), annotation_name_gp = grid::gpar(fontsize = 8),
    annotation_name_side = "left", annotation_label = c("Phased-read m6A", "cCREs"))
  title = paste0(region$annotation, " | ", nrow(a), " full-window phased molecules | ",
    if (split_by_cluster) cluster_label else "no clustering",
    "\nTop: fraction of displayed phased molecules with m6A; Unassigned reads excluded\n",
    region$chr, ":", region$analysis_start, "-", region$analysis_end)
  heatmap_options = list(top_annotation = top, column_title = title,
    column_title_gp = grid::gpar(fontsize = 10))
  ordering = order(result$assignments$cluster, a$start, a$RID)
  heatmap_options$left_annotation = pair_read_annotations(
    result$assignments[ordering, , drop = FALSE], region, include_cluster = split_by_cluster)
  if (!split_by_cluster) {
    heatmap_options["row_split"] = list(NULL)
    heatmap_options["row_title"] = list(NULL)
  }
  h = plot_read_heatmap(result, region, sample_palette(a$sample_name),
    include_haplotype = FALSE, show_cluster_profiles = FALSE,
    heatmap_options = heatmap_options)
  list(heatmap = h, aggregate = data.table(position = positions, m6a_fraction = aggregate,
    n_molecules = nrow(a)),
    cre_left = (c(pair$cre1_start, pair$cre2_start) - region$start) / region$width,
    cre_right = (c(pair$cre1_end, pair$cre2_end) - region$start) / region$width,
    snp_x = (region$focal_pos - region$analysis_start + .5) / region$width)
}

draw_pair_m6a_heatmap = function(x) {
  grid::grid.newpage()
  h = ComplexHeatmap::draw(x$heatmap, newpage = FALSE,
    heatmap_legend_side = "right", annotation_legend_side = "right")
  highlight = function() {
    grid::grid.rect(x = x$cre_left, y = 0, width = x$cre_right - x$cre_left,
      height = 1, just = c("left", "bottom"),
      gp = grid::gpar(fill = grDevices::adjustcolor("#FFFF00", alpha.f = .18), col = NA))
    grid::grid.lines(x = c(x$snp_x, x$snp_x), y = c(0, 1),
      gp = grid::gpar(col = "#984EA3", lty = 2, lwd = .7))
  }
  order = ComplexHeatmap::row_order(h)
  n_slices = if (is.list(order)) length(order) else 1L
  for (i in seq_len(n_slices))
    ComplexHeatmap::decorate_heatmap_body("m6A", highlight(), slice = i)
  ComplexHeatmap::decorate_annotation("Phased_reads_m6A", highlight())
  invisible(h)
}

# cre-pair-plot-function
plot_cre_pair = function(cre_pair, label = NULL, output_prefix = NULL,
    pooled_pairs = pooled, allele_pairs = pair_allele, config = params,
    output_dir = examples_dir, region_description = NULL) {
  # Accept a complete CRE_pair ID or the two full CRE IDs in either order.
  if (!is.character(cre_pair) || !length(cre_pair) %in% c(1L, 2L) ||
      anyNA(cre_pair) || any(!nzchar(cre_pair)))
    stop("cre_pair must be one CRE_pair ID or a vector containing two full CRE IDs.")
  candidates = if (length(cre_pair) == 1L) cre_pair else
    c(paste(cre_pair, collapse = "-"), paste(rev(cre_pair), collapse = "-"))
  pair = copy(as.data.table(pooled_pairs)[CRE_pair %in% candidates])
  if (nrow(pair) != 1L) stop("Expected exactly one matching row in pooled_pairs; found ", nrow(pair), ".")
  pair_id = pair$CRE_pair
  pair_allele_row = copy(as.data.table(allele_pairs)[CRE_pair == pair_id])
  if (nrow(pair_allele_row) != 1L || is.na(pair_allele_row$variant_id) ||
      !nzchar(pair_allele_row$variant_id))
    stop("The selected pair must have exactly one allele-test row with a focal SNV.")
  short_ids = sub("^[^.]+\\.([^.]+)\\..*$", "\\1", c(pair$CRE1, pair$CRE2))
  if (is.null(label) && "gene_names" %in% names(pair) && !is.na(pair$gene_names)) {
    symbols = trimws(strsplit(pair$gene_names, ";", fixed = TRUE)[[1]])
    symbols = unique(symbols[nzchar(symbols) & !grepl("^ENSG[0-9]+", symbols)])
    if (length(symbols)) label = paste(symbols, collapse = "; ")
  }
  folder_label = if (is.null(label)) "CRE_pair" else label
  if (is.null(label)) label = paste(short_ids, collapse = " / ")
  if (is.null(output_prefix)) output_prefix = paste0("LCL_hsca_", paste(short_ids, collapse = "_"))
  if (length(label) != 1L || is.na(label) || !nzchar(label)) stop("label must be one nonempty string.")
  if (length(output_prefix) != 1L || is.na(output_prefix) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", output_prefix))
    stop("output_prefix must be a filename stem containing letters, numbers, dots, underscores or hyphens.")
  stopifnot(length(config$locus_flank_bp) == 1L, is.finite(config$locus_flank_bp),
    config$locus_flank_bp >= 0)
  folder_label = gsub("[^A-Za-z0-9_.-]+", "_", folder_label)
  output_dir = file.path(output_dir, paste0(folder_label, "_", paste(short_ids, collapse = "__")))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  molecule_manifest = fread(file.path(config$root, "universe", "sample_manifest.tsv"))
  stopifnot(!anyDuplicated(molecule_manifest$sample_name))
  snp = strsplit(pair_allele_row$variant_id, ":", fixed = TRUE)[[1]]
  if (length(snp) != 4L || snp[1] != pair$chrom ||
      !all(snp[3:4] %in% c("A", "C", "G", "T")) || snp[3] == snp[4] ||
      is.na(suppressWarnings(as.integer(snp[2]))))
    stop("The selected pair requires a valid focal SNV (chromosome:position:REF:ALT).")
  pair_region = data.frame(region_id = paste0(output_prefix, "_coaccessible_pair"), region_type = "coaccess_pair",
    annotation = label, chr = pair$chrom,
    start = max(0L, min(pair$cre1_start, pair$cre2_start) - config$locus_flank_bp),
    end = max(pair$cre1_end, pair$cre2_end) + config$locus_flank_bp,
    focal_pos = as.integer(snp[2]), focal_snp = paste0(snp[1], ":", snp[2]),
    ref = snp[3], alt = snp[4])
  pair_region$analysis_start = pair_region$start + 1L
  pair_region$analysis_end = pair_region$end
  pair_region$width = pair_region$end - pair_region$start
  if (is.null(region_description)) region_description = paste0("CRE pair at the ", label, " locus")
  if (length(region_description) != 1L || is.na(region_description) || !nzchar(region_description))
    stop("region_description must be one nonempty string.")
  pair_region$gene_label = label
  pair_region$region_description = region_description
  # Titles use the complete displayed window, in 1-based inclusive hg38 coordinates.
  pair_region$annotation = paste0(label, " | ", region_description, "\n",
    pair_region$chr, ":", format(pair_region$analysis_start, big.mark = ",", scientific = FALSE, trim = TRUE),
    "-", format(pair_region$analysis_end, big.mark = ",", scientific = FALSE, trim = TRUE), " (hg38)")
  phased_input = lcl_build_phased_m6a_input(pair_region, as.data.frame(molecule_manifest),
    config$fire_root, config$phasing_root)
  window = list(chrom = pair_region$chr, start = pair_region$start, end = pair_region$end)
  pair_res = load_region(window, molecule_manifest$sample_name, config$root, config$fire_root,
    tabix_bin = config$tabix, span_source = "aligned_blocks")
  labels_all = label_reads(pair_res,
    list(start = pair$cre1_start, end = pair$cre1_end),
    list(start = pair$cre2_start, end = pair$cre2_end), fire_overlap_fraction = .5)
  original_counts = coaccess_2x2(labels_all)
  stopifnot(identical(as.numeric(unlist(original_counts[c("neither", "cre1_only", "cre2_only", "both", "n_shared")])),
    as.numeric(unlist(pair[, .(co_closed, CRE1_access, CRE2_access, co_access, n_shared_reads)]))))
  a = as.data.table(phased_input$rids_df)
  a[, key := paste(sample_name, original_RID)]
  full_primary = pair_res$spans[start <= window$start & end >= window$end, key]
  a[, plot_retained := key %in% full_primary & key %in% labels_all[shared == TRUE, key]]
  # Keep the phasing audit, including full-span molecules excluded by phasing.
  audit = as.data.table(phased_input$read_filter_audit)
  audit[, plot_retained := RID %in% a[plot_retained == TRUE, RID]]
  a = a[plot_retained == TRUE]
  stopifnot(nrow(a) > 0L, !anyDuplicated(a$key), all(a$haplotype %in% c("HP1", "HP2")),
    all(a$focal_genotype %in% c("0|1", "1|0")), all(a$allele_status == "phased_focal_genotype"))
  allele_bases = sub(".*: ", "", a$allele_label)
  stopifnot(all(allele_bases %in% c(snp[3], snp[4])))
  a[, allele_group := factor(ifelse(allele_bases == snp[3], paste0("REF-", snp[3]), paste0("ALT-", snp[4])),
    levels = c(paste0("REF-", snp[3]), paste0("ALT-", snp[4])))]
  a[, cluster := allele_group] # Renderer facet variable; no clustering is performed.
  pair_labels = merge(a, labels_all[, .(key, shared, acc1, acc2, config)], by = "key", sort = FALSE)
  pair_labels[, `:=`(CRE_pair = pair_id,
    donor = molecule_manifest$donor[match(sample_name, molecule_manifest$sample_name)],
    variant_id = pair_allele_row$variant_id,
    primary_start = pair_res$spans$start[match(key, pair_res$spans$key)],
    primary_end = pair_res$spans$end[match(key, pair_res$spans$key)])]
  pair_contributions = pair_labels[, .(n_shared = .N, both = sum(acc1 & acc2),
    cre1_only = sum(acc1 & !acc2), cre2_only = sum(!acc1 & acc2), neither = sum(!acc1 & !acc2)),
    by = .(sample_name, donor, allele_group)]
  pair_result = list(region = pair_region, assignments = as.data.frame(a),
    site_met_mat = phased_input$met_mat[a$RID, , drop = FALSE])
  # Footprints keep the autocorrelation figure's size ranges and colors.
  # Apply the size filter before clipping calls at aligned-block gaps.
  display_res = pair_res
  display_res$blocks = NULL
  display_res = load_ft_tracks(display_res, keys = a$key, ft_root = config$fire_root,
    tabix_bin = config$tabix, nuc_source = "ft", tf_source = "recalled_tf", tf_root = config$tf_root)
  nuc = copy(display_res$nuc)[end - start >= 130 & end - start <= 160]
  nuc[, track := "ft_nuc_130-160bp"]
  nuc = clip_to_aligned_blocks(nuc, pair_res$blocks)
  tf = clip_to_aligned_blocks(copy(display_res$size_fps), pair_res$blocks)
  tf[, track := paste0("FiberHMM_", sub("^size", "", class), "bp")]
  footprints = rbindlist(list(nuc[, .(key, start, end, track)], tf[, .(key, start, end, track)]))
  footprints[, RID := a$RID[match(key, a$key)]]
  stopifnot(!anyNA(footprints$RID))
  records = rbindlist(list(as.data.table(m6a_intervals(pair_result)), footprints), fill = TRUE)
  pair_plot = plot_phased_pair_tracks(pair_result, as.data.frame(records), pair)
  bar_data = pair_labels[, .N, by = .(sample_name, allele_group, config)]
  pair_plot$bars = ggplot(bar_data, aes(x = sub("_.*$", "", sample_name), y = N, fill = config)) +
    geom_col(position = "fill") + facet_grid(allele_group ~ .) + scale_fill_manual(values = CONFIG_COLS, drop = FALSE) +
    scale_y_continuous(labels = scales::percent) + cowplot::theme_cowplot(font_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), legend.position = "bottom") +
    labs(x = "Sample", y = "Fraction of retained molecules", fill = "Configuration",
      title = paste0(pair_region$annotation, " | full-window phased molecules by SNP allele"))
  ggsave(file.path(output_dir, paste0(output_prefix, "_molecule_tracks.pdf")), pair_plot$panels$combined,
    width = 12, height = pair_plot$height, limitsize = FALSE)
  ggsave(file.path(output_dir, paste0(output_prefix, "_configurations.pdf")), pair_plot$bars,
    width = 12, height = 6)
  pair_all_input = assemble_region_m6a(as.data.frame(molecule_manifest), pair_region,
    full_span = TRUE, matrix_dir = NULL, reuse = FALSE)
  input = list(region = pair_region, assignments = as.data.frame(pair_all_input$rids_df),
    site_met_mat = pair_all_input$met_mat)
  all_keys = paste(input$assignments$sample_name, input$assignments$original_RID)
  input$assignments$in_phased_plot = all_keys %in% pair_labels$key
  input$assignments$coaccess_eligible = all_keys %in% labels_all[shared == TRUE, key]
  input$assignments$n_m6a_calls = as.numeric(Matrix::rowSums(input$site_met_mat))
  stopifnot(nrow(input$assignments) == sum(pair_all_input$qc$full_span_reads),
    identical(input$assignments$RID, rownames(input$site_met_mat)),
    sum(input$assignments$in_phased_plot) == nrow(pair_labels))
  input = annotate_pair_alleles(input, audit)
  pair_clustered = cluster_phased_pair(input, audit, config)
  pair_heatmap_unclustered = pair_m6a_heatmap(pair_clustered, pair, split_by_cluster = FALSE)
  pair_heatmap_leiden = pair_m6a_heatmap(pair_clustered, pair, split_by_cluster = TRUE)
  pair_heatmap_height = max(8, .012 * nrow(pair_clustered$assignments) + 3)
  for (kind in c("phased", "leiden")) {
    heatmap = if (kind == "phased") pair_heatmap_unclustered else pair_heatmap_leiden
    save_figure(heatmap, file.path(output_dir, paste0(output_prefix, "_m6a_heatmap_", kind, ".pdf")),
      width = 14, height = pair_heatmap_height, draw = draw_pair_m6a_heatmap)
  }
  # Dense profiles retain zero-call positions, with the denominator fixed within each group.
  positions = pair_heatmap_unclustered$aggregate$position
  groups = levels(pair_clustered$assignments$cluster)
  pair_cluster_profiles = rbindlist(lapply(c("All phased molecules", groups), function(group) {
    selected = if (group == "All phased molecules") rep(TRUE, nrow(pair_clustered$assignments)) else
      pair_clustered$assignments$cluster == group
    fraction = numeric(length(positions))
    fraction[match(as.integer(colnames(pair_clustered$site_met_mat)), positions)] =
      as.numeric(Matrix::colMeans(pair_clustered$site_met_mat[selected, , drop = FALSE]))
    data.table(group = group, position = positions, fraction = fraction, n_reads = sum(selected))
  }))
  pair_cluster_profiles[, group := factor(group, levels = c("All phased molecules", groups))]
  cre_bands = data.frame(left = c(pair$cre1_start, pair$cre2_start) + 1L,
    right = c(pair$cre1_end, pair$cre2_end))
  pair_profile_plot = ggplot(pair_cluster_profiles, aes(position, fraction, color = group)) +
    geom_rect(data = cre_bands, aes(xmin = left, xmax = right, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "#FFFF00", alpha = .22) + geom_line(linewidth = .35) +
    facet_grid(group ~ .) + geom_vline(xintercept = pair_region$focal_pos,
      color = "#984EA3", linetype = "dashed") +
    scale_color_manual(values = c("All phased molecules" = "#333333", cluster_id_palette(groups))) +
    scale_y_continuous(limits = c(0, 1)) + cowplot::theme_cowplot(font_size = 10) +
    theme(legend.position = "none") +
    labs(title = paste0(pair_region$annotation, " | aggregate and Leiden-cluster m6A profiles"),
      x = paste0(pair_region$chr, " position (1-based)"), y = "Fraction of group molecules with m6A")
  pair_knn_plot = plot_knn_graph(pair_clustered)
  ggsave(file.path(output_dir, paste0(output_prefix, "_leiden_m6a_profiles.pdf")), pair_profile_plot,
    width = 12, height = 1.3 * (length(groups) + 1) + 1)
  ggsave(file.path(output_dir, paste0(output_prefix, "_manhattan_leiden_graph.pdf")), pair_knn_plot,
    width = 10, height = 8)
  acf_plots = write_pair_allele_acf_plots(pair_clustered, pair, config, output_dir, output_prefix)
  pdf_suffixes = c("molecule_tracks", "configurations", "m6a_heatmap_phased",
    "m6a_heatmap_leiden", "leiden_m6a_profiles", "manhattan_leiden_graph")
  files = setNames(file.path(output_dir, paste0(output_prefix, "_", pdf_suffixes, ".pdf")), pdf_suffixes)
  stopifnot(all(file.exists(c(files, acf_plots$files))))
  invisible(list(pair_id = pair_id, pair = pair, allele = pair_allele_row, region = pair_region,
    phased = pair_result, read_labels = pair_labels, sample_counts = pair_contributions,
    phasing_audit = audit, records = records, phased_plot = pair_plot,
    phased_molecules = pair_clustered, heatmap_unclustered = pair_heatmap_unclustered, heatmap_leiden = pair_heatmap_leiden,
    heatmap_height = pair_heatmap_height, cluster_profiles = pair_cluster_profiles,
    profile_plot = pair_profile_plot, graph_plot = pair_knn_plot,
    autocorrelation = acf_plots, output_dir = output_dir, files = c(files, acf_plots$files)))
}

# domain-ctcf-functions
# Exon rows (BED coordinates) of each protein-coding and lncRNA gene's Ensembl_canonical transcript.
load_canonical_gene_models = function(gtf) {
  awk = '$3 == "exon" && $9 ~ /Ensembl_canonical/ && $9 ~ /gene_type "(protein_coding|lncRNA)"/ {print $1, $4 - 1, $5, $7, $9}'
  exons = fread(cmd = paste("zcat", shQuote(gtf), "| awk -F'\\t' -v OFS='\\t'", shQuote(awk)),
    header = FALSE, sep = "\t", quote = "", col.names = c("chrom", "start", "end", "strand", "attr"))
  exons[, `:=`(gene_name = sub('.*gene_name "([^"]+)".*', "\\1", attr),
    gene_type = sub('.*gene_type "([^"]+)".*', "\\1", attr),
    transcript_id = sub('.*transcript_id "([^"]+)".*', "\\1", attr), attr = NULL)]
  exons[]
}

# Rao et al. 2014 GM12878 HiCCUPS loops and Arrowhead domains lifted to hg38
# (rao2014_GM12878/liftover_to_hg38.sh). loopN = row N of the lifted loop list. Rao's motif
# columns are 1-based starts (sequence = hg38[motif_x1 - 1, motif_x2 - 1)); stored here as BED.
load_rao_hg38 = function(dir) {
  read = function(name) fread(file.path(dir, paste0("GSE63525_GM12878_primary+replicate_", name, ".hg38.txt.gz")))
  key = c("chr1", "x1", "x2", "y1", "y2")
  loops = read("HiCCUPS_looplist")[, ..key][, loop_id := paste0("loop", .I)]
  motifs = read("HiCCUPS_looplist_with_motifs")[, c(key, "motif_x1", "motif_x2", "orientation_1",
    "motif_y1", "motif_y2", "orientation_2"), with = FALSE]
  stopifnot(!anyDuplicated(loops[, ..key]), !anyDuplicated(motifs[, ..key]))
  loops = merge(loops, motifs, by = key, all.x = TRUE)[order(as.integer(sub("loop", "", loop_id)))]
  loops[, convergent := fifelse(is.na(orientation_1) | is.na(orientation_2), NA,
    orientation_1 == "p" & orientation_2 == "n")]
  strand = c(p = "+", n = "-")
  side = function(s, e, ms, me, o, os, oe, label) loops[, .(loop_id, side = label, chrom = chr1,
    start = get(s), end = get(e), motif_start = get(ms) - 1L, motif_end = get(me) - 1L,
    motif_strand = unname(strand[get(o)]), partner = paste0(chr1, ":", get(os), "-", get(oe)),
    loop_size = y1 - x1, convergent)]
  anchors = rbind(side("x1", "x2", "motif_x1", "motif_x2", "orientation_1", "y1", "y2", "A"),
    side("y1", "y2", "motif_y1", "motif_y2", "orientation_2", "x1", "x2", "B"))
  domains = read("Arrowhead_domainlist")[, .(chrom = chr1, start = x1, end = x2)]
  edges = unique(rbind(domains[, .(chrom, pos = start)], domains[, .(chrom, pos = end)]))
  list(loops = loops, anchors = anchors, domains = domains, edges = edges)
}

# Shared look of the browser tracks: no y axis; each track sets its own x range.
domain_ctcf_theme = function() cowplot::theme_cowplot(font_size = 9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = .5, size = 8), plot.margin = margin(2, 5, 2, 5))

# Canonical transcripts overlapping win_chrom:[win_start, win_end) (BED), packed into rows: intron
# line with strand chevrons, exon boxes, gene name. With region, x is SNP-relative bp
# (plot_positions()); otherwise genomic Mb. A zoom window without transcripts names the nearest
# protein-coding TSS. Returns list(plot, rows).
plot_gene_track = function(models, win_chrom, win_start, win_end, xlim, region = NULL, focus = NULL) {
  to_x = function(pos) if (is.null(region)) pos / 1e6 else plot_positions(pos, region)
  tx = models[chrom == win_chrom, .(start = min(start), end = max(end)),
    by = .(gene_name, gene_type, strand, transcript_id)]
  shown = tx[start < win_end & end > win_start][order(start)]
  p = ggplot() + scale_x_continuous(expand = c(0, 0)) + coord_cartesian(xlim = xlim) +
    labs(x = NULL, y = "Genes") + domain_ctcf_theme()
  if (!is.null(focus)) p = p + annotate("rect", xmin = to_x(focus[1] + 1), xmax = to_x(focus[2]),
    ymin = -Inf, ymax = Inf, fill = "#FFFF00", alpha = .35)
  if (!nrow(shown)) {
    label = "No canonical transcript in window"
    if (!is.null(region)) {
      tss = tx[gene_type == "protein_coding", .(gene_name, tss = fifelse(strand == "+", start + 1L, end))]
      near = tss[which.min(abs(tss - region$focal_pos))]
      label = sprintf("%s; nearest protein-coding TSS: %s at %s:%s (%+.0f kb from SNP)", label,
        near$gene_name, win_chrom, format(near$tss, big.mark = ","), (near$tss - region$focal_pos) / 1e3)
    }
    return(list(plot = p + annotate("text", x = mean(xlim), y = 0, size = 3, label = label) +
      scale_y_continuous(limits = c(-1, 1)), rows = 1L))
  }
  # greedy packing; the gap leaves room for the name
  gap = .12 * (win_end - win_start)
  row_end = numeric()
  shown[, row := 0L]
  for (i in seq_len(nrow(shown))) {
    r = which(row_end + gap < shown$start[i])[1]
    if (is.na(r)) r = length(row_end) + 1L
    row_end[r] = shown$end[i]
    set(shown, i, "row", r)
  }
  shown[, `:=`(left = to_x(pmax(start, win_start) + 1), right = to_x(pmin(end, win_end)))]
  exons = merge(models[chrom == win_chrom & start < win_end & end > win_start,
    .(transcript_id, xmin = to_x(pmax(start, win_start) + 1), xmax = to_x(pmin(end, win_end)))],
    shown[, .(transcript_id, row)], by = "transcript_id")
  orientation = if (is.null(region)) 1 else plot_anchor(region)$direction
  step = diff(xlim) / 40
  chevrons = shown[, .(x = seq(min(left, right), max(left, right), by = step)), by = .(transcript_id, row, strand)]
  chevrons[, xend := x + fifelse(strand == "+", 1, -1) * orientation * .35 * step]
  name = if (is.null(region)) shown$gene_name else paste0(shown$gene_name, " (", shown$strand, ")")
  p = p + geom_segment(data = shown, aes(x = left, xend = right, y = -row, yend = -row),
      color = "grey30", linewidth = .3) +
    geom_segment(data = chevrons, aes(x = x, xend = xend, y = -row, yend = -row), color = "grey30",
      linewidth = .3, arrow = arrow(length = unit(1, "mm"), type = "open")) +
    geom_rect(data = exons, aes(xmin = xmin, xmax = xmax, ymin = -row - .22, ymax = -row + .22),
      fill = "#1F3A93") +
    annotate("text", x = (shown$left + shown$right) / 2, y = -shown$row + .45, label = name,
      size = 2.6, fontface = "italic") +
    scale_y_continuous(limits = c(-max(shown$row) - .5, -.3))
  list(plot = p, rows = max(shown$row))
}

# Context tracks over chrom:[start, end): HiCCUPS loops with both anchors in view as arcs between
# anchor midpoints (loops with an anchor in focus coloured and labelled), anchors and CTCF motif
# arrows at the base; Arrowhead domains as triangles with edge ticks.
plot_rao_context = function(rao, win_chrom, win_start, win_end, focus) {
  xlim = c(win_start, win_end) / 1e6
  band = annotate("rect", xmin = (focus[1] + 1) / 1e6, xmax = focus[2] / 1e6, ymin = -Inf, ymax = Inf,
    fill = "#FFFF00", alpha = .35)
  loops = rao$loops[chr1 == win_chrom & x1 >= win_start & y2 <= win_end]
  loops[, focal := (x1 < focus[2] & x2 > focus[1]) | (y1 < focus[2] & y2 > focus[1])]
  tallest = max(c(loops$y1 - loops$x1, 1))
  arc_t = seq(0, pi, length.out = 60)
  arcs = loops[, .(x = ((x1 + x2 + y1 + y2) / 4 - (y1 + y2 - x1 - x2) / 4 * cos(arc_t)) / 1e6,
    y = (y1 - x1) / tallest * sin(arc_t)), by = .(loop_id, focal)]
  anchors = rao$anchors[loop_id %in% loops$loop_id]
  anchors[, focal := loop_id %in% loops[focal == TRUE, loop_id]]
  motifs = unique(anchors[!is.na(motif_start), .(motif_start, motif_end, motif_strand)])
  arrow_bp = .012 * (win_end - win_start)
  labels = loops[focal == TRUE, .(x = (x1 + x2 + y1 + y2) / 4e6, y = (y1 - x1) / tallest + .1,
    label = sprintf("%s (%.0f kb%s)", loop_id, (y1 - x1) / 1e3,
      fifelse(is.na(convergent), "", fifelse(convergent, ", convergent CTCF", ", not convergent"))))]
  colors = c("TRUE" = "#B2182B", "FALSE" = "grey70")
  p_loops = ggplot() + band +
    geom_path(data = arcs, aes(x, y, group = loop_id, color = as.character(focal)), linewidth = .4) +
    geom_rect(data = anchors, aes(xmin = (start + 1) / 1e6, xmax = end / 1e6, ymin = -.12, ymax = -.02,
      fill = as.character(focal))) +
    geom_segment(data = motifs, aes(x = (motif_start + motif_end) / 2e6,
      xend = ((motif_start + motif_end) / 2 + fifelse(motif_strand == "+", arrow_bp, -arrow_bp)) / 1e6,
      y = -.22, yend = -.22, color = motif_strand), linewidth = .5,
      arrow = arrow(length = unit(1.3, "mm"), type = "closed")) +
    geom_text(data = labels, aes(x, y, label = label), size = 2.5, color = colors[["TRUE"]]) +
    scale_color_manual(values = c(colors, "+" = "#D6604D", "-" = "#4393C3"), guide = "none") +
    scale_fill_manual(values = colors, guide = "none") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(-.3, 1.25)) +
    coord_cartesian(xlim = xlim) + labs(x = NULL, y = "HiCCUPS\nloops\n(CTCF motif\narrows)") +
    domain_ctcf_theme()
  if (!nrow(loops)) p_loops = p_loops + annotate("text", x = mean(xlim), y = .5, size = 3,
    label = "No HiCCUPS loop with both anchors in this window")
  domains = rao$domains[chrom == win_chrom & start < win_end & end > win_start]
  widest = max(c(domains$end - domains$start, 1))
  triangles = domains[, .(x = c(start, (start + end) / 2, end) / 1e6, y = c(0, (end - start) / widest, 0)),
    by = .(id = paste(start, end))]
  edges = rao$edges[chrom == win_chrom & pos >= win_start & pos <= win_end]
  p_domains = ggplot() + band +
    geom_polygon(data = triangles, aes(x, y, group = id), fill = "grey60", alpha = .25,
      color = "grey35", linewidth = .25) +
    geom_segment(data = edges, aes(x = pos / 1e6, xend = pos / 1e6, y = -.12, yend = 0),
      color = "grey20", linewidth = .3) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(-.15, 1.05)) +
    coord_cartesian(xlim = xlim) +
    labs(x = paste0(win_chrom, " position (Mb, hg38)"), y = "Arrowhead\ndomains") + domain_ctcf_theme()
  list(loops = p_loops, domains = p_domains)
}

# Zoom track over the pair window: loop anchors overlapping it (clipped, one row per anchor interval,
# labelled with every loop using it), CTCF motif arrows and domain edges; SNP-relative x.
# Returns list(plot, rows).
plot_rao_zoom = function(rao, region, xlim) {
  s = region$start
  e = region$end
  to_x = function(pos) plot_positions(pos, region)
  anchors = rao$anchors[chrom == region$chr & start < e & end > s]
  rows = anchors[, .(label = paste(sprintf("%s %s anchor -> partner %s (%.0f kb loop%s)", loop_id, side,
      partner, loop_size / 1e3, fifelse(is.na(convergent), "", fifelse(convergent, ", convergent CTCF",
      ", not convergent"))), collapse = "\n")), by = .(start, end)][order(start, end)]
  rows[, row := .I]
  motifs = unique(anchors[!is.na(motif_start) & motif_start < e & motif_end > s,
    .(motif_start, motif_end, motif_strand)])
  edges = rao$edges[chrom == region$chr & pos > s & pos <= e]
  p = ggplot() + scale_x_continuous(expand = c(0, 0)) + coord_cartesian(xlim = xlim) +
    labs(x = NULL, y = "Rao 2014") + domain_ctcf_theme()
  if (!nrow(rows) && !nrow(edges)) {
    distance = function(pos) min(abs(pos - region$focal_pos)) / 1e3
    others = rao$anchors[chrom == region$chr]
    return(list(plot = p + scale_y_continuous(limits = c(-1, 1)) + annotate("text", x = mean(xlim), y = 0,
      size = 3, label = sprintf("No loop anchor or domain edge in window; nearest anchor %.0f kb, nearest domain edge %.0f kb from SNP",
        distance(c(others$start, others$end)), distance(rao$edges[chrom == region$chr, pos]))), rows = 1L))
  }
  if (nrow(edges)) p = p + geom_vline(xintercept = to_x(edges$pos), linetype = "dashed",
    color = "grey20", linewidth = .4) +
    annotate("text", x = to_x(edges$pos), y = .6, label = " domain edge", size = 2.4, hjust = 0)
  if (nrow(rows)) p = p + geom_rect(data = rows, aes(xmin = to_x(pmax(start, s) + 1),
      xmax = to_x(pmin(end, e)), ymin = -row - .35, ymax = -row + .05), fill = "#B2182B", alpha = .3) +
    annotate("text", x = xlim[1] + .01 * diff(xlim), y = -rows$row - .15, label = rows$label,
      hjust = 0, size = 2.3, lineheight = .9)
  if (nrow(motifs)) p = p +
    geom_rect(data = motifs, aes(xmin = to_x(motif_start + 1), xmax = to_x(motif_end), ymin = .15,
      ymax = .35, color = motif_strand), fill = NA, linewidth = .3) +
    geom_segment(data = motifs, aes(x = to_x((motif_start + motif_end) / 2),
      xend = to_x((motif_start + motif_end) / 2) + fifelse(motif_strand == "+", 1, -1) * .02 * diff(xlim),
      y = .25, yend = .25, color = motif_strand), linewidth = .6,
      arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
    annotate("text", x = to_x((motifs$motif_start + motifs$motif_end) / 2), y = .5, size = 2.4,
      label = paste0("CTCF motif (", motifs$motif_strand, ")")) +
    scale_color_manual(values = c("+" = "#D6604D", "-" = "#4393C3"), guide = "none")
  n = max(nrow(rows), 1L)
  list(plot = p + scale_y_continuous(limits = c(-n - .6, .75)), rows = n)
}

# Genotype of every sample at the focal SNV from each read's own base (hsca_functions.read_alleles,
# the rule of the pooled allele test), so homozygous donors and samples with wrong or missing VCFs
# are labelled too. A sample needs >= 5 allele-informative reads; it is het when the minor allele
# has >= 2 reads and >= 15% of them (a lone discordant HiFi read, or 2 of 48, is an error, not a
# het), else homozygous for the majority allele.
pair_genotype_calls = function(region, manifest, python, hsca_functions) {
  work = tempfile("pair_genotypes_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  fwrite(manifest[, .(sample_name, cram_path)], file.path(work, "manifest.tsv"), sep = "\t")
  runner = r"---(
import importlib.util, sys
import pandas as pd
source, manifest, chrom, pos, ref, alt, out = sys.argv[1:8]
spec = importlib.util.spec_from_file_location("hsca_functions", source)
hsca = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hsca)
variant = pd.DataFrame(dict(variant_id=["focal"], pos=[int(pos)], ref=[ref], alt=[alt]))
rows = [(s.sample_name, read, allele)
        for s in pd.read_csv(manifest, sep="\t").itertuples(index=False)
        for read, allele in hsca.read_alleles(s.cram_path, chrom, variant)["focal"].items()]
pd.DataFrame(rows, columns=["sample_name", "read_name", "allele"]).to_csv(out, sep="\t", index=False)
)---"
  writeLines(runner, file.path(work, "run.py"))
  log = file.path(work, "run.log")
  args = c(file.path(work, "run.py"), hsca_functions, file.path(work, "manifest.tsv"), region$chr,
    region$focal_pos, region$ref, region$alt, file.path(work, "alleles.tsv"))
  status = system2(python, shQuote(args), stdout = log, stderr = log, env = "PYTHONDONTWRITEBYTECODE=1")
  if (status != 0L) stop("Focal SNV pileup failed:\n", paste(readLines(log), collapse = "\n"))
  alleles = fread(file.path(work, "alleles.tsv"))
  calls = merge(manifest[, .(sample_name, donor)],
    alleles[, .(n_ref = sum(allele == "ref"), n_alt = sum(allele == "alt")), by = sample_name],
    by = "sample_name", all.x = TRUE)
  calls[is.na(n_ref), `:=`(n_ref = 0L, n_alt = 0L)]
  calls[, alt_fraction := n_alt / (n_ref + n_alt)]
  calls[, genotype := fcase(n_ref + n_alt < 5L, "uncalled", pmin(n_ref, n_alt) >= 2L &
    pmin(n_ref, n_alt) >= .15 * (n_ref + n_alt), "het",
    n_ref > n_alt, "hom_ref", default = "hom_alt")]
  list(calls = calls[order(as.integer(sub("^AL-?([0-9]+).*$", "\\1", sample_name)))], alleles = alleles)
}

# One PDF per pair: context (genes, HiCCUPS loops with CTCF motifs, Arrowhead domains) over the pair
# window (genes, Rao features, cCREs, per-base m6A by genotype group, m6A molecules). Reads: the
# full-span molecules of all samples that are full-window primary alignments covering both cCREs
# (plot_cre_pair()'s molecule-track universe without its phasing requirement). Groups: hom-REF
# donors, het donors' REF-base and ALT-base reads, hom-ALT donors. No clustering.
plot_cre_pair_domain_ctcf = function(cre_pair, label, region_description = NULL,
    context_flank_bp = 250000, output_dir = file.path(hsca_dir, "CRE_domain_CTCF_anno"),
    config = params, models = gene_models, rao = rao_hg38,
    python = "/project/spott/cshan/envs/Jupyter-notebook/bin/python3") {
  candidates = if (length(cre_pair) == 1L) cre_pair else
    c(paste(cre_pair, collapse = "-"), paste(rev(cre_pair), collapse = "-"))
  pair = pooled[CRE_pair %in% candidates]
  stopifnot(nrow(pair) == 1L)
  allele_row = pair_allele[CRE_pair == pair$CRE_pair]
  stopifnot(nrow(allele_row) == 1L, !is.na(allele_row$variant_id))
  snp = strsplit(allele_row$variant_id, ":", fixed = TRUE)[[1]]
  short_ids = sub("^[^.]+\\.([^.]+)\\..*$", "\\1", c(pair$CRE1, pair$CRE2))
  if (is.null(region_description)) region_description = paste0("CRE pair at the ", label, " locus")
  region = data.frame(region_id = paste0(label, "_domain_ctcf"), region_type = "coaccess_pair",
    chr = pair$chrom, start = max(0L, min(pair$cre1_start, pair$cre2_start) - config$locus_flank_bp),
    end = max(pair$cre1_end, pair$cre2_end) + config$locus_flank_bp,
    focal_pos = as.integer(snp[2]), focal_snp = paste0(snp[1], ":", snp[2]), ref = snp[3], alt = snp[4])
  region$analysis_start = region$start + 1L
  region$analysis_end = region$end
  region$width = region$end - region$start
  region$annotation = paste0(label, " | ", region_description, "\n", region$chr, ":",
    format(region$analysis_start, big.mark = ","), "-", format(region$analysis_end, big.mark = ","), " (hg38)")

  # reads, with the pooled 2x2 reproduced as a check of the read rule
  manifest = fread(file.path(config$root, "universe", "sample_manifest.tsv"))
  window = list(chrom = region$chr, start = region$start, end = region$end)
  spans = load_region(window, manifest$sample_name, config$root, config$fire_root,
    tabix_bin = config$tabix, span_source = "aligned_blocks")
  labels_all = label_reads(spans, list(start = pair$cre1_start, end = pair$cre1_end),
    list(start = pair$cre2_start, end = pair$cre2_end), fire_overlap_fraction = .5)
  counts = coaccess_2x2(labels_all)
  stopifnot(identical(as.numeric(unlist(counts[c("neither", "cre1_only", "cre2_only", "both", "n_shared")])),
    as.numeric(unlist(pair[, .(co_closed, CRE1_access, CRE2_access, co_access, n_shared_reads)]))))
  m6a = assemble_region_m6a(as.data.frame(manifest), region, full_span = TRUE, matrix_dir = NULL, reuse = FALSE)
  a = as.data.table(m6a$rids_df)
  a[, key := paste(sample_name, original_RID)]
  a = a[key %in% spans$spans[start <= window$start & end >= window$end, key] &
    key %in% labels_all[shared == TRUE, key]]
  a[, donor := manifest$donor[match(sample_name, manifest$sample_name)]]

  # genotype groups
  genotypes = pair_genotype_calls(region, manifest, python,
    file.path(config$code_root, "co-accessibility/LCL/hsca/hsca_functions.py"))
  group_labels = c(hom_ref = paste0(region$ref, region$ref),
    het_ref = paste0(region$ref, region$alt, ": ", region$ref, " reads"),
    het_alt = paste0(region$ref, region$alt, ": ", region$alt, " reads"),
    hom_alt = paste0(region$alt, region$alt))
  group_colors = setNames(c("#E6550D", "#FDBE85", "#9ECAE1", "#3182BD"), group_labels)
  bases = genotypes$alleles[, setNames(allele, paste(sample_name, read_name))]
  a[, genotype := genotypes$calls$genotype[match(sample_name, genotypes$calls$sample_name)]]
  a[, base := unname(bases[key])]
  a[, group := fcase(genotype == "hom_ref", "hom_ref", genotype == "hom_alt", "hom_alt",
    genotype == "het" & base == "ref", "het_ref", genotype == "het" & base == "alt", "het_alt")]
  het_base_rate = a[genotype == "het", mean(!is.na(base))]
  n_excluded = a[is.na(group), .N]
  a = a[!is.na(group)]
  a[, cluster := droplevels(factor(group_labels[group], levels = group_labels))]
  group_sizes = a[, .(molecules = .N, samples = uniqueN(sample_name), donors = uniqueN(donor)),
    by = cluster][order(cluster)]
  missing = setdiff(group_labels, levels(a$cluster))
  result = list(region = region, assignments = as.data.frame(a),
    site_met_mat = m6a$met_mat[a$RID, , drop = FALSE])

  # zoom x range = plot_smf_reads() limits without the haplotype bar
  anchor = plot_anchor(region)
  bounds = c(anchor$left - .085 * max(1, anchor$right - anchor$left), anchor$right + .5)
  ticks = pretty(c(anchor$left, anchor$right), n = 6)
  zoom_x = function(p) suppressMessages(p + scale_x_continuous(expand = c(0, 0),
    breaks = ticks[ticks >= anchor$left & ticks <= anchor$right]))
  classes = sub(".*\\.", "", c(pair$CRE1, pair$CRE2))
  cre = data.frame(left = plot_positions(c(pair$cre1_start, pair$cre2_start) + 1L, region),
    right = plot_positions(c(pair$cre1_end, pair$cre2_end), region),
    label = paste0(c("CRE1 (", "CRE2 ("), classes, ")"),
    color = ifelse(is.na(CCRE_COLORS[classes]), "grey40", CCRE_COLORS[classes]))
  bands = function() geom_rect(data = cre, aes(xmin = left, xmax = right, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#FFFF00", alpha = .22)
  snp_line = function() geom_vline(xintercept = 0, color = "#984EA3", linetype = "dashed", linewidth = .4)

  # context: pair midpoint +- flank, widened to whole loops with an anchor in the pair window
  mid = (region$start + region$end) %/% 2
  touching = rao$loops[chr1 == region$chr & ((x1 < region$end & x2 > region$start) |
    (y1 < region$end & y2 > region$start))]
  context = c(max(0, min(mid - context_flank_bp, touching$x1 - 20000)),
    max(mid + context_flank_bp, touching$y2 + 20000))
  focus = c(region$start, region$end)
  context_genes = plot_gene_track(models, region$chr, context[1], context[2], xlim = context / 1e6, focus = focus)
  context_genes$plot = context_genes$plot + labs(title = region$annotation, subtitle = sprintf(
    "SNP %s %s>%s | top: context, yellow band = pair window | bottom: pair window, dashed line at the SNP%s",
    format(region$focal_pos, big.mark = ","), region$ref, region$alt,
    if (length(missing)) paste0(" | no molecules: ", paste(missing, collapse = ", ")) else ""))
  context_rao = plot_rao_context(rao, region$chr, context[1], context[2], focus)

  # zoom tracks
  zoom_genes = plot_gene_track(models, region$chr, region$start, region$end, xlim = bounds, region = region)
  zoom_rao = plot_rao_zoom(rao, region, bounds)
  cre_strip = ggplot(cre) + bands() +
    geom_rect(aes(xmin = left, xmax = right, ymin = .3, ymax = .6, fill = color)) +
    geom_text(aes(x = (left + right) / 2, y = .8, label = label), size = 3) + snp_line() +
    annotate("point", x = 0, y = .1, shape = 18, size = 3, color = "#984EA3") +
    scale_fill_identity() + scale_x_continuous(expand = c(0, 0)) + coord_cartesian(xlim = bounds) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + labs(x = NULL, y = "cCREs") + domain_ctcf_theme()
  profile = as.data.table(cluster_site_profiles(result, result$site_met_mat))
  strip_labels = group_sizes[, setNames(sprintf("%s\n%d molecules\n%d donors", cluster, molecules, donors),
    as.character(cluster))]
  profile[, strip := factor(strip_labels[cluster], levels = strip_labels)]
  profiles = plot_group_profile(profile, group_col = "strip", x_col = "pos", y_col = "met",
    facet_col = "strip", facet_layout = "grid", facet_scales = "fixed",
    colors = setNames(group_colors[names(strip_labels)], strip_labels), coordinates = "relative",
    region = region, layers = list(list(geom = "line", linewidth = .3)))
  profiles$layers = c(list(bands()), profiles$layers)
  profiles = profiles + snp_line() + scale_x_continuous(expand = c(0, 0)) + coord_cartesian(xlim = bounds) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    labs(x = NULL, y = "Per-base\nfraction of\nmolecules\nwith m6A") + domain_ctcf_theme() +
    theme(legend.position = "none", strip.text.y = element_text(angle = 0, size = 7),
      axis.text.y = element_text(size = 7), axis.ticks.y = element_line(), axis.line.y = element_line())

  # molecules: plot_smf_reads() with the group palette on its first bar (private copy)
  renderer = plot_smf_reads
  renderer_env = new.env(parent = environment(renderer))
  renderer_env$cluster_id_palette = function(groups) group_colors[groups]
  environment(renderer) = renderer_env
  molecules = renderer(result, m6a_intervals(result), sample_palette(a$sample_name), tracks = "m6A")
  molecules$layers = c(list(bands()), molecules$layers)
  molecules = suppressMessages(molecules + snp_line() +
    geom_vline(xintercept = c(cre$left, cre$right), color = "grey60", linetype = "dotted", linewidth = .25) +
    labs(title = NULL, subtitle = NULL, y = "m6A: one row per molecule", caption = paste0(
      "Left bars: genotype group, sample. Genotype from each sample's reads at the SNP ",
      "(>= 5 reads; het if the minor allele has >= 2 reads and >= 15%, else the majority allele). ",
      n_excluded, " molecules without a group excluded.")) +
    guides(fill = guide_legend(ncol = 11, byrow = TRUE)))

  zoom_genes$plot = zoom_x(zoom_genes$plot)
  zoom_rao$plot = zoom_x(zoom_rao$plot)
  cre_strip = zoom_x(cre_strip)
  profiles = zoom_x(profiles)
  heights = c(.9 + .3 * context_genes$rows, 1.8, 1.1, .5 + .3 * zoom_genes$rows, .6 + .35 * zoom_rao$rows,
    .7, .8 * nlevels(a$cluster), max(3.5, .035 * nrow(a) + 1.8))
  combined = cowplot::plot_grid(context_genes$plot, context_rao$loops, context_rao$domains, zoom_genes$plot,
    zoom_rao$plot, cre_strip, profiles, molecules, ncol = 1, align = "v", axis = "lr", rel_heights = heights)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  file = file.path(output_dir, paste0("LCL_hsca_", gsub("[^A-Za-z0-9_.-]+", "_", label), "_",
    paste(short_ids, collapse = "__"), "_domain_CTCF.pdf"))
  ggsave(file, combined, width = 12, height = sum(heights), limitsize = FALSE)
  invisible(list(pair = pair, allele = allele_row, region = region, genotypes = genotypes$calls,
    group_sizes = group_sizes, het_base_rate = het_base_rate, n_excluded = n_excluded,
    assignments = a, result = result, profile = profile, context = context, combined = combined, file = file))
}
