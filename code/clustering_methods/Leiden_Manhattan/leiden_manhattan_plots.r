source("/project/spott/cshan/fiber-seq/code/parsing_functions/plotting_functions.r", local = TRUE)

# leiden_manhattan_plots.r
#
# Plot helpers for Leiden + Manhattan single-molecule clustering, shared by
# leiden_LCL.Rmd, marcophage_leiden_manhattan.Rmd, the Snakemake plot rule and
# the fire_frequency / autocorrelation notebooks. Generic builders (palettes,
# save_figure(), plot_stacked_proportion(), plot_group_profile(),
# plot_read_heatmap(), ...) live in code/parsing_functions/plotting_functions.r,
# sourced above. LCL-only report orchestration stays in leiden_LCL.Rmd.
#
#   report helpers          embed_report_png(), as_timepoint_factor()
#   plot inputs             cluster_composition_inputs(), met_fraction_inputs(),
#                           signal_profile_inputs(), category_composition_inputs()
#   LCL Fiber-seq figures   plot_smf_reads(), plot_fiberseq_profiles(), plot_umap(),
#                           plot_knn_graph(), plot_cluster_pie(),
#                           methylation_by_cluster(), save_fiberseq_plots() and
#                           their palettes, legends and tables
#   plot coordinates        plot_anchor(), plot_positions(), prepare_read_tracks()
#   footprints              extract_nucleosomes(), read_footprint_region(),
#                           cache_footprint_tracks(), region_footprints(),
#                           footprint_profiles(), signal_profiles()
#
# The bar panels are drawn on a fixed 0-1 axis so clusters and genes stay
# comparable; the line plot uses the data range, since that is the plot for
# reading profile shape.

suppressMessages({
  requireNamespace("ComplexHeatmap")
})

# ---------------------------------------------------------------------------
# Embed a saved PNG in a knitted HTML report. The PNG is copied into
# image_dir, which retains the preview until Pandoc embeds it in the HTML, and
# an HTML <figure> with the image and caption is printed. Call it from a
# results = "asis" chunk, e.g. wrapped as the save_figure() preview$embed
# callback.
#
# Inputs:
#   path      - existing, non-empty PNG file
#   label     - caption text, also the image alt text (HTML-escaped)
#   image_dir - folder that keeps the copied PNG; created if needed
# Output:
#   invisible NULL; copies the PNG to <image_dir>/plot-<random>.png and prints
#   the <figure> HTML to stdout
# ---------------------------------------------------------------------------
embed_report_png <- function(path, label, image_dir) {
  stopifnot(file.exists(path), file.info(path)$size > 0)
  dir.create(image_dir, recursive = TRUE, showWarnings = FALSE)
  image_path <- tempfile("plot-", tmpdir = image_dir, fileext = ".png")
  stopifnot(file.copy(path, image_path))
  label_attribute <- as.character(htmltools::htmlEscape(label, attribute = TRUE))
  label <- as.character(htmltools::htmlEscape(label))
  image_path <- as.character(htmltools::htmlEscape(normalizePath(image_path), attribute = TRUE))
  cat('\n\n<figure><img src="', image_path,
      '" alt="', label_attribute, '" style="max-width:100%;height:auto;" />',
      '<figcaption>', label, '</figcaption></figure>\n\n', sep = "")
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Turn sample_name / timepoint labels into a factor whose levels follow the
# palette order (LPS_0, LPS_5, ...), whatever the caller passed in, so plots
# and legends list timepoints in time order. Labels missing from the palette
# come after, in order of appearance; a factor is returned unchanged.
#
# Inputs:
#   x              - timepoint labels (character or factor)
#   timepoint_cols - named colour vector; its names give the level order
# Output:
#   factor the same length as x
# ---------------------------------------------------------------------------
as_timepoint_factor <- function(x, timepoint_cols = LEIDEN_TIMEPOINT_COLORS) {
  if (is.factor(x)) return(x)
  lv <- intersect(names(timepoint_cols), unique(as.character(x)))
  factor(as.character(x), levels = c(lv, setdiff(unique(as.character(x)), lv)))
}


# ---------------------------------------------------------------------------
# Tables for the cluster-composition-by-timepoint plots: cluster proportions
# within each timepoint (stacked bars drawn from the per-read table), and each
# timepoint's reads spread over the clusters (the proportions table). Needs
# dplyr's %>%.
#
# Inputs:
#   res            - clustering result; res$assignments needs cluster and
#                    sample_name (the timepoint)
#   timepoint_cols - named timepoint colours giving the timepoint order
#                    (as_timepoint_factor())
# Output:
#   list(assignments = res$assignments with sample_name as a factor in
#        timepoint order,
#        proportions = one row per cluster x timepoint: cluster, sample_name,
#        n (reads) and proportion (n / reads of that timepoint))
# ---------------------------------------------------------------------------
cluster_composition_inputs <- function(res, timepoint_cols = LEIDEN_TIMEPOINT_COLORS) {
  df <- res$assignments
  df$sample_name <- as_timepoint_factor(df$sample_name, timepoint_cols)
  prop <- df %>%
    dplyr::count(cluster, sample_name) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(proportion = n / sum(n)) %>%
    dplyr::ungroup()
  list(assignments = df, proportions = prop)
}


# ---- Fiber-seq read, allele and composition plots ----

# read haplotype colours (HP1 / HP2, unphased, and reads pooled without phasing)
LCL_HAPLOTYPE_COLORS <- c(HP1 = "#ADD8E6", HP2 = "#FFF2AE", unphased = "#999999", pooled = "#BBBBBB")


# colour of each read-level track: m6A, nucleosome and FiberHMM TF footprints
LCL_TRACK_COLORS <- c(m6A = "#800080", "ft_nuc_130-160bp" = "#4d4d4d",
  "FiberHMM_10-30bp" = "#fdae6b", "FiberHMM_40-60bp" = "#f16913",
  "FiberHMM_60-80bp" = "#a63603")
# legend label of each track
LCL_TRACK_LABELS <- c(m6A = "m6A", "ft_nuc_130-160bp" = "nucleosome footprint (130-160 bp)",
  "FiberHMM_10-30bp" = "TF footprint 10-30 bp", "FiberHMM_40-60bp" = "TF footprint 40-60 bp",
  "FiberHMM_60-80bp" = "TF footprint 60-80 bp")
# ---------------------------------------------------------------------------
# Look up the plot colour of each track, stopping if a track has none.
#
# Inputs:
#   tracks - track names, e.g. c("m6A", "ft_nuc_130-160bp")
#   colors - named colour vector (default LCL_TRACK_COLORS)
# Output:
#   named character vector of colours in the order of tracks, names = tracks
# ---------------------------------------------------------------------------
feature_colors <- function(tracks, colors = LCL_TRACK_COLORS) {
  stopifnot(all(tracks %in% names(colors)))
  colors[tracks]
}

# ---------------------------------------------------------------------------
# Single-molecule footprint plot (as TNF Fig 1): one row per read, faceted by
# cluster, with the read span as a grey line and its footprints (and m6A
# calls, if "m6A" is in tracks) as coloured rectangles. Colour bars left of
# the reads mark each read's cluster, sample and, with include_haplotype, its
# focal SNP allele ("top_asfire_het" regions) or haplotype. x is relative to
# the region anchor (plot_anchor(): TSS, focal SNP or region centre); read
# rows and coordinates come from prepare_read_tracks().
#
# Inputs:
#   result              - clustering result: assignments (RID, cluster factor,
#                         start, end, sample_name, optional haplotype /
#                         allele_display), region (incl. annotation,
#                         region_type), optional group_id (added to the
#                         title), and allele_display_levels for
#                         "top_asfire_het" regions
#   records             - per-read intervals with RID, start (0-based BED),
#                         end and track, e.g. signal_profiles()$records;
#                         tracks not in `tracks` are dropped
#   sample_colors       - colours named by sample label (sample_name up to its
#                         first "_", or the sample_label_column values)
#   tracks              - tracks to draw, bottom to top (default
#                         LCL_FOOTPRINT_TRACKS: footprints only, no m6A)
#   include_haplotype   - add the allele / haplotype bar and its legend
#   track_colors, track_labels - colour and legend label per track
#   sample_label_column - optional assignments column holding sample labels
# Output:
#   ggplot object
# ---------------------------------------------------------------------------
plot_smf_reads <- function(result, records, sample_colors, tracks = LCL_FOOTPRINT_TRACKS,
                          include_haplotype = FALSE, track_colors = LCL_TRACK_COLORS,
                          track_labels = LCL_TRACK_LABELS, sample_label_column = NULL) {
  records <- records[records$track %in% tracks, , drop = FALSE]
  prepared <- prepare_read_tracks(result, records, sample_colors, sample_label_column)
  reads <- prepared$reads
  features <- prepared$features[prepared$features$track %in% tracks, , drop = FALSE]
  anchor <- prepared$anchor
  width <- max(1, anchor$right - anchor$left)
  feature_colors <- feature_colors(tracks, track_colors)
  # Draw nucleosomes first and TF footprints on top, as in TNF Fig 1.
  features <- features[order(match(features$track, tracks)), ]
  features$ymin <- features$row - 0.45
  features$ymax <- features$row + 0.45
  features$fill <- unname(feature_colors[features$track])
  cluster_colors <- cluster_id_palette(levels(reads$cluster))
  # colour bars left of the window: cluster, sample, then (optional) allele
  make_bar <- function(offset, colors) data.frame(cluster = reads$cluster, row = reads$row,
    left = anchor$left - offset * width, right = anchor$left - (offset - 0.025) * width,
    fill = unname(colors))
  bars <- rbind(make_bar(0.075, cluster_colors[as.character(reads$cluster)]),
                 make_bar(0.040, sample_colors[reads$sample_label]))
  allele_colors <- if (identical(result$region$region_type, "top_asfire_het"))
    fiberseq_category_palette(result$allele_display_levels) else LCL_HAPLOTYPE_COLORS
  allele <- if (identical(result$region$region_type, "top_asfire_het")) as.character(reads$allele_display) else reads$haplotype
  if (include_haplotype) bars <- rbind(bars, make_bar(0.110, allele_colors[allele]))
  legend_colors <- c(feature_colors, sample_colors,
    if (include_haplotype) allele_colors)
  legend_labels <- c(unname(track_labels[tracks]), names(sample_colors),
    if (include_haplotype) names(allele_colors))
  counts <- table(reads$cluster)
  cluster_labels <- setNames(paste0(names(counts), " (n=", as.integer(counts), ")"), names(counts))
  ticks <- pretty(c(anchor$left, anchor$right), n = 6)
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = reads, ggplot2::aes(x = left, xend = right, y = row, yend = row),
      color = "grey80", linewidth = 0.25) +
    ggplot2::geom_rect(data = features,
      ggplot2::aes(xmin = left, xmax = right, ymin = ymin, ymax = ymax, fill = fill)) +
    ggplot2::geom_rect(data = bars,
      ggplot2::aes(xmin = left, xmax = right, ymin = row - 0.5, ymax = row + 0.5, fill = fill)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
    ggplot2::scale_fill_identity(guide = "legend", name = NULL,
      breaks = unname(legend_colors), labels = legend_labels) +
    ggplot2::scale_x_continuous(limits = c(anchor$left - (if (include_haplotype) 0.12 * width else 0.085 * width),
      anchor$right + 0.5), breaks = ticks[ticks >= anchor$left & ticks <= anchor$right], expand = c(0, 0)) +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::facet_grid(cluster ~ ., scales = "free_y", space = "free_y",
      labeller = ggplot2::as_labeller(cluster_labels)) +
    ggplot2::labs(x = anchor$x_label, y = "SMF reads",
      title = paste(c(result$region$annotation, result$group_id), collapse = " | "),
      subtitle = "One row per retained molecule; saved m6A-defined Leiden clusters") +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4, byrow = TRUE)) +
    cowplot::theme_cowplot(font_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(angle = 0),
      panel.spacing.y = grid::unit(2, "mm"), legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8), legend.key.size = grid::unit(3, "mm"))
}

# ---------------------------------------------------------------------------
# Long table of the per-cluster mean m6A of every clustering feature, for the
# line plot with all features of every cluster overlaid in one panel.
# Promoter positions become strand-oriented positions from the TSS
# (plot_positions()); an optional centred rolling mean smooths each cluster.
#
# Inputs:
#   res      - clustering result: profiles (clusters x features matrix,
#              rownames = clusters, colnames = genomic feature positions) and
#              assignments$cluster (for the read counts)
#   smooth_k - rolling-mean width in features; <= 1 means no smoothing
#   region   - region row used by plot_anchor() (default res$region)
# Output:
#   data.frame with cluster (factor, levels relabelled "<cluster> (n=<reads>)"),
#   pos and value; rows with NA value (e.g. smoothing edges) dropped
# ---------------------------------------------------------------------------
met_fraction_inputs <- function(res, smooth_k = 1, region = res$region) {
  P <- res$profiles
  pos <- as.numeric(colnames(P))
  if (plot_anchor(region)$promoter) pos <- plot_positions(pos, region)
  smooth_row <- function(v) {
    if (smooth_k <= 1) return(v)
    as.numeric(stats::filter(v, rep(1 / smooth_k, smooth_k), sides = 2))
  }

  df <- do.call(rbind, lapply(rownames(P), function(cl)
    data.frame(cluster = cl, pos = pos, value = smooth_row(P[cl, ]))))
  df$cluster <- factor(df$cluster, levels = rownames(P))
  n <- table(res$assignments$cluster)
  levels(df$cluster) <- sprintf("%s (n=%d)", rownames(P), as.integer(n[rownames(P)]))

  df[!is.na(df$value), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Prepare signal_profiles()$profiles for the per-cluster occupancy plot: keep
# the tracks that have a colour, add anchor-relative positions and set the
# cluster and track factor orders.
#
# Inputs:
#   profiles     - long table with cluster, pos, track, fraction, n_reads
#   region       - region row for plot_positions()
#   clusters     - cluster levels in display order
#   track_colors - named track colours; tracks not named are dropped and the
#                  order of the names sets the track levels
# Output:
#   profiles with relative_pos added, cluster and track as factors, sorted by
#   cluster, track and relative_pos
# ---------------------------------------------------------------------------
signal_profile_inputs <- function(profiles, region, clusters, track_colors = LCL_TRACK_COLORS) {
  profiles <- profiles[profiles$track %in% names(track_colors), , drop = FALSE]
  profiles$relative_pos <- plot_positions(profiles$pos, region)
  profiles$cluster <- factor(profiles$cluster, levels = clusters)
  profiles$track <- factor(profiles$track, levels = names(track_colors)[names(track_colors) %in% profiles$track])
  profiles <- profiles[order(profiles$cluster, profiles$track, profiles$relative_pos), ]
  profiles
}

# ---------------------------------------------------------------------------
# One "Dynamic" HCL colour per LCL sample label, ordered by the number in the
# label (AL2 before AL10). Labels without an "AL<N>" / "AL-<N>" number go last
# (with an as.integer() NA warning).
#
# Inputs:
#   sample_names  - sample names, e.g. "AL10_..." (repeats allowed); only used
#                   for the default sample_labels
#   sample_labels - labels to colour (default: sample_names up to the first "_")
# Output:
#   named character vector of colours, names = unique sample labels
# ---------------------------------------------------------------------------
sample_palette <- function(sample_names, sample_labels = sub("_.*$", "", sample_names)) {
  labels <- sample_labels
  labels <- unique(labels[order(as.integer(sub("^AL-?([0-9]+).*$", "\\1", labels)))])
  setNames(grDevices::hcl.colors(length(labels), "Dynamic"), labels)
}



# ---------------------------------------------------------------------------
# Per-read table for a "category composition within cluster" stacked bar,
# e.g. plot_stacked_proportion(x = "cluster", fill = "category"): adds each
# read's category, optionally relabelled, as a factor in palette order.
#
# Inputs:
#   res       - clustering result; res$assignments must contain `column`
#   column    - assignments column holding the category (e.g. "sample_name")
#   colors    - colours named by category; every (relabelled) category needs
#               one, and their order sets the factor levels
#   label_fun - function applied to the category values first (default
#               identity), e.g. to shorten sample names
# Output:
#   res$assignments with an added `category` factor column (one row per read)
# ---------------------------------------------------------------------------
category_composition_inputs <- function(res, column, colors, label_fun = identity) {
  assignments <- res$assignments
  stopifnot(column %in% names(assignments))
  categories <- label_fun(as.character(assignments[[column]]))
  stopifnot(length(categories) == nrow(assignments),
            !anyNA(categories), all(categories %in% names(colors)))
  assignments$category <- factor(categories, levels = names(colors))
  assignments
}



# ---------------------------------------------------------------------------
# Colours for the sample / allele categories of the Fiber-seq panels.
# Categories are sorted, with "unphased" and "Unknown allele" moved last and
# coloured dark grey; focal SNP alleles ("rs123: A") get fixed colour-blind-
# safe colours (up to four); all others get "Dynamic" HCL colours.
#
# Inputs:
#   categories - category labels (repeats allowed)
# Output:
#   named character vector of colours, names = sorted unique categories
# ---------------------------------------------------------------------------
fiberseq_category_palette <- function(categories) {
  categories <- sort(unique(as.character(categories)))
  categories <- c(setdiff(categories, c("unphased", "Unknown allele")),
    intersect(c("unphased", "Unknown allele"), categories))
  colors <- setNames(grDevices::hcl.colors(length(categories), "Dynamic"), categories)
  focal <- categories[grepl("^rs[0-9]+: [ACGT]$", categories)]
  colors[focal] <- c("#0072B2", "#D55E00", "#CC79A7", "#009E73")[seq_along(focal)]
  unknown <- intersect(c("unphased", "Unknown allele"), categories)
  colors[unknown] <- "#555555"
  colors
}

# ---------------------------------------------------------------------------
# Standalone legend (coloured points) to stack under cowplot figures, taken
# from a throwaway ggplot with cowplot::get_legend().
#
# Inputs:
#   colors  - named colours, one legend key per name in that order
#   title   - legend title
#   columns - number of legend columns (filled by row)
#   labels  - key labels (default names(colors))
# Output:
#   legend grob (gtable) for cowplot::plot_grid()
# ---------------------------------------------------------------------------
fiberseq_legend <- function(colors, title, columns = 4L, labels = names(colors)) {
  d <- data.frame(category = factor(names(colors), levels = names(colors)), x = 1, y = 1)
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y, color = category)) +
    ggplot2::geom_point() + ggplot2::scale_color_manual(values = colors, labels = labels, drop = FALSE) +
    ggplot2::guides(color = ggplot2::guide_legend(title = title, title.position = "top", ncol = columns, byrow = TRUE,
      override.aes = list(size = 3))) + ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8))
  cowplot::get_legend(p)
}

# ---------------------------------------------------------------------------
# Height of a fiberseq_legend() grob plus 0.2 in of padding, used as its row
# height when stacking figure rows.
#
# Inputs:
#   legend - legend grob from fiberseq_legend()
# Output:
#   numeric height in inches
# ---------------------------------------------------------------------------
fiberseq_legend_height <- function(legend) {
  grid::convertHeight(sum(legend$heights), "in", valueOnly = TRUE) + .2
}

# ---------------------------------------------------------------------------
# Per-cluster profile + composition figure: one row per cluster with, on the
# left, the m6A call fraction (cluster colour) over the nucleosome occupancy
# (grey line and fill) and, on the right, horizontal bars of the sample and
# allele composition within the cluster. A header gives the title, colour key
# and focal markers; sample and allele legends go underneath. Promoters are
# drawn relative to the TSS (dashed line at 0), other regions in genomic
# coordinates; markers are dotted lines. Uses plot_group_profile() and
# plot_stacked_proportion() from plotting_functions.r.
#
# Inputs:
#   tables           - fiberseq_tables() output (groups, counts, profiles,
#                      sample, allele)
#   region           - region row (chr, analysis_start, analysis_end,
#                      region_type; strand and tss for promoters)
#   markers          - data.frame with pos (1-based genomic) and label, e.g.
#                      focal_markers(); zero rows draws none
#   cluster_colors   - colours named by cluster
#   sample_colors    - colours named by the sample categories of tables$sample
#   allele_colors    - colours named by the allele categories of tables$allele
#   title            - header title (default region$annotation)
#   nucleosome_track - profile track(s) drawn as nucleosome occupancy
#   nucleosome_label - its name in the header text
#   nucleosome_color, nucleosome_fill - line and ribbon colour of that track
#   sample_labels    - sample legend labels (default: names(sample_colors) up
#                      to the first "_")
# Output:
#   list(plot = cowplot figure, height = suggested height in inches: 0.85
#        header + 1.35 per cluster + the legend heights)
# ---------------------------------------------------------------------------
plot_fiberseq_profiles <- function(tables, region, markers, cluster_colors, sample_colors, allele_colors,
                                  title = region$annotation, nucleosome_track = "ft_nuc_130-160bp",
                                  nucleosome_label = "130-160 bp nucleosome occupancy",
                                  nucleosome_color = "#666666", nucleosome_fill = "#808080",
                                  sample_labels = sub("_.*$", "", names(sample_colors))) {
  anchor <- plot_anchor(region)
  limits <- if (anchor$promoter) c(anchor$left, anchor$right) else c(region$analysis_start, region$analysis_end)
  x_label <- if (anchor$promoter) anchor$x_label else paste0(region$chr, " (1-based bp)")
  if (anchor$promoter) markers$pos <- plot_positions(markers$pos, region)
  allele_heading <- if (identical(region$region_type, "top_asfire_het"))
    "Focal SNP allele" else "Allele / local haplotype"
  rows <- lapply(tables$groups, function(cluster) {
    d <- tables$profiles[tables$profiles$cluster == cluster, ]
    if (anchor$promoter) d$pos <- plot_positions(d$pos, region)
    n <- tables$counts$n_reads[match(cluster, tables$counts$cluster)]
    nuc <- d[d$track %in% nucleosome_track, , drop = FALSE]
    met <- d[d$track == "m6A", , drop = FALSE]
    # left panel: nucleosome ribbon + line, then the m6A line on top
    p <- plot_group_profile(d, style = "ribbon", mapping = ggplot2::aes(pos, fraction),
      layers = list(list(geom = "ribbon", data = nuc,
        mapping = ggplot2::aes(ymin = 0, ymax = fraction),
        fill = nucleosome_fill, alpha = .18, color = NA),
        list(geom = "line", data = nuc, color = nucleosome_color, linewidth = .55),
        list(geom = "line", data = met, color = cluster_colors[[cluster]], linewidth = .45))) +
      ggplot2::scale_x_continuous(limits = limits,
        labels = function(x) format(x, scientific = FALSE, trim = TRUE), expand = ggplot2::expansion(mult = 0)) +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
      ggplot2::labs(title = paste0(cluster, " (n = ", n, ")"), x = x_label, y = "Read fraction") +
      ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "none",
        plot.title = ggplot2::element_text(color = cluster_colors[[cluster]], face = "bold"))
    if (anchor$promoter) p <- p + ggplot2::geom_vline(xintercept = 0,
      color = "grey40", linetype = "dashed", linewidth = .35)
    if (nrow(markers)) {
      p <- p +
        ggplot2::geom_vline(data = markers, ggplot2::aes(xintercept = pos),
          inherit.aes = FALSE, color = "#666666", linetype = "dotted", linewidth = .35)
    }
    # right panels: sample and allele composition of this cluster
    composition_plots <- lapply(c("sample", "allele"), function(kind) {
      tab <- tables[[kind]]
      tab <- tab[tab$cluster == cluster, ]
      colors <- if (kind == "sample") sample_colors else allele_colors
      heading <- if (kind == "sample") "LCL sample" else allele_heading
      tab$category <- factor(tab$category, levels = names(colors))
      tab$bar <- 1
      plot_stacked_proportion(tab, "bar", "category", colors, y = "fraction",
        position = "stack", reverse = TRUE, horizontal = TRUE, width = .6,
        legend = list(drop = FALSE),
        scales = list(ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1),
          labels = scales::percent, expand = ggplot2::expansion(mult = 0)),
          ggplot2::scale_x_continuous(breaks = NULL)),
        labels = list(title = heading, subtitle = paste0("n = ", n), x = NULL, y = "Within cluster"),
        theme = ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "none",
          panel.grid = ggplot2::element_blank(), plot.margin = ggplot2::margin(5.5, 16, 5.5, 5.5)))
    })
    cowplot::plot_grid(p, composition_plots[[1]], composition_plots[[2]], nrow = 1,
      rel_widths = c(2.7, 1, 1.5), align = "h", axis = "tb")
  })
  marker_text <- if (nrow(markers)) paste(paste0(markers$label, ": ", markers$pos),
    collapse = "; ") else "No focal annotation resolved"
  header <- cowplot::ggdraw() + cowplot::draw_label(paste0(title,
    "\nCluster color: m6A calls; grey line / light grey fill: ", nucleosome_label, "\n",
    paste(strwrap(marker_text, width = 155), collapse = "\n")), size = 10)
  sample_legend <- fiberseq_legend(sample_colors, "LCL sample", 6, sample_labels)
  allele_legend <- fiberseq_legend(allele_colors, allele_heading, 3)
  heights <- c(.85, rep(1.35, length(rows)), fiberseq_legend_height(sample_legend), fiberseq_legend_height(allele_legend))
  plot <- cowplot::plot_grid(plotlist = c(list(header), rows, list(sample_legend, allele_legend)),
                            ncol = 1, rel_heights = heights)
  list(plot = plot, height = sum(heights))
}

# ---------------------------------------------------------------------------
# UMAP scatter of the reads coloured by one column (cluster, sample, allele).
#
# Inputs:
#   embedding   - one row per read with UMAP1, UMAP2 and `column`
#   column      - column mapped to colour
#   colors      - colours named by the levels of `column`; also the legend
#                 order (unused levels kept)
#   title       - plot title
#   labels      - legend labels (default names(colors))
#   show_legend - TRUE shows the legend on the right
# Output:
#   ggplot object
# ---------------------------------------------------------------------------
plot_umap <- function(embedding, column, colors, title, labels = names(colors),
                           show_legend = FALSE) {
  ggplot2::ggplot(embedding, ggplot2::aes(UMAP1, UMAP2, color = .data[[column]])) +
    ggplot2::geom_point(size = .9, alpha = .75) +
    ggplot2::scale_color_manual(values = colors, labels = labels, breaks = names(colors), drop = FALSE) +
    ggplot2::coord_equal() + ggplot2::labs(title = title, color = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = if (show_legend) "right" else "none")
}

# ---------------------------------------------------------------------------
# Draw the saved Manhattan KNN graph behind the Leiden clustering with a
# Fruchterman-Reingold layout: reads are points coloured by Leiden cluster or
# by focal SNP allele (knn_allele_summary()), edges are grey segments. The
# layout is seeded and the caller's random seed is restored afterwards.
#
# Inputs:
#   result      - clustering result: graph (undirected igraph, vertex names
#                 = the assignments RIDs, optional positive edge weights),
#                 assignments (RID, cluster; allele_label for "allele"),
#                 region$annotation, and params (seed; k_eff or k_neighbors
#                 and resolution, shown in the subtitle when present)
#   main        - title (default: region annotation + what the colour shows)
#   seed        - layout seed (default result$params$seed, else 1)
#   vertex_size - point size
#   edge_width, edge_alpha - edge line width and transparency
#   color_by    - "cluster" or "allele" (focal-SNP regions only)
# Output:
#   ggplot object; legend labels give the read count of each colour group
# ---------------------------------------------------------------------------
plot_knn_graph <- function(result, main = NULL, seed = NULL, vertex_size = 1.5,
                           edge_width = 0.3, edge_alpha = 0.2, color_by = c("cluster", "allele")) {
  color_by <- match.arg(color_by)
  graph <- result$graph
  assignments <- result$assignments
  if (!igraph::is_igraph(graph)) stop("The clustering result must contain its saved KNN graph")
  stopifnot(all(c("RID", "cluster") %in% names(assignments)),
    !igraph::is_directed(graph), igraph::vcount(graph) > 0L,
    !anyNA(assignments$RID), !anyDuplicated(assignments$RID), !anyNA(assignments$cluster))
  read_ids <- igraph::V(graph)$name
  stopifnot(length(read_ids) == igraph::vcount(graph), !anyNA(read_ids),
    !anyDuplicated(read_ids), setequal(read_ids, as.character(assignments$RID)))
  assignments <- assignments[match(read_ids, assignments$RID), , drop = FALSE]
  clusters <- if (is.factor(assignments$cluster)) levels(droplevels(assignments$cluster)) else
    unique(as.character(assignments$cluster))
  colors <- cluster_id_palette(clusters)
  weights <- igraph::edge_attr(graph, "weight")
  if (length(weights)) stopifnot(all(is.finite(weights)), all(weights > 0))
  if (is.null(seed)) seed <- if (is.null(result$params$seed)) 1L else result$params$seed
  stopifnot(length(seed) == 1L, is.finite(seed))
  # save the global RNG state and restore it on exit
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  coordinates <- igraph::layout_with_fr(graph, weights = weights)
  nodes <- data.frame(RID = read_ids, layout_x = coordinates[, 1], layout_y = coordinates[, 2],
    cluster = factor(assignments$cluster, levels = clusters))
  nodes$color_group <- nodes$cluster
  color_levels <- clusters
  legend_title <- "Leiden cluster"
  if (color_by == "allele") {
    summary <- knn_allele_summary(result)
    color_levels <- summary$allele_levels
    nodes$color_group <- factor(assignments$allele_label, levels = color_levels)
    colors <- fiberseq_category_palette(color_levels)
    legend_title <- "Focal SNP allele"
  }
  endpoints <- igraph::as_edgelist(graph, names = FALSE)
  edges <- data.frame(from_x = coordinates[endpoints[, 1], 1],
    from_y = coordinates[endpoints[, 1], 2], to_x = coordinates[endpoints[, 2], 1],
    to_y = coordinates[endpoints[, 2], 2])
  counts <- table(nodes$color_group)
  labels <- paste0(color_levels, " (n=", as.integer(counts[color_levels]), ")")
  if (is.null(main)) main <- paste(c(result$region$annotation,
    paste("Manhattan KNN graph colored by", if (color_by == "allele") "focal SNP allele" else "Leiden cluster")),
    collapse = "\n")
  neighbors <- if (is.null(result$params$k_eff)) result$params$k_neighbors else result$params$k_eff
  settings <- c(paste(nrow(nodes), "reads"), paste(nrow(edges), "edges"),
    if (!is.null(neighbors)) paste0("k = ", neighbors),
    if (!is.null(result$params$resolution)) paste0("resolution = ", result$params$resolution))
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges,
      ggplot2::aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
      color = "grey60", linewidth = edge_width, alpha = edge_alpha) +
    ggplot2::geom_point(data = nodes, ggplot2::aes(layout_x, layout_y, color = color_group),
      size = vertex_size) +
    ggplot2::scale_color_manual(values = colors, breaks = color_levels, labels = labels, drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = main, subtitle = paste(settings, collapse = " | "), color = legend_title) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      plot.margin = ggplot2::margin(12, 12, 12, 12))
}

# ---------------------------------------------------------------------------
# Focal-SNP allele mixing in the saved KNN graph: the allele composition of
# all reads and of each cluster, and for every read how many of its graph
# neighbours carry the same or the opposite allele. Every read must have a
# resolved phased focal allele. Used by plot_knn_graph(color_by = "allele")
# and plot_knn_allele_connectivity().
#
# Inputs:
#   result - clustering result: graph (undirected, simple igraph with vertex
#            names = assignments RIDs), focal (one row: variant_id (rsID),
#            ref, alt) and assignments (RID, cluster, allele_label =
#            "<variant_id>: <base>", allele_status all
#            "phased_focal_genotype")
# Output:
#   list of
#     allele_levels - the two allele labels, sorted
#     composition   - group ("All reads", then clusters), allele, n_reads,
#                     total_reads (of the group) and fraction
#     reads         - per read: RID, cluster, allele, degree,
#                     same_allele_neighbors, opposite_allele_neighbors,
#                     same_fraction, opposite_fraction (NA without neighbours)
#     connectivity  - per allele: n_reads, n_reads_with_neighbors,
#                     n_isolated_reads, mean same_fraction and
#                     opposite_fraction over connected reads, and
#                     expected_same_fraction ((n_allele - 1) / (n - 1), the
#                     random-mixing baseline)
# ---------------------------------------------------------------------------
knn_allele_summary <- function(result) {
  graph <- result$graph
  assignments <- result$assignments
  focal <- result$focal
  stopifnot(igraph::is_igraph(graph), !igraph::is_directed(graph), igraph::is_simple(graph),
    nrow(assignments) > 1L, !is.null(focal), nrow(focal) == 1L,
    all(c("RID", "cluster", "allele_label", "allele_status") %in% names(assignments)),
    !anyNA(assignments$RID), !anyDuplicated(assignments$RID), !anyNA(assignments$cluster),
    !anyNA(assignments$allele_label), all(assignments$allele_status == "phased_focal_genotype"))
  alleles <- sort(paste0(focal$variant_id, ": ", c(focal$ref, focal$alt)))
  stopifnot(length(alleles) == 2L, length(unique(alleles)) == 2L,
    all(grepl("^rs[0-9]+: [ACGT]$", alleles)), all(assignments$allele_label %in% alleles))
  read_ids <- igraph::V(graph)$name
  stopifnot(length(read_ids) == igraph::vcount(graph), !anyNA(read_ids),
    !anyDuplicated(read_ids), setequal(read_ids, as.character(assignments$RID)))
  assignments <- assignments[match(read_ids, assignments$RID), , drop = FALSE]
  clusters <- if (is.factor(assignments$cluster)) levels(droplevels(assignments$cluster)) else
    unique(as.character(assignments$cluster))
  stopifnot(!"All reads" %in% clusters)
  allele <- factor(assignments$allele_label, levels = alleles)
  composition <- as.data.frame(table(group = factor(assignments$cluster, levels = clusters),
    allele = allele), stringsAsFactors = FALSE)
  names(composition)[3] <- "n_reads"
  overall <- data.frame(group = "All reads", allele = alleles, n_reads = as.integer(table(allele)))
  composition <- dplyr::bind_rows(overall, composition)
  composition$group <- factor(composition$group, levels = c("All reads", clusters))
  composition$allele <- factor(composition$allele, levels = alleles)
  composition$total_reads <- ave(composition$n_reads, composition$group, FUN = sum)
  composition$fraction <- composition$n_reads / composition$total_reads
  endpoints <- igraph::as_edgelist(graph, names = FALSE)
  # count every undirected edge from both of its ends
  source_nodes <- c(endpoints[, 1], endpoints[, 2])
  target_nodes <- c(endpoints[, 2], endpoints[, 1])
  degree <- tabulate(source_nodes, nbins = length(read_ids))
  same <- tabulate(source_nodes[allele[source_nodes] == allele[target_nodes]], nbins = length(read_ids))
  reads <- data.frame(RID = read_ids, cluster = assignments$cluster, allele = allele,
    degree = degree, same_allele_neighbors = same, opposite_allele_neighbors = degree - same)
  reads$same_fraction <- ifelse(degree > 0L, same / degree, NA_real_)
  reads$opposite_fraction <- ifelse(degree > 0L, (degree - same) / degree, NA_real_)
  connectivity <- dplyr::bind_rows(lapply(alleles, function(label) {
    selected <- reads$allele == label
    connected <- selected & reads$degree > 0L
    data.frame(allele = label, n_reads = sum(selected), n_reads_with_neighbors = sum(connected),
      n_isolated_reads = sum(selected & reads$degree == 0L),
      same_fraction = if (any(connected)) mean(reads$same_fraction[connected]) else NA_real_,
      opposite_fraction = if (any(connected)) mean(reads$opposite_fraction[connected]) else NA_real_,
      expected_same_fraction = if (any(selected)) (sum(selected) - 1L) / (nrow(reads) - 1L) else NA_real_)
  }))
  connectivity$allele <- factor(connectivity$allele, levels = alleles)
  list(allele_levels = alleles, composition = composition, reads = reads, connectivity = connectivity)
}


# ---------------------------------------------------------------------------
# Pie chart of the share of reads in each cluster, with percentages on slices
# of at least 3% (white or black text by slice brightness) and legend labels
# "<cluster>\nn=<reads> (<pct>%)"; matches the original save_fiberseq_plots()
# pie. Callers with an existing count summary (e.g. fiberseq_tables()$counts)
# can pass it in.
#
# Inputs:
#   result         - clustering result; assignments$cluster (factor) gives the
#                    counts when counts is NULL, and the total in the title
#   counts         - optional data.frame with cluster, n_reads and proportion;
#                    its row order sets the slice order
#   cluster_colors - colours named by cluster (default cluster_id_palette())
# Output:
#   ggplot object
# ---------------------------------------------------------------------------
plot_cluster_pie <- function(result, counts = NULL, cluster_colors = NULL) {
  groups <- if (is.null(counts)) levels(result$assignments$cluster) else as.character(counts$cluster)
  if (is.null(counts)) {
    counts <- data.frame(cluster = factor(groups, levels = groups),
      n_reads = as.integer(table(result$assignments$cluster)))
    counts$proportion <- counts$n_reads / sum(counts$n_reads)
  }
  counts$cluster <- factor(counts$cluster, levels = groups)
  counts$label <- sprintf("%s\nn=%d (%.1f%%)", counts$cluster, counts$n_reads, 100 * counts$proportion)
  if (is.null(cluster_colors)) cluster_colors <- cluster_id_palette(groups)
  # white text on dark slices (luma < 0.5), black otherwise
  rgb <- grDevices::col2rgb(cluster_colors) / 255
  text_colors <- setNames(ifelse(colSums(rgb * c(.299, .587, .114)) < .5, "white", "black"), names(cluster_colors))
  ggplot2::ggplot(counts, ggplot2::aes(x = "", y = proportion, fill = cluster)) +
    ggplot2::geom_col(width = 1, color = "white", position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(proportion >= .03, sprintf("%.1f%%", 100 * proportion), ""), color = cluster),
      position = ggplot2::position_stack(vjust = .5, reverse = TRUE), size = 3) +
    ggplot2::scale_color_manual(values = text_colors, guide = "none") +
    ggplot2::coord_polar(theta = "y") + ggplot2::scale_fill_manual(values = cluster_colors, labels = counts$label) +
    ggplot2::labs(title = paste0("Cluster proportions: ", nrow(result$assignments), " reads"), fill = NULL) +
    ggplot2::theme_void() + ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size = 8))
}

# ---------------------------------------------------------------------------
# Mean m6A call proportion of each cluster: m6A calls / (reads x observed m6A
# sites) over the cluster's reads in result$site_met_mat.
#
# Inputs:
#   result - clustering result: assignments (RID, cluster factor) and
#            site_met_mat (reads x observed m6A sites)
# Output:
#   data.frame, one row per cluster: cluster, n_reads, n_sites, n_m6a_calls,
#   denominator (n_reads x n_sites) and methylation_proportion (NA when the
#   denominator is 0)
# ---------------------------------------------------------------------------
methylation_by_cluster <- function(result) {
  dplyr::bind_rows(lapply(levels(result$assignments$cluster), function(cluster) {
    ids <- result$assignments$RID[result$assignments$cluster == cluster]
    m <- result$site_met_mat[ids, , drop = FALSE]
    n_calls <- sum(m > 0)
    denominator <- nrow(m) * ncol(m)
    data.frame(cluster = cluster, n_reads = nrow(m), n_sites = ncol(m),
      n_m6a_calls = n_calls, denominator = denominator,
      methylation_proportion = if (denominator) n_calls / denominator else NA_real_)
  }))
}

# ---------------------------------------------------------------------------
# Bar chart of the methylation_by_cluster() proportions: one bar per cluster
# on a fixed 0-1 axis, coloured with cluster_id_palette().
#
# Inputs:
#   summary - table from methylation_by_cluster()
#   region  - region row; region$annotation is used as the title
# Output:
#   ggplot object
# ---------------------------------------------------------------------------
plot_methylation_by_cluster <- function(summary, region) {
  summary$cluster <- factor(summary$cluster, levels = summary$cluster)
  ggplot2::ggplot(summary, ggplot2::aes(cluster, methylation_proportion, fill = cluster)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = cluster_id_palette(levels(summary$cluster)), guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.03))) +
    ggplot2::labs(x = "m6A-defined cluster", y = "Mean m6A call proportion", title = region$annotation,
      subtitle = "Mean across observed m6A sites; all cluster reads included") +
    cowplot::theme_cowplot(font_size = 10) + cowplot::panel_border()
}

# ---------------------------------------------------------------------------
# Build and save the Fiber-seq report figures for one region's clustering:
#   fiberseq_profiles_composition - plot_fiberseq_profiles()
#   fiberseq_cluster_proportions  - plot_cluster_pie()
#   fiberseq_umap_overview        - UMAPs by cluster, sample and allele, the
#                                   pie and legends
#   fiberseq_haplotype_example    - the overview above the profiles, one page
# The result first goes through fiberseq_display_result() and the tables come
# from fiberseq_tables(). include_read_panels also writes the LCL detail
# panels, using methylation_by_cluster() / plot_methylation_by_cluster()
# (above) and plot_read_heatmap() (plotting_functions.r).
#
# Inputs:
#   result              - clustering result: assignments, region, focal /
#                         variants, site_met_mat, feat_mat; for the detail
#                         panels also graph, params and n_clusters
#   footprints          - footprint intervals (RID, start, end, track), e.g.
#                         dplyr::bind_rows() of extract_nucleosomes() and
#                         region_footprints()$records
#   output_dir          - figures go to <output_dir>/plots/ (always created)
#   sample_colors       - colours named by full sample_name (every
#                         assignments$sample_name needs one)
#   embedding           - optional precomputed UMAP (RID in the same order as
#                         result$assignments, UMAP1, UMAP2); NULL runs
#                         fiberseq_umap()
#   example_only        - TRUE saves only fiberseq_haplotype_example
#   plot_writer         - optional function(plot, name, width, height) called
#                         instead of writing the four figures as PDFs
#   include_read_panels - TRUE also writes the detail panels listed below
# Output:
#   invisible list(tables = fiberseq_tables() output, embedding = assignments
#   + UMAP1, UMAP2, display_cluster). Writes <output_dir>/plots/<name>.pdf for
#   each figure above (unless plot_writer is given); include_read_panels adds
#   knn_graph.pdf, methylation_proportion_by_cluster.pdf,
#   heatmap_m6a_footprints.pdf, fig1_read_footprints.pdf,
#   fig2_occupancy_by_cluster.pdf and cluster_sample_composition.pdf in the
#   same folder (always with save_figure(), even with plot_writer)
# ---------------------------------------------------------------------------
save_fiberseq_plots <- function(result, footprints, output_dir, sample_colors,
                               embedding = NULL, example_only = FALSE,
                               plot_writer = NULL, include_read_panels = FALSE) {
  result <- fiberseq_display_result(result)
  allele_heading <- if (identical(result$region$region_type, "top_asfire_het"))
    "Focal SNP allele" else "Allele / sample-local haplotype"
  plot_dir <- file.path(output_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- fiberseq_tables(result, footprints)
  cluster_colors <- cluster_id_palette(tables$groups)
  allele_colors <- fiberseq_category_palette(result$allele_display_levels)
  stopifnot(all(result$assignments$sample_name %in% names(sample_colors)))
  profiles <- plot_fiberseq_profiles(tables, result$region, result$markers,
    cluster_colors, sample_colors, allele_colors,
    result$region$annotation)
  figures <- list(fiberseq_profiles_composition = list(plot = profiles$plot, width = 16, height = profiles$height))
  if (is.null(embedding)) embedding <- fiberseq_umap(result)
  stopifnot(identical(embedding$RID, result$assignments$RID))
  embedding <- data.frame(result$assignments, embedding[, c("UMAP1", "UMAP2")])
  embedding$display_cluster <- factor(result$assignments$cluster, levels = tables$groups)
  cluster_labels <- paste0(tables$groups, " (n=", tables$counts$n_reads, ")")
  pc <- plot_umap(embedding, "display_cluster", cluster_colors, "UMAP: cluster")
  ps <- plot_umap(embedding, "sample_name", sample_colors, "UMAP: LCL sample")
  pa <- plot_umap(embedding, "allele_display", allele_colors, paste0("UMAP: ", allele_heading))
  pie <- plot_cluster_pie(result, counts = tables$counts, cluster_colors = cluster_colors)
  figures$fiberseq_cluster_proportions <- list(plot = pie, width = 8, height = 5)
  legends <- list(fiberseq_legend(cluster_colors, "Cluster", 5, cluster_labels),
    fiberseq_legend(sample_colors, "LCL sample", 6, sub("_.*$", "", names(sample_colors))),
    fiberseq_legend(allele_colors, allele_heading, 3))
  lh <- vapply(legends, fiberseq_legend_height, numeric(1))
  heading <- cowplot::ggdraw() + cowplot::draw_label(result$region$annotation, size = 12)
  overview <- cowplot::plot_grid(plotlist = c(list(heading, cowplot::plot_grid(pc, ps, pa, pie, ncol = 2)), legends),
                                ncol = 1, rel_heights = c(.45, 8, lh))
  overview_height <- 8.45 + sum(lh)
  figures$fiberseq_umap_overview <- list(plot = overview, width = 16, height = overview_height)
  # One PDF page combines the same UMAP coordinates, cluster pie, aggregate
  # accessibility/nucleosomes, and within-cluster sample/actual-allele bars.
  example <- cowplot::plot_grid(overview, profiles$plot, ncol = 1,
    rel_heights = c(overview_height, profiles$height))
  figures$fiberseq_haplotype_example <- list(plot = example, width = 16,
    height = overview_height + profiles$height)
  for (name in names(figures)) {
    if (example_only && name != "fiberseq_haplotype_example") next
    figure <- figures[[name]]
    if (!is.null(plot_writer)) {
      plot_writer(figure$plot, name, figure$width, figure$height)
    } else {
      save_figure(figure$plot, file.path(plot_dir, paste0(name, ".pdf")),
        width = figure$width, height = figure$height, limitsize = FALSE, bg = "white")
    }
  }
  if (include_read_panels) {
    # Optional LCL detail panels previously dispatched by the notebook.
    # detail panels name samples by their short label (sample_name before "_")
    detail_colors <- setNames(unname(sample_colors), sub("_.*$", "", names(sample_colors)))
    save_figure(plot_knn_graph(result), file.path(plot_dir, "knn_graph.pdf"),
      width = 12, height = 9, bg = "white")
    signals <- signal_profiles(result, footprints)
    profiles <- signals$profiles
    records <- signals$records
    matched <- match(records$RID, result$assignments$RID)
    for (column in intersect(c("cluster", "sample_name", "haplotype", "phase_set", "group_id"),
                             names(result$assignments))) records[[column]] <- result$assignments[[column]][matched]
    summary <- methylation_by_cluster(result)
    save_figure(plot_methylation_by_cluster(summary, result$region), file.path(plot_dir, "methylation_proportion_by_cluster.pdf"),
      width = 8, height = 4.5)
    heatmap <- plot_read_heatmap(result, result$region, detail_colors,
      include_haplotype = TRUE, variants = result$variants,
      show_cluster_profiles = TRUE)
    save_figure(heatmap, file.path(plot_dir, "heatmap_m6a_footprints.pdf"), width = 14,
      height = 11 + 0.55 * result$n_clusters,
      draw = function(x) ComplexHeatmap::draw(x, newpage = FALSE))
    save_figure(plot_smf_reads(result, records, detail_colors, include_haplotype = TRUE), file.path(plot_dir, "fig1_read_footprints.pdf"),
      width = 12, height = max(7, 0.03 * nrow(result$assignments) + 3.5 + 0.3 * result$n_clusters), limitsize = FALSE)
    # fig2: m6A and footprint occupancy profiles, one facet per cluster
    save_figure(local({
      .profile_args <- list(profiles = profiles,
      region = result$region,
      clusters = levels(result$assignments$cluster))
      profiles <- .profile_args$profiles
      region <- .profile_args$region
      clusters <- .profile_args$clusters
      title <- region$annotation
      track_colors <- LCL_TRACK_COLORS
      track_labels <- LCL_TRACK_LABELS
      nucleosome_track <- "ft_nuc_130-160bp"
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
    }), file.path(plot_dir, "fig2_occupancy_by_cluster.pdf"),
      width = 10, height = 1.4 * result$n_clusters + 2.5, limitsize = FALSE)
    # sample composition within each cluster (filled bars)
    composition <- category_composition_inputs(result, "sample_name", detail_colors,
      label_fun = function(x) sub("_.*$", "", x))
    save_figure(plot_stacked_proportion(composition, "cluster", "category", detail_colors,
      position = "fill", width = .75, legend = list(drop = FALSE),
      scales = list(ggplot2::scale_x_discrete(labels = function(x) paste0(x, "\n(n=", table(composition$cluster)[x], ")")),
        ggplot2::scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, .25),
          expand = ggplot2::expansion(mult = c(0, .02)))),
      labels = list(x = "Cluster", y = "Sample composition within cluster", fill = "LCL sample",
        title = result$region$annotation, subtitle = "Each bar sums to 100% of retained reads in that cluster"),
      theme = ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "right")), file.path(plot_dir, "cluster_sample_composition.pdf"),
      width = 12, height = 5)
  }
  invisible(list(tables = tables, embedding = embedding))
}


# ---------------------------------------------------------------------------
# Per-cluster mean m6A call at every bp of the region. Columns of met_mat are
# the positions with at least one call; every other bp is 0 because all reads
# are full-span and have no call there.
#
# Inputs:
#   res     - clustering result: assignments (RID, cluster factor) and the
#             region bounds, region$analysis_start / analysis_end (else
#             params$region_start / region_end)
#   met_mat - reads x positions m6A call matrix (rownames = RIDs, colnames =
#             1-based positions inside the region), e.g. res$site_met_mat
# Output:
#   data.frame, one row per cluster x bp: cluster, pos, met (mean call over
#   the cluster's reads, 0-1) and n_reads (cluster reads found in met_mat)
# ---------------------------------------------------------------------------
cluster_site_profiles <- function(res, met_mat) {
  M  <- as.matrix(met_mat)
  bounds <- if (!is.null(res$region$analysis_start)) c(res$region$analysis_start, res$region$analysis_end) else
    c(res$params$region_start, res$params$region_end)
  positions <- seq.int(bounds[1], bounds[2])
  at <- match(as.integer(colnames(M)), positions)
  stopifnot(!anyNA(at))
  df <- res$assignments
  df <- df[df$RID %in% rownames(M), ]
  do.call(rbind, lapply(levels(df$cluster), function(cl) {
    rids <- df$RID[df$cluster == cl]
    site_met <- colMeans(M[rids, , drop = FALSE], na.rm = TRUE)
    site_met[is.nan(site_met)] <- NA
    met <- numeric(length(positions))
    met[at] <- site_met
    data.frame(cluster = cl, pos = as.numeric(positions), met = met,
               n_reads = length(rids), row.names = NULL)
  }))
}

# ---------------------------------------------------------------------------
# Path of each sample's fibertools `ft extract` BED file for one feature and
# chromosome.
#
# Inputs:
#   sample_table - samples with fire_dir and sample_name columns
#   chromosome   - chromosome name, e.g. "chr5"
#   feature      - extracted feature, e.g. "nuc"
# Output:
#   character vector, one path per sample:
#   <fire_dir>/extracted_results/<feature>_by_chr/
#     <sample_name>.ft_extracted_<feature>.<chromosome>.bed.gz
# ---------------------------------------------------------------------------
extracted_path <- function(sample_table, chromosome, feature) {
  file.path(sample_table$fire_dir, "extracted_results", paste0(feature, "_by_chr"),
            paste0(sample_table$sample_name, ".ft_extracted_", feature, ".", chromosome, ".bed.gz"))
}

# ---------------------------------------------------------------------------
# Collapse sample-level variants to exact in-window SNP positions and REF/ALT
# bases, for the SNP marks under plot_read_heatmap(). Only a single-base REF
# with single-base ALT allele(s) counts; labels of several variants at one
# position are joined with "; ".
#
# Inputs:
#   variants - variant table with chr, pos, ref, alt and variant_id ("." or NA
#              is shown as "SNP"); NULL or empty gives no SNPs
#   region   - region row with chr, analysis_start and analysis_end
# Output:
#   data.frame with pos (integer) and label ("<id> <pos> <ref>><alt>"), one
#   row per SNP position; zero rows when there is none
# ---------------------------------------------------------------------------
window_snps <- function(variants, region) {
  empty <- data.frame(pos = integer(), label = character())
  if (is.null(variants) || !nrow(variants)) return(empty)
  snps <- variants[variants$chr == region$chr & variants$pos >= region$analysis_start &
    variants$pos <= region$analysis_end & nchar(variants$ref) == 1L &
    grepl("^[ACGT](,[ACGT])*$", variants$alt), , drop = FALSE]
  if (!nrow(snps)) return(empty)
  snps$id <- ifelse(is.na(snps$variant_id) | snps$variant_id == ".", "SNP", snps$variant_id)
  snps$label <- paste0(snps$id, " ", snps$pos, " ", snps$ref, ">", snps$alt)
  labels <- tapply(snps$label, snps$pos, function(x) paste(unique(x), collapse = "; "))
  data.frame(pos = as.integer(names(labels)), label = unname(labels))
}


# ---------------------------------------------------------------------------
# Position to mark on profile plots: only the observed focal SNP, else the
# region's focal SNP, else its annotated TSS.
#
# Inputs:
#   result - clustering result: optional focal (pos, variant_id, ref, alt) and
#            region (optional focal_pos, focal_snp, tss, gene)
# Output:
#   data.frame with pos, end (= pos) and label ("<variant_id> <ref>><alt>",
#   the focal_snp, or "<gene> TSS"); zero rows when nothing is annotated
# ---------------------------------------------------------------------------
focal_markers <- function(result) {
  focal <- result$focal
  if (!is.null(focal) && nrow(focal))
    return(data.frame(pos = focal$pos, end = focal$pos,
      label = paste0(focal$variant_id, " ", focal$ref, ">", focal$alt)))
  region <- result$region
  if (!is.null(region$focal_pos) && !is.na(region$focal_pos))
    return(data.frame(pos = region$focal_pos, end = region$focal_pos, label = region$focal_snp))
  if (!is.null(region$tss) && !is.na(region$tss))
    return(data.frame(pos = region$tss, end = region$tss, label = paste(region$gene, "TSS")))
  data.frame(pos = integer(), end = integer(), label = character())
}

# ---------------------------------------------------------------------------
# Add the display fields used by the Fiber-seq plots. Detailed phasing
# provenance is kept, but unknown focal alleles are combined for display: in a
# focal-SNP region (region$focal_pos set) every read without a resolved
# phased focal genotype becomes "Unknown allele", and the levels are the focal
# SNP alleles ("rs123: A", observed or REF/ALT in result$variants) plus
# "Unknown allele". "top_asfire_het" regions are strict: every read must be a
# resolved heterozygote (0|1 / 1|0) and the levels are exactly the region's
# REF and ALT alleles. Other regions show allele_label as is.
#
# Inputs:
#   result - clustering result: assignments (allele_label; in focal-SNP
#            regions also allele_status, and focal_genotype for
#            "top_asfire_het"), region, and optional focal, variants and
#            allele_levels
# Output:
#   result with markers (focal_markers()), allele_display_levels and
#   assignments$allele_display (factor) added; stops if a read's display
#   allele is not among the levels
# ---------------------------------------------------------------------------
fiberseq_display_result <- function(result) {
  result$markers <- focal_markers(result)
  a <- result$assignments
  a$allele_display <- a$allele_label
  anchored <- !is.null(result$region$focal_pos) && !is.na(result$region$focal_pos)
  if (anchored) {
    known <- a$allele_status == "phased_focal_genotype"
    strict <- identical(result$region$region_type, "top_asfire_het")
    if (strict && (!all(known) || !all(a$focal_genotype %in% c("0|1", "1|0"))))
      stop("Top AS-FIRE plots require only resolved phased heterozygous reads")
    a$allele_display[!known] <- "Unknown allele"
    bases <- a$allele_label[known]
    v <- result$variants
    if (!is.null(v) && nrow(v)) {
      v <- v[v$chr == result$region$chr & v$pos == result$region$focal_pos, , drop = FALSE]
      bases <- c(bases, paste0(if (is.null(result$region$focal_snp)) result$region$region_id else result$region$focal_snp, ": ",
        unique(c(v$ref, unlist(strsplit(v$alt, ",", fixed = TRUE))))))
    }
    bases <- sort(unique(bases[grepl("^rs[0-9]+: [ACGT]$", bases)]))
    if (strict) {
      bases <- sort(paste0(result$region$focal_snp, ": ", c(result$region$ref, result$region$alt)))
      stopifnot(length(unique(bases)) == 2L, all(a$allele_label %in% bases))
    }
    result$allele_display_levels <- if (strict) bases else c(bases, "Unknown allele")
  } else {
    result$allele_display_levels <- if (is.null(result$allele_levels))
      sort(unique(a$allele_label)) else result$allele_levels
  }
  a$allele_display <- factor(a$allele_display, levels = result$allele_display_levels)
  stopifnot(!anyNA(a$allele_display))
  result$assignments <- a
  result
}

# ---------------------------------------------------------------------------
# Summary tables behind the Fiber-seq figures: read counts per cluster, the
# sample and allele composition within each cluster (every composition
# denominator is the cluster read count), and per-bp profiles of the m6A call
# fraction (cluster_site_profiles()) and nucleosome occupancy
# (footprint_profiles()). Applies fiberseq_display_result() first.
#
# Inputs:
#   result           - clustering result: assignments (RID, cluster,
#                      sample_column, allele fields), site_met_mat and region
#                      (chr, start, end, analysis_start, analysis_end, width)
#   footprints       - footprint intervals with RID, start, end, track
#   nucleosome_track - footprint track(s) profiled as nucleosome occupancy
#   sample_column    - assignments column used for the sample composition
# Output:
#   list of
#     counts         - cluster, n_reads, total_reads, proportion
#     sample, allele - cluster, category, n_reads, cluster_reads, fraction
#                      (fractions sum to 1 within each cluster)
#     profiles       - cluster, pos (1-based), track ("m6A" or the nucleosome
#                      track), fraction, n_reads, chr, coordinate_system
#     groups         - cluster labels in display order (unused levels dropped)
# ---------------------------------------------------------------------------
fiberseq_tables <- function(result, footprints, nucleosome_track = "ft_nuc_130-160bp",
                            sample_column = "sample_name") {
  result <- fiberseq_display_result(result)
  a <- result$assignments
  groups <- if (is.factor(a$cluster)) levels(droplevels(a$cluster)) else unique(a$cluster)
  group <- factor(a$cluster, levels = groups)
  counts <- data.frame(cluster = groups, n_reads = as.integer(table(group)), total_reads = nrow(a))
  counts$proportion <- counts$n_reads / counts$total_reads
  composition <- function(column) {
    tab <- as.data.frame(table(cluster = group, category = a[[column]]), stringsAsFactors = FALSE)
    names(tab)[3] <- "n_reads"
    tab$cluster <- factor(tab$cluster, levels = groups)
    tab$cluster_reads <- counts$n_reads[match(tab$cluster, counts$cluster)]
    tab$fraction <- tab$n_reads / tab$cluster_reads
    stopifnot(all(abs(tapply(tab$fraction, tab$cluster, sum) - 1) < 1e-10))
    tab
  }
  view <- result
  view$assignments$cluster <- group
  nuc <- footprint_profiles(footprints, nucleosome_track, view$assignments, view$region)
  met <- cluster_site_profiles(view, view$site_met_mat)
  profiles <- dplyr::bind_rows(data.frame(cluster = met$cluster, pos = met$pos,
    track = "m6A", fraction = met$met, n_reads = met$n_reads), nuc)
  profiles$cluster <- factor(profiles$cluster, levels = groups)
  profiles$chr <- result$region$chr
  profiles$coordinate_system <- "1-based inclusive"
  stopifnot(all(is.finite(profiles$fraction)), all(profiles$fraction >= 0 & profiles$fraction <= 1))
  list(counts = counts, sample = composition(sample_column), allele = composition("allele_display"),
       profiles = profiles, groups = groups)
}

# ---------------------------------------------------------------------------
# 2-D UMAP of the reads from their clustering features with the Manhattan
# metric (as the Leiden KNN graph), random init and single threads so the
# seed reproduces it. Needs the uwot package; sets the global random seed.
#
# Inputs:
#   result      - clustering result: assignments$RID and feat_mat (reads x
#                 features, no NA, at least 3 reads)
#   seed        - set.seed() value
#   n_neighbors - UMAP neighbours (capped at the number of reads - 1)
#   min_dist    - UMAP min_dist
# Output:
#   result$assignments with UMAP1 and UMAP2 columns added
# ---------------------------------------------------------------------------
fiberseq_umap <- function(result, seed = 1L, n_neighbors = 15L, min_dist = 0.1) {
  stopifnot(requireNamespace("uwot", quietly = TRUE))
  a <- result$assignments
  x <- as.matrix(result$feat_mat[a$RID, , drop = FALSE])
  stopifnot(nrow(x) >= 3L, !anyNA(x), !anyDuplicated(a$RID))
  set.seed(seed)
  xy <- uwot::umap(x, metric = "manhattan", n_neighbors = min(n_neighbors, nrow(x) - 1L),
    min_dist = min_dist, n_threads = 1L, n_sgd_threads = 1L, init = "random", verbose = FALSE)
  data.frame(a, UMAP1 = xy[, 1], UMAP2 = xy[, 2])
}


# ---- Single-molecule footprint processing ----
# Process LCL nucleosome/TF footprints on fixed m6A-defined read clusters.
# BED12 parsing is provided by parsing_footprints_functions.r; display helpers follow.
# Footprint display tracks only; the clustering input remains the m6A matrix.
LCL_FOOTPRINT_TRACKS <- c("ft_nuc_130-160bp", "FiberHMM_10-30bp",
                          "FiberHMM_40-60bp", "FiberHMM_60-80bp")

# ---------------------------------------------------------------------------
# Turn every m6A call in site_met_mat into a 1-bp BED-style interval, in the
# same layout as the footprint records so both can be drawn together
# (signal_profiles(), plot_smf_reads()).
#
# Inputs:
#   result - clustering result: site_met_mat (reads x 1-based positions,
#            rownames = RIDs), assignments (RID, original_RID), region$chr
# Output:
#   data.frame, one row per call: RID, original_RID, chr, start (0-based),
#   end, size (1), track ("m6A")
# ---------------------------------------------------------------------------
m6a_intervals <- function(result) {
  sites <- Matrix::summary(as(result$site_met_mat, "dgCMatrix"))
  sites <- sites[sites$x > 0, , drop = FALSE]
  positions <- as.integer(colnames(result$site_met_mat)[sites$j])
  identifiers <- rownames(result$site_met_mat)[sites$i]
  data.frame(RID = identifiers,
    original_RID = result$assignments$original_RID[match(identifiers, result$assignments$RID)],
    chr = rep(result$region$chr, length(identifiers)), start = positions - 1L,
    end = positions, size = rep(1L, length(identifiers)), track = rep("m6A", length(identifiers)))
}

# ---------------------------------------------------------------------------
# Choose the x-axis origin and orientation for a region. A focal-SNP region is
# centred on the SNP and always retains genome order, even when promoter
# metadata is present. Only a promoter with a known +/- strand and TSS is
# TSS-anchored, flipped on the minus strand so upstream is negative. Any other
# region is centred on the middle of its window.
#
# Inputs:
#   region - region row: analysis_start / analysis_end (else start / end),
#            and optionally focal_pos, focal_snp, region_type, tss, strand
# Output:
#   list(anchor = origin position, direction = 1, or -1 for a minus-strand
#   promoter, left / right = window ends relative to the anchor, promoter =
#   TRUE if TSS-anchored, x_label = matching axis title)
# ---------------------------------------------------------------------------
plot_anchor <- function(region) {
  focal <- !is.null(region$focal_pos) && !is.na(region$focal_pos)
  promoter <- !focal && isTRUE(region$region_type == "promoter") &&
    !is.null(region$tss) && !is.na(region$tss) && isTRUE(region$strand %in% c("+", "-"))
  start <- if (!is.null(region$analysis_start)) region$analysis_start else region$start
  end <- if (!is.null(region$analysis_end)) region$analysis_end else region$end
  anchor <- if (focal) region$focal_pos else if (promoter) region$tss else if (length(c(start, end))) floor(mean(c(start, end))) else NA_real_
  direction <- if (promoter && isTRUE(region$strand == "-")) -1L else 1L
  bounds <- sort(direction * (c(start, end) - anchor))
  list(anchor = anchor, direction = direction, left = bounds[1], right = bounds[2],
    promoter = promoter,
    x_label = if (focal) paste0("Position relative to ", region$focal_snp, " (bp)") else if (promoter) "Position relative to canonical TSS (bp); upstream < 0" else "Position relative to region centre (bp)")
}

# ---------------------------------------------------------------------------
# Convert genomic positions to plot positions relative to the region anchor
# (plot_anchor()), strand-oriented for promoters. Positions are 1-based
# inclusive: convert BED starts with +1 before calling; BED exclusive ends
# already equal the last included 1-based position.
#
# Inputs:
#   pos    - 1-based genomic positions
#   region - region row (see plot_anchor())
# Output:
#   numeric vector of bp from the anchor (upstream < 0 for promoters)
# ---------------------------------------------------------------------------
plot_positions <- function(pos, region) {
  anchor <- plot_anchor(region)
  anchor$direction * (pos - anchor$anchor)
}

# ---------------------------------------------------------------------------
# Lay out reads and their intervals for plot_smf_reads(): reads are ordered by
# cluster, start and RID and numbered within their cluster, and reads and
# intervals get anchor-relative left / right ends clipped to the window
# (intervals widened by 0.5 bp on each side so 1-bp calls stay visible).
#
# Inputs:
#   result              - clustering result: assignments (RID, cluster, start,
#                         end, sample_name, optional haplotype) and region
#   records             - per-read intervals with RID, start (0-based BED),
#                         end; intervals of reads not in assignments dropped
#   sample_colors       - colours named by sample label; every read's label
#                         must be present
#   sample_label_column - optional assignments column with sample labels
#                         (default: sample_name up to its first "_")
# Output:
#   list(reads = assignments plus row, sample_label, haplotype ("pooled" if
#   absent), left, right; features = records plus row, cluster, left, right;
#   anchor = plot_anchor() of the region; region)
# ---------------------------------------------------------------------------
prepare_read_tracks <- function(result, records, sample_colors, sample_label_column = NULL) {
  anchor <- plot_anchor(result$region)
  reads <- result$assignments
  reads <- reads[order(reads$cluster, reads$start, reads$RID), , drop = FALSE]
  reads$row <- as.integer(ave(seq_len(nrow(reads)), reads$cluster, FUN = seq_along))
  reads$sample_label <- if (is.null(sample_label_column)) sub("_.*$", "", reads$sample_name) else
    reads[[sample_label_column]]
  stopifnot(length(reads$sample_label) == nrow(reads))
  stopifnot(all(reads$sample_label %in% names(sample_colors)))
  if (!"haplotype" %in% names(reads)) reads$haplotype <- "pooled"
  relative_start <- plot_positions(reads$start, result$region)
  relative_end <- plot_positions(reads$end, result$region)
  reads$left <- pmax(pmin(relative_start, relative_end), anchor$left)
  reads$right <- pmin(pmax(relative_start, relative_end), anchor$right)
  matched <- match(records$RID, reads$RID)
  records <- records[!is.na(matched), , drop = FALSE]
  matched <- matched[!is.na(matched)]
  records$row <- reads$row[matched]
  records$cluster <- reads$cluster[matched]
  relative_start <- plot_positions(records$start + 1L, result$region)
  relative_end <- plot_positions(records$end, result$region)
  records$left <- pmax(pmin(relative_start, relative_end) - 0.5, anchor$left - 0.5)
  records$right <- pmin(pmax(relative_start, relative_end) + 0.5, anchor$right + 0.5)
  list(reads = reads, features = records, anchor = anchor, region = result$region)
}


# ---------------------------------------------------------------------------
# Read the records of selected reads overlapping a region from an unindexed
# footprint BED. gzip/awk scans with bounded memory and emits only relevant
# records; every scanned record is checked against the column count of the
# explicitly configured format. Used by read_footprint_region().
#
# Inputs:
#   path         - footprint BED file (.gz files are decompressed with gzip)
#   format       - configured format name (footprint_format_columns())
#   region       - region row with chr, start (0-based) and end
#   original_ids - raw read names (BED column 4) to keep
# Output:
#   data.frame of the matching records with all-character columns V1, V2, ...;
#   an empty data.frame() when none match. Stops on a column-count mismatch
#   or read failure
# ---------------------------------------------------------------------------
stream_footprint_region <- function(path, format, region, original_ids) {
  expected <- footprint_format_columns(format)
  id_file <- tempfile("footprint_reads_", fileext = ".txt")
  out_file <- tempfile("footprint_region_", fileext = ".bed")
  err_file <- tempfile("footprint_stderr_", fileext = ".txt")
  on.exit(unlink(c(id_file, out_file, err_file)), add = TRUE)
  writeLines(unique(as.character(original_ids)), id_file)
  # awk: first file = wanted read names, then filter the BED arriving on stdin
  awk <- paste(
    'FILENAME == ARGV[1] { wanted[$0] = 1; next }',
    'NF != expected { printf "Footprint file %s (configured format %s): expected %d columns, found %d\\n", source_path, format_name, expected, NF > "/dev/stderr"; exit 23 }',
    '$1 == chrom && $2 < right && $3 > left && ($4 in wanted) { print }')
  reader <- if (grepl("\\.gz$", path)) paste("gzip -cd --", shQuote(path)) else paste("cat --", shQuote(path))
  command <- paste(reader, "| awk -F", shQuote("\t"),
    "-v", shQuote(paste0("expected=", expected)),
    "-v", shQuote(paste0("source_path=", path)), "-v", shQuote(paste0("format_name=", format)),
    "-v", shQuote(paste0("chrom=", region$chr)),
    "-v", shQuote(paste0("left=", format(region$start, scientific = FALSE, trim = TRUE))),
    "-v", shQuote(paste0("right=", format(region$end, scientific = FALSE, trim = TRUE))),
    shQuote(awk), shQuote(id_file), "-")
  status <- system2("bash", c("-o", "pipefail", "-c", shQuote(command)), stdout = out_file, stderr = err_file)
  if (status != 0L) stop("Unable to read ", path, " (configured format '", format, "'): ",
                        paste(readLines(err_file, warn = FALSE), collapse = "\n"), call. = FALSE)
  if (file.info(out_file)$size == 0) return(data.frame())
  read.delim(out_file, header = FALSE, quote = "", comment.char = "", colClasses = "character",
             stringsAsFactors = FALSE, check.names = FALSE)
}

# ---------------------------------------------------------------------------
# Read one sample's footprint records for a region: a tabix query when
# <path>.tbi exists (read_tabix_region()), otherwise a streamed scan
# (stream_footprint_region()). Checks the column count and BED coordinates,
# then keeps the records on region$chr that overlap the region and belong to
# the requested reads.
#
# Inputs:
#   path         - footprint BED (.bed.gz with optional .tbi index)
#   format       - configured format name (footprint_format_columns())
#   region       - region row with chr, start (0-based), end, analysis_start
#                  and analysis_end
#   original_ids - raw read names (BED column 4) to keep; empty returns none
# Output:
#   data.frame of raw BED columns (columns 2-3 integer, column 4 character);
#   an empty data.frame() when nothing matches. Stops on a missing file, a
#   wrong column count or invalid coordinates
# ---------------------------------------------------------------------------
read_footprint_region <- function(path, format, region, original_ids) {
  expected <- footprint_format_columns(format)
  if (!file.exists(path)) stop("Missing footprint file ", path, " (configured format '", format, "')")
  if (!length(original_ids)) return(data.frame())
  if (file.exists(paste0(path, ".tbi"))) {
    query <- GenomicRanges::GRanges(region$chr,
      IRanges::IRanges(region$analysis_start, region$analysis_end))
    records <- tryCatch(read_tabix_region(path, query), error = function(e)
      stop("Unable to read ", path, " (configured format '", format, "'): ", conditionMessage(e), call. = FALSE))
  } else {
    records <- stream_footprint_region(path, format, region, original_ids)
  }
  if (!nrow(records)) return(data.frame())
  if (ncol(records) != expected) footprint_column_error(path, format, ncol(records))
  start <- suppressWarnings(as.numeric(records[[2]]))
  end <- suppressWarnings(as.numeric(records[[3]]))
  if (any(!is.finite(start)) || any(!is.finite(end)) || any(start < 0) ||
      any(start != as.integer(start)) || any(end != as.integer(end)) || any(end < start))
    stop("Invalid BED coordinates in ", path, " (configured format '", format, "')")
  records[[2]] <- as.integer(start)
  records[[3]] <- as.integer(end)
  records[[4]] <- as.character(records[[4]])
  selected <- records[[1]] == region$chr & start < region$end & end > region$start &
    records[[4]] %in% original_ids
  records[selected, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Nucleosome footprints of the clustered reads in a region, one row per
# nucleosome block, filtered by size. input_path = NULL retains the original
# LCL fibertools path/reader: each sample's nuc BED12 (extracted_path(), tabix
# query, longest alignment per read), with reads matched as
# "<sample_name>::<read name>" RIDs. Otherwise input_path may be
# function(sample, chromosome), where sample is one row of sample_table, or a
# template with {sample_name}, {chr}, {chromosome}, {sample} and other
# sample-table columns; those files are read with read_footprint_region().
# Custom sources match their raw read names through original_RID when
# supplied. BED12 parsing comes from parsing_footprints_functions.r.
#
# Inputs:
#   sample_table       - one row per sample (sample_name; fire_dir for the
#                        default path)
#   region             - region row: chr, start (0-based), end,
#                        analysis_start, analysis_end
#   assignments        - clustered reads: RID, original_RID, optional
#                        sample_name
#   min_size, max_size - nucleosome size range in bp (inclusive); NULL drops
#                        that bound
#   input_path         - NULL, a function or a path template (see above)
#   format             - "bed12_fibertools" or "bed13_fiberhmm" (the latter
#                        needs input_path)
#   strict_blocks      - TRUE validates the blocks of bed13_fiberhmm files
# Output:
#   data.frame with RID, original_RID, start, end (0-based BED), size, chr and
#   track: "ft_nuc_<min>-<max>bp" (fibertools) or "nuc_<min>-<max>bp"
#   (FiberHMM) when both bounds are set, else "nuc_all", "nuc_gt<min-1>bp" or
#   "nuc_le<max>bp"
# ---------------------------------------------------------------------------
extract_nucleosomes <- function(sample_table, region, assignments, min_size = 130L, max_size = 160L,
                                input_path = NULL, format = "bed12_fibertools", strict_blocks = FALSE) {
  if (length(format) != 1L || !format %in% c("bed12_fibertools", "bed13_fiberhmm"))
    stop("Nucleosomes require bed12_fibertools or bed13_fiberhmm format")
  if (is.null(input_path) && format != "bed12_fibertools")
    stop("Supply input_path for bed13_fiberhmm nucleosomes")
  if (is.null(input_path)) paths <- extracted_path(sample_table, region$chr, "nuc")
  result <- lapply(seq_len(nrow(sample_table)), function(sample_index) {
    if (is.null(input_path)) {
      query <- GenomicRanges::GRanges(region$chr, IRanges::IRanges(region$analysis_start, region$analysis_end))
      bed <- read_ft_bed12(paths[sample_index], query, longest_alignment = TRUE)
      if (!nrow(bed)) return(NULL)
      bed$RID <- paste(sample_table$sample_name[sample_index], bed$RID, sep = "::")
      bed <- bed[bed$RID %in% assignments$RID, , drop = FALSE]
      if (!nrow(bed)) return(NULL)
      blocks <- convert_ft_bed12_to_bed6(bed,
        format = "bed12_fibertools", source = paths[sample_index])
      blocks$original_RID <- assignments$original_RID[match(blocks$RID, assignments$RID)]
    } else {
      sample <- as.data.frame(sample_table)[sample_index, , drop = FALSE]
      reads <- as.data.frame(assignments)
      if ("sample_name" %in% names(reads))
        reads <- reads[reads$sample_name == sample$sample_name, , drop = FALSE]
      if (!nrow(reads)) return(NULL)
      original_ids <- if ("original_RID" %in% names(reads)) as.character(reads$original_RID) else as.character(reads$RID)
      path <- if (is.function(input_path)) input_path(sample, region$chr) else {
        if (!is.character(input_path) || length(input_path) != 1L || is.na(input_path))
          stop("input_path must be a function or one path template")
        path <- input_path
        values <- c(as.list(sample), list(chr = region$chr, chromosome = region$chr,
          sample = sample$sample_name))
        for (field in names(values))
          path <- gsub(paste0("{", field, "}"), as.character(values[[field]]), path, fixed = TRUE)
        if (grepl("\\{[^{}]+\\}", path)) stop("Unresolved input_path template: ", path)
        path
      }
      if (length(path) != 1L || is.na(path)) stop("input_path must resolve to one file per sample and chromosome")
      bed <- read_footprint_region(path, format, region, original_ids)
      if (!nrow(bed)) return(NULL)
      blocks <- convert_ft_bed12_to_bed6(bed, format = format,
        longest_alignment = (format == "bed12_fibertools"),
        validate_blocks = strict_blocks && format == "bed13_fiberhmm", source = path)
      if (!nrow(blocks)) return(NULL)
      matched <- match(as.character(blocks$RID), original_ids)
      stopifnot(!anyNA(matched))
      blocks$original_RID <- as.character(blocks$RID)
      blocks$RID <- as.character(reads$RID[matched])
    }
    blocks$size <- blocks$end - blocks$start
    selected <- blocks$start < region$end & blocks$end > region$start
    if (!is.null(min_size)) selected <- selected & blocks$size >= min_size
    if (!is.null(max_size)) selected <- selected & blocks$size <= max_size
    blocks <- blocks[selected, , drop = FALSE]
    blocks[, c("RID", "original_RID", "start", "end", "size"), drop = FALSE]
  })
  result <- dplyr::bind_rows(result)
  if (!ncol(result)) result <- data.frame(RID = character(), original_RID = character(),
    start = integer(), end = integer(), size = integer())
  result$chr <- rep(region$chr, nrow(result))
  track <- if (is.null(min_size) && is.null(max_size)) "nuc_all" else
    if (is.null(max_size)) paste0("nuc_gt", min_size - 1L, "bp") else
    if (is.null(min_size)) paste0("nuc_le", max_size, "bp") else
    paste0(if (format == "bed12_fibertools") "ft_nuc_" else "nuc_", min_size, "-", max_size, "bp")
  result$track <- rep(track, nrow(result))
  result
}

# ---------------------------------------------------------------------------
# Pre-extract the pooled FiberHMM TF footprints (10-30, 40-60 and 60-80 bp
# BED4 files) that overlap any of the regions into one RDS cache per file, so
# per-region lookups (region_footprints()) avoid rescanning the large files.
# Files are scanned with gzip | awk in parallel; an existing cache is reused
# when it was built for the same regions and file.
#
# Inputs:
#   regions    - regions with chr, start (0-based) and end
#   tracks_dir - folder with one subfolder per chromosome holding
#                combined_<chr>_<size>bp_fps.bed.gz files
#   cache_dir  - output folder for the caches (created if needed)
#   workers    - parallel::mclapply cores
#   reuse      - FALSE always re-extracts
# Output:
#   character vector of cache paths, one per footprint file:
#   <cache_dir>/<footprint file name>.rds, each holding list(records =
#   data.frame(chr, start, end, original_RID, size, track =
#   "FiberHMM_<size>bp"), signature = list(version, regions, path))
# ---------------------------------------------------------------------------
cache_footprint_tracks <- function(regions, tracks_dir, cache_dir, workers = 2L, reuse = TRUE) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- unlist(lapply(unique(regions$chr), function(chromosome) {
    files <- list.files(file.path(tracks_dir, chromosome), pattern = "_(10-30|40-60|60-80)bp_fps\\.bed\\.gz$", full.names = TRUE)
    if (!length(files)) stop("No footprint tracks for ", chromosome)
    files
  }), use.names = FALSE)
  jobs <- lapply(paths, function(path) {
    chromosome <- basename(dirname(path))
    list(path = path, regions = regions[regions$chr == chromosome, c("chr", "start", "end")])
  })
  extract_job <- function(job) {
    signature <- list(version = 1L, regions = job$regions, path = job$path)
    cache <- file.path(cache_dir, paste0(basename(job$path), ".rds"))
    if (reuse && file.exists(cache)) {
      previous <- readRDS(cache)$signature
      previous_path <- if (!is.null(previous$path)) previous$path else previous$source$path
      if (identical(previous$version, signature$version) &&
          identical(previous$regions, signature$regions) &&
          identical(previous_path, signature$path)) return(cache)
    }
    conditions <- sprintf("($1 == \"%s\" && $2 < %.0f && $3 > %.0f)",
                            job$regions$chr, job$regions$end, job$regions$start)
    program <- paste0("BEGIN { FS = OFS = \"\\t\" } (", paste(conditions, collapse = " || "), ") { print }")
    temporary <- tempfile(tmpdir = cache_dir, fileext = ".bed")
    on.exit(unlink(temporary), add = TRUE)
    command <- paste("set -o pipefail; gzip -cd", shQuote(job$path), "| awk", shQuote(program), ">", shQuote(temporary))
    message("Extracting regional footprints: ", basename(job$path))
    status <- system2("/bin/bash", c("-c", shQuote(command)))
    if (status != 0L) stop("Footprint extraction failed: ", job$path)
    records <- data.frame(chr = character(), start = integer(), end = integer(), original_RID = character())
    if (file.info(temporary)$size > 0) {
      records <- data.table::fread(temporary, header = FALSE, data.table = FALSE)
      if (ncol(records) != 4L) stop("Expected BED4 footprints: ", job$path)
      names(records) <- c("chr", "start", "end", "original_RID")
    }
    records$size <- records$end - records$start
    records$track <- rep(sub("^combined_chr[^_]+_(.*)_fps\\.bed\\.gz$", "FiberHMM_\\1", basename(job$path)), nrow(records))
    saveRDS(list(records = records, signature = signature), cache)
    cache
  }
  caches <- parallel::mclapply(jobs, extract_job, mc.cores = workers, mc.preschedule = FALSE)
  if (any(vapply(caches, inherits, logical(1), "try-error"))) stop("One or more footprint extractions failed")
  unlist(caches, use.names = FALSE)
}

# ---------------------------------------------------------------------------
# One region's TF footprints from the cache_footprint_tracks() caches, mapped
# to the clustered reads. The pooled BED4 files carry no sample ID, so the
# raw read names must be unique across samples.
#
# Inputs:
#   caches      - cache paths from cache_footprint_tracks()
#   region      - region row with region_id, chr, start (0-based) and end
#   assignments - clustered reads with RID and original_RID (unique)
# Output:
#   list(records = footprints of the clustered reads overlapping the region:
#   chr, start, end, original_RID, size, track, RID; tracks = track name of
#   each cache of this chromosome, named by cache path)
# ---------------------------------------------------------------------------
region_footprints <- function(caches, region, assignments) {
  if (anyDuplicated(assignments$original_RID)) {
    stop("Combined BED4 lacks sample IDs: ambiguous original read names in ", region$region_id)
  }
  caches <- caches[startsWith(basename(caches), paste0("combined_", region$chr, "_"))]
  pieces <- lapply(caches, function(cache) {
    records <- readRDS(cache)$records
    selected <- records$start < region$end & records$end > region$start
    records <- records[selected, , drop = FALSE]
    matched <- match(records$original_RID, assignments$original_RID)
    records$RID <- assignments$RID[matched]
    records[!is.na(matched), , drop = FALSE]
  })
  list(records = dplyr::bind_rows(pieces),
       tracks = vapply(caches, function(cache) {
         sub("^combined_chr[^_]+_(.*)_fps\\.bed\\.gz.rds$", "FiberHMM_\\1", basename(cache))
       }, character(1)))
}

# ---------------------------------------------------------------------------
# Per-read, per-bp footprint occupancy: 1 where one of the read's intervals
# covers the bp, else 0. Intervals are clipped to the region.
#
# Inputs:
#   records     - intervals with RID, start (0-based BED) and end; reads not
#                 in assignments are skipped
#   assignments - reads (RID) giving the row order
#   region      - region row: start (0-based, = analysis_start - 1), end,
#                 analysis_start, analysis_end, width
# Output:
#   integer matrix, reads x bp (rownames = RIDs, colnames = 1-based positions
#   analysis_start..analysis_end)
# ---------------------------------------------------------------------------
occupancy_matrix <- function(records, assignments, region) {
  occupancy <- matrix(0L, nrow(assignments), region$width,
                       dimnames = list(assignments$RID, seq.int(region$analysis_start, region$analysis_end)))
  if (!nrow(records)) return(occupancy)
  for (record_index in seq_len(nrow(records))) {
    read_index <- match(records$RID[record_index], assignments$RID)
    if (is.na(read_index)) next
    # BED [start, end) -> window columns start + 1 .. end
    left <- max(records$start[record_index], region$start) - region$start + 1L
    right <- min(records$end[record_index], region$end) - region$start
    if (left <= right) occupancy[read_index, seq.int(left, right)] <- 1L
  }
  occupancy
}

# ---------------------------------------------------------------------------
# Per-cluster footprint occupancy profiles: for each track and cluster, the
# fraction of the cluster's reads covered at each bp (occupancy_matrix()).
#
# Inputs:
#   records     - footprint intervals with RID, start, end and track
#   tracks      - tracks to profile
#   assignments - reads with RID and cluster (factor; its levels are the
#                 clusters)
#   region      - region row (see occupancy_matrix())
# Output:
#   data.frame, one row per track x cluster x bp: cluster, pos (1-based),
#   track, fraction, n_reads
# ---------------------------------------------------------------------------
footprint_profiles <- function(records, tracks, assignments, region) {
  dplyr::bind_rows(lapply(tracks, function(track) {
    occupancy <- occupancy_matrix(records[records$track == track, ], assignments, region)
    dplyr::bind_rows(lapply(levels(assignments$cluster), function(cluster) {
      selected <- assignments$cluster == cluster
      data.frame(cluster = cluster, pos = as.integer(colnames(occupancy)), track = track,
                 fraction = colMeans(occupancy[selected, , drop = FALSE]), n_reads = sum(selected))
    }))
  }))
}


# ---------------------------------------------------------------------------
# Read-level records and per-cluster profiles of the m6A and footprint tracks
# of one region. Footprints use per-base occupancy (footprint_profiles()); m6A
# profiles are per-bp cluster means, 0 where no read has a call
# (cluster_site_profiles()).
#
# Inputs:
#   result     - clustering result: assignments (RID, original_RID, cluster),
#                site_met_mat and region
#   footprints - footprint intervals with RID, start, end, track
#   tracks     - footprint tracks to keep (default LCL_FOOTPRINT_TRACKS)
# Output:
#   list(records = m6a_intervals() plus the kept footprints, for
#   plot_smf_reads(); tracks = c("m6A", tracks); profiles = cluster, pos,
#   track, fraction, n_reads)
# ---------------------------------------------------------------------------
signal_profiles <- function(result, footprints, tracks = LCL_FOOTPRINT_TRACKS) {
  footprints <- footprints[footprints$track %in% tracks, , drop = FALSE]
  records <- dplyr::bind_rows(m6a_intervals(result), footprints)
  fp <- footprint_profiles(footprints, tracks, result$assignments, result$region)
  met <- cluster_site_profiles(result, result$site_met_mat)
  met <- data.frame(cluster = met$cluster, pos = met$pos, track = "m6A",
    fraction = met$met, n_reads = met$n_reads)
  list(records = records, tracks = c("m6A", tracks),
       profiles = dplyr::bind_rows(met, fp))
}
