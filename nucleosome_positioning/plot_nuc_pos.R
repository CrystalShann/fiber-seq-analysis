#' Plots for the Leiden-cluster nucleosome-positioning analysis.
#'
#' Every function takes the tables built by process_nuc_pos.R filtered to ONE
#' region, plus that region's row of the anchor table (`meta`, which carries
#' the window in relative coordinates), and returns a ggplot:
#'   plot_read_footprints()       Fig 1  read-level footprint tracks, rows grouped
#'                                       by Leiden cluster, then timepoint
#'   plot_occupancy_by_cluster()  Fig 2  nucleosome occupancy per cluster, one
#'                                       panel per cluster, timepoints pooled
#'   plot_position_summary()      supporting: -1 / +1 nucleosome position and
#'                                NFR width per read, by cluster and timepoint
#'
#' Cluster and timepoint colours mirror leiden_manhattan_plots.r
#' (LEIDEN_CLUSTER_COLORS / LEIDEN_TIMEPOINT_COLORS), so clusterN is coloured
#' exactly as in the accessibility heatmaps. TF footprint size classes share
#' one warm ramp (small -> large); nucleosomes are dark grey.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

# 25-colour cluster palette of the Leiden / topic-model notebooks
CLUSTER_COLORS <- c(
  "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1",
  "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
  "green1", "yellow4", "yellow3", "darkorange4", "brown"
)

# "clusterN" -> colour N of CLUSTER_COLORS (position in `levels` otherwise)
cluster_palette <- function(levels) {
  idx <- suppressWarnings(as.integer(sub("^cluster", "", levels)))
  if (anyNA(idx)) idx <- seq_along(levels)
  cols <- if (max(idx) <= length(CLUSTER_COLORS)) CLUSTER_COLORS[idx]
          else grDevices::colorRampPalette(CLUSTER_COLORS)(max(idx))[idx]
  setNames(cols, levels)
}

# timepoint colours keyed by LPS minutes (sequential, as in the Leiden notebook)
TP_COLORS <- c("0" = "#bdbdbd", "5" = "#6baed6", "10" = "#2171b5", "15" = "#08306b")
TP_LABELS <- setNames(paste0("LPS ", names(TP_COLORS), " min"), names(TP_COLORS))

# footprint classes drawn in the read tracks (NUC_ prefix: the topic-model
# helpers define their own FP_CLASS_COLORS)
NUC_FP_COLORS <- c(nucleosome = "#4d4d4d", "tf_10-30" = "#fdae6b",
                   "tf_40-60" = "#f16913", "tf_60-80" = "#a63603")
NUC_FP_LABELS <- c(nucleosome = "nucleosome footprint (FiberHMM)",
                   "tf_10-30" = "TF footprint 10-30 bp", "tf_40-60" = "TF footprint 40-60 bp",
                   "tf_60-80" = "TF footprint 60-80 bp")

# ---------------------------------------------------------------------------
# Shared pieces
# ---------------------------------------------------------------------------

position_label <- function(meta) {
  if (meta$region_type == "promoter") "Position relative to canonical TSS (bp)"
  else "Position relative to region centre (bp)"
}

region_title <- function(meta) {
  if (meta$region_type == "promoter")
    sprintf("%s (%s)   %s:%s-%s, %s strand", meta$gene, meta$gencode_name, meta$chr,
            format(meta$view_start, big.mark = ","), format(meta$view_end, big.mark = ","),
            meta$strand)
  else
    sprintf("%s   (%s:%s-%s)", meta$region_id, meta$chr,
            format(meta$view_start, big.mark = ","), format(meta$view_end, big.mark = ","))
}

# dashed line at the anchor (TSS / region centre)
anchor_line <- function() {
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.4)
}

# x axis = the region's window in relative coordinates; extra_left leaves room
# for the annotation bars of the read tracks
position_scale <- function(meta, extra_left = 0) {
  lo <- meta$rel_view_start; hi <- meta$rel_view_end
  br <- pretty(c(lo, hi), n = 6)
  scale_x_continuous(limits = c(lo - extra_left, hi), breaks = br[br >= lo & br <= hi],
                     expand = c(0, 0))
}

# profile plots: room for the outermost tick labels
profile_theme <- function() {
  theme_cowplot(font_size = 10) +
    theme(plot.margin = margin(5.5, 16, 5.5, 5.5))
}

# "cluster1" -> "cluster1 (n=45)" for strips and legends
cluster_labels <- function(reads_r) {
  n <- table(droplevels(reads_r$cluster))
  setNames(sprintf("%s (n=%d)", names(n), as.integer(n)), names(n))
}

# ---------------------------------------------------------------------------
# Fig 1. Read-level footprint tracks
# ---------------------------------------------------------------------------
# One row per read, footprints drawn at their full width and clipped to the
# window. Rows are grouped by Leiden cluster (facets, with a gap between
# clusters), then timepoint, then the position of the nearest downstream (+1)
# nucleosome. The two bars left of the window give each read's cluster and
# timepoint.
plot_read_footprints <- function(fp_r, reads_r, pos_r, meta) {
  lo <- meta$rel_view_start; hi <- meta$rel_view_end; w <- hi - lo
  ord <- merge(reads_r[, .(RID, cluster, timepoint, rel_read_start, rel_read_end)],
               pos_r[, .(RID, plus1_mid)], by = "RID", all.x = TRUE)
  ord[, cluster := droplevels(cluster)]
  setorder(ord, cluster, timepoint, plus1_mid, RID, na.last = TRUE)
  ord[, row := seq_len(.N), by = cluster]
  pal <- cluster_palette(levels(ord$cluster))

  spans <- ord[, .(cluster, row, x0 = pmax(rel_read_start, lo),
                   x1 = pmin(rel_read_end, hi))][x0 <= x1]
  fps <- merge(fp_r[rel_end >= lo & rel_start <= hi, .(RID, fp_class, rel_start, rel_end)],
               ord[, .(RID, cluster, row)], by = "RID")
  fps[, `:=`(xmin = pmax(rel_start - 0.5, lo), xmax = pmin(rel_end + 0.5, hi))]
  fps[, fp_class := factor(fp_class, levels = names(NUC_FP_COLORS))]
  setorder(fps, fp_class)                         # nucleosomes first, TF on top
  fps[, fill := NUC_FP_COLORS[as.character(fp_class)]]

  bars <- rbind(
    ord[, .(cluster, row, xmin = lo - 0.075 * w, xmax = lo - 0.045 * w,
            fill = unname(pal[as.character(cluster)]))],
    ord[, .(cluster, row, xmin = lo - 0.040 * w, xmax = lo - 0.010 * w,
            fill = unname(TP_COLORS[as.character(timepoint)]))])

  classes <- levels(droplevels(fps$fp_class))
  breaks  <- c(NUC_FP_COLORS[classes], TP_COLORS)
  labels  <- c(NUC_FP_LABELS[classes], TP_LABELS)

  ggplot() +
    geom_segment(data = spans, aes(x = x0, xend = x1, y = row, yend = row),
                 colour = "grey80", linewidth = 0.25) +
    geom_rect(data = fps, aes(xmin = xmin, xmax = xmax, ymin = row - 0.45, ymax = row + 0.45,
                              fill = fill)) +
    geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = row - 0.5, ymax = row + 0.5,
                               fill = fill)) +
    anchor_line() +
    scale_fill_identity(guide = "legend", name = NULL,
                        breaks = unname(breaks), labels = unname(labels)) +
    position_scale(meta, extra_left = 0.08 * w) +
    scale_y_reverse(expand = expansion(add = 0.6)) +
    facet_grid(cluster ~ ., scales = "free_y", space = "free_y",
               labeller = as_labeller(cluster_labels(reads_r))) +
    labs(x = position_label(meta), y = "reads", title = region_title(meta)) +
    theme_cowplot(font_size = 10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), strip.text.y = element_text(angle = 0),
          panel.spacing.y = unit(2, "mm"), legend.position = "bottom",
          legend.text = element_text(size = 8), legend.key.size = unit(3.5, "mm"))
}

# ---------------------------------------------------------------------------
# Fig 2. Nucleosome occupancy per cluster, one panel per cluster
# ---------------------------------------------------------------------------
plot_occupancy_by_cluster <- function(occ_cl, reads_r, meta) {
  d <- copy(occ_cl)[, cluster := droplevels(cluster)]
  ggplot(d, aes(pos, fraction, colour = cluster)) +
    anchor_line() +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = cluster_palette(levels(d$cluster)), guide = "none") +
    position_scale(meta) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       expand = expansion(mult = c(0, 0.02))) +
    facet_wrap(~ cluster, ncol = 1, labeller = as_labeller(cluster_labels(reads_r))) +
    labs(x = position_label(meta), y = "Fraction of reads with a nucleosome",
         title = "Nucleosome occupancy per Leiden cluster",
         subtitle = paste0(region_title(meta), ", timepoints pooled")) +
    profile_theme() +
    panel_border()
}

# ---------------------------------------------------------------------------
# Supporting: -1 / +1 nucleosome position and NFR width per read
# ---------------------------------------------------------------------------
plot_position_summary <- function(pos_r, reads_r, meta) {
  long <- melt(pos_r[, .(RID, cluster = droplevels(cluster), timepoint,
                         minus1_mid, plus1_mid, nfr_width)],
               id.vars = c("RID", "cluster", "timepoint"),
               variable.name = "metric", value.name = "value")[!is.na(value)]
  metric_labels <- c(minus1_mid = "-1 nucleosome midpoint (bp)",
                     plus1_mid  = "+1 nucleosome midpoint (bp)",
                     nfr_width  = "NFR width (bp)")
  ggplot(long, aes(timepoint, value)) +
    geom_boxplot(aes(fill = timepoint), outlier.shape = NA, width = 0.6, alpha = 0.7,
                 linewidth = 0.3) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.5, colour = "grey20") +
    facet_grid(metric ~ cluster, scales = "free_y",
               labeller = labeller(metric = metric_labels, cluster = cluster_labels(reads_r))) +
    scale_fill_manual(values = TP_COLORS, guide = "none") +
    labs(x = "LPS (min)", y = NULL,
         title = "Nearest upstream (-1) and downstream (+1) nucleosome per read",
         subtitle = region_title(meta)) +
    theme_cowplot(font_size = 9) +
    panel_border()
}
