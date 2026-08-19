#' Plot FiberHMM footprint distributions around CAGE TSS across LPS timepoints.
#'
#' Reads the TSVs written by process_nuc_pos.R and produces:
#'   * genome-wide footprint metaprofile around TSS, timepoints overlaid
#'     (raw and per-TSS normalized);
#'   * per selected-gene metaprofiles, faceted, timepoints overlaid;
#'   * nearest upstream/downstream footprint position distributions across
#'     timepoints (the TSS-relative shift);
#'   * nucleosome-to-nucleosome spacing (pairwise autocorrelation) over DHSs.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
})

RES_DIR  <- "/project/spott/cshan/fiber-seq/macrophage_project/nucleosome_positioning/results/tables"
PLOT_DIR <- "/project/spott/cshan/fiber-seq/macrophage_project/nucleosome_positioning/results/plots"
WINDOW_BP <- 2000L
TP_LEVELS <- c("0", "5", "10", "15")

# Blue -> red gradient over the LPS timeline
TP_COLORS <- c("0" = "#2166ac", "5" = "#67a9cf", "10" = "#ef8a62", "15" = "#b2182b")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

as_tp_factor <- function(x) factor(as.character(x), levels = TP_LEVELS)

kb_axis <- function() {
  scale_x_continuous(
    breaks = seq(-WINDOW_BP, WINDOW_BP, 500),
    labels = function(x) paste0(x / 1000, " kb")
  )
}

# ---------------------------------------------------------------------------
# 1. Genome-wide metaprofiles, timepoints overlaid
# ---------------------------------------------------------------------------
plot_genomewide_metaprofile <- function(meta_dt,
                                        which_mode = c("aggregated", "per_read"),
                                        which_norm = c("raw", "per_tss", "per_read")) {
  which_mode <- match.arg(which_mode)
  which_norm <- match.arg(which_norm)
  d <- meta_dt[mode == which_mode & norm == which_norm]
  d[, timepoint := as_tp_factor(timepoint)]

  ylab <- switch(which_norm,
    raw      = "Footprint midpoint count",
    per_tss  = "Footprint midpoints per TSS",
    per_read = "Mean footprints per read-TSS pair")
  ttl <- sprintf("FiberHMM footprints (130-160bp) around CAGE TSS [%s, %s]",
                 which_mode, which_norm)

  ggplot(d, aes(x = bin, y = N, color = timepoint)) +
    geom_line(linewidth = 0.4, alpha = 0.5) +
    geom_smooth(method = "loess", span = 0.05, se = FALSE, linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    scale_color_manual(values = TP_COLORS, name = "LPS (min)") +
    kb_axis() +
    labs(x = "Distance to TSS", y = ylab, title = ttl) +
    theme_bw(base_size = 12)
}

# Combined selected-gene metaprofile (all selected genes pooled), one panel per
# timepoint (same faceted style as the per-gene plots).
plot_combined_selected_metaprofile <- function(meta_comb_dt,
                                               which_mode = c("aggregated", "per_read"),
                                               which_norm = c("raw", "per_tss", "per_read")) {
  which_mode <- match.arg(which_mode)
  which_norm <- match.arg(which_norm)
  d <- meta_comb_dt[mode == which_mode & norm == which_norm]
  d[, timepoint := as_tp_factor(timepoint)]

  ylab <- switch(which_norm,
    raw      = "Footprint midpoint count",
    per_tss  = "Footprint midpoints per TSS",
    per_read = "Mean footprints per read-TSS pair")
  ttl <- sprintf("FiberHMM footprints (130-160bp) around combined selected-gene TSS [%s, %s]",
                 which_mode, which_norm)

  ggplot(d, aes(x = bin, y = N, color = timepoint)) +
    geom_line(linewidth = 0.4, alpha = 0.5) +
    geom_smooth(method = "loess", span = 0.05, se = FALSE, linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    scale_color_manual(values = TP_COLORS, name = "LPS (min)", guide = "none") +
    kb_axis() +
    facet_wrap(~ timepoint, nrow = 1,
               labeller = labeller(timepoint = function(x) paste0("LPS ", x, " min"))) +
    labs(x = "Distance to TSS", y = ylab, title = ttl) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---------------------------------------------------------------------------
# 2. Per-gene metaprofile, one panel per timepoint (same style as genome-wide)
# ---------------------------------------------------------------------------
plot_gene_metaprofile <- function(meta_sel_dt, which_gene,
                                  which_mode = c("aggregated", "per_read"),
                                  which_norm = c("raw", "per_tss", "per_read")) {
  which_mode <- match.arg(which_mode)
  which_norm <- match.arg(which_norm)
  d <- meta_sel_dt[gene == which_gene & mode == which_mode & norm == which_norm]
  d[, timepoint := as_tp_factor(timepoint)]

  ylab <- switch(which_norm,
    raw      = "Footprint midpoint count",
    per_tss  = "Footprint midpoints per TSS",
    per_read = "Mean footprints per read-TSS pair")

  ggplot(d, aes(x = bin, y = N, color = timepoint)) +
    geom_line(linewidth = 0.4, alpha = 0.5) +
    geom_smooth(method = "loess", span = 0.05, se = FALSE, linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    scale_color_manual(values = TP_COLORS, name = "LPS (min)", guide = "none") +
    kb_axis() +
    facet_wrap(~ timepoint, nrow = 1,
               labeller = labeller(timepoint = function(x) paste0("LPS ", x, " min"))) +
    labs(x = "Distance to TSS", y = ylab,
         title = sprintf("FiberHMM footprints around %s TSS", which_gene)) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---------------------------------------------------------------------------
# 3. Nearest up/down footprint distributions across timepoints
# ---------------------------------------------------------------------------
plot_nearest_updown <- function(nearest_dt) {
  d <- copy(nearest_dt)
  d[, timepoint := as_tp_factor(timepoint)]
  ggplot(d, aes(x = rel_pos, color = timepoint)) +
    geom_density(linewidth = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    scale_color_manual(values = TP_COLORS, name = "LPS (min)") +
    facet_wrap(~ direction, scales = "free") +
    labs(x = "Nearest footprint distance to TSS (bp, strand-aware)",
         y = "Density",
         title = "Nearest up/downstream footprint position vs TSS") +
    theme_bw(base_size = 12)
}

# Histogram of nearest-footprint sizes, split by up/downstream (count on y-axis).
# Timepoints overlaid as outlines so the four distributions are comparable.
plot_footprint_size_hist <- function(nearest_dt, binwidth = 10L, max_size = 1000L) {
  d <- copy(nearest_dt)
  d[, timepoint := as_tp_factor(timepoint)]
  d <- d[fp_size <= max_size]
  d[, direction := factor(direction, levels = c("upstream", "downstream"))]
  ggplot(d, aes(x = fp_size, color = timepoint)) +
    geom_histogram(binwidth = binwidth, fill = NA, position = "identity") +
    scale_color_manual(values = TP_COLORS, name = "LPS (min)") +
    facet_wrap(~ direction) +
    labs(x = "Footprint size (bp)", y = "Count",
         title = "Nearest up/downstream footprint size distribution") +
    theme_bw(base_size = 12)
}

# Median-shift point range across timepoints (compact view of the shift)
plot_nearest_shift <- function(shift_summary) {
  d <- copy(shift_summary)
  d[, timepoint := as_tp_factor(timepoint)]
  ggplot(d, aes(x = timepoint, y = median_rel_bp, color = direction, group = direction)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    labs(x = "LPS (min)", y = "Median nearest-footprint distance to TSS (bp)",
         title = "Nearest-footprint TSS shift across LPS timepoints",
         color = "Direction") +
    theme_bw(base_size = 12)
}

# ---------------------------------------------------------------------------
# 3b. Footprint size distribution (all footprints, NOT size-filtered)
# ---------------------------------------------------------------------------

# Histogram of FiberHMM footprint sizes, one panel per timepoint stacked
# vertically. Input is footprint_size_distribution.tsv.gz from
# process_nuc_pos.R::main_footprint_size_table(), which tallies every footprint
# with no size filter. The shaded band marks the 130-160bp mono-nucleosome
# window that every other analysis here filters to.
plot_footprint_size_distribution <- function(size_dt, binwidth = 5L,
                                             max_size = 1000L,
                                             normalize = FALSE) {
  d <- copy(size_dt)[fp_size <= max_size]
  d[, bin := (fp_size %/% binwidth) * binwidth]
  db <- d[, .(N = sum(N)), by = .(timepoint, bin)]
  db[, timepoint := as_tp_factor(timepoint)]
  db[, N := as.numeric(N)]   # counts come back from fread as integer; keep fractions exact
  if (normalize) db[, N := N / sum(N), by = timepoint]

  ylab <- if (normalize) "Fraction of footprints" else "Footprint count"
  ggplot(db, aes(x = bin, y = N, fill = timepoint)) +
    annotate("rect", xmin = FP_SIZE_MIN, xmax = FP_SIZE_MAX,
             ymin = -Inf, ymax = Inf, fill = "grey60", alpha = 0.3) +
    geom_col(width = binwidth) +
    scale_fill_manual(values = TP_COLORS, guide = "none") +
    facet_wrap(~ timepoint, ncol = 1, scales = "free_y",
               labeller = labeller(timepoint = function(x) paste0("LPS ", x, " min"))) +
    labs(x = sprintf("Footprint size (bp, %dbp bins)", binwidth), y = ylab,
         title = "FiberHMM footprint size distribution by timepoint",
         subtitle = sprintf("shaded = %d-%dbp mono-nucleosome window used elsewhere",
                            FP_SIZE_MIN, FP_SIZE_MAX)) +
    theme_bw(base_size = 12)
}

# ---------------------------------------------------------------------------
# 4. Nucleosome-to-nucleosome spacing (autocorrelation) over DHSs
# ---------------------------------------------------------------------------

# value: "N" (raw pair counts) or "N_norm" (fraction of pairs, depth-normalized)
plot_nuc_spacing <- function(spacing_dt, value = c("N_norm", "N")) {
  value <- match.arg(value)
  d <- copy(spacing_dt)
  d[, timepoint := as_tp_factor(timepoint)]
  d[, y := get(value)]

  ylab <- if (value == "N_norm") "Fraction of nucleosome pairs" else "Nucleosome pair count"

  ggplot(d, aes(x = dist_bin, y = y, color = timepoint)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = TP_COLORS, name = "LPS (min)") +
    scale_x_continuous(breaks = seq(0, max(d$dist_bin), 200)) +
    labs(x = "Nucleosome-to-nucleosome distance (bp)", y = ylab,
         title = "Nucleosome spacing over DHSs (pairwise autocorrelation)",
         subtitle = "Peaks at the nucleosome repeat length and its multiples = phased arrays") +
    theme_bw(base_size = 12)
}

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
main_plots <- function() {
  meta_all  <- fread(file.path(RES_DIR, "metaprofile_genomewide_by_timepoint.tsv"))
  meta_sel  <- fread(file.path(RES_DIR, "metaprofile_selected_genes_by_timepoint.tsv"))
  meta_comb <- fread(file.path(RES_DIR, "metaprofile_combined_selected_by_timepoint.tsv"))
  nearest   <- fread(file.path(RES_DIR, "nearest_updown_footprint_by_timepoint.tsv.gz"))
  shift     <- fread(file.path(RES_DIR, "nearest_updown_shift_summary.tsv"))
  spacing   <- fread(file.path(RES_DIR, "nuc_spacing_over_DHS_by_timepoint.tsv"))

  # Genome-wide: aggregated (raw + per-TSS) and per-read (raw + per-read-TSS).
  ggsave(file.path(PLOT_DIR, "metaprofile_genomewide_aggregated_raw.pdf"),
         plot_genomewide_metaprofile(meta_all, "aggregated", "raw"), width = 7, height = 4.5)
  ggsave(file.path(PLOT_DIR, "metaprofile_genomewide_aggregated_perTSS.pdf"),
         plot_genomewide_metaprofile(meta_all, "aggregated", "per_tss"), width = 7, height = 4.5)
  ggsave(file.path(PLOT_DIR, "metaprofile_genomewide_perRead_raw.pdf"),
         plot_genomewide_metaprofile(meta_all, "per_read", "raw"), width = 7, height = 4.5)
  ggsave(file.path(PLOT_DIR, "metaprofile_genomewide_perRead_perReadTSS.pdf"),
         plot_genomewide_metaprofile(meta_all, "per_read", "per_read"), width = 7, height = 4.5)

  # Combined selected genes (all pooled), one panel per timepoint.
  ggsave(file.path(PLOT_DIR, "metaprofile_combined_selected_aggregated_raw.pdf"),
         plot_combined_selected_metaprofile(meta_comb, "aggregated", "raw"), width = 11, height = 3.5)
  ggsave(file.path(PLOT_DIR, "metaprofile_combined_selected_aggregated_perTSS.pdf"),
         plot_combined_selected_metaprofile(meta_comb, "aggregated", "per_tss"), width = 11, height = 3.5)
  ggsave(file.path(PLOT_DIR, "metaprofile_combined_selected_perRead_raw.pdf"),
         plot_combined_selected_metaprofile(meta_comb, "per_read", "raw"), width = 11, height = 3.5)
  ggsave(file.path(PLOT_DIR, "metaprofile_combined_selected_perRead_perReadTSS.pdf"),
         plot_combined_selected_metaprofile(meta_comb, "per_read", "per_read"), width = 11, height = 3.5)

  ggsave(file.path(PLOT_DIR, "nearest_updown_density.pdf"),   plot_nearest_updown(nearest), width = 9, height = 4.5)
  ggsave(file.path(PLOT_DIR, "nearest_updown_shift.pdf"),     plot_nearest_shift(shift),    width = 6, height = 4.5)
  ggsave(file.path(PLOT_DIR, "nearest_updown_size_hist.pdf"), plot_footprint_size_hist(nearest), width = 9, height = 4.5)

  # All-footprint size distribution (unfiltered), one panel per timepoint.
  size_tsv <- file.path(RES_DIR, "footprint_size_distribution.tsv.gz")
  if (file.exists(size_tsv)) {
    fp_sizes <- fread(size_tsv)
    ggsave(file.path(PLOT_DIR, "footprint_size_distribution_by_timepoint.pdf"),
           plot_footprint_size_distribution(fp_sizes), width = 7, height = 9)
    # fraction version: timepoints directly comparable despite differing depths
    ggsave(file.path(PLOT_DIR, "footprint_size_distribution_by_timepoint_fraction.pdf"),
           plot_footprint_size_distribution(fp_sizes, normalize = TRUE),
           width = 7, height = 9)
  } else {
    message("Skipping footprint size distribution; run ",
            "process_nuc_pos.R::main_footprint_size_table() first.")
  }

  # One PDF per selected gene, timepoints as side-by-side panels:
  # aggregated raw, and per-read raw.
  gene_dir <- file.path(PLOT_DIR, "selected_genes")
  dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)
  for (g in sort(unique(meta_sel$gene))) {
    ggsave(file.path(gene_dir, sprintf("metaprofile_%s_aggregated_raw.pdf", g)),
           plot_gene_metaprofile(meta_sel, g, "aggregated", "raw"), width = 11, height = 3.5)
    ggsave(file.path(gene_dir, sprintf("metaprofile_%s_perRead_raw.pdf", g)),
           plot_gene_metaprofile(meta_sel, g, "per_read", "raw"), width = 11, height = 3.5)
  }

  # Nucleosome spacing over DHSs (depth-normalized and raw).
  ggsave(file.path(PLOT_DIR, "nuc_spacing_over_DHS_norm.pdf"),
         plot_nuc_spacing(spacing, "N_norm"), width = 8, height = 4.5)
  ggsave(file.path(PLOT_DIR, "nuc_spacing_over_DHS_raw.pdf"),
         plot_nuc_spacing(spacing, "N"), width = 8, height = 4.5)

  message("Wrote plots to ", PLOT_DIR)
}

# ===========================================================================
# Nucleosome occupancy + m6A topic heatmap, side by side, per promoter gene
# ---------------------------------------------------------------------------
# Combines, at a matched scale (same reads, same vertical order, same genomic
# range = the promoter interval the topic-model m6A heatmap covers):
#   * the m6A cluster heatmap (metA_mat, reads x A-positions, 0/1) reproduced
#     with the same darkgray/blue colours and order(clusters) row order as
#     1_m6a_promoter_topic_modelling.Rmd;
#   * a nucleosome-occupancy heatmap at single-base resolution (reads x bp,
#     binary: 1 if a 130-160bp FiberHMM footprint covers that base), same rows;
#   * a per-read nucleosome track (rectangles, cluster-coloured), same rows -
#     the plot_nuc_by_topic_cluster.py view rendered at the heatmap scale.
# Two figures per gene: track|m6A heatmap, and occupancy|m6A heatmap.
# ===========================================================================

FP_ROOT   <- "/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_footprint"
TOPIC_DIR <- "/project/spott/cshan/fiber-seq/macrophage_project/topic_model"
HM_DIR    <- file.path(PLOT_DIR, "nuc_topic_heatmap")

FP_SIZE_MIN <- 130L
FP_SIZE_MAX <- 160L

# colors_25[1:6] from the topic-model Rmd, the RowSideColors of the m6A heatmap.
CLUSTER_COLORS <- c("#1C86EE", "#E31A1C", "#008B00", "#6A3D9A", "#FF7F00", "#000000")

NUC_PROMOTERS <- data.table(
  gene  = c("IL1B","TNF","IL1A","PTGS2","IL8","CXCL2","IL6","RELB","rs2836882"),
  chr   = c("chr2","chr6","chr2","chr1","chr4","chr4","chr7","chr19","chr21"),
  start = c(112836323L,31575265L,112793233L,186679985L,73740006L,74098863L,22726376L,45001011L,39094468L),
  end   = c(112837385L,31575824L,112794969L,186681163L,73741332L,74099651L,22727094L,45002417L,39094712L)
)

# Read the topic-assignment TSV; add a heatmap row order = (cluster, table row).
# The TSV row order is the metA_mat / fit$L order, so ordering by (cluster, row)
# reproduces order(clusters) - the exact order the m6A heatmap stacks reads.
load_topic_order <- function(outname) {
  tsv <- file.path(TOPIC_DIR, outname, sprintf("read_topic_assignments_%s.tsv", outname))
  dt <- fread(tsv, colClasses = list(character = "RID"))
  dt[, table_row := .I]
  dt[, cluster_num := as.integer(sub("cluster", "", cluster))]
  dt <- dt[order(cluster_num, table_row)]
  dt[, .(RID, sample_name, cluster, cluster_num, table_row)]
}

# Filtered metA_mat (reads x A-positions, 0/1) ordered to match `rids`. Subset the
# saved region_met_data metA_mat to the topic reads (= the filter_met_mat rows) and
# drop all-zero columns, reproducing the heatmap matrix without refiltering.
load_metA_ordered <- function(outname, rids) {
  rd <- readRDS(file.path(TOPIC_DIR, outname, sprintf("region_met_data_%s.rds", outname)))
  m <- rd$metA_mat[rids, , drop = FALSE]
  m <- m[, Matrix::colSums(m) > 0, drop = FALSE]
  as.matrix(m)
}

# Per-read 130-160bp footprints (fp_start, fp_end) whose midpoint falls in the
# region, from each read's own sample BED. Returns a named list RID -> matrix(,2).
extract_read_footprints <- function(topic_dt, chr, start, end) {
  fp_by_rid <- list()
  for (sn in unique(topic_dt$sample_name)) {
    sdir <- gsub("_", "", sn)
    bed  <- file.path(FP_ROOT, sdir,
                      sprintf("%s_hmm_extracted_footprint_%s.bed.gz", sdir, chr))
    if (!file.exists(bed)) { message("  [WARN] missing ", bed); next }
    txt <- tryCatch(
      Rsamtools::scanTabix(
        Rsamtools::TabixFile(bed),
        param = GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))[[1]],
      error = function(e) character(0))
    if (!length(txt)) next
    rows <- data.table::fread(text = paste(txt, collapse = "\n"),
                              header = FALSE, select = c(2, 4, 11, 12),
                              col.names = c("V2", "V4", "V11", "V12"),
                              colClasses = list(character = c("V4", "V11", "V12")))
    want <- topic_dt[sample_name == sn, RID]
    rows <- rows[V4 %in% want]
    if (!nrow(rows)) next

    sizes_l  <- strsplit(rows$V11, ",", fixed = TRUE)
    starts_l <- strsplit(rows$V12, ",", fixed = TRUE)
    nblk <- lengths(sizes_l)
    sizes  <- as.integer(unlist(sizes_l,  use.names = FALSE))
    bstart <- as.integer(unlist(starts_l, use.names = FALSE))
    cstart <- rep(rows$V2, nblk)
    rid    <- rep(rows$V4, nblk)
    fp_start <- cstart + bstart
    fp_end   <- fp_start + sizes
    fp_mid   <- fp_start + sizes %/% 2L
    keep <- !is.na(sizes) & sizes >= FP_SIZE_MIN & sizes < FP_SIZE_MAX &
            fp_mid >= start & fp_mid <= end
    if (!any(keep)) next
    fp <- data.table(RID = rid[keep], fp_start = fp_start[keep], fp_end = fp_end[keep])
    for (r in unique(fp$RID))
      fp_by_rid[[r]] <- as.matrix(fp[RID == r, .(fp_start, fp_end)])
  }
  fp_by_rid
}

# Own-cluster footprint pile midpoints, written by main_own_pile_midpoints() in
# process_nuc_pos.R. Each read cluster is scanned on its own footprints alone, so
# a line drawn from this table describes only the band it is drawn in.
load_own_pile_midpoints <- function(gene_sym) {
  f <- file.path(RES_DIR, "nuc_own_pile_midpoints_by_topic_cluster.tsv")
  if (!file.exists(f)) {
    message("  [WARN] no ", basename(f), " -- run main_own_pile_midpoints(); ",
            "drawing the track without midpoint lines")
    return(data.table())
  }
  d <- fread(f)
  if (!nrow(d)) return(d)
  d[gene == gene_sym, .(cluster_num = as.integer(sub("cluster", "", read_cluster)),
                        mid_mean, n_fp, fold)]
}

# --- panel drawing (base graphics; row 1 of the ordered matrix at the top) ---

draw_sidebar <- function(cluster_num) {
  n <- length(cluster_num)
  image(z = matrix(cluster_num[n:1], nrow = 1), col = CLUSTER_COLORS,
        breaks = seq(0.5, length(CLUSTER_COLORS) + 0.5, by = 1),
        axes = FALSE, xlab = "", ylab = "")
  mtext("cluster", side = 1, line = 0.3, cex = 0.6)
}

draw_heat <- function(M, one_col, chr, start, end, title) {
  n <- nrow(M)
  image(z = t(M[n:1, , drop = FALSE]), col = c("darkgray", one_col),
        breaks = c(-0.5, 0.5, 1.5), axes = FALSE, xlab = "", ylab = "",
        useRaster = TRUE)
  axis(1, at = c(0, 1), labels = format(c(start, end), big.mark = ","), cex.axis = 0.6)
  mtext(sprintf("%s  (%s:%s-%s)", title, chr, format(start, big.mark = ","),
                format(end, big.mark = ",")), side = 3, line = 0.4, cex = 0.8)
}

draw_track <- function(rids, cluster_num, fp_by_rid, chr, start, end, title,
                       midlines = NULL) {
  n <- length(rids)
  plot(NA, xlim = c(start, end), ylim = c(0.5, n + 0.5), xlab = "", ylab = "",
       axes = FALSE, xaxs = "i", yaxs = "i")
  for (i in seq_along(rids)) {
    y <- n - i + 1L                       # read 1 at the top
    col <- CLUSTER_COLORS[cluster_num[i]]
    fps <- fp_by_rid[[rids[i]]]
    if (is.null(fps)) next
    rect(fps[, 1], y - 0.45, fps[, 2], y + 0.45, col = col, border = NA)
  }

  # Midpoint of each footprint pile, spanning only the rows of the read cluster
  # it was computed from. Reads are ordered by (cluster, table row) and drawn top
  # down, so a cluster's rows are contiguous and its y-range is just the range of
  # n - i + 1 over that cluster's row indices.
  if (!is.null(midlines) && nrow(midlines)) {
    for (k in unique(midlines$cluster_num)) {
      rows <- which(cluster_num == k)
      if (!length(rows)) next
      ys <- n - range(rows) + 1L
      for (x in midlines[cluster_num == k, mid_mean]) {
        # solid white halo first so the dashed line stays legible over the
        # footprint colour (a dashed halo would show the colour through the gaps)
        segments(x, min(ys) - 0.5, x, max(ys) + 0.5, col = "white", lwd = 3.2)
        segments(x, min(ys) - 0.5, x, max(ys) + 0.5, col = "black", lwd = 1.4,
                 lty = 2)
      }
    }
  }

  # White rule between heatmap clusters, drawn last so it cuts cleanly through
  # the footprints and the ends of the midpoint lines. Reads are ordered by
  # (cluster, table row), so each cluster is a contiguous block and a boundary
  # sits between rows i and i+1 wherever the cluster number changes.
  brk <- which(diff(cluster_num) != 0L)
  if (length(brk)) abline(h = n - brk + 0.5, col = "white", lwd = 3)

  axis(1, at = c(start, end), labels = format(c(start, end), big.mark = ","),
       cex.axis = 0.6)
  mtext(sprintf("%s  (%s:%s-%s)", title, chr, format(start, big.mark = ","),
                format(end, big.mark = ",")), side = 3, line = 0.4, cex = 0.8)
}

# One gene: per-read nucleosome track (with each read cluster's own footprint
# pile midpoints marked) beside the m6A cluster heatmap.
plot_gene_nuc_topic <- function(gene, chr, start, end) {
  outname <- sprintf("%s_%s_%d_%d", gene, chr, start, end)
  topic <- load_topic_order(outname)
  rids  <- topic$RID
  metA  <- load_metA_ordered(outname, rids)
  fp    <- extract_read_footprints(topic, chr, start, end)
  mids  <- load_own_pile_midpoints(gene)
  n     <- length(rids)
  with_nuc <- sum(vapply(rids, function(r) !is.null(fp[[r]]), logical(1)))
  cat(sprintf("%-10s %d reads (%d with a nucleosome), metA %d x %d, %d pile midpoints\n",
              gene, n, with_nuc, nrow(metA), ncol(metA), nrow(mids)))

  # panels: cluster sidebar | track | spacer | m6A heatmap
  ttl <- sprintf("%s  |  reads by topic cluster (heatmap order, n=%d)", gene, n)
  pdf(file.path(HM_DIR, sprintf("%s_track_vs_m6a.pdf", gene)),
      width = 11, height = max(4, 0.03 * n + 2))
  layout(matrix(1:4, nrow = 1), widths = c(0.45, 4, 0.35, 4))
  par(oma = c(2, 1, 3, 1))
  par(mar = c(3, 0.4, 2, 0.2)); draw_sidebar(topic$cluster_num)
  par(mar = c(3, 0.4, 2, 0.5))
  draw_track(rids, topic$cluster_num, fp, chr, start, end,
             "nucleosome footprints (130-160bp), line = pile midpoint", midlines = mids)
  par(mar = c(3, 0, 2, 0)); plot.new()
  par(mar = c(3, 0.5, 2, 0.5))
  draw_heat(metA, "blue", chr, start, end, "m6A (per-A methylation)")
  mtext(ttl, outer = TRUE, cex = 1)
  dev.off()
}

main_nuc_topic_heatmaps <- function() {
  dir.create(HM_DIR, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(NUC_PROMOTERS)))
    plot_gene_nuc_topic(NUC_PROMOTERS$gene[i], NUC_PROMOTERS$chr[i],
                        NUC_PROMOTERS$start[i], NUC_PROMOTERS$end[i])
  message("Wrote nucleosome/m6A side-by-side heatmaps to ", HM_DIR)
}


# ===========================================================================
# Distribution-level comparison of read clusters
# ===========================================================================
# Reads the tables written by main_cluster_distributions() in process_nuc_pos.R.
#
#   plot_cluster_intensity_profiles()  one panel per read cluster, that cluster's
#     footprint intensity in colour against the other clusters in grey. The curve
#     is footprints per READ per kb, not a probability density, so its HEIGHT is
#     occupancy and its POSITION is positioning -- both questions in one panel,
#     with no second y axis.
#   plot_cluster_signal_quadrant()  one point per gene, evidence for a positional
#     difference against evidence for an occupancy difference, so it is visible
#     that the two mostly pick out different genes.
#
# Six read clusters cannot be told apart by colour alone: CLUSTER_COLORS is the
# topic model's legacy palette and its green/red pair fails CVD separation
# (deutan dE 5.1). So identity lives in the panel title and colour only
# reinforces it -- never six overlaid curves distinguished by hue.

OCC_DIR <- "/project/spott/cshan/fiber-seq/macrophage_project/nucleosome_positioning/results/nuc_occpancy_by_topic_cluster"

# cluster name -> its colour in the m6A heatmap RowSideColors
cluster_color_map <- function(clusters) {
  cl <- sort(unique(grep("^cluster[0-9]+$", clusters, value = TRUE)))
  n  <- as.integer(sub("cluster", "", cl))
  setNames(CLUSTER_COLORS[((n - 1L) %% length(CLUSTER_COLORS)) + 1L], cl)
}

# ---------------------------------------------------------------------------
# 1. Footprint intensity per read cluster (spotlight small multiples)
# ---------------------------------------------------------------------------
plot_cluster_intensity_profiles <- function(prof_dt, which_gene, stats_dt = NULL) {
  d <- prof_dt[gene == which_gene]
  if (!nrow(d)) return(NULL)
  d <- copy(d)[, intensity_kb := intensity * 1000]        # per read per kb
  # background = every cluster's curve with the facet column dropped, so it is
  # redrawn in all panels
  bg <- d[, .(pos, intensity_kb, grp = read_cluster)]
  lab <- d[, .(lab = sprintf("%s  (n=%d reads, %.2f footprints/read)",
                             read_cluster[1], n_reads[1], fp_per_read[1])),
           by = read_cluster]
  d <- merge(d, lab, by = "read_cluster")
  d[, lab := factor(lab, levels = lab[order(read_cluster)][!duplicated(lab[order(read_cluster)])])]

  sub <- sprintf("%s:%s-%s", d$chrom[1],
                 format(d$view_start[1], big.mark = ","),
                 format(d$view_end[1], big.mark = ","))
  if (!is.null(stats_dt) && nrow(stats_dt[gene == which_gene]))
    sub <- paste0(sub, sprintf("   |   position q = %.3f, occupancy q = %.3f",
                               stats_dt[gene == which_gene, q_position],
                               stats_dt[gene == which_gene, q_occupancy]))

  ggplot(d, aes(pos, intensity_kb)) +
    geom_line(data = bg, aes(group = grp), color = "grey85", linewidth = 0.35) +
    geom_line(aes(color = read_cluster), linewidth = 0.9) +
    scale_color_manual(values = cluster_color_map(d$read_cluster), guide = "none") +
    scale_x_continuous(labels = function(x) format(round(x), big.mark = ",")) +
    facet_wrap(~ lab, ncol = 2) +
    labs(x = sprintf("Genomic position (%s)", d$chrom[1]),
         y = "Footprints per read per kb",
         title = sprintf("%s: nucleosome footprint intensity by topic-model read cluster",
                         which_gene),
         subtitle = paste0(sub,
                           "\nCurve height = occupancy, curve position = positioning. ",
                           "Grey = the other clusters.")) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(size = 8.5),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5))
}

# ---------------------------------------------------------------------------
# 2. Which genes carry which signal
# ---------------------------------------------------------------------------
plot_cluster_signal_quadrant <- function(stats_dt, alpha = 0.05) {
  d <- copy(stats_dt)
  d[, `:=`(x = -log10(q_position), y = -log10(q_occupancy))]
  ref <- -log10(alpha)
  xmax <- max(d$x, ref) * 1.18; ymax <- max(d$y, ref) * 1.22
  ins  <- 0.02   # inset so the corner labels do not sit on the panel border

  ggplot(d, aes(x, y)) +
    geom_hline(yintercept = ref, linetype = "dashed", color = "grey55", linewidth = 0.4) +
    geom_vline(xintercept = ref, linetype = "dashed", color = "grey55", linewidth = 0.4) +
    annotate("text", x = xmax * (1 - ins), y = ymax * (1 - ins),
             hjust = 1, vjust = 1, size = 3, color = "grey45", label = "both") +
    annotate("text", x = xmax * ins, y = ymax * (1 - ins),
             hjust = 0, vjust = 1, size = 3, color = "grey45", label = "occupancy differs") +
    annotate("text", x = xmax * (1 - ins), y = ymax * ins,
             hjust = 1, vjust = 0, size = 3, color = "grey45", label = "position differs") +
    geom_point(size = 2.6, color = "#2166ac") +
    ggrepel::geom_text_repel(aes(label = gene), size = 3.2, color = "grey15",
                             min.segment.length = 0.2, seed = 1) +
    scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0.02, 0))) +
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0.02, 0))) +
    labs(x = expression(Position~signal~~-log[10]~(q)~","~~permuted~Wasserstein),
         y = expression(Occupancy~signal~~-log[10]~(q)~","~~Fisher~exact),
         title = "Read clusters differ in position and in occupancy, but at different genes",
         subtitle = sprintf("Dashed lines mark q = %.2f. BH across the %d promoters, each axis separately.",
                            alpha, nrow(d))) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

main_cluster_distribution_plots <- function() {
  dir.create(OCC_DIR, recursive = TRUE, showWarnings = FALSE)
  prof  <- fread(file.path(RES_DIR, "cluster_intensity_profiles.tsv.gz"))
  stats <- fread(file.path(RES_DIR, "cluster_distribution_tests.tsv"))
  if (!nrow(prof)) { message("No intensity profiles to plot"); return(invisible(NULL)) }

  for (g in sort(unique(prof$gene))) {
    p <- plot_cluster_intensity_profiles(prof, g, stats)
    if (is.null(p)) next
    n_cl <- uniqueN(prof[gene == g, read_cluster])
    ggsave(file.path(OCC_DIR, sprintf("%s_intensity_by_cluster.pdf", g)), p,
           width = 10, height = 1.9 * ceiling(n_cl / 2) + 2)
  }
  ggsave(file.path(OCC_DIR, "gene_signal_quadrant.pdf"),
         plot_cluster_signal_quadrant(stats), width = 7.5, height = 6)
  message("Wrote cluster distribution plots to ", OCC_DIR)
}

if (sys.nframe() == 0L) {
  main_plots()
  main_nuc_topic_heatmaps()
  main_cluster_distribution_plots()
}
