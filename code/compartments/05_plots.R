source("/project/spott/cshan/fiber-seq/code/compartments/compartments_functions.r")
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })
PLOTS <- file.path(OUT, "plots"); dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)
read_saved <- function(name) fread(file.path(OUT, "tables", name))
domains <- read_saved("rao_domains.tsv")
s <- read_saved("boundary_summary.tsv")
s[, eligible := all(n_pairs >= 1000), by = .(signal, lag)]
s[eligible == FALSE, c("mean", "q25", "q75", "n_pairs", "within_minus_across") := NA_real_]
s[, distance_kb := distance_bp / 1000]
last_distance <- if (any(s$eligible)) max(s[eligible == TRUE, distance_bp]) / 1000 else 100
if (!is.finite(last_distance)) last_distance <- 100
distance_scale <- function() scale_x_continuous(limits = c(10, last_distance + 10))
shade <- function() annotate("rect", xmin = 10, xmax = 50, ymin = -Inf, ymax = Inf, fill = "grey70", alpha = .3)
colours <- c(within = "#2166ac", across = "#b2182b")
p <- ggplot(s, aes(distance_kb, mean, colour = class, fill = class)) + shade() +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = .18, colour = NA, na.rm = TRUE) + geom_line(na.rm = TRUE) +
  facet_wrap(~signal, nrow = 1) + scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
  labs(x = NULL, y = "Across-sample Pearson correlation") + theme_fiberseq() + distance_scale()
c <- ggplot(s, aes(distance_kb, n_pairs, fill = class)) + shade() +
  geom_col(position = "dodge", width = 16, na.rm = TRUE) + facet_wrap(~signal, nrow = 1) +
  scale_fill_manual(values = colours) + labs(x = "Distance (kb)", y = "Pairs") + theme_fiberseq() + distance_scale()
save_figure(plot_grid(p, c, ncol = 1, rel_heights = c(3, 1)), file.path(PLOTS, "within_vs_across_boundary.pdf"), 12, 6)
p <- ggplot(unique(s[, .(signal, distance_kb, within_minus_across)]), aes(distance_kb, within_minus_across)) +
  shade() + geom_hline(yintercept = 0, colour = "grey50") + geom_line(na.rm = TRUE) + facet_wrap(~signal, nrow = 1) +
  labs(x = "Distance (kb)", y = "Within minus across mean correlation") + theme_fiberseq() + distance_scale()
save_figure(p, file.path(PLOTS, "within_minus_across_boundary.pdf"), 12, 4)
pooled <- list(m6a = read_saved("bin_m6a_rate_pooled.tsv.gz"), nucleosome = read_saved("bin_nucleosome_fraction_pooled.tsv.gz"))
pooled$m6a[, value := rate]; pooled$nucleosome[, value := fraction]
domain_plot <- function(ch, limits) {
  d <- copy(domains[chrom == ch & start < limits[2] * 1e6 & end > limits[1] * 1e6])
  # Greedy lane packing keeps nested/overlapping domains visible.
  setorder(d, start, end); lane_end <- numeric(); d[, lane := 0L]
  for (i in seq_len(nrow(d))) {
    lane <- which(lane_end <= d$start[i])[1]
    if (is.na(lane)) lane <- length(lane_end) + 1L
    lane_end[lane] <- d$end[i]; d$lane[i] <- lane
  }
  ggplot(d) + geom_segment(aes(start / 1e6, lane, xend = end / 1e6, yend = lane), linewidth = 1.5) +
    coord_cartesian(xlim = limits) + labs(x = NULL, y = "Domains") + theme_fiberseq() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
signal_plot <- function(ch, signal, limits) ggplot(pooled[[signal]][chr == ch], aes(midpoint / 1e6, value)) +
  geom_line(linewidth = .25, na.rm = TRUE) + coord_cartesian(xlim = limits) +
  labs(x = "Position (Mb)", y = signal) + theme_fiberseq()
grDevices::pdf(file.path(PLOTS, "bin_signals_by_chromosome.pdf"), width = 12, height = 7)
for (ch in CHROMS) {
  limits <- c(0, max(pooled$m6a[chr == ch, end]) / 1e6)
  print(plot_grid(signal_plot(ch, "m6a", limits) + ggtitle(ch), signal_plot(ch, "nucleosome", limits),
    domain_plot(ch, limits), ncol = 1, rel_heights = c(2, 2, 1)))
}
dev.off()
# Base rasterImage embeds actual bitmap images in PDFs, keeping hundreds of pages small.
palette <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(256)
to_colours <- function(z, limits) {
  idx <- 1L + floor(255 * pmax(0, pmin(1, (z - limits[1]) / diff(limits))))
  ans <- matrix(palette[idx], nrow(z), ncol(z)); ans[is.na(z)] <- "grey85"; ans
}
heatmap_file <- file.path(PLOTS, "bin_signals_per_sample_heatmap.pdf")
grDevices::pdf(heatmap_file, width = 14, height = 8)
for (signal in names(pooled)) {
  x <- read_saved(paste0("normalized_", signal, ".tsv.gz"))
  sample_cols <- setdiff(names(x), c("bin", "chr", "start", "end", "midpoint"))
  lim <- range(as.matrix(x[, ..sample_cols]), finite = TRUE)
  if (diff(lim) == 0) lim <- lim + c(-.5, .5)
  for (ch in CHROMS) {
    b <- x[chr == ch]; setorder(b, start); limits <- range(c(b$start, b$end)) / 1e6
    layout(matrix(c(1, 2), 2, 1), heights = c(1, 5))
    par(mar = c(0, 14, 3, 3)); plot(NA, xlim = limits, ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "", main = paste(ch, signal))
    d <- domains[chrom == ch]; if (nrow(d)) segments(d$start / 1e6, .5, d$end / 1e6, .5, lwd = 3)
    par(mar = c(4, 14, 0, 3)); n <- length(sample_cols)
    plot(NA, xlim = limits, ylim = c(.5, n + .5), xaxs = "i", yaxs = "i", axes = FALSE, xlab = "Position (Mb)", ylab = "")
    rasterImage(to_colours(t(as.matrix(b[, ..sample_cols])), lim), limits[1], .5, limits[2], n + .5, interpolate = FALSE)
    axis(1); axis(2, at = n:1, labels = sample_cols, las = 2, cex.axis = .65)
    mtext(sprintf("QN log value: blue %.2f; white %.2f; red %.2f; missing grey", lim[1], mean(lim), lim[2]), side = 1, line = 2.6, cex = .7)
  }
}
dev.off()
grDevices::pdf(file.path(PLOTS, "covariation_maps.pdf"), width = 12, height = 7)
for (signal in names(pooled)) for (ch in CHROMS) {
  b <- pooled[[signal]][chr == ch]; setorder(b, start)
  pairs <- read_saved(file.path("pair_correlations", signal, paste0(ch, ".tsv.gz")))
  pairs[, left := b$midpoint[match(bin, b$bin)]]
  pairs[, center := left + lag * 10000]
  for (ws in seq(0, max(b$end) - 1, by = 10000000)) {
    we <- min(ws + 10000000, max(b$end)); limits <- c(ws, we) / 1e6
    # Rotated bin pairs lie on a 10 kb x 20 kb lattice. Duplicate the half-cell
    # offsets horizontally to avoid an artificial checkerboard in the raster.
    p <- pairs[center >= ws & center < we]
    z <- matrix(NA_real_, 100, ceiling((we - ws) / 10000))
    col <- floor((p$center - ws) / 10000) + 1L
    z[cbind(101L - p$lag, col)] <- p$r
    valid <- col > 1L
    z[cbind(101L - p$lag[valid], col[valid] - 1L)] <- p$r[valid]
    layout(matrix(1:3, 3, 1), heights = c(4, 1.2, .8))
    par(mar = c(0, 4, 3, 2))
    plot(NA, xlim = limits, ylim = c(0, 2), xaxs = "i", yaxs = "i", xlab = "", ylab = "Distance (Mb)", xaxt = "n",
      main = paste(signal, ch, sprintf("%.1f-%.1f Mb", limits[1], limits[2]), " | r: blue -1, white 0, red 1"))
    rasterImage(to_colours(z, c(-1, 1)), limits[1], 0, limits[2], 2, interpolate = FALSE)
    d <- domains[chrom == ch & start < we & end > ws]
    for (k in seq_len(nrow(d))) lines(c(d$start[k], (d$start[k] + d$end[k]) / 2, d$end[k]) / 1e6,
      c(0, (d$end[k] - d$start[k]) / 1e6, 0), col = "black", lwd = .6)
    rect(limits[1], 0, limits[2], .05, col = adjustcolor("grey60", alpha.f = .25), border = NA)
    par(mar = c(0, 4, 0, 2)); bb <- b[midpoint >= ws & midpoint < we]
    plot(bb$midpoint / 1e6, bb$value, type = "l", xlim = limits, ylim = c(0, 1), xaxs = "i", xlab = "", ylab = "Pooled", xaxt = "n")
    par(mar = c(4, 4, 0, 2)); plot(NA, xlim = limits, ylim = c(0, 1), xaxs = "i", xlab = "Position (Mb)", ylab = "Domains", yaxt = "n")
    if (nrow(d)) segments(d$start / 1e6, .5, d$end / 1e6, .5, lwd = 3)
  }
}
dev.off()
