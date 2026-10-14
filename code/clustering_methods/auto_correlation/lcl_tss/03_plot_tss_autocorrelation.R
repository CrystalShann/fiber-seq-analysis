#' Figures for the high- vs low-pausing TSS autocorrelation comparison.
#'
#' PDF only. Reads the tables and per-region acf.rds written by
#' 02_tss_autocorrelation.R and does no computation beyond summarising for display.
#'
#' Figures, under <OUT_DIR>/results/plots/ :
#'   01_group_mean_acf.pdf        the headline: group mean ACF, high- vs low-PI
#'   02_per_pair_mean_acf.pdf     the same per matched pair, to show consistency
#'   03_acf_heatmap.pdf           per-molecule ACF, one row per molecule, by group
#'   04_m6a_call_fraction.pdf     the density confound, per gene and per group
#'   05_nrl_per_molecule.pdf      per-molecule repeat_peak lag (UNCALIBRATED)
#'   06_paired_comparison.pdf     the PRIMARY statistic: paired per-gene peak height
#'                                and lag, one line per matched pair
#'
#' Curves start at lag 25 bp: the lag-0 shoulder dominates any shared axis and carries
#' no nucleosome information.

suppressPackageStartupMessages({ library(data.table) })

OUT_DIR  <- "/project/spott/cshan/fiber-seq/LCL_project/auto_correlation/lcl_tss"
TAB      <- file.path(OUT_DIR, "tables")
PLOTS    <- file.path(OUT_DIR, "results", "plots")

LAG_FROM <- 25L     # first lag shown
LAG_TO   <- 1000L   # last lag shown; beyond this the single-molecule ACF is noise
NRL_BAND <- c(140L, 250L)
BAND_FILL <- rgb(0, 0, 0, 0.05)   # shading for the repeat_peak search window
COL      <- c(high_pi = "#D55E00", low_pi = "#2878B5")
LAB      <- c(high_pi = "high pausing index", low_pi = "low pausing index")

dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)
grp   <- fread(file.path(TAB, "acf_group_summary.tsv"))
gene  <- fread(file.path(TAB, "acf_per_gene.tsv"))
qc    <- fread(file.path(TAB, "read_qc.tsv"))
peaks <- fread(file.path(TAB, "group_peaks.tsv"))
mol   <- fread(file.path(TAB, "nrl_per_molecule.tsv"))
regions <- fread(file.path(TAB, "tss_regions.tsv"))
stopifnot(nrow(grp) > 0, nrow(gene) > 0, uniqueN(grp$pi_group) == 2L)

# ---- 01 group mean ACF ------------------------------------------------------
pdf(file.path(PLOTS, "01_group_mean_acf.pdf"), width = 9, height = 6, useDingbats = FALSE)
d <- grp[lag_bp >= LAG_FROM & lag_bp <= LAG_TO]
plot(NA, xlim = range(d$lag_bp), ylim = range(c(d$mean_acf - d$sem_acf,
                                                d$mean_acf + d$sem_acf)),
     xlab = "Lag (bp)", ylab = "Mean m6A autocorrelation",
     main = "Single-molecule m6A autocorrelation at the TSS", las = 1)
abline(h = 0, col = "grey70")
rect(NRL_BAND[1], par("usr")[3], NRL_BAND[2], par("usr")[4],
     col = BAND_FILL, border = NA)
for (g in names(COL)) {
  dd <- d[pi_group == g]
  polygon(c(dd$lag_bp, rev(dd$lag_bp)),
          c(dd$mean_acf - dd$sem_acf, rev(dd$mean_acf + dd$sem_acf)),
          col = adjustcolor(COL[[g]], 0.25), border = NA)
  lines(dd$lag_bp, dd$mean_acf, col = COL[[g]], lwd = 2)
  pk <- peaks[level == "pi_group" & label == g]
  if (nrow(pk) && !is.na(pk$peak_lag_bp))
    points(pk$peak_lag_bp, dd[lag_bp == pk$peak_lag_bp]$mean_acf,
           pch = 8, col = COL[[g]], cex = 1.4)
}
legend("topright", bty = "n", lwd = 2, col = unname(COL[names(COL)]),
       legend = vapply(names(COL), function(g) {
         pk <- peaks[level == "pi_group" & label == g]
         n  <- grp[pi_group == g]$n_reads[1]
         sprintf("%s  (n=%d, peak %s bp)", LAB[[g]], n,
                 if (nrow(pk)) pk$peak_lag_bp else NA)
       }, character(1)))
mtext(sprintf("shaded band = repeat_peak search window %d-%d bp; ribbons = SEM; lags from %d bp",
              NRL_BAND[1], NRL_BAND[2], LAG_FROM), side = 3, line = 0.2, cex = 0.75)
dev.off()

# ---- 02 per matched pair ----------------------------------------------------
pdf(file.path(PLOTS, "02_per_pair_mean_acf.pdf"), width = 10, height = 7, useDingbats = FALSE)
op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
for (p in sort(unique(gene$pair))) {
  dd <- gene[pair == p & lag_bp >= LAG_FROM & lag_bp <= LAG_TO]
  plot(NA, xlim = range(dd$lag_bp), ylim = range(dd$mean_acf),
       xlab = "Lag (bp)", ylab = "Mean ACF", las = 1,
       main = paste0("pair ", p))
  abline(h = 0, col = "grey70")
  for (g in names(COL)) {
    x <- dd[pi_group == g]
    if (!nrow(x)) next
    lines(x$lag_bp, x$mean_acf, col = COL[[g]], lwd = 1.8)
  }
  lg <- regions[pair == p][order(pi_group)]
  legend("topright", bty = "n", cex = 0.75, lwd = 1.8,
         col = COL[lg$pi_group],
         legend = sprintf("%s  PI=%.2f  %.1f TPM", lg$gene_name, lg$pi, lg$mean_tpm))
}
par(op)
mtext("Matched pairs: same expression, ~1000-fold different pausing index",
      outer = TRUE, cex = 0.95)
dev.off()

# ---- 03 per-molecule ACF heatmap -------------------------------------------
pdf(file.path(PLOTS, "03_acf_heatmap.pdf"), width = 10, height = 7, useDingbats = FALSE)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
lag_keep <- LAG_FROM:LAG_TO
for (g in names(COL)) {
  rs <- regions[pi_group == g]$region_id
  m  <- do.call(rbind, lapply(rs, function(r) {
    x  <- readRDS(file.path(OUT_DIR, "outputs", r, "acf.rds"))
    ok <- x$records$acf_valid
    a  <- x$acf[ok, lag_keep + 1L, drop = FALSE]
    a[order(-x$records$m6a_count[ok]), , drop = FALSE]   # densest molecules on top
  }))
  lim <- quantile(abs(m), 0.99, na.rm = TRUE)
  image(x = lag_keep, y = seq_len(nrow(m)), z = t(pmax(pmin(m, lim), -lim)),
        col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(255),
        zlim = c(-lim, lim), useRaster = TRUE,
        xlab = "Lag (bp)", ylab = "Molecule", main = LAB[[g]], las = 1)
  abline(v = NRL_BAND, col = "grey30", lty = 2)
  mtext(sprintf("n=%d; colour saturates at the 99th pct of |ACF| (%.3f)", nrow(m), lim),
        side = 3, line = 0.2, cex = 0.7)
}
par(op)
dev.off()

# ---- 04 m6A call fraction: the density confound ----------------------------
pdf(file.path(PLOTS, "04_m6a_call_fraction.pdf"), width = 9, height = 6, useDingbats = FALSE)
op <- par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))
b <- split(mol$m6a_call_fraction, factor(mol$pi_group, levels = names(COL)))
boxplot(b, col = adjustcolor(unname(COL[names(COL)]), 0.4), las = 2, outline = FALSE,
        names = LAB[names(COL)], ylab = "m6A calls / window width",
        main = "m6A density by group")
pv <- wilcox.test(m6a_call_fraction ~ pi_group, data = mol)$p.value
mtext(sprintf("Wilcoxon p = %.2g", pv), side = 3, line = 0.2, cex = 0.75)
qq <- qc[order(pi_group, -mean_m6a_call_fraction)]
barplot(qq$mean_m6a_call_fraction, names.arg = qq$gene_name, las = 2,
        col = COL[qq$pi_group], ylab = "mean m6A calls / window width",
        main = "m6A density by gene")
par(op)
dev.off()

# ---- 05 per-molecule NRL (uncalibrated) ------------------------------------
pdf(file.path(PLOTS, "05_nrl_per_molecule.pdf"), width = 9, height = 6, useDingbats = FALSE)
op <- par(mar = c(5, 4, 4, 1))
h <- split(mol[!is.na(nrl_bp)]$nrl_bp, mol[!is.na(nrl_bp)]$pi_group)
plot(NA, xlim = NRL_BAND, ylim = c(0, 0.035), xlab = "Per-molecule repeat_peak lag (bp)",
     ylab = "Density", main = "Per-molecule nucleosome repeat length", las = 1)
for (g in names(COL)) {
  if (is.null(h[[g]])) next
  d <- density(h[[g]], from = NRL_BAND[1], to = NRL_BAND[2])
  lines(d$x, d$y, col = COL[[g]], lwd = 2)
  abline(v = median(h[[g]]), col = COL[[g]], lty = 2)
}
legend("topright", bty = "n", lwd = 2, col = unname(COL[names(COL)]),
       legend = sprintf("%s  median %.0f bp", LAB[names(COL)],
                        vapply(names(COL), function(g) median(h[[g]]), numeric(1))))
mtext(paste("UNCALIBRATED: repeat_peak() applies no prominence constraint, so single",
            "molecules can return a noise bump. Group curves are the primary readout."),
      side = 3, line = 0.2, cex = 0.7)
par(op)
dev.off()

# ---- 06 the primary statistic: paired, per gene ---------------------------
# One line per matched pair. This is the inferential figure: molecules within a gene
# are not independent, so the gene is the experimental unit and the pairs are the
# replicates. The pooled curves in figure 01 are descriptive only.
wide   <- fread(file.path(TAB, "pair_comparison.tsv"))
paired <- fread(file.path(TAB, "paired_tests.tsv"))
pdf(file.path(PLOTS, "06_paired_comparison.pdf"), width = 9, height = 5, useDingbats = FALSE)
op <- par(mfrow = c(1, 2), mar = c(5, 4.5, 4, 1))
panel <- function(hi, lo, ylab, main, stat) {
  keep <- is.finite(hi) & is.finite(lo)
  plot(NA, xlim = c(0.7, 2.3), ylim = range(c(hi, lo), na.rm = TRUE), xaxt = "n",
       xlab = "", ylab = ylab, main = main, las = 1)
  axis(1, at = 1:2, labels = c("high PI", "low PI"))
  segments(1, hi[keep], 2, lo[keep], col = "grey55")
  points(rep(1, sum(keep)), hi[keep], pch = 19, col = COL[["high_pi"]], cex = 1.3)
  points(rep(2, sum(keep)), lo[keep], pch = 19, col = COL[["low_pi"]], cex = 1.3)
  text(rep(1, sum(keep)), hi[keep], wide$gene_name_high_pi[keep], pos = 2, cex = 0.6)
  text(rep(2, sum(keep)), lo[keep], wide$gene_name_low_pi[keep],  pos = 4, cex = 0.6)
  st <- paired[statistic == stat]
  mtext(sprintf("paired Wilcoxon p = %.3f (n = %d pairs); %d/%d favour high PI",
                st$wilcoxon_p, st$n_pairs, st$n_high_greater, st$n_pairs),
        side = 3, line = 0.2, cex = 0.72)
}
panel(wide$peak_height_high_pi, wide$peak_height_low_pi,
      "ACF peak height (140-250 bp)", "Phasing amplitude", "acf_peak_height")
panel(wide$peak_lag_bp_high_pi, wide$peak_lag_bp_low_pi,
      "ACF peak lag (bp)", "Nucleosome repeat length", "acf_peak_lag_bp")
par(op)
mtext("Primary statistic: gene is the unit, pairs are matched on expression",
      outer = TRUE, line = -1.3, cex = 0.9)
dev.off()

message("Wrote 6 PDFs to ", PLOTS)
print(list.files(PLOTS))
