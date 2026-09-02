# Shared definitions for the SAM-seq-style accessibility heterogeneity
# analysis (Leduque et al. 2024, NAR: 4-bin Shannon entropy of per-read m6A
# fractions, reads pooled across the LPS timepoints).
# Sourced by 01_entropy_windows.R, 03_entropy_violin.R and 05_read_plots.Rmd.

ENTROPY_DIR  <- "/project/spott/cshan/fiber-seq/macrophage_project/clustering_TSS_methods/entropy"
TAB_DIR      <- file.path(ENTROPY_DIR, "tables")
PLOT_DIR     <- file.path(ENTROPY_DIR, "plots")

SAMPLES <- c("LPS_0", "LPS_5", "LPS_10", "LPS_15")

MIN_READS   <- 20                      # pooled spanning reads needed for entropy
FRAC_BREAKS <- c(0, 0.25, 0.5, 0.75, 1)
MAX_ENTROPY <- log(4)
PLOT_HALF   <- 1000L                   # region plots span center +/- PLOT_HALF

#' 4-bin Shannon entropy of per-read m6A fractions (natural log, max ln 4).
#' @param fracs per-read m6A fractions of one window, pooled over timepoints.
#' @return one-row data.frame: n_reads, n_bin1..4, mean_frac, entropy.
shannon_entropy_4bin <- function(fracs) {
  stopifnot(all(fracs >= 0 & fracs <= 1))
  n <- tabulate(cut(fracs, breaks = FRAC_BREAKS, include.lowest = TRUE,
                    labels = FALSE), nbins = 4L)
  p <- n / sum(n)
  data.frame(n_reads = length(fracs),
             n_bin1 = n[1], n_bin2 = n[2], n_bin3 = n[3], n_bin4 = n[4],
             mean_frac = mean(fracs),
             entropy = -sum(ifelse(p > 0, p * log(p), 0)))
}
