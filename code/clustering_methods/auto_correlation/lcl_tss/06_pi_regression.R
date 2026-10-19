#' Does pausing index predict TSS nucleosome phasing at fixed expression?
#'
#' One row per gene from 05_tss_acf_all_genes.R; the gene is the unit, not the molecule.
#' For each outcome (ACF peak height and peak lag, 140-250 bp) an ordinary least-squares fit
#'     outcome ~ log10(pi) + log10(mean_tpm)
#' over the genes with at least MIN_READS ACF-valid reads and a peak. Expression is the
#' only covariate: it takes the place of the pilot's expression matching (stage 01), since
#' log10 PI and log10 TPM are correlated.
#'
#' Input:  <OUT_DIR>/tables/gene_acf_summary.tsv   (stage 05)
#' Output: <OUT_DIR>/tables/pi_regression.tsv      per outcome and term: estimate, SE,
#'                                                 95% CI, t, p, n genes, R^2

suppressPackageStartupMessages(library(data.table))

OUT_DIR  <- "/project/spott/cshan/fiber-seq/LCL_project/auto_correlation/lcl_tss/continuous_pi"
GENES    <- file.path(OUT_DIR, "tables", "gene_acf_summary.tsv")
OUTCOMES <- c("peak_height", "peak_lag_bp")
MIN_READS <- 20L   # ACF-valid reads behind a gene's mean ACF

stopifnot(file.exists(GENES))
genes <- fread(GENES)
stopifnot(nrow(genes) > 1000L, genes[, all(pi > 0 & mean_tpm > 0)])
n_all <- nrow(genes)
genes <- genes[n_acf_valid >= MIN_READS]
message(n_all - nrow(genes), " of ", n_all, " genes dropped with < ", MIN_READS,
        " ACF-valid reads")
genes[, `:=`(log10_pi = log10(pi), log10_tpm = log10(mean_tpm))]

fits <- rbindlist(lapply(OUTCOMES, function(y) {
  d   <- genes[is.finite(get(y))]   # genes with no 140-250 bp peak drop out
  fit <- lm(reformulate(c("log10_pi", "log10_tpm"), y), data = d)
  co  <- coef(summary(fit)); ci <- confint(fit)
  data.table(outcome = y, term = rownames(co), estimate = co[, 1], std_error = co[, 2],
             conf_low = ci[, 1], conf_high = ci[, 2], t_value = co[, 3], p_value = co[, 4],
             n_genes = nrow(d), r_squared = summary(fit)$r.squared)
}))
fwrite(fits, file.path(OUT_DIR, "tables", "pi_regression.tsv"), sep = "\t")
print(fits)
message("Wrote ", file.path(OUT_DIR, "tables", "pi_regression.tsv"))
