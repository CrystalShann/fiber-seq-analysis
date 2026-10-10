#!/usr/bin/env Rscript
# Step 4: summary of per-read FIRE frequency by region set (and by cCRE class within active_cre),
# ECDF by set, and between-sample correlation.
# Usage: Rscript summarize.R <run_dir>   reads <run_dir>/fire_freq.<set>.tsv, writes <run_dir>/summary/
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

run_dir <- commandArgs(trailingOnly = TRUE)[1]
out_dir <- file.path(run_dir, "summary")
dir.create(out_dir, showWarnings = FALSE)
sets <- c("active_cre", "low_dnase_cre", "no_dnase_cre_region")
min_reads <- 10

tabs <- setNames(lapply(sets, function(s) fread(file.path(run_dir, sprintf("fire_freq.%s.tsv", s)))), sets)
samples <- sub("^n_reads_", "", grep("^n_reads_", names(tabs[[1]]), value = TRUE))
cat(sprintf("Samples (%d): %s\n\n", length(samples), paste(samples, collapse = ", ")))

for (s in sets) {
  x <- tabs[[s]]
  bad <- x[n_fire > n_reads, .N] +
    sum(sapply(samples, function(k) sum(x[[paste0("n_fire_", k)]] > x[[paste0("n_reads_", k)]])))
  cat(sprintf("%-20s n_fire > n_reads violations (pooled + per sample): %d\n", s, bad))
}

# median and fraction > 0.1 use regions with >= min_reads spanning reads; pooled uses every region
stats <- function(n_reads, n_fire, fire_freq) {
  ok <- n_reads >= min_reads
  list(n_regions = length(n_reads), n_reads_ge10 = sum(ok), median_n_reads = as.numeric(median(n_reads)),
       median_fire_freq = median(fire_freq[ok]), pooled_fire_freq = sum(n_fire) / sum(n_reads),
       frac_fire_freq_gt_0.1 = mean(fire_freq[ok] > 0.1))
}
by_set <- rbindlist(lapply(sets, function(s)
  tabs[[s]][, c(list(set = s, class = "all"), stats(n_reads, n_fire, fire_freq))]))
by_class <- tabs$active_cre[, c(list(set = "active_cre"), stats(n_reads, n_fire, fire_freq)), by = class]
summary_dt <- rbind(by_set, by_class[order(-n_regions)], use.names = TRUE)
fwrite(summary_dt, file.path(out_dir, "fire_freq_summary.tsv"), sep = "\t")
cat("\nFIRE frequency summary (median / frac > 0.1 over regions with n_reads >= 10):\n")
print(summary_dt[, lapply(.SD, function(v) if (is.double(v)) signif(v, 3) else v)], row.names = FALSE)

ecdf_dt <- rbindlist(lapply(sets, function(s) tabs[[s]][n_reads >= min_reads, .(set = s, fire_freq)]))
n_ok <- ecdf_dt[, .N, by = set][match(sets, set), N]
ecdf_dt[, set := factor(set, levels = sets,
                        labels = sprintf("%s (n = %s)", sets, format(n_ok, big.mark = ",", trim = TRUE)))]
p <- ggplot(ecdf_dt, aes(fire_freq, colour = set)) +
  stat_ecdf(linewidth = 0.7) +
  scale_colour_manual(values = c("#2a78d6", "#eb6834", "#1baf7a"), name = NULL) +
  labs(x = "Per-read FIRE frequency (n_fire / n_reads)", y = "Fraction of regions",
       title = "Per-read FIRE frequency by region set",
       subtitle = sprintf("Regions with >= %d spanning reads, pooled over %d sample(s)", min_reads, length(samples))) +
  theme_classic(base_size = 11) +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.04), legend.justification = c(1, 0))
ggsave(file.path(out_dir, "fire_freq_ecdf.pdf"), p, width = 6.5, height = 4.5)

if (length(samples) > 1) {
  sample_pairs <- CJ(a = samples, b = samples)[a < b]
  cors <- rbindlist(lapply(sets, function(s) {
    x <- tabs[[s]]
    sample_pairs[, {
      ra <- x[[paste0("n_reads_", a)]]; rb <- x[[paste0("n_reads_", b)]]
      ok <- ra >= min_reads & rb >= min_reads
      fa <- x[[paste0("n_fire_", a)]][ok] / ra[ok]
      fb <- x[[paste0("n_fire_", b)]][ok] / rb[ok]
      list(set = s, n_regions = sum(ok), pearson = suppressWarnings(cor(fa, fb)),
           spearman = suppressWarnings(cor(fa, fb, method = "spearman")))
    }, by = .(a, b)]
  }))
  fwrite(cors, file.path(out_dir, "fire_freq_sample_correlation.tsv"), sep = "\t")
  cat(sprintf("\nBetween-sample correlation of fire_freq (regions with n_reads >= %d in both):\n", min_reads))
  print(cors[, .(n_pairs = .N, median_n_regions = median(n_regions),
                 pearson_median = median(pearson, na.rm = TRUE), pearson_min = min(pearson, na.rm = TRUE),
                 pearson_max = max(pearson, na.rm = TRUE), spearman_median = median(spearman, na.rm = TRUE)),
             by = set][, lapply(.SD, function(v) if (is.double(v)) signif(v, 3) else v)], row.names = FALSE)
}
cat(sprintf("\nWrote %s/{fire_freq_summary.tsv, fire_freq_ecdf.pdf%s}\n", out_dir,
            if (length(samples) > 1) ", fire_freq_sample_correlation.tsv" else ""))
