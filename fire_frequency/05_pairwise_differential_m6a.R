#!/usr/bin/env Rscript
# Pairwise differential m6A level across the LPS timecourse.
#
# For each region of the union peak universe and each of the 6 timepoint pairs,
# a two-sided Wilcoxon rank-sum test on the per-read m6A fractions
# (read's called m6A / read's covered A/T sites, within the region) at
# timepoint A vs timepoint B. Reads are the independent unit, so unlike a
# Fisher test on pooled read x site events these p-values are not inflated by
# within-read correlation. Normal approximation (exact = FALSE): per-read
# fractions are heavily tied and group sizes are ~25-80.
#
# Effect size (delta_m6a) stays the pooled region-level difference
# (n_mod / n_cov, from 04), matching the quantity shown in the bin profiles.
# Coverage filter uses n_reads with the same windows logic as
# 02_pairwise_differential.R.

suppressMessages({
  library(data.table)
  library(readr)
})


region_file <- "/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/methylation/m6a_region_long.tsv.gz"
reads_file  <- "/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/methylation/m6a_reads_long.tsv.gz"
out_file    <- "/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/methylation/pairwise_wilcoxon_m6a.tsv.gz"
timepoints <- c("LPS_0", "LPS_5", "LPS_10", "LPS_15")

fire_peak_universe <- "/project/spott/cshan/fiber-seq/macrophage_project/co-accessibility/universe/fire_peaks_union.bed"
fire_peak_universe_df <- read.csv(fire_peak_universe, sep = "\t", header = FALSE)
n_regions_expected <- nrow(fire_peak_universe_df)

# "chrom" "start" "end" "region_id" "timepoint" "n_sites" "n_cov" "n_mod" "n_reads" "m6a"
long <- fread(region_file)
stopifnot(nrow(long) == n_regions_expected * length(timepoints))

# "chrom" "region_id" "timepoint" "n_sites_cov" "n_mod" "frac"
reads <- fread(reads_file)
stopifnot(nrow(reads) == long[, sum(n_reads)])
# a read's fraction is defined only when it covers >= 1 A/T site in the region
reads <- reads[!is.na(frac)]

# one row per region, per timepoint n_cov / n_mod / n_reads columns
wide <- dcast(long, chrom + start + end + region_id ~ timepoint,
              value.var = c("n_cov", "n_mod", "n_reads"))
stopifnot(nrow(wide) == n_regions_expected)
setkey(wide, chrom, start, end)

# read-depth window per timepoint, same formula as 02_pairwise_differential.R
cov_lo <- cov_hi <- setNames(numeric(length(timepoints)), timepoints)
for (s in timepoints) {
  m <- median(long[timepoint == s, n_reads])
  cov_lo[s] <- max(m - 3 * sqrt(m), 10)
  cov_hi[s] <- m + 3 * sqrt(m)
  cat(sprintf("%s: median n_reads %.1f, coverage window [%.1f, %.1f]\n",
              s, m, cov_lo[s], cov_hi[s]))
}

pairs <- t(combn(timepoints, 2))   # A earlier, B later (timepoints is in time order)

res_list <- vector("list", nrow(pairs))
for (k in seq_len(nrow(pairs))) {
  A <- pairs[k, 1]; B <- pairs[k, 2]
  comparison <- paste0(A, "_vs_", B)

  d <- wide[, .(chrom, start, end, region_id,
                n_cov_A   = get(paste0("n_cov_",   A)),
                n_mod_A   = get(paste0("n_mod_",   A)),
                n_reads_A = get(paste0("n_reads_", A)),
                n_cov_B   = get(paste0("n_cov_",   B)),
                n_mod_B   = get(paste0("n_mod_",   B)),
                n_reads_B = get(paste0("n_reads_", B)))]
  d[, comparison := comparison]

  d[, m6a_A := ifelse(n_cov_A > 0, n_mod_A / n_cov_A, NA_real_)]
  d[, m6a_B := ifelse(n_cov_B > 0, n_mod_B / n_cov_B, NA_real_)]

  # positive = higher m6A (more accessible) at the later timepoint
  d[, delta_m6a := m6a_B - m6a_A]

  # apply coverage filter (read depth at both timepoints of the pair)
  d[, pass_coverage := n_reads_A >= cov_lo[A] & n_reads_A <= cov_hi[A] &
                       n_reads_B >= cov_lo[B] & n_reads_B <= cov_hi[B]]

  ################################
  # Wilcoxon rank-sum test on per-read fractions
  ################################
  pw <- reads[timepoint %in% c(A, B),
    {
      fa <- frac[timepoint == A]
      fb <- frac[timepoint == B]
      p <- if (length(fa) == 0 || length(fb) == 0) {
        NA_real_                              # one side unobserved: untestable
      } else if (uniqueN(c(fa, fb)) <= 1) {
        1                                     # all values tied: no evidence
      } else {
        wilcox.test(fa, fb, exact = FALSE)$p.value
      }
      .(n_frac_A = length(fa), n_frac_B = length(fb), p_value = p)
    },
    by = region_id]
  d <- pw[d, on = "region_id"]
  d[is.na(n_frac_A), `:=`(n_frac_A = 0L, n_frac_B = 0L)]
  d[, p_value := pmin(p_value, 1)]

  # rank by raw p among passing, testable regions only
  d[, rank_p := NA_integer_]
  d[pass_coverage == TRUE & !is.na(p_value),
    rank_p := frank(p_value, ties.method = "min")]

  top1000 <- d[pass_coverage == TRUE & !is.na(p_value)][order(p_value)][1:1000]
  cat(sprintf("%s: %d/%d pass coverage (%d untestable); p<0.05: %d; p<1e-3: %d; median |delta_m6a| in top 1000: %.4f\n",
              comparison, sum(d$pass_coverage), nrow(d),
              d[pass_coverage == TRUE, sum(is.na(p_value))],
              d[pass_coverage == TRUE, sum(p_value < 0.05, na.rm = TRUE)],
              d[pass_coverage == TRUE, sum(p_value < 1e-3, na.rm = TRUE)],
              top1000[, median(abs(delta_m6a))]))
  cat(sprintf("  p histogram (passing, 10 bins 0-1): %s\n",
              paste(d[pass_coverage == TRUE & !is.na(p_value),
                      table(cut(p_value, seq(0, 1, 0.1), include.lowest = TRUE))],
                    collapse = " ")))

  res_list[[k]] <- d[, .(chrom, start, end, region_id, comparison,
                         n_cov_A, n_mod_A, m6a_A, n_cov_B, n_mod_B, m6a_B,
                         n_reads_A, n_reads_B, n_frac_A, n_frac_B,
                         delta_m6a, p_value, rank_p, pass_coverage)]
}

res <- rbindlist(res_list)
stopifnot(nrow(res) == n_regions_expected * nrow(pairs))
stopifnot(res[!is.na(p_value), all(p_value >= 0 & p_value <= 1)])
stopifnot(res[, all(is.na(rank_p) == (!pass_coverage | is.na(p_value)))])
stopifnot(res[!is.na(delta_m6a), all(abs(delta_m6a) <= 1)])

fwrite(res, out_file, sep = "\t")
cat(sprintf("wrote %s (%d rows)\n", out_file, nrow(res)))
