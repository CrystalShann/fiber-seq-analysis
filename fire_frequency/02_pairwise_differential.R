#!/usr/bin/env Rscript
# Pairwise differential FIRE frequency across the LPS timecourse.
#
# For each region of the union peak universe and each of the 6 timepoint pairs,
# a two-sided Fisher exact test on the 2x2 of (FIRE, non-FIRE) read counts:
#
#                timepoint A     timepoint B
#   FIRE         n_fire_A        n_fire_B
#   non-FIRE     n_reads_A       n_reads_B 



suppressMessages({
  library(data.table)
  library(readr)
})


in_file  <- "/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/fire_frequency_long.tsv.gz"
out_file <- "/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/pairwise_fisher.tsv.gz"
timepoints <- c("LPS_0", "LPS_5", "LPS_10", "LPS_15")

fire_peak_universe <- "/project/spott/cshan/fiber-seq/macrophage_project/co-accessibility/universe/fire_peaks_union.bed"
fire_peak_universe_df <- read.csv(fire_peak_universe, sep = "\t", header = FALSE)
n_regions_expected <- nrow(fire_peak_universe_df)

# [1] "chrom"     "start"     "end"       "region_id" "timepoint" "n_reads"   "n_fire"    "freq"     
long <- fread(in_file)
stopifnot(nrow(long) == n_regions_expected * length(timepoints))

# one row per region, per timepoint n_reads / n_fire columns
# convert long to wide format

# [1] "chrom"          "start"          "end"            "region_id"      "n_reads_LPS_0"  "n_reads_LPS_10" "n_reads_LPS_15"
# [8] "n_reads_LPS_5"  "n_fire_LPS_0"   "n_fire_LPS_10"  "n_fire_LPS_15"  "n_fire_LPS_5"  

wide <- dcast(long, chrom + start + end + region_id ~ timepoint,
              value.var = c("n_reads", "n_fire"))
stopifnot(nrow(wide) == n_regions_expected)
setkey(wide, chrom, start, end)

# define lower and upper read depth limit for each time point
cov_lo <- cov_hi <- setNames(numeric(length(timepoints)), timepoints)

# calculate the median number of reads covering a FIRE universe region
    # lower threshold median read coverage of the timepoint minus a buffer 
    # for expected statistical fluctuation, with a hard floor of 10 reads
  
    # upper threshold: median read coverage plus the same buffer for expected fluctuation

for (s in timepoints) {
  m <- median(long[timepoint == s, n_reads])
  cov_lo[s] <- max(m - 3 * sqrt(m), 10)
  cov_hi[s] <- m + 3 * sqrt(m)
  cat(sprintf("%s: median n_reads %.1f, coverage window [%.1f, %.1f]\n",
              s, m, cov_lo[s], cov_hi[s]))
}

# generate all pairwise timepoint comparisons
pairs <- t(combn(timepoints, 2))   # A earlier, B later (timepoints is in time order)

# loop over each timepoint
  # A = earlier time point
  # B = later time point
res_list <- vector("list", nrow(pairs))

# extract counts needed for time point A and B
for (k in seq_len(nrow(pairs))) {
  A <- pairs[k, 1]; B <- pairs[k, 2]
  comparison <- paste0(A, "_vs_", B)

  d <- wide[, .(chrom, start, end, region_id,
                n_reads_A = get(paste0("n_reads_", A)),
                n_fire_A  = get(paste0("n_fire_",  A)),
                n_reads_B = get(paste0("n_reads_", B)),
                n_fire_B  = get(paste0("n_fire_",  B)))]
  d[, comparison := comparison]
  
  # Recalculate FIRE frequencies for A and B
    # frequency A = N fire, A / N reads, A
  d[, freq_A := ifelse(n_reads_A > 0, n_fire_A / n_reads_A, NA_real_)]
  d[, freq_B := ifelse(n_reads_B > 0, n_fire_B / n_reads_B, NA_real_)]
  
  # calculate change in FIRE frequency
  d[, delta_freq := freq_B - freq_A]
  
  ################################
  # Haldane +0.5 log2 odds ratio 
  ################################
  d[, log2_or := log2(((n_fire_B + 0.5) * (n_reads_A - n_fire_A + 0.5)) /
                      ((n_fire_A + 0.5) * (n_reads_B - n_fire_B + 0.5)))]
  # apply coverage filter
  d[, pass_coverage := n_reads_A >= cov_lo[A] & n_reads_A <= cov_hi[A] &
                       n_reads_B >= cov_lo[B] & n_reads_B <= cov_hi[B]]
  
  ################################
  # Fisher test
  ################################
  
  # keep unique count combinations
  uniq <- unique(d[, .(n_fire_A, n_reads_A, n_fire_B, n_reads_B)])
  # run fisher test
    # fa = FIRE reads at A
    # ra = total reads at A
    # fb = FIRE reads at B
    # rb = total reads at B
  uniq[, p_value := mapply(function(fa, ra, fb, rb) {
    fisher.test(matrix(c(fa, ra - fa, fb, rb - fb), nrow = 2))$p.value
  }, n_fire_A, n_reads_A, n_fire_B, n_reads_B)]
  uniq[, p_value := pmin(p_value, 1)]
  # join fisher results back to every region
  d <- uniq[d, on = .(n_fire_A, n_reads_A, n_fire_B, n_reads_B)]

  # rank by raw p among passing regions only 
  d[, rank_p := NA_integer_]
  d[pass_coverage == TRUE, rank_p := frank(p_value, ties.method = "min")]

  # sort regions by p value
  top1000 <- d[pass_coverage == TRUE][order(p_value)][1:1000]
  cat(sprintf("%s: %d distinct tables; %d/%d pass coverage; p<1e-3: %d; p<1e-5: %d; median |delta_freq| in top 1000: %.3f\n",
              comparison, nrow(uniq), sum(d$pass_coverage), nrow(d),
              d[pass_coverage == TRUE, sum(p_value < 1e-3)],
              d[pass_coverage == TRUE, sum(p_value < 1e-5)],
              top1000[, median(abs(delta_freq))]))
  cat(sprintf("  p histogram (passing, 10 bins 0-1): %s\n",
              paste(d[pass_coverage == TRUE, table(cut(p_value, seq(0, 1, 0.1), include.lowest = TRUE))],
                    collapse = " ")))

  res_list[[k]] <- d[, .(chrom, start, end, region_id, comparison,
                         n_reads_A, n_fire_A, freq_A, n_reads_B, n_fire_B, freq_B,
                         delta_freq, log2_or, p_value, rank_p, pass_coverage)]
}

res <- rbindlist(res_list)
stopifnot(nrow(res) == n_regions_expected * nrow(pairs))
stopifnot(res[, all(p_value >= 0 & p_value <= 1)])
stopifnot(res[, all(is.na(rank_p) == !pass_coverage)])
stopifnot(res[!is.na(delta_freq), all(abs(delta_freq) <= 1)])

fwrite(res, out_file, sep = "\t")
cat(sprintf("wrote %s (%d rows)\n", out_file, nrow(res)))
