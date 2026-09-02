#!/usr/bin/env Rscript
# 1000 bp windows for the accessibility entropy analysis: TSS +/- 500 for
# every gene in tss_expression_bins.tsv, plus the 10 FIRE-frequency regions --
# among coverage-passing pairwise Fisher results, per region the min-p
# comparison; top 5 by p with delta_freq > 0 (gain over the LPS timecourse)
# and top 5 with delta_freq < 0 (loss) -- at region midpoint +/- 500.

suppressPackageStartupMessages(library(data.table))

source("/project/spott/cshan/fiber-seq/code/clustering_TSS_methods/entropy/entropy_functions.r")

HALF <- 500L

BINS_TSV   <- "/project/spott/cshan/fiber-seq/macrophage_project/expr_access/tables/tss_expression_bins.tsv"
FISHER_TSV <- "/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/pairwise_fisher.tsv.gz"

dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

## TSS windows ---------------------------------------------------------------
genes <- fread(BINS_TSV)
genes[, `:=`(win_start = as.integer(tss - HALF),
             win_end   = as.integer(tss + HALF))]
stopifnot(genes[, all(win_end - win_start == 2L * HALF)])
setorder(genes, chrom, win_start)

out_tss <- file.path(TAB_DIR, "entropy_windows_tss.tsv")
fwrite(genes, out_tss, sep = "\t")
cat("wrote", out_tss, ":", nrow(genes), "TSS windows,",
    uniqueN(genes[, .(chrom, win_start)]), "distinct loci\n")

## FIRE windows --------------------------------------------------------------
fisher <- fread(FISHER_TSV)
fisher <- fisher[pass_coverage == TRUE]
stopifnot(nrow(fisher) > 0, !anyNA(fisher$p_value))

COMPARISONS <- c("LPS_0_vs_LPS_5", "LPS_0_vs_LPS_10", "LPS_0_vs_LPS_15",
                 "LPS_5_vs_LPS_10", "LPS_5_vs_LPS_15", "LPS_10_vs_LPS_15")
stopifnot(all(fisher$comparison %in% COMPARISONS))

# per region the min-p row; deterministic tie-break: p, then |delta_freq|
# desc, then |log2_or| desc, then fixed comparison order
fisher[, `:=`(abs_delta  = abs(delta_freq),
              abs_or     = abs(log2_or),
              comp_order = match(comparison, COMPARISONS))]
setorder(fisher, p_value, -abs_delta, -abs_or, comp_order)
best <- fisher[, .SD[1L], by = region_id]

setorder(best, p_value)
top <- rbind(best[delta_freq > 0][1:5][, direction := "gain"],
             best[delta_freq < 0][1:5][, direction := "loss"])
stopifnot(nrow(top) == 10L, !anyNA(top$region_id))

top[, `:=`(region_width = end - start,
           midpoint     = (start + end) %/% 2L)]
top[, `:=`(win_start = as.integer(midpoint - HALF),
           win_end   = as.integer(midpoint + HALF))]
wide <- top[region_width > 2L * HALF]
if (nrow(wide) > 0)
  cat("note:", nrow(wide), "selected regions wider than", 2L * HALF, "bp;",
      "the window will not cover them fully:",
      paste(wide$region_id, collapse = ", "), "\n")

out <- top[, .(region_id, chrom, region_start = start, region_end = end,
               region_width, midpoint, win_start, win_end, direction,
               best_comparison = comparison, min_p = p_value,
               delta_freq_best = delta_freq, log2_or_best = log2_or,
               n_reads_A, n_fire_A, freq_A, n_reads_B, n_fire_B, freq_B)]
out_fire <- file.path(TAB_DIR, "entropy_windows_fire.tsv")
fwrite(out, out_fire, sep = "\t")
cat("wrote", out_fire, ":", nrow(out), "FIRE windows\n")
print(out[, .(region_id, direction, best_comparison, min_p, delta_freq_best,
              region_width)])
