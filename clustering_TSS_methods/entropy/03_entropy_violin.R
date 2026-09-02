#!/usr/bin/env Rscript
# Pooled-timepoint 4-bin Shannon entropy per 1000 bp window (Leduque et al.
# 2024 accessibility heterogeneity): per-gene and FIRE-region entropy tables,
# and selected_regions.tsv (top 5 highest-entropy genes per quartile + the 10
# FIRE regions). Tables only; the Leiden-based entropy analysis lives in
# 04_leiden_entropy.R.

suppressPackageStartupMessages(library(data.table))

source("/project/spott/cshan/fiber-seq/code/clustering_TSS_methods/entropy/entropy_functions.r")

QUARTILES <- c("Q1_low", "Q2", "Q3", "Q4_high")

COUNT_COLS <- c("n_reads", "n_bin1", "n_bin2", "n_bin3", "n_bin4")

# entropy per window, reads pooled across timepoints (timepoint is ignored)
entropy_by_window <- function(fr) {
  if (nrow(fr) == 0)
    return(data.table(window_id = character(), n_reads = integer(),
                      n_bin1 = integer(), n_bin2 = integer(),
                      n_bin3 = integer(), n_bin4 = integer(),
                      mean_frac = numeric(), entropy = numeric()))
  stopifnot(fr[, all(frac >= 0 & frac <= 1)])
  ent <- fr[, shannon_entropy_4bin(frac), by = window_id]
  stopifnot(ent[, all(n_bin1 + n_bin2 + n_bin3 + n_bin4 == n_reads)],
            ent[, all(entropy <= MAX_ENTROPY + 1e-9)])
  ent
}

## TSS entropy ---------------------------------------------------------------
win_tss <- fread(file.path(TAB_DIR, "entropy_windows_tss.tsv"))
fr_tss  <- fread(file.path(TAB_DIR, "read_fractions_tss.tsv.gz"))

tss <- merge(win_tss, entropy_by_window(fr_tss),
             by.x = "gene_id", by.y = "window_id", all.x = TRUE)
for (col in COUNT_COLS)
  set(tss, which(is.na(tss[[col]])), col, 0L)
tss[n_reads < MIN_READS, entropy := NA_real_]

cat("entropy defined for", tss[!is.na(entropy), .N], "of", nrow(tss),
    "genes (MIN_READS =", MIN_READS, ")\n")
print(tss[, .(n = .N, n_entropy = sum(!is.na(entropy)),
              median_entropy = median(entropy, na.rm = TRUE)), by = expr_bin])

setorder(tss, chrom, win_start)
out_tss <- file.path(TAB_DIR, "tss_entropy.tsv.gz")
fwrite(tss, out_tss, sep = "\t")
cat("wrote", out_tss, "\n")

## FIRE entropy --------------------------------------------------------------
win_fire <- fread(file.path(TAB_DIR, "entropy_windows_fire.tsv"))
fr_fire  <- fread(file.path(TAB_DIR, "read_fractions_fire.tsv.gz"))

fire <- merge(win_fire, entropy_by_window(fr_fire),
              by.x = "region_id", by.y = "window_id", all.x = TRUE)
for (col in COUNT_COLS)
  set(fire, which(is.na(fire[[col]])), col, 0L)
fire[n_reads < MIN_READS, entropy := NA_real_]

out_fire <- file.path(TAB_DIR, "fire_entropy.tsv")
fwrite(fire, out_fire, sep = "\t")
cat("wrote", out_fire, "\n")
print(fire[, .(region_id, direction, min_p, delta_freq_best, n_reads, entropy)])

## Top 5 highest-entropy genes per quartile ----------------------------------
sel <- tss[expr_bin %in% QUARTILES & !is.na(entropy)]
# genes sharing a TSS give identical windows; keep the higher-TPM gene per locus
setorder(sel, chrom, win_start, -mean_tpm)
sel <- unique(sel, by = c("chrom", "win_start"))
sel[, expr_bin := factor(expr_bin, levels = QUARTILES)]
setorder(sel, expr_bin, -entropy, -n_reads, gene_id)
top5 <- sel[, head(.SD, 5L), by = expr_bin]
if (nrow(top5) != 20L)
  warning("selected ", nrow(top5), " genes, not 20 (partial-genome test run?)")

top5[, outname := gene_name]
dup <- duplicated(top5$outname) | duplicated(top5$outname, fromLast = TRUE)
top5[dup, outname := paste0(outname, "_", sub("\\..*$", "", gene_id))]

sel_tss <- top5[, .(outname,
                    region_type = "tss",
                    plot_group = as.character(expr_bin),
                    window_id = gene_id,
                    label = paste0(gene_name, " (", expr_bin, ", TPM ",
                                   signif(mean_tpm, 3), ")"),
                    chrom, win_start, win_end,
                    plot_start = as.integer(tss - PLOT_HALF),
                    plot_end   = as.integer(tss + PLOT_HALF),
                    highlight_start = as.integer(tss - 25L),
                    highlight_end   = as.integer(tss + 25L),
                    entropy, n_reads, n_bin1, n_bin2, n_bin3, n_bin4,
                    expr_bin)]
setorder(sel_tss, expr_bin, -entropy)

## 30-region driver table ----------------------------------------------------
sel_fire <- fire[, .(outname = paste(chrom, region_start, region_end,
                                     sep = "_"),
                     region_type = paste0("fire_", direction),
                     plot_group = paste0("fire_", direction),
                     window_id = region_id,
                     label = paste0(region_id, " (FIRE ", direction,
                                    ", p=", signif(min_p, 3), ")"),
                     chrom, win_start, win_end,
                     plot_start = as.integer(midpoint - PLOT_HALF),
                     plot_end   = as.integer(midpoint + PLOT_HALF),
                     highlight_start = region_start,
                     highlight_end   = region_end,
                     entropy, n_reads, n_bin1, n_bin2, n_bin3, n_bin4,
                     expr_bin = NA_character_)]
setorder(sel_fire, region_type, -entropy, na.last = TRUE)

selected <- rbind(sel_tss, sel_fire)
stopifnot(!anyDuplicated(selected$outname))
out_sel <- file.path(TAB_DIR, "selected_regions.tsv")
fwrite(selected, out_sel, sep = "\t")
cat("wrote", out_sel, ":", nrow(selected), "regions\n")
print(selected[, .(outname, region_type, label, entropy = round(entropy, 3),
                   n_reads)])
