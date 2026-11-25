# Chromatin states of intergenic Pol II consensus peaks vs Pol II promoter (TSS) peaks in LCLs.
#
# Usage (or submit all three steps with run_intergenic_vs_tss_pol2_chromatin_states.sh):
#   Rscript intergenic_vs_tss_pol2_chromatin_states.r states      # regions, cCRE, ChromHMM tables
#   Rscript intergenic_vs_tss_pol2_chromatin_states.r m6a chr21   # m6A for one chromosome
#   Rscript intergenic_vs_tss_pol2_chromatin_states.r summarize   # m6A summaries (after all chromosomes)
#
# Outputs in tables/:
#   states     chromatin_regions.tsv.gz, ccre_summary.tsv, chromhmm_center_summary.tsv,
#              chromhmm_bp_summary.tsv, rank_bin_summary.tsv
#   m6a        m6a/<chrom>.tsv.gz
#   summarize  m6a_summary.tsv, m6a_by_state_summary.tsv, fire_fp_grid_summary.tsv
#
# No hard footprint classes: accessible intergenic peaks (>= 50% of the peak inside the LCL
# FIRE union) are ranked by FIRE frequency and, separately, by 30-60 bp footprint frequency
# (pol2_consensus_fire_footprints.py, counted in the same 500-bp windows). Summaries are per
# axis (FIRE / footprint) x bin: promoter, D1 (lowest) ... D10 (highest decile), closed.
# Pol II ChIP (a ChIP summit within +/- HALF_WINDOW bp of the centroid) is reported per bin,
# not used for grouping.
#
# Every region is a window of CENTER +/- HALF_WINDOW bp around
#   intergenic  Langen et al. 2023 consensus centroid (mean Pol II summit)
#   promoter    ENCFF886PSD Pol II IDR peak summit (strongest per peak), for peaks whose
#               summit is within +/- TSS_SUMMIT_DIST bp of a canonical GENCODE v46 TSS
# Annotations per region:
#   ENCODE cCRE  class at the centre (Langen et al. rule) and class with the largest
#                overlap inside the window
#   E116 18-state ChromHMM  state at the centre and bp fraction per state in the window
#                (hg38 liftover; bases with no state are NA, not Quies)
#   m6A          pooled 31-sample m6A proportion over full-span fibers
#                (assemble_region_m6a from parsing_functions: calls / (fibers x positions
#                with >= 1 call), no A/T denominator)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

CODE <- "/project/spott/cshan/fiber-seq/code"
ANNOT <- "/project/spott/cshan/annotations"
OUT <- "/project/spott/cshan/fiber-seq/LCL_project/Pol2/intergenic_Pol2/tables"
POOLED <- file.path(OUT, "pooled_per_peak.tsv.gz")
POL2_CHIP <- file.path(ANNOT, "LCL_public_data/Pol2_chip/ENCFF886PSD.bed.gz")
TSS <- file.path(ANNOT, "gencodev46_Ensembl_canonical_TSS.bed")
CCRE <- file.path(ANNOT, "GRCh38-cCREs.bed")
CHROMHMM <- file.path(ANNOT, "LCL_public_data/roadmap_E116_chromHMM/hg38lift/E116_18_core_K27ac_hg38lift_dense.bed.gz")
SAMPLES <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv"
CHROMS <- paste0("chr", c(1:22, "X", "Y"))
HALF_WINDOW <- 250L
TSS_SUMMIT_DIST <- 250L   # promoter: Pol II ChIP summit within this distance of a canonical TSS
N_CLOSED <- 5000L          # closed intergenic peaks sampled for m6A
BINS <- c("promoter", paste0("D", 1:10), "closed")
CCRE_CLASSES <- c("PLS", "pELS", "dELS", "CA-H3K4me3", "CA-CTCF", "CA-TF", "CA", "TF", "none")

as_gr <- function(x) GRanges(x$chrom, IRanges(x$start + 1L, x$end))

#################################################
# region sets
#################################################
intergenic_regions <- function() {
  d <- fread(POOLED)
  setorder(d, cid)
  d[, group := fifelse(accessible, "intergenic accessible", "intergenic closed")]
  # percentile rank (ties averaged) and decile among accessible peaks with covering fibers
  acc <- d$accessible & !is.na(d$fire_freq)
  d[acc, `:=`(fire_pct = frank(fire_freq, ties.method = "average") / .N,
              fp_pct = frank(fp_freq, ties.method = "average") / .N)]
  d[, `:=`(fire_decile = as.integer(ceiling(10 * fire_pct)), fp_decile = as.integer(ceiling(10 * fp_pct)))]
  d[, .(region_id, chrom, center, group, category, fire_freq, fp_freq, in_pol2_chip,
        fire_pct, fp_pct, fire_decile, fp_decile)]
}

promoter_regions <- function() {
  pk <- fread(POL2_CHIP, header = FALSE)[, .(chrom = V1, start = V2, end = V3, signal = V7, summit = V10)]
  # narrowPeak repeats a peak once per summit: keep the strongest summit per peak
  pk <- pk[chrom %in% CHROMS][order(-signal)]
  pk <- pk[!duplicated(pk, by = c("chrom", "start", "end"))]
  pk[, center := start + summit]
  tss <- fread(TSS, header = FALSE)[, .(chrom = V1, start = pmax(V2 - TSS_SUMMIT_DIST, 0L),
                                        end = V2 + TSS_SUMMIT_DIST + 1L)]
  at_tss <- overlapsAny(as_gr(pk[, .(chrom, start = center, end = center + 1L)]), as_gr(tss))
  pk <- pk[at_tss]
  pk[, .(region_id = sprintf("promoter_%s:%d-%d", chrom, start, end), chrom,
         center, group = "promoter Pol II", category = "promoter")]
}

# one row per region and axis: bin = promoter, D1..D10 (accessible, by rank) or closed
rank_bins <- function(regions) {
  b <- rbindlist(lapply(c(FIRE = "fire_decile", footprint = "fp_decile"), function(col)
    regions[, .(region_id, bin = fcase(group == "promoter Pol II", "promoter",
                                       group == "intergenic closed", "closed",
                                       !is.na(get(col)), paste0("D", get(col))))]), idcol = "axis")
  b[!is.na(bin)]
}

region_windows <- function(regions) {
  regions[, `:=`(start = pmax(center - HALF_WINDOW, 0L), end = center + HALF_WINDOW)]
  regions[]
}

#################################################
# ENCODE cCREs
#################################################
ccre_annotation <- function(regions) {
  ccre <- fread(CCRE, header = FALSE, col.names = c("chrom", "start", "end", "dhs", "ccre_id", "class"))
  cgr <- as_gr(ccre)
  centre <- GRanges(regions$chrom, IRanges(regions$center + 1L, width = 1L))
  hit <- findOverlaps(centre, cgr, select = "first")
  regions[, ccre_center := fifelse(is.na(hit), "none", ccre$class[hit])]
  ov <- findOverlaps(as_gr(regions), cgr)
  w <- width(pintersect(as_gr(regions)[queryHits(ov)], cgr[subjectHits(ov)]))
  best <- data.table(q = queryHits(ov), class = ccre$class[subjectHits(ov)], w = w)[
    order(q, -w)][, .SD[1], by = q]
  regions[, ccre_best := "none"][best$q, ccre_best := best$class]
  regions[]
}

#################################################
# E116 18-state ChromHMM
#################################################
read_chromhmm <- function() {
  fread(CHROMHMM, skip = 1, header = FALSE, select = 1:4, col.names = c("chrom", "start", "end", "state"))
}

chromhmm_annotation <- function(regions, hmm = read_chromhmm()) {
  hgr <- as_gr(hmm)
  centre <- GRanges(regions$chrom, IRanges(regions$center + 1L, width = 1L))
  hit <- findOverlaps(centre, hgr, select = "first")
  regions[, state_center := hmm$state[hit]]            # NA = not covered by the liftover
  ov <- findOverlaps(as_gr(regions), hgr)
  bp <- data.table(region_id = regions$region_id[queryHits(ov)], state = hmm$state[subjectHits(ov)],
                   bp = width(pintersect(as_gr(regions)[queryHits(ov)], hgr[subjectHits(ov)])))
  bp <- bp[, .(bp = sum(bp)), by = .(region_id, state)]
  list(regions = regions[], state_bp = bp)
}

#################################################
# m6A proportion (parsing_functions)
#################################################
m6a_proportion <- function(regions, sample_table) {
  rbindlist(lapply(seq_len(nrow(regions)), function(i) {
    r <- data.frame(chr = regions$chrom[i], start = regions$start[i], end = regions$end[i],
                    region_id = regions$region_id[i])
    r$analysis_start <- r$start + 1L
    r$analysis_end <- r$end
    a <- tryCatch(assemble_region_m6a(sample_table, r, modality = "m6a", full_span = TRUE, reuse = FALSE),
                  error = function(e) NULL)
    if (is.null(a)) return(data.table(region_id = r$region_id, n_fibers = NA_integer_,
                                      n_positions = NA_integer_, m6a_prop = NA_real_, m6a_prop_fiber_median = NA_real_))
    m <- a$met_mat
    per_fiber <- Matrix::rowSums(m) / ncol(m)
    data.table(region_id = r$region_id, n_fibers = nrow(m), n_positions = ncol(m),
               m6a_prop = sum(m) / (nrow(m) * ncol(m)), m6a_prop_fiber_median = median(per_fiber))
  }))
}

#################################################
# step "states": region sets, ENCODE cCRE class, E116 18-state ChromHMM
#################################################
run_states <- function() {
  regions <- region_windows(rbind(promoter_regions(), intergenic_regions(), fill = TRUE))
  stopifnot(!anyDuplicated(regions$region_id))
  regions <- ccre_annotation(regions)
  hmm <- chromhmm_annotation(regions)
  regions <- hmm$regions

  # m6A is run on every promoter and accessible intergenic peak, plus N_CLOSED closed peaks
  set.seed(1)
  closed <- regions[group == "intergenic closed", sample(region_id, N_CLOSED)]
  regions[, m6a_selected := group != "intergenic closed" | region_id %in% closed]
  fwrite(regions, file.path(OUT, "chromatin_regions.tsv.gz"), sep = "\t")

  b <- merge(rank_bins(regions), regions, by = "region_id")
  b[, bin := factor(bin, BINS)]
  pct <- function(x, levels) {
    t <- dcast(x[, .N, by = .(axis, bin, value)][, pct := 100 * N / sum(N), by = .(axis, bin)],
               axis + bin ~ value, value.var = "pct", fill = 0)
    setcolorder(t, c("axis", "bin", intersect(levels, names(t))))
    merge(t, x[, .(n = .N), by = .(axis, bin)], by = c("axis", "bin"))[order(axis, bin)]
  }

  ccre <- rbind(
    cbind(rule = "center", pct(b[, .(axis, bin, value = ccre_center)], CCRE_CLASSES)),
    cbind(rule = "largest_overlap", pct(b[, .(axis, bin, value = ccre_best)], CCRE_CLASSES)), fill = TRUE)
  fwrite(ccre, file.path(OUT, "ccre_summary.tsv"), sep = "\t")

  state_levels <- read_chromhmm()[, unique(state)]
  state_levels <- state_levels[order(as.integer(sub("_.*", "", state_levels)))]
  center <- pct(b[, .(axis, bin, value = fifelse(is.na(state_center), "no_state", state_center))],
                c(state_levels, "no_state"))
  fwrite(center, file.path(OUT, "chromhmm_center_summary.tsv"), sep = "\t")

  # bp-weighted: share of window bases in each state (denominator = all window bases, so the
  # remainder is bases with no lifted state)
  bp <- merge(hmm$state_bp, b[, .(region_id, axis, bin)], by = "region_id", allow.cartesian = TRUE)
  tot <- b[, .(total = sum(end - start)), by = .(axis, bin)]
  bp <- merge(bp[, .(bp = sum(bp)), by = .(axis, bin, state)], tot, by = c("axis", "bin"))
  bp <- dcast(bp[, pct := 100 * bp / total], axis + bin ~ state, value.var = "pct", fill = 0)
  setcolorder(bp, c("axis", "bin", intersect(state_levels, names(bp))))
  fwrite(bp[order(axis, bin)], file.path(OUT, "chromhmm_bp_summary.tsv"), sep = "\t")

  # Pol II ChIP overlap per bin (promoters are ChIP peaks by definition: NA)
  ranks <- b[, .(n = .N, median_fire_freq = median(fire_freq), median_fp_freq = median(fp_freq),
                 frac_pol2_chip = mean(in_pol2_chip)), by = .(axis, bin)][order(axis, bin)]
  fwrite(ranks, file.path(OUT, "rank_bin_summary.tsv"), sep = "\t")

  options(width = 250)
  print(ranks, digits = 3)
  cat(sprintf("%d regions; %d selected for m6A\n", nrow(regions), sum(regions$m6a_selected)))
}

#################################################
# per chromosome
#################################################
run_m6a_chrom <- function(target_chrom) {
  source(file.path(CODE, "parsing_functions/parsing_footprints_functions.r"))
  regions <- fread(file.path(OUT, "chromatin_regions.tsv.gz"))[m6a_selected == TRUE & chrom == target_chrom]
  sample_table <- fread(SAMPLES)[, .(sample_name, fire_dir)]
  stopifnot(nrow(sample_table) == 31L)
  res <- m6a_proportion(regions, sample_table)
  dir.create(file.path(OUT, "m6a"), showWarnings = FALSE)
  fwrite(res, file.path(OUT, "m6a", paste0(target_chrom, ".tsv.gz")), sep = "\t")
  cat(sprintf("%s: %d regions, %d with m6A\n", target_chrom, nrow(res), sum(!is.na(res$m6a_prop))))
}

#################################################
# summaries: by group and by group x E116 centre state
#################################################
summarize_m6a <- function() {
  regions <- fread(file.path(OUT, "chromatin_regions.tsv.gz"), na.strings = "")[m6a_selected == TRUE]
  files <- file.path(OUT, "m6a", paste0(unique(regions$chrom), ".tsv.gz"))
  stopifnot(all(file.exists(files)))
  m <- merge(regions, rbindlist(lapply(files, fread)), by = "region_id")
  stopifnot(nrow(m) == nrow(regions))
  b <- merge(rank_bins(m), m, by = "region_id")
  b[, bin := factor(bin, BINS)]

  by_bin <- b[, .(n = .N, n_with_m6a = sum(!is.na(m6a_prop)),
                  median_fibers = as.double(median(n_fibers, na.rm = TRUE)),
                  mean_m6a_prop = mean(m6a_prop, na.rm = TRUE),
                  median_m6a_prop = median(m6a_prop, na.rm = TRUE),
                  q25 = quantile(m6a_prop, 0.25, na.rm = TRUE),
                  q75 = quantile(m6a_prop, 0.75, na.rm = TRUE)), by = .(axis, bin)][order(axis, bin)]
  fwrite(by_bin, file.path(OUT, "m6a_summary.tsv"), sep = "\t")

  # long format: axis x bin x E116 centre state, strata with >= 20 valid regions
  by_state <- b[!is.na(state_center), .(n_valid = sum(!is.na(m6a_prop)),
                                        median_m6a_prop = median(m6a_prop, na.rm = TRUE)),
                by = .(axis, bin, state_center)][n_valid >= 20]
  by_state <- by_state[order(axis, bin, as.integer(sub("_.*", "", state_center)))]
  fwrite(by_state, file.path(OUT, "m6a_by_state_summary.tsv"), sep = "\t")

  # FIRE decile x footprint decile among accessible peaks
  grid <- m[!is.na(fire_decile) & !is.na(fp_decile),
            .(n = .N, frac_pol2_chip = mean(in_pol2_chip), median_m6a_prop = median(m6a_prop, na.rm = TRUE)),
            by = .(fire_decile, fp_decile)][order(fire_decile, fp_decile)]
  fwrite(grid, file.path(OUT, "fire_fp_grid_summary.tsv"), sep = "\t")

  options(width = 250)
  print(by_bin, digits = 3)
}

#################################################
# command line
#################################################
args <- commandArgs(trailingOnly = TRUE)
step <- args[1]
if (identical(step, "states")) {
  run_states()
} else if (identical(step, "m6a")) {
  stopifnot(length(args) == 2L, args[2] %in% CHROMS)
  run_m6a_chrom(args[2])
} else if (identical(step, "summarize")) {
  summarize_m6a()
} else {
  stop("usage: Rscript intergenic_vs_tss_pol2_chromatin_states.r states | m6a <chrom> | summarize")
}
