# Chromatin states of intergenic Pol II consensus peaks vs Pol II promoter (TSS) peaks in LCLs.
#
# Usage (or submit all three steps with run_intergenic_vs_tss_pol2_chromatin_states.sh):
#   Rscript intergenic_vs_tss_pol2_chromatin_states.r states      # regions, cCRE, ChromHMM tables
#   Rscript intergenic_vs_tss_pol2_chromatin_states.r m6a chr21   # m6A for one chromosome
#   Rscript intergenic_vs_tss_pol2_chromatin_states.r summarize   # m6A summaries (after all chromosomes)
#
# Outputs in tables/:
#   states     chromatin_regions.tsv.gz, ccre_summary.tsv, chromhmm_center_summary.tsv,
#              chromhmm_bp_summary.tsv
#   m6a        m6a/<chrom>.tsv.gz
#   summarize  m6a_summary.tsv, m6a_by_state_summary.tsv
#
# Every region is a window of CENTER +/- HALF_WINDOW bp around
#   intergenic  Langen et al. 2023 consensus centroid (mean Pol II summit)
#   promoter    ENCFF886PSD Pol II IDR peak summit, for peaks within +/-1 kb of a
#               canonical GENCODE v46 TSS
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
OUT <- "/project/spott/cshan/fiber-seq/LCL_project/pausing/intergenic_Pol2/tables"
POOLED <- file.path(OUT, "pooled_per_peak.tsv.gz")
CONSENSUS <- file.path(ANNOT, "LCL_public_data/intergenic_pol2_Langen_2023/pol2_consensuses/consensuses.bed")
POL2_CHIP <- file.path(ANNOT, "LCL_public_data/Pol2_chip/ENCFF886PSD.bed.gz")
TSS <- file.path(ANNOT, "gencodev46_Ensembl_canonical_TSS.bed")
CCRE <- file.path(ANNOT, "GRCh38-cCREs.bed")
CHROMHMM <- file.path(ANNOT, "LCL_public_data/roadmap_E116_chromHMM/hg38lift/E116_18_core_K27ac_hg38lift_dense.bed.gz")
SAMPLES <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv"
CHROMS <- paste0("chr", c(1:22, "X", "Y"))
HALF_WINDOW <- 250L
TSS_FLANK <- 1000L
N_CLOSED <- 5000L          # closed intergenic peaks sampled for m6A
GROUPS <- c("promoter Pol II", "intergenic acc+fp+ChIP", "intergenic acc+fp, no ChIP",
            "intergenic acc, no fp", "intergenic closed")
CCRE_CLASSES <- c("PLS", "pELS", "dELS", "CA-H3K4me3", "CA-CTCF", "CA-TF", "CA", "TF", "none")

as_gr <- function(x) GRanges(x$chrom, IRanges(x$start + 1L, x$end))

#################################################
# region sets
#################################################
intergenic_regions <- function() {
  d <- fread(POOLED)
  cen <- fread(CONSENSUS, select = c(4, 7), col.names = c("cid", "center"))
  d <- merge(d, cen, by = "cid")
  d[, group := fcase(accessible & has_fp & in_pol2_chip, "intergenic acc+fp+ChIP",
                     accessible & has_fp, "intergenic acc+fp, no ChIP",
                     accessible, "intergenic acc, no fp",
                     default = "intergenic closed")]
  d[, .(region_id, chrom, center, group, category)]
}

promoter_regions <- function() {
  pk <- fread(POL2_CHIP, header = FALSE)[, .(chrom = V1, start = V2, end = V3, signal = V7, summit = V10)]
  # narrowPeak repeats a peak once per summit: keep the strongest summit per peak
  pk <- pk[chrom %in% CHROMS][order(-signal)]
  pk <- pk[!duplicated(pk, by = c("chrom", "start", "end"))]
  tss <- fread(TSS, header = FALSE)[, .(chrom = V1, start = pmax(V2 - TSS_FLANK, 0L), end = V2 + TSS_FLANK + 1L)]
  at_tss <- overlapsAny(as_gr(pk), as_gr(tss))
  pk <- pk[at_tss]
  pk[, .(region_id = sprintf("promoter_%s:%d-%d", chrom, start, end), chrom,
         center = start + summit, group = "promoter Pol II", category = "promoter")]
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
  regions <- region_windows(rbind(promoter_regions(), intergenic_regions()))
  stopifnot(!anyDuplicated(regions$region_id))
  regions <- ccre_annotation(regions)
  hmm <- chromhmm_annotation(regions)
  regions <- hmm$regions

  # m6A is run on every promoter and accessible intergenic peak, plus N_CLOSED closed peaks
  set.seed(1)
  closed <- regions[group == "intergenic closed", sample(region_id, N_CLOSED)]
  regions[, m6a_selected := group != "intergenic closed" | region_id %in% closed]
  fwrite(regions, file.path(OUT, "chromatin_regions.tsv.gz"), sep = "\t")

  regions[, group := factor(group, GROUPS)]
  pct <- function(x, levels) {
    t <- dcast(x[, .N, by = .(group, value)][, pct := 100 * N / sum(N), by = group],
               group ~ value, value.var = "pct", fill = 0)
    setcolorder(t, c("group", intersect(levels, names(t))))
    t[, n := as.integer(table(x$group)[as.character(group)])][]
  }

  ccre <- rbind(
    cbind(rule = "center", pct(regions[, .(group, value = ccre_center)], CCRE_CLASSES)),
    cbind(rule = "largest_overlap", pct(regions[, .(group, value = ccre_best)], CCRE_CLASSES)), fill = TRUE)
  fwrite(ccre, file.path(OUT, "ccre_summary.tsv"), sep = "\t")

  state_levels <- read_chromhmm()[, unique(state)]
  state_levels <- state_levels[order(as.integer(sub("_.*", "", state_levels)))]
  center <- pct(regions[, .(group, value = fifelse(is.na(state_center), "no_state", state_center))],
                c(state_levels, "no_state"))
  fwrite(center, file.path(OUT, "chromhmm_center_summary.tsv"), sep = "\t")

  # bp-weighted: share of window bases in each state (denominator = all window bases, so the
  # remainder is bases with no lifted state)
  bp <- merge(hmm$state_bp, regions[, .(region_id, group, w = end - start)], by = "region_id")
  tot <- regions[, .(total = sum(end - start)), by = group]
  bp <- merge(bp[, .(bp = sum(bp)), by = .(group, state)], tot, by = "group")
  bp <- dcast(bp[, pct := 100 * bp / total], group ~ state, value.var = "pct", fill = 0)
  setcolorder(bp, c("group", intersect(state_levels, names(bp))))
  fwrite(bp, file.path(OUT, "chromhmm_bp_summary.tsv"), sep = "\t")

  options(width = 250)
  print(ccre, digits = 3)
  print(center, digits = 3)
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
  m[, group := factor(group, GROUPS)]

  by_group <- m[, .(n = .N, n_with_m6a = sum(!is.na(m6a_prop)),
                    median_fibers = as.double(median(n_fibers, na.rm = TRUE)),
                    mean_m6a_prop = mean(m6a_prop, na.rm = TRUE),
                    median_m6a_prop = median(m6a_prop, na.rm = TRUE),
                    q25 = quantile(m6a_prop, 0.25, na.rm = TRUE),
                    q75 = quantile(m6a_prop, 0.75, na.rm = TRUE)), by = group][order(group)]
  fwrite(by_group, file.path(OUT, "m6a_summary.tsv"), sep = "\t")

  by_state <- m[!is.na(state_center), .(n = .N, median_m6a_prop = median(m6a_prop, na.rm = TRUE)),
                by = .(group, state_center)]
  by_state <- dcast(by_state[n >= 20], state_center ~ group, value.var = "median_m6a_prop")
  by_state <- by_state[order(as.integer(sub("_.*", "", state_center)))]
  fwrite(by_state, file.path(OUT, "m6a_by_state_summary.tsv"), sep = "\t")

  options(width = 250)
  print(by_group, digits = 3)
  print(by_state, digits = 3)
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
