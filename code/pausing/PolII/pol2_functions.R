# Paths, parameters and helpers for the PolII example genes (gene_pausing_example.rmd) and
# their plots (polII_pausing_plots.r).
#
# Table coordinates are 1-based inclusive. Read-level positions are offsets from
# the pause_index.tsv TSS (RefSeq isoform with the strongest Pol II TSSR) in the
# direction of transcription (minus-strand genes are flipped), so offset 0 is the
# TSS base and positive offsets are downstream.

suppressPackageStartupMessages({
  library(GenomicRanges); library(data.table)   # data.table last: its shift() masks IRanges'
})

PROJECT <- "/project/spott/cshan/fiber-seq"
ANN     <- "/project/spott/cshan/annotations"
LCL_PUB <- file.path(ANN, "LCL_public_data")
SHARED  <- "/project/spott/1_Shared_projects/LCL_Fiber_seq"

PLOTTING     <- file.path(PROJECT, "code/parsing_functions/plotting_functions.r")
PARSING      <- file.path(PROJECT, "code/parsing_functions/parsing_footprints_functions.r")
PI_TABLE     <- file.path(PROJECT, "LCL_project/pausing/pause_index.tsv")   # code/pausing/cal_pausing_index/pause_index.R
OUT_DIR      <- file.path(PROJECT, "LCL_project/pausing/PolII")
DOMINANT_TSS <- file.path(ANN, "fantom5_LCL/fantom5.hg38.LCL.dominant_TSS.bed.gz")
SAMPLES      <- file.path(SHARED, "Data/LCL_sample_metatable_merged_samples_31samples.csv")
# FiberHMM v1 footprints, all samples combined: closed size bins, together 10-80 bp (a footprint on
# a bin edge is in two files); plain gzip sorted by start, no index. The FiberHMM_v2 combined
# tracks are not usable here: one row per read (record span), not per footprint.
COMBINED_ROOT <- file.path(SHARED, "FiberHMM/merged/combined/joint_trained_tracks")
COMBINED_BINS <- c("10-30bp", "30-60bp", "60-80bp")
# fibertools nucleosome calls per sample (ft extract --nuc, BED12 with 1-bp sentinel blocks, tabix)
nuc_file <- function(sample, chrom)
  file.path(SHARED, "FIRE/results", sample, "extracted_results/nuc_by_chr",
            sprintf("%s.ft_extracted_nuc.%s.bed.gz", sample, chrom))

CHAIN  <- file.path(ANN, "hg19ToHg38.over.chain")
GROSEQ <- c(`+` = file.path(LCL_PUB, "GRO-seq/GSM1480326_GM12878_GROseq_plus.bigWig"),     # hg19
            `-` = file.path(LCL_PUB, "GRO-seq/GSM1480326_GM12878_GROseq_minus.bigWig"))
DNASE  <- file.path(LCL_PUB, "DNase.ENCFF743ULW.bigWig")    # GM12878 ENCSR000EMT, GRCh38
CAGE   <- file.path(ANN, "fantom5_LCL/fantom5.hg38.LCL.consensus.TSS.bed.gz")

AUTOSOMES    <- paste0("chr", 1:22)   # FiberHMM covers chr1-22 only
FLANK        <- 1000L                 # read window = TSS +/- FLANK
MIN_READS    <- 50L                   # example genes: fibers with a nucleosome record in the window
OTHER_TSS_BP <- 1000L                 # example genes: no other dominant TSS this close
N_EXAMPLES   <- 5L                    # example genes per PI bin
PI_BINS      <- c("low", "mid", "high")
CONTROL_GENES <- c("GADD45B", "HSPA1A", "HSPA1B", "HSPA8", "HSP90AA1", "DNAJB1",
                   "FOS", "JUNB", "EGR1", "IER2", "MYC")   # textbook paused genes
PLUS1_OFFSETS <- 0:299                # +1 nucleosome = first nucleosome-call base in here
SMOOTH_BP    <- 25L                   # running mean of the per-position fractions

# FiberHMM footprint size classes: assign_size_class() bins of parsing_footprints_functions.r
# matching COMBINED_BINS. A size on a shared break goes to the upper bin (assign_size_class()
# assigns bins in order, so 30 bp is 30-60 and 60 bp is 60-80, although its comment says lower).
# PPP candidates are the 30-60 bp class, PIC candidates the 60-80 bp class.
TF_BREAKS  <- c(10L, 30L, 60L, 80L)
CALL_CLASS <- c(PPP = "30-60 bp", PIC = "60-80 bp")

# PIC green / PPP orange (categorical slots 3 and 2; PPP vs PIC CVD dE 9.2); other footprints use
# MOLECULE_TRACK_COLORS (plotting_functions.r), all-pairs CVD dE >= 9.2 against these two. A pure
# green (#008300) is CVD dE 3.2 from this orange, so it is not used
CLASS_COLORS  <- c(PIC = "#1baf7a", PPP = "#eb6834")
INK <-c(primary = "#0b0b0b", secondary = "#52514e", muted = "#898781", grid = "#e1e0d9")


# TSS-relative first (lo) and last (hi) base of BED intervals [start0, end0)
to_relative <- function(start0, end0, tss, strand) {
  plus <- strand == "+"
  list(lo = as.integer(ifelse(plus, start0 + 1 - tss, tss - end0)),
       hi = as.integer(ifelse(plus, end0 - tss, tss - start0 - 1)))
}

# GM12878 GRO-seq of one chromosome lifted hg19 -> hg38:
# list(`+`, `-`) of GRanges with per-base read counts (minus stored negative)
lift_groseq <- function(chrom, chain) {
  lapply(GROSEQ, function(f) {
    len <- GenomeInfoDb::seqlengths(rtracklayer::BigWigFile(f))[[chrom]]
    gr <- rtracklayer::import(f, which = GRanges(chrom, IRanges(1L, len)))
    gr$score <- abs(gr$score)
    lo <- unlist(rtracklayer::liftOver(gr, chain))
    lo[as.character(seqnames(lo)) == chrom]
  })
}

# fibertools nucleosome calls of every LCL sample in chrom:start..end (1-based), pooled. One
# alignment per read (longest; palindromic pairs collapse), sentinel blocks dropped; a read in
# several samples is kept in the first.
# Output: list(reads = RID, sample, start0, end0 (read span); nucs = RID, start0, end0)
window_nucleosomes <- function(chrom, start, end) {
  region <- GRanges(chrom, IRanges(start, end))
  per_sample <- lapply(fread(SAMPLES)$sample_name, function(s) {
    f <- nuc_file(s, chrom)
    d <- read_tabix_region(f, region)
    if (!nrow(d)) return(NULL)
    names(d)[seq_along(FT_BED12_COLS)] <- FT_BED12_COLS
    d <- keep_longest_alignment(d)
    b <- convert_ft_bed12_to_bed6(d, format = "bed12_fibertools", source = f)
    list(reads = data.table(RID = d$RID, sample = s, start0 = d$start, end0 = d$end),
         nucs = data.table(RID = b$RID, sample = s, start0 = b$start, end0 = b$end))
  })
  reads <- rbindlist(lapply(per_sample, `[[`, "reads"))
  reads <- reads[!duplicated(RID)]
  nucs <- rbindlist(lapply(per_sample, `[[`, "nucs"))[reads[, .(RID, sample)], on = .(RID, sample), nomatch = 0L]
  list(reads = reads, nucs = nucs[, .(RID, start0, end0)])
}

# FiberHMM combined footprints (COMBINED_BINS, 10-80 bp) overlapping chrom:start..end (1-based).
# The files have no index, so each is streamed until its records start past `end`, in parallel.
# Output: RID, start0, end0 (unique)
window_tf_footprints <- function(chrom, start, end) {
  files <- file.path(COMBINED_ROOT, chrom, sprintf("combined_%s_%s_fps.bed.gz", chrom, COMBINED_BINS))
  cmds <- sprintf("zcat %s 2>/dev/null | awk -F'\\t' -v s=%d -v e=%d '$2 >= e {exit} $3 >= s'",
                  shQuote(files), as.integer(start), as.integer(end))
  out <- parallel::mclapply(cmds, function(x) {
    d <- fread(cmd = x, header = FALSE, sep = "\t")
    if (nrow(d)) setnames(d, c("chrom", "start0", "end0", "RID"))[, .(RID, start0, end0)]
  }, mc.cores = length(cmds))
  unique(rbindlist(out))
}

# Per-read matrix in the Tullius format: one row per read (the order of `reads`),
# one column per offset -FLANK..FLANK. A footprint base holds the footprint's
# length, any other base 0.
#   reads - RID; fps - RID, lo, hi, size (one gene)
footprint_matrix <- function(reads, fps) {
  m <- matrix(0L, nrow(reads), 2L * FLANK + 1L,
              dimnames = list(reads$RID, as.character(-FLANK:FLANK)))
  fill <- function(rid, lo, hi, value) {
    lo <- pmax(lo, -FLANK); hi <- pmin(hi, FLANK)
    keep <- lo <= hi & rid %in% reads$RID
    if (!any(keep)) return(invisible())
    len <- hi[keep] - lo[keep] + 1L
    cells <- cbind(rep(match(rid[keep], reads$RID), len), sequence(len, from = lo[keep] + FLANK + 1L))
    m[cells] <<- rep(value[keep], len)
  }
  fill(fps$RID, fps$lo, fps$hi, as.integer(fps$size))
  m
}

# Fraction of the reads covering each offset (cover > 0) that have a call there (calls > 0);
# both footprint_matrix() results over the same reads
call_fraction <- function(calls, cover) {
  colSums(calls > 0 & cover > 0) / colSums(cover > 0)
}

# +1 nucleosome: first offset in PLUS1_OFFSETS covered by a nucleosome call (m = footprint_matrix()
# of the nucleosome calls), NA when there is none
plus1_offset <- function(m) {
  sub <- (m[, as.character(PLUS1_OFFSETS), drop = FALSE] > 0) * 1L
  first <- max.col(sub, ties.method = "first")
  ifelse(rowSums(sub) > 0, PLUS1_OFFSETS[first], NA_integer_)
}

# Centred running mean over SMOOTH_BP positions (NA at the edges)
running_mean <- function(x, k = SMOOTH_BP) as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))


# --- example genes (gene_pausing_example.rmd) -------------------------------------

# Nearest other gene's CAGE dominant TSS (DOMINANT_TSS, any biotype, either strand) to each gene's
# TSS; the gene itself is matched by ENSG or symbol.
#   genes - gene_id, ensg, gene_name, chrom, tss, strand
# Output: gene_id, nearest_tss_bp, nearest_tss_gene, nearest_tss_orientation
nearest_other_tss <- function(genes) {
  dom <- fread(DOMINANT_TSS)[, .(dom_ensg = sub("\\..*", "", gene_id), dom_gene = gene_name,
                                 chrom = `#chrom`, pos = end, dom_strand = strand)]
  near <- merge(genes[, .(gene_id, ensg, gene_name, chrom, tss, strand)], dom, by = "chrom",
                allow.cartesian = TRUE)
  near <- near[(is.na(ensg) | dom_ensg != ensg) & dom_gene != gene_name]
  near[, dist := abs(pos - tss)]
  setorder(near, gene_id, dist)
  near[!duplicated(gene_id), .(gene_id, nearest_tss_bp = dist, nearest_tss_gene = dom_gene,
                               nearest_tss_orientation = fifelse(dom_strand == strand,
                                                                 "same strand", "opposite strand"))]
}

# Fibers (31 samples pooled, one alignment each) with an ft nucleosome record overlapping the
# gene's TSS +/- FLANK; g = one pause_index.tsv row
window_fiber_count <- function(g) {
  nrow(window_nucleosomes(g$chrom, g$tss - FLANK, g$tss + FLANK)$reads)
}

# Analysis-set genes with their nearest other TSS and matching scores:
#   candidate     no other gene's dominant TSS within OTHER_TSS_BP
#   in_fpkm_band  FPKM in the middle 50% of the analysis set
#   distance      |PI percentile - percentile of the tertile median| + |FPKM percentile - 0.5|
#   pi_tab - pause_index.tsv
example_candidates <- function(pi_tab) {
  an <- pi_tab[analysis_set == TRUE]
  an <- merge(an, nearest_other_tss(an), by = "gene_id", all.x = TRUE)
  an[, candidate := is.na(nearest_tss_bp) | nearest_tss_bp > OTHER_TSS_BP]
  band <- an[, quantile(fpkm, c(0.25, 0.75))]
  an[, in_fpkm_band := fpkm >= band[1] & fpkm <= band[2]]
  pi_ecdf <- ecdf(an$pi)
  an[, `:=`(pct_pi = frank(pi) / .N, pct_fpkm = frank(fpkm) / .N)]
  an[, target_pct_pi := pi_ecdf(median(pi)), by = pi_bin]
  an[, distance := abs(pct_pi - target_pct_pi) + abs(pct_fpkm - 0.5)]
  an[, pi_bin := factor(pi_bin, levels = PI_BINS)]
  setorder(an, pi_bin, distance)
  an[]
}

# Per PI tertile, walk the candidates in the FPKM band from the nearest to the matching target and
# keep the first N_EXAMPLES with >= MIN_READS fibers (window_fiber_count(), counted only for the
# genes looked at); adds n_reads
pick_examples <- function(an) {
  pool <- an[candidate & in_fpkm_band]
  picks <- rbindlist(lapply(PI_BINS, function(b) {
    kept <- list()
    for (i in which(pool$pi_bin == b)) {
      n <- window_fiber_count(pool[i])
      if (n >= MIN_READS) kept[[length(kept) + 1L]] <- pool[i][, n_reads := n]
      if (length(kept) == N_EXAMPLES) break
    }
    rbindlist(kept)
  }))
  picks[, group := as.character(pi_bin)][]
}

# Status of each CONTROL_GENES gene; the highest-PI one that "passes" (high tertile, no other TSS
# within OTHER_TSS_BP, >= MIN_READS fibers; FPKM band not required) is the positive control.
# Fibers are counted only for genes that pass the other rules.
control_status <- function(pi_tab, an) {
  controls <- merge(data.table(gene_name = CONTROL_GENES),
                    pi_tab[, .(gene_name, gene_id, chrom, tss, pi, fpkm, near_other_gene, short,
                               analysis_set, pi_bin)],
                    by = "gene_name", all.x = TRUE)
  controls <- merge(controls, an[, .(gene_id, nearest_tss_bp, candidate)], by = "gene_id", all.x = TRUE)
  controls[, status := fcase(is.na(gene_id), "not a RefSeq gene in pause_index.tsv",
                             !analysis_set & near_other_gene, "another gene within 3 kb",
                             !analysis_set & short, "shorter than 1 kb",
                             !analysis_set, "not in analysis set (expression, biotype, chrom or PI)",
                             pi_bin != "high", paste(pi_bin, "PI"),
                             !candidate, "other dominant TSS within 1 kb",
                             default = "passes")]
  controls[, n_reads := NA_integer_]
  for (i in which(controls$status == "passes")) controls[i, n_reads := window_fiber_count(controls[i])]
  controls[status == "passes" & n_reads < MIN_READS, status := paste("fewer than", MIN_READS, "fibers")]
  setorder(controls, -pi)
  controls[]
}

# picks plus the positive control (first "passes" row of control_status()), if any
add_control <- function(picks, an, controls) {
  control <- controls[status == "passes"]
  if (!nrow(control)) return(picks)
  row <- an[gene_id == control$gene_id[1]][, `:=`(n_reads = control$n_reads[1], group = "control")]
  rbind(picks, row, fill = TRUE)
}

# Write the chosen genes to tables/example_genes.tsv
write_example_table <- function(picks) {
  cols <- c("group", "gene_id", "gene_name", "chrom", "tss", "strand", "tes", "tx_id", "pi", "pi_bin",
            "fpkm", "tssr_dens", "body_dens", "n_reads",
            "nearest_tss_bp", "nearest_tss_gene", "nearest_tss_orientation",
            "candidate", "in_fpkm_band", "pct_pi", "pct_fpkm", "distance")
  path <- file.path(OUT_DIR, "tables/example_genes.tsv")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(picks[, ..cols], path, sep = "\t")
  invisible(path)
}
