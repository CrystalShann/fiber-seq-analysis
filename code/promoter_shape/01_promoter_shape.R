#!/usr/bin/env Rscript
# Promoter shape of the PolII example genes from FANTOM5 GM12878 CAGE (README Steps 1-8, 10, 11):
# paths, parameters and helpers, then the pipeline run_promoter_shape().
#
# Run: bash code/promoter_shape/run_promoter_shape.sh   (login node is fine: ~1 min, < 4 GB)
# source() this file to get the helpers only (02_promoter_shape_figures.Rmd, test_promoter_shape.R);
# the pipeline runs only when the file is executed with Rscript.
#
# Shape index (SI) of Hoskins et al. 2011 (Genome Res 21:182), cut at -1 as in Hoskins (peaked /
# broad) and Schor et al. 2017 (Nat Genet 49:550; narrow / broad). Method notes and the published
# vs adapted steps are in README.md.
#
# Coordinates are 1-based inclusive: promoter start = FANTOM5 BED start + 1, end = BED end, and a
# CTSS or TSS is the BED end of its 1-bp record (as dominant_TSS in the PolII tables). Offsets run
# in the direction of transcription (minus-strand promoters are flipped): positive = downstream.

suppressPackageStartupMessages({
  library(GenomicRanges); library(data.table)   # data.table last: its shift() masks IRanges'
})

PROJECT <- "/project/spott/cshan/fiber-seq"
ANN     <- "/project/spott/cshan/annotations"
F5      <- file.path(ANN, "fantom5_LCL")
CODE_DIR <- file.path(PROJECT, "code/promoter_shape")
OUT_DIR  <- file.path(PROJECT, "LCL_project/promoter_shape")

EXAMPLE_GENES <- file.path(PROJECT, "code/pausing/PolII/example_gene_list.tsv")
DOMINANT_TSS  <- file.path(F5, "fantom5.hg38.LCL.dominant_TSS.bed.gz")
CANONICAL_TSS <- file.path(ANN, "gencodev46_Ensembl_canonical_TSS.bed")
PEAKS      <- file.path(F5, "fantom5.hg38.CAGE_peaks.bed.gz")   # hg38_fair+new_CAGE_peaks_phase1and2
PEAK_ANNOT <- file.path(F5, "hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt.gz")
PEAK_NAMES <- file.path(F5, "human_phase1and2_CAGE_Peak_name.txt.gz")
HGNC       <- file.path(F5, "hgnc_complete_set.txt")
# FANTOM5 GM12878 hCAGE biological replicates, hg38 CTSS (download links in fantom5_LCL/README.md)
SAMPLES <- c(CNhs12331 = "rep1", CNhs12332 = "rep2", CNhs12333 = "rep3")
ctss_file <- function(sample) file.path(F5, sprintf("fantom5.hg38.LCL.%s.CTSS.bed.gz", SAMPLES[[sample]]))

# not_expressed genes are excluded (user request 2026-11-13): almost none have CAGE tags
GROUPS <- c("housekeeping", "cell_type_specific")
GROUP_COLORS <- setNames(c("#2878B5", "#D77A21"), GROUPS)   # as the PolII summaries

SI_CUT          <- -1         # peaked / narrow: SI > -1; broad: SI <= -1 (both papers)
MIN_TAGS        <- 20L        # fewer aggregated tags = low coverage, left unclassified
N_RESAMPLE      <- 500L       # bootstrap resamples per promoter, each of the observed tag total
SEED            <- 20261113L  # promoter seed = SEED + FANTOM5 peak number (hg_<number>.1)
STABILITY_MIN   <- 0.9        # class kept when >= 90% of resamples fall on the same side of SI_CUT
MAX_DEPTH_RATIO <- 1.5        # replicate library sizes must agree within this ratio to sum raw counts

strand_sign <- function(strand) fifelse(strand == "+", 1L, -1L)

#####################################
# Inputs
#####################################

# Example genes with their dominant_TSS, looked up exactly as resolve_example_tss() in
# code/pausing/PolII/pol2_pausing_functions.r (unchanged, not recalculated)
#
# 1. FANTOM5 LCL dominant TSS BED by gene name; more than one hit is an error
# 2. only when absent, the GENCODE v46 canonical one-base BED, which must have exactly one hit
# 3. TSS = BED end (1-based); the source is kept as dominant_TSS_source
# Only genes of the analysed GROUPS are kept.
example_genes <- function() {
  genes <- fread(EXAMPLE_GENES)[group %in% GROUPS]
  dominant <- fread(DOMINANT_TSS)
  canonical <- fread(CANONICAL_TSS, header = FALSE,
    col.names = c("chrom", "start0", "end", "annotation", "score", "strand"))
  canonical[, c("gene_id", "tx_id", "gene_name") := tstrsplit(annotation, ";", fixed = TRUE, keep = 1:3)]
  tss <- rbindlist(lapply(genes$gene_name, function(g) {
    hit <- dominant[gene_name == g]
    if (nrow(hit) > 1L) stop("multiple dominant TSS found for ", g, call. = FALSE)
    if (nrow(hit)) return(hit[, .(tss_gene_id = gene_id, tss_chrom = `#chrom`, strand,
      dominant_TSS = as.integer(end), dominant_TSS_source = "FANTOM5 dominant TSS")])
    hit <- canonical[gene_name == g]
    if (nrow(hit) != 1L) stop(nrow(hit), " GENCODE canonical TSSs found for ", g, call. = FALSE)
    stopifnot(hit$end - hit$start0 == 1L)
    hit[, .(tss_gene_id = gene_id, tss_chrom = chrom, strand,
            dominant_TSS = as.integer(end), dominant_TSS_source = "GENCODE v46 canonical TSS")]
  }))
  out <- cbind(genes[, .(gene_group = group, gene_name, gene_id, chr = chrom)], tss)
  stopifnot(out$gene_group %in% GROUPS, out$chr == out$tss_chrom,
            out$gene_id == sub("\\..*", "", out$tss_gene_id), !anyDuplicated(out$gene_id))
  out[, .(gene_group, gene_name, gene_id, chr, strand, dominant_TSS, dominant_TSS_source)]
}

# FANTOM5 hg38 CAGE peaks with their promoter name and the genes FANTOM5 assigns them to
#
# 1. peak BED (0-based start -> 1-based start); promoter_id = the hg38 peak ID (hg_<n>.1)
# 2. promoter_name: the p<n>@<gene> name of the peak-name file (New_TSS_name, else Old_TSS_name)
# 3. links: peak -> Ensembl gene through the HGNC IDs of the hg38 peak annotation file
#    (one row per peak x gene; a peak can belong to several genes)
fantom5_peaks <- function() {
  peaks <- fread(PEAKS, header = FALSE, select = c(1:4, 6),
                 col.names = c("chr", "start0", "end", "peak_full", "strand"))
  peaks[, `:=`(promoter_id = sub(".*;", "", peak_full), start = start0 + 1L)]
  stopifnot(!anyDuplicated(peaks$promoter_id), grepl("^hg_[0-9]+\\.[0-9]+$", peaks$promoter_id))
  pnames <- fread(PEAK_NAMES, sep = "\t", quote = "", col.names = c("peak_full", "new_name", "old_name"))
  pnames[, promoter_name := fifelse(grepl("@", new_name), new_name,
                                    fifelse(grepl("@", old_name), old_name, NA_character_))]
  peaks <- merge(peaks, pnames[, .(peak_full, promoter_name)], by = "peak_full", all.x = TRUE)

  annot <- fread(PEAK_ANNOT, sep = "\t", quote = "", select = c(1L, 5L), col.names = c("peak_full", "hgnc"))
  hgnc <- fread(HGNC, select = c("hgnc_id", "symbol", "ensembl_gene_id"), quote = "")
  links <- annot[hgnc != "", .(hgnc_id = unlist(strsplit(hgnc, "[ ,]+"))), by = peak_full]
  links <- merge(links, hgnc, by = "hgnc_id")[ensembl_gene_id != ""]
  links <- merge(links, peaks[, .(peak_full, promoter_id)], by = "peak_full")
  list(peaks = peaks[, .(promoter_id, promoter_name, chr, start, end, strand)],
       links = unique(links[, .(promoter_id, gene_id = ensembl_gene_id, symbol)]))
}

# One GM12878 replicate: 1-based position, tag count and strand of every CTSS
read_ctss <- function(sample) {
  x <- fread(ctss_file(sample), header = FALSE, select = c(1:3, 5:6),
             col.names = c("chr", "start0", "pos", "count", "strand"))
  stopifnot(x$pos - x$start0 == 1L, x$count > 0L, x$strand %in% c("+", "-"),
            !anyDuplicated(x[, .(chr, pos, strand)]))
  x[, start0 := NULL]
  x
}

#####################################
# Step 1. Candidate promoters of each gene
#####################################

# Gene x peak candidates: same-strand peaks FANTOM5 assigns to the gene, plus the peak that
# contains the gene's dominant_TSS on its strand (even when FANTOM5 assigns it elsewhere)
#
# 1. FANTOM5 links of the gene (annotated_to_gene)
# 2. the peak that overlaps dominant_TSS on the gene's strand (contains_dominant_TSS); FANTOM5
#    peaks on one strand do not overlap, so a TSS lies in at most one peak
# 3. same_strand: peak on the gene's chromosome and strand
promoter_candidates <- function(genes, peaks, links) {
  annotated <- links[gene_id %in% genes$gene_id, .(gene_id, promoter_id, annotated_to_gene = TRUE)]
  tss <- GRanges(genes$chr, IRanges(genes$dominant_TSS, width = 1L), strand = genes$strand)
  hits <- findOverlaps(tss, GRanges(peaks$chr, IRanges(peaks$start, peaks$end), strand = peaks$strand))
  stopifnot(!anyDuplicated(queryHits(hits)))
  containing <- data.table(gene_id = genes$gene_id[queryHits(hits)],
                           promoter_id = peaks$promoter_id[subjectHits(hits)], contains_dominant_TSS = TRUE)
  cand <- merge(unique(annotated), containing, by = c("gene_id", "promoter_id"), all = TRUE)
  cand[is.na(annotated_to_gene), annotated_to_gene := FALSE]
  cand[is.na(contains_dominant_TSS), contains_dominant_TSS := FALSE]
  cand <- merge(cand, peaks, by = "promoter_id")
  cand <- merge(cand, genes[, .(gene_id, gene_name, gene_group, gene_chr = chr, gene_strand = strand,
                                dominant_TSS, dominant_TSS_source)], by = "gene_id")
  cand[, same_strand := chr == gene_chr & strand == gene_strand]
  cand[, c("gene_chr", "gene_strand") := NULL]
  cand[]
}

#####################################
# Step 2. Nucleotide-resolution CTSS profiles
#####################################

# Tag count of every sample at every nucleotide of every interval, on the interval's strand;
# positions without a CTSS are 0
nucleotide_profiles <- function(intervals, ctss) {
  base <- intervals[, .(pos = seq.int(start, end)), by = .(promoter_id, chr, strand)]
  rbindlist(lapply(names(ctss), function(s) {
    x <- copy(base)
    x[ctss[[s]], on = .(chr, strand, pos), count := i.count]
    x[is.na(count), count := 0L]
    x[, sample := s]
  }))
}

#####################################
# Steps 4-8. Shape features
#####################################

# Hoskins et al. 2011 shape index: SI = 2 + sum over positions with >= 1 tag of p_i log2 p_i
shape_index <- function(count) {
  n <- sum(count)
  if (n == 0) return(NA_real_)
  p <- count[count > 0] / n
  2 + sum(p * log2(p))
}

# 10th-90th percentile initiation width: positions in transcription order, q10 / q90 = the first
# position at which the cumulative tag count reaches 10% / 90% of the total; width = |q90 - q10|
# bp (0 when both quantiles fall on one nucleotide)
tss_width_10_90 <- function(pos, count, strand) {
  o <- order(pos * strand_sign(strand))
  pos <- pos[o]; cum <- cumsum(count[o]); n <- cum[length(cum)]
  if (n == 0) return(NA_integer_)
  as.integer(abs(pos[which(10 * cum >= 9 * n)[1L]] - pos[which(10 * cum >= n)[1L]]))
}

# Features of one aggregated promoter profile (all NA except the counts when there are no tags)
#
# 1. strongest_CAGE_TSS: the highest-count nucleotide; ties go to the tied position closest to
#    dominant_TSS, then the 5'-most one (n_strongest_tied counts the tied positions)
# 2. dominant_TSS_fraction: share of the tags at strongest_CAGE_TSS
# 3. dominant_TSS_shift: strongest_CAGE_TSS - dominant_TSS, sign flipped on the minus strand,
#    so positive = strongest_CAGE_TSS downstream of dominant_TSS
# 4. dominant_TSS_tags: tags at dominant_TSS when it lies in the promoter, else NA
shape_features <- function(pos, count, strand, dominant_tss) {
  total <- sum(count)
  at_dominant <- count[pos == dominant_tss]
  out <- list(total_CAGE_tags = total, n_active_CTSS = sum(count > 0L),
              dominant_TSS_tags = if (length(at_dominant)) at_dominant else NA_integer_,
              shape_index = NA_real_, strongest_CAGE_TSS = NA_integer_, n_strongest_tied = NA_integer_,
              dominant_TSS_shift = NA_integer_, dominant_TSS_fraction = NA_real_, TSS_width_10_90 = NA_integer_)
  if (total == 0L) return(out)
  top <- pos[count == max(count)]
  near <- top[abs(top - dominant_tss) == min(abs(top - dominant_tss))]
  strongest <- if (strand == "+") min(near) else max(near)
  out$shape_index <- shape_index(count)
  out$strongest_CAGE_TSS <- as.integer(strongest)
  out$n_strongest_tied <- length(top)
  out$dominant_TSS_shift <- as.integer((strongest - dominant_tss) * strand_sign(strand))
  out$dominant_TSS_fraction <- max(count) / total
  out$TSS_width_10_90 <- tss_width_10_90(pos, count, strand)
  out
}

# Deterministic per-promoter seed, independent of which other promoters are analysed
promoter_seed <- function(promoter_id) SEED + as.integer(sub("^hg_([0-9]+)\\..*$", "\\1", promoter_id))

# Bootstrap of the observed initiation distribution: N_RESAMPLE multinomial draws of the observed
# tag total from the observed per-nucleotide proportions, SI of each draw, and the fraction of
# draws on the same side of SI_CUT as the observed SI
resample_stability <- function(count, seed, n_resample = N_RESAMPLE) {
  count <- count[count > 0]
  si <- shape_index(count)
  set.seed(seed)
  p <- rmultinom(n_resample, sum(count), count) / sum(count)
  plogp <- p * log2(p)
  plogp[p == 0] <- 0
  resampled <- 2 + colSums(plogp)
  list(classification_stability = mean((resampled > SI_CUT) == (si > SI_CUT)))
}

#####################################
# Pipeline
#####################################
#
# 1. promoters: the FANTOM5 hg38 peak containing each gene's dominant_TSS (kept even without
#    tags), plus the other same-strand peaks FANTOM5 assigns to the gene that have >= 1 tag in
#    every GM12878 replicate
# 2. tag count of each replicate at every nucleotide of every promoter, on its strand
# 3. replicate depth check; aggregated profile = sum of the raw replicate counts
# 4-8. strongest_CAGE_TSS, shape index, Hoskins / Schor classes, bootstrap stability, features
#
# Outputs (LCL_project/promoter_shape):
#   promoter_shape_annotations.tsv               one row per gene x promoter
#   CTSS_profiles/sample_CTSS_profiles.tsv.gz    every nucleotide x replicate: count, CPM
#   CTSS_profiles/aggregated_CTSS_profiles.tsv.gz every nucleotide: summed count, fraction
run_promoter_shape <- function() {
  for (d in c("CTSS_profiles", "figures")) dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)
  ids <- names(SAMPLES)

  #####################################
  # Step 1. Promoter regions
  #####################################
  genes <- example_genes()
  fantom5 <- fantom5_peaks()
  cand <- promoter_candidates(genes, fantom5$peaks, fantom5$links)

  #####################################
  # Step 2. Per-nucleotide CTSS counts of every GM12878 replicate
  #####################################
  ctss <- setNames(lapply(ids, read_ctss), ids)
  depth <- data.table(sample = ids, replicate = unname(SAMPLES),
                      library_tags = vapply(ctss, function(x) sum(as.numeric(x$count)), numeric(1)),
                      ctss_positions = vapply(ctss, nrow, integer(1)))

  profiles <- nucleotide_profiles(unique(cand[same_strand == TRUE, .(promoter_id, chr, start, end, strand)]), ctss)
  tags <- dcast(profiles[, .(tags = sum(count)), by = .(promoter_id, sample)], promoter_id ~ sample, value.var = "tags")
  setnames(tags, ids, paste0("tags_", ids))
  cand <- merge(cand, tags, by = "promoter_id", all.x = TRUE)
  cand[, detected_in_all_samples := same_strand & rowSums(.SD > 0) == length(ids), .SDcols = paste0("tags_", ids)]
  cand[, retained := same_strand & (contains_dominant_TSS | detected_in_all_samples)]
  promoters <- cand[retained == TRUE, .(gene_id, gene_name, gene_group, promoter_id, chr, start, end, strand,
                                        dominant_TSS, dominant_TSS_source, contains_dominant_TSS)]

  #####################################
  # Step 3. Replicate depth and aggregation
  #####################################
  # FANTOM5 CTSS keep only reads with mapping quality above 20 (single-mapped) and >= 85% identity,
  # so no further read filtering is possible or done here
  depth_ratio <- max(depth$library_tags) / min(depth$library_tags)
  if (depth_ratio > MAX_DEPTH_RATIO)
    stop(sprintf("Replicate library sizes differ %.2f-fold (> %.1f): normalise before summing", depth_ratio, MAX_DEPTH_RATIO))

  sample_profiles <- profiles[promoter_id %in% promoters$promoter_id]
  sample_profiles[depth, on = "sample", cpm := count / i.library_tags * 1e6]
  aggregated <- sample_profiles[, .(count = sum(count)), by = .(promoter_id, chr, strand, pos)]
  aggregated <- merge(promoters[, .(gene_id, gene_name, promoter_id, dominant_TSS)], aggregated,
                      by = "promoter_id", allow.cartesian = TRUE)
  aggregated[, `:=`(offset_from_dominant_TSS = (pos - dominant_TSS) * strand_sign(strand),
                    fraction = if (sum(count) > 0L) count / sum(count) else NA_real_), by = .(gene_id, promoter_id)]

  #####################################
  # Steps 4-8. Shape features, classes and bootstrap stability
  #####################################
  features <- aggregated[, shape_features(pos, count, strand[1L], dominant_TSS[1L]), by = .(gene_id, promoter_id)]
  stability <- aggregated[, if (sum(count) >= MIN_TAGS) resample_stability(count, promoter_seed(promoter_id[1L])),
                          by = .(gene_id, promoter_id)]
  res <- merge(promoters, features, by = c("gene_id", "promoter_id"))
  res <- merge(res, stability, by = c("gene_id", "promoter_id"), all.x = TRUE)
  res[, QC_status := fcase(total_CAGE_tags == 0L, "no_CAGE_tags",
                           total_CAGE_tags < MIN_TAGS, "low_coverage",
                           classification_stability < STABILITY_MIN, "unstable_classification",
                           default = "pass")]
  res[, `:=`(Hoskins_shape_class = fifelse(shape_index > SI_CUT, "peaked", "broad"),
             Schor_shape_class = fifelse(shape_index > SI_CUT, "narrow", "broad"))]
  res[QC_status != "pass", c("Hoskins_shape_class", "Schor_shape_class") := "unclassified"]

  #####################################
  # Step 10. Outputs
  #####################################
  res[, gene_group := factor(gene_group, levels = GROUPS)]
  res[, gene_order := match(gene_id, genes$gene_id)]
  setorder(res, gene_group, gene_order, -contains_dominant_TSS, start)
  res[, gene_group := as.character(gene_group)]
  annotation_cols <- c("gene_id", "gene_name", "gene_group", "promoter_id", "chr", "start", "end", "strand",
    "dominant_TSS", "dominant_TSS_source", "contains_dominant_TSS", "strongest_CAGE_TSS", "dominant_TSS_shift",
    "total_CAGE_tags", "n_active_CTSS", "shape_index", "dominant_TSS_fraction", "TSS_width_10_90",
    "Hoskins_shape_class", "Schor_shape_class", "classification_stability", "QC_status")
  fwrite(res[, ..annotation_cols], file.path(OUT_DIR, "promoter_shape_annotations.tsv"), sep = "\t", na = "NA", quote = FALSE)

  key <- res[, .(gene_id, promoter_id, row = .I)]
  sample_out <- merge(promoters[, .(gene_id, gene_name, promoter_id, dominant_TSS)], sample_profiles,
                      by = "promoter_id", allow.cartesian = TRUE)
  sample_out <- merge(sample_out, key, by = c("gene_id", "promoter_id"))
  sample_out[, `:=`(offset_from_dominant_TSS = (pos - dominant_TSS) * strand_sign(strand),
                    replicate = unname(SAMPLES[sample]))]
  setorder(sample_out, row, pos, sample)
  fwrite(sample_out[, .(gene_id, gene_name, promoter_id, chr, pos, strand, offset_from_dominant_TSS,
                        sample, replicate, count, cpm)],
         file.path(OUT_DIR, "CTSS_profiles/sample_CTSS_profiles.tsv.gz"), sep = "\t")
  aggregated <- merge(aggregated, key, by = c("gene_id", "promoter_id"))
  setorder(aggregated, row, pos)
  fwrite(aggregated[, .(gene_id, gene_name, promoter_id, chr, pos, strand, offset_from_dominant_TSS, count, fraction)],
         file.path(OUT_DIR, "CTSS_profiles/aggregated_CTSS_profiles.tsv.gz"), sep = "\t", na = "NA", quote = FALSE)

  message(nrow(res), " promoters for ", uniqueN(res$gene_id), " of ", nrow(genes), " genes; QC: ",
          paste(names(table(res$QC_status)), table(res$QC_status), sep = "=", collapse = ", "))
  message("Wrote ", file.path(OUT_DIR, "promoter_shape_annotations.tsv"))
}

if (sys.nframe() == 0L) run_promoter_shape()
