#' Processing for nucleosome positioning around CAGE TSS and nucleosome-to-
#' nucleosome spacing over DHSs.
#'
#' FiberHMM output:
#'   * the BED12 files carry 13 columns (last col is an all-zero block field) and
#'     are read as 12-column BED12;
#'   * FiberHMM has NO sentinel first/last blocks, so every block is a real
#'     footprint (drop_terminal_blocks = FALSE).
#'
#' Features profiled: FiberHMM footprint midpoints (per read), strand-aware
#' relative to CAGE TSS. Produces (a) binned metaprofiles per timepoint and
#' (b) nearest upstream/downstream footprint per (TSS, read) for the shift
#' analysis.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
FP_ROOT   <- "/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_footprint"
# Dominant CTSS per gene (one anchor per gene). Per-CTSS windows are far too
# dense (median gap ~1bp) for tractable genome-wide metaprofiles.
CAGE_PATH <- "/project/spott/cshan/fiber-seq/macrophage_project/annotations/fantom5_thp1/THP1.consensus.CAGE_peaks.dominantPerGene.bed.gz"
# All DHS summits genome-wide (ENCFF359IAQ), anchors for the spacing analysis.
DHS_PATH  <- "/project/spott/cshan/fiber-seq/macrophage_project/nucleosome_positioning/results/annotations/ENCFF359IAQ.summits.bed.gz"
OUT_DIR   <- "/project/spott/cshan/fiber-seq/macrophage_project/nucleosome_positioning/results/tables"

# LPS timepoints: directory name -> label (0 = control)
TIMEPOINTS <- c(LPS0 = "0", LPS5 = "5", LPS10 = "10", LPS15 = "15")

WINDOW_BP <- 2000L
BIN_SIZE  <- 10L  # TSS-offset bin for metaprofiles
N_RANKS   <- 3L   # nearest N footprints per direction per (TSS, read)

# Spacing analysis: DHS summit +/- window, pairwise-distance histogram.
# SPACING_BIN_SIZE bins a distance axis, not a TSS offset -- kept separate from
# BIN_SIZE so the two can diverge.
SPACING_WINDOW_BP <- 2000L
SPACING_BIN_SIZE  <- 10L
SPACING_MAX_DIST  <- 1500L   # cap pairwise distance (~5 nucleosomes)

# Restrict to mono-nucleosome-sized footprints [FP_SIZE_MIN, FP_SIZE_MAX).
# Applied at extraction, so every downstream profile (aggregated, per-read,
# nearest-updown) inherits the filter.
FP_SIZE_MIN <- 130L
FP_SIZE_MAX <- 160L

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# IO helpers
# ---------------------------------------------------------------------------

# FiberHMM per-chr footprint BED (13 cols); read the standard 12 BED12 cols.
read_fiberhmm_bed12 <- function(path) {
  fread(
    path, header = FALSE, select = 1:12,
    col.names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
                  "thickStart", "thickEnd", "itemRgb",
                  "blockCount", "blockSizes", "blockStarts")
  )
}

list_sample_beds <- function(sample_dir) {
  sort(list.files(
    sample_dir,
    pattern = "_hmm_extracted_footprint_chr.*\\.bed\\.gz$",
    full.names = TRUE
  ))
}

# CAGE TSS: BED9 with gene_id in col4, CAGE strand in col6.
load_cage_tss <- function(path) {
  cols <- c("chrom", "start", "end", "gene_id", "score", "strand",
            "thickStart", "thickEnd", "itemRgb")
  dt <- fread(path, header = FALSE, col.names = cols)
  GRanges(
    seqnames = dt$chrom,
    ranges   = IRanges(dt$start + 1L, dt$end),  # BED 0-based, half-open
    strand   = dt$strand,
    gene_id  = dt$gene_id,
    score    = dt$score
  )
}

# DHS summit +/- window, merged so overlapping windows count a region once.
load_dhs_windows <- function(path, window_bp = SPACING_WINDOW_BP) {
  cols <- c("chrom", "start", "end", "id", "score", "strand")
  dt <- fread(path, header = FALSE, col.names = cols)
  gr <- GRanges(seqnames = dt$chrom,
                ranges = IRanges(pmax(dt$end - window_bp, 1L),
                                 dt$end + window_bp))
  reduce(gr)  # merge overlapping windows
}

# ---------------------------------------------------------------------------
# Per-read footprint extraction 
# ---------------------------------------------------------------------------

# One row per footprint block, preserving the read id (RID).
# Vectorized: split every block list at once and expand read metadata with rep(),
# which is far faster and lighter than a per-read `by` group over millions of reads.

extract_footprints_per_read <- function(bed12) {
  # result is a list, because every read may contain a different number of blocks
  sizes_l  <- strsplit(bed12$blockSizes, ",", fixed = TRUE) # fixed = Treats the comma as an ordinary literal character rather than a regular expression
  starts_l <- strsplit(bed12$blockStarts, ",", fixed = TRUE)
  
  nblocks  <- lengths(sizes_l)

  # combine all list entries into one long vector
  # size and starts are flattend in the same read order
  sizes  <- as.integer(unlist(sizes_l,  use.names = FALSE))
  starts <- as.integer(unlist(starts_l, use.names = FALSE))
  
  # repeat each read's chr according to the number of blocks in that read
  chrom  <- rep(bed12$chrom,      nblocks)
  cstart <- rep(bed12$chromStart, nblocks)
  rid    <- rep(bed12$name,       nblocks)

  # construct one row per footprint and calculate midpoints
  dt <- data.table(chrom = chrom, RID = rid,
                   fp_size = sizes,
                   fp_mid = as.integer(cstart + starts + sizes / 2))
  # Keep only mono-nucleosome-sized footprints (filter at the single choke point)
  dt[!is.na(fp_size) & !is.na(fp_mid) &
       fp_size >= FP_SIZE_MIN & fp_size < FP_SIZE_MAX]
}

# Strand-aware per-read footprint offsets to TSS, ranked per direction.
# Returns one row per (TSS, read, footprint) within the window.
#
# Scalable design: overlap footprint MIDPOINTS (each carrying its RID) directly
# with the TSS +/- window via findOverlaps. This is an interval join, a footprint is counted
# for a TSS only when its midpoint lands inside that TSS window - so the per-read structure (RID) is preserved without
# joining on read spans.
find_footprints_near_tss_per_read <- function(bed12, tss_gr,
                                              window_bp = WINDOW_BP,
                                              n_ranks = N_RANKS) {
  
  # creates a +/- 2kb window near CAGE TSS
  tss_windows <- promoters(tss_gr, upstream = window_bp, downstream = window_bp + 1L)

  # Pre-filter to reads whose span overlaps a TSS window before parsing blocks.
  # Only these reads can contribute footprints near a TSS, and this shrinks the
  # per-chromosome block table ~10x (peak memory), since most reads are far from
  # any anchor
  
  # Convert fiberHMM bed12 to Granges, each row is one footprint interval
  reads_gr <- GRanges(
    seqnames = bed12$chrom,
    ranges   = IRanges(bed12$chromStart + 1L, bed12$chromEnd)
  )
  
  # queryHits(...): extracts quary indices - read row numbers
  ## Keeps only the BED12 rows whose indices appear in keep
  keep <- unique(queryHits(findOverlaps(reads_gr, tss_windows)))
  if (!length(keep)) return(data.table())
  bed12 <- bed12[keep]

  # extract individual footprints
  fp_read <- extract_footprints_per_read(bed12)
  if (!nrow(fp_read)) return(data.table())

  # convert footprints to GRanges
  fp_gr <- GRanges(
    seqnames = fp_read$chrom,
    ranges   = IRanges(fp_read$fp_mid + 1L, fp_read$fp_mid + 1L),
    strand   = "*"
  )

  # finds every footprint midpoint that lies inside a window
  hits <- findOverlaps(fp_gr, tss_windows)
  if (!length(hits)) return(data.table())

  # query = footprint midpoint
  # subject = TSS window
  qi <- queryHits(hits) # extract footprint indices
  si <- subjectHits(hits) # corresponding TSS indices

  # creates a table with one row per footprint-TSS overlap
  rtp <- data.table(
    tss_idx    = si,
    RID        = fp_read$RID[qi],
    fp_mid     = fp_read$fp_mid[qi],
    fp_size    = fp_read$fp_size[qi],
    tss_pos    = start(tss_gr)[si],
    tss_strand = as.character(strand(tss_gr))[si],
    gene_id    = as.character(mcols(tss_gr)$gene_id)[si]
  )

  # add new column with relative position to TSS
  rtp[, rel_pos := fifelse(tss_strand == "+", fp_mid - tss_pos, tss_pos - fp_mid)]
  rtp[, abs_dist := abs(rel_pos)]
  
  # negative positive = upstream
  rtp[, direction := fifelse(rel_pos < 0L, "upstream", "downstream")]
  setorder(rtp, tss_idx, RID, direction, abs_dist)
  
  # for each TSS, read, direction combination, sort by ascending order
  rtp[, rank := seq_len(.N), by = .(tss_idx, RID, direction)]
  # assign ranking in each group
  rtp <- rtp[rank <= n_ranks]
  # keep only the nearest 3 footprints in each direction
  rtp[, signed_rank := fifelse(direction == "upstream", -rank, rank)]
  rtp[]
}

# Aggregated (non-per-read) strand-aware footprint offsets to TSS.
# Every footprint midpoint in a TSS window is counted once -- no read grouping,
# no ranking Returns one row per (footprint midpoint, TSS) overlap.
find_footprints_near_tss_aggregated <- function(bed12, tss_gr,
                                                window_bp = WINDOW_BP) {
  
  # create windows around each TSS
  tss_windows <- promoters(tss_gr, upstream = window_bp, downstream = window_bp + 1L)

  # Same read pre-filter as the per-read path, then overlap all
  # surviving footprint midpoints against the TSS windows
  
  # convert each full read into a genomic range
  reads_gr <- GRanges(
    seqnames = bed12$chrom,
    ranges   = IRanges(bed12$chromStart + 1L, bed12$chromEnd)
  )
  keep <- unique(queryHits(findOverlaps(reads_gr, tss_windows)))
  if (!length(keep)) return(data.table())
  fp <- extract_footprints_per_read(bed12[keep])
  if (!nrow(fp)) return(data.table())

  fp_gr <- GRanges(
    seqnames = fp$chrom,
    ranges   = IRanges(fp$fp_mid + 1L, fp$fp_mid + 1L),
    strand   = "*"
  )
  hits <- findOverlaps(fp_gr, tss_windows)
  if (!length(hits)) return(data.table())

  qi <- queryHits(hits)
  si <- subjectHits(hits)
  dt <- data.table(
    tss_idx    = si,
    fp_mid     = fp$fp_mid[qi],
    fp_size    = fp$fp_size[qi],
    tss_pos    = start(tss_gr)[si],
    tss_strand = as.character(strand(tss_gr))[si],
    gene_id    = as.character(mcols(tss_gr)$gene_id)[si]
  )
  dt[, rel_pos := fifelse(tss_strand == "+", fp_mid - tss_pos, tss_pos - fp_mid)]
  dt[]
}

# ---------------------------------------------------------------------------
# Aggregation across chromosomes for one sample
# ---------------------------------------------------------------------------

# Bin footprint offsets into a zero-filled metaprofile.
# normalize:
#   "none"     raw footprint-midpoint counts per bin
#   "tss"      per COVERED TSS -- denom = n_covered_tss from process_sample() =
#              number of TSS whose +/- window overlaps >= 1 read (annotated TSS
#              with no overlapping read are excluded)
#   "read_tss" per (read, TSS) pair that carries a footprint in the window --
#              denominator = unique (tss_idx, RID) in `hits` (self-computed)
bin_metaprofile <- function(hits, denom = NULL, window_bp = WINDOW_BP,
                            bin_size = BIN_SIZE, normalize = "none") {

  # creates a table containing every bin position from upstream to downstream of the TSS
  bins <- data.table(bin = seq(-window_bp, window_bp, by = bin_size))
  if (!nrow(hits)) {
    # adds a column named N and assigns zero to every bin
    bins[, N := 0]
    return(bins[])
  }

  # Assign each footprint to a bin
  hits[, bin := round(rel_pos / bin_size) * bin_size]

  # count footprint midpoints in each observed bin
  meta <- hits[, .N, by = bin]

  # add zero count binds
  meta <- merge(bins, meta, by = "bin", all.x = TRUE, sort = TRUE)
  # replace missing counts with 0
  meta[is.na(N), N := 0]

  if (normalize == "tss") {
    # per covered TSS (denominator passed in from process_sample coverage)
    meta[, N := if (is.null(denom) || denom == 0) 0 else N / denom]
  } else if (normalize == "read_tss") {
    # per (read, TSS) pair that carries a footprint in the window
    total_read_tss_pairs <- nrow(unique(hits[, .(tss_idx, RID)]))
    meta[, N := if (total_read_tss_pairs == 0) 0 else N / total_read_tss_pairs]
  }
  meta[order(bin)]
}

# Per-gene metaprofiles with a mode (aggregated | per_read) x norm layout.
#   aggregated: raw, per_tss (divide by n TSS anchors for that gene)
#   per_read:   raw, per_read (divide by that gene's read-TSS pair count)
# Returns long: gene, bin, N, mode, norm.
selected_gene_metaprofiles <- function(sel_hits, sel_agg, sel_genes,
                                        window_bp = WINDOW_BP, bin_size = BIN_SIZE) {
  bins <- CJ(gene = sel_genes,
             bin  = seq(-window_bp, window_bp, by = bin_size))
  out <- list()

  # aggregated
  if (nrow(sel_agg)) {
    a <- copy(sel_agg)
    setnames(a, "gene_id", "gene")
    a[, bin := round(rel_pos / bin_size) * bin_size]
    n_tss <- a[, .(n_tss = uniqueN(tss_idx)), by = gene]
    ma <- a[, .N, by = .(gene, bin)]
    ma <- merge(bins, ma, by = c("gene", "bin"), all.x = TRUE)
    ma[is.na(N), N := 0]
    ma <- merge(ma, n_tss, by = "gene", all.x = TRUE)
    out$agg_raw <- ma[, .(gene, bin, N, mode = "aggregated", norm = "raw")]
    out$agg_tss <- ma[, .(gene, bin,
                          N = fifelse(is.na(n_tss) | n_tss == 0, 0, N / n_tss),
                          mode = "aggregated", norm = "per_tss")]
  }

  # per-read
  if (nrow(sel_hits)) {
    h <- copy(sel_hits)
    setnames(h, "gene_id", "gene")
    h[, bin := round(rel_pos / bin_size) * bin_size]
    denom <- h[, .(pairs = uniqueN(paste(tss_idx, RID))), by = gene]
    mh <- h[, .N, by = .(gene, bin)]
    mh <- merge(bins, mh, by = c("gene", "bin"), all.x = TRUE)
    mh[is.na(N), N := 0]
    mh <- merge(mh, denom, by = "gene", all.x = TRUE)
    out$pr_raw <- mh[, .(gene, bin, N, mode = "per_read", norm = "raw")]
    out$pr_pr  <- mh[, .(gene, bin,
                         N = fifelse(is.na(pairs) | pairs == 0, 0, N / pairs),
                         mode = "per_read", norm = "per_read")]
  }

  rbindlist(out, use.names = TRUE, fill = TRUE)
}

# Collect nearest up/down footprint (rank 1) per (TSS, read) for the shift and
# size analyses (fp_size = footprint block width).
nearest_updown <- function(hits) {
  if (!nrow(hits)) return(data.table())
  hits[rank == 1L,
       .(tss_idx, RID, gene_id, direction, rel_pos, abs_dist, fp_size)]
}

# Process one timepoint: read the per-chr beds, tally hits against the TSS set.
# Only chromosomes that carry at least one TSS are read (skips extracting the
# many blocks on chromosomes with no anchors -- important for the small
# selected-gene anchor set).
process_sample <- function(sample_dir, tss_gr,
                           window_bp = WINDOW_BP, n_ranks = N_RANKS,
                           mode = c("per_read", "aggregated")) {
  mode <- match.arg(mode)
  bed_files <- list_sample_beds(sample_dir)
  if (!length(bed_files)) stop("No footprint BEDs in ", sample_dir)

  tss_chroms <- unique(as.character(seqnames(tss_gr)))
  file_chrom <- sub(".*_footprint_(chr[^.]+)\\.bed\\.gz$", "\\1", basename(bed_files))
  bed_files  <- bed_files[file_chrom %in% tss_chroms]
  if (!length(bed_files)) {
    empty <- data.table()
    setattr(empty, "n_covered_tss", 0)
    return(empty)
  }

  # +/- window per TSS. Used to count read coverage for the "tss" normalization:
  # a TSS is "covered" if its window overlaps >= 1 read span, independent of
  # whether the read carries an extracted footprint. TSS indices are disjoint
  # across chromosomes, so per-chrom covered-TSS counts sum to the sample total.
  tss_windows <- promoters(tss_gr, upstream = window_bp, downstream = window_bp + 1L)

  res_list <- lapply(bed_files, function(f) {
    bed12 <- read_fiberhmm_bed12(f)
    reads_gr <- GRanges(bed12$chrom,
                        IRanges(bed12$chromStart + 1L, bed12$chromEnd))
    ov <- findOverlaps(reads_gr, tss_windows)
    hits <- if (mode == "aggregated")
      find_footprints_near_tss_aggregated(bed12, tss_gr, window_bp)
    else
      find_footprints_near_tss_per_read(bed12, tss_gr, window_bp, n_ranks)
    list(hits  = hits,
         n_tss = as.numeric(uniqueN(subjectHits(ov))))
  })

  res <- rbindlist(lapply(res_list, `[[`, "hits"), use.names = TRUE, fill = TRUE)
  setattr(res, "n_covered_tss", sum(vapply(res_list, function(x) x$n_tss, numeric(1))))
  res
}

# ---------------------------------------------------------------------------
# Nucleosome-to-nucleosome spacing (point-process autocorrelation) over DHSs
# ---------------------------------------------------------------------------
#
# For each timepoint, pool all pairwise distances between mono-nucleosome
# footprint midpoints that fall inside a DHS window, on the same read. The
# histogram of these distances is the pair-correlation of nucleosome positions:
# peaks at the nucleosome repeat length and its multiples (2x, 3x...) indicate
# phased arrays. Anchors are ALL DHS regions genome-wide -- no TSS restriction.

# For one chromosome BED, return a binned distance histogram (raw counts) plus
# the number of reads contributing (for normalization). Distances are computed
# per (read, DHS window): only midpoints on the same read AND inside the same
# merged window are paired, so pairs are counted once per genomic region.
spacing_hist_one_chrom <- function(bed12, dhs_win,
                                   bin_size = SPACING_BIN_SIZE,
                                   max_dist = SPACING_MAX_DIST) {
  empty <- list(hist = data.table(dist_bin = integer(), N = numeric()),
                n_reads = 0L)

  # Pre-filter to reads whose span overlaps a DHS window (memory).
  reads_gr <- GRanges(bed12$chrom,
                      IRanges(bed12$chromStart + 1L, bed12$chromEnd))
  keep <- unique(queryHits(findOverlaps(reads_gr, dhs_win)))
  if (!length(keep)) return(empty)
  fp <- extract_footprints_per_read(bed12[keep])[, .(chrom, RID, fp_mid)]
  if (!nrow(fp)) return(empty)

  # Assign each midpoint to the DHS window it falls in (1 bp point overlap).
  fp_gr <- GRanges(fp$chrom, IRanges(fp$fp_mid + 1L, fp$fp_mid + 1L))
  ov <- findOverlaps(fp_gr, dhs_win)
  if (!length(ov)) return(empty)
  m <- fp[queryHits(ov)]
  m[, win_id := subjectHits(ov)]

  # Pairwise distances within each (window, read) group of >= 2 midpoints.
  m <- m[, if (.N >= 2L) .(fp_mid = fp_mid) else NULL, by = .(win_id, RID)]
  if (!nrow(m)) return(empty)

  pair_dists <- m[, {
    if (.N < 2L) NULL else {
      p <- combn(sort(fp_mid), 2L)
      list(d = p[2L, ] - p[1L, ])
    }
  }, by = .(win_id, RID)]

  pair_dists <- pair_dists[d > 0L & d <= max_dist]
  if (!nrow(pair_dists)) return(empty)

  pair_dists[, dist_bin := round(d / bin_size) * bin_size]
  hist <- pair_dists[, .(N = .N), by = dist_bin]
  list(hist = hist, n_reads = uniqueN(m$RID))
}

# Process one timepoint across all chromosomes -> binned distance histogram.
spacing_hist_sample <- function(sample_dir, dhs_win_by_chrom,
                                bin_size = SPACING_BIN_SIZE,
                                max_dist = SPACING_MAX_DIST) {
  bed_files <- list_sample_beds(sample_dir)
  if (!length(bed_files)) stop("No footprint BEDs in ", sample_dir)
  file_chrom <- sub(".*_footprint_(chr[^.]+)\\.bed\\.gz$", "\\1", basename(bed_files))

  hist_list <- list(); total_reads <- 0L
  for (i in seq_along(bed_files)) {
    ch <- file_chrom[i]
    dhs_win <- dhs_win_by_chrom[[ch]]
    if (is.null(dhs_win) || !length(dhs_win)) next
    res <- spacing_hist_one_chrom(read_fiberhmm_bed12(bed_files[i]), dhs_win,
                                  bin_size, max_dist)
    if (nrow(res$hist)) hist_list[[ch]] <- res$hist
    total_reads <- total_reads + res$n_reads
  }
  if (!length(hist_list)) {
    return(list(hist = data.table(dist_bin = integer(), N = numeric()),
                n_reads = 0L))
  }
  hist <- rbindlist(hist_list)[, .(N = sum(N)), by = dist_bin][order(dist_bin)]
  list(hist = hist, n_reads = total_reads)
}

# ---------------------------------------------------------------------------
# Footprint size distribution (UNFILTERED by size)
# ---------------------------------------------------------------------------
#
# Every other path in this file keeps the mono-nucleosome filter
# [FP_SIZE_MIN, FP_SIZE_MAX) applied inside extract_footprints_per_read(). This
# section deliberately does NOT filter by size -- the point is to show the whole
# distribution and where the 130-160bp nucleosome window sits inside it.
#
# Sizes are tallied per chromosome with tabulate() and accumulated into one count
# vector per sample, so no large per-footprint vector is ever held.

MAX_FP_SIZE_TALLY <- 20000L   # sizes above this are ignored in the tally

# Count footprint block sizes for one sample. Returns data.table(fp_size, N).
footprint_size_counts <- function(sample_dir, chroms = NULL) {
  bed_files <- list_sample_beds(sample_dir)
  if (!length(bed_files)) stop("No footprint BEDs in ", sample_dir)
  if (!is.null(chroms)) {
    fc <- sub(".*_footprint_(chr[^.]+)\\.bed\\.gz$", "\\1", basename(bed_files))
    bed_files <- bed_files[fc %in% chroms]
  }

  counts <- numeric(MAX_FP_SIZE_TALLY)
  for (f in bed_files) {
    dt <- fread(f, header = FALSE, select = 11, colClasses = "character",
                col.names = "blockSizes")
    sizes <- as.integer(unlist(strsplit(dt$blockSizes, ",", fixed = TRUE),
                               use.names = FALSE))
    rm(dt)
    sizes <- sizes[!is.na(sizes) & sizes >= 1L & sizes <= MAX_FP_SIZE_TALLY]
    if (length(sizes)) counts <- counts + tabulate(sizes, nbins = MAX_FP_SIZE_TALLY)
    rm(sizes); gc(verbose = FALSE)
  }
  data.table(fp_size = seq_len(MAX_FP_SIZE_TALLY), N = counts)[N > 0]
}

# Build the size table for every timepoint and save it. Plotting lives in
# plot_nuc_pos.R (plot_footprint_size_distribution()).
main_footprint_size_table <- function(chroms = NULL) {
  out <- list()
  for (i in seq_along(TIMEPOINTS)) {
    sample <- names(TIMEPOINTS)[i]; tp <- TIMEPOINTS[i]
    message("Tallying footprint sizes for ", sample, " (LPS ", tp, ") ...")
    d <- footprint_size_counts(file.path(FP_ROOT, sample), chroms = chroms)
    d[, `:=`(sample = sample, timepoint = tp)]
    message("  ", format(sum(d$N), big.mark = ","), " footprints; median ",
            d[, fp_size[which.max(cumsum(N) >= sum(N) / 2)]], " bp; ",
            sprintf("%.1f%% in %d-%dbp",
                    100 * d[fp_size >= FP_SIZE_MIN & fp_size < FP_SIZE_MAX, sum(N)] /
                      sum(d$N), FP_SIZE_MIN, FP_SIZE_MAX))
    out[[sample]] <- d
  }
  size_dt <- rbindlist(out)
  fwrite(size_dt, file.path(OUT_DIR, "footprint_size_distribution.tsv.gz"), sep = "\t")
  message("Wrote ", file.path(OUT_DIR, "footprint_size_distribution.tsv.gz"))
  invisible(size_dt)
}

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

# Selected promoter genes from the rmd. TSS subset = CAGE TSS whose gene_id
# base matches these genes' Ensembl ids. We resolve by gene SYMBOL via the
# rmd coordinates instead: subset CAGE TSS overlapping each promoter interval.
SELECTED_PROMOTERS <- data.table(
  gene = c("IL1B","TNF","IL1A","PTGS2","IL8","CXCL2","IL6","RELB","rs2836882"),
  chr  = c("chr2","chr6","chr2","chr1","chr4","chr4","chr7","chr19","chr21"),
  start = c(112836323L,31575265L,112793233L,186679985L,73740006L,74098863L,22726376L,45001011L,39094468L),
  end   = c(112837385L,31575824L,112794969L,186681163L,73741332L,74099651L,22727094L,45002417L,39094712L)
)

# gencode gene TSS (strand-aware single-bp point) used as a fallback anchor for
# selected genes with no THP-1 CAGE TSS in their promoter window.
GENCODE_TSS_PATH <- "/project/spott/cshan/annotations/TSS.gencode.v49.bed"

# Symbol -> gencode ENSG for the selected genes (IL8 = CXCL8). rs2836882 is a SNP
# and has no gencode gene TSS, so it is intentionally absent.
SELECTED_ENSG <- c(
  IL1B  = "ENSG00000125538", TNF   = "ENSG00000232810", IL1A = "ENSG00000115008",
  PTGS2 = "ENSG00000073756", IL8   = "ENSG00000169429", CXCL2 = "ENSG00000081041",
  IL6   = "ENSG00000136244", RELB  = "ENSG00000104856"
)

# Build the selected-gene TSS anchor set: one GRanges of TSS points tagged with
# gene symbol and source. Uses CAGE CTSS inside the promoter window where any
# exist, otherwise falls back to the gencode gene TSS.
build_selected_tss <- function(cage_gr, promoters_dt) {
  prom_gr <- GRanges(
    seqnames = promoters_dt$chr,
    ranges   = IRanges(promoters_dt$start + 1L, promoters_dt$end),
    gene     = promoters_dt$gene
  )
  ov <- findOverlaps(cage_gr, prom_gr)
  cage_sel <- cage_gr[queryHits(ov)]
  cage_gene <- prom_gr$gene[subjectHits(ov)]
  genes_with_cage <- unique(cage_gene)

  parts <- list()
  if (length(cage_sel)) {
    parts$cage <- GRanges(
      seqnames = seqnames(cage_sel), ranges = ranges(cage_sel),
      strand = strand(cage_sel),
      gene = cage_gene, source = "CAGE"
    )
  }

  # gencode fallback for selected genes lacking CAGE
  need_fallback <- setdiff(names(SELECTED_ENSG), genes_with_cage)
  if (length(need_fallback)) {
    gc <- fread(GENCODE_TSS_PATH, header = FALSE,
                col.names = c("chrom","start","end","gene_id","score","strand"))
    gc[, ensg := sub("\\..*$", "", gene_id)]
    want <- data.table(gene = need_fallback, ensg = SELECTED_ENSG[need_fallback])
    gc_sel <- merge(gc, want, by = "ensg")
    if (nrow(gc_sel)) {
      parts$gencode <- GRanges(
        seqnames = gc_sel$chrom,
        ranges   = IRanges(gc_sel$start + 1L, gc_sel$end),
        strand   = gc_sel$strand,
        gene     = gc_sel$gene, source = "gencode"
      )
    }
    dropped <- setdiff(need_fallback, if (nrow(gc_sel)) gc_sel$gene else character())
    if (length(dropped))
      message("  selected genes with neither CAGE nor gencode TSS: ",
              paste(dropped, collapse = ", "))
  }

  sel_gr <- do.call(c, unname(parts))
  message("  selected-gene anchors: ", length(sel_gr), " TSS across ",
          uniqueN(sel_gr$gene), " genes (",
          paste(sort(unique(sel_gr$gene)), collapse = ", "), ")")
  sel_gr
}

main_nuc_pos <- function() {
  message("Loading CAGE TSS: ", CAGE_PATH)
  tss_gr <- load_cage_tss(CAGE_PATH)
  tss_gr <- tss_gr[!is.na(mcols(tss_gr)$gene_id) & mcols(tss_gr)$gene_id != "NA"]
  message("  ", length(tss_gr), " CAGE TSS")

  # Selected-gene TSS anchors: CAGE where available, gencode TSS as fallback.
  # gene symbol is carried in gene_id so it propagates through the per-read finder.
  sel_gr <- build_selected_tss(tss_gr, SELECTED_PROMOTERS)
  mcols(sel_gr)$gene_id <- sel_gr$gene
  sel_genes <- sort(unique(sel_gr$gene))

  meta_all       <- list()  # genome-wide metaprofiles
  meta_selected  <- list()  # per selected-gene metaprofiles
  meta_sel_comb  <- list()  # selected genes pooled into one metaprofile
  nearest_all    <- list()  # nearest up/down summary (genome-wide)

  for (i in seq_along(TIMEPOINTS)) {
    sample <- names(TIMEPOINTS)[i]
    tp     <- TIMEPOINTS[i]
    message("Processing ", sample, " (LPS ", tp, ") ...")
    sample_dir <- file.path(FP_ROOT, sample)

    # --- Genome-wide ---
    # Per-read hits (ranked; used for nearest-updown + per-read metaprofiles).
    hits <- process_sample(sample_dir, tss_gr, mode = "per_read")
    # Covered-TSS denominator for the "tss" normalization (captured before the
    # in-place :=, which may drop attributes).
    n_cov_tss <- attr(hits, "n_covered_tss")
    hits[, `:=`(sample = sample, timepoint = tp)]
    # Aggregated hits (every footprint midpoint counted once).
    agg <- process_sample(sample_dir, tss_gr, mode = "aggregated")

    meta_all[[sample]] <- rbind(
      # aggregated: raw + per-covered-TSS
      bin_metaprofile(agg,  normalize = "none")[, `:=`(mode = "aggregated", norm = "raw")],
      bin_metaprofile(agg,  denom = n_cov_tss, normalize = "tss")[,  `:=`(mode = "aggregated", norm = "per_tss")],
      # per-read: raw + per-covered-TSS + per-footprint-read-TSS-pair
      bin_metaprofile(hits, normalize = "none")[,     `:=`(mode = "per_read", norm = "raw")],
      bin_metaprofile(hits, denom = n_cov_tss, normalize = "tss")[,  `:=`(mode = "per_read", norm = "per_tss")],
      bin_metaprofile(hits, normalize = "read_tss")[, `:=`(mode = "per_read", norm = "per_read")]
    )[, `:=`(sample = sample, timepoint = tp)]

    # --- Selected genes ---
    sel_hits <- process_sample(sample_dir, sel_gr, mode = "per_read")
    sel_cov_tss <- attr(sel_hits, "n_covered_tss")
    sel_agg  <- process_sample(sample_dir, sel_gr, mode = "aggregated")
    meta_selected[[sample]] <-
      selected_gene_metaprofiles(sel_hits, sel_agg, sel_genes)[, `:=`(sample = sample, timepoint = tp)]

    # Combined: pool all selected genes into one metaprofile (gains power lost by
    # splitting per gene). Identical treatment to the genome-wide path, but with
    # the selected-gene anchor set as the TSS reference.
    meta_sel_comb[[sample]] <- rbind(
      bin_metaprofile(sel_agg,  normalize = "none")[, `:=`(mode = "aggregated", norm = "raw")],
      bin_metaprofile(sel_agg,  denom = sel_cov_tss, normalize = "tss")[,  `:=`(mode = "aggregated", norm = "per_tss")],
      bin_metaprofile(sel_hits, normalize = "none")[,     `:=`(mode = "per_read", norm = "raw")],
      bin_metaprofile(sel_hits, denom = sel_cov_tss, normalize = "tss")[,  `:=`(mode = "per_read", norm = "per_tss")],
      bin_metaprofile(sel_hits, normalize = "read_tss")[, `:=`(mode = "per_read", norm = "per_read")]
    )[, `:=`(sample = sample, timepoint = tp)]

    # --- Nearest up/down (rank-1) footprints, genome-wide (per-read) ---
    nu <- nearest_updown(hits)
    if (nrow(nu)) nu[, `:=`(sample = sample, timepoint = tp)]
    nearest_all[[sample]] <- nu

    # Persist the raw per-read hit table per timepoint (can be large)
    fwrite(hits, file.path(OUT_DIR, sprintf("%s_footprint_tss_hits.tsv.gz", sample)),
           sep = "\t")
  }

  tp_levels <- unname(TIMEPOINTS)

  meta_all_dt <- rbindlist(meta_all, use.names = TRUE, fill = TRUE)
  meta_all_dt[, timepoint := factor(timepoint, levels = tp_levels)]
  fwrite(meta_all_dt, file.path(OUT_DIR, "metaprofile_genomewide_by_timepoint.tsv"), sep = "\t")

  meta_sel_dt <- rbindlist(meta_selected, use.names = TRUE, fill = TRUE)
  if (nrow(meta_sel_dt)) meta_sel_dt[, timepoint := factor(timepoint, levels = tp_levels)]
  fwrite(meta_sel_dt, file.path(OUT_DIR, "metaprofile_selected_genes_by_timepoint.tsv"), sep = "\t")

  meta_sel_comb_dt <- rbindlist(meta_sel_comb, use.names = TRUE, fill = TRUE)
  if (nrow(meta_sel_comb_dt)) meta_sel_comb_dt[, timepoint := factor(timepoint, levels = tp_levels)]
  fwrite(meta_sel_comb_dt, file.path(OUT_DIR, "metaprofile_combined_selected_by_timepoint.tsv"), sep = "\t")

  nearest_dt <- rbindlist(nearest_all, use.names = TRUE, fill = TRUE)
  nearest_dt[, timepoint := factor(timepoint, levels = tp_levels)]
  fwrite(nearest_dt, file.path(OUT_DIR, "nearest_updown_footprint_by_timepoint.tsv.gz"), sep = "\t")

  # Nearest up/down shift summary (median signed distance per timepoint/direction)
  shift_summary <- nearest_dt[, .(
    n            = .N,
    mean_rel_bp  = mean(rel_pos),
    median_rel_bp = as.numeric(median(rel_pos)),
    sd_rel_bp    = sd(rel_pos)
  ), by = .(timepoint, direction)][order(direction, timepoint)]
  fwrite(shift_summary, file.path(OUT_DIR, "nearest_updown_shift_summary.tsv"), sep = "\t")

  message("Done. Outputs in ", OUT_DIR)
  print(shift_summary)
}

main_nuc_spacing <- function() {
  message("Loading DHS windows: ", DHS_PATH)
  dhs_win <- load_dhs_windows(DHS_PATH)
  message("  ", length(dhs_win), " merged DHS windows (+/-", SPACING_WINDOW_BP, "bp)")
  dhs_by_chrom <- split(dhs_win, as.character(seqnames(dhs_win)))

  full_bins <- data.table(dist_bin = seq(0L, SPACING_MAX_DIST, by = SPACING_BIN_SIZE))
  out <- list()
  for (i in seq_along(TIMEPOINTS)) {
    sample <- names(TIMEPOINTS)[i]; tp <- TIMEPOINTS[i]
    message("Processing ", sample, " (LPS ", tp, ") ...")
    res <- spacing_hist_sample(file.path(FP_ROOT, sample), dhs_by_chrom)

    d <- merge(full_bins, res$hist, by = "dist_bin", all.x = TRUE)
    d[is.na(N), N := 0]
    total <- sum(d$N)
    d[, `:=`(
      N_norm   = if (total > 0) N / total else 0,   # fraction of all pairs
      sample   = sample,
      timepoint = tp,
      n_reads  = res$n_reads
    )]
    out[[sample]] <- d
    message("  ", sample, ": ", format(total, big.mark = ","),
            " pairs from ", format(res$n_reads, big.mark = ","), " reads")
  }

  spacing <- rbindlist(out)
  spacing[, timepoint := factor(timepoint, levels = unname(TIMEPOINTS))]
  fwrite(spacing, file.path(OUT_DIR, "nuc_spacing_over_DHS_by_timepoint.tsv"),
         sep = "\t")
  message("Done. Wrote ", file.path(OUT_DIR, "nuc_spacing_over_DHS_by_timepoint.tsv"))
}

# ---------------------------------------------------------------------------
# Nucleosome footprint clusters per topic-model read cluster
# ---------------------------------------------------------------------------
#
# The topic model (code/topic_model/1_m6a_promoter_topic_modelling.Rmd) k-means
# clusters the reads over each promoter by their m6A topic loadings, and
# plot_nuc_pos.R (plot_gene_nuc_topic) draws those read clusters side by side
# with their nucleosome footprints. This section turns that picture into
# numbers: inside each read cluster, find where footprints PILE UP along the
# promoter and report the midpoint of each pile, so that the same nucleosome can
# be compared across read clusters.
#
# Calling piles -- call_footprint_piles():
#   A pile is a window that is busier than usual FOR THAT READ CLUSTER. There is
#   no null model and no p-value: the line drawn on the track is a description of
#   footprints the reader can see underneath it, not a claim of significance.
#   * spread a read cluster's M footprints evenly over the region to get what a
#     PILE_WIN_BP window would hold on average: expected = M * win / region_len;
#   * slide that window across the region and keep any window holding at least
#     PILE_FOLD times its own cluster's expected count, AND at least PILE_MIN_OBS
#     footprints outright. The fold part is what makes clusters comparable -- a
#     flat count would just track how many footprints a cluster has (at these
#     settings pile count correlates with footprint count at rho 0.93 for a flat
#     count vs 0.77 here).
#     The absolute floor stops a sparse cluster gaming the ratio: one lone
#     footprint against an expectation of 0.3 is "3x enriched" but is not a pile;
#   * a real pile lights up a run of overlapping windows as the scan crosses it,
#     so calls are taken greedily, best fold first, keeping calls
#     PILE_MIN_SEP_BP apart -- one nucleosome's width, since two nucleosomes
#     cannot sit closer than that. One pile gives one call.
#
# Midpoints -- own_pile_midpoints():
#   each read cluster is scanned alone, and a pile's midpoint is the mean of the
#   footprint midpoints inside its winning window. No read cluster's numbers
#   depend on any other, so a midpoint describes only the band it is drawn in.
#   The midpoints are marked on the per-read track by plot_gene_nuc_topic() in
#   plot_nuc_pos.R.

TOPIC_DIR <- "/project/spott/cshan/fiber-seq/macrophage_project/topic_model"

# 50bp roughly matches the width of a real pile: a positioned nucleosome's
# midpoints scatter with SD ~25-30bp. Narrower windows find more piles (128 at
# 40bp vs 116 here vs 94 at 75bp) because the expected count scales with the
# window, so a tight pile is measured against a smaller bar -- but below ~40bp
# the window is narrower than the pile and PILE_MIN_OBS starts discarding real
# ones. Step barely matters (the min-separation rule collapses runs of windows
# anyway) and costs nothing, since the fold rule applies no multiplicity penalty.
PILE_WIN_BP     <- 50L   # sliding window
PILE_STEP_BP    <- 10L   # window step
PILE_FOLD       <- 2.5   # window must hold this many x the cluster's own average
PILE_MIN_OBS    <- 3L    # ...and at least this many footprints outright
PILE_MIN_SEP_BP <- 147L  # min spacing between piles (one nucleosome)

# Topic-model output directory name for a region (as written by the Rmd).
topic_outname <- function(gene, chrom, start, end)
  sprintf("%s_%s_%d_%d", gene, chrom, start, end)

# Per-read k-means cluster + alignment span from the topic-model TSV.
load_topic_read_clusters <- function(outname) {
  tsv <- file.path(TOPIC_DIR, outname, sprintf("read_topic_assignments_%s.tsv", outname))
  dt <- fread(tsv, colClasses = list(character = "RID"))
  dt[, .(RID, sample_name, cluster,
         cluster_num = as.integer(sub("cluster", "", cluster)),
         r_start = start, r_end = end)]
}

# Mono-nucleosome footprints for the topic reads over one region, read through
# the tabix index (the per-chromosome BEDs are far too big to fread whole).
# Same midpoint-in-region rule as extract_read_footprints() in plot_nuc_pos.R,
# so these are exactly the footprints drawn in the side-by-side figures.
read_topic_footprints <- function(topic_dt, chrom, view_start, view_end) {
  out <- list()
  for (sn in unique(topic_dt$sample_name)) {
    sdir <- gsub("_", "", sn)
    bed  <- file.path(FP_ROOT, sdir,
                      sprintf("%s_hmm_extracted_footprint_%s.bed.gz", sdir, chrom))
    if (!file.exists(bed)) { message("  [WARN] missing ", bed); next }
    txt <- tryCatch(
      Rsamtools::scanTabix(Rsamtools::TabixFile(bed),
                           param = GRanges(chrom, IRanges(view_start, view_end)))[[1]],
      error = function(e) character(0))
    if (!length(txt)) next
    rows <- fread(text = paste(txt, collapse = "\n"), header = FALSE,
                  select = c(2, 4, 11, 12), col.names = c("V2", "V4", "V11", "V12"),
                  colClasses = list(character = c("V4", "V11", "V12")))
    rows <- rows[V4 %in% topic_dt[sample_name == sn, RID]]
    if (!nrow(rows)) next

    sizes_l  <- strsplit(rows$V11, ",", fixed = TRUE)
    starts_l <- strsplit(rows$V12, ",", fixed = TRUE)
    nblk     <- lengths(sizes_l)
    sizes    <- as.integer(unlist(sizes_l, use.names = FALSE))
    fp_start <- rep(rows$V2, nblk) + as.integer(unlist(starts_l, use.names = FALSE))
    d <- data.table(RID = rep(rows$V4, nblk), fp_start = fp_start,
                    fp_end = fp_start + sizes, fp_size = sizes,
                    fp_mid = fp_start + sizes %/% 2L)
    out[[sn]] <- d[!is.na(fp_size) & fp_size >= FP_SIZE_MIN & fp_size < FP_SIZE_MAX &
                     fp_mid >= view_start & fp_mid <= view_end]
  }
  rbindlist(out)
}

# Sliding-window scan for footprint piles, relative to the read cluster's OWN
# average density. Returns one row per pile, empty data.table when none pass.
call_footprint_piles <- function(fp_mid, view_start, view_end,
                                 win = PILE_WIN_BP, step = PILE_STEP_BP,
                                 min_fold = PILE_FOLD, min_obs = PILE_MIN_OBS,
                                 min_sep = PILE_MIN_SEP_BP) {
  M <- length(fp_mid)
  if (M == 0L) return(data.table())

  half    <- win %/% 2L
  centers <- seq(view_start, view_end, by = step)
  wlo     <- pmax(centers - half, view_start)   # windows are clipped at the
  whi     <- pmin(centers + half, view_end)     # edges, and so is `expected`
  # what this cluster's own footprints would give if spread evenly
  expected <- M * (whi - wlo + 1L) / (view_end - view_start + 1L)

  srt <- sort(fp_mid)
  obs <- findInterval(whi, srt) - findInterval(wlo - 1L, srt)

  w <- data.table(peak_center = centers, win_start = wlo, win_end = whi,
                  obs = obs, expected = expected, fold = obs / expected)
  hit <- w[obs >= min_obs & fold >= min_fold][order(-fold, -obs, peak_center)]
  if (!nrow(hit)) return(data.table())

  # greedy: take the busiest window, then only windows >= min_sep from every pile
  # already kept, so one pile of footprints yields exactly one call.
  keep <- integer(0)
  for (k in seq_len(nrow(hit)))
    if (!length(keep) ||
        all(abs(hit$peak_center[k] - hit$peak_center[keep]) >= min_sep))
      keep <- c(keep, k)

  hit[keep][order(peak_center)]
}

# Per-read-cluster footprint piles and their midpoints, with NO pooling across
# read clusters: each cluster is scanned on its own footprints against its own
# average density, and each pile's midpoint is the mean of the footprint
# midpoints inside its winning window. Nothing from any other read cluster enters
# a cluster's numbers, so these midpoints are exactly what is visible in that
# cluster's band of the per-read track.
own_pile_midpoints <- function(gene_sym, chrom_id, view_start, view_end) {
  topic <- load_topic_read_clusters(topic_outname(gene_sym, chrom_id, view_start, view_end))
  fp    <- read_topic_footprints(topic, chrom_id, view_start, view_end)
  if (!nrow(fp)) { message("  no footprints in region"); return(data.table()) }
  fp <- merge(fp, topic[, .(RID, cluster)], by = "RID")

  out <- list()
  for (grp in sort(unique(topic$cluster))) {
    mids <- fp[cluster == grp, fp_mid]
    cc   <- call_footprint_piles(mids, view_start, view_end)
    if (!nrow(cc)) next

    # midpoint = mean of the footprints inside the pile's own window
    m <- cc[, {
      x <- mids[mids >= win_start & mids <= win_end]
      .(n_fp = length(x), mid_mean = mean(x), mid_median = as.numeric(median(x)),
        mid_sd = sd(x))
    }, by = .(peak_center, win_start, win_end, obs, expected, fold)]
    m[, read_cluster := grp]
    setorder(m, peak_center)
    m[, pile_id := seq_len(.N)]
    out[[grp]] <- m
  }
  res <- rbindlist(out, use.names = TRUE, fill = TRUE)
  if (!nrow(res)) return(res)
  res[, `:=`(gene = gene_sym, chrom = chrom_id,
             view_start = view_start, view_end = view_end,
             mid_se = mid_sd / sqrt(n_fp))]
  res[, .(gene, chrom, view_start, view_end, read_cluster, pile_id,
          peak_center, win_start, win_end, n_fp, obs, expected, fold,
          mid_mean, mid_median, mid_sd, mid_se)]
}

main_own_pile_midpoints <- function(regions = SELECTED_PROMOTERS) {
  out <- list()
  for (i in seq_len(nrow(regions))) {
    r <- regions[i]
    message("Own-cluster footprint piles for ", r$gene, " ...")
    d <- own_pile_midpoints(r$gene, r$chr, r$start, r$end)
    if (!nrow(d)) next
    out[[r$gene]] <- d
    message("  ", d[, .N], " piles across ", uniqueN(d$read_cluster), " read clusters")
  }
  res <- rbindlist(out, use.names = TRUE, fill = TRUE)
  fwrite(res, file.path(OUT_DIR, "nuc_own_pile_midpoints_by_topic_cluster.tsv"), sep = "\t")
  message("Done. Wrote nuc_own_pile_midpoints_by_topic_cluster.tsv to ", OUT_DIR)
  invisible(res)
}

# ---------------------------------------------------------------------------
# Distribution-level comparison of read clusters (position and occupancy)
# ---------------------------------------------------------------------------
#
# Comparing read clusters pile-by-pile needs a correspondence between piles that
# the data does not supply -- a cluster with one pile and a cluster with three
# cannot be lined up by rank without mismatching different nucleosomes. This
# section compares the WHOLE footprint distribution of each read cluster instead,
# so no matching is required.
#
# The profile is an INTENSITY, footprints per read per bp, not a probability
# density. Normalising to a density would divide occupancy away; keeping the
# per-read scale means one curve carries both questions -- its height is how
# often a nucleosome is there, its position is where.
#
# Two statistics per region, testing the two axes separately:
#   * POSITION -- 1-D Wasserstein-1 distance between each pair of read clusters,
#     the area between their ECDFs, in bp ("how far mass has to move to turn one
#     profile into the other"). Averaged over all pairs for one omnibus number.
#     Wasserstein compares shape only, so it is blind to occupancy by design.
#   * OCCUPANCY -- the fraction of a cluster's reads carrying any footprint in the
#     region, tested across clusters with Fisher's exact test on the 2 x k table.
#     Each read is counted once, so reads are independent and no permutation is
#     needed here.
#
# The null for the position statistic PERMUTES CLUSTER LABELS ACROSS WHOLE READS
# (a read's footprints travel together, cluster sizes preserved). Footprints
# within a read are strongly dependent -- nucleosomes are spaced, not scattered --
# so treating footprints as independent observations is badly anticonservative.
# Permuting whole reads leaves all within-read structure intact and destroys only
# the association with m6A cluster identity, which is the hypothesis.
#
# Note on scale: Wasserstein distance in bp grows with region width and sample
# size, so the raw statistic is NOT comparable between genes -- the null median
# runs from ~45bp on the smallest region to ~160bp on the largest. Compare each
# gene against its own null (which is what the permutation p-value does); never
# pool or average the raw bp across genes.

PROFILE_BW_BP   <- 30L    # KDE bandwidth, ~ nucleosome positional jitter
PROFILE_STEP_BP <- 5L     # profile grid step
PERM_B          <- 2000L  # label permutations for the position null

# 1-D Wasserstein-1 distance in bp: the area between two ECDFs.
wasserstein1 <- function(x, y) {
  if (!length(x) || !length(y)) return(NA_real_)
  pts <- sort(unique(c(x, y)))
  if (length(pts) < 2L) return(0)
  Fx <- ecdf(x)(pts); Fy <- ecdf(y)(pts)
  sum(abs(Fx[-length(pts)] - Fy[-length(pts)]) * diff(pts))
}

# All pairwise W1 for one labelling of the footprints.
pairwise_w1 <- function(mids, lab, levs) {
  s  <- split(mids, factor(lab, levels = levs))
  cb <- combn(length(levs), 2L)
  vapply(seq_len(ncol(cb)), function(j) wasserstein1(s[[cb[1, j]]], s[[cb[2, j]]]),
         numeric(1))
}

# Footprints per read per bp for each read cluster, on a common grid.
cluster_intensity_profiles <- function(gene_sym, chrom_id, view_start, view_end,
                                       bw = PROFILE_BW_BP, step = PROFILE_STEP_BP) {
  topic <- load_topic_read_clusters(topic_outname(gene_sym, chrom_id, view_start, view_end))
  fp    <- read_topic_footprints(topic, chrom_id, view_start, view_end)
  if (!nrow(fp)) return(data.table())
  fp <- merge(fp, topic[, .(RID, cluster)], by = "RID")

  npts <- as.integer((view_end - view_start) %/% step) + 1L
  out  <- list()
  for (grp in sort(unique(topic$cluster))) {
    m  <- fp[cluster == grp, fp_mid]
    nr <- topic[cluster == grp, .N]
    if (length(m) >= 2L) {
      d   <- density(m, bw = bw, from = view_start, to = view_end, n = npts)
      pos <- d$x
      # scale the unit-area density to footprints per read per bp
      y   <- d$y * length(m) / nr
    } else {
      pos <- seq(view_start, view_end, length.out = npts)
      y   <- numeric(npts)
    }
    out[[grp]] <- data.table(read_cluster = grp, pos = pos, intensity = y,
                             n_reads = nr, n_fp = length(m),
                             fp_per_read = length(m) / nr)
  }
  rbindlist(out)[, `:=`(gene = gene_sym, chrom = chrom_id,
                        view_start = view_start, view_end = view_end)][]
}

# Position (permuted W1) and occupancy (Fisher) for one region.
cluster_distribution_tests <- function(gene_sym, chrom_id, view_start, view_end,
                                       B = PERM_B) {
  topic <- load_topic_read_clusters(topic_outname(gene_sym, chrom_id, view_start, view_end))
  fp    <- read_topic_footprints(topic, chrom_id, view_start, view_end)
  if (!nrow(fp)) return(NULL)
  levs    <- sort(unique(topic$cluster))
  ridx    <- match(fp$RID, topic$RID)   # footprint -> its read
  read_cl <- topic$cluster              # one label per READ
  mids    <- fp$fp_mid

  obs <- mean(pairwise_w1(mids, read_cl[ridx], levs), na.rm = TRUE)
  nul <- vapply(seq_len(B), function(b)
    mean(pairwise_w1(mids, sample(read_cl)[ridx], levs), na.rm = TRUE), numeric(1))

  # occupancy: reads carrying any footprint, by cluster
  pr <- merge(topic[, .(RID, cluster)], fp[, .(n = .N), by = RID], by = "RID", all.x = TRUE)
  pr[is.na(n), n := 0L]
  occ <- pr[, .(n_reads = .N, with_fp = sum(n > 0)), by = cluster][order(cluster)]
  occ[, pct := 100 * with_fp / n_reads]
  ft <- fisher.test(as.matrix(occ[, .(with_fp, no_fp = n_reads - with_fp)]),
                    simulate.p.value = TRUE, B = 20000L)

  list(
    gene_stats = data.table(
      gene = gene_sym, chrom = chrom_id,
      view_start = view_start, view_end = view_end,
      n_reads = nrow(topic), n_fp = nrow(fp), n_clusters = length(levs),
      w1_obs = obs, w1_null_med = median(nul), w1_null_p95 = quantile(nul, 0.95),
      w1_fold = obs / median(nul),
      p_position = (sum(nul >= obs) + 1) / (B + 1),
      occ_min_pct = min(occ$pct), occ_max_pct = max(occ$pct),
      p_occupancy = ft$p.value),
    occupancy = occ[, .(gene = gene_sym, read_cluster = cluster, n_reads,
                        with_fp, pct)])
}

# `seed` is not cosmetic: the position p-value is a Monte-Carlo estimate, and a
# gene sitting near q = 0.05 (PTGS2 does) flips either side of it between
# unseeded runs. Fix the seed so the tables are reproducible, and read a
# borderline gene as borderline rather than as significant.
main_cluster_distributions <- function(regions = SELECTED_PROMOTERS, seed = 1L) {
  set.seed(seed)
  prof <- list(); stats <- list(); occ <- list()
  for (i in seq_len(nrow(regions))) {
    r <- regions[i]
    message("Distribution comparison for ", r$gene, " ...")
    p <- cluster_intensity_profiles(r$gene, r$chr, r$start, r$end)
    if (!nrow(p)) { message("  no footprints"); next }
    prof[[r$gene]] <- p
    tt <- cluster_distribution_tests(r$gene, r$chr, r$start, r$end)
    if (is.null(tt)) next
    stats[[r$gene]] <- tt$gene_stats
    occ[[r$gene]]   <- tt$occupancy
    message(sprintf("  position: W1 %.0fbp vs null %.0fbp, p=%.4f | occupancy: %.0f-%.0f%%, p=%.4f",
                    tt$gene_stats$w1_obs, tt$gene_stats$w1_null_med,
                    tt$gene_stats$p_position, tt$gene_stats$occ_min_pct,
                    tt$gene_stats$occ_max_pct, tt$gene_stats$p_occupancy))
  }
  prof_dt  <- rbindlist(prof,  use.names = TRUE, fill = TRUE)
  stats_dt <- rbindlist(stats, use.names = TRUE, fill = TRUE)
  occ_dt   <- rbindlist(occ,   use.names = TRUE, fill = TRUE)
  # BH across genes, separately for each axis
  if (nrow(stats_dt)) stats_dt[, `:=`(q_position  = p.adjust(p_position,  "BH"),
                                      q_occupancy = p.adjust(p_occupancy, "BH"))]

  fwrite(prof_dt,  file.path(OUT_DIR, "cluster_intensity_profiles.tsv.gz"), sep = "\t")
  fwrite(stats_dt, file.path(OUT_DIR, "cluster_distribution_tests.tsv"), sep = "\t")
  fwrite(occ_dt,   file.path(OUT_DIR, "cluster_occupancy_by_gene.tsv"), sep = "\t")
  message("Done. Wrote cluster_intensity_profiles.tsv.gz, cluster_distribution_tests.tsv ",
          "and cluster_occupancy_by_gene.tsv to ", OUT_DIR)
  invisible(list(profiles = prof_dt, stats = stats_dt, occupancy = occ_dt))
}

main <- function() {
  main_nuc_pos()
  main_nuc_spacing()
  main_own_pile_midpoints()
  main_cluster_distributions()
}

if (sys.nframe() == 0L) main()
