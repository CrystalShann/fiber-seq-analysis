#' Nucleosome positioning around GENCODE canonical TSS, per Leiden
#' accessibility cluster, across the LPS time course (0, 5, 10, 15 min).
#'
#' Data preparation only. Plotting lives in plot_nuc_pos.R; region selection,
#' function calls and interpretation live in nuc_pos.rmd. The cluster
#' assignments come from
#' code/clustering_TSS_methods/Leiden_Manhattan/leiden_manhattan.Rmd, whose
#' promoter windows are the same canonical TSS +/- VIEW_HALF_WIDTH used here;
#' the loader accepts only the run computed on exactly the plotted window.
#'
#' Inputs
#'   * FiberHMM footprint BED12, one tabix-indexed file per sample x chromosome.
#'     The files carry 13 columns (the last is an all-zero block field); only
#'     the standard BED12 columns are read. FiberHMM writes NO sentinel
#'     first/last blocks, so every block is a real nucleosome footprint and all
#'     of them are used (no size filter).
#'   * FiberHMM TF footprints split by size class: BED6 with the footprint size
#'     in the score column, one tabix-indexed file per size class x sample x
#'     chromosome.
#'   * GENCODE v46 Ensembl-canonical TSS BED: one 20-bp interval per gene with
#'     the annotation "gene_id;transcript_id;gene_name;gene_type;tags".
#'   * The promoter table written by make_promoter_table.R (genes, ENSG,
#'     strand; its own v49 windows are not used) and the Leiden read-cluster
#'     assignments.
#'
#' Coordinates
#'   * BED is 0-based half-open; fp_start/fp_end keep that convention. Read
#'     spans from the Leiden assignments are 1-based closed.
#'   * Relative positions are closed integer intervals in bp from the anchor
#'     (anchor base = 0), oriented so transcription runs left to right:
#'     upstream < 0 < downstream, mirrored for minus-strand genes.
#'   * Promoter anchor = GENCODE canonical TSS of the promoter-table ENSG;
#'     window = anchor +/- VIEW_HALF_WIDTH. Custom regions have no biological
#'     anchor: their anchor is the region centre on the plus strand ("position
#'     relative to region centre") and their window is the region itself.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(Rsamtools)
})

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
# FiberHMM footprint BED12: <FP_ROOT>/LPS0/LPS0_hmm_extracted_footprint_<chr>.bed.gz
FP_ROOT <- "/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_footprint"
# TF footprints by size class: <TF_ROOT>/size10-30/LPS0/LPS0_tf_size10-30_<chr>.bed.gz
TF_ROOT <- "/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_tf/ft_by_size"
TF_SIZE_CLASSES <- c("size10-30", "size40-60", "size60-80")
# GENCODE v46 Ensembl-canonical TSS (the TSS annotation; replaces THP-1 CAGE)
GENCODE_TSS_PATH <- "/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed"

TIMEPOINT_LEVELS <- c("0", "5", "10", "15")   # LPS minutes, 0 = control

# Promoter window = canonical TSS +/- VIEW_HALF_WIDTH: the 2-kb window
# clustered in leiden_manhattan.Rmd.
VIEW_HALF_WIDTH <- 1000L
# Positional bin (bp) of the occupancy profiles. Visualisation only -- the
# Leiden clustering itself runs unbinned (window_size = 0).
PROFILE_BIN_BP <- 10L

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Leiden sample_name "LPS_0" -> FiberHMM directory / file prefix "LPS0"
sample_dir <- function(sample_name) gsub("_", "", sample_name, fixed = TRUE)

# "ENSG00000125538.12" -> "ENSG00000125538"
strip_ensembl_version <- function(x) sub("\\.[0-9]+$", "", x)

# Lines of a tabix-indexed BED overlapping chr:start-end (1-based, closed).
tabix_lines <- function(bed, chr, start, end) {
  if (!file.exists(bed)) stop("missing BED: ", bed)
  tryCatch(
    scanTabix(TabixFile(bed), param = GRanges(chr, IRanges(start, end)))[[1]],
    error = function(e) character(0))
}

# 1-based closed genomic interval -> closed interval in bp relative to the
# anchor, transcription to the right (mirrored on the minus strand).
rel_interval <- function(start1, end1, anchor, strand) {
  minus <- rep_len(strand == "-", length(start1))
  list(rel_start = fifelse(minus, anchor - end1,   start1 - anchor),
       rel_end   = fifelse(minus, anchor - start1, end1 - anchor))
}

# ---------------------------------------------------------------------------
# GENCODE canonical TSS
# ---------------------------------------------------------------------------

# One row per gene: 1-based TSS, strand and the parsed annotation field
# (gene_id, transcript_id, gene_name, gene_type, transcript_tags) plus the
# version-stripped `ensg`, so gene_id -> gene_name and ENSG -> TSS lookups are
# a merge on `ensg`. The 20-bp BED interval is [TSS - 10, TSS + 10) in 0-based
# coordinates (it reproduces the single-base GENCODE v49 TSS of TNF and IL1A
# exactly), so the 1-based TSS base is start + 11.
load_gencode_tss <- function(path = GENCODE_TSS_PATH) {
  dt <- fread(path, header = FALSE,
              col.names = c("chrom", "start", "end", "annotation", "score", "strand"))
  dt[, c("gene_id", "transcript_id", "gene_name", "gene_type", "transcript_tags") :=
       tstrsplit(annotation, ";", fixed = TRUE)]
  dt[, ensg := strip_ensembl_version(gene_id)]
  dt[, tss  := start + 11L]
  dt[, .(chrom, tss, strand, gene_id, ensg, transcript_id, gene_name, gene_type,
         transcript_tags)]
}

# ---------------------------------------------------------------------------
# Regions
# ---------------------------------------------------------------------------

# Promoter table from make_promoter_table.R. `gene` is the project name (IL8),
# `gencode_name` the GENCODE symbol (CXCL8); `ensg_base` drops the version so
# it can be matched to the canonical-TSS table.
load_promoter_table <- function(path) {
  p <- fread(path)
  need <- c("gene", "gencode_name", "ensg", "chr", "start", "end", "tss",
            "gene_start", "gene_end", "gene_strand")
  missing <- setdiff(need, names(p))
  if (length(missing)) stop("promoter table lacks: ", paste(missing, collapse = ", "))
  p[, ensg_base := strip_ensembl_version(ensg)]
  p[]
}

# One row per region with the anchor of its relative coordinate system and the
# window (view_start-view_end, 1-based closed) on which the reads were
# clustered and are plotted, plus that window in relative coordinates.
#   promoters (region_id = gene): anchor = GENCODE canonical TSS of the
#     promoter-table ENSG, strand = canonical strand, window = anchor +/-
#     half_width. The promoter table's own TSS (v49 gene 5' end) is kept as
#     promoter_tss with its signed offset in the direction of transcription
#     (tss_offset = canonical - promoter TSS).
#   custom regions: anchor = region centre, plus strand, window = the region.
build_region_anchors <- function(regions, promoters, gencode_tss,
                                 half_width = VIEW_HALF_WIDTH) {
  reg <- as.data.table(regions)[, .(region_id = as.character(region_id),
                                     chr = as.character(chr),
                                     start = as.integer(start), end = as.integer(end))]
  prom <- promoters[, .(region_id = gene, gene, gencode_name, ensg = ensg_base,
                        gene_strand, promoter_tss = tss)]
  prom <- merge(prom,
                gencode_tss[, .(ensg, canonical_chrom = chrom, canonical_tss = tss,
                                canonical_strand = strand)],
                by = "ensg", all.x = TRUE)
  if (anyNA(prom$canonical_tss))
    stop("no GENCODE canonical TSS for: ",
         paste(prom[is.na(canonical_tss), gene], collapse = ", "))

  a <- merge(reg, prom, by = "region_id", all.x = TRUE)
  a[, region_type := fifelse(is.na(gene), "custom", "promoter")]
  a[region_type == "promoter", `:=`(
    anchor = canonical_tss, strand = canonical_strand,
    anchor_source = "GENCODE v46 canonical TSS",
    tss_offset = fifelse(canonical_strand == "+", canonical_tss - promoter_tss,
                         promoter_tss - canonical_tss),
    view_start = canonical_tss - half_width, view_end = canonical_tss + half_width)]
  a[region_type == "custom", `:=`(
    anchor = (start + end) %/% 2L, strand = "+", anchor_source = "region centre",
    view_start = start, view_end = end)]
  a[, c("rel_view_start", "rel_view_end") := rel_interval(view_start, view_end, anchor, strand)]

  bad <- a[region_type == "promoter" &
             (chr != canonical_chrom | gene_strand != canonical_strand)]
  if (nrow(bad))
    warning("promoter table and canonical TSS disagree on chr/strand for: ",
            paste(bad$gene, collapse = ", "))

  a[, c("canonical_chrom", "canonical_strand", "gene_strand") := NULL]
  setcolorder(a, c("region_id", "region_type", "chr", "start", "end", "gene",
                   "gencode_name", "ensg", "strand", "anchor", "anchor_source",
                   "view_start", "view_end", "rel_view_start", "rel_view_end",
                   "promoter_tss", "canonical_tss", "tss_offset"))
  a[match(reg$region_id, region_id)]
}

# ---------------------------------------------------------------------------
# Leiden cluster assignments (leiden_manhattan.Rmd outputs)
# ---------------------------------------------------------------------------

# TRUE when run_parameters_<region>.tsv (written by leiden_manhattan.Rmd)
# carries every value of `expected`, compared numerically.
leiden_run_matches <- function(run_file, expected) {
  if (!file.exists(run_file)) return(FALSE)
  p <- fread(run_file, colClasses = "character")
  v <- setNames(p$value, p$parameter)
  all(vapply(names(expected), function(n)
    n %in% names(v) && isTRUE(all.equal(as.numeric(v[[n]]), as.numeric(expected[[n]]))),
    logical(1)))
}

# Assignments of leiden_manhattan.Rmd for each region, one file per region:
# <leiden_dir>/<region_id>/bin<window_size>/k<n>/read_cluster_assignments_<region_id>.tsv
# (RID, sample_name, timepoint, cluster, read span start/end 1-based closed).
# The k<n> directory is the run whose parameters match `params` AND whose
# region_start/region_end equal the region's window in `anchors`, so the
# clusters used here were computed on exactly the plotted window; runs on
# other windows are ignored. Reads of all four timepoints were pooled and
# clustered once per region, so cluster labels are comparable across
# timepoints. region_id is added from the file location.
load_leiden_assignments <- function(leiden_dir, anchors, params) {
  out <- list()
  for (i in seq_len(nrow(anchors))) {
    m <- anchors[i]
    expected <- c(params, list(region_start = m$view_start, region_end = m$view_end))
    bin_dir  <- file.path(leiden_dir, m$region_id, paste0("bin", params$window_size))
    k_dirs   <- list.files(bin_dir, pattern = "^k[0-9]+$", full.names = TRUE)
    hit <- k_dirs[vapply(k_dirs, function(d)
      leiden_run_matches(file.path(d, sprintf("run_parameters_%s.tsv", m$region_id)), expected),
      logical(1))]
    if (length(hit) != 1)
      stop(length(hit), " Leiden runs for ", m$region_id, " on ", m$chr, ":", m$view_start,
           "-", m$view_end, " with ", paste(names(params), params, sep = "=", collapse = ", "),
           " under ", bin_dir, "; rerun leiden_manhattan.Rmd on this window")
    a <- fread(file.path(hit, sprintf("read_cluster_assignments_%s.tsv", m$region_id)),
               colClasses = list(character = c("RID", "sample_name", "cluster")))
    a[, region_id := m$region_id]
    out[[m$region_id]] <- a
  }
  a <- rbindlist(out, use.names = TRUE, fill = TRUE)
  need <- c("RID", "sample_name", "timepoint", "cluster", "region_id", "start", "end")
  missing <- setdiff(need, names(a))
  if (length(missing)) stop("Leiden assignments lack: ", paste(missing, collapse = ", "))
  setnames(a, c("start", "end"), c("read_start", "read_end"))
  a[, timepoint := factor(as.character(timepoint), levels = TIMEPOINT_LEVELS)]
  if (anyNA(a$timepoint)) stop("timepoints outside ", paste(TIMEPOINT_LEVELS, collapse = "/"))
  k <- max(as.integer(sub("^cluster", "", a$cluster)))
  a[, cluster := factor(cluster, levels = paste0("cluster", seq_len(k)))]
  if (anyDuplicated(a[, .(region_id, RID)])) stop("duplicated RID within a region")
  a[, .(region_id, RID, sample_name, timepoint, cluster, read_start, read_end)]
}

# Reads with their spans in the relative coordinates of their region. Every
# clustered read is kept, with or without footprints: the reads carry the
# denominators of the occupancy profiles and the baselines of the read tracks.
build_read_table <- function(assignments, anchors) {
  r <- merge(assignments,
             anchors[, .(region_id, region_type, gene, ensg, anchor, strand)],
             by = "region_id")
  r[, c("rel_read_start", "rel_read_end") :=
        rel_interval(read_start, read_end, anchor, strand)]
  setorder(r, region_id, cluster, timepoint, RID)
  r[]
}

# ---------------------------------------------------------------------------
# Footprint extraction (tabix, per region, RID preserved)
# ---------------------------------------------------------------------------

# BED12 block lists -> one row per block (0-based half-open), RID preserved.
expand_bed12_blocks <- function(rows) {
  sizes_l  <- strsplit(rows$blockSizes,  ",", fixed = TRUE)
  starts_l <- strsplit(rows$blockStarts, ",", fixed = TRUE)
  n        <- lengths(sizes_l)
  sizes    <- as.integer(unlist(sizes_l,  use.names = FALSE))
  starts   <- as.integer(unlist(starts_l, use.names = FALSE))
  stopifnot(length(sizes) == length(starts))
  fp_start <- rep(rows$chromStart, n) + starts
  data.table(RID = rep(rows$RID, n), fp_start = fp_start,
             fp_end = fp_start + sizes, fp_size = sizes)[!is.na(fp_size)]
}

# FiberHMM nucleosome footprints (every block, no size filter) of `reads`
# (RID, sample_name) that overlap chr:view_start-view_end, read through tabix
# from each read's own sample BED12.
extract_nuc_footprints <- function(reads, chr, view_start, view_end) {
  out <- list()
  for (sn in unique(reads$sample_name)) {
    sd  <- sample_dir(sn)
    bed <- file.path(FP_ROOT, sd, sprintf("%s_hmm_extracted_footprint_%s.bed.gz", sd, chr))
    txt <- tabix_lines(bed, chr, view_start, view_end)
    if (!length(txt)) next
    rows <- fread(text = txt, header = FALSE, select = c(2, 4, 11, 12),
                  col.names = c("chromStart", "RID", "blockSizes", "blockStarts"),
                  colClasses = list(character = c(4, 11, 12)))
    rows <- rows[RID %in% reads[sample_name == sn, RID]]
    if (!nrow(rows)) next
    fp <- expand_bed12_blocks(rows)
    fp <- fp[fp_start < view_end & fp_end > view_start - 1L]
    if (nrow(fp)) out[[sn]] <- fp[, `:=`(sample_name = sn, fp_class = "nucleosome")]
  }
  rbindlist(out)
}

# TF footprints of `reads` by size class (BED6, size in column 5) overlapping
# the view window; tabix returns exactly the overlapping intervals.
extract_tf_footprints <- function(reads, chr, view_start, view_end,
                                  size_classes = TF_SIZE_CLASSES) {
  out <- list()
  for (cls in size_classes) for (sn in unique(reads$sample_name)) {
    sd  <- sample_dir(sn)
    bed <- file.path(TF_ROOT, cls, sd, sprintf("%s_tf_%s_%s.bed.gz", sd, cls, chr))
    txt <- tabix_lines(bed, chr, view_start, view_end)
    if (!length(txt)) next
    fp <- fread(text = txt, header = FALSE, select = 2:5,
                col.names = c("fp_start", "fp_end", "RID", "fp_size"),
                colClasses = list(character = 4))
    fp <- fp[RID %in% reads[sample_name == sn, RID]]
    if (nrow(fp))
      out[[paste(cls, sn)]] <- fp[, `:=`(sample_name = sn, fp_class = sub("^size", "tf_", cls))]
  }
  rbindlist(out)
}

# Footprints of every clustered read over its region's window, in genomic
# (BED) and anchor-relative coordinates. One row per footprint; join to the
# read table on (region_id, RID, sample_name) to attach cluster and timepoint.
extract_region_footprints <- function(reads, anchors) {
  out <- list()
  for (i in seq_len(nrow(anchors))) {
    m  <- anchors[i]
    rr <- reads[region_id == m$region_id, .(RID, sample_name)]
    if (!nrow(rr)) next
    fp <- rbindlist(list(
      extract_nuc_footprints(rr, m$chr, m$view_start, m$view_end),
      extract_tf_footprints(rr, m$chr, m$view_start, m$view_end)),
      use.names = TRUE)
    if (!nrow(fp)) next
    fp[, c("rel_start", "rel_end") := rel_interval(fp_start + 1L, fp_end, m$anchor, m$strand)]
    fp[, `:=`(fp_mid = (fp_start + 1L + fp_end) / 2, rel_mid = (rel_start + rel_end) / 2,
              region_id = m$region_id)]
    message(sprintf("%-26s %4d reads  %5d nucleosome  %5d TF footprints",
                    m$region_id, nrow(rr), fp[fp_class == "nucleosome", .N],
                    fp[fp_class != "nucleosome", .N]))
    out[[m$region_id]] <- fp
  }
  fp <- rbindlist(out)
  setcolorder(fp, c("region_id", "RID", "sample_name", "fp_class", "fp_start", "fp_end",
                    "fp_mid", "rel_start", "rel_end", "rel_mid", "fp_size"))
  fp[]
}

# ---------------------------------------------------------------------------
# Summaries by region x cluster (x timepoint)
# ---------------------------------------------------------------------------

# Grid positions (bin centres) covered by each closed interval
# [start_col, end_col]: the input rows repeated once per covered position,
# with a `pos` column.
covered_positions <- function(dt, start_col, end_col, positions) {
  step <- positions[2] - positions[1]
  i0 <- pmax(1L, as.integer(ceiling((dt[[start_col]] - positions[1]) / step)) + 1L)
  i1 <- pmin(length(positions), as.integer(floor((dt[[end_col]] - positions[1]) / step)) + 1L)
  n  <- pmax(0L, i1 - i0 + 1L)
  out <- dt[rep(seq_len(nrow(dt)), n)]
  out[, pos := positions[sequence(n, from = i0)]]
  out[]
}

# Fraction of reads with a nucleosome covering each position of the region's
# window (bin centres on a grid through 0), per group:
#   numerator   reads of the group with a nucleosome footprint covering pos
#   denominator reads of the group whose alignment spans pos
# `by` must include region_id (relative coordinates are per region).
occupancy_profile <- function(footprints, reads, anchors, by = c("region_id", "cluster"),
                              bin_size = PROFILE_BIN_BP) {
  half <- ceiling(max(abs(c(anchors$rel_view_start, anchors$rel_view_end))) / bin_size) * bin_size
  positions <- seq(-half, half, by = bin_size)
  nuc <- footprints[fp_class == "nucleosome", .(region_id, RID, rel_start, rel_end)]
  covered <- unique(covered_positions(nuc, "rel_start", "rel_end", positions)[
    , .(region_id, RID, pos)])
  cols <- c(unique(c("region_id", "RID", by)), "rel_read_start", "rel_read_end")
  span <- covered_positions(reads[, cols, with = FALSE],
                            "rel_read_start", "rel_read_end", positions)
  span[, occupied := FALSE]
  span[covered, on = .(region_id, RID, pos), occupied := TRUE]
  out <- span[, .(n_reads = .N, n_occupied = sum(occupied)), by = c(by, "pos")]
  out[, fraction := n_occupied / n_reads]
  out <- merge(out, anchors[, .(region_id, rel_view_start, rel_view_end)], by = "region_id")
  out <- out[pos >= rel_view_start & pos <= rel_view_end]
  out[, c("rel_view_start", "rel_view_end") := NULL]
  setorderv(out, c(by, "pos"))
  out[]
}

# Per-read positional summary of the nucleosome footprints around the anchor
# (TSS for promoters, region centre for custom regions):
#   minus1_mid / plus1_mid  midpoint of the nearest upstream (-1) / downstream
#                           (+1) nucleosome, by the side of the anchor its
#                           midpoint falls on
#   nfr_width               bases between the -1 footprint end and the +1
#                           footprint start (NA unless both exist)
#   anchor_occupied         a nucleosome footprint covers the anchor base
#   n_nuc                   nucleosome footprints overlapping the window
read_position_summary <- function(footprints, reads) {
  nuc  <- footprints[fp_class == "nucleosome"]
  up   <- nuc[rel_mid < 0][order(-rel_mid),
              .(minus1_mid = rel_mid[1], minus1_end = rel_end[1]), by = .(region_id, RID)]
  down <- nuc[rel_mid >= 0][order(rel_mid),
              .(plus1_mid = rel_mid[1], plus1_start = rel_start[1]), by = .(region_id, RID)]
  occ  <- nuc[, .(anchor_occupied = any(rel_start <= 0 & rel_end >= 0), n_nuc = .N),
              by = .(region_id, RID)]
  s <- reads[, .(region_id, gene, RID, sample_name, timepoint, cluster)]
  for (x in list(up, down, occ)) s <- merge(s, x, by = c("region_id", "RID"), all.x = TRUE)
  s[is.na(n_nuc), `:=`(n_nuc = 0L, anchor_occupied = FALSE)]
  s[, nfr_width := as.numeric(plus1_start - minus1_end - 1L)]
  setorder(s, region_id, cluster, timepoint, RID)
  s[]
}

# Median positional summaries per (region, cluster, timepoint).
position_summary_by_group <- function(read_positions) {
  read_positions[, .(
    n_reads        = .N,
    n_minus1       = sum(!is.na(minus1_mid)),
    minus1_median  = as.numeric(median(minus1_mid, na.rm = TRUE)),
    n_plus1        = sum(!is.na(plus1_mid)),
    plus1_median   = as.numeric(median(plus1_mid, na.rm = TRUE)),
    nfr_median     = as.numeric(median(nfr_width, na.rm = TRUE)),
    frac_anchor_occupied = mean(anchor_occupied)
  ), by = .(region_id, cluster, timepoint)][order(region_id, cluster, timepoint)]
}
