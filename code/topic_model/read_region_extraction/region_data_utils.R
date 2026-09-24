# Readers and feature helpers for read-level region plots.
#
# Self-contained: nothing here sources or depends on the fiberhub code. Every input is
# a small file already sliced to the region by extract_region_result_macrophage.sh, so
# these are plain reads of plain files - no tabix querying from R.
#
# Requires: data.table, dplyr, GenomicRanges, IRanges

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(IRanges)
})

FT_BED12_COLS <- c("chr", "start", "end", "RID", "score", "strand",
                   "read_start", "read_end", "rgb", "blockCount",
                   "blockSizes", "blockStarts")
FIRE_BED_COLS <- c("chr", "start", "end", "RID", "score", "strand",
                   "thickStart", "thickEnd", "rgb", "fire_score", "tag")

# FIRE colours the per-read segmentation by class; everything that is not grey
# (nucleosome) or purple (linker) is a FIRE element.
RGB_NUCLEOSOME <- "169,169,169"
RGB_LINKER     <- "147,112,219"


#' Read a headerless bed-like file, returning an empty data.frame if it has no rows.
read_bed <- function(file, col_names = NULL) {
  if (!file.exists(file))
    stop("file not found: ", file)
  if (file.size(file) == 0)
    return(data.frame())
  df <- data.table::fread(file, header = FALSE, sep = "\t", data.table = FALSE)
  if (!is.null(col_names)) {
    if (ncol(df) < length(col_names))
      stop(file, " has ", ncol(df), " columns, expected at least ", length(col_names))
    colnames(df)[seq_along(col_names)] <- col_names
  }
  df
}


#' Expand the blocks of a BED12(+) frame into one row per block.
#'
#' @param df BED12 frame with blockCount/blockSizes/blockStarts (and optionally
#'   blockScores) columns.
#' @param drop_flanking drop the first and last block of every record. `ft extract`
#'   brackets each read with zero-length sentinel blocks at the alignment ends;
#'   FiberHMM's footprint/tf beds do not, so this is TRUE for the former only.
#' @return data.frame with one row per block: chr, RID, strand, read_start, read_end,
#'   block_start (0-based), block_end, size, and block_score when available.
expand_bed12_blocks <- function(df, drop_flanking = FALSE) {
  if (nrow(df) == 0) return(data.frame())

  sizes  <- strsplit(as.character(df$blockSizes), ",", fixed = TRUE)
  starts <- strsplit(as.character(df$blockStarts), ",", fixed = TRUE)
  n_blk  <- lengths(starts)

  if (any(n_blk != lengths(sizes)))
    stop("blockSizes and blockStarts have different lengths")

  has_scores <- "blockScores" %in% colnames(df)
  if (has_scores) {
    scores <- strsplit(as.character(df$blockScores), ",", fixed = TRUE)
    if (any(lengths(scores) != n_blk))
      stop("blockScores length does not match blockStarts")
  }

  keep_row <- n_blk > 0
  if (drop_flanking) keep_row <- n_blk > 2
  if (!any(keep_row)) return(data.frame())

  row_idx <- rep(which(keep_row), n_blk[keep_row])
  blk_start <- as.integer(unlist(starts[keep_row], use.names = FALSE))
  blk_size  <- as.integer(unlist(sizes[keep_row], use.names = FALSE))

  out <- data.frame(
    chr         = df$chr[row_idx],
    RID         = as.character(df$RID)[row_idx],
    strand      = df$strand[row_idx],
    read_start  = df$read_start[row_idx],
    read_end    = df$read_end[row_idx],
    block_start = df$start[row_idx] + blk_start,
    size        = blk_size,
    stringsAsFactors = FALSE)
  out$block_end <- out$block_start + out$size
  if (has_scores)
    out$block_score <- as.numeric(unlist(scores[keep_row], use.names = FALSE))

  if (drop_flanking) {
    # position of each block within its record
    pos_in_rec <- sequence(n_blk[keep_row])
    last_of_rec <- rep(n_blk[keep_row], n_blk[keep_row])
    out <- out[pos_in_rec != 1L & pos_in_rec != last_of_rec, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}


#' Read `ft extract --m6a`/`--cpg` output for a region into one row per modified base.
#'
#' @return data.frame: chr, RID, strand, read_start (1-based), read_end, pos (1-based
#'   position of the modified base).
read_ft_mod_region <- function(file, region_start, region_end) {
  df <- read_bed(file, FT_BED12_COLS)
  if (nrow(df) == 0) return(data.frame())

  # one alignment per read: keep the longest if a read is split
  if (any(duplicated(df$RID))) {
    df <- df %>%
      dplyr::group_by(RID) %>%
      dplyr::slice_max(end - start, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      as.data.frame()
  }

  blocks <- expand_bed12_blocks(df, drop_flanking = TRUE)
  if (nrow(blocks) == 0) return(data.frame())

  # blocks are 0-based half-open and 1 bp wide for a modification call, so the
  # block end is the 1-based coordinate of the modified base
  out <- data.frame(
    chr        = blocks$chr,
    RID        = blocks$RID,
    strand     = blocks$strand,
    read_start = blocks$read_start + 1L,
    read_end   = blocks$read_end,
    pos        = blocks$block_end,
    stringsAsFactors = FALSE)

  out[out$pos >= region_start & out$pos <= region_end, , drop = FALSE]
}


#' Read `ft fire --extract` output (fire.bed) and label each segment.
read_fire_region <- function(file) {
  df <- read_bed(file, FIRE_BED_COLS)
  if (nrow(df) == 0)
    stop("fire.bed is empty: ", file)
  df$RID <- as.character(df$RID)
  df$class <- factor(
    ifelse(df$rgb == RGB_NUCLEOSOME, "nucleosome",
           ifelse(df$rgb == RGB_LINKER, "linker", "FIRE")),
    levels = c("nucleosome", "linker", "FIRE"))
  df$fire_score <- pmin(df$fire_score, 1)
  df
}


#' Read the FIRE peak calls for the region.
#'
#' The v0.1 peak bed has 29 columns; only the peak interval and the FDR are needed.
#' logFDR is stored x10, as in the FIRE track hub.
read_fire_peaks_region <- function(file) {
  df <- read_bed(file)
  if (nrow(df) == 0) return(NULL)
  if (ncol(df) != 29)
    stop(file, " has ", ncol(df), " columns, expected 29 (FIRE v0.1 peaks)")
  data.frame(chr    = df[[1]],
             start  = df[[2]],
             end    = df[[3]],
             FDR    = df[[21]],
             logFDR = df[[22]] / 10,
             stringsAsFactors = FALSE)
}


#' Bin footprint sizes into the classes used for colouring.
#'
#' Sizes outside [min(breaks), max(breaks)] get NA and are dropped when plotting -
#' this is what keeps nucleosome-sized calls out of the footprint track.
assign_size_class <- function(size, breaks) {
  breaks <- sort(breaks)
  if (length(breaks) < 2)
    stop("need at least two size breaks")
  labs <- paste0(head(breaks, -1), "-", tail(breaks, -1), " bp")
  idx <- rep(NA_integer_, length(size))
  for (i in seq_along(labs))
    idx[size >= breaks[i] & size <= breaks[i + 1]] <- i
  factor(labs[idx], levels = labs)
}


#' Read FiberHMM footprint calls for the region and expand them into footprints.
#'
#' @param feature_cols column names for the file: 15 for `tf`
#'   (blockScores/blockEdgeLeft/blockEdgeRight), 13 for `footprint`/`msp` (blockScores).
read_fiberhmm_region <- function(file,
                                 min_score = 0,
                                 size_breaks = c(10, 30, 60, 80)) {
  df <- read_bed(file)
  if (nrow(df) == 0) return(data.frame())

  cols <- switch(as.character(ncol(df)),
    "12" = FT_BED12_COLS,
    "13" = c(FT_BED12_COLS, "blockScores"),
    "15" = c(FT_BED12_COLS, "blockScores", "blockEdgeLeft", "blockEdgeRight"),
    stop(file, " has ", ncol(df), " columns; expected 12, 13 or 15"))
  colnames(df) <- cols

  # FiberHMM beds have no sentinel blocks, unlike `ft extract` output
  blocks <- expand_bed12_blocks(df, drop_flanking = FALSE)
  if (nrow(blocks) == 0) return(data.frame())

  if ("block_score" %in% colnames(blocks) && min_score > 0)
    blocks <- blocks[blocks$block_score >= min_score, , drop = FALSE]

  data.frame(chr        = blocks$chr,
             start      = blocks$block_start,
             end        = blocks$block_end,
             RID        = blocks$RID,
             size       = blocks$size,
             score      = if ("block_score" %in% colnames(blocks)) blocks$block_score else NA_real_,
             read_start = blocks$read_start,
             read_end   = blocks$read_end,
             class      = assign_size_class(blocks$size, size_breaks),
             stringsAsFactors = FALSE)
}


#' Keep only footprints that fall entirely inside a FIRE element of the same read.
#'
#' Uses the read ID as the GRanges seqname, so overlaps are automatically confined to
#' the read a footprint came from.
subset_footprints_in_fire <- function(fps, fire) {
  if (nrow(fps) == 0) return(fps)
  fire_el <- fire[fire$class == "FIRE", , drop = FALSE]
  if (nrow(fire_el) == 0) return(fps[0, , drop = FALSE])

  fps.gr  <- GRanges(seqnames = fps$RID,
                     ranges = IRanges(start = fps$start + 1L, end = fps$end))
  fire.gr <- GRanges(seqnames = fire_el$RID,
                     ranges = IRanges(start = fire_el$start + 1L, end = fire_el$end))
  hits <- findOverlaps(fps.gr, fire.gr, type = "within")
  fps[sort(unique(queryHits(hits))), , drop = FALSE]
}


#' Per-read summary: alignment span, strand, haplotype tag, arrow direction.
build_rids_df <- function(fire) {
  fire %>%
    dplyr::group_by(RID) %>%
    dplyr::summarise(sample_name = unique(sample_name),
                     chr    = unique(chr),
                     start  = min(start),
                     end    = max(end),
                     strand = unique(strand)[1],
                     tag    = unique(tag)[1],
                     .groups = "drop") %>%
    dplyr::mutate(arrow_end = ifelse(strand == "+", "last", "first")) %>%
    as.data.frame()
}


#' Read x position modification matrix.
#'
#' 1 = modified, 0 = covered but unmodified, NA = position not covered by the read.
#' `ft extract` reports the read's full alignment span, so "covered but unmodified"
#' is a real observation, not missing data.
build_met_mat <- function(reads, rids_df, positions) {
  rid_levels <- as.character(rids_df$RID)
  m <- matrix(NA_real_, nrow = length(rid_levels), ncol = length(positions),
              dimnames = list(rid_levels, positions))

  covered <- outer(rids_df$start, positions, "<=") &
             outer(rids_df$end,   positions, ">=")
  m[covered] <- 0

  ri <- match(as.character(reads$RID), rid_levels)
  ci <- match(reads$pos, positions)
  ok <- !is.na(ri) & !is.na(ci)
  m[cbind(ri[ok], ci[ok])] <- 1
  m
}


#' Sliding-window mean over positions (window of +/- n/2 bp).
smooth_mean_slidewindow <- function(x, pos, n = 20) {
  sapply(pos, function(i) mean(x[which(pos >= (i - n / 2) & pos <= (i + n / 2))]))
}


#' Modification pileup over the region, optionally split into groups of reads.
#'
#' The denominator is the number of reads in the group, not the number covering each
#' position, so a position a read does not cover counts as unmethylated.
pileup_reads <- function(reads, rids_df, window_start, window_end,
                         split_by = NULL, smooth = TRUE, window_n = 10) {

  if (!is.null(split_by) && !split_by %in% colnames(rids_df))
    stop("column '", split_by, "' not found in rids_df")

  out <- list()
  for (b in unique(reads$base)) {
    positions <- sort(unique(reads$pos[reads$base == b]))
    positions <- positions[positions >= window_start & positions <= window_end]
    if (length(positions) == 0) next

    # positions are restricted to this base, but the matrix is filled from every read
    met_mat <- build_met_mat(reads, rids_df, positions)

    groups <- if (is.null(split_by)) list(All = seq_len(nrow(rids_df))) else
      split(seq_len(nrow(rids_df)), rids_df[[split_by]])

    for (g in names(groups)) {
      mm  <- met_mat[groups[[g]], , drop = FALSE]
      cov <- nrow(mm)
      met <- colSums(mm, na.rm = TRUE)
      df <- data.frame(pos = positions, base = b, group = g,
                       cov = cov, met = met, frac = met / cov,
                       stringsAsFactors = FALSE)
      if (smooth) {
        df$smooth_cov  <- smooth_mean_slidewindow(df$cov, df$pos, window_n)
        df$smooth_frac <- smooth_mean_slidewindow(df$frac, df$pos, window_n)
      } else {
        df$smooth_cov  <- df$cov
        df$smooth_frac <- df$frac
      }
      out[[length(out) + 1]] <- df
    }
  }
  if (length(out) == 0)
    stop("no modified positions inside the region")
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}


#' Pairwise Euclidean distance over the positions where both rows are non-NA.
#'
#' Pairs sharing no non-NA position keep the initial distance of 100.
calculate_dist <- function(x) {
  d_mx <- matrix(100, nrow(x), nrow(x))
  dimnames(d_mx) <- list(rownames(x), rownames(x))
  for (i in 1:nrow(x)) {
    for (j in i:nrow(x)) {
      ix <- which(!is.na(x[i, ]) & !is.na(x[j, ]))
      if (length(ix) > 0) {
        tmp <- as.matrix(dist(x[c(i, j), ix]))
        d_mx[rownames(tmp), rownames(tmp)] <- tmp
      }
    }
  }
  as.dist(d_mx)
}


#' Hierarchically cluster reads on their footprint calls inside the region.
#'
#' Intersects the footprints with the region, expands them to one row per base pair,
#' and builds a read x position matrix holding the footprint size class as an integer.
#' Positions with no footprint on a read are NA, so two reads are compared only over
#' positions where both carry a size-classed footprint.
#'
#' @return an hclust object, or NULL when no footprint falls in the region.
hclust_reads_by_footprints <- function(fps, region_chr, region_start, region_end,
                                       size_levels) {
  if (is.null(fps) || nrow(fps) == 0) return(NULL)

  fps.gr <- GRanges(seqnames = fps$chr,
                    ranges = IRanges(start = fps$start + 1L, end = fps$end),
                    RID = as.character(fps$RID),
                    class = factor(as.character(fps$class), levels = size_levels))
  fps.gr <- sort(fps.gr)

  region.gr <- GRanges(seqnames = region_chr,
                       ranges = IRanges(start = region_start, end = region_end))

  hits <- findOverlaps(fps.gr, region.gr)
  if (length(hits) == 0) return(NULL)
  inter <- pintersect(fps.gr[queryHits(hits)], region.gr[subjectHits(hits)])

  # one row per base pair covered by a footprint
  n_bp <- width(inter)
  df <- data.frame(
    RID = rep(mcols(inter)$RID, n_bp),
    pos = sequence(n_bp, from = start(inter)),
    class_value = rep(as.numeric(mcols(inter)$class), n_bp),
    stringsAsFactors = FALSE)

  df <- df[order(df$RID, df$pos), , drop = FALSE]
  df <- df[!duplicated(paste0(df$RID, ".", df$pos)), , drop = FALSE]
  if (length(unique(df$RID)) < 2) return(NULL)

  rid_f <- factor(df$RID, levels = unique(df$RID))
  pos_f <- factor(df$pos, levels = unique(df$pos))
  m <- matrix(NA_real_, nlevels(rid_f), nlevels(pos_f),
              dimnames = list(levels(rid_f), levels(pos_f)))
  m[cbind(as.integer(rid_f), as.integer(pos_f))] <- df$class_value

  hclust(calculate_dist(m))
}


#' Read order for plotting: reads without footprints first, then in clustering order.
order_reads_by_footprints <- function(fps, rid_levels, region_chr,
                                      region_start, region_end, size_levels) {
  hc <- hclust_reads_by_footprints(fps, region_chr, region_start, region_end, size_levels)
  if (is.null(hc)) return(rid_levels)
  clustered <- hc$labels[hc$order]
  clustered <- clustered[clustered %in% rid_levels]
  c(setdiff(rid_levels, clustered), clustered)
}
