#' Utility functions for nucleosome distribution plots
suppressPackageStartupMessages({
 library(GenomicRanges)
 library(rtracklayer)
 library(ggplot2)
 library(ggpubr)
 library(dplyr)
 library(data.table)
 library(Gviz)
})

# bed12 -> nucleosome midpoints
calculate_nucleosome_midpoints <- function(bed12) {
  bed12[, {
    parsed <- parse_bed12_blocks(blockSizes, blockStarts,
                                 drop_terminal_blocks = TRUE)
    sizes <- parsed$sizes
    starts <- parsed$starts
    if (!length(sizes)) return(NULL)
    abs_start <- chromStart + starts
    abs_end   <- chromStart + starts + sizes
    midpoint  <- as.integer((abs_start + abs_end) / 2)
    .(chrom = chrom, nuc_mid = midpoint, nuc_size = sizes)
  }, by = .(chrom, chromStart, chromEnd, name, strand)]
}

# load a six-column TSS BED into GRanges
load_tss_annotations <- function(tss_path) {
  tss <- fread(tss_path, header = FALSE,
               col.names = c("chrom", "start", "end", "gene_id", "score", "strand"))
  GRanges(seqnames = tss$chrom,
          ranges   = IRanges(tss$start + 1L, tss$end),
          strand   = tss$strand,
          gene_id  = tss$gene_id)
}

# convert midpoint summary table -> GRanges
nuc_midpoints_to_granges <- function(nucs) {
  GRanges(seqnames = nucs$chrom,
          ranges   = IRanges(nucs$nuc_mid + 1L, nucs$nuc_mid + 1L),
          strand   = "*",
          nuc_size = nucs$nuc_size)
}

# find nucleosome midpoints within a symmetric window around TSS
find_nucs_near_tss <- function(nuc_gr, tss_gr, window_bp = 2000L) {
  tss_windows <- promoters(tss_gr, upstream = window_bp, downstream = window_bp + 1L)
  hits <- findOverlaps(nuc_gr, tss_windows) 
  nuc_pos <- start(nuc_gr)[queryHits(hits)]
  tss_pos <- start(tss_gr)[subjectHits(hits)]
  tss_str <- as.character(strand(tss_gr))[subjectHits(hits)]
  data.table(
    seqnames = as.character(seqnames(tss_gr))[subjectHits(hits)],
    strand   = tss_str,
    nuc_pos  = nuc_pos,
    tss_pos  = tss_pos,
    rel_pos  = ifelse(tss_str == "+", nuc_pos - tss_pos, tss_pos - nuc_pos)
  )
}

# create a metaprofile line plot of nucleosome midpoints around the TSS
create_meta_profile_plot <- function(plot_dt, window_bp = 2000L) {
  bin_size <- 10L
  plot_dt[, bin := round(rel_pos / bin_size) * bin_size]
  meta <- plot_dt[, .N, by = bin][order(bin)]
  ggplot(meta, aes(x = bin, y = N)) +
    geom_line(linewidth = 0.6) +
    geom_smooth(method = "loess", span = 0.05, se = FALSE,
                color = "steelblue", linewidth = 1) +
    scale_x_continuous(breaks = seq(-window_bp, window_bp, 500),
                       labels = function(x) paste0(x / 1000, " kb")) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "red", linewidth = 0.5) +
    labs(x = "Distance to TSS", y = "Nucleosome midpoint count",
         title = "Nucleosome midpoint distribution around TSS") +
    theme_bw(base_size = 12)
}

# build normalized metaprofile data.table for a set of TSS
make_meta_normalized_TSS <- function(nuc_gr, tss_gr, window_bp = 2000L,
                                     bin_size = 10L, label = NULL) {
  plot_dt <- find_nucs_near_tss(nuc_gr, tss_gr, window_bp = window_bp)
  plot_dt[, bin := round(rel_pos / bin_size) * bin_size]
  meta <- plot_dt[, .N, by = bin][order(bin)]
  meta[, N := N / length(tss_gr)]
  if (!is.null(label)) meta[, source := label]
  meta
}

# raw counts per bin around TSS
make_meta_counts <- function(nuc_gr, tss_gr, window_bp = 2000L,
                             bin_size = 10L, label = NULL) {
  plot_dt <- find_nucs_near_tss(nuc_gr, tss_gr, window_bp = window_bp)
  plot_dt[, bin := round(rel_pos / bin_size) * bin_size]
  meta <- plot_dt[, .N, by = bin][order(bin)]
  if (!is.null(label)) meta[, source := label]
  meta
}

# read one or more BED12 files and return a GRanges of nucleosome midpoints
load_nucleosome_annotations <- function(bed_files) {
  beds <- lapply(bed_files, fread, header = FALSE,
                 col.names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
                               "thickStart", "thickEnd", "itemRgb", "blockCount",
                               "blockSizes", "blockStarts"))
  mids <- rbindlist(lapply(beds, calculate_nucleosome_midpoints))
  nuc_midpoints_to_granges(mids)
}

# tidy up blockSizes into a long table of nucleosome sizes
extract_nuc_sizes <- function(bed12) {
  bed12[, blockSizes := strsplit(blockSizes, ",")]
  bed12 <- bed12[, .(chrom, blockSizes)]
  bed12 <- bed12[, .(blockSizes = unlist(blockSizes)), by = chrom]
  bed12[, blockSizes := as.numeric(blockSizes)]
  bed12[, .(chrom, nuc_size = blockSizes)]
}

# bin nucleosome sizes into canonical classes
bin_nuc_sizes <- function(dt) {
  dt[, size_class := fcase(
    nuc_size <  100,                   "< 100 bp",
    nuc_size >= 100 & nuc_size < 130,  "100-130 bp",
    nuc_size >= 130 & nuc_size < 160,  "130-160 bp",
    nuc_size >= 160 & nuc_size < 190,  "160-190 bp",
    nuc_size >= 190,                   "> 190 bp"
  )]
  dt[, size_class := factor(size_class, levels = c(
    "< 100 bp",
    "100-130 bp",
    "130-160 bp",
    "160-190 bp",
    "> 190 bp"
  ))]
  dt
}

process_chr_bed <- function(bed_file) {
  bed12 <- fread(bed_file, header = FALSE,
                 col.names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
                               "thickStart", "thickEnd", "itemRgb",
                               "blockCount", "blockSizes", "blockStarts"))
  extract_nuc_sizes(bed12)
}

# plot nucleosome size distribution for a sample chromosome
plot_nuc_size_dist <- function(bed_files, chr = "chr1", sample_name = "sample") {
  pattern <- paste0("\\.", chr, "\\.")
  matched <- bed_files[grepl(pattern, bed_files)]
  if (length(matched) == 0L) {
    stop("No BED file matched chromosome: ", chr)
  }
  nucs <- rbindlist(lapply(matched, process_chr_bed))
  nucs <- bin_nuc_sizes(nucs)

  size_class_colors <- c(
    "< 100 bp"   = "#d62728",
    "100-130 bp" = "#ff7f0e",
    "130-160 bp" = "#2ca02c",
    "160-190 bp" = "#1f77b4",
    "> 190 bp"   = "#9467bd"
  )

  counts <- nucs[, .N, by = size_class]

  ggplot(counts, aes(x = size_class, y = N, fill = size_class)) +
    geom_col() +
    geom_text(aes(label = format(N, big.mark = ",")), vjust = -0.4, size = 3) +
    scale_fill_manual(values = size_class_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                       labels = scales::comma) +
    labs(
      x = "Size class",
      y = "Number of nucleosomes",
      title = sprintf("Nucleosome count per size class on %s", chr),
      subtitle = sprintf("Sample: %s | total = %s", sample_name, format(nrow(nucs), big.mark = ","))
    ) +
    theme_bw(base_size = 12)
}


# BED9 (chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb)
load_tss_annotations_bed9 <- function(tss_path) {
  cols <- c("chrom","start","end","name","score","strand",
            "thickStart","thickEnd","itemRgb")
  tss <- data.table::fread(tss_path, header = FALSE, col.names = cols)
  GenomicRanges::GRanges(
    seqnames = tss$chrom,
    ranges   = IRanges::IRanges(tss$start + 1L, tss$end),  # BED is 0‑based, half‑open
    strand   = tss$strand,
    tss_id   = tss$name,
    score    = tss$score
  )
}


# function used to load CAGE TSS annotations with gene_id in column 4 instead of name
load_tss_annotations_bed9_with_gene_id <- function(path) {
  cols <- c("chrom","start","end","gene_id","score","strand","thickStart","thickEnd","itemRgb")
  dt <- data.table::fread(path, header = FALSE, col.names = cols)
  GenomicRanges::GRanges(
    seqnames = dt$chrom,
    ranges   = IRanges::IRanges(dt$start + 1L, dt$end),
    strand   = dt$strand,
    gene_id  = dt$gene_id,
    score    = dt$score
  )
}


# read a generic BED12 file produced by ft extract
read_fiber_bed12 <- function(path) {
  data.table::fread(
    path,
    header = FALSE,
    col.names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
                  "thickStart", "thickEnd", "itemRgb", "blockCount",
                  "blockSizes", "blockStarts")
  )
}

# parse BED12 block columns and optionally drop the first/last sentinel blocks
parse_bed12_blocks <- function(block_sizes, block_starts, drop_terminal_blocks = TRUE) {
  sizes <- as.integer(strsplit(block_sizes, ",")[[1]])
  starts <- as.integer(strsplit(block_starts, ",")[[1]])
  keep <- !is.na(sizes) & !is.na(starts)
  sizes <- sizes[keep]
  starts <- starts[keep]

  if (drop_terminal_blocks && length(sizes) >= 2L) {
    idx <- seq_along(sizes)
    idx <- idx[-c(1L, length(idx))]
    sizes <- sizes[idx]
    starts <- starts[idx]
  }

  list(sizes = sizes, starts = starts)
}

# expand BED12 features into one row per biological feature block
expand_fiber_bed12_features <- function(bed12, feature_name = "feature",
                                        drop_terminal_blocks = TRUE) {
  bed12[, {
    parsed <- parse_bed12_blocks(blockSizes, blockStarts,
                                 drop_terminal_blocks = drop_terminal_blocks)
    sizes <- parsed$sizes
    starts <- parsed$starts
    if (!length(sizes)) return(NULL)

    abs_start <- chromStart + starts
    abs_end <- chromStart + starts + sizes
    midpoint <- as.integer((abs_start + abs_end) / 2)

    .(
      chrom = chrom,
      feature_start = abs_start,
      feature_end = abs_end,
      feature_pos = midpoint,
      feature_size = sizes,
      feature_name = feature_name
    )
  }, by = .(chrom, chromStart, chromEnd, name, strand)]
}

# convert expanded feature positions into 1 bp GRanges for metaprofiles
feature_positions_to_granges <- function(features) {
  GRanges(
    seqnames = features$chrom,
    ranges = IRanges(features$feature_pos + 1L, features$feature_pos + 1L),
    strand = "*",
    feature_size = features$feature_size,
    feature_name = features$feature_name
  )
}

# load ft extract BED12 files and return midpoint-like feature positions
load_fiber_feature_positions <- function(bed_files, feature_name = "feature",
                                         drop_terminal_blocks = TRUE) {
  beds <- lapply(bed_files, read_fiber_bed12)
  features <- rbindlist(lapply(
    beds,
    function(bed12) {
      expand_fiber_bed12_features(
        bed12,
        feature_name = feature_name,
        drop_terminal_blocks = drop_terminal_blocks
      )
    }
  ), use.names = TRUE, fill = TRUE)

  if (!nrow(features)) {
    return(GRanges())
  }

  feature_positions_to_granges(features)
}

# build a zero-filled metaprofile for any point-like GRanges around TSS
make_position_meta <- function(position_gr, tss_gr, window_bp = 2000L,
                               bin_size = 10L, normalize_by_tss = FALSE,
                               label = NULL) {
  bins <- data.table(bin = seq(-window_bp, window_bp, by = bin_size))
  if (!length(position_gr) || !length(tss_gr)) {
    bins[, N := 0]
    if (!is.null(label)) bins[, source := label]
    return(bins[])
  }

  plot_dt <- find_nucs_near_tss(position_gr, tss_gr, window_bp = window_bp)
  if (!nrow(plot_dt)) {
    bins[, N := 0]
    if (!is.null(label)) bins[, source := label]
    return(bins[])
  }

  plot_dt[, bin := round(rel_pos / bin_size) * bin_size]
  meta <- plot_dt[, .N, by = bin]
  meta <- merge(bins, meta, by = "bin", all.x = TRUE, sort = TRUE)
  meta[is.na(N), N := 0]
  if (normalize_by_tss) meta[, N := N / length(tss_gr)]
  if (!is.null(label)) meta[, source := label]
  meta[order(bin)]
}



# ─── PER-READ AWARE FUNCTIONS  ────────────────────────────────

# Step 1: Build GRanges of reads (one row per sequencing molecule) with RID
build_reads_granges <- function(bed12) {
  GRanges(
    seqnames = bed12$chrom,
    ranges   = IRanges::IRanges(bed12$chromStart + 1L, bed12$chromEnd),
    strand   = "*",
    RID      = bed12$name
  )
}

# Step 2: Extract nucleosomes per read, preserving RID
extract_nucleosomes_per_read <- function(bed12) {
  nucs <- bed12[, {
    sizes  <- as.integer(strsplit(blockSizes, ",")[[1]])
    starts <- as.integer(strsplit(blockStarts, ",")[[1]])
    keep <- !is.na(sizes) & !is.na(starts)
    sizes <- sizes[keep]
    starts <- starts[keep]
    if (length(sizes) >= 2L) {
      idx <- seq_along(sizes)
      idx <- idx[-c(1L, length(idx))]  # Remove first and last block
      sizes <- sizes[idx]
      starts <- starts[idx]
    }
    if (!length(sizes)) return(NULL)
    abs_start <- chromStart + starts
    abs_end   <- chromStart + starts + sizes
    nuc_mid <- as.integer((abs_start + abs_end) / 2)
    .(RID = name, nuc_mid = nuc_mid, nuc_size = sizes)
  }, by = .(chrom, chromStart, chromEnd, name, strand)]
  nucs
}

# Main per-read function: Find nucleosomes on reads that overlap TSS windows
find_nucs_near_tss_per_read <- function(bed12, tss_gr, window_bp = 2000L, n_ranks = 5L) {
  # Step 1: Build reads GRanges
  reads_gr <- build_reads_granges(bed12)
  
  # Step 2: Extract nucleosomes per read
  nuc_read_cols <- extract_nucleosomes_per_read(bed12)
  
  # Step 3: Create TSS lookup table
  tss_dt <- data.table(
    tss_idx    = seq_along(tss_gr),
    tss_pos    = start(tss_gr),
    tss_strand = as.character(strand(tss_gr))
  )
  
  # Step 4: Find reads overlapping TSS windows
  tss_windows_gr <- promoters(tss_gr, upstream = window_bp, downstream = window_bp + 1L)
  hits_tss <- findOverlaps(reads_gr, tss_windows_gr)
  
  read_tss_pairs <- data.table(
    RID     = reads_gr$RID[S4Vectors::queryHits(hits_tss)],
    tss_idx = S4Vectors::subjectHits(hits_tss)
  )
  
  # Step 5: Join reads, TSS, and nucleosomes
  rtp_nuc <- merge(read_tss_pairs, nuc_read_cols, by = "RID", allow.cartesian = TRUE)
  rtp_nuc <- merge(rtp_nuc, tss_dt, by = "tss_idx", allow.cartesian = TRUE)
  
  # Step 6: Compute strand-aware offsets
  rtp_nuc[, rel_pos  := fifelse(tss_strand == "+",
                                (nuc_mid + 1L) - tss_pos,
                                tss_pos - (nuc_mid + 1L))]
  rtp_nuc[, abs_dist := abs(rel_pos)]
  rtp_nuc <- rtp_nuc[abs_dist <= window_bp]
  rtp_nuc[, direction := fifelse(rel_pos < 0L, "upstream", "downstream")]

  # Step 7: Rank nucleosomes within each (TSS, read, direction)
  setorder(rtp_nuc, tss_idx, RID, direction, abs_dist)
  rtp_nuc[, rank := seq_len(.N), by = .(tss_idx, RID, direction)]
  rtp_nuc <- rtp_nuc[rank <= n_ranks]
  
  rtp_nuc[, signed_rank := fifelse(direction == "upstream", -rank, rank)]
  
  rtp_nuc
}

# Per-read version with gene information
find_nucs_near_tss_with_gene_per_read <- function(bed12, tss_gr, window_bp = 2000L, n_ranks = 5L) {
  # Build everything per-read first
  rtp_nuc <- find_nucs_near_tss_per_read(bed12, tss_gr, window_bp, n_ranks)
  
  # Add gene_id from TSS GRanges
  gene_dt <- data.table(
    tss_idx = seq_along(tss_gr),
    gene_id = as.character(mcols(tss_gr)$gene_id)
  )
  rtp_nuc <- merge(rtp_nuc, gene_dt, by = "tss_idx", all.x = TRUE)
  
  rtp_nuc
}

# ─── META-PROFILE FUNCTIONS FOR PER-READ DATA ──────────────────────────────

# Raw meta-profile: counts per bin (per-read aware)
# Returns: bin, N (count), n_reads (number of unique reads per bin)
make_meta_per_read_raw <- function(bed12, tss_gr, window_bp = 2000L, bin_size = 10L, label = NULL) {
  plot_dt <- find_nucs_near_tss_per_read(bed12, tss_gr, window_bp = window_bp)
  plot_dt[, bin := round(rel_pos / bin_size) * bin_size]
  
  meta <- plot_dt[, .(
    N = .N,                          # Total nucleosome observations
    n_reads = uniqueN(RID),          # Number of unique reads per bin
    n_tss = uniqueN(tss_idx)         # Number of unique TSS per bin
  ), by = bin][order(bin)]
  
  if (!is.null(label)) meta[, source := label]
  meta
}

# Per-read normalized meta-profile: mean nucleosomes per read
# Accounts for varying numbers of reads per TSS
make_meta_per_read_normalized <- function(bed12, tss_gr, window_bp = 2000L, bin_size = 10L, label = NULL) {
  plot_dt <- find_nucs_near_tss_per_read(bed12, tss_gr, window_bp = window_bp)
  plot_dt[, bin := round(rel_pos / bin_size) * bin_size]
  
  # Count unique (read, TSS) pairs first to normalize
  read_tss_pairs <- plot_dt[, .(n_reads = uniqueN(RID)), by = tss_idx]
  total_read_tss_pairs <- sum(read_tss_pairs$n_reads)
  
  meta <- plot_dt[, .(
    N = .N / total_read_tss_pairs,      # Normalize by total read-TSS pairs
    n_reads = uniqueN(RID),
    n_tss = uniqueN(tss_idx)
  ), by = bin][order(bin)]
  
  if (!is.null(label)) meta[, source := label]
  meta
}

# Generic meta-profile function: flexible normalization
# normalize_by: "none" (raw counts), "tss" (per TSS), "read_tss" (per read-TSS pair)
make_position_meta_per_read <- function(bed12, tss_gr, window_bp = 2000L, bin_size = 10L,
                                        normalize_by = "none", label = NULL) {
  plot_dt <- find_nucs_near_tss_per_read(bed12, tss_gr, window_bp = window_bp)
  plot_dt[, bin := round(rel_pos / bin_size) * bin_size]
  
  meta <- plot_dt[, .(
    N_raw = .N,
    n_reads = uniqueN(RID),
    n_tss = uniqueN(tss_idx)
  ), by = bin]
  
  # Apply normalization
  if (normalize_by == "tss") {
    meta[, N := N_raw / length(tss_gr)]
  } else if (normalize_by == "read_tss") {
    # Count unique (read, TSS) pairs
    read_tss_pairs <- plot_dt[, .(n_reads = uniqueN(RID)), by = tss_idx]
    total_read_tss_pairs <- sum(read_tss_pairs$n_reads)
    meta[, N := N_raw / total_read_tss_pairs]
  } else {
    meta[, N := N_raw]
  }
  
  setorder(meta, bin)
  meta[, N_raw := NULL]  # Drop intermediate column
  
  if (!is.null(label)) meta[, source := label]
  meta
}

# Gene-level per-read summary
summarize_gene_tss_per_read <- function(rtp_nuc_with_gene) {
  if (!nrow(rtp_nuc_with_gene) || !"gene_id" %in% names(rtp_nuc_with_gene)) {
    return(data.table(
      gene_id = character(),
      mean_rel_pos_bp = numeric(),
      median_rel_pos_bp = numeric(),
      mean_abs_rel_pos_bp = numeric(),
      median_abs_rel_pos_bp = numeric(),
      sd_rel_pos_bp = numeric(),
      n_nuc = integer(),
      n_reads = integer(),
      n_tss = integer()
    ))
  }
  
  rtp_nuc_with_gene[, .(
    mean_rel_pos_bp = mean(rel_pos, na.rm = TRUE),
    median_rel_pos_bp = as.integer(median(rel_pos, na.rm = TRUE)),
    mean_abs_rel_pos_bp = mean(abs_dist, na.rm = TRUE),
    median_abs_rel_pos_bp = as.integer(median(abs_dist, na.rm = TRUE)),
    sd_rel_pos_bp = sd(rel_pos, na.rm = TRUE),
    n_nuc = .N,
    n_reads = uniqueN(RID),
    n_tss = uniqueN(tss_idx)
  ), by = gene_id][order(mean_abs_rel_pos_bp, -n_nuc)]
}

