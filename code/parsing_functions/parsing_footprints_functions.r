# parsing_footprints_functions.r
#
# Shared Fiber-seq m6A/CpG and footprint parsing, read metadata, and methylation
# matrix functions

# ---------------------------------------------------------------------------
# Query a bgzipped, tabix-indexed BED file for the records overlapping a region.
# Rsamtools::scanTabix returns the matching lines as text, which read.table
# then splits on tabs into columns.
#
# Inputs:
#   ft_extracted_file - path to a .bed.gz file with its .tbi index alongside
#   region_gr         - GRanges of the region(s) to query
# Output:
#   data.frame, one row per overlapping record, unnamed columns V1, V2, ...;
#   an empty data.frame() when no record overlaps
# ---------------------------------------------------------------------------
read_tabix_region <- function(ft_extracted_file, region_gr) {
  tabix_index <- Rsamtools::TabixFile(ft_extracted_file)
  compressed_records <- Rsamtools::scanTabix(tabix_index, param = region_gr)
  # scanTabix gives one character vector of lines per range; pool them
  raw_text <- unlist(compressed_records, use.names = FALSE)

  if (length(raw_text) == 0) {
    return(data.frame())
  }

  read.table(
    text = raw_text,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Expected column count of each supported footprint / modification format:
#   bed12_fibertools   12  fibertools `ft extract` BED12 (m6A, CpG, nuc, msp);
#                          first and last blocks are 1-bp sentinels
#   bed13_fiberhmm     13  FiberHMM BED12 + per-block scores; no sentinels
#   bed15_fiberhmm_tf  15  FiberHMM TF BED12 + per-block tq, left edge, right edge
#   bed6_per_sample     6  one interval per row, one file per sample
#   bed4_pooled         4  one interval per row, samples pooled (no sample ID)
#
# Inputs:
#   format - a single format name from the table above
# Output:
#   integer column count; stops for an unknown or non-scalar format
# ---------------------------------------------------------------------------
footprint_format_columns <- function(format) {
  columns <- c(bed12_fibertools = 12L, bed13_fiberhmm = 13L,
               bed15_fiberhmm_tf = 15L, bed6_per_sample = 6L, bed4_pooled = 4L)
  if (length(format) != 1L || is.na(format) || !format %in% names(columns))
    stop("Unsupported footprint format: ", paste(format, collapse = ", "))
  unname(columns[[format]])
}

# ---------------------------------------------------------------------------
# Stop with a message saying which file has the wrong number of columns for
# its configured format.
#
# Inputs:
#   path   - file the table was read from (only used in the message)
#   format - configured format name (see footprint_format_columns())
#   actual - number of columns actually found
# Output:
#   none; always raises an error
# ---------------------------------------------------------------------------
footprint_column_error <- function(path, format, actual) {
  stop("Footprint file ", path, " (configured format '", format, "'): expected ",
       footprint_format_columns(format), " columns, found ", actual, call. = FALSE)
}


# ---------------------------------------------------------------------------
# Keep one alignment per read. A read with supplementary alignments appears on
# several rows with the same RID; sorting by RID and then by decreasing
# alignment length (end - start) puts the longest alignment first, and
# duplicated() drops the rest.
#
# Inputs:
#   df - table with RID, start and end columns (one row per alignment)
# Output:
#   data.frame with one row per RID (its longest alignment), sorted by RID
# ---------------------------------------------------------------------------
keep_longest_alignment <- function(df) {
  df <- as.data.frame(df)
  if (!nrow(df)) return(df)
  df <- df[order(df$RID, -(df$end - df$start)), , drop = FALSE]
  df[!duplicated(df$RID), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Read a fibertools BED12 file (`ft extract` output), either the whole file or
# only the reads overlapping a region (tabix query), and name its 12 columns.
#
# Inputs:
#   bed_file          - path to the BED12 (.bed.gz; needs a .tbi index when
#                       region is given)
#   region            - NULL to read the whole file, or a GRanges to query
#   longest_alignment - TRUE keeps only the longest alignment of each read
#                       (keep_longest_alignment())
# Output:
#   one row per read alignment with columns chr, start, end (0-based BED
#   coordinates of the alignment), RID, score, strand, x1, x2 (thickStart /
#   thickEnd), rgb, blockCount, blockSizes, blockStarts (comma-separated
#   strings). A data.table when the whole file is read and
#   longest_alignment = FALSE, otherwise a data.frame
# ---------------------------------------------------------------------------
read_ft_bed12 <- function(bed_file, region = NULL, longest_alignment = FALSE) {
  if (is.null(region)) {
    df <- data.table::fread(bed_file, header = FALSE)
  } else {
    if (!file.exists(paste0(bed_file, ".tbi"))) stop("Missing tabix index for ", bed_file)
    df <- read_tabix_region(bed_file, region)
  }
  if (!nrow(df)) return(as.data.frame(df))
  if (ncol(df) != 12L) stop("Expected fibertools BED12 in ", bed_file)
  names(df) <- c('chr','start','end','RID','score','strand','x1','x2','rgb','blockCount','blockSizes','blockStarts')
  if (longest_alignment) {
    df <- keep_longest_alignment(df)
  }
  return(df)
}


# ---------------------------------------------------------------------------
# Expand a BED12 table (one row per read) into one row per block. Each block
# is one feature on the read: a 1-bp m6A / CpG call, or a nucleosome, MSP or
# footprint interval. Block i starts at chromStart + blockStarts[i] and is
# blockSizes[i] long.
#
# fibertools BED12 files put a 1-bp sentinel block at each end of every read;
# drop_sentinels = TRUE removes the first and last block of every RID.
#
# Inputs:
#   bed12_df               - BED12 table with 12 (fibertools), 13 (FiberHMM,
#                            + blockScores) or 15 (FiberHMM TF, + blockScores,
#                            blockEdgeLeft, blockEdgeRight) columns
#   include_read_start_end - TRUE adds read_start / read_end (BED columns 7-8,
#                            unchanged) to every block
#   drop_sentinels         - TRUE drops the first and last block of each read
#   keep_block_scores      - TRUE adds block_score (BED column 13) when present
#   format                 - optional format name (footprint_format_columns()).
#                            When given, the column count is checked and
#                            drop_sentinels is set by the format (TRUE only for
#                            bed12_fibertools); a conflicting explicit
#                            drop_sentinels is an error. When NULL, a 12-column
#                            table is required unless sentinels are kept or
#                            block scores requested (then 12 or 13 columns)
#   source                 - file name used in error messages only
#   longest_alignment      - TRUE first keeps one alignment per read
#                            (keep_longest_alignment())
#   validate_blocks        - TRUE stops on empty blocks or blocks starting
#                            before the read
# Output:
#   data.frame, one row per block: chr, start, end (0-based, half-open), RID,
#   score, strand, plus block_score and/or read_start, read_end when
#   requested. With drop_sentinels the rows are grouped by RID; an empty
#   data.frame() for an empty input. Example (m6A calls):
#
#      chr start   end                                   RID score strand
#   1 chr5 10368 10369 m84241_260613_040640_s4/177673938/ccs    29      +
#   2 chr5 10674 10675 m84241_260613_040640_s4/177673938/ccs    29      +
#   3 chr5 11296 11297 m84241_260613_040640_s4/177673938/ccs    29      +
# ---------------------------------------------------------------------------
convert_ft_bed12_to_bed6 <- function(bed12_df, include_read_start_end = FALSE,
                                   drop_sentinels = TRUE, keep_block_scores = FALSE,
                                   format = NULL, source = NULL,
                                   longest_alignment = FALSE, validate_blocks = FALSE) {
  # an explicit format fixes the expected column count and the sentinel rule
  if (!is.null(format)) {
    expected <- footprint_format_columns(format)
    if (!format %in% c("bed12_fibertools", "bed13_fiberhmm", "bed15_fiberhmm_tf"))
      stop("Not a BED12 block format: ", format)
    configured_sentinels <- format == "bed12_fibertools"
    if (!missing(drop_sentinels) && !identical(drop_sentinels, configured_sentinels))
      stop("drop_sentinels conflicts with configured format '", format, "'")
    drop_sentinels <- configured_sentinels
    if (ncol(bed12_df) != expected)
      footprint_column_error(if (is.null(source)) "<unknown>" else source,
                             format, ncol(bed12_df))
  }
  if (nrow(bed12_df) == 0) {
    return(data.frame())
  }
  # Existing calls retain the BED12-only contract. FiberHMM callers opt in by
  # keeping all blocks or requesting its per-block scores.
  if (is.null(format)) {
    legacy <- isTRUE(drop_sentinels) && !isTRUE(keep_block_scores)
    if (legacy && ncol(bed12_df) != 12L) stop("Expected 12 BED columns")
    if (!ncol(bed12_df) %in% c(12L, 13L)) stop("Expected 12 or 13 BED columns")
  }
  has_block_scores <- if (is.null(format)) ncol(bed12_df) == 13L else format != "bed12_fibertools"
  colnames(bed12_df) <- c('chr','start','end','RID','score','strand','read_start','read_end','rgb','blockCount','blockSizes','blockStarts',
                         if (has_block_scores) 'blockScores',
                         if (identical(format, "bed15_fiberhmm_tf")) c('blockEdgeLeft', 'blockEdgeRight'))
  if (identical(format, "bed12_fibertools")) {
    bed12_df$blockCount <- as.integer(bed12_df$blockCount)
    bed12_df$blockSizes <- as.character(bed12_df$blockSizes)
    bed12_df$blockStarts <- as.character(bed12_df$blockStarts)
  }
  if (longest_alignment) bed12_df <- keep_longest_alignment(bed12_df)

  # expand BED12 into one row per block
  block_sizes_list  <- strsplit(sub(",$", "", bed12_df$blockSizes),  ",", fixed = TRUE)
  block_starts_list <- strsplit(sub(",$", "", bed12_df$blockStarts), ",", fixed = TRUE)
  n_blocks <- lengths(block_sizes_list)
  offsets <- suppressWarnings(as.integer(unlist(block_starts_list)))
  sizes <- suppressWarnings(as.integer(unlist(block_sizes_list)))
  # blockCount, blockSizes and blockStarts must agree and be integers
  if (anyNA(bed12_df$blockCount) || any(n_blocks != bed12_df$blockCount) ||
      any(lengths(block_starts_list) != n_blocks) || anyNA(offsets) || anyNA(sizes))
    stop("Invalid BED12 blocks")
  # blockStarts are offsets from chromStart; repeat each read's fields per block
  starts <- rep(bed12_df$start, n_blocks) + offsets
  bed6_df <- data.frame(
    chr    = rep(bed12_df$chr,    n_blocks),
    start  = starts,
    end    = starts + sizes,
    RID    = rep(bed12_df$RID,    n_blocks),
    score  = rep(bed12_df$score,  n_blocks),
    strand = rep(bed12_df$strand, n_blocks),
    stringsAsFactors = FALSE
  )

  if (validate_blocks && (any(bed6_df$end <= bed6_df$start) ||
      any(bed6_df$start < rep(bed12_df$start, n_blocks))))
    stop("Invalid blocks in ", if (is.null(source)) "<unknown>" else source,
         " (configured format '", if (is.null(format)) "<unspecified>" else format, "')")

  if (keep_block_scores && has_block_scores) {
    block_scores_list <- strsplit(sub(",$", "", bed12_df$blockScores), ",", fixed = TRUE)
    if (any(lengths(block_scores_list) != n_blocks))
      stop("blockScores length does not match blockStarts")
    bed6_df$block_score <- as.numeric(unlist(block_scores_list, use.names = FALSE))
  }

  # Preserve the original RID grouping and ordering for every existing LCL
  # call. Keeping all blocks instead preserves alignment/block order.
  if (drop_sentinels) {
    bed6_df <- dplyr::ungroup(dplyr::slice(dplyr::group_by(bed6_df, RID),
                                         -c(1, dplyr::n())))
  }

  if (include_read_start_end) {
    if (drop_sentinels) {
      read_start_end_df <- dplyr::select(bed12_df, RID, read_start, read_end)
      bed6_df <- dplyr::left_join(bed6_df, read_start_end_df, by = "RID")
    } else {
      bed6_df$read_start <- rep(bed12_df$read_start, n_blocks)
      bed6_df$read_end <- rep(bed12_df$read_end, n_blocks)
    }
  }
  return(as.data.frame(bed6_df))
}


# ---------------------------------------------------------------------------
# Read the modification calls (m6A or CpG) of every read overlapping one
# genomic region from a fibertools `ft extract` BED12. Uses a tabix query when
# a .tbi index exists, otherwise reads the whole file and filters by overlap.
# Reads with several alignments keep only their longest one; sentinel blocks
# are dropped and every call is reported as a 1-based position.
#
# Inputs:
#   ft_extracted_file       - fibertools BED12 (.bed.gz) of one sample
#   region                  - a single region as GRanges, or anything
#                             as(, "GRanges") accepts (e.g. "chr1:100-200")
#   keep_pos_in_region_only - TRUE keeps only calls inside the region; FALSE
#                             keeps every call on the overlapping reads
#   verbose                 - unused
# Output:
#   data.frame, one row per call, sorted by RID: chrom, pos (1-based call
#   position), read_start (1-based), read_end, RID, score, strand; an empty
#   data.frame() when no read overlaps the region
# ---------------------------------------------------------------------------
extract_ft_region_reads <- function(ft_extracted_file, region, keep_pos_in_region_only = TRUE, verbose = FALSE) {
  
  # convert to GRanges
  if (!inherits(region, "GRanges")) {
    region_gr <- as(region, "GRanges")
  } else {
    region_gr <- region
  }
  if (length(region_gr) != 1) {
    stop("region_gr must contain exactly one region!")
  }
  region_chr <- as.character(seqnames(region_gr))
  region_start <- start(region_gr)
  region_end <- end(region_gr)

  if (!file.exists(ft_extracted_file)) {
    stop("ft_extracted_file does not exist!")
  }

  if (file.exists(paste0(ft_extracted_file, ".tbi"))) {
    extracted_data_bed12 <- read_tabix_region(ft_extracted_file, region_gr)
  } else {
    extracted_data_bed12 <- read_ft_bed12(ft_extracted_file)
    extracted_data_bed12 <- extracted_data_bed12 %>%
      dplyr::filter(chr == region_chr, start <= region_end, end >= region_start)
  }

  if (nrow(extracted_data_bed12) == 0) {
    return(data.frame())
  }
  colnames(extracted_data_bed12) <- c('chr','start','end','RID','score','strand','read_start','read_end','rgb','blockCount','blockSizes','blockStarts')

  # if duplicated RIDs are found, keep the longest alignment per RID
  if (any(duplicated(extracted_data_bed12$RID))) {
    extracted_data_bed12 <- keep_longest_alignment(extracted_data_bed12)
    extracted_data_bed12 <- extracted_data_bed12[
      order(extracted_data_bed12$RID, method = "radix"), , drop = FALSE]
  }

  # one row per call (1-bp block); sentinel blocks dropped
  reads <- convert_ft_bed12_to_bed6(extracted_data_bed12, include_read_start_end = TRUE,
    format = "bed12_fibertools", source = ft_extracted_file)
  reads <- reads %>% dplyr::rename(chrom = chr)
  # use 1-based positions
  reads <- reads %>% dplyr::mutate(pos = end, read_start = read_start + 1)
  reads <- reads %>% dplyr::select(chrom, pos, read_start, read_end, RID, score, strand)

  if (keep_pos_in_region_only) {
    reads <- reads %>% dplyr::filter(pos >= region_start & pos <= region_end)
  }

  reads <- reads %>% dplyr::arrange(RID)
  return(reads)
}

# ---------------------------------------------------------------------------
# Create one summary row per read, either
#   1. from an existing calls table produced by extract_ft_region_reads()
#      (pass `reads`): start / end = min(read_start) / max(read_end) of the
#      read's rows, other fields taken from its first row, or
#   2. directly from the fibertools BED12 file for one genomic region (pass
#      ft_extracted_file and region): one row per overlapping record, start =
#      BED start + 1. Reads with several alignments are not collapsed here.
#
# Inputs:
#   reads             - calls table with RID, chr or chrom, read_start,
#                       read_end, strand (mode 1)
#   ft_extracted_file - fibertools BED12 (.bed.gz) path (mode 2)
#   region            - single region as GRanges or coercible to it (mode 2)
#   keep_columns      - extra columns to carry over (e.g. sample_name, score);
#                       names missing from the input are skipped with a message
# Output:
#   data.frame, one row per read, sorted by RID: RID, chr, start (1-based),
#   end, strand, plus keep_columns; an empty data.frame() when there are no
#   reads
# ---------------------------------------------------------------------------
extract_ft_read_info <- function(reads, ft_extracted_file, region, keep_columns = NULL) {
  # mode 1: summarise an existing calls table
  if (!missing(reads)) {
    if (nrow(reads) == 0) {
      return(data.frame())
    }
    if (!"read_start" %in% colnames(reads) | !"read_end" %in% colnames(reads)) {
      stop("read_start and read_end columns are required in reads!")
    }
    if ("chrom" %in% colnames(reads))
      reads <- reads %>% dplyr::rename(chr = chrom)

    if (!is.null(keep_columns)) {
      if (!all(keep_columns %in% colnames(reads))) {
        cat("Some keep_columns are not found in reads! Only use those in reads.\n")
        keep_columns <- keep_columns[keep_columns %in% colnames(reads)]
      }
      reads <- reads %>% dplyr::select(RID, chr, read_start, read_end, strand, all_of(keep_columns))
    } else {
      reads <- reads %>% dplyr::select(RID, chr, read_start, read_end, strand)
    }

    rids_df <- reads %>%
      dplyr::group_by(RID) %>%
      dplyr::summarise(chr = chr[1], start = min(read_start), end = max(read_end),
                       strand = strand[1], across(all_of(keep_columns), ~ .x[1])) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(RID)
  } else {
    # mode 2: read the overlapping records from the BED12 file
    if (!inherits(region, "GRanges")) {
      region_gr <- as(region, "GRanges")
    } else {
      region_gr <- region
    }
    if (length(region_gr) != 1) {
      stop("region_gr must contain exactly one region!")
    }
    region_chr <- as.character(seqnames(region_gr))
    region_start <- start(region_gr)
    region_end <- end(region_gr)

    if (!file.exists(ft_extracted_file)) {
      stop("ft_extracted_file does not exist!")
    }

    if (file.exists(paste0(ft_extracted_file, ".tbi"))) {
      extracted_data_bed12 <- read_tabix_region(ft_extracted_file, region_gr)
    } else {
      extracted_data_bed12 <- read_ft_bed12(ft_extracted_file)
      extracted_data_bed12 <- extracted_data_bed12 %>%
        dplyr::filter(chr == region_chr, start <= region_end, end >= region_start)
    }

    if (nrow(extracted_data_bed12) == 0) {
      return(data.frame())
    }
    colnames(extracted_data_bed12) <- c('chr','start','end','RID','score','strand','read_start','read_end','rgb','blockCount','blockSizes','blockStarts')

    if (!is.null(keep_columns)) {
      if (!all(keep_columns %in% colnames(extracted_data_bed12))) {
        cat("Columns not included: ", keep_columns[!keep_columns %in% colnames(extracted_data_bed12)], "\n")
        keep_columns <- keep_columns[keep_columns %in% colnames(extracted_data_bed12)]
      }
      rids_df <- extracted_data_bed12 %>% dplyr::select(RID, chr, start, end, strand, all_of(keep_columns))
    } else {
      rids_df <- extracted_data_bed12 %>% dplyr::select(RID, chr, start, end, strand)
    }

    # BED start is 0-based; convert to 1-based
    rids_df <- rids_df %>% dplyr::mutate(start = start + 1)
    rids_df <- rids_df %>% dplyr::arrange(RID)
  }

  return(as.data.frame(rids_df))
}

# ---------------------------------------------------------------------------
# Pooled read x position matrix for one region (m6A by default, or CpG),
# tagging each read with its sample of origin. For every sample the region is
# read from its fibertools BED12 with tabix, each read keeps its longest
# alignment, sentinel blocks are dropped and the calls inside the analysis
# window become the 1 entries. Each sample's file is
#   <fire_dir>/extracted_results/<modality>_by_chr/
#     <sample>.ft_extracted_<modality>.<chr>.bed.gz
# unless m6a_paths is given.
#
# Inputs:
#   sample_table - one row per sample with unique sample_name and fire_dir
#   region       - one row: region_id, chr, start / end (BED, 0-based start)
#                  and analysis_start / analysis_end (1-based inclusive;
#                  analysis_start = start + 1, analysis_end = end)
#   modality     - "m6a" reads m6a_by_chr, "cpg" reads cpg_by_chr
#   full_span    - TRUE keeps reads whose longest alignment spans the whole
#                  analysis window (1 = call, 0 = no call); FALSE keeps every
#                  overlapping read and sets NA where the read does not cover
#                  a position
#   positions    - matrix columns; NULL = positions with at least one call in
#                  the window, or e.g. seq(analysis_start, analysis_end) for
#                  every base
#   matrix_dir   - optional cache folder; the result is saved to
#                  <matrix_dir>/<region_id>/<region_id>_<modality>_matrix.rds
#                  with per-sample read counts in <region_id>_read_qc.tsv
#                  (<region_id>_cpg_read_qc.tsv for CpG)
#   reuse        - load the cached result when its signature (region,
#                  samples, paths and options) is identical
#   m6a_paths    - optional BED12 paths, one per sample in sample_table order,
#                  replacing the default path (for either modality)
# Output:
#   list(met_mat = reads x positions dgCMatrix, rownames = RID, colnames =
#   position; rids_df = one row per read: RID ("<sample>::<read>"),
#   original_RID, chr, start (1-based), end, strand, sample_name, score;
#   qc = per sample overlapping_reads and full_span_reads; signature = cache
#   key). Stops with fewer than 3 reads or no calls in the window
# ---------------------------------------------------------------------------
assemble_region_m6a <- function(sample_table, region, modality = c("m6a", "cpg"),
                                full_span = FALSE, positions = NULL,
                                matrix_dir = NULL, reuse = TRUE, m6a_paths = NULL) {
  modality <- match.arg(modality)
  label <- c(m6a = "m6A", cpg = "CpG")[[modality]]
  stopifnot(nrow(region) == 1L,
    all(c("sample_name", "fire_dir") %in% names(sample_table)),
    !anyDuplicated(sample_table$sample_name),
    region$analysis_start == region$start + 1L, region$analysis_end == region$end,
    region$analysis_start <= region$analysis_end)
  if (!is.null(positions))
    stopifnot(!anyNA(positions), !anyDuplicated(positions),
      all(positions >= region$analysis_start & positions <= region$analysis_end))
  if (is.null(m6a_paths)) {
    paths <- file.path(sample_table$fire_dir, "extracted_results", paste0(modality, "_by_chr"),
      paste0(sample_table$sample_name, ".ft_extracted_", modality, ".", region$chr, ".bed.gz"))
  } else {
    if (!is.character(m6a_paths) || length(m6a_paths) != nrow(sample_table) ||
        anyNA(m6a_paths) || any(!nzchar(m6a_paths)))
      stop("m6a_paths must contain one nonempty path per sample, in sample_table order")
    paths <- unname(m6a_paths)
  }
  required <- c(paths, paste0(paths, ".tbi"))
  if (!all(file.exists(required))) stop("Missing input: ", paste(required[!file.exists(required)], collapse = ", "))
  signature <- list(version = 1L, region = region, samples = sample_table)
  if (!is.null(m6a_paths)) signature$m6a_paths <- paths
  if (modality != "m6a") signature$modality <- modality
  if (!full_span) signature$full_span <- FALSE
  if (!is.null(positions)) signature$positions <- positions
  cache <- NULL
  if (!is.null(matrix_dir)) {
    region_matrix_dir <- file.path(matrix_dir, region$region_id)
    dir.create(region_matrix_dir, recursive = TRUE, showWarnings = FALSE)
    cache <- file.path(region_matrix_dir, paste0(region$region_id, "_", modality, "_matrix.rds"))
    if (reuse && file.exists(cache)) {
      previous <- readRDS(cache)
      if (identical(previous$signature[names(signature)], signature)) {
        stopifnot(!anyDuplicated(previous$rids_df$RID),
          !full_span || all(previous$rids_df$start <= region$analysis_start),
          !full_span || all(previous$rids_df$end >= region$analysis_end),
          identical(rownames(previous$met_mat), as.character(previous$rids_df$RID)))
        return(previous)
      }
    }
  }
  region_gr <- GenomicRanges::GRanges(region$chr,
    IRanges::IRanges(region$analysis_start, region$analysis_end))
  reads_list <- metadata_list <- qc_list <- vector("list", nrow(sample_table))
  for (sample_index in seq_len(nrow(sample_table))) {
    sample_name <- sample_table$sample_name[sample_index]
    bed <- read_ft_bed12(paths[sample_index], region_gr, longest_alignment = TRUE)
    overlapping <- nrow(bed)
    spanning <- if (overlapping) bed$start <= region$start & bed$end >= region$end else logical()
    if (full_span) bed <- bed[spanning, , drop = FALSE]
    qc_list[[sample_index]] <- data.frame(sample_name = sample_name,
      overlapping_reads = overlapping, full_span_reads = sum(spanning))
    if (!nrow(bed)) next
    original_ids <- as.character(bed$RID)
    bed$RID <- paste(sample_name, original_ids, sep = "::")
    metadata_list[[sample_index]] <- data.frame(
      RID = bed$RID, original_RID = original_ids, chr = bed$chr,
      start = bed$start + 1L, end = bed$end, strand = bed$strand,
      sample_name = sample_name, score = bed$score)
    blocks <- convert_ft_bed12_to_bed6(bed,
      format = "bed12_fibertools", source = paths[sample_index])
    if (any(blocks$end - blocks$start != 1L)) stop("Non-single-base ", label, " block in ", paths[sample_index])
    blocks <- blocks[blocks$end >= region$analysis_start & blocks$end <= region$analysis_end, , drop = FALSE]
    reads_list[[sample_index]] <- data.frame(RID = as.character(blocks$RID), pos = blocks$end)
  }
  rids_df <- dplyr::bind_rows(metadata_list)
  reads <- dplyr::bind_rows(reads_list)
  if (nrow(rids_df) < 3L || !nrow(reads))
    stop("Insufficient ", if (full_span) "full-span ", "reads/", label, " sites: ", region$region_id)

  # read x position matrix: 1 = call, 0 = covered without a call
  if (is.null(positions)) positions <- sort(unique(reads$pos))
  pairs <- unique(data.frame(read_index = match(reads$RID, rids_df$RID),
    site_index = match(reads$pos, positions)))
  pairs <- pairs[!is.na(pairs$site_index), , drop = FALSE]
  met_mat <- Matrix::sparseMatrix(i = pairs$read_index, j = pairs$site_index,
    x = rep(1, nrow(pairs)), dims = c(nrow(rids_df), length(positions)),
    dimnames = list(rids_df$RID, as.character(positions)))
  if (!full_span) {
    # NA = the read does not cover that position
    met_mat <- as.matrix(met_mat)
    met_mat[!(outer(rids_df$start, positions, "<=") & outer(rids_df$end, positions, ">="))] <- NA
    met_mat <- as(met_mat, "dgCMatrix")
  }
  stopifnot(!full_span || !anyNA(met_mat), !anyDuplicated(rownames(met_mat)))
  result <- list(met_mat = met_mat, rids_df = rids_df, qc = dplyr::bind_rows(qc_list),
    signature = signature)
  if (!is.null(cache)) {
    saveRDS(result, cache)
    data.table::fwrite(result$qc, file.path(region_matrix_dir,
      paste0(region$region_id, if (modality == "cpg") "_cpg", "_read_qc.tsv")), sep = "\t")
  }
  result
}

# ---------------------------------------------------------------------------
# Read one modification modality for each sample, preserving sample order.
# Each sample's per-chromosome fibertools file is
#   <ft_result_dir>/<sample>/extracted_results/<modality>_by_chr/
#     <sample>.ft_extracted_<modality>.<chr>.bed.gz
# and is parsed with extract_ft_region_reads().
#
# Inputs:
#   sample_names            - sample folder names under ft_result_dir
#   region_gr               - GRanges of one region
#   ft_result_dir           - root of the per-sample fibertools results
#   modality                - "m6a" or "cpg"
#   keep_pos_in_region_only - passed to extract_ft_region_reads()
#   verbose                 - passed to extract_ft_region_reads()
# Output:
#   list with one element per sample: that sample's calls table with a
#   leading sample_name column, or NULL when the sample has no reads in the
#   region. Callers drop the NULLs before combining the tables.
# ---------------------------------------------------------------------------
read_sample_region_reads <- function(sample_names, region_gr, ft_result_dir,
                                     modality = c("m6a", "cpg"),
                                     keep_pos_in_region_only, verbose) {
  modality <- match.arg(modality)
  region_chr <- as.character(seqnames(region_gr))
  lapply(sample_names, function(sample_name) {
    extracted_dir <- file.path(ft_result_dir, sample_name, "extracted_results",
                               paste0(modality, "_by_chr"))
    extracted_file <- file.path(extracted_dir,
      paste0(sample_name, ".ft_extracted_", modality, ".", region_chr, ".bed.gz"))
    sample_reads <- extract_ft_region_reads(extracted_file, region_gr,
      keep_pos_in_region_only = keep_pos_in_region_only, verbose = verbose)
    if (nrow(sample_reads) == 0) return(NULL)
    dplyr::mutate(sample_reads, sample_name = sample_name, .before = 1)
  })
}


# ===========================================================================
# Read-level region plots: per-region extraction and loading
#
# Read the per-sample folders written by
# extract_region_result_macrophage.sh or extract_region_result_lcl.sh (both in
# code/parsing_functions/, run through extract_region_results()):
#   <out_root>/<outname>/<sample>/parsed/
# The figure is drawn by plot_region_panels() / plot_region_example() in
# plotting_functions.r. Package calls are qualified, so sourcing attaches nothing.
# ===========================================================================

# Column names of a fibertools `ft extract` BED12 and of a FIRE `ft fire --extract` BED
FT_BED12_COLS <- c("chr", "start", "end", "RID", "score", "strand",
                   "read_start", "read_end", "rgb", "blockCount",
                   "blockSizes", "blockStarts")
FIRE_BED_COLS <- c("chr", "start", "end", "RID", "score", "strand",
                   "thickStart", "thickEnd", "rgb", "fire_score", "tag")

# FIRE colours the per-read segmentation by class; everything that is not grey
# (nucleosome) or purple (linker) is a FIRE element.
RGB_NUCLEOSOME <- "169,169,169"
RGB_LINKER     <- "147,112,219"


# ---------------------------------------------------------------------------
# Read a headerless BED-like file, naming its leading columns.
#
# Inputs:
#   file      - path (plain or .gz); must exist
#   col_names - optional names for the first length(col_names) columns
# Output:
#   data.frame; an empty data.frame() when the file has size 0. Stops when the
#   file is missing or has fewer columns than col_names
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Read `ft extract --m6a` / `--cpg` output for a region into one row per
# modified base. A read with several alignments keeps its longest; sentinel
# blocks are dropped.
#
# Inputs:
#   file                     - fibertools BED12 (e.g. parsed/extracted.m6a.bed.gz)
#   region_start, region_end - 1-based inclusive window; calls outside are dropped
# Output:
#   data.frame chr, RID, strand, read_start (1-based), read_end, pos (1-based
#   position of the modified base); an empty data.frame() for an empty file
# ---------------------------------------------------------------------------
read_ft_mod_region <- function(file, region_start, region_end) {
  df <- read_bed(file, FT_BED12_COLS)
  if (nrow(df) == 0) return(data.frame())

  # one alignment per read: keep the longest if a read is split
  if (any(duplicated(df$RID))) {
    df <- keep_longest_alignment(df)
    df <- df[order(df$RID, method = "radix"), , drop = FALSE]
  }

  blocks <- convert_ft_bed12_to_bed6(df, include_read_start_end = TRUE,
                                    format = "bed12_fibertools", source = file)
  if (nrow(blocks) == 0) return(data.frame())
  blocks <- blocks[order(match(blocks$RID, df$RID)), , drop = FALSE]

  # blocks are 0-based half-open and 1 bp wide for a modification call, so the
  # block end is the 1-based coordinate of the modified base
  out <- data.frame(
    chr        = blocks$chr,
    RID        = as.character(blocks$RID),
    strand     = blocks$strand,
    read_start = blocks$read_start + 1L,
    read_end   = blocks$read_end,
    pos        = blocks$end,
    stringsAsFactors = FALSE)

  out[out$pos >= region_start & out$pos <= region_end, , drop = FALSE]
}


# ---------------------------------------------------------------------------
# Read `ft fire --extract` output (parsed/fire.bed) and label each segment as
# nucleosome (grey), linker (purple) or FIRE (any other colour). FIRE scores
# above 1 are capped at 1.
#
# Inputs:
#   file - FIRE per-read segmentation BED (11 columns, FIRE_BED_COLS)
# Output:
#   data.frame with FIRE_BED_COLS plus class (factor nucleosome / linker /
#   FIRE); stops when the file is empty
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Read the FIRE peak calls for the region. The v0.1 peak BED has 29 columns;
# only the peak interval and the FDR are kept. logFDR is stored x10, as in the
# FIRE track hub.
#
# Inputs:
#   file - parsed/fire_peaks.bed (may be empty)
# Output:
#   data.frame chr, start, end, FDR, logFDR; NULL when the file is empty.
#   Stops when the file does not have 29 columns
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Bin footprint sizes into the classes used for colouring. Bins are closed
# [lo, hi]; a size on a shared break goes to the lower bin. Sizes outside
# [min(breaks), max(breaks)] get NA and are dropped when plotting - this is
# what keeps nucleosome-sized calls out of the footprint track.
#
# Inputs:
#   size   - footprint sizes (bp)
#   breaks - at least two size breaks, e.g. c(10, 30, 60, 80)
# Output:
#   factor with levels "<lo>-<hi> bp", one value per size (NA outside)
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Read FiberHMM footprint calls for the region and expand them into one row
# per footprint (FiberHMM has no sentinel blocks; the TF edge columns of the
# 15-column format are ignored).
#
# Inputs:
#   file        - parsed/region.fiberhmm_tf.bed or region.fiberhmm_footprint.bed
#   format      - "bed15_fiberhmm_tf" for `tf`, "bed13_fiberhmm" for
#                 `footprint` / `msp`
#   min_score   - minimum per-block score (0 keeps every call)
#   size_breaks - size classes, see assign_size_class()
# Output:
#   data.frame chr, start, end, RID, size, score, read_start, read_end, class;
#   an empty data.frame() for an empty file
# ---------------------------------------------------------------------------
read_fiberhmm_region <- function(file, format,
                                 min_score = 0,
                                 size_breaks = c(10, 30, 60, 80)) {
  if (length(format) != 1L || is.na(format) ||
      !format %in% c("bed13_fiberhmm", "bed15_fiberhmm_tf"))
    stop("FiberHMM requires bed13_fiberhmm or bed15_fiberhmm_tf format")
  df <- read_bed(file)
  if (nrow(df) == 0) return(data.frame())

  # FiberHMM has no sentinels; the configured TF edge columns are ignored.
  blocks <- convert_ft_bed12_to_bed6(df, format = format, source = file,
                                    include_read_start_end = TRUE,
                                    keep_block_scores = TRUE)
  if (nrow(blocks) == 0) return(data.frame())

  if ("block_score" %in% colnames(blocks) && min_score > 0)
    blocks <- blocks[blocks$block_score >= min_score, , drop = FALSE]

  data.frame(chr        = blocks$chr,
             start      = blocks$start,
             end        = blocks$end,
             RID        = as.character(blocks$RID),
             size       = as.integer(blocks$end - blocks$start),
             score      = if ("block_score" %in% colnames(blocks)) blocks$block_score else NA_real_,
             read_start = blocks$read_start,
             read_end   = blocks$read_end,
             class      = assign_size_class(blocks$end - blocks$start, size_breaks),
             stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------------
# Read fibertools nucleosome calls (`ft extract --nuc` BED12, the LCL
# parsed/extracted.nuc.bed.gz) for the region, one row per nucleosome. A read
# with several alignments keeps its longest; sentinel blocks are dropped.
#
# Inputs:
#   file        - fibertools nucleosome BED12 (plain or .gz)
#   size_breaks - nucleosome size window, see assign_size_class()
# Output:
#   data.frame with the columns of read_fiberhmm_region() (score is NA);
#   an empty data.frame() for an empty file
# ---------------------------------------------------------------------------
read_ft_nuc_region <- function(file, size_breaks = c(130, 160)) {
  df <- read_bed(file, FT_BED12_COLS)
  if (nrow(df) == 0) return(data.frame())
  if (any(duplicated(df$RID))) {
    df <- keep_longest_alignment(df)
    df <- df[order(df$RID, method = "radix"), , drop = FALSE]
  }
  blocks <- convert_ft_bed12_to_bed6(df, include_read_start_end = TRUE,
                                    format = "bed12_fibertools", source = file)
  if (nrow(blocks) == 0) return(data.frame())
  data.frame(chr        = blocks$chr,
             start      = blocks$start,
             end        = blocks$end,
             RID        = as.character(blocks$RID),
             size       = as.integer(blocks$end - blocks$start),
             score      = NA_real_,
             read_start = blocks$read_start,
             read_end   = blocks$read_end,
             class      = assign_size_class(blocks$end - blocks$start, size_breaks),
             stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------------
# Keep only footprints that fall entirely inside a FIRE element of the same
# read. The read ID is used as the GRanges seqname, so overlaps are confined to
# the read a footprint came from.
#
# Inputs:
#   fps  - footprints (read_fiberhmm_region() rows with RID, start, end)
#   fire - read_fire_region() rows; only class == "FIRE" is used
# Output:
#   the subset of fps, in its original order
# ---------------------------------------------------------------------------
subset_footprints_in_fire <- function(fps, fire) {
  if (nrow(fps) == 0) return(fps)
  fire_el <- fire[fire$class == "FIRE", , drop = FALSE]
  if (nrow(fire_el) == 0) return(fps[0, , drop = FALSE])

  fps.gr  <- GenomicRanges::GRanges(seqnames = fps$RID,
                                    ranges = IRanges::IRanges(start = fps$start + 1L, end = fps$end))
  fire.gr <- GenomicRanges::GRanges(seqnames = fire_el$RID,
                                    ranges = IRanges::IRanges(start = fire_el$start + 1L, end = fire_el$end))
  hits <- IRanges::findOverlaps(fps.gr, fire.gr, type = "within")
  fps[sort(unique(S4Vectors::queryHits(hits))), , drop = FALSE]
}


# ---------------------------------------------------------------------------
# One row per read from its FIRE segmentation: alignment span, strand,
# haplotype tag and arrow direction.
#
# Inputs:
#   fire - read_fire_region() rows with a sample_name column
# Output:
#   data.frame RID, sample_name, chr, start (min), end (max), strand, tag,
#   arrow_end ("last" for + strand, "first" otherwise)
# ---------------------------------------------------------------------------
build_rids_df <- function(fire) {
  rids <- dplyr::group_by(fire, RID)
  rids <- dplyr::summarise(rids, sample_name = unique(sample_name),
                           chr    = unique(chr),
                           start  = min(start),
                           end    = max(end),
                           strand = unique(strand)[1],
                           tag    = unique(tag)[1],
                           .groups = "drop")
  rids <- dplyr::mutate(rids, arrow_end = ifelse(strand == "+", "last", "first"))
  as.data.frame(rids)
}


# ---------------------------------------------------------------------------
# Read x position modification matrix. `ft extract` reports the read's full
# alignment span, so "covered but unmodified" is a real observation.
#
# Inputs:
#   reads     - calls with RID and pos (read_ft_mod_region())
#   rids_df   - one row per read with RID, start, end (build_rids_df())
#   positions - matrix columns (1-based positions)
# Output:
#   matrix reads x positions: 1 = modified, 0 = covered but unmodified,
#   NA = not covered by the read
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Sliding-window mean over positions (window of +/- n/2 bp).
#
# Inputs:
#   x   - values, one per position
#   pos - positions of x
#   n   - window width in bp
# Output:
#   numeric vector, the window mean at each position
# ---------------------------------------------------------------------------
smooth_mean_slidewindow <- function(x, pos, n = 20) {
  sapply(pos, function(i) mean(x[which(pos >= (i - n / 2) & pos <= (i + n / 2))]))
}


# ---------------------------------------------------------------------------
# Modification pileup over the region, optionally split into groups of reads.
# The denominator is the number of reads in the group, not the number covering
# each position, so a position a read does not cover counts as unmethylated.
#
# Inputs:
#   reads, rids_df           - calls and per-read table (load_region_results())
#   window_start, window_end - 1-based inclusive window
#   split_by                 - optional rids_df column defining the groups
#   smooth                   - TRUE adds the sliding-window mean
#   window_n                 - smoothing window (bp)
# Output:
#   data.frame pos, base, group, cov, met, frac, smooth_cov, smooth_frac;
#   stops when no modified position falls inside the window
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Pairwise Euclidean distance over the positions where both rows are non-NA.
# Pairs sharing no non-NA position keep the initial distance of 100.
#
# Inputs:
#   x - matrix, rows = reads (rownames = RIDs)
# Output:
#   dist object
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Hierarchically cluster reads on their footprint calls inside the region.
# Footprints are intersected with the region and expanded to one row per base
# pair; the read x position matrix holds the footprint size class as an
# integer, NA where a read has no footprint, so two reads are compared only
# over positions where both carry a size-classed footprint.
#
# Inputs:
#   fps                                  - footprints with chr, start, end, RID, class
#   region_chr, region_start, region_end - the region (1-based start)
#   size_levels                          - class levels, in size order
# Output:
#   hclust object, or NULL when fewer than two reads have a footprint in the region
# ---------------------------------------------------------------------------
hclust_reads_by_footprints <- function(fps, region_chr, region_start, region_end,
                                       size_levels) {
  if (is.null(fps) || nrow(fps) == 0) return(NULL)

  fps.gr <- GenomicRanges::GRanges(seqnames = fps$chr,
                                   ranges = IRanges::IRanges(start = fps$start + 1L, end = fps$end),
                                   RID = as.character(fps$RID),
                                   class = factor(as.character(fps$class), levels = size_levels))
  fps.gr <- BiocGenerics::sort(fps.gr)

  region.gr <- GenomicRanges::GRanges(seqnames = region_chr,
                                      ranges = IRanges::IRanges(start = region_start, end = region_end))

  hits <- IRanges::findOverlaps(fps.gr, region.gr)
  if (length(hits) == 0) return(NULL)
  inter <- IRanges::pintersect(fps.gr[S4Vectors::queryHits(hits)],
                               region.gr[S4Vectors::subjectHits(hits)])

  # one row per base pair covered by a footprint
  n_bp <- BiocGenerics::width(inter)
  df <- data.frame(
    RID = rep(S4Vectors::mcols(inter)$RID, n_bp),
    pos = sequence(n_bp, from = BiocGenerics::start(inter)),
    class_value = rep(as.numeric(S4Vectors::mcols(inter)$class), n_bp),
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


# ---------------------------------------------------------------------------
# Read order for plotting: reads without footprints first (in their given
# order), then the clustered reads in hclust order.
#
# Inputs:
#   fps                                  - footprints, see hclust_reads_by_footprints()
#   rid_levels                           - every read ID, in the fallback order
#   region_chr, region_start, region_end - the region
#   size_levels                          - footprint class levels
# Output:
#   character vector of read IDs
# ---------------------------------------------------------------------------
order_reads_by_footprints <- function(fps, rid_levels, region_chr,
                                      region_start, region_end, size_levels) {
  hc <- hclust_reads_by_footprints(fps, region_chr, region_start, region_end, size_levels)
  if (is.null(hc)) return(rid_levels)
  clustered <- hc$labels[hc$order]
  clustered <- clustered[clustered %in% rid_levels]
  c(setdiff(rid_levels, clustered), clustered)
}


# ---------------------------------------------------------------------------
# Run a per-region extraction script once per sample and return the parsed/
# folders that load_region_results() reads. A sample is skipped when its
# parsed/fire.bed exists, unless regenerate = TRUE or <sample>/region.txt
# (written after each extraction) records a different region. Script output
# goes to <out_root>/<outname>/<sample>/extract.log.
#
# Inputs:
#   sample_names - sample labels, the script's first argument
#   region       - one row with chr, start, end, passed as "chr:start-end" (the
#                  same numbers load_region_results() is given)
#   outname      - region folder name under out_root
#   out_root     - output root; results go to <out_root>/<outname>/<sample>/parsed
#   script       - extract_region_result_macrophage.sh or extract_region_result_lcl.sh
#   regenerate   - TRUE re-runs the script even when the outputs exist
#   bash         - bash executable
# Output:
#   character vector of parsed/ folders, in sample_names order. Stops when the
#   script exits non-zero, naming the log
# ---------------------------------------------------------------------------
extract_region_results <- function(sample_names, region, outname, out_root, script,
                                   regenerate = FALSE, bash = "bash") {
  stopifnot(is.data.frame(region), nrow(region) == 1L, file.exists(script))
  region_str <- sprintf("%s:%d-%d", region$chr, as.integer(region$start), as.integer(region$end))
  vapply(sample_names, function(s) {
    sample_dir <- file.path(out_root, outname, s)
    parsed <- file.path(sample_dir, "parsed")
    stamp <- file.path(sample_dir, "region.txt")
    same_region <- !file.exists(stamp) || identical(readLines(stamp, warn = FALSE), region_str)
    if (!regenerate && same_region && file.exists(file.path(parsed, "fire.bed"))) return(parsed)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    log <- file.path(sample_dir, "extract.log")
    rc <- system2(bash, shQuote(c(script, s, region_str, outname, out_root)),
                  stdout = log, stderr = log)
    if (rc != 0) stop("extraction failed for ", outname, " ", s, " (exit ", rc, "); see ", log)
    writeLines(region_str, stamp)
    parsed
  }, character(1), USE.NAMES = FALSE)
}


# ---------------------------------------------------------------------------
# Load one region's extracted results across samples into the tables that
# plot_region_panels() draws. Reads present in both the m6A slice and fire.bed
# are kept and ordered by footprint similarity (order_reads_by_footprints()).
#
# Inputs:
#   region                 - one row with chr, start, end (1-based inclusive window)
#   sample_names           - sample labels, e.g. c("LPS_0", "LPS_5", ...)
#   result_dirs            - the matching parsed/ folders, same order
#   include_cpg            - also load CpG calls (bars in the methylation panel,
#                            red tiles in the read panel)
#   fiberHMM_feature       - FiberHMM file for the TF footprints ("tf")
#   min_fiberHMM_fp_score  - minimum FiberHMM per-block score (`tf` calls are
#                            already floored at 50)
#   fiberHMM_size_breaks   - footprint size bins; calls outside get no class
#   nucleosome_feature     - FiberHMM file for nucleosomes ("footprint"), used
#                            when nucleosome_source = "fiberhmm"
#   nucleosome_size_breaks - nucleosome size window
#   nucleosome_label       - legend label of the nucleosome class
#   nucleosome_source      - "fiberhmm" (region.fiberhmm_<feature>.bed, macrophage)
#                            or "ft" (extracted.nuc.bed.gz from `ft extract --nuc`,
#                            LCL). With "ft" the nucleosome rows of fire.bed (an
#                            `ft fire --extract --all` slice) are dropped after the
#                            read spans are taken, so nucleosomes are drawn once
#   window_n, smooth_pileup - pileup smoothing, see pileup_reads()
#   verbose                - print progress and counts
# Output:
#   list region, sample_names, reads, rids_df, fire, fire_peaks, fps,
#   fps_infire, nucs, pileup, size_levels, nuc_label, window_n, smooth_pileup,
#   group_col (NULL)
# ---------------------------------------------------------------------------
load_region_results <- function(region,
                                sample_names,
                                result_dirs,
                                include_cpg = FALSE,
                                fiberHMM_feature = "tf",
                                min_fiberHMM_fp_score = 50,
                                fiberHMM_size_breaks = c(10, 30, 60, 80),
                                nucleosome_feature = "footprint",
                                nucleosome_size_breaks = c(130, 160),
                                nucleosome_label = "130-160 bp nucleosome",
                                nucleosome_source = c("fiberhmm", "ft"),
                                window_n = 10,
                                smooth_pileup = TRUE,
                                verbose = TRUE) {

  stopifnot(is.data.frame(region), nrow(region) == 1)
  stopifnot(length(sample_names) == length(result_dirs))
  nucleosome_source <- match.arg(nucleosome_source)

  region <- list(chr = as.character(region$chr),
                 start = as.numeric(region$start),
                 end = as.numeric(region$end))

  fiberhmm_formats <- c(tf = "bed15_fiberhmm_tf", footprint = "bed13_fiberhmm",
                        msp = "bed13_fiberhmm")
  for (feature in list(fiberHMM_feature, nucleosome_feature)) {
    if (length(feature) != 1L || is.na(feature) || !feature %in% names(fiberhmm_formats))
      stop("Unsupported FiberHMM feature: ", paste(feature, collapse = ", "))
  }

  reads_l <- list(); fire_l <- list(); peaks_l <- list(); fps_l <- list(); nuc_l <- list()

  for (i in seq_along(sample_names)) {
    s <- sample_names[i]
    d <- result_dirs[i]
    if (!dir.exists(d)) stop("result dir not found: ", d)
    if (verbose) cat("Loading", s, "from", d, "\n")

    mods <- read_ft_mod_region(file.path(d, "extracted.m6a.bed.gz"),
                               region$start, region$end)
    if (nrow(mods) == 0) stop("no m6A calls in the region for ", s)
    mods$base <- "A"
    if (include_cpg) {
      cg <- read_ft_mod_region(file.path(d, "extracted.cpg.bed.gz"),
                               region$start, region$end)
      if (nrow(cg) > 0) { cg$base <- "CG"; mods <- rbind(mods, cg) }
    }
    mods$sample_name <- s
    reads_l[[i]] <- mods

    fr <- read_fire_region(file.path(d, "fire.bed"))
    fr$sample_name <- s
    fire_l[[i]] <- fr

    peaks_l[[i]] <- read_fire_peaks_region(file.path(d, "fire_peaks.bed"))

    fp <- read_fiberhmm_region(
      file.path(d, paste0("region.fiberhmm_", fiberHMM_feature, ".bed")),
      format = unname(fiberhmm_formats[[fiberHMM_feature]]),
      min_score = min_fiberHMM_fp_score,
      size_breaks = fiberHMM_size_breaks)
    if (nrow(fp) > 0) fp$sample_name <- s
    fps_l[[i]] <- fp

    # Nucleosomes. The FiberHMM `footprint` bed carries no per-block score (all
    # zero), so no score filter here; ft nucleosomes have no score at all.
    nuc <- if (nucleosome_source == "ft") {
      read_ft_nuc_region(file.path(d, "extracted.nuc.bed.gz"),
                         size_breaks = nucleosome_size_breaks)
    } else {
      read_fiberhmm_region(
        file.path(d, paste0("region.fiberhmm_", nucleosome_feature, ".bed")),
        format = unname(fiberhmm_formats[[nucleosome_feature]]),
        min_score = 0,
        size_breaks = nucleosome_size_breaks)
    }
    if (nrow(nuc) > 0) {
      nuc <- nuc[!is.na(nuc$class), , drop = FALSE]
      nuc$class <- factor(nucleosome_label, levels = nucleosome_label)
      nuc$sample_name <- s
    }
    nuc_l[[i]] <- nuc
  }

  reads <- do.call(rbind, reads_l)
  fire  <- do.call(rbind, fire_l)
  fps   <- do.call(rbind, Filter(function(x) nrow(x) > 0, fps_l))
  if (is.null(fps)) fps <- data.frame()
  nucs  <- do.call(rbind, Filter(function(x) nrow(x) > 0, nuc_l))
  if (is.null(nucs)) nucs <- data.frame()
  fire_peaks <- do.call(rbind, Filter(Negate(is.null), peaks_l))

  # A read name must not appear in two samples, or per-read joins become ambiguous.
  rid_sample <- unique(fire[, c("RID", "sample_name")])
  if (any(duplicated(rid_sample$RID)))
    stop(sum(duplicated(rid_sample$RID)), " read IDs appear in more than one sample")

  # keep reads that have both a modification call and a FIRE segmentation
  keep <- intersect(unique(reads$RID), unique(fire$RID))
  if (length(keep) == 0) stop("no reads shared between the m6A and FIRE data")
  reads <- reads[reads$RID %in% keep, , drop = FALSE]
  fire  <- fire[fire$RID %in% keep, , drop = FALSE]
  if (nrow(fps) > 0)  fps  <- fps[fps$RID %in% keep, , drop = FALSE]
  if (nrow(nucs) > 0) nucs <- nucs[nucs$RID %in% keep, , drop = FALSE]

  rids_df <- build_rids_df(fire)
  if (verbose)
    cat(nrow(rids_df), "reads in the region across", length(sample_names), "samples\n")

  # an `ft fire --extract --all` fire.bed also carries the ft nucleosomes, which
  # are already in nucs; keep them only for the read spans above
  if (nucleosome_source == "ft")
    fire <- fire[fire$class != "nucleosome", , drop = FALSE]

  size_levels <- levels(assign_size_class(numeric(0), fiberHMM_size_breaks))

  # Order reads by footprint similarity so they stay sorted within each facet. Reads
  # with no footprint in the region are not clustered and keep their initial
  # (alphabetical) order at the top. The pairwise distance is an O(n^2) loop, so this
  # is the slow step - roughly a minute at a few hundred reads.
  rid_levels <- order_reads_by_footprints(fps,
                                          sort(unique(as.character(rids_df$RID))),
                                          region$chr, region$start, region$end,
                                          size_levels)

  # put a table's rows in the footprint read order, RID as a factor
  as_ordered <- function(df) {
    df$RID <- factor(as.character(df$RID), levels = rid_levels)
    df[order(df$RID), , drop = FALSE]
  }
  reads   <- as_ordered(reads)
  fire    <- as_ordered(fire)
  rids_df <- as_ordered(rids_df)
  if (nrow(fps) > 0)  fps  <- as_ordered(fps)
  if (nrow(nucs) > 0) nucs <- as_ordered(nucs)

  reads$base <- factor(reads$base, levels = c("A", "CG"))

  # Only the TF footprints are FIRE-polished. Nucleosomes are protected DNA and FIRE
  # elements are accessible patches, so intersecting the two would discard them all.
  fps_infire <- if (nrow(fps) > 0) subset_footprints_in_fire(fps, fire) else fps
  if (verbose) {
    cat(nrow(fps), "TF footprints,", nrow(fps_infire), "of them inside a FIRE element\n")
    cat(nrow(nucs), "nucleosome calls in",
        paste0(nucleosome_size_breaks, collapse = "-"), "bp\n")
  }

  pileup <- pileup_reads(reads, rids_df, region$start, region$end,
                         split_by = NULL, smooth = smooth_pileup, window_n = window_n)

  list(region = region,
       sample_names = sample_names,
       reads = reads,
       rids_df = rids_df,
       fire = fire,
       fire_peaks = fire_peaks,
       fps = fps,
       fps_infire = fps_infire,
       nucs = nucs,
       pileup = pileup,
       size_levels = size_levels,
       nuc_label = nucleosome_label,
       window_n = window_n,
       smooth_pileup = smooth_pileup,
       group_col = NULL)
}


# ---------------------------------------------------------------------------
# Attach a per-read grouping to a region result and recompute the pileup per
# group, so the methylation panel is faceted like the read panels.
#
# Inputs:
#   res             - load_region_results() result
#   labels          - named character vector, read ID -> group label
#   group_col       - name of the attached column
#   group_levels    - facet order; defaults to the sorted unique labels
#   drop_unassigned - drop reads absent from labels
#   verbose         - print counts and the group x sample table
# Output:
#   res with group_col added to reads, rids_df, fire, fps, fps_infire, nucs,
#   a per-group pileup and res$group_col set
# ---------------------------------------------------------------------------
add_read_groups <- function(res, labels, group_col, group_levels = NULL,
                            drop_unassigned = TRUE, verbose = TRUE) {

  if (is.null(names(labels)))
    stop("labels must be a named vector of read ID -> group")
  labels <- labels[!is.na(labels)]
  if (is.null(group_levels))
    group_levels <- sort(unique(unname(labels)))

  keep <- levels(res$rids_df$RID)            # preserves the footprint ordering
  n_before <- length(keep)
  if (drop_unassigned) keep <- keep[keep %in% names(labels)]
  if (length(keep) == 0)
    stop("no reads left after grouping by ", group_col)

  if (verbose)
    cat(n_before, "reads in region result;", length(keep), "grouped;",
        n_before - length(keep), "dropped.\n")

  # restrict one table to the kept reads and add the group column
  attach_group <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)
    df <- df[as.character(df$RID) %in% keep, , drop = FALSE]
    df$RID <- factor(as.character(df$RID), levels = keep)
    df[[group_col]] <- factor(unname(labels[as.character(df$RID)]), levels = group_levels)
    df
  }
  for (nm in c("reads", "rids_df", "fire", "fps", "fps_infire", "nucs"))
    res[[nm]] <- attach_group(res[[nm]])

  # the pileup from load_region_results() is ungrouped; recompute it per group or the
  # methylation panel would show one pooled curve while the read panels are faceted
  res$pileup <- pileup_reads(res$reads, res$rids_df,
                             res$region$start, res$region$end,
                             split_by = group_col,
                             smooth = res$smooth_pileup,
                             window_n = res$window_n)
  names(res$pileup)[names(res$pileup) == "group"] <- group_col
  res$pileup[[group_col]] <- factor(res$pileup[[group_col]], levels = group_levels)

  if (verbose) print(table(res$rids_df[[group_col]], res$rids_df$sample_name))

  res$group_col <- group_col
  res
}


# ---------------------------------------------------------------------------
# Join per-read topic-model clusters onto a region result, from the
# read_topic_assignments_<outname>.tsv written by the topic-model pipeline.
#
# Inputs:
#   res              - load_region_results() result
#   assignments_file - path to read_topic_assignments_<outname>.tsv
#   assignment_col   - column holding the label ("cluster" or "dominant_topic")
#   group_col        - name of the attached column
#   drop_unassigned  - drop reads absent from the file (TRUE keeps the plotted
#                      read set identical to the topic model's)
#   verbose          - print counts
# Output:
#   add_read_groups() result; stops when the file is missing, lacks the
#   columns or has duplicated RIDs
# ---------------------------------------------------------------------------
add_topic_clusters <- function(res,
                               assignments_file,
                               assignment_col = "cluster",
                               group_col = "cluster",
                               drop_unassigned = TRUE,
                               verbose = TRUE) {

  if (!file.exists(assignments_file))
    stop("assignments file not found: ", assignments_file)

  assign_df <- data.table::fread(assignments_file, header = TRUE, data.table = FALSE)
  if (!all(c("RID", assignment_col) %in% colnames(assign_df)))
    stop("assignments file must have columns RID and ", assignment_col)
  if (any(duplicated(assign_df$RID)))
    stop("duplicated RIDs in ", assignments_file)

  if (verbose)
    cat(sum(!assign_df$RID %in% levels(res$rids_df$RID)),
        "assigned reads are absent from the region result.\n")

  add_read_groups(res,
                  labels = setNames(as.character(assign_df[[assignment_col]]), assign_df$RID),
                  group_col = group_col,
                  drop_unassigned = drop_unassigned,
                  verbose = verbose)
}


# ---------------------------------------------------------------------------
# Split a region result by haplotype, using the HP tag carried by the reads in
# fire.bed. Useful at heterozygous sites, where the two alleles can differ.
#
# Inputs:
#   res          - load_region_results() result
#   group_col    - name of the attached column
#   drop_unknown - drop reads whose haplotype is UNK (unphased)
#   verbose      - print the tag table and counts
# Output:
#   add_read_groups() result; stops when no tag column or no phased read
# ---------------------------------------------------------------------------
add_haplotype_groups <- function(res, group_col = "haplotype",
                                 drop_unknown = TRUE, verbose = TRUE) {

  if (!"tag" %in% colnames(res$rids_df))
    stop("no haplotype tag in rids_df")

  tags <- setNames(as.character(res$rids_df$tag), as.character(res$rids_df$RID))
  if (verbose) {
    cat("haplotype tags:\n"); print(table(tags, useNA = "ifany"))
  }
  if (drop_unknown) tags <- tags[tags %in% c("H1", "H2")]
  if (length(tags) == 0)
    stop("no phased reads in this region")

  add_read_groups(res, labels = tags, group_col = group_col,
                  group_levels = sort(unique(unname(tags))),
                  drop_unassigned = TRUE, verbose = verbose)
}


# ===========================================================================
# cCRE-pair co-accessibility: fibers, FIRE elements and display tracks
#
# Shared by code/co-accessibility/macrophage/coaccess_macrophage.Rmd and
# code/co-accessibility/LCL/LCL_co-access.Rmd. Every path is an argument. The
# accessibility call is always the FIRE elements; m6A, nucleosome and TF
# footprint tracks are display only. Figures: plot_coaccess_pair() in
# plotting_functions.r.
# ===========================================================================

# FiberHMM TF footprint size bins of the co-accessibility figures, half-open
# [lo, hi) as in macrophage FiberHMM/extract/code/split_footprints_by_size.sh
FP_SIZE_BINS <- c("size10-30", "size40-60", "size60-80")

# read configuration at the pair, in the order they are stacked (matches the 2x2)
CONFIG_LEVELS <- c("both accessible", "CRE1 only", "CRE2 only", "neither")

# ENCODE SCREEN v4 cCRE classes (last field of a CRE_ID accession1.accession2.class)
CCRE_CLASSES <- c("PLS", "pELS", "dELS", "CA-H3K4me3", "CA-CTCF", "CA-TF", "CA", "TF")


# ---------------------------------------------------------------------------
# Tabix a BED region into a data.table with the tabix command line tool. The
# query is 1-based inclusive (BED start + 1 to end).
#
# Inputs:
#   path       - bgzipped, tabix-indexed BED
#   chrom      - chromosome
#   start, end - BED (0-based start) coordinates of the region
#   col_names  - names for the leading columns
#   tabix_bin  - tabix executable
# Output:
#   data.table with col_names; an empty table with those columns when nothing
#   overlaps. Stops when the file is missing
# ---------------------------------------------------------------------------
tabix_region <- function(path, chrom, start, end, col_names, tabix_bin) {
  if (!file.exists(path)) stop("file not found: ", path)
  q <- sprintf("%s:%d-%d", chrom, start + 1L, end)     # tabix is 1-based inclusive
  txt <- suppressWarnings(system2(tabix_bin, c(shQuote(path), shQuote(q)),
                                  stdout = TRUE, stderr = FALSE))
  status <- attr(txt, "status")
  if (!is.null(status) && status != 0L)
    stop("tabix failed (status ", status, ") for ", path, " at ", q)
  if (length(txt) == 0)
    return(data.table::data.table(matrix(character(0), ncol = length(col_names),
                                         dimnames = list(NULL, col_names))))
  dt <- data.table::fread(text = paste(txt, collapse = "\n"), header = FALSE,
                          sep = "\t", showProgress = FALSE)
  data.table::setnames(dt, seq_along(col_names), col_names)
  dt
}


# ---------------------------------------------------------------------------
# Load fibers and their FIRE elements for one window, across samples. Fibers
# are keyed "<sample> <read name>", because the same read name can occur in two
# samples.
#
# Inputs:
#   region      - list(chrom, start, end) in BED coordinates
#   samples     - sample names
#   root        - co-accessibility output root holding <s>/<s>.read_spans.bed.gz
#                 (span_source = "read_spans"; unused otherwise)
#   fire_root   - FIRE results root: <s>/additional-outputs-<fire_ver>/fire-peaks/
#                 <s>-<fire_ver>-fire-elements.bed.gz, and for span_source =
#                 "fire_all" <s>/extracted_results/<s>.fire_all.bed.gz
#   fire_ver    - FIRE version in the file names
#   tabix_bin   - tabix executable
#   span_source - "read_spans": one primary alignment per row (both datasets,
#                 their respective *_read_spans.sh); "fire_all": `ft fire --extract --all` rows
#                 collapsed to one span per read (min start, max end of the rows
#                 overlapping the window) with the read's HP tag (LCL)
# Output:
#   list(region, samples, spans, elements); spans and elements carry
#   sample_name and key. Stops when no fiber overlaps the window
# ---------------------------------------------------------------------------
load_region <- function(region, samples, root, fire_root, fire_ver = "v0.1",
                        tabix_bin, span_source = c("read_spans", "fire_all")) {
  span_source <- match.arg(span_source)
  span_cols <- c("chrom", "start", "end", "RID", "mapq", "strand")
  elem_cols <- c("chrom", "start", "end", "RID", "score", "strand",
                 "tstart", "tend", "rgb", "fdr", "HP")

  spans <- data.table::rbindlist(lapply(samples, function(s) {
    if (span_source == "read_spans") {
      d <- tabix_region(file.path(root, s, paste0(s, ".read_spans.bed.gz")),
                        region$chrom, region$start, region$end, span_cols, tabix_bin)
    } else {
      d <- tabix_region(file.path(fire_root, s, "extracted_results", paste0(s, ".fire_all.bed.gz")),
                        region$chrom, region$start, region$end, elem_cols, tabix_bin)
      if (nrow(d)) d <- d[, .(chrom = chrom[1], start = min(start), end = max(end),
                              strand = strand[1], HP = HP[1]), by = RID]
    }
    if (nrow(d)) d[, sample_name := s]
    d
  }), fill = TRUE)

  elements <- data.table::rbindlist(lapply(samples, function(s) {
    d <- tabix_region(file.path(fire_root, s, paste0("additional-outputs-", fire_ver),
                                "fire-peaks",
                                paste0(s, "-", fire_ver, "-fire-elements.bed.gz")),
                      region$chrom, region$start, region$end, elem_cols, tabix_bin)
    if (nrow(d)) d[, sample_name := s]
    d
  }), fill = TRUE)

  if (nrow(spans) == 0) stop("no reads in ", region$chrom, ":", region$start, "-", region$end)

  # a read name is unique within a sample (spans are primary alignments only), but
  # the same name can occur in two samples - key on both
  spans[, key := paste(sample_name, RID)]
  if (nrow(elements)) elements[, key := paste(sample_name, RID)]

  list(region = region, samples = samples, spans = spans, elements = elements)
}


# ---------------------------------------------------------------------------
# Attach the display tracks to a load_region() result: the per-fiber m6A
# calls (`ft extract`), the nucleosome calls and the size-binned TF
# footprints. The accessibility call stays the FIRE elements; a missing track
# file gives a warning and that track is left out.
#
# Inputs:
#   res          - load_region() result
#   keys         - restrict to these fibers (the ones being plotted; the m6A
#                  track is ~240 blocks per fiber)
#   ft_root      - `ft extract` root: <s>/extracted_results/<kind>_by_chr/
#                  <s>.ft_extracted_<kind>.<chr>.bed.gz (BED12, sentinels)
#   hmm_root     - FiberHMM extract root (macrophage): firehmm_footprint/ and
#                  firehmm_tf/ft_by_size/; needed for nuc_source = "fiberhmm"
#                  and tf_source = "by_size"
#   tabix_bin    - tabix executable (size-binned footprints)
#   nuc_source   - "fiberhmm": firehmm_footprint BED13 (no sentinels); "ft":
#                  nuc_by_chr from ft_root
#   tf_source    - "by_size": firehmm_tf/ft_by_size/<bin> BED6 files;
#                  "recalled_tf": FiberHMM v2 <tf_root>/<s>/<s>.recalled_tf.bed.gz
#                  (BED15) binned into FP_SIZE_BINS; "none"
#   tf_root      - FiberHMM v2 results root, for tf_source = "recalled_tf"
#   hmm_label    - FiberHMM file label of a sample (LPS_0 -> LPS0)
#   min_tf_score - minimum per-block score of recalled_tf calls
# Output:
#   res with m6a, nuc and size_fps tables (key, RID, start, end, sample_name;
#   size_fps also class)
# ---------------------------------------------------------------------------
load_ft_tracks <- function(res, keys = NULL, ft_root, hmm_root = NULL, tabix_bin,
                           nuc_source = c("fiberhmm", "ft"),
                           tf_source = c("by_size", "recalled_tf", "none"),
                           tf_root = NULL,
                           hmm_label = function(s) gsub("_", "", s, fixed = TRUE),
                           min_tf_score = 50) {
  nuc_source <- match.arg(nuc_source)
  tf_source <- match.arg(tf_source)
  if ((nuc_source == "fiberhmm" || tf_source == "by_size") && is.null(hmm_root))
    stop("hmm_root is required for nuc_source = 'fiberhmm' or tf_source = 'by_size'")
  if (tf_source == "recalled_tf" && is.null(tf_root))
    stop("tf_root is required for tf_source = 'recalled_tf'")
  region <- res$region
  region_gr <- GenomicRanges::GRanges(region$chrom, IRanges::IRanges(region$start + 1L, region$end))

  # one BED12 track (ft extract, or FiberHMM with fiberhmm = TRUE) for every sample
  grab <- function(kind, fiberhmm = FALSE) {
    data.table::rbindlist(lapply(res$samples, function(s) {
      if (fiberhmm) {
        # Only the file names use the FiberHMM label; read keys retain sample_name.
        s2 <- hmm_label(s)
        p <- file.path(hmm_root, paste0("firehmm_", kind), s2,
                       sprintf("%s_hmm_extracted_%s_%s.bed.gz", s2, kind, region$chrom))
      } else {
        p <- file.path(ft_root, s, "extracted_results", paste0(kind, "_by_chr"),
                       sprintf("%s.ft_extracted_%s.%s.bed.gz", s, kind, region$chrom))
      }
      if (!file.exists(p)) {
        warning(if (fiberhmm) "missing FiberHMM track: " else "missing ft extract track: ",
                p, call. = FALSE)
        return(NULL)
      }
      rows <- read_tabix_region(p, region_gr)
      if (nrow(rows) && !is.null(keys))
        rows <- rows[paste(s, rows[[4]]) %in% keys, , drop = FALSE]
      if (!nrow(rows)) {
        d <- data.table::data.table(RID = character(0), start = numeric(0), end = numeric(0))
        d[, key := character(0)]
        return(d)
      }
      d <- data.table::as.data.table(convert_ft_bed12_to_bed6(rows,
        format = if (fiberhmm) "bed13_fiberhmm" else "bed12_fibertools", source = p))
      # LCL sentinel removal groups by RID; retain this track's original read
      # order, and its rule excluding any remaining zero-width blocks.
      d <- d[order(match(RID, rows[[4]]))]
      d <- d[end > start & end > region$start & start < region$end,
             .(RID, start = as.numeric(start), end = as.numeric(end))]
      d[, key := paste(s, RID)]
      data.table::setcolorder(d, c("key", "RID", "start", "end"))
      if (nrow(d)) d[, sample_name := s]
      d
    }), fill = TRUE)
  }
  # size-binned FiberHMM footprints: plain BED6, one row per footprint
  grab_size <- function(bin) {
    data.table::rbindlist(lapply(res$samples, function(s) {
      s2 <- hmm_label(s)
      p <- file.path(hmm_root, "firehmm_tf", "ft_by_size", bin, s2,
                     sprintf("%s_tf_%s_%s.bed.gz", s2, bin, region$chrom))
      if (!file.exists(p)) {
        warning("missing FiberHMM size track: ", p, call. = FALSE)
        return(NULL)
      }
      d <- tabix_region(p, region$chrom, region$start, region$end,
                        c("chrom", "start", "end", "RID", "size", "strand"), tabix_bin)
      if (nrow(d) == 0) return(NULL)
      d[, key := paste(s, RID)]
      if (!is.null(keys)) d <- d[key %in% keys]
      if (nrow(d) == 0) return(NULL)
      d[, `:=`(sample_name = s, class = bin)]
      d[, .(key, RID, start, end, sample_name, class)]
    }), fill = TRUE)
  }
  # FiberHMM v2 recalled_tf (BED15): every block is a call; bin sizes into FP_SIZE_BINS
  grab_recalled_tf <- function() {
    lo <- as.integer(sub("^size([0-9]+)-([0-9]+)$", "\\1", FP_SIZE_BINS))
    hi <- as.integer(sub("^size([0-9]+)-([0-9]+)$", "\\2", FP_SIZE_BINS))
    data.table::rbindlist(lapply(res$samples, function(s) {
      p <- file.path(tf_root, s, paste0(s, ".recalled_tf.bed.gz"))
      if (!file.exists(p)) {
        warning("missing FiberHMM recalled_tf track: ", p, call. = FALSE)
        return(NULL)
      }
      rows <- read_tabix_region(p, region_gr)
      if (nrow(rows) && !is.null(keys))
        rows <- rows[paste(s, rows[[4]]) %in% keys, , drop = FALSE]
      if (!nrow(rows)) return(NULL)
      b <- data.table::as.data.table(convert_ft_bed12_to_bed6(rows, format = "bed15_fiberhmm_tf",
        keep_block_scores = TRUE, source = p))
      b <- b[block_score >= min_tf_score & end > region$start & start < region$end]
      if (!nrow(b)) return(NULL)
      size <- b$end - b$start
      bin <- vapply(size, function(x) { i <- which(x >= lo & x < hi); if (length(i)) i[1] else NA_integer_ },
                    integer(1))
      b <- b[!is.na(bin)][, class := FP_SIZE_BINS[bin[!is.na(bin)]]]
      if (!nrow(b)) return(NULL)
      b[, .(key = paste(s, RID), RID, start = as.numeric(start), end = as.numeric(end),
            sample_name = s, class)]
    }), fill = TRUE)
  }
  res$m6a <- grab("m6a")
  res$nuc <- if (nuc_source == "fiberhmm") grab("footprint", fiberhmm = TRUE) else grab("nuc")
  res$size_fps <- switch(tf_source,
    by_size     = data.table::rbindlist(lapply(FP_SIZE_BINS, grab_size), fill = TRUE),
    recalled_tf = grab_recalled_tf(),
    none        = data.table::data.table())
  res
}


# ---------------------------------------------------------------------------
# Per-fiber configuration at a cCRE pair. A fiber is accessible at a cCRE
# when ONE of its own FIRE elements overlaps at least fire_overlap_fraction
# of the cCRE length. Separate elements are not summed and the overlap is not
# reciprocal. Both this fraction and read_rule must match the statistics.
#
# Inputs:
#   res        - load_region() result
#   cre1, cre2 - list(start, end) in BED coordinates
#   read_rule  - "any" (a fiber counts if it overlaps the cCRE at all, Kevin's
#                rule) or "contain" (it must span the whole cCRE)
#   fire_overlap_fraction - fraction of cCRE width covered by one FIRE element;
#                0 preserves legacy >= 1 bp behavior; current notebooks use 0.5
# Output:
#   data.table key, sample_name, shared, acc1, acc2, config (factor with
#   CONFIG_LEVELS; NA for fibers that do not cover both cCREs)
# ---------------------------------------------------------------------------
label_reads <- function(res, cre1, cre2, read_rule = c("any", "contain"),
                        fire_overlap_fraction = 0) {
  read_rule <- match.arg(read_rule)
  if (length(fire_overlap_fraction) != 1L || !is.finite(fire_overlap_fraction) ||
      fire_overlap_fraction < 0 || fire_overlap_fraction > 1)
    stop("fire_overlap_fraction must be one number in [0, 1]")
  if (cre1$end <= cre1$start || cre2$end <= cre2$start)
    stop("cCREs must have positive width")
  sp <- res$spans
  el <- res$elements

  # which fibers cover the interval s-e under read_rule
  covers <- function(s, e) {
    if (read_rule == "contain") sp$start <= s & sp$end >= e
    else sp$start < e & sp$end > s
  }
  cov1 <- covers(cre1$start, cre1$end)
  cov2 <- covers(cre2$start, cre2$end)

  # Apply the threshold to each single element, then take distinct fiber keys.
  hit <- function(s, e) {
    if (nrow(el) == 0) return(character(0))
    overlap <- pmin(el$end, e) - pmax(el$start, s)
    unique(el$key[overlap > 0 & overlap >= fire_overlap_fraction * (e - s)])
  }
  a1 <- sp$key %in% hit(cre1$start, cre1$end)
  a2 <- sp$key %in% hit(cre2$start, cre2$end)

  shared <- cov1 & cov2
  cfg <- ifelse(!shared, NA_character_,
         ifelse( a1 &  a2, "both accessible",
         ifelse( a1 & !a2, "CRE1 only",
         ifelse(!a1 &  a2, "CRE2 only", "neither"))))

  # `key` is a reserved argument of data.table(), so assign that column afterwards
  out <- data.table::data.table(sample_name = sp$sample_name,
                                shared = shared, acc1 = a1, acc2 = a2,
                                config = factor(cfg, levels = CONFIG_LEVELS))
  out[, key := sp$key]
  data.table::setcolorder(out, "key")[]
}


# ---------------------------------------------------------------------------
# Order fibers for the raster: by sample, then configuration, then start -
# the visual counterpart of the 2x2 table (Kevin's cluster = "fire_configs").
#
# Inputs:
#   res     - load_region() result
#   labels  - label_reads() result
#   samples - sample (facet) order
# Output:
#   character vector of fiber keys
# ---------------------------------------------------------------------------
order_reads <- function(res, labels, samples) {
  d <- merge(res$spans[, .(key, sample_name, start)], labels[, .(key, config)],
             by = "key", all.x = TRUE)
  d[, samp_rank := match(sample_name, samples)]
  d[, cfg_rank := ifelse(is.na(config), length(CONFIG_LEVELS) + 1L, as.integer(config))]
  data.table::setorder(d, samp_rank, cfg_rank, start)
  d$key
}


# ---------------------------------------------------------------------------
# Fraction of covering fibers methylated at each observed A position, per
# sample, unbinned: distinct methylated fibers / fibers whose span covers the
# position.
#
# Inputs:
#   m       - m6A calls (key, start, sample_name), already restricted to the drawn fibers
#   sp      - fiber spans (key, start, end, sample_name)
#   samples - sample (facet) order
# Output:
#   data.table pos, sample_name (factor), frac; empty when there are no calls
# ---------------------------------------------------------------------------
coaccess_m6a_fraction <- function(m, sp, samples) {
  if (is.null(m) || nrow(m) == 0) return(data.table::data.table())
  prop <- data.table::rbindlist(lapply(samples, function(s) {
    ms <- m[sample_name == s]
    ss <- sp[sample_name == s]
    if (nrow(ms) == 0 || nrow(ss) == 0) return(NULL)
    d <- ms[, .(met_n = data.table::uniqueN(key)), by = .(pos = start)]
    d[, cov_n := vapply(pos, function(x) sum(ss$start <= x & ss$end > x),
                        numeric(1))]
    d[, .(pos, sample_name = s, frac = met_n / cov_n)]
  }))
  if (nrow(prop) > 0) prop[, sample_name := factor(sample_name, levels = samples)]
  prop
}


# ---------------------------------------------------------------------------
# The pair's 2x2 read off the fibers that cover both cCREs, with the Fisher
# test of the co-accessibility tables (fisher.test(table + 1)).
#
# Inputs:
#   labels - label_reads() result
# Output:
#   one-row data.frame both, cre1_only, cre2_only, neither, n_shared,
#   fisher_or (conditional MLE), fisher_p
# ---------------------------------------------------------------------------
coaccess_2x2 <- function(labels) {
  n <- table(factor(labels$config[!is.na(labels$config)], levels = CONFIG_LEVELS))
  # rows CRE1 FALSE/TRUE, columns CRE2 FALSE/TRUE, as table(fire_region1, fire_region2)
  tab <- matrix(c(n[["neither"]], n[["CRE2 only"]], n[["CRE1 only"]], n[["both accessible"]]),
                2, byrow = TRUE)
  ft <- stats::fisher.test(tab + 1)
  data.frame(both = n[["both accessible"]], cre1_only = n[["CRE1 only"]],
             cre2_only = n[["CRE2 only"]], neither = n[["neither"]], n_shared = sum(n),
             fisher_or = unname(ft$estimate), fisher_p = ft$p.value)
}


# ---------------------------------------------------------------------------
# cCRE class of a CRE_ID (accession1.accession2.class), e.g. "CA-CTCF".
#
# Inputs:
#   ids - character vector of CRE_IDs
# Output:
#   character vector of classes
# ---------------------------------------------------------------------------
cre_class_from_id <- function(ids) {
  sub("^[^.]*\\.[^.]*\\.", "", ids)
}


# ---------------------------------------------------------------------------
# Keep the cCRE pairs whose member classes match cre_types. Classes are taken
# from the CRE_IDs (exact match, any of CCRE_CLASSES) and the pair is treated
# as unordered, since which member is CRE1 is arbitrary between tables.
#
# Inputs:
#   pairs      - pair table with CRE_ID columns cre1_col and cre2_col
#   cre_types  - classes to accept; NULL keeps every pair (any cCRE type)
#   cre_match  - "either": at least one member in cre_types; "both": both
#                members in cre_types; "pair": one member in cre_types and the
#                other in cre_types2
#   cre_types2 - second class set, for cre_match = "pair"
#   cre1_col, cre2_col - columns holding the two CRE_IDs
# Output:
#   the matching rows of pairs (same class as the input)
# ---------------------------------------------------------------------------
select_cre_pairs <- function(pairs, cre_types = NULL,
                             cre_match = c("either", "both", "pair"),
                             cre_types2 = NULL, cre1_col = "CRE1", cre2_col = "CRE2") {
  cre_match <- match.arg(cre_match)
  if (is.null(cre_types)) return(pairs)
  unknown <- setdiff(c(cre_types, cre_types2), CCRE_CLASSES)
  if (length(unknown))
    stop("unknown cCRE class: ", paste(unknown, collapse = ", "),
         "; expected one of ", paste(CCRE_CLASSES, collapse = ", "))
  if (cre_match == "pair" && is.null(cre_types2)) stop("cre_match = 'pair' needs cre_types2")
  c1 <- cre_class_from_id(pairs[[cre1_col]])
  c2 <- cre_class_from_id(pairs[[cre2_col]])
  keep <- switch(cre_match,
    either = c1 %in% cre_types | c2 %in% cre_types,
    both   = c1 %in% cre_types & c2 %in% cre_types,
    pair   = (c1 %in% cre_types & c2 %in% cre_types2) | (c1 %in% cre_types2 & c2 %in% cre_types))
  if (data.table::is.data.table(pairs)) pairs[keep] else pairs[keep, , drop = FALSE]
}
