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
