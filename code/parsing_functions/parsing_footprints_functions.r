# parsing_footprints_functions.r
#
# Shared Fiber-seq m6A/CpG and footprint parsing, read metadata, and methylation
# matrix functions. Adapted from Kevin Luo's `process_fiberseq_data.R` and
# `topic_model_utils.R`; topic_modelling_functions.r sources this file.
# These functions rely on dplyr, GenomicRanges, IRanges, Rsamtools, data.table,
# and Matrix being available as they are now. Callers load the packages.

read_tabix_region <- function(ft_extracted_file, region_gr) {
  tabix_index <- Rsamtools::TabixFile(ft_extracted_file)
  compressed_records <- Rsamtools::scanTabix(tabix_index, param = region_gr)
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
    df <- as.data.frame(df)
    df <- df[order(df$RID, -(df$end - df$start)), , drop = FALSE]
    df <- df[!duplicated(df$RID), , drop = FALSE]
  }
  return(df)
}


################################################
# this function retunrs one row per bed12 block

# output
# chr start   end                                          RID score strand
# 1 chr5 10368 10369 m84241_260613_040640_s4/177673938/ccs      29      +
#   2 chr5 10674 10675 m84241_260613_040640_s4/177673938/ccs      29      +
#   3 chr5 11296 11297 m84241_260613_040640_s4/177673938/ccs      29      +

################################################

convert_ft_bed12_to_bed6 <- function(bed12_df, include_read_start_end = FALSE,
                                   drop_sentinels = TRUE, keep_block_scores = FALSE) {
  if (nrow(bed12_df) == 0) {
    return(data.frame())
  }
  # Existing calls retain the BED12-only contract. FiberHMM callers opt in by
  # keeping all blocks or requesting its per-block scores.
  legacy <- isTRUE(drop_sentinels) && !isTRUE(keep_block_scores)
  if (legacy && ncol(bed12_df) != 12L) stop("Expected 12 BED columns")
  if (!ncol(bed12_df) %in% c(12L, 13L)) stop("Expected 12 or 13 BED columns")
  has_block_scores <- ncol(bed12_df) == 13L
  colnames(bed12_df) <- c('chr','start','end','RID','score','strand','read_start','read_end','rgb','blockCount','blockSizes','blockStarts',
                         if (has_block_scores) 'blockScores')

  # expand BED12 into one row per block
  block_sizes_list  <- strsplit(sub(",$", "", bed12_df$blockSizes),  ",", fixed = TRUE)
  block_starts_list <- strsplit(sub(",$", "", bed12_df$blockStarts), ",", fixed = TRUE)
  n_blocks <- lengths(block_sizes_list)
  offsets <- suppressWarnings(as.integer(unlist(block_starts_list)))
  sizes <- suppressWarnings(as.integer(unlist(block_sizes_list)))
  if (anyNA(bed12_df$blockCount) || any(n_blocks != bed12_df$blockCount) ||
      any(lengths(block_starts_list) != n_blocks) || anyNA(offsets) || anyNA(sizes))
    stop("Invalid BED12 blocks")
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


################################################
# takes ft extracted bed12 file and find a single genomic region 
################################################
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
    extracted_data_bed12 <- extracted_data_bed12 %>%
      dplyr::group_by(RID) %>%
      dplyr::slice_max(end - start, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
  }

  reads <- convert_ft_bed12_to_bed6(extracted_data_bed12, include_read_start_end = TRUE)
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

###################################
# creates one summary row per read
# 1. From an existing reads table produced by extract_ft_region_reads().
# 2. Directly from the original fibertools BED12 file for a requested genomic region
# output = read level dataframe: RID    chr    start    end    strand
###################################
extract_ft_read_info <- function(reads, ft_extracted_file, region, keep_columns = NULL) {
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

    rids_df <- rids_df %>% dplyr::mutate(start = start + 1)
    rids_df <- rids_df %>% dplyr::arrange(RID)
  }

  return(as.data.frame(rids_df))
}

################################################
# convert modification calls into a read by position matrix

# each row is one read
# each column is one genomic A or CpG position

# 1 = modified call
# 0 = no modified call
# NA = the read does not cover that position
################################################
get_sparse_met_mat <- function(reads, rids_df, window_start = NULL, window_end = NULL, all_met_pos = NULL,
                               base = c("A", "CG"), require_full_span = FALSE) {
  base <- match.arg(base)

  if ("base" %in% colnames(reads))
    reads <- reads %>% dplyr::filter(.data$base == !!base)

  if ("chrom" %in% colnames(reads))
    reads <- reads %>% dplyr::rename(chr = chrom)
  if ("ref_position" %in% colnames(reads))
    reads <- reads %>% dplyr::rename(pos = ref_position)

  if (missing(rids_df)) {
    rids_df <- extract_ft_read_info(reads)
  }

  if (is.null(all_met_pos)) {
    all_met_pos <- sort(unique(reads$pos))
  }

  if (!is.null(window_start) & !is.null(window_end)) {
    all_met_pos <- all_met_pos[all_met_pos >= window_start & all_met_pos <= window_end]
  }

  n_rids <- nrow(rids_df)
  n_pos <- length(all_met_pos)
  if (require_full_span) {
    stopifnot(length(window_start) == 1L, length(window_end) == 1L,
      !anyDuplicated(rids_df$RID), !anyNA(rids_df$RID),
      all(rids_df$start <= window_start), all(rids_df$end >= window_end))
    pairs <- unique(data.frame(
      read_index = match(as.character(reads$RID), as.character(rids_df$RID)),
      site_index = match(reads$pos, all_met_pos)))
    pairs <- pairs[!is.na(pairs$read_index) & !is.na(pairs$site_index), , drop = FALSE]
    return(Matrix::sparseMatrix(i = pairs$read_index, j = pairs$site_index,
      x = rep(1, nrow(pairs)), dims = c(n_rids, n_pos),
      dimnames = list(as.character(rids_df$RID), as.character(all_met_pos))))
  }
  met_mat <- matrix(NA, nrow = n_rids, ncol = n_pos)
  rownames(met_mat) <- as.character(rids_df$RID)
  colnames(met_mat) <- all_met_pos

  if (n_rids == 0 | n_pos == 0) {
    return(met_mat)
  }

  # mark covered positions (0 = covered but not methylated, NA = uncovered)
  in_cov <- outer(rids_df$start, all_met_pos, "<=") & outer(rids_df$end, all_met_pos, ">=")
  met_mat[in_cov] <- 0

  # set methylated positions to 1 for matched RID/position pairs
  rid_idx <- match(as.character(reads$RID), as.character(rids_df$RID))
  pos_idx <- match(reads$pos, all_met_pos)
  keep <- !is.na(rid_idx) & !is.na(pos_idx)
  if (any(keep)) {
    met_mat[cbind(rid_idx[keep], pos_idx[keep])] <- 1
  }

  met_mat <- as(met_mat, "dgCMatrix")
  return(met_mat)
}

# Read one modification modality for each sample, preserving sample order.
# Empty samples return NULL; callers drop them before combining the tables.
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
