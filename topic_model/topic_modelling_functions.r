# topic_modelling_functions.r
#
# Functions for topic modelling on Fiber-seq promoter data, sourced by
# cluster_footprints_notebooks/m6a_promoter_topic_modelling.Rmd.
#
# The m6A / methylation functions are adapted from Kevin Luo's
# `process_fiberseq_data.R` and `topic_model_utils.R`. The footprint functions
# at the bottom build the same multi-class footprint occupancy matrix used in
# `cluster_by_footprint_class.ipynb`, for clustering reads via topic modelling.


# =============================================================================
# m6A / methylation functions
# =============================================================================

read_tabix_region <- function(ft_extracted_file, region_gr) {
  tabix_index <- TabixFile(ft_extracted_file)
  compressed_records <- scanTabix(tabix_index, param = region_gr)
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

read_ft_bed12 <- function(bed_file) {
  df <- data.table::fread(bed_file)
  names(df) <- c('chr','start','end','RID','score','strand','x1','x2','rgb','blockCount','blockSizes','blockStarts')
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

convert_ft_bed12_to_bed6 <- function(bed12_df, include_read_start_end = FALSE) {
  if (nrow(bed12_df) == 0) {
    return(data.frame())
  }
  colnames(bed12_df) <- c('chr','start','end','RID','score','strand','read_start','read_end','rgb','blockCount','blockSizes','blockStarts')

  # expand BED12 into one row per block
  block_sizes_list  <- strsplit(sub(",$", "", bed12_df$blockSizes),  ",", fixed = TRUE)
  block_starts_list <- strsplit(sub(",$", "", bed12_df$blockStarts), ",", fixed = TRUE)
  n_blocks <- lengths(block_sizes_list)
  starts <- rep(bed12_df$start, n_blocks) + as.integer(unlist(block_starts_list))
  sizes  <- as.integer(unlist(block_sizes_list))
  bed6_df <- data.frame(
    chr    = rep(bed12_df$chr,    n_blocks),
    start  = starts,
    end    = starts + sizes,
    RID    = rep(bed12_df$RID,    n_blocks),
    score  = rep(bed12_df$score,  n_blocks),
    strand = rep(bed12_df$strand, n_blocks),
    stringsAsFactors = FALSE
  )

  # remove the first and last rows for each RID (they are added by ft to make the bed12 format work)
  bed6_df <- bed6_df %>% dplyr::group_by(RID) %>% dplyr::slice(-c(1, n())) %>% dplyr::ungroup()

  if (include_read_start_end) {
    read_start_end_df <- bed12_df %>% dplyr::select(RID, read_start, read_end)
    bed6_df <- bed6_df %>% dplyr::left_join(read_start_end_df, by = "RID")
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
get_sparse_met_mat <- function(reads, rids_df, window_start = NULL, window_end = NULL, all_met_pos = NULL, base = c("A", "CG")) {
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

filter_met_mat <- function(met_mat, verbose = FALSE) {
  if (nrow(met_mat) == 0 | ncol(met_mat) == 0) {
    if (verbose) cat("Met matrix is empty. \n")
    return(met_mat)
  }
  # remove rows with NAs
  row_idx_na <- which(rowSums(is.na(met_mat)) > 0)
  if (length(row_idx_na) > 0) {
    if (verbose) cat("Removing", length(row_idx_na), "rows with NAs in the met matrix.\n")
    met_mat_cleaned <- met_mat[-row_idx_na, ]
  } else {
    met_mat_cleaned <- met_mat
  }
  # remove columns with all zeros
  col_idx_zero <- which(colSums(met_mat_cleaned) == 0)
  if (length(col_idx_zero) > 0) {
    if (verbose) cat("Removing", length(col_idx_zero), "columns with all zeros in the met matrix.\n")
    met_mat_cleaned <- met_mat_cleaned[, -col_idx_zero]
  }
  return(met_mat_cleaned)
}

################################################
# builds a combined read by position methylation for one genomic region across multiple samples
################################################
assemble_region_met_data <- function(sample_names, region_chr, region_start, region_end,
                                     ft_result_dir = "/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results",
                                     filter_NAs = FALSE,
                                     verbose = TRUE) {

  region_gr <- GRanges(seqnames = region_chr, ranges = IRanges(start = region_start, end = region_end))
  if (length(region_gr) != 1)
    stop("region_gr should be a single region!")

  if (verbose)
    cat("Assembling region data for", length(sample_names), "samples in",
        paste0(region_chr, ":", region_start, "-", region_end), "...\n")

  # assemble data for metA
  if (verbose) cat("Assembling data for m6A...\n")
  reads_list <- lapply(sample_names, function(sample_name) {
    extracted_dir <- file.path(ft_result_dir, sample_name, "extracted_results", "m6a_by_chr")
    extracted_file <- file.path(extracted_dir, paste0(sample_name, ".ft_extracted_m6a.", region_chr, ".bed.gz"))
    sample_reads <- extract_ft_region_reads(extracted_file, region_gr, keep_pos_in_region_only = FALSE, verbose = verbose)
    sample_reads <- sample_reads %>% dplyr::mutate(sample_name = sample_name, .before = 1)
    return(sample_reads)
  })
  combined_reads_metA <- do.call(rbind, reads_list)
  combined_reads_metA <- combined_reads_metA %>% dplyr::mutate(base = "A")

  combined_rids_metA_df <- extract_ft_read_info(combined_reads_metA, keep_columns = c("sample_name", "score"))
  cat("Number of RIDs in metA data:", nrow(combined_rids_metA_df), "\n")

  # assemble data for metCG
  if (verbose) cat("Assembling data for metCG...\n")
  reads_list <- lapply(sample_names, function(sample_name) {
    extracted_dir <- file.path(ft_result_dir, sample_name, "extracted_results", "cpg_by_chr")
    extracted_file <- file.path(extracted_dir, paste0(sample_name, ".ft_extracted_cpg.", region_chr, ".bed.gz"))
    sample_reads <- extract_ft_region_reads(extracted_file, region_gr, keep_pos_in_region_only = FALSE, verbose = verbose)
    sample_reads <- sample_reads %>% dplyr::mutate(sample_name = sample_name, .before = 1)
    return(sample_reads)
  })
  combined_reads_metCG <- do.call(rbind, reads_list)
  combined_reads_metCG <- combined_reads_metCG %>% dplyr::mutate(base = "CG")

  combined_rids_metCG_df <- extract_ft_read_info(combined_reads_metCG, keep_columns = c("sample_name", "score"))
  cat("Number of RIDs in metCG data:", nrow(combined_rids_metCG_df), "\n")

  combined_rids_df <- unique(rbind(combined_rids_metA_df, combined_rids_metCG_df))
  rownames(combined_rids_df) <- NULL
  cat("Number of RIDs in combined data:", nrow(combined_rids_df), "\n")

  combined_reads <- rbind(dplyr::mutate(combined_reads_metA, base = "A"),
                          dplyr::mutate(combined_reads_metCG, base = "CG"))

  metA_mat <- get_sparse_met_mat(combined_reads, combined_rids_df, window_start = region_start, window_end = region_end, base = "A")
  if (nrow(metA_mat) == 0 | ncol(metA_mat) == 0) {
    if (verbose) cat("metA matrix is empty. \n")
  } else {
    colnames(metA_mat) <- paste0("A_", colnames(metA_mat))
    if (verbose) cat("Dimensions of the metA matrix:", nrow(metA_mat), "x", ncol(metA_mat), "\n")
  }

  metCG_mat <- get_sparse_met_mat(combined_reads, combined_rids_df, window_start = region_start, window_end = region_end, base = "CG")
  if (nrow(metCG_mat) == 0 | ncol(metCG_mat) == 0) {
    if (verbose) cat("metCG matrix is empty. \n")
  } else {
    colnames(metCG_mat) <- paste0("CG_", colnames(metCG_mat))
    if (verbose) cat("Dimensions of the metCG matrix:", nrow(metCG_mat), "x", ncol(metCG_mat), "\n")
  }

  if (!all(row.names(metA_mat) == row.names(metCG_mat))) {
    stop("RIDs in metA and metCG matrices do not match!")
  }
  combined_met_mat <- cbind(metA_mat, metCG_mat)
  if (nrow(combined_met_mat) == 0 | ncol(combined_met_mat) == 0) {
    if (verbose) cat("Combined met matrix is empty. \n")
  } else {
    if (verbose) cat("Dimensions of the combined met matrix:", nrow(combined_met_mat), "x", ncol(combined_met_mat), "\n")
  }

  if (filter_NAs) {
    combined_met_mat <- filter_met_mat(combined_met_mat, verbose = verbose)
    metA_mat <- combined_met_mat[, grep("^A_", colnames(combined_met_mat))]
    metCG_mat <- combined_met_mat[, grep("^CG_", colnames(combined_met_mat))]
    if (verbose) {
      cat("Dimensions of the metA matrix after filtering:", nrow(metA_mat), "x", ncol(metA_mat), "\n")
      cat("Dimensions of the metCG matrix after filtering:", nrow(metCG_mat), "x", ncol(metCG_mat), "\n")
      cat("Dimensions of the combined met matrix after filtering:", nrow(combined_met_mat), "x", ncol(combined_met_mat), "\n")
    }
  }

  return(list(chr = region_chr,
              start = region_start,
              end = region_end,
              sample_names = sample_names,
              rids_df = combined_rids_df,
              reads = combined_reads,
              metA_mat = metA_mat,
              metCG_mat = metCG_mat,
              combined_met_mat = combined_met_mat))
}


################################################
# fits a binomial topic model to a binary matrix
################################################

fit_binomial_topic_model <- function(X, k, numem = 100) {
  fit_pois <- init_poisson_nmf(X, k = k)
  fit_pois <- fit_poisson_nmf(X, fit0 = fit_pois,
                              control = list(extrapolate = TRUE),
                              verbose = "none")
  fit_binom <- poisson2binom(X, fit_pois, numem = numem)
  return(fit_binom)
}

################################################
# visualize the fitted topic model

# per read topic memberships from fit$L: structure stacked bar plot
## each bar represents one read
## each color represents one topic

# topic-specific genomic profiles from fit$F: one vertical bar per topic
## x = genomic position
## y = topic specific probability of observing a modification at that position
################################################

plot_topic_model_result <- function(fit,
                                    grouping,
                                    region = NULL,
                                    topic_colors,
                                    color_by = c("topic", "base"),
                                    plot = c("all", "structure", "topics"),
                                    title = NULL,
                                    rel_heights = c(1, 2)) {
  color_by <- match.arg(color_by)
  plot <- match.arg(plot)

  n_topics <- ncol(fit$L)
  if (length(topic_colors) < n_topics) {
    stop("Not enough colors provided for the number of topics.")
  }

  if (plot == "structure" || plot == "all") {
    gg_structure <- structure_plot(fit$L, colors = topic_colors,
                                   topics = 1:n_topics,
                                   grouping = grouping) +
                    theme_cowplot(font_size = 10) +
                    ggtitle(title) +
                    theme(plot.title = element_text(hjust = 0.5),
                          axis.text.x  = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x = element_blank())
  } else {
    gg_structure <- NULL
  }

  if (plot == "topics" || plot == "all") {
    df <- as.data.frame(fit$F)
    if (grepl("^A|^CG", rownames(df)[1])) {
      df <- df %>% dplyr::mutate(base_pos = rownames(df)) %>%
        tidyr::separate(base_pos, into = c("base", "pos"), sep = "_") %>%
        dplyr::mutate(pos = as.numeric(pos))
    } else {
      df <- df %>% dplyr::mutate(pos = as.numeric(rownames(df)))
    }
    df2 <- df %>% pivot_longer(cols = starts_with("k"), names_to = "topic")
    p_list <- vector("list", length = n_topics)
    topics <- unique(df2$topic)
    p_list <- lapply(1:length(p_list), function(i) {
      topic <- topics[i]
      gg_tmp <- ggplot(df2[df2$topic == topic, ], aes(x = pos, y = value, fill = base)) +
        geom_col() +
        xlab(NULL) + ylab(paste0("k", i)) +
        theme_cowplot(font_size = 10) +
        theme(legend.position = "none")

      if (color_by == "topic") {
        gg_tmp <- gg_tmp + scale_fill_manual(values = topic_colors[i])
      } else if (color_by == "base") {
        gg_tmp <- gg_tmp + scale_fill_manual(values = c("A" = "blue", "CG" = "red"))
      }

      if (!is.null(region)) {
        gg_tmp <- gg_tmp + scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0)) +
          coord_cartesian(xlim = c(region$start, region$end), expand = FALSE)
      }
      if (i < n_topics) {
        gg_tmp <- gg_tmp + theme(
          axis.text.x  = element_text(color = NA),
          axis.ticks.x = element_line(color = NA),
          axis.title.x = element_text(color = NA),
          axis.line.x  = element_line(color = NA),
          panel.grid.minor = element_blank()
        )
      }
      return(gg_tmp)
    })
    names(p_list) <- paste0("k", 1:n_topics)
    gg_topics <- cowplot::plot_grid(plotlist = p_list, ncol = 1, align = "v", axis = "tblr",
                                    rel_heights = rep(1, length(p_list)))
  } else {
    gg_topics <- NULL
  }

  if (plot == "all") {
    gg_plots <- cowplot::plot_grid(gg_structure, gg_topics, ncol = 1, rel_heights = rel_heights)
  } else if (plot == "structure") {
    gg_plots <- gg_structure
  } else if (plot == "topics") {
    gg_plots <- gg_topics
  } else {
    gg_plots <- NULL
  }

  return(gg_plots)
}

# =============================================================================
# Footprint functions
# =============================================================================
#
# These build the same multi-class footprint occupancy matrix as
# `cluster_by_footprint_class.ipynb`, reading the tabix-indexed BED files at
# FP_BED_DIR/<class_subdir>/<chrom>_footprints.bed.gz
# (columns: chrom, fp_start, fp_end, read_name, fp_class, fp_mid).

FP_BED_DIR <- "/project/spott/cshan/fiber-seq/results/PolII/footprint_summary_beds"

# Merged 3-class definition (matches the single-class / TSS-proximal sections of
# cluster_by_footprint_class.ipynb): each class maps to one or more BED subdirs.
FP_CLASS_SOURCES <- list(
  PPP_unknown     = c("PPP", "unknown"),
  PIC             = c("PIC"),
  FIRE_nucleosome = c("FIRE_nucleosome")
)

FP_CLASS_COLORS <- c(
  PPP_unknown     = "#FF69B4",
  PIC             = "#4169E1",
  FIRE_nucleosome = "#FFA500"
)

# Read one class subdir's tabix-indexed BED for a region; returns read_name,
# fp_start, fp_end (genomic, BED half-open coordinates).
read_footprint_class_region <- function(class_subdir, region_gr, fp_bed_dir = FP_BED_DIR) {
  region_chr <- as.character(seqnames(region_gr))
  bed_file <- file.path(fp_bed_dir, class_subdir, paste0(region_chr, "_footprints.bed.gz"))
  if (!file.exists(bed_file)) {
    warning("Missing footprint BED: ", bed_file)
    return(data.frame(read_name = character(), fp_start = integer(), fp_end = integer()))
  }
  df <- read_tabix_region(bed_file, region_gr)
  if (nrow(df) == 0) {
    return(data.frame(read_name = character(), fp_start = integer(), fp_end = integer()))
  }
  # BED cols: chrom, fp_start, fp_end, read_name, fp_class, fp_mid
  data.frame(read_name = as.character(df[[4]]),
             fp_start  = as.integer(df[[2]]),
             fp_end    = as.integer(df[[3]]),
             stringsAsFactors = FALSE)
}

# Build a per-read, full-resolution binary footprint occupancy matrix over a
# strand-aware TSS-relative window (-hw .. +hw). For each merged class the
# per-class binary vector is concatenated into a single combined feature vector,
# giving a (reads x [n_class * (2*hw+1)]) sparse matrix. Columns are named
# "<class>@<relpos>". Reads are the union of all reads with >=1 footprint (of any
# merged class) in the window. All-zero rows and columns are dropped, as
# required by fastTopics.
assemble_footprint_mat <- function(gene, region_chr, tss, strand, hw = 300,
                                   class_sources = FP_CLASS_SOURCES,
                                   fp_bed_dir = FP_BED_DIR,
                                   verbose = TRUE) {
  view_start <- max(1, tss - hw)
  view_end   <- tss + hw
  region_gr  <- GRanges(seqnames = region_chr,
                        ranges = IRanges(start = view_start, end = view_end))
  n_pos   <- 2 * hw + 1
  rel_pos <- -hw:hw
  classes <- names(class_sources)

  # per-class footprints, merging source subdirs, keyed by read_name
  per_class_reads <- lapply(classes, function(cls) {
    parts <- lapply(class_sources[[cls]],
                    function(src) read_footprint_class_region(src, region_gr, fp_bed_dir))
    parts <- parts[vapply(parts, nrow, integer(1)) > 0]
    if (length(parts) == 0) {
      return(data.frame(read_name = character(), fp_start = integer(), fp_end = integer()))
    }
    do.call(rbind, parts)
  })
  names(per_class_reads) <- classes

  # union of reads with any footprint in the window
  read_ids <- sort(unique(unlist(lapply(per_class_reads, function(d) d$read_name),
                                 use.names = FALSE)))
  n_reads <- length(read_ids)
  if (verbose) cat(gene, ": n_reads =", n_reads, "\n")
  if (n_reads == 0) return(NULL)
  rid_index <- setNames(seq_along(read_ids), read_ids)

  # collect sparse-matrix triplets (row, col) of covered positions per class
  i_list <- vector("list", length(classes))
  j_list <- vector("list", length(classes))
  for (ci in seq_along(classes)) {
    d <- per_class_reads[[ci]]
    col_off <- (ci - 1) * n_pos
    if (nrow(d) == 0) next

    # strand-aware TSS-relative bounds for each footprint (BED half-open -> inclusive end)
    g_end_incl <- d$fp_end - 1
    if (strand == "+") {
      lo <- d$fp_start - tss
      hi <- g_end_incl - tss
    } else {
      lo <- tss - g_end_incl
      hi <- tss - d$fp_start
    }
    lo <- pmax(lo, -hw)
    hi <- pmin(hi, hw)
    keep <- lo <= hi
    if (!any(keep)) next

    rows   <- rid_index[d$read_name]
    col_lo <- lo + hw + 1 + col_off
    col_hi <- hi + hw + 1 + col_off
    lens   <- hi - lo + 1

    j_list[[ci]] <- unlist(Map(seq.int, col_lo[keep], col_hi[keep]), use.names = FALSE)
    i_list[[ci]] <- rep(rows[keep], lens[keep])
  }

  occ_mat <- Matrix::sparseMatrix(
    i = unlist(i_list), j = unlist(j_list), x = 1,
    dims = c(n_reads, n_pos * length(classes))
  )
  occ_mat@x <- rep(1, length(occ_mat@x))   # binary: collapse any overlap counts to 1
  colnames(occ_mat) <- unlist(lapply(classes, function(cls) paste0(cls, "@", rel_pos)))
  rownames(occ_mat) <- read_ids

  # drop all-zero rows and columns (fastTopics requires nonzero rows/cols)
  row_keep <- Matrix::rowSums(occ_mat) > 0
  occ_mat  <- occ_mat[row_keep, , drop = FALSE]
  col_keep <- Matrix::colSums(occ_mat) > 0
  occ_mat  <- occ_mat[, col_keep, drop = FALSE]
  if (verbose)
    cat("  footprint matrix after filtering:", nrow(occ_mat), "reads x",
        ncol(occ_mat), "features\n")

  list(gene = gene, chr = region_chr, tss = tss, strand = strand, hw = hw,
       classes = classes, rel_pos = rel_pos,
       read_ids = read_ids[row_keep], occ_mat = occ_mat)
}

# Plot the per-topic footprint factor profiles. Feature names "<class>@<pos>" are
# split into class + TSS-relative position; one panel per topic, bars colored by
# footprint class.
plot_footprint_topic_factors <- function(fit, class_colors = FP_CLASS_COLORS,
                                         region = NULL) {
  df <- as.data.frame(fit$F)
  feat <- rownames(df)
  df$class <- factor(sub("@.*$", "", feat), levels = names(class_colors))
  df$pos   <- as.numeric(sub("^.*@", "", feat))

  n_topics <- sum(grepl("^k", colnames(df)))
  df2 <- df %>% pivot_longer(cols = starts_with("k"), names_to = "topic")
  topics <- unique(df2$topic)

  p_list <- lapply(seq_along(topics), function(i) {
    gg_tmp <- ggplot(df2[df2$topic == topics[i], ],
                     aes(x = pos, y = value, fill = class)) +
      geom_col(width = 1) +
      scale_fill_manual(values = class_colors, drop = FALSE) +
      xlab(NULL) + ylab(paste0("k", i)) +
      theme_cowplot(font_size = 10) +
      theme(legend.position = if (i == 1) "top" else "none")
    if (!is.null(region)) {
      gg_tmp <- gg_tmp +
        coord_cartesian(xlim = c(region$start, region$end), expand = FALSE)
    }
    if (i < n_topics) {
      gg_tmp <- gg_tmp + theme(axis.text.x = element_blank(),
                               axis.ticks.x = element_blank())
    }
    gg_tmp
  })
  cowplot::plot_grid(plotlist = p_list, ncol = 1, align = "v", axis = "tblr")
}

# Heatmap of the footprint occupancy matrix with reads ordered by cluster. Each
# cell shows the covering footprint class (later classes in `classes` overwrite
# earlier ones); uncovered positions are white. Rows are annotated with their
# cluster color. Drawn to the active graphics device (call inside pdf() to save).
plot_footprint_cluster_heatmap <- function(occ_mat, clusters, classes, class_colors,
                                           cluster_colors, main = NULL) {
  feat   <- colnames(occ_mat)
  fclass <- sub("@.*$", "", feat)
  fpos   <- as.numeric(sub("^.*@", "", feat))
  positions <- sort(unique(fpos))
  M <- as.matrix(occ_mat)

  # integer class code per (read, position): 0 = uncovered, ci = class i
  code <- matrix(0L, nrow(M), length(positions))
  for (ci in seq_along(classes)) {
    cols <- which(fclass == classes[ci])
    if (!length(cols)) next
    full <- matrix(0L, nrow(M), length(positions))
    full[, match(fpos[cols], positions)] <- M[, cols, drop = FALSE]
    code[full > 0] <- ci
  }

  ord      <- order(clusters)
  code_ord <- code[ord, , drop = FALSE]
  row_side <- cluster_colors[as.integer(clusters)][ord]

  K <- length(classes)
  gplots::heatmap.2(
    x = code_ord, Rowv = FALSE, Colv = FALSE, dendrogram = "none", scale = "none",
    col = c("white", class_colors[classes]),
    breaks = seq(-0.5, K + 0.5, by = 1),
    RowSideColors = row_side,
    labRow = FALSE, labCol = FALSE, main = main,
    margins = c(3, 4), trace = "none", keysize = 0.5, key = FALSE
  )
}

# Relabel a strand-aware TSS-relative occupancy matrix ("<class>@<relpos>") into
# absolute genomic coordinates ("<class>@<genomic_pos>"), matching the coordinate
# system of the m6A methylation heatmap (genomic, not strand-flipped). Only the
# column names change; plot_footprint_cluster_heatmap re-sorts columns by the
# parsed position, so for minus-strand genes the genomic view is mirror-flipped
# relative to the TSS-relative view.
relabel_occ_to_genomic <- function(occ_mat, tss, strand) {
  feat <- colnames(occ_mat)
  cls  <- sub("@.*$", "", feat)
  rel  <- as.numeric(sub("^.*@", "", feat))
  gpos <- if (strand == "+") tss + rel else tss - rel
  colnames(occ_mat) <- paste0(cls, "@", gpos)
  occ_mat
}

# Assign each read to a cluster by k-means on the PCA of the topic-proportions
# matrix L (same recipe as the inline clustering in the notebook). Returns a
# factor "cluster1..clusterk" aligned to the rows of fit$L.
cluster_reads_by_topics <- function(fit, k, seed = 1) {
  set.seed(seed)
  pca <- prcomp(fit$L)$x
  cl  <- kmeans(pca, centers = k, iter.max = 100)$cluster
  cl  <- factor(cl)
  levels(cl) <- paste0("cluster", seq_along(levels(cl)))
  cl
}

# Per-read table combining the topic proportions (fit$L), the k-means cluster from
# cluster_reads_by_topics(), and the read metadata in rids_df. rownames(fit$L) are
# the post-filter RIDs, so rids_df is matched against them (no refiltering needed).
# Pass the same fit that was used for clustering in process_region_m6a() (the
# binomial one) to reproduce the saved cluster assignments.
read_topic_assignments <- function(fit, rids_df, seed = 1) {
  L <- as.matrix(fit$L)
  clusters <- cluster_reads_by_topics(fit, ncol(L), seed = seed)

  meta_cols <- intersect(c("sample_name", "chr", "start", "end", "strand"), colnames(rids_df))
  meta <- rids_df[match(rownames(L), rids_df$RID), meta_cols, drop = FALSE]

  df <- data.frame(
    RID            = rownames(L),
    meta,
    cluster        = as.character(clusters),
    dominant_topic = colnames(L)[max.col(L, ties.method = "first")],
    max_loading    = apply(L, 1, max),
    L,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  stopifnot(nrow(df) == nrow(L))
  df
}

# Topic profiles fit$F in long format: one row per (position, topic). Rownames of F
# are "<base>_<genomic_pos>" (e.g. "A_112836326"), split into `base` and `pos`.
topic_profiles_long <- function(fit) {
  Fm   <- as.matrix(fit$F)
  feat <- rep(rownames(Fm), times = ncol(Fm))
  data.frame(
    feature = feat,
    base    = sub("_.*$", "", feat),
    pos     = as.numeric(sub("^.*_", "", feat)),
    topic   = rep(colnames(Fm), each = nrow(Fm)),
    prob    = as.vector(Fm),
    stringsAsFactors = FALSE
  )
}

# Per-cluster aggregate footprint profiles: mean per-class occupancy at each
# TSS-relative position, one facet per cluster, one line per footprint class.
# `clusters` must be a factor aligned to the rows of occ_mat.
plot_footprint_cluster_profiles <- function(occ_mat, clusters, classes, class_colors,
                                            title = NULL,
                                            xlab = "Distance to TSS (bp, strand-oriented)") {
  feat   <- colnames(occ_mat)
  fclass <- sub("@.*$", "", feat)
  fpos   <- as.numeric(sub("^.*@", "", feat))
  M <- as.matrix(occ_mat)

  cl_levels <- levels(droplevels(clusters))
  parts <- list()
  for (cl in cl_levels) {
    mask <- which(clusters == cl)
    if (!length(mask)) next
    for (ci in seq_along(classes)) {
      cols <- which(fclass == classes[ci])
      if (!length(cols)) next
      parts[[length(parts) + 1]] <- data.frame(
        cluster = paste0(cl, " (n=", length(mask), ")"),
        class   = classes[ci],
        pos     = fpos[cols],
        value   = colMeans(M[mask, cols, drop = FALSE]),
        stringsAsFactors = FALSE
      )
    }
  }
  df <- do.call(rbind, parts)
  df$class <- factor(df$class, levels = names(class_colors))

  ggplot(df, aes(x = pos, y = value, color = class)) +
    geom_line() +
    scale_color_manual(values = class_colors, drop = FALSE) +
    facet_wrap(~ cluster, ncol = 1, strip.position = "right") +
    labs(x = xlab, y = "mean occupancy", title = title) +
    theme_cowplot(font_size = 10)
}
