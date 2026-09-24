library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)
library(SummarizedExperiment)
library(BiocParallel)
library(Rsamtools)
library(footprintR)
library(ggplot2)
library(patchwork)
library(scales)

assert_file_exists <- function(path, label = "file") {
  if (!file.exists(path)) {
    stop(sprintf("Missing %s: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

write_gz_tsv <- function(x, path) {
  data.table::fwrite(x, file = path, sep = "\t")
  invisible(path)
}

coerce_iranges_coord <- function(x, label = "coordinate") {
  if (inherits(x, "integer64")) {
    x <- as.integer(x)
  } else if (is.character(x) || is.factor(x)) {
    x <- suppressWarnings(as.integer(as.character(x)))
  } else if (!is.integer(x)) {
    x <- as.integer(x)
  }

  if (anyNA(x)) {
    stop(sprintf("Non-numeric %s values cannot be used to build IRanges.", label), call. = FALSE)
  }

  x
}

read_bam_seqinfo <- function(bam_path) {
  targets <- Rsamtools::scanBamHeader(bam_path)[[1]]$targets
  GenomeInfoDb::Seqinfo(
    seqnames = names(targets),
    seqlengths = as.numeric(targets)
  )
}

read_gencode_tss <- function(path, chrom = NULL, window_bp = 2000L, seqinfo = NULL) {
  dt <- data.table::fread(
    path,
    header = FALSE,
    col.names = c("chromosome", "start0", "end0", "tss_id", "score", "tss_strand")
  )

  if (!is.null(chrom)) {
    dt <- dt[chromosome == chrom]
  }

  if (!nrow(dt)) {
    stop("No TSS rows remained after chromosome filtering.", call. = FALSE)
  }

  dt[, `:=`(
    chromosome = as.character(chromosome),
    start0 = coerce_iranges_coord(start0, "TSS start"),
    end0 = coerce_iranges_coord(end0, "TSS end"),
    tss_strand = as.character(tss_strand)
  )]

  dt[, row_id := .I]
  dt[, tss_coordinate := coerce_iranges_coord(start0 + 1L, "TSS coordinate")]
  dt[, gene_name_or_tss_id := tss_id]
  dt[, tss_uid := paste(chromosome, tss_coordinate, tss_strand, tss_id, row_id, sep = "|")]

  tss_gr <- GenomicRanges::GRanges(
    seqnames = dt$chromosome,
    ranges = IRanges::IRanges(
      start = coerce_iranges_coord(dt$tss_coordinate, "TSS coordinate"),
      end = coerce_iranges_coord(dt$tss_coordinate, "TSS coordinate")
    ),
    strand = dt$tss_strand
  )

  if (!is.null(seqinfo)) {
    GenomeInfoDb::seqinfo(tss_gr) <- seqinfo[unique(dt$chromosome)]
  }

  mcols(tss_gr)$tss_uid <- dt$tss_uid
  mcols(tss_gr)$tss_id <- dt$tss_id
  mcols(tss_gr)$gene_name_or_tss_id <- dt$gene_name_or_tss_id

  promoter_windows <- GenomicRanges::promoters(
    tss_gr,
    upstream = window_bp,
    downstream = window_bp + 1L
  )

  promoter_windows <- GenomicRanges::trim(promoter_windows)

  promoter_meta <- data.table::data.table(
    tss_uid = dt$tss_uid,
    chromosome = dt$chromosome,
    tss_coordinate = dt$tss_coordinate,
    tss_strand = dt$tss_strand,
    gene_name_or_tss_id = dt$gene_name_or_tss_id,
    tss_id = dt$tss_id,
    promoter_start = start(promoter_windows),
    promoter_end = end(promoter_windows)
  )

  tss_table <- copy(promoter_meta)
  tss_table[, strand := tss_strand]
  tss_table[, tss_strand := NULL]
  data.table::setcolorder(
    tss_table,
    c(
      "chromosome",
      "tss_coordinate",
      "strand",
      "gene_name_or_tss_id",
      "tss_id",
      "tss_uid",
      "promoter_start",
      "promoter_end"
    )
  )

  list(
    tss_table = tss_table,
    tss_gr = tss_gr,
    promoter_windows = promoter_windows,
    promoter_meta = promoter_meta
  )
}

read_fiber_bed12 <- function(path) {
  data.table::fread(
    path,
    header = FALSE,
    col.names = c(
      "chrom",
      "chromStart",
      "chromEnd",
      "name",
      "score",
      "strand",
      "thickStart",
      "thickEnd",
      "itemRgb",
      "blockCount",
      "blockSizes",
      "blockStarts"
    )
  )
}

parse_bed12_blocks <- function(block_sizes, block_starts, drop_terminal_blocks = TRUE) {
  sizes <- as.integer(strsplit(block_sizes, ",", fixed = TRUE)[[1]])
  starts <- as.integer(strsplit(block_starts, ",", fixed = TRUE)[[1]])

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

expand_ft_bed12_features <- function(bed12, feature_name = "feature", drop_terminal_blocks = TRUE) {
  feature_name_local <- feature_name
  
  # { ... } = run once per row/read
  # for each group, expand into multiple rows and return a new table using .()
  bed12[, {
    parsed <- parse_bed12_blocks(
      # block size = length of 1
      block_sizes = blockSizes, # a compressed string representation of multiple blocks for a single read
      block_starts = blockStarts,
      drop_terminal_blocks = drop_terminal_blocks
    )

    sizes <- parsed$sizes
    starts0 <- parsed$starts

    if (!length(sizes)) {
      return(NULL)
    }

    feature_start0 <- chromStart + starts0
    feature_end0 <- chromStart + starts0 + sizes
    feature_mid0 <- as.integer((feature_start0 + feature_end0) / 2)
 
    # .() is shorthand for list()
    .(
      feature_chromosome = chrom,
      read_id = name,
      read_strand = strand,
      feature_start0 = feature_start0,
      feature_end0 = feature_end0,
      feature_mid0 = feature_mid0,
      feature_pos = feature_mid0 + 1L,
      feature_size = sizes,
      feature_name = feature_name_local
    )
  }, by = .(chrom, chromStart, chromEnd, name, strand)]
}

# Converts feature table into GRanges
feature_blocks_to_granges <- function(feature_dt) {
  GenomicRanges::GRanges(
    seqnames = feature_dt$feature_chromosome,
    ranges = IRanges::IRanges(
      start = coerce_iranges_coord(feature_dt$feature_start0 + 1L, "feature start"),
      end = coerce_iranges_coord(feature_dt$feature_end0, "feature end")
    ),
    strand = "*"
  )
}


# Finds which features overlap TSS promoter windows and computes relative positions

## feature_dt: a data.table containing feature information, like chromosome, read ID, feature position, feature size
## feature_gr: the same features, but stored as a GRanges object for overlap operations
## promoter_windows: genomic promoter intervals around TSSs, also as GRanges
## promoter_meta: a table containing metadata for those promoter windows, such as TSS coordinate, strand, gene/TSS ID

align_feature_blocks_to_tss <- function(feature_dt, feature_gr, promoter_windows, promoter_meta) {
  if (!nrow(feature_dt) || !length(feature_gr) || !length(promoter_windows)) {
    # if any of the input is empty, return an empty data table instead of giving errors
    return(data.table::data.table())
  }

  # which feature block overlap which promoter windows?
  hits <- GenomicRanges::findOverlaps(
    # set query ranges to be the feature ranges --> intervals want to test for overlap
    query = feature_gr, 
    # subject ranges = promoter windows --> reference interval
    subject = promoter_windows,
    # overlap is only based on genommic coordinate
    ignore.strand = TRUE
  )

  if (!length(hits)) {
    return(data.table::data.table())
  }

  aligned <- cbind(
    # subset feature meta table to only include overlapping features
    feature_dt[S4Vectors::queryHits(hits)], 
    # subset promoter meta table to only include overlapping promoter windows
    promoter_meta[S4Vectors::subjectHits(hits)]
  )

  # Create a relative-position column based on strand
  ## add a new column to indicate relative position to TSS
  aligned[, rel_pos := ifelse(
    tss_strand == "+",
    # if + strand, relative position = feature - TSS coordinate --> positive number = donwstream
    feature_pos - tss_coordinate,
    tss_coordinate - feature_pos
  )]
  # return aligned data table after modifications
  aligned[]
}


# summarize nucleosome count, size and distance to TSS for each unique TSS/read
summarize_nucleosomes_by_tss <- function(nuc_aligned, window_bp = 2000L) {
  # takes output from  <align_feature_blocks_to_tss> to construct nuc_aligned
  if (!nrow(nuc_aligned)) {
    return(data.table::data.table(tss_uid = character()))
  }

  nuc_aligned[, {
    # .N = number of rows in the group
    n_nuc <- .N
    .(
      # n_nuc = nucleosome count per TSS
      structure_nucleosome_count = n_nuc,
      # average nucleosome size in this TSS group
      structure_mean_nucleosome_size = mean(feature_size),
      # Computes the median nucleosome size per TSS
      structure_median_nucleosome_size = as.double(median(feature_size)),
      # How many nucleosomes per kilobase are present around this TSS
      structure_nucleosome_density_per_kb = n_nuc / ((2L * window_bp + 1L) / 1000),
      # Counts how many nucleosomes have midpoint positions within 150 bp of the TSS
      structure_core_nucleosome_count = sum(abs(rel_pos) <= 150L),
      # Computes the average signed midpoint offset from the TSS.
      ## On average, are nucleosome midpoints shifted upstream or downstream of the TSS?
      structure_mean_midpoint_offset = mean(rel_pos),
      # How far away, on average, are nucleosome midpoints from the TSS?
      structure_mean_abs_midpoint_offset = mean(abs(rel_pos))
    )
  }, by = tss_uid]
}

## modbase_code == "a" ~ "m6A",
## modbase_code == "m" ~ "5mC"

# function that incorporates <footprintR::readModBam> to count how many bases were observed and how many were modified at each genomic position for this sample, then return the results as a clean table
read_modbam_summary_positions <- function(
  bamfiles,
  sample_name,
  promoter_windows,
  modbase_code,
  seqinfo, # A chromosome metadata object containing chromosome names and lengths
  bpparam,
  modProbThreshold = 0.9
) {
  se_summary <- footprintR::readModBam(
    bamfiles = bamfiles,
    regions = promoter_windows,
    modbase = setNames(modbase_code, sample_name), # Creates a named vector where the modification code is assigned to the sample name
    # aggregated counts at each genomic position
    ## count total valid bases per position/strand
    ## count modified bases per position/strand
    ## return assays Nvalid, Nmod, and FracM
    level = "summary",
    modProbThreshold = modProbThreshold,
    seqinfo = seqinfo,
    trim = TRUE,
    BPPARAM = bpparam,
    verbose = TRUE
  )
  
  # rows = genomic positions
  # columns = samples
  # assays = count summaries for those positions
  rr_gr <- SummarizedExperiment::rowRanges(se_summary)
  rr <- data.table::data.table(
    call_chromosome = as.character(GenomeInfoDb::seqnames(rr_gr)),
    position = coerce_iranges_coord(IRanges::start(rr_gr), "modBAM position"),
    call_strand = as.character(BiocGenerics::strand(rr_gr))
  )
  
  # Pulls the "Nvalid" assay from the SummarizedExperiment
  ## Nvalid is the total number of valid/callable bases at each position for each sample
  total_calls <- as.integer(SummarizedExperiment::assay(se_summary, "Nvalid")[, sample_name])
  
  ## Nmod is the number of bases classified as modified at each position.
  modified_calls <- as.integer(SummarizedExperiment::assay(se_summary, "Nmod")[, sample_name])

  rr[, sample_name := sample_name]
  rr[, modbase_code := modbase_code]
  rr[, total_calls := total_calls]
  rr[, modified_calls := modified_calls]
  # Computes how many calls were not modified
  rr[, unmodified_calls := total_calls - modified_calls]
  rr[, frac_modified := data.table::fifelse(total_calls > 0, modified_calls / total_calls, NA_real_)]

  # only rows where there was at least one valid call
  rr[total_calls > 0]
}


# converts per-position modified vs unmodified counts into a long-format table with a status label and n_calls column
make_modified_status_table <- function(position_dt) {
  if (!nrow(position_dt)) {
    return(data.table::data.table())
  }

  modified_dt <- copy(position_dt)
  modified_dt[, `:=`(status = "modified", n_calls = modified_calls)]

  unmodified_dt <- copy(position_dt)
  unmodified_dt[, `:=`(status = "not_modified", n_calls = unmodified_calls)]

  data.table::rbindlist(
    list(modified_dt, unmodified_dt),
    use.names = TRUE,
    fill = TRUE
  )[n_calls > 0]
}

align_positions_to_tss <- function(position_dt, promoter_windows, promoter_meta) {
  if (!nrow(position_dt) || !length(promoter_windows)) {
    return(data.table::data.table())
  }

  pos_gr <- GenomicRanges::GRanges(
    seqnames = position_dt$call_chromosome,
    ranges = IRanges::IRanges(
      start = coerce_iranges_coord(position_dt$position, "modBAM position"),
      end = coerce_iranges_coord(position_dt$position, "modBAM position")
    ),
    strand = "*"
  )

  hits <- GenomicRanges::findOverlaps(
    query = pos_gr,
    subject = promoter_windows,
    ignore.strand = TRUE
  )

  if (!length(hits)) {
    return(data.table::data.table())
  }

  aligned <- cbind(
    position_dt[S4Vectors::queryHits(hits)],
    promoter_meta[S4Vectors::subjectHits(hits)]
  )

  aligned[, rel_pos := ifelse(
    tss_strand == "+",
    position - tss_coordinate,
    tss_coordinate - position
  )]

  aligned[]
}

summarize_modality_by_tss <- function(aligned_dt, prefix = c("accessibility", "methylation")) {
  prefix <- match.arg(prefix)

  if (!nrow(aligned_dt)) {
    return(data.table::data.table(tss_uid = character()))
  }

  out <- aligned_dt[, {
    total <- sum(total_calls)
    modified <- sum(modified_calls)
    unmodified <- sum(unmodified_calls)
    .(
      modified_calls = modified,
      unmodified_calls = unmodified,
      total_calls = total,
      positions_with_coverage = sum(total_calls > 0),
      fraction_modified = if (total > 0) modified / total else NA_real_
    )
  }, by = tss_uid]

  data.table::setnames(
    out,
    old = c("modified_calls", "unmodified_calls", "total_calls", "positions_with_coverage", "fraction_modified"),
    new = c(
      sprintf("%s_modified_calls", prefix),
      sprintf("%s_unmodified_calls", prefix),
      sprintf("%s_total_calls", prefix),
      sprintf("%s_positions_with_coverage", prefix),
      sprintf("%s_fraction_modified", prefix)
    )
  )

  out[]
}

make_nucleosome_meta <- function(nuc_aligned, num_tss, window_bp = 2000L, bin_size = 10L) {
  bins <- data.table::data.table(bin = seq(-window_bp, window_bp, by = bin_size))

  if (!nrow(nuc_aligned)) {
    bins[, nucleosome_midpoints_per_tss := 0]
    return(bins[])
  }

  tmp <- copy(nuc_aligned)
  tmp[, bin := round(rel_pos / bin_size) * bin_size]

  meta <- tmp[, .N, by = bin]
  meta <- merge(bins, meta, by = "bin", all.x = TRUE, sort = TRUE)
  meta[is.na(N), N := 0]
  meta[, nucleosome_midpoints_per_tss := N / num_tss]
  meta[]
}

# takes modification calls already aligned relative to TSSs, groups them into position bins, summarizes modified and total calls per bin, computes fraction modified and modified calls per TSS

make_modification_meta <- function(aligned_dt, num_tss, label, window_bp = 2000L, bin_size = 10L) {
  bins <- data.table::data.table(bin = seq(-window_bp, window_bp, by = bin_size))
  label_value <- label

  if (!nrow(aligned_dt)) {
    bins[, `:=`(
      modified_calls = 0L,
      total_calls = 0L,
      fraction_modified = NA_real_,
      modified_calls_per_tss = 0,
      label = label_value
    )]
    return(bins[])
  }

  tmp <- copy(aligned_dt)
  # converts each relative position rel_pos into the nearest bin value
  tmp[, bin := round(rel_pos / bin_size) * bin_size]

  meta <- tmp[, .(
    modified_calls = sum(modified_calls),
    total_calls = sum(total_calls)
  ), by = bin]

  meta <- merge(bins, meta, by = "bin", all.x = TRUE, sort = TRUE)
  meta[is.na(modified_calls), modified_calls := 0L]
  meta[is.na(total_calls), total_calls := 0L]
  meta[, fraction_modified := data.table::fifelse(total_calls > 0, modified_calls / total_calls, NA_real_)]
  meta[, modified_calls_per_tss := modified_calls / num_tss]
  meta[, label := label_value]
  meta[]
}

plot_structure_meta <- function(
  meta_dt,
  sample_name,
  chrom,
  window_bp        = 5000L,   # x-axis range (+/- bp from TSS)
  x_break_interval = 500,     # spacing of x-axis tick marks (bp)
  loess_span       = 0.06,    # LOESS smoothing bandwidth (smaller = tighter fit)
  ylim             = NULL,    # c(lo, hi) to fix y-axis, NULL for auto
  line_color       = "#4D4D4D",
  smooth_color     = "#1f77b4",
  linewidth        = 0.6,     # raw data line width
  smooth_linewidth = 0.9,     # LOESS line width
  base_size        = 12       # ggplot theme base font size
) {
  p <- ggplot2::ggplot(meta_dt, ggplot2::aes(x = bin, y = nucleosome_midpoints_per_tss)) +
    ggplot2::geom_line(color = line_color, linewidth = linewidth) +
    ggplot2::geom_smooth(
      method = "loess",
      span   = loess_span,
      se     = FALSE,
      color  = smooth_color,
      linewidth = smooth_linewidth
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(-window_bp, window_bp, x_break_interval),
      labels = function(x) paste0(x / 1000, " kb")
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
    ggplot2::labs(
      x     = "Distance to TSS",
      y     = "Nucleosome midpoints per TSS",
      title = sprintf("Promoter structure: %s %s", sample_name, chrom)
    ) +
    ggplot2::theme_bw(base_size = base_size)

  if (!is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = ylim)
  }
  p
}

plot_modification_meta <- function(
  meta_dt,
  title,
  color,
  window_bp        = 5000L,   # x-axis range (+/- bp from TSS)
  x_break_interval = 500,     # spacing of x-axis tick marks (bp)
  loess_span       = 0.06,    # LOESS smoothing bandwidth (smaller = tighter fit)
  ylim             = NULL,    # c(lo, hi) to fix y-axis, NULL for auto
  linewidth        = 0.6,     # raw data line width
  smooth_linewidth = 0.9,     # LOESS line width
  y_label          = "Modified fraction",
  base_size        = 12       # ggplot theme base font size
) {
  p <- ggplot2::ggplot(meta_dt, ggplot2::aes(x = bin, y = fraction_modified)) +
    ggplot2::geom_line(color = color, linewidth = linewidth, na.rm = TRUE) +
    ggplot2::geom_smooth(
      method = "loess",
      span   = loess_span,
      se     = FALSE,
      color  = color,
      linewidth = smooth_linewidth,
      na.rm  = TRUE
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(-window_bp, window_bp, x_break_interval),
      labels = function(x) paste0(x / 1000, " kb")
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
    ggplot2::labs(
      x     = "Distance to TSS",
      y     = y_label,
      title = title
    ) +
    ggplot2::theme_bw(base_size = base_size)

  if (!is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = ylim)
  }
  p
}

# Adapted from /project/spott/cshan/code/fiberseqFunctions.R

## takes a SummarizedExperiment containing read-level modification probabilities, converts those probabilities into 0/1 modified calls using a threshold, reshapes the data into a long table with one row per genomic position per read, and annotates each row with chromosome, position, strand, sample, and modification label
getBinarizedModificationsLongTable <- function(footprintRObject, sampleName, modProbThreshold = 0.9) {
  if (!inherits(footprintRObject, "SummarizedExperiment")) {
    stop("footprintRObject must be a SummarizedExperiment.", call. = FALSE)
  }

  # default is at per read level
  ## Extracts the assay named "mod_prob" from the SummarizedExperiment
  mod_prob_assay <- SummarizedExperiment::assay(footprintRObject, "mod_prob")[[sampleName]]
  mod_prob_matrix <- as(mod_prob_assay, "matrix")

  # 1 = probability greater than threshold = modified 
  mod_prob_matrix[mod_prob_matrix > modProbThreshold] <- 1
  mod_prob_matrix[mod_prob_matrix <= modProbThreshold] <- 0

  rownames(mod_prob_matrix) <- rownames(footprintRObject)
  mod_wide <- tibble::as_tibble(mod_prob_matrix, rownames = "genomic_location")

  # Reshapes from wide to long format
  mod_long <- tidyr::pivot_longer(
    mod_wide,
    cols = -genomic_location,
    names_to = "read_id",
    values_to = "mod_status"
  )

  # Parse genomic location and annotate call state
  mod_long <- dplyr::mutate(
    mod_long,
    call_chromosome = stringr::str_extract(genomic_location, "^[^:]+"),
    position = as.integer(stringr::str_extract(genomic_location, "(?<=:)[0-9]+")),
    call_strand = stringr::str_extract(genomic_location, "[+-]$"),
    sample_name = sampleName,
    call_state = ifelse(mod_status == 1, "modified", "not_modified")
  )

  modbase_code <- SummarizedExperiment::colData(footprintRObject)$modbase[[sampleName]]

  mod_long <- dplyr::mutate(
    mod_long,
    modbase_code = modbase_code,
    modbase_label = dplyr::case_when(
      modbase_code == "a" ~ "m6A",
      modbase_code == "m" ~ "5mC",
      TRUE ~ modbase_code
    )
  )

  dplyr::select(
    mod_long,
    sample_name,
    modbase_code,
    modbase_label,
    read_id,
    call_chromosome,
    position,
    call_strand,
    mod_status,
    call_state
  )
}

# Make a vector with 20 positive weights, then 120 negative weights, then 20 positive weights
build_nucleosome_weight_vector <- function() {
  c(
    rep(0.375, 20),
    rep(-0.125, 120),
    rep(0.375, 20)
  )
}

make_chapter5_example_plot <- function(
  bamfiles,
  sample_name,
  seqinfo,
  promoter_windows,
  promoter_meta,
  mod_prob_threshold,
  example_tss_index,
  output_dir,
  chrom
) {
  if (length(promoter_windows) < example_tss_index) {
    return(NULL)
  }

  example_region <- promoter_windows[example_tss_index]
  example_tss_uid <- promoter_meta$tss_uid[example_tss_index]

  se_example <- footprintR::readModBam(
    bamfiles = bamfiles,
    regions = example_region,
    modbase = setNames("a", sample_name),
    level = "quickread",
    seqinfo = seqinfo,
    trim = TRUE,
    BPPARAM = BiocParallel::SerialParam(),
    verbose = TRUE
  )

  se_example <- footprintR::flattenReadLevelAssay(se_example)
  wgt <- build_nucleosome_weight_vector()
  se_example <- footprintR::addFootprints(
    se = se_example,
    wgt = wgt,
    thresh = 0.02,
    name = "nucl"
  )

  example_binary_calls <- getBinarizedModificationsLongTable(
    footprintRObject = se_example,
    sampleName = sample_name,
    modProbThreshold = mod_prob_threshold
  )

  example_binary_calls <- data.table::as.data.table(example_binary_calls)
  example_binary_calls[, example_tss_uid := example_tss_uid]
  write_gz_tsv(
    example_binary_calls,
    file.path(
      output_dir,
      sprintf("%s_%s_exampleTSS_read_level_m6A_binary_calls.tsv.gz", sample_name, chrom)
    )
  )

  p <- footprintR::plotRegion(
    se_example,
    region = example_region,
    minCoveredFraction = 0,
    tracks = list(
      list(
        trackData = "mod_prob",
        trackType = "Heatmap",
        interpolate = TRUE,
        footprintColumns = "nucl",
        arglistFootprints = list(
          nucl = list(color = scales::alpha("firebrick", 0.5))
        ),
        orderReads = "squish"
      ),
      list(
        trackData = "FracMod",
        trackType = "Smooth",
        smoothMethod = "rollingMean",
        windowSize = 41
      )
    )
  ) + patchwork::plot_layout(heights = c(2, 0.6))

  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      sprintf("%s_%s_exampleTSS_footprintR_nucleosome_plot.pdf", sample_name, chrom)
    ),
    plot = p,
    width = 10,
    height = 6
  )

  list(
    plot = p,
    binary_calls = example_binary_calls,
    example_tss_uid = example_tss_uid
  )
}
