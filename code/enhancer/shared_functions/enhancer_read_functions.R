# Shared full-span enhancer fiber extraction at one column per genomic base.
# Callers load data.table, dplyr, GenomicRanges, Rsamtools and Matrix, and source
# code/topic_model/topic_modelling_functions.r before using this helper.

collect_enhancer_read_inputs <- function(regions, samples, ft_dir, reference) {
  # BED windows are [midpoint0 - 500, midpoint0 + 500); columns are individual bp.
  regions <- data.table::as.data.table(regions)
  stopifnot(all(c("enhancer_id", "enhancer_class", "chr", "midpoint0") %in% names(regions)),
            !anyNA(regions[, .(enhancer_id, enhancer_class, chr, midpoint0)]),
            !anyDuplicated(regions$enhancer_id), !anyDuplicated(samples))
  reference_lengths <- GenomeInfoDb::seqlengths(
    Rsamtools::scanFaIndex(Rsamtools::FaFile(reference)))
  qc <- data.table::data.table(
    enhancer_id = as.character(regions$enhancer_id),
    enhancer_class = as.character(regions$enhancer_class),
    status = rep("no_full_span_reads", nrow(regions)), n_reads = integer(nrow(regions)))
  metadata <- list()
  row_parts <- column_parts <- list()
  n_reads <- 0L
  part <- 0L
  for (i in seq_len(nrow(regions))) {
    if (i == 1L || i %% 250L == 0L || i == nrow(regions))
      message("Extracting enhancer ", i, "/", nrow(regions), " (", regions$enhancer_class[i], ")")
    region <- regions[i]
    left0 <- as.integer(region$midpoint0 - 500L)
    right0 <- as.integer(region$midpoint0 + 500L)
    if (left0 < 0L || !region$chr %in% names(reference_lengths) ||
        right0 > reference_lengths[[region$chr]]) {
      qc$status[i] <- "window_outside_reference"
      next
    }
    gr <- GenomicRanges::GRanges(region$chr, IRanges::IRanges(left0 + 1L, right0))
    for (sample_id in samples) {
      path <- file.path(ft_dir, sample_id, "extracted_results", "m6a_by_chr",
                        paste0(sample_id, ".ft_extracted_m6a.", region$chr, ".bed.gz"))
      if (!all(file.exists(c(path, paste0(path, ".tbi")))))
        stop("Missing m6A BED12 or tabix index: ", path)
      bed <- read_ft_bed12(path, region = gr, longest_alignment = TRUE)
      if (!nrow(bed)) next
      bed <- bed[bed$start <= left0 & bed$end >= right0, , drop = FALSE]
      if (!nrow(bed)) next
      original_rid <- as.character(bed$RID)
      bed$RID <- paste(region$enhancer_id, sample_id, original_rid, sep = "::")
      # The shared converter removes fibertools' first/last sentinel blocks.
      # Keep BED rows even when no genuine m6A calls remain after conversion.
      blocks <- convert_ft_bed12_to_bed6(bed)
      part <- part + 1L
      row_parts[[part]] <- column_parts[[part]] <- integer()
      if (nrow(blocks)) {
        if (any(blocks$end - blocks$start != 1L)) stop("Expected single-base m6A blocks")
        blocks <- blocks[blocks$start >= left0 & blocks$end <= right0, , drop = FALSE]
        blocks <- unique(blocks[, c("RID", "start"), drop = FALSE])
        if (nrow(blocks)) {
          row_index <- match(blocks$RID, bed$RID)
          if (anyNA(row_index)) stop("m6A call has no matching full-span read")
          row_parts[[part]] <- n_reads + row_index
          column_parts[[part]] <- as.integer(blocks$start - left0 + 1L)
        }
      }
      metadata[[part]] <- data.table::data.table(
        RID = as.character(bed$RID), read_id = original_rid,
        enhancer_id = as.character(region$enhancer_id),
        enhancer_class = as.character(region$enhancer_class),
        sample_name = as.character(sample_id), chr = as.character(bed$chr),
        start = as.integer(bed$start + 1L), end = as.integer(bed$end),
        strand = as.character(bed$strand))
      n_reads <- n_reads + nrow(bed)
      qc$n_reads[i] <- qc$n_reads[i] + nrow(bed)
      qc$status[i] <- "retained"
    }
  }
  metadata <- if (length(metadata)) data.table::rbindlist(metadata) else
    data.table::data.table(RID = character(), read_id = character(),
      enhancer_id = character(), enhancer_class = character(), sample_name = character(),
      chr = character(), start = integer(), end = integer(), strand = character())
  row_index <- as.integer(unlist(row_parts, use.names = FALSE))
  column_index <- as.integer(unlist(column_parts, use.names = FALSE))
  # Absence of a call is zero at each covered genomic base, including G/C.
  # No genomic positions are compressed, averaged, or removed.
  mat <- Matrix::sparseMatrix(i = row_index, j = column_index,
    x = rep(1, length(row_index)), dims = c(n_reads, 1000L),
    dimnames = list(metadata$RID, as.character(-500:499)))
  stopifnot(!anyDuplicated(metadata$RID), nrow(mat) == nrow(metadata),
            all(mat@x == 1), !anyNA(mat))
  list(mat = mat, metadata = metadata, qc = qc)
}

