# Shared input for the LCL Fourier and autocorrelation workflows.
# Source parsing_footprints_functions.r, leiden_manhattan_functions.r and
# LCL_phasing.r first. The binary representation is the original Fourier one;
# the read-selection rule is the existing LCL Leiden focal-heterozygote filter.

lcl_build_phased_m6a_input <- function(region, sample_table, ft_result_dir,
    phasing_root, existing_phase_dir = NULL) {
  required <- c("region_id", "chr", "start", "end", "analysis_start", "analysis_end",
    "focal_pos", "focal_snp", "ref", "alt")
  stopifnot(nrow(region) == 1L, all(required %in% names(region)),
    !anyNA(region[, required]), region$analysis_start >= 1L,
    region$focal_pos >= region$analysis_start, region$focal_pos <= region$analysis_end,
    region$ref %in% c("A", "C", "G", "T"), region$alt %in% c("A", "C", "G", "T"),
    region$ref != region$alt, !anyDuplicated(sample_table$sample_name))
  samples <- lcl_region_samples(region, sample_table)
  stopifnot(nrow(samples) > 0L, !anyNA(samples$sample_name))
  positions <- seq.int(region$analysis_start, region$analysis_end)
  paths <- file.path(ft_result_dir, samples$sample_name, "extracted_results/m6a_by_chr",
    paste0(samples$sample_name, ".ft_extracted_m6a.", region$chr, ".bed.gz"))
  dat <- assemble_region_m6a(sample_table = samples, region = region, full_span = TRUE,
    positions = positions, matrix_dir = NULL, m6a_paths = paths)
  # Keep full-span reads with at least one m6A call in the window, ordered by
  # sample and then read ID, with sample__RID row names.
  info <- dat$rids_df
  keep <- which(Matrix::rowSums(dat$met_mat) > 0)
  keep <- keep[order(match(info$sample_name[keep], samples$sample_name),
    info$original_RID[keep], method = "radix")]
  mat <- dat$met_mat[keep, , drop = FALSE]
  reads <- info[keep, c("chr", "start", "end", "strand", "original_RID", "sample_name")]
  reads <- data.frame(RID = paste(reads$sample_name, reads$original_RID, sep = "__"), reads)
  rownames(reads) <- NULL
  rownames(mat) <- reads$RID
  stopifnot(!anyNA(mat), all(mat@x %in% c(0, 1)), !anyDuplicated(reads$RID),
    identical(rownames(mat), reads$RID), identical(as.integer(colnames(mat)), positions))

  # Reuse a read-only cache only if it covers this actual region's input reads.
  phases <- setNames(lapply(samples$sample_name, function(sample_name) {
    path <- if (!is.null(existing_phase_dir))
      file.path(existing_phase_dir, paste0(sample_name, "_haplotags.rds")) else ""
    if (nzchar(path) && file.exists(path)) readRDS(path) else NULL
  }), samples$sample_name)
  refresh <- samples$sample_name[vapply(samples$sample_name, function(sample_name) {
    phase <- phases[[sample_name]]
    wanted <- reads$original_RID[reads$sample_name == sample_name]
    is.null(phase) || !isTRUE(phase$available) ||
      !all(wanted %in% phase$signature$reads) ||
      !identical(phase$signature$cell_line, samples$cell_line[match(sample_name, samples$sample_name)])
  }, logical(1))]
  if (length(refresh)) {
    temporary <- tempfile(pattern = "lcl_phased_input_")
    dir.create(temporary)
    on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
    phases[refresh] <- cache_lcl_haplotags(samples[samples$sample_name %in% refresh, , drop = FALSE],
      list(list(rids_df = reads)), phasing_root, temporary, reuse = FALSE)[refresh]
  }
  filtered <- lcl_filter_focal_heterozygotes(list(rids_df = reads, met_mat = mat),
    region, samples, phases)
  retained <- filtered$rids_df
  # The shared Leiden filter drops empty columns for distance clustering.
  # Fourier/ACF must preserve every original base, including all-zero columns.
  filtered$met_mat <- mat[retained$RID, , drop = FALSE]
  stopifnot(ncol(filtered$met_mat) == length(positions), !anyNA(filtered$met_mat),
    all(retained$haplotype %in% c("HP1", "HP2")),
    all(retained$focal_genotype %in% c("0|1", "1|0")),
    all(retained$allele_status == "phased_focal_genotype"))
  retained$row_id <- retained$RID
  retained$RID <- retained$original_RID
  retained$read_start <- retained$start
  retained$read_end <- retained$end
  filtered$read_meta <- retained
  filtered$phasing_sources <- lapply(phases, function(x) x[c("available", "reason", "paths", "signature")])
  message(region$region_id, ": retained ", nrow(retained), " / ", nrow(reads),
    " full-span reads after focal heterozygote phasing; ", length(positions), " bases")
  filtered
}
