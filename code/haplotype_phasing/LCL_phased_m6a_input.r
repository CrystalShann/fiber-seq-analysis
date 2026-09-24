# Shared input for the LCL Fourier and autocorrelation workflows.
# Source topic_modelling_functions.r, leiden_manhattan_functions.r and
# LCL_phasing.r first. The binary representation is the original Fourier one;
# the read-selection rule is the existing LCL Leiden focal-heterozygote filter.

lcl_build_phased_m6a_input <- function(region, sample_table, ft_result_dir,
    phasing_root, existing_phase_dir = NULL) {
  required <- c("region_id", "chr", "analysis_start", "analysis_end",
    "focal_pos", "focal_snp", "ref", "alt")
  stopifnot(nrow(region) == 1L, all(required %in% names(region)),
    !anyNA(region[, required]), region$analysis_start >= 1L,
    region$focal_pos >= region$analysis_start, region$focal_pos <= region$analysis_end,
    region$ref %in% c("A", "C", "G", "T"), region$alt %in% c("A", "C", "G", "T"),
    region$ref != region$alt, !anyDuplicated(sample_table$sample_name))
  samples <- lcl_region_samples(region, sample_table)
  stopifnot(nrow(samples) > 0L, !anyNA(samples$sample_name))
  positions <- seq.int(region$analysis_start, region$analysis_end)
  query <- GenomicRanges::GRanges(region$chr,
    IRanges::IRanges(region$analysis_start, region$analysis_end))
  matrices <- metadata <- list()
  for (sample_name in samples$sample_name) {
    path <- file.path(ft_result_dir, sample_name, "extracted_results/m6a_by_chr",
      paste0(sample_name, ".ft_extracted_m6a.", region$chr, ".bed.gz"))
    if (!all(file.exists(c(path, paste0(path, ".tbi"))))) stop("Missing indexed m6A BED: ", path)
    reads <- extract_ft_region_reads(path, query)
    if (is.null(reads) || !nrow(reads)) next
    info <- extract_ft_read_info(reads)
    if (is.null(info) || !nrow(info)) next
    mat <- get_sparse_met_mat(reads = reads, rids_df = info,
      window_start = region$analysis_start, window_end = region$analysis_end,
      all_met_pos = positions, base = "A")
    if (is.null(mat) || !nrow(mat)) next
    mat <- mat[rowSums(is.na(mat)) == 0L, , drop = FALSE]
    if (!nrow(mat)) next
    info <- info[match(rownames(mat), info$RID), , drop = FALSE]
    info$original_RID <- info$RID
    info$sample_name <- sample_name
    info$RID <- paste(sample_name, info$original_RID, sep = "__")
    rownames(mat) <- info$RID
    stopifnot(!anyNA(info$RID), all(info$start <= region$analysis_start),
      all(info$end >= region$analysis_end))
    matrices[[sample_name]] <- mat
    metadata[[sample_name]] <- info
  }
  if (!length(matrices)) stop("No fully spanning reads: ", region$region_id)
  mat <- do.call(rbind, matrices)
  reads <- dplyr::bind_rows(metadata)
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
    # Give the unchanged LCL cache loader a guard scoped to this disposable
    # metadata directory. No methylation/feature matrices are written.
    lcl_output_path <- function(path) {
      path <- normalizePath(path, mustWork = FALSE)
      root <- normalizePath(temporary)
      if (!(identical(path, root) || startsWith(path, paste0(root, "/")))) stop("Invalid temporary phasing path")
      path
    }
    loader <- cache_lcl_haplotags
    environment(loader) <- environment()
    phases[refresh] <- loader(samples[samples$sample_name %in% refresh, , drop = FALSE],
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
