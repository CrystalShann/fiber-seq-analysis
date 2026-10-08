# Preserve the existing selector's complete region and focal-SNP metadata.
select_asfire_top_regions <- function(spec) {
  required <- c("asfire_file", "top_n", "width_bp", "q_cutoff", "min_coverage", "min_fire_freq")
  if (!all(required %in% names(spec))) stop("AS-FIRE configuration needs: ", paste(required, collapse = ", "))
  result <- select_lcl_top_asfire_regions(readRDS(spec$asfire_file),
    top_n = spec$top_n, width_bp = spec$width_bp, q_cutoff = spec$q_cutoff,
    min_coverage = spec$min_coverage, min_fire_freq = spec$min_fire_freq)
  mode <- if (is.null(spec$samples)) "contributing" else spec$samples
  if (!mode %in% c("all", "contributing")) stop("AS-FIRE samples must be all or contributing")
  result$samples_mode <- mode
  result
}
