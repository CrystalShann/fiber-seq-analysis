# Configured custom coordinates are 1-based inclusive.
select_custom_regions <- function(spec) {
  if (!length(spec)) stop("regions.custom has no intervals")
  rows <- lapply(spec, function(x) {
    required <- c("region_id", "chr", "start", "end")
    if (!all(required %in% names(x))) stop("Each custom region needs: ", paste(required, collapse = ", "))
    bounds <- as.numeric(c(x$start, x$end))
    if (length(bounds) != 2L || any(!is.finite(bounds)) || any(bounds != as.integer(bounds)) ||
        bounds[1] < 1L || bounds[2] < bounds[1]) stop("Invalid custom coordinates for ", x$region_id)
    mode <- if (is.null(x$samples)) "all" else x$samples
    if (!identical(mode, "all")) stop("Custom regions require samples: all")
    data.frame(region_id = x$region_id, chr = x$chr, start = as.integer(bounds[1]) - 1L,
               end = as.integer(bounds[2]), analysis_start = as.integer(bounds[1]),
               analysis_end = as.integer(bounds[2]), width = as.integer(diff(bounds)) + 1L,
               region_type = "custom", strand = "*", tss = NA_integer_,
               annotation = if (is.null(x$annotation)) x$region_id else x$annotation,
               anchor = mean(bounds), samples_mode = mode,
               coordinate_system = "BED: 0-based, half-open; analysis: 1-based inclusive",
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
