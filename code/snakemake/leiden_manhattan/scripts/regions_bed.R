# Custom regions supplied as a BED file. BED is 0-based, half-open; the
# workflow stores analysis coordinates as 1-based, inclusive.
select_bed_regions <- function(spec) {
  path <- spec$bed
  if (is.null(path) || length(path) != 1L || !nzchar(path) || !file.exists(path))
    stop("regions.bed requires an existing BED file")
  bed <- read.delim(path, header = FALSE, stringsAsFactors = FALSE,
                    quote = "", comment.char = "", check.names = FALSE)
  if (ncol(bed) < 3L) stop("BED file needs at least 3 columns: ", path)
  names(bed)[1:3] <- c("chr", "start0", "end0")
  start0 <- suppressWarnings(as.numeric(bed$start0))
  end0 <- suppressWarnings(as.numeric(bed$end0))
  if (any(!is.finite(start0)) || any(!is.finite(end0)) ||
      any(start0 != as.integer(start0)) || any(end0 != as.integer(end0)) ||
      any(start0 < 0L) || any(end0 <= start0))
    stop("Invalid BED coordinates: ", path)
  ids <- if (ncol(bed) >= 4L) as.character(bed[[4]]) else rep(NA_character_, nrow(bed))
  missing_id <- is.na(ids) | !nzchar(ids) | ids == "."
  ids[missing_id] <- paste(as.character(bed$chr[missing_id]),
                           as.integer(start0[missing_id]) + 1L,
                           as.integer(end0[missing_id]), sep = "_")
  if (!nrow(bed) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("BED region IDs must be nonempty and unique: ", path)
  annotation <- if (ncol(bed) >= 5L) as.character(bed[[5]]) else ids
  data.frame(region_id = ids, chr = as.character(bed$chr),
    start = as.integer(start0), end = as.integer(end0),
    analysis_start = as.integer(start0) + 1L, analysis_end = as.integer(end0),
    width = as.integer(end0 - start0), region_type = "bed",
    strand = if (ncol(bed) >= 6L) as.character(bed[[6]]) else "*",
    tss = NA_integer_, annotation = annotation, anchor = (start0 + end0) / 2,
    samples_mode = "all",
    coordinate_system = "BED: 0-based, half-open; analysis: 1-based inclusive",
    stringsAsFactors = FALSE, row.names = NULL)
}
