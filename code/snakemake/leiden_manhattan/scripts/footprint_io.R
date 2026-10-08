# Config-driven footprint I/O. Coordinates in records remain BED half-open.
empty_footprints <- function() {
  data.frame(RID = character(), original_RID = character(), chr = character(),
             start = integer(), end = integer(), size = integer(), track = character(),
             sample_name = character(), stringsAsFactors = FALSE)
}

footprint_size_bound <- function(spec, defaults, field) {
  value <- if (field %in% names(spec)) spec[[field]] else defaults[[field]]
  if (is.null(value)) return(NULL)
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
      value < 0 || value != as.integer(value)) stop("footprints.nucleosome.", field, " must be null or a nonnegative integer")
  as.integer(value)
}

footprint_track_spec <- function(ds, cfg) {
  tracks <- character()
  colors <- c(m6A = "#800080")
  labels <- c(m6A = "m6A")
  nuc_track <- NULL
  minimum <- maximum <- NULL
  if (!is.null(ds$footprints$nucleosome)) {
    nuc <- ds$footprints$nucleosome
    minimum <- footprint_size_bound(nuc, cfg$nucleosome_size_defaults, "min_size")
    maximum <- footprint_size_bound(nuc, cfg$nucleosome_size_defaults, "max_size")
    if (!is.null(minimum) && !is.null(maximum) && minimum > maximum)
      stop("Nucleosome min_size exceeds max_size")
    nuc_track <- if (is.null(minimum) && is.null(maximum)) "nuc_all" else
      if (is.null(maximum)) paste0("nuc_gt", minimum - 1L, "bp") else
      if (is.null(minimum)) paste0("nuc_le", maximum, "bp") else
        paste0("nuc_", minimum, "-", maximum, "bp")
    tracks <- nuc_track
    colors[nuc_track] <- "#4d4d4d"
    labels[nuc_track] <- nuc_track
  }
  sizes <- as.character(unlist(ds$footprints$tf$sizes, use.names = FALSE))
  if (length(sizes)) {
    if (anyNA(sizes) || any(!grepl("^[0-9]+-[0-9]+$", sizes)) || anyDuplicated(sizes))
      stop("footprints.tf.sizes must contain unique min-max size strings")
    tf_tracks <- paste0("TF_", sizes, "bp")
    tf_colors <- c("#fdae6b", "#f16913", "#a63603")
    if (length(sizes) > length(tf_colors))
      tf_colors <- grDevices::colorRampPalette(tf_colors)(length(sizes))
    tracks <- c(tracks, tf_tracks)
    colors[tf_tracks] <- tf_colors[seq_along(sizes)]
    labels[tf_tracks] <- tf_tracks
  }
  list(tracks = tracks, colors = colors, labels = labels, nucleosome_track = nuc_track,
       min_size = minimum, max_size = maximum)
}

collect_footprints <- function(ds, cfg, sample_table, region, assignments) {
  track_spec <- footprint_track_spec(ds, cfg)
  needed <- c("RID", "original_RID", "sample_name")
  if (!all(needed %in% names(assignments))) stop("Footprint read table lacks: ", paste(setdiff(needed, names(assignments)), collapse = ", "))
  if (anyNA(assignments[, needed, drop = FALSE]) || anyDuplicated(assignments$RID))
    stop("Footprint read table has missing IDs or duplicate RID values")
  pieces <- tasks <- list()
  if (!is.null(ds$footprints$nucleosome)) {
    spec <- ds$footprints$nucleosome
    records <- extract_nucleosomes(sample_table, region, assignments,
      min_size = track_spec$min_size, max_size = track_spec$max_size,
      input_path = function(sample, chromosome) {
        values <- as.list(sample)
        for (field in names(region)) values[[field]] <- region[[field]][1L]
        values$chr <- chromosome
        expand_template(spec$source, values)
      }, format = spec$format, strict_blocks = TRUE)
    if (nrow(records)) {
      records$sample_name <- as.character(assignments$sample_name[match(records$RID, assignments$RID)])
      records$track <- track_spec$nucleosome_track
      pieces[[length(pieces) + 1L]] <- records[, names(empty_footprints()), drop = FALSE]
    }
  }
  for (size in as.character(unlist(ds$footprints$tf$sizes, use.names = FALSE)))
    tasks[[length(tasks) + 1L]] <- list(spec = ds$footprints$tf,
      track = paste0("TF_", size, "bp"), size = size)
  for (task in tasks) {
    spec <- task$spec
    footprint_format_columns(spec$format)
    if (is.null(spec$source) || length(spec$source) != 1L) stop("Footprint source must be one path template")
    pooled <- spec$format == "bed4_pooled"
    if (pooled && anyDuplicated(assignments$original_RID))
      stop("Combined BED4 lacks sample IDs: ambiguous original read names in ", region$region_id)
    indices <- if (pooled) 0L else seq_len(nrow(sample_table))
    for (sample_index in indices) {
      reads <- if (pooled) assignments else assignments[assignments$sample_name == sample_table$sample_name[sample_index], , drop = FALSE]
      if (!nrow(reads)) next
      values <- if (pooled) list() else as.list(sample_table[sample_index, , drop = FALSE])
      for (field in names(region)) values[[field]] <- region[[field]][1L]
      if (!is.null(task$size)) values$size <- task$size
      path <- expand_template(spec$source, values)
      bed <- read_footprint_region(path, spec$format, region, reads$original_RID)
      if (!nrow(bed)) next
      if (spec$format %in% c("bed12_fibertools", "bed13_fiberhmm")) {
        blocks <- convert_ft_bed12_to_bed6(bed, format = spec$format,
          longest_alignment = (spec$format == "bed12_fibertools"),
          validate_blocks = (spec$format == "bed13_fiberhmm"), source = path)
      } else {
        blocks <- as.data.frame(bed[, 1:4, drop = FALSE], stringsAsFactors = FALSE)
        names(blocks) <- c("chr", "start", "end", "RID")
      }
      if (!nrow(blocks)) next
      blocks$size <- blocks$end - blocks$start
      selected <- blocks$start < region$end & blocks$end > region$start
      blocks <- blocks[selected, , drop = FALSE]
      if (!nrow(blocks)) next
      matched <- match(as.character(blocks$RID), as.character(reads$original_RID))
      stopifnot(!anyNA(matched))
      blocks$original_RID <- as.character(blocks$RID)
      blocks$RID <- as.character(reads$RID[matched])
      blocks$sample_name <- as.character(reads$sample_name[matched])
      blocks$track <- task$track
      pieces[[length(pieces) + 1L]] <- blocks[, names(empty_footprints()), drop = FALSE]
    }
  }
  records <- if (length(pieces)) dplyr::bind_rows(pieces) else empty_footprints()
  c(list(records = records), track_spec)
}
