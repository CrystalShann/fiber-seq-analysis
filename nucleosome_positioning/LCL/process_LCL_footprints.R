# Process LCL nucleosome/TF footprints on fixed m6A-defined read clusters.
# Source leiden_LCL_functions.R first for BED12 and shared data utilities.
# Footprint display tracks only; the clustering input remains the m6A matrix.
LCL_FOOTPRINT_TRACKS <- c("ft_nuc_130-160bp", "FiberHMM_10-30bp",
                          "FiberHMM_40-60bp", "FiberHMM_60-80bp")

lcl_m6a_intervals <- function(result) {
  sites <- Matrix::summary(as(result$site_met_mat, "dgCMatrix"))
  sites <- sites[sites$x > 0, , drop = FALSE]
  positions <- as.integer(colnames(result$site_met_mat)[sites$j])
  identifiers <- rownames(result$site_met_mat)[sites$i]
  data.frame(RID = identifiers,
    original_RID = result$assignments$original_RID[match(identifiers, result$assignments$RID)],
    chr = rep(result$region$chr, length(identifiers)), start = positions - 1L,
    end = positions, size = rep(1L, length(identifiers)), track = rep("m6A", length(identifiers)))
}

lcl_plot_anchor <- function(region) {
  promoter <- region$region_type == "promoter" && !is.na(region$tss)
  anchor <- if (promoter) region$tss else floor(mean(c(region$analysis_start, region$analysis_end)))
  direction <- if (promoter && region$strand == "-") -1L else 1L
  bounds <- sort(direction * (c(region$analysis_start, region$analysis_end) - anchor))
  list(anchor = anchor, direction = direction, left = bounds[1], right = bounds[2],
    x_label = if (promoter) "Position relative to canonical TSS (bp)" else "Position relative to region centre (bp)")
}

prepare_lcl_read_tracks <- function(result, records, sample_colors) {
  anchor <- lcl_plot_anchor(result$region)
  reads <- result$assignments
  reads <- reads[order(reads$cluster, reads$start, reads$RID), , drop = FALSE]
  reads$row <- as.integer(ave(seq_len(nrow(reads)), reads$cluster, FUN = seq_along))
  reads$sample_label <- sub("_.*$", "", reads$sample_name)
  stopifnot(all(reads$sample_label %in% names(sample_colors)))
  if (!"haplotype" %in% names(reads)) reads$haplotype <- "pooled"
  relative_start <- anchor$direction * (reads$start - anchor$anchor)
  relative_end <- anchor$direction * (reads$end - anchor$anchor)
  reads$left <- pmax(pmin(relative_start, relative_end), anchor$left)
  reads$right <- pmin(pmax(relative_start, relative_end), anchor$right)
  matched <- match(records$RID, reads$RID)
  records <- records[!is.na(matched), , drop = FALSE]
  matched <- matched[!is.na(matched)]
  records$row <- reads$row[matched]
  records$cluster <- reads$cluster[matched]
  relative_start <- anchor$direction * (records$start + 1L - anchor$anchor)
  relative_end <- anchor$direction * (records$end - anchor$anchor)
  records$left <- pmax(pmin(relative_start, relative_end) - 0.5, anchor$left - 0.5)
  records$right <- pmin(pmax(relative_start, relative_end) + 0.5, anchor$right + 0.5)
  list(reads = reads, features = records, anchor = anchor, region = result$region)
}


extract_lcl_nucleosomes <- function(sample_table, region, assignments, min_size = 130L, max_size = 160L) {
  paths <- lcl_extracted_path(sample_table, region$chr, "nuc")
  result <- lapply(seq_len(nrow(sample_table)), function(sample_index) {
    bed <- read_lcl_bed12(paths[sample_index], region)
    if (!nrow(bed)) return(NULL)
    bed$RID <- paste(sample_table$sample_name[sample_index], bed$original_RID, sep = "::")
    bed <- bed[bed$RID %in% assignments$RID, , drop = FALSE]
    blocks <- lcl_bed12_blocks(bed)
    blocks[blocks$size >= min_size & blocks$size <= max_size &
             blocks$start < region$end & blocks$end > region$start, , drop = FALSE]
  })
  result <- dplyr::bind_rows(result)
  if (!ncol(result)) result <- lcl_bed12_blocks(data.frame())
  result$chr <- rep(region$chr, nrow(result))
  result$track <- rep(paste0("ft_nuc_", min_size, "-", max_size, "bp"), nrow(result))
  result
}

cache_lcl_footprint_tracks <- function(regions, tracks_dir, cache_dir, workers = 2L, reuse = TRUE) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- unlist(lapply(unique(regions$chr), function(chromosome) {
    files <- list.files(file.path(tracks_dir, chromosome), pattern = "_(10-30|40-60|60-80)bp_fps\\.bed\\.gz$", full.names = TRUE)
    if (!length(files)) stop("No footprint tracks for ", chromosome)
    files
  }), use.names = FALSE)
  jobs <- lapply(paths, function(path) {
    chromosome <- basename(dirname(path))
    list(path = path, regions = regions[regions$chr == chromosome, c("chr", "start", "end")])
  })
  extract_job <- function(job) {
    signature <- list(version = 1L, regions = job$regions, source = lcl_source_signature(job$path))
    cache <- file.path(cache_dir, paste0(basename(job$path), ".rds"))
    if (reuse && file.exists(cache) && identical(readRDS(cache)$signature, signature)) return(cache)
    conditions <- sprintf("($1 == \"%s\" && $2 < %.0f && $3 > %.0f)",
                            job$regions$chr, job$regions$end, job$regions$start)
    program <- paste0("BEGIN { FS = OFS = \"\\t\" } (", paste(conditions, collapse = " || "), ") { print }")
    temporary <- tempfile(tmpdir = cache_dir, fileext = ".bed")
    on.exit(unlink(temporary), add = TRUE)
    command <- paste("set -o pipefail; gzip -cd", shQuote(job$path), "| awk", shQuote(program), ">", shQuote(temporary))
    message("Extracting regional footprints: ", basename(job$path))
    status <- system2("/bin/bash", c("-c", shQuote(command)))
    if (status != 0L) stop("Footprint extraction failed: ", job$path)
    records <- data.frame(chr = character(), start = integer(), end = integer(), original_RID = character())
    if (file.info(temporary)$size > 0) {
      records <- data.table::fread(temporary, header = FALSE, data.table = FALSE)
      if (ncol(records) != 4L) stop("Expected BED4 footprints: ", job$path)
      names(records) <- c("chr", "start", "end", "original_RID")
    }
    records$size <- records$end - records$start
    records$track <- rep(sub("^combined_chr[^_]+_(.*)_fps\\.bed\\.gz$", "FiberHMM_\\1", basename(job$path)), nrow(records))
    saveRDS(list(records = records, signature = signature), cache)
    cache
  }
  caches <- parallel::mclapply(jobs, extract_job, mc.cores = workers, mc.preschedule = FALSE)
  if (any(vapply(caches, inherits, logical(1), "try-error"))) stop("One or more footprint extractions failed")
  unlist(caches, use.names = FALSE)
}

lcl_region_footprints <- function(caches, region, assignments) {
  if (anyDuplicated(assignments$original_RID)) {
    stop("Combined BED4 lacks sample IDs: ambiguous original read names in ", region$region_id)
  }
  caches <- caches[startsWith(basename(caches), paste0("combined_", region$chr, "_"))]
  pieces <- lapply(caches, function(cache) {
    records <- readRDS(cache)$records
    selected <- records$start < region$end & records$end > region$start
    records <- records[selected, , drop = FALSE]
    matched <- match(records$original_RID, assignments$original_RID)
    records$RID <- assignments$RID[matched]
    records[!is.na(matched), , drop = FALSE]
  })
  list(records = dplyr::bind_rows(pieces),
       tracks = vapply(caches, function(cache) {
         sub("^combined_chr[^_]+_(.*)_fps\\.bed\\.gz.rds$", "FiberHMM_\\1", basename(cache))
       }, character(1)))
}

lcl_occupancy_matrix <- function(records, assignments, region) {
  occupancy <- matrix(0L, nrow(assignments), region$width,
                       dimnames = list(assignments$RID, seq.int(region$analysis_start, region$analysis_end)))
  if (!nrow(records)) return(occupancy)
  for (record_index in seq_len(nrow(records))) {
    read_index <- match(records$RID[record_index], assignments$RID)
    if (is.na(read_index)) next
    left <- max(records$start[record_index], region$start) - region$start + 1L
    right <- min(records$end[record_index], region$end) - region$start
    if (left <= right) occupancy[read_index, seq.int(left, right)] <- 1L
  }
  occupancy
}

lcl_footprint_profiles <- function(records, tracks, assignments, region) {
  dplyr::bind_rows(lapply(tracks, function(track) {
    occupancy <- lcl_occupancy_matrix(records[records$track == track, ], assignments, region)
    dplyr::bind_rows(lapply(levels(assignments$cluster), function(cluster) {
      selected <- assignments$cluster == cluster
      data.frame(cluster = cluster, pos = as.integer(colnames(occupancy)), track = track,
                 fraction = colMeans(occupancy[selected, , drop = FALSE]), n_reads = sum(selected))
    }))
  }))
}


# Footprints use per-base occupancy; m6A profiles use the observed-site matrix.
# Do not insert zeros at unobserved genomic bases when drawing the m6A line.
lcl_signal_profiles <- function(result, footprints) {
  footprints <- footprints[footprints$track %in% LCL_FOOTPRINT_TRACKS, , drop = FALSE]
  records <- dplyr::bind_rows(lcl_m6a_intervals(result), footprints)
  fp <- lcl_footprint_profiles(footprints, LCL_FOOTPRINT_TRACKS, result$assignments, result$region)
  met <- cluster_site_profiles(result, result$site_met_mat)
  met <- data.frame(cluster = met$cluster, pos = met$pos, track = "m6A",
    fraction = met$met, n_reads = met$n_reads)
  list(records = records, tracks = c("m6A", LCL_FOOTPRINT_TRACKS),
       profiles = dplyr::bind_rows(met, fp))
}

process_lcl_footprints <- function(result_paths, sample_table, regions, tracks_dir,
                                   output_dir, workers = 2L, reuse = TRUE) {
  dir.create(file.path(output_dir, "summary tables"), recursive = TRUE, showWarnings = FALSE)
  caches <- cache_lcl_footprint_tracks(regions, tracks_dir, file.path(output_dir, "footprint_cache"),
    workers = workers, reuse = reuse)
  data.table::fwrite(dplyr::bind_rows(lapply(caches, function(path) readRDS(path)$signature$source)),
    file.path(output_dir, "summary tables", "footprint_source_manifest.tsv"), sep = "\t")
  for (path in result_paths) {
    result <- readRDS(path)
    message("Processing footprints: ", result$region$region_id)
    nuc <- extract_lcl_nucleosomes(sample_table, result$region, result$assignments)
    tf <- lcl_region_footprints(caches, result$region, result$assignments)$records
    records <- dplyr::bind_rows(nuc, tf)
    records <- records[records$track %in% LCL_FOOTPRINT_TRACKS, , drop = FALSE]
    matched <- match(records$RID, result$assignments$RID)
    records$cluster <- result$assignments$cluster[matched]
    records$sample_name <- result$assignments$sample_name[matched]
    records$region_id <- rep(result$region$region_id, nrow(records))
    data.table::fwrite(records, file.path(dirname(path), "footprint_positions.tsv.gz"), sep = "\t")
  }
  invisible(result_paths)
}
