# Shared footprint extraction for the fixed Manhattan/ACF cohort.
collect_enhancer_footprints <- function(assignments, regions, extract_root) {
  a <- data.table::as.data.table(assignments)
  regions <- data.table::as.data.table(regions)
  empty <- data.table::data.table(RID = character(), track = character(),
    start = integer(), end = integer(), original_size = integer())
  out <- list(); part <- 0L
  groups <- unique(a[, .(enhancer_id, sample_name)])
  for (i in seq_len(nrow(groups))) {
    if (i == 1L || i %% 250L == 0L || i == nrow(groups))
      message("Footprint group ", i, "/", nrow(groups))
    enhancer <- groups$enhancer_id[i]; sample <- groups$sample_name[i]
    selected <- a[enhancer_id == enhancer & sample_name == sample]
    region <- regions[match(enhancer, enhancer_id)]
    stopifnot(nrow(region) == 1L, !is.na(region$midpoint0))
    left0 <- region$midpoint0 - 500L; right0 <- region$midpoint0 + 500L
    query <- GenomicRanges::GRanges(region$chr, IRanges::IRanges(left0 + 1L, right0))
    compact <- gsub("_", "", sample, fixed = TRUE)
    for (kind in c("footprint", "tf")) {
      path <- file.path(extract_root, paste0("firehmm_", kind), compact,
        paste0(compact, "_hmm_extracted_", kind, "_", region$chr, ".bed.gz"))
      if (!all(file.exists(c(path, paste0(path, ".tbi"))))) stop("Missing footprint BED/index: ", path)
      lines <- unique(unlist(Rsamtools::scanTabix(path, param = query), use.names = FALSE))
      if (!length(lines)) next
      bed <- data.table::fread(text = paste(lines, collapse = "\n"), header = FALSE,
                              sep = "\t", showProgress = FALSE)
      if (ncol(bed) < 12L) stop("Expected FiberHMM BED12: ", path)
      bed <- bed[as.character(bed[[4L]]) %in% selected$read_id]
      if (!nrow(bed)) next
      counts <- as.integer(bed[[10L]])
      stopifnot(!anyNA(counts), all(counts >= 0L))
      bed <- bed[counts > 0L]
      if (!nrow(bed)) next
      # FiberHMM blocks are all real calls: unlike ft m6A, do not drop sentinels.
      sizes <- lapply(strsplit(sub(",$", "", as.character(bed[[11L]])), ",", fixed = TRUE), as.integer)
      offsets <- lapply(strsplit(sub(",$", "", as.character(bed[[12L]])), ",", fixed = TRUE), as.integer)
      counts <- as.integer(bed[[10L]])
      stopifnot(!anyNA(counts), all(lengths(sizes) == counts), all(lengths(offsets) == counts))
      parent <- rep(seq_len(nrow(bed)), counts)
      size <- unlist(sizes, use.names = FALSE); offset <- unlist(offsets, use.names = FALSE)
      start0 <- as.integer(bed[[2L]][parent]) + offset
      stopifnot(!anyNA(size), !anyNA(offset), all(size > 0L), all(offset >= 0L),
                all(start0 + size <= as.integer(bed[[3L]][parent])))
      keep <- start0 < right0 & start0 + size > left0 &
        if (kind == "tf") size < 60L else size > 90L
      if (!any(keep)) next
      part <- part + 1L
      out[[part]] <- data.table::data.table(
        RID = selected$RID[match(as.character(bed[[4L]][parent[keep]]), selected$read_id)],
        track = if (kind == "tf") "TF <60 bp" else "Nucleosome >90 bp",
        start = as.integer(pmax(start0[keep], left0) - region$midpoint0),
        end = as.integer(pmin(start0[keep] + size[keep], right0) - region$midpoint0),
        original_size = as.integer(size[keep]))
    }
  }
  if (!length(out)) return(empty)
  records <- unique(data.table::rbindlist(out))
  stopifnot(!anyNA(records), all(records$RID %in% a$RID),
            all(records$start >= -500L), all(records$end <= 500L), all(records$end > records$start))
  records
}

# Match the macrophage reference's plotting contract, with relative coordinates.
