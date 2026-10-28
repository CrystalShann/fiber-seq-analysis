#!/usr/bin/env Rscript
# Rebuild consensus from existing inputs; source URLs are in the annotation README.
suppressPackageStartupMessages(library(data.table))

validate_gzip <- function(path) {
    if (file.info(path)$size <= 0)
        stop("Empty input: ", path)
    con <- gzfile(path, "rb")
    on.exit(close(con))
    withCallingHandlers(while (length(readBin(con, "raw", 1048576L))) {
    }, warning = function(w) stop(conditionMessage(w)))
}

read_gzip_table <- function(path, ...) {
    con <- gzfile(path, "rt")
    on.exit(close(con))
    fread(text = paste(readLines(con, warn = FALSE), collapse = "\n"), sep = "\t", quote = "", na.strings = "",
        ...)
}

read_replicate_counts <- function(path, chunk_size = 200000L) {
    con <- gzfile(path, "rt")
    on.exit(close(con))
    parts <- list()
    withCallingHandlers(repeat {
        lines <- readLines(con, n = chunk_size, warn = FALSE)
        if (!length(lines))
            break
        tab <- fread(text = paste0(paste(lines, collapse = "\n"), "\n"), header = FALSE, sep = "\t", quote = "",
            select = c(1L, 2L, 3L, 5L, 6L))
        setnames(tab, c("chrom", "start0", "end0", "count", "strand"))
        if (!is.numeric(tab$count) || anyNA(tab) || any(!is.finite(tab$count)) || any(tab$count <= 0 | tab$count%%1 !=
            0) || any(tab$end0 != tab$start0 + 1) || any(tab$start0 < 0 | tab$start0%%1 != 0) || any(!tab$strand %in%
            c("+", "-")))
            stop("Invalid single-base CTSS record: ", path)
        parts[[length(parts) + 1L]] <- tab
    }, warning = function(w) stop(conditionMessage(w)))
    x <- rbindlist(parts)
    if (!nrow(x))
        stop("Empty CTSS file: ", path)
    if (anyDuplicated(x[, .(chrom, start0, strand)]))
        stop("Duplicate CTSS within replicate: ", path)
    x
}

consensus_bed <- function(common, minimum_count = 1) {
    stopifnot(length(minimum_count) == 1L, is.finite(minimum_count), minimum_count >= 1, minimum_count%%1 == 0)
    common[minimum_replicate_count >= minimum_count, .(chrom, start0, end0, name = paste0(chrom, ":", start0, "..",
        end0, ",", strand), score = summed_count, strand, thickStart = start0, thickEnd = end0, itemRgb = "0,0,0")]
}

intersect_replicate_counts <- function(paths) {
  stopifnot(length(paths) >= 2L, !anyDuplicated(paths))
  common <- NULL
  for (i in seq_along(paths)) {
    message("Reading replicate ", i)
    x <- read_replicate_counts(paths[i])
    setnames(x, "count", paste0("rep", i, "_count"))
    common <- if (is.null(common)) x else merge(common, x, by = c("chrom", "start0", "end0", "strand"), sort = FALSE)
    rm(x); gc(verbose = FALSE)
  }
  count_cols <- paste0("rep", seq_along(paths), "_count")
  common[, summed_count := rowSums(.SD), .SDcols = count_cols]
  common[, minimum_replicate_count := do.call(pmin, .SD), .SDcols = count_cols]
  setorder(common, chrom, start0, strand)
  common
}

setDTthreads(2L)
cfg <- list(input_dir = "/project/spott/cshan/annotations/fantom5_LCL",
  output_file = "fantom5.hg38.LCL.consensus.TSS.bed.gz", minimum_count = 1)
for (arg in commandArgs(trailingOnly = TRUE)) {
  if (!grepl("^--[a-z_-]+=", arg)) stop("Use --option=value")
  name <- gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", arg))
  if (!name %in% names(cfg)) stop("Unknown option: ", name)
  value <- sub("^--[^=]+=", "", arg)
  cfg[[name]] <- if (is.numeric(cfg[[name]])) as.numeric(value) else value
}
stopifnot(basename(cfg$output_file) == cfg$output_file, endsWith(cfg$output_file, ".gz"))
paths <- file.path(cfg$input_dir, paste0("fantom5.hg38.LCL.rep", 1:3, ".CTSS.bed.gz"))
if (!all(file.exists(paths))) stop("Missing CTSS inputs; see the annotation README for source URLs.")
bed <- consensus_bed(intersect_replicate_counts(paths), cfg$minimum_count)
if (!nrow(bed)) stop("No positions meet the consensus criteria")
output <- file.path(cfg$input_dir, cfg$output_file)
reuse <- FALSE
if (file.exists(output)) {
  old <- read_gzip_table(output, header = FALSE)
  if (ncol(old) == ncol(bed)) {
    setnames(old, names(bed)); setorder(old, chrom, start0, strand)
    reuse <- isTRUE(all.equal(as.data.frame(old), as.data.frame(bed), check.attributes = FALSE))
  }
}
if (reuse) {
  message("Reusing identical consensus: ", output)
} else {
  part <- paste0(output, ".part.gz")
  con <- gzfile(part, "wt")
  tryCatch(writeLines(do.call(paste, c(as.list(bed), sep = "\t")), con), finally = close(con))
  validate_gzip(part)
  stopifnot(file.rename(part, output))
}
message("Consensus positions: ", nrow(bed))
