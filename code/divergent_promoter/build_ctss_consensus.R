#!/usr/bin/env Rscript
# Rebuild consensus from existing inputs; source URLs are in the annotation README.
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1L])
source(file.path(dirname(normalizePath(script)), "divergent_promoter_functions.R"))
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
