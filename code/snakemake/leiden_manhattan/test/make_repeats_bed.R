#!/usr/bin/env Rscript
# Convert repeats_subset (UCSC RepeatMasker rmsk table subset) into a BED file
# usable as a workflow regions.bed (see ../README.md "Custom regions"):
# 0-based, half-open; column 4 = unique region_id.
#
#   Rscript make_repeats_bed.R [repeats_subset.csv] [output.bed]
#
# genoStart/genoEnd in the rmsk table are already 0-based, half-open, same as BED.

args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[[1]] else "repeats_subset"
output <- if (length(args) >= 2) args[[2]] else "repeats_regions.bed"

x <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("genoName", "genoStart", "genoEnd", "repName", "strand")
missing <- setdiff(required, names(x))
if (length(missing)) stop("repeats table lacks columns: ", paste(missing, collapse = ", "))

rep_name <- gsub("[^A-Za-z0-9_.-]", "_", x$repName)
region_id <- paste(rep_name, x$genoName, x$genoStart + 1L, x$genoEnd, sep = "_")
if (anyDuplicated(region_id)) stop("Generated region IDs are not unique")

bed <- data.frame(chrom = x$genoName, start = x$genoStart, end = x$genoEnd,
  region_id = region_id, repName = x$repName, strand = x$strand,
  stringsAsFactors = FALSE)
bed <- bed[order(bed$chrom, bed$start), ]

write.table(bed, output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
message("Wrote ", nrow(bed), " regions to ", normalizePath(output))
