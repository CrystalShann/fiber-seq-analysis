#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr); library(Matrix); library(GenomicRanges); library(Rsamtools)
})

##########################################
# load and process m6a regions for input to lcl_build_phased_m6a_input
##########################################

# Resize each window around the lead SNP for each FIRE-SNP and returns a region table 

# 1. load `/project/spott/cshan/fiber-seq/LCL_project/Leiden_manhattan/summary tables/selection_signature.rds`
# 2. calculate the analysis_start and analysis_end for each region based on the window_size

load_m6a_regions <- function(selection, region_table = NULL, window_size = 2000L) {
  regions <- selection$regions
  # The RDS is the region authority; an external table is an optional cross-check.
  if (!is.null(region_table)) {
    saved <- read.delim(region_table, check.names = FALSE,
      colClasses = c(ref = "character", alt = "character"))
    required <- c("region_id", "chr", "analysis_start", "analysis_end")
    if (!isTRUE(all.equal(regions[, required], saved[, required],
        check.attributes = FALSE, tolerance = 0)))
      stop("Leiden selection signature and region table disagree")
  }
  stopifnot(nrow(regions) > 0L, !anyDuplicated(regions$region_id),
    all(regions$analysis_start == regions$start + 1L),
    all(regions$analysis_end == regions$end),
    all(regions$analysis_end - regions$analysis_start + 1L == regions$width),
    length(window_size) == 1L, is.finite(window_size),
    window_size == as.integer(window_size), window_size >= 3L)
  regions$width <- as.integer(window_size)
  regions$analysis_start <- as.integer(regions$focal_pos - window_size %/% 2L)
  regions$analysis_end <- regions$analysis_start + window_size - 1L
  regions$start <- regions$analysis_start - 1L
  regions$end <- regions$analysis_end
  if (anyNA(regions$analysis_start) || any(regions$analysis_start < 1L))
    stop("Requested window extends before chromosome position 1 or has missing coordinates")
  regions
}

build_m6a_input <- function(region, sample_table, ft_result_dir, project) {
  source(file.path(project, "code/topic_model/topic_modelling_functions.r"), local = TRUE)
  source(file.path(project, "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r"), local = TRUE)
  source(file.path(project, "code/haplotype_phasing/LCL_phasing.r"), local = TRUE)
  source(file.path(project, "code/haplotype_phasing/LCL_phased_m6a_input.r"), local = TRUE)
  lcl_build_phased_m6a_input(region, sample_table, ft_result_dir,
    "/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples",
    file.path(project, "LCL_project/Leiden_manhattan/phase_summary"))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %in% c(3L, 4L) && args[1] == "--regions") {
    region_table <- if (length(args) == 4L) args[3] else NULL
    regions <- load_m6a_regions(readRDS(args[2]), region_table, as.numeric(tail(args, 1L)))
    write.table(regions, stdout(), sep = "\t", quote = TRUE, row.names = FALSE)
    return(invisible(NULL))
  }
  if (!length(args) %in% c(4L, 5L))
    stop("Usage: 02_build_m6a_input.r functions.r selection_signature.rds region_id ft_result_dir [window_size=2000]\n",
      "   or: 02_build_m6a_input.r --regions selection_signature.rds [selected_regions.tsv] window_size")
  project <- dirname(dirname(dirname(normalizePath(args[1]))))
  selection <- readRDS(args[2])
  window_size <- if (length(args) >= 5L) as.numeric(args[5]) else 2000L
  regions <- load_m6a_regions(selection, window_size = window_size)
  region <- regions[regions$region_id == args[3], , drop = FALSE]
  stopifnot(nrow(region) == 1L)
  message(R.version.string)
  z <- build_m6a_input(region, selection$sample_table, args[4], project)
  calls <- summary(z$met_mat)
  offsets <- split(calls$j[calls$x == 1] - 1L, calls$i[calls$x == 1])
  z$read_meta$call_offsets <- vapply(seq_len(nrow(z$met_mat)), function(i)
    paste(offsets[[as.character(i)]], collapse = ","), "")
  # Pipe sparse call offsets to Python; never write a matrix or input cache.
  write.table(z$read_meta, stdout(), sep = "\t", quote = TRUE, row.names = FALSE)
}
if (sys.nframe() == 0L) main()
