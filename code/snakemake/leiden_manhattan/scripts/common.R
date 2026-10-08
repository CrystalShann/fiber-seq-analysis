# Shared I/O only; clustering and plotting methods are sourced from their owners.
load_shared <- function(cfg, plots = FALSE, footprints = FALSE) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(Matrix)
    library(GenomicRanges)
    library(data.table)
    library(Rsamtools)
    library(igraph)
  })
  source(file.path(cfg$project_root, cfg$parsing_functions), local = .GlobalEnv)
  source(file.path(cfg$project_root, cfg$plotting_functions), local = .GlobalEnv)
  source(file.path(cfg$code_dir, "clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r"), local = .GlobalEnv)
  source(file.path(cfg$code_dir, "haplotype_phasing/LCL_phasing.r"), local = .GlobalEnv)
  if (plots || footprints) {
    source(file.path(cfg$code_dir, "clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r"), local = .GlobalEnv)
  }
  if (plots) {
    suppressPackageStartupMessages({library(ggplot2); library(cowplot); library(grid)})
    source(file.path(cfg$code_dir, "clustering_methods/Leiden_Manhattan/diff_avg_m6a.R"), local = .GlobalEnv)
  }
  data.table::setDTthreads(1L)
  options(scipen = 999)
  invisible(NULL)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

read_samples <- function(path) {
  x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sample_name", "sample_label", "fire_dir")
  if (!all(required %in% names(x))) stop("Sample sheet lacks: ", paste(setdiff(required, names(x)), collapse = ", "))
  if (!nrow(x) || anyNA(x[, required]) || any(!nzchar(x$sample_name)) ||
      anyDuplicated(x$sample_name) || anyDuplicated(x$sample_label))
    stop("Sample names and labels must be nonempty and unique: ", path)
  x
}

expand_template <- function(template, values) {
  fields <- regmatches(template, gregexpr("\\{[^{}]+\\}", template))[[1]]
  for (field in unique(fields)) {
    key <- substring(field, 2L, nchar(field) - 1L)
    value <- values[[key]]
    if (is.null(value) || length(value) != 1L || is.na(value)) stop("Missing template value: ", key)
    template <- gsub(field, as.character(value), template, fixed = TRUE)
  }
  template
}

read_region <- function(path, id) {
  regions <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  region <- regions[regions$region_id == id, , drop = FALSE]
  if (nrow(region) != 1L) stop("Expected one region: ", id)
  region
}

region_samples <- function(region, samples) {
  mode <- if ("samples_mode" %in% names(region)) region$samples_mode else "all"
  if (identical(mode, "contributing")) return(lcl_region_samples(region, samples))
  if (!identical(mode, "all")) stop("Unknown region samples mode: ", mode)
  samples
}

save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path)
}
