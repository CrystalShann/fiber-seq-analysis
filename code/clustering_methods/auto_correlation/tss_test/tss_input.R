#!/usr/bin/env Rscript
# Read-only bridge: JSON travels through process pipes; no data files are saved.
suppressPackageStartupMessages({
  library(dplyr); library(Matrix); library(GenomicRanges); library(Rsamtools)
  library(data.table); library(jsonlite)
})
data.table::setDTthreads(1L)

expression_in_memory <- function(project) {
  # Source the exact expression calculation in an isolated environment. Mask
  # its two output functions locally, without changing the original script.
  e <- new.env(parent = globalenv())
  e$dir.create <- function(...) invisible(FALSE)
  e$fwrite <- function(x, ...) { e$captured <- data.table::copy(x); invisible(NULL) }
  e$message <- function(...) {
    txt <- paste0(..., collapse = "")
    if (!startsWith(txt, "Wrote ")) base::message(txt)
  }
  # fread(cmd=...) normally spools a temporary table. Read the command through
  # a pipe instead, so the GTF-derived exon table also stays in memory.
  e$fread <- function(..., cmd = NULL) {
    if (is.null(cmd)) return(data.table::fread(...))
    con <- pipe(cmd, open = "r")
    on.exit(close(con))
    data.table::fread(text = paste(readLines(con, warn = FALSE), collapse = "\n"), ...)
  }
  sink(stderr())
  on.exit(sink(), add = TRUE)
  sys.source(file.path(project, "code/accessibility/expr_access/01_expression_bins.R"), envir = e)
  stopifnot(!is.null(e$captured), identical(e$USE_TIMES, c(0L, 5L, 10L, 15L)), e$TPM_EXPR == 1)
  list(mapping = as.data.frame(e$captured), quartiles = unname(e$qs),
       rna_samples = as.character(e$sm$sample))
}

extract_tss <- function(region, samples, ft_root, max_reads = NULL, seed = 0L) {
  start <- as.integer(region$analysis_start)
  end <- as.integer(region$analysis_end)
  stopifnot(end - start + 1L == 2000L, start >= 1L)
  query <- GenomicRanges::GRanges(region$chrom, IRanges::IRanges(start, end))
  records <- audit <- list()
  for (sample in samples) {
    path <- file.path(ft_root, sample, "extracted_results/m6a_by_chr",
                      paste0(sample, ".ft_extracted_m6a.", region$chrom, ".bed.gz"))
    if (!all(file.exists(c(path, paste0(path, ".tbi")))))
      stop("Missing indexed m6A BED: ", path)
    # Identical call/sentinel/longest-alignment logic to 02_build_m6a_input.r's
    # lcl_build_phased_m6a_input, before the focal SNP genotype filter.
    reads <- extract_ft_region_reads(path, query)
    info <- extract_ft_read_info(reads)
    overlapping <- nrow(info)
    if (overlapping) info <- info[info$start <= start & info$end >= end, , drop = FALSE]
    full_span <- nrow(info)
    if (!is.null(max_reads) && full_span > max_reads) {
      set.seed(seed)
      info <- info[sort(sample.int(full_span, max_reads)), , drop = FALSE]
    }
    audit[[sample]] <- data.frame(region_id = region$region_id, sample_name = sample,
      overlapping_with_calls = overlapping, full_span = full_span, retained = nrow(info))
    if (!nrow(info)) next
    reads <- reads[reads$RID %in% info$RID, , drop = FALSE]
    mat <- get_sparse_met_mat(reads, info, start, end, seq.int(start, end), base = "A")
    stopifnot(ncol(mat) == 2000L, !anyNA(mat), all(mat@x %in% c(0, 1)),
              identical(rownames(mat), info$RID))
    calls <- summary(mat)
    offsets <- split(calls$j[calls$x == 1] - 1L, calls$i[calls$x == 1])
    info$call_offsets <- lapply(seq_len(nrow(info)), function(i)
      as.integer(offsets[[as.character(i)]]))
    info$sample_name <- sample
    info$region_id <- region$region_id
    info$row_id <- paste(region$region_id, sample, info$RID, sep = "::")
    names(info)[names(info) == "start"] <- "read_start"
    names(info)[names(info) == "end"] <- "read_end"
    records[[sample]] <- info
  }
  list(records = dplyr::bind_rows(records), audit = dplyr::bind_rows(audit))
}

main <- function() {
  project <- commandArgs(trailingOnly = TRUE)[1]
  source(file.path(project, "code/topic_model/topic_modelling_functions.r"), local = TRUE)
  # Resolve extraction functions in this source environment, without globals.
  environment(extract_tss) <- environment()
  input <- file("stdin", open = "r")
  on.exit(close(input))
  repeat {
    line <- readLines(input, n = 1L, warn = FALSE)
    if (!length(line)) break
    req <- jsonlite::fromJSON(line)
    answer <- switch(req$mode,
      expression = expression_in_memory(project),
      extract = extract_tss(req$region, req$samples, req$ft_root, req$max_reads, req$seed),
      stop("Unknown request mode"))
    cat(jsonlite::toJSON(answer, dataframe = "rows", auto_unbox = TRUE,
                         digits = NA, na = "null", null = "null"), "\n", sep = "")
    flush(stdout())
  }
}
if (sys.nframe() == 0L) main()
