#!/usr/bin/env Rscript
# Shared helpers and batch worker; sourcing this file does not run the analysis.
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(Rsamtools)
})

bed_ranges <- function(x) {
  stopifnot(all(x$start0 >= 0), all(x$end0 > x$start0))
  GenomicRanges::GRanges(x$chr, IRanges::IRanges(x$start0 + 1L, x$end0))
}

tabix_bed <- function(path, regions) {
  if (!file.exists(path) || !file.exists(paste0(path, ".tbi")))
    stop("Missing BED or tabix index: ", path)
  if (!length(regions)) return(data.table::data.table())
  tf <- Rsamtools::TabixFile(path)
  open(tf)
  on.exit(close(tf))
  # A long fiber may overlap several query regions: remove duplicate lines.
  lines <- unique(unlist(Rsamtools::scanTabix(tf, param = GenomicRanges::reduce(regions)),
                         use.names = FALSE))
  if (!length(lines)) return(data.table::data.table())
  data.table::fread(text = paste(lines, collapse = "\n"), header = FALSE,
                    sep = "\t", showProgress = FALSE)
}

expand_hmm_blocks <- function(raw) {
  empty <- data.table(chr = character(), start0 = integer(), end0 = integer(),
                      read_id = character(), strand = character(), length = integer())
  if (!nrow(raw)) return(empty)
  if (ncol(raw) < 12L) stop("Expected FiberHMM BED12 or BED12+ input")
  sizes <- strsplit(sub(",$", "", as.character(raw[[11L]])), ",", fixed = TRUE)
  offsets <- strsplit(sub(",$", "", as.character(raw[[12L]])), ",", fixed = TRUE)
  counts <- as.integer(raw[[10L]])
  if (anyNA(counts) || any(lengths(sizes) != counts) || any(lengths(offsets) != counts))
    stop("Invalid FiberHMM block counts")
  parent <- rep(seq_len(nrow(raw)), counts)
  size <- as.integer(unlist(sizes, use.names = FALSE))
  offset <- as.integer(unlist(offsets, use.names = FALSE))
  start0 <- as.integer(raw[[2L]][parent]) + offset
  if (anyNA(size) || anyNA(offset) || any(size <= 0L) || any(offset < 0L) ||
      any(start0 + size > as.integer(raw[[3L]][parent]))) stop("Invalid FiberHMM blocks")
  unique(data.table(chr = as.character(raw[[1L]][parent]), start0 = start0,
                    end0 = start0 + size, read_id = as.character(raw[[4L]][parent]),
                    strand = as.character(raw[[6L]][parent]), length = size))
}

footprint_length_bin <- function(length_bp) {
  # Lower bin edges: 0 = 1–9 bp, 10 = 10–19 bp, ..., 500 = >=500 bp.
  pmin(500L, as.integer(length_bp %/% 10L) * 10L)
}

stream_footprint_lengths <- function(path, footprint_type, chunk_lines = 5000L) {
  stopifnot(footprint_type %in% c("TF", "nucleosome"), chunk_lines > 0L)
  if (!file.exists(path)) stop("Missing footprint input: ", path)
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  totals <- numeric(51L)
  repeat {
    lines <- readLines(con, n = chunk_lines, warn = FALSE)
    if (!length(lines)) break
    raw <- data.table::fread(text = paste(lines, collapse = "\n"),
                            header = FALSE, sep = "\t", select = 11L,
                            colClasses = "character", showProgress = FALSE)
    sizes <- as.integer(unlist(strsplit(raw[[1L]], ",", fixed = TRUE),
                               use.names = FALSE))
    if (anyNA(sizes)) stop("Invalid footprint block size in ", path)
    sizes <- if (footprint_type == "TF") sizes[sizes > 0L & sizes < 60L] else
      sizes[sizes > 90L]
    if (length(sizes)) {
      bins <- footprint_length_bin(sizes)
      totals <- totals + tabulate(bins %/% 10L + 1L, nbins = 51L)
    }
  }
  data.table::data.table(length_bin = seq.int(0L, 500L, 10L), n = totals)
}

expected_footprint_lengths <- function(sample_table, chromosomes, extract_root) {
  result <- vector("list", nrow(sample_table))
  for (s in seq_len(nrow(sample_table))) {
    sample_name <- sample_table$sample[s]
    compact_name <- gsub("_", "", sample_name, fixed = TRUE)
    pieces <- list()
    for (chromosome in chromosomes) {
      for (kind in c("tf", "footprint")) {
        path <- file.path(extract_root, paste0("firehmm_", kind), compact_name,
                          paste0(compact_name, "_hmm_extracted_", kind, "_",
                                 chromosome, ".bed.gz"))
        pieces[[length(pieces) + 1L]] <- stream_footprint_lengths(
          path, if (kind == "tf") "TF" else "nucleosome")
      }
    }
    combined <- data.table::rbindlist(pieces)[, .(n = sum(n)), by = length_bin]
    combined[, `:=`(sample = sample_name, time_min = sample_table$time_min[s])]
    result[[s]] <- combined[, .(sample, time_min, length_bin, n)]
  }
  data.table::rbindlist(result)
}

count_enhancer_footprints <- function(analysis_enhancers, sample_table,
                                      extract_root, ft_result_dir,
                                      batch_size = 100L) {
  stopifnot(!anyDuplicated(analysis_enhancers$enhancer_id), batch_size > 0L)
  enhancer_parts <- list()
  length_parts <- list()
  for (s in seq_len(nrow(sample_table))) {
    sample_name <- sample_table$sample[s]
    compact_name <- gsub("_", "", sample_name, fixed = TRUE)
    for (chromosome in unique(analysis_enhancers$chr)) {
      enh <- data.table::copy(analysis_enhancers[chr == chromosome])
      m6a_path <- file.path(ft_result_dir, sample_name, "extracted_results",
                           "m6a_by_chr", paste0(sample_name,
                             ".ft_extracted_m6a.", chromosome, ".bed.gz"))
      tf_path <- file.path(extract_root, "firehmm_tf", compact_name,
                          paste0(compact_name, "_hmm_extracted_tf_",
                                 chromosome, ".bed.gz"))
      nuc_path <- file.path(extract_root, "firehmm_footprint", compact_name,
                           paste0(compact_name, "_hmm_extracted_footprint_",
                                  chromosome, ".bed.gz"))
      for (first in seq.int(1L, nrow(enh), by = batch_size)) {
        last <- min(nrow(enh), first + batch_size - 1L)
        current <- enh[first:last]
        enh_gr <- bed_ranges(current)
        query_gr <- GenomicRanges::reduce(enh_gr)
        raw_spans <- tabix_bed(m6a_path, query_gr)
        coverage <- data.table::data.table(enhancer_id = character(),
                                          read_id = character())
        if (nrow(raw_spans)) {
          spans <- unique(data.table::data.table(
            chr = as.character(raw_spans[[1L]]), start0 = as.integer(raw_spans[[2L]]),
            end0 = as.integer(raw_spans[[3L]]), read_id = as.character(raw_spans[[4L]])))
          spans[, alignment_width := end0 - start0]
          data.table::setorderv(spans, c("read_id", "alignment_width"), c(1L, -1L))
          spans <- unique(spans, by = "read_id")
          hits <- GenomicRanges::findOverlaps(enh_gr, bed_ranges(spans),
                                               type = "within", ignore.strand = TRUE)
          coverage <- unique(data.table::data.table(
            enhancer_id = current$enhancer_id[S4Vectors::queryHits(hits)],
            read_id = spans$read_id[S4Vectors::subjectHits(hits)]))
        }
        observed <- list()
        for (kind in c("TF", "nucleosome")) {
          raw <- tabix_bed(if (kind == "TF") tf_path else nuc_path, query_gr)
          fp <- expand_hmm_blocks(raw)
          fp <- if (kind == "TF") fp[length > 0L & length < 60L] else fp[length > 90L]
          if (!nrow(fp) || !nrow(coverage)) next
          fp <- unique(fp, by = c("chr", "read_id", "start0", "end0"))
          # Convert BED midpoint to a 1-based point; length remains unclipped.
          midpoint <- floor((fp$start0 + fp$end0) / 2) + 1L
          mid_gr <- GenomicRanges::GRanges(fp$chr, IRanges::IRanges(midpoint, midpoint))
          hits <- GenomicRanges::findOverlaps(mid_gr, enh_gr, ignore.strand = TRUE)
          if (!length(hits)) next
          assigned <- data.table::data.table(
            enhancer_id = current$enhancer_id[S4Vectors::subjectHits(hits)],
            read_id = fp$read_id[S4Vectors::queryHits(hits)],
            length = fp$length[S4Vectors::queryHits(hits)])
          assigned <- merge(assigned, coverage, by = c("enhancer_id", "read_id"),
                            all = FALSE, sort = FALSE)
          if (!nrow(assigned)) next
          assigned[, footprint_type := kind]
          observed[[kind]] <- assigned
        }
        summary <- current[, .(enhancer_id, enhancer_class)]
        summary[, `:=`(n_reads = 0L, n_cooccupied = 0L, n_tf = 0L, n_nuc = 0L)]
        if (nrow(coverage)) {
          covered_counts <- coverage[, .(n_reads = .N), by = enhancer_id]
          summary[covered_counts, n_reads := i.n_reads, on = "enhancer_id"]
        }
        if (length(observed)) {
          calls <- data.table::rbindlist(observed)
          tf_counts <- calls[footprint_type == "TF", .N, by = .(enhancer_id, read_id)]
          if (nrow(tf_counts)) {
            tf_summary <- tf_counts[, .(n_tf = sum(N), n_cooccupied = sum(N >= 2L)),
                                     by = enhancer_id]
            summary[tf_summary, `:=`(n_tf = i.n_tf, n_cooccupied = i.n_cooccupied),
                    on = "enhancer_id"]
          }
          nuc_counts <- calls[footprint_type == "nucleosome", .(n_nuc = .N),
                              by = enhancer_id]
          if (nrow(nuc_counts)) summary[nuc_counts, n_nuc := i.n_nuc, on = "enhancer_id"]
          calls <- merge(calls, current[, .(enhancer_id, enhancer_class)],
                         by = "enhancer_id", all.x = TRUE, sort = FALSE)
          calls[, length_bin := footprint_length_bin(length)]
          lengths <- calls[, .(n = .N), by = .(enhancer_class, length_bin)]
          lengths[, `:=`(sample = sample_name, time_min = sample_table$time_min[s])]
          length_parts[[length(length_parts) + 1L]] <- lengths
        }
        summary[, `:=`(sample = sample_name, time_min = sample_table$time_min[s])]
        stopifnot(all(summary$n_cooccupied <= summary$n_reads),
                  all(summary$n_tf >= 2L * summary$n_cooccupied))
        enhancer_parts[[length(enhancer_parts) + 1L]] <- summary[
          , .(enhancer_id, sample, time_min, enhancer_class,
              n_reads, n_cooccupied, n_tf, n_nuc)]
      }
    }
  }
  counts <- data.table::rbindlist(enhancer_parts)
  if (length(length_parts)) {
    lengths <- data.table::rbindlist(length_parts)[
      , .(n = sum(n)), by = .(sample, time_min, enhancer_class, length_bin)]
  } else {
    lengths <- data.table::data.table(sample = character(), time_min = numeric(),
                                     enhancer_class = character(),
                                     length_bin = integer(), n = numeric())
  }
  list(footprint_counts = counts, length_counts = lengths)
}

footprint_job_inputs <- function(enhancers, samples, chromosomes, extract_root,
                                 ft_result_dir, background) {
  stopifnot(background %in% c("genomewide", "inactive"), nrow(enhancers) > 0L,
            !anyDuplicated(enhancers$enhancer_id), nrow(samples) > 0L,
            !anyDuplicated(samples$sample), length(chromosomes) > 0L)
  list(enhancers = as.data.frame(enhancers[, .(enhancer_id, chr, start0, end0, enhancer_class)]),
       samples = as.data.frame(samples[, .(sample, time_min)]),
       chromosomes = chromosomes, extract_root = extract_root,
       ft_result_dir = ft_result_dir, expected_background = background)
}

run_enhancer_footprints <- function(config_file) {
  config <- readRDS(config_file)
  inputs <- config$inputs
  enhancers <- as.data.table(inputs$enhancers)
  samples <- as.data.table(inputs$samples)
  stopifnot(inputs$expected_background %in% c("genomewide", "inactive"),
            dir.exists(dirname(config$result_file)), dir.exists(dirname(config$counts_file)))
  data.table::setDTthreads(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
  message("Counting footprints for ", nrow(enhancers), " enhancers and ", nrow(samples),
          " samples at ", Sys.time())
  result <- count_enhancer_footprints(enhancers, samples,
                                     inputs$extract_root, inputs$ft_result_dir)
  message("Calculating ", inputs$expected_background, " background at ", Sys.time())
  if (inputs$expected_background == "genomewide") {
    background_counts <- expected_footprint_lengths(samples, inputs$chromosomes,
                                                     inputs$extract_root)
  } else {
    background_counts <- merge(
      CJ(sample = samples$sample, length_bin = seq.int(0L, 500L, 10L)),
      result$length_counts[enhancer_class == "inactive", .(sample, length_bin, n)],
      by = c("sample", "length_bin"), all.x = TRUE)
    background_counts[is.na(n), n := 0]
    background_counts[, time_min := samples$time_min[match(sample, samples$sample)]]
  }
  result$expected_length_counts <- background_counts
  result$inputs <- inputs
  result$job_id <- Sys.getenv("SLURM_JOB_ID", "interactive")
  result$completed_at <- Sys.time()
  # Publish the RDS last so the notebook never loads a partially written result.
  temporary_rds <- tempfile(".footprints_", tmpdir = dirname(config$result_file))
  temporary_tsv <- tempfile(".footprints_", tmpdir = dirname(config$counts_file))
  on.exit(unlink(c(temporary_rds, temporary_tsv)), add = TRUE)
  saveRDS(result, temporary_rds)
  fwrite(result$footprint_counts, temporary_tsv, sep = "\t", na = "NA")
  if (!file.rename(temporary_tsv, config$counts_file) ||
      !file.rename(temporary_rds, config$result_file)) stop("Could not publish footprint results")
  message("Saved ", config$result_file, " at ", Sys.time())
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Usage: Rscript run_enhancer_footprints.R INPUT_CONFIG.rds")
  run_enhancer_footprints(args[[1L]])
}
