#!/usr/bin/env Rscript
# Regression for preserving Analysis 3's batch-dependent alignment selection.
# Run from the repository root. All BED queries and outputs stay in memory.

read_notebook_chunk <- function(lines, label) {
  first <- which(startsWith(lines, paste0("```{r ", label, ",")) |
                 startsWith(lines, paste0("```{r ", label, "}")))
  stopifnot(length(first) == 1L)
  last <- first + which(lines[(first + 1L):length(lines)] == "```")[1L]
  stopifnot(!is.na(last))
  parse(text = lines[(first + 1L):(last - 1L)])
}

project_root <- Sys.getenv("ENHANCER_PROJECT_ROOT", getwd())
notebook <- readLines(file.path(project_root, "code/enhancer/enhancer_co-occ.Rmd"))
fixture <- new.env(parent = globalenv())
sys.source(file.path(project_root, "code/enhancer/co_occ_functions/run_enhancer_footprints.R"), envir = fixture)
eval(read_notebook_chunk(notebook, "site-pair-input-functions"), envir = fixture)
eval(read_notebook_chunk(notebook, "site-pair-statistic-functions"), envir = fixture)
calculation <- read_notebook_chunk(notebook, "calculate-site-pair-cooccupancy")

eval(quote({
  # Original batch 1: active target + 98 empty loci + an inactive competitor.
  # Original batch 2: active target + a low-coverage active competitor.
  # Filtering either competitor out before batching changes the winning
  # alignment for its split read, even though all input files are identical.
  target_one <- data.table::data.table(enhancer_id = "active_target_1", chr = "chr1",
    start0 = 100L, end0 = 300L, enhancer_class = "active")
  fillers <- data.table::data.table(enhancer_id = paste0("filler_", seq_len(98L)),
    chr = "chr1", start0 = 10000L + seq_len(98L) * 1000L,
    end0 = 10100L + seq_len(98L) * 1000L, enhancer_class = "inactive")
  inactive_competitor <- data.table::data.table(enhancer_id = "inactive_competitor", chr = "chr1",
    start0 = 1100L, end0 = 1300L, enhancer_class = "inactive")
  target_two <- data.table::data.table(enhancer_id = "active_target_2", chr = "chr1",
    start0 = 2100L, end0 = 2300L, enhancer_class = "inactive")
  active_competitor <- data.table::data.table(enhancer_id = "active_low_coverage", chr = "chr1",
    start0 = 3100L, end0 = 3300L, enhancer_class = "active")
  original_enhancers <- data.table::rbindlist(list(target_one, fillers,
    inactive_competitor, target_two, active_competitor))
  stopifnot(nrow(original_enhancers) == 102L)
  sample_table <- data.table::data.table(sample = paste0("LPS_", c(0L, 5L, 10L, 15L)),
                                       time_min = c(0L, 5L, 10L, 15L))
  raw_spans <- data.table::rbindlist(lapply(0:1, function(index) {
    shift <- 2000L * index
    data.table::data.table(chr = "chr1", start0 = c(rep(90L + shift, 11L), 1000L + shift),
      end0 = c(rep(310L + shift, 11L), 1500L + shift),
      read_id = c(paste0("target", index + 1L, "_stable", seq_len(10L)),
                  rep(paste0("target", index + 1L, "_split"), 2L)))
  }))
  raw_tf <- data.table::rbindlist(lapply(0:1, function(index) {
    shift <- 2000L * index
    data.table::rbindlist(lapply(c(seq_len(8L), 11L), function(read_index) {
      # Among the 10 stable fibers: 4 both-bound, 2 site-1-only,
      # 2 site-2-only and 2 with no TF call. Split fibers have both calls.
      offsets <- if (read_index <= 4L || read_index == 11L) c(55L, 145L) else
        if (read_index <= 6L) 55L else 145L
      rid <- if (read_index == 11L) paste0("target", index + 1L, "_split") else
        paste0("target", index + 1L, "_stable", read_index)
      data.table::data.table(V1 = "chr1", V2 = 90L + shift, V3 = 310L + shift,
        V4 = rid, V5 = 0L, V6 = "+", V7 = 90L + shift, V8 = 310L + shift,
        V9 = "0", V10 = length(offsets), V11 = paste(rep(10L, length(offsets)), collapse = ","),
        V12 = paste(offsets, collapse = ","))
    }))
  }))

  extract_root <- "/in-memory-fixture"
  ft_result_dir <- "/in-memory-fixture/m6a"
  hmm_file <- function(sample, chr, kind)
    paste0("/in-memory-fixture/firehmm_", kind, "/", sample, "/", chr)
  tabix_bed <- function(path, regions) {
    raw <- if (grepl("m6a_by_chr", path, fixed = TRUE)) raw_spans else
      if (grepl("firehmm_tf", path, fixed = TRUE)) raw_tf else data.table::data.table()
    if (!nrow(raw)) return(raw)
    spans <- data.table::data.table(chr = raw[[1L]], start0 = raw[[2L]], end0 = raw[[3L]])
    keep <- unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(
      bed_ranges(spans), regions, ignore.strand = TRUE)))
    data.table::copy(raw[keep])
  }
  saved_tables <- list()
  save_table <- function(value, name) saved_tables[[name]] <<- data.table::copy(value)

  footprint_result <- count_enhancer_footprints(original_enhancers, sample_table,
                                               extract_root, ft_result_dir)
  footprint_result$inputs <- footprint_job_inputs(original_enhancers, sample_table,
    "chr1", extract_root, ft_result_dir, "genomewide")
  footprint_counts <- data.table::copy(footprint_result$footprint_counts)
  eligible_enhancers <- footprint_counts[, .(keep = .N == nrow(sample_table) &&
                                             all(n_reads >= 10L)), by = enhancer_id][keep == TRUE, enhancer_id]
  stopifnot(setequal(eligible_enhancers, c("active_target_1", "active_target_2")),
    all(footprint_counts[enhancer_id %in% eligible_enhancers, n_reads] == 10L))
  # An interactive reorder must not replace the completed job's batch order.
  analysis_enhancers <- data.table::copy(original_enhancers[rev(seq_len(.N))])
  active_targets <- analysis_enhancers[enhancer_id %in% eligible_enhancers]
  old_coverage <- load_site_pair_fibers(active_targets, sample_table)$coverage
  old_totals <- old_coverage[, .(loaded_reads = .N), by = .(enhancer_id, sample)]
  stopifnot(nrow(old_totals) == 8L, all(old_totals$loaded_reads == 11L))
  cat("Reproduced prior bug: active/eligible rebatching gives 11 fibers; saved Analysis 3 gives 10.\n")

  site_window_bp <- 50L
  min_site_observations <- 5L
  min_pair_reads <- 10L
  cobinding_alpha <- 0.01
  cobinding_significance <- "q_value"
}), envir = fixture)

eval(calculation, envir = fixture)
eval(quote({
  tested <- site_pair_cooccupancy[status == "tested"]
  stopifnot(nrow(tested) == 8L, all(tested$total == 10L),
    setequal(tested$enhancer_id, c("active_target_1", "active_target_2")),
    setequal(tested$enhancer_class,c("active","inactive")),
    all(tested$n_tf_tf == 4L), all(tested$n_tf_naked == 2L),
    all(tested$n_naked_tf == 2L), all(tested$n_naked_naked == 2L),
    all(is.finite(tested$p_value)), all(is.finite(tested$q_value)),
    identical(saved_tables[["05_site_pair_cooccupancy.tsv"]], site_pair_cooccupancy))
}), envir = fixture)

# A real reference/input discrepancy must still stop; preserving old batches
# must not turn the audit into an unconditional pass or overwrite the reference.
fixture$footprint_counts[enhancer_id == "active_target_1" & sample == "LPS_0", n_reads := n_reads + 1L]
failure <- tryCatch({eval(calculation, envir = fixture); NULL}, error = identity)
stopifnot(inherits(failure, "error"),
          grepl("spanning-read|Analysis 3|denominator", conditionMessage(failure), ignore.case = TRUE))

# A missing reference row must be caught by the expected enhancer/sample grid.
fixture$footprint_counts <- data.table::copy(fixture$footprint_result$footprint_counts)[
  !(enhancer_id == "active_target_1" & sample == "LPS_0")]
missing_reference <- tryCatch({eval(calculation, envir = fixture); NULL}, error = identity)
stopifnot(inherits(missing_reference, "error"),
  grepl("spanning-read|Analysis 3|denominator", conditionMessage(missing_reference), ignore.case = TRUE))

# No eligible active targets is a valid empty analysis, with a typed result.
fixture$footprint_counts <- data.table::copy(fixture$footprint_result$footprint_counts)
fixture$eligible_enhancers <- character()
eval(calculation, envir = fixture)
stopifnot(nrow(fixture$site_pair_cooccupancy) == 0L,
  all(names(fixture$empty_site_pair_results()) %in% names(fixture$site_pair_cooccupancy)),
  nrow(fixture$saved_tables[["05_site_pair_cooccupancy.tsv"]]) == 0L)
cat("PASS: original saved batches preserve pair-state totals; mismatched/missing references stop; empty targets succeed.\n")
