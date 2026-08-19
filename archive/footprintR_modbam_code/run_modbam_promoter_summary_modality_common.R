source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/modbam_footprintR_functions.R")

default_sample <- Sys.getenv("SAMPLE", unset = "AL10_bc2178_19130")
default_chrom <- Sys.getenv("CHROM", unset = "chr1")
window_bp <- as.integer(Sys.getenv("WINDOW_BP", unset = "10000"))
bin_size <- as.integer(Sys.getenv("BIN_SIZE", unset = "100"))
mod_prob_threshold <- as.numeric(Sys.getenv("MOD_PROB_THRESHOLD", unset = "0.9"))
modbam_tss_chunk_size <- as.integer(Sys.getenv("MODBAM_TSS_CHUNK_SIZE", unset = "250"))

if (is.na(modbam_tss_chunk_size) || modbam_tss_chunk_size < 1L) {
  modbam_tss_chunk_size <- 250L
}

results_root_dir <- "/project/spott/cshan/fiber-seq/results/promoter_summary/10kb_TSS"
threshold_dir <- file.path(results_root_dir, "threshold_0.9")

bam_path <- file.path(
  "/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples",
  sprintf("%s.5mC.6mA.aligned.phased.bam", default_sample)
)
tss_path <- "/project/spott/cshan/annotations/TSS.gencode.v49.bed"

output_dir <- file.path(threshold_dir, default_sample)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bamfiles <- setNames(bam_path, default_sample)
bpparam_modbam <- BiocParallel::SerialParam()

tss_table_path <- file.path(output_dir, sprintf("%s_%s_tss_table.tsv.gz", default_sample, default_chrom))

assert_file_exists(bam_path, "BAM")
assert_file_exists(paste0(bam_path, ".bai"), "BAM index")
assert_file_exists(tss_path, "TSS BED")
assert_file_exists(tss_table_path, "part1 TSS table")

bam_seqinfo <- read_bam_seqinfo(bam_path)
tss_objects <- read_gencode_tss(tss_path, chrom = default_chrom, window_bp = window_bp, seqinfo = bam_seqinfo)
tss_table <- data.table::fread(tss_table_path)
promoter_windows <- tss_objects$promoter_windows
promoter_meta <- tss_objects$promoter_meta

append_tsv <- function(dt, path, wrote_once) {
  if (!nrow(dt)) {
    return(wrote_once)
  }

  data.table::fwrite(
    dt,
    file = path,
    sep = "\t",
    append = wrote_once,
    col.names = !wrote_once
  )
  TRUE
}

write_empty_if_missing <- function(path, empty_dt) {
  if (!file.exists(path)) {
    write_gz_tsv(empty_dt, path)
  }
}

run_modality_chunked <- function(modbase_code, summary_prefix, mod_label, meta_label) {
  paths <- list(
    positions = file.path(output_dir, sprintf("%s_%s_%s_position_calls.tsv.gz", default_sample, default_chrom, mod_label)),
    status = file.path(output_dir, sprintf("%s_%s_%s_modified_not_modified.tsv.gz", default_sample, default_chrom, mod_label)),
    aligned = file.path(output_dir, sprintf("%s_%s_%s_tss_aligned_positions.tsv.gz", default_sample, default_chrom, mod_label)),
    aligned_status = file.path(output_dir, sprintf("%s_%s_%s_tss_aligned_modified_not_modified.tsv.gz", default_sample, default_chrom, mod_label)),
    summary = file.path(output_dir, sprintf("%s_%s_%s.tsv.gz", default_sample, default_chrom, summary_prefix)),
    meta = file.path(output_dir, sprintf("%s_%s_%s_meta.tsv.gz", default_sample, default_chrom, mod_label))
  )

  all_paths <- unlist(paths, use.names = FALSE)

  if (all(file.exists(all_paths))) {
    message(mod_label, " outputs already exist; skipping.")
    return(invisible(paths))
  }

  existing <- all_paths[file.exists(all_paths)]
  if (length(existing)) {
    file.remove(existing)
  }

  n_tss <- nrow(tss_table)
  chunk_starts <- seq.int(1L, n_tss, by = modbam_tss_chunk_size)

  message("Processing ", mod_label, " in ", length(chunk_starts), " chunks (", modbam_tss_chunk_size, " TSS/chunk).")

  summary_chunks <- vector("list", length(chunk_starts))
  meta_acc <- data.table::data.table(bin = seq(-window_bp, window_bp, by = bin_size))
  meta_acc[, `:=`(modified_calls = 0, total_calls = 0)]

  wrote_positions <- FALSE
  wrote_status <- FALSE
  wrote_aligned <- FALSE
  wrote_aligned_status <- FALSE

  for (i in seq_along(chunk_starts)) {
    from <- chunk_starts[i]
    to <- min(from + modbam_tss_chunk_size - 1L, n_tss)

    message(sprintf("%s chunk %d/%d: TSS %d-%d", mod_label, i, length(chunk_starts), from, to))

    chunk_windows <- promoter_windows[from:to]
    chunk_meta <- promoter_meta[from:to]

    chunk_positions <- read_modbam_summary_positions(
      bamfiles = bamfiles,
      sample_name = default_sample,
      promoter_windows = chunk_windows,
      modbase_code = modbase_code,
      seqinfo = bam_seqinfo,
      bpparam = bpparam_modbam,
      modProbThreshold = mod_prob_threshold
    )

    chunk_status <- make_modified_status_table(chunk_positions)

    chunk_aligned <- align_positions_to_tss(
      position_dt = chunk_positions,
      promoter_windows = chunk_windows,
      promoter_meta = chunk_meta
    )

    chunk_aligned_status <- make_modified_status_table(chunk_aligned)

    summary_chunks[[i]] <- summarize_modality_by_tss(
      aligned_dt = chunk_aligned,
      prefix = summary_prefix
    )

    if (nrow(chunk_aligned)) {
      bin_chunk <- copy(chunk_aligned)
      bin_chunk[, bin := round(rel_pos / bin_size) * bin_size]
      bin_chunk <- bin_chunk[, .(
        modified_calls = sum(modified_calls),
        total_calls = sum(total_calls)
      ), by = bin]

      meta_acc[bin_chunk, on = "bin", `:=`(
        modified_calls = modified_calls + i.modified_calls,
        total_calls = total_calls + i.total_calls
      )]
    }

    wrote_positions <- append_tsv(chunk_positions, paths$positions, wrote_positions)
    wrote_status <- append_tsv(chunk_status, paths$status, wrote_status)
    wrote_aligned <- append_tsv(chunk_aligned, paths$aligned, wrote_aligned)
    wrote_aligned_status <- append_tsv(chunk_aligned_status, paths$aligned_status, wrote_aligned_status)

    rm(chunk_positions, chunk_status, chunk_aligned, chunk_aligned_status)
    gc()
  }

  summary_dt <- data.table::rbindlist(summary_chunks, use.names = TRUE, fill = TRUE)

  modified_col <- sprintf("%s_modified_calls", summary_prefix)
  unmodified_col <- sprintf("%s_unmodified_calls", summary_prefix)
  total_col <- sprintf("%s_total_calls", summary_prefix)
  coverage_col <- sprintf("%s_positions_with_coverage", summary_prefix)
  frac_col <- sprintf("%s_fraction_modified", summary_prefix)

  if (!nrow(summary_dt)) {
    summary_dt <- data.table::data.table(tss_uid = character())
    summary_dt[, (modified_col) := integer()]
    summary_dt[, (unmodified_col) := integer()]
    summary_dt[, (total_col) := integer()]
    summary_dt[, (coverage_col) := integer()]
    summary_dt[, (frac_col) := numeric()]
  } else {
    summary_dt <- summary_dt[, .(
      modified_calls = sum(get(modified_col)),
      unmodified_calls = sum(get(unmodified_col)),
      total_calls = sum(get(total_col)),
      positions_with_coverage = sum(get(coverage_col))
    ), by = tss_uid]

    summary_dt[, fraction_modified := data.table::fifelse(total_calls > 0, modified_calls / total_calls, NA_real_)]

    data.table::setnames(
      summary_dt,
      old = c("modified_calls", "unmodified_calls", "total_calls", "positions_with_coverage", "fraction_modified"),
      new = c(modified_col, unmodified_col, total_col, coverage_col, frac_col)
    )
  }

  meta_dt <- copy(meta_acc)
  meta_dt[, fraction_modified := data.table::fifelse(total_calls > 0, modified_calls / total_calls, NA_real_)]
  meta_dt[, modified_calls_per_tss := modified_calls / n_tss]
  meta_dt[, label := meta_label]

  write_empty_if_missing(
    paths$positions,
    data.table::data.table(
      call_chromosome = character(),
      position = integer(),
      call_strand = character(),
      sample_name = character(),
      modbase_code = character(),
      total_calls = integer(),
      modified_calls = integer(),
      unmodified_calls = integer(),
      frac_modified = numeric()
    )
  )

  write_empty_if_missing(
    paths$status,
    data.table::data.table(
      call_chromosome = character(),
      position = integer(),
      call_strand = character(),
      sample_name = character(),
      modbase_code = character(),
      total_calls = integer(),
      modified_calls = integer(),
      unmodified_calls = integer(),
      frac_modified = numeric(),
      status = character(),
      n_calls = integer()
    )
  )

  write_empty_if_missing(
    paths$aligned,
    data.table::data.table(
      call_chromosome = character(),
      position = integer(),
      call_strand = character(),
      sample_name = character(),
      modbase_code = character(),
      total_calls = integer(),
      modified_calls = integer(),
      unmodified_calls = integer(),
      frac_modified = numeric(),
      tss_uid = character(),
      chromosome = character(),
      tss_coordinate = integer(),
      tss_strand = character(),
      gene_name_or_tss_id = character(),
      tss_id = character(),
      promoter_start = integer(),
      promoter_end = integer(),
      rel_pos = integer()
    )
  )

  write_empty_if_missing(
    paths$aligned_status,
    data.table::data.table(
      call_chromosome = character(),
      position = integer(),
      call_strand = character(),
      sample_name = character(),
      modbase_code = character(),
      total_calls = integer(),
      modified_calls = integer(),
      unmodified_calls = integer(),
      frac_modified = numeric(),
      tss_uid = character(),
      chromosome = character(),
      tss_coordinate = integer(),
      tss_strand = character(),
      gene_name_or_tss_id = character(),
      tss_id = character(),
      promoter_start = integer(),
      promoter_end = integer(),
      rel_pos = integer(),
      status = character(),
      n_calls = integer()
    )
  )

  write_gz_tsv(summary_dt, paths$summary)
  write_gz_tsv(meta_dt, paths$meta)

  invisible(paths)
}
