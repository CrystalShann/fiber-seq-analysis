# ---------------------------------------------------------------------------
# Part 2: m6A accessibility aligned to silencer centers (m6A only)
#
# Parameters (override via env vars):
#   SAMPLE             [AL10_bc2178_19130]
#   CHROM              [chr1]
#   WINDOW_BP          [300]
#   BIN_SIZE           [1]
#   MOD_PROB_THRESHOLD [0.9]
#   SILENCER_CHUNK_SIZE [500]
#   SILENCER_CLASS     one of REST-Enhancers, REST-Silencers,
#                      STARR-Silencers.Robust, STARR-Silencers.Stringent,
#                      or ALL [ALL]
# ---------------------------------------------------------------------------

source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/modbam_footprintR_functions.R")

default_sample      <- Sys.getenv("SAMPLE",              unset = "AL10_bc2178_19130")
default_chrom       <- Sys.getenv("CHROM",               unset = "chr1")
window_bp           <- as.integer(Sys.getenv("WINDOW_BP",           unset = "300"))
bin_size            <- as.integer(Sys.getenv("BIN_SIZE",            unset = "1"))
mod_prob_threshold  <- as.numeric(Sys.getenv("MOD_PROB_THRESHOLD",  unset = "0.9"))
silencer_chunk_size <- as.integer(Sys.getenv("SILENCER_CHUNK_SIZE", unset = "500"))
silencer_class      <- Sys.getenv("SILENCER_CLASS",      unset = "ALL")
output_scope        <- Sys.getenv("OUTPUT_SCOPE",        unset = "sample_chr")
if (is.na(silencer_chunk_size) || silencer_chunk_size < 1L) silencer_chunk_size <- 500L

if (silencer_class == "ALL") {
  class_tag <- "all4classes"
} else {
  class_tag <- gsub("[^A-Za-z0-9._-]", "_", silencer_class)
}

# ── Paths ─────────────────────────────────────────────────────────────────────
results_root_dir <- file.path(
  "/project/spott/cshan/fiber-seq/results/cCRE_summary/silencer_summary/m6a_modbam",
  output_scope
)
output_dir       <- file.path(results_root_dir, "threshold_0.9", default_sample, class_tag)

bam_path <- file.path(
  "/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples",
  sprintf("%s.5mC.6mA.aligned.phased.bam", default_sample)
)

out_silencer_table <- file.path(output_dir, sprintf("%s_%s_%s_silencer_table.tsv.gz", default_sample, default_chrom, class_tag))
out_m6a_summary    <- file.path(output_dir, sprintf("%s_%s_%s_silencer_accessibility.tsv.gz", default_sample, default_chrom, class_tag))
out_m6a_meta       <- file.path(output_dir, sprintf("%s_%s_%s_m6A_meta.tsv.gz", default_sample, default_chrom, class_tag))
out_multiomic      <- file.path(output_dir, sprintf("%s_%s_%s_silencer_structure_accessibility.tsv.gz", default_sample, default_chrom, class_tag))

assert_file_exists(bam_path, "BAM")
assert_file_exists(paste0(bam_path, ".bai"), "BAM index")
assert_file_exists(out_silencer_table, "Part 1 silencer table — run part1 first")

message("Loading silencer table from Part 1 ...")
silencer_table <- data.table::fread(out_silencer_table)

bam_seqinfo    <- read_bam_seqinfo(bam_path)
bamfiles       <- setNames(bam_path, default_sample)
bpparam_modbam <- BiocParallel::SerialParam()

silencer_wins <- GenomicRanges::GRanges(
  seqnames = silencer_table$chromosome,
  ranges   = IRanges::IRanges(silencer_table$win_start + 1L, silencer_table$win_end),
  strand   = "*"
)
GenomeInfoDb::seqinfo(silencer_wins) <- bam_seqinfo[unique(silencer_table$chromosome)]

align_positions_to_silencer <- function(position_dt, chunk_wins, chunk_table) {
  if (!nrow(position_dt)) return(data.table::data.table())

  pos_gr <- GenomicRanges::GRanges(
    seqnames = position_dt$call_chromosome,
    ranges   = IRanges::IRanges(position_dt$position, position_dt$position),
    strand   = "*"
  )

  hits <- GenomicRanges::findOverlaps(pos_gr, chunk_wins, ignore.strand = TRUE)
  if (!length(hits)) return(data.table::data.table())

  aligned <- cbind(
    position_dt[S4Vectors::queryHits(hits)],
    chunk_table[S4Vectors::subjectHits(hits)]
  )
  aligned[, rel_pos := position - center]
  aligned[]
}

summarize_modality_by_silencer <- function(aligned_dt, prefix) {
  if (!nrow(aligned_dt)) return(data.table::data.table(silencer_uid = character()))

  out <- aligned_dt[, {
    total <- sum(total_calls)
    modified <- sum(modified_calls)
    .(
      modified_calls          = modified,
      unmodified_calls        = sum(unmodified_calls),
      total_calls             = total,
      positions_with_coverage = sum(total_calls > 0),
      fraction_modified       = if (total > 0) modified / total else NA_real_
    )
  }, by = silencer_uid]

  data.table::setnames(
    out,
    old = c("modified_calls", "unmodified_calls", "total_calls", "positions_with_coverage", "fraction_modified"),
    new = sprintf(c("%s_modified_calls", "%s_unmodified_calls", "%s_total_calls", "%s_positions_with_coverage", "%s_fraction_modified"), prefix)
  )
  out[]
}

build_meta_from_acc <- function(meta_acc, silencer_table, label, window_bp, bin_size) {
  all_bins    <- seq(-window_bp, window_bp, by = bin_size)
  all_classes <- unique(silencer_table$silencer_class)
  full_grid   <- data.table::CJ(silencer_class = all_classes, bin = all_bins)
  n_per_class <- silencer_table[, .N, by = silencer_class]

  meta <- merge(full_grid, meta_acc, by = c("silencer_class", "bin"), all.x = TRUE, sort = TRUE)
  meta[is.na(modified_calls), modified_calls := 0L]
  meta[is.na(total_calls), total_calls := 0L]
  meta[, fraction_modified := data.table::fifelse(total_calls > 0, modified_calls / total_calls, NA_real_)]
  meta <- merge(meta, n_per_class, by = "silencer_class", all.x = TRUE)
  meta[, modified_calls_per_silencer := modified_calls / N]
  meta[, N := NULL]
  meta[, label := label]
  meta[]
}

process_m6a_chunked <- function() {
  if (all(file.exists(c(out_m6a_summary, out_m6a_meta)))) {
    message("m6A outputs already exist; skipping.")
    return(invisible(NULL))
  }
  if (file.exists(out_m6a_summary)) file.remove(out_m6a_summary)
  if (file.exists(out_m6a_meta)) file.remove(out_m6a_meta)

  n_silencer   <- nrow(silencer_table)
  chunk_starts <- seq.int(1L, n_silencer, by = silencer_chunk_size)

  all_bins    <- seq(-window_bp, window_bp, by = bin_size)
  all_classes <- unique(silencer_table$silencer_class)
  meta_acc    <- data.table::CJ(silencer_class = all_classes, bin = all_bins)
  meta_acc[, `:=`(modified_calls = 0L, total_calls = 0L)]

  summary_list <- vector("list", length(chunk_starts))

  for (i in seq_along(chunk_starts)) {
    from <- chunk_starts[i]
    to   <- min(from + silencer_chunk_size - 1L, n_silencer)

    chunk_table <- silencer_table[from:to]
    chunk_wins  <- silencer_wins[from:to]

    message(sprintf("  m6A chunk %d/%d (silencers %d-%d)", i, length(chunk_starts), from, to))

    chunk_pos <- read_modbam_summary_positions(
      bamfiles = bamfiles,
      sample_name = default_sample,
      promoter_windows = chunk_wins,
      modbase_code = "a",
      seqinfo = bam_seqinfo,
      bpparam = bpparam_modbam,
      modProbThreshold = mod_prob_threshold
    )

    if (!nrow(chunk_pos)) {
      summary_list[[i]] <- data.table::data.table(silencer_uid = character())
      next
    }

    chunk_aligned <- align_positions_to_silencer(chunk_pos, chunk_wins, chunk_table)
    summary_list[[i]] <- summarize_modality_by_silencer(chunk_aligned, "accessibility")

    if (nrow(chunk_aligned)) {
      bin_chunk <- data.table::copy(chunk_aligned)
      bin_chunk[, bin := as.integer(round(rel_pos / bin_size) * bin_size)]
      bin_chunk <- bin_chunk[bin >= -window_bp & bin <= window_bp,
        .(modified_calls = sum(modified_calls), total_calls = sum(total_calls)),
        by = .(silencer_class, bin)]

      meta_acc[bin_chunk, on = c("silencer_class", "bin"), `:=`(
        modified_calls = modified_calls + i.modified_calls,
        total_calls = total_calls + i.total_calls
      )]
    }

    rm(chunk_pos, chunk_aligned)
    gc()
  }

  summary_dt <- data.table::rbindlist(summary_list, use.names = TRUE, fill = TRUE)
  if (nrow(summary_dt)) {
    summary_dt <- summary_dt[, .(
      modified_calls = sum(accessibility_modified_calls, na.rm = TRUE),
      unmodified_calls = sum(accessibility_unmodified_calls, na.rm = TRUE),
      total_calls = sum(accessibility_total_calls, na.rm = TRUE),
      positions_with_coverage = sum(accessibility_positions_with_coverage, na.rm = TRUE)
    ), by = silencer_uid]
    summary_dt[, fraction_modified := data.table::fifelse(total_calls > 0, modified_calls / total_calls, NA_real_)]

    data.table::setnames(
      summary_dt,
      old = c("modified_calls", "unmodified_calls", "total_calls", "positions_with_coverage", "fraction_modified"),
      new = c(
        "accessibility_modified_calls", "accessibility_unmodified_calls", "accessibility_total_calls",
        "accessibility_positions_with_coverage", "accessibility_fraction_modified"
      )
    )
  }

  meta_dt <- build_meta_from_acc(meta_acc, silencer_table, "Accessibility (m6A)", window_bp, bin_size)

  write_gz_tsv(summary_dt, out_m6a_summary)
  write_gz_tsv(meta_dt, out_m6a_meta)
}

process_m6a_chunked()

if (!file.exists(out_multiomic)) {
  message("Building silencer structure + accessibility summary ...")

  out_nuc_summary <- file.path(output_dir, sprintf("%s_%s_%s_silencer_structure.tsv.gz", default_sample, default_chrom, class_tag))
  m6a_summary <- data.table::fread(out_m6a_summary)

  key_cols <- intersect(
    c("silencer_uid", "chromosome", "center", "silencer_class", "region_id", "win_start", "win_end"),
    names(silencer_table)
  )
  silencer_multiomic <- data.table::copy(silencer_table[, ..key_cols])
  silencer_multiomic <- merge(silencer_multiomic, m6a_summary, by = "silencer_uid", all.x = TRUE)

  if (file.exists(out_nuc_summary)) {
    nuc_summary <- data.table::fread(out_nuc_summary)
    silencer_multiomic <- merge(silencer_multiomic, nuc_summary, by = "silencer_uid", all.x = TRUE)
  }

  zero_cols <- intersect(
    c(
      "structure_nucleosome_count", "structure_nucleosome_density_per_kb", "structure_core_nucleosome_count",
      "accessibility_modified_calls", "accessibility_unmodified_calls", "accessibility_total_calls",
      "accessibility_positions_with_coverage"
    ),
    names(silencer_multiomic)
  )
  for (col in zero_cols) silencer_multiomic[is.na(get(col)), (col) := 0]

  na_cols <- intersect(
    c(
      "structure_mean_nucleosome_size", "structure_median_nucleosome_size",
      "structure_mean_midpoint_offset", "structure_mean_abs_midpoint_offset",
      "accessibility_fraction_modified"
    ),
    names(silencer_multiomic)
  )
  for (col in na_cols) {
    silencer_multiomic[is.infinite(get(col)) | is.nan(get(col)), (col) := NA_real_]
  }

  write_gz_tsv(silencer_multiomic, out_multiomic)
}

message("Part 2 complete.")