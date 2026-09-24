# ---------------------------------------------------------------------------
# Part 2: m6A accessibility + 5mC methylation aligned to cCRE centers
#
# Reads Part 1 output (ccre_table TSV), queries the BAM for m6A and 5mC
# calls in chunks, aligns each position to its cCRE center
# (rel_pos = position - center, no strand flip), accumulates bin-level
# counts, and writes per-class metaprofiles.
#
# Parameters (override via env vars):
#   SAMPLE             [AL10_bc2178_19130]
#   CHROM              [chr1]
#   WINDOW_BP          flanking bp each side of cCRE center  [300]
#   BIN_SIZE           bp per metaprofile bin                [1]
#   MOD_PROB_THRESHOLD probability cutoff for modified call  [0.9]
#   CCRE_CHUNK_SIZE    cCRE rows per BAM query chunk         [500]
# ---------------------------------------------------------------------------

source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/modbam_footprintR_functions.R")

default_sample     <- Sys.getenv("SAMPLE",             unset = "AL10_bc2178_19130")
default_chrom      <- Sys.getenv("CHROM",              unset = "chr1")
window_bp          <- as.integer(Sys.getenv("WINDOW_BP",          unset = "300"))
bin_size           <- as.integer(Sys.getenv("BIN_SIZE",           unset = "1"))
mod_prob_threshold <- as.numeric(Sys.getenv("MOD_PROB_THRESHOLD", unset = "0.9"))
ccre_chunk_size    <- as.integer(Sys.getenv("CCRE_CHUNK_SIZE",    unset = "500"))
if (is.na(ccre_chunk_size) || ccre_chunk_size < 1L) ccre_chunk_size <- 500L

# ── Paths ─────────────────────────────────────────────────────────────────────
results_root_dir <- "/project/spott/cshan/fiber-seq/results/cCRE_summary/300bp_center"
output_dir       <- file.path(results_root_dir, "threshold_0.9", default_sample)

bam_path  <- file.path(
  "/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples",
  sprintf("%s.5mC.6mA.aligned.phased.bam", default_sample)
)

out_ccre_table  <- file.path(output_dir, sprintf("%s_%s_ccre_table.tsv.gz",                       default_sample, default_chrom))
out_m6a_summary <- file.path(output_dir, sprintf("%s_%s_ccre_accessibility.tsv.gz",                default_sample, default_chrom))
out_m6a_meta    <- file.path(output_dir, sprintf("%s_%s_m6A_meta.tsv.gz",                          default_sample, default_chrom))
out_5mc_summary <- file.path(output_dir, sprintf("%s_%s_ccre_methylation.tsv.gz",                  default_sample, default_chrom))
out_5mc_meta    <- file.path(output_dir, sprintf("%s_%s_5mC_meta.tsv.gz",                          default_sample, default_chrom))
out_multiomic   <- file.path(output_dir, sprintf("%s_%s_ccre_structure_accessibility_methylation.tsv.gz", default_sample, default_chrom))

assert_file_exists(bam_path,        "BAM")
assert_file_exists(paste0(bam_path, ".bai"), "BAM index")
assert_file_exists(out_ccre_table,  "Part 1 cCRE table — run part1 first")

# ── Load Part 1 cCRE table ────────────────────────────────────────────────────
message("Loading cCRE table from Part 1 ...")
ccre_table <- data.table::fread(out_ccre_table)

bam_seqinfo    <- read_bam_seqinfo(bam_path)
bamfiles       <- setNames(bam_path, default_sample)
bpparam_modbam <- BiocParallel::SerialParam()

# Rebuild GRanges from saved (trimmed) window coordinates
ccre_wins <- GenomicRanges::GRanges(
  seqnames = ccre_table$chromosome,
  ranges   = IRanges::IRanges(ccre_table$win_start + 1L, ccre_table$win_end),
  strand   = "*"
)
GenomeInfoDb::seqinfo(ccre_wins) <- bam_seqinfo[unique(ccre_table$chromosome)]

# ── Helper: align per-position BAM calls to cCRE center (no strand flip) ──────
align_positions_to_ccre <- function(position_dt, chunk_wins, chunk_table) {
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
  # Unstranded: distance from cCRE center; positive = downstream
  aligned[, rel_pos := position - center]
  aligned[]
}

# ── Helper: per-cCRE modification fraction summary ────────────────────────────
summarize_modality_by_ccre <- function(aligned_dt, prefix) {
  if (!nrow(aligned_dt)) return(data.table::data.table(ccre_uid = character()))

  out <- aligned_dt[, {
    total    <- sum(total_calls)
    modified <- sum(modified_calls)
    .(
      modified_calls          = modified,
      unmodified_calls        = sum(unmodified_calls),
      total_calls             = total,
      positions_with_coverage = sum(total_calls > 0),
      fraction_modified       = if (total > 0) modified / total else NA_real_
    )
  }, by = ccre_uid]

  data.table::setnames(out,
    old = c("modified_calls","unmodified_calls","total_calls",
            "positions_with_coverage","fraction_modified"),
    new = sprintf(c("%s_modified_calls","%s_unmodified_calls","%s_total_calls",
                    "%s_positions_with_coverage","%s_fraction_modified"), prefix)
  )
  out[]
}

# ── Helper: build final metaprofile from accumulated bin counts ────────────────
build_meta_from_acc <- function(meta_acc, ccre_table, label, window_bp, bin_size) {
  all_bins    <- seq(-window_bp, window_bp, by = bin_size)
  all_classes <- unique(ccre_table$ccre_class)
  full_grid   <- data.table::CJ(ccre_class = all_classes, bin = all_bins)
  n_per_class <- ccre_table[, .N, by = ccre_class]

  meta <- merge(full_grid, meta_acc, by = c("ccre_class", "bin"), all.x = TRUE, sort = TRUE)
  meta[is.na(modified_calls), modified_calls := 0L]
  meta[is.na(total_calls),    total_calls    := 0L]
  meta[, fraction_modified    := data.table::fifelse(total_calls > 0,
                                    modified_calls / total_calls, NA_real_)]
  meta <- merge(meta, n_per_class, by = "ccre_class", all.x = TRUE)
  meta[, modified_calls_per_ccre := modified_calls / N]
  meta[, c("N") := NULL]
  meta[, label := label]
  meta[]
}

# ── Generic chunked BAM processor ─────────────────────────────────────────────
process_modality_chunked <- function(modbase_code, summary_prefix,
                                     mod_label, meta_label,
                                     out_summary, out_meta) {

  if (all(file.exists(c(out_summary, out_meta)))) {
    message(mod_label, " outputs already exist; skipping.")
    return(invisible(NULL))
  }
  if (file.exists(out_summary)) file.remove(out_summary)
  if (file.exists(out_meta))    file.remove(out_meta)

  n_ccre       <- nrow(ccre_table)
  chunk_starts <- seq.int(1L, n_ccre, by = ccre_chunk_size)
  message(sprintf("Processing %s: %d chunks of %d cCREs each ...",
                  mod_label, length(chunk_starts), ccre_chunk_size))

  # Accumulator: bin counts per (ccre_class, bin)
  all_bins    <- seq(-window_bp, window_bp, by = bin_size)
  all_classes <- unique(ccre_table$ccre_class)
  meta_acc    <- data.table::CJ(ccre_class = all_classes, bin = all_bins)
  meta_acc[, `:=`(modified_calls = 0L, total_calls = 0L)]

  summary_list <- vector("list", length(chunk_starts))

  for (i in seq_along(chunk_starts)) {
    from <- chunk_starts[i]
    to   <- min(from + ccre_chunk_size - 1L, n_ccre)

    chunk_table <- ccre_table[from:to]
    chunk_wins  <- ccre_wins[from:to]

    message(sprintf("  %s chunk %d/%d (cCREs %d-%d)",
                    mod_label, i, length(chunk_starts), from, to))

    chunk_pos <- read_modbam_summary_positions(
      bamfiles         = bamfiles,
      sample_name      = default_sample,
      promoter_windows = chunk_wins,
      modbase_code     = modbase_code,
      seqinfo          = bam_seqinfo,
      bpparam          = bpparam_modbam,
      modProbThreshold = mod_prob_threshold
    )

    if (!nrow(chunk_pos)) {
      summary_list[[i]] <- data.table::data.table(ccre_uid = character())
      next
    }

    chunk_aligned <- align_positions_to_ccre(chunk_pos, chunk_wins, chunk_table)

    summary_list[[i]] <- summarize_modality_by_ccre(chunk_aligned, summary_prefix)

    # Accumulate bin counts across chunks
    if (nrow(chunk_aligned)) {
      bin_chunk <- data.table::copy(chunk_aligned)
      bin_chunk[, bin := as.integer(round(rel_pos / bin_size) * bin_size)]
      bin_chunk <- bin_chunk[bin >= -window_bp & bin <= window_bp,
                             .(modified_calls = sum(modified_calls),
                               total_calls    = sum(total_calls)),
                             by = .(ccre_class, bin)]

      meta_acc[bin_chunk, on = c("ccre_class", "bin"), `:=`(
        modified_calls = modified_calls + i.modified_calls,
        total_calls    = total_calls    + i.total_calls
      )]
    }

    rm(chunk_pos, chunk_aligned)
    gc()
  }

  # ── Collapse per-cCRE summary ─────────────────────────────────────────────
  summary_dt <- data.table::rbindlist(summary_list, use.names = TRUE, fill = TRUE)

  mod_col     <- sprintf("%s_modified_calls",          summary_prefix)
  unmod_col   <- sprintf("%s_unmodified_calls",        summary_prefix)
  total_col   <- sprintf("%s_total_calls",             summary_prefix)
  cov_col     <- sprintf("%s_positions_with_coverage", summary_prefix)
  frac_col    <- sprintf("%s_fraction_modified",       summary_prefix)

  if (nrow(summary_dt)) {
    summary_dt <- summary_dt[, .(
      modified_calls          = sum(get(mod_col),   na.rm = TRUE),
      unmodified_calls        = sum(get(unmod_col), na.rm = TRUE),
      total_calls             = sum(get(total_col), na.rm = TRUE),
      positions_with_coverage = sum(get(cov_col),   na.rm = TRUE)
    ), by = ccre_uid]
    summary_dt[, fraction_modified := data.table::fifelse(
      total_calls > 0, modified_calls / total_calls, NA_real_
    )]
    data.table::setnames(summary_dt,
      old = c("modified_calls","unmodified_calls","total_calls",
              "positions_with_coverage","fraction_modified"),
      new = c(mod_col, unmod_col, total_col, cov_col, frac_col)
    )
  }

  meta_dt <- build_meta_from_acc(meta_acc, ccre_table, meta_label, window_bp, bin_size)

  write_gz_tsv(summary_dt, out_summary)
  write_gz_tsv(meta_dt,    out_meta)
  message(mod_label, " done.")
  invisible(NULL)
}

# ── Run both modalities ────────────────────────────────────────────────────────
process_modality_chunked("a", "accessibility", "m6A", "Accessibility (m6A)",
                         out_m6a_summary, out_m6a_meta)

process_modality_chunked("m", "methylation",   "5mC", "Methylation (5mC)",
                         out_5mc_summary, out_5mc_meta)

# ── Join per-cCRE multi-omic summary (structure + accessibility + methylation) ─
if (!file.exists(out_multiomic)) {
  message("Building cCRE multi-omic summary ...")

  out_nuc_summary <- file.path(output_dir,
    sprintf("%s_%s_ccre_structure.tsv.gz", default_sample, default_chrom))

  m6a_summary <- data.table::fread(out_m6a_summary)
  m5c_summary <- data.table::fread(out_5mc_summary)

  ccre_key_cols <- intersect(
    c("ccre_uid","chromosome","center","ccre_class","ccre_id","win_start","win_end"),
    names(ccre_table)
  )
  ccre_multiomic <- data.table::copy(ccre_table[, ..ccre_key_cols])
  ccre_multiomic <- merge(ccre_multiomic, m6a_summary, by = "ccre_uid", all.x = TRUE)
  ccre_multiomic <- merge(ccre_multiomic, m5c_summary, by = "ccre_uid", all.x = TRUE)

  # Optionally join nucleosome structure if Part 1 produced it
  if (file.exists(out_nuc_summary)) {
    nuc_summary <- data.table::fread(out_nuc_summary)
    ccre_multiomic <- merge(ccre_multiomic, nuc_summary, by = "ccre_uid", all.x = TRUE)
  }

  # Fill missing numeric columns
  zero_cols <- intersect(
    c("structure_nucleosome_count","structure_nucleosome_density_per_kb",
      "structure_core_nucleosome_count",
      "accessibility_modified_calls","accessibility_unmodified_calls",
      "accessibility_total_calls","accessibility_positions_with_coverage",
      "methylation_modified_calls","methylation_unmodified_calls",
      "methylation_total_calls","methylation_positions_with_coverage"),
    names(ccre_multiomic)
  )
  for (col in zero_cols) ccre_multiomic[is.na(get(col)), (col) := 0]

  na_cols <- intersect(
    c("structure_mean_nucleosome_size","structure_median_nucleosome_size",
      "structure_mean_midpoint_offset","structure_mean_abs_midpoint_offset",
      "accessibility_fraction_modified","methylation_fraction_modified"),
    names(ccre_multiomic)
  )
  for (col in na_cols) {
    ccre_multiomic[is.infinite(get(col)) | is.nan(get(col)), (col) := NA_real_]
  }

  write_gz_tsv(ccre_multiomic, out_multiomic)
  message("Multi-omic summary written to: ", out_multiomic)
}

message("Part 2 complete.")
