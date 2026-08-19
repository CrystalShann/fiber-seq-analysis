# ---------------------------------------------------------------------------
# Part 1: cCRE setup, window building, nucleosome alignment
#
# Analogous to run_modbam_promoter_summary_part1.R but anchored on cCRE
# centers instead of TSS coordinates.  No strand flip: cCREs are unstranded,
# so rel_pos = feature_pos - center (positive = downstream of center).
#
# Parameters (override via env vars for batch submission):
#   SAMPLE            sample name   [AL10_bc2178_19130]
#   CHROM             chromosome    [chr1]
#   WINDOW_BP         flanking bp each side of cCRE center  [300]
#   BIN_SIZE          bp per metaprofile bin                [1]
# ---------------------------------------------------------------------------

source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/modbam_footprintR_functions.R")

default_sample <- Sys.getenv("SAMPLE",    unset = "AL10_bc2178_19130")
default_chrom  <- Sys.getenv("CHROM",     unset = "chr1")
window_bp      <- as.integer(Sys.getenv("WINDOW_BP", unset = "300"))
bin_size       <- as.integer(Sys.getenv("BIN_SIZE",  unset = "1"))

# ── Paths ─────────────────────────────────────────────────────────────────────
results_root_dir <- "/project/spott/cshan/fiber-seq/results/cCRE_summary/300bp_center"
threshold_dir    <- file.path(results_root_dir, "threshold_0.9")
output_dir       <- file.path(threshold_dir, default_sample)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bam_path <- file.path(
  "/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples",
  sprintf("%s.5mC.6mA.aligned.phased.bam", default_sample)
)

sample_extract_dir <- file.path(
  "/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results",
  default_sample,
  "extracted_results"
)

nuc_path <- file.path(
  sample_extract_dir,
  "nuc_by_chr",
  sprintf("%s.ft_extracted_nuc.%s.bed.gz", default_sample, default_chrom)
)

ccre_path <- "/project/spott/cshan/annotations/GRCh38-cCREs.bed"

# ── cCRE-specific helpers ─────────────────────────────────────────────────────

# All ENCODE Registry V4 classes (+ synthetic Enhancers = pELS + dELS)
ccre_classes_base <- c("PLS", "pELS", "dELS", "CA-TF", "CA-CTCF", "CA-H3K4me3", "CA", "TF")

read_ccre_table <- function(path, chrom = NULL, seqinfo = NULL, window_bp = 300L) {
  dt <- data.table::fread(
    path,
    header    = FALSE,
    col.names = c("chromosome", "start0", "end0", "ccre_id", "ccre_element_id", "ccre_class")
  )

  if (!is.null(chrom)) {
    dt <- dt[chromosome == chrom & ccre_class %in% ccre_classes_base]
  }

  if (!nrow(dt)) {
    stop("No cCRE rows remained after chromosome/class filtering.", call. = FALSE)
  }

  # Synthetic "Enhancers" = pELS + dELS relabelled
  enhancer_rows <- data.table::copy(dt[ccre_class %in% c("pELS", "dELS")])
  enhancer_rows[, ccre_class := "Enhancers"]
  dt_all <- data.table::rbindlist(list(dt, enhancer_rows), use.names = TRUE, fill = TRUE)

  dt_all[, row_id   := .I]
  dt_all[, center   := as.integer((start0 + end0) / 2L)]
  dt_all[, ccre_uid := paste(chromosome, center, ccre_class, ccre_id, row_id, sep = "|")]

  # Build ±window_bp windows around each center
  dt_all[, win_start := pmax(0L, center - window_bp)]
  dt_all[, win_end   := center + window_bp]

  # GRanges for the center windows (strand = "*"; cCREs are unstranded)
  wins_gr <- GenomicRanges::GRanges(
    seqnames = dt_all$chromosome,
    ranges   = IRanges::IRanges(dt_all$win_start + 1L, dt_all$win_end),
    strand   = "*"
  )
  if (!is.null(seqinfo)) {
    GenomeInfoDb::seqinfo(wins_gr) <- seqinfo[unique(dt_all$chromosome)]
  }
  wins_gr <- GenomicRanges::trim(wins_gr)

  # Update win_start / win_end after trimming (for edge-of-chromosome cCREs)
  dt_all[, win_start := as.integer(GenomicRanges::start(wins_gr)) - 1L]
  dt_all[, win_end   := as.integer(GenomicRanges::end(wins_gr))]

  list(
    ccre_table = dt_all,
    ccre_wins  = wins_gr
  )
}

# Align nucleosome feature blocks to cCRE center windows.
# No strand flip: rel_pos = feature_pos - center.
align_feature_blocks_to_ccre <- function(feature_dt, feature_gr, ccre_wins, ccre_table) {
  if (!nrow(feature_dt) || !length(feature_gr) || !length(ccre_wins)) {
    return(data.table::data.table())
  }

  hits <- GenomicRanges::findOverlaps(
    query   = feature_gr,
    subject = ccre_wins,
    ignore.strand = TRUE
  )

  if (!length(hits)) {
    return(data.table::data.table())
  }

  aligned <- cbind(
    feature_dt[S4Vectors::queryHits(hits)],
    ccre_table[S4Vectors::subjectHits(hits)]
  )

  # Unstranded: rel_pos = feature midpoint - cCRE center
  aligned[, rel_pos := feature_pos - center]
  aligned[]
}

# Nucleosome count summary per cCRE (analogous to summarize_nucleosomes_by_tss)
summarize_nucleosomes_by_ccre <- function(nuc_aligned, window_bp = 300L) {
  if (!nrow(nuc_aligned)) {
    return(data.table::data.table(ccre_uid = character()))
  }

  nuc_aligned[, {
    n_nuc <- .N
    .(
      structure_nucleosome_count             = n_nuc,
      structure_mean_nucleosome_size         = mean(feature_size),
      structure_median_nucleosome_size       = as.double(median(feature_size)),
      structure_nucleosome_density_per_kb    = n_nuc / ((2L * window_bp + 1L) / 1000),
      structure_core_nucleosome_count        = sum(abs(rel_pos) <= 75L),
      structure_mean_midpoint_offset         = mean(rel_pos),
      structure_mean_abs_midpoint_offset     = mean(abs(rel_pos))
    )
  }, by = ccre_uid]
}

# Nucleosome metaprofile per cCRE class (bin counts / n_ccre_in_class)
make_nucleosome_meta_by_class <- function(nuc_aligned, ccre_table, window_bp = 300L, bin_size = 1L) {
  bins        <- seq(-window_bp, window_bp, by = bin_size)
  classes     <- unique(ccre_table$ccre_class)
  full_grid   <- data.table::CJ(ccre_class = classes, bin = bins)

  if (!nrow(nuc_aligned)) {
    full_grid[, nucleosome_midpoints_per_ccre := 0]
    return(full_grid[])
  }

  n_per_class <- ccre_table[, .N, by = ccre_class]

  tmp <- data.table::copy(nuc_aligned)
  tmp[, bin       := as.integer(round(rel_pos / bin_size) * bin_size)]
  tmp[, ccre_class := gsub("\\|.*", "", gsub("^[^|]+\\|[^|]+\\|", "", ccre_uid))]
  # ccre_class is already in ccre_table columns; re-extract from aligned table
  if ("ccre_class" %in% names(nuc_aligned)) {
    tmp[, ccre_class := nuc_aligned$ccre_class[.I]]
  }

  meta <- tmp[, .N, by = .(ccre_class, bin)]
  meta <- merge(full_grid, meta, by = c("ccre_class", "bin"), all.x = TRUE, sort = TRUE)
  meta[is.na(N), N := 0L]
  meta <- merge(meta, n_per_class, by = "ccre_class", all.x = TRUE)
  meta[, nucleosome_midpoints_per_ccre := N / N.y]
  meta[, c("N", "N.y") := NULL]
  meta[]
}

# Modification metaprofile per cCRE class
make_modification_meta_by_class <- function(aligned_dt, ccre_table, label,
                                            window_bp = 300L, bin_size = 1L) {
  bins      <- seq(-window_bp, window_bp, by = bin_size)
  classes   <- unique(ccre_table$ccre_class)
  full_grid <- data.table::CJ(ccre_class = classes, bin = bins)

  label_value <- label

  if (!nrow(aligned_dt)) {
    full_grid[, `:=`(
      modified_calls   = 0L,
      total_calls      = 0L,
      fraction_modified = NA_real_,
      modified_calls_per_ccre = 0,
      label            = label_value
    )]
    return(full_grid[])
  }

  n_per_class <- ccre_table[, .N, by = ccre_class]

  tmp <- data.table::copy(aligned_dt)
  tmp[, bin := as.integer(round(rel_pos / bin_size) * bin_size)]

  meta <- tmp[, .(
    modified_calls = sum(modified_calls),
    total_calls    = sum(total_calls)
  ), by = .(ccre_class, bin)]

  meta <- merge(full_grid, meta, by = c("ccre_class", "bin"), all.x = TRUE, sort = TRUE)
  meta[is.na(modified_calls), modified_calls := 0L]
  meta[is.na(total_calls),    total_calls    := 0L]
  meta[, fraction_modified     := data.table::fifelse(total_calls > 0, modified_calls / total_calls, NA_real_)]
  meta <- merge(meta, n_per_class, by = "ccre_class", all.x = TRUE)
  meta[, modified_calls_per_ccre := modified_calls / N]
  meta[, N     := NULL]
  meta[, label := label_value]
  meta[]
}

# ── Check existing outputs ─────────────────────────────────────────────────────
out_ccre_table   <- file.path(output_dir, sprintf("%s_%s_ccre_table.tsv.gz",            default_sample, default_chrom))
out_nuc_aligned  <- file.path(output_dir, sprintf("%s_%s_nucleosome_ccre_aligned.tsv.gz", default_sample, default_chrom))
out_nuc_summary  <- file.path(output_dir, sprintf("%s_%s_ccre_structure.tsv.gz",          default_sample, default_chrom))
out_nuc_meta     <- file.path(output_dir, sprintf("%s_%s_nucleosome_meta.tsv.gz",          default_sample, default_chrom))

if (all(file.exists(c(out_ccre_table, out_nuc_aligned, out_nuc_summary, out_nuc_meta)))) {
  message("Part 1 outputs already exist; skipping nucleosome stage.")
  bam_seqinfo <- read_bam_seqinfo(bam_path)
  ccre_objects <- read_ccre_table(ccre_path, chrom = default_chrom, seqinfo = bam_seqinfo, window_bp = window_bp)
  ccre_table   <- ccre_objects$ccre_table
  ccre_wins    <- ccre_objects$ccre_wins
} else {
  assert_file_exists(bam_path,  "BAM")
  assert_file_exists(paste0(bam_path, ".bai"), "BAM index")
  assert_file_exists(nuc_path,  "nucleosome BED12")
  assert_file_exists(ccre_path, "cCRE BED")

  # ── Setup ──────────────────────────────────────────────────────────────────
  message("Reading BAM seqinfo ...")
  bam_seqinfo <- read_bam_seqinfo(bam_path)

  message("Reading cCRE annotations ...")
  ccre_objects <- read_ccre_table(ccre_path, chrom = default_chrom, seqinfo = bam_seqinfo, window_bp = window_bp)
  ccre_table   <- ccre_objects$ccre_table
  ccre_wins    <- ccre_objects$ccre_wins

  message(sprintf(
    "cCRE table: %d rows across %d classes on %s (including synthetic Enhancers)",
    nrow(ccre_table), length(unique(ccre_table$ccre_class)), default_chrom
  ))

  write_gz_tsv(ccre_table, out_ccre_table)

  # ── Nucleosome alignment ───────────────────────────────────────────────────
  message("Reading nucleosome BED12 ...")
  nuc_bed12 <- read_fiber_bed12(nuc_path)

  message("Expanding nucleosome feature blocks ...")
  nuc_features <- expand_ft_bed12_features(
    bed12                = nuc_bed12,
    feature_name         = "nucleosome",
    drop_terminal_blocks = TRUE
  )
  nuc_gr <- feature_blocks_to_granges(nuc_features)

  message("Aligning nucleosomes to cCRE windows ...")
  nuc_ccre_aligned <- align_feature_blocks_to_ccre(
    feature_dt = nuc_features,
    feature_gr = nuc_gr,
    ccre_wins  = ccre_wins,
    ccre_table = ccre_table
  )

  message(sprintf("Aligned %d nucleosome–cCRE pairs.", nrow(nuc_ccre_aligned)))

  nuc_summary <- summarize_nucleosomes_by_ccre(nuc_ccre_aligned, window_bp = window_bp)

  # ccre_class is already a column in nuc_ccre_aligned (from ccre_table merge)
  nuc_meta <- make_nucleosome_meta_by_class(
    nuc_aligned = nuc_ccre_aligned,
    ccre_table  = ccre_table,
    window_bp   = window_bp,
    bin_size    = bin_size
  )

  write_gz_tsv(nuc_ccre_aligned, out_nuc_aligned)
  write_gz_tsv(nuc_summary,      out_nuc_summary)
  write_gz_tsv(nuc_meta,         out_nuc_meta)

  rm(nuc_bed12, nuc_features, nuc_gr, nuc_meta)
  gc()

  message("Part 1 complete.")
}
