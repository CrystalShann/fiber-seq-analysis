# ---------------------------------------------------------------------------
# Part 1: Silencer setup, window building, nucleosome alignment
#
# Parameters (override via env vars for batch submission):
#   SAMPLE            sample name                              [AL10_bc2178_19130]
#   CHROM             chromosome                               [chr1]
#   WINDOW_BP         flanking bp each side of silencer center [300]
#   BIN_SIZE          bp per metaprofile bin                   [1]
#   SILENCER_CLASS    one of REST-Enhancers, REST-Silencers,
#                     STARR-Silencers.Robust, STARR-Silencers.Stringent,
#                     or ALL                                   [ALL]
# ---------------------------------------------------------------------------

source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/modbam_footprintR_functions.R")

default_sample  <- Sys.getenv("SAMPLE",         unset = "AL10_bc2178_19130")
default_chrom   <- Sys.getenv("CHROM",          unset = "chr1")
window_bp       <- as.integer(Sys.getenv("WINDOW_BP",      unset = "300"))
bin_size        <- as.integer(Sys.getenv("BIN_SIZE",       unset = "1"))
silencer_class  <- Sys.getenv("SILENCER_CLASS", unset = "ALL")
output_scope    <- Sys.getenv("OUTPUT_SCOPE",   unset = "sample_chr")

silencer_files <- c(
  "REST-Enhancers"            = "/project/spott/cshan/annotations/human_silencers/REST-Enhancers.bed",
  "REST-Silencers"            = "/project/spott/cshan/annotations/human_silencers/REST-Silencers.bed",
  "STARR-Silencers.Robust"    = "/project/spott/cshan/annotations/human_silencers/STARR-Silencers.Robust.bed",
  "STARR-Silencers.Stringent" = "/project/spott/cshan/annotations/human_silencers/STARR-Silencers.Stringent.bed"
)

if (silencer_class == "ALL") {
  selected_classes <- names(silencer_files)
  class_tag <- "all4classes"
} else {
  if (!silencer_class %in% names(silencer_files)) {
    stop("Invalid SILENCER_CLASS value.", call. = FALSE)
  }
  selected_classes <- silencer_class
  class_tag <- gsub("[^A-Za-z0-9._-]", "_", silencer_class)
}

# ── Paths ─────────────────────────────────────────────────────────────────────
results_root_dir <- file.path(
  "/project/spott/cshan/fiber-seq/results/cCRE_summary/silencer_summary/m6a_modbam",
  output_scope
)
threshold_dir    <- file.path(results_root_dir, "threshold_0.9")
output_dir       <- file.path(threshold_dir, default_sample, class_tag)
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

# ── Silencer helpers ─────────────────────────────────────────────────────────
read_silencer_table <- function(file_map, selected_classes, chrom = NULL, seqinfo = NULL, window_bp = 300L) {
  dt_list <- lapply(selected_classes, function(cls) {
    path <- file_map[[cls]]
    data.table::fread(
      path,
      header = FALSE,
      col.names = c("chromosome", "start0", "end0", "region_id", "element_id", "annotation_class")
    )[, silencer_class := cls]
  })

  dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)

  if (!is.null(chrom)) {
    dt <- dt[chromosome == chrom]
  }

  if (!nrow(dt)) {
    stop("No silencer rows remained after chromosome/class filtering.", call. = FALSE)
  }

  dt[, row_id := .I]
  dt[, center := as.integer((start0 + end0) / 2L)]
  dt[, silencer_uid := paste(chromosome, center, silencer_class, region_id, row_id, sep = "|")]

  dt[, win_start := pmax(0L, center - window_bp)]
  dt[, win_end   := center + window_bp]

  wins_gr <- GenomicRanges::GRanges(
    seqnames = dt$chromosome,
    ranges   = IRanges::IRanges(dt$win_start + 1L, dt$win_end),
    strand   = "*"
  )
  if (!is.null(seqinfo)) {
    GenomeInfoDb::seqinfo(wins_gr) <- seqinfo[unique(dt$chromosome)]
  }
  wins_gr <- GenomicRanges::trim(wins_gr)

  dt[, win_start := as.integer(GenomicRanges::start(wins_gr)) - 1L]
  dt[, win_end   := as.integer(GenomicRanges::end(wins_gr))]

  list(silencer_table = dt, silencer_wins = wins_gr)
}

align_feature_blocks_to_silencer <- function(feature_dt, feature_gr, silencer_wins, silencer_table) {
  if (!nrow(feature_dt) || !length(feature_gr) || !length(silencer_wins)) {
    return(data.table::data.table())
  }

  hits <- GenomicRanges::findOverlaps(
    query = feature_gr,
    subject = silencer_wins,
    ignore.strand = TRUE
  )
  if (!length(hits)) {
    return(data.table::data.table())
  }

  aligned <- cbind(
    feature_dt[S4Vectors::queryHits(hits)],
    silencer_table[S4Vectors::subjectHits(hits)]
  )
  aligned[, rel_pos := feature_pos - center]
  aligned[]
}

summarize_nucleosomes_by_silencer <- function(nuc_aligned, window_bp = 300L) {
  if (!nrow(nuc_aligned)) {
    return(data.table::data.table(silencer_uid = character()))
  }

  nuc_aligned[, {
    n_nuc <- .N
    .(
      structure_nucleosome_count         = n_nuc,
      structure_mean_nucleosome_size     = mean(feature_size),
      structure_median_nucleosome_size   = as.double(median(feature_size)),
      structure_nucleosome_density_per_kb = n_nuc / ((2L * window_bp + 1L) / 1000),
      structure_core_nucleosome_count    = sum(abs(rel_pos) <= 75L),
      structure_mean_midpoint_offset     = mean(rel_pos),
      structure_mean_abs_midpoint_offset = mean(abs(rel_pos))
    )
  }, by = silencer_uid]
}

make_nucleosome_meta_by_class <- function(nuc_aligned, silencer_table, window_bp = 300L, bin_size = 1L) {
  bins      <- seq(-window_bp, window_bp, by = bin_size)
  classes   <- unique(silencer_table$silencer_class)
  full_grid <- data.table::CJ(silencer_class = classes, bin = bins)

  if (!nrow(nuc_aligned)) {
    full_grid[, nucleosome_midpoints_per_silencer := 0]
    return(full_grid[])
  }

  n_per_class <- silencer_table[, .(n_silencers = .N), by = silencer_class]

  tmp <- data.table::copy(nuc_aligned)
  tmp[, bin := as.integer(round(rel_pos / bin_size) * bin_size)]

  meta <- tmp[, .(n_midpoints = .N), by = .(silencer_class, bin)]
  meta <- merge(full_grid, meta, by = c("silencer_class", "bin"), all.x = TRUE, sort = TRUE)
  meta[is.na(n_midpoints), n_midpoints := 0L]
  meta <- merge(meta, n_per_class, by = "silencer_class", all.x = TRUE)
  meta[, nucleosome_midpoints_per_silencer := data.table::fifelse(n_silencers > 0, n_midpoints / n_silencers, NA_real_)]
  meta[, c("n_midpoints", "n_silencers") := NULL]
  meta[]
}

# ── Check existing outputs ─────────────────────────────────────────────────────
out_silencer_table <- file.path(output_dir, sprintf("%s_%s_%s_silencer_table.tsv.gz", default_sample, default_chrom, class_tag))
out_nuc_aligned    <- file.path(output_dir, sprintf("%s_%s_%s_nucleosome_silencer_aligned.tsv.gz", default_sample, default_chrom, class_tag))
out_nuc_summary    <- file.path(output_dir, sprintf("%s_%s_%s_silencer_structure.tsv.gz", default_sample, default_chrom, class_tag))
out_nuc_meta       <- file.path(output_dir, sprintf("%s_%s_%s_nucleosome_meta.tsv.gz", default_sample, default_chrom, class_tag))

if (all(file.exists(c(out_silencer_table, out_nuc_aligned, out_nuc_summary, out_nuc_meta)))) {
  message("Part 1 outputs already exist; skipping nucleosome stage.")
} else {
  assert_file_exists(bam_path, "BAM")
  assert_file_exists(paste0(bam_path, ".bai"), "BAM index")
  assert_file_exists(nuc_path, "nucleosome BED12")
  for (p in silencer_files[selected_classes]) {
    assert_file_exists(p, "silencer BED")
  }

  message("Reading BAM seqinfo ...")
  bam_seqinfo <- read_bam_seqinfo(bam_path)

  message("Reading silencer annotations ...")
  silencer_objects <- read_silencer_table(
    silencer_files,
    selected_classes = selected_classes,
    chrom = default_chrom,
    seqinfo = bam_seqinfo,
    window_bp = window_bp
  )
  silencer_table <- silencer_objects$silencer_table
  silencer_wins  <- silencer_objects$silencer_wins

  write_gz_tsv(silencer_table, out_silencer_table)

  message("Reading nucleosome BED12 ...")
  nuc_bed12 <- read_fiber_bed12(nuc_path)

  message("Expanding nucleosome feature blocks ...")
  nuc_features <- expand_ft_bed12_features(
    bed12 = nuc_bed12,
    feature_name = "nucleosome",
    drop_terminal_blocks = TRUE
  )
  nuc_gr <- feature_blocks_to_granges(nuc_features)

  message("Aligning nucleosomes to silencer windows ...")
  nuc_silencer_aligned <- align_feature_blocks_to_silencer(
    feature_dt = nuc_features,
    feature_gr = nuc_gr,
    silencer_wins = silencer_wins,
    silencer_table = silencer_table
  )

  nuc_summary <- summarize_nucleosomes_by_silencer(nuc_silencer_aligned, window_bp = window_bp)
  nuc_meta <- make_nucleosome_meta_by_class(
    nuc_aligned = nuc_silencer_aligned,
    silencer_table = silencer_table,
    window_bp = window_bp,
    bin_size = bin_size
  )

  write_gz_tsv(nuc_silencer_aligned, out_nuc_aligned)
  write_gz_tsv(nuc_summary, out_nuc_summary)
  write_gz_tsv(nuc_meta, out_nuc_meta)

  message("Part 1 complete.")
}