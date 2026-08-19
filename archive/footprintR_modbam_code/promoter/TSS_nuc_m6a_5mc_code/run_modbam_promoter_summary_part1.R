# ---------------------------------------------------------------------------
# Part 1: validation, TSS setup, nucleosome summaries
# ---------------------------------------------------------------------------
source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/modbam_footprintR_functions.R")

default_sample <- Sys.getenv("SAMPLE", unset = "AL10_bc2178_19130")
default_chrom <- Sys.getenv("CHROM", unset = "chr1")
window_bp <- as.integer(Sys.getenv("WINDOW_BP", unset = "10000"))
bin_size <- as.integer(Sys.getenv("BIN_SIZE", unset = "100"))

results_root_dir <- "/project/spott/cshan/fiber-seq/results/promoter_summary/10kb_TSS"
threshold_dir <- file.path(results_root_dir, "threshold_0.9")
ccre_output_dir <- file.path(results_root_dir, "ccre")

bam_path <- file.path(
  "/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples",
  sprintf("%s.5mC.6mA.aligned.phased.bam", default_sample)
)
tss_path <- "/project/spott/cshan/annotations/TSS.gencode.v49.bed"

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

m6a_ft_path <- file.path(
  sample_extract_dir,
  "m6a_by_chr",
  sprintf("%s.ft_extracted_m6a.%s.bed.gz", default_sample, default_chrom)
)

cpg_ft_path <- file.path(
  sample_extract_dir,
  "cpg_by_chr",
  sprintf("%s.ft_extracted_cpg.%s.bed.gz", default_sample, default_chrom)
)

fire_peak_path <- file.path(
  "/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results",
  default_sample,
  sprintf("%s-fire-v0.1-peaks.bed.gz", default_sample)
)

fire_element_path <- file.path(
  "/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results",
  default_sample,
  "additional-outputs-v0.1",
  "fire-peaks",
  sprintf("%s-v0.1-fire-elements.bed.gz", default_sample)
)

output_dir <- file.path(threshold_dir, default_sample)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ccre_output_dir, recursive = TRUE, showWarnings = FALSE)

if (all(file.exists(c(file.path(output_dir, sprintf("%s_%s_tss_table.tsv.gz", default_sample, default_chrom)), file.path(output_dir, sprintf("%s_%s_nucleosome_tss_aligned.tsv.gz", default_sample, default_chrom)), file.path(output_dir, sprintf("%s_%s_promoter_structure.tsv.gz", default_sample, default_chrom)), file.path(output_dir, sprintf("%s_%s_nucleosome_meta.tsv.gz", default_sample, default_chrom)))))) {
  message("Part 1 outputs already exist; skipping nucleosome stage.")
  bam_seqinfo <- read_bam_seqinfo(bam_path)
  tss_objects <- read_gencode_tss(tss_path, chrom = default_chrom, window_bp = window_bp, seqinfo = bam_seqinfo)
  tss_table <- tss_objects$tss_table
  tss_gr <- tss_objects$tss_gr
  promoter_windows <- tss_objects$promoter_windows
  promoter_meta <- tss_objects$promoter_meta
  nuc_summary <- data.table::fread(file.path(output_dir, sprintf("%s_%s_promoter_structure.tsv.gz", default_sample, default_chrom)))
} else {
  assert_file_exists(bam_path, "BAM")
  assert_file_exists(paste0(bam_path, ".bai"), "BAM index")
  assert_file_exists(tss_path, "TSS BED")
  assert_file_exists(nuc_path, "nucleosome BED12")

  message("m6A ft extract BED present: ", file.exists(m6a_ft_path))
  message("CpG ft extract BED present: ", file.exists(cpg_ft_path))
  message("FIRE peak BED present:      ", file.exists(fire_peak_path))
  message("FIRE element BED present:   ", file.exists(fire_element_path))

# ---------------------------------------------------------------------------
# Setup: seqinfo, TSS objects
# ---------------------------------------------------------------------------
  message("Reading BAM seqinfo ...")
  bam_seqinfo <- read_bam_seqinfo(bam_path)

  message("Reading TSS annotations ...")
  tss_objects <- read_gencode_tss(tss_path, chrom = default_chrom, window_bp = window_bp, seqinfo = bam_seqinfo)
  tss_table <- tss_objects$tss_table
  tss_gr <- tss_objects$tss_gr
  promoter_windows <- tss_objects$promoter_windows
  promoter_meta <- tss_objects$promoter_meta

  write_gz_tsv(
    tss_table,
    file.path(output_dir, sprintf("%s_%s_tss_table.tsv.gz", default_sample, default_chrom))
  )

  nuc_bed12 <- read_fiber_bed12(nuc_path)

  nuc_features <- expand_ft_bed12_features(
    bed12 = nuc_bed12,
    feature_name = "nucleosome",
    drop_terminal_blocks = TRUE
  )

  nuc_gr <- feature_blocks_to_granges(nuc_features)

  nuc_tss_aligned <- align_feature_blocks_to_tss(
    feature_dt = nuc_features,
    feature_gr = nuc_gr,
    promoter_windows = promoter_windows,
    promoter_meta = promoter_meta
  )

  nuc_summary <- summarize_nucleosomes_by_tss(
    nuc_aligned = nuc_tss_aligned,
    window_bp = window_bp
  )

  nuc_meta <- make_nucleosome_meta(
    nuc_aligned = nuc_tss_aligned,
    num_tss = nrow(tss_table),
    window_bp = window_bp,
    bin_size = bin_size
  )

  write_gz_tsv(
    nuc_tss_aligned,
    file.path(output_dir, sprintf("%s_%s_nucleosome_tss_aligned.tsv.gz", default_sample, default_chrom))
  )

  write_gz_tsv(
    nuc_summary,
    file.path(output_dir, sprintf("%s_%s_promoter_structure.tsv.gz", default_sample, default_chrom))
  )

  write_gz_tsv(
    nuc_meta,
    file.path(output_dir, sprintf("%s_%s_nucleosome_meta.tsv.gz", default_sample, default_chrom))
  )

  rm(nuc_bed12, nuc_features, nuc_gr, nuc_meta)
  gc()
}
