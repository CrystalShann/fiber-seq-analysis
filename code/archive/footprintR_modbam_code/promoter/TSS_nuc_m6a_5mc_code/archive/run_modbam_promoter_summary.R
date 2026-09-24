library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)
library(SummarizedExperiment)
library(BiocParallel)
library(Rsamtools)
library(footprintR)
library(ggplot2)
library(patchwork)
library(scales)



# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
default_sample <- "AL10_bc2178_19130"
default_chrom  <- "chr1"
window_bp          <- 5000L
bin_size           <- 10L
mod_prob_threshold <- 0.9
example_tss_index  <- 1L
results_root_dir   <- "/project/spott/cshan/fiber-seq/results/promoter_summary/5kb_TSS"
threshold_dir      <- file.path(results_root_dir, "threshold_0.9")
ccre_output_dir    <- file.path(results_root_dir, "ccre")

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

bamfiles <- setNames(bam_path, default_sample)
bpparam_modbam <- BiocParallel::SerialParam()

# ---------------------------------------------------------------------------
# Split execution stages
# ---------------------------------------------------------------------------
source(file.path("/project/spott/cshan/fiber-seq/footprintR_modbam_code", "modbam_footprintR_functions.R"))
source(file.path("/project/spott/cshan/fiber-seq/footprintR_modbam_code", "run_modbam_promoter_summary_part1.R"))
source(file.path("/project/spott/cshan/fiber-seq/footprintR_modbam_code", "run_modbam_promoter_summary_m6A.R"))
source(file.path("/project/spott/cshan/fiber-seq/footprintR_modbam_code", "run_modbam_promoter_summary_5mC.R"))
source(file.path("/project/spott/cshan/fiber-seq/footprintR_modbam_code", "run_modbam_promoter_summary_join.R"))
