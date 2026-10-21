#!/usr/bin/env Rscript
# Rscript code/divergent_promoter/run_divergent_promoters.R [--option=value]
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1L])
script_dir <- dirname(normalizePath(script))
source(file.path(script_dir, "divergent_promoter_functions.R"))
setDTthreads(2L)
cfg <- list(input_dir = "/project/spott/cshan/annotations/fantom5_LCL",
            ctss_file = "fantom5.hg38.LCL.consensus.TSS.bed.gz",
            tpm_file = "hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz",
            output_dir = "/project/spott/cshan/fiber-seq/LCL_project/divergent_tss",
            gencode = "/project/spott/cshan/annotations/gencode.v49.primary_assembly.annotation.gtf",
            score_threshold = 1, tpm_threshold = 1, chunk_size = 500L,
            min_distance = 1L, max_distance = 500L)
for (arg in args) {
  if (!grepl("^--[a-z_-]+=", arg)) stop("Use --option=value; options: ", paste(names(cfg), collapse = ", "))
  name <- gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", arg))
  if (!name %in% names(cfg)) stop("Unknown option: ", name)
  value <- sub("^--[^=]+=", "", arg)
  cfg[[name]] <- if (is.numeric(cfg[[name]])) as.numeric(value) else value
}
stopifnot(is.finite(cfg$score_threshold), cfg$score_threshold > 0,
          is.finite(cfg$tpm_threshold), cfg$tpm_threshold > 0,
          is.finite(cfg$chunk_size), cfg$chunk_size >= 1, cfg$chunk_size %% 1 == 0,
          is.finite(cfg$min_distance), cfg$min_distance >= 1, cfg$min_distance %% 1 == 0,
          is.finite(cfg$max_distance), cfg$max_distance >= cfg$min_distance, cfg$max_distance %% 1 == 0)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
run_ctss_workflow(cfg)
