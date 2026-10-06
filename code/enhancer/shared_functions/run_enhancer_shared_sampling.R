#!/usr/bin/env Rscript
# Prepare ONE capped fiber cohort; both pooled analyses consume this saved matrix.
cap_per_timepoint <- 3L
cap_per_enhancer <- 10L
n_fibers_per_class <- 10000L
seed <- 1L

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(Rsamtools)
  library(Matrix)
})

run_enhancer_shared_sampling <- function() {
  project_root <- Sys.getenv("ENHANCER_PROJECT_ROOT", "/project/spott/cshan/fiber-seq")
  table_dir <- file.path(project_root, "macrophage_project/enhancer/TF_co-occ/tables")
  source(file.path(project_root, "code/enhancer/shared_functions/enhancer_fiber_sampling.R"), local = TRUE)
  source(file.path(project_root, "code/enhancer/shared_functions/enhancer_shared_sampling.R"), local = TRUE)
  parameter <- function(name, default, minimum = 1L)
    validate_sampling_integer(as.numeric(Sys.getenv(name, as.character(default))), name, minimum)
  cap_per_timepoint <- parameter("ENHANCER_CAP_PER_TIMEPOINT", cap_per_timepoint)
  cap_per_enhancer <- parameter("ENHANCER_CAP_PER_ENHANCER", cap_per_enhancer)
  n_fibers_per_class <- parameter("ENHANCER_N_FIBERS_PER_CLASS", n_fibers_per_class)
  sampling_seed <- parameter("ENHANCER_SAMPLING_SEED", seed, 0L)
  enhancer_file <- Sys.getenv("ENHANCER_SHARED_ENHANCERS",
    file.path(table_dir, "02_sampled_enhancer_classes.tsv"))
  output <- Sys.getenv("ENHANCER_SHARED_FIBERS",
    file.path(table_dir, "enhancer_shared_fibers_pooled_capped.rds"))
  ft_result_dir <- Sys.getenv("ENHANCER_SHARED_FT_ROOT",
    file.path(project_root, "macrophage_project/FiberHMM/extract/ft_result_dir"))
  reference <- Sys.getenv("ENHANCER_SHARED_REFERENCE", "/project/spott/reference/human/GRCh38/hg38.fa")
  samples <- strsplit(Sys.getenv("ENHANCER_SHARED_SAMPLES", "LPS_0,LPS_5,LPS_10,LPS_15"), ",", fixed = TRUE)[[1L]]
  regions <- fread(enhancer_file)
  stopifnot(all(c("enhancer_id", "enhancer_class", "chr", "midpoint0") %in% names(regions)),
    !anyNA(regions[, .(enhancer_id, enhancer_class, chr, midpoint0)]),
    !anyDuplicated(regions$enhancer_id), setequal(regions$enhancer_class, c("active", "inactive")),
    !anyDuplicated(samples), all(nzchar(samples)),
    file.exists(reference), file.exists(paste0(reference, ".fai")))
  source(file.path(project_root, "code/topic_model/topic_modelling_functions.r"), local = TRUE)
  source(file.path(project_root, "code/enhancer/shared_functions/enhancer_read_functions.R"), local = TRUE)
  source(file.path(project_root, "code/enhancer/shared_functions/enhancer_footprint_functions.R"), local = TRUE)
  parts <- list()
  for (class_id in c("active", "inactive")) {
    message("Preparing shared capped sample for ", class_id, " at ", Sys.time())
    inputs <- collect_enhancer_read_inputs(regions[enhancer_class == class_id],
      samples, ft_result_dir, reference)
    selection <- sample_capped_enhancer_fibers(inputs$metadata,
      cap_per_timepoint, cap_per_enhancer, n_fibers_per_class, sampling_seed, samples, class_id)
    inputs$mat <- inputs$mat[selection$selected, , drop = FALSE]
    inputs$metadata <- inputs$metadata[selection$selected]
    print(selection$diagnostics)
    info <- as.data.table(selection$info)
    info[, enhancer_class := class_id]
    parts[[class_id]] <- list(mat = inputs$mat, metadata = inputs$metadata,
      footprints = collect_enhancer_footprints(inputs$metadata, regions, dirname(ft_result_dir)),
      qc = inputs$qc, fiber_sampling = info, sampling_diagnostics = selection$diagnostics)
    rm(inputs, selection)
    invisible(gc())
  }
  input_regions <- as.data.frame(regions[, .(enhancer_id = as.character(enhancer_id),
    enhancer_class = as.character(enhancer_class), chr = as.character(chr),
    midpoint0 = as.integer(midpoint0))])
  input_regions <- input_regions[order(input_regions$enhancer_id), , drop = FALSE]
  rownames(input_regions) <- NULL
  shared <- list(mat = do.call(rbind, lapply(parts, `[[`, "mat")),
    metadata = rbindlist(lapply(parts, `[[`, "metadata")),
    footprints = rbindlist(lapply(parts, `[[`, "footprints")),
    qc = rbindlist(lapply(parts, `[[`, "qc")),
    fiber_sampling = rbindlist(lapply(parts, `[[`, "fiber_sampling")),
    sampling_diagnostics = rbindlist(lapply(parts, `[[`, "sampling_diagnostics")),
    inputs = list(regions = input_regions, samples = samples,
      enhancer_file = normalizePath(enhancer_file), ft_result_dir = ft_result_dir, reference = reference),
    parameters = list(cap_per_timepoint = cap_per_timepoint, cap_per_enhancer = cap_per_enhancer,
      n_fibers_per_class = n_fibers_per_class, sampling_seed = sampling_seed,
      sampling_method = "capped_enhancer_timepoint_round_robin", position_bp = -500:499,
      signal = "raw_binary_m6a"),
    sample_id = paste0("pooled_capped_", Sys.getenv("SLURM_JOB_ID", "interactive"), "_",
      format(Sys.time(), "%Y%m%dT%H%M%OS6", tz = "UTC")),
    job_id = Sys.getenv("SLURM_JOB_ID", "interactive"), completed_at = Sys.time())
  validate_shared_enhancer_fibers(shared)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  write_table <- function(x, stem) fwrite(x,
    file.path(dirname(output), paste0(stem, "_pooled_capped.tsv")), sep = "\t")
  manifest <- copy(shared$metadata)
  manifest[, sample_id := shared$sample_id]
  write_table(manifest, "shared_fiber_manifest")
  write_table(shared$fiber_sampling, "shared_fiber_sampling")
  write_table(shared$sampling_diagnostics, "shared_sampling_diagnostics")
  write_table(shared$qc, "shared_extraction_qc")
  atomic_save_enhancer_rds(shared, output)
  message("Saved shared sample: ", nrow(shared$mat), " fibers (active and inactive pooled): ", output)
}

if (sys.nframe() == 0L) run_enhancer_shared_sampling()
