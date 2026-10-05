#!/usr/bin/env Rscript
# Run through run_enhancer_acf.sh; preserve the notebook's extraction and ACF methods.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(Rsamtools)
  library(Matrix)
})

run_enhancer_acf <- function() {
  project_root <- Sys.getenv("ENHANCER_PROJECT_ROOT", "/project/spott/cshan/fiber-seq")
  table_dir <- file.path(project_root, "macrophage_project/enhancer/TF_co-occ/tables")
  enhancer_file <- Sys.getenv("ENHANCER_ACF_ENHANCERS",
                             file.path(table_dir, "02_sampled_enhancer_classes.tsv"))
  result_file <- Sys.getenv("ENHANCER_ACF_OUTPUT",
                           file.path(table_dir, "enhancer_acf_results.rds"))
  ft_result_dir <- Sys.getenv("ENHANCER_ACF_FT_ROOT",
    file.path(project_root, "macrophage_project/FiberHMM/extract/ft_result_dir"))
  sample_table <- data.table(sample = strsplit(Sys.getenv("ENHANCER_ACF_SAMPLES",
    "LPS_0,LPS_5,LPS_10,LPS_15"), ",", fixed = TRUE)[[1L]])
  analysis_enhancers <- fread(enhancer_file)
  required_columns <- c("enhancer_id", "enhancer_class", "chr", "midpoint0")
  stopifnot(all(required_columns %in% names(analysis_enhancers)),
            nrow(analysis_enhancers) > 0L,
            !anyDuplicated(analysis_enhancers$enhancer_id),
            !anyNA(analysis_enhancers[, ..required_columns]),
            setequal(analysis_enhancers$enhancer_class, c("active", "inactive")),
            !anyDuplicated(sample_table$sample), all(nzchar(sample_table$sample)))
  source(file.path(project_root, "code/topic_model/topic_modelling_functions.r"), local = TRUE)
  source(file.path(project_root, "code/enhancer/enhancer_read_functions.R"), local = TRUE)
  acf_reference <- Sys.getenv("ENHANCER_ACF_REFERENCE", "/project/spott/reference/human/GRCh38/hg38.fa") # chromosome lengths only
  acf_python_path <- Sys.getenv("ENHANCER_ACF_PYTHON", "/project/spott/cshan/envs/Jupyter-notebook/bin/python")
  acf_source_dir <- file.path(project_root, "code/clustering_methods/auto_correlation")
  acf_n_features <- 1000L  # every lag 0:999 for the 1000-bp window
  acf_n_pcs <- 50L
  acf_n_neighbors <- 10L
  acf_resolution <- 0.4
  acf_seed <- 0L
  acf_batch_size <- 512L
  # Same pooled-class fiber sample as run_enhancer_manhattan.R: identical
  # extraction order, max_fibers and seed select the identical fibers.
  acf_max_fibers <- 10000L
  acf_sampling_seed <- 1L
  stopifnot(file.exists(acf_reference), file.exists(paste0(acf_reference, ".fai")),
            file.exists(acf_python_path), acf_n_features == 1000L)

  if (!requireNamespace("reticulate", quietly = TRUE)) stop("R package reticulate is required")
  Sys.setenv(PYTHONDONTWRITEBYTECODE = "1",
             NUMBA_CACHE_DIR = file.path(tempdir(), "enhancer-acf-numba"),
             MPLCONFIGDIR = file.path(tempdir(), "enhancer-acf-matplotlib"))
  dir.create(Sys.getenv("NUMBA_CACHE_DIR"), recursive = TRUE, showWarnings = FALSE)
  dir.create(Sys.getenv("MPLCONFIGDIR"), recursive = TRUE, showWarnings = FALSE)
  reticulate::use_python(acf_python_path, required = TRUE)
  acf_sys <- reticulate::import("sys", convert = FALSE)
  acf_sys$dont_write_bytecode <- TRUE
  # Check clustering dependencies before extracting the full enhancer universe.
  invisible(reticulate::import("scanpy", convert = FALSE))
  invisible(reticulate::import("leidenalg", convert = FALSE))
  acf_compute <- reticulate::import_from_path("03_compute_autocorrelations",
                                              path = acf_source_dir, convert = TRUE)
  acf_cluster <- reticulate::import_from_path("04_cluster_autocorrelations",
                                              path = acf_source_dir, convert = TRUE)

  cluster_enhancer_acf <- function(binary, metadata, n_features = 1000L, n_pcs = 50L,
                                   n_neighbors = 10L, resolution = 0.4, seed = 0L,
                                   batch_size = 512L) {
    stopifnot(identical(as.character(rownames(binary)), as.character(metadata$RID)),
              ncol(binary) == 1000L, n_features == ncol(binary), batch_size >= 1L)
    profiles <- matrix(NA_real_, nrow(binary), n_features)
    valid <- logical(nrow(binary))
    if (nrow(binary)) {
      for (first in seq.int(1L, nrow(binary), by = batch_size)) {
        rows <- seq.int(first, min(nrow(binary), first + batch_size - 1L))
        # Densify only this batch for the original per-fiber ACF function.
        batch <- acf_compute$autocorrelations(as.matrix(binary[rows, , drop = FALSE]),
                                              n_features = as.integer(n_features))
        profiles[rows, ] <- batch[[1L]]
        valid[rows] <- as.logical(batch[[2L]])
      }
    }
    stopifnot(all(is.finite(profiles[valid, , drop = FALSE])),
              all(abs(profiles[valid, 1L] - 1) < 1e-10))
    message("Clustering ", sum(valid), " valid fibers out of ", length(valid))
    fit <- acf_cluster$cluster_profiles(
      profiles, as.array(valid), n_pcs = as.integer(n_pcs),
      n_neighbors = as.integer(n_neighbors), resolution = resolution, seed = as.integer(seed))
    assignments <- data.table::copy(metadata)
    assignments[, `:=`(acf_status = as.character(fit[[2L]]),
                        cluster = NA_character_, UMAP1 = fit[[3L]][, 1L],
                        UMAP2 = fit[[3L]][, 2L])]
    clustered <- which(assignments$acf_status == "clustered")
    assignments$cluster[clustered] <- paste0("cluster", as.integer(fit[[1L]][clustered]) + 1L)
    cluster_names <- unique(assignments$cluster[clustered])
    if (length(cluster_names)) {
      cluster_names <- cluster_names[order(as.integer(sub("cluster", "", cluster_names)))]
    }
    mean_acf <- matrix(NA_real_, length(cluster_names), n_features,
                       dimnames = list(cluster_names, as.character(seq_len(n_features) - 1L)))
    mean_m6a <- matrix(NA_real_, length(cluster_names), ncol(binary),
                       dimnames = list(cluster_names, colnames(binary)))
    for (i in seq_along(cluster_names)) {
      rows <- which(assignments$cluster == cluster_names[i])
      mean_acf[i, ] <- colMeans(profiles[rows, , drop = FALSE])
      mean_m6a[i, ] <- Matrix::colMeans(binary[rows, , drop = FALSE])
    }
    # Keep per-fiber metadata in memory for UMAP/composition; export compact
    # summaries only. Dense ACF matrices and graph edges are not saved.
    list(assignments = assignments, acf_profiles = mean_acf,
         m6a_profiles = mean_m6a, info = fit[[4L]])
  }

  enhancer_acf_results <- list()
  for (class_id in unique(analysis_enhancers$enhancer_class)) {
    message("Starting ", class_id, " enhancer extraction at ", Sys.time())
    inputs <- collect_enhancer_read_inputs(
      analysis_enhancers[enhancer_class == class_id], as.character(sample_table$sample),
      ft_result_dir, reference = acf_reference)
    n_eligible <- nrow(inputs$mat)
    set.seed(acf_sampling_seed)
    selected <- if (n_eligible > acf_max_fibers) sort(sample.int(n_eligible, acf_max_fibers)) else seq_len(n_eligible)
    inputs$mat <- inputs$mat[selected, , drop = FALSE]
    inputs$metadata <- inputs$metadata[selected]
    invisible(gc())
    message("Fiber sampling: ", length(selected), " of ", n_eligible,
            " eligible fibers; seed = ", acf_sampling_seed)
    message("Computing ACF for ", nrow(inputs$mat), " fibers in ", class_id)
    result <- cluster_enhancer_acf(
      inputs$mat, inputs$metadata, n_features = acf_n_features, n_pcs = acf_n_pcs,
      n_neighbors = acf_n_neighbors, resolution = acf_resolution, seed = acf_seed,
      batch_size = acf_batch_size)
    result$extraction_qc <- inputs$qc
    result$fiber_sampling <- list(n_eligible = n_eligible, n_selected = length(selected),
      max_fibers = acf_max_fibers, seed = acf_sampling_seed, method = "uniform_without_replacement")
    enhancer_acf_results[[class_id]] <- result
    message("Finished ", class_id, " at ", Sys.time())
    rm(inputs, result)
    invisible(gc())
  }

  input_regions <- as.data.frame(analysis_enhancers[, .(
    enhancer_id = as.character(enhancer_id), enhancer_class = as.character(enhancer_class),
    chr = as.character(chr), midpoint0 = as.integer(midpoint0))])
  input_regions <- input_regions[order(input_regions$enhancer_id), , drop = FALSE]
  rownames(input_regions) <- NULL
  saved <- list(
    results = enhancer_acf_results,
    inputs = list(regions = input_regions, samples = as.character(sample_table$sample),
                  enhancer_file = normalizePath(enhancer_file), ft_result_dir = ft_result_dir,
                  reference = acf_reference),
    parameters = list(n_features = acf_n_features, n_pcs = acf_n_pcs,
      n_neighbors = acf_n_neighbors, resolution = acf_resolution,
      seed = acf_seed, batch_size = acf_batch_size,
      max_fibers = acf_max_fibers, sampling_seed = acf_sampling_seed),
    job_id = Sys.getenv("SLURM_JOB_ID", "interactive"), completed_at = Sys.time())
  dir.create(dirname(result_file), recursive = TRUE, showWarnings = FALSE)
  partial <- tempfile(".enhancer_acf_results_", tmpdir = dirname(result_file), fileext = ".rds")
  on.exit(unlink(partial), add = TRUE)
  saveRDS(saved, partial)
  if (!file.rename(partial, result_file)) stop("Could not publish ACF results: ", result_file)
  message("Saved completed ACF results: ", result_file)
}

run_enhancer_acf()
