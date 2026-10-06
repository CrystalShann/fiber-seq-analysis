#!/usr/bin/env Rscript
# Pool both enhancer classes from the shared capped cohort; preserve ACF methods.
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
  result_file <- Sys.getenv("ENHANCER_ACF_OUTPUT",
                           file.path(table_dir, "enhancer_acf_results_pooled_capped.rds"))
  shared_file <- Sys.getenv("ENHANCER_SHARED_FIBERS",
    file.path(table_dir, "enhancer_shared_fibers_pooled_capped.rds"))
  source(file.path(project_root, "code/enhancer/shared_functions/enhancer_shared_sampling.R"), local = TRUE)
  source(file.path(project_root, "code/enhancer/shared_functions/enhancer_fiber_sampling.R"), local = TRUE)
  source(file.path(project_root, "code/enhancer/acf/enhancer_acf_plots.R"), local = TRUE)
  shared <- load_shared_enhancer_fibers(shared_file)
  stopifnot(identical(as.character(rownames(shared$mat)), as.character(shared$metadata$RID)),
            !anyDuplicated(shared$metadata$RID), ncol(shared$mat) == 1000L,
            all(shared$metadata$enhancer_class %in% c("active", "inactive")))
  acf_python_path <- Sys.getenv("ENHANCER_ACF_PYTHON", "/project/spott/cshan/envs/Jupyter-notebook/bin/python")
  acf_source_dir <- file.path(project_root, "code/clustering_methods/auto_correlation")
  acf_n_features <- 1000L  # every lag 0:999 for the 1000-bp window
  acf_n_pcs <- 50L
  acf_n_neighbors <- 10L
  acf_resolution <- 0.4
  acf_seed <- 0L
  acf_batch_size <- 512L
  stopifnot(file.exists(acf_python_path), acf_n_features == 1000L)

  if (!requireNamespace("reticulate", quietly = TRUE)) stop("R package reticulate is required")
  Sys.setenv(PYTHONDONTWRITEBYTECODE = "1",
             NUMBA_CACHE_DIR = file.path(tempdir(), "enhancer-acf-numba"),
             MPLCONFIGDIR = file.path(tempdir(), "enhancer-acf-matplotlib"))
  dir.create(Sys.getenv("NUMBA_CACHE_DIR"), recursive = TRUE, showWarnings = FALSE)
  dir.create(Sys.getenv("MPLCONFIGDIR"), recursive = TRUE, showWarnings = FALSE)
  reticulate::use_python(acf_python_path, required = TRUE)
  acf_sys <- reticulate::import("sys", convert = FALSE)
  acf_sys$dont_write_bytecode <- TRUE
  # Both methods consume this saved cohort; neither extracts or samples again.
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

  message("Computing pooled ACF for ", nrow(shared$mat),
          " shared fibers from active and inactive enhancers at ", Sys.time())
  result <- cluster_enhancer_acf(
    shared$mat, shared$metadata, n_features = acf_n_features, n_pcs = acf_n_pcs,
    n_neighbors = acf_n_neighbors, resolution = acf_resolution, seed = acf_seed,
    batch_size = acf_batch_size)
  stopifnot(identical(as.character(result$assignments$RID), as.character(shared$metadata$RID)))
  result$extraction_qc <- shared$qc
  result$fiber_sampling <- shared$fiber_sampling
  result$sampling_diagnostics <- shared$sampling_diagnostics
  result$cluster_diagnostics <- enhancer_cluster_diagnostics(
    result$assignments, "pooled", shared$inputs$samples)
  message("Pooled ACF cluster contributions (constant fibers remain unclustered):")
  print(result$cluster_diagnostics)
  if (any(result$cluster_diagnostics$flagged))
    warning("Pooled ACF contribution flags: ",
      paste(result$cluster_diagnostics[flagged == TRUE, cluster], collapse = ", "), call. = FALSE)
  saved <- list(
    results = list(pooled = result), inputs = shared$inputs,
    shared_sample = list(path = normalizePath(shared_file), md5 = shared$md5,
                         sample_id = shared$sample_id),
    parameters = utils::modifyList(shared$parameters, list(
      pooling = "active_inactive_combined", n_features = acf_n_features, n_pcs = acf_n_pcs,
      n_neighbors = acf_n_neighbors, resolution = acf_resolution,
      seed = acf_seed, batch_size = acf_batch_size)),
    job_id = Sys.getenv("SLURM_JOB_ID", "interactive"), completed_at = Sys.time())
  dir.create(dirname(result_file), recursive = TRUE, showWarnings = FALSE)
  partial <- tempfile(".enhancer_acf_results_", tmpdir = dirname(result_file), fileext = ".rds")
  on.exit(unlink(partial), add = TRUE)
  saveRDS(saved, partial)
  if (!file.rename(partial, result_file)) stop("Could not publish ACF results: ", result_file)
  message("Saved completed ACF results: ", result_file)
  plot_dir <- Sys.getenv("ENHANCER_ACF_PLOT_DIR", file.path(dirname(dirname(result_file)), "plots"))
  save_enhancer_acf_outputs(result, plot_dir, dirname(result_file),
    samples = shared$inputs$samples, suffix = "_pooled_capped")
  message("Saved pooled ACF tables and PDFs: ", plot_dir)
}

if (sys.nframe() == 0L) run_enhancer_acf()
