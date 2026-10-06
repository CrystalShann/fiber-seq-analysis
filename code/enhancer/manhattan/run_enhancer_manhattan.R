#!/usr/bin/env Rscript
# Pool active and inactive fibers from the shared capped sample for Manhattan + Leiden.
# Extraction and sampling are performed once by prepare_enhancer_shared_fibers.R.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(Rsamtools)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(ComplexHeatmap)
})

cluster_enhancer_manhattan <- function(binary, metadata, neighbor_builder, partitioner,
                                       k_neighbors = 100L, resolution = 1,
                                       sigma = NULL, seed = 1L,
                                       neighbor_method = "approximate", n_jobs = 4L) {
  stopifnot(ncol(binary) == 1000L,
    identical(colnames(binary), as.character(-500:499)),
    identical(as.character(rownames(binary)), as.character(metadata$RID)),
    !anyDuplicated(metadata$RID), !anyNA(binary), all(binary@x == 1))
  assignments <- data.table::copy(as.data.table(metadata))
  assignments[, `:=`(cluster = NA_character_, clustering_status = "insufficient_reads")]
  profiles <- matrix(numeric(), nrow = 0L, ncol = ncol(binary),
                     dimnames = list(character(), colnames(binary)))
  info <- list(n_reads = nrow(binary), n_features = ncol(binary),
    n_zero_call_reads = sum(Matrix::rowSums(binary) == 0),
    neighbor_method = neighbor_method, k_neighbors = k_neighbors,
    resolution = resolution, seed = seed, window_size = 1L,
    impute_missing = FALSE, n_clusters = 0L)
  if (nrow(binary) < 3L)
    return(list(assignments = assignments, m6a_profiles = profiles, info = info))

  message("Finding ", neighbor_method, " Manhattan neighbors for ", nrow(binary), " fibers")
  nn <- neighbor_builder$build_manhattan_knn(binary, k_neighbors = as.integer(k_neighbors),
    method = neighbor_method, seed = as.integer(seed), n_jobs = as.integer(n_jobs))
  k_eff <- as.integer(nn$k_eff)
  stopifnot(identical(as.integer(dim(nn$indices)), c(nrow(binary), k_eff)),
    identical(dim(nn$indices), dim(nn$distances)),
    all(is.finite(nn$distances)), all(nn$distances >= 0),
    all(nn$indices >= 0 & nn$indices < nrow(binary)))
  mean_knn_dist <- mean(nn$distances)
  if (is.null(sigma)) sigma <- mean_knn_dist
  if (!is.finite(sigma) || sigma <= 0) {
    warning("All retained Manhattan distances are zero; using sigma = 1")
    sigma <- 1
  }
  # Same exponential weighting and undirected UNION as manhattan_knn_graph().
  # Numeric vertices avoid storing each long enhancer/fiber ID on every edge.
  edges <- data.table(from = rep(seq_len(nrow(binary)), each = k_eff),
    to = as.integer(t(nn$indices)) + 1L, distance = as.numeric(t(nn$distances)))
  stopifnot(all(edges$from != edges$to))
  edges[, `:=`(a = pmin(from, to), b = pmax(from, to))]
  edges <- unique(edges[, .(a, b, distance)], by = c("a", "b"))
  graph <- igraph::graph_from_edgelist(as.matrix(edges[, .(a, b)]), directed = FALSE)
  graph <- igraph::set_edge_attr(graph, "weight", value = exp(-edges$distance / sigma))
  graph <- igraph::set_edge_attr(graph, "manhattan_distance", value = edges$distance)
  stopifnot(igraph::vcount(graph) == nrow(binary))
  # Original reference's modularity objective, weights and convergence settings.
  part <- partitioner(graph, resolution = resolution, n_iterations = -1L, seed = seed)
  membership <- as.integer(part$membership)
  stopifnot(length(membership) == nrow(binary), !anyNA(membership))
  cluster_order <- names(sort(table(membership), decreasing = TRUE))
  cluster_levels <- paste0("cluster", seq_along(cluster_order))
  assignments[, `:=`(cluster = paste0("cluster", match(as.character(membership), cluster_order)),
                     clustering_status = "clustered")]
  profiles <- matrix(NA_real_, length(cluster_levels), ncol(binary),
                     dimnames = list(cluster_levels, colnames(binary)))
  for (i in seq_along(cluster_levels))
    profiles[i, ] <- Matrix::colMeans(binary[assignments$cluster == cluster_levels[i], , drop = FALSE])
  info$n_clusters <- length(cluster_levels)
  info$k_eff <- k_eff
  info$sigma <- sigma
  info$mean_knn_dist <- mean_knn_dist
  info$n_edges <- igraph::ecount(graph)
  info$quality <- part$quality
  info$neighbor_details <- nn[setdiff(names(nn), c("indices", "distances"))]
  message("Leiden: ", info$n_clusters, " clusters, ", info$n_edges, " edges; sigma = ", signif(sigma, 5))
  # Store assignments and 1-bp means; no dense pairwise distances or graph in the RDS.
  list(assignments = assignments, m6a_profiles = profiles, info = info)
}

run_enhancer_manhattan <- function() {
  project_root <- Sys.getenv("ENHANCER_PROJECT_ROOT", "/project/spott/cshan/fiber-seq")
  source(file.path(project_root, "code/enhancer/shared_functions/enhancer_fiber_sampling.R"), local = TRUE)
  source(file.path(project_root, "code/enhancer/shared_functions/enhancer_shared_sampling.R"), local = TRUE)
  table_dir <- file.path(project_root, "macrophage_project/enhancer/TF_co-occ/tables")
  shared_path <- Sys.getenv("ENHANCER_SHARED_FIBERS",
    file.path(table_dir, "enhancer_shared_fibers_pooled_capped.rds"))
  shared <- load_shared_enhancer_fibers(shared_path)
  result_file <- Sys.getenv("ENHANCER_MANHATTAN_OUTPUT",
    file.path(table_dir, "enhancer_manhattan_results_pooled_capped.rds"))
  dir.create(dirname(result_file), recursive = TRUE, showWarnings = FALSE)
  diagnostic_file <- function(kind) file.path(dirname(result_file),
    paste0("manhattan_", kind, "_pooled_capped.tsv"))
  python_path <- Sys.getenv("ENHANCER_MANHATTAN_PYTHON", "/project/spott/cshan/envs/Jupyter-notebook/bin/python")
  samples <- as.character(shared$inputs$samples)
  neighbor_method <- match.arg(Sys.getenv("ENHANCER_MANHATTAN_NEIGHBORS", "approximate"),
                               c("approximate", "exact"))
  # Manhattan graph method follows the macrophage reference; requested k=100.
  k_neighbors <- as.integer(Sys.getenv("ENHANCER_MANHATTAN_K", "100"))
  stopifnot(!is.na(k_neighbors), k_neighbors >= 1L)
  resolution <- 1
  kernel_sigma <- NULL  # mean of directed retained neighbor distances
  leiden_seed <- 1L
  n_jobs <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
  stopifnot(file.exists(python_path), is.finite(n_jobs), n_jobs >= 1L)
  source(file.path(project_root, "code/enhancer/manhattan/enhancer_manhattan_plots.R"), local = TRUE)
  source(file.path(project_root, "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r"), local = TRUE)
  source(file.path(project_root, "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r"), local = TRUE)
  if (!requireNamespace("igraph", quietly = TRUE) || !requireNamespace("reticulate", quietly = TRUE))
    stop("R packages igraph and reticulate are required")
  Sys.setenv(PYTHONDONTWRITEBYTECODE = "1", NUMBA_CACHE_DIR = file.path(tempdir(), "enhancer-manhattan-numba"))
  dir.create(Sys.getenv("NUMBA_CACHE_DIR"), recursive = TRUE, showWarnings = FALSE)
  reticulate::use_python(python_path, required = TRUE)
  py_sys <- reticulate::import("sys", convert = FALSE)
  py_sys$dont_write_bytecode <- TRUE
  neighbor_builder <- reticulate::import_from_path("enhancer_manhattan_neighbors",
    path = file.path(project_root, "code/enhancer/manhattan"), convert = TRUE)
  if (neighbor_method == "approximate") invisible(reticulate::import("pynndescent", convert = FALSE))

  message("Clustering the shared active/inactive sample: ", nrow(shared$metadata), " fibers")
  message("Shared sample: ", shared$sample_id, "; MD5: ", shared$md5)
  print(shared$fiber_sampling)
  fwrite(shared$sampling_diagnostics, diagnostic_file("sampling_diagnostics"), sep = "\t")
  fwrite(shared$fiber_sampling, diagnostic_file("fiber_sampling"), sep = "\t")
  result <- cluster_enhancer_manhattan(shared$mat, shared$metadata, neighbor_builder, leiden_partition,
    k_neighbors = k_neighbors, resolution = resolution, sigma = kernel_sigma,
    seed = leiden_seed, neighbor_method = neighbor_method, n_jobs = n_jobs)
  result$sampling_diagnostics <- shared$sampling_diagnostics
  result$cluster_diagnostics <- enhancer_cluster_diagnostics(result$assignments, "pooled", samples)
  message("Cluster contributions (>10% from one enhancer or <50 enhancers are flagged):")
  print(result$cluster_diagnostics)
  fwrite(result$cluster_diagnostics, diagnostic_file("cluster_diagnostics"), sep = "\t")
  if (any(result$cluster_diagnostics$flagged))
    warning("Pooled Manhattan: flagged clusters: ",
      paste(result$cluster_diagnostics[flagged == TRUE, cluster], collapse = ", "), call. = FALSE)
  result$feat_mat <- shared$mat
  result$footprints <- shared$footprints
  result$extraction_qc <- shared$qc
  result$fiber_sampling <- shared$fiber_sampling
  results <- list(pooled = result)
  saved <- list(results = results, inputs = shared$inputs,
    shared_sample = list(path = normalizePath(shared_path), md5 = shared$md5, sample_id = shared$sample_id),
    parameters = c(shared$parameters,
      list(n_features = 1000L, neighbor_method = neighbor_method,
        k_neighbors = k_neighbors, resolution = resolution, sigma = kernel_sigma,
        seed = leiden_seed, n_jobs = n_jobs, igraph_version = as.character(packageVersion("igraph")))),
    job_id = Sys.getenv("SLURM_JOB_ID", "interactive"),
    enhancer_class = "pooled", completed_at = Sys.time())
  partial <- tempfile(".enhancer_manhattan_results_", tmpdir = dirname(result_file), fileext = ".rds")
  on.exit(unlink(partial), add = TRUE)
  saveRDS(saved, partial)
  if (!file.rename(partial, result_file)) stop("Could not publish Manhattan results: ", result_file)
  message("Saved completed pooled Manhattan results: ", result_file)
  plot_dir <- Sys.getenv("ENHANCER_MANHATTAN_PLOT_DIR",
    file.path(dirname(dirname(result_file)), "plots"))
  save_enhancer_manhattan_pdfs(result, "pooled", plot_dir,
                              dirname(result_file), expected_k = k_neighbors, suffix = "_capped")
  message("Saved pooled Manhattan heatmap and footprint PDFs: ", plot_dir)
}

if (sys.nframe() == 0L) run_enhancer_manhattan()
