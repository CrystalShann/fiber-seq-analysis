#!/usr/bin/env Rscript
# Pool raw 1-bp m6A patterns within each enhancer class, then Manhattan + Leiden.
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
  table_dir <- file.path(project_root, "macrophage_project/enhancer/TF_co-occ/tables")
  enhancer_file <- Sys.getenv("ENHANCER_MANHATTAN_ENHANCERS",
    file.path(table_dir, "02_sampled_enhancer_classes.tsv"))
  class_id <- Sys.getenv("ENHANCER_MANHATTAN_CLASS")
  if (!class_id %in% c("active", "inactive"))
    stop("ENHANCER_MANHATTAN_CLASS must be active or inactive; submit the array script")
  result_file <- Sys.getenv("ENHANCER_MANHATTAN_OUTPUT",
    file.path(table_dir, paste0("enhancer_manhattan_", class_id, "_results.rds")))
  ft_result_dir <- Sys.getenv("ENHANCER_MANHATTAN_FT_ROOT",
    file.path(project_root, "macrophage_project/FiberHMM/extract/ft_result_dir"))
  reference <- Sys.getenv("ENHANCER_MANHATTAN_REFERENCE", "/project/spott/reference/human/GRCh38/hg38.fa")
  python_path <- Sys.getenv("ENHANCER_MANHATTAN_PYTHON", "/project/spott/cshan/envs/Jupyter-notebook/bin/python")
  samples <- strsplit(Sys.getenv("ENHANCER_MANHATTAN_SAMPLES", "LPS_0,LPS_5,LPS_10,LPS_15"), ",", fixed = TRUE)[[1L]]
  neighbor_method <- match.arg(Sys.getenv("ENHANCER_MANHATTAN_NEIGHBORS", "approximate"),
                               c("approximate", "exact"))
  # Manhattan graph method follows the macrophage reference; requested k=100.
  k_neighbors <- as.integer(Sys.getenv("ENHANCER_MANHATTAN_K", "100"))
  stopifnot(!is.na(k_neighbors), k_neighbors >= 1L)
  resolution <- 1
  kernel_sigma <- NULL  # mean of directed retained neighbor distances
  seed <- 1L
  max_fibers <- as.integer(Sys.getenv("ENHANCER_MANHATTAN_MAX_FIBERS", "10000"))
  stopifnot(length(max_fibers) == 1L, !is.na(max_fibers), max_fibers >= 3L)
  n_jobs <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
  stopifnot(file.exists(reference), file.exists(paste0(reference, ".fai")),
    file.exists(python_path), !anyDuplicated(samples), length(samples) > 0L,
    all(nzchar(samples)), is.finite(n_jobs), n_jobs >= 1L)
  analysis_enhancers <- fread(enhancer_file)
  required_columns <- c("enhancer_id", "enhancer_class", "chr", "midpoint0")
  stopifnot(all(required_columns %in% names(analysis_enhancers)), nrow(analysis_enhancers) > 0L,
    !anyDuplicated(analysis_enhancers$enhancer_id),
    !anyNA(analysis_enhancers[, ..required_columns]),
    setequal(analysis_enhancers$enhancer_class, c("active", "inactive")))
  analysis_enhancers <- analysis_enhancers[enhancer_class == class_id]
  message("Array class: ", class_id, "; sampled enhancers: ", nrow(analysis_enhancers))
  source(file.path(project_root, "code/topic_model/topic_modelling_functions.r"), local = TRUE)
  source(file.path(project_root, "code/enhancer/enhancer_read_functions.R"), local = TRUE)
  source(file.path(project_root, "code/enhancer/enhancer_manhattan_plots.R"), local = TRUE)
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
    path = file.path(project_root, "code/enhancer"), convert = TRUE)
  # Check dependencies before spending time on extraction.
  if (neighbor_method == "approximate") invisible(reticulate::import("pynndescent", convert = FALSE))

  results <- list()
  for (class_id in unique(analysis_enhancers$enhancer_class)) {
    message("Starting ", class_id, " enhancer extraction at ", Sys.time())
    inputs <- collect_enhancer_read_inputs(analysis_enhancers[enhancer_class == class_id],
                                           samples, ft_result_dir, reference)
    n_eligible <- nrow(inputs$mat)
    set.seed(seed)
    selected <- if (n_eligible > max_fibers) sort(sample.int(n_eligible, max_fibers)) else seq_len(n_eligible)
    inputs$mat <- inputs$mat[selected, , drop = FALSE]
    inputs$metadata <- inputs$metadata[selected]
    invisible(gc())
    message("Fiber sampling: ", length(selected), " of ", n_eligible,
            " eligible fibers; seed = ", seed)
    result <- cluster_enhancer_manhattan(inputs$mat, inputs$metadata, neighbor_builder, leiden_partition,
      k_neighbors = k_neighbors, resolution = resolution, sigma = kernel_sigma,
      seed = seed, neighbor_method = neighbor_method, n_jobs = n_jobs)
    result$feat_mat <- inputs$mat
    result$footprints <- collect_enhancer_footprints(result$assignments,
      analysis_enhancers, dirname(ft_result_dir))
    result$extraction_qc <- inputs$qc
    result$fiber_sampling <- list(n_eligible = n_eligible, n_selected = length(selected),
      max_fibers = max_fibers, seed = seed, method = "uniform_without_replacement")
    results[[class_id]] <- result
    message("Finished ", class_id, " at ", Sys.time())
    rm(inputs, result)
    invisible(gc())
  }
  input_regions <- as.data.frame(analysis_enhancers[, .(
    enhancer_id = as.character(enhancer_id), enhancer_class = as.character(enhancer_class),
    chr = as.character(chr), midpoint0 = as.integer(midpoint0))])
  input_regions <- input_regions[order(input_regions$enhancer_id), , drop = FALSE]
  rownames(input_regions) <- NULL
  saved <- list(results = results,
    inputs = list(regions = input_regions, samples = samples,
      enhancer_file = normalizePath(enhancer_file), ft_result_dir = ft_result_dir, reference = reference),
    parameters = list(n_features = 1000L, position_bp = -500:499,
      signal = "raw_binary_m6a", neighbor_method = neighbor_method,
      max_fibers = max_fibers, sampling_seed = seed,
      k_neighbors = k_neighbors, resolution = resolution, sigma = kernel_sigma,
      seed = seed, n_jobs = n_jobs, igraph_version = as.character(packageVersion("igraph"))),
    job_id = Sys.getenv("SLURM_JOB_ID", "interactive"),
    array_job_id = Sys.getenv("SLURM_ARRAY_JOB_ID", "interactive"),
    array_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", ""),
    enhancer_class = class_id, completed_at = Sys.time())
  dir.create(dirname(result_file), recursive = TRUE, showWarnings = FALSE)
  partial <- tempfile(".enhancer_manhattan_results_", tmpdir = dirname(result_file), fileext = ".rds")
  on.exit(unlink(partial), add = TRUE)
  saveRDS(saved, partial)
  if (!file.rename(partial, result_file)) stop("Could not publish Manhattan results: ", result_file)
  message("Saved completed Manhattan results: ", result_file)
  plot_dir <- Sys.getenv("ENHANCER_MANHATTAN_PLOT_DIR",
    file.path(dirname(dirname(result_file)), "plots"))
  save_enhancer_manhattan_pdfs(results[[class_id]], class_id, plot_dir,
                              dirname(result_file), expected_k = k_neighbors)
  message("Saved Manhattan heatmap and footprint PDFs: ", plot_dir)
}

if (sys.nframe() == 0L) run_enhancer_manhattan()
