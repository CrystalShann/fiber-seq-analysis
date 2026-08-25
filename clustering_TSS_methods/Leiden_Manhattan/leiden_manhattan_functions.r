# leiden_manhattan_functions.r
#
# Single-molecule Leiden clustering of Fiber-seq m6A reads, replicating the
# nano-NOMe-seq read-clustering procedure of Raviram, Jiang, Schippke, Cova &
# Skok (2026), "Cohesin collisions maintain ordered nucleosome architecture at
# boundaries and promoters" (bioRxiv 2026.05.22.727261), Methods / "Clustering
# reads", with their GpC (GCH) accessibility calls replaced by Fiber-seq m6A
# calls (1 = m6A / accessible, 0 = covered but unmodified, NA = not covered).
#
# The paper's procedure, step by step:
#   1. bin the per-read accessibility signal into 50-bp windows across the 2-kb
#      region; a read's bin value is the mean methylation call in the bin;
#   2. impute the bins still missing after binning by KNN imputation
#      (k = 5, scikit-learn);
#   3. balance conditions - sample an equal number of reads per condition -
#      and pool the sampled reads before clustering;
#   4. read-read similarity = Manhattan distance over the bins,
#      d(i,j) = sum_b |x_ib - x_jb|;
#   5. KNN graph over reads (k = 50 nearest neighbours, scikit-learn);
#   6. edge weights = affinities from an exponential kernel on those distances;
#   7. Leiden community detection (leidenalg, RBConfigurationVertexPartition)
#      at a given resolution;
#   8. cluster profiles = mean accessibility over the reads of each cluster.
#
# Deviations from the paper (all switched in the notebook's Settings chunk):
#   * Data type: Fiber-seq m6A instead of nano-NOMe GpC methylation.
#   * Reads: the full-span read set (filter_met_mat(), the topic model's row
#     filter, also used by the hamming and spearman clusterings), so all four
#     clusterings share one read universe. All four LPS timepoints are pooled
#     and clustered jointly, as in the topic model.
#   * No condition balancing (step 3): every full-span read enters the
#     clustering, and each read keeps its timepoint tag (sample_name) alongside
#     its cluster, so a read's cluster and timepoint are both recoverable from
#     the assignment table. The paper subsamples to equalise conditions; here
#     the timepoints are deliberately left at their observed depths.
#   * Window-size sweep instead of the paper's fixed 50 bp: window size in bp,
#     or 0 for no windowing at all (each m6A site is its own feature, the
#     default). At window size 0 the read x feature matrix is the 0/1 site
#     matrix itself.
#   * KNN imputation (step 2) is optional and off by default. With full-span
#     reads it is a no-op: every read covers every m6A site, so the only
#     missing bins are windows containing no m6A site at all, and those are
#     dropped as features (they carry no information for any read) rather than
#     imputed.
#   * No meta-clustering: the paper groups Leiden clusters into meta-clusters
#     by FFT + hierarchical clustering for display; the Leiden clusters are
#     used directly here.
#
# Implementation notes (details the paper's Methods leave open):
#   * Exponential kernel (step 6): affinity = exp(-d / sigma). sigma defaults
#     to the mean of the retained KNN distances, which keeps the weight
#     contrast comparable across window sizes (both d and the number of
#     features change with the window size). Passing sigma = ncol(feature
#     matrix) instead reproduces scikit-learn's laplacian_kernel default
#     (gamma = 1 / n_features), the most literal reading of "exponential
#     kernel" applied to Manhattan distances.
#   * The KNN graph is directed (read i's k nearest neighbours), so it is
#     symmetrised by the union of the two directions before clustering.
#     Manhattan distance is symmetric, so both directions of a kept pair carry
#     the same affinity and the union is just the set of unique unordered
#     pairs.
#   * igraph::cluster_leiden(objective_function = "modularity", resolution = g)
#     is the RBConfigurationVertexPartition of leidenalg: modularity with a
#     resolution parameter and vertex weights fixed to the (weighted) degree.
#     Vertex weights are therefore not passed explicitly (igraph warns that
#     doing so contradicts the modularity objective). n_iterations = -1 runs
#     Leiden to convergence.
#   * Leiden is randomised; every run is seeded (see the `seed` argument).
#   * Clusters are renumbered largest-first as cluster1..clusterN, so cluster
#     numbering is deterministic. The paper's cluster order comes from its
#     meta-clustering step, which is not used here.
#
# Input matrices come from get_sparse_met_mat() in
# code/topic_model/topic_modelling_functions.r (reads x m6A-site positions);
# source that file before this one.
#
# igraph needs libglpk at load time; if library(igraph) fails with
# "libglpk.so.40: cannot open shared object file", add glpk's lib directory
# (on this cluster /software/glpk-5.0-el8-x86_64/lib) to LD_LIBRARY_PATH
# before starting R.

suppressMessages({
  requireNamespace("igraph")
  requireNamespace("Matrix")
})


# ---------------------------------------------------------------------------
# Pooled m6A read x position matrix for one region, tagging each read with its
# sample of origin (identical to the hamming and spearman notebooks, so the
# four clusterings start from the same reads).
# ---------------------------------------------------------------------------
assemble_region_m6a <- function(sample_names, region_chr, region_start, region_end,
                                ft_result_dir, verbose = TRUE) {
  region_gr <- GRanges(region_chr, IRanges(region_start, region_end))
  reads_list <- lapply(sample_names, function(sample_name) {
    extracted_file <- file.path(
      ft_result_dir, sample_name, "extracted_results", "m6a_by_chr",
      paste0(sample_name, ".ft_extracted_m6a.", region_chr, ".bed.gz"))
    sample_reads <- extract_ft_region_reads(extracted_file, region_gr,
                                            keep_pos_in_region_only = TRUE,
                                            verbose = verbose)
    if (nrow(sample_reads) == 0) return(NULL)
    dplyr::mutate(sample_reads, sample_name = sample_name, .before = 1)
  })
  reads <- dplyr::bind_rows(reads_list)
  if (nrow(reads) == 0) stop("no reads in region")
  rids_df <- extract_ft_read_info(reads, keep_columns = c("sample_name", "score"))
  met_mat <- get_sparse_met_mat(reads, rids_df,
                                window_start = region_start,
                                window_end = region_end, base = "A")
  list(reads = reads, rids_df = rids_df, met_mat = met_mat)
}


# ---------------------------------------------------------------------------
# Step 1: read x feature matrix.
#
# window_size > 0: consecutive windows tiling [region_start, region_end] (the
# last window keeps the remainder); a read's feature value is the mean of its
# m6A site calls in the window. Windows containing no m6A site carry no
# information for any read and are dropped, and are recorded in the returned
# window annotation with n_sites = 0.
# window_size = 0: no windowing - each m6A site is its own feature and the
# matrix is the 0/1 site matrix itself.
#
# Feature matrix column names are the genomic midpoint of the window (the site
# position when window_size = 0), so every downstream plot can read genomic
# coordinates straight off the columns.
# ---------------------------------------------------------------------------
bin_read_matrix <- function(met_mat, region_start, region_end, window_size = 0) {
  M   <- as.matrix(met_mat)
  pos <- as.integer(colnames(M))

  if (window_size == 0) {
    anno <- data.frame(feature = seq_along(pos), win_start = pos, win_end = pos,
                       mid = pos, n_sites = 1L)
    colnames(M) <- as.character(pos)
    return(list(mat = M, window_anno = anno, window_size = 0))
  }

  starts <- seq(region_start, region_end, by = window_size)
  anno <- data.frame(feature   = seq_along(starts),
                     win_start = starts,
                     win_end   = pmin(starts + window_size - 1, region_end))
  anno$mid <- as.integer(floor((anno$win_start + anno$win_end) / 2))

  win_of_pos <- findInterval(pos, anno$win_start)
  stopifnot(all(win_of_pos >= 1), all(win_of_pos <= nrow(anno)))
  anno$n_sites <- tabulate(win_of_pos, nbins = nrow(anno))

  # sum of each read's site calls per occupied window, then divide by the
  # number of sites in that window -> mean call per window
  win_sums <- t(rowsum(t(M), group = win_of_pos))
  occupied <- as.integer(colnames(win_sums))
  mat <- sweep(win_sums, 2, anno$n_sites[occupied], "/")
  colnames(mat) <- as.character(anno$mid[occupied])
  rownames(mat) <- rownames(M)

  list(mat = mat, window_anno = anno, window_size = window_size)
}


# ---------------------------------------------------------------------------
# Step 2 (optional): KNN imputation of the bins still missing after binning,
# following scikit-learn's KNNImputer - donors for a missing feature are the
# rows that observe it, ranked by the nan-euclidean distance
#   d(i,j) = sqrt(n_features / n_co-observed * sum_co (x_i - x_j)^2),
# and the imputed value is the unweighted mean of the k nearest donors.
#
# Off by default: with full-span reads there is nothing to impute (see the
# header note).
# ---------------------------------------------------------------------------
knn_impute <- function(mat, k = 5, verbose = TRUE) {
  na_idx <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(na_idx) == 0) {
    if (verbose) cat("KNN imputation: no missing values, nothing to impute\n")
    return(mat)
  }
  if (verbose)
    cat(sprintf("KNN imputation (k = %d): %d missing values in %d of %d reads\n",
                k, nrow(na_idx), length(unique(na_idx[, "row"])), nrow(mat)))

  obs <- !is.na(mat)
  X0  <- mat; X0[!obs] <- 0
  n_feat <- ncol(mat)

  # nan-euclidean over co-observed features, by matrix algebra:
  # sum_co (x-y)^2 = sum_co x^2 + sum_co y^2 - 2 sum_co xy
  n_co  <- obs %*% t(obs)
  D2    <- (X0^2) %*% t(obs) + obs %*% t(X0^2) - 2 * (X0 %*% t(X0))
  D     <- sqrt(pmax(D2, 0) * n_feat / pmax(n_co, 1))
  D[n_co == 0] <- Inf
  diag(D) <- Inf

  for (col in unique(na_idx[, "col"])) {
    donors <- which(obs[, col])
    if (length(donors) == 0) next          # no read observes this feature
    rows <- na_idx[na_idx[, "col"] == col, "row"]
    for (r in rows) {
      d <- D[r, donors]
      use <- donors[order(d)][seq_len(min(k, sum(is.finite(d))))]
      if (length(use) > 0) mat[r, col] <- mean(mat[use, col])
    }
  }
  mat
}


# ---------------------------------------------------------------------------
# Steps 4-6: Manhattan distances between reads, KNN graph, exponential-kernel
# affinities as edge weights.
#
# k_neighbours is capped at nrow(mat) - 1. sigma = NULL uses the mean of the
# retained KNN distances; pass sigma = ncol(mat) for scikit-learn's
# laplacian_kernel default (gamma = 1 / n_features).
# ---------------------------------------------------------------------------
manhattan_knn_graph <- function(mat, k_neighbors = 50, sigma = NULL, verbose = TRUE) {
  n <- nrow(mat)
  if (n < 3) stop("fewer than 3 reads to cluster")
  k_eff <- min(k_neighbors, n - 1)
  if (verbose && k_eff < k_neighbors)
    cat(sprintf("KNN graph: k reduced from %d to %d (only %d reads)\n",
                k_neighbors, k_eff, n))

  D <- as.matrix(dist(mat, method = "manhattan"))
  diag(D) <- Inf

  # k nearest neighbours of every read
  nn <- t(apply(D, 1, function(d) order(d)[seq_len(k_eff)]))
  from <- rep(seq_len(n), each = k_eff)
  to   <- as.vector(t(nn))
  d_nn <- D[cbind(from, to)]

  if (is.null(sigma)) sigma <- mean(d_nn)
  if (!is.finite(sigma) || sigma <= 0) {
    warning("all KNN distances are 0 (identical reads); using sigma = 1")
    sigma <- 1
  }
  affinity <- exp(-d_nn / sigma)

  # union symmetrisation: Manhattan distance is symmetric, so the two
  # directions of a kept pair carry the same affinity and the union is the set
  # of unique unordered pairs
  a <- pmin(from, to); b <- pmax(from, to)
  keep <- !duplicated(a * n + b)

  edges <- data.frame(from = rownames(mat)[a[keep]],
                      to   = rownames(mat)[b[keep]],
                      weight = affinity[keep], stringsAsFactors = FALSE)
  g <- igraph::graph_from_data_frame(
    edges, directed = FALSE,
    vertices = data.frame(name = rownames(mat), stringsAsFactors = FALSE))

  if (verbose)
    cat(sprintf("KNN graph: %d reads, %d edges, sigma = %.4g, mean KNN distance = %.4g\n",
                n, igraph::ecount(g), sigma, mean(d_nn)))

  list(graph = g, sigma = sigma, k_eff = k_eff,
       mean_knn_dist = mean(d_nn), dist = D)
}


# ---------------------------------------------------------------------------
# Step 7: Leiden community detection with the RBConfigurationVertexPartition
# quality function (= igraph's modularity objective with a resolution
# parameter). Returns the membership vector named by read.
# ---------------------------------------------------------------------------
leiden_partition <- function(graph, resolution = 1, n_iterations = -1, seed = 1) {
  set.seed(seed)
  cl <- igraph::cluster_leiden(graph,
                               objective_function = "modularity",
                               weights    = igraph::E(graph)$weight,
                               resolution = resolution,
                               n_iterations = n_iterations)
  list(membership = igraph::membership(cl), quality = cl$quality)
}


# ---------------------------------------------------------------------------
# Steps 1-8 orchestrator.
#
# met_mat:  reads x m6A-site matrix, already reduced to the full-span read set
#           (filter_met_mat()) by the caller.
# rids_df:  read info (RID, chr, start, end, strand, sample_name, ...); extra
#           columns are carried into the assignment table, so every read keeps
#           its timepoint next to its cluster.
# window_size / k_neighbors / resolution / impute_missing / sigma: see above.
#
# Returns the feature matrix, the per-read assignments, the per-cluster mean
# profiles, the KNN graph and the parameters of the run.
# ---------------------------------------------------------------------------
leiden_manhattan_cluster <- function(met_mat, rids_df,
                                     region_start, region_end,
                                     window_size = 0,
                                     k_neighbors = 50,
                                     resolution = 1,
                                     impute_missing = FALSE, impute_k = 5,
                                     sigma = NULL, seed = 1, verbose = TRUE) {
  rids_df <- rids_df[match(rownames(met_mat), as.character(rids_df$RID)), ]
  stopifnot(!any(is.na(rids_df$RID)))

  ## step 1: read x feature matrix
  binned <- bin_read_matrix(met_mat, region_start, region_end, window_size)
  feat   <- binned$mat
  if (verbose) {
    n_empty <- sum(binned$window_anno$n_sites == 0)
    cat(sprintf("features: %d %s%s\n", ncol(feat),
                if (window_size == 0) "m6A sites (no windowing)"
                else paste0(window_size, "-bp windows"),
                if (n_empty > 0)
                  sprintf(" (%d window(s) with no m6A site dropped)", n_empty) else ""))
  }
  if (ncol(feat) < 2) stop("fewer than 2 informative features")

  ## step 2: optional KNN imputation
  if (impute_missing) feat <- knn_impute(feat, k = impute_k, verbose = verbose)
  n_na <- sum(is.na(feat))
  if (n_na > 0)
    warning(sprintf(paste("%d missing feature value(s) remain; dist() will",
                          "rescale the Manhattan distance over co-observed",
                          "features. Use full-span reads or impute_missing = TRUE."),
                    n_na))

  ## steps 4-6: Manhattan KNN graph with exponential-kernel affinities
  knn <- manhattan_knn_graph(feat, k_neighbors = k_neighbors, sigma = sigma,
                             verbose = verbose)

  ## step 7: Leiden
  part <- leiden_partition(knn$graph, resolution = resolution, seed = seed)
  memb <- part$membership[rownames(feat)]      # order by read, not by vertex id

  # renumber largest cluster first, so labels are deterministic
  ord <- names(sort(table(memb), decreasing = TRUE))
  cluster <- factor(paste0("cluster", match(as.character(memb), ord)),
                    levels = paste0("cluster", seq_along(ord)))
  if (verbose) {
    cat(sprintf("Leiden (resolution = %g): %d clusters\n", resolution, nlevels(cluster)))
    print(table(cluster))
  }

  ## per-read table: cluster + timepoint + coordinates for every read
  meta_cols <- intersect(c("sample_name", "chr", "start", "end", "strand"),
                         colnames(rids_df))
  assignments <- data.frame(RID = rownames(feat), cluster = cluster,
                            stringsAsFactors = FALSE)
  assignments <- cbind(assignments,
                       rids_df[match(assignments$RID, as.character(rids_df$RID)),
                               meta_cols, drop = FALSE])
  rownames(assignments) <- NULL

  ## step 8: cluster profiles = mean feature value per cluster
  profiles <- t(sapply(levels(cluster), function(cl)
    colMeans(feat[cluster == cl, , drop = FALSE], na.rm = TRUE)))
  colnames(profiles) <- colnames(feat)

  list(feat_mat    = feat,
       window_anno = binned$window_anno,
       assignments = assignments,
       profiles    = profiles,
       graph       = knn$graph,
       n_clusters  = nlevels(cluster),
       params      = list(window_size = window_size,
                          k_neighbors = k_neighbors, k_eff = knn$k_eff,
                          resolution = resolution, sigma = knn$sigma,
                          mean_knn_dist = knn$mean_knn_dist,
                          quality = part$quality,
                          impute_missing = impute_missing, impute_k = impute_k,
                          seed = seed,
                          region_start = region_start, region_end = region_end,
                          n_reads = nrow(feat), n_features = ncol(feat)))
}


# ---------------------------------------------------------------------------
# Per-cluster mean m6A call at each m6A site, computed at site resolution
# straight from met_mat (only the cluster labels are feature-derived). This is
# the topic model's cluster_met_profiles quantity, so the clusterings can be
# compared panel for panel at any window size.
# ---------------------------------------------------------------------------
cluster_site_profiles <- function(res, met_mat) {
  M  <- as.matrix(met_mat)
  df <- res$assignments
  df <- df[df$RID %in% rownames(M), ]
  do.call(rbind, lapply(levels(df$cluster), function(cl) {
    rids <- df$RID[df$cluster == cl]
    met  <- colMeans(M[rids, , drop = FALSE], na.rm = TRUE)
    met[is.nan(met)] <- NA
    data.frame(cluster = cl, pos = as.numeric(colnames(M)), met = met,
               n_reads = length(rids), row.names = NULL)
  }))
}
