# hamming_leiden_functions.r
#
# Single-molecule Leiden clustering of Fiber-seq m6A reads on masked Hamming
# distances - a hybrid of the two existing clusterings:
#
#   * the distance is the hamming_hierarchical analysis's masked Hamming
#     (masked_hamming() in phasing_functions.r; Bellini, Do & Skok 2025),
#     computed on the same binary-bin feature matrices that analysis clusters
#     (fractional bins compress all masked Hamming distances into one blob -
#     see the hamming_hierarchical notebook - so bins are binarized: 1 = any
#     m6A call in the bin);
#   * the graph construction, community detection and number-of-clusters
#     selection are the Leiden_Manhattan analysis's (Raviram et al. 2026,
#     steps 5-7): KNN graph over reads, exponential-kernel affinities
#     exp(-d / sigma), Leiden with the RBConfigurationVertexPartition quality
#     at a set resolution. The number of clusters is whatever Leiden returns
#     at that resolution - not the hamming analysis's dispersion elbow.
#
# No KNN imputation step: masking is the masked Hamming distance's own way of
# handling missing features (each pair is compared over its co-observed
# features only), so nothing needs imputing even off the full-span read set.
#
# Note on full-span reads: with no missing feature values the masked Hamming
# distance is exactly the Manhattan distance divided by the number of
# features, and the default kernel width (mean retained KNN distance) is
# scale-free - so with binarize_bins = FALSE this pipeline reproduces
# Leiden_Manhattan exactly. The binary bins (and the masking, whenever
# partial-coverage reads are ever clustered) are what make it different.
#
# Source order (the notebook does this): topic_modelling_functions.r, then
# phasing_functions.r (masked_hamming()), then leiden_manhattan_functions.r
# (bin_read_matrix(), leiden_partition(), assemble_region_m6a(),
# cluster_site_profiles()), then leiden_manhattan_plots.r, then this file.
# phasing_functions.r must come BEFORE leiden_manhattan_plots.r: both define
# plot_cluster_met_profiles() / plot_cluster_composition(), and the Leiden
# versions (which take this file's result object) must win.

suppressMessages({
  requireNamespace("igraph")
  requireNamespace("Matrix")
})


# ---------------------------------------------------------------------------
# Masked-Hamming KNN graph with exponential-kernel affinities - the
# Leiden_Manhattan manhattan_knn_graph() with dist(method = "manhattan")
# replaced by masked_hamming(). sigma = NULL uses the mean of the retained
# KNN distances.
#
# masked_hamming() returns d = 0 for a pair sharing no observed feature (the
# Bellini paper's default for profile assignment); in a KNN graph that would
# link such reads maximally, so those pairs are set to Inf instead. Full-span
# reads never hit this.
# ---------------------------------------------------------------------------
hamming_knn_graph <- function(mat, k_neighbors = 50, sigma = NULL, verbose = TRUE) {
  n <- nrow(mat)
  if (n < 3) stop("fewer than 3 reads to cluster")
  k_eff <- min(k_neighbors, n - 1)
  if (verbose && k_eff < k_neighbors)
    cat(sprintf("KNN graph: k reduced from %d to %d (only %d reads)\n",
                k_neighbors, k_eff, n))

  D <- masked_hamming(mat)
  D[attr(D, "n_obs") == 0] <- Inf   # never link reads sharing no feature
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

  # union symmetrisation: masked Hamming distance is symmetric, so the two
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
# Orchestrator - the Leiden_Manhattan leiden_manhattan_cluster() with the
# Manhattan KNN graph replaced by the masked-Hamming one and the hamming
# analysis's binary bins.
#
# met_mat:  reads x m6A-site matrix, already reduced to the full-span read set
#           (filter_met_mat()) by the caller.
# rids_df:  read info (RID, chr, start, end, strand, sample_name, ...); extra
#           columns are carried into the assignment table.
# bin_size: bp per bin (bin_read_matrix()'s windowing: genomic-midpoint column
#           names, windows with no m6A site dropped); 0 = no binning, the 0/1
#           site matrix itself.
# binarize_bins: a read's bin value is 1 if the read has any m6A call in the
#           bin (the hamming analysis's setting); no effect at bin_size 0.
# k_neighbors / resolution / sigma / seed: see hamming_knn_graph() and
#           leiden_partition().
#
# Returns the same object shape as leiden_manhattan_cluster(), so every
# leiden_manhattan_plots.r function works on it unchanged (params keeps the
# name window_size for that reason).
# ---------------------------------------------------------------------------
hamming_leiden_cluster <- function(met_mat, rids_df,
                                   region_start, region_end,
                                   bin_size = 0,
                                   binarize_bins = TRUE,
                                   k_neighbors = 50,
                                   resolution = 1,
                                   sigma = NULL, seed = 1, verbose = TRUE) {
  rids_df <- rids_df[match(rownames(met_mat), as.character(rids_df$RID)), ]
  stopifnot(!any(is.na(rids_df$RID)))

  ## read x feature matrix
  binned <- bin_read_matrix(met_mat, region_start, region_end, bin_size)
  feat   <- binned$mat
  if (binarize_bins && bin_size > 0) feat[which(feat > 0)] <- 1
  if (verbose) {
    n_empty <- sum(binned$window_anno$n_sites == 0)
    cat(sprintf("features: %d %s%s%s\n", ncol(feat),
                if (bin_size == 0) "m6A sites (no binning)"
                else paste0(bin_size, "-bp bins"),
                if (binarize_bins && bin_size > 0) ", binarized" else "",
                if (n_empty > 0)
                  sprintf(" (%d bin(s) with no m6A site dropped)", n_empty) else ""))
  }
  if (ncol(feat) < 2) stop("fewer than 2 informative features")

  ## masked-Hamming KNN graph with exponential-kernel affinities
  knn <- hamming_knn_graph(feat, k_neighbors = k_neighbors, sigma = sigma,
                           verbose = verbose)

  ## Leiden (leiden_partition() from leiden_manhattan_functions.r)
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

  ## cluster profiles = mean feature value per cluster (with binary bins: the
  ## fraction of the cluster's reads with any m6A call in the bin)
  profiles <- t(sapply(levels(cluster), function(cl)
    colMeans(feat[cluster == cl, , drop = FALSE], na.rm = TRUE)))
  colnames(profiles) <- colnames(feat)

  list(feat_mat    = feat,
       window_anno = binned$window_anno,
       assignments = assignments,
       profiles    = profiles,
       graph       = knn$graph,
       n_clusters  = nlevels(cluster),
       params      = list(distance = "masked_hamming",
                          window_size = bin_size,
                          binarize_bins = binarize_bins,
                          k_neighbors = k_neighbors, k_eff = knn$k_eff,
                          resolution = resolution, sigma = knn$sigma,
                          mean_knn_dist = knn$mean_knn_dist,
                          quality = part$quality,
                          seed = seed,
                          region_start = region_start, region_end = region_end,
                          n_reads = nrow(feat), n_features = ncol(feat)))
}
