# Centered PCA -> exact correlation neighbors -> fuzzy graph -> Leiden.
# UMAP is for display only. All computation uses native R packages.
acf_correlation_neighbors <- function(scores, k, block_size = 256L) {
  centered <- scores - rowMeans(scores)
  unit <- centered / sqrt(rowSums(centered^2))
  n <- nrow(unit)
  stopifnot(all(is.finite(unit)), k >= 2L, k < n)
  indices <- matrix(0L, n, k)
  distances <- matrix(0, n, k)
  # Block multiplication avoids allocating a dense n x n distance matrix.
  for (first in seq.int(1L, n, by = block_size)) {
    rows <- seq.int(first, min(n, first + block_size - 1L))
    d <- 1 - tcrossprod(unit[rows, , drop = FALSE], unit)
    d[] <- pmax(0, pmin(2, d))
    for (j in seq_along(rows)) {
      i <- rows[j]
      d[j, i] <- Inf
      other <- order(d[j, ], seq_len(n))[seq_len(k - 1L)]
      indices[i, ] <- c(i, other)  # uwot counts self among k neighbors
      distances[i, ] <- c(0, d[j, other])
    }
  }
  list(idx = indices, dist = distances)
}

cluster_profiles <- function(profiles, valid, n_pcs = 50L, n_neighbors = 10L,
                             resolution = 0.4, seed = 0L) {
  stopifnot(is.matrix(profiles), is.logical(valid), length(valid) == nrow(profiles),
    !anyNA(valid), all(is.finite(profiles[valid, , drop = FALSE])),
    length(n_pcs) == 1L, is.finite(n_pcs), n_pcs >= 2L, n_pcs == floor(n_pcs),
    length(n_neighbors) == 1L, is.finite(n_neighbors), n_neighbors >= 2L,
    n_neighbors == floor(n_neighbors), length(resolution) == 1L,
    is.finite(resolution), resolution > 0, length(seed) == 1L, is.finite(seed))
  n <- nrow(profiles)
  result <- list(labels = rep("", n),
    status = ifelse(valid, "insufficient_reads", "zero_variance"),
    embedding = matrix(NA_real_, n, 2L),
    edges = data.frame(source = integer(), target = integer(), weight = numeric()),
    info = list(method = "R PCA -> correlation neighbors -> uwot fuzzy graph -> igraph Leiden; UMAP for display",
      requested_pcs = n_pcs, requested_neighbors = n_neighbors,
      resolution = resolution, seed = seed, n_clustered = 0L, n_clusters = 0L))
  indices <- which(valid)
  if (length(indices) < 4L) return(result)
  values <- profiles[indices, , drop = FALSE]
  if (all(abs(sweep(values, 2L, values[1L, ], "-")) <= 1e-12)) {
    result$status[indices] <- "identical_profiles"
    return(result)
  }
  pcs <- min(n_pcs, nrow(values) - 1L, ncol(values) - 1L)
  if (pcs < 2L) stop("At least three lag features are required for clustering")
  # Retain lag zero, center columns, and do not variance-scale the features.
  set.seed(seed)
  pca <- if (pcs < min(dim(values))/2) {
    irlba::prcomp_irlba(values, n = pcs, center = TRUE, scale. = FALSE, tol = 1e-8)
  } else {
    stats::prcomp(values, center = TRUE, scale. = FALSE, rank. = pcs)
  }
  total_variance <- sum(apply(values, 2L, stats::var))
  scores <- pca$x[, seq_len(pcs), drop = FALSE]
  # Correlation between score rows depends on component signs. Use a stable
  # loading-based convention: each component's largest absolute loading is positive.
  signs <- vapply(seq_len(pcs), function(j)
    sign(pca$rotation[which.max(abs(pca$rotation[,j])),j]), numeric(1))
  scores <- sweep(scores, 2L, signs, "*")
  usable <- apply(scores, 1L, function(x) all(is.finite(x)) && stats::sd(x) > 1e-12)
  result$status[indices[!usable]] <- "undefined_pca_correlation"
  indices <- indices[usable]
  scores <- scores[usable, , drop = FALSE]
  if (length(indices) < 4L) return(result)
  k <- min(n_neighbors, length(indices) - 1L)
  nn <- acf_correlation_neighbors(scores, k)
  set.seed(seed)
  embedding <- uwot::umap(scores, nn_method = nn, n_neighbors = k,
    n_components = 2L, init = "random", min_dist = 0.5, spread = 1,
    n_threads = 1L, n_sgd_threads = 1L, fast_sgd = FALSE,
    ret_extra = "fgraph", verbose = FALSE)
  adjacency <- embedding$fgraph
  stopifnot(all(is.finite(adjacency@x)), all(is.finite(embedding$embedding)))
  graph <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected",
    weighted = TRUE, diag = FALSE)
  set.seed(seed)
  partition <- igraph::cluster_leiden(graph, objective_function = "modularity",
    resolution = resolution, weights = igraph::E(graph)$weight, n_iterations = -1L)
  membership <- igraph::membership(partition)
  # Stable labels by decreasing size; report loader adds 1 for display.
  sizes <- table(membership)
  groups <- names(sizes)[order(-as.integer(sizes), as.integer(names(sizes)))]
  result$labels[indices] <- as.character(match(as.character(membership), groups) - 1L)
  result$status[indices] <- "clustered"
  result$embedding[indices, ] <- embedding$embedding
  edges <- as.data.frame(igraph::as_data_frame(graph, what = "edges"))
  result$edges <- data.frame(source = indices[as.integer(edges$from)],
    target = indices[as.integer(edges$to)], weight = edges$weight)
  result$info <- c(result$info[!names(result$info) %in% c("n_clustered", "n_clusters")],
    list(actual_pcs = pcs, actual_neighbors = k, n_clustered = length(indices),
      n_clusters = length(groups), pca_variance_explained = sum(pca$sdev[seq_len(pcs)]^2)/total_variance,
      neighbor_metric = "correlation", neighbor_representation = "PCA scores",
      leiden_flavor = "igraph", leiden_directed = FALSE,
      package_versions = c(igraph = as.character(utils::packageVersion("igraph")),
                           uwot = as.character(utils::packageVersion("uwot")))))
  result
}
