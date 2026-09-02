# leiden_manhattan_functions.r
#
# Single-molecule Leiden clustering of Fiber-seq m6A reads, replicating the
# nano-NOMe-seq read-clustering procedure of Raviram, Jiang, Schippke, Cova &
# Skok (2026), "Cohesin collisions maintain ordered nucleosome architecture at
# boundaries and promoters" (bioRxiv 2026.05.22.727261)

# Paper's clustering method:

#   1. bin the per-read accessibility signal into indows across the 2-kb
#      region; a read's bin value is the mean methylation call in the bin;
#   2. impute the bins still missing after binning by KNN imputation - optional
#   3. balance conditions - sample an equal number of reads per condition -
#      and pool the sampled reads before clustering 
#   4. read-read similarity = Manhattan distance over the bins
#   5. KNN graph over reads (k = 50 nearest neighbours, scikit-learn)
#   6. edge weights = affinities from an exponential kernel on those distances
#   7. Leiden community detection (leidenalg, RBConfigurationVertexPartition)
#      at a given resolution;
#   8. cluster profiles = mean accessibility over the reads of each cluster
#
# Method:
#   * Reads: the full-span read set (filter_met_mat(), the topic model's row
#     filter. All four LPS timepoints are pooled
#   * No condition balancing (step 3)


suppressMessages({
  requireNamespace("igraph")
  requireNamespace("Matrix")
})


# ---------------------------------------------------------------------------
# Pooled m6A read x position matrix for one region, tagging each read with its
# sample of origin
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
# Convenience wrapper: run the full Leiden + Manhattan workflow across a table
# of explicit regions. Regions must contain chr/start/end, and may optionally
# carry a region_id column; otherwise coordinate-based IDs are created.
# ---------------------------------------------------------------------------
leiden_manhattan_cluster_regions <- function(sample_names, regions, ft_result_dir,
                                             window_size = 0,
                                             k_neighbors = 10,
                                             resolution = 1,
                                             impute_missing = FALSE,
                                             impute_k = 5,
                                             sigma = NULL,
                                             seed = 1,
                                             sample_table = NULL,
                                             verbose = TRUE) {
  required_cols <- c("chr", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(regions))
  if (length(missing_cols) > 0) {
    stop("regions is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  if (nrow(regions) == 0) stop("regions has no rows")

  if (length(window_size) != 1) {
    stop("window_size must be a single non-negative integer")
  }
  window_size <- as.integer(window_size)
  if (is.na(window_size) || window_size < 0) {
    stop("window_size must be a single non-negative integer")
  }

  regions <- as.data.frame(regions, stringsAsFactors = FALSE)
  regions$chr <- as.character(regions$chr)
  regions$start <- as.integer(regions$start)
  regions$end <- as.integer(regions$end)

  if (any(is.na(regions$start)) || any(is.na(regions$end))) {
    stop("regions$start and regions$end must be integer-like")
  }
  if (any(regions$start > regions$end)) {
    stop("each region must satisfy start <= end")
  }

  if (!"region_id" %in% colnames(regions)) {
    regions$region_id <- paste(regions$chr, regions$start, regions$end, sep = "_")
  }
  regions$region_id <- as.character(regions$region_id)
  if (anyDuplicated(regions$region_id)) {
    stop("region_id values must be unique")
  }

  region_results <- setNames(vector("list", nrow(regions)), regions$region_id)

  for (i in seq_len(nrow(regions))) {
    region <- regions[i, , drop = FALSE]
    region_id <- region$region_id[[1]]
    region_chr <- region$chr[[1]]
    region_start <- region$start[[1]]
    region_end <- region$end[[1]]

    if (verbose) {
      cat("\n=====", region_id, "=====\n")
      cat("Region:", region_chr, region_start, region_end, "\n")
    }

    dat <- assemble_region_m6a(
      sample_names = sample_names,
      region_chr = region_chr,
      region_start = region_start,
      region_end = region_end,
      ft_result_dir = ft_result_dir,
      verbose = verbose
    )
    if (verbose) {
      cat(region_id, ": ", nrow(dat$rids_df), " reads x ", ncol(dat$met_mat),
          " m6A sites before filtering\n", sep = "")
    }

    met_mat <- filter_met_mat(dat$met_mat, verbose = verbose)
    rids_df <- dat$rids_df[match(rownames(met_mat), as.character(dat$rids_df$RID)), ,
                           drop = FALSE]
    if (verbose) {
      cat(region_id, ": ", nrow(met_mat), " full-span reads x ", ncol(met_mat),
          " m6A sites\n", sep = "")
    }

    res <- leiden_manhattan_cluster(
      met_mat = met_mat,
      rids_df = rids_df,
      region_start = region_start,
      region_end = region_end,
      window_size = window_size,
      k_neighbors = k_neighbors,
      resolution = resolution,
      impute_missing = impute_missing,
      impute_k = impute_k,
      sigma = sigma,
      seed = seed,
      verbose = verbose
    )

    res$assignments$sample_name <- factor(res$assignments$sample_name,
                                          levels = sample_names)
    if (!is.null(sample_table) &&
        all(c("sample_name", "timepoint") %in% colnames(sample_table))) {
      res$assignments$timepoint <- sample_table$timepoint[
        match(res$assignments$sample_name, sample_table$sample_name)
      ]
    }

    res$region <- region
    res$met_mat <- met_mat
    res$site_met_mat <- met_mat
    res$params$region_chr <- region_chr
    res$params$region_id <- region_id
    region_results[[region_id]] <- res
  }

  region_results
}


# ---------------------------------------------------------------------------
# Step 1: read x feature matrix.
#
# window_size > 0: consecutive windows tiling [region_start, region_end]
# a read's feature value is the mean of its
# m6A site calls in the window. Windows containing no m6A site carry no
# information for any read and are dropped, and are recorded in the returned
# window annotation with n_sites = 0

# window_size = 0: no windowing - each m6A site is its own feature and the
# matrix is the 0/1 site matrix itself.
#
# Feature matrix column names are the genomic midpoint of the window (the site
# position when window_size = 0)
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
  # Counts how many m6A sites fall into each window
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
# knn_impute <- function(mat, k = 5, verbose = TRUE) {
#   na_idx <- which(is.na(mat), arr.ind = TRUE)
#   if (nrow(na_idx) == 0) {
#     if (verbose) cat("KNN imputation: no missing values, nothing to impute\n")
#     return(mat)
#   }
#   if (verbose)
#     cat(sprintf("KNN imputation (k = %d): %d missing values in %d of %d reads\n",
#                 k, nrow(na_idx), length(unique(na_idx[, "row"])), nrow(mat)))
# 
#   obs <- !is.na(mat)
#   X0  <- mat; X0[!obs] <- 0
#   n_feat <- ncol(mat)
# 
#   # nan-euclidean over co-observed features, by matrix algebra:
#   # sum_co (x-y)^2 = sum_co x^2 + sum_co y^2 - 2 sum_co xy
#   n_co  <- obs %*% t(obs)
#   D2    <- (X0^2) %*% t(obs) + obs %*% t(X0^2) - 2 * (X0 %*% t(X0))
#   D     <- sqrt(pmax(D2, 0) * n_feat / pmax(n_co, 1))
#   D[n_co == 0] <- Inf
#   diag(D) <- Inf
# 
#   for (col in unique(na_idx[, "col"])) {
#     donors <- which(obs[, col])
#     if (length(donors) == 0) next          # no read observes this feature
#     rows <- na_idx[na_idx[, "col"] == col, "row"]
#     for (r in rows) {
#       d <- D[r, donors]
#       use <- donors[order(d)][seq_len(min(k, sum(is.finite(d))))]
#       if (length(use) > 0) mat[r, col] <- mean(mat[use, col])
#     }
#   }
#   mat
# }


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
  # determine how many (k) neighbors can be used
  k_eff <- min(k_neighbors, n - 1)
  if (verbose && k_eff < k_neighbors)
    cat(sprintf("KNN graph: k reduced from %d to %d (only %d reads)\n",
                k_neighbors, k_eff, n))

  # calculate the Manhattan distance between every pair of reads
  D <- as.matrix(dist(mat, method = "manhattan"))
  # diagonal contains each read compared with itself, but a read can not count itself as its nearest neighbor
  # change 0 to inf
  diag(D) <- Inf

  # k nearest neighbours of every read
  nn <- t(apply(D, 1, function(d) order(d)[seq_len(k_eff)]))
  # turn the neighbor table into pairs
  from <- rep(seq_len(n), each = k_eff)
  # gets the corresponding nearest neighbor read ID
  to   <- as.vector(t(nn))
  # get distance for each pair
  d_nn <- D[cbind(from, to)]

  
  # convert distance into similarity - use the average distance among nearest neighbor pairs
  if (is.null(sigma)) sigma <- mean(d_nn)
  if (!is.finite(sigma) || sigma <= 0) {
    warning("all KNN distances are 0 (identical reads); using sigma = 1")
    sigma <- 1
  }
  affinity <- exp(-d_nn / sigma)

  # remove duplicate direction edges
  # ex: 1-->3 and 3-->1 is taken as the same representation
  # output is unique undirected pairs --> union KNN graph 
    # if A considers B a neighbor or B considers A a neighbor, A and B get connected
    # does not require both reads to choose each other
  a <- pmin(from, to); b <- pmax(from, to)
  keep <- !duplicated(a * n + b)

  # create the edge table using read ID
  edges <- data.frame(from = rownames(mat)[a[keep]],
                      to   = rownames(mat)[b[keep]],
                      weight = affinity[keep], stringsAsFactors = FALSE)
  # convert edge table into igraph object
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
# parameter). Returns the membership vector named by read
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
# straight from met_mat 
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
