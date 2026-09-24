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
assemble_region_m6a <- function(sample_names = NULL, region_chr = NULL, region_start = NULL, region_end = NULL,
                                ft_result_dir = NULL, verbose = TRUE, full_span = FALSE,
                                sample_table = NULL, region = NULL, matrix_dir = NULL, reuse = TRUE) {
  if (full_span) {
    stopifnot(!is.null(sample_table), !is.null(region), nrow(region) == 1L,
      all(c("sample_name", "fire_dir") %in% names(sample_table)),
      !anyDuplicated(sample_table$sample_name),
      region$analysis_start == region$start + 1L, region$analysis_end == region$end,
      region$analysis_start <= region$analysis_end)
    paths <- file.path(sample_table$fire_dir, "extracted_results", "m6a_by_chr",
      paste0(sample_table$sample_name, ".ft_extracted_m6a.", region$chr, ".bed.gz"))
    required <- c(paths, paste0(paths, ".tbi"))
    if (!all(file.exists(required))) stop("Missing input: ", paste(required[!file.exists(required)], collapse = ", "))
    signature <- list(version = 1L, region = region, samples = sample_table)
    cache <- NULL
    if (!is.null(matrix_dir)) {
      region_matrix_dir <- file.path(matrix_dir, region$region_id)
      dir.create(region_matrix_dir, recursive = TRUE, showWarnings = FALSE)
      cache <- file.path(region_matrix_dir, paste0(region$region_id, "_m6a_matrix.rds"))
      if (reuse && file.exists(cache)) {
        previous <- readRDS(cache)
        if (identical(previous$signature[names(signature)], signature)) {
          stopifnot(!anyDuplicated(previous$rids_df$RID),
            all(previous$rids_df$start <= region$analysis_start),
            all(previous$rids_df$end >= region$analysis_end),
            identical(rownames(previous$met_mat), as.character(previous$rids_df$RID)))
          return(previous)
        }
      }
    }
    region_gr <- GenomicRanges::GRanges(region$chr,
      IRanges::IRanges(region$analysis_start, region$analysis_end))
    reads_list <- metadata_list <- qc_list <- vector("list", nrow(sample_table))
    for (sample_index in seq_len(nrow(sample_table))) {
      sample_name <- sample_table$sample_name[sample_index]
      bed <- read_ft_bed12(paths[sample_index], region_gr, longest_alignment = TRUE)
      overlapping <- nrow(bed)
      if (overlapping) bed <- bed[bed$start <= region$start & bed$end >= region$end, , drop = FALSE]
      qc_list[[sample_index]] <- data.frame(sample_name = sample_name,
        overlapping_reads = overlapping, full_span_reads = nrow(bed))
      if (!nrow(bed)) next
      original_ids <- as.character(bed$RID)
      bed$RID <- paste(sample_name, original_ids, sep = "::")
      metadata_list[[sample_index]] <- data.frame(
        RID = bed$RID, original_RID = original_ids, chr = bed$chr,
        start = bed$start + 1L, end = bed$end, strand = bed$strand,
        sample_name = sample_name, score = bed$score)
      blocks <- convert_ft_bed12_to_bed6(bed)
      if (any(blocks$end - blocks$start != 1L)) stop("Non-single-base m6A block in ", paths[sample_index])
      blocks <- blocks[blocks$end >= region$analysis_start & blocks$end <= region$analysis_end, , drop = FALSE]
      reads_list[[sample_index]] <- data.frame(RID = as.character(blocks$RID), pos = blocks$end)
    }
    rids_df <- dplyr::bind_rows(metadata_list)
    reads <- dplyr::bind_rows(reads_list)
    if (nrow(rids_df) < 3L || !nrow(reads)) stop("Insufficient full-span reads/m6A sites: ", region$region_id)
    met_mat <- get_sparse_met_mat(reads, rids_df, region$analysis_start, region$analysis_end,
      base = "A", require_full_span = TRUE)
    stopifnot(!anyNA(met_mat), !anyDuplicated(rownames(met_mat)))
    result <- list(met_mat = met_mat, rids_df = rids_df, qc = dplyr::bind_rows(qc_list),
      signature = signature)
    if (!is.null(cache)) {
      saveRDS(result, cache)
      data.table::fwrite(result$qc, file.path(region_matrix_dir, paste0(region$region_id, "_read_qc.tsv")), sep = "\t")
    }
    return(result)
  }
  
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
# run the full Leiden + Manhattan workflow across a table
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
# donors for a missing feature are the rows that observe it, 
# ranked by the nan-euclidean distance
# d(i,j) = sqrt(n_features / n_co-observed * sum_co (x_i - x_j)^2),
# and the imputed value is the unweighted mean of the k nearest donors
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



# ---------------------------------------------------------------------------------
# Region selection and haplotype phasing functions for LCL 

# LCL top-AS-FIRE selection, m6A matrices, and Leiden clustering
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# select siginficant AS-FIRE regions using combined_AS_fire_freq_merged_peaks_window2000_firefreq0.1_het_only.rds
# ---------------------------------------------------------------------------------

# /project/spott/kevinluo/Fiber_seq/results/QTL/fireQTL/AS_fire_freq_results/combined_results/combined_AS_fire_freq_merged_peaks_window2000_firefreq0.1_het_only.rds

select_lcl_top_asfire_regions <- function(asfire_results, top_n = 10L, width_bp = 2000L,
                                          q_cutoff = 0.1, min_coverage = 20,
                                          min_fire_freq = 0.1) {
  stopifnot(top_n >= 1L, top_n == as.integer(top_n), width_bp == 2000L)
  x <- as.data.frame(asfire_results)
  if (!"chr" %in% names(x) && "seqnames" %in% names(x)) x$chr <- as.character(x$seqnames)
  required <- c("chr", "start", "end", "peak_name", "rsID", "snp_pos", "ref", "alt",
                "fisher_pvalue", "coverage", "coverage_ref", "coverage_alt", "sample_name")
  if (!all(required %in% names(x))) stop("AS-FIRE input lacks: ", paste(setdiff(required, names(x)), collapse = ", "))
  if (!"fire_freq" %in% names(x)) {
    if (!"fire_coverage" %in% names(x)) stop("AS-FIRE input needs fire_freq or fire_coverage")
    x$fire_freq <- x$fire_coverage / x$coverage
  }
  # Match the coverage/FIRE-frequency filters before computing missing q-values
  valid <- is.finite(x$fisher_pvalue) & x$fisher_pvalue >= 0 &
    is.finite(x$coverage) & x$coverage > min_coverage &
    is.finite(x$coverage_ref) & x$coverage_ref > 0 &
    is.finite(x$coverage_alt) & x$coverage_alt > 0 &
    is.finite(x$fire_freq) & x$fire_freq > min_fire_freq
  x <- x[which(valid), , drop = FALSE]
  if (!nrow(x)) stop("No AS-FIRE pairs pass coverage and FIRE-frequency filters")
  x$fisher_pvalue <- pmin(x$fisher_pvalue, 1)
  q_source <- "input fisher_qvalue"
  if (!"fisher_qvalue" %in% names(x)) {
    if (!requireNamespace("qvalue", quietly = TRUE)) stop("Install qvalue to compute missing fisher_qvalue")
    x$fisher_qvalue <- qvalue::qvalue(x$fisher_pvalue)$qvalues
    q_source <- "qvalue over all coverage/FIRE-frequency-passing input pairs"
  }
  
  # supplying sig_ASfire_het_res.gr, keep its existing q-values
    # filter for SNPs with q < 0.01
  x <- x[which(is.finite(x$fisher_qvalue) & x$fisher_qvalue >= 0 &
                 x$fisher_qvalue < q_cutoff), , drop = FALSE]
  valid_snp <- !is.na(x$rsID) & grepl("^rs[0-9]+$", x$rsID) &
    !is.na(x$ref) & grepl("^[ACGT]$", x$ref) & !is.na(x$alt) & grepl("^[ACGT]$", x$alt) &
    x$ref != x$alt & is.finite(x$snp_pos) & x$snp_pos == as.integer(x$snp_pos) &
    x$snp_pos > width_bp / 2 & !is.na(x$peak_name) & nzchar(x$peak_name)
  x <- x[which(valid_snp), , drop = FALSE]
  
  # select one lead SNP for each FIRE peak used on the smallest p value
  x <- x[order(x$peak_name, x$fisher_pvalue, x$snp_pos, x$rsID), , drop = FALSE]
  x <- x[!duplicated(x$peak_name), , drop = FALSE]
  
  # rank the selected lead SNP by p value
  x <- x[order(x$fisher_pvalue, x$peak_name, x$snp_pos, x$rsID), , drop = FALSE]
  if (nrow(x) < top_n) stop("Only ", nrow(x), " eligible distinct peaks; requested ", top_n)
  x <- x[seq_len(top_n), , drop = FALSE]
  x$rank <- seq_len(nrow(x))
  x$qvalue_source <- q_source
  
  # keep the AS FIRE peak coordinates ex: chr1_11161_11494
  x$peak_start <- x$start
  x$peak_end <- x$end
  
  # get all the samples contain the SNP
  x$contributing_samples <- vapply(as.character(x$sample_name), function(z) {
    if (is.na(z)) stop("Missing contributing sample list")
    ids <- unique(trimws(strsplit(z, ",", fixed = TRUE)[[1]]))
    ids <- ids[nzchar(ids)]
    if (!length(ids)) stop("Empty contributing sample list")
    paste(ids, collapse = ",")
  }, character(1))
  x$n_listed_samples <- lengths(strsplit(x$contributing_samples, ",", fixed = TRUE))
  
  # construct 2000 bp window based on the coordinate of lead SNP
  x$analysis_start <- as.integer(x$snp_pos - width_bp / 2)
  x$analysis_end <- as.integer(x$analysis_start + width_bp - 1L)
  x$start <- x$analysis_start - 1L
  x$end <- x$analysis_end
  x$width <- width_bp
  
  # make region id
  x$region_id <- paste0("rank", sprintf("%02d", x$rank), "_", x$rsID, "_",
    gsub("[^A-Za-z0-9_-]", "_", x$peak_name))
  x$focal_snp <- x$rsID
  x$focal_pos <- as.integer(x$snp_pos)
  x$annotation <- paste0("Top AS-FIRE #", x$rank, ": ", x$rsID, " | ", x$peak_name)
  x$region_type <- "top_asfire_het"
  x$gene <- x$ensg <- NA_character_
  x$tss <- NA_integer_
  x$strand <- "*"
  x$coordinate_system <- "BED: 0-based, half-open; analysis: 1-based inclusive"
  x$overlapping_selected_windows <- vapply(seq_len(nrow(x)), function(i) {
    hit <- which(seq_len(nrow(x)) != i & x$chr == x$chr[i] &
      x$analysis_start <= x$analysis_end[i] & x$analysis_end >= x$analysis_start[i])
    paste(x$region_id[hit], collapse = ",")
  }, character(1))
  x$overlapping_selected_peaks <- vapply(seq_len(nrow(x)), function(i) {
    hit <- which(seq_len(nrow(x)) != i & x$chr == x$chr[i] &
      x$peak_start <= x$peak_end[i] & x$peak_end >= x$peak_start[i])
    paste(x$peak_name[hit], collapse = ",")
  }, character(1))
  
  # for each selected SNP, calculate its distance to the closest other selected SNP on the same chromosome
  x$nearest_selected_snp_distance_bp <- vapply(seq_len(nrow(x)), function(i) {
    other <- which(seq_len(nrow(x)) != i & x$chr == x$chr[i])
    if (!length(other)) return(NA_real_)
    min(abs(x$focal_pos[other] - x$focal_pos[i]))
  }, numeric(1))
  
  # if two different FIRE peak select the same lead SNP, if so, drop both calls
  x$same_focal_snp_as_another_peak <- duplicated(paste(x$chr, x$focal_pos)) |
    duplicated(paste(x$chr, x$focal_pos), fromLast = TRUE)
  rownames(x) <- NULL
  x
}

# ---------------------------------------------------------------------------------
# filter for LCL samples containing AS-FIRE SNP
# ---------------------------------------------------------------------------------

# creates a sample table for each AS-FIRE region centered on the SNP

lcl_region_samples <- function(region, sample_table) {
  if (is.null(region$contributing_samples)) return(sample_table)
  ids <- strsplit(region$contributing_samples, ",", fixed = TRUE)[[1]]
  missing <- setdiff(ids, sample_table$sample_name)
  if (length(missing)) stop("Contributing LCLs missing from sample table: ", paste(missing, collapse = ", "))
  sample_table[match(ids, sample_table$sample_name), , drop = FALSE]
}

# ---------------------------------------------------------------------------------
# Build m6a binary matrix, run leiden clustering using Manhattan distance
# ---------------------------------------------------------------------------------

# build m6a matrix for reads covering each complete region, then retains reads with lead SNP
# from heterozygous samples. Cluster both alleles together with Manhattan distance (bin = 0,
# k =10, resolution = 1)

run_lcl_clustering <- function(regions, sample_table, matrix_dir, output_dir,
                                window_size = 0L, k_neighbors = 10L,
                                leiden_resolution = 1, leiden_seed = 1L,
                                kernel_sigma = NULL, reuse_cache = TRUE, phase_cache = NULL) {
  output_dir <- lcl_output_path(output_dir)
  matrix_dir <- lcl_output_path(matrix_dir)
  dir.create(file.path(output_dir, "summary tables"), recursive = TRUE, showWarnings = FALSE)
  run_summary <- list()
  result_paths <- setNames(character(nrow(regions)), regions$region_id)
  for (i in seq_len(nrow(regions))) {
    region <- regions[i, , drop = FALSE]
    selected_samples <- lcl_region_samples(region, sample_table)
    dat <- assemble_region_m6a(sample_table = selected_samples, region = region,
      matrix_dir = matrix_dir, reuse = reuse_cache, full_span = TRUE)
    if (region$region_type != "top_asfire_het" || is.null(phase_cache))
      stop("This LCL workflow requires top AS-FIRE regions and focal-SNP phasing")
    dat <- lcl_filter_focal_heterozygotes(dat, region, selected_samples, phase_cache)
    signature <- list(version = "top10_asfire_het_2kb_v1", region = region,
      matrix = dat$signature, retained_reads = dat$rids_df,
      params = list(window_size, k_neighbors, leiden_resolution, leiden_seed, kernel_sigma))
    out <- lcl_output_path(file.path(output_dir, region$region_id, paste0("bin", window_size),
      paste0("k", k_neighbors), paste0("resolution", leiden_resolution), "tables"))
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    result_paths[i] <- file.path(out, "clustering.rds")
    if (file.exists(result_paths[i])) {
      result <- readRDS(result_paths[i])
      if (!identical(result$analysis_signature, signature))
        stop("Saved clustering inputs differ; use a new analysis output directory: ", out)
    } else {
      result <- leiden_manhattan_cluster(dat$met_mat, dat$rids_df, region$analysis_start,
        region$analysis_end, window_size = window_size, k_neighbors = k_neighbors,
        resolution = leiden_resolution, sigma = kernel_sigma, seed = leiden_seed, impute_missing = FALSE)
      meta <- dat$rids_df[match(result$assignments$RID, dat$rids_df$RID), , drop = FALSE]
      for (key in setdiff(names(meta), names(result$assignments))) result$assignments[[key]] <- meta[[key]]
      result$site_met_mat <- dat$met_mat
      result$region <- region
      result$variants <- dat$variants
      result$focal <- dat$focal
      result$analysis_signature <- signature
      result$params$region_chr <- region$chr
      result$params$coordinate_system <- "1-based inclusive"
      result$assignments$sample_label <- sub("_.*$", "", result$assignments$sample_name)
      result$assignments$region_id <- region$region_id
      result$assignments$annotation <- region$annotation
      saveRDS(result, result_paths[i])
    }
    run_summary[[i]] <- data.frame(region_id = region$region_id, rank = region$rank,
      peak_name = region$peak_name, focal_snp = region$focal_snp, fisher_pvalue = region$fisher_pvalue,
      fisher_qvalue = region$fisher_qvalue, asfire_coverage = region$coverage,
      asfire_coverage_ref = region$coverage_ref, asfire_coverage_alt = region$coverage_alt,
      n_listed_samples = region$n_listed_samples,
      n_retained_samples = length(unique(result$assignments$sample_name)),
      n_full_span_reads = nrow(dat$read_filter_audit), n_reads = nrow(result$assignments),
      n_excluded_reads = sum(!dat$read_filter_audit$retained_for_analysis),
      n_clusters = result$n_clusters, n_alleles = length(unique(result$assignments$allele_label)))
    data.table::fwrite(dplyr::bind_rows(run_summary), file.path(output_dir, "summary tables", "run_summary.tsv"), sep = "\t")
  }
  saveRDS(result_paths, file.path(output_dir, "summary tables", "result_paths.rds"))
  result_paths
}
