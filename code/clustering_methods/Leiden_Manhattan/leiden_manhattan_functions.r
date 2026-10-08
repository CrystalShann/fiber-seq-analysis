# leiden_manhattan_functions.r
#
# Single-molecule Leiden clustering of Fiber-seq m6A reads, replicating the
# nano-NOMe-seq read-clustering procedure of Raviram, Jiang, Schippke, Cova &
# Skok (2026), "Cohesin collisions maintain ordered nucleosome architecture at
# boundaries and promoters" (bioRxiv 2026.05.22.727261)



#   1. bin the per-read accessibility signal into indows across the 2-kb
#      region; a read's bin value is the mean methylation call in the bin;
#   2. balance conditions - sample an equal number of reads per condition -
#      and pool the sampled reads before clustering 
#   3. read-read similarity = Manhattan distance over the bins
#   4. KNN graph over reads (k = 10 nearest neighbours)
#   5. edge weights = affinities from an exponential kernel on those distances
#   6. Leiden community detection (leidenalg, RBConfigurationVertexPartition)
#      at a given resolution;
#   7. cluster profiles = mean accessibility over the reads of each cluster



suppressMessages({
  requireNamespace("igraph")
  requireNamespace("Matrix")
})


# ---------------------------------------------------------------------------
# Run the full Leiden + Manhattan workflow on a table of regions: for each
# region, assemble the full-span m6A matrix of all samples with
# assemble_region_m6a() (parsing_footprints_functions.r) and cluster it with
# leiden_manhattan_cluster(). Each sample's files are read from
# <ft_result_dir>/<sample>/extracted_results/m6a_by_chr/.
#
# Inputs:
#   sample_names  - sample folder names under ft_result_dir; also the factor
#                   levels of assignments$sample_name
#   regions       - data.frame with 1-based inclusive chr, start, end and an
#                   optional unique region_id (default "<chr>_<start>_<end>")
#   ft_result_dir - root of the per-sample fibertools results
#   window_size, k_neighbors, resolution, sigma, seed - passed to
#                   leiden_manhattan_cluster()
#   sample_table  - optional table with sample_name and timepoint; adds a
#                   timepoint column to the assignments
#   verbose       - print progress
# Output:
#   list named by region_id of leiden_manhattan_cluster() results, each also
#   holding region, met_mat and site_met_mat (the full-span read x m6A-site
#   matrix), params$region_chr / region_id, and assignments with
#   original_RID (and timepoint)
# ---------------------------------------------------------------------------
leiden_manhattan_cluster_regions <- function(sample_names, regions, ft_result_dir,
                                             window_size = 0,
                                             k_neighbors = 10,
                                             resolution = 1,
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

  assembly_samples <- data.frame(
    sample_name = as.character(sample_names),
    fire_dir = file.path(ft_result_dir, sample_names),
    stringsAsFactors = FALSE
  )
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

    assembly_region <- data.frame(
      region_id = region_id, chr = region_chr,
      start = region_start - 1L, end = region_end,
      analysis_start = region_start, analysis_end = region_end,
      stringsAsFactors = FALSE
    )
    dat <- assemble_region_m6a(
      full_span = TRUE, sample_table = assembly_samples, region = assembly_region,
      matrix_dir = NULL, reuse = FALSE
    )
    met_mat <- dat$met_mat
    rids_df <- dat$rids_df
    if (verbose) {
      cat(region_id, ": ", nrow(met_mat), " full-span reads x ", ncol(met_mat),
          " m6A sites\n", sep = "")
    }

    res <- leiden_manhattan_cluster(
      met_mat = met_mat,
      rids_df = rids_df,
      region_start = assembly_region$analysis_start,
      region_end = assembly_region$analysis_end,
      window_size = window_size,
      k_neighbors = k_neighbors,
      resolution = resolution,
      sigma = sigma,
      seed = seed,
      verbose = verbose
    )

    res$assignments$original_RID <- rids_df$original_RID[
      match(res$assignments$RID, rids_df$RID)
    ]
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
# window_size > 0: consecutive windows tiling [region_start, region_end]; a
# read's feature value is the mean of its m6A site calls in the window.
# Windows containing no m6A site carry no information for any read and are
# dropped; they stay in the returned window annotation with n_sites = 0.
# window_size = 0: no windowing - each m6A site is its own feature and the
# matrix is the 0/1 site matrix itself.
# Feature column names are the genomic midpoint of the window (the site
# position when window_size = 0).
#
# Inputs:
#   met_mat                  - reads x m6A sites matrix (colnames = positions)
#   region_start, region_end - 1-based window tiled when window_size > 0
#   window_size              - window width in bp; 0 = one feature per site
# Output:
#   list(mat = reads x features matrix, window_anno = data.frame(feature,
#   win_start, win_end, mid, n_sites) with one row per window, window_size)
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
# Steps 4-6: Manhattan distances between reads, KNN graph, exponential-kernel
# affinities exp(-distance / sigma) as edge weights. The graph joins two reads
# if either is among the other's k nearest neighbours. Each edge also keeps
# its original Manhattan distance (before the kernel) as the
# manhattan_distance attribute. k_neighbors is capped at nrow(mat) - 1.
#
# Inputs:
#   mat         - reads x features matrix (rownames = read IDs)
#   k_neighbors - neighbours per read
#   sigma       - kernel width; NULL = mean of the retained KNN distances.
#                 sigma = ncol(mat) gives scikit-learn's laplacian_kernel
#                 default (gamma = 1 / n_features). Falls back to 1 when all
#                 KNN distances are 0
#   verbose     - print the graph size and sigma
# Output:
#   list(graph = undirected igraph, vertex names = read IDs, edge attributes
#   weight and manhattan_distance; sigma; k_eff = k used; mean_knn_dist;
#   dist = full reads x reads distance matrix with an Inf diagonal)
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
                      weight = affinity[keep],
                      manhattan_distance = d_nn[keep], stringsAsFactors = FALSE)
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
# parameter), weighted by the edge affinities.
#
# Inputs:
#   graph        - igraph from manhattan_knn_graph() (edge attribute weight)
#   resolution   - higher gives more, smaller clusters
#   n_iterations - Leiden iterations; negative = until the partition is stable
#   seed         - random seed set before clustering
# Output:
#   list(membership = cluster number per read, named by read ID;
#   quality = quality of the partition)
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
# Rename clusters by size: the largest becomes "cluster1", the next
# "cluster2", and so on, so labels do not depend on Leiden's arbitrary
# numbering. Keep the existing size ordering (including table()'s order for
# ties).
#
# Inputs:
#   membership - cluster number per read
# Output:
#   factor of "cluster<N>" labels in the order of membership, levels
#   cluster1 ... clusterK
# ---------------------------------------------------------------------------
relabel_clusters_by_size <- function(membership) {
  ord <- names(sort(table(membership), decreasing = TRUE))
  factor(paste0("cluster", match(as.character(membership), ord)),
         levels = paste0("cluster", seq_along(ord)))
}


# ---------------------------------------------------------------------------
# Step 8: mean feature value (m6A fraction) of each cluster at every feature.
# Feature rows and cluster labels must already be in the same read order.
#
# Inputs:
#   feat_mat - reads x features matrix
#   cluster  - factor of cluster labels, one per row of feat_mat
#   na.rm    - ignore NA values in the means
# Output:
#   clusters x features matrix (rownames = cluster levels, colnames =
#   colnames(feat_mat))
# ---------------------------------------------------------------------------
cluster_feature_profiles <- function(feat_mat, cluster, na.rm = TRUE) {
  profiles <- t(sapply(levels(cluster), function(cl)
    colMeans(feat_mat[cluster == cl, , drop = FALSE], na.rm = na.rm)))
  colnames(profiles) <- colnames(feat_mat)
  profiles
}


# ---------------------------------------------------------------------------
# Steps 1-8 orchestrator for one region: bin the read x site matrix into
# features (bin_read_matrix), build the Manhattan KNN graph
# (manhattan_knn_graph), run Leiden (leiden_partition), name clusters by size
# (relabel_clusters_by_size) and average each cluster's features
# (cluster_feature_profiles). Warns when NA features remain.
#
# Inputs:
#   met_mat     - reads x m6A-site matrix from assemble_region_m6a(full_span =
#                 TRUE): longest alignments physically span the entire
#                 analysis window; reads with zero m6A calls in the window
#                 remain in the matrix
#   rids_df     - read info with RID (= rownames of met_mat), and optionally
#                 sample_name, chr, start, end, strand, which are copied into
#                 the assignment table
#   region_start, region_end - 1-based analysis window (used when
#                 window_size > 0)
#   window_size - feature window in bp; 0 = one feature per m6A site
#   k_neighbors, sigma - see manhattan_knn_graph()
#   resolution, seed   - see leiden_partition()
#   verbose     - print feature, graph and cluster summaries
# Output:
#   list(feat_mat, window_anno (from bin_read_matrix), assignments =
#   data.frame(RID, cluster, sample_name, chr, start, end, strand), profiles =
#   clusters x features means, graph, n_clusters, params = the settings plus
#   k_eff, sigma, mean_knn_dist, Leiden quality, region_start / region_end,
#   n_reads, n_features)
# ---------------------------------------------------------------------------
leiden_manhattan_cluster <- function(met_mat, rids_df,
                                     region_start, region_end,
                                     window_size = 0,
                                     k_neighbors = 50,
                                     resolution = 1,
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

  n_na <- sum(is.na(feat))
  if (n_na > 0)
    warning(sprintf(paste("%d missing feature value(s) remain; dist() will",
                          "rescale the Manhattan distance over co-observed",
                          "features. Use full-span reads."),
                    n_na))

  ## steps 4-6: Manhattan KNN graph with exponential-kernel affinities
  knn <- manhattan_knn_graph(feat, k_neighbors = k_neighbors, sigma = sigma,
                             verbose = verbose)

  ## step 7: Leiden
  part <- leiden_partition(knn$graph, resolution = resolution, seed = seed)
  memb <- part$membership[rownames(feat)]      # order by read, not by vertex id

  # renumber largest cluster first, so labels are deterministic
  cluster <- relabel_clusters_by_size(memb)
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
  profiles <- cluster_feature_profiles(feat, cluster)

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
                          seed = seed,
                          region_start = region_start, region_end = region_end,
                          n_reads = nrow(feat), n_features = ncol(feat)))
}



# ---------------------------------------------------------------------------------
# Region selection and haplotype phasing functions for LCL 

# LCL top-AS-FIRE selection, m6A matrices, and Leiden clustering
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Select the top LCL allele-specific FIRE (AS-FIRE) regions: filter peak x SNP
# tests on coverage and FIRE frequency, compute Fisher q-values when missing
# (qvalue package), keep q < q_cutoff with a valid rsID and SNV, take one lead
# SNP per FIRE peak (smallest p-value), rank peaks by p-value and keep the
# first top_n. Each region is a width_bp window centred on its lead SNP.
#
# The LCL workflow reads
#   /project/spott/kevinluo/Fiber_seq/results/QTL/fireQTL/AS_fire_freq_results/combined_results/combined_AS_fire_freq_merged_peaks_window2000_firefreq0.1_het_only.rds
#
# Inputs:
#   asfire_results - table of AS-FIRE tests with chr (or seqnames), start,
#                    end, peak_name, rsID, snp_pos, ref, alt, fisher_pvalue,
#                    coverage, coverage_ref, coverage_alt, sample_name
#                    (comma-separated samples carrying the SNP), and
#                    fire_freq or fire_coverage; fisher_qvalue is optional
#   top_n          - number of regions to keep
#   width_bp       - window width; must be 2000
#   q_cutoff       - keep fisher_qvalue < q_cutoff
#   min_coverage   - keep coverage > min_coverage (both alleles must be > 0)
#   min_fire_freq  - keep fire_freq > min_fire_freq
# Output:
#   data.frame, one row per region in rank order: the input columns plus
#   rank, region_id ("rank01_<rsID>_<peak>"), analysis_start / analysis_end
#   (1-based), start / end (BED), width, focal_snp, focal_pos,
#   contributing_samples, n_listed_samples, annotation, region_type =
#   "top_asfire_het", qvalue_source, peak_start / peak_end, and diagnostics
#   (overlapping windows / peaks, nearest selected SNP distance,
#   same_focal_snp_as_another_peak)
# ---------------------------------------------------------------------------
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
  
  # flag regions whose lead SNP is also the lead SNP of another selected peak
  # (flagged only, not dropped)
  x$same_focal_snp_as_another_peak <- duplicated(paste(x$chr, x$focal_pos)) |
    duplicated(paste(x$chr, x$focal_pos), fromLast = TRUE)
  rownames(x) <- NULL
  x
}

# ---------------------------------------------------------------------------
# Sample table for one AS-FIRE region: the rows of sample_table listed in the
# region's comma-separated contributing_samples (the LCLs carrying the SNP),
# in that order. Without a contributing_samples column the whole table is
# returned.
#
# Inputs:
#   region       - one region row, optionally with contributing_samples
#   sample_table - all samples (sample_name, fire_dir, ...)
# Output:
#   data.frame subset of sample_table; stops if a listed sample is missing
# ---------------------------------------------------------------------------
lcl_region_samples <- function(region, sample_table) {
  if (is.null(region$contributing_samples)) return(sample_table)
  ids <- strsplit(region$contributing_samples, ",", fixed = TRUE)[[1]]
  missing <- setdiff(ids, sample_table$sample_name)
  if (length(missing)) stop("Contributing LCLs missing from sample table: ", paste(missing, collapse = ", "))
  sample_table[match(ids, sample_table$sample_name), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# LCL AS-FIRE clustering: for each region, build the m6A matrix of reads
# covering the complete region (assemble_region_m6a, cached under
# matrix_dir), keep reads phased at the heterozygous lead SNP
# (lcl_filter_focal_heterozygotes() from LCL_phasing.r) and cluster both
# alleles together with leiden_manhattan_cluster() (defaults: bin = 0,
# k = 10, resolution = 1). A saved clustering.rds is reused only when its
# inputs and parameters are identical; otherwise the function stops. Needs
# lcl_output_path() from leiden_LCL.Rmd.
#
# Inputs:
#   regions      - rows from select_lcl_top_asfire_regions() (region_type
#                  "top_asfire_het")
#   sample_table - LCL samples (sample_name, fire_dir, phasing columns)
#   matrix_dir   - cache folder for assemble_region_m6a()
#   output_dir   - results folder
#   window_size, k_neighbors, leiden_resolution, leiden_seed, kernel_sigma -
#                  clustering parameters (see leiden_manhattan_cluster())
#   reuse_cache  - reuse matching m6A matrix caches
#   phase_cache  - per-sample haplotags from cache_lcl_haplotags()
# Output:
#   character vector of clustering.rds paths named by region_id. Writes
#     <output_dir>/<region_id>/bin<W>/k<K>/resolution<R>/tables/clustering.rds
#     <output_dir>/summary tables/run_summary.tsv
#     <output_dir>/summary tables/result_paths.rds
#     <matrix_dir>/<region_id>/<region_id>_m6a_matrix.rds and _read_qc.tsv
# ---------------------------------------------------------------------------
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
        resolution = leiden_resolution, sigma = kernel_sigma, seed = leiden_seed)
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
