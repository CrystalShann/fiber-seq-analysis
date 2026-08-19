# phasing_functions.r
#
# Single-molecule cluster-based phasing of Fiber-seq m6A reads, replicating the
# procedure of Bellini, Do & Skok (2025), "Phasing single-molecule nano-NOMe-seq
# reveals chromatin state heterogeneity in the context of transcription and
# long-range interactions" (bioRxiv 2025.09.08.674887), with the GpC methylation
# calls replaced by m6A calls (1 = m6A / accessible, 0 = covered but unmodified,
# NA = position not covered by the read). All parameters are the paper's; by
# default t is selected from the data instead of predefined (see notes), and
# two optional deviations are off by default: full_span_only (keep only reads
# covering every m6A site column - the topic model's filter_met_mat row
# filter) and binarize_bins (bin value = any m6A call in the bin instead of
# the mean site call):
#
#   50 bp bins; reads < 10 bp or with < 10 methylated sites removed; fully
#   methylated / fully unmethylated reads pre-assigned and excluded from
#   clustering; sliding windows of 60 valid bins with >= 20 overlapping bins;
#   reads with < 30 valid bins per window excluded; masked Hamming distance;
#   average-linkage agglomerative clustering to t clusters per window; Hungarian
#   matching of cluster means between adjacent windows; global reassignment of
#   every read to its closest stitched cluster profile; held-out long-read
#   validation accuracy A(r) = 1 - min masked Hamming distance (>= 50 kb reads).
#
# Input matrices come from get_sparse_met_mat() in
# code/topic_model/topic_modelling_functions.r (reads x m6A-site positions).
#
# Implementation notes (unstated details in the paper, chosen once and used
# everywhere):
#   * A read's value in a 50-bp bin is the mean of its observed site calls in
#     that bin (NA if the read observes no site there). The paper's masked
#     Hamming formula sum(m_i * |u_i - v_i|) / sum(m_i) is applied unchanged;
#     for binary vectors it is exactly the masked Hamming distance, and it is
#     the same formula the paper already applies to the (fractional) cluster
#     mean profiles during stitching.
#   * "Fully methylated" = every observed site on the read is methylated;
#     "fully unmethylated" = no observed site is methylated.
#   * The paper stops merging at a predefined t but never states how t was
#     chosen. Here t is not predefined: choose_t() cuts every window's
#     average-linkage tree at each candidate t (2..t_max) and picks the t
#     maximizing the read-weighted mean silhouette width across windows. A
#     single t is used in every window because the paper's one-to-one
#     Hungarian stitching requires equal cluster counts in adjacent windows.

suppressMessages({
  requireNamespace("clue")
  requireNamespace("cluster")
})


# ---------------------------------------------------------------------------
# Step 1: per-site coverage (mean coverage at each methylation site)
# ---------------------------------------------------------------------------
site_coverage_stats <- function(met_mat) {
  M <- as.matrix(met_mat)
  data.frame(
    pos       = as.numeric(colnames(M)),
    n_covered = colSums(!is.na(M)),
    n_meth    = colSums(M == 1, na.rm = TRUE),
    row.names = NULL
  )
}


# ---------------------------------------------------------------------------
# Step 2: remove low-information reads
# (read length < min_read_len bp, or < min_meth_sites methylated sites)
# ---------------------------------------------------------------------------
filter_informative_reads <- function(met_mat, rids_df,
                                     min_read_len = 10, min_meth_sites = 10) {
  M <- as.matrix(met_mat)
  rids_df <- rids_df[match(rownames(M), rids_df$RID), ]
  read_len <- rids_df$end - rids_df$start + 1
  n_meth   <- rowSums(M == 1, na.rm = TRUE)

  too_short <- read_len < min_read_len
  low_meth  <- !too_short & n_meth < min_meth_sites
  keep      <- !too_short & !low_meth

  list(keep_rids  = rownames(M)[keep],
       status     = ifelse(too_short, "removed_short_read",
                    ifelse(low_meth, "removed_low_meth", "kept")),
       n_meth     = n_meth,
       read_len   = read_len)
}


# ---------------------------------------------------------------------------
# Step 5: pre-assign fully methylated / fully unmethylated reads
# ---------------------------------------------------------------------------
classify_extreme_reads <- function(met_mat) {
  M <- as.matrix(met_mat)
  n_obs  <- rowSums(!is.na(M))
  n_meth <- rowSums(M == 1, na.rm = TRUE)
  list(fully_meth   = rownames(M)[n_obs > 0 & n_meth == n_obs],
       fully_unmeth = rownames(M)[n_obs > 0 & n_meth == 0])
}


# ---------------------------------------------------------------------------
# Steps 3-4: 50-bp binned ("smoothed") read x bin matrix
#
# Bins tile from the first read start to the last read end; only bins that
# contain at least one m6A site column can carry information, so only those
# appear as columns (absent bins are invalid bins by construction). Column
# names are the genomic start position of each bin.
# ---------------------------------------------------------------------------
bin_met_matrix <- function(met_mat, rids_df, bin_size = 50) {
  M <- as.matrix(met_mat)
  rids_df <- rids_df[match(rownames(M), rids_df$RID), ]
  site_pos <- as.numeric(colnames(M))

  grid_start <- min(rids_df$start)
  bin_idx    <- (site_pos - grid_start) %/% bin_size

  bin_ids <- sort(unique(bin_idx))
  bin_mat <- sapply(bin_ids, function(b) {
    v <- rowMeans(M[, bin_idx == b, drop = FALSE], na.rm = TRUE)
    v[is.nan(v)] <- NA
    v
  })
  if (is.null(dim(bin_mat))) bin_mat <- matrix(bin_mat, ncol = length(bin_ids))
  rownames(bin_mat) <- rownames(M)
  colnames(bin_mat) <- as.character(grid_start + bin_ids * bin_size)

  list(bin_mat    = bin_mat,
       bin_size   = bin_size,
       grid_start = grid_start,
       n_sites    = as.integer(table(factor(bin_idx, levels = bin_ids))))
}


# ---------------------------------------------------------------------------
# Step 8 / 15 / 19: masked Hamming distance between all row pairs of A and B
# (paper's formula: sum over co-observed positions of |u_i - v_i|, divided by
# the number of co-observed positions; 0 if there are none). The matrix of
# co-observed counts is returned as attribute "n_obs".
# ---------------------------------------------------------------------------
masked_hamming <- function(A, B = NULL) {
  if (is.null(dim(A))) A <- matrix(A, nrow = 1)
  if (is.null(B)) B <- A
  if (is.null(dim(B))) B <- matrix(B, nrow = 1)
  stopifnot(ncol(A) == ncol(B))

  num <- matrix(0, nrow(A), nrow(B))
  den <- matrix(0, nrow(A), nrow(B))
  for (j in seq_len(ncol(A))) {
    a <- A[, j]; b <- B[, j]
    oa <- !is.na(a); ob <- !is.na(b)
    if (!any(oa) || !any(ob)) next
    m  <- tcrossprod(as.numeric(oa), as.numeric(ob))
    a0 <- ifelse(oa, a, 0); b0 <- ifelse(ob, b, 0)
    num <- num + abs(outer(a0, b0, "-")) * m
    den <- den + m
  }
  D <- num / den
  D[den == 0] <- 0
  dimnames(D) <- list(rownames(A), rownames(B))
  attr(D, "n_obs") <- den
  D
}


# ---------------------------------------------------------------------------
# Step 6: overlapping sliding windows of `window_valid_bins` valid bins with
# >= `min_overlap_bins` shared bins between adjacent windows. `valid_cols` is
# the ordered vector of valid bin column indices; a list of index vectors is
# returned. A locus with <= window_valid_bins valid bins gives one window.
# ---------------------------------------------------------------------------
make_windows <- function(valid_cols, window_valid_bins = 60, min_overlap_bins = 20) {
  n <- length(valid_cols)
  if (n <= window_valid_bins) return(list(valid_cols))

  step   <- window_valid_bins - min_overlap_bins
  starts <- seq(1, n - window_valid_bins + 1, by = step)
  if (tail(starts, 1) != n - window_valid_bins + 1)
    starts <- c(starts, n - window_valid_bins + 1)
  lapply(starts, function(s) valid_cols[s:(s + window_valid_bins - 1)])
}


# ---------------------------------------------------------------------------
# Data-driven choice of t (the paper predefines t; here it is found instead).
# For every window, reads are filtered and clustered exactly as in
# cluster_window(); each candidate t in 2..t_max cuts the average-linkage tree
# at t and is scored by the average silhouette width on the window's masked
# Hamming distances. The t maximizing the read-weighted mean score across
# windows is returned and used in every window (the one-to-one Hungarian
# stitching requires the same t in adjacent windows).
# ---------------------------------------------------------------------------
choose_t <- function(bin_mat, windows, min_valid_bins = 30, t_max = 10) {
  per_win <- lapply(windows, function(w) {
    sub  <- bin_mat[, w, drop = FALSE]
    keep <- rowSums(!is.na(sub)) >= min_valid_bins
    if (sum(keep) < 3) return(NULL)
    sub <- sub[keep, , drop = FALSE]
    D   <- masked_hamming(sub)
    hc  <- hclust(as.dist(D), method = "average")
    ts  <- 2:min(t_max, nrow(sub) - 1)
    sil <- vapply(ts, function(t) {
      cl <- cutree(hc, k = t)
      mean(cluster::silhouette(cl, dmatrix = D)[, "sil_width"])
    }, numeric(1))
    list(n = nrow(sub), ts = ts, sil = sil)
  })
  per_win <- per_win[!vapply(per_win, is.null, logical(1))]
  if (length(per_win) == 0) stop("no window has enough reads to choose t")

  ts_all <- sort(unique(unlist(lapply(per_win, `[[`, "ts"))))
  score  <- vapply(ts_all, function(t) {
    sil <- vapply(per_win, function(x) x$sil[match(t, x$ts)], numeric(1))
    n   <- vapply(per_win, `[[`, numeric(1), "n")
    ok  <- !is.na(sil)
    sum(sil[ok] * n[ok]) / sum(n[ok])
  }, numeric(1))
  list(t = ts_all[which.max(score)], t_candidates = ts_all, score = score)
}


# ---------------------------------------------------------------------------
# Steps 7-14: cluster one window. Reads with < min_valid_bins observed bins in
# the window are excluded; the rest are clustered by average-linkage
# agglomerative clustering (hclust "average" implements exactly the
# size-weighted Lance-Williams update in the paper) on the masked Hamming
# distance, cut at t clusters. Returns the per-read labels and the cluster
# mean profiles over the window's bins.
# ---------------------------------------------------------------------------
cluster_window <- function(bin_mat, window_cols, t, min_valid_bins = 30) {
  sub  <- bin_mat[, window_cols, drop = FALSE]
  keep <- rowSums(!is.na(sub)) >= min_valid_bins
  if (sum(keep) < 2) return(NULL)
  sub <- sub[keep, , drop = FALSE]

  t_eff  <- min(t, nrow(sub))
  D      <- masked_hamming(sub)
  hc     <- hclust(as.dist(D), method = "average")
  labels <- cutree(hc, k = t_eff)

  means <- t(sapply(seq_len(t_eff), function(cl) {
    v <- colMeans(sub[labels == cl, , drop = FALSE], na.rm = TRUE)
    v[is.nan(v)] <- NA
    v
  }))
  rownames(means) <- as.character(seq_len(t_eff))
  colnames(means) <- colnames(sub)

  list(cols = window_cols, reads = rownames(sub), labels = labels, means = means)
}


# ---------------------------------------------------------------------------
# Steps 15-18: stitch adjacent windows by Hungarian matching of their cluster
# mean profiles (masked Hamming distance over the union of the two windows'
# bins; only co-observed - i.e. overlapping - bins contribute). Returns, per
# window, the local -> global cluster label map. Window 1's local labels are
# the global labels.
# ---------------------------------------------------------------------------
stitch_windows <- function(win_res) {
  M <- length(win_res)
  maps <- vector("list", M)
  maps[[1]] <- setNames(seq_len(nrow(win_res[[1]]$means)),
                        rownames(win_res[[1]]$means))
  if (M == 1) return(maps)

  for (N in 2:M) {
    A <- win_res[[N - 1]]$means
    rownames(A) <- as.character(maps[[N - 1]][rownames(A)])  # global labels
    B <- win_res[[N]]$means

    bins_u <- union(colnames(A), colnames(B))
    Au <- matrix(NA_real_, nrow(A), length(bins_u),
                 dimnames = list(rownames(A), bins_u))
    Au[, colnames(A)] <- A
    Bu <- matrix(NA_real_, nrow(B), length(bins_u),
                 dimnames = list(rownames(B), bins_u))
    Bu[, colnames(B)] <- B

    D <- masked_hamming(Au, Bu)   # paper's d = 0 default for zero-overlap pairs

    # pad to square so solve_LSAP always gets a complete assignment problem
    # (pad entries share one constant, which cannot change the optimal
    # matching among the real clusters)
    tA <- nrow(D); tB <- ncol(D); tt <- max(tA, tB)
    Dsq <- matrix(max(D) + 1, tt, tt)
    Dsq[seq_len(tA), seq_len(tB)] <- D
    sol <- clue::solve_LSAP(Dsq)

    map_cur <- integer(tB)
    for (i in seq_len(tt)) {
      j <- sol[i]
      if (i <= tA && j <= tB) map_cur[j] <- as.integer(rownames(D)[i])
    }
    if (any(map_cur == 0)) {   # clusters matched to a padded row: new labels
      used <- setdiff(seq_len(tt + max(unlist(maps[[N - 1]]))), map_cur)
      map_cur[map_cur == 0] <- used[seq_len(sum(map_cur == 0))]
    }
    maps[[N]] <- setNames(map_cur, rownames(win_res[[N]]$means))
  }
  maps
}


# ---------------------------------------------------------------------------
# Step 18 (result): global stitched cluster profiles across all bins. Where
# two overlapping windows both provide a mean for a (cluster, bin) pair, the
# means are averaged.
# ---------------------------------------------------------------------------
global_cluster_profiles <- function(win_res, maps) {
  all_bins <- sort(unique(as.numeric(unlist(
    lapply(win_res, function(w) colnames(w$means))))))
  t_total <- max(unlist(maps))

  P   <- matrix(0, t_total, length(all_bins),
                dimnames = list(paste0("cluster", seq_len(t_total)),
                                as.character(all_bins)))
  Cnt <- matrix(0, t_total, length(all_bins))

  for (N in seq_along(win_res)) {
    mu <- win_res[[N]]$means
    for (r in seq_len(nrow(mu))) {
      g  <- maps[[N]][rownames(mu)[r]]
      v  <- mu[r, ]
      ok <- !is.na(v)
      if (!any(ok)) next
      cols <- match(colnames(mu)[ok], colnames(P))
      P[g, cols]   <- P[g, cols] + v[ok]
      Cnt[g, cols] <- Cnt[g, cols] + 1
    }
  }
  out <- P / Cnt
  out[Cnt == 0] <- NA
  out
}


# ---------------------------------------------------------------------------
# Step 19: reassign every read to the closest global stitched profile by
# masked Hamming distance, with the paper's d = 0 default when a read and a
# profile share no observed bin (such ties resolve to the first profile).
# ---------------------------------------------------------------------------
assign_reads_to_profiles <- function(bin_mat, profiles) {
  common <- intersect(colnames(bin_mat), colnames(profiles))
  D <- masked_hamming(bin_mat[, common, drop = FALSE],
                      profiles[, common, drop = FALSE])
  idx <- apply(D, 1, which.min)
  data.frame(
    RID           = rownames(bin_mat),
    final_cluster = rownames(profiles)[idx],
    min_dist      = D[cbind(seq_len(nrow(D)), idx)],
    stringsAsFactors = FALSE, row.names = NULL
  )
}


# ---------------------------------------------------------------------------
# Steps 2-19 orchestrator.
#
# met_mat: reads x m6A-site-position matrix from get_sparse_met_mat()
# rids_df: read info (RID, chr, start, end, strand, ...) from
#          extract_ft_read_info(); extra columns (sample_name, ...) are carried
#          into the returned assignment table.
# t:       number of clusters per window; NULL (the default) selects t from
#          the data via choose_t().
# full_span_only: drop reads that do not cover every m6A site column (the
#          topic model's filter_met_mat row filter); dropped reads are
#          omitted from the output, as in the topic model.
# binarize_bins: a read's 50-bp bin value is 1 if any observed site in the
#          bin is methylated, 0 otherwise (instead of the mean site call).
# ---------------------------------------------------------------------------
phase_reads <- function(met_mat, rids_df, t = NULL,
                        bin_size = 50,
                        window_valid_bins = 60, min_overlap_bins = 20,
                        min_valid_bins_per_read = 30,
                        min_read_len = 10, min_meth_sites = 10,
                        full_span_only = FALSE, binarize_bins = FALSE,
                        verbose = TRUE) {
  M <- as.matrix(met_mat)
  rids_df <- rids_df[match(rownames(M), rids_df$RID), ]

  # -- optional: full-span reads only -------------------------------------
  if (full_span_only) {
    keep <- rowSums(is.na(M)) == 0
    if (verbose)
      cat(sprintf("full-span filter: %d of %d reads cover all sites, %d removed\n",
                  sum(keep), nrow(M), sum(!keep)))
    M <- M[keep, , drop = FALSE]
    rids_df <- rids_df[keep, , drop = FALSE]
    if (nrow(M) < 2) stop("fewer than 2 full-span reads")
  }

  site_cov <- site_coverage_stats(M)

  # -- step 2: low-information reads --------------------------------------
  filt <- filter_informative_reads(M, rids_df, min_read_len, min_meth_sites)
  M_f  <- M[filt$keep_rids, , drop = FALSE]
  if (verbose)
    cat(sprintf("reads: %d total, %d removed (<%d bp), %d removed (<%d m6A sites), %d kept\n",
                nrow(M), sum(filt$status == "removed_short_read"), min_read_len,
                sum(filt$status == "removed_low_meth"), min_meth_sites, nrow(M_f)))
  if (nrow(M_f) < 2) stop("fewer than 2 informative reads")

  # -- step 5: fully methylated / unmethylated reads ----------------------
  extreme <- classify_extreme_reads(M_f)
  cluster_rids <- setdiff(rownames(M_f),
                          c(extreme$fully_meth, extreme$fully_unmeth))
  if (verbose)
    cat(sprintf("pre-assigned: %d fully methylated, %d fully unmethylated; %d reads enter clustering\n",
                length(extreme$fully_meth), length(extreme$fully_unmeth),
                length(cluster_rids)))
  if (length(cluster_rids) < 2) stop("fewer than 2 heterogeneous reads to cluster")

  # -- steps 3-4: binned matrix (grid over all informative reads) ---------
  binned  <- bin_met_matrix(M_f, rids_df, bin_size)
  if (binarize_bins)
    binned$bin_mat[which(binned$bin_mat > 0)] <- 1
  bin_mat <- binned$bin_mat[cluster_rids, , drop = FALSE]

  # -- step 6: windows over valid bins ------------------------------------
  valid_cols <- which(colSums(!is.na(bin_mat)) > 0)
  if (verbose)
    cat(sprintf("bins: %d with sites, %d valid\n", ncol(bin_mat), length(valid_cols)))
  windows <- make_windows(valid_cols, window_valid_bins, min_overlap_bins)

  # -- t: selected from the data unless supplied --------------------------
  t_selection <- NULL
  if (is.null(t)) {
    t_selection <- choose_t(bin_mat, windows, min_valid_bins_per_read)
    t <- t_selection$t
    if (verbose)
      cat(sprintf("t = %d selected by silhouette (candidates %d..%d)\n",
                  t, min(t_selection$t_candidates),
                  max(t_selection$t_candidates)))
  }

  # -- steps 7-14: cluster each window ------------------------------------
  win_res <- lapply(windows, function(w)
    cluster_window(bin_mat, w, t, min_valid_bins_per_read))
  skipped <- vapply(win_res, is.null, logical(1))
  if (verbose && any(skipped))
    cat(sprintf("windows: %d of %d skipped (fewer than 2 reads with >= %d valid bins)\n",
                sum(skipped), length(windows), min_valid_bins_per_read))
  win_res <- win_res[!skipped]
  if (length(win_res) == 0) stop("no window could be clustered")

  # -- steps 15-18: stitch and build global profiles ----------------------
  maps     <- stitch_windows(win_res)
  profiles <- global_cluster_profiles(win_res, maps)

  # -- step 19: global reassignment of all clustering-eligible reads ------
  assign_df <- assign_reads_to_profiles(bin_mat, profiles)

  # -- assemble the full per-read table -----------------------------------
  full <- data.frame(RID = rownames(M), status = filt$status,
                     n_meth_sites = filt$n_meth, read_len = filt$read_len,
                     stringsAsFactors = FALSE)
  full$final_cluster <- full$status
  full$final_cluster[full$status == "kept"] <- NA
  full$final_cluster[full$RID %in% extreme$fully_meth]   <- "fully_methylated"
  full$final_cluster[full$RID %in% extreme$fully_unmeth] <- "fully_unmethylated"
  i <- match(assign_df$RID, full$RID)
  full$final_cluster[i] <- assign_df$final_cluster
  full$min_dist <- NA_real_
  full$min_dist[i] <- assign_df$min_dist
  full$n_valid_bins <- NA_integer_
  full$n_valid_bins[match(rownames(bin_mat), full$RID)] <-
    rowSums(!is.na(bin_mat))

  meta_cols <- intersect(c("sample_name", "chr", "start", "end", "strand"),
                         colnames(rids_df))
  full <- cbind(full, rids_df[match(full$RID, rids_df$RID), meta_cols,
                              drop = FALSE])
  rownames(full) <- NULL

  list(assignments = full,
       profiles    = profiles,
       bin_mat     = binned$bin_mat,      # all informative reads x bins
       bin_size    = bin_size,
       site_cov    = site_cov,
       windows     = windows,
       win_res     = win_res,
       stitch_maps = maps,
       extreme     = extreme,
       t_selection = t_selection,
       params      = list(t = t, bin_size = bin_size,
                          window_valid_bins = window_valid_bins,
                          min_overlap_bins = min_overlap_bins,
                          min_valid_bins_per_read = min_valid_bins_per_read,
                          min_read_len = min_read_len,
                          min_meth_sites = min_meth_sites,
                          full_span_only = full_span_only,
                          binarize_bins = binarize_bins))
}


# ---------------------------------------------------------------------------
# Steps 21-25: held-out long-read validation. Reads >= long_read_min_len are
# removed, the phasing is rebuilt without them, and each held-out read is
# compared to the resulting global profiles: A(r) = 1 - min masked Hamming
# distance. Returns NULL (with a message) if no read passes the length
# threshold - with the paper's 50 kb threshold that is the expected outcome
# for HiFi Fiber-seq libraries.
# ---------------------------------------------------------------------------
holdout_long_read_accuracy <- function(met_mat, rids_df, t = NULL,
                                       long_read_min_len = 50000,
                                       binarize_bins = FALSE, ...) {
  M <- as.matrix(met_mat)
  rids_df  <- rids_df[match(rownames(M), rids_df$RID), ]
  read_len <- rids_df$end - rids_df$start + 1
  long_ids <- rownames(M)[read_len >= long_read_min_len]
  if (length(long_ids) == 0) {
    message(sprintf("no read >= %d bp; long-read validation not applicable",
                    long_read_min_len))
    return(NULL)
  }

  res <- phase_reads(M[!rownames(M) %in% long_ids, , drop = FALSE],
                     rids_df[!rids_df$RID %in% long_ids, ], t = t,
                     binarize_bins = binarize_bins, ...)

  binned_long <- bin_met_matrix(M[long_ids, , drop = FALSE],
                                rids_df[rids_df$RID %in% long_ids, ],
                                res$bin_size)
  if (binarize_bins)
    binned_long$bin_mat[which(binned_long$bin_mat > 0)] <- 1
  # align the held-out reads' bins to the profile grid on genomic bin start
  common <- intersect(colnames(binned_long$bin_mat), colnames(res$profiles))
  if (length(common) == 0) {
    message("held-out reads share no bins with the fitted profiles")
    return(NULL)
  }
  D <- masked_hamming(binned_long$bin_mat[, common, drop = FALSE],
                      res$profiles[, common, drop = FALSE])
  acc <- 1 - apply(D, 1, min)

  data.frame(RID = long_ids, read_len = read_len[match(long_ids, rownames(M))],
             accuracy = acc, row.names = NULL)
}


# ===========================================================================
# Plotting
# ===========================================================================

PHASE_CLUSTER_COLORS <- c(
  "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
  "gold1", "skyblue2", "palegreen2", "#CAB2D6", "maroon",
  "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
  "green1", "yellow4", "yellow3", "darkorange4", "brown"
)

# numeric sort of "cluster<i>" labels (lexicographic sort puts cluster10
# before cluster2)
sort_cluster_labels <- function(labels) {
  cl <- unique(grep("^cluster", labels, value = TRUE))
  cl[order(as.integer(sub("^cluster", "", cl)))]
}

# color per final_cluster level, matching cluster1..t + the special groups
phase_cluster_palette <- function(levels) {
  cl  <- sort_cluster_labels(levels)
  pal <- setNames(PHASE_CLUSTER_COLORS[seq_along(cl)], cl)
  c(pal, fully_methylated = "grey20", fully_unmethylated = "grey60")[levels]
}

# order reads for heatmaps: cluster1..t, fully methylated, fully unmethylated;
# within a cluster by genomic start, then distance to its profile. Reads
# shorter than the locus otherwise interleave their uncovered flanks and the
# heatmap reads as scattered NA-grey lines.
order_reads_by_cluster <- function(assign_df) {
  lv <- c(sort_cluster_labels(unique(assign_df$final_cluster)),
          "fully_methylated", "fully_unmethylated")
  o <- order(match(assign_df$final_cluster, lv), assign_df$start,
             assign_df$min_dist)
  assign_df[o[assign_df$final_cluster[o] %in% lv], ]
}

# Binned (50-bp) heatmap of reads ordered by final cluster, NA grey. With
# binarize_bins a cell is 1 if the read has any m6A call in the bin; without it
# a cell is the read's mean site call in the bin (white = protected, dark =
# accessible). The ramp saturates at zmax: fractional bins rarely exceed ~0.3,
# so a 0-1 ramp would leave almost all real signal near-white.
plot_phasing_bin_heatmap <- function(res, main = NULL, zmax = 0.5) {
  df <- order_reads_by_cluster(res$assignments)
  mat <- res$bin_mat[df$RID, , drop = FALSE]
  binary <- isTRUE(res$params$binarize_bins)
  gplots::heatmap.2(
    pmin(mat, zmax), Rowv = FALSE, Colv = FALSE, dendrogram = "none", scale = "none",
    col = colorRampPalette(c("white", "#08306B"))(50),
    breaks = seq(0, zmax, length.out = 51), na.color = "grey85",
    RowSideColors = phase_cluster_palette(df$final_cluster),
    labRow = FALSE, labCol = FALSE, main = main,
    margins = c(3, 4), trace = "none", keysize = 0.9, key = TRUE,
    key.title = NA, density.info = "none",
    key.xlab = if (binary) "any m6A in bin (0/1)" else "mean m6A / bin"
  )
}

# Site-level heatmap (1 = m6A call black, 0 = protected white, NA grey).
plot_phasing_site_heatmap <- function(res, met_mat, main = NULL) {
  df  <- order_reads_by_cluster(res$assignments)
  df  <- df[df$RID %in% rownames(met_mat), ]
  mat <- as.matrix(met_mat)[df$RID, , drop = FALSE]
  gplots::heatmap.2(
    mat, Rowv = FALSE, Colv = FALSE, dendrogram = "none", scale = "none",
    col = c("white", "black"), breaks = c(-0.5, 0.5, 1.5), na.color = "grey92",
    RowSideColors = phase_cluster_palette(df$final_cluster),
    labRow = FALSE, labCol = FALSE, main = main,
    margins = c(3, 4), trace = "none", keysize = 0.5, key = FALSE
  )
}

# Global stitched cluster profiles, one facet per cluster.
plot_global_profiles <- function(res, tss = NULL, main = NULL) {
  P <- res$profiles
  df <- data.frame(
    cluster = rep(rownames(P), times = ncol(P)),
    pos     = rep(as.numeric(colnames(P)) + res$bin_size / 2, each = nrow(P)),
    value   = as.vector(P)
  )
  n <- table(res$assignments$final_cluster)
  df$cluster <- factor(df$cluster, levels = rownames(P),
                       labels = sprintf("%s (n=%d)", rownames(P),
                                        ifelse(is.na(n[rownames(P)]), 0L,
                                               as.integer(n[rownames(P)]))))
  ylab <- if (isTRUE(res$params$binarize_bins))
    sprintf("fraction of reads with m6A per %d-bp bin", res$bin_size)
  else
    sprintf("mean m6A per %d-bp bin", res$bin_size)
  gg <- ggplot(df[!is.na(df$value), ], aes(x = pos, y = value)) +
    geom_line(color = "#08306B") +
    facet_wrap(~ cluster, ncol = 1, strip.position = "right") +
    labs(x = "genomic position", y = ylab, title = main) +
    theme_cowplot(font_size = 10)
  if (!is.null(tss))
    gg <- gg + geom_vline(xintercept = tss, linetype = "dashed",
                          color = "grey40")
  gg
}

# Per-cluster m6A methylation proportion at each site, as in the topic model's
# cluster_met_profiles plot: one panel per cluster, bar height = mean call over
# that cluster's reads at that m6A site. Unlike plot_global_profiles (which
# draws the stitched cluster x bin profiles) this is computed at site
# resolution straight from met_mat; only the cluster labels are bin-derived.
# Panels follow the heatmap blocks, so the pre-assigned fully methylated /
# unmethylated groups get a panel each.
plot_cluster_met_profiles <- function(res, met_mat, tss = NULL, main = NULL) {
  M  <- as.matrix(met_mat)
  df <- res$assignments
  lv <- c(sort_cluster_labels(unique(df$final_cluster)),
          "fully_methylated", "fully_unmethylated")
  lv <- lv[lv %in% df$final_cluster]
  df <- df[df$final_cluster %in% lv & df$RID %in% rownames(M), ]
  if (nrow(df) == 0) stop("no assigned read is present in met_mat")

  pos <- as.numeric(colnames(M))
  pal <- phase_cluster_palette(lv)

  p_list <- lapply(seq_along(lv), function(i) {
    rids <- df$RID[df$final_cluster == lv[i]]
    met  <- colMeans(M[rids, , drop = FALSE], na.rm = TRUE)
    met[is.nan(met)] <- NA          # sites no read in the cluster covers
    gg <- ggplot(data.frame(pos = pos, met = met), aes(x = pos, y = met)) +
      geom_col(fill = pal[i]) +
      ylim(0, 1) +
      ylab("met prop.") + xlab("pos") +
      ggtitle(sprintf("%s (n=%d)", lv[i], length(rids))) +
      theme_cowplot(font_size = 10) +
      theme(plot.title = element_text(hjust = 0.5))
    if (!is.null(tss))
      gg <- gg + geom_vline(xintercept = tss, linetype = "dashed",
                            color = "grey40")
    gg
  })
  cowplot::plot_grid(plotlist = p_list, ncol = 1)
}

# Cluster composition by sample (LPS timepoint): stacked proportions.
plot_cluster_composition <- function(res, main = NULL) {
  df <- res$assignments
  df <- df[grepl("^cluster|^fully", df$final_cluster) & !is.na(df$sample_name), ]
  lv <- c(sort_cluster_labels(unique(df$final_cluster)),
          "fully_methylated", "fully_unmethylated")
  lv <- lv[lv %in% df$final_cluster]
  df$final_cluster <- factor(df$final_cluster, levels = lv)
  ggplot(df, aes(x = sample_name, fill = final_cluster)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = phase_cluster_palette(lv)) +
    labs(x = NULL, y = "fraction of reads", fill = "cluster", title = main) +
    theme_cowplot(font_size = 10)
}
