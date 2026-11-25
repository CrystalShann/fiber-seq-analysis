#########
# Abundance-weighted Manhattan distance Q between cluster m6A profiles
#######
# 1. Check the binary read x base matrix and cluster labels
# 2. Compute cluster frequencies p_k
# 3. Average each cluster's reads into a mean m6A profile
# 4. Sum p_k * p_l * Manhattan distance between profiles, divide by window width
# 5. Return fiber count, cluster count and Q
# Main gene-level score: abundance-weighted Manhattan profile distance Q.
# Distances use the full genomic width, including all-zero bases.
manhattan_diversity <- function(mat, cluster) {
  mat <- as.matrix(mat)
  stopifnot(ncol(mat) > 0L, length(cluster) == nrow(mat),
            !anyNA(mat), all(mat %in% c(0, 1)), !anyNA(cluster),
            all(nzchar(as.character(cluster))))
  n <- nrow(mat); width <- ncol(mat)
  if (!n) return(data.frame(n_fibers = 0L, n_clusters = 0L,
    manhattan_between = NA_real_))
  cluster <- droplevels(factor(cluster))
  p <- as.numeric(table(cluster)) / n
  mu <- do.call(rbind, lapply(levels(cluster), function(k)
    colMeans(mat[cluster == k, , drop = FALSE])))
  q <- sum(outer(p, p) * as.matrix(dist(mu, method = "manhattan"))) / width
  data.frame(n_fibers = n, n_clusters = length(p), manhattan_between = q)
}

#########
# Self-test of the diversity functions on toy inputs
#######
# 1. Score small matrices with known Q values
# 2. Check that an extra all-zero base halves Q and one cluster gives Q = 0
# 3. Check the Q -> 1 / (1 - Q) conversion
# Known values exercise abundance weights, full-width scaling and the conversion.
validate_diversity <- function() {
  x <- matrix(c(0, 0, 1, 1, 0, 1, 1, 1), nrow = 4)
  r <- manhattan_diversity(x, c('a', 'a', 'a', 'b'))
  stopifnot(abs(r$manhattan_between - 0.1875) < 1e-12)
  # Adding an all-zero base halves the score for a one-base matrix.
  x <- matrix(c(0, 0, 0, 1), ncol = 1)
  r <- manhattan_diversity(x, c('a', 'a', 'b', 'c'))
  stopifnot(abs(r$manhattan_between - 0.375) < 1e-12,
    abs(manhattan_diversity(cbind(x, 0), c('a','a','b','c'))$manhattan_between - 0.1875) < 1e-12,
    manhattan_diversity(x, rep('a', 4))$manhattan_between == 0)
  z <- add_quadratic_diversity(data.frame(manhattan_between = c(0, .375, NA),
    status = c('valid', 'valid', 'upstream_failure')))
  stopifnot(identical(z$cluster_quadratic_diversity, c(1, 1.6, NA_real_)))
  invisible(TRUE)
}

#########
# Per-group summary statistics for each diversity metric
#######
# 1. For each metric and gene group, keep finite values
# 2. Report counts, mean, sd, median, quartiles, min and max
summarize_diversity <- function(d, metrics, groups) {
  do.call(rbind, lapply(metrics, function(m) do.call(rbind, lapply(groups, function(g) {
    values <- d[[m]][d$gene_group == g]; x <- values[is.finite(values)]
    stat <- function(f) if (length(x)) as.numeric(f(x)) else NA_real_
    data.frame(gene_group = g, metric = m, n_genes = length(values),
      n_defined = length(x), n_undefined = sum(!is.finite(values)),
      mean = stat(mean), sd = stat(sd), median = stat(median),
      q25 = stat(function(x) quantile(x, .25)), q75 = stat(function(x) quantile(x, .75)),
      min = stat(min), max = stat(max))
  }))))
}

#########
# Convert Manhattan Q into diversity 1 / (1 - Q)
#######
# 1. Check Q is finite and in [0, 1) for valid genes
# 2. Add cluster_quadratic_diversity = 1 / (1 - Q) for valid genes, NA otherwise
# Smouse, Banks & Peakall (2017), doi:10.1371/journal.pone.0185499,
# equations 1 and 3, adapted with cluster identity as allele identity.
# Q is the frequency-weighted Manhattan distance of cluster means
# divided by the fixed genomic width (manhattan_between).
# This is not the paper's hierarchical between-stratum beta diversity.
add_quadratic_diversity <- function(d) {
  d <- as.data.frame(d)
  stopifnot(all(c('manhattan_between', 'status') %in% names(d)))
  valid <- d$status == 'valid'
  stopifnot(!anyNA(valid), all(is.finite(d$manhattan_between[valid])),
    all(d$manhattan_between[valid] >= 0 & d$manhattan_between[valid] < 1))
  d$cluster_quadratic_diversity <- NA_real_
  d$cluster_quadratic_diversity[valid] <- 1 / (1 - d$manhattan_between[valid])
  d
}
