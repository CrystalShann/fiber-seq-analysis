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

#########
# GO biological-process terms ordered by where their genes sit in the diversity ranking
#######
# 1. Check the ranking: unique genes, ranks 1..N (1 = highest cluster_quadratic_diversity)
# 2. Annotate panel genes with GO BP terms (GOALL, ancestor-propagated, ENSEMBL keys)
# 3. Keep terms covering >= min_genes panel genes but not every panel gene
# 4. Collapse terms with identical member sets, keeping the one smallest genome-wide
# 5. Per term: members' median rank, shift from the panel median rank, rank percentile, median diversity
# Descriptive only: the universe is the ranked panel, pooled across gene groups, no tests.
# rank_shift < 0 means the term's genes sit toward high diversity, > 0 toward low.
go_term_rank_shift <- function(d, min_genes = 5L) {
  stopifnot(requireNamespace('org.Hs.eg.db', quietly = TRUE),
    requireNamespace('GO.db', quietly = TRUE), requireNamespace('AnnotationDbi', quietly = TRUE))
  d <- as.data.frame(d)
  N <- nrow(d)
  stopifnot(all(c('gene_id', 'gene_name', 'rank_highest_first', 'cluster_quadratic_diversity') %in% names(d)),
    !anyDuplicated(d$gene_id), setequal(d$rank_highest_first, seq_len(N)))
  org <- org.Hs.eg.db::org.Hs.eg.db
  ann <- suppressMessages(AnnotationDbi::select(org, keys = intersect(d$gene_id,
    AnnotationDbi::keys(org, keytype = 'ENSEMBL')), keytype = 'ENSEMBL', columns = c('GOALL', 'ONTOLOGYALL')))
  ann <- unique(ann[!is.na(ann$GOALL) & ann$ONTOLOGYALL == 'BP', c('ENSEMBL', 'GOALL')])
  sets <- lapply(split(ann$ENSEMBL, ann$GOALL), sort)
  sets <- sets[lengths(sets) >= min_genes & lengths(sets) < N]
  stopifnot(length(sets) > 0L)
  # Genome-wide BP size per term, used to keep the most specific of identical sets.
  genome_size <- lengths(AnnotationDbi::mget(names(sets), org.Hs.eg.db::org.Hs.egGO2ALLEGS))
  key <- vapply(sets, paste, character(1), collapse = ';')
  keep <- names(sets)[order(genome_size, names(sets))]
  keep <- keep[!duplicated(key[keep])]
  terms <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = names(sets),
    keytype = 'GOID', columns = 'TERM'))
  term_of <- setNames(terms$TERM, terms$GOID)
  panel_median <- stats::median(seq_len(N))
  go <- do.call(rbind, lapply(keep, function(id) {
    i <- match(sets[[id]], d$gene_id); r <- d$rank_highest_first[i]
    same <- setdiff(names(key)[key == key[[id]]], id)
    data.frame(go_id = id, term = term_of[[id]], n_genes = length(i),
      median_rank = stats::median(r), rank_shift = stats::median(r) - panel_median,
      rank_percentile = (stats::median(r) - 0.5) / N,
      median_diversity = stats::median(d$cluster_quadratic_diversity[i]),
      gene_names = paste(d$gene_name[i][order(r)], collapse = ';'),
      collapsed_terms = paste(term_of[same], collapse = ';'))
  }))
  go <- go[order(go$rank_shift, go$go_id), , drop = FALSE]
  rownames(go) <- NULL
  go
}

#########
# Member-gene ranks for the GO terms shifted most toward high and low diversity
#######
# 1. Take the n_each most negative and n_each most positive rank_shift terms
# 2. One row per term, one point per member gene at its diversity rank
# 3. Mark each term's median rank (black diamond) and the panel median (dashed line)
# Returns the ggplot; saving is left to the caller.
plot_go_term_rank_shift <- function(go, d, n_each = 15L) {
  d <- as.data.frame(d)
  top <- unique(rbind(head(go, n_each), tail(go, n_each)))
  pts <- do.call(rbind, lapply(seq_len(nrow(top)), function(j) {
    g <- strsplit(top$gene_names[j], ';', fixed = TRUE)[[1]]
    data.frame(label = sprintf('%s (n=%d)', top$term[j], top$n_genes[j]), gene_name = g,
      rank = d$rank_highest_first[match(g, d$gene_name)])
  }))
  lev <- sprintf('%s (n=%d)', top$term, top$n_genes)
  pts$label <- factor(pts$label, levels = rev(lev))
  med <- data.frame(label = factor(lev, levels = rev(lev)), rank = top$median_rank)
  ggplot2::ggplot(pts, ggplot2::aes(rank, label)) +
    ggplot2::geom_vline(xintercept = stats::median(seq_len(nrow(d))), linetype = 'dashed', colour = 'grey50') +
    ggplot2::geom_point(alpha = 0.5, size = 1.2, colour = 'grey30') +
    ggplot2::geom_point(data = med, shape = 18, size = 3.5, colour = 'black') +
    ggplot2::scale_x_continuous(limits = c(1, nrow(d)), breaks = c(1, 50, 100, nrow(d))) +
    ggplot2::labs(x = 'Diversity rank (1 = highest cluster_quadratic_diversity)', y = NULL) +
    ggplot2::theme_bw(base_size = 9)
}
