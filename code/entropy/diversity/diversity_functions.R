# Exact empirical Rao diversity: ordered pairs, including self-pairs (N^2).
# For binary X, |Xi-Xj| = (Xi-Xj)^2 at each position.
# All distances are divided by the full genomic width, including all-zero bases.
rao_diversity <- function(mat, cluster) {
  mat <- as.matrix(mat)
  stopifnot(ncol(mat) > 0L, length(cluster) == nrow(mat),
            !anyNA(mat), all(mat %in% c(0, 1)), !anyNA(cluster),
            all(nzchar(as.character(cluster))))
  n <- nrow(mat); width <- ncol(mat)
  if (!n) return(data.frame(n_fibers = 0L, n_clusters = 0L,
    rao_total = NA_real_, rao_between = NA_real_,
    between_cluster_share = NA_real_, q_over_g = NA_real_))
  cluster <- droplevels(factor(cluster))
  p <- as.numeric(table(cluster)) / n
  mu <- do.call(rbind, lapply(levels(cluster), function(k)
    colMeans(mat[cluster == k, , drop = FALSE])))
  overall <- colMeans(mat)
  total <- 2 * sum(overall * (1 - overall)) / width
  # Compute B independently from pairwise squared distances of cluster means.
  between <- sum(outer(p, p) * as.matrix(dist(mu))^2) / width
  # Within is sum_k p_k Q_k, not sum_k p_k^2 Q_k.
  within <- sum(p * (2 * rowSums(mu * (1 - mu)))) / width
  stopifnot(abs(total - within - between) < 1e-12,
            between >= 0, between <= total + 1e-12)
  g <- 1 - sum(p^2)
  data.frame(n_fibers = n, n_clusters = length(p), rao_total = total,
    rao_between = between,
    between_cluster_share = if (total > 0) between / total else NA_real_,
    q_over_g = if (g > 0) between / g else NA_real_)
}

# Deterministic small fixtures verify the N^2 convention against literal distances,
# unequal cluster weights, constant profiles, one cluster, and singleton clusters.
validate_rao <- function() {
  x <- as.matrix(expand.grid(a = 0:1, b = 0:1, c = 0:1))
  for (cl in list(rep('a', 8), as.character(1:8), c(rep('a', 5), rep('b', 3)))) {
    r <- rao_diversity(x, cl)
    stopifnot(abs(r$rao_total - sum(as.matrix(dist(x, method = 'manhattan'))) /
      (nrow(x)^2 * ncol(x))) < 1e-12)
    p <- as.numeric(table(cl)) / length(cl)
    within <- sum(p * vapply(sort(unique(cl)), function(k) {
      z <- x[cl == k, , drop = FALSE]
      sum(as.matrix(dist(z, method = 'manhattan'))) / (nrow(z)^2 * ncol(z))
    }, numeric(1)))
    stopifnot(abs(r$rao_total - within - r$rao_between) < 1e-12)
  }
  r <- rao_diversity(matrix(0, 3, 4), rep('a', 3))
  stopifnot(r$rao_total == 0, is.na(r$between_cluster_share), is.na(r$q_over_g))
  r <- rao_diversity(matrix(c(0, 1), 2, 1), c('a', 'b'))
  stopifnot(r$rao_total == .5, r$rao_between == .5, r$q_over_g == 1)
  invisible(TRUE)
}

compare_diversity_groups <- function(d, metrics, groups) {
  omnibus <- do.call(rbind, lapply(metrics, function(m) {
    z <- d[is.finite(d[[m]]), , drop = FALSE]
    p <- if (length(unique(z[[m]])) < 2L) 1 else
      kruskal.test(z[[m]], factor(z$gene_group))$p.value
    data.frame(metric = m, kruskal_p = p)
  }))
  omnibus$kruskal_p_BH <- p.adjust(omnibus$kruskal_p, 'BH')
  pairs <- do.call(rbind, lapply(metrics, function(m) {
    do.call(rbind, lapply(combn(groups, 2, simplify = FALSE), function(g) {
      a <- d[[m]][d$gene_group == g[1]]; a <- a[is.finite(a)]
      b <- d[[m]][d$gene_group == g[2]]; b <- b[is.finite(b)]
      p <- if (!length(a) || !length(b)) NA_real_ else if (length(unique(c(a,b))) == 1) 1 else
        wilcox.test(a, b, exact = FALSE)$p.value
      data.frame(metric = m, group1 = g[1], group2 = g[2], p_value = p)
    }))
  }))
  pairs$p_BH <- p.adjust(pairs$p_value, 'BH') # All 15 pairwise comparisons.
  list(omnibus = omnibus, pairwise = pairs)
}

summarize_diversity <- function(d, metrics, groups, tests) {
  do.call(rbind, lapply(metrics, function(m) do.call(rbind, lapply(groups, function(g) {
    values <- d[[m]][d$gene_group == g]; x <- values[is.finite(values)]
    stat <- function(f) if (length(x)) as.numeric(f(x)) else NA_real_
    data.frame(gene_group = g, metric = m, n_genes = length(values),
      n_defined = length(x), n_undefined = sum(!is.finite(values)),
      mean = stat(mean), sd = stat(sd), median = stat(median),
      q25 = stat(function(x) quantile(x, .25)), q75 = stat(function(x) quantile(x, .75)),
      min = stat(min), max = stat(max),
      kruskal_p = tests$kruskal_p[match(m, tests$metric)],
      kruskal_p_BH = tests$kruskal_p_BH[match(m, tests$metric)])
  }))))
}

# GO Biological Process ORA within the preselected scored gene panel.
# Exact upper-tail hypergeometric tests; no permutations or resampling.
rao_go_enrichment <- function(d) {
  stopifnot(requireNamespace('org.Hs.eg.db', quietly = TRUE),
    requireNamespace('GO.db', quietly = TRUE), requireNamespace('AnnotationDbi', quietly = TRUE))
  d <- as.data.frame(d)
  d <- d[is.finite(d$rao_total), , drop = FALSE]
  d <- d[order(d$rao_total, d$gene_id), , drop = FALSE]
  stopifnot(nrow(d) >= 3L, !anyDuplicated(d$gene_id))
  # Fixed pooled thirds (148 genes: 49 low, 50 intermediate, 49 high). No gene-group split.
  n_end <- round(nrow(d) / 3)
  d$rao_category <- rep(c('low', 'intermediate', 'high'), c(n_end, nrow(d) - 2L * n_end, n_end))
  d$rao_rank <- seq_len(nrow(d))
  stopifnot(d$rao_total[n_end] < d$rao_total[n_end + 1L],
            d$rao_total[nrow(d) - n_end] < d$rao_total[nrow(d) - n_end + 1L])
  org <- org.Hs.eg.db::org.Hs.eg.db
  known <- intersect(d$gene_id, AnnotationDbi::keys(org, keytype = 'ENSEMBL'))
  ann <- AnnotationDbi::select(org, keys = known, keytype = 'ENSEMBL',
    columns = c('GOALL', 'ONTOLOGYALL'))
  # GOALL propagates annotations to ancestor terms; unique Ensembl genes are the units.
  ann <- unique(ann[!is.na(ann$GOALL) & ann$ONTOLOGYALL == 'BP', c('ENSEMBL', 'GOALL')])
  d$has_bp_annotation <- d$gene_id %in% ann$ENSEMBL
  sets <- split(ann$ENSEMBL, ann$GOALL)
  # Avoid single-gene anecdotes; test every term covering >=3 panel genes,
  # including terms with zero foreground hits. Exclude universal terms.
  sets <- sets[lengths(sets) >= 3L & lengths(sets) < nrow(d)]
  stopifnot(length(sets) > 0L)
  descriptions <- AnnotationDbi::select(GO.db::GO.db, keys = names(sets),
    keytype = 'GOID', columns = 'TERM')
  results <- do.call(rbind, lapply(c('low','intermediate','high'), function(category) {
    ids <- d$gene_id[d$rao_category == category]
    do.call(rbind, lapply(names(sets), function(term) {
      members <- sets[[term]]; hits <- intersect(ids, members)
      N <- nrow(d); n <- length(ids); K <- length(members); k <- length(hits)
      p <- phyper(k - 1L, K, N - K, n, lower.tail = FALSE)
      # Independent check of the exact test against Fisher's test.
      contingency <- matrix(c(k, n-k, K-k, N-n-K+k), nrow = 2L, byrow = TRUE)
      stopifnot(all(contingency >= 0),
        abs(p - fisher.test(contingency, alternative = 'greater')$p.value) < 1e-10)
      data.frame(rao_category = category, go_id = term,
        description = descriptions$TERM[match(term, descriptions$GOID)],
        n_hits = k, n_category = n, n_background_hits = K, n_background = N,
        fold_enrichment = (k/n)/(K/N), p_value = p,
        gene_ids = paste(hits, collapse = ';'),
        gene_names = paste(d$gene_name[match(hits, d$gene_id)], collapse = ';'))
    }))
  }))
  results$p_adjusted_BH <- p.adjust(results$p_value, method = 'BH')
  # Sensitivity to missing annotation: retain the same category boundaries and
  # test family, but restrict denominators to the BP-annotated genes.
  annotated_N <- sum(d$has_bp_annotation)
  annotated_n <- vapply(results$rao_category, function(category)
    sum(d$rao_category == category & d$has_bp_annotation), integer(1))
  results$n_category_annotated <- annotated_n
  results$n_background_annotated <- annotated_N
  results$fold_enrichment_annotated <- (results$n_hits / annotated_n) /
    (results$n_background_hits / annotated_N)
  results$p_value_annotated <- phyper(results$n_hits - 1L,
    results$n_background_hits, annotated_N - results$n_background_hits,
    annotated_n, lower.tail = FALSE)
  results$p_adjusted_BH_annotated <- p.adjust(results$p_value_annotated, 'BH')
  results <- results[order(results$p_adjusted_BH, results$p_value, results$go_id), ]
  rownames(results) <- NULL
  meta <- AnnotationDbi::metadata(org)
  go_date <- paste(meta$value[grepl('GO.*SOURCEDATE', meta$name) &
    !is.na(meta$value) & nzchar(meta$value)], collapse = ';')
  if (!nzchar(go_date)) go_date <- NA_character_
  summary <- do.call(rbind, lapply(c('low','intermediate','high'), function(category) {
    genes <- d[d$rao_category == category, ]; z <- results[results$rao_category == category, ]
    data.frame(rao_category = category, n_genes = nrow(genes),
      rao_min = min(genes$rao_total), rao_median = median(genes$rao_total), rao_max = max(genes$rao_total),
      n_genes_with_bp_annotation = sum(genes$has_bp_annotation), n_background = nrow(d),
      n_background_with_bp_annotation = sum(d$has_bp_annotation),
      n_terms_tested = nrow(z), n_tests_in_BH_family = nrow(results),
      n_significant_BH_0_05 = sum(z$p_adjusted_BH < .05),
      n_significant_annotated_background_BH_0_05 = sum(z$p_adjusted_BH_annotated < .05),
      min_adjusted_p = min(z$p_adjusted_BH),
      org_Hs_eg_db_version = as.character(utils::packageVersion('org.Hs.eg.db')),
      GO_db_version = as.character(utils::packageVersion('GO.db')), go_source_date = go_date,
      method = 'One-sided exact hypergeometric; BH across all three categories',
      universe = sprintf('All %d scored genes, including genes without BP annotations', nrow(d)),
      term_filter = sprintf('At least 3 and fewer than %d background genes; GOALL BP ancestors', nrow(d)),
      interpretation = 'Exploratory within the preselected gene panel; no genome-wide inference')
  }))
  list(assignments = d[, c('gene_id','gene_name','rao_total','rao_rank','rao_category','has_bp_annotation')],
    results = results, summary = summary)
}
