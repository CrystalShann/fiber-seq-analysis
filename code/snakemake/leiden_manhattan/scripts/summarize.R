source(file.path(snakemake@params[["scripts_dir"]], "common.R"))
load_shared(snakemake@config)
runs <- clusters <- list()
for (path in snakemake@input[["results"]]) {
  result <- readRDS(path)
  a <- result$assignments
  p <- result$params
  row <- data.frame(dataset = snakemake@wildcards[["dataset"]], region_id = result$region$region_id,
    window_size = p$window_size, k_neighbors = p$k_neighbors, resolution = p$resolution,
    seed = p$seed, n_clusters = result$n_clusters, n_reads = nrow(a),
    n_features = p$n_features, k_eff = p$k_eff, sigma = p$sigma,
    mean_knn_dist = p$mean_knn_dist, clustering_path = path)
  for (key in intersect(c("rank", "focal_snp", "fisher_pvalue", "fisher_qvalue"), names(result$region)))
    row[[key]] <- result$region[[key]]
  runs[[path]] <- row
  count <- as.data.frame(table(cluster = a$cluster), stringsAsFactors = FALSE)
  names(count)[2L] <- "n_reads"
  base <- row[rep(1L, nrow(count)), c("dataset", "region_id", "window_size", "k_neighbors", "resolution", "seed", "n_clusters"), drop = FALSE]
  clusters[[path]] <- data.frame(base, count, fraction = count$n_reads / nrow(a), row.names = NULL)
}
write_tsv(dplyr::bind_rows(runs), snakemake@output[["run_summary"]])
write_tsv(dplyr::bind_rows(clusters), snakemake@output[["cluster_summary"]])
