# Column-permutation null for the Leiden / Manhattan clusters of the PolII example genes
# (code/pausing/PolII/example_gene_list.tsv).
#
# Per gene, the TSS +/- 1 kb window, read cohort and clustering settings are the saved ones of
# code/pausing/PolII/01_run_gene_examples.R (LCL_project/pausing/PolII/recurrent_patterns/tables/
# <group>/<gene>_clustering.rds). Shuffling each column (base) of the reads x bases m6A matrix
# across reads keeps every base's m6A fraction but breaks which calls sit on the same molecule.
# If the clusters reflect molecule-level m6A patterns, the permuted matrices should give a weaker
# partition (lower modularity) and more distant nearest neighbours (higher mean KNN Manhattan
# distance) than the real one. 10,000 permutations per gene, run in parallel within each gene.
#
# Outputs under LCL_project/Leiden_manhattan/permute_clusters:
#   tables/permutations/<group>/<gene>.tsv  observed row (permutation = 0) + one row per permutation
#   tables/permutation_runs.tsv             all genes' rows of the files above
#   tables/permutation_summary.tsv          one row per gene: observed vs permuted, empirical p
#   plots/modularity_observed_vs_permuted.pdf, plots/n_clusters_observed_vs_permuted.pdf,
#   plots/mean_knn_dist_observed_vs_permuted.pdf
# A gene whose tables/permutations file already exists is read back instead of rerun, so a job
# that hits the time limit can be resubmitted with the same command.
#
# sbatch --account=pi-spott --partition=caslake --cpus-per-task=48 --mem=128G --time=36:00:00 \
#   --job-name=permute_clusters \
#   --output=/project/spott/cshan/fiber-seq/LCL_project/Leiden_manhattan/permute_clusters/run.log \
#   --export=ALL,PERMUTE_WORKERS=48 /project/spott/cshan/fiber-seq/code/pausing/PolII/run_gene_examples.sh \
#   /project/spott/cshan/fiber-seq/code/clustering_methods/Leiden_Manhattan/permute_clusters.R

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
PROJECT <- "/project/spott/cshan/fiber-seq"
source(file.path(PROJECT, "code/parsing_functions/parsing_footprints_functions.r"))   # assemble_region_m6a()
source(file.path(PROJECT, "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r"))   # leiden_manhattan_cluster()

EXAMPLE_GENES <- file.path(PROJECT, "code/pausing/PolII/example_gene_list.tsv")
SAVED_DIR <- file.path(PROJECT, "LCL_project/pausing/PolII/recurrent_patterns/tables")
SAMPLES <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv"
OUT_DIR <- file.path(PROJECT, "LCL_project/Leiden_manhattan/permute_clusters")
N_PERM <- 10000L
WORKERS <- as.integer(Sys.getenv("PERMUTE_WORKERS", "1"))
data.table::setDTthreads(1L)

# Partition statistics of one leiden_manhattan_cluster() result
#
# 1. modularity of the Leiden partition on the weighted KNN graph
# 2. number of clusters, fraction of reads in the largest one and mean KNN Manhattan distance
cluster_stats <- function(result) {
  memb <- as.integer(result$assignments$cluster)[match(igraph::V(result$graph)$name, result$assignments$RID)]
  data.table(n_clusters = result$n_clusters,
             modularity = igraph::modularity(result$graph, memb, weights = igraph::E(result$graph)$weight),
             largest_cluster_fraction = max(table(result$assignments$cluster)) / nrow(result$assignments),
             mean_knn_dist = result$params$mean_knn_dist)
}

# Observed and N_PERM column-permuted clusterings of one gene
#
# 1. read the saved region and clustering of 01_run_gene_examples.R
# 2. rebuild the full-span m6A matrix over every base of the window with assemble_region_m6a()
# 3. re-cluster with the saved settings; flag whether the assignments match the saved ones
# 4. per permutation p (in parallel over WORKERS): set.seed(p), shuffle every column across
#    reads, cluster again
# 5. write and return one row per clustering (permutation 0 = observed)
permute_gene <- function(group, gene, samples) {
  out_file <- file.path(OUT_DIR, "tables/permutations", group, paste0(gene, ".tsv"))
  if (file.exists(out_file)) return(fread(out_file))
  saved <- readRDS(file.path(SAVED_DIR, group, paste0(gene, "_clustering.rds")))
  region <- saved$region
  p <- saved$params
  dat <- assemble_region_m6a(samples, region, full_span = TRUE,
    positions = seq(region$analysis_start, region$analysis_end))
  run <- function(mat) leiden_manhattan_cluster(mat, dat$rids_df, region$analysis_start, region$analysis_end,
    window_size = p$window_size, k_neighbors = p$k_neighbors, resolution = p$resolution,
    seed = p$seed, verbose = FALSE)
  observed <- run(dat$met_mat)
  same <- identical(observed$assignments$RID, saved$assignments$RID) &&
    identical(as.character(observed$assignments$cluster), as.character(saved$assignments$cluster))
  mat <- as.matrix(dat$met_mat)
  permuted <- parallel::mclapply(seq_len(N_PERM), function(i) {
    set.seed(i)
    shuffled <- apply(mat, 2, sample)
    dimnames(shuffled) <- dimnames(mat)
    cbind(permutation = i, cluster_stats(run(shuffled)))
  }, mc.cores = WORKERS)
  failed <- vapply(permuted, inherits, logical(1), "try-error")
  if (any(failed)) stop(gene, ": ", sum(failed), " permutations failed\n", as.character(permuted[[which(failed)[1]]]))
  out <- cbind(group = group, gene_name = gene, n_reads = nrow(mat), matches_saved_clusters = same,
               rbindlist(c(list(cbind(permutation = 0L, cluster_stats(observed))), permuted)))
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  fwrite(out, out_file, sep = "\t")
  out
}

genes <- fread(EXAMPLE_GENES)
genes <- genes[file.exists(file.path(SAVED_DIR, group, paste0(gene_name, "_clustering.rds")))]
samples <- as.data.frame(fread(SAMPLES)[, .(sample_name, fire_dir)])
runs <- rbindlist(lapply(seq_len(nrow(genes)), function(i) {
  message(format(Sys.time()), " gene ", i, "/", nrow(genes), ": ", genes$gene_name[i])
  permute_gene(genes$group[i], genes$gene_name[i], samples)
}))
fwrite(runs, file.path(OUT_DIR, "tables/permutation_runs.tsv"), sep = "\t")

# Observed vs permuted per gene; p = (1 + permutations at least as structured) / (1 + N_PERM)
summary <- runs[, {
  obs <- .SD[permutation == 0L]; null <- .SD[permutation > 0L]
  .(observed_n_clusters = obs$n_clusters, permuted_n_clusters_mean = mean(null$n_clusters),
    observed_modularity = obs$modularity, permuted_modularity_mean = mean(null$modularity),
    permuted_modularity_sd = sd(null$modularity),
    z_modularity = (obs$modularity - mean(null$modularity)) / sd(null$modularity),
    p_modularity = (1 + sum(null$modularity >= obs$modularity)) / (1 + nrow(null)),
    observed_mean_knn_dist = obs$mean_knn_dist, permuted_mean_knn_dist_mean = mean(null$mean_knn_dist),
    permuted_mean_knn_dist_sd = sd(null$mean_knn_dist),
    z_mean_knn_dist = (obs$mean_knn_dist - mean(null$mean_knn_dist)) / sd(null$mean_knn_dist),
    p_mean_knn_dist = (1 + sum(null$mean_knn_dist <= obs$mean_knn_dist)) / (1 + nrow(null)))
}, by = .(group, gene_name, n_reads, matches_saved_clusters)]
fwrite(summary, file.path(OUT_DIR, "tables/permutation_summary.tsv"), sep = "\t")

dir.create(file.path(OUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
runs[, gene_name := factor(gene_name, levels = summary[order(observed_modularity), gene_name])]
plot_metric <- function(metric, label, file) {
  fig <- ggplot(runs[permutation > 0L], aes(.data[[metric]], gene_name)) +
    geom_boxplot(colour = "grey50", outlier.size = 0.3, linewidth = 0.3) +
    geom_point(data = runs[permutation == 0L], colour = "#D55E00", size = 1.2) +
    facet_wrap(~ group, scales = "free_y") +
    labs(x = label, y = NULL,
         title = sprintf("Observed (orange) vs %d column-permuted m6A matrices (grey)", N_PERM)) +
    theme_bw(base_size = 8)
  ggsave(file.path(OUT_DIR, "plots", file), fig, width = 12, height = 9)
}
plot_metric("modularity", "Leiden partition modularity", "modularity_observed_vs_permuted.pdf")
plot_metric("n_clusters", "number of Leiden clusters", "n_clusters_observed_vs_permuted.pdf")
plot_metric("mean_knn_dist", "mean KNN Manhattan distance", "mean_knn_dist_observed_vs_permuted.pdf")
