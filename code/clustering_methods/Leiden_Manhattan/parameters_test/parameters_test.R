# Seed / k / resolution sweep of the Leiden / Manhattan clusters of the PolII example genes
# (code/Pol2/PolII/example_gene_list.tsv): how do the cluster number and sizes change?
#
# Per gene, the TSS +/- 1 kb window and read cohort are the saved ones of
# code/Pol2/PolII/01_run_gene_examples.R (LCL_project/Pol2/PolII/recurrent_patterns/tables/
# <group>/<gene>_clustering.rds; saved settings k = 10, resolution = 1, seed = 1, window_size = 0).
# The full-span m6A matrix is rebuilt once per gene, the Manhattan KNN graph once per k, and Leiden
# is rerun for every resolution x seed.
#
# Outputs under LCL_project/Leiden_manhattan/parameters_test:
#   tables/per_gene/<group>/<gene>.tsv  one row per k x resolution x seed
#   tables/parameter_runs.tsv           all genes' rows of the files above
#   tables/parameter_summary.tsv        one row per k x resolution x seed: medians over genes
#   plots/n_clusters_heatmap.pdf
# A gene whose tables/per_gene file already exists is read back instead of rerun.
#
# Submit with code/clustering_methods/Leiden_Manhattan/parameters_test/parameters_test.sh

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
PROJECT <- "/project/spott/cshan/fiber-seq"
source(file.path(PROJECT, "code/parsing_functions/parsing_footprints_functions.r"))   # assemble_region_m6a()
source(file.path(PROJECT, "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r"))

EXAMPLE_GENES <- file.path(PROJECT, "code/Pol2/PolII/example_gene_list.tsv")
SAVED_DIR <- file.path(PROJECT, "LCL_project/Pol2/PolII/recurrent_patterns/tables")
SAMPLES <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv"
OUT_DIR <- file.path(PROJECT, "LCL_project/Leiden_manhattan/parameters_test")
K_VALUES <- c(5L, 10L, 15L, 20L, 30L, 50L)
RESOLUTIONS <- c(0.25, 0.5, 1, 1.5, 2)
SEEDS <- 1:5
WORKERS <- as.integer(Sys.getenv("PARAM_WORKERS", "1"))
data.table::setDTthreads(1L)

# All k x resolution x seed clusterings of one gene
#
# 1. read the saved region and clustering of 01_run_gene_examples.R
# 2. rebuild the full-span m6A matrix over every base of the window with assemble_region_m6a()
# 3. per k: Manhattan KNN graph; per resolution x seed: Leiden, clusters renamed by size
# 4. cluster-size statistics; check that k = 10, resolution = 1, seed = 1 reproduces the saved
#    clusters (same grouping of reads, whatever the labels)
# 5. write and return one row per clustering
sweep_gene <- function(group, gene, samples) {
  out_file <- file.path(OUT_DIR, "tables/per_gene", group, paste0(gene, ".tsv"))
  if (file.exists(out_file)) return(fread(out_file))
  saved <- readRDS(file.path(SAVED_DIR, group, paste0(gene, "_clustering.rds")))
  region <- saved$region
  stopifnot(saved$params$window_size == 0L)
  dat <- assemble_region_m6a(samples, region, full_span = TRUE,
    positions = seq(region$analysis_start, region$analysis_end))
  feat <- bin_read_matrix(dat$met_mat, region$analysis_start, region$analysis_end, 0)$mat
  saved_memb <- as.integer(saved$assignments$cluster[match(rownames(feat), saved$assignments$RID)])
  rows <- list()
  matches_saved <- FALSE
  for (k in K_VALUES) {
    knn <- manhattan_knn_graph(feat, k_neighbors = k, verbose = FALSE)
    for (res in RESOLUTIONS) {
      for (seed in SEEDS) {
        part <- leiden_partition(knn$graph, resolution = res, seed = seed)
        cluster <- relabel_clusters_by_size(part$membership[rownames(feat)])
        if (k == 10L && res == 1 && seed == 1L && !anyNA(saved_memb))
          matches_saved <- nrow(unique(data.frame(cluster, saved_memb))) == uniqueN(cluster) &&
            uniqueN(cluster) == uniqueN(saved_memb)
        sizes <- as.integer(table(cluster))
        rows[[length(rows) + 1L]] <- data.table(k = k, k_eff = knn$k_eff, resolution = res, seed = seed,
          n_clusters = length(sizes), largest_cluster_fraction = max(sizes) / sum(sizes),
          smallest_cluster_fraction = min(sizes) / sum(sizes),
          cluster_sizes = paste(sizes, collapse = ","), quality = part$quality)
      }
    }
  }
  out <- cbind(group = group, gene_name = gene, n_reads = nrow(feat), rbindlist(rows))
  out[, matches_saved_clusters := matches_saved]
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  fwrite(out, out_file, sep = "\t")
  out
}

genes <- fread(EXAMPLE_GENES)
genes <- genes[file.exists(file.path(SAVED_DIR, group, paste0(gene_name, "_clustering.rds")))]
samples <- as.data.frame(fread(SAMPLES)[, .(sample_name, fire_dir)])
runs <- parallel::mclapply(seq_len(nrow(genes)), function(i) {
  message(format(Sys.time()), " gene ", i, "/", nrow(genes), ": ", genes$gene_name[i])
  sweep_gene(genes$group[i], genes$gene_name[i], samples)
}, mc.cores = WORKERS)
failed <- vapply(runs, inherits, logical(1), "try-error")
if (any(failed)) stop(sum(failed), " genes failed\n", as.character(runs[[which(failed)[1]]]))
runs <- rbindlist(runs)
fwrite(runs, file.path(OUT_DIR, "tables/parameter_runs.tsv"), sep = "\t")
message(runs[, sum(!matches_saved_clusters[1]), by = gene_name][, sum(V1)],
        " genes do not reproduce the saved clusters at k = 10, resolution = 1, seed = 1")

summary <- runs[, .(n_genes = .N, median_n_clusters = median(n_clusters),
  single_cluster_genes = sum(n_clusters == 1L),
  median_largest_cluster_fraction = median(largest_cluster_fraction),
  median_smallest_cluster_fraction = median(smallest_cluster_fraction)),
  by = .(k, resolution, seed)]
fwrite(summary, file.path(OUT_DIR, "tables/parameter_summary.tsv"), sep = "\t")

# Heatmap of k x resolution, one panel per seed (medians over genes)
dir.create(file.path(OUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
plot_heatmap <- function(metric, label, file) {
  fig <- ggplot(summary, aes(factor(resolution), factor(k), fill = .data[[metric]])) +
    geom_tile(colour = "white") +
    geom_text(aes(label = signif(.data[[metric]], 2)), size = 2.5) +
    scale_fill_gradient(low = "#eff3ff", high = "#3182bd") +
    facet_wrap(~ seed, nrow = 1, labeller = label_both) +
    labs(x = "Leiden resolution", y = "KNN neighbours (k)", fill = label,
         title = sprintf("%s over %d PolII example genes", label, uniqueN(runs$gene_name))) +
    theme_minimal(base_size = 9) + theme(panel.grid = element_blank())
  ggsave(file.path(OUT_DIR, "plots", file), fig, width = 14, height = 3.5)
}
plot_heatmap("median_n_clusters", "median number of clusters", "n_clusters_heatmap.pdf")
