# Column-permutation null for the Leiden / Manhattan clusters of the PolII example genes
# (code/Pol2/PolII/example_gene_list.tsv).
#
# Per gene, the TSS +/- 1 kb window, read cohort and clustering settings are the saved ones of
# code/Pol2/PolII/01_run_gene_examples.R (LCL_project/Pol2/PolII/recurrent_patterns/tables/
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
# A gene whose tables/permutations file already exists is read back instead of rerun, so a job
# that hits the time limit can be resubmitted with the same command.
#
# The last stage draws, per gene, the observed and permutation-1 m6A heatmaps side by side (each
# matrix clustered on its own) and checks their statistics against the permutation tables.
# Outputs under LCL_project/Leiden_manhattan/parameters_test:
#   plots/permutation/<gene>.pdf             observed | permutation 1 heatmaps
#   tables/permutation/<gene>.rds            both matrices, cluster assignments and statistics
#   tables/permutation/genes.tsv, runs.tsv   plotted genes; one row per gene x panel
# PERMUTE_PLOT_GENE=<gene> plots one gene only (<gene>_genes.tsv, <gene>_runs.tsv);
# PERMUTE_PLOT_WORKERS sets the plotting cores.
#
# sbatch --account=pi-spott --partition=caslake --cpus-per-task=48 --mem=128G --time=36:00:00 \
#   --job-name=permute_clusters \
#   --output=/project/spott/cshan/fiber-seq/LCL_project/Leiden_manhattan/permute_clusters/run.log \
#   --export=ALL,PERMUTE_WORKERS=48 /project/spott/cshan/fiber-seq/code/Pol2/PolII/run_gene_examples.sh \
#   /project/spott/cshan/fiber-seq/code/clustering_methods/Leiden_Manhattan/permutation/permute_clusters.R
# Once every gene's permutation table exists, permute_clusters_side_by_side.sh reruns the
# script cheaply (tables are read back) to redraw the heatmaps.

dyn.load("/project/spott/cshan/envs/dimelo/lib/libcurl.so.4", local = FALSE)
# Defines PROJECT, EXAMPLE_GENES, SAMPLES and the heatmap helpers; also sources
# parsing_footprints_functions.r (assemble_region_m6a()) and leiden_manhattan_functions.r
source("/project/spott/cshan/fiber-seq/code/Pol2/PolII/pol2_pausing_functions.r")

SAVED_DIR <- file.path(PROJECT, "LCL_project/Pol2/PolII/recurrent_patterns/tables")
OUT_DIR <- file.path(PROJECT, "LCL_project/Leiden_manhattan/permute_clusters")
PLOT_DIR <- file.path(PROJECT, "LCL_project/Leiden_manhattan/parameters_test/plots/permutation")
PLOT_TABLE_DIR <- file.path(PROJECT, "LCL_project/Leiden_manhattan/parameters_test/tables/permutation")
N_PERM <- 10000L
WORKERS <- as.integer(Sys.getenv("PERMUTE_WORKERS", "1"))
data.table::setDTthreads(1L)

#########
# Partition statistics of one leiden_manhattan_cluster() result
#######
# 1. Modularity of the Leiden partition on the weighted KNN graph
# 2. Number of clusters, fraction of reads in the largest one and mean KNN Manhattan distance
cluster_stats <- function(result) {
  memb <- as.integer(result$assignments$cluster)[match(igraph::V(result$graph)$name, result$assignments$RID)]
  data.table(n_clusters = result$n_clusters,
             modularity = igraph::modularity(result$graph, memb, weights = igraph::E(result$graph)$weight),
             largest_cluster_fraction = max(table(result$assignments$cluster)) / nrow(result$assignments),
             mean_knn_dist = result$params$mean_knn_dist)
}

#########
# Observed and N_PERM column-permuted clusterings of one gene
#######
# 1. Read the saved region and clustering of 01_run_gene_examples.R (cached table: read it back)
# 2. Rebuild the full-span m6A matrix over every base of the window with assemble_region_m6a()
# 3. Re-cluster with the saved settings; flag whether the assignments match the saved ones
# 4. Per permutation p (in parallel over WORKERS): set.seed(p), shuffle every column across
#    reads, cluster again
# 5. Write and return one row per clustering (permutation 0 = observed)
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

# Side-by-side observed / permutation 1 heatmaps of every gene in the summary
plot_genes <- summary[, .(group, gene_name)]
test_gene <- Sys.getenv("PERMUTE_PLOT_GENE", "")
if (nzchar(test_gene)) {
  plot_genes <- plot_genes[gene_name == test_gene]
  stopifnot(nrow(plot_genes) == 1L)
}
ComplexHeatmap::ht_opt(message = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(!anyDuplicated(plot_genes$gene_name))
fwrite(plot_genes, file.path(PLOT_TABLE_DIR, if (nzchar(test_gene)) paste0(test_gene, "_genes.tsv") else "genes.tsv"), sep = "\t")

#########
# Observed and permutation-1 heatmaps of one gene, side by side
#######
# 1. Rebuild the gene's full-span m6A matrix from its saved region
# 2. Cluster the observed matrix and the permutation-1 matrix (set.seed(1), columns shuffled)
# 3. Check clusters, modularity and mean KNN distance against the permutation table
# 4. Draw both clustered heatmaps in one PDF
# 5. Save both matrices and assignments; return one statistics row per panel
plot_gene <- function(group, gene) {
  message(format(Sys.time()), " plotting ", gene)
  saved <- readRDS(file.path(SAVED_DIR, group, paste0(gene, "_clustering.rds")))
  region <- saved$region
  p <- saved$params
  dat <- assemble_region_m6a(samples, region, full_span = TRUE,
    positions = seq(region$analysis_start, region$analysis_end))
  mat <- as.matrix(dat$met_mat)
  meta <- dat$rids_df
  plain <- meta$RID[order(match(meta$sample_name, samples$sample_name), meta$start, meta$RID)]
  reference <- runs[gene_name == gene]
  panels <- list()
  matrices <- list()
  assignments <- list()
  records <- lapply(0:1, function(i) {
    if (i == 0L) {
      current <- mat
    } else {
      set.seed(i)
      current <- apply(mat, 2, sample)
      dimnames(current) <- dimnames(mat)
    }
    stopifnot(identical(colSums(current), colSums(mat)))
    result <- leiden_manhattan_cluster(current, meta, region$analysis_start, region$analysis_end,
      window_size = p$window_size, k_neighbors = p$k_neighbors,
      resolution = p$resolution, seed = p$seed, verbose = FALSE)
    a <- result$assignments
    membership <- as.integer(a$cluster)[match(igraph::V(result$graph)$name, a$RID)]
    modularity <- igraph::modularity(result$graph, membership, weights = igraph::E(result$graph)$weight)
    expected <- reference[permutation == i]
    stopifnot(nrow(expected) == 1L, result$n_clusters == expected$n_clusters,
      isTRUE(all.equal(modularity, expected$modularity, tolerance = 1e-8)),
      isTRUE(all.equal(result$params$mean_knn_dist, expected$mean_knn_dist, tolerance = 1e-8)))
    ordered <- a$RID[order(a$cluster, match(a$RID, plain))]
    rows <- example_rows(ordered, meta, a$cluster[match(ordered, a$RID)])
    label <- if (i == 0L) "observed" else sprintf("permutation %d", i)
    g <- data.table(group = group, gene_name = gene, chrom = region$chr,
      tss = region$tss, strand = region$strand, tss_source = region$tss_source)
    current_sparse <- Matrix::Matrix(current, sparse = TRUE)
    title <- sprintf("%s | %s\n%s | %d reads, %d clusters | Q = %.3f\nk = %d, resolution = %g | TSS-relative bp",
      gene, group, label, nrow(mat), result$n_clusters, modularity, p$k_neighbors, p$resolution)
    panels[[i + 1L]] <<- plot_read_heatmap(
      list(assignments = rows, site_met_mat = current_sparse), example_region(g),
      example_sample_colors(rows$sample_name), show_cluster_profiles = TRUE,
      heatmap_options = list(column_title = title, column_title_gp = grid::gpar(fontsize = 10),
        left_annotation = ComplexHeatmap::rowAnnotation(cluster = rows$cluster,
          col = list(cluster = cluster_id_palette(levels(rows$cluster))))))
    matrices[[i + 1L]] <<- current_sparse
    assignments[[i + 1L]] <<- a
    output <- file.path(PLOT_DIR, paste0(gene, ".pdf"))
    sizes <- table(a$cluster)
    data.table(group = group, gene_name = gene, permutation = i, n_reads = nrow(mat),
      n_clusters = result$n_clusters, cluster_sizes = paste(as.integer(sizes), collapse = ","),
      largest_cluster_fraction = max(sizes) / nrow(mat), modularity = modularity,
      mean_knn_dist = result$params$mean_knn_dist, matches_saved_statistics = TRUE, plot = output)
  })
  records <- rbindlist(records)
  save_figure(panels, records$plot[1], width = 20,
    height = 6 + 0.55 * max(records$n_clusters), draw = function(x) {
      grid::grid.newpage()
      for (j in 1:2) {
        grid::pushViewport(grid::viewport(x = (j - 0.5) / 2, width = 0.5))
        ComplexHeatmap::draw(x[[j]], newpage = FALSE, merge_legend = TRUE,
          padding = grid::unit(c(5, 5, 5, 22), "mm"))
        grid::upViewport(0)
      }
    })
  artifact <- file.path(PLOT_TABLE_DIR, paste0(gene, ".rds"))
  dir.create(dirname(artifact), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(region = region, params = p, permutation = 1L,
    matrices = setNames(matrices, c("observed", "permuted")),
    assignments = setNames(assignments, c("observed", "permuted")), statistics = records), artifact)
  records
}
plot_workers <- as.integer(Sys.getenv("PERMUTE_PLOT_WORKERS", "1"))
plot_runs <- parallel::mclapply(seq_len(nrow(plot_genes)), function(i)
  plot_gene(plot_genes$group[i], plot_genes$gene_name[i]), mc.cores = plot_workers, mc.preschedule = FALSE)
failed <- vapply(plot_runs, inherits, logical(1), "try-error")
if (any(failed)) stop(paste(unlist(plot_runs[failed]), collapse = "\n"))
fwrite(rbindlist(plot_runs), file.path(PLOT_TABLE_DIR, if (nzchar(test_gene)) paste0(test_gene, "_runs.tsv") else "runs.tsv"), sep = "\t")
