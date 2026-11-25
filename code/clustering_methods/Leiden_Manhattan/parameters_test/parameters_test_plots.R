# Consensus and non-consensus clustered m6A heatmaps of the parameters_test.R sweep for every PolII
# example gene of tables/parameter_runs.tsv (149 genes; SLC44A4 has no saved clustering): one figure per
# gene x k x resolution x seed replicate x mode. Non-consensus uses seeds 1-5. Replicate s is leiden_consensus_partition(consensus =
# TRUE, tau = 0.5) over Leiden seeds (s - 1) * 100 + 1 ... s * 100, so the replicates of the 5 tested
# seeds use disjoint seed sets.
#
# Per gene the window, read cohort and m6A matrix are rebuilt as in parameters_test.R; the figure
# is plot_example_m6a(clustered = TRUE) of code/Pol2/PolII/pol2_pausing_functions.r, rows
# ordered as in plot_curated_gene() (cluster, then sample, start and read ID).
#
# Outputs under LCL_project/Leiden_manhattan/parameters_test:
#   tables/plotted_genes.tsv        the genes drawn
#   tables/consensus_comparison_runs.tsv  both modes, with actual seed and plot path
#   tables/non_consensus_runs.tsv   non-consensus runs
#   tables/parameter_test_consensus/ per-run read cluster assignments
#   tables/consensus_runs.tsv       one row per gene x k x resolution x replicate: consensus clusters,
#                                   rounds, converged
#   plots/parameter_test_consensus/<group>/<gene>/<gene>_k<k>_res<resolution>_<mode>_seed<s>.pdf
#
# Submit with code/clustering_methods/Leiden_Manhattan/parameters_test/parameters_test_plots.sh

dyn.load("/project/spott/cshan/envs/dimelo/lib/libcurl.so.4", local = FALSE)
source("/project/spott/cshan/fiber-seq/code/Pol2/PolII/pol2_pausing_functions.r")   # also sources
                   # parsing_footprints_functions.r and leiden_manhattan_functions.r
SAVED_DIR <- file.path(PROJECT, "LCL_project/Pol2/PolII/recurrent_patterns/tables")
PARAM_DIR <- file.path(PROJECT, "LCL_project/Leiden_manhattan/parameters_test")
K_VALUES <- c(5L, 10L, 15L, 20L, 30L, 50L)
RESOLUTIONS <- c(0.25, 0.5, 1, 1.5, 2)
SEEDS <- 1:5
CONSENSUS_PARTITIONS <- 100L
CONSENSUS_TAU <- 0.5
WORKERS <- as.integer(Sys.getenv("PARAM_WORKERS", "1"))
data.table::setDTthreads(1L)

# Consensus and non-consensus figures of one gene at one k
#
# 1. read the saved region of 01_run_gene_examples.R and rebuild the full-span m6A matrix
# 2. Manhattan KNN graph; per resolution x replicate: consensus Leiden, clusters renamed by size
# 3. order the rows by cluster, then sample, start and read ID, and draw plot_example_m6a();
#    the parameters go into the figure title through the gene label
# 4. return one row per resolution x replicate x mode
plot_gene_k <- function(group, gene, k, samples) {
  saved <- readRDS(file.path(SAVED_DIR, group, paste0(gene, "_clustering.rds")))
  region <- saved$region
  g <- data.table(group = group, gene_name = gene, chrom = region$chr, tss = region$tss,
                  strand = region$strand, tss_source = region$tss_source)
  dat <- assemble_region_m6a(samples, region, full_span = TRUE,
    positions = seq(region$analysis_start, region$analysis_end))
  meta <- dat$rids_df
  plain <- meta$RID[order(match(meta$sample_name, samples$sample_name), meta$start, meta$RID)]
  feat <- bin_read_matrix(dat$met_mat, region$analysis_start, region$analysis_end, 0)$mat
  consensus_dir <- file.path(PARAM_DIR, "plots/parameter_test_consensus", group, gene)
  knn <- manhattan_knn_graph(feat, k_neighbors = k, verbose = FALSE)
  rows <- list()
  for (res in RESOLUTIONS) {
    for (s in SEEDS) {
      for (use_consensus in c(FALSE, TRUE)) {
        mode <- if (use_consensus) "consensus" else "non_consensus"
        first <- if (use_consensus) (s - 1L) * CONSENSUS_PARTITIONS + 1L else s
        part <- leiden_consensus_partition(knn$graph, resolution = res, seed = first,
          consensus = use_consensus, n_partitions = CONSENSUS_PARTITIONS, tau = CONSENSUS_TAU)
        cluster <- relabel_clusters_by_size(part$membership[rownames(feat)])
        clustered <- rownames(feat)[order(cluster, match(rownames(feat), plain))]
        seed_label <- if (use_consensus)
          sprintf("consensus replicate %d: seeds %d-%d, tau = %g", s, first,
            first + CONSENSUS_PARTITIONS - 1L, CONSENSUS_TAU) else sprintf("non-consensus seed %d", s)
        label <- copy(g)[, gene_name := sprintf("%s (k = %d, resolution = %g, %s)",
          gene, k, res, seed_label)]
        output <- file.path(consensus_dir, sprintf("%s_k%d_res%g_%s_seed%d.pdf", gene, k, res, mode, s))
        plot_example_m6a(label, dat, example_rows(clustered, meta, cluster[match(clustered, rownames(feat))]),
          output, clustered = TRUE)
        sizes <- as.integer(table(cluster))
        assignment_file <- file.path(PARAM_DIR, "tables/parameter_test_consensus", group, gene,
          sprintf("%s_k%d_res%g_%s_seed%d.tsv", gene, k, res, mode, s))
        dir.create(dirname(assignment_file), recursive = TRUE, showWarnings = FALSE)
        fwrite(data.table(RID = rownames(feat), cluster = as.character(cluster)), assignment_file, sep = "\t")
        rows[[length(rows) + 1L]] <- data.table(group = group, gene_name = gene, n_reads = nrow(feat),
          k = k, k_eff = knn$k_eff, resolution = res, mode = mode, seed = s,
          consensus_seed = if (use_consensus) s else NA_integer_, first_leiden_seed = first,
          n_partitions = if (use_consensus) CONSENSUS_PARTITIONS else 1L,
          n_clusters = length(sizes), cluster_sizes = paste(sizes, collapse = ","), quality = part$quality,
          n_rounds = if (use_consensus) part$n_rounds else NA_integer_,
          converged = if (use_consensus) part$converged else NA, plot = output)
      }
    }
  }
  rbindlist(rows)
}

genes <- unique(fread(file.path(PARAM_DIR, "tables/parameter_runs.tsv"))[, .(group, gene_name, n_reads)])
dir.create(file.path(PARAM_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
fwrite(genes, file.path(PARAM_DIR, "tables/plotted_genes.tsv"), sep = "\t")
samples <- as.data.frame(fread(SAMPLES)[, .(sample_name, fire_dir)])
tasks <- CJ(gene_row = seq_len(nrow(genes)), k = K_VALUES)
runs <- parallel::mclapply(seq_len(nrow(tasks)), function(i) {
  gene <- genes[tasks$gene_row[i]]
  message(format(Sys.time()), " task ", i, "/", nrow(tasks), ": ", gene$gene_name, " k = ", tasks$k[i])
  plot_gene_k(gene$group, gene$gene_name, tasks$k[i], samples)
}, mc.cores = WORKERS, mc.preschedule = FALSE)
failed <- vapply(runs, inherits, logical(1), "try-error")
if (any(failed)) stop(sum(failed), " tasks failed\n", as.character(runs[[which(failed)[1]]]))
runs <- rbindlist(runs)
fwrite(runs, file.path(PARAM_DIR, "tables/consensus_comparison_runs.tsv"), sep = "\t")
fwrite(runs[mode == "consensus"], file.path(PARAM_DIR, "tables/consensus_runs.tsv"), sep = "\t")
fwrite(runs[mode == "non_consensus"], file.path(PARAM_DIR, "tables/non_consensus_runs.tsv"), sep = "\t")
