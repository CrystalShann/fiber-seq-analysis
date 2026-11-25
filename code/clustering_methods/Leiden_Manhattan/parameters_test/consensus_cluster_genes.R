# Override settings with CONSENSUS_K and CONSENSUS_RESOLUTION; non-default runs use
# plots/consensus_cluster_k<K>_n<R> and tables/consensus_cluster_k<K>_n<R>.
# Consensus Leiden / Manhattan clusters of every PolII example gene (code/Pol2/PolII/example_gene_list.tsv)
# at the settings of code/Pol2/PolII/01_run_gene_examples.R (k = 10, resolution = 1,
# window_size = 0): leiden_consensus_partition(consensus = TRUE), seeds 1-100, tau = 0.5.
#
# Per gene the window and read cohort are the saved ones of 01_run_gene_examples.R
# (LCL_project/Pol2/PolII/recurrent_patterns/tables/<group>/<gene>_clustering.rds; genes without
# one, e.g. SLC44A4 with too few full-span reads, are skipped). The figure is the clustered m6A heatmap
# of plot_curated_gene() (code/Pol2/PolII/pol2_pausing_functions.r) with the consensus clusters,
# rows ordered by cluster, then sample, start and read ID.
#
# Outputs under LCL_project/Leiden_manhattan/parameters_test:
#   plots/consensus_cluster_k10_n1_original/<group>/<gene>/<gene>_consensus_clustered_m6a.pdf
#   tables/consensus_cluster/<group>/<gene>_read_clusters.tsv   RID, cluster, sample, coordinates
#   tables/consensus_cluster/consensus_summary.tsv              one row per gene
#
# Submit with code/clustering_methods/Leiden_Manhattan/parameters_test/consensus_cluster_genes.sh

dyn.load("/project/spott/cshan/envs/dimelo/lib/libcurl.so.4", local = FALSE)
source("/project/spott/cshan/fiber-seq/code/Pol2/PolII/pol2_pausing_functions.r")   # also sources
                   # parsing_footprints_functions.r and leiden_manhattan_functions.r
SAVED_DIR <- file.path(PROJECT, "LCL_project/Pol2/PolII/recurrent_patterns/tables")
PARAM_DIR <- file.path(PROJECT, "LCL_project/Leiden_manhattan/parameters_test")
K_NEIGHBORS <- as.integer(Sys.getenv("CONSENSUS_K", "10"))
RESOLUTION <- as.numeric(Sys.getenv("CONSENSUS_RESOLUTION", "1"))
stopifnot(is.finite(K_NEIGHBORS), K_NEIGHBORS >= 1L, is.finite(RESOLUTION), RESOLUTION > 0)
RUN_NAME <- sprintf("consensus_cluster_k%d_n%g", K_NEIGHBORS, RESOLUTION)
PLOT_DIR <- file.path(PARAM_DIR, "plots", if (K_NEIGHBORS == 10L && RESOLUTION == 1)
  "consensus_cluster_k10_n1_original" else RUN_NAME)
TABLE_DIR <- file.path(PARAM_DIR, "tables", if (K_NEIGHBORS == 10L && RESOLUTION == 1)
  "consensus_cluster" else RUN_NAME)
CONSENSUS_PARTITIONS <- 100L
CONSENSUS_TAU <- 0.5
WORKERS <- as.integer(Sys.getenv("PARAM_WORKERS", "1"))
data.table::setDTthreads(1L)

#########
# Consensus Leiden clusters and clustered m6A heatmap of one gene
#######
# 1. Read the saved region and rebuild the full-span m6A matrix
# 2. Consensus Leiden clusters with leiden_manhattan_consensus_cluster() (100 seeds, tau = 0.5)
# 3. Order rows by cluster, then sample, start and read ID
# 4. Write the per-read cluster table
# 5. Draw the clustered m6A heatmap, settings in the figure title
# 6. Return the gene's summary row
consensus_gene <- function(group, gene, samples) {
  saved <- readRDS(file.path(SAVED_DIR, group, paste0(gene, "_clustering.rds")))
  region <- saved$region
  g <- data.table(group = group, gene_name = gene, chrom = region$chr, tss = region$tss,
                  strand = region$strand, tss_source = region$tss_source)
  dat <- assemble_region_m6a(samples, region, full_span = TRUE,
    positions = seq(region$analysis_start, region$analysis_end))
  meta <- dat$rids_df
  plain <- meta$RID[order(match(meta$sample_name, samples$sample_name), meta$start, meta$RID)]
  res <- leiden_manhattan_consensus_cluster(dat$met_mat, meta, region$analysis_start, region$analysis_end,
    k_neighbors = K_NEIGHBORS, resolution = RESOLUTION, n_partitions = CONSENSUS_PARTITIONS,
    tau = CONSENSUS_TAU, verbose = FALSE)
  a <- res$assignments
  clustered <- a$RID[order(a$cluster, match(a$RID, plain))]
  rows <- example_rows(clustered, meta, a$cluster[match(clustered, a$RID)])

  a$clustered_row <- match(a$RID, clustered)
  tables <- file.path(TABLE_DIR, group)
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)
  fwrite(a, file.path(tables, paste0(gene, "_read_clusters.tsv")), sep = "\t")

  label <- copy(g)[, gene_name := sprintf("%s (consensus of seeds 1-%d, tau = %g, k = %d, resolution = %g)",
                                         gene, CONSENSUS_PARTITIONS, CONSENSUS_TAU, K_NEIGHBORS, RESOLUTION)]
  path <- function(suffix) file.path(PLOT_DIR, group, gene,
                                     paste0(gene, "_", suffix, ".pdf"))
  plot_example_m6a(label, dat, rows, path("consensus_clustered_m6a"), clustered = TRUE)

  sizes <- as.integer(table(a$cluster))
  p <- res$params
  data.table(group = group, gene_name = gene, n_reads = p$n_reads, k = K_NEIGHBORS, resolution = RESOLUTION,
             n_partitions = CONSENSUS_PARTITIONS, tau = CONSENSUS_TAU, k_eff = p$k_eff,
             n_clusters = length(sizes), cluster_sizes = paste(sizes, collapse = ","),
             quality = p$quality, n_rounds = p$n_rounds, converged = p$converged)
}

genes <- curated_example_genes()
genes <- genes[file.exists(file.path(SAVED_DIR, group, paste0(gene_name, "_clustering.rds")))]
samples <- as.data.frame(fread(SAMPLES)[, .(sample_name, fire_dir)])
runs <- parallel::mclapply(seq_len(nrow(genes)), function(i) {
  message(format(Sys.time()), " gene ", i, "/", nrow(genes), ": ", genes$gene_name[i])
  consensus_gene(genes$group[i], genes$gene_name[i], samples)
}, mc.cores = WORKERS, mc.preschedule = FALSE)
failed <- vapply(runs, inherits, logical(1), "try-error")
if (any(failed)) stop(sum(failed), " genes failed\n", as.character(runs[[which(failed)[1]]]))
fwrite(rbindlist(runs), file.path(TABLE_DIR, "consensus_summary.tsv"), sep = "\t")
