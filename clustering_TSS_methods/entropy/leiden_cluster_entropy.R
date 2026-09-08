#!/usr/bin/env Rscript
# Cluster-level positional Shannon entropy of the Leiden + Manhattan clusters


#  every region's saved result object is loaded from the notebook's
#  output folder, and for each
# Leiden cluster the per-site accessible fraction, the binary positional
# entropy and the mean over sites (H_c) are computed on the full-span filtered
# site matrix used for clustering
#
# Outputs (entropy folder):
#   tables/leiden_cluster_entropy.tsv      region, cluster, n_reads, n_sites, cluster_entropy
#   tables/leiden_positional_entropy.tsv   region, cluster, position, p_accessible, positional_entropy
#   plots/leiden_cluster_entropy_heatmap_<region>.pdf   the Leiden_Manhattan read
#         heatmap with a continuous cluster-entropy track
#   plots/leiden_cluster_entropy_summary.pdf            H_c per region x cluster

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(ComplexHeatmap)
  library(circlize)
})

source("/project/spott/cshan/fiber-seq/code/clustering_TSS_methods/Leiden_Manhattan/leiden_manhattan_plots.r")
source("/project/spott/cshan/fiber-seq/code/clustering_TSS_methods/entropy/entropy_functions.r")

LEIDEN_DIR <- "/project/spott/cshan/fiber-seq/macrophage_project/clustering_TSS_methods/Leiden_Manhattan"

dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

runs <- read.delim(file.path(LEIDEN_DIR, "run_summary_all_regions.tsv"))
cat(nrow(runs), "regions in run_summary_all_regions.tsv\n")

clusters_l <- list(); positions_l <- list()

for (i in seq_len(nrow(runs))) {
  region_id <- runs$region_id[i]
  k_dir <- file.path(LEIDEN_DIR, region_id, runs$bins[i],
                     paste0("k", runs$n_clusters[i]))
  res <- readRDS(file.path(k_dir, paste0("leiden_manhattan_", region_id, ".rds")))

  # full-span filtered site matrix (filter_met_mat()); at window size 0 the
  # feature matrix is that matrix itself
  met_mat <- res$site_met_mat
  if (is.null(met_mat)) {
    stopifnot(res$params$window_size == 0)
    met_mat <- res$feat_mat
  }

  ent <- cluster_positional_entropy(met_mat, res$assignments)
  clusters_l[[region_id]]  <- cbind(region = region_id, ent$clusters)
  positions_l[[region_id]] <- cbind(region = region_id, ent$positions)
  cat("== ", region_id, ": ", nrow(res$assignments), " reads, ", ncol(met_mat),
      " sites, H_c = ", paste(sprintf("%.3f", ent$clusters$cluster_entropy),
                              collapse = " "), "\n", sep = "")

  pdf(file.path(PLOT_DIR, paste0("leiden_cluster_entropy_heatmap_", region_id, ".pdf")),
      width = 8, height = 9)
  ComplexHeatmap::draw(plot_cluster_entropy_heatmap(
    res, ent$clusters,
    main = sprintf("%s (%s, k=%d, res=%g)", region_id, runs$bins[i],
                   res$params$k_neighbors, res$params$resolution)))
  dev.off()
}

cluster_entropy    <- do.call(rbind, clusters_l)
positional_entropy <- do.call(rbind, positions_l)
rownames(cluster_entropy) <- rownames(positional_entropy) <- NULL

stopifnot(all(cluster_entropy$cluster_entropy >= 0 & cluster_entropy$cluster_entropy <= 1),
          all(positional_entropy$p_accessible >= 0 & positional_entropy$p_accessible <= 1))

write.table(cluster_entropy, file.path(TAB_DIR, "leiden_cluster_entropy.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(positional_entropy, file.path(TAB_DIR, "leiden_positional_entropy.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(cluster_entropy)

gg <- plot_cluster_entropy_summary(
  cluster_entropy,
  main = "Cluster entropy of the Leiden + Manhattan clusters (timepoints pooled)")
ggsave(file.path(PLOT_DIR, "leiden_cluster_entropy_summary.pdf"), gg,
       width = 9, height = 4)
cat("wrote leiden_{cluster,positional}_entropy.tsv to", TAB_DIR,
    "and leiden_cluster_entropy_{heatmap_<region>,summary}.pdf to", PLOT_DIR, "\n")
