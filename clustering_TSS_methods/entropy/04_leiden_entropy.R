#!/usr/bin/env Rscript
# Read-level Shannon entropy within Leiden clusters at the 8 early-response
# promoters (TSS +/- 500), reads POOLED across the four LPS timepoints before
# clustering so Leiden clusters represent the same molecular states across
# time.
#
# Clustering (reusing the Leiden_Manhattan pipeline):
#   reads -> 10 bp binned accessibility profiles -> Manhattan kNN graph ->
#   igraph Leiden -> assignments split back by timepoint.
#
# Entropy: per READ, from the raw calls over the region's called-position
# universe (get_sparse_met_mat's default; assigned reads are full-span so
# they cover every site). With p = n_methylated / n_sites_covered per read:
#
#   H_read = -p log2 p - (1 - p) log2 (1 - p)     (0 log 0 = 0; H in [0,1])
#
# H_read ~ 0 = the read is uniformly closed (or uniformly methylated) along
# the window; H_read ~ 1 = the read mixes accessible and inaccessible sites
# about half and half. Read-level entropy is then compared (Kruskal-Wallis):
#   (a) among Leiden clusters within each region (timepoints pooled);
#   (b) across 0/5/10/15 min within each region x cluster.
# Read-level heatmaps order reads timepoint -> cluster -> H_read and annotate
# every read with H_read and p_accessible.

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(GenomicRanges)
  library(data.table)
  library(Rsamtools)
  library(igraph)
  library(cowplot)
  library(ComplexHeatmap)
  library(circlize)
})

source("/project/spott/cshan/fiber-seq/code/topic_model/topic_modelling_functions.r")
source("/project/spott/cshan/fiber-seq/code/clustering_TSS_methods/Leiden_Manhattan/leiden_manhattan_functions.r")
source("/project/spott/cshan/fiber-seq/code/clustering_TSS_methods/Leiden_Manhattan/leiden_manhattan_plots.r")
source("/project/spott/cshan/fiber-seq/code/clustering_TSS_methods/entropy/entropy_functions.r")

options(scipen = 999)
set.seed(123)

HALF <- 500L            # cluster on the 1 kb entropy window (TSS +/- 500)
WINDOW_SIZE <- 10L      # bin size of the CLUSTERING features (entropy is unbinned)
K_NEIGHBORS <- 10L      # as leiden_manhattan.Rmd (paper's 50 over-connects here)
RESOLUTION <- 1
LEIDEN_SEED <- 1

FT_RESULT_DIR <- "/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/ft_result_dir"
PROMOTER_TABLE <- "/project/spott/cshan/fiber-seq/macrophage_project/clustering_TSS_methods/hamming_hierarchical/promoter_regions_gencode_v49.tsv"

dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

sample_table <- data.table(sample_name = SAMPLES, time = c(0, 5, 10, 15))

binary_entropy <- function(p) {
  ifelse(is.na(p), NA_real_,
         -ifelse(p > 0, p * log2(p), 0) - ifelse(p < 1, (1 - p) * log2(1 - p), 0))
}

kw_p <- function(x, g) {
  g <- droplevels(as.factor(g))
  if (nlevels(g) < 2) return(NA_real_)
  kruskal.test(x, g)$p.value
}

promoters <- fread(PROMOTER_TABLE)
promoters[, `:=`(win_start = as.integer(tss - HALF),
                 win_end   = as.integer(tss + HALF))]

assign_l <- list(); props_l <- list(); read_l <- list()

for (i in seq_len(nrow(promoters))) {
  gene <- promoters$gene[i]
  cat("== ", gene, " ==\n", sep = "")

  dat <- assemble_region_m6a(SAMPLES, promoters$chr[i],
                             promoters$win_start[i], promoters$win_end[i],
                             FT_RESULT_DIR, verbose = FALSE)
  met_mat <- filter_met_mat(dat$met_mat, verbose = TRUE)   # full-span reads
  rids_df <- dat$rids_df[match(rownames(met_mat),
                               as.character(dat$rids_df$RID)), ]

  res <- leiden_manhattan_cluster(met_mat, rids_df,
                                  region_start = promoters$win_start[i],
                                  region_end   = promoters$win_end[i],
                                  window_size  = WINDOW_SIZE,
                                  k_neighbors  = K_NEIGHBORS,
                                  resolution   = RESOLUTION,
                                  sigma        = NULL,
                                  seed         = LEIDEN_SEED,
                                  verbose      = FALSE)
  K <- res$n_clusters
  cat(nrow(res$assignments), "pooled full-span reads,", K, "Leiden clusters\n")

  asg <- as.data.table(res$assignments)
  asg[, sample_name := factor(sample_name, levels = SAMPLES)]
  assign_l[[gene]] <- data.table(region = gene, asg)

  ## cluster proportions per timepoint (composition, for the stacked bars) ----
  grid <- CJ(cluster = levels(asg$cluster), sample_name = SAMPLES)
  counts <- merge(grid, asg[, .N, by = .(cluster, sample_name)],
                  by = c("cluster", "sample_name"), all.x = TRUE)
  counts[is.na(N), N := 0L]
  counts[, n_total := sum(N), by = sample_name]
  counts[, prop := fifelse(n_total > 0, N / n_total, NA_real_)]
  counts <- merge(counts, sample_table, by = "sample_name")
  props_l[[gene]] <- counts[order(sample_name, cluster),
                            .(region = gene, sample_name, time, cluster,
                              n_reads = N, n_total, prop)]

  ## read-level entropy over the called-position (default) site universe -----
  M <- as.matrix(dat$met_mat)[as.character(asg$RID), , drop = FALSE]
  stopifnot(identical(rownames(M), as.character(asg$RID)))
  rd <- data.table(region = gene,
                   RID = as.character(asg$RID),
                   sample_name = asg$sample_name,
                   cluster = asg$cluster,
                   n_sites_covered = rowSums(!is.na(M)),
                   n_methylated = rowSums(M == 1, na.rm = TRUE))
  rd[, p_accessible := n_methylated / n_sites_covered]
  rd[, entropy := binary_entropy(p_accessible)]
  rd <- merge(rd, sample_table, by = "sample_name")
  rd[, sample_name := factor(sample_name, levels = SAMPLES)]  # merge drops the factor
  read_l[[gene]] <- rd[, .(region, RID, sample_name, time, cluster,
                           n_sites_covered, n_methylated, p_accessible,
                           entropy)]

  ## read-level heatmap: reads x sites, rows timepoint -> cluster -> H_read --
  o <- rd[, order(sample_name, cluster, -entropy)]
  rdo <- rd[o]
  hm_mat <- M[rdo$RID, , drop = FALSE]
  hm_mat <- matrix(ifelse(is.na(hm_mat), NA,
                          ifelse(hm_mat > 0, "m6A", "no m6A")),
                   nrow(hm_mat), ncol(hm_mat), dimnames = dimnames(hm_mat))
  ha <- ComplexHeatmap::rowAnnotation(
    timepoint    = rdo$sample_name,
    cluster      = rdo$cluster,
    H_read       = rdo$entropy,
    p_accessible = rdo$p_accessible,
    col = list(
      timepoint    = LEIDEN_TIMEPOINT_COLORS[levels(rdo$sample_name)],
      cluster      = cluster_palette(levels(rdo$cluster)),
      H_read       = colorRamp2(c(0, 0.5, 1), c("#440154", "#21918c", "#fde725")),
      p_accessible = colorRamp2(c(0, 1), c("white", "darkred"))),
    annotation_name_gp = grid::gpar(fontsize = 8))
  hm <- ComplexHeatmap::Heatmap(
    hm_mat,
    name = "m6A call",
    col  = c("m6A" = "black", "no m6A" = "white"),
    na_col = "grey85",
    show_row_names = FALSE, show_column_names = FALSE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = rdo$sample_name, row_gap = unit(0.8, "mm"),
    row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 8),
    width = unit(11, "cm"), height = unit(14, "cm"),
    use_raster = TRUE,
    column_title = paste0(gene, ": reads ordered timepoint -> cluster -> H_read"),
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    left_annotation = ha)
  pdf(file.path(PLOT_DIR, paste0("leiden_read_heatmap_", gene, ".pdf")),
      width = 8.5, height = 8)
  ComplexHeatmap::draw(hm)
  dev.off()

  ## per-gene entropy-vs-time violins, one facet per cluster ------------------
  kw_time <- rd[, .(p = kw_p(entropy, sample_name)), by = cluster][order(cluster)]
  p_time <- ggplot(rd, aes(sample_name, entropy, fill = sample_name)) +
    geom_violin(scale = "width", linewidth = 0.3) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 linewidth = 0.3) +
    facet_wrap(~ cluster, nrow = 1) +
    scale_fill_manual(values = LEIDEN_TIMEPOINT_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "read-level entropy (bits)",
         title = paste0(gene, ": read-level entropy across the LPS timecourse"),
         subtitle = paste0("Kruskal-Wallis across timepoints: ",
                           paste(sprintf("%s p=%.3g", kw_time$cluster,
                                         kw_time$p), collapse = ", "))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(PLOT_DIR,
                   paste0("leiden_read_entropy_by_time_", gene, ".pdf")),
         p_time, width = 1.5 + 1.8 * K, height = 4, limitsize = FALSE)
}

## tables ---------------------------------------------------------------------
assigns <- rbindlist(assign_l)
props   <- rbindlist(props_l)
reads   <- rbindlist(read_l)

stopifnot(reads[, all(n_methylated <= n_sites_covered)],
          reads[, all(p_accessible >= 0 & p_accessible <= 1)],
          reads[, all(entropy >= 0 & entropy <= 1)],
          nrow(reads) == nrow(assigns))

fwrite(assigns, file.path(TAB_DIR, "leiden_read_assignments.tsv.gz"), sep = "\t")
fwrite(props,   file.path(TAB_DIR, "leiden_cluster_props.tsv"), sep = "\t")
fwrite(reads,   file.path(TAB_DIR, "leiden_read_entropy.tsv.gz"), sep = "\t")

## Kruskal-Wallis comparisons --------------------------------------------------
med_str <- function(x, g) {
  m <- tapply(x, droplevels(as.factor(g)), median)
  paste(sprintf("%s:%.3f", names(m), m), collapse = ";")
}
tests_a <- reads[, .(test = "among_clusters", cluster = NA_character_,
                     n_reads = .N, n_groups = uniqueN(cluster),
                     kw_p = kw_p(entropy, cluster),
                     median_by_group = med_str(entropy, cluster)),
                 by = region]
tests_b <- reads[, .(test = "across_timepoints", n_reads = .N,
                     n_groups = uniqueN(sample_name),
                     kw_p = kw_p(entropy, sample_name),
                     median_by_group = med_str(entropy, sample_name)),
                 by = .(region, cluster)]
tests <- rbind(tests_a,
               tests_b[, .(region, test, cluster = as.character(cluster),
                           n_reads, n_groups, kw_p, median_by_group)])
setorder(tests, region, test, cluster, na.last = FALSE)
fwrite(tests, file.path(TAB_DIR, "leiden_read_entropy_tests.tsv"), sep = "\t")
cat("wrote leiden_{read_assignments,cluster_props,read_entropy,read_entropy_tests} to",
    TAB_DIR, "\n")
print(tests[test == "among_clusters"])

## entropy by cluster (timepoints pooled), one facet per region ----------------
pal <- cluster_palette(sort(unique(reads$cluster)))
kw_lab <- tests_a[, .(region, lab = sprintf("KW p = %.2g", kw_p))]
p_cl <- ggplot(reads, aes(cluster, entropy, fill = cluster)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
               linewidth = 0.3) +
  geom_text(data = kw_lab, aes(x = Inf, y = 0.04, label = lab),
            inherit.aes = FALSE, hjust = 1.1, size = 2.8, colour = "grey30") +
  facet_wrap(~ region, ncol = 4, scales = "free_x") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "read-level entropy (bits)",
       title = "Read-level entropy by Leiden cluster (timepoints pooled)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(PLOT_DIR, "leiden_read_entropy_by_cluster.pdf"), p_cl,
       width = 11, height = 5.5)

## cluster composition across time ---------------------------------------------
p_props <- ggplot(props[n_total > 0],
                  aes(factor(time, levels = c(0, 5, 10, 15)), prop,
                      fill = cluster)) +
  geom_col() +
  facet_wrap(~ region, ncol = 4) +
  scale_fill_manual(values = pal) +
  labs(x = "minutes LPS", y = "fraction of reads",
       title = "Leiden cluster proportions across the LPS timecourse",
       subtitle = "clusters from pooled-timepoint Manhattan-kNN Leiden clustering") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "leiden_cluster_props.pdf"), p_props,
       width = 11, height = 5.5)
cat("wrote leiden_read_entropy_by_cluster.pdf, leiden_cluster_props.pdf and",
    "per-gene leiden_read_{heatmap,entropy_by_time}_<gene>.pdf to", PLOT_DIR, "\n")
