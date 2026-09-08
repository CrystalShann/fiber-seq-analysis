# entropy_functions.r
#
# Cluster-level positional Shannon entropy of the Leiden + Manhattan single
# molecule clusters (Leiden_Manhattan/leiden_manhattan.Rmd). The clustering is
# a fixed input: entropy quantifies how internally homogeneous each Leiden
# cluster is and never modifies, splits or reruns it.
#
# For each Leiden cluster c and m6A site j of the full-span filtered site
# matrix (filter_met_mat() output, as in the topic model):
#   p_{c,j} = n_accessible,c,j / n_reads,c
#   H_{c,j} = -p log2 p - (1 - p) log2 (1 - p)      (0 log 0 = 0; H in [0,1])
#   H_c     = (1/L) sum_j H_{c,j}                    (L = number of sites)
# H_c ~ 0: the cluster's reads agree at nearly every site; H_c ~ 1: many sites
# are split about half and half among the reads.
#
# Sourced by leiden_cluster_entropy.R. The plot functions need
# Leiden_Manhattan/leiden_manhattan_plots.r (cluster_palette,
# LEIDEN_TIMEPOINT_COLORS, as_timepoint_factor) sourced first.

ENTROPY_DIR <- "/project/spott/cshan/fiber-seq/macrophage_project/clustering_TSS_methods/entropy"
TAB_DIR     <- file.path(ENTROPY_DIR, "tables")
PLOT_DIR    <- file.path(ENTROPY_DIR, "plots")

# viridis gradient for entropy tracks (0 = dark purple, 1 = yellow)
entropy_col <- function() circlize::colorRamp2(c(0, 0.5, 1),
                                               c("#440154", "#21918c", "#fde725"))


# ---------------------------------------------------------------------------
# Binary Shannon entropy in bits, 0 log 0 = 0.
# ---------------------------------------------------------------------------
binary_entropy <- function(p) {
  h  <- numeric(length(p))
  ok <- p > 0 & p < 1
  h[ok] <- -p[ok] * log2(p[ok]) - (1 - p[ok]) * log2(1 - p[ok])
  h
}


# ---------------------------------------------------------------------------
# met_mat:     full-span filtered reads x m6A-site 0/1 matrix used for the
#              clustering (no NA expected)
# assignments: res$assignments of leiden_manhattan_cluster() (RID, cluster)
# Returns clusters (cluster, n_reads, n_sites, cluster_entropy) and positions
# (cluster, position, p_accessible, positional_entropy). No read threshold.
# ---------------------------------------------------------------------------
cluster_positional_entropy <- function(met_mat, assignments) {
  M <- as.matrix(met_mat)[as.character(assignments$RID), , drop = FALSE]
  stopifnot(!anyNA(M))
  lv  <- levels(assignments$cluster)
  pos <- as.integer(colnames(M))

  positions <- do.call(rbind, lapply(lv, function(cl) {
    rows <- assignments$cluster == cl
    p <- colSums(M[rows, , drop = FALSE] == 1) / sum(rows)
    data.frame(cluster = cl, position = pos, p_accessible = as.numeric(p),
               positional_entropy = binary_entropy(as.numeric(p)),
               row.names = NULL)
  }))
  positions$cluster <- factor(positions$cluster, levels = lv)

  clusters <- data.frame(
    cluster = factor(lv, levels = lv),
    n_reads = as.integer(table(assignments$cluster)[lv]),
    n_sites = ncol(M),
    cluster_entropy = as.numeric(tapply(positions$positional_entropy,
                                        positions$cluster, mean)[lv]),
    row.names = NULL)

  list(clusters = clusters, positions = positions)
}


# ---------------------------------------------------------------------------
# The Leiden_Manhattan read x feature heatmap (plot_cluster_heatmap(): rows =
# reads split by cluster and ordered by read start, black = any m6A call,
# cluster + timepoint row tracks) with one extra continuous row track: each
# read carries its cluster's H_c.
# cluster_entropy = the `clusters` table of cluster_positional_entropy().
# ---------------------------------------------------------------------------
plot_cluster_entropy_heatmap <- function(res, cluster_entropy, main = NULL,
                                         timepoint_cols = LEIDEN_TIMEPOINT_COLORS) {
  df <- res$assignments
  df$sample_name <- as_timepoint_factor(df$sample_name, timepoint_cols)
  df <- df[order(df$cluster, df$start), ]
  mat <- res$feat_mat[df$RID, , drop = FALSE]
  mat <- matrix(ifelse(is.na(mat), NA, ifelse(mat > 0, "m6A", "no m6A")),
                nrow(mat), ncol(mat), dimnames = dimnames(mat))

  h_c <- cluster_entropy$cluster_entropy[
    match(as.character(df$cluster), as.character(cluster_entropy$cluster))]
  ha <- ComplexHeatmap::rowAnnotation(
    cluster         = df$cluster,
    timepoint       = df$sample_name,
    cluster_entropy = h_c,
    col = list(cluster         = cluster_palette(levels(df$cluster)),
               timepoint       = timepoint_cols[levels(df$sample_name)],
               cluster_entropy = entropy_col()))

  ComplexHeatmap::Heatmap(
    mat,
    name = "m6A call",
    col  = c("m6A" = "black", "no m6A" = "white"),
    na_col = "grey85",
    show_row_names = FALSE, show_column_names = FALSE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = df$cluster, row_gap = grid::unit(0.6, "mm"),
    row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 8),
    width = grid::unit(11, "cm"), height = grid::unit(14, "cm"),
    use_raster = TRUE,
    column_title = main,
    column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
    left_annotation = ha)
}


# ---------------------------------------------------------------------------
# Cluster entropy summary: one bar per region x cluster (H_c), labelled with
# the cluster's read count. `tab` = the bound `clusters` tables with a
# `region` column.
# ---------------------------------------------------------------------------
plot_cluster_entropy_summary <- function(tab, main = NULL) {
  tab$cluster <- factor(as.character(tab$cluster),
                        levels = unique(as.character(tab$cluster)))
  tab$region  <- factor(as.character(tab$region),
                        levels = unique(as.character(tab$region)))
  ggplot(tab, aes(x = region, y = cluster_entropy, fill = cluster)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = n_reads, y = cluster_entropy + 0.02),
              position = position_dodge(width = 0.8), size = 2.2, vjust = 0) +
    scale_fill_manual(values = cluster_palette(levels(tab$cluster))) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = "cluster entropy (mean positional H, bits)",
         title = main) +
    cowplot::theme_cowplot(font_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
