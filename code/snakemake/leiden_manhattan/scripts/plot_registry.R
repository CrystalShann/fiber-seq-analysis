# The single panel-name -> PDF filename/function registry for every dataset.
# A renderer returns plot_panel(); expensive shared summaries are cached in context.
panel_entry <- function(filename, render, enabled = NULL) {
  if (is.character(filename)) {
    name <- filename
    filename <- function(region_id) name
  }
  list(filename = filename, render = render, enabled = enabled)
}

panel_registry <- list(
  cluster_met_profiles = panel_entry(function(id) paste0("cluster_met_profiles_", id, ".pdf"), function(c) {
    r <- context_profile_result(c)
    p <- plot_cluster_met_profiles(r, r$site_met_mat, main = r$region$region_id,
                                   palette = cluster_id_palette, x_label = plot_anchor(c$result$region)$x_label)
    plot_panel(p, width = 6, height = 1.4 * r$n_clusters + 1.4)
  }),
  cluster_composition = panel_entry(function(id) paste0("cluster_composition_", id, ".pdf"), function(c) {
    plot_panel(plot_cluster_composition(c$result, main = c$result$region$region_id,
      timepoint_cols = c$timepoint_colors, palette = cluster_id_palette), width = 7, height = 8)
  }),
  heatmap = panel_entry(function(id) paste0("heatmap_", id, ".pdf"), function(c) {
    plot_panel(plot_cluster_heatmap(c$result, main = c$result$region$region_id,
      timepoint_cols = c$timepoint_colors, palette = cluster_id_palette), width = 8, height = 9,
      draw = function(x) ComplexHeatmap::draw(x, newpage = FALSE))
  }),
  met_fraction_lines = panel_entry(function(id) paste0("met_fraction_lines_", id, ".pdf"), function(c) {
    r <- context_profile_result(c)
    smooth <- c$ds$plot_options$line_smooth_k
    if (is.null(smooth)) smooth <- 1L
    p <- plot_met_fraction_lines(r, main = r$region$region_id, smooth_k = smooth,
                                 palette = cluster_id_palette, x_label = plot_anchor(c$result$region)$x_label)
    plot_panel(p, width = 8, height = 1.4 * r$n_clusters + 1)
  }),
  fiberseq_profiles_composition = panel_entry("fiberseq_profiles_composition.pdf", function(c) {
    p <- context_fiberseq_profiles(c)
    plot_panel(p$plot, width = 16, height = p$height)
  }),
  cluster_proportions = panel_entry("fiberseq_cluster_proportions.pdf", function(c) {
    plot_panel(plot_cluster_pie(c$result), width = 8, height = 5)
  }),
  umap_overview = panel_entry("fiberseq_umap_overview.pdf", context_umap_panel),
  haplotype_example = panel_entry("fiberseq_haplotype_example.pdf", context_haplotype_panel),
  knn_graph = panel_entry("knn_graph.pdf", function(c) {
    plot_panel(plot_knn_graph(c$result), width = 12, height = 9)
  }),
  knn_graph_by_allele = panel_entry("knn_graph_by_allele.pdf", function(c) {
    context_allele_summary(c)
    plot_panel(plot_knn_graph(c$result, color_by = "allele"), width = 12, height = 9)
  }),
  knn_allele_composition = panel_entry("knn_allele_composition.pdf", function(c) {
    summary <- context_allele_summary(c)
    p <- plot_knn_allele_composition(summary,
      paste(c$result$region$region_id, "Focal SNP allele proportions", sep = "\n"))
    plot_panel(p, width = max(10, .7 * nlevels(summary$composition$group) + 3), height = 6)
  }),
  knn_allele_connectivity = panel_entry("knn_allele_connectivity.pdf", function(c) {
    p <- plot_knn_allele_connectivity(context_allele_summary(c),
      paste(c$result$region$region_id, "Allele mixing in the KNN graph", sep = "\n"))
    plot_panel(p, width = 11, height = 7)
  }),
  methylation_proportion_by_cluster = panel_entry("methylation_proportion_by_cluster.pdf", function(c) {
    summary <- methylation_by_cluster(c$result)
    plot_table(c, summary, "methylation_proportion_by_cluster.tsv")
    plot_panel(plot_methylation_by_cluster(summary, c$result$region), width = 8, height = 4.5)
  }),
  heatmap_m6a_footprints = panel_entry("heatmap_m6a_footprints.pdf", context_genomic_heatmap),
  read_footprints = panel_entry("fig1_read_footprints.pdf", context_read_footprints),
  occupancy_by_cluster = panel_entry("fig2_occupancy_by_cluster.pdf", context_occupancy),
  cluster_sample_composition = panel_entry("cluster_sample_composition.pdf", function(c) {
    p <- plot_sample_composition(c$result, c$sample_colors, c$result$region$annotation,
                                 sample_column = "sample_label", label_fun = identity)
    plot_panel(p, width = 12, height = 5)
  }),
  pairwise_m6a_delta = panel_entry("pairwise_m6a_delta.pdf", context_pairwise_delta,
    enabled = function(c) {
      limit <- c$ds$plot_options$delta_top_n
      rank <- c$result$region$rank
      is.null(limit) || is.null(rank) || (!is.na(rank) && rank <= limit)
    })
)
