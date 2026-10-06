# Tables and plots for one ACF fit across the shared active/inactive cohort.

ACF_CLASS_COLORS <- c(active = "#D55E00", inactive = "#0072B2")
ACF_TIMEPOINT_COLORS <- c(LPS_0 = "#CC79A7", LPS_5 = "#56B4E9",
                         LPS_10 = "#009E73", LPS_15 = "#E69F00")

acf_ordered_clusters <- function(clusters) {
  clusters <- unique(as.character(clusters[!is.na(clusters)]))
  numbered <- setdiff(clusters, "unclustered")
  numbered <- numbered[order(as.integer(sub("^cluster", "", numbered)))]
  c(numbered, if ("unclustered" %in% clusters) "unclustered")
}

acf_cluster_palette <- function(clusters) {
  numbered <- setdiff(clusters, "unclustered")
  colors <- stats::setNames(grDevices::hcl.colors(length(numbered), "Dark 3"), numbered)
  if ("unclustered" %in% clusters) colors <- c(colors, unclustered = "grey65")
  colors
}

pooled_acf_cluster_counts <- function(assignments) {
  a <- data.table::copy(data.table::as.data.table(assignments))
  a[, cluster := as.character(cluster)]
  a[is.na(cluster), cluster := "unclustered"]
  out <- a[, .(n_reads = .N, n_enhancers = data.table::uniqueN(enhancer_id)),
           by = .(enhancer_class, sample_name, cluster)]
  out[, fraction := n_reads / sum(n_reads), by = .(enhancer_class, sample_name)]
  out
}

pooled_acf_qc <- function(result) {
  extraction <- data.table::as.data.table(result$extraction_qc)[,
    .(n_enhancers = data.table::uniqueN(enhancer_id), n_reads = sum(n_reads)),
    by = .(enhancer_class, status)]
  extraction[, stage := "extraction"]
  clustering <- data.table::as.data.table(result$assignments)[,
    .(n_enhancers = data.table::uniqueN(enhancer_id), n_reads = .N),
    by = .(enhancer_class, status = acf_status)]
  clustering[, stage := "acf_clustering"]
  out <- data.table::rbindlist(list(extraction, clustering), use.names = TRUE)
  out[, .(enhancer_class, stage, status, n_enhancers, n_reads)]
}

pooled_acf_profile_table <- function(mat, position_name, value_name) {
  out <- data.table::data.table(cluster = character(), position = integer(), value = numeric())
  if (!is.null(mat) && nrow(mat))
    out <- data.table::data.table(cluster = rep(rownames(mat), each = ncol(mat)),
      position = rep(as.integer(colnames(mat)), nrow(mat)), value = as.vector(t(mat)))
  data.table::setnames(out, c("position", "value"), c(position_name, value_name))
  out
}

plot_pooled_acf_profile <- function(profiles, profile_type = c("acf", "m6a")) {
  profile_type <- match.arg(profile_type)
  dt <- data.table::copy(data.table::as.data.table(profiles))
  if (!nrow(dt)) return(NULL)
  if (profile_type == "acf") {
    data.table::setnames(dt, c("lag_bp", "mean_acf"), c("position", "mean_value"))
    dt <- dt[position > 0L]
  } else {
    data.table::setnames(dt, c("position_bp", "mean_m6a"), c("position", "mean_value"))
  }
  dt[, cluster := factor(cluster, levels = acf_ordered_clusters(cluster))]
  p <- ggplot2::ggplot(dt, ggplot2::aes(position, mean_value)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
    ggplot2::geom_line(linewidth = 0.3, colour = "#176B87") +
    ggplot2::facet_wrap(~cluster, ncol = 1) + ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = "Pooled active and inactive enhancers: ACF clusters",
      subtitle = "Shared cluster assignments across both enhancer classes; unsmoothed means",
      x = if (profile_type == "acf") "Lag (bp)" else "Position relative to enhancer midpoint (bp)",
      y = if (profile_type == "acf") "Mean autocorrelation" else "Fraction of fibers with m6A")
  if (profile_type == "m6a")
    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = 0.3)
  p
}

plot_pooled_acf_composition <- function(counts, samples) {
  if (!nrow(counts)) return(NULL)
  dt <- data.table::copy(data.table::as.data.table(counts))
  clusters <- acf_ordered_clusters(dt$cluster)
  dt[, `:=`(cluster = factor(cluster, levels = clusters),
    enhancer_class = factor(enhancer_class, levels = names(ACF_CLASS_COLORS)),
    sample_name = factor(sample_name, levels = samples))]
  ggplot2::ggplot(dt, ggplot2::aes(sample_name, fraction, fill = cluster)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::facet_wrap(~enhancer_class, nrow = 1, drop = FALSE) +
    ggplot2::scale_fill_manual(values = acf_cluster_palette(clusters), drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE, labels = function(x) sub("^LPS_", "", x)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(100 * x, "%")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = "ACF cluster composition: shared pooled clustering",
      subtitle = "Denominator: all selected fibers within each enhancer class and timepoint",
      x = "Time after LPS (min)", y = "Fraction of selected fibers", fill = "Cluster")
}

plot_pooled_acf_umap <- function(assignments, samples, colour_by = c("cluster", "sample_name", "enhancer_class")) {
  colour_by <- match.arg(colour_by)
  dt <- data.table::copy(data.table::as.data.table(assignments))[
    acf_status == "clustered" & is.finite(UMAP1) & is.finite(UMAP2)]
  if (!nrow(dt)) return(NULL)
  groups <- switch(colour_by, cluster = acf_ordered_clusters(dt$cluster),
                   sample_name = samples, enhancer_class = names(ACF_CLASS_COLORS))
  dt[, plot_group := factor(get(colour_by), levels = groups)]
  p <- ggplot2::ggplot(dt, ggplot2::aes(UMAP1, UMAP2, colour = plot_group)) +
    ggplot2::geom_point(size = 0.35, alpha = 0.45, stroke = 0) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
    ggplot2::coord_equal() + ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = "Pooled active and inactive enhancers: ACF UMAP",
      subtitle = "One embedding for the shared cohort; fibers without coordinates are omitted",
      colour = switch(colour_by, cluster = "Cluster", sample_name = "Time after LPS",
                      enhancer_class = "Enhancer class"))
  colors <- switch(colour_by, cluster = acf_cluster_palette(groups),
                    sample_name = ACF_TIMEPOINT_COLORS, enhancer_class = ACF_CLASS_COLORS)
  if (all(groups %in% names(colors)))
    p <- p + ggplot2::scale_colour_manual(values = colors, drop = FALSE)
  p
}

save_enhancer_acf_outputs <- function(result, plot_dir, table_dir,
                                     samples = c("LPS_0", "LPS_5", "LPS_10", "LPS_15"),
                                     suffix = "_pooled_capped") {
  stopifnot(!anyDuplicated(result$assignments$RID),
    all(result$assignments$enhancer_class %in% c("active", "inactive")),
    all(result$assignments$sample_name %in% samples))
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  save_table <- function(x, stem)
    data.table::fwrite(x, file.path(table_dir, paste0("acf_", stem, suffix, ".tsv")), sep = "\t")
  counts <- pooled_acf_cluster_counts(result$assignments)
  qc <- pooled_acf_qc(result)
  acf_profiles <- pooled_acf_profile_table(result$acf_profiles, "lag_bp", "mean_acf")
  m6a_profiles <- pooled_acf_profile_table(result$m6a_profiles, "position_bp", "mean_m6a")
  save_table(result$assignments, "assignments")
  save_table(counts, "cluster_counts")
  save_table(qc, "qc")
  save_table(acf_profiles, "cluster_profiles")
  save_table(m6a_profiles, "m6a_base_cluster_profiles")
  save_table(result$fiber_sampling, "fiber_sampling")
  save_table(result$sampling_diagnostics, "sampling_diagnostics")
  save_table(result$cluster_diagnostics, "cluster_diagnostics")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for ACF PDF export")
  save_plot <- function(plot, stem, width = 8, height = 6) {
    if (!is.null(plot)) ggplot2::ggsave(
      file.path(plot_dir, paste0("acf_", stem, suffix, ".pdf")),
      plot = plot, width = width, height = height, limitsize = FALSE)
  }
  n_clusters <- length(rownames(result$acf_profiles))
  save_plot(plot_pooled_acf_profile(acf_profiles, "acf"), "cluster_profiles",
            height = max(4, 1.5 * n_clusters + 1))
  save_plot(plot_pooled_acf_profile(m6a_profiles, "m6a"), "m6a_base_cluster_profiles",
            height = max(4, 1.5 * n_clusters + 1))
  save_plot(plot_pooled_acf_composition(counts, samples), "cluster_composition", width = 10)
  for (colour_by in c("cluster", "sample_name", "enhancer_class"))
    save_plot(plot_pooled_acf_umap(result$assignments, samples, colour_by), paste0("umap_", colour_by))
  if (!n_clusters)
    warning("No clustered ACF fibers: assignments, diagnostics and available composition outputs were saved.",
            call. = FALSE)
  invisible(list(counts = counts, qc = qc, acf_profiles = acf_profiles, m6a_profiles = m6a_profiles))
}
