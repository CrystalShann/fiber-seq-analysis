# Selected-fiber signals and plots for pooled enhancer Manhattan clusters.
# Coordinates in records are midpoint-relative, 0-based half-open [-500, 500).
ENHANCER_CLASS_COLORS <- c(active = "#E69F00", inactive = "#009E73")

enhancer_manhattan_plot_view <- function(result) {
  if (is.null(result$feat_mat) || is.null(result$footprints))
    stop("This result predates saved plot data. Rerun the pooled Manhattan job with k=100.")
  a <- data.table::copy(data.table::as.data.table(result$assignments))
  levels <- rownames(result$m6a_profiles)
  a[, cluster := factor(cluster, levels = levels)]
  a[, sample_name := factor(sample_name, levels = names(LEIDEN_TIMEPOINT_COLORS))]
  if ("enhancer_class" %in% names(a))
    a[, enhancer_class := factor(enhancer_class, levels = names(ENHANCER_CLASS_COLORS))]
  data.table::setorder(a, cluster, sample_name, enhancer_id, read_id)
  # The reference orders by start; equal relative starts preserve the order above.
  a[, start := -500L]
  stopifnot(setequal(rownames(result$feat_mat), a$RID), ncol(result$feat_mat) == 1000L)
  list(assignments = as.data.frame(a), feat_mat = as.matrix(result$feat_mat[a$RID, , drop = FALSE]),
       profiles = result$m6a_profiles, params = list(window_size = 1L),
       n_clusters = length(levels))
}

plot_enhancer_manhattan_heatmap <- function(view, main = NULL) {
  ht <- plot_cluster_heatmap(view, main = main)
  if (!"enhancer_class" %in% names(view$assignments)) return(ht)
  a <- view$assignments
  a <- a[order(a$cluster, a$start), , drop = FALSE]
  ht + ComplexHeatmap::rowAnnotation(
    enhancer_class = a$enhancer_class,
    col = list(enhancer_class = ENHANCER_CLASS_COLORS))
}

plot_enhancer_manhattan_composition <- function(view, main = NULL) {
  original <- plot_cluster_composition(view, main = main)
  if (!"enhancer_class" %in% names(view$assignments)) return(original)
  a <- view$assignments
  by_class <- ggplot2::ggplot(a, ggplot2::aes(sample_name, fill = cluster)) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::facet_wrap(~enhancer_class, nrow = 1) +
    ggplot2::scale_fill_manual(values = cluster_palette(levels(a$cluster))) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::labs(x = "Timepoint", y = "Fraction within class and timepoint",
      title = "Shared clusters within active and inactive enhancers")
  class_fraction <- ggplot2::ggplot(a, ggplot2::aes(cluster, fill = enhancer_class)) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::scale_fill_manual(values = ENHANCER_CLASS_COLORS, drop = FALSE) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = "Shared cluster", y = "Fraction of cluster fibers",
      fill = "Enhancer class", title = "Active/inactive composition of each cluster")
  cowplot::plot_grid(original, by_class, class_fraction, ncol = 1, rel_heights = c(2, 1, 1))
}

enhancer_manhattan_signal_profiles <- function(result) {
  a <- data.table::as.data.table(result$assignments)
  records <- data.table::as.data.table(result$footprints)
  profiles <- list(); part <- 0L
  for (cl in rownames(result$m6a_profiles)) {
    ids <- a[cluster == cl, RID]
    part <- part + 1L
    profiles[[part]] <- data.table::data.table(cluster = cl, track = "m6A",
      position_bp = -500:499, fraction = as.numeric(result$m6a_profiles[cl, ]), n_reads = length(ids))
    for (track_id in c("Nucleosome >90 bp", "TF <60 bp")) {
      r <- records[RID %in% ids & track == track_id]
      fraction <- numeric(1000L)
      if (nrow(r)) {
        ranges <- IRanges::IRanges(r$start + 501L, r$end + 500L)
        # Union calls within a fiber before computing per-base occupancy.
        merged <- unlist(IRanges::reduce(split(ranges, r$RID)), use.names = FALSE)
        fraction <- as.numeric(IRanges::coverage(merged, width = 1000L)) / length(ids)
      }
      part <- part + 1L
      profiles[[part]] <- data.table::data.table(cluster = cl, track = track_id,
        position_bp = -500:499, fraction = fraction, n_reads = length(ids))
    }
  }
  data.table::rbindlist(profiles)
}

plot_enhancer_manhattan_signals <- function(profiles, class_id) {
  levels <- unique(profiles$cluster)
  profiles <- data.table::copy(profiles)
  profiles[, cluster := factor(cluster, levels = levels)]
  ggplot2::ggplot(profiles, ggplot2::aes(position_bp, fraction, colour = track)) +
    ggplot2::geom_line(linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    ggplot2::facet_wrap(~cluster, ncol = 1) +
    ggplot2::scale_colour_manual(values = c("m6A" = "black", "Nucleosome >90 bp" = "#6488A3", "TF <60 bp" = "#CC3377")) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(title = paste(class_id, "enhancers: Manhattan-Leiden"),
      x = "Position relative to enhancer midpoint (bp)", y = "Fraction of sampled fibers", colour = NULL)
}

enhancer_manhattan_footprint_height <- function(result) {
  # Keep a single vector PDF page within the usual 200-inch page-size limit.
  # This limits physical page height, never the number of fibers displayed.
  min(190, max(10, 2.5 + 0.018 * nrow(result$assignments) +
                 0.25 * nrow(result$m6a_profiles)))
}

plot_enhancer_manhattan_footprints <- function(result, class_id) {
  a <- data.table::copy(data.table::as.data.table(result$assignments))
  cluster_levels <- rownames(result$m6a_profiles)
  a <- a[!is.na(cluster) & cluster %in% cluster_levels]
  if (!nrow(a)) return(NULL)
  a[, cluster := factor(cluster, levels = cluster_levels)]
  a[, sample_name := factor(sample_name, levels = names(LEIDEN_TIMEPOINT_COLORS))]
  if (anyNA(a$sample_name)) stop("Unknown timepoint in saved Manhattan assignments")
  data.table::setorder(a, cluster, sample_name, enhancer_id, read_id)
  a[, row := seq_len(.N), by = cluster]
  counts <- a[, .(n = .N), by = cluster]
  cluster_labels <- stats::setNames(
    paste0(counts$cluster, "\n", counts$n, " displayed / ", counts$n, " total"),
    as.character(counts$cluster))
  cluster_cols <- cluster_palette(cluster_levels)
  time_cols <- LEIDEN_TIMEPOINT_COLORS
  has_class <- "enhancer_class" %in% names(a)
  class_cols <- if (has_class) ENHANCER_CLASS_COLORS else character()
  if (has_class && anyNA(match(a$enhancer_class, names(class_cols))))
    stop("Unknown enhancer class in saved Manhattan assignments")
  track_cols <- c("Nucleosome >90 bp" = "#6488A3", "TF <60 bp" = "#CC3377")
  r <- data.table::copy(data.table::as.data.table(result$footprints))[RID %in% a$RID]
  r[, `:=`(row = a$row[match(RID, a$RID)], cluster = a$cluster[match(RID, a$RID)])]
  r[, track := factor(track, levels = c("Nucleosome >90 bp", "TF <60 bp"))]
  data.table::setorder(r, track)
  calls <- summary(result$feat_mat[a$RID, , drop = FALSE])
  m <- data.table::data.table(position = calls$j - 501L + 0.5,
    row = a$row[calls$i], cluster = a$cluster[calls$i])
  bars <- data.table::rbindlist(list(
    a[, .(cluster, row, xmin = -556, xmax = -534, value = as.character(cluster))],
    a[, .(cluster, row, xmin = -528, xmax = -506, value = as.character(sample_name))]))
  if (has_class)
    bars <- data.table::rbindlist(list(
      a[, .(cluster, row, xmin = -584, xmax = -562, value = as.character(enhancer_class))], bars))
  # Scale the baseline stroke down if the page-height limit compresses rows.
  density_scale <- min(1, enhancer_manhattan_footprint_height(result) /
                         (2.5 + 0.018 * nrow(a) + 0.25 * length(cluster_levels)))
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = a, ggplot2::aes(x = -500, xend = 500, y = row, yend = row),
      colour = "grey80", linewidth = 0.15 * density_scale) +
    ggplot2::geom_rect(data = r, ggplot2::aes(xmin = start, xmax = end, ymin = row - .38, ymax = row + .38, fill = track)) +
    ggplot2::geom_segment(data = m, ggplot2::aes(x = position, xend = position, y = row - .3, yend = row + .3), colour = "black", linewidth = 0.2) +
    ggplot2::geom_rect(data = bars, ggplot2::aes(xmin = xmin, xmax = xmax,
      ymin = row - .5, ymax = row + .5, fill = value)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = .25) +
    ggplot2::facet_grid(cluster ~ ., scales = "free_y", space = "free_y",
      labeller = ggplot2::as_labeller(cluster_labels)) +
    ggplot2::scale_fill_manual(values = c(track_cols, cluster_cols, time_cols, class_cols),
      breaks = c(names(track_cols), names(time_cols), names(class_cols)), drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq(-500, 500, 250), expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(if (has_class) -590 else -562, 500), expand = FALSE) +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(add = .6)) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(), legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text.y = ggplot2::element_text(angle = 0), panel.spacing.y = grid::unit(1.5, "mm")) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    ggplot2::labs(title = paste(class_id, "enhancers | Manhattan-Leiden k =", result$info$k_neighbors),
      subtitle = paste0("All ", nrow(a), " selected fibers in ", length(cluster_levels),
        " clusters; left bars: ", if (has_class) "class, " else "",
        "cluster, timepoint; black marks: m6A\n",
        "Within each cluster: timepoint, enhancer, read; dashed line: enhancer midpoint"),
      x = "Position relative to enhancer midpoint (bp)", y = "Fibers", fill = NULL)
}

save_enhancer_manhattan_pdfs <- function(result, class_id, plot_dir, table_dir,
                                       expected_k = 100L, suffix = "_capped") {
  if (is.null(result$m6a_profiles) || !nrow(result$m6a_profiles)) {
    warning("No clustered fibers to plot for ", class_id, "; saved diagnostics remain available.", call. = FALSE)
    return(invisible(NULL))
  }
  if (!identical(as.integer(result$info$k_neighbors), as.integer(expected_k)))
    stop("Saved Manhattan clusters use k=", result$info$k_neighbors, "; rerun the array for k=", expected_k)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  view <- enhancer_manhattan_plot_view(result)
  title <- paste(class_id, "enhancers | Manhattan-Leiden k =", expected_k)
  pdf_path <- function(stem) file.path(plot_dir, paste0("manhattan_", stem, "_", class_id, suffix, ".pdf"))
  grDevices::pdf(pdf_path("m6a_heatmap"), width = 8, height = 10)
  tryCatch(ComplexHeatmap::draw(plot_enhancer_manhattan_heatmap(view, main = title)), finally = grDevices::dev.off())
  ggplot2::ggsave(pdf_path("cluster_composition"), plot_enhancer_manhattan_composition(view, main = title),
    width = 10, height = if ("enhancer_class" %in% names(view$assignments)) 14 else 8)
  ggplot2::ggsave(pdf_path("cluster_structure"), plot_cluster_structure(view, main = title), width = 9, height = 3)
  p <- plot_met_fraction_lines(view, tss = 0, main = title, smooth_k = 1) +
    ggplot2::labs(x = "Position relative to enhancer midpoint (bp)")
  ggplot2::ggsave(pdf_path("m6a_base_cluster_profiles"), p, width = 8, height = max(4, 1.4 * view$n_clusters + 1), limitsize = FALSE)
  profiles <- enhancer_manhattan_signal_profiles(result)
  data.table::fwrite(profiles, file.path(table_dir, paste0("manhattan_signal_profiles_", class_id, suffix, ".tsv")), sep = "\t")
  ggplot2::ggsave(pdf_path("footprint_profiles"), plot_enhancer_manhattan_signals(profiles, class_id),
    width = 9, height = max(4, 1.5 * view$n_clusters + 1), limitsize = FALSE)
  ggplot2::ggsave(pdf_path("single_fiber_footprints"),
    plot_enhancer_manhattan_footprints(result, class_id), width = 11,
    height = enhancer_manhattan_footprint_height(result), limitsize = FALSE)
  invisible(profiles)
}
