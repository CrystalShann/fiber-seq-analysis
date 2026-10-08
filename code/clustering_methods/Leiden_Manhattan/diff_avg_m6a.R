source("/project/spott/cshan/fiber-seq/code/parsing_functions/plotting_functions.r", local = TRUE)

# diff_avg_m6a.R
#
# Pairwise per-bp differences in mean m6A between Leiden clusters of one
# region result. Requires leiden_manhattan_plots.r to be sourced first for
# cluster_site_profiles() and plot_anchor().
#
# Sign convention: delta = met(A) - met(B); positive = A more methylated than B.
# Pairs are labelled "A − B" (U+2212 minus), e.g. "cluster1 − cluster2".
#
#   m6a_delta_pairs()   clusters with reads, all unique (A, B) pairs or
#                       baseline vs each other cluster
#   m6a_delta_table()   long table: pair, cluster_A, cluster_B, n_A, n_B, pos,
#                       delta, delta_smooth (centred rolling mean)
#   plot_m6a_delta()    one panel per pair, raw + smoothed line, optional shading
#   run_m6a_delta()     one saved clustering.rds -> TSV + PDF
#   run_all_m6a_delta() several saved clustering.rds files

M6A_DELTA_MINUS <- "−"

m6a_delta_pairs <- function(res, baseline = NULL) {
  cluster <- res$assignments$cluster
  clusters <- if (is.factor(cluster)) levels(cluster) else unique(as.character(cluster))
  counts <- table(factor(as.character(cluster), levels = clusters))
  clusters <- clusters[counts > 0]
  if (length(clusters) < 2L) {
    message("Fewer than 2 clusters with reads: skipping ", res$region$region_id)
    return(NULL)
  }
  if (is.null(baseline)) {
    pairs <- t(utils::combn(clusters, 2L))
  } else {
    stopifnot(length(baseline) == 1L, baseline %in% clusters)
    pairs <- cbind(baseline, setdiff(clusters, baseline))
  }
  data.frame(cluster_A = pairs[, 1], cluster_B = pairs[, 2],
             n_A = as.integer(counts[pairs[, 1]]), n_B = as.integer(counts[pairs[, 2]]),
             stringsAsFactors = FALSE)
}

m6a_delta_table <- function(res, baseline = NULL, smooth_bp = 25) {
  stopifnot(length(smooth_bp) == 1L, smooth_bp >= 1, smooth_bp == round(smooth_bp))
  pairs <- m6a_delta_pairs(res, baseline)
  if (is.null(pairs)) return(NULL)
  prof <- cluster_site_profiles(res, res$site_met_mat)
  out <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    a <- prof[prof$cluster == pairs$cluster_A[i], ]
    b <- prof[prof$cluster == pairs$cluster_B[i], ]
    a <- a[order(a$pos), ]
    b <- b[order(b$pos), ]
    if (!identical(a$pos, b$pos))
      stop("Position sets differ between ", pairs$cluster_A[i], " and ", pairs$cluster_B[i])
    delta <- a$met - b$met
    smooth <- if (smooth_bp > 1) as.numeric(stats::filter(delta, rep(1 / smooth_bp, smooth_bp), sides = 2)) else delta
    data.frame(pair = paste(pairs$cluster_A[i], M6A_DELTA_MINUS, pairs$cluster_B[i]),
               cluster_A = pairs$cluster_A[i], cluster_B = pairs$cluster_B[i],
               n_A = pairs$n_A[i], n_B = pairs$n_B[i],
               pos = a$pos, delta = delta, delta_smooth = smooth,
               stringsAsFactors = FALSE)
  }))
  stopifnot(all(is.finite(out$delta)), all(abs(out$delta) <= 1))
  out$pair <- factor(out$pair, levels = unique(out$pair))
  out
}

# Insert y = 0 rows where the smoothed line crosses zero so the above/below
# ribbons meet exactly at the crossing instead of leaving a gap.
m6a_delta_zero_crossings <- function(x, y) {
  ok <- !is.na(y[-length(y)]) & !is.na(y[-1])
  cross <- which(ok & y[-length(y)] * y[-1] < 0)
  if (!length(cross)) return(data.frame(x = x, y = y))
  x0 <- x[cross] + (x[cross + 1] - x[cross]) * y[cross] / (y[cross] - y[cross + 1])
  d <- rbind(data.frame(x = x, y = y), data.frame(x = x0, y = 0))
  d[order(d$x), ]
}

plot_m6a_delta <- function(delta, region = NULL, ylim = c(-1, 1), shade = TRUE, smooth_bp = 25,
                           above_color = "#C0392B", below_color = "#2874A6") {
  stopifnot(length(ylim) == 2L, ylim[1] < 0, ylim[2] > 0)
  if (!is.null(region)) {
    anchor <- plot_anchor(region)
    delta$x <- plot_positions(delta$pos, region)
    x_label <- anchor$x_label
    xlim <- c(anchor$left - 0.5, anchor$right + 0.5)
    title <- paste0(region$annotation, " | ", region$chr, ":", region$analysis_start, "-", region$analysis_end)
  } else {
    delta$x <- delta$pos
    x_label <- "Genomic position (bp)"
    xlim <- range(delta$pos) + c(-0.5, 0.5)
    title <- "Pairwise m6A difference"
  }
  pairs <- unique(delta[, c("pair", "cluster_A", "cluster_B", "n_A", "n_B")])
  pairs <- pairs[order(pairs$pair), ]
  labels <- setNames(paste0(pairs$cluster_A, " (n=", pairs$n_A, ") ", M6A_DELTA_MINUS, " ",
                            pairs$cluster_B, " (n=", pairs$n_B, ")"), as.character(pairs$pair))
  delta <- delta[order(delta$pair, delta$x), ]
  # raw per-bp line first, then the fill on top of it so dense spikes do not hide the shading
  p <- ggplot2::ggplot(delta, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = delta), color = "grey75", linewidth = .25)
  if (shade) {
    fill <- do.call(rbind, lapply(split(delta, delta$pair), function(d) {
      z <- m6a_delta_zero_crossings(d$x, d$delta_smooth)
      data.frame(pair = d$pair[1], x = z$x, y = z$y)
    }))
    p <- p +
      ggplot2::geom_ribbon(data = fill, ggplot2::aes(x = x, ymin = 0, ymax = pmax(y, 0)),
        fill = above_color, alpha = .45, na.rm = TRUE) +
      ggplot2::geom_ribbon(data = fill, ggplot2::aes(x = x, ymin = pmin(y, 0), ymax = 0),
        fill = below_color, alpha = .45, na.rm = TRUE)
  }
  p <- p +
    ggplot2::geom_hline(yintercept = 0, color = "grey40", linewidth = .4) +
    ggplot2::geom_line(ggplot2::aes(y = delta_smooth), color = "black", linewidth = .7, na.rm = TRUE) +
    ggplot2::facet_wrap(~pair, ncol = 1, labeller = ggplot2::as_labeller(labels)) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::scale_y_continuous(breaks = c(ylim[1], 0, ylim[2])) +
    ggplot2::labs(x = x_label, y = paste0("delta m6A (A ", M6A_DELTA_MINUS, " B)"), title = title,
      subtitle = paste0("delta = A ", M6A_DELTA_MINUS, " B; positive = A more methylated. ",
        "Grey: per-bp delta; black: centred ", smooth_bp, " bp rolling mean")) +
    cowplot::theme_cowplot(font_size = 10) + cowplot::panel_border() +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0))
  if (!is.null(region))
    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = .4)
  list(plot = p, height = 1.6 * nrow(pairs) + 1.5)
}

run_m6a_delta <- function(result_path, output_root, baseline = NULL, smooth_bp = 25,
                          ylim = c(-1, 1), shade = TRUE, width = 10) {
  res <- readRDS(result_path)
  delta <- m6a_delta_table(res, baseline, smooth_bp)
  if (is.null(delta)) return(invisible(NULL))
  region_id <- if (!is.null(res$region$region_id)) res$region$region_id else basename(dirname(dirname(result_path)))
  table_dir <- file.path(output_root, region_id, "tables")
  plot_dir <- file.path(output_root, region_id, "plots")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  tsv <- file.path(table_dir, "pairwise_m6a_delta.tsv")
  pdf <- file.path(plot_dir, "pairwise_m6a_delta.pdf")
  utils::write.table(delta, tsv, sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  p <- plot_m6a_delta(delta, res$region, ylim = ylim, shade = shade, smooth_bp = smooth_bp)
  ggplot2::ggsave(pdf, p$plot, width = width, height = p$height, device = grDevices::cairo_pdf,
                  limitsize = FALSE, bg = "white")
  invisible(list(table = tsv, plot = pdf, delta = delta))
}

run_all_m6a_delta <- function(result_paths, output_root, ...) {
  stopifnot(all(file.exists(result_paths)))
  out <- lapply(result_paths, run_m6a_delta, output_root = output_root, ...)
  names(out) <- result_paths
  invisible(out)
}
