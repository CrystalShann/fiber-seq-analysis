# Selected-fiber signals and plots for pooled enhancer Manhattan clusters.
# Coordinates in records are midpoint-relative, 0-based half-open [-500, 500).

collect_enhancer_footprints <- function(assignments, regions, extract_root) {
  a <- data.table::as.data.table(assignments)
  regions <- data.table::as.data.table(regions)
  empty <- data.table::data.table(RID = character(), track = character(),
    start = integer(), end = integer(), original_size = integer())
  out <- list(); part <- 0L
  groups <- unique(a[, .(enhancer_id, sample_name)])
  for (i in seq_len(nrow(groups))) {
    if (i == 1L || i %% 250L == 0L || i == nrow(groups))
      message("Footprint group ", i, "/", nrow(groups))
    enhancer <- groups$enhancer_id[i]; sample <- groups$sample_name[i]
    selected <- a[enhancer_id == enhancer & sample_name == sample]
    region <- regions[match(enhancer, enhancer_id)]
    stopifnot(nrow(region) == 1L, !is.na(region$midpoint0))
    left0 <- region$midpoint0 - 500L; right0 <- region$midpoint0 + 500L
    query <- GenomicRanges::GRanges(region$chr, IRanges::IRanges(left0 + 1L, right0))
    compact <- gsub("_", "", sample, fixed = TRUE)
    for (kind in c("footprint", "tf")) {
      path <- file.path(extract_root, paste0("firehmm_", kind), compact,
        paste0(compact, "_hmm_extracted_", kind, "_", region$chr, ".bed.gz"))
      if (!all(file.exists(c(path, paste0(path, ".tbi"))))) stop("Missing footprint BED/index: ", path)
      lines <- unique(unlist(Rsamtools::scanTabix(path, param = query), use.names = FALSE))
      if (!length(lines)) next
      bed <- data.table::fread(text = paste(lines, collapse = "\n"), header = FALSE,
                              sep = "\t", showProgress = FALSE)
      if (ncol(bed) < 12L) stop("Expected FiberHMM BED12: ", path)
      bed <- bed[as.character(bed[[4L]]) %in% selected$read_id]
      if (!nrow(bed)) next
      counts <- as.integer(bed[[10L]])
      stopifnot(!anyNA(counts), all(counts >= 0L))
      bed <- bed[counts > 0L]
      if (!nrow(bed)) next
      # FiberHMM blocks are all real calls: unlike ft m6A, do not drop sentinels.
      sizes <- lapply(strsplit(sub(",$", "", as.character(bed[[11L]])), ",", fixed = TRUE), as.integer)
      offsets <- lapply(strsplit(sub(",$", "", as.character(bed[[12L]])), ",", fixed = TRUE), as.integer)
      counts <- as.integer(bed[[10L]])
      stopifnot(!anyNA(counts), all(lengths(sizes) == counts), all(lengths(offsets) == counts))
      parent <- rep(seq_len(nrow(bed)), counts)
      size <- unlist(sizes, use.names = FALSE); offset <- unlist(offsets, use.names = FALSE)
      start0 <- as.integer(bed[[2L]][parent]) + offset
      stopifnot(!anyNA(size), !anyNA(offset), all(size > 0L), all(offset >= 0L),
                all(start0 + size <= as.integer(bed[[3L]][parent])))
      keep <- start0 < right0 & start0 + size > left0 &
        if (kind == "tf") size < 60L else size > 90L
      if (!any(keep)) next
      part <- part + 1L
      out[[part]] <- data.table::data.table(
        RID = selected$RID[match(as.character(bed[[4L]][parent[keep]]), selected$read_id)],
        track = if (kind == "tf") "TF <60 bp" else "Nucleosome >90 bp",
        start = as.integer(pmax(start0[keep], left0) - region$midpoint0),
        end = as.integer(pmin(start0[keep] + size[keep], right0) - region$midpoint0),
        original_size = as.integer(size[keep]))
    }
  }
  if (!length(out)) return(empty)
  records <- unique(data.table::rbindlist(out))
  stopifnot(!anyNA(records), all(records$RID %in% a$RID),
            all(records$start >= -500L), all(records$end <= 500L), all(records$end > records$start))
  records
}

# Match the macrophage reference's plotting contract, with relative coordinates.
enhancer_manhattan_plot_view <- function(result) {
  if (is.null(result$feat_mat) || is.null(result$footprints))
    stop("This result predates saved plot data. Rerun the Manhattan array with k=100.")
  a <- data.table::copy(data.table::as.data.table(result$assignments))
  levels <- rownames(result$m6a_profiles)
  a[, cluster := factor(cluster, levels = levels)]
  a[, sample_name := factor(sample_name, levels = names(LEIDEN_TIMEPOINT_COLORS))]
  data.table::setorder(a, cluster, sample_name, enhancer_id, read_id)
  # The reference orders by start; equal relative starts preserve the order above.
  a[, start := -500L]
  stopifnot(setequal(rownames(result$feat_mat), a$RID), ncol(result$feat_mat) == 1000L)
  list(assignments = as.data.frame(a), feat_mat = as.matrix(result$feat_mat[a$RID, , drop = FALSE]),
       profiles = result$m6a_profiles, params = list(window_size = 1L),
       n_clusters = length(levels))
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

plot_enhancer_manhattan_footprints <- function(result, cluster_id, class_id, max_reads = 100L, seed = 1L) {
  a <- data.table::copy(data.table::as.data.table(result$assignments))[cluster == cluster_id]
  n_total <- nrow(a)
  data.table::setorder(a, sample_name, enhancer_id, read_id)
  set.seed(seed)
  if (nrow(a) > max_reads) a <- a[sort(sample.int(nrow(a), max_reads))]
  a[, sample_name := factor(sample_name, levels = names(LEIDEN_TIMEPOINT_COLORS))]
  data.table::setorder(a, sample_name, enhancer_id, read_id)
  a[, row := seq_len(.N), by = sample_name]
  r <- data.table::copy(data.table::as.data.table(result$footprints))[RID %in% a$RID]
  r[, `:=`(row = a$row[match(RID, a$RID)], sample_name = a$sample_name[match(RID, a$RID)])]
  r[, track := factor(track, levels = c("Nucleosome >90 bp", "TF <60 bp"))]
  data.table::setorder(r, track)
  calls <- summary(result$feat_mat[a$RID, , drop = FALSE])
  m <- data.table::data.table(position = calls$j - 501L + 0.5,
    row = a$row[calls$i], sample_name = a$sample_name[calls$i])
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = a, ggplot2::aes(x = -500, xend = 500, y = row, yend = row), colour = "grey80", linewidth = 0.25) +
    ggplot2::geom_rect(data = r, ggplot2::aes(xmin = start, xmax = end, ymin = row - .38, ymax = row + .38, fill = track)) +
    ggplot2::geom_segment(data = m, ggplot2::aes(x = position, xend = position, y = row - .3, yend = row + .3), colour = "black", linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = .25) +
    ggplot2::facet_grid(sample_name ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_manual(values = c("Nucleosome >90 bp" = "#6488A3", "TF <60 bp" = "#CC3377"), drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(-500, 500), expand = c(0, 0)) +
    ggplot2::scale_y_reverse() + ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), legend.position = "bottom") +
    ggplot2::labs(title = paste(class_id, cluster_id, "| Manhattan-Leiden k =", result$info$k_neighbors),
      subtitle = paste0(nrow(a), " displayed / ", n_total, " sampled fibers; black marks: m6A"),
      x = "Position relative to enhancer midpoint (bp)", y = "Fibers", fill = NULL)
}

save_enhancer_manhattan_pdfs <- function(result, class_id, plot_dir, table_dir, expected_k = 100L) {
  if (!identical(as.integer(result$info$k_neighbors), as.integer(expected_k)))
    stop("Saved Manhattan clusters use k=", result$info$k_neighbors, "; rerun the array for k=", expected_k)
  if (!nrow(result$m6a_profiles)) stop("No clustered fibers to plot for ", class_id)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  view <- enhancer_manhattan_plot_view(result)
  title <- paste(class_id, "enhancers | Manhattan-Leiden k =", expected_k)
  pdf_path <- function(stem) file.path(plot_dir, paste0("manhattan_", stem, "_", class_id, ".pdf"))
  grDevices::pdf(pdf_path("m6a_heatmap"), width = 8, height = 10)
  tryCatch(ComplexHeatmap::draw(plot_cluster_heatmap(view, main = title)), finally = grDevices::dev.off())
  ggplot2::ggsave(pdf_path("cluster_composition"), plot_cluster_composition(view, main = title), width = 8, height = 8)
  ggplot2::ggsave(pdf_path("cluster_structure"), plot_cluster_structure(view, main = title), width = 9, height = 3)
  p <- plot_met_fraction_lines(view, tss = 0, main = title, smooth_k = 1) +
    ggplot2::labs(x = "Position relative to enhancer midpoint (bp)")
  ggplot2::ggsave(pdf_path("m6a_base_cluster_profiles"), p, width = 8, height = max(4, 1.4 * view$n_clusters + 1), limitsize = FALSE)
  profiles <- enhancer_manhattan_signal_profiles(result)
  data.table::fwrite(profiles, file.path(table_dir, paste0("manhattan_signal_profiles_", class_id, ".tsv")), sep = "\t")
  ggplot2::ggsave(pdf_path("footprint_profiles"), plot_enhancer_manhattan_signals(profiles, class_id),
    width = 9, height = max(4, 1.5 * view$n_clusters + 1), limitsize = FALSE)
  grDevices::pdf(pdf_path("single_fiber_footprints"), width = 11, height = 10)
  tryCatch(for (cl in rownames(result$m6a_profiles))
    print(plot_enhancer_manhattan_footprints(result, cl, class_id)), finally = grDevices::dev.off())
  invisible(profiles)
}
