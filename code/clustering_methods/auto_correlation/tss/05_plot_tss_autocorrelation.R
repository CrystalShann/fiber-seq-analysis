#!/usr/bin/env Rscript

# Base R only. Rasterization changes rendering, never the per-base ACF matrix.
args <- commandArgs(trailingOnly = TRUE)
out_dir <- "/project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss"
write_png <- FALSE
force <- FALSE
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--out-dir") {
    i <- i + 1L
    if (i > length(args)) stop("--out-dir requires a path")
    out_dir <- args[i]
  } else if (args[i] == "--png") {
    write_png <- TRUE
  } else if (args[i] == "--force") {
    force <- TRUE
  } else if (args[i] %in% c("--help", "-h")) {
    cat("Usage: Rscript 05_plot_tss_autocorrelation.R [--out-dir PATH] [--png] [--force]\n")
    quit(save = "no", status = 0L)
  } else {
    stop("Unknown argument: ", args[i])
  }
  i <- i + 1L
}
out_dir <- normalizePath(out_dir, mustWork = TRUE)
table_dir <- file.path(out_dir, "tables")
plot_dir <- file.path(out_dir, "plots")
validation_dir <- file.path(out_dir, "validation")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
input_paths <- file.path(table_dir, c("plot_metadata.tsv", "acf_heatmap.tsv.gz",
  "cluster_acf.tsv", "expression_acf.tsv", "cluster_composition.tsv"))
if (!all(file.exists(input_paths))) stop("Missing inputs: ",
  paste(input_paths[!file.exists(input_paths)], collapse = ", "))
signature <- list(input_md5 = tools::md5sum(c(input_paths, script_path)),
  png = write_png, R_version = R.version.string)
manifest_path <- file.path(plot_dir, "plot_manifest.rds")
validation_path <- file.path(validation_dir, "plot_validation.tsv")
old_manifest <- if (file.exists(manifest_path)) {
  tryCatch(readRDS(manifest_path), error = function(e) NULL)
} else NULL
if (!force && !is.null(old_manifest) && identical(signature, old_manifest$signature) &&
    length(old_manifest$output_md5) > 0L &&
    all(file.exists(names(old_manifest$output_md5))) &&
    identical(tools::md5sum(names(old_manifest$output_md5)), old_manifest$output_md5)) {
  message("[", format(Sys.time(), "%F %T"), "] Plot outputs and checksums are valid; skipping.")
  quit(save = "no", status = 0L)
}
message("[", format(Sys.time(), "%F %T"), "] Reading and validating plot tables.")
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("NA", "NaN", "nan", ""), quote = "", comment.char = "")
metadata <- read_tsv(input_paths[1])
cluster_acf <- read_tsv(input_paths[3])
expression_acf <- read_tsv(input_paths[4])
composition <- read_tsv(input_paths[5])
required <- c("row_index", "read_id", "sample", "timepoint", "gene_id", "gene_name",
  "chrom", "tss", "strand", "mean_tpm", "expr_bin", "m6a_count", "cluster",
  "status", "acf_valid", "umap1", "umap2", "heatmap_rank", "window_start", "window_end")
if (!all(required %in% names(metadata))) stop("Missing metadata columns: ",
  paste(setdiff(required, names(metadata)), collapse = ", "))
checks <- data.frame(check = character(), passed = logical(), details = character())
check <- function(name, condition, details) {
  if (!isTRUE(condition)) stop("Plot validation failed [", name, "]: ", details)
  checks <<- rbind(checks, data.frame(check = name, passed = TRUE, details = details))
}
bins <- c("Q1_low", "Q2", "Q3", "Q4_high")
n <- nrow(metadata)
check("nonempty_metadata", n > 0L, paste(n, "selected molecules"))
window_width <- unique(metadata$window_end - metadata$window_start)
check("single_window_width", length(window_width) == 1L && is.finite(window_width) && window_width >= 3,
  paste(window_width, "bp window around the TSS"))
check("unique_physical_reads", !anyNA(metadata$read_id) && !anyDuplicated(metadata$read_id),
  "Each physical read ID appears once")
check("original_metadata_order", identical(as.integer(metadata$row_index), seq_len(n) - 1L),
  "Metadata remains in original zero-based matrix row order")
check("expressed_bins_only", !anyNA(metadata$expr_bin) && all(metadata$expr_bin %in% bins),
  "Only Q1_low, Q2, Q3, Q4_high are included")
bin_counts <- table(factor(metadata$expr_bin, levels = bins))
check("balanced_bins", all(bin_counts > 0L) && length(unique(as.integer(bin_counts))) == 1L,
  paste(names(bin_counts), as.integer(bin_counts), collapse = "; "))
check("valid_m6a_counts", all(is.finite(metadata$m6a_count)) &&
  all(metadata$m6a_count >= 0 & metadata$m6a_count == floor(metadata$m6a_count)),
  "m6A counts are finite nonnegative integers")
valid_text <- tolower(as.character(metadata$acf_valid))
check("valid_acf_flags", all(valid_text %in% c("true", "false", "1", "0")),
  "Every molecule has an explicit ACF validity flag")
valid <- valid_text %in% c("true", "1")
clustered <- !is.na(metadata$status) & metadata$status == "clustered"
metadata$cluster <- as.character(metadata$cluster)
check("invalid_reads_unclustered", all(!is.na(metadata$cluster)) &&
  all(metadata$cluster[!valid] == "Unclustered") &&
  all(metadata$cluster[!clustered] == "Unclustered") &&
  all(metadata$cluster[clustered] != "Unclustered") && all(valid[clustered]),
  paste(sum(!valid), "zero-variance molecules retained as Unclustered"))
cluster_ids <- unique(metadata$cluster[clustered])
cluster_numbers <- suppressWarnings(as.integer(cluster_ids))
check("global_numeric_cluster_ids", length(cluster_ids) > 0L && !anyNA(cluster_numbers) &&
  all(as.character(cluster_numbers) == cluster_ids), "Global Leiden labels are integer IDs")
cluster_levels <- cluster_ids[order(cluster_numbers)]
if (any(!valid) || "Unclustered" %in% composition$cluster)
  cluster_levels <- c(cluster_levels, "Unclustered")
cluster_acf$cluster <- as.character(cluster_acf$cluster)
composition$cluster <- as.character(composition$cluster)
check("heatmap_rank_permutation", !anyNA(metadata$heatmap_rank) &&
  identical(sort(as.integer(metadata$heatmap_rank)), seq_len(n)),
  "Display ranks are a permutation of all selected molecules")
heatmap_order <- order(metadata$heatmap_rank)
cluster_rank <- match(metadata$cluster[heatmap_order], cluster_levels)
check("heatmap_global_cluster_order", all(diff(cluster_rank) >= 0L),
  "Numeric global Leiden labels precede Unclustered")
check("heatmap_m6a_rank", all(vapply(split(metadata$m6a_count[heatmap_order],
  metadata$cluster[heatmap_order]), function(x) all(diff(x) <= 0), logical(1))),
  "m6A counts descend within each global cluster")

# scan() avoids a separate R object for each of up to 2,000 ACF columns.
connection <- gzfile(input_paths[2], open = "rt")
header <- strsplit(readLines(connection, n = 1L), "\t", fixed = TRUE)[[1]]
flat <- tryCatch(scan(connection, what = double(), sep = "\t", quiet = TRUE,
  na.strings = c("NA", "NaN", "nan", "")), finally = close(connection))
check("acf_table_shape", length(header) >= 3L && length(flat) == n * length(header),
  paste(n, "rows and", length(header) - 1L, "lag columns"))
lags <- suppressWarnings(as.integer(sub("^lag_", "", header[-1])))
check("contiguous_unbinned_lags", header[1] == "row_index" && !anyNA(lags) &&
  identical(lags, seq_along(lags) - 1L) && all(header[-1] == paste0("lag_", lags)),
  "Columns are consecutive, unbinned base-pair lags starting at zero")
acf_with_index <- matrix(flat, nrow = n, byrow = TRUE)
rm(flat)
check("acf_metadata_row_order", identical(as.integer(acf_with_index[, 1]),
  as.integer(metadata$row_index)), "ACF and metadata rows match exactly")
acf <- acf_with_index[, -1, drop = FALSE]
rm(acf_with_index)
check("valid_acf_finite", all(is.finite(acf[valid, , drop = FALSE])),
  "All lags are finite for valid molecules")
check("zero_variance_acf_missing", all(is.na(acf[!valid, , drop = FALSE])),
  "Undefined zero-variance ACF rows are missing, shown in gray")
check("acf_lag_zero", all(abs(acf[valid, 1] - 1) < 1e-7),
  "All valid lag-zero ACF values equal one within 1e-7")
check("umap_validity", all(is.finite(metadata$umap1[clustered])) &&
  all(is.finite(metadata$umap2[clustered])), "All clustered molecules have finite UMAP coordinates")
for (item in list(list(data = cluster_acf, group = "cluster", levels = cluster_ids),
                  list(data = expression_acf, group = "expr_bin", levels = bins))) {
  tab <- item$data
  check(paste0(item$group, "_curve_columns"),
    all(c(item$group, "lag_bp", "mean_acf", "median_acf", "n_reads") %in% names(tab)),
    "Aggregate curve columns are present")
  check(paste0(item$group, "_curve_labels"),
    setequal(as.character(tab[[item$group]]), item$levels), "Aggregate curves retain global groups")
  check(paste0(item$group, "_curve_lags"),
    all(vapply(split(tab$lag_bp, tab[[item$group]]), function(x)
      identical(sort(as.integer(x)), lags), logical(1))), "Each curve includes every original lag")
  check(paste0(item$group, "_curve_finite"),
    all(is.finite(tab$n_reads) & tab$n_reads >= 0) &&
    all(is.finite(tab$mean_acf[tab$n_reads > 0])) &&
    all(is.finite(tab$median_acf[tab$n_reads > 0])) &&
    all(is.na(tab$mean_acf[tab$n_reads == 0])) &&
    all(is.na(tab$median_acf[tab$n_reads == 0])),
    "Nonempty aggregate curves are finite; empty groups retain undefined curves")
}
check("composition_columns", all(c("cluster", "expr_bin", "n_reads", "n_bin_total",
  "fraction_within_bin", "fraction_within_cluster") %in% names(composition)),
  "Both directions of expression/cluster composition are available")
check("composition_labels", all(composition$cluster %in% cluster_levels) &&
  all(composition$expr_bin %in% bins) &&
  !anyDuplicated(paste(composition$cluster, composition$expr_bin)),
  "Composition has unique pairs of global Leiden labels and expression bins")
count_matrix <- matrix(0, nrow = length(cluster_levels), ncol = length(bins),
  dimnames = list(cluster_levels, bins))
count_matrix[cbind(match(composition$cluster, cluster_levels), match(composition$expr_bin, bins))] <-
  composition$n_reads
expected_counts <- table(factor(metadata$cluster, levels = cluster_levels),
  factor(metadata$expr_bin, levels = bins))
check("composition_counts", all(count_matrix == expected_counts),
  "Composition includes every selected molecule, including Unclustered")
within_bin <- sweep(count_matrix, 2, colSums(count_matrix), "/")
within_cluster <- sweep(count_matrix, 1, rowSums(count_matrix), "/")
composition_index <- cbind(match(composition$cluster, cluster_levels),
  match(composition$expr_bin, bins))
same_numeric <- function(observed, expected, tolerance = 1e-6) {
  identical(is.na(observed), is.na(expected)) &&
    all(abs(observed[!is.na(expected)] - expected[!is.na(expected)]) < tolerance)
}
check("composition_fractions", same_numeric(composition$fraction_within_bin,
  within_bin[composition_index]) && same_numeric(composition$fraction_within_cluster,
  within_cluster[composition_index]) &&
  all(composition$n_bin_total == as.integer(bin_counts[composition$expr_bin])),
  "Fractions use all selected molecules as the appropriate denominators")

# Plots show clustered molecules only. The tables above keep the Unclustered
# category, so the composition matrix is renormalized over clustered reads here.
plot_levels <- cluster_levels[cluster_levels != "Unclustered"]
cluster_colors <- setNames(grDevices::hcl.colors(length(plot_levels), "Dark 3"), plot_levels)
heatmap_order <- heatmap_order[clustered[heatmap_order]]
plot_counts <- count_matrix[plot_levels, , drop = FALSE]
within_bin <- sweep(plot_counts, 2, colSums(plot_counts), "/")
within_cluster <- within_cluster[plot_levels, , drop = FALSE]
bin_colors <- setNames(c("#4477AA", "#66CCEE", "#EEAA33", "#CC6677"), bins)
acf_palette <- grDevices::colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(257)
acf_limit <- max(0.01, unname(quantile(abs(acf[valid, -1, drop = FALSE]),
  probs = 0.99, na.rm = TRUE)))
output_paths <- character()
save_plot <- function(name, draw, width = 11, height = 8) {
  for (format in c("pdf", if (write_png) "png")) {
    path <- file.path(plot_dir, paste0(name, ".", format))
    tmp <- paste0(path, ".tmp")
    if (format == "pdf") {
      grDevices::pdf(tmp, width = width, height = height, onefile = TRUE,
        useDingbats = FALSE, compress = TRUE)
    } else {
      grDevices::png(tmp, width = width, height = height, units = "in", res = 180,
        type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
    }
    tryCatch(draw(), finally = grDevices::dev.off())
    if (!file.rename(tmp, path)) stop("Could not finalize plot: ", path)
    output_paths <<- c(output_paths, path)
  }
  message("[", format(Sys.time(), "%F %T"), "] Wrote ", name)
}
heatmap_panel <- function(rows, title) {
  par(mar = c(4.4, 6.3, 3.7, 1.1), mgp = c(2.6, 0.7, 0))
  nr <- length(rows)
  values <- acf[rows, , drop = FALSE]
  indices <- as.integer(round((pmax(-acf_limit, pmin(acf_limit, values)) /
    acf_limit + 1) * 128)) + 1L
  colors <- acf_palette[indices]
  raster <- as.raster(matrix(colors, nrow = nr, ncol = length(lags)))
  rm(values, indices, colors)
  plot.new()
  plot.window(xlim = c(-0.5, max(lags) + 0.5), ylim = c(0, nr), xaxs = "i", yaxs = "i")
  rasterImage(raster, -0.5, 0, max(lags) + 0.5, nr, interpolate = FALSE)
  groups <- rle(metadata$cluster[rows])
  bounds <- cumsum(groups$lengths)
  centers <- nr - bounds + groups$lengths / 2
  labels <- paste0("C", groups$values, " (", format(groups$lengths, big.mark = ","), ")")
  axis(1, at = pretty(range(lags), n = 5), cex.axis = 0.8)
  axis(2, at = centers, labels = labels, las = 1, tick = FALSE, cex.axis = 0.7)
  if (length(bounds) > 1L) abline(h = nr - bounds[-length(bounds)], col = "#333333", lwd = 0.45)
  strip_width <- max(1, max(lags)) * 0.015
  rect(-0.5 - strip_width, nr - bounds, -0.5,
    nr - bounds + groups$lengths, col = cluster_colors[groups$values], border = NA, xpd = NA)
  box()
  title(main = title, xlab = "Lag (bp)", cex.main = 1)
  mtext(paste0(format(nr, big.mark = ","), " clustered molecules; descending m6A within each cluster"),
    side = 3, line = 0.25, cex = 0.65)
}
acf_key <- function() {
  par(mar = c(5, 0.5, 4.5, 4.3), mgp = c(2, 0.5, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(-acf_limit, acf_limit), xaxs = "i", yaxs = "i")
  rasterImage(as.raster(matrix(rev(acf_palette), ncol = 1)), 0, -acf_limit,
    0.8, acf_limit, interpolate = FALSE)
  axis(4, at = c(-acf_limit, 0, acf_limit),
    labels = formatC(c(-acf_limit, 0, acf_limit), digits = 3, format = "fg"), las = 1, cex.axis = 0.7)
  title(main = "ACF", cex.main = 0.9)
  mtext("Colors saturate at 99th percentile of |ACF| at lags > 0", side = 4, line = 2.7, cex = 0.6)
}
save_plot("01_acf_heatmap_all", function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 0.16))
  heatmap_panel(heatmap_order, "TSS autocorrelograms: global Leiden classes")
  acf_key()
}, width = 12, height = 11)
save_plot("03_acf_heatmaps_by_expression", function() {
  layout(matrix(c(1, 2, 5, 3, 4, 5), nrow = 2, byrow = TRUE), widths = c(1, 1, 0.15))
  for (bin in bins) heatmap_panel(heatmap_order[metadata$expr_bin[heatmap_order] == bin], bin)
  acf_key()
}, width = 16, height = 12)

curve_plot <- function(tab, grouping, levels, colors, title, min_lag = 0) {
  # min_lag = 0: two panels (full ACF, then lags > 0). min_lag > 0: one panel from that lag.
  lower_bounds <- if (min_lag > 0) min_lag else c(0, 1)
  k <- length(lower_bounds)
  layout(matrix(c(seq_len(k), rep(k + 1, k)), nrow = 2, byrow = TRUE), heights = c(1, 0.25))
  for (lower in lower_bounds) {
    par(mar = c(4.3, 4.5, 3.5, 1), mgp = c(2.7, 0.7, 0))
    use <- tab$lag_bp >= lower
    yrange <- range(c(0, tab$mean_acf[use]), finite = TRUE)
    if (diff(yrange) == 0) yrange <- yrange + c(-0.01, 0.01)
    plot(NA, xlim = range(tab$lag_bp[use]), ylim = yrange,
      xlab = "Lag (bp)", ylab = "Mean ACF", main = title, cex.main = 1)
    abline(h = 0, col = "#BBBBBB", lty = 3)
    for (group in levels) {
      sub <- tab[as.character(tab[[grouping]]) == group & use, , drop = FALSE]
      sub <- sub[order(sub$lag_bp), , drop = FALSE]
      lines(sub$lag_bp, sub$mean_acf, col = colors[group], lwd = 1.25)
    }
    mtext(if (lower == 0) "Full ACF, including lag 0" else
      paste0("Lags from ", lower, " bp; shorter lags omitted from this panel"),
      side = 3, line = 0.25, cex = 0.7)
  }
  par(mar = rep(0, 4))
  plot.new()
  labels <- vapply(levels, function(group) {
    number <- unique(tab$n_reads[as.character(tab[[grouping]]) == group])
    paste0(if (grouping == "cluster") paste0("C", group) else group,
      " (n=", paste(number, collapse = ","), ")")
  }, character(1))
  legend("center", legend = labels, col = colors[levels], lwd = 2,
    ncol = min(4L, length(levels)), bty = "n", cex = 0.8)
  mtext("Means use valid ACF rows only; raw per-base curves, without smoothing", side = 1, line = -0.9, cex = 0.7)
}
save_plot("02_mean_acf_by_cluster", function() curve_plot(cluster_acf, "cluster",
  plot_levels, cluster_colors, "Mean ACF by global Leiden class"))
save_plot("07_mean_acf_by_expression", function() curve_plot(expression_acf, "expr_bin",
  bins, bin_colors, "Mean ACF by expression bin", min_lag = 25))

save_plot("04_cluster_fraction_within_expression", function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 0.35))
  par(mar = c(4.5, 5, 3, 1))
  barplot(within_bin, col = cluster_colors[rownames(within_bin)], border = NA,
    ylim = c(0, 1), ylab = "Fraction of clustered molecules within expression bin",
    main = "Global Leiden-class composition across expression bins")
  par(mar = rep(0, 4)); plot.new()
  legend("center", legend = paste0("C", plot_levels), fill = cluster_colors[plot_levels], bty = "n", cex = 0.9)
})
save_plot("04b_expression_fraction_within_cluster", function() {
  par(mar = c(5.5, 5, 4, 1))
  nonempty <- rowSums(plot_counts) > 0L
  barplot(t(within_cluster[nonempty, , drop = FALSE]), col = bin_colors, border = NA, ylim = c(0, 1),
    ylab = "Fraction of molecules within Leiden class", xlab = "Global Leiden class",
    main = "Expression-bin composition of global Leiden classes", las = 2)
  legend("top", legend = bins, fill = bin_colors, horiz = TRUE, bty = "n", inset = c(0, -0.12), xpd = NA)
})
boxplot_panel <- function(rows, title) {
  rows <- rows[clustered[rows]]
  values <- lapply(plot_levels, function(group) metadata$m6a_count[rows][metadata$cluster[rows] == group])
  names(values) <- paste0("C", plot_levels)
  boxplot(values, col = cluster_colors[plot_levels], border = "#333333", las = 2,
    ylab = paste0("m6A count per ", window_width, "-bp window"), main = title, outline = TRUE,
    pch = 16, cex = 0.35, ylim = range(c(0, metadata$m6a_count[clustered])), cex.axis = 0.75)
}
save_plot("05_m6a_by_cluster_and_expression", function() {
  par(mfrow = c(2, 2), mar = c(6.3, 4.6, 2.8, 0.8), mgp = c(2.8, 0.7, 0))
  for (bin in bins) boxplot_panel(which(metadata$expr_bin == bin), bin)
}, width = max(12, length(plot_levels) * 0.75), height = 10)
save_plot("05b_m6a_by_cluster", function() {
  par(mar = c(6.3, 4.6, 3, 1))
  boxplot_panel(seq_len(n), "m6A counts by global Leiden class: all expression bins")
}, width = max(10, length(plot_levels) * 0.6), height = 7)

set.seed(0L) # Reproducible drawing order only; does not affect analysis or row order.
umap_rows <- sample(which(clustered))
umap_panel <- function(point_colors, title) {
  par(mar = c(4.4, 4.4, 4, 1), mgp = c(2.7, 0.7, 0))
  plot(metadata$umap1[umap_rows], metadata$umap2[umap_rows],
    col = adjustcolor(point_colors[umap_rows], alpha.f = 0.65), pch = 16,
    cex = if (n > 2000) 0.4 else 0.7, xlab = "UMAP 1", ylab = "UMAP 2", main = title)
  mtext(paste0(sum(clustered), " clustered molecules"), side = 3, line = 0.3, cex = 0.7)
}
save_plot("06_umap_by_cluster", function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 0.3))
  umap_panel(cluster_colors[metadata$cluster], "UMAP: global Leiden classes")
  par(mar = rep(0, 4)); plot.new()
  legend("center", legend = paste0("C", plot_levels), col = cluster_colors[plot_levels],
    pch = 16, bty = "n", cex = 0.85)
})
save_plot("06b_umap_by_expression", function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 0.3))
  umap_panel(bin_colors[metadata$expr_bin], "UMAP: expression bins")
  par(mar = rep(0, 4)); plot.new()
  legend("center", legend = bins, col = bin_colors, pch = 16, bty = "n")
})
save_plot("06c_umap_by_m6a", function() {
  palette <- grDevices::hcl.colors(256, "YlOrRd", rev = TRUE)
  limits <- range(metadata$m6a_count[clustered])
  denominator <- max(1, diff(limits))
  indices <- pmax(1L, pmin(256L, as.integer(1 + 255 * (metadata$m6a_count - limits[1]) / denominator)))
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 0.18))
  umap_panel(palette[indices], "UMAP: m6A count")
  par(mar = c(5, 0.5, 5, 4))
  plot.new()
  display_limits <- if (diff(limits) > 0) limits else limits + c(-0.5, 0.5)
  plot.window(xlim = c(0, 1), ylim = display_limits, xaxs = "i", yaxs = "i")
  rasterImage(as.raster(matrix(rev(palette), ncol = 1)), 0, display_limits[1], 0.8,
    display_limits[2], interpolate = FALSE)
  axis(4, at = pretty(limits), las = 1, cex.axis = 0.8)
  title(main = "m6A count", cex.main = 0.8)
})
check("plot_completion", length(output_paths) == 11L * (if (write_png) 2L else 1L) &&
  all(file.exists(output_paths)) && all(file.info(output_paths)$size > 0L),
  paste(length(output_paths), "plots exported successfully"))
write.table(checks, paste0(validation_path, ".tmp"), sep = "\t", row.names = FALSE, quote = FALSE)
if (!file.rename(paste0(validation_path, ".tmp"), validation_path)) stop("Could not finalize plot validation")
manifest <- list(signature = signature, output_md5 = tools::md5sum(c(output_paths, validation_path)),
  completed_at = format(Sys.time(), "%F %T %Z"), n_molecules = n,
  n_valid = sum(valid), n_clustered = sum(clustered), n_zero_variance = sum(!valid),
  acf_color_limit = acf_limit)
saveRDS(manifest, paste0(manifest_path, ".tmp"))
if (!file.rename(paste0(manifest_path, ".tmp"), manifest_path)) stop("Could not finalize plot manifest")
message("[", format(Sys.time(), "%F %T"), "] R plotting complete: ", plot_dir)
