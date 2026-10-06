#!/usr/bin/env Rscript

# Cross-region sliding-window plots and tests (base R only). Reads the
# sliding tables of the three regions written by 04c_sliding_windows.py and
# writes plots/ and tables/ under <root>/combined[/test]. Runs on its own or
# as the dependent job submitted by 06_run_tss_autocorrelation.sh.
args <- commandArgs(trailingOnly = TRUE)
root <- "/project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss"
out_dir <- NULL
test_mode <- FALSE
n_boot <- 1000L
seed <- 0L
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--root") { i <- i + 1L; root <- args[i] }
  else if (args[i] == "--out-dir") { i <- i + 1L; out_dir <- args[i] }
  else if (args[i] == "--test") test_mode <- TRUE
  else if (args[i] == "--boot") { i <- i + 1L; n_boot <- as.integer(args[i]) }
  else if (args[i] == "--seed") { i <- i + 1L; seed <- as.integer(args[i]) }
  else if (args[i] %in% c("--help", "-h")) {
    cat("Usage: Rscript 08_plot_sliding_combined.R [--root PATH] [--test] [--out-dir PATH] [--boot N] [--seed N]\n")
    quit(save = "no", status = 0L)
  } else stop("Unknown argument: ", args[i])
  i <- i + 1L
}
if (is.na(n_boot) || n_boot < 10L) stop("--boot must be at least 10")
regions <- c(span_2kb = "2000_tss", upstream = "upstream_1000_tss", downstream = "downstream_1000_tss")
suffix <- if (test_mode) "test" else ""
region_dirs <- setNames(file.path(root, regions, suffix), names(regions))
region_dirs <- sub("/$", "", region_dirs)
if (is.null(out_dir)) out_dir <- sub("/$", "", file.path(root, "combined", suffix))
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "validation"), recursive = TRUE, showWarnings = FALSE)
msg <- function(...) message("[", format(Sys.time(), "%F %T"), "] ", ...)
read_tsv <- function(path) {
  tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("NA", "NaN", "nan", ""), quote = "", comment.char = "")
  # Python writes True/False; make the flag a proper logical.
  if ("overlaps_ndr_core" %in% names(tab))
    tab$overlaps_ndr_core <- tolower(as.character(tab$overlaps_ndr_core)) %in% c("true", "1")
  tab
}
checks <- data.frame(check = character(), passed = logical(), details = character())
check <- function(name, condition, details) {
  if (!isTRUE(condition)) stop("Combined validation failed [", name, "]: ", details)
  checks <<- rbind(checks, data.frame(check = name, passed = TRUE, details = details))
}
bins <- c("Q1_low", "Q2", "Q3", "Q4_high")
bin_colors <- setNames(c("#4477AA", "#66CCEE", "#EEAA33", "#CC6677"), bins)
timepoints <- c("LPS_0", "LPS_5", "LPS_10", "LPS_15")
ndr_core <- c(-100, 100)

msg("Reading sliding tables from ", paste(region_dirs, collapse = ", "))
long <- list(); mean_acf <- list(); at_control <- list(); expected_windows <- c(span_2kb = 16L, upstream = 5L, downstream = 5L)
for (r in names(regions)) {
  tab <- read_tsv(file.path(region_dirs[[r]], "tables", "sliding_nrl_long.tsv.gz"))
  sampled <- read_tsv(file.path(region_dirs[[r]], "tables", "sampled_molecules.tsv"))
  check(paste0(r, "_region_label"), all(tab$region == r), paste(r, "rows carry their region label"))
  check(paste0(r, "_read_ids_match_sampled"), setequal(tab$read_id, sampled$read_id) &&
    all(table(tab$read_id) == length(unique(tab$win_start))),
    paste(length(unique(tab$read_id)), "molecules, each in every window"))
  span <- c(unique(sampled$window_offset_start), unique(sampled$window_offset_end))
  wins <- unique(tab[, c("win_start", "win_end", "win_center")])
  width <- unique(wins$win_end - wins$win_start)
  check(paste0(r, "_windows_inside_span"), length(width) == 1L && all(wins$win_start >= span[1]) &&
    all(wins$win_end <= span[2]), paste(nrow(wins), "windows of", width, "bp inside", sprintf("[%d, %d)", span[1], span[2])))
  if (width == 500 && length(unique(diff(sort(wins$win_start)))) == 1L && unique(diff(sort(wins$win_start))) == 100)
    check(paste0(r, "_window_count"), nrow(wins) == expected_windows[[r]],
      paste(nrow(wins), "windows at width 500 / step 100"))
  if (r != "span_2kb") check(paste0(r, "_avoids_ndr_core"), all(wins$win_end <= ndr_core[1] | wins$win_start >= ndr_core[2]),
    "one-sided windows never touch [-100, +100)")
  tab$run <- if (r == "span_2kb") "span_2kb" else "one_sided"
  tab$has_peak <- !is.na(tab$nrl_bp)
  tab$m6a_density <- tab$m6a_count_win / width
  long[[r]] <- tab
  mean_acf[[r]] <- read_tsv(file.path(region_dirs[[r]], "tables", "sliding_mean_acf.tsv.gz"))
  at_control[[r]] <- read_tsv(file.path(region_dirs[[r]], "tables", "sliding_at_control.tsv"))
}
long <- do.call(rbind, long)
mean_acf <- do.call(rbind, mean_acf)
at_control <- do.call(rbind, at_control)
win_width <- unique(long$win_end - long$win_start)
check("single_window_width", length(win_width) == 1L, paste(win_width, "bp sliding windows"))
check("expression_bins", setequal(long$expr_bin, bins) && all(long$timepoint %in% timepoints),
  "all four bins and LPS timepoints present")
check("nrl_lags_within_band", all(long$nrl_bp[long$has_peak] >= 1), "NRL lags are positive")

# Metrics: read-level statistics; CIs come from resampling genes (not reads).
metric_defs <- list(
  nrl = list(label = "Median NRL (bp)", column = "nrl_bp", fun = "median", control = "median_nrl_bp"),
  peak_fraction = list(label = "Fraction of molecules with an NRL peak", column = "has_peak", fun = "mean",
    control = "fraction_with_nrl_peak"),
  decay = list(label = "Median decay length (bp)", column = "decay_length_bp", fun = "median",
    control = "median_decay_length_bp"),
  damping = list(label = "Median damping lag (bp)", column = "damping_lag_bp", fun = "median",
    control = "median_damping_lag_bp"),
  m6a = list(label = "Mean m6A density (calls per bp)", column = "m6a_density", fun = "mean", control = NA))
stat_of <- function(values, fun) {
  values <- values[!is.na(values)]
  if (!length(values)) return(NA_real_)
  if (fun == "median") median(values) else mean(values)
}
gene_bootstrap <- function(d, B) {
  idx <- split(seq_len(nrow(d)), d$gene_id)
  G <- length(idx)
  columns <- lapply(metric_defs, function(m) as.numeric(d[[m$column]]))
  point <- vapply(names(metric_defs), function(k) stat_of(columns[[k]], metric_defs[[k]]$fun), numeric(1))
  draws <- matrix(NA_real_, nrow = B, ncol = length(metric_defs), dimnames = list(NULL, names(metric_defs)))
  for (b in seq_len(B)) {
    rows <- unlist(idx[sample.int(G, G, replace = TRUE)], use.names = FALSE)
    for (k in names(metric_defs)) draws[b, k] <- stat_of(columns[[k]][rows], metric_defs[[k]]$fun)
  }
  ci <- apply(draws, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  data.frame(metric = names(metric_defs), estimate = unname(point), ci_low = ci[1, ], ci_high = ci[2, ],
    n_reads = nrow(d), n_genes = G, n_with_peak = sum(d$has_peak), row.names = NULL)
}
set.seed(seed)
msg("Bootstrapping ", n_boot, " gene resamples per bin x window (pooled and per timepoint)")
groups <- unique(long[, c("run", "region", "expr_bin", "win_start", "win_end", "win_center", "overlaps_ndr_core")])
ci_rows <- list()
for (g in seq_len(nrow(groups))) {
  key <- groups[g, ]
  sel <- long$region == key$region & long$expr_bin == key$expr_bin & long$win_start == key$win_start
  for (tp in c("all", timepoints)) {
    d <- long[sel & (tp == "all" | long$timepoint == tp), , drop = FALSE]
    if (!nrow(d)) next
    ci_rows[[length(ci_rows) + 1L]] <- cbind(key, timepoint = tp, gene_bootstrap(d, n_boot), row.names = NULL)
  }
}
ci <- do.call(rbind, ci_rows)
write.table(ci, file.path(out_dir, "tables", "sliding_bootstrap_ci.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
msg("Bootstrap table: ", nrow(ci), " rows")

# Plot helpers -----------------------------------------------------------------
save_plot <- function(name, draw, width = 11, height = 8) {
  path <- file.path(out_dir, "plots", paste0(name, ".pdf"))
  grDevices::pdf(paste0(path, ".tmp"), width = width, height = height, onefile = TRUE, useDingbats = FALSE)
  tryCatch(draw(), finally = grDevices::dev.off())
  if (!file.rename(paste0(path, ".tmp"), path)) stop("Could not finalize ", path)
  msg("Wrote ", name)
  invisible(path)
}
control_curve <- function(run, metric) {
  if (is.na(metric_defs[[metric]]$control)) return(NULL)
  sub <- at_control[at_control$expr_bin == "all" & (if (run == "span_2kb") at_control$region == "span_2kb"
    else at_control$region != "span_2kb"), ]
  sub <- sub[order(sub$win_center), ]
  data.frame(x = sub$win_center, y = sub[[metric_defs[[metric]]$control]], region = sub$region)
}
draw_segments <- function(x, y, ...) {
  # Lines broken at the one-sided gap: separate upstream (x < 0) and downstream (x > 0) runs.
  for (side in list(x < 0, x > 0)) if (any(side)) lines(x[side], y[side], ...)
}
metric_panel <- function(run, metric, tp, title, log_y = FALSE) {
  sub <- ci[ci$run == run & ci$metric == metric & ci$timepoint == tp, ]
  sub <- sub[order(sub$win_center), ]
  control <- control_curve(run, metric)
  ys <- c(sub$estimate, sub$ci_low, sub$ci_high, control$y)
  ys <- ys[is.finite(ys)]
  if (log_y) ys <- ys[ys > 0]
  par(mar = c(4.3, 4.6, 3.2, 0.8), mgp = c(2.7, 0.7, 0))
  if (!length(ys)) { plot.new(); title(main = title, cex.main = 0.9); text(0.5, 0.5, "No estimates"); return(invisible()) }
  ylim <- range(ys); if (diff(ylim) == 0) ylim <- ylim * c(0.9, 1.1) + c(-1e-6, 1e-6)
  plot(NA, xlim = c(-1000, 1000), ylim = ylim, log = if (log_y) "y" else "", xlab = "Window centre relative to TSS (bp)",
    ylab = metric_defs[[metric]]$label, main = title, cex.main = 0.9, xaxt = "n")
  axis(1, at = seq(-1000, 1000, 250))
  usr <- par("usr"); ybottom <- if (log_y) 10^usr[3] else usr[3]; ytop <- if (log_y) 10^usr[4] else usr[4]
  if (run == "span_2kb") {
    ndr <- sub[sub$overlaps_ndr_core, ]
    if (nrow(ndr)) rect(min(ndr$win_center) - win_width / 2, ybottom, max(ndr$win_center) + win_width / 2, ytop,
      col = "#F2DCDB", border = NA)
  } else {
    rect(ndr_core[1], ybottom, ndr_core[2], ytop, col = "#E6E6E6", border = NA)
  }
  abline(v = 0, lty = 3, col = "#888888")
  for (b in bins) {
    s <- sub[sub$expr_bin == b, ]
    ok <- is.finite(s$ci_low) & is.finite(s$ci_high) & (!log_y | (s$ci_low > 0))
    for (side in list(s$win_center < 0, s$win_center > 0)) {
      use <- ok & side
      if (sum(use) >= 2) polygon(c(s$win_center[use], rev(s$win_center[use])), c(s$ci_low[use], rev(s$ci_high[use])),
        col = adjustcolor(bin_colors[b], 0.18), border = NA)
    }
    if (run == "span_2kb") lines(s$win_center, s$estimate, col = bin_colors[b], lwd = 1.6)
    else draw_segments(s$win_center, s$estimate, col = bin_colors[b], lwd = 1.6)
    points(s$win_center, s$estimate, col = bin_colors[b], pch = 16, cex = 0.7)
  }
  if (!is.null(control) && nrow(control)) {
    if (run == "span_2kb") lines(control$x, control$y, col = "#777777", lwd = 1.4, lty = 2)
    else draw_segments(control$x, control$y, col = "#777777", lwd = 1.4, lty = 2)
  }
  box()
}
metric_legend <- function(with_control) {
  par(mar = rep(0, 4)); plot.new()
  labels <- c(bins, if (with_control) "A/T reference control (all bins)")
  legend("center", legend = labels, col = c(bin_colors, if (with_control) "#777777"), lwd = 2,
    lty = c(rep(1, 4), if (with_control) 2), ncol = length(labels), bty = "n", cex = 0.85)
  mtext(paste0("Points: read-level statistic per ", win_width, " bp window; bands: bootstrap 95% CI over genes. ",
    "Pink: windows overlapping the NDR core [-100, +100); grey: the gap of the one-sided regions."),
    side = 1, line = -1.2, cex = 0.65)
}
metric_figure <- function(metrics, titles, log_flags, by_timepoint = FALSE) {
  function() {
    tps <- if (by_timepoint) timepoints else "all"
    k <- length(metrics)
    panels <- matrix(seq_len(2 * k * length(tps)), nrow = 2 * k, byrow = TRUE)
    layout(rbind(panels, rep(max(panels) + 1L, length(tps))), heights = c(rep(1, 2 * k), 0.22))
    for (m in seq_len(k)) for (run in c("span_2kb", "one_sided")) for (tp in tps) {
      label <- if (run == "span_2kb") "2 kb span [-1000, +1000)" else "one-sided regions, upstream | downstream"
      metric_panel(run, metrics[m], tp, paste0(titles[m], ": ", label, if (by_timepoint) paste0(" (", tp, ")")),
        log_flags[m])
    }
    metric_legend(any(!is.na(unlist(lapply(metrics, function(m) metric_defs[[m]]$control)))))
  }
}
figures <- list(
  list(name = "01_sliding_nrl", metrics = "nrl", titles = "Median NRL", log = FALSE),
  list(name = "02_sliding_peak_fraction", metrics = "peak_fraction", titles = "Fraction with an NRL peak", log = FALSE),
  list(name = "03_sliding_decay_damping", metrics = c("decay", "damping"),
    titles = c("Median decay length", "Median damping lag"), log = c(TRUE, FALSE)),
  list(name = "04_sliding_m6a_density", metrics = "m6a", titles = "m6A density", log = FALSE))
outputs <- character()
for (f in figures) {
  k <- length(f$metrics)
  outputs <- c(outputs, save_plot(f$name, metric_figure(f$metrics, f$titles, f$log), width = 10, height = 4.2 * 2 * k + 1))
  outputs <- c(outputs, save_plot(paste0(f$name, "_by_timepoint"), metric_figure(f$metrics, f$titles, f$log, TRUE),
    width = 18, height = 3.8 * 2 * k + 1))
}

# Mean ACF heatmaps: lags 0..win_width-1 (rows) x window centre (columns) per bin.
acf_palette <- grDevices::colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(257)
heatmap_figure <- function(run) {
  function() {
    sub <- mean_acf[mean_acf$signal == "m6a" & (if (run == "span_2kb") mean_acf$region == "span_2kb" else mean_acf$region != "span_2kb"), ]
    limit <- max(0.005, quantile(abs(sub$mean_acf[sub$lag_bp > 0 & is.finite(sub$mean_acf)]), 0.99))
    layout(matrix(c(1, 2, 5, 3, 4, 5), nrow = 2, byrow = TRUE), widths = c(1, 1, 0.16))
    lags <- sort(unique(sub$lag_bp))
    for (b in bins) {
      s <- sub[sub$expr_bin == b, ]
      centers <- sort(unique(s$win_center))
      par(mar = c(4.4, 4.6, 3.2, 0.8), mgp = c(2.7, 0.7, 0))
      plot(NA, xlim = c(-1000, 1000), ylim = c(-0.5, max(lags) + 0.5), xaxs = "i", yaxs = "i",
        xlab = "Window centre relative to TSS (bp)", ylab = "Lag (bp)",
        main = paste0(b, ": mean ACF per ", win_width, " bp window (", if (run == "span_2kb") "2 kb span" else "one-sided", ")"),
        cex.main = 0.85)
      for (cc in centers) {
        v <- s[s$win_center == cc, ]
        v <- v[order(v$lag_bp), ]
        idx <- as.integer(round((pmax(-limit, pmin(limit, v$mean_acf)) / limit + 1) * 128)) + 1L
        idx[!is.finite(v$mean_acf)] <- NA
        colors <- ifelse(is.na(idx), "#BBBBBB", acf_palette[pmax(1L, pmin(257L, idx))])
        # Column raster: lag 0 at the bottom. The step is the spacing between neighbouring window centres.
        rasterImage(as.raster(matrix(rev(colors), ncol = 1)), cc - 50, -0.5, cc + 50, max(lags) + 0.5, interpolate = FALSE)
      }
      if (run == "span_2kb") {
        ndr <- unique(s$win_center[s$overlaps_ndr_core])
        if (length(ndr)) rect(min(ndr) - 50, -0.5, max(ndr) + 50, max(lags) + 0.5, border = "#B2182B", lwd = 1.2, lty = 2)
      } else rect(ndr_core[1], -0.5, ndr_core[2], max(lags) + 0.5, col = "#E6E6E6", border = NA)
      abline(v = 0, lty = 3, col = "#444444")
      box()
    }
    par(mar = c(5, 0.5, 4.5, 4.3)); plot.new()
    plot.window(xlim = c(0, 1), ylim = c(-limit, limit), xaxs = "i", yaxs = "i")
    rasterImage(as.raster(matrix(rev(acf_palette), ncol = 1)), 0, -limit, 0.8, limit, interpolate = FALSE)
    axis(4, at = c(-limit, 0, limit), labels = formatC(c(-limit, 0, limit), digits = 3, format = "fg"), las = 1, cex.axis = 0.7)
    title(main = "Mean ACF", cex.main = 0.9)
    mtext("Colours saturate at the 99th percentile of |mean ACF| at lags > 0; 100 bp columns per window centre",
      side = 4, line = 2.7, cex = 0.55)
  }
}
outputs <- c(outputs, save_plot("05_sliding_acf_heatmap_2kb", heatmap_figure("span_2kb"), width = 15, height = 11))
outputs <- c(outputs, save_plot("05b_sliding_acf_heatmap_one_sided", heatmap_figure("one_sided"), width = 15, height = 11))

# Statistics -----------------------------------------------------------------
# (a) Upstream vs downstream: different reads, so unpaired Wilcoxon on gene-level medians per bin.
# (b) Within the 2 kb span: mirror windows share molecules, so paired (signed-rank) by read.
gene_level <- function(d) {
  out <- lapply(split(d, d$gene_id), function(g) vapply(names(metric_defs), function(k)
    stat_of(as.numeric(g[[metric_defs[[k]]$column]]), metric_defs[[k]]$fun), numeric(1)))
  as.data.frame(do.call(rbind, out))
}
tests <- list()
safe_wilcox <- function(a, b, paired = FALSE) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (paired) { keep <- is.finite(a) & is.finite(b); a <- a[keep]; b <- b[keep]; if (length(a) < 3 || all(a == b)) return(c(NA, NA, length(a), length(a))) }
  else { a <- a[is.finite(a)]; b <- b[is.finite(b)]; if (length(a) < 3 || length(b) < 3) return(c(NA, NA, length(a), length(b))) }
  w <- suppressWarnings(wilcox.test(a, b, paired = paired, exact = FALSE))
  c(unname(w$statistic), w$p.value, length(a), length(b))
}
for (b in bins) {
  up <- gene_level(long[long$region == "upstream" & long$expr_bin == b, ])
  down <- gene_level(long[long$region == "downstream" & long$expr_bin == b, ])
  for (k in names(metric_defs)) {
    res <- safe_wilcox(up[[k]], down[[k]])
    tests[[length(tests) + 1L]] <- data.frame(test_type = "upstream_vs_downstream_unpaired_gene_medians",
      comparison = "upstream [-1000,-100) vs downstream [+100,+1000)", expr_bin = b, metric = k, unit = "gene",
      n_a = res[3], n_b = res[4], median_a = median(up[[k]], na.rm = TRUE), median_b = median(down[[k]], na.rm = TRUE),
      statistic = res[1], p_value = res[2])
  }
}
span <- long[long$region == "span_2kb", ]
centers <- sort(unique(span$win_center))
mirror <- centers[centers < 0 & -centers %in% centers]
mirror <- mirror[!span$overlaps_ndr_core[match(mirror, span$win_center)] & !span$overlaps_ndr_core[match(-mirror, span$win_center)]]
for (cc in mirror) for (b in bins) {
  a <- span[span$win_center == cc & span$expr_bin == b, ]
  d <- span[span$win_center == -cc & span$expr_bin == b, ]
  d <- d[match(a$read_id, d$read_id), ]
  stopifnot(identical(a$read_id, d$read_id))
  for (k in names(metric_defs)) {
    col <- metric_defs[[k]]$column
    res <- safe_wilcox(a[[col]], d[[col]], paired = TRUE)
    tests[[length(tests) + 1L]] <- data.frame(test_type = "span_2kb_mirror_windows_paired_by_read",
      comparison = sprintf("centre %+d vs %+d", as.integer(cc), as.integer(-cc)), expr_bin = b, metric = k, unit = "read",
      n_a = res[3], n_b = res[4], median_a = median(as.numeric(a[[col]]), na.rm = TRUE),
      median_b = median(as.numeric(d[[col]]), na.rm = TRUE), statistic = res[1], p_value = res[2])
  }
}
tests <- do.call(rbind, tests)
tests$p_adj_bh <- NA_real_
for (family in unique(tests$test_type)) {
  sel <- tests$test_type == family & is.finite(tests$p_value)
  tests$p_adj_bh[sel] <- p.adjust(tests$p_value[sel], method = "BH")
}
write.table(tests, file.path(out_dir, "tables", "sliding_tests.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
summary_all <- ci[ci$timepoint == "all", ]
write.table(summary_all, file.path(out_dir, "tables", "sliding_combined_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
check("plot_completion", length(outputs) == 10L && all(file.exists(outputs)) && all(file.info(outputs)$size > 0),
  paste(length(outputs), "combined plots written"))
check("tests_written", nrow(tests) > 0 && all(c("p_value", "p_adj_bh") %in% names(tests)),
  paste(nrow(tests), "test rows;", sum(is.finite(tests$p_value)), "with a p-value"))
write.table(checks, file.path(out_dir, "validation", "combined_validation.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
msg("Combined sliding-window plots complete: ", out_dir)
