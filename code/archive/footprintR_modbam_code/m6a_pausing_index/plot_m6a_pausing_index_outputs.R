#!/usr/bin/env Rscript

# Make extra plots from an existing m6A pausing-index output directory.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

get_env_chr <- function(name, default) {
  value <- trimws(Sys.getenv(name, unset = default))
  if (!nzchar(value)) default else value
}

get_env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value)) default else value
}

assert_file_exists <- function(path, label = "file") {
  if (!file.exists(path)) stop(sprintf("Missing %s: %s", label, path), call. = FALSE)
  invisible(path)
}

infer_latest_output_dir <- function(output_root, sample_name) {
  candidates <- list.dirs(output_root, recursive = FALSE, full.names = TRUE)
  candidates <- candidates[grepl(paste0("^", sample_name, "_"), basename(candidates))]
  if (!length(candidates)) {
    stop("No output directories found under: ", output_root, call. = FALSE)
  }
  candidates[which.max(file.info(candidates)$mtime)]
}

summarize_gene_bin_m6a_from_aligned <- function(aligned_dt, gene_summary, window_bp, bin_size) {
  bins <- data.table::data.table(bin = seq(-window_bp, window_bp, by = bin_size))
  genes <- unique(gene_summary[, .(
    tss_uid, gene_id, gene_id_base, gene_name, chromosome, tss_coordinate,
    tss_strand, tss_source, PI, pausing_group
  )])
  genes[, join_key := 1L]
  bins[, join_key := 1L]
  gene_bins <- merge(genes, bins, by = "join_key", allow.cartesian = TRUE)
  gene_bins[, join_key := NULL]

  tmp <- copy(aligned_dt)
  tmp[, bin := round(rel_pos / bin_size) * bin_size]
  tmp <- tmp[bin >= -window_bp & bin <= window_bp]

  observed <- tmp[, {
    total <- sum(total_calls)
    modified <- sum(modified_calls)
    .(
      m6A_modified_calls = modified,
      m6A_unmodified_calls = sum(unmodified_calls),
      m6A_total_calls = total,
      m6A_fraction_modified = if (total > 0) modified / total else NA_real_
    )
  }, by = .(
    tss_uid, gene_id, gene_id_base, gene_name, chromosome, tss_coordinate,
    tss_strand, tss_source, PI, pausing_group, bin
  )]

  out <- merge(
    gene_bins,
    observed,
    by = c(
      "tss_uid", "gene_id", "gene_id_base", "gene_name", "chromosome",
      "tss_coordinate", "tss_strand", "tss_source", "PI", "pausing_group", "bin"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  out[is.na(m6A_modified_calls), m6A_modified_calls := 0L]
  out[is.na(m6A_unmodified_calls), m6A_unmodified_calls := 0L]
  out[is.na(m6A_total_calls), m6A_total_calls := 0L]
  out[]
}

summarize_gene_bin_histogram <- function(gene_bin_dt, max_modified_bin = 10L) {
  hist_dt <- copy(gene_bin_dt)
  hist_dt[, modified_call_bin := data.table::fifelse(
    m6A_modified_calls >= max_modified_bin,
    paste0(max_modified_bin, "+"),
    as.character(m6A_modified_calls)
  )]
  hist_dt[, modified_call_bin := factor(
    modified_call_bin,
    levels = c(as.character(seq.int(0L, max_modified_bin - 1L)), paste0(max_modified_bin, "+"))
  )]
  hist_dt[, .(
    n_genes = uniqueN(gene_id_base),
    n_tss = uniqueN(tss_uid)
  ), by = .(pausing_group, bin, modified_call_bin)]
}

plot_gene_m6a_boxplot <- function(gene_summary, output_pdf, sample_name) {
  plot_dt <- copy(gene_summary[!is.na(m6A_fraction_modified)])
  plot_dt[, pausing_group := factor(pausing_group, levels = c("Low PI", "High PI"))]

  p_value <- NA_real_
  if (length(unique(plot_dt$pausing_group)) == 2L) {
    p_value <- stats::wilcox.test(m6A_fraction_modified ~ pausing_group, data = plot_dt)$p.value
  }
  subtitle <- if (is.na(p_value)) NULL else sprintf("Wilcoxon rank-sum p = %.3g", p_value)

  label_dt <- plot_dt[, .(
    y = max(m6A_fraction_modified, na.rm = TRUE),
    label = sprintf("n=%s", uniqueN(gene_id_base))
  ), by = pausing_group]

  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = pausing_group, y = m6A_fraction_modified, fill = pausing_group)) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
    ggplot2::geom_jitter(width = 0.12, height = 0, size = 0.6, alpha = 0.35) +
    ggplot2::geom_text(data = label_dt, ggplot2::aes(y = y, label = label), vjust = -0.6, size = 3.2, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c("Low PI" = "#0072B2", "High PI" = "#D55E00")) +
    ggplot2::labs(
      x = NULL,
      y = "Gene-level m6A fraction modified",
      title = sprintf("%s gene-level m6A by Pol II pausing index", sample_name),
      subtitle = subtitle
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none")

  ggplot2::ggsave(output_pdf, p, width = 5.5, height = 5)
}

plot_gene_bin_histogram <- function(hist_dt, output_pdf, sample_name, window_bp) {
  plot_dt <- copy(hist_dt)
  plot_dt[, pausing_group := factor(pausing_group, levels = c("Low PI", "High PI"))]

  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = bin, y = modified_call_bin, fill = n_genes)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~pausing_group, ncol = 1) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.35) +
    ggplot2::scale_fill_viridis_c(option = "magma", direction = -1) +
    ggplot2::scale_x_continuous(
      breaks = seq(-window_bp, window_bp, by = max(500L, window_bp / 5L)),
      labels = function(x) paste0(x / 1000, " kb")
    ) +
    ggplot2::labs(
      x = "Distance to TSS",
      y = "Modified m6A calls per gene per bin",
      fill = "Genes",
      title = sprintf("%s m6A modified-call distribution across TSS window", sample_name)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  ggplot2::ggsave(output_pdf, p, width = 9, height = 7)
}

sample_name <- get_env_chr("SAMPLE", "AL10_bc2178_19130")
output_root <- get_env_chr("OUTPUT_ROOT", "/project/spott/cshan/fiber-seq/results/PolII/m6a_pausing_index")
output_dir <- get_env_chr("OUTPUT_DIR", infer_latest_output_dir(output_root, sample_name))
window_bp <- get_env_int("WINDOW_BP", 1000L)
bin_size <- get_env_int("BIN_SIZE", 10L)
max_modified_hist_bin <- get_env_int("MAX_MODIFIED_HIST_BIN", 10L)

table_dir <- file.path(output_dir, "tables")
plot_dir <- file.path(output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

gene_summary_path <- file.path(table_dir, sprintf("%s_gene_m6A_summary_high_low_PI.tsv", sample_name))
aligned_path <- file.path(table_dir, sprintf("%s_m6A_tss_aligned_positions.tsv.gz", sample_name))
gene_bin_summary_path <- file.path(table_dir, sprintf("%s_gene_bin_m6A_summary_high_low_PI.tsv.gz", sample_name))
histogram_path <- file.path(table_dir, sprintf("%s_gene_bin_m6A_modified_call_histogram_high_low_PI.tsv", sample_name))

assert_file_exists(gene_summary_path, "gene summary table")
gene_summary <- data.table::fread(gene_summary_path)

if (file.exists(gene_bin_summary_path)) {
  gene_bin_summary <- data.table::fread(gene_bin_summary_path)
} else {
  assert_file_exists(aligned_path, "aligned m6A table")
  aligned_dt <- data.table::fread(aligned_path)
  gene_bin_summary <- summarize_gene_bin_m6a_from_aligned(aligned_dt, gene_summary, window_bp, bin_size)
  data.table::fwrite(gene_bin_summary[order(pausing_group, gene_id_base, bin)], gene_bin_summary_path, sep = "\t")
}

gene_bin_histogram <- summarize_gene_bin_histogram(gene_bin_summary, max_modified_hist_bin)
data.table::fwrite(gene_bin_histogram[order(pausing_group, bin, modified_call_bin)], histogram_path, sep = "\t")

boxplot_pdf <- file.path(plot_dir, sprintf("%s_gene_m6A_fraction_boxplot_high_low_PI.pdf", sample_name))
histogram_pdf <- file.path(plot_dir, sprintf("%s_gene_bin_m6A_modified_call_histogram_high_low_PI.pdf", sample_name))

plot_gene_m6a_boxplot(gene_summary, boxplot_pdf, sample_name)
plot_gene_bin_histogram(gene_bin_histogram, histogram_pdf, sample_name, window_bp)

message("Read output directory: ", output_dir)
message("Wrote boxplot: ", boxplot_pdf)
message("Wrote gene-bin summary: ", gene_bin_summary_path)
message("Wrote histogram table: ", histogram_path)
message("Wrote histogram plot: ", histogram_pdf)
