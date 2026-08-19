#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Aggregate all nucleosome midpoint calls in the TSS window for the Pol II
# pausing-index analysis. This mirrors the current nearest-nucleosome workflow's
# gene groups and output naming, but uses every nucleosome assigned to each
# strand-aware TSS-relative bin.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit)) sub(prefix, "", hit[[length(hit)]], fixed = TRUE) else default
}

as_bool <- function(x) {
  if (is.null(x)) return(FALSE)
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

sample_name <- get_arg("sample", "AL10_bc2178_19130")
nuc_source <- toupper(get_arg("nuc-source", "CAGE"))
feature_name <- get_arg("feature", "Nucleosome")
window_bp <- as.integer(get_arg("window-bp", "2000"))
bin_size <- as.integer(get_arg("bin-size", "10"))
low_q <- as.numeric(get_arg("low-q", "0.25"))
high_q <- as.numeric(get_arg("high-q", "0.75"))
force_rebuild <- as_bool(get_arg("force-rebuild", "false"))

repo_dir <- get_arg("repo-dir", "/project/spott/cshan/fiber-seq")
polii_annotation_dir <- get_arg(
  "polii-annotation-dir",
  file.path(repo_dir, "results/PolII/annotations")
)
aggregate_nuc_dir <- get_arg(
  "aggregate-nuc-dir",
  file.path(repo_dir, "results/nuc/all_chr_aggregated_FIRE_features", sample_name, "all_chr")
)
out_dir <- get_arg(
  "out-dir",
  file.path(repo_dir, "results/PolII/promoter_pausing_nucleosome_positioning")
)
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (nuc_source %in% c("", "ALL", "ANY", "NONE", "NULL")) {
  nuc_source_filter <- NULL
  nuc_source_label <- "ALL"
} else {
  nuc_source_filter <- nuc_source
  nuc_source_label <- nuc_source
}

message("Sample: ", sample_name)
message("Selected TSS source: ", nuc_source_label)
message("Output table directory: ", table_dir)
message("Output plot directory: ", plot_dir)

cohen_d <- function(high, low) {
  high <- as.numeric(na.omit(high))
  low <- as.numeric(na.omit(low))
  n_h <- length(high)
  n_l <- length(low)
  if (n_h < 2L || n_l < 2L) return(NA_real_)
  sd_h <- stats::sd(high)
  sd_l <- stats::sd(low)
  pooled <- sqrt(((n_h - 1) * sd_h^2 + (n_l - 1) * sd_l^2) / (n_h + n_l - 2))
  if (is.na(pooled) || pooled == 0) NA_real_ else (mean(high) - mean(low)) / pooled
}

se_mean <- function(x) {
  x <- as.numeric(x)
  stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
}

base_gene_id <- function(x) sub("\\.\\d+$", "", as.character(x))

parse_bed12_nucleosomes <- function(fields) {
  chrom_start <- as.integer(fields[[2]])
  read_name <- fields[[4]]
  block_sizes <- as.integer(strsplit(sub(",$", "", fields[[11]]), ",", fixed = TRUE)[[1]])
  block_starts <- as.integer(strsplit(sub(",$", "", fields[[12]]), ",", fixed = TRUE)[[1]])
  keep <- !is.na(block_sizes) & !is.na(block_starts)
  block_sizes <- block_sizes[keep]
  block_starts <- block_starts[keep]
  if (length(block_sizes) >= 2L) {
    idx <- seq_along(block_sizes)[-c(1L, length(block_sizes))]
    block_sizes <- block_sizes[idx]
    block_starts <- block_starts[idx]
  }
  if (!length(block_sizes)) {
    return(data.table(RID = character(), nuc_start = integer(), nuc_end = integer(),
                      nuc_mid = integer(), nuc_size = integer()))
  }
  nuc_start <- chrom_start + block_starts
  nuc_end <- nuc_start + block_sizes
  data.table(
    RID = read_name,
    nuc_start = nuc_start,
    nuc_end = nuc_end,
    nuc_mid = as.integer(round((nuc_start + nuc_end) / 2)),
    nuc_size = block_sizes
  )
}

read_fiber_bed12 <- function(path) {
  fread(
    path,
    header = FALSE,
    col.names = c(
      "chrom", "chromStart", "chromEnd", "name", "score", "strand",
      "thickStart", "thickEnd", "itemRgb", "blockCount",
      "blockSizes", "blockStarts"
    )
  )
}

expand_bed12_nucleosome_midpoints <- function(bed12) {
  bed12[, {
    fields <- c(
      chrom, chromStart, chromEnd, name, score, strand,
      thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
    )
    parsed <- parse_bed12_nucleosomes(fields)
    if (!nrow(parsed)) return(NULL)
    parsed[, .(RID, nuc_start, nuc_end, nuc_mid, nuc_size)]
  }, by = .(chrom, chromStart, chromEnd, name, strand)]
}

prepare_pausing_groups <- function(pausing_path) {
  pausing <- fread(pausing_path)
  pausing[, PI := as.numeric(PI)]
  pausing[, gene_id_base := base_gene_id(gene_id)]
  pausing <- pausing[
    !is.na(PI) &
      !is.na(gene_id_base) &
      grepl("^chr([0-9]+|X|Y)$", as.character(chrom))
  ]
  setorder(pausing, -PI)
  pausing <- unique(pausing, by = "gene_id_base")

  if (!"cage_tss" %in% names(pausing)) pausing[, cage_tss := NA_real_]
  if (!"tss_gencode" %in% names(pausing)) pausing[, tss_gencode := tss]
  pausing[, cage_tss := as.numeric(cage_tss)]
  pausing[, tss_gencode := as.numeric(tss_gencode)]
  pausing[, tss := as.numeric(tss)]
  pausing[, tss_source := fifelse(!is.na(cage_tss), "CAGE", "GENCODE")]

  if (identical(nuc_source_filter, "CAGE")) {
    pausing <- pausing[!is.na(cage_tss)]
    pausing[, tss_for_nearest_nuc := as.integer(cage_tss)]
    pausing[, tss_source := "CAGE"]
  } else if (identical(nuc_source_filter, "GENCODE")) {
    pausing <- pausing[!is.na(tss_gencode)]
    pausing[, tss_for_nearest_nuc := as.integer(tss_gencode)]
    pausing[, tss_source := "GENCODE"]
  } else {
    pausing <- pausing[!is.na(tss)]
    pausing[, tss_for_nearest_nuc := as.integer(tss)]
  }

  low_cut <- stats::quantile(pausing$PI, low_q, na.rm = TRUE)
  high_cut <- stats::quantile(pausing$PI, high_q, na.rm = TRUE)
  pausing[, pausing_group := fifelse(
    PI >= high_cut, "High PI",
    fifelse(PI <= low_cut, "Low PI", "Middle PI")
  )]
  pausing_groups <- pausing[pausing_group %in% c("High PI", "Low PI")]

  message(sprintf("PI low cutoff (%.0f%%): %.3f", low_q * 100, low_cut))
  message(sprintf("PI high cutoff (%.0f%%): %.3f", high_q * 100, high_cut))
  print(pausing_groups[, .N, by = pausing_group])
  print(pausing_groups[, .N, by = tss_source])

  group_path <- file.path(table_dir, sprintf(
    "%s_%s_high_low_pausing_gene_groups.tsv",
    sample_name, nuc_source_label
  ))
  fwrite(
    pausing_groups[, .(
      gene_name, gene_id, gene_id_base, chrom, strand, tss_for_nearest_nuc,
      tss_source, cage_tss, tss_gencode, PI, pausing_group
    )],
    group_path,
    sep = "\t"
  )
  message("Saved: ", group_path)
  pausing_groups
}

build_all_nucleosome_gene_bin_signal <- function(pausing_groups) {
  gene_bin_path <- file.path(
    aggregate_nuc_dir,
    sprintf("%s_all_chr_gene_bin_signal.tsv.gz", sample_name)
  )
  gene_tss_summary_path <- file.path(
    aggregate_nuc_dir,
    sprintf("%s_all_chr_gene_tss_summary.tsv", sample_name)
  )
  if (!file.exists(gene_bin_path)) stop("Missing aggregate nucleosome table: ", gene_bin_path)
  if (!file.exists(gene_tss_summary_path)) stop("Missing gene TSS summary: ", gene_tss_summary_path)

  gene_tss <- fread(gene_tss_summary_path)
  if (is.null(nuc_source_filter)) {
    gene_tss <- gene_tss[sample == sample_name]
  } else {
    gene_tss <- gene_tss[sample == sample_name & toupper(as.character(source)) == nuc_source_filter]
  }
  gene_tss[, gene_id_base := as.character(gene_id_base)]

  nuc <- fread(gene_bin_path)
  nuc <- nuc[sample == sample_name]
  if (!is.null(nuc_source_filter)) {
    nuc <- nuc[toupper(as.character(source)) == nuc_source_filter]
  }
  nuc <- nuc[feature == feature_name]
  nuc[, gene_id_base := as.character(gene_id_base)]
  nuc[, bin := as.integer(bin)]
  nuc[, signal := fifelse(is.na(as.numeric(signal)), 0, as.numeric(signal))]

  nuc <- merge(
    nuc,
    unique(gene_tss[, .(gene_id_base, n_tss_gene)], by = "gene_id_base"),
    by = "gene_id_base",
    all.x = TRUE
  )
  nuc[is.na(n_tss_gene), n_tss_gene := 1]
  nuc[, signal_per_tss := signal / pmax(as.numeric(n_tss_gene), 1)]

  group_cols <- unique(pausing_groups[, .(gene_id_base, gene_name, PI, pausing_group)], by = "gene_id_base")
  nuc <- merge(
    nuc,
    group_cols,
    by = "gene_id_base",
    all = FALSE,
    allow.cartesian = TRUE
  )

  bins <- seq(-window_bp, window_bp, by = bin_size)
  group_genes <- unique(nuc[, .(gene_id_base, gene_name, PI, pausing_group)], by = "gene_id_base")
  message("Genes with aggregate nucleosome profiles after join:")
  print(group_genes[, .N, by = pausing_group])

  gene_index <- CJ(gene_id_base = group_genes$gene_id_base, bin = bins, unique = TRUE)
  gene_index <- merge(gene_index, group_genes, by = "gene_id_base", all.x = TRUE)
  nuc_signal <- nuc[, .(signal_per_tss = sum(signal_per_tss, na.rm = TRUE)), by = .(gene_id_base, bin)]
  nuc_complete <- merge(gene_index, nuc_signal, by = c("gene_id_base", "bin"), all.x = TRUE)
  nuc_complete[is.na(signal_per_tss), signal_per_tss := 0]
  setcolorder(nuc_complete, c("gene_id_base", "gene_name", "PI", "pausing_group", "bin", "signal_per_tss"))
  setorder(nuc_complete, pausing_group, gene_id_base, bin)
  nuc_complete
}

write_metagene_outputs <- function(nuc_complete) {
  gene_bin_out <- file.path(table_dir, sprintf(
    "%s_%s_nucleosome_gene_bin_signal_high_low_PI.tsv.gz",
    sample_name, nuc_source_label
  ))
  fwrite(nuc_complete, gene_bin_out, sep = "\t")
  message("Saved: ", gene_bin_out)

  nuc_meta <- nuc_complete[, .(
    mean_signal_per_tss = mean(signal_per_tss, na.rm = TRUE),
    median_signal_per_tss = stats::median(signal_per_tss, na.rm = TRUE),
    sem_signal_per_tss = se_mean(signal_per_tss),
    n_genes = uniqueN(gene_id_base)
  ), by = .(pausing_group, bin)]
  setorder(nuc_meta, pausing_group, bin)

  meta_out <- file.path(table_dir, sprintf(
    "%s_%s_nucleosome_metagene_high_low_PI.tsv",
    sample_name, nuc_source_label
  ))
  fwrite(nuc_meta, meta_out, sep = "\t")
  message("Saved: ", meta_out)

  effect <- nuc_complete[, {
    high <- signal_per_tss[pausing_group == "High PI"]
    low <- signal_per_tss[pausing_group == "Low PI"]
    .(
      high_mean = mean(high, na.rm = TRUE),
      low_mean = mean(low, na.rm = TRUE),
      mean_diff_high_minus_low = mean(high, na.rm = TRUE) - mean(low, na.rm = TRUE),
      cohen_d_high_minus_low = cohen_d(high, low),
      n_high = sum(pausing_group == "High PI"),
      n_low = sum(pausing_group == "Low PI")
    )
  }, by = bin]
  setorder(effect, bin)

  effect_out <- file.path(table_dir, sprintf(
    "%s_%s_nucleosome_effect_size_high_low_PI.tsv",
    sample_name, nuc_source_label
  ))
  fwrite(effect, effect_out, sep = "\t")
  message("Saved: ", effect_out)

  list(nuc_meta = nuc_meta, effect = effect)
}

plot_metagene <- function(nuc_meta) {
  colors <- c("High PI" = "#D55E00", "Low PI" = "#0072B2")
  p <- ggplot(nuc_meta, aes(x = bin, y = mean_signal_per_tss, color = pausing_group, fill = pausing_group)) +
    geom_ribbon(
      aes(ymin = mean_signal_per_tss - sem_signal_per_tss,
          ymax = mean_signal_per_tss + sem_signal_per_tss),
      alpha = 0.15,
      linewidth = 0,
      color = NA
    ) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "black") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    coord_cartesian(xlim = c(-1000, 1000)) +
    labs(
      x = "Distance to TSS (bp, transcript strand-aware)",
      y = "Mean nucleosome midpoint count per TSS",
      color = NULL,
      fill = NULL,
      title = sprintf("All nucleosome midpoints in TSS window (%s)", nuc_source_label)
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

  out <- file.path(plot_dir, sprintf(
    "%s_%s_nucleosome_metagene_high_low_PI.pdf",
    sample_name, nuc_source_label
  ))
  ggsave(out, p, width = 8.5, height = 4.5, units = "in")
  message("Saved: ", out)
}

plot_effect <- function(effect) {
  effect_long <- rbindlist(list(
    effect[, .(bin, metric = "Mean difference", value = mean_diff_high_minus_low)],
    effect[, .(bin, metric = "Cohen's d", value = cohen_d_high_minus_low)]
  ))
  effect_long[, metric := factor(metric, levels = c("Mean difference", "Cohen's d"))]

  p <- ggplot(effect_long, aes(x = bin, y = value)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "black") +
    geom_line(color = "#3B3B3B", linewidth = 0.8) +
    facet_wrap(~ metric, ncol = 1, scales = "free_y") +
    coord_cartesian(xlim = c(-1000, 1000)) +
    labs(
      x = "Distance to TSS (bp, transcript strand-aware)",
      y = NULL,
      title = "High PI - Low PI all-nucleosome signal"
    ) +
    theme_bw(base_size = 11)

  out <- file.path(plot_dir, sprintf(
    "%s_%s_nucleosome_effect_size_high_low_PI.pdf",
    sample_name, nuc_source_label
  ))
  ggsave(out, p, width = 8.5, height = 5.5, units = "in")
  message("Saved: ", out)
}

pausing_path <- file.path(polii_annotation_dir, "pausing_index_principal_with_CAGE_TSS_all_genes.tsv")
if (!file.exists(pausing_path)) stop("Missing pausing table: ", pausing_path)

gene_bin_out <- file.path(table_dir, sprintf(
  "%s_%s_nucleosome_gene_bin_signal_high_low_PI.tsv.gz",
  sample_name, nuc_source_label
))

pausing_groups <- prepare_pausing_groups(pausing_path)

if (file.exists(gene_bin_out) && !force_rebuild) {
  message("Using existing all-nucleosome gene-bin table from results/PolII: ", gene_bin_out)
  nuc_complete <- fread(gene_bin_out)
  nuc_complete[, gene_id_base := as.character(gene_id_base)]
  nuc_complete[, bin := as.integer(bin)]
  nuc_complete[, signal_per_tss := as.numeric(signal_per_tss)]
  keep_cols <- setdiff(names(nuc_complete), c("gene_name", "PI", "pausing_group"))
  nuc_complete <- merge(
    nuc_complete[, ..keep_cols],
    unique(pausing_groups[, .(gene_id_base, gene_name, PI, pausing_group)], by = "gene_id_base"),
    by = "gene_id_base",
    all = FALSE
  )
  setcolorder(nuc_complete, c("gene_id_base", "gene_name", "PI", "pausing_group", "bin", "signal_per_tss"))
} else {
  message("Building all-nucleosome gene-bin table from aggregate FIRE feature output.")
  nuc_complete <- build_all_nucleosome_gene_bin_signal(pausing_groups)
}

outputs <- write_metagene_outputs(nuc_complete)
plot_metagene(outputs$nuc_meta)
plot_effect(outputs$effect)

message("Done.")
