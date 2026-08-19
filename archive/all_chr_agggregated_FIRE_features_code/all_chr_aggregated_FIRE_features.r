source("/project/spott/cshan/fiber-seq/plot_nuc_distribution_functions.R")

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
})

cage_tss_path <- "/project/spott/cshan/annotations/fantom5/fantom5.hg38.LCL.consensus.CAGE_peaks.withGene.bed.gz"
sample_root <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results"
default_sample <- "AL10_bc2178_19130"
default_window_bp <- 2000L
default_bin_size <- 10L

plot_dir <- "/project/spott/cshan/fiber-seq/results/plots/nuc/all_chr_aggregated_FIRE_features"
table_dir <- "/project/spott/cshan/fiber-seq/results/nuc/all_chr_aggregated_FIRE_features"

feature_specs <- list(
  Nucleosome = list(subdir = "nuc_by_chr", stem = "nuc"),
  m6A = list(subdir = "m6a_by_chr", stem = "m6a"),
  CpG = list(subdir = "cpg_by_chr", stem = "cpg"),
  MSP = list(subdir = "msp_by_chr", stem = "msp")
)

feature_levels <- names(feature_specs)
source_levels <- c("CAGE")
feature_colors <- c(
  "Nucleosome" = "#4D4D4D",
  "m6A" = "#7B3294",
  "CpG" = "#8C510A",
  "MSP" = "#1B9E77"
)
source_colors <- c("CAGE" = "#be1c54ff")

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) y else x
}

parse_named_args <- function(args) {
  if (!length(args)) return(list())
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) next
    key <- sub("=.*$", "", arg)
    value <- sub("^[^=]*=", "", arg)
    out[[key]] <- value
  }
  out
}

resolve_run_config <- function() {
  args <- parse_named_args(commandArgs(trailingOnly = TRUE))
  sample_name <- args[["sample"]] %||% Sys.getenv("FIBER_SAMPLE", default_sample)
  window_value <- as.integer(args[["window_bp"]] %||% Sys.getenv("FIBER_WINDOW_BP", as.character(default_window_bp)))
  bin_value <- as.integer(args[["bin_size"]] %||% Sys.getenv("FIBER_BIN_SIZE", as.character(default_bin_size)))
  mode_value <- tolower(args[["mode"]] %||% "full")
  chrom_value <- args[["chrom"]] %||% Sys.getenv("FIBER_CHROM", "")
  chroms_value <- args[["chroms"]] %||% Sys.getenv("FIBER_CHROMS", "")
  chrom_set_value <- tolower(args[["chrom_set"]] %||% Sys.getenv("FIBER_CHROM_SET", "autosomes_xy"))

  default_chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
  chroms_value <- trimws(chroms_value)
  chroms <- if (nzchar(chroms_value)) {
    unique(trimws(unlist(strsplit(chroms_value, ",", fixed = TRUE))))
  } else if (chrom_set_value %in% c("autosomes_xy", "default", "human")) {
    default_chroms
  } else if (chrom_set_value == "all") {
    character()
  } else {
    default_chroms
  }

  list(
    sample = sample_name,
    window_bp = window_value,
    bin_size = bin_value,
    mode = mode_value,
    chrom = chrom_value,
    chroms = chroms,
    chrom_set = chrom_set_value
  )
}

extract_chrom_from_feature_path <- function(path) {
  sub("^.*\\.(chr[^.]+)\\.bed\\.gz$", "\\1", basename(path))
}

annotate_gene_id_base <- function(gr) {
  mcols(gr)$gene_id_base <- sub("\\..*$", "", mcols(gr)$gene_id)
  gr
}

discover_feature_manifest <- function(sample_name) {
  sample_extract_dir <- file.path(sample_root, sample_name, "extracted_results")

  manifest <- rbindlist(lapply(names(feature_specs), function(feature_nm) {
    spec <- feature_specs[[feature_nm]]
    feature_dir <- file.path(sample_extract_dir, spec$subdir)
    if (!dir.exists(feature_dir)) {
      stop("Missing extracted feature directory: ", feature_dir)
    }

    pattern <- sprintf("^%s\\.ft_extracted_%s\\.chr[^.]+\\.bed\\.gz$", sample_name, spec$stem)
    paths <- sort(list.files(feature_dir, pattern = pattern, full.names = TRUE))
    if (!length(paths)) {
      stop("No chromosome-level BED files found for feature ", feature_nm, " in ", feature_dir)
    }

    data.table(
      sample = sample_name,
      feature = feature_nm,
      chromosome = vapply(paths, extract_chrom_from_feature_path, character(1)),
      feature_file = basename(paths),
      feature_path = paths
    )
  }), use.names = TRUE, fill = TRUE)

  manifest[]
}

load_cage_tss_annotations <- function() {
  tss_gr <- load_tss_annotations_bed9_with_gene_id(cage_tss_path)
  tss_gr <- tss_gr[!is.na(mcols(tss_gr)$gene_id) & mcols(tss_gr)$gene_id != "NA"]
  annotate_gene_id_base(tss_gr)
}

split_tss_by_chromosome <- function(tss_gr) {
  split(tss_gr, as.character(seqnames(tss_gr)))
}

expected_chunk_columns <- list(
  feature_manifest = c("sample", "feature", "chromosome", "feature_file", "feature_path"),
  raw_profiles = c("bin", "N", "sample", "chrom_scope", "source", "feature"),
  gene_bin_signal = c(
    "gene_id", "gene_id_base", "bin", "signal", "n_tss_with_signal",
    "sample", "chrom_scope", "source", "feature"
  ),
  gene_tss_summary = c(
    "sample", "chrom_scope", "source", "gene_id", "gene_id_base",
    "n_tss_gene", "n_gene_chromosomes"
  )
)

chunk_output_paths <- function(sample_name, chrom, chunk_dir) {
  chrom_dir <- file.path(chunk_dir, chrom)
  list(
    chrom_dir = chrom_dir,
    feature_manifest = file.path(chrom_dir, sprintf("%s_%s_feature_manifest.tsv", sample_name, chrom)),
    raw_profiles = file.path(chrom_dir, sprintf("%s_%s_raw_profiles.tsv.gz", sample_name, chrom)),
    gene_bin_signal = file.path(chrom_dir, sprintf("%s_%s_gene_bin_signal.tsv.gz", sample_name, chrom)),
    gene_tss_summary = file.path(chrom_dir, sprintf("%s_%s_gene_tss_summary.tsv", sample_name, chrom)),
    partial_rds = file.path(chrom_dir, sprintf("%s_%s_partial_objects.rds", sample_name, chrom))
  )
}

file_is_nonempty <- function(path) {
  file.exists(path) && isTRUE(file.info(path)$size > 0)
}

tsv_has_expected_columns <- function(path, expected_cols) {
  if (!file_is_nonempty(path)) return(FALSE)

  header <- tryCatch(
    names(fread(path, nrows = 0, showProgress = FALSE)),
    error = function(e) character()
  )

  length(header) > 0 && all(expected_cols %in% header)
}

chunk_outputs_are_complete <- function(sample_name, chrom, chunk_dir) {
  paths <- chunk_output_paths(sample_name, chrom, chunk_dir)
  file_is_nonempty(paths$partial_rds) &&
    tsv_has_expected_columns(paths$feature_manifest, expected_chunk_columns$feature_manifest) &&
    tsv_has_expected_columns(paths$raw_profiles, expected_chunk_columns$raw_profiles) &&
    tsv_has_expected_columns(paths$gene_bin_signal, expected_chunk_columns$gene_bin_signal) &&
    tsv_has_expected_columns(paths$gene_tss_summary, expected_chunk_columns$gene_tss_summary)
}

load_cached_chromosome_result <- function(sample_name, chrom, chunk_dir) {
  paths <- chunk_output_paths(sample_name, chrom, chunk_dir)
  readRDS(paths$partial_rds)
}

select_analysis_chroms <- function(feature_manifest, config) {
  available <- sort(unique(feature_manifest$chromosome))
  requested <- config$chroms

  if (identical(config$chrom_set, "all")) {
    requested <- available
  }

  if (is.null(requested) || !length(requested)) {
    requested <- available
  }

  requested <- intersect(requested, available)
  if (!length(requested)) {
    stop("No requested chromosomes were found in the feature manifest.")
  }

  requested
}

aggregate_partial_chrom_results <- function(chrom_results, total_tss_count) {
  meta_raw <- rbindlist(lapply(chrom_results, `[[`, "meta_raw"), use.names = TRUE, fill = TRUE)
  gene_bin_signal <- rbindlist(lapply(chrom_results, `[[`, "gene_bin_signal"), use.names = TRUE, fill = TRUE)

  meta_raw <- meta_raw[, .(N = sum(N)), by = .(sample, source, feature, bin)]
  meta_raw[, chrom_scope := "all_chr"]
  setcolorder(meta_raw, c("sample", "chrom_scope", "source", "feature", "bin", "N"))

  if (total_tss_count > 0) {
    meta_norm <- copy(meta_raw)
    meta_norm[, N := N / total_tss_count]
  } else {
    meta_norm <- copy(meta_raw)
    meta_norm[, N := 0]
  }

  gene_bin_signal <- gene_bin_signal[, .(
    signal = sum(signal),
    n_tss_with_signal = sum(n_tss_with_signal)
  ), by = .(sample, source, feature, gene_id, gene_id_base, bin)]
  gene_bin_signal[, chrom_scope := "all_chr"]
  setcolorder(gene_bin_signal, c("sample", "chrom_scope", "source", "feature", "gene_id", "gene_id_base", "bin", "signal", "n_tss_with_signal"))

  list(
    meta_raw = meta_raw,
    meta_norm = meta_norm,
    gene_bin_signal = gene_bin_signal
  )
}

sum_count_table <- function(existing, incoming, group_cols, value_col) {
  if (is.null(existing) || !nrow(existing)) {
    return(copy(incoming))
  }

  rbindlist(list(existing, incoming), use.names = TRUE, fill = TRUE)[
    , .(value = sum(get(value_col), na.rm = TRUE)),
    by = group_cols
  ][, (value_col) := value][, value := NULL]
}

aggregate_partial_chrom_results_from_paths <- function(partial_paths, total_tss_count) {
  meta_raw_sum <- NULL
  gene_bin_signal_sum <- NULL

  for (path in partial_paths) {
    chrom_result <- readRDS(path)

    meta_raw <- chrom_result$meta_raw[, .(N = sum(N)), by = .(sample, source, feature, bin)]
    gene_bin_signal <- chrom_result$gene_bin_signal[, .(
      signal = sum(signal),
      n_tss_with_signal = sum(n_tss_with_signal)
    ), by = .(sample, source, feature, gene_id, gene_id_base, bin)]

    meta_raw_sum <- sum_count_table(
      meta_raw_sum,
      meta_raw,
      group_cols = c("sample", "source", "feature", "bin"),
      value_col = "N"
    )
    gene_bin_signal_sum <- rbindlist(list(gene_bin_signal_sum, gene_bin_signal), use.names = TRUE, fill = TRUE)[
      , .(
        signal = sum(signal),
        n_tss_with_signal = sum(n_tss_with_signal)
      ),
      by = .(sample, source, feature, gene_id, gene_id_base, bin)
    ]

    rm(chrom_result, meta_raw, gene_bin_signal)
    gc(verbose = FALSE)
  }

  meta_raw_sum[, chrom_scope := "all_chr"]
  setcolorder(meta_raw_sum, c("sample", "chrom_scope", "source", "feature", "bin", "N"))

  if (total_tss_count > 0) {
    meta_norm <- copy(meta_raw_sum)
    meta_norm[, N := N / total_tss_count]
  } else {
    meta_norm <- copy(meta_raw_sum)
    meta_norm[, N := 0]
  }

  gene_bin_signal_sum[, chrom_scope := "all_chr"]
  setcolorder(gene_bin_signal_sum, c("sample", "chrom_scope", "source", "feature", "gene_id", "gene_id_base", "bin", "signal", "n_tss_with_signal"))

  list(
    meta_raw = meta_raw_sum,
    meta_norm = meta_norm,
    gene_bin_signal = gene_bin_signal_sum
  )
}

get_feature_paths_for_chrom <- function(feature_manifest, chrom) {
  paths <- setNames(vector("list", length(feature_levels)), feature_levels)
  for (feature_nm in feature_levels) {
    path <- feature_manifest[chromosome == chrom & feature == feature_nm, feature_path]
    if (length(path) != 1L) {
      stop("Expected exactly one file for ", feature_nm, " on ", chrom, " but found ", length(path))
    }
    paths[[feature_nm]] <- path
  }
  paths
}

read_feature_granges <- function(feature_paths) {
  setNames(
    lapply(names(feature_paths), function(feature_nm) {
      load_fiber_feature_positions(feature_paths[[feature_nm]], feature_name = feature_nm)
    }),
    names(feature_paths)
  )
}

empty_feature_hit_table <- function() {
  data.table(
    gene_id = character(),
    gene_id_base = character(),
    seqnames = character(),
    strand = character(),
    tss_index = integer(),
    feature_pos = integer(),
    tss_pos = integer(),
    rel_pos = integer(),
    bin = integer()
  )
}

find_feature_hits_near_tss_with_gene <- function(position_gr, tss_gr,
                                                 window_bp = 2000L,
                                                 bin_size = 10L) {
  if (!length(position_gr) || !length(tss_gr)) return(empty_feature_hit_table())

  tss_windows <- promoters(tss_gr, upstream = window_bp, downstream = window_bp + 1L)
  hits <- findOverlaps(position_gr, tss_windows)
  if (!length(hits)) return(empty_feature_hit_table())

  query_idx <- queryHits(hits)
  subject_idx <- subjectHits(hits)
  feature_pos <- start(position_gr)[query_idx]
  tss_pos <- start(tss_gr)[subject_idx]
  tss_str <- as.character(strand(tss_gr))[subject_idx]
  rel_pos <- ifelse(tss_str == "+", feature_pos - tss_pos, tss_pos - feature_pos)

  data.table(
    gene_id = as.character(mcols(tss_gr)$gene_id[subject_idx]),
    gene_id_base = as.character(mcols(tss_gr)$gene_id_base[subject_idx]),
    seqnames = as.character(seqnames(tss_gr))[subject_idx],
    strand = tss_str,
    tss_index = subject_idx,
    feature_pos = as.integer(feature_pos),
    tss_pos = as.integer(tss_pos),
    rel_pos = as.integer(rel_pos),
    bin = as.integer(round(rel_pos / bin_size) * bin_size)
  )[bin >= -window_bp & bin <= window_bp]
}

summarize_tss_by_gene <- function(tss_gr, sample_name, scope_label) {
  if (!length(tss_gr)) {
    return(data.table(
      sample = character(),
      chrom_scope = character(),
      source = character(),
      gene_id = character(),
      gene_id_base = character(),
      n_tss_gene = integer(),
      n_gene_chromosomes = integer()
    ))
  }

  data.table(
    sample = sample_name,
    chrom_scope = scope_label,
    source = "CAGE",
    gene_id = as.character(mcols(tss_gr)$gene_id),
    gene_id_base = as.character(mcols(tss_gr)$gene_id_base),
    seqnames = as.character(seqnames(tss_gr))
  )[, .(
    n_tss_gene = .N,
    n_gene_chromosomes = uniqueN(seqnames)
  ), by = .(sample, chrom_scope, source, gene_id, gene_id_base)][]
}

summarize_meta_from_hits <- function(hit_dt, total_tss_count,
                                     sample_name, scope_label,
                                     feature_name, window_bp = 2000L, bin_size = 10L) {
  template <- data.table(bin = seq(-window_bp, window_bp, by = bin_size))
  counts <- if (nrow(hit_dt)) {
    hit_dt[, .(N = .N), by = bin]
  } else {
    data.table(bin = integer(), N = integer())
  }

  meta_raw <- merge(template, counts, by = "bin", all.x = TRUE, sort = TRUE)
  meta_raw[is.na(N), N := 0L]
  meta_raw[, `:=`(
    sample = sample_name,
    chrom_scope = scope_label,
    source = "CAGE",
    feature = feature_name
  )]

  meta_norm <- copy(meta_raw)
  if (total_tss_count > 0) {
    meta_norm[, N := N / total_tss_count]
  } else {
    meta_norm[, N := 0]
  }

  list(raw = meta_raw[], norm = meta_norm[])
}

summarize_gene_bin_signal <- function(hit_dt, sample_name, scope_label, feature_name) {
  if (!nrow(hit_dt)) {
    return(data.table(
      sample = character(),
      chrom_scope = character(),
      source = character(),
      feature = character(),
      gene_id = character(),
      gene_id_base = character(),
      bin = integer(),
      signal = integer(),
      n_tss_with_signal = integer()
    ))
  }

  hit_dt[, .(
    signal = .N,
    n_tss_with_signal = uniqueN(tss_index)
  ), by = .(gene_id, gene_id_base, bin)][
    , `:=`(
      sample = sample_name,
      chrom_scope = scope_label,
      source = "CAGE",
      feature = feature_name
    )
  ][]
}

summarize_gene_signal <- function(gene_bin_signal, gene_tss_summary) {
  if (!nrow(gene_bin_signal)) {
    return(data.table(
      sample = character(),
      chrom_scope = character(),
      source = character(),
      feature = character(),
      gene_id = character(),
      gene_id_base = character(),
      total_signal = numeric(),
      central_signal_250bp = numeric(),
      peak_bin = integer(),
      peak_bin_signal = numeric(),
      bins_with_signal = integer(),
      n_tss_with_signal = integer(),
      n_tss_gene = integer(),
      n_gene_chromosomes = integer(),
      signal_per_tss_gene = numeric(),
      central_signal_fraction = numeric(),
      signal_fraction = numeric(),
      gene_rank = integer(),
      cumulative_signal_fraction = numeric()
    ))
  }

  gene_summary <- gene_bin_signal[
    order(bin),
    .(
      total_signal = sum(signal),
      central_signal_250bp = sum(signal[abs(bin) <= 250L]),
      peak_bin = bin[which.max(signal)][1L],
      peak_bin_signal = max(signal),
      bins_with_signal = .N,
      n_tss_with_signal = max(n_tss_with_signal)
    ),
    by = .(sample, chrom_scope, source, feature, gene_id, gene_id_base)
  ]

  gene_summary <- merge(
    gene_summary,
    gene_tss_summary,
    by = c("sample", "chrom_scope", "source", "gene_id", "gene_id_base"),
    all.x = TRUE,
    sort = FALSE
  )

  gene_summary[is.na(n_tss_gene), n_tss_gene := 0L]
  gene_summary[is.na(n_gene_chromosomes), n_gene_chromosomes := 0L]
  gene_summary[, signal_per_tss_gene := fifelse(n_tss_gene > 0, total_signal / n_tss_gene, 0)]
  gene_summary[, central_signal_fraction := fifelse(total_signal > 0, central_signal_250bp / total_signal, 0)]
  gene_summary <- gene_summary[order(feature, -total_signal, gene_id)]
  gene_summary[, group_signal := sum(total_signal), by = .(sample, chrom_scope, feature, source)]
  gene_summary[, signal_fraction := fifelse(group_signal > 0, total_signal / group_signal, 0)]
  gene_summary[, group_signal := NULL]
  gene_summary[, gene_rank := seq_len(.N), by = .(sample, chrom_scope, feature, source)]
  gene_summary[, cumulative_signal_fraction := cumsum(signal_fraction), by = .(sample, chrom_scope, feature, source)]
  gene_summary[]
}

threshold_rank <- function(cum_frac, threshold_value) {
  idx <- which(cum_frac >= threshold_value)
  if (length(idx)) idx[1L] else NA_integer_
}

summarize_gene_dominance <- function(gene_signal_summary) {
  if (!nrow(gene_signal_summary)) {
    return(data.table(
      sample = character(),
      chrom_scope = character(),
      feature = character(),
      source = character(),
      n_genes = integer(),
      total_signal = numeric(),
      top_gene_id = character(),
      top1_fraction = numeric(),
      top5_fraction = numeric(),
      top10_fraction = numeric(),
      genes_for_25pct = integer(),
      genes_for_50pct = integer(),
      genes_for_75pct = integer(),
      genes_for_90pct = integer(),
      effective_gene_number = numeric()
    ))
  }

  gene_signal_summary[
    order(-total_signal, gene_id),
    {
      shares <- if (sum(total_signal) > 0) total_signal / sum(total_signal) else rep(0, .N)
      cum_frac <- cumsum(shares)
      data.table(
        n_genes = .N,
        total_signal = sum(total_signal),
        top_gene_id = gene_id[1L],
        top1_fraction = sum(head(shares, 1L)),
        top5_fraction = sum(head(shares, 5L)),
        top10_fraction = sum(head(shares, 10L)),
        genes_for_25pct = threshold_rank(cum_frac, 0.25),
        genes_for_50pct = threshold_rank(cum_frac, 0.50),
        genes_for_75pct = threshold_rank(cum_frac, 0.75),
        genes_for_90pct = threshold_rank(cum_frac, 0.90),
        effective_gene_number = if (sum(shares^2) > 0) 1 / sum(shares^2) else NA_real_
      )
    },
    by = .(sample, chrom_scope, feature, source)
  ]
}

build_top_gene_table <- function(gene_signal_summary, top_n = 50L) {
  if (!nrow(gene_signal_summary)) return(gene_signal_summary)
  gene_signal_summary[
    order(feature, source, gene_rank),
    head(.SD, top_n),
    by = .(sample, chrom_scope, feature, source)
  ]
}

build_dominance_curve <- function(gene_signal_summary) {
  if (!nrow(gene_signal_summary)) {
    return(data.table(
      sample = character(),
      chrom_scope = character(),
      feature = character(),
      source = character(),
      gene_id = character(),
      gene_rank = integer(),
      gene_rank_fraction = numeric(),
      cumulative_signal_fraction = numeric(),
      signal_fraction = numeric()
    ))
  }

  curve_dt <- copy(gene_signal_summary)
  curve_dt[, gene_rank_fraction := gene_rank / .N,
           by = .(sample, chrom_scope, feature, source)]
  curve_dt[, .(
    sample,
    chrom_scope,
    feature,
    source,
    gene_id,
    gene_rank,
    gene_rank_fraction,
    cumulative_signal_fraction,
    signal_fraction
  )]
}

set_profile_factors <- function(dt) {
  if (!nrow(dt)) return(dt)
  if ("source" %in% names(dt)) {
    dt[, source := factor(source, levels = source_levels)]
  }
  if ("feature" %in% names(dt)) {
    dt[, feature := factor(feature, levels = feature_levels)]
  }
  dt
}

compute_profile_correlation_table <- function(meta_norm) {
  rbindlist(lapply(setdiff(feature_levels, "Nucleosome"), function(feature_nm) {
    merged <- merge(
      meta_norm[feature == "Nucleosome", .(bin, source, nuc_signal = N)],
      meta_norm[feature == feature_nm, .(bin, source, feature_signal = N)],
      by = c("bin", "source"),
      all = FALSE
    )
    merged[, .(
      pearson = if (length(unique(nuc_signal)) < 2L || length(unique(feature_signal)) < 2L) NA_real_ else suppressWarnings(cor(nuc_signal, feature_signal, method = "pearson")),
      spearman = if (length(unique(nuc_signal)) < 2L || length(unique(feature_signal)) < 2L) NA_real_ else suppressWarnings(cor(nuc_signal, feature_signal, method = "spearman"))
    ), by = source][, feature := feature_nm]
  }), use.names = TRUE, fill = TRUE)
}

build_scope_plots <- function(meta_raw, meta_norm, profile_corr, gene_dominance,
                              dominance_curve, sample_name, scope_label,
                              window_bp = 2000L) {
  x_scale <- scale_x_continuous(
    breaks = seq(-window_bp, window_bp, 500),
    labels = function(x) paste0(x / 1000, " kb")
  )
  scope_subtitle <- sprintf("Sample: %s | %s", sample_name, scope_label)

  raw_plot <- ggplot(meta_raw, aes(x = bin, y = N, color = feature)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.08, se = FALSE, linewidth = 0.9) +
    facet_wrap(~source, nrow = 1) +
    scale_color_manual(values = feature_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x = "Distance to TSS",
      y = "Feature midpoints per 10 bp bin",
      color = "Feature",
      title = "Fiber-seq features around CAGE TSS: raw counts",
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  norm_plot <- ggplot(meta_norm, aes(x = bin, y = N, color = feature)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.08, se = FALSE, linewidth = 0.9) +
    facet_wrap(~source, nrow = 1) +
    scale_color_manual(values = feature_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x = "Distance to TSS",
      y = "Average feature midpoints per TSS (10 bp bins)",
      color = "Feature",
      title = "Fiber-seq features around CAGE TSS: per-TSS normalized profiles",
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  meta_z <- copy(meta_norm)
  meta_z[, z := if (sd(N) == 0) 0 else as.numeric(scale(N)), by = .(feature, source)]

  overlay_plot <- ggplot(meta_z, aes(x = bin, y = z, color = feature)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~source, nrow = 1) +
    scale_color_manual(values = feature_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x = "Distance to TSS",
      y = "Z-scored normalized signal",
      color = "Feature",
      title = "Relative positioning of nucleosome, m6A, CpG, and MSP signals",
      subtitle = "Z-scores are computed within each feature across all chromosomes"
    ) +
    theme_bw(base_size = 12)

  corr_plot <- ggplot(profile_corr, aes(x = source, y = feature, fill = spearman)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(spearman), "NA", sprintf("%.2f", spearman))), size = 3.5) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, limits = c(-1, 1), na.value = "grey90"
    ) +
    labs(
      x = "TSS source",
      y = "Feature",
      fill = "Spearman\nrho",
      title = "Correlation between nucleosome and feature profiles across TSS bins",
      subtitle = "Computed on the all-chromosome aggregated profile"
    ) +
    theme_bw(base_size = 12)

  dominance_plot <- ggplot(
    dominance_curve,
    aes(x = gene_rank_fraction, y = cumulative_signal_fraction)
  ) +
    geom_line(linewidth = 0.8, color = "#1f77b4") +
    facet_wrap(~feature, nrow = 2) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Fraction of genes ranked by descending signal",
      y = "Cumulative fraction of signal explained",
      title = "Gene dominance curves for CAGE-centered fiber feature signal",
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  list(
    raw_plot = raw_plot,
    norm_plot = norm_plot,
    overlay_plot = overlay_plot,
    corr_plot = corr_plot,
    dominance_plot = dominance_plot
  )
}

write_partial_outputs <- function(sample_name, chrom, feature_manifest, chrom_result, chunk_dir) {
  chrom_dir <- file.path(chunk_dir, chrom)
  dir.create(chrom_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(
    feature_manifest[chromosome == chrom],
    file.path(chrom_dir, sprintf("%s_%s_feature_manifest.tsv", sample_name, chrom)),
    sep = "\t"
  )
  fwrite(
    chrom_result$meta_raw,
    file.path(chrom_dir, sprintf("%s_%s_raw_profiles.tsv.gz", sample_name, chrom)),
    sep = "\t"
  )
  fwrite(
    chrom_result$gene_bin_signal,
    file.path(chrom_dir, sprintf("%s_%s_gene_bin_signal.tsv.gz", sample_name, chrom)),
    sep = "\t"
  )
  fwrite(
    chrom_result$gene_tss_summary,
    file.path(chrom_dir, sprintf("%s_%s_gene_tss_summary.tsv", sample_name, chrom)),
    sep = "\t"
  )
  saveRDS(
    chrom_result,
    file.path(chrom_dir, sprintf("%s_%s_partial_objects.rds", sample_name, chrom))
  )
}

process_chromosome <- function(sample_name, chrom, feature_manifest, tss_by_chr,
                               total_tss_count, window_bp = 2000L, bin_size = 10L) {
  feature_paths <- get_feature_paths_for_chrom(feature_manifest, chrom)
  feature_grs <- read_feature_granges(feature_paths)
  tss_chr <- tss_by_chr[[chrom]]
  if (is.null(tss_chr)) {
    tss_chr <- GRanges()
  }

  meta_raw_list <- list()
  meta_norm_list <- list()
  gene_bin_signal_list <- list()

  for (feature_nm in feature_levels) {
    hit_dt <- find_feature_hits_near_tss_with_gene(
      feature_grs[[feature_nm]],
      tss_chr,
      window_bp = window_bp,
      bin_size = bin_size
    )

    meta_pair <- summarize_meta_from_hits(
      hit_dt,
      total_tss_count = total_tss_count,
      sample_name = sample_name,
      scope_label = chrom,
      feature_name = feature_nm,
      window_bp = window_bp,
      bin_size = bin_size
    )
    meta_raw_list[[length(meta_raw_list) + 1L]] <- meta_pair$raw
    meta_norm_list[[length(meta_norm_list) + 1L]] <- meta_pair$norm

    gene_bin_signal_list[[length(gene_bin_signal_list) + 1L]] <- summarize_gene_bin_signal(
      hit_dt,
      sample_name = sample_name,
      scope_label = chrom,
      feature_name = feature_nm
    )
  }

  gene_tss_summary <- summarize_tss_by_gene(tss_chr, sample_name, chrom)

  list(
    meta_raw = set_profile_factors(rbindlist(meta_raw_list, use.names = TRUE, fill = TRUE)),
    meta_norm = set_profile_factors(rbindlist(meta_norm_list, use.names = TRUE, fill = TRUE)),
    gene_bin_signal = set_profile_factors(rbindlist(gene_bin_signal_list, use.names = TRUE, fill = TRUE)),
    gene_tss_summary = gene_tss_summary
  )
}

run_chunk_analysis <- function(sample_name, chrom, config) {
  feature_manifest <- discover_feature_manifest(sample_name)
  tss_gr <- load_cage_tss_annotations()
  tss_by_chr <- split_tss_by_chromosome(tss_gr)
  allowed_chroms <- select_analysis_chroms(feature_manifest, config)

  if (!chrom %in% allowed_chroms) {
    stop("Chromosome ", chrom, " is not in the selected chromosome set.")
  }

  if (!chrom %in% unique(feature_manifest$chromosome)) {
    message("Skipping ", chrom, " because no feature files were discovered.")
    return(invisible(NULL))
  }
  if (!chrom %in% names(tss_by_chr)) {
    message("No CAGE TSS annotations found for ", chrom, "; writing an empty chunk result.")
  }

  total_tss_count <- length(tss_gr)
  chrom_result <- process_chromosome(
    sample_name = sample_name,
    chrom = chrom,
    feature_manifest = feature_manifest,
    tss_by_chr = tss_by_chr,
    total_tss_count = total_tss_count,
    window_bp = config$window_bp,
    bin_size = config$bin_size
  )

  chunk_dir <- file.path(table_dir, sample_name, "chrom_chunks")
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
  write_partial_outputs(
    sample_name = sample_name,
    chrom = chrom,
    feature_manifest = feature_manifest,
    chrom_result = chrom_result,
    chunk_dir = chunk_dir
  )

  invisible(chrom_result)
}

process_and_write_chromosome <- function(sample_name, chrom, feature_manifest, tss_by_chr,
                                         total_tss_count, config, chunk_dir) {
  if (chunk_outputs_are_complete(sample_name, chrom, chunk_dir)) {
    message("Reusing cached chromosome: ", chrom)
    return(load_cached_chromosome_result(sample_name, chrom, chunk_dir))
  }

  chrom_result <- process_chromosome(
    sample_name = sample_name,
    chrom = chrom,
    feature_manifest = feature_manifest,
    tss_by_chr = tss_by_chr,
    total_tss_count = total_tss_count,
    window_bp = config$window_bp,
    bin_size = config$bin_size
  )

  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
  write_partial_outputs(
    sample_name = sample_name,
    chrom = chrom,
    feature_manifest = feature_manifest,
    chrom_result = chrom_result,
    chunk_dir = chunk_dir
  )

  chrom_result
}

run_full_analysis <- function(sample_name, config) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  feature_manifest <- discover_feature_manifest(sample_name)
  tss_gr <- load_cage_tss_annotations()
  tss_by_chr <- split_tss_by_chromosome(tss_gr)
  total_tss_count <- length(tss_gr)
  chroms <- select_analysis_chroms(feature_manifest, config)
  if (!length(chroms)) {
    stop("No chromosome-level feature files were discovered for sample ", sample_name)
  }

  chunk_dir <- file.path(table_dir, sample_name, "chrom_chunks")
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

  for (chrom in chroms) {
    message("Processing chromosome: ", chrom)
    chrom_result <- process_and_write_chromosome(
      sample_name = sample_name,
      chrom = chrom,
      feature_manifest = feature_manifest,
      tss_by_chr = tss_by_chr,
      total_tss_count = total_tss_count,
      config = config,
      chunk_dir = chunk_dir
    )
    rm(chrom_result)
    gc(verbose = FALSE)
  }

  combine_chunk_results(sample_name, config)
}

combine_chunk_results <- function(sample_name, config) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  message(
    "Combining partial chromosome results for sample ", sample_name,
    " | window_bp=", config$window_bp,
    " | bin_size=", config$bin_size,
    " | CAGE TSS only"
  )

  feature_manifest <- discover_feature_manifest(sample_name)

  chunk_dir <- file.path(table_dir, sample_name, "chrom_chunks")
  final_dir <- file.path(table_dir, sample_name, "all_chr")
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

  chroms <- select_analysis_chroms(feature_manifest, config)
  if (!length(chroms)) {
    stop("No chromosome-level feature files were discovered for sample ", config$sample)
  }

  partial_paths <- file.path(
    chunk_dir,
    chroms,
    sprintf("%s_%s_partial_objects.rds", sample_name, chroms)
  )
  missing_partial <- partial_paths[!file.exists(partial_paths)]
  if (length(missing_partial)) {
    stop(
      "Missing partial result files. Expected chromosome chunks at:\n",
      paste(missing_partial, collapse = "\n")
    )
  }

  chrom_manifest <- data.table(
    sample = sample_name,
    chrom_scope = chroms,
    partial_dir = file.path(chunk_dir, chroms)
  )

  tss_gr <- load_cage_tss_annotations()
  total_tss_count <- length(tss_gr)

  aggregated <- aggregate_partial_chrom_results_from_paths(partial_paths, total_tss_count)
  meta_raw <- aggregated$meta_raw
  meta_norm <- aggregated$meta_norm
  gene_bin_signal <- aggregated$gene_bin_signal

  gene_tss_summary <- summarize_tss_by_gene(tss_gr, sample_name, "all_chr")
  gene_signal_summary <- summarize_gene_signal(gene_bin_signal, gene_tss_summary)
  gene_dominance <- summarize_gene_dominance(gene_signal_summary)
  top_gene_contributors <- build_top_gene_table(gene_signal_summary, top_n = 50L)
  dominance_curve <- build_dominance_curve(gene_signal_summary)
  profile_corr <- compute_profile_correlation_table(meta_norm)

  meta_raw <- set_profile_factors(meta_raw)
  meta_norm <- set_profile_factors(meta_norm)
  gene_bin_signal <- set_profile_factors(gene_bin_signal)
  gene_signal_summary <- set_profile_factors(gene_signal_summary)
  gene_dominance <- set_profile_factors(gene_dominance)
  top_gene_contributors <- set_profile_factors(top_gene_contributors)
  dominance_curve <- set_profile_factors(dominance_curve)
  profile_corr <- set_profile_factors(profile_corr)

  plots <- build_scope_plots(
    meta_raw = meta_raw,
    meta_norm = meta_norm,
    profile_corr = profile_corr,
    gene_dominance = gene_dominance,
    dominance_curve = dominance_curve,
    sample_name = config$sample,
    scope_label = "all_chr",
    window_bp = config$window_bp
  )

  print(plots$raw_plot)
  print(plots$norm_plot)
  print(plots$overlay_plot)
  print(plots$corr_plot)
  print(plots$dominance_plot)

  ggsave(
    file.path(plot_dir, sprintf("%s_all_chr_FIRE_features_raw.pdf", config$sample)),
    plot = plots$raw_plot,
    width = 9,
    height = 6
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_all_chr_FIRE_features_norm.pdf", config$sample)),
    plot = plots$norm_plot,
    width = 9,
    height = 6
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_all_chr_FIRE_features_overlay.pdf", config$sample)),
    plot = plots$overlay_plot,
    width = 8,
    height = 5
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_all_chr_FIRE_features_correlation_heatmap.pdf", config$sample)),
    plot = plots$corr_plot,
    width = 6.5,
    height = 4
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_all_chr_FIRE_features_gene_dominance_curves.pdf", config$sample)),
    plot = plots$dominance_plot,
    width = 8.5,
    height = 7
  )

  analysis_config <- data.table(
    sample = sample_name,
    chrom_scope = "all_chr",
    window_bp = config$window_bp,
    bin_size = config$bin_size,
    tss_source = "CAGE",
    feature_levels = paste(feature_levels, collapse = ","),
    total_tss_count = total_tss_count,
    n_chromosomes = length(chroms)
  )

  fwrite(
    analysis_config,
    file.path(final_dir, sprintf("%s_all_chr_analysis_config.tsv", sample_name)),
    sep = "\t"
  )
  fwrite(
    feature_manifest,
    file.path(final_dir, sprintf("%s_all_chr_feature_manifest.tsv", sample_name)),
    sep = "\t"
  )
  fwrite(
    meta_raw,
    file.path(final_dir, sprintf("%s_all_chr_fiber_features_raw_profiles.tsv.gz", sample_name)),
    sep = "\t"
  )
  fwrite(
    meta_norm,
    file.path(final_dir, sprintf("%s_all_chr_fiber_features_norm_profiles.tsv.gz", sample_name)),
    sep = "\t"
  )
  fwrite(
    profile_corr,
    file.path(final_dir, sprintf("%s_all_chr_fiber_features_vs_nucleosome_correlations.tsv", sample_name)),
    sep = "\t"
  )
  fwrite(
    gene_tss_summary,
    file.path(final_dir, sprintf("%s_all_chr_gene_tss_summary.tsv", sample_name)),
    sep = "\t"
  )
  fwrite(
    gene_bin_signal,
    file.path(final_dir, sprintf("%s_all_chr_gene_bin_signal.tsv.gz", sample_name)),
    sep = "\t"
  )
  fwrite(
    gene_signal_summary,
    file.path(final_dir, sprintf("%s_all_chr_gene_signal_summary.tsv.gz", sample_name)),
    sep = "\t"
  )
  fwrite(
    gene_dominance,
    file.path(final_dir, sprintf("%s_all_chr_gene_dominance.tsv", sample_name)),
    sep = "\t"
  )
  fwrite(
    top_gene_contributors,
    file.path(final_dir, sprintf("%s_all_chr_top_gene_contributors.tsv", sample_name)),
    sep = "\t"
  )
  fwrite(
    dominance_curve,
    file.path(final_dir, sprintf("%s_all_chr_gene_dominance_curves.tsv.gz", sample_name)),
    sep = "\t"
  )
  fwrite(
    chrom_manifest,
    file.path(final_dir, sprintf("%s_all_chr_partial_output_manifest.tsv", sample_name)),
    sep = "\t"
  )

  saveRDS(
    list(
      analysis_config = analysis_config,
      feature_manifest = feature_manifest,
      meta_raw = meta_raw,
      meta_norm = meta_norm,
      profile_corr = profile_corr,
      gene_tss_summary = gene_tss_summary,
      gene_bin_signal = gene_bin_signal,
      gene_signal_summary = gene_signal_summary,
      gene_dominance = gene_dominance,
      top_gene_contributors = top_gene_contributors,
      dominance_curve = dominance_curve
    ),
    file.path(final_dir, sprintf("%s_all_chr_FIRE_features_analysis_objects.rds", sample_name))
  )

  invisible(list(
    sample = sample_name,
    feature_manifest = feature_manifest,
    meta_raw = meta_raw,
    meta_norm = meta_norm,
    profile_corr = profile_corr,
    gene_tss_summary = gene_tss_summary,
    gene_bin_signal = gene_bin_signal,
    gene_signal_summary = gene_signal_summary,
    gene_dominance = gene_dominance
  ))
}

run_all_chr_aggregated_FIRE_features <- function() {
  config <- resolve_run_config()
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  if (config$mode == "chunk") {
    if (!nzchar(config$chrom)) {
      stop("mode=chunk requires chrom=<chromosome>.")
    }
    message(
      "Running chunk for sample ", config$sample,
      " | chrom=", config$chrom,
      " | window_bp=", config$window_bp,
      " | bin_size=", config$bin_size
    )
    return(run_chunk_analysis(config$sample, config$chrom, config))
  }

  if (config$mode == "full") {
    message(
      "Running full sequential chromosome processing for sample ", config$sample,
      " | window_bp=", config$window_bp,
      " | bin_size=", config$bin_size
    )
    return(run_full_analysis(config$sample, config))
  }

  if (config$mode == "combine") {
    message(
      "Combining partial chromosome results for sample ", config$sample,
      " | window_bp=", config$window_bp,
      " | bin_size=", config$bin_size
    )
    return(combine_chunk_results(config$sample, config))
  }

  stop("Unsupported mode: ", config$mode, ". Use mode=chunk, mode=full, or mode=combine.")
}

run_all_chr_aggregated_FIRE_features()

# Example:
# Rscript /project/spott/cshan/fiber-seq/code/all_chr_aggregated_FIRE_features.r sample=AL10_bc2178_19130 mode=chunk chrom=chrM
# Rscript /project/spott/cshan/fiber-seq/code/all_chr_aggregated_FIRE_features.r sample=AL10_bc2178_19130 mode=combine
