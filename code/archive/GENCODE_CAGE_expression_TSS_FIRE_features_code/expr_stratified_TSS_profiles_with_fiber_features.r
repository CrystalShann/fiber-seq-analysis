source("/project/spott/cshan/fiber-seq/plot_nuc_distribution_functions.R")

gencode_tss_path <- "/project/spott/cshan/annotations/TSS.gencode.v49.bed"
cage_tss_path <- "/project/spott/cshan/annotations/fantom5/fantom5.hg38.LCL.consensus.CAGE_peaks.withGene.bed.gz"
sample_root <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results"
expr_path <- "/project/spott/cshan/annotations/ENCFF873VWU.gene_quan.tsv"
default_sample <- "AL10_bc2178_19130"
default_chrom <- "chr1"
window_bp <- 2000L
bin_size <- 10L

plot_dir <- "/project/spott/cshan/fiber-seq/results/plots/nuc"
table_dir <- "/project/spott/cshan/fiber-seq/results/nuc/expression_tables_TSS_profiles"
ft_extract_feature_dir <- "/project/spott/cshan/fiber-seq/results/ft_extract_features"

source_colors <- c("GENCODE" = "#1f77b4", "CAGE" = "#be1c54ff")
expr_colors <- c("high" = "#FF7F0E", "low" = "#2CA02C")
feature_colors <- c(
  "Nucleosome" = "#4D4D4D",
  "m6A" = "#7B3294",
  "CpG" = "#8C510A",
  "MSP" = "#1B9E77"
)
feature_levels <- c("Nucleosome", "m6A", "CpG", "MSP")
expr_group_levels <- c("high", "low")
tss_source_levels <- c("GENCODE", "CAGE")

annotate_gene_id_base <- function(gr) {
  mcols(gr)$gene_id_base <- sub("\\..*$", "", mcols(gr)$gene_id)
  gr
}

expr <- fread(expr_path)
expr <- expr[grepl("^ENSG", gene_id)]
expr[, gene_id_base := sub("\\..*$", "", gene_id)]
qs <- quantile(expr$TPM, probs = c(0.25, 0.75), na.rm = TRUE)
expr[, expr_group := fcase(
  TPM >= qs[2], "high",
  TPM <= qs[1], "low",
  default = NA_character_
)]
expr <- expr[!is.na(expr_group), .(gene_id, gene_id_base, TPM, expr_group)]

tss_sets <- list(
  GENCODE = annotate_gene_id_base(load_tss_annotations(gencode_tss_path)),
  CAGE = annotate_gene_id_base(load_tss_annotations_bed9_with_gene_id(cage_tss_path))
)

align_tss <- function(gr, source) {
  idx_full <- match(mcols(gr)$gene_id, expr$gene_id)
  idx_base <- match(mcols(gr)$gene_id_base, expr$gene_id_base)
  idx <- ifelse(!is.na(idx_full), idx_full, idx_base)
  valid <- !is.na(idx)
  if (!any(valid)) return(NULL)
  gr2 <- gr[valid]
  mcols(gr2)$expr_group <- expr$expr_group[idx[valid]]
  mcols(gr2)$TPM <- expr$TPM[idx[valid]]
  list(source = source, gr = gr2)
}

tss_by_group <- Filter(
  Negate(is.null),
  lapply(names(tss_sets), function(source) align_tss(tss_sets[[source]], source))
)

matched_sources <- vapply(tss_by_group, function(x) x$source, character(1))
if (!all(c("GENCODE", "CAGE") %in% matched_sources)) {
  stop("One or more TSS sources did not match the expression table; check CAGE gene_id annotation.")
}

safe_cor <- function(x, y, method = "pearson") {
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) return(NA_real_)
  suppressWarnings(cor(x, y, method = method))
}

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

parse_bool_flag <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) return(default)
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

is_all_chrom_scope <- function(chrom_value) {
  tolower(chrom_value) %in% c("all", "all_chr", "all_chrom", "genome", "*")
}

normalize_scope_label <- function(chrom_value) {
  if (is_all_chrom_scope(chrom_value)) "all_chr" else chrom_value
}

scope_subtitle_label <- function(scope_label) {
  if (identical(scope_label, "all_chr")) "All chromosomes" else paste("Chromosome:", scope_label)
}

resolve_run_config <- function() {
  args <- parse_named_args(commandArgs(trailingOnly = TRUE))
  sample_name <- args[["sample"]] %||% Sys.getenv("FIBER_SAMPLE", default_sample)
  chrom_value <- args[["chrom"]] %||% Sys.getenv("FIBER_CHROM", default_chrom)
  window_value <- as.integer(args[["window_bp"]] %||% Sys.getenv("FIBER_WINDOW_BP", as.character(window_bp)))
  bin_value <- as.integer(args[["bin_size"]] %||% Sys.getenv("FIBER_BIN_SIZE", as.character(bin_size)))
  run_all_chr <- parse_bool_flag(
    args[["run_all_chr"]] %||% Sys.getenv("FIBER_RUN_ALL_CHR", ""),
    default = TRUE
  )

  if (is_all_chrom_scope(chrom_value)) {
    scope_labels <- "all_chr"
  } else if (run_all_chr) {
    scope_labels <- unique(c(chrom_value, "all_chr"))
  } else {
    scope_labels <- chrom_value
  }

  list(
    sample = sample_name,
    scopes = scope_labels,
    window_bp = window_value,
    bin_size = bin_value,
    run_all_chr = run_all_chr
  )
}

feature_specs <- list(
  Nucleosome = list(subdir = "nuc_by_chr", stem = "nuc"),
  m6A = list(subdir = "m6a_by_chr", stem = "m6a"),
  CpG = list(subdir = "cpg_by_chr", stem = "cpg"),
  MSP = list(subdir = "msp_by_chr", stem = "msp")
)

extract_chrom_from_feature_path <- function(path) {
  sub("^.*\\.(chr[^.]+)\\.bed\\.gz$", "\\1", basename(path))
}

collect_feature_files_for_scope <- function(sample_name, scope_label) {
  sample_extract_dir <- file.path(sample_root, sample_name, "extracted_results")

  feature_files <- lapply(names(feature_specs), function(feature_nm) {
    spec <- feature_specs[[feature_nm]]
    feature_dir <- file.path(sample_extract_dir, spec$subdir)
    if (!dir.exists(feature_dir)) {
      stop("Missing extracted feature directory: ", feature_dir)
    }

    if (identical(scope_label, "all_chr")) {
      pattern <- sprintf("^%s\\.ft_extracted_%s\\.chr[^.]+\\.bed\\.gz$", sample_name, spec$stem)
      paths <- sort(list.files(feature_dir, pattern = pattern, full.names = TRUE))
      if (!length(paths)) {
        stop("No feature BED files matched ", feature_nm, " for ", sample_name, " in ", feature_dir)
      }
      paths
    } else {
      path <- file.path(
        feature_dir,
        sprintf("%s.ft_extracted_%s.%s.bed.gz", sample_name, spec$stem, scope_label)
      )
      if (!file.exists(path)) {
        stop("Missing feature BED12 file for ", feature_nm, ": ", path)
      }
      path
    }
  })

  names(feature_files) <- names(feature_specs)
  feature_files
}

build_feature_file_manifest <- function(sample_name, scope_label, feature_files) {
  rbindlist(lapply(names(feature_files), function(feature_nm) {
    paths <- feature_files[[feature_nm]]
    data.table(
      sample = sample_name,
      chrom_scope = scope_label,
      feature = feature_nm,
      chromosome = vapply(paths, extract_chrom_from_feature_path, character(1)),
      feature_file = basename(paths),
      feature_path = paths
    )
  }), use.names = TRUE, fill = TRUE)
}

load_feature_granges_for_scope <- function(feature_files) {
  setNames(
    lapply(names(feature_files), function(feature_nm) {
      load_fiber_feature_positions(feature_files[[feature_nm]], feature_name = feature_nm)
    }),
    names(feature_files)
  )
}

make_tss_group_counts <- function(tss_gr) {
  counts <- data.table(
    expr_group = expr_group_levels,
    n_tss = 0L
  )
  observed <- data.table(expr_group = as.character(mcols(tss_gr)$expr_group))[
    , .(n_tss = .N), by = expr_group
  ]
  counts <- merge(counts, observed, by = "expr_group", all.x = TRUE, suffixes = c("", "_obs"))
  counts[, n_tss := fifelse(is.na(n_tss_obs), n_tss, n_tss_obs)]
  counts[, n_tss_obs := NULL]
  counts[]
}

summarize_tss_by_gene <- function(tss_gr, source_name, sample_name, scope_label) {
  dt <- data.table(
    sample = sample_name,
    chrom_scope = scope_label,
    source = source_name,
    gene_id = as.character(mcols(tss_gr)$gene_id),
    gene_id_base = as.character(mcols(tss_gr)$gene_id_base),
    expr_group = as.character(mcols(tss_gr)$expr_group),
    seqnames = as.character(seqnames(tss_gr))
  )

  dt[, .(
    n_tss_gene = .N,
    n_gene_chromosomes = uniqueN(seqnames)
  ), by = .(sample, chrom_scope, source, gene_id, gene_id_base, expr_group)]
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
    bin = integer(),
    expr_group = character(),
    TPM = numeric()
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
    bin = as.integer(round(rel_pos / bin_size) * bin_size),
    expr_group = as.character(mcols(tss_gr)$expr_group[subject_idx]),
    TPM = as.numeric(mcols(tss_gr)$TPM[subject_idx])
  )[bin >= -window_bp & bin <= window_bp]
}

summarize_meta_from_hits <- function(hit_dt, tss_group_counts,
                                     sample_name, scope_label,
                                     source_name, feature_name,
                                     window_bp = 2000L, bin_size = 10L) {
  template <- as.data.table(expand.grid(
    expr_group = expr_group_levels,
    bin = seq(-window_bp, window_bp, by = bin_size),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ))

  counts <- if (nrow(hit_dt)) {
    hit_dt[, .(N = .N), by = .(expr_group, bin)]
  } else {
    data.table(expr_group = character(), bin = integer(), N = integer())
  }

  meta_raw <- merge(template, counts, by = c("expr_group", "bin"), all.x = TRUE, sort = TRUE)
  meta_raw[is.na(N), N := 0L]
  meta_raw[, `:=`(
    sample = sample_name,
    chrom_scope = scope_label,
    source = source_name,
    feature = feature_name
  )]

  meta_norm <- merge(meta_raw, tss_group_counts, by = "expr_group", all.x = TRUE, sort = FALSE)
  meta_norm[is.na(n_tss), n_tss := 0L]
  meta_norm[, N := fifelse(n_tss > 0, N / n_tss, 0)]
  meta_norm[, n_tss := NULL]

  list(raw = meta_raw[], norm = meta_norm[])
}

summarize_gene_bin_signal <- function(hit_dt, sample_name, scope_label,
                                      source_name, feature_name) {
  if (!nrow(hit_dt)) {
    return(data.table(
      sample = character(),
      chrom_scope = character(),
      source = character(),
      feature = character(),
      expr_group = character(),
      gene_id = character(),
      gene_id_base = character(),
      TPM = numeric(),
      bin = integer(),
      signal = integer(),
      n_tss_with_signal = integer()
    ))
  }

  hit_dt[, .(
    signal = .N,
    n_tss_with_signal = uniqueN(tss_index)
  ), by = .(expr_group, gene_id, gene_id_base, TPM, bin)][
    , `:=`(
      sample = sample_name,
      chrom_scope = scope_label,
      source = source_name,
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
      expr_group = character(),
      gene_id = character(),
      gene_id_base = character(),
      TPM = numeric(),
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
    by = .(sample, chrom_scope, source, feature, expr_group, gene_id, gene_id_base, TPM)
  ]

  gene_summary <- merge(
    gene_summary,
    gene_tss_summary,
    by = c("sample", "chrom_scope", "source", "gene_id", "gene_id_base", "expr_group"),
    all.x = TRUE,
    sort = FALSE
  )

  gene_summary[is.na(n_tss_gene), n_tss_gene := 0L]
  gene_summary[is.na(n_gene_chromosomes), n_gene_chromosomes := 0L]
  gene_summary[, signal_per_tss_gene := fifelse(n_tss_gene > 0, total_signal / n_tss_gene, 0)]
  gene_summary[, central_signal_fraction := fifelse(total_signal > 0, central_signal_250bp / total_signal, 0)]
  gene_summary <- gene_summary[order(feature, source, expr_group, -total_signal, gene_id)]
  gene_summary[, group_signal := sum(total_signal), by = .(sample, chrom_scope, feature, source, expr_group)]
  gene_summary[, signal_fraction := fifelse(group_signal > 0, total_signal / group_signal, 0)]
  gene_summary[, group_signal := NULL]
  gene_summary[, gene_rank := seq_len(.N), by = .(sample, chrom_scope, feature, source, expr_group)]
  gene_summary[, cumulative_signal_fraction := cumsum(signal_fraction),
               by = .(sample, chrom_scope, feature, source, expr_group)]
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
      expr_group = character(),
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
    by = .(sample, chrom_scope, feature, source, expr_group)
  ]
}

build_top_gene_table <- function(gene_signal_summary, top_n = 50L) {
  if (!nrow(gene_signal_summary)) return(gene_signal_summary)
  gene_signal_summary[
    order(feature, source, expr_group, gene_rank),
    head(.SD, top_n),
    by = .(sample, chrom_scope, feature, source, expr_group)
  ]
}

build_dominance_curve <- function(gene_signal_summary) {
  if (!nrow(gene_signal_summary)) {
    return(data.table(
      sample = character(),
      chrom_scope = character(),
      feature = character(),
      source = character(),
      expr_group = character(),
      gene_id = character(),
      gene_rank = integer(),
      gene_rank_fraction = numeric(),
      cumulative_signal_fraction = numeric(),
      signal_fraction = numeric()
    ))
  }

  curve_dt <- copy(gene_signal_summary)
  curve_dt[, gene_rank_fraction := fifelse(.N > 0, gene_rank / .N, 0),
           by = .(sample, chrom_scope, feature, source, expr_group)]
  curve_dt[, .(
    sample,
    chrom_scope,
    feature,
    source,
    expr_group,
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
    dt[, source := factor(source, levels = tss_source_levels)]
  }
  if ("expr_group" %in% names(dt)) {
    dt[, expr_group := factor(expr_group, levels = expr_group_levels)]
  }
  if ("feature" %in% names(dt)) {
    dt[, feature := factor(feature, levels = feature_levels)]
  }
  dt
}

compute_profile_correlation_table <- function(meta_norm) {
  profile_corr <- rbindlist(lapply(setdiff(feature_levels, "Nucleosome"), function(feature_nm) {
    merged <- merge(
      meta_norm[feature == "Nucleosome", .(bin, source, expr_group, nuc_signal = N)],
      meta_norm[feature == feature_nm, .(bin, source, expr_group, feature_signal = N)],
      by = c("bin", "source", "expr_group"),
      all = FALSE
    )
    merged[, .(
      pearson = safe_cor(nuc_signal, feature_signal, method = "pearson"),
      spearman = safe_cor(nuc_signal, feature_signal, method = "spearman")
    ), by = .(source, expr_group)][, feature := feature_nm]
  }), use.names = TRUE, fill = TRUE)

  profile_corr[, feature := factor(feature, levels = c("m6A", "CpG", "MSP"))]
  profile_corr
}

compute_expr_enrichment_table <- function(meta_norm) {
  feature_window_signal <- meta_norm[, .(signal_per_tss = sum(N)),
                                     by = .(feature, source, expr_group)]
  expr_enrichment <- dcast(
    feature_window_signal,
    feature + source ~ expr_group,
    value.var = "signal_per_tss",
    fill = 0
  )
  expr_enrichment[, log2_high_vs_low := log2((high + 1e-6) / (low + 1e-6))]
  expr_enrichment[, feature := factor(feature, levels = feature_levels)]
  expr_enrichment
}

build_scope_plots <- function(meta_raw, meta_norm, profile_corr, expr_enrichment,
                              dominance_curve, sample_name, scope_label,
                              window_bp = 2000L) {
  x_scale <- scale_x_continuous(
    breaks = seq(-window_bp, window_bp, 500),
    labels = function(x) paste0(x / 1000, " kb")
  )
  scope_subtitle <- sprintf("Sample: %s | %s", sample_name, scope_subtitle_label(scope_label))

  raw_source_plot <- ggplot(meta_raw, aes(x = bin, y = N, color = source)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.06, se = FALSE, linewidth = 0.9) +
    facet_grid(feature ~ expr_group, scales = "free_y") +
    scale_color_manual(values = source_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x = "Distance to TSS",
      y = "Feature midpoints per 10 bp bin",
      color = "TSS source",
      title = "Fiber-seq features around TSS: raw counts by expression group",
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  norm_source_plot <- ggplot(meta_norm, aes(x = bin, y = N, color = source)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.06, se = FALSE, linewidth = 0.9) +
    facet_grid(feature ~ expr_group, scales = "free_y") +
    scale_color_manual(values = source_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x = "Distance to TSS",
      y = "Average feature midpoints per TSS (10 bp bins)",
      color = "TSS source",
      title = "Fiber-seq features around TSS: per-TSS normalized profiles",
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  raw_expr_plot <- ggplot(meta_raw, aes(x = bin, y = N, color = expr_group)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.06, se = FALSE, linewidth = 0.9) +
    facet_grid(feature ~ source, scales = "free_y") +
    scale_color_manual(values = expr_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x = "Distance to TSS",
      y = "Feature midpoints per 10 bp bin",
      color = "Expression bin",
      title = "Fiber-seq features around TSS: raw counts by TSS source",
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  norm_expr_plot <- ggplot(meta_norm, aes(x = bin, y = N, color = expr_group)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.06, se = FALSE, linewidth = 0.9) +
    facet_grid(feature ~ source, scales = "free_y") +
    scale_color_manual(values = expr_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x = "Distance to TSS",
      y = "Average feature midpoints per TSS (10 bp bins)",
      color = "Expression bin",
      title = "Fiber-seq features around TSS: normalized by TSS source",
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  meta_z <- copy(meta_norm)
  meta_z[, z := if (sd(N) == 0) 0 else as.numeric(scale(N)),
         by = .(feature, source, expr_group)]

  shape_overlay_plot <- ggplot(meta_z, aes(x = bin, y = z, color = feature)) +
    geom_line(linewidth = 0.8) +
    facet_grid(expr_group ~ source) +
    scale_color_manual(values = feature_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x = "Distance to TSS",
      y = "Z-scored normalized signal",
      color = "Feature",
      title = "Relative positioning of nucleosome, m6A, CpG, and MSP signals",
      subtitle = "Z-scores are computed within each feature, TSS source, and expression group"
    ) +
    theme_bw(base_size = 12)

  corr_plot <- ggplot(profile_corr, aes(x = source, y = feature, fill = spearman)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", spearman)), size = 3.5) +
    facet_wrap(~expr_group, nrow = 1) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, limits = c(-1, 1), na.value = "grey90"
    ) +
    labs(
      x = "TSS source",
      y = "Feature",
      fill = "Spearman\nrho",
      title = "Correlation between nucleosome and feature profiles across TSS bins",
      subtitle = "Positive values indicate co-enrichment in the same TSS-relative bins"
    ) +
    theme_bw(base_size = 12)

  expr_enrichment_plot <- ggplot(expr_enrichment, aes(x = source, y = feature, fill = log2_high_vs_low)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", log2_high_vs_low)), size = 3.5) +
    scale_fill_gradient2(low = "#2CA02C", mid = "white", high = "#FF7F0E", midpoint = 0) +
    labs(
      x = "TSS source",
      y = "Feature",
      fill = "log2(high/low)",
      title = "Expression association of fiber-seq features within +/-2 kb of TSS",
      subtitle = "Computed from per-TSS normalized signal summed across bins"
    ) +
    theme_bw(base_size = 12)

  dominance_plot <- ggplot(
    dominance_curve,
    aes(x = gene_rank_fraction, y = cumulative_signal_fraction, color = expr_group)
  ) +
    geom_line(linewidth = 0.8) +
    facet_grid(feature ~ source) +
    scale_color_manual(values = expr_colors, drop = FALSE) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Fraction of genes ranked by descending signal",
      y = "Cumulative fraction of signal explained",
      color = "Expression bin",
      title = "Gene dominance curves for TSS-centered fiber feature signal",
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  list(
    raw_source_plot = raw_source_plot,
    norm_source_plot = norm_source_plot,
    raw_expr_plot = raw_expr_plot,
    norm_expr_plot = norm_expr_plot,
    shape_overlay_plot = shape_overlay_plot,
    corr_plot = corr_plot,
    expr_enrichment_plot = expr_enrichment_plot,
    dominance_plot = dominance_plot
  )
}

write_scope_outputs <- function(sample_name, scope_label,
                                feature_manifest, analysis_config,
                                meta_raw, meta_norm, profile_corr,
                                expr_enrichment, gene_bin_signal,
                                gene_signal_summary, gene_dominance,
                                top_gene_contributors, dominance_curve,
                                plots) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  scope_out_dir <- file.path(ft_extract_feature_dir, sample_name, scope_label)
  dir.create(scope_out_dir, recursive = TRUE, showWarnings = FALSE)

  print(plots$raw_source_plot)
  print(plots$norm_source_plot)
  print(plots$raw_expr_plot)
  print(plots$norm_expr_plot)
  print(plots$shape_overlay_plot)
  print(plots$corr_plot)
  print(plots$expr_enrichment_plot)
  print(plots$dominance_plot)

  ggsave(
    file.path(plot_dir, sprintf("%s_%s_fiber_features_raw_by_source_expr.pdf", sample_name, scope_label)),
    plot = plots$raw_source_plot, width = 10, height = 10
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_%s_fiber_features_norm_by_source_expr.pdf", sample_name, scope_label)),
    plot = plots$norm_source_plot, width = 10, height = 10
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_%s_fiber_features_raw_by_expr_source.pdf", sample_name, scope_label)),
    plot = plots$raw_expr_plot, width = 10, height = 10
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_%s_fiber_features_norm_by_expr_source.pdf", sample_name, scope_label)),
    plot = plots$norm_expr_plot, width = 10, height = 10
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_%s_fiber_features_zscore_overlay.pdf", sample_name, scope_label)),
    plot = plots$shape_overlay_plot, width = 8, height = 5
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_%s_fiber_features_vs_nucleosome_correlation_heatmap.pdf", sample_name, scope_label)),
    plot = plots$corr_plot, width = 7, height = 4
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_%s_fiber_features_expression_heatmap.pdf", sample_name, scope_label)),
    plot = plots$expr_enrichment_plot, width = 7, height = 4.5
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_%s_fiber_features_gene_dominance_curves.pdf", sample_name, scope_label)),
    plot = plots$dominance_plot, width = 8.5, height = 7
  )

  fwrite(
    meta_raw,
    file.path(table_dir, sprintf("%s_%s_fiber_features_raw_profiles.tsv", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    meta_norm,
    file.path(table_dir, sprintf("%s_%s_fiber_features_norm_profiles.tsv", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    profile_corr,
    file.path(table_dir, sprintf("%s_%s_fiber_features_vs_nucleosome_correlations.tsv", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    expr_enrichment,
    file.path(table_dir, sprintf("%s_%s_fiber_features_high_vs_low_expression.tsv", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    gene_dominance,
    file.path(table_dir, sprintf("%s_%s_fiber_features_gene_dominance.tsv", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    top_gene_contributors,
    file.path(table_dir, sprintf("%s_%s_fiber_features_top_gene_contributors.tsv", sample_name, scope_label)),
    sep = "\t"
  )

  fwrite(
    analysis_config,
    file.path(scope_out_dir, sprintf("%s_%s_analysis_config.tsv", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    feature_manifest,
    file.path(scope_out_dir, sprintf("%s_%s_feature_file_manifest.tsv", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    meta_raw,
    file.path(scope_out_dir, sprintf("%s_%s_fiber_features_raw_profiles.tsv.gz", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    meta_norm,
    file.path(scope_out_dir, sprintf("%s_%s_fiber_features_norm_profiles.tsv.gz", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    profile_corr,
    file.path(scope_out_dir, sprintf("%s_%s_fiber_features_vs_nucleosome_correlations.tsv", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    expr_enrichment,
    file.path(scope_out_dir, sprintf("%s_%s_fiber_features_high_vs_low_expression.tsv", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    gene_bin_signal,
    file.path(scope_out_dir, sprintf("%s_%s_fiber_features_gene_bin_signal.tsv.gz", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    gene_signal_summary,
    file.path(scope_out_dir, sprintf("%s_%s_fiber_features_gene_signal_summary.tsv.gz", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    gene_dominance,
    file.path(scope_out_dir, sprintf("%s_%s_fiber_features_gene_dominance.tsv", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    top_gene_contributors,
    file.path(scope_out_dir, sprintf("%s_%s_fiber_features_top_gene_contributors.tsv", sample_name, scope_label)),
    sep = "\t"
  )
  fwrite(
    dominance_curve,
    file.path(scope_out_dir, sprintf("%s_%s_fiber_features_gene_dominance_curves.tsv.gz", sample_name, scope_label)),
    sep = "\t"
  )

  saveRDS(
    list(
      analysis_config = analysis_config,
      feature_manifest = feature_manifest,
      meta_raw = meta_raw,
      meta_norm = meta_norm,
      profile_corr = profile_corr,
      expr_enrichment = expr_enrichment,
      gene_bin_signal = gene_bin_signal,
      gene_signal_summary = gene_signal_summary,
      gene_dominance = gene_dominance,
      top_gene_contributors = top_gene_contributors,
      dominance_curve = dominance_curve
    ),
    file.path(scope_out_dir, sprintf("%s_%s_fiber_features_analysis_objects.rds", sample_name, scope_label))
  )
}

run_scope_analysis <- function(sample_name, scope_label,
                               window_bp = 2000L, bin_size = 10L) {
  feature_files <- collect_feature_files_for_scope(sample_name, scope_label)
  feature_manifest <- build_feature_file_manifest(sample_name, scope_label, feature_files)

  meta_raw_list <- list()
  meta_norm_list <- list()
  gene_bin_signal_list <- list()
  gene_signal_summary_list <- list()

  for (entry in tss_by_group) {
    tss_group_counts <- make_tss_group_counts(entry$gr)
    gene_tss_summary <- summarize_tss_by_gene(entry$gr, entry$source, sample_name, scope_label)

    for (feature_nm in names(feature_files)) {
      feature_gr <- load_fiber_feature_positions(feature_files[[feature_nm]], feature_name = feature_nm)
      hit_dt <- find_feature_hits_near_tss_with_gene(
        feature_gr,
        entry$gr,
        window_bp = window_bp,
        bin_size = bin_size
      )

      meta_pair <- summarize_meta_from_hits(
        hit_dt,
        tss_group_counts,
        sample_name = sample_name,
        scope_label = scope_label,
        source_name = entry$source,
        feature_name = feature_nm,
        window_bp = window_bp,
        bin_size = bin_size
      )
      meta_raw_list[[length(meta_raw_list) + 1L]] <- meta_pair$raw
      meta_norm_list[[length(meta_norm_list) + 1L]] <- meta_pair$norm

      gene_bin_signal <- summarize_gene_bin_signal(
        hit_dt,
        sample_name = sample_name,
        scope_label = scope_label,
        source_name = entry$source,
        feature_name = feature_nm
      )
      gene_signal_summary <- summarize_gene_signal(gene_bin_signal, gene_tss_summary)

      gene_bin_signal_list[[length(gene_bin_signal_list) + 1L]] <- gene_bin_signal
      gene_signal_summary_list[[length(gene_signal_summary_list) + 1L]] <- gene_signal_summary
      rm(feature_gr)
      gc()
    }
  }

  meta_raw <- rbindlist(meta_raw_list, use.names = TRUE, fill = TRUE)
  meta_norm <- rbindlist(meta_norm_list, use.names = TRUE, fill = TRUE)
  gene_bin_signal <- rbindlist(gene_bin_signal_list, use.names = TRUE, fill = TRUE)
  gene_signal_summary <- rbindlist(gene_signal_summary_list, use.names = TRUE, fill = TRUE)

  meta_raw <- set_profile_factors(meta_raw)
  meta_norm <- set_profile_factors(meta_norm)
  gene_bin_signal <- set_profile_factors(gene_bin_signal)
  gene_signal_summary <- set_profile_factors(gene_signal_summary)

  profile_corr <- compute_profile_correlation_table(meta_norm)
  expr_enrichment <- compute_expr_enrichment_table(meta_norm)
  gene_dominance <- summarize_gene_dominance(gene_signal_summary)
  top_gene_contributors <- build_top_gene_table(gene_signal_summary, top_n = 50L)
  dominance_curve <- build_dominance_curve(gene_signal_summary)

  plots <- build_scope_plots(
    meta_raw = meta_raw,
    meta_norm = meta_norm,
    profile_corr = profile_corr,
    expr_enrichment = expr_enrichment,
    dominance_curve = dominance_curve,
    sample_name = sample_name,
    scope_label = scope_label,
    window_bp = window_bp
  )

  analysis_config <- data.table(
    sample = sample_name,
    chrom_scope = scope_label,
    window_bp = window_bp,
    bin_size = bin_size,
    tss_sources = paste(tss_source_levels, collapse = ","),
    expression_groups = paste(expr_group_levels, collapse = ","),
    feature_levels = paste(feature_levels, collapse = ",")
  )

  write_scope_outputs(
    sample_name = sample_name,
    scope_label = scope_label,
    feature_manifest = feature_manifest,
    analysis_config = analysis_config,
    meta_raw = meta_raw,
    meta_norm = meta_norm,
    profile_corr = profile_corr,
    expr_enrichment = expr_enrichment,
    gene_bin_signal = gene_bin_signal,
    gene_signal_summary = gene_signal_summary,
    gene_dominance = gene_dominance,
    top_gene_contributors = top_gene_contributors,
    dominance_curve = dominance_curve,
    plots = plots
  )

  invisible(list(
    sample = sample_name,
    chrom_scope = scope_label,
    feature_manifest = feature_manifest,
    meta_raw = meta_raw,
    meta_norm = meta_norm,
    gene_bin_signal = gene_bin_signal,
    gene_signal_summary = gene_signal_summary,
    gene_dominance = gene_dominance
  ))
}

run_expr_stratified_TSS_profiles <- function() {
  config <- resolve_run_config()
  message(
    "Running sample ", config$sample,
    " for scopes: ", paste(config$scopes, collapse = ", "),
    " | window_bp=", config$window_bp,
    " | bin_size=", config$bin_size
  )

  lapply(config$scopes, function(scope_label) {
    message("Starting scope: ", scope_label)
    run_scope_analysis(
      sample_name = config$sample,
      scope_label = scope_label,
      window_bp = config$window_bp,
      bin_size = config$bin_size
    )
  })
}

run_expr_stratified_TSS_profiles()


# Rscript /project/spott/cshan/fiber-seq/code/expr_stratified_TSS_profiles_with_fiber_features.r sample=AL10_bc2178_19130 chrom=all
# Rscript /project/spott/cshan/fiber-seq/code/expr_stratified_TSS_profiles_with_fiber_features.r sample=AL10_bc2178_19130 chrom=chr5 run_all_chr=false
