source("/project/spott/cshan/fiber-seq/plot_nuc_distribution_functions.R")

# ── Paths ─────────────────────────────────────────────────────────────────────
cage_tss_path <- "/project/spott/cshan/annotations/fantom5/fantom5.hg38.LCL.consensus.CAGE_peaks.withGene.bed.gz"

PEAKS_DIR <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged/combined/joint_trained_peaks/"
FP_FILES <- list(
  "10-30"   = file.path(PEAKS_DIR, "combined_all_chrs_10-30bp_fps.bed.gz"),
  "30-45"   = file.path(PEAKS_DIR, "combined_all_chrs_30-45bp_fps.bed.gz"),
  "45-60"   = file.path(PEAKS_DIR, "combined_all_chrs_45-60bp_fps.bed.gz"),
  "60-80"   = file.path(PEAKS_DIR, "combined_all_chrs_60-80bp_fps.bed.gz"),
  "140-160" = file.path(PEAKS_DIR, "combined_all_chrs_140-160bp_fps.bed.gz")
)

default_sample <- "AL10_bc2178_19130"
default_chrom  <- "chr1"
window_bp      <- 2000L
bin_size       <- 10L

out_dir        <- "/project/spott/cshan/fiber-seq/results/PolII"

fp_size_levels <- c("10-30", "30-45", "45-60", "60-80", "140-160")
fp_colors <- c(
  "10-30"   = "#1f77b4",
  "30-45"   = "#ff7f0e",
  "45-60"   = "#2ca02c",
  "60-80"   = "#d62728",
  "140-160" = "#9467bd"
)

# ── Argument parsing helpers (same pattern as example script) ─────────────────
`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) y else x
}

parse_named_args <- function(args) {
  if (!length(args)) return(list())
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) next
    key   <- sub("=.*$", "", arg)
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
  args        <- parse_named_args(commandArgs(trailingOnly = TRUE))
  sample_name <- args[["sample"]]     %||% Sys.getenv("FIBER_SAMPLE",    default_sample)
  chrom_value <- args[["chrom"]]      %||% Sys.getenv("FIBER_CHROM",     default_chrom)
  window_value <- as.integer(args[["window_bp"]] %||% Sys.getenv("FIBER_WINDOW_BP", as.character(window_bp)))
  bin_value    <- as.integer(args[["bin_size"]]  %||% Sys.getenv("FIBER_BIN_SIZE",  as.character(bin_size)))
  run_all_chr  <- parse_bool_flag(
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
    sample     = sample_name,
    scopes     = scope_labels,
    window_bp  = window_value,
    bin_size   = bin_value,
    run_all_chr = run_all_chr
  )
}

# ── Load CAGE TSS (no expression annotation) ─────────────────────────────────
load_cage_tss <- function() {
  load_tss_annotations_bed9_with_gene_id(cage_tss_path)
}

# ── Load footprint peaks from a BED file, return midpoint GRanges ─────────────
# Supports BED3+ (chrom, start, end, ...) – uses peak midpoint as position.
load_footprint_midpoints <- function(fp_path, size_label, chrom_filter = NULL) {
  dt <- fread(fp_path, header = FALSE, select = 1:3,
              col.names = c("chrom", "start", "end"))
  if (!is.null(chrom_filter)) {
    dt <- dt[chrom == chrom_filter]
  }
  if (!nrow(dt)) return(GRanges())
  dt[, midpoint := as.integer((start + end) / 2L)]
  GRanges(
    seqnames = dt$chrom,
    ranges   = IRanges(dt$midpoint + 1L, dt$midpoint + 1L),
    strand   = "*",
    fp_size  = size_label
  )
}

# ── Find footprint midpoints near TSS ────────────────────────────────────────
find_fp_hits_near_tss <- function(fp_gr, tss_gr, window_bp = 2000L, bin_size = 10L) {
  if (!length(fp_gr) || !length(tss_gr)) {
    return(data.table(
      tss_index = integer(), fp_pos = integer(), tss_pos = integer(),
      rel_pos = integer(), bin = integer()
    ))
  }
  tss_windows <- promoters(tss_gr, upstream = window_bp, downstream = window_bp + 1L)
  hits        <- findOverlaps(fp_gr, tss_windows)
  if (!length(hits)) {
    return(data.table(
      tss_index = integer(), fp_pos = integer(), tss_pos = integer(),
      rel_pos = integer(), bin = integer()
    ))
  }
  query_idx <- queryHits(hits)
  subj_idx  <- subjectHits(hits)
  fp_pos    <- start(fp_gr)[query_idx]
  tss_pos   <- start(tss_gr)[subj_idx]
  tss_str   <- as.character(strand(tss_gr))[subj_idx]
  rel_pos   <- ifelse(tss_str == "+", fp_pos - tss_pos, tss_pos - fp_pos)
  bin_col   <- as.integer(round(rel_pos / bin_size) * bin_size)

  dt <- data.table(
    tss_index = subj_idx,
    fp_pos    = as.integer(fp_pos),
    tss_pos   = as.integer(tss_pos),
    rel_pos   = as.integer(rel_pos),
    bin       = bin_col
  )
  dt[bin >= -window_bp & bin <= window_bp]
}

# ── Summarize hits into meta-profile ─────────────────────────────────────────
summarize_fp_meta <- function(hit_dt, n_tss, size_label, sample_name, scope_label,
                               window_bp = 2000L, bin_size = 10L) {
  template <- data.table(bin = seq(-window_bp, window_bp, by = bin_size))

  counts <- if (nrow(hit_dt)) {
    hit_dt[, .(N = .N), by = bin]
  } else {
    data.table(bin = integer(), N = integer())
  }

  meta_raw <- merge(template, counts, by = "bin", all.x = TRUE, sort = TRUE)
  meta_raw[is.na(N), N := 0L]
  meta_raw[, `:=`(sample = sample_name, chrom_scope = scope_label, fp_size = size_label)]

  meta_norm <- copy(meta_raw)
  meta_norm[, N := if (n_tss > 0) N / n_tss else 0]

  list(raw = meta_raw[], norm = meta_norm[])
}

# ── Set factor levels ─────────────────────────────────────────────────────────
set_profile_factors <- function(dt) {
  if (!nrow(dt)) return(dt)
  if ("fp_size" %in% names(dt)) {
    dt[, fp_size := factor(fp_size, levels = fp_size_levels)]
  }
  dt
}

# ── Resolve chromosome filter from scope ─────────────────────────────────────
# Returns NULL for all_chr (load everything), or the chromosome string.
scope_to_chrom_filter <- function(scope_label) {
  if (identical(scope_label, "all_chr")) NULL else scope_label
}

# ── Build plots ──────────────────────────────────────────────────────────────
build_scope_plots <- function(meta_raw, meta_norm, sample_name, scope_label,
                               window_bp = 2000L) {
  x_scale <- scale_x_continuous(
    breaks = seq(-window_bp, window_bp, 500),
    labels = function(x) paste0(x / 1000, " kb")
  )
  scope_subtitle <- sprintf("Sample: %s | %s", sample_name, scope_subtitle_label(scope_label))

  raw_plot <- ggplot(meta_raw, aes(x = bin, y = N, color = fp_size)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.06, se = FALSE, linewidth = 0.9) +
    facet_wrap(~fp_size, scales = "free_y", ncol = 1) +
    scale_color_manual(values = fp_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x     = "Distance to TSS",
      y     = "Footprint midpoints per 10 bp bin",
      color = "Footprint size",
      title = "Footprint size profiles around CAGE TSS: raw counts",
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  norm_plot <- ggplot(meta_norm, aes(x = bin, y = N, color = fp_size)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.06, se = FALSE, linewidth = 0.9) +
    facet_wrap(~fp_size, scales = "free_y", ncol = 1) +
    scale_color_manual(values = fp_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x     = "Distance to TSS",
      y     = "Average footprint midpoints per TSS (10 bp bins)",
      color = "Footprint size",
      title = "Footprint size profiles around CAGE TSS: per-TSS normalized",
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  # Overlay all size classes on one panel (normalized)
  overlay_plot <- ggplot(meta_norm, aes(x = bin, y = N, color = fp_size)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = fp_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x     = "Distance to TSS",
      y     = "Average footprint midpoints per TSS (10 bp bins)",
      color = "Footprint size",
      title = "Footprint size profiles around CAGE TSS: all sizes overlaid",
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  # Z-score overlay – highlights relative shape differences
  meta_z <- copy(meta_norm)
  meta_z[, z := if (sd(N) == 0) 0 else as.numeric(scale(N)),
         by = fp_size]

  zscore_plot <- ggplot(meta_z, aes(x = bin, y = z, color = fp_size)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = fp_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x     = "Distance to TSS",
      y     = "Z-scored normalized signal",
      color = "Footprint size",
      title = "Footprint size profiles around CAGE TSS: Z-score overlay",
      subtitle = "Z-scores computed within each footprint size class"
    ) +
    theme_bw(base_size = 12)

  list(
    raw_plot     = raw_plot,
    norm_plot    = norm_plot,
    overlay_plot = overlay_plot,
    zscore_plot  = zscore_plot
  )
}

# ── Write outputs ─────────────────────────────────────────────────────────────
write_scope_outputs <- function(sample_name, scope_label, n_tss,
                                 meta_raw, meta_norm, analysis_config, plots) {
  plot_dir  <- file.path(out_dir, "plots")
  table_dir <- file.path(out_dir, "tables")
  dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  print(plots$raw_plot)
  print(plots$norm_plot)
  print(plots$overlay_plot)
  print(plots$zscore_plot)

  ggsave(
    file.path(plot_dir, sprintf("%s_%s_footprint_profiles_raw.pdf",     sample_name, scope_label)),
    plot = plots$raw_plot,     width = 8, height = 14
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_%s_footprint_profiles_norm.pdf",    sample_name, scope_label)),
    plot = plots$norm_plot,    width = 8, height = 14
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_%s_footprint_profiles_overlay.pdf", sample_name, scope_label)),
    plot = plots$overlay_plot, width = 8, height = 5
  )
  ggsave(
    file.path(plot_dir, sprintf("%s_%s_footprint_profiles_zscore.pdf",  sample_name, scope_label)),
    plot = plots$zscore_plot,  width = 8, height = 5
  )

  fwrite(meta_raw,  file.path(table_dir, sprintf("%s_%s_footprint_profiles_raw.tsv.gz",  sample_name, scope_label)), sep = "\t")
  fwrite(meta_norm, file.path(table_dir, sprintf("%s_%s_footprint_profiles_norm.tsv.gz", sample_name, scope_label)), sep = "\t")
  fwrite(analysis_config, file.path(table_dir, sprintf("%s_%s_footprint_analysis_config.tsv", sample_name, scope_label)), sep = "\t")

  message(sprintf("Outputs written to %s", out_dir))
}

# ── Per-scope analysis ────────────────────────────────────────────────────────
run_scope_analysis <- function(sample_name, scope_label,
                                window_bp = 2000L, bin_size = 10L) {
  chrom_filter <- scope_to_chrom_filter(scope_label)

  # Load CAGE TSS, optionally filter to chromosome
  tss_gr <- load_cage_tss()
  if (!is.null(chrom_filter)) {
    tss_gr <- tss_gr[as.character(seqnames(tss_gr)) == chrom_filter]
  }
  n_tss <- length(tss_gr)
  message(sprintf("  CAGE TSS loaded: %d sites", n_tss))

  meta_raw_list  <- list()
  meta_norm_list <- list()

  for (size_label in fp_size_levels) {
    fp_path <- FP_FILES[[size_label]]
    message(sprintf("  Loading footprint peaks [%s]: %s", size_label, basename(fp_path)))

    fp_gr <- load_footprint_midpoints(fp_path, size_label, chrom_filter = chrom_filter)
    message(sprintf("    %d footprint midpoints loaded", length(fp_gr)))

    hit_dt <- find_fp_hits_near_tss(fp_gr, tss_gr, window_bp = window_bp, bin_size = bin_size)
    message(sprintf("    %d overlapping hits within +/-%d bp of TSS", nrow(hit_dt), window_bp))

    meta_pair <- summarize_fp_meta(
      hit_dt, n_tss, size_label,
      sample_name = sample_name, scope_label = scope_label,
      window_bp = window_bp, bin_size = bin_size
    )
    meta_raw_list[[size_label]]  <- meta_pair$raw
    meta_norm_list[[size_label]] <- meta_pair$norm
    rm(fp_gr, hit_dt)
    gc()
  }

  meta_raw  <- set_profile_factors(rbindlist(meta_raw_list,  use.names = TRUE, fill = TRUE))
  meta_norm <- set_profile_factors(rbindlist(meta_norm_list, use.names = TRUE, fill = TRUE))

  plots <- build_scope_plots(
    meta_raw   = meta_raw,
    meta_norm  = meta_norm,
    sample_name = sample_name,
    scope_label = scope_label,
    window_bp   = window_bp
  )

  analysis_config <- data.table(
    sample      = sample_name,
    chrom_scope = scope_label,
    n_tss       = n_tss,
    tss_source  = "CAGE",
    window_bp   = window_bp,
    bin_size    = bin_size,
    fp_sizes    = paste(fp_size_levels, collapse = ",")
  )

  write_scope_outputs(
    sample_name     = sample_name,
    scope_label     = scope_label,
    n_tss           = n_tss,
    meta_raw        = meta_raw,
    meta_norm       = meta_norm,
    analysis_config = analysis_config,
    plots           = plots
  )

  invisible(list(
    sample      = sample_name,
    chrom_scope = scope_label,
    n_tss       = n_tss,
    meta_raw    = meta_raw,
    meta_norm   = meta_norm
  ))
}

# ── Entry point ───────────────────────────────────────────────────────────────
run_footprint_TSS_profiles <- function() {
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
      window_bp   = config$window_bp,
      bin_size    = config$bin_size
    )
  })
}

run_footprint_TSS_profiles()

# Rscript /project/spott/cshan/fiber-seq/code/footprints_TSS_profiles_with_fiber_features.r sample=AL10_bc2178_19130 chrom=all


#!/bin/bash
#SBATCH --job-name=ft_FIRE_tss
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=3
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/ft_FIRE_tss%A_%a.err

# export LANG=C
# export LC_ALL=C

# module load R/4.4.1
# Rscript /project/spott/cshan/fiber-seq/code/footprints_TSS_profiles_with_fiber_features.r sample=AL10_bc2178_19130 chrom=chr1 run_all_chr=false
