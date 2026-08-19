source("/project/spott/cshan/fiber-seq/plot_nuc_distribution_functions.R")

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
})

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

  list(
    sample = args[["sample"]] %||% Sys.getenv("FIBER_SAMPLE", "AL10_bc2178_19130"),
    window_bp = as.integer(args[["window_bp"]] %||% Sys.getenv("FIBER_WINDOW_BP", "2000")),
    bin_size = as.integer(args[["bin_size"]] %||% Sys.getenv("FIBER_BIN_SIZE", "10")),
    mode = tolower(args[["mode"]] %||% Sys.getenv("FIBER_MODE", "combine")),
    chrom = args[["chrom"]] %||% Sys.getenv("FIBER_CHROM", ""),
    chroms = args[["chroms"]] %||% Sys.getenv("FIBER_CHROMS", ""),
    chrom_set = tolower(args[["chrom_set"]] %||% Sys.getenv("FIBER_CHROM_SET", "autosomes_xy"))
  )
}

extract_chrom_from_nuc_path <- function(path) {
  # Example: <sample>.ft_extracted_nuc.chr1.bed.gz -> chr1
  sub("^.*\\.ft_extracted_nuc\\.(chr[^.]+)\\.bed\\.gz$", "\\1", basename(path))
}

discover_nuc_bed_by_chr <- function(sample_root, sample_name) {
  nuc_by_chr_dir <- file.path(sample_root, sample_name, "extracted_results", "nuc_by_chr")
  bed_files <- sort(list.files(
    nuc_by_chr_dir,
    pattern = "ft_extracted_nuc\\.chr.*\\.bed\\.gz$",
    full.names = TRUE
  ))

  if (!length(bed_files)) {
    stop("No chromosome-level nucleosome BED files found in ", nuc_by_chr_dir)
  }

  chroms <- vapply(bed_files, extract_chrom_from_nuc_path, character(1))
  if (any(!nzchar(chroms))) {
    stop("Failed to parse chromosome names from some nucleosome BED filenames under ", nuc_by_chr_dir)
  }

  dt <- data.table(chrom = chroms, bed_path = bed_files)
  dt <- dt[!duplicated(chrom)]
  dt[order(chrom)]
}

read_nuc_bed12 <- function(path) {
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

make_meta_per_read_coverage_normalized <- function(bed12, tss_gr,
                                                   window_bp = 2000L,
                                                   bin_size = 10L,
                                                   label = NULL) {
  plot_dt <- find_nucs_near_tss_per_read(bed12, tss_gr, window_bp = window_bp)
  plot_dt[, bin := round(rel_pos / bin_size) * bin_size]

  template <- data.table(bin = seq(-window_bp, window_bp, by = bin_size))
  meta <- if (nrow(plot_dt)) {
    plot_dt[, .(
      N = .N,
      n_reads = uniqueN(RID),
      n_tss = uniqueN(tss_idx)
    ), by = bin]
  } else {
    data.table(bin = integer(), N = integer(), n_reads = integer(), n_tss = integer())
  }

  meta <- merge(template, meta, by = "bin", all.x = TRUE, sort = TRUE)
  meta[is.na(N), `:=`(N = 0L, n_reads = 0L, n_tss = 0L)]
  meta[, N_per_read := fifelse(n_reads > 0L, N / n_reads, 0)]

  if (!is.null(label)) meta[, source := label]
  meta[]
}

build_gene_summary <- function(rtp_nuc_with_gene) {
  if (!nrow(rtp_nuc_with_gene)) {
    return(data.table(
      gene_id = character(),
      mean_rel_pos_bp = numeric(),
      median_rel_pos_bp = numeric(),
      mean_abs_rel_pos_bp = numeric(),
      median_abs_rel_pos_bp = numeric(),
      sd_rel_pos_bp = numeric(),
      n_nuc = integer(),
      n_reads = integer(),
      n_tss = integer()
    ))
  }

  rtp_nuc_with_gene[, .(
    mean_rel_pos_bp = mean(rel_pos, na.rm = TRUE),
    median_rel_pos_bp = as.integer(median(rel_pos, na.rm = TRUE)),
    mean_abs_rel_pos_bp = mean(abs_dist, na.rm = TRUE),
    median_abs_rel_pos_bp = as.integer(median(abs_dist, na.rm = TRUE)),
    sd_rel_pos_bp = sd(rel_pos, na.rm = TRUE),
    n_nuc = .N,
    n_reads = uniqueN(RID),
    n_tss = uniqueN(tss_idx)
  ), by = gene_id][order(mean_abs_rel_pos_bp, -n_nuc)]
}

build_plots <- function(meta_dt, gene_summary, sample_name, window_bp) {
  raw_plot <- ggplot(meta_dt, aes(x = bin, y = N)) +
    geom_line(linewidth = 0.6, color = "#2ca02c") +
    geom_smooth(method = "loess", span = 0.05, se = FALSE, color = "#2ca02c", linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_x_continuous(
      breaks = seq(-window_bp, window_bp, 500),
      labels = function(x) paste0(x / 1000, " kb")
    ) +
    labs(
      x = "Distance to TSS",
      y = "Nucleosome count",
      title = "Per-read raw counts around CAGE TSS",
      subtitle = sprintf("Sample: %s | all chromosomes | each read-TSS pair counted", sample_name)
    ) +
    theme_bw(base_size = 12)

  corrected_plot <- ggplot(meta_dt, aes(x = bin, y = N_per_read)) +
    geom_line(linewidth = 0.6, color = "#d62728") +
    geom_smooth(method = "loess", span = 0.05, se = FALSE, color = "#d62728", linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_x_continuous(
      breaks = seq(-window_bp, window_bp, 500),
      labels = function(x) paste0(x / 1000, " kb")
    ) +
    labs(
      x = "Distance to TSS",
      y = "Mean nucleosomes per contributing read",
      title = "Per-read coverage-normalised CAGE TSS meta-profile",
      subtitle = sprintf("Sample: %s | all chromosomes | N / n_reads per bin", sample_name)
    ) +
    theme_bw(base_size = 12)

  compare_dt <- rbindlist(list(
    meta_dt[, .(
      bin,
      value = fifelse(max(N, na.rm = TRUE) > 0, N / max(N, na.rm = TRUE), 0),
      method = "Uncorrected (N)"
    )],
    meta_dt[, .(
      bin,
      value = fifelse(max(N_per_read, na.rm = TRUE) > 0, N_per_read / max(N_per_read, na.rm = TRUE), 0),
      method = "Corrected (N / n_reads)"
    )]
  ))

  compare_plot <- ggplot(compare_dt, aes(x = bin, y = value, color = method)) +
    geom_line(linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_color_manual(values = c(
      "Uncorrected (N)" = "#2ca02c",
      "Corrected (N / n_reads)" = "#d62728"
    )) +
    scale_x_continuous(
      breaks = seq(-window_bp, window_bp, 500),
      labels = function(x) paste0(x / 1000, " kb")
    ) +
    labs(
      x = "Distance to TSS",
      y = "Scaled signal (max = 1)",
      color = NULL,
      title = "Effect of coverage normalisation on per-read meta-profile",
      subtitle = "Both traces are scaled to max = 1 to compare shape only"
    ) +
    theme_bw(base_size = 12)

  coverage_plot <- ggplot(meta_dt, aes(x = bin, y = n_reads)) +
    geom_line(linewidth = 0.7, color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_x_continuous(
      breaks = seq(-window_bp, window_bp, 500),
      labels = function(x) paste0(x / 1000, " kb")
    ) +
    labs(
      x = "Distance to TSS",
      y = "n_reads (unique reads contributing to bin)",
      title = "Read coverage by distance",
      subtitle = "This is the denominator removed by N / n_reads"
    ) +
    theme_bw(base_size = 12)

  gene_plot <- ggplot(gene_summary, aes(x = mean_abs_rel_pos_bp)) +
    geom_histogram(binwidth = 25, fill = "#be1c54ff", color = "white") +
    labs(
      x = "Mean absolute distance to TSS per gene (bp)",
      y = "Number of genes",
      title = "Per-gene average nucleosome distance to CAGE TSS",
      subtitle = sprintf("Sample: %s | all chromosomes | genes = %s",
                         sample_name, format(nrow(gene_summary), big.mark = ","))
    ) +
    theme_bw(base_size = 12)

  reads_vs_distance_plot <- ggplot(
    gene_summary,
    aes(x = n_reads, y = mean_abs_rel_pos_bp, size = n_nuc, color = n_tss)
  ) +
    geom_point(alpha = 0.6) +
    scale_size_continuous(name = "n nucleosomes") +
    scale_color_continuous(name = "n TSS", type = "viridis") +
    labs(
      x = "Number of reads covering gene TSS",
      y = "Mean absolute distance to TSS (bp)",
      title = "Per-gene nucleosome positioning by read coverage",
      subtitle = "Size = total nucleosomes; Color = number of TSS"
    ) +
    theme_bw(base_size = 12)

  list(
    raw_plot = raw_plot,
    corrected_plot = corrected_plot,
    compare_plot = compare_plot,
    coverage_plot = coverage_plot,
    gene_plot = gene_plot,
    reads_vs_distance_plot = reads_vs_distance_plot
  )
}

select_analysis_chroms <- function(available_chroms, config) {
  available_chroms <- sort(unique(available_chroms))
  default_chroms <- c(paste0("chr", 1:22), "chrX", "chrY")

  chroms_value <- trimws(config$chroms %||% "")
  requested <- if (nzchar(chroms_value)) {
    unique(trimws(unlist(strsplit(chroms_value, ",", fixed = TRUE))))
  } else if (identical(config$chrom_set, "all")) {
    available_chroms
  } else if (config$chrom_set %in% c("autosomes_xy", "default", "human")) {
    default_chroms
  } else {
    default_chroms
  }

  requested <- intersect(requested, available_chroms)
  if (!length(requested)) stop("No requested chromosomes were found in the nucleosome manifest.")
  requested
}

chunk_output_paths <- function(base_dir, sample_name, chrom) {
  chrom_dir <- file.path(base_dir, sample_name, "chrom_chunks", chrom)
  list(
    chrom_dir = chrom_dir,
    meta = file.path(chrom_dir, sprintf("per_read_%s_%s_meta.tsv.gz", sample_name, chrom)),
    bed_manifest = file.path(chrom_dir, sprintf("per_read_%s_%s_input_bed.tsv", sample_name, chrom)),
    analysis_config = file.path(chrom_dir, sprintf("per_read_%s_%s_config.tsv", sample_name, chrom)),
    partial_rds = file.path(chrom_dir, sprintf("per_read_%s_%s_partial_objects.rds", sample_name, chrom))
  )
}

file_is_nonempty <- function(path) file.exists(path) && isTRUE(file.info(path)$size > 0)

chunk_is_complete <- function(base_dir, sample_name, chrom) {
  paths <- chunk_output_paths(base_dir, sample_name, chrom)
  file_is_nonempty(paths$partial_rds) && file_is_nonempty(paths$meta)
}

run_chunk <- function(sample_root, tss_gr, nuc_manifest, config, out_base) {
  if (!nzchar(config$chrom)) stop("mode=chunk requires chrom=<chromosome>.")
  chrom <- config$chrom

  if (!chrom %in% nuc_manifest$chrom) {
    stop("Chromosome ", chrom, " not found in nucleosome manifest for sample ", config$sample)
  }

  paths <- chunk_output_paths(out_base, config$sample, chrom)
  dir.create(paths$chrom_dir, recursive = TRUE, showWarnings = FALSE)

  if (chunk_is_complete(out_base, config$sample, chrom)) {
    message("Reusing cached chunk: ", chrom)
    return(invisible(readRDS(paths$partial_rds)))
  }

  chrom_id <- chrom
  bed_path <- nuc_manifest[chrom == chrom_id, bed_path][1]
  bed12 <- read_nuc_bed12(bed_path)

  meta_dt <- make_meta_per_read_coverage_normalized(
    bed12,
    tss_gr,
    window_bp = config$window_bp,
    bin_size = config$bin_size,
    label = "CAGE TSS"
  )

  fwrite(
    meta_dt,
    paths$meta,
    sep = "\t"
  )
  fwrite(
    data.table(bed_path = bed_path),
    paths$bed_manifest,
    sep = "\t"
  )
  fwrite(
    data.table(
      sample = config$sample,
      chrom = chrom,
      window_bp = config$window_bp,
      bin_size = config$bin_size,
      tss_source = "CAGE",
      total_tss = length(tss_gr)
    ),
    paths$analysis_config,
    sep = "\t"
  )

  obj <- list(
    sample = config$sample,
    chrom = chrom,
    meta = meta_dt
  )
  saveRDS(obj, paths$partial_rds)
  invisible(obj)
}

combine_chunks <- function(config, nuc_manifest, out_base, plot_dir) {
  chroms <- select_analysis_chroms(nuc_manifest$chrom, config)
  final_dir <- file.path(out_base, config$sample, "all_chr")
  dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  partial_paths <- file.path(
    out_base, config$sample, "chrom_chunks", chroms,
    sprintf("per_read_%s_%s_partial_objects.rds", config$sample, chroms)
  )
  missing <- partial_paths[!file.exists(partial_paths)]
  if (length(missing)) {
    stop(
      "Missing chunk files. Run mode=chunk for these chromosomes first:\n",
      paste(missing, collapse = "\n")
    )
  }

  metas <- rbindlist(lapply(partial_paths, function(p) readRDS(p)$meta), use.names = TRUE, fill = TRUE)
  meta_sum <- metas[, .(
    N = sum(N, na.rm = TRUE),
    n_reads = sum(n_reads, na.rm = TRUE),
    n_tss = sum(n_tss, na.rm = TRUE)
  ), by = .(bin, source)][order(bin)]
  meta_sum[, N_per_read := fifelse(n_reads > 0L, N / n_reads, 0)]

  # No gene-level combine here; chunking is primarily for metaprofile + coverage correction.
  gene_summary <- data.table(
    gene_id = character(),
    mean_rel_pos_bp = numeric(),
    median_rel_pos_bp = numeric(),
    mean_abs_rel_pos_bp = numeric(),
    median_abs_rel_pos_bp = numeric(),
    sd_rel_pos_bp = numeric(),
    n_nuc = integer(),
    n_reads = integer(),
    n_tss = integer()
  )

  plots <- build_plots(
    meta_dt = meta_sum,
    gene_summary = gene_summary,
    sample_name = config$sample,
    window_bp = config$window_bp
  )

  fwrite(
    data.table(
      sample = config$sample,
      chrom_scope = "all_chr",
      window_bp = config$window_bp,
      bin_size = config$bin_size,
      tss_source = "CAGE",
      n_chromosomes = length(chroms)
    ),
    file.path(final_dir, sprintf("per_read_%s_all_chr_config.tsv", config$sample)),
    sep = "\t"
  )
  fwrite(
    data.table(chrom = chroms),
    file.path(final_dir, sprintf("per_read_%s_all_chr_chunk_manifest.tsv", config$sample)),
    sep = "\t"
  )
  fwrite(
    meta_sum,
    file.path(final_dir, sprintf("per_read_%s_all_chr_meta.tsv.gz", config$sample)),
    sep = "\t"
  )

  ggsave(
    file.path(plot_dir, sprintf("per_read_%s_all_chr_raw.pdf", config$sample)),
    plot = plots$raw_plot,
    width = 7,
    height = 5
  )
  ggsave(
    file.path(plot_dir, sprintf("per_read_%s_all_chr_corrected.pdf", config$sample)),
    plot = plots$corrected_plot,
    width = 7,
    height = 5
  )
  ggsave(
    file.path(plot_dir, sprintf("per_read_%s_all_chr_shape_compare.pdf", config$sample)),
    plot = plots$compare_plot,
    width = 7,
    height = 4.5
  )
  ggsave(
    file.path(plot_dir, sprintf("per_read_%s_all_chr_coverage_curve.pdf", config$sample)),
    plot = plots$coverage_plot,
    width = 7,
    height = 4.5
  )

  saveRDS(
    list(
      analysis_config = data.table(
        sample = config$sample,
        chrom_scope = "all_chr",
        window_bp = config$window_bp,
        bin_size = config$bin_size,
        tss_source = "CAGE"
      ),
      meta = meta_sum
    ),
    file.path(final_dir, sprintf("per_read_%s_all_chr_objects.rds", config$sample))
  )

  invisible(meta_sum)
}

run_all_chr_per_read_CAGE_TSS <- function() {
  config <- resolve_run_config()

  sample_root <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results"
  tss_path <- "/project/spott/cshan/annotations/fantom5/fantom5.hg38.LCL.consensus.CAGE_peaks.withGene.bed.gz"
  plot_dir <- "/project/spott/cshan/fiber-seq/results/nuc/per_read_CAGE_TSS_all_chr/plots"
  out_base <- "/project/spott/cshan/fiber-seq/results/nuc/per_read_CAGE_TSS_all_chr"

  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  message(
    "Running all-chromosome per-read CAGE TSS profile for sample ", config$sample,
    " | window_bp=", config$window_bp,
    " | bin_size=", config$bin_size
  )

  tss_gr <- load_tss_annotations_bed9_with_gene_id(tss_path)
  tss_gr <- tss_gr[!is.na(mcols(tss_gr)$gene_id) & mcols(tss_gr)$gene_id != "NA"]

  nuc_manifest <- discover_nuc_bed_by_chr(sample_root, config$sample)

  if (config$mode == "chunk") {
    return(run_chunk(sample_root, tss_gr, nuc_manifest, config, out_base))
  }

  if (config$mode == "combine") {
    return(combine_chunks(config, nuc_manifest, out_base, plot_dir))
  }

  if (config$mode == "full") {
    chroms <- select_analysis_chroms(nuc_manifest$chrom, config)
    for (chrom in chroms) {
      config_one <- config
      config_one$chrom <- chrom
      run_chunk(sample_root, tss_gr, nuc_manifest, config_one, out_base)
      gc(verbose = FALSE)
    }
    return(combine_chunks(config, nuc_manifest, out_base, plot_dir))
  }

  stop("Unsupported mode: ", config$mode, ". Use mode=chunk, mode=combine, or mode=full.")
}

run_all_chr_per_read_CAGE_TSS()

# Example:
# Rscript /project/spott/cshan/fiber-seq/code/CAGE_TSS_per_read_all_chr/all_chr_per_read_CAGE_TSS_nucleosome_meta_profile.r sample=AL10_bc2178_19130 mode=chunk chrom=chr1
# Rscript /project/spott/cshan/fiber-seq/code/CAGE_TSS_per_read_all_chr/all_chr_per_read_CAGE_TSS_nucleosome_meta_profile.r sample=AL10_bc2178_19130 mode=combine
