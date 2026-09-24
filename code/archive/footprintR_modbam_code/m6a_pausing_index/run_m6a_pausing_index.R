#!/usr/bin/env Rscript

# Compare m6A levels around TSSs for genes with high and low Pol II pausing index.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
  library(BiocParallel)
  library(Rsamtools)
  library(ggplot2)
})

source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/modbam_footprintR_functions.R")

get_env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = default)
  value <- trimws(value)
  if (!nzchar(value)) default else value
}

get_env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value)) default else value
}

get_env_num <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value)) default else value
}

standard_chromosomes <- function(x) {
  x[grepl("^chr([0-9]+|X|Y)$", x)]
}

strip_gene_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

make_pausing_groups <- function(
  pausing_path,
  tss_source_filter = c("ALL", "CAGE", "GENCODE"),
  low_q = 0.25,
  high_q = 0.75
) {
  tss_source_filter <- match.arg(tss_source_filter)
  pausing <- data.table::fread(pausing_path)

  required_cols <- c("PI", "gene_id", "gene_name", "chrom", "strand", "tss")
  missing_cols <- setdiff(required_cols, names(pausing))
  if (length(missing_cols)) {
    stop("Missing required pausing-index columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  pausing[, PI := suppressWarnings(as.numeric(PI))]
  pausing[, gene_id_base := strip_gene_version(gene_id)]
  pausing <- pausing[
    !is.na(PI) &
      !is.na(gene_id_base) &
      chrom %chin% standard_chromosomes(unique(chrom))
  ]

  pausing <- pausing[order(-PI)]
  pausing <- pausing[!duplicated(gene_id_base)]

  if (!"cage_tss" %in% names(pausing)) {
    pausing[, cage_tss := NA_integer_]
  }
  pausing[, cage_tss := suppressWarnings(as.integer(cage_tss))]
  if (!"tss_gencode" %in% names(pausing)) {
    pausing[, tss_gencode := tss]
  }
  pausing[, tss_gencode := suppressWarnings(as.integer(tss_gencode))]
  pausing[, tss := suppressWarnings(as.integer(tss))]

  if (tss_source_filter == "CAGE") {
    pausing <- pausing[!is.na(cage_tss)]
    pausing[, `:=`(tss_anchor = cage_tss, tss_source = "CAGE")]
  } else if (tss_source_filter == "GENCODE") {
    pausing <- pausing[!is.na(tss_gencode)]
    pausing[, `:=`(tss_anchor = tss_gencode, tss_source = "GENCODE")]
  } else {
    pausing[, `:=`(
      tss_anchor = data.table::fifelse(!is.na(cage_tss), cage_tss, tss),
      tss_source = data.table::fifelse(!is.na(cage_tss), "CAGE", "GENCODE")
    )]
  }

  pausing <- pausing[!is.na(tss_anchor)]

  low_cut <- as.numeric(stats::quantile(pausing$PI, low_q, na.rm = TRUE))
  high_cut <- as.numeric(stats::quantile(pausing$PI, high_q, na.rm = TRUE))

  pausing[, pausing_group := data.table::fifelse(
    PI >= high_cut,
    "High PI",
    data.table::fifelse(PI <= low_cut, "Low PI", "Middle PI")
  )]

  pausing_groups <- pausing[pausing_group %chin% c("High PI", "Low PI")]
  pausing_groups[, `:=`(
    tss_anchor = as.integer(tss_anchor),
    low_q = low_q,
    high_q = high_q,
    low_cut = low_cut,
    high_cut = high_cut
  )]

  pausing_groups[]
}

make_tss_windows <- function(pausing_groups, chrom_name, window_bp, seqinfo = NULL) {
  dt <- copy(pausing_groups[chrom == chrom_name])
  if (!nrow(dt)) {
    return(list(promoter_windows = GRanges(), promoter_meta = data.table()))
  }

  dt[, row_id := .I]
  dt[, tss_uid := paste(chrom, tss_anchor, strand, gene_id_base, row_id, sep = "|")]

  tss_gr <- GenomicRanges::GRanges(
    seqnames = dt$chrom,
    ranges = IRanges::IRanges(start = dt$tss_anchor, end = dt$tss_anchor),
    strand = dt$strand
  )

  if (!is.null(seqinfo)) {
    GenomeInfoDb::seqinfo(tss_gr) <- seqinfo[unique(dt$chrom)]
  }

  mcols(tss_gr)$tss_uid <- dt$tss_uid
  mcols(tss_gr)$gene_id_base <- dt$gene_id_base
  mcols(tss_gr)$gene_name <- dt$gene_name
  mcols(tss_gr)$PI <- dt$PI
  mcols(tss_gr)$pausing_group <- dt$pausing_group

  promoter_windows <- GenomicRanges::promoters(
    tss_gr,
    upstream = window_bp,
    downstream = window_bp + 1L
  )
  promoter_windows <- GenomicRanges::trim(promoter_windows)

  promoter_meta <- data.table::data.table(
    tss_uid = dt$tss_uid,
    chromosome = dt$chrom,
    tss_coordinate = dt$tss_anchor,
    tss_strand = dt$strand,
    gene_name_or_tss_id = dt$gene_name,
    tss_id = dt$gene_id_base,
    gene_id = dt$gene_id,
    gene_id_base = dt$gene_id_base,
    gene_name = dt$gene_name,
    PI = dt$PI,
    pausing_group = dt$pausing_group,
    tss_source = dt$tss_source,
    promoter_start = start(promoter_windows),
    promoter_end = end(promoter_windows)
  )

  list(promoter_windows = promoter_windows, promoter_meta = promoter_meta)
}

append_tsv <- function(dt, path, wrote_once) {
  if (!nrow(dt)) return(wrote_once)
  data.table::fwrite(dt, file = path, sep = "\t", append = wrote_once, col.names = !wrote_once)
  TRUE
}

summarize_group_meta <- function(aligned_dt, groups_dt, window_bp, bin_size) {
  bins <- data.table::CJ(
    pausing_group = c("Low PI", "High PI"),
    bin = seq(-window_bp, window_bp, by = bin_size)
  )

  n_by_group <- groups_dt[, .(n_tss = uniqueN(tss_uid)), by = pausing_group]

  if (!nrow(aligned_dt)) {
    out <- merge(bins, n_by_group, by = "pausing_group", all.x = TRUE)
    out[, `:=`(
      modified_calls = 0L,
      total_calls = 0L,
      fraction_modified = NA_real_,
      modified_calls_per_tss = 0
    )]
    return(out[])
  }

  tmp <- copy(aligned_dt)
  tmp[, bin := round(rel_pos / bin_size) * bin_size]
  meta <- tmp[, .(
    modified_calls = sum(modified_calls),
    total_calls = sum(total_calls)
  ), by = .(pausing_group, bin)]

  out <- merge(bins, meta, by = c("pausing_group", "bin"), all.x = TRUE, sort = TRUE)
  out <- merge(out, n_by_group, by = "pausing_group", all.x = TRUE)
  out[is.na(modified_calls), modified_calls := 0L]
  out[is.na(total_calls), total_calls := 0L]
  out[, fraction_modified := data.table::fifelse(total_calls > 0, modified_calls / total_calls, NA_real_)]
  out[, modified_calls_per_tss := data.table::fifelse(n_tss > 0, modified_calls / n_tss, NA_real_)]
  out[]
}

summarize_gene_m6a <- function(aligned_dt, groups_dt) {
  if (!nrow(aligned_dt)) {
    out <- unique(groups_dt[, .(
      tss_uid, gene_id, gene_id_base, gene_name, chromosome, tss_coordinate,
      tss_strand, tss_source, PI, pausing_group
    )])
    out[, `:=`(
      m6A_modified_calls = 0L,
      m6A_unmodified_calls = 0L,
      m6A_total_calls = 0L,
      m6A_positions_with_coverage = 0L,
      m6A_fraction_modified = NA_real_
    )]
    return(out[])
  }

  gene_summary <- aligned_dt[, {
    total <- sum(total_calls)
    modified <- sum(modified_calls)
    .(
      m6A_modified_calls = modified,
      m6A_unmodified_calls = sum(unmodified_calls),
      m6A_total_calls = total,
      m6A_positions_with_coverage = sum(total_calls > 0),
      m6A_fraction_modified = if (total > 0) modified / total else NA_real_
    )
  }, by = .(
    tss_uid, gene_id, gene_id_base, gene_name, chromosome, tss_coordinate,
    tss_strand, tss_source, PI, pausing_group
  )]

  missing <- groups_dt[!gene_summary, on = "tss_uid"]
  if (nrow(missing)) {
    missing <- missing[, .(
      tss_uid, gene_id, gene_id_base, gene_name, chromosome, tss_coordinate,
      tss_strand, tss_source, PI, pausing_group
    )]
    missing[, `:=`(
      m6A_modified_calls = 0L,
      m6A_unmodified_calls = 0L,
      m6A_total_calls = 0L,
      m6A_positions_with_coverage = 0L,
      m6A_fraction_modified = NA_real_
    )]
    gene_summary <- data.table::rbindlist(list(gene_summary, missing), use.names = TRUE)
  }

  gene_summary[order(pausing_group, -PI)]
}

summarize_gene_bin_m6a <- function(aligned_dt, groups_dt, window_bp, bin_size) {
  bins <- data.table::data.table(bin = seq(-window_bp, window_bp, by = bin_size))
  tss_dt <- unique(groups_dt[, .(
    tss_uid, gene_id, gene_id_base, gene_name, chromosome, tss_coordinate,
    tss_strand, tss_source, PI, pausing_group
  )])
  tss_dt[, join_key := 1L]
  bins[, join_key := 1L]
  tss_bins <- merge(tss_dt, bins, by = "join_key", allow.cartesian = TRUE)
  tss_bins[, join_key := NULL]

  if (!nrow(aligned_dt)) {
    tss_bins[, `:=`(
      m6A_modified_calls = 0L,
      m6A_unmodified_calls = 0L,
      m6A_total_calls = 0L,
      m6A_fraction_modified = NA_real_
    )]
    return(tss_bins[])
  }

  tmp <- copy(aligned_dt)
  tmp[, bin := round(rel_pos / bin_size) * bin_size]
  tmp <- tmp[bin >= -window_bp & bin <= window_bp]

  gene_bin <- tmp[, {
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

  gene_bin <- merge(
    tss_bins,
    gene_bin,
    by = c(
      "tss_uid", "gene_id", "gene_id_base", "gene_name", "chromosome",
      "tss_coordinate", "tss_strand", "tss_source", "PI", "pausing_group", "bin"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  gene_bin[is.na(m6A_modified_calls), m6A_modified_calls := 0L]
  gene_bin[is.na(m6A_unmodified_calls), m6A_unmodified_calls := 0L]
  gene_bin[is.na(m6A_total_calls), m6A_total_calls := 0L]
  gene_bin[is.na(m6A_fraction_modified) & m6A_total_calls > 0,
    m6A_fraction_modified := m6A_modified_calls / m6A_total_calls
  ]
  gene_bin[]
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

plot_m6a_pausing_meta <- function(meta_dt, output_pdf, sample_name, window_bp) {
  plot_dt <- copy(meta_dt)
  plot_dt[, pausing_group := factor(pausing_group, levels = c("Low PI", "High PI"))]

  p <- ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(x = bin, y = fraction_modified, color = pausing_group)
  ) +
    ggplot2::geom_line(linewidth = 0.45, alpha = 0.55, na.rm = TRUE) +
    ggplot2::geom_smooth(method = "loess", span = 0.08, se = FALSE, linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.35) +
    ggplot2::scale_color_manual(values = c("Low PI" = "#0072B2", "High PI" = "#D55E00")) +
    ggplot2::scale_x_continuous(
      breaks = seq(-window_bp, window_bp, by = max(1000L, window_bp / 5L)),
      labels = function(x) paste0(x / 1000, " kb")
    ) +
    ggplot2::labs(
      x = "Distance to TSS",
      y = "m6A fraction modified",
      color = NULL,
      title = sprintf("%s m6A around TSS by Pol II pausing index", sample_name)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "top")

  ggplot2::ggsave(output_pdf, p, width = 8, height = 5)
  invisible(p)
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
  invisible(p)
}

plot_gene_bin_histogram <- function(hist_dt, output_pdf, sample_name, window_bp) {
  plot_dt <- copy(hist_dt)
  plot_dt[, pausing_group := factor(pausing_group, levels = c("Low PI", "High PI"))]

  p <- ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(x = bin, y = modified_call_bin, fill = n_genes)
  ) +
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
  invisible(p)
}

sample_name <- get_env_chr("SAMPLE", "AL10_bc2178_19130")
window_bp <- get_env_int("WINDOW_BP", 10000L)
bin_size <- get_env_int("BIN_SIZE", 100L)
mod_prob_threshold <- get_env_num("MOD_PROB_THRESHOLD", 0.9)
chunk_size <- get_env_int("MODBAM_TSS_CHUNK_SIZE", 250L)
low_q <- get_env_num("LOW_Q", 0.25)
high_q <- get_env_num("HIGH_Q", 0.75)
max_modified_hist_bin <- get_env_int("MAX_MODIFIED_HIST_BIN", 10L)
tss_source_filter <- toupper(get_env_chr("TSS_SOURCE_FILTER", "ALL"))
chroms_env <- get_env_chr("CHROMS", "AUTO")

polii_root <- get_env_chr("POLII_ROOT", "/project/spott/cshan/fiber-seq/results/PolII")
pausing_path <- get_env_chr(
  "PAUSING_PATH",
  file.path(polii_root, "annotations", "pausing_index_principal_with_CAGE_TSS_all_genes.tsv")
)
output_root <- get_env_chr(
  "OUTPUT_ROOT",
  file.path("/project/spott/cshan/fiber-seq/results/PolII", "m6a_pausing_index")
)
bam_path <- get_env_chr(
  "BAM",
  file.path(
    "/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples",
    sprintf("%s.5mC.6mA.aligned.phased.bam", sample_name)
  )
)

if (chunk_size < 1L) chunk_size <- 250L
if (low_q <= 0 || low_q >= 1 || high_q <= 0 || high_q >= 1 || low_q >= high_q) {
  stop("LOW_Q and HIGH_Q must be between 0 and 1 with LOW_Q < HIGH_Q.", call. = FALSE)
}

assert_file_exists(pausing_path, "pausing index table")
assert_file_exists(bam_path, "BAM")
assert_file_exists(paste0(bam_path, ".bai"), "BAM index")

output_dir <- file.path(
  output_root,
  sprintf(
    "%s_%skb_bin%s_modthresh%s_%s",
    sample_name,
    window_bp / 1000,
    bin_size,
    mod_prob_threshold,
    tolower(tss_source_filter)
  )
)
table_dir <- file.path(output_dir, "tables")
plot_dir <- file.path(output_dir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading pausing groups from: ", pausing_path)
pausing_groups <- make_pausing_groups(
  pausing_path = pausing_path,
  tss_source_filter = tss_source_filter,
  low_q = low_q,
  high_q = high_q
)

if (!nrow(pausing_groups)) {
  stop("No high/low pausing-index genes remained after filtering.", call. = FALSE)
}

if (chroms_env == "AUTO") {
  chroms <- standard_chromosomes(sort(unique(pausing_groups$chrom)))
} else {
  chroms <- trimws(strsplit(chroms_env, ",", fixed = TRUE)[[1]])
  chroms <- chroms[nzchar(chroms)]
}

pausing_groups <- pausing_groups[chrom %chin% chroms]
if (!nrow(pausing_groups)) {
  stop("No high/low pausing-index genes remained on selected chromosomes.", call. = FALSE)
}

data.table::fwrite(
  pausing_groups[, .(
    gene_name, gene_id, gene_id_base, chrom, strand, tss_anchor, tss_source,
    cage_tss, tss_gencode, PI, pausing_group, low_q, high_q, low_cut, high_cut
  )],
  file.path(table_dir, sprintf("%s_high_low_pausing_gene_groups.tsv", sample_name)),
  sep = "\t"
)

message("Selected chromosomes: ", paste(chroms, collapse = ", "))
message("High PI genes: ", nrow(pausing_groups[pausing_group == "High PI"]))
message("Low PI genes: ", nrow(pausing_groups[pausing_group == "Low PI"]))
message("PI low cutoff: ", unique(pausing_groups$low_cut)[1])
message("PI high cutoff: ", unique(pausing_groups$high_cut)[1])

bam_seqinfo <- read_bam_seqinfo(bam_path)
bamfiles <- setNames(bam_path, sample_name)
bpparam_modbam <- BiocParallel::SerialParam()

all_gene_summary <- list()
all_gene_bin_summary <- list()
all_meta <- list()
aligned_path <- file.path(table_dir, sprintf("%s_m6A_tss_aligned_positions.tsv.gz", sample_name))
if (file.exists(aligned_path)) file.remove(aligned_path)
wrote_aligned <- FALSE

for (chrom in chroms) {
  message("Processing ", chrom)
  tss_objects <- make_tss_windows(pausing_groups, chrom, window_bp, seqinfo = bam_seqinfo)
  promoter_windows <- tss_objects$promoter_windows
  promoter_meta <- tss_objects$promoter_meta

  if (!nrow(promoter_meta)) {
    next
  }

  chrom_groups <- copy(promoter_meta)
  n_tss <- nrow(promoter_meta)
  chunk_starts <- seq.int(1L, n_tss, by = chunk_size)
  chrom_aligned <- vector("list", length(chunk_starts))

  for (i in seq_along(chunk_starts)) {
    from <- chunk_starts[i]
    to <- min(from + chunk_size - 1L, n_tss)
    message(sprintf("%s chunk %d/%d: TSS %d-%d", chrom, i, length(chunk_starts), from, to))

    chunk_windows <- promoter_windows[from:to]
    chunk_meta <- promoter_meta[from:to]

    chunk_positions <- read_modbam_summary_positions(
      bamfiles = bamfiles,
      sample_name = sample_name,
      promoter_windows = chunk_windows,
      modbase_code = "a",
      seqinfo = bam_seqinfo,
      bpparam = bpparam_modbam,
      modProbThreshold = mod_prob_threshold
    )

    chunk_aligned <- align_positions_to_tss(
      position_dt = chunk_positions,
      promoter_windows = chunk_windows,
      promoter_meta = chunk_meta
    )

    chrom_aligned[[i]] <- chunk_aligned
    wrote_aligned <- append_tsv(chunk_aligned, aligned_path, wrote_aligned)

    rm(chunk_positions, chunk_aligned)
    gc()
  }

  chrom_aligned_dt <- data.table::rbindlist(chrom_aligned, use.names = TRUE, fill = TRUE)
  all_gene_summary[[chrom]] <- summarize_gene_m6a(chrom_aligned_dt, chrom_groups)
  all_gene_bin_summary[[chrom]] <- summarize_gene_bin_m6a(chrom_aligned_dt, chrom_groups, window_bp, bin_size)
  all_meta[[chrom]] <- summarize_group_meta(chrom_aligned_dt, chrom_groups, window_bp, bin_size)
}

gene_summary <- data.table::rbindlist(all_gene_summary, use.names = TRUE, fill = TRUE)
gene_bin_summary <- data.table::rbindlist(all_gene_bin_summary, use.names = TRUE, fill = TRUE)
gene_bin_histogram <- summarize_gene_bin_histogram(gene_bin_summary, max_modified_hist_bin)
group_meta <- data.table::rbindlist(all_meta, use.names = TRUE, fill = TRUE)
group_n_tss <- gene_summary[, .(n_tss = uniqueN(tss_uid)), by = pausing_group]

group_meta <- group_meta[, .(
  modified_calls = sum(modified_calls),
  total_calls = sum(total_calls)
), by = .(pausing_group, bin)]
group_meta <- merge(group_meta, group_n_tss, by = "pausing_group", all.x = TRUE)
group_meta[, fraction_modified := data.table::fifelse(total_calls > 0, modified_calls / total_calls, NA_real_)]
group_meta[, modified_calls_per_tss := data.table::fifelse(n_tss > 0, modified_calls / n_tss, NA_real_)]

gene_summary_path <- file.path(table_dir, sprintf("%s_gene_m6A_summary_high_low_PI.tsv", sample_name))
gene_bin_summary_path <- file.path(table_dir, sprintf("%s_gene_bin_m6A_summary_high_low_PI.tsv.gz", sample_name))
gene_bin_histogram_path <- file.path(table_dir, sprintf("%s_gene_bin_m6A_modified_call_histogram_high_low_PI.tsv", sample_name))
meta_path <- file.path(table_dir, sprintf("%s_m6A_metaprofile_high_low_PI.tsv", sample_name))
group_summary_path <- file.path(table_dir, sprintf("%s_m6A_group_summary_high_low_PI.tsv", sample_name))
plot_path <- file.path(plot_dir, sprintf("%s_m6A_metaprofile_high_low_PI.pdf", sample_name))
boxplot_path <- file.path(plot_dir, sprintf("%s_gene_m6A_fraction_boxplot_high_low_PI.pdf", sample_name))
histogram_path <- file.path(plot_dir, sprintf("%s_gene_bin_m6A_modified_call_histogram_high_low_PI.pdf", sample_name))

data.table::fwrite(gene_summary, gene_summary_path, sep = "\t")
data.table::fwrite(gene_bin_summary[order(pausing_group, gene_id_base, bin)], gene_bin_summary_path, sep = "\t")
data.table::fwrite(gene_bin_histogram[order(pausing_group, bin, modified_call_bin)], gene_bin_histogram_path, sep = "\t")
data.table::fwrite(group_meta[order(pausing_group, bin)], meta_path, sep = "\t")

group_summary <- gene_summary[, .(
  n_genes = uniqueN(gene_id_base),
  n_tss = uniqueN(tss_uid),
  median_PI = stats::median(PI, na.rm = TRUE),
  mean_m6A_fraction = mean(m6A_fraction_modified, na.rm = TRUE),
  median_m6A_fraction = stats::median(m6A_fraction_modified, na.rm = TRUE),
  mean_total_calls = mean(m6A_total_calls, na.rm = TRUE),
  median_total_calls = stats::median(m6A_total_calls, na.rm = TRUE)
), by = pausing_group]

data.table::fwrite(group_summary, group_summary_path, sep = "\t")
plot_m6a_pausing_meta(group_meta, plot_path, sample_name, window_bp)
plot_gene_m6a_boxplot(gene_summary, boxplot_path, sample_name)
plot_gene_bin_histogram(gene_bin_histogram, histogram_path, sample_name, window_bp)

message("Wrote gene summary: ", gene_summary_path)
message("Wrote gene-bin summary: ", gene_bin_summary_path)
message("Wrote gene-bin histogram: ", gene_bin_histogram_path)
message("Wrote metaprofile table: ", meta_path)
message("Wrote group summary: ", group_summary_path)
message("Wrote plot: ", plot_path)
message("Wrote boxplot: ", boxplot_path)
message("Wrote histogram: ", histogram_path)
