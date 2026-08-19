#!/usr/bin/env Rscript

# m6A around TSS across all 4 Pol II pausing index quartiles.
#
# Differences from run_m6a_pausing_index.R:
#   1. All genes are kept and divided into 4 PI quartiles (Q1–Q4).
#      No middle 50% is discarded.
#   2. All summaries and plots show all 4 quartile groups.
#   3. Per-read TSS accessibility is computed within a small window
#      (TSS_ACCESS_WINDOW bp, default ±100 bp) using quickread-level
#      modBAM data. Each read-TSS pair is labelled tss_accessible=TRUE/FALSE
#      and written to a separate table for downstream footprint analysis.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
  library(BiocParallel)
  library(Rsamtools)
  library(ggplot2)
})
# footprintR, SummarizedExperiment, and their dependencies are loaded by
# modbam_footprintR_functions.R below. Loading footprintR a second time
# before that source triggers a segfault in its .onLoad hook.

source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/modbam_footprintR_functions.R")


# ── Environment helpers ────────────────────────────────────────────────────────

get_env_chr <- function(name, default) {
  v <- trimws(Sys.getenv(name, unset = default))
  if (!nzchar(v)) default else v
}
get_env_int <- function(name, default) {
  v <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(v)) default else v
}
get_env_num <- function(name, default) {
  v <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (is.na(v)) default else v
}

standard_chromosomes <- function(x) x[grepl("^chr([0-9]+|X|Y)$", x)]
strip_gene_version   <- function(x) sub("\\.[0-9]+$", "", as.character(x))


# ── Quartile palette ───────────────────────────────────────────────────────────

QUARTILE_LEVELS <- c("Q1 (Low PI)", "Q2", "Q3", "Q4 (High PI)")

QUARTILE_COLORS <- c(
  "Q1 (Low PI)"  = "#0072B2",
  "Q2"           = "#56B4E9",
  "Q3"           = "#E69F00",
  "Q4 (High PI)" = "#D55E00"
)


# ── make_pausing_quartiles ─────────────────────────────────────────────────────


make_pausing_quartiles <- function(
  pausing_path,
  tss_source_filter = c("ALL", "CAGE", "GENCODE")
) {
  tss_source_filter <- match.arg(tss_source_filter)
  pausing <- data.table::fread(pausing_path)

  required_cols <- c("PI", "gene_id", "gene_name", "chrom", "strand", "tss")
  missing_cols  <- setdiff(required_cols, names(pausing))
  if (length(missing_cols)) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  pausing[, PI          := suppressWarnings(as.numeric(PI))]
  pausing[, gene_id_base := strip_gene_version(gene_id)]
  pausing <- pausing[
    !is.na(PI) & !is.na(gene_id_base) &
      chrom %chin% standard_chromosomes(unique(chrom))
  ]
  pausing <- pausing[order(-PI)][!duplicated(gene_id_base)]

  if (!"cage_tss"    %in% names(pausing)) pausing[, cage_tss    := NA_integer_]
  if (!"tss_gencode" %in% names(pausing)) pausing[, tss_gencode := tss]
  pausing[, cage_tss    := suppressWarnings(as.integer(cage_tss))]
  pausing[, tss_gencode := suppressWarnings(as.integer(tss_gencode))]
  pausing[, tss         := suppressWarnings(as.integer(tss))]

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

  cuts   <- stats::quantile(pausing$PI, c(0.25, 0.5, 0.75), na.rm = TRUE)
  q1_cut <- as.numeric(cuts[1])
  q2_cut <- as.numeric(cuts[2])
  q3_cut <- as.numeric(cuts[3])

  pausing[, pausing_group := data.table::fcase(
    PI <= q1_cut,                "Q1 (Low PI)",
    PI >  q1_cut & PI <= q2_cut, "Q2",
    PI >  q2_cut & PI <= q3_cut, "Q3",
    PI >  q3_cut,                "Q4 (High PI)"
  )]
  pausing[, pausing_group := factor(pausing_group, levels = QUARTILE_LEVELS)]
  pausing[, `:=`(
    tss_anchor = as.integer(tss_anchor),
    q1_cut = q1_cut, q2_cut = q2_cut, q3_cut = q3_cut
  )]
  pausing[]
}


# ── make_tss_windows ───────────────────────────────────────────────────────────

make_tss_windows <- function(pausing_groups, chrom_name, window_bp, seqinfo = NULL) {
  dt <- data.table::copy(pausing_groups[chrom == chrom_name])
  if (!nrow(dt)) {
    return(list(
      promoter_windows = GenomicRanges::GRanges(),
      promoter_meta    = data.table::data.table()
    ))
  }

  dt[, row_id  := .I]
  dt[, tss_uid := paste(chrom, tss_anchor, strand, gene_id_base, row_id, sep = "|")]

  tss_gr <- GenomicRanges::GRanges(
    seqnames = dt$chrom,
    ranges   = IRanges::IRanges(start = dt$tss_anchor, end = dt$tss_anchor),
    strand   = dt$strand
  )
  if (!is.null(seqinfo)) {
    GenomeInfoDb::seqinfo(tss_gr) <- seqinfo[unique(dt$chrom)]
  }
  GenomicRanges::mcols(tss_gr)$tss_uid       <- dt$tss_uid
  GenomicRanges::mcols(tss_gr)$gene_id_base  <- dt$gene_id_base
  GenomicRanges::mcols(tss_gr)$gene_name     <- dt$gene_name
  GenomicRanges::mcols(tss_gr)$PI            <- dt$PI
  GenomicRanges::mcols(tss_gr)$pausing_group <- as.character(dt$pausing_group)

  promoter_windows <- GenomicRanges::promoters(tss_gr, upstream = window_bp, downstream = window_bp + 1L)
  promoter_windows <- GenomicRanges::trim(promoter_windows)

  promoter_meta <- data.table::data.table(
    tss_uid             = dt$tss_uid,
    chromosome          = dt$chrom,
    tss_coordinate      = dt$tss_anchor,
    tss_strand          = dt$strand,
    gene_name_or_tss_id = dt$gene_name,
    tss_id              = dt$gene_id_base,
    gene_id             = dt$gene_id,
    gene_id_base        = dt$gene_id_base,
    gene_name           = dt$gene_name,
    PI                  = dt$PI,
    pausing_group       = as.character(dt$pausing_group),
    tss_source          = dt$tss_source,
    promoter_start      = start(promoter_windows),
    promoter_end        = end(promoter_windows)
  )

  list(promoter_windows = promoter_windows, promoter_meta = promoter_meta)
}


# ── append_tsv ────────────────────────────────────────────────────────────────

append_tsv <- function(dt, path, wrote_once) {
  if (!nrow(dt)) return(wrote_once)
  data.table::fwrite(dt, file = path, sep = "\t", append = wrote_once, col.names = !wrote_once)
  TRUE
}


# ── summarize_group_meta (4 groups) ───────────────────────────────────────────

summarize_group_meta <- function(aligned_dt, groups_dt, window_bp, bin_size) {
  bins       <- data.table::CJ(pausing_group = QUARTILE_LEVELS,
                                bin           = seq(-window_bp, window_bp, by = bin_size))
  n_by_group <- groups_dt[, .(n_tss = data.table::uniqueN(tss_uid)), by = pausing_group]

  if (!nrow(aligned_dt)) {
    out <- merge(bins, n_by_group, by = "pausing_group", all.x = TRUE)
    out[, `:=`(modified_calls = 0L, total_calls = 0L,
               fraction_modified = NA_real_, modified_calls_per_tss = 0)]
    return(out[])
  }

  tmp  <- data.table::copy(aligned_dt)
  tmp[, bin := round(rel_pos / bin_size) * bin_size]
  meta <- tmp[, .(modified_calls = sum(modified_calls), total_calls = sum(total_calls)),
              by = .(pausing_group, bin)]

  out  <- merge(bins, meta, by = c("pausing_group", "bin"), all.x = TRUE, sort = TRUE)
  out  <- merge(out, n_by_group, by = "pausing_group", all.x = TRUE)
  out[is.na(modified_calls), modified_calls := 0L]
  out[is.na(total_calls),    total_calls    := 0L]
  out[, fraction_modified      := data.table::fifelse(total_calls > 0, modified_calls / total_calls, NA_real_)]
  out[, modified_calls_per_tss := data.table::fifelse(n_tss > 0, modified_calls / n_tss, NA_real_)]
  out[]
}


# ── summarize_gene_m6a (4 groups) ─────────────────────────────────────────────

.KEY_COLS <- c("tss_uid", "gene_id", "gene_id_base", "gene_name", "chromosome",
               "tss_coordinate", "tss_strand", "tss_source", "PI", "pausing_group")

summarize_gene_m6a <- function(aligned_dt, groups_dt) {
  if (!nrow(aligned_dt)) {
    out <- unique(groups_dt[, .KEY_COLS, with = FALSE])
    out[, `:=`(m6A_modified_calls = 0L, m6A_unmodified_calls = 0L,
               m6A_total_calls = 0L, m6A_positions_with_coverage = 0L,
               m6A_fraction_modified = NA_real_)]
    return(out[])
  }

  gene_summary <- aligned_dt[, {
    total    <- sum(total_calls)
    modified <- sum(modified_calls)
    .(
      m6A_modified_calls          = modified,
      m6A_unmodified_calls        = sum(unmodified_calls),
      m6A_total_calls             = total,
      m6A_positions_with_coverage = sum(total_calls > 0),
      m6A_fraction_modified       = if (total > 0) modified / total else NA_real_
    )
  }, by = .KEY_COLS]

  missing <- groups_dt[!tss_uid %chin% gene_summary$tss_uid][, .KEY_COLS, with = FALSE]
  if (nrow(missing)) {
    missing[, `:=`(m6A_modified_calls = 0L, m6A_unmodified_calls = 0L,
                   m6A_total_calls = 0L, m6A_positions_with_coverage = 0L,
                   m6A_fraction_modified = NA_real_)]
    gene_summary <- data.table::rbindlist(list(gene_summary, missing), use.names = TRUE)
  }
  gene_summary[order(pausing_group, -PI)]
}


# ── summarize_gene_bin_m6a (4 groups) ─────────────────────────────────────────

summarize_gene_bin_m6a <- function(aligned_dt, groups_dt, window_bp, bin_size) {
  bins     <- data.table::data.table(bin = seq(-window_bp, window_bp, by = bin_size))
  tss_dt   <- unique(groups_dt[, .KEY_COLS, with = FALSE])
  tss_dt[, join_key := 1L]
  bins[,   join_key := 1L]
  tss_bins <- merge(tss_dt, bins, by = "join_key", allow.cartesian = TRUE)
  tss_bins[, join_key := NULL]

  if (!nrow(aligned_dt)) {
    tss_bins[, `:=`(m6A_modified_calls = 0L, m6A_unmodified_calls = 0L,
                    m6A_total_calls = 0L, m6A_fraction_modified = NA_real_)]
    return(tss_bins[])
  }

  tmp <- data.table::copy(aligned_dt)
  tmp[, bin := round(rel_pos / bin_size) * bin_size]
  tmp <- tmp[bin >= -window_bp & bin <= window_bp]

  gene_bin <- tmp[, {
    total    <- sum(total_calls)
    modified <- sum(modified_calls)
    .(m6A_modified_calls   = modified,
      m6A_unmodified_calls = sum(unmodified_calls),
      m6A_total_calls      = total,
      m6A_fraction_modified = if (total > 0) modified / total else NA_real_)
  }, by = c(.KEY_COLS, "bin")]

  gene_bin <- merge(tss_bins, gene_bin, by = c(.KEY_COLS, "bin"), all.x = TRUE, sort = FALSE)
  gene_bin[is.na(m6A_modified_calls),   m6A_modified_calls   := 0L]
  gene_bin[is.na(m6A_unmodified_calls), m6A_unmodified_calls := 0L]
  gene_bin[is.na(m6A_total_calls),      m6A_total_calls      := 0L]
  gene_bin[is.na(m6A_fraction_modified) & m6A_total_calls > 0,
           m6A_fraction_modified := m6A_modified_calls / m6A_total_calls]
  gene_bin[]
}


# ── summarize_gene_bin_histogram ──────────────────────────────────────────────

summarize_gene_bin_histogram <- function(gene_bin_dt, max_modified_bin = 10L) {
  hist_dt <- data.table::copy(gene_bin_dt)
  hist_dt[, modified_call_bin := data.table::fifelse(
    m6A_modified_calls >= max_modified_bin,
    paste0(max_modified_bin, "+"), as.character(m6A_modified_calls)
  )]
  hist_dt[, modified_call_bin := factor(
    modified_call_bin,
    levels = c(as.character(seq.int(0L, max_modified_bin - 1L)), paste0(max_modified_bin, "+"))
  )]
  hist_dt[, .(n_genes = data.table::uniqueN(gene_id_base), n_tss = data.table::uniqueN(tss_uid)),
          by = .(pausing_group, bin, modified_call_bin)]
}


# ── label_read_tss_accessibility ──────────────────────────────────────────────
# For each read that overlaps the TSS core window, labels whether the read
# carries m6A signal (tss_accessible = TRUE) or not (tss_accessible = FALSE).
# Uses quickread-level modBAM data so individual read IDs are preserved.

label_read_tss_accessibility <- function(
  bamfiles,
  sample_name,
  access_windows,   # GRanges: TSS ± tss_access_window for this chunk
  promoter_meta,    # matching metadata data.table
  seqinfo,
  bpparam,
  modProbThreshold = 0.9
) {
  se_reads <- tryCatch(
    footprintR::readModBam(
      bamfiles         = bamfiles,
      regions          = access_windows,
      modbase          = setNames("a", sample_name),
      level            = "quickread",
      seqinfo          = seqinfo,
      trim             = TRUE,
      BPPARAM          = bpparam,
      verbose          = FALSE
    ),
    error = function(e) { message("  [WARN] quickread failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(se_reads)) return(data.table::data.table())

  se_flat <- tryCatch(
    footprintR::flattenReadLevelAssay(se_reads),
    error = function(e) { message("  [WARN] flattenReadLevelAssay failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(se_flat)) return(data.table::data.table())

  binary_calls <- tryCatch(
    getBinarizedModificationsLongTable(se_flat, sample_name, modProbThreshold),
    error = function(e) { message("  [WARN] getBinarizedModificationsLongTable failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(binary_calls) || nrow(binary_calls) == 0L) return(data.table::data.table())

  calls_dt <- data.table::as.data.table(binary_calls)

  pos_gr <- GenomicRanges::GRanges(
    seqnames = calls_dt$call_chromosome,
    ranges   = IRanges::IRanges(
      start = coerce_iranges_coord(calls_dt$position, "read position"),
      end   = coerce_iranges_coord(calls_dt$position, "read position")
    ),
    strand = "*"
  )

  hits <- GenomicRanges::findOverlaps(pos_gr, access_windows, ignore.strand = TRUE)
  if (!length(hits)) return(data.table::data.table())

  aligned <- cbind(
    calls_dt[S4Vectors::queryHits(hits)],
    promoter_meta[S4Vectors::subjectHits(hits)]
  )
  aligned[, rel_pos := ifelse(
    tss_strand == "+",
    position - tss_coordinate,
    tss_coordinate - position
  )]

  # Per (read_id, tss_uid): count modified and unmodified adenines in window.
  # tss_accessible = TRUE  if >= 1 m6A call (chromatin open, methyltransferase accessible)
  # tss_accessible = FALSE if read is covered but 0 m6A calls (protected/bound)
  read_access <- aligned[, .(
    n_positions_covered = .N,
    n_m6a_modified      = sum(mod_status == 1L, na.rm = TRUE),
    n_m6a_unmodified    = sum(mod_status == 0L, na.rm = TRUE)
  ), by = .(
    read_id, tss_uid, gene_id, gene_id_base, gene_name,
    chromosome, tss_coordinate, tss_strand, tss_source, PI, pausing_group
  )]

  read_access[, tss_accessible := n_m6a_modified > 0L]
  read_access[]
}


# ── Plots ──────────────────────────────────────────────────────────────────────

plot_m6a_quartile_meta <- function(meta_dt, output_pdf, sample_name, window_bp) {
  plot_dt <- data.table::copy(meta_dt)
  plot_dt[, pausing_group := factor(pausing_group, levels = QUARTILE_LEVELS)]

  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = bin, y = fraction_modified, color = pausing_group)) +
    ggplot2::geom_line(linewidth = 0.45, alpha = 0.55, na.rm = TRUE) +
    ggplot2::geom_smooth(method = "loess", span = 0.08, se = FALSE, linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.35) +
    ggplot2::scale_color_manual(values = QUARTILE_COLORS) +
    ggplot2::scale_x_continuous(
      breaks = seq(-window_bp, window_bp, by = max(1000L, window_bp / 5L)),
      labels = function(x) paste0(x / 1000, " kb")
    ) +
    ggplot2::labs(
      x     = "Distance to TSS",
      y     = "m6A fraction modified",
      color = NULL,
      title = sprintf("%s  m6A around TSS  –  4 PI quartiles", sample_name)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "top")

  ggplot2::ggsave(output_pdf, p, width = 8, height = 5)
  invisible(p)
}

plot_gene_m6a_boxplot_quartiles <- function(gene_summary, output_pdf, sample_name) {
  plot_dt <- data.table::copy(gene_summary[!is.na(m6A_fraction_modified)])
  plot_dt[, pausing_group := factor(pausing_group, levels = QUARTILE_LEVELS)]

  kw_p <- tryCatch(
    stats::kruskal.test(m6A_fraction_modified ~ pausing_group, data = plot_dt)$p.value,
    error = function(e) NA_real_
  )
  subtitle <- if (is.na(kw_p)) NULL else sprintf("Kruskal-Wallis p = %.3g", kw_p)

  label_dt <- plot_dt[, .(
    y     = max(m6A_fraction_modified, na.rm = TRUE),
    label = sprintf("n=%s", data.table::uniqueN(gene_id_base))
  ), by = pausing_group]

  p <- ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(x = pausing_group, y = m6A_fraction_modified, fill = pausing_group)
  ) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
    ggplot2::geom_jitter(width = 0.12, height = 0, size = 0.6, alpha = 0.35) +
    ggplot2::geom_text(
      data = label_dt,
      ggplot2::aes(y = y, label = label),
      vjust = -0.6, size = 3.2, show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = QUARTILE_COLORS) +
    ggplot2::labs(
      x        = NULL,
      y        = "Gene-level m6A fraction modified",
      title    = sprintf("%s  gene-level m6A  –  4 PI quartiles", sample_name),
      subtitle = subtitle
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none")

  ggplot2::ggsave(output_pdf, p, width = 7, height = 5)
  invisible(p)
}

plot_gene_bin_histogram_quartiles <- function(hist_dt, output_pdf, sample_name, window_bp) {
  plot_dt <- data.table::copy(hist_dt)
  plot_dt[, pausing_group := factor(pausing_group, levels = QUARTILE_LEVELS)]

  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = bin, y = modified_call_bin, fill = n_genes)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~pausing_group, ncol = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.35) +
    ggplot2::scale_fill_viridis_c(option = "magma", direction = -1) +
    ggplot2::scale_x_continuous(
      breaks = seq(-window_bp, window_bp, by = max(500L, window_bp / 5L)),
      labels = function(x) paste0(x / 1000, " kb")
    ) +
    ggplot2::labs(
      x     = "Distance to TSS",
      y     = "Modified m6A calls per gene per bin",
      fill  = "Genes",
      title = sprintf("%s  m6A call distribution  –  4 PI quartiles", sample_name)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  ggplot2::ggsave(output_pdf, p, width = 12, height = 8)
  invisible(p)
}


# ── Parameters ────────────────────────────────────────────────────────────────

sample_name           <- get_env_chr("SAMPLE",                "AL10_bc2178_19130")
window_bp             <- get_env_int("WINDOW_BP",             10000L)
bin_size              <- get_env_int("BIN_SIZE",              100L)
mod_prob_threshold    <- get_env_num("MOD_PROB_THRESHOLD",    0.9)
chunk_size            <- get_env_int("MODBAM_TSS_CHUNK_SIZE", 250L)
max_modified_hist_bin <- get_env_int("MAX_MODIFIED_HIST_BIN", 10L)
tss_access_window     <- get_env_int("TSS_ACCESS_WINDOW",     100L)  # bp either side of TSS
tss_source_filter     <- toupper(get_env_chr("TSS_SOURCE_FILTER", "ALL"))
chroms_env            <- get_env_chr("CHROMS",                "AUTO")
do_position_scan      <- get_env_int("DO_POSITION_SCAN",      1L) > 0L  # FALSE = skip m6A position scan and the 4 derived summary tables/plots; only the per-read accessibility table is produced. Use when re-running for a different TSS_ACCESS_WINDOW after the position-level outputs have already been written.

polii_root   <- get_env_chr("POLII_ROOT", "/project/spott/cshan/fiber-seq/results/PolII")
pausing_path <- get_env_chr(
  "PAUSING_PATH",
  file.path(polii_root, "annotations", "pausing_index_principal_with_CAGE_TSS_all_genes.tsv")
)
output_root  <- get_env_chr(
  "OUTPUT_ROOT",
  file.path(polii_root, "m6a_pausing_quartiles")
)
bam_path <- get_env_chr(
  "BAM",
  file.path(
    "/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples",
    sprintf("%s.5mC.6mA.aligned.phased.bam", sample_name)
  )
)

if (chunk_size < 1L)        chunk_size        <- 250L
if (tss_access_window < 1L) tss_access_window <- 100L

assert_file_exists(pausing_path,            "pausing index table")
assert_file_exists(bam_path,                "BAM")
assert_file_exists(paste0(bam_path, ".bai"), "BAM index")

output_dir <- file.path(
  output_root,
  sprintf("%s_%skb_bin%s_modthresh%s_%s",
          sample_name, window_bp / 1000, bin_size, mod_prob_threshold,
          tolower(tss_source_filter))
)
table_dir <- file.path(output_dir, "tables")
plot_dir  <- file.path(output_dir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)


# ── Load & partition genes ────────────────────────────────────────────────────

message("Loading pausing quartiles from: ", pausing_path)
pausing_groups <- make_pausing_quartiles(pausing_path, tss_source_filter)

if (!nrow(pausing_groups)) stop("No genes remained after filtering.", call. = FALSE)

if (chroms_env == "AUTO") {
  chroms <- standard_chromosomes(sort(unique(pausing_groups$chrom)))
} else {
  chroms <- trimws(strsplit(chroms_env, ",", fixed = TRUE)[[1]])
  chroms <- chroms[nzchar(chroms)]
}

pausing_groups <- pausing_groups[chrom %chin% chroms]
if (!nrow(pausing_groups)) stop("No genes on selected chromosomes.", call. = FALSE)

data.table::fwrite(
  pausing_groups[, .(
    gene_name, gene_id, gene_id_base, chrom, strand, tss_anchor, tss_source,
    cage_tss, tss_gencode, PI, pausing_group, q1_cut, q2_cut, q3_cut
  )],
  file.path(table_dir, sprintf("%s_quartile_gene_groups.tsv", sample_name)),
  sep = "\t"
)

message(sprintf("TSS_ACCESS_WINDOW = +/- %d bp  |  DO_POSITION_SCAN = %s", tss_access_window, do_position_scan))
message("Selected chromosomes: ", paste(chroms, collapse = ", "))
for (grp in QUARTILE_LEVELS) {
  message(sprintf("  %-16s: %d genes", grp, nrow(pausing_groups[pausing_group == grp])))
}
message(sprintf(
  "PI cuts:  Q1 <= %.3g  |  Q2 <= %.3g  |  Q3 <= %.3g  |  Q4 > %.3g",
  unique(pausing_groups$q1_cut)[1], unique(pausing_groups$q2_cut)[1],
  unique(pausing_groups$q3_cut)[1], unique(pausing_groups$q3_cut)[1]
))


# ── Main loop ──────────────────────────────────────────────────────────────────

bam_seqinfo    <- read_bam_seqinfo(bam_path)
bamfiles       <- setNames(bam_path, sample_name)
bpparam_serial <- BiocParallel::SerialParam()

all_gene_summary     <- list()
all_gene_bin_summary <- list()
all_meta             <- list()

aligned_path <- file.path(table_dir, sprintf("%s_m6A_tss_aligned_positions.tsv.gz", sample_name))
access_path  <- file.path(
  table_dir,
  sprintf("%s_tss_read_accessibility_aw%dbp.tsv.gz", sample_name, tss_access_window)
)
if (do_position_scan && file.exists(aligned_path)) file.remove(aligned_path)
if (file.exists(access_path))  file.remove(access_path)
wrote_aligned <- FALSE
wrote_access  <- FALSE

for (chrom in chroms) {
  message("Processing ", chrom)
  tss_objects      <- make_tss_windows(pausing_groups, chrom, window_bp, seqinfo = bam_seqinfo)
  promoter_windows <- tss_objects$promoter_windows
  promoter_meta    <- tss_objects$promoter_meta
  if (!nrow(promoter_meta)) next

  chrom_groups <- data.table::copy(promoter_meta)
  n_tss        <- nrow(promoter_meta)
  chunk_starts <- seq.int(1L, n_tss, by = chunk_size)
  chrom_aligned <- vector("list", length(chunk_starts))

  for (i in seq_along(chunk_starts)) {
    from <- chunk_starts[i]
    to   <- min(from + chunk_size - 1L, n_tss)
    message(sprintf("  %s chunk %d/%d: TSS %d-%d", chrom, i, length(chunk_starts), from, to))

    chunk_windows <- promoter_windows[from:to]
    chunk_meta    <- promoter_meta[from:to]

    if (do_position_scan) {
      # ── Summary-level m6A across the full promoter window ───────────────
      chunk_positions <- read_modbam_summary_positions(
        bamfiles         = bamfiles,
        sample_name      = sample_name,
        promoter_windows = chunk_windows,
        modbase_code     = "a",
        seqinfo          = bam_seqinfo,
        bpparam          = bpparam_serial,
        modProbThreshold = mod_prob_threshold
      )
      chunk_aligned <- align_positions_to_tss(
        position_dt      = chunk_positions,
        promoter_windows = chunk_windows,
        promoter_meta    = chunk_meta
      )
      chrom_aligned[[i]] <- chunk_aligned
      wrote_aligned <- append_tsv(chunk_aligned, aligned_path, wrote_aligned)
    }

    # ── Per-read accessibility in a small TSS-centred window ──────────────
    # Build GRanges for ±tss_access_window around each TSS in this chunk.
    tss_gr_chunk <- GenomicRanges::GRanges(
      seqnames = chunk_meta$chromosome,
      ranges   = IRanges::IRanges(
        start = chunk_meta$tss_coordinate,
        end   = chunk_meta$tss_coordinate
      ),
      strand = chunk_meta$tss_strand
    )
    GenomeInfoDb::seqinfo(tss_gr_chunk) <- bam_seqinfo[unique(chunk_meta$chromosome)]

    access_windows_chunk <- GenomicRanges::promoters(
      tss_gr_chunk,
      upstream   = tss_access_window,
      downstream = tss_access_window + 1L
    )
    access_windows_chunk <- GenomicRanges::trim(access_windows_chunk)

    access_meta <- data.table::copy(chunk_meta)
    access_meta[, promoter_start := start(access_windows_chunk)]
    access_meta[, promoter_end   := end(access_windows_chunk)]

    chunk_access <- label_read_tss_accessibility(
      bamfiles         = bamfiles,
      sample_name      = sample_name,
      access_windows   = access_windows_chunk,
      promoter_meta    = access_meta,
      seqinfo          = bam_seqinfo,
      bpparam          = bpparam_serial,
      modProbThreshold = mod_prob_threshold
    )
    if (nrow(chunk_access)) {
      wrote_access <- append_tsv(chunk_access, access_path, wrote_access)
    }

    if (exists("chunk_positions")) rm(chunk_positions)
    if (exists("chunk_aligned"))   rm(chunk_aligned)
    rm(chunk_access)
    gc()
  }

  if (do_position_scan) {
    chrom_aligned_dt <- data.table::rbindlist(chrom_aligned, use.names = TRUE, fill = TRUE)
    all_gene_summary[[chrom]]     <- summarize_gene_m6a(chrom_aligned_dt, chrom_groups)
    all_gene_bin_summary[[chrom]] <- summarize_gene_bin_m6a(chrom_aligned_dt, chrom_groups, window_bp, bin_size)
    all_meta[[chrom]]             <- summarize_group_meta(chrom_aligned_dt, chrom_groups, window_bp, bin_size)
  }
}


# ── Aggregate and write ────────────────────────────────────────────────────────

if (do_position_scan) {
  gene_summary     <- data.table::rbindlist(all_gene_summary,     use.names = TRUE, fill = TRUE)
  gene_bin_summary <- data.table::rbindlist(all_gene_bin_summary, use.names = TRUE, fill = TRUE)
  gene_bin_hist    <- summarize_gene_bin_histogram(gene_bin_summary, max_modified_hist_bin)
  group_meta_raw   <- data.table::rbindlist(all_meta,             use.names = TRUE, fill = TRUE)

  group_n_tss <- gene_summary[, .(n_tss = data.table::uniqueN(tss_uid)), by = pausing_group]
  group_meta  <- group_meta_raw[, .(modified_calls = sum(modified_calls),
                                     total_calls    = sum(total_calls)),
                                 by = .(pausing_group, bin)]
  group_meta  <- merge(group_meta, group_n_tss, by = "pausing_group", all.x = TRUE)
  group_meta[, fraction_modified      := data.table::fifelse(total_calls > 0, modified_calls / total_calls, NA_real_)]
  group_meta[, modified_calls_per_tss := data.table::fifelse(n_tss > 0, modified_calls / n_tss, NA_real_)]

  pfx  <- function(s) file.path(table_dir, sprintf("%s_%s", sample_name, s))
  ppfx <- function(s) file.path(plot_dir,  sprintf("%s_%s", sample_name, s))

  gene_summary_path     <- pfx("gene_m6A_summary_quartiles.tsv")
  gene_bin_summary_path <- pfx("gene_bin_m6A_summary_quartiles.tsv.gz")
  histogram_table_path  <- pfx("gene_bin_m6A_histogram_quartiles.tsv")
  meta_path             <- pfx("m6A_metaprofile_quartiles.tsv")
  group_summary_path    <- pfx("m6A_group_summary_quartiles.tsv")

  data.table::fwrite(gene_summary,                                                  gene_summary_path,     sep = "\t")
  data.table::fwrite(gene_bin_summary[order(pausing_group, gene_id_base, bin)],     gene_bin_summary_path, sep = "\t")
  data.table::fwrite(gene_bin_hist[order(pausing_group, bin, modified_call_bin)],   histogram_table_path,  sep = "\t")
  data.table::fwrite(group_meta[order(pausing_group, bin)],                         meta_path,             sep = "\t")

  group_summary <- gene_summary[, .(
    n_genes             = data.table::uniqueN(gene_id_base),
    n_tss               = data.table::uniqueN(tss_uid),
    median_PI           = stats::median(PI, na.rm = TRUE),
    mean_m6A_fraction   = mean(m6A_fraction_modified, na.rm = TRUE),
    median_m6A_fraction = stats::median(m6A_fraction_modified, na.rm = TRUE),
    mean_total_calls    = mean(m6A_total_calls, na.rm = TRUE),
    median_total_calls  = as.numeric(stats::median(m6A_total_calls, na.rm = TRUE))
  ), by = pausing_group]
  data.table::fwrite(group_summary, group_summary_path, sep = "\t")

  meta_plot_path      <- ppfx("m6A_metaprofile_quartiles.pdf")
  boxplot_path        <- ppfx("gene_m6A_fraction_boxplot_quartiles.pdf")
  histogram_plot_path <- ppfx("gene_bin_m6A_histogram_quartiles.pdf")

  plot_m6a_quartile_meta(group_meta,     meta_plot_path,      sample_name, window_bp)
  plot_gene_m6a_boxplot_quartiles(gene_summary, boxplot_path,  sample_name)
  plot_gene_bin_histogram_quartiles(gene_bin_hist, histogram_plot_path, sample_name, window_bp)

  message("Wrote gene summary:        ", gene_summary_path)
  message("Wrote gene-bin summary:    ", gene_bin_summary_path)
  message("Wrote histogram table:     ", histogram_table_path)
  message("Wrote metaprofile table:   ", meta_path)
  message("Wrote group summary:       ", group_summary_path)
  message("Wrote aligned positions:   ", aligned_path)
  message("Wrote metaprofile plot:    ", meta_plot_path)
  message("Wrote boxplot:             ", boxplot_path)
  message("Wrote histogram plot:      ", histogram_plot_path)
} else {
  message("DO_POSITION_SCAN=0 -> skipped m6A position scan and the position-level summary tables/plots.")
}

message("Wrote read accessibility:  ", access_path)
