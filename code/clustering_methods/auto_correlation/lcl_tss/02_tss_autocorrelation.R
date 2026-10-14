#' Per-molecule m6A autocorrelation at the selected TSS windows.
#'
#' Nothing scientific is reimplemented here. The two functions that do the work are
#' sourced unchanged from the parent pipeline through reticulate:
#'
#'   autocorrelations()  ../03_compute_autocorrelations.py
#'                       ACF(k) = sum((x[t]-mean x)(x[t+k]-mean x)) / (N * var x),
#'                       per read, on the unbinned per-base 0/1 signal, every integer
#'                       lag 0..width-1. Zero-variance rows come back all-NaN with a
#'                       false validity flag.
#'   repeat_peak()       ../05_summarize_autocorrelations.py
#'                       strongest positive local maximum of an ACF in 140-250 bp.
#'
#' Read extraction is also reused unchanged: assemble_region_m6a() from
#' code/parsing_functions/parsing_footprints_functions.r, with full_span = TRUE and an
#' explicit `positions` vector. The explicit positions are REQUIRED -- the default
#' collapses the matrix to only those positions carrying a call, which would destroy
#' the base-pair lag axis.
#'
#' Reads are pooled across all 31 LCL samples and are unphased: the pausing groups are
#' a property of the gene, not of a haplotype, so no phasing step is involved.
#'
#' Orientation. Each row is built in genomic order and then reversed for '-' strand
#' genes, so column j always holds TSS-relative offset j - 1001 in the direction of
#' transcription (negative = upstream). The ACF is invariant to this reversal; the
#' orientation matters for the heatmaps and the mean oriented m6A profile.
#'
#' Outputs under <OUT_DIR>/
#'   outputs/<region_id>/acf.rds            per-read ACF matrix + metadata
#'   tables/acf_group_summary.tsv           mean/median ACF per PI group per lag
#'   tables/acf_per_gene.tsv                mean ACF per gene per lag
#'   tables/nrl_per_molecule.tsv            per-molecule repeat_peak (uncalibrated)
#'   tables/group_peaks.tsv                 repeat_peak of each group mean ACF
#'   tables/read_qc.tsv                     per-region read counts and m6A density

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(reticulate)
})

PROJECT  <- "/project/spott/cshan/fiber-seq"
HERE     <- file.path(PROJECT, "code/clustering_methods/auto_correlation/lcl_tss")
PARENT   <- file.path(PROJECT, "code/clustering_methods/auto_correlation")
FUNCTIONS <- file.path(PROJECT, "code/parsing_functions/parsing_footprints_functions.r")
ACF_PY   <- file.path(PARENT, "03_compute_autocorrelations.py")
SUMM_PY  <- file.path(PARENT, "05_summarize_autocorrelations.py")
SAMPLES  <- paste0("/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/",
                   "LCL_sample_metatable_merged_samples_31samples.csv")
OUT_DIR  <- file.path(PROJECT, "LCL_project/auto_correlation/lcl_tss")
REGIONS  <- file.path(OUT_DIR, "tables", "tss_regions.tsv")

MIN_READS <- 100L   # a region with fewer spanning reads is reported and skipped

# Corrected inputs that shadow the shared tree. One shared m6A extract,
# AL6_bc2162_19093 chr1, was re-compressed with plain gzip instead of bgzip, so it is
# not BGZF and cannot be tabix-indexed: htslib returns nothing and reports it on stderr
# only, silently dropping that sample from every chr1 region. The shared tree is
# read-only here, so a lossless BGZF recompression (identical decompressed md5,
# 221,378 records) lives under FIXED_INPUTS and is substituted per sample/chromosome
# through assemble_region_m6a()'s own `m6a_paths` argument. Remove the override once the
# shared file is itself re-compressed.
FIXED_INPUTS <- file.path(PROJECT, "LCL_project/fixed_inputs")
NRL_BAND  <- c(140L, 250L)   # the band repeat_peak() searches, for the record

stopifnot(all(file.exists(c(FUNCTIONS, ACF_PY, SUMM_PY, SAMPLES, REGIONS))))
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

source(FUNCTIONS)
source_python(ACF_PY)
source_python(SUMM_PY)
stopifnot(is.function(autocorrelations), is.function(repeat_peak))
message("sourced autocorrelations() and repeat_peak() from the parent pipeline")

samples <- fread(SAMPLES)
stopifnot(all(c("sample_name", "fire_dir") %in% names(samples)),
          !anyDuplicated(samples$sample_name))
regions <- fread(REGIONS)
message(nrow(regions), " regions; ", nrow(samples), " LCL samples")

acc  <- list()   # per-region results kept for the pooled summaries
qc   <- list()
sqc  <- list()   # per-sample contribution, to catch silent tabix dropouts

for (i in seq_len(nrow(regions))) {
  region <- as.data.frame(regions[i])
  positions <- seq.int(region$analysis_start, region$analysis_end)

  # default paths, then swap in any corrected file that exists for this chromosome
  paths <- file.path(samples$fire_dir, "extracted_results", "m6a_by_chr",
                     paste0(samples$sample_name, ".ft_extracted_m6a.", region$chr, ".bed.gz"))
  fixed <- file.path(FIXED_INPUTS, samples$sample_name, "extracted_results", "m6a_by_chr",
                     paste0(samples$sample_name, ".ft_extracted_m6a.", region$chr, ".bed.gz"))
  use_fixed <- file.exists(fixed) & file.exists(paste0(fixed, ".tbi"))
  if (any(use_fixed)) {
    paths[use_fixed] <- fixed[use_fixed]
    message("  corrected input for ", paste(samples$sample_name[use_fixed], collapse = ", "),
            " ", region$chr)
  }
  dat <- assemble_region_m6a(sample_table = samples, region = region,
                             full_span = TRUE, positions = positions,
                             matrix_dir = NULL, m6a_paths = paths)
  mat <- dat$met_mat
  stopifnot(ncol(mat) == region$width, !anyNA(mat),
            identical(as.integer(colnames(mat)), positions),
            identical(rownames(mat), as.character(dat$rids_df$RID)))

  # orient to the direction of transcription
  oriented_positions <- positions
  if (region$strand == "-") {
    mat <- mat[, rev(seq_len(ncol(mat))), drop = FALSE]
    oriented_positions <- rev(positions)
  }
  # column tss_col must be the TSS base on either strand
  stopifnot(oriented_positions[region$tss_col] == region$tss + 1L)

  dense <- as.matrix(mat)
  storage.mode(dense) <- "integer"
  res   <- autocorrelations(dense)
  acf   <- res[[1]]; valid <- as.logical(res[[2]])
  stopifnot(ncol(acf) == region$width, nrow(acf) == nrow(dense),
            all(abs(acf[valid, 1] - 1) < 1e-10))

  m6a_frac <- Matrix::rowSums(mat) / ncol(mat)
  records <- data.table(
    region_id = region$region_id, pair = region$pair, pi_group = region$pi_group,
    gene_name = region$gene_name, pi = region$pi, mean_tpm = region$mean_tpm,
    strand = region$strand, RID = dat$rids_df$RID,
    sample_name = dat$rids_df$sample_name,
    row_index = seq_len(nrow(dense)) - 1L,
    m6a_count = Matrix::rowSums(mat), m6a_call_fraction = m6a_frac,
    acf_valid = valid)

  saveRDS(list(region = region, records = records, acf = acf,
               oriented_positions = oriented_positions,
               orientation = "transcriptional"),
          { d <- file.path(OUT_DIR, "outputs", region$region_id)
            dir.create(d, recursive = TRUE, showWarnings = FALSE)
            file.path(d, "acf.rds") })

  qc[[i]] <- data.table(region_id = region$region_id, pi_group = region$pi_group,
                        gene_name = region$gene_name, pi = region$pi,
                        mean_tpm = region$mean_tpm,
                        n_reads = nrow(dense), n_acf_valid = sum(valid),
                        n_zero_variance = sum(!valid),
                        mean_m6a_call_fraction = mean(m6a_frac),
                        enough_reads = nrow(dense) >= MIN_READS)
  sqc[[i]] <- data.table(region_id = region$region_id, chr = region$chr,
                         pi_group = region$pi_group, dat$qc)
  acc[[i]] <- list(records = records, acf = acf)
  message(sprintf("  %-8s %-7s n=%4d valid=%4d m6A frac=%.4f",
                  region$region_id, region$pi_group, nrow(dense), sum(valid),
                  mean(m6a_frac)))
}

qc <- rbindlist(qc)
fwrite(qc, file.path(OUT_DIR, "tables", "read_qc.tsv"), sep = "\t")
if (any(!qc$enough_reads))
  warning("regions below ", MIN_READS, " reads: ",
          paste(qc[!enough_reads]$region_id, collapse = ", "))

# Per-sample contributions. A sample that returns zero overlapping reads for every
# region on a chromosome is not biology: it is a stale tabix index whose offsets point
# past the end of its .bed.gz, and htslib reports it only on stderr. Known instance:
# AL6_bc2162_19093 chr1, whose index was built against a larger file than the one on
# disk, so that sample silently contributes nothing to any chr1 region.
sqc <- rbindlist(sqc)
fwrite(sqc, file.path(OUT_DIR, "tables", "sample_contribution.tsv"), sep = "\t")
dropout <- sqc[, .(regions = .N, zero = sum(overlapping_reads == 0L)),
               by = .(sample_name, chr)][regions == zero]
if (nrow(dropout)) {
  msg <- dropout[, paste0(sample_name, " ", chr, " (0 reads in all ", regions,
                          " region(s))", collapse = "; ")]
  warning("sample/chromosome dropouts, likely stale tabix indexes: ", msg)
  message("DROPOUT: ", msg)
} else message("no sample/chromosome dropouts")

records <- rbindlist(lapply(acc, `[[`, "records"))
acf_all <- do.call(rbind, lapply(acc, `[[`, "acf"))
stopifnot(nrow(records) == nrow(acf_all))
lags <- seq_len(ncol(acf_all)) - 1L

#' Mean, median and SEM across a set of rows, ignoring invalid ones.
curve_of <- function(rows, label_cols) {
  m <- acf_all[rows & records$acf_valid, , drop = FALSE]
  if (!nrow(m)) return(NULL)
  data.table(label_cols, lag_bp = lags,
             mean_acf = colMeans(m), median_acf = apply(m, 2, median),
             sem_acf = apply(m, 2, sd) / sqrt(nrow(m)), n_reads = nrow(m))
}

group_curves <- rbindlist(lapply(sort(unique(records$pi_group)), function(g)
  curve_of(records$pi_group == g, data.table(pi_group = g))))
fwrite(group_curves, file.path(OUT_DIR, "tables", "acf_group_summary.tsv"), sep = "\t")

gene_curves <- rbindlist(lapply(unique(records$region_id), function(r) {
  rr <- records$region_id == r
  curve_of(rr, data.table(region_id = r, pi_group = records$pi_group[rr][1],
                          pair = records$pair[rr][1]))
}))
fwrite(gene_curves, file.path(OUT_DIR, "tables", "acf_per_gene.tsv"), sep = "\t")

# --- repeat_peak() on the group and per-gene mean curves ---------------------
peak_of <- function(curve) {
  p <- repeat_peak(np_array(curve))
  list(lag = if (is.null(p[[1]])) NA_integer_ else as.integer(p[[1]]),
       height = if (is.null(p[[2]])) NA_real_ else as.numeric(p[[2]]))
}
group_peaks <- rbindlist(lapply(split(group_curves, by = "pi_group"), function(d) {
  p <- peak_of(d$mean_acf)
  data.table(level = "pi_group", label = d$pi_group[1], n_reads = d$n_reads[1],
             peak_lag_bp = p$lag, peak_height = p$height)
}))
gene_peaks <- rbindlist(lapply(split(gene_curves, by = "region_id"), function(d) {
  p <- peak_of(d$mean_acf)
  data.table(level = "gene", label = d$region_id[1], n_reads = d$n_reads[1],
             peak_lag_bp = p$lag, peak_height = p$height,
             pi_group = d$pi_group[1], pair = d$pair[1])
}))
fwrite(rbind(group_peaks, gene_peaks, fill = TRUE),
       file.path(OUT_DIR, "tables", "group_peaks.tsv"), sep = "\t")

# --- the primary statistic: paired, per gene --------------------------------
# Molecules within a gene share one promoter and are NOT independent, so a
# molecule-level test is anticonservative by thousands of degrees of freedom. The
# experimental unit is the gene, and the design is paired on expression, so the test
# is a paired comparison over the matched pairs. This also protects against the
# pooled curve being dominated by whichever genes happen to contribute most molecules
# (283-669 here).
gp  <- merge(gene_peaks[, .(region_id = label, n_reads, peak_lag_bp, peak_height)],
             regions[, .(region_id, pair, pi_group, gene_name, pi, mean_tpm)],
             by = "region_id")
wide <- dcast(gp, pair ~ pi_group,
              value.var = c("gene_name", "n_reads", "peak_lag_bp", "peak_height", "pi"))
wide[, `:=`(height_diff = peak_height_high_pi - peak_height_low_pi,
            lag_diff    = peak_lag_bp_high_pi - peak_lag_bp_low_pi)]
fwrite(wide, file.path(OUT_DIR, "tables", "pair_comparison.tsv"), sep = "\t")

paired_test <- function(a, b, what) {
  keep <- is.finite(a) & is.finite(b)
  if (sum(keep) < 3L)
    return(data.table(statistic = what, n_pairs = sum(keep), n_high_greater = NA_integer_,
                      median_diff = NA_real_, wilcoxon_p = NA_real_, t_p = NA_real_))
  data.table(statistic = what, n_pairs = sum(keep),
             n_high_greater = sum(a[keep] > b[keep]),
             median_diff = median(a[keep] - b[keep]),
             wilcoxon_p = wilcox.test(a[keep], b[keep], paired = TRUE)$p.value,
             t_p = t.test(a[keep], b[keep], paired = TRUE)$p.value)
}
paired <- rbind(
  paired_test(wide$peak_height_high_pi, wide$peak_height_low_pi, "acf_peak_height"),
  paired_test(wide$peak_lag_bp_high_pi, wide$peak_lag_bp_low_pi, "acf_peak_lag_bp"))
fwrite(paired, file.path(OUT_DIR, "tables", "paired_tests.tsv"), sep = "\t")

# --- per-molecule NRL: secondary, explicitly uncalibrated --------------------
# repeat_peak() applies no prominence constraint, so on a single molecule it will
# return a noise bump inside the band (per-base single-molecule ACF noise is about
# 1/sqrt(width)). Reported for completeness; the group curves are the primary readout.
nrl <- vapply(seq_len(nrow(acf_all)), function(k) {
  if (!records$acf_valid[k]) return(c(NA_real_, NA_real_))
  p <- peak_of(acf_all[k, ])
  c(as.numeric(p$lag), p$height)
}, numeric(2))
per_mol <- cbind(records, nrl_bp = nrl[1, ], nrl_height = nrl[2, ])
per_mol[, nrl_band_lo := NRL_BAND[1]][, nrl_band_hi := NRL_BAND[2]]
per_mol[, nrl_calibrated := FALSE]
fwrite(per_mol, file.path(OUT_DIR, "tables", "nrl_per_molecule.tsv"), sep = "\t")

message("\nmolecules: ", nrow(records), " (",
        paste(sprintf("%s %d", names(table(records$pi_group)), table(records$pi_group)),
              collapse = ", "), ")")
message("\npooled group curves (molecule-weighted, NOT the inferential statistic):")
print(group_peaks)
message("\nprimary statistic -- paired over the ", nrow(wide), " matched gene pairs:")
print(paired)
if (all(is.na(paired$wilcoxon_p)) || min(paired$wilcoxon_p, na.rm = TRUE) >= 0.05)
  message("\nNo significant paired difference in phasing amplitude or NRL at this ",
          "pilot size. Any gap between the pooled group curves reflects between-gene ",
          "heterogeneity and unequal molecule counts, not pausing.")
message("Wrote tables to ", file.path(OUT_DIR, "tables"))
