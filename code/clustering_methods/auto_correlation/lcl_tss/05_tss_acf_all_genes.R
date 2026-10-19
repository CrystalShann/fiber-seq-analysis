#' Per-gene m6A autocorrelation at every TSS region selected by 04_select_all_genes.R.

#' Usage (through run_continuous_pi.sh, which sets the library environment):
#'   Rscript 05_tss_acf_all_genes.R            all regions, then the tables
#'   Rscript 05_tss_acf_all_genes.R validate   the 12 pilot genes vs the pilot tables
#' Set LCL_TSS_KEEP_MOLECULES="GENE1,GENE2" to also save the per-molecule ACF of those
#' genes, in stage 02's acf.rds layout.
#'
#' Outputs under <OUT_DIR>/
#'   outputs/acf_all_genes.rds              every gene: scalars, curves, m6A profile, qc
#'                                          (or its error)
#'   outputs/molecules/<region_id>/acf.rds  per-molecule ACF, only for KEEP_MOLECULES
#'   tables/gene_acf_summary.tsv            one row per gene, scalar fields
#'   tables/acf_per_gene.tsv.gz             long format mean/median/SEM ACF
#'   tables/sample_contribution.tsv         assemble_region_m6a() per-sample qc
#'   tables/failed_regions.tsv              regions that failed, and why
#'   tables/pilot_validation.tsv            (validate) the 12 pilot genes, new vs pilot

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(reticulate)
})

PROJECT  <- "/project/spott/cshan/fiber-seq"
PARENT   <- file.path(PROJECT, "code/clustering_methods/auto_correlation")
FUNCTIONS <- file.path(PROJECT, "code/parsing_functions/parsing_footprints_functions.r")
ACF_PY   <- file.path(PARENT, "03_compute_autocorrelations.py")
SUMM_PY  <- file.path(PARENT, "05_summarize_autocorrelations.py")
SAMPLES  <- paste0("/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/",
                   "LCL_sample_metatable_merged_samples_31samples.csv")
PILOT_DIR <- file.path(PROJECT, "LCL_project/auto_correlation/lcl_tss")
OUT_DIR  <- file.path(PILOT_DIR, "continuous_pi")
REGIONS  <- file.path(OUT_DIR, "tables", "tss_regions.tsv")
RDS_OUT  <- file.path(OUT_DIR, "outputs", "acf_all_genes.rds")
CHECKPOINT_EVERY <- 500L   # genes between rewrites of RDS_OUT

# Corrected inputs that shadow the shared tree (see 02_tss_autocorrelation.R): the
# shared AL6_bc2162_19093 chr1 m6A extract is plain gzip, not BGZF, so tabix returns
# nothing for it. A lossless BGZF recompression lives under FIXED_INPUTS.
FIXED_INPUTS <- file.path(PROJECT, "LCL_project/fixed_inputs")
KEEP_MOLECULES <- strsplit(Sys.getenv("LCL_TSS_KEEP_MOLECULES", ""), ",")[[1]]
META_COLS <- c("region_id", "gene_name", "ensg", "chr", "strand", "tss", "pi",
               "mean_tpm", "tpm1", "tpm2", "gene_len")

stopifnot(all(file.exists(c(FUNCTIONS, ACF_PY, SUMM_PY, SAMPLES, REGIONS))))
args <- commandArgs(trailingOnly = TRUE)
regions <- fread(REGIONS)
stopifnot(!anyDuplicated(regions$region_id))

if (identical(args, "validate")) {
  # ------------------------------------------- validation against the pilot ---
  genes  <- fread(file.path(OUT_DIR, "tables", "gene_acf_summary.tsv"))
  curves <- fread(file.path(OUT_DIR, "tables", "acf_per_gene.tsv.gz"))
  pq <- fread(file.path(PILOT_DIR, "tables", "read_qc.tsv"))
  pp <- fread(file.path(PILOT_DIR, "tables", "group_peaks.tsv"))[level == "gene"]
  pc <- fread(file.path(PILOT_DIR, "tables", "acf_per_gene.tsv"))
  stopifnot(all(pq$region_id %in% genes$region_id))
  v <- merge(genes[region_id %in% pq$region_id,
                   .(region_id, n_reads, n_acf_valid, peak_lag_bp, peak_height)],
             merge(pq[, .(region_id, qc_n_reads = n_reads, qc_n_acf_valid = n_acf_valid)],
                   pp[, .(region_id = label, pk_n_reads = n_reads,
                          pilot_peak_lag_bp = peak_lag_bp, pilot_peak_height = peak_height)],
                   by = "region_id"), by = "region_id")
  cc <- merge(curves[region_id %in% pq$region_id, .(region_id, lag_bp, mean_acf)],
              pc[, .(region_id, lag_bp, pilot_mean_acf = mean_acf)], by = c("region_id", "lag_bp"))
  v <- merge(v, cc[, .(n_lags = .N, max_abs_diff_mean_acf = max(abs(mean_acf - pilot_mean_acf))),
                   by = region_id], by = "region_id")
  v[, height_abs_diff := abs(peak_height - pilot_peak_height)]
  fwrite(v, file.path(OUT_DIR, "tables", "pilot_validation.tsv"), sep = "\t")
  print(v[, .(region_id, n_reads, qc_n_reads, n_acf_valid, pk_n_reads, peak_lag_bp,
              pilot_peak_lag_bp, height_abs_diff, max_abs_diff_mean_acf)])
  ok <- c(
    n_reads      = v[, all(n_reads == qc_n_reads)],
    n_acf_valid  = v[, all(n_acf_valid == qc_n_acf_valid & n_acf_valid == pk_n_reads)],
    peak_lag_bp  = v[, identical(peak_lag_bp, as.integer(pilot_peak_lag_bp))],
    peak_height  = v[, identical(is.na(peak_height), is.na(pilot_peak_height)) &&
                       all(height_abs_diff <= 1e-15, na.rm = TRUE)],
    mean_acf     = v[, all(n_lags == 2000L) && all(max_abs_diff_mean_acf < 1e-10)])
  print(ok)
  if (!all(ok)) stop("stage 05 does not reproduce the pilot for: ",
                     paste(names(ok)[!ok], collapse = ", "))
  message("pilot validation passed: ", nrow(v), " genes; max |peak height diff| = ",
          signif(max(v$height_abs_diff, na.rm = TRUE), 3),
          ", max |mean ACF diff| = ", signif(max(v$max_abs_diff_mean_acf), 3))
  quit(save = "no")
}

# ------------------------------------------------------------ all regions ---
stopifnot(length(args) == 0L)
dir.create(dirname(RDS_OUT), recursive = TRUE, showWarnings = FALSE)

# resume from the checkpoint, which must come from this same region table
regions_md5 <- unname(tools::md5sum(REGIONS))
results <- list()
if (file.exists(RDS_OUT)) {
  saved <- readRDS(RDS_OUT)
  if (!identical(saved$regions_md5, regions_md5))
    stop(RDS_OUT, " was built from a different tss_regions.tsv; delete it to start over")
  results <- saved$genes
}
done <- names(results)[vapply(results, function(r) is.null(r$error), logical(1))]
todo <- which(!regions$region_id %in% done)
message(nrow(regions) - length(todo), " of ", nrow(regions), " regions already done, ",
        length(todo), " to run")

if (length(todo)) {
  source(FUNCTIONS)
  source_python(ACF_PY)
  source_python(SUMM_PY)
  stopifnot(is.function(autocorrelations), is.function(summarize))
  samples <- fread(SAMPLES)
  stopifnot(all(c("sample_name", "fire_dir") %in% names(samples)),
            !anyDuplicated(samples$sample_name))
}

for (k in seq_along(todo)) {
  region <- as.data.frame(regions[todo[k]])
  t0 <- Sys.time()
  results[[region$region_id]] <- tryCatch({
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
    dat <-assemble_region_m6a(sample_table = samples, region = region,
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
    records <- data.frame(region_id = region$region_id, RID = dat$rids_df$RID,
                          sample_name = dat$rids_df$sample_name,
                          m6a_call_fraction = m6a_frac, cluster = 0L,
                          status = ifelse(valid, "clustered", "zero_variance"))

    # parent summarize(): all valid reads of the gene as one cluster
    lag_bp <- seq_len(ncol(acf)) - 1L
    curves <- data.table(lag_bp = lag_bp, mean_acf = NA_real_, median_acf = NA_real_,
                         sem_acf = NA_real_)
    peak <- list(lag = NA_integer_, height = NA_real_)
    if (any(valid)) {
      s <- summarize(acf, records)
      avg <- s[[1]]; st <- s[[2]]
      stopifnot(nrow(st) == 1L, st$n_reads == sum(valid), all(avg$lag_bp == lag_bp))
      curves[, `:=`(mean_acf = avg$mean_acf,
                    median_acf = apply(acf[valid, , drop = FALSE], 2, median),
                    sem_acf = avg$sem_acf)]
      if (is.finite(st$positive_local_peak_140_250_bp))
        peak <- list(lag = as.integer(st$positive_local_peak_140_250_bp),
                     height = as.numeric(st$peak_acf))
    }

    gene <- data.table(region[META_COLS],
                       n_reads = nrow(dense), n_acf_valid = sum(valid),
                       n_zero_variance = sum(!valid),
                       mean_m6a_call_fraction = mean(m6a_frac),
                       peak_lag_bp = peak$lag, peak_height = peak$height)
    # oriented: column j is TSS-relative offset j - tss_col
    m6a_profile <- data.table(offset_bp = seq_len(ncol(mat)) - region$tss_col,
                              mean_m6a = Matrix::colMeans(mat))
    qc <- data.table(region_id = region$region_id, chr = region$chr, dat$qc)

    if (region$region_id %in% KEEP_MOLECULES) {
      d <- file.path(OUT_DIR, "outputs", "molecules", region$region_id)
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      saveRDS(list(region = region, records = records, acf = acf,
                   oriented_positions = oriented_positions,
                   orientation = "transcriptional"), file.path(d, "acf.rds"))
    }
    message(sprintf("  %-10s n=%4d valid=%4d peak=%s  %.1fs", region$region_id,
                    nrow(dense), sum(valid), peak$lag,
                    as.numeric(Sys.time() - t0, units = "secs")))
    list(gene = gene, curves = curves, m6a_profile = m6a_profile, qc = qc)
  }, error = function(e) {
    message("  FAILED ", region$region_id, ": ", conditionMessage(e))
    list(error = conditionMessage(e))
  })

  if (k %% CHECKPOINT_EVERY == 0L || k == length(todo)) {
    # atomic: a job killed mid-save must not leave a truncated checkpoint
    tmp <- paste0(RDS_OUT, ".tmp.", Sys.getpid())
    saveRDS(list(regions_md5 = regions_md5, genes = results), tmp)
    stopifnot(file.rename(tmp, RDS_OUT))
    message("checkpoint: ", k, " of ", length(todo), " run; ", length(results),
            " genes in ", RDS_OUT)
  }
}

# ------------------------------------------------------------------ tables ---
results <- results[intersect(regions$region_id, names(results))]
ok  <- vapply(results, function(r) is.null(r$error), logical(1))
genes <- rbindlist(lapply(results[ok], `[[`, "gene"))
fwrite(genes, file.path(OUT_DIR, "tables", "gene_acf_summary.tsv"), sep = "\t")
curves <- rbindlist(lapply(results[ok], function(r)
  data.table(region_id = r$gene$region_id, r$curves, n_reads = r$gene$n_acf_valid)))
fwrite(curves, file.path(OUT_DIR, "tables", "acf_per_gene.tsv.gz"), sep = "\t")
failed <- data.table(region_id = names(results)[!ok],
                     error = vapply(results[!ok], `[[`, character(1), "error"))
fwrite(failed, file.path(OUT_DIR, "tables", "failed_regions.tsv"), sep = "\t")
message(sum(ok), " / ", nrow(regions), " regions done; ", nrow(failed), " failed")
if (nrow(failed)) print(failed[, .N, by = error][order(-N)])

# Per-sample contributions, and stage 02's dropout detector: a sample that returns zero
# overlapping reads for every region on a chromosome is an unreadable input, not biology.
sqc <- rbindlist(lapply(results[ok], `[[`, "qc"))
fwrite(sqc, file.path(OUT_DIR, "tables", "sample_contribution.tsv"), sep = "\t")
dropout <- sqc[, .(regions = .N, zero = sum(overlapping_reads == 0L)),
               by = .(sample_name, chr)][regions == zero]
if (nrow(dropout)) {
  msg <- dropout[, paste0(sample_name, " ", chr, " (0 reads in all ", regions,
                          " region(s))", collapse = "; ")]
  warning("sample/chromosome dropouts, likely stale tabix indexes: ", msg)
  message("DROPOUT: ", msg)
} else message("no sample/chromosome dropouts")
message("Wrote ", RDS_OUT, " and tables to ", file.path(OUT_DIR, "tables"))
