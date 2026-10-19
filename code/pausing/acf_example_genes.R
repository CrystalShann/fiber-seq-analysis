#' m6A autocorrelation from TSS +50 to +2500 bp for the 15 example genes picked from
#' pause_index.tsv, and whether the three pausing / expression groups differ.
#'
#' Groups (FPKM and PI thirds over expressed genes; genes nearest each group's median
#' PI and FPKM, protein-coding, >= 100 GRO-seq reads in both windows, no other canonical
#' TSS within 2 kb):
#'   non-paused, high   lowest PI third, highest FPKM third
#'   paused, high       highest PI third, highest FPKM third
#'   paused, moderate   highest PI third, middle FPKM third
#'
#' Reads are pooled across the 31 LCL samples, unphased, and must span the whole window
#' (assemble_region_m6a(full_span = TRUE) with every base as a column). The per-read ACF is
#' the parent autocorrelations() on the per-base 0/1 m6A signal; the ACF is invariant to
#' reversing a read, so minus-strand windows are not flipped. Per-gene and per-group mean
#' curves and their 140-250 bp repeat peaks come from the parent summarize(); group curves
#' weight genes equally.
#'
#' Tests, gene = unit, for non-paused high vs paused high (pausing at matched expression)
#' and paused high vs paused moderate (expression at matched pausing):
#'   profile   squared difference of the two group mean curves summed over lags 25-1000,
#'             exact permutation over all 252 splits of the 10 genes
#'   peaks     Wilcoxon rank-sum on per-gene peak height and lag
#'   m6A       Wilcoxon rank-sum on per-gene m6A call fraction (ACF amplitude depends on it)
#'
#' Outputs under <OUT_DIR>/: tables/{genes,acf_per_gene,acf_by_group,tests,sample_contribution}.tsv
#' and plots/{01_group_mean_acf,02_per_gene_acf}.pdf.
#'
#' Run through run_acf_example_genes.sh, which sets the library paths and RETICULATE_PYTHON.

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(reticulate); library(ggplot2)
})

PROJECT   <- "/project/spott/cshan/fiber-seq"
PARENT    <- file.path(PROJECT, "code/clustering_methods/auto_correlation")
FUNCTIONS <- file.path(PROJECT, "code/parsing_functions/parsing_footprints_functions.r")
PLOTTING  <- file.path(PROJECT, "code/parsing_functions/plotting_functions.r")
SAMPLES   <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv"
PI_TSV    <- file.path(PROJECT, "LCL_project/pausing/pause_index.tsv")
OUT_DIR   <- file.path(PROJECT, "LCL_project/pausing/acf_tss_downstream")

GENES <- list(`non-paused, high` = c("FBXO7", "SLC7A1", "COASY", "RSL24D1", "SLC15A3"),
              `paused, high`     = c("POLR2K", "RPS6KB2", "MRPL39", "CCNC", "HEATR1"),
              `paused, moderate` = c("MLH3", "PDLIM5", "MLXIP", "GLMP", "RLIG1"))
WINDOW <- c(50L, 2500L)        # TSS-relative, direction of transcription, inclusive
LAGS   <- 25:1000              # lags for the profile test and the figures
COMPARISONS <- list(c("non-paused, high", "paused, high"), c("paused, high", "paused, moderate"))

stopifnot(all(file.exists(c(FUNCTIONS, PLOTTING, SAMPLES, PI_TSV))))
for (d in c("tables", "plots")) dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)

source(FUNCTIONS)
source(PLOTTING)
source_python(file.path(PARENT, "03_compute_autocorrelations.py"))
source_python(file.path(PARENT, "05_summarize_autocorrelations.py"))
# report_function_code() reads the UTF-8 parent Rmd in the locale encoding; SLURM nodes only
# have the C locale and embedded Python ignores PYTHONUTF8, so run it in a -X utf8 subprocess
report_code <- system2(Sys.getenv("RETICULATE_PYTHON"), c("-X", "utf8", "-c", shQuote(sprintf(
  "import sys; sys.path.insert(0, %s); import importlib; print(importlib.import_module(\"01_run_autocorrelation\").report_function_code())",
  shQuote(PARENT, "sh")))), stdout = TRUE)
stopifnot(is.null(attr(report_code, "status")))
eval(parse(text = report_code))

samples <- fread(SAMPLES)
pi_tab  <- fread(PI_TSV)
genes <- rbindlist(lapply(names(GENES), function(g)
  data.table(group = g, gene_name = GENES[[g]])))
genes[, group_id := match(group, names(GENES))]
genes <- merge(genes, pi_tab[, .(gene_name, chrom, tss, strand, pi, fpkm)], by = "gene_name", sort = FALSE)
stopifnot(nrow(genes) == length(unlist(GENES)))
genes[, analysis_start := ifelse(strand == "+", tss + WINDOW[1], tss - WINDOW[2])]
genes[, analysis_end   := ifelse(strand == "+", tss + WINDOW[2], tss - WINDOW[1])]

# --- per-read ACF, summarised per gene ---------------------------------------
per_gene <- lapply(seq_len(nrow(genes)), function(i) {
  g <- genes[i]
  region <- data.frame(region_id = g$gene_name, chr = g$chrom, start = g$analysis_start - 1L,
                       end = g$analysis_end, analysis_start = g$analysis_start,
                       analysis_end = g$analysis_end)
  positions <- seq.int(g$analysis_start, g$analysis_end)
  dat <- assemble_region_m6a(sample_table = samples, region = region, full_span = TRUE,
                             positions = positions)
  dense <- as.matrix(dat$met_mat)
  storage.mode(dense) <- "integer"
  stopifnot(ncol(dense) == length(positions), !anyNA(dense))
  res   <- autocorrelations(dense)
  valid <- as.logical(res[[2]])
  records <- data.frame(region_id = g$gene_name, RID = dat$rids_df$RID[valid],
                        sample_name = dat$rids_df$sample_name[valid],
                        m6a_call_fraction = rowMeans(dense)[valid], cluster = 1L,
                        status = "clustered")
  s <- summarize(res[[1]][valid, , drop = FALSE], records)
  message(sprintf("  %-8s %-17s reads=%4d valid=%4d m6A=%.4f peak=%s bp",
                  g$gene_name, g$group, nrow(dense), sum(valid), mean(records$m6a_call_fraction),
                  s[[2]]$positive_local_peak_140_250_bp))
  list(curve = as.data.table(s[[1]]), stats = as.data.table(s[[2]]),
       n_reads = nrow(dense), qc = data.table(gene_name = g$gene_name, chr = g$chrom, dat$qc))
})

genes[, n_reads    := vapply(per_gene, `[[`, integer(1), "n_reads")]
genes[, n_acf_valid := vapply(per_gene, function(r) r$stats$n_reads, numeric(1))]
genes[, mean_m6a_call_fraction := vapply(per_gene, function(r) r$stats$mean_m6a_call_fraction, numeric(1))]
genes[, peak_lag_bp := vapply(per_gene, function(r) as.numeric(r$stats$positive_local_peak_140_250_bp), numeric(1))]
genes[, peak_acf    := vapply(per_gene, function(r) r$stats$peak_acf, numeric(1))]
curves <- rbindlist(lapply(seq_len(nrow(genes)), function(i)
  cbind(genes[i, .(group, gene_name)], per_gene[[i]]$curve[, .(lag_bp, mean_acf, sem_acf, n_reads)])))

# a non-BGZF or stale-indexed extract makes tabix return nothing with only a stderr
# message: report any sample/chromosome that gave no reads in every region
sqc <- rbindlist(lapply(per_gene, `[[`, "qc"))
dropout <- sqc[, .(regions = .N, zero = sum(overlapping_reads == 0L)), by = .(sample_name, chr)][regions == zero]
if (nrow(dropout)) message("DROPOUT (0 reads in every region on the chromosome): ",
                           dropout[, paste(sample_name, chr, collapse = "; ")])

# --- group mean curves, genes weighted equally --------------------------------
gene_mat <- do.call(rbind, lapply(per_gene, function(r) r$curve$mean_acf))
gene_records <- data.frame(region_id = "example_genes", RID = genes$gene_name, sample_name = "all",
                           m6a_call_fraction = genes$mean_m6a_call_fraction,
                           cluster = genes$group_id, status = "clustered")
s <- summarize(gene_mat, gene_records)
group_curves <- as.data.table(s[[1]])[, group := names(GENES)[cluster]]
group_peaks  <- as.data.table(s[[2]])[, group := names(GENES)[cluster]]

# --- tests ---------------------------------------------------------------------
lag_cols <- LAGS + 1L
profile_stat <- function(in_a, m) sum((colMeans(m[in_a, , drop = FALSE]) - colMeans(m[!in_a, , drop = FALSE]))^2)
tests <- rbindlist(lapply(COMPARISONS, function(cmp) {
  rows <- which(genes$group %in% cmp)
  m    <- gene_mat[rows, lag_cols, drop = FALSE]
  in_a <- genes$group[rows] == cmp[1]
  obs  <- profile_stat(in_a, m)
  null <- apply(combn(length(rows), sum(in_a)), 2, function(k) profile_stat(seq_along(rows) %in% k, m))
  wil  <- function(col) {
    a <- genes[rows][in_a, get(col)];  b <- genes[rows][!in_a, get(col)]
    a <- a[is.finite(a)];  b <- b[is.finite(b)]
    c(median_a = median(a), median_b = median(b),
      p = if (length(a) >= 2 && length(b) >= 2) wilcox.test(a, b, exact = TRUE)$p.value else NA_real_)
  }
  rbind(data.table(a = cmp[1], b = cmp[2], statistic = sprintf("profile SSD, lags %d-%d", min(LAGS), max(LAGS)),
                   value_a = NA_real_, value_b = NA_real_, observed = obs,
                   p_value = mean(null >= obs - 1e-12), n_permutations = length(null)),
        rbindlist(lapply(c("peak_acf", "peak_lag_bp", "mean_m6a_call_fraction"), function(col) {
          w <- wil(col)
          data.table(a = cmp[1], b = cmp[2], statistic = paste("median", col),
                     value_a = w[["median_a"]], value_b = w[["median_b"]], observed = NA_real_,
                     p_value = w[["p"]], n_permutations = NA_integer_)
        })))
}))

fwrite(genes, file.path(OUT_DIR, "tables", "genes.tsv"), sep = "\t")
fwrite(curves, file.path(OUT_DIR, "tables", "acf_per_gene.tsv"), sep = "\t")
fwrite(group_curves, file.path(OUT_DIR, "tables", "acf_by_group.tsv"), sep = "\t")
fwrite(tests, file.path(OUT_DIR, "tables", "tests.tsv"), sep = "\t")
fwrite(sqc, file.path(OUT_DIR, "tables", "sample_contribution.tsv"), sep = "\t")

# --- figures -------------------------------------------------------------------
region_data <- list(
  region = list(region_id = sprintf("TSS +%d to +%d bp, 5 genes per group, genes weighted equally",
                                    WINDOW[1], WINDOW[2]),
                focal_snp = "groups", ref = "C1 non-paused high, C2 paused high",
                alt = "C3 paused moderate"),
  tables = list(cluster_acf_summary = as.data.frame(group_curves)))
p <- plot_report_cluster_acf(region_data, lower = min(LAGS), upper = max(LAGS),
                             title = "Mean ACF per group; no smoothing")
save_figure(p, file.path(OUT_DIR, "plots", "01_group_mean_acf.pdf"), width = 9, height = 7)

d <- curves[lag_bp %in% LAGS]
d[, group := factor(group, levels = names(GENES))]
m <- group_curves[lag_bp %in% LAGS][, group := factor(group, levels = names(GENES))]
p <- ggplot(d, aes(lag_bp, mean_acf)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = .3) +
  geom_line(aes(colour = gene_name), linewidth = .3) +
  geom_line(data = m, colour = "black", linewidth = .6) +
  facet_wrap(~ group, ncol = 1) +
  labs(x = "Lag (bp)", y = "Mean ACF", colour = "Gene",
       title = sprintf("Per-gene mean m6A ACF, TSS +%d to +%d bp; black = group mean; no smoothing",
                       WINDOW[1], WINDOW[2])) +
  theme_fiberseq("bw")
save_figure(p, file.path(OUT_DIR, "plots", "02_per_gene_acf.pdf"), width = 9, height = 8)

message("\ngroup repeat peaks (group mean curves):")
print(group_peaks[, .(group, genes = n_reads, peak_lag_bp = positive_local_peak_140_250_bp,
                      peak_acf = signif(peak_acf, 3))])
message("\ntests:")
print(tests[, .(a, b, statistic, value_a = signif(value_a, 3), value_b = signif(value_b, 3),
                observed = signif(observed, 3), p_value = signif(p_value, 3))])
message("Wrote ", OUT_DIR)
