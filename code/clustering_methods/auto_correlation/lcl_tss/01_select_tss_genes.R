#' Select the TSS regions for the high- vs low-pausing autocorrelation comparison.
#'
#' Reads the pausing-index table from 00_pause_index.R and picks example genes from the
#' two tails of PI, 1:1 matched on expression. The matching is the point: PI is weakly
#' ANTI-correlated with expression in this data (Spearman log10 PI vs log10 TPM about
#' -0.27), so unmatched tails would differ in transcription rate as well as pausing and
#' the contrast would be uninterpretable.
#'
#' Filters, in order:
#'   expressed            both replicates > 0 and mean TPM >= 1   (PI is undefined for
#'                        silent genes: both windows are empty, giving 0/Inf/NA)
#'   isolated             no other canonical TSS within +-2 kb, so the 2 kb analysis
#'                        window is not contaminated by a neighbouring promoter
#'   long enough          gene_len >= 4 kb, so the body window is a real measurement
#'   covered              tssr_sum >= 50 and body_sum >= 20
#'   named                a real gene symbol, not an ENSG placeholder
#'   TPM band 5-60        pilot genes only. Below ~5 TPM both PI windows are sparse and
#'                        PI is dominated by counting noise; above ~60 TPM the pool is
#'                        mostly ribosomal-protein and translation genes with their own
#'                        promoter architecture. Without this band the greedy
#'                        strongest-PI-first matching drifts to 1-4 TPM genes.
#'
#' Output: <OUT_DIR>/tables/tss_regions.tsv, in exactly the schema
#' assemble_region_m6a() requires (region_id, chr, start, end, analysis_start,
#' analysis_end, width) plus the gene and group annotation.
#'
#' Window convention. Offsets are TSS-relative and transcriptional, [-1000, +1000),
#' so after the minus-strand reversal in 02 every matrix column j is offset j - 1001
#' and column 1001 is the TSS base. That requires the genomic window to be shifted by
#' one base on the minus strand:
#'   + strand: 0-based [tss - 1000, tss +  999]  ->  1-based [tss - 999, tss + 1000]
#'   - strand: 0-based [tss -  999, tss + 1000]  ->  1-based [tss - 998, tss + 1001]
#' both 2000 bp wide, both with offset 0 landing on column 1001.

suppressPackageStartupMessages(library(data.table))

OUT_DIR <- "/project/spott/cshan/fiber-seq/LCL_project/auto_correlation/lcl_tss"
PI_TSV  <- file.path(OUT_DIR, "tables", "pause_index.tsv")

WINDOW      <- 2000L   # analysis window width (bp)
HALF        <- WINDOW %/% 2L
TSS_COL     <- HALF + 1L   # oriented column holding TSS-relative offset 0
N_PAIRS     <- 6L
PI_TAIL     <- 0.02   # top/bottom fraction of PI to draw from
TPM_TOL_DEX <- 0.15   # maximum |log10 TPM| difference within a matched pair
MIN_TPM     <- 1
TPM_BAND    <- c(5, 60)   # pilot genes must sit in this TPM band (see note below)
MIN_LEN     <- 4000L
MIN_TSSR    <- 50
MIN_BODY    <- 20
SEED        <- 0L

stopifnot(file.exists(PI_TSV))
set.seed(SEED)
pit <- fread(PI_TSV)

n0 <- nrow(pit)
pit <- pit[tpm1 > 0 & tpm2 > 0 & mean_tpm >= MIN_TPM]; n1 <- nrow(pit)
pit <- pit[n_near == 0L];                              n2 <- nrow(pit)
pit <- pit[gene_len >= MIN_LEN];                       n3 <- nrow(pit)
pit <- pit[is.finite(pi) & tssr_sum >= MIN_TSSR & body_sum >= MIN_BODY]; n4 <- nrow(pit)
pit <- pit[grepl("^[A-Z][A-Z0-9-]*$", gene_name)];     n5 <- nrow(pit)
pit <- pit[mean_tpm >= TPM_BAND[1] & mean_tpm <= TPM_BAND[2]]; n6 <- nrow(pit)
message(sprintf(paste("%d genes -> expressed %d -> isolated %d -> >=%dkb %d ->",
                      "covered %d -> named %d -> TPM %g-%g %d"),
                n0, n1, n2, MIN_LEN %/% 1000L, n3, n4, n5, TPM_BAND[1], TPM_BAND[2], n6))
stopifnot(n6 > 200L)

pi_hi_cut <- quantile(pit$pi, 1 - PI_TAIL)
pi_lo_cut <- quantile(pit$pi,     PI_TAIL)
hi <- pit[pi >= pi_hi_cut][order(-pi)]
lo <- pit[pi <= pi_lo_cut][order(pi)]
message("PI tails: high n=", nrow(hi), " (", round(min(hi$pi), 1), "-", round(max(hi$pi), 1),
        "), low n=", nrow(lo), " (", signif(min(lo$pi), 2), "-", signif(max(lo$pi), 2), ")")

# greedy 1:1 match on log10 expression, strongest-PI genes first
used  <- character(0)
pairs <- list()
for (i in seq_len(nrow(hi))) {
  cand <- lo[!gene_name %in% used]
  if (!nrow(cand)) break
  d <- abs(log10(cand$mean_tpm) - log10(hi$mean_tpm[i]))
  j <- which.min(d)
  if (d[j] > TPM_TOL_DEX) next
  used <- c(used, cand$gene_name[j])
  pairs[[length(pairs) + 1L]] <- rbind(
    cbind(hi[i],      pi_group = "high_pi", pair = length(pairs) + 1L),
    cbind(cand[j],    pi_group = "low_pi",  pair = length(pairs) + 1L))
  if (length(pairs) == N_PAIRS) break
}
stopifnot(length(pairs) == N_PAIRS)
sel <- rbindlist(pairs)

# --- region schema -----------------------------------------------------------
plus <- sel$strand == "+"
sel[, start          := ifelse(plus, tss - HALF,      tss - HALF + 1L)]
sel[, end            := ifelse(plus, tss + HALF,      tss + HALF + 1L)]
sel[, analysis_start := start + 1L]
sel[, analysis_end   := end]
sel[, width          := analysis_end - analysis_start + 1L]
sel[, region_id      := gene_name]
sel[, chr            := chrom]
sel[, tss_col        := TSS_COL]

# --- checks ------------------------------------------------------------------
stopifnot(nrow(sel) == 2L * N_PAIRS,
          sel[, all(width == WINDOW)],
          sel[, all(analysis_start == start + 1L)],
          sel[, all(analysis_end == end)],
          sel[, all(analysis_start <= tss & tss <= analysis_end)],
          sel[, uniqueN(region_id) == .N],
          sel[, .N, by = pi_group][, all(N == N_PAIRS)])
# offset 0 must land on column TSS_COL on both strands
off0 <- sel[, ifelse(strand == "+", analysis_start + TSS_COL - 1L,
                     analysis_end - TSS_COL + 1L)]
stopifnot(all(off0 == sel$tss + 1L))          # 1-based TSS base

w <- wilcox.test(mean_tpm ~ pi_group, data = sel)$p.value
message("matched pairs: ", N_PAIRS,
        " | TPM median high=", round(median(sel[pi_group == "high_pi"]$mean_tpm), 1),
        " low=", round(median(sel[pi_group == "low_pi"]$mean_tpm), 1),
        " | Wilcoxon p=", round(w, 3))
if (w < 0.05) warning("expression is not matched between the PI groups (p < 0.05)")

cols <- c("region_id", "chr", "start", "end", "analysis_start", "analysis_end", "width",
          "pair", "pi_group", "gene_name", "ensg", "tx_id", "tss", "tss_col", "strand",
          "tes", "gene_len", "tpm1", "tpm2", "mean_tpm", "tssr_sum", "body_sum", "pi")
setcolorder(sel, cols)
fwrite(sel[, ..cols], file.path(OUT_DIR, "tables", "tss_regions.tsv"), sep = "\t")
print(sel[order(pair, pi_group), .(pair, pi_group, gene_name, chr, tss, strand,
                                   mean_tpm = round(mean_tpm, 1), pi = round(pi, 2))])
message("Wrote ", file.path(OUT_DIR, "tables", "tss_regions.tsv"))
