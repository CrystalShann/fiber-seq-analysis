#' Select every TSS region for the continuous pausing-index analysis.
#'
#' The continuous counterpart of 01_select_tss_genes.R. It applies stage 01's filters 1-5
#' in the same order with the same constants, and then STOPS: no 5-60 TPM band, no PI
#' tails, no matching. Every gene that survives is selected; PI enters stage 06 as a
#' continuous predictor and expression as a regression covariate instead of a matching
#' variable. `pi_group` and `pair` are kept as NA columns so the schema matches stage 01.
#'
#' Filters, in order (see 01_select_tss_genes.R for the rationale of each):
#'   expressed     both replicates > 0 and mean TPM >= 1
#'   isolated      no other canonical TSS within +-2 kb
#'   long enough   gene_len >= 4 kb
#'   covered       tssr_sum >= 50 and body_sum >= 20 (and a finite PI)
#'   named         a real gene symbol, not an ENSG placeholder
#'
#' The window convention is stage 01's, unchanged: [-1000, +1000) TSS-relative in the
#' direction of transcription, shifted one base on the minus strand so that after the
#' reversal in stage 05 column 1001 is the TSS base on both strands.
#'
#' Input:  <PILOT_DIR>/tables/pause_index.tsv     (stage 00; not rerun here)
#' Output: <OUT_DIR>/tables/tss_regions.tsv        every selected region, stage 01 schema
#'         <OUT_DIR>/tables/filter_funnel.tsv      genes remaining after each filter
#'         <OUT_DIR>/tables/selection_summary.tsv  PI / TPM quantiles, Spearman rho

suppressPackageStartupMessages(library(data.table))

PILOT_DIR <- "/project/spott/cshan/fiber-seq/LCL_project/auto_correlation/lcl_tss"
OUT_DIR   <- file.path(PILOT_DIR, "continuous_pi")
PI_TSV    <- file.path(PILOT_DIR, "tables", "pause_index.tsv")

WINDOW   <- 2000L   # analysis window width (bp)
HALF     <- WINDOW %/% 2L
TSS_COL  <- HALF + 1L   # oriented column holding TSS-relative offset 0
MIN_TPM  <- 1
MIN_LEN  <- 4000L
MIN_TSSR <- 50
MIN_BODY <- 20

stopifnot(file.exists(PI_TSV))
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
pit <- fread(PI_TSV)

n0 <- nrow(pit)
pit <- pit[tpm1 > 0 & tpm2 > 0 & mean_tpm >= MIN_TPM]; n1 <- nrow(pit)
pit <- pit[n_near == 0L];                              n2 <- nrow(pit)
pit <- pit[gene_len >= MIN_LEN];                       n3 <- nrow(pit)
pit <- pit[is.finite(pi) & tssr_sum >= MIN_TSSR & body_sum >= MIN_BODY]; n4 <- nrow(pit)
pit <- pit[grepl("^[A-Z][A-Z0-9-]*$", gene_name)];     n5 <- nrow(pit)

funnel <- data.table(
  step   = 0:5,
  filter = c("pause_index.tsv", "expressed", "isolated", "long_enough", "pi_measurable", "named"),
  rule   = c("protein-coding canonical, chr1-22",
             sprintf("tpm1 > 0 & tpm2 > 0 & mean_tpm >= %g", MIN_TPM),
             "n_near == 0",
             sprintf("gene_len >= %d", MIN_LEN),
             sprintf("is.finite(pi) & tssr_sum >= %g & body_sum >= %g", MIN_TSSR, MIN_BODY),
             "grepl('^[A-Z][A-Z0-9-]*$', gene_name)"),
  n_genes = c(n0, n1, n2, n3, n4, n5))
funnel[, n_removed := c(NA_integer_, -diff(n_genes))]
fwrite(funnel, file.path(OUT_DIR, "tables", "filter_funnel.tsv"), sep = "\t")
print(funnel)
stopifnot(n5 > 1000L)

sel <- copy(pit)
sel[, pi_group := NA_character_]
sel[, pair     := NA_integer_]

# --- region schema (stage 01, unchanged) -------------------------------------
plus <- sel$strand == "+"
sel[, start          := ifelse(plus, tss - HALF,      tss - HALF + 1L)]
sel[, end            := ifelse(plus, tss + HALF,      tss + HALF + 1L)]
sel[, analysis_start := start + 1L]
sel[, analysis_end   := end]
sel[, width          := analysis_end - analysis_start + 1L]
sel[, region_id      := gene_name]
sel[, chr            := chrom]
sel[, tss_col        := TSS_COL]

# --- checks (stage 01's, minus the ones that assume pairs or groups) ---------
stopifnot(nrow(sel) == n5,
          sel[, all(width == WINDOW)],
          sel[, all(analysis_start == start + 1L)],
          sel[, all(analysis_end == end)],
          sel[, all(analysis_start <= tss & tss <= analysis_end)],
          sel[, uniqueN(region_id) == .N],
          sel[, all(pi > 0)])   # log10(pi) is the stage 06 predictor
# offset 0 must land on column TSS_COL on both strands
off0 <- sel[, ifelse(strand == "+", analysis_start + TSS_COL - 1L,
                     analysis_end - TSS_COL + 1L)]
stopifnot(all(off0 == sel$tss + 1L))          # 1-based TSS base

# the pilot genes were drawn from this same filtered pool, so must all be present
pilot <- fread(file.path(PILOT_DIR, "tables", "tss_regions.tsv"))
stopifnot(all(pilot$region_id %in% sel$region_id))
chk <- merge(pilot[, .(region_id, start, end, analysis_start, analysis_end, tss, strand)],
             sel[,   .(region_id, start, end, analysis_start, analysis_end, tss, strand)],
             by = "region_id", suffixes = c(".pilot", ".new"))
stopifnot(chk[, all(start.pilot == start.new & end.pilot == end.new &
                    analysis_start.pilot == analysis_start.new &
                    analysis_end.pilot == analysis_end.new &
                    tss.pilot == tss.new & strand.pilot == strand.new)])
message("all ", nrow(pilot), " pilot regions present with identical coordinates")

# --- selected-set description ------------------------------------------------
probs <- c(0, .05, .25, .5, .75, .95, 1)
rho   <- cor(log10(sel$pi), log10(sel$mean_tpm), method = "spearman")
summ  <- rbind(
  data.table(statistic = "pi_quantile",       prob = probs, value = quantile(sel$pi, probs)),
  data.table(statistic = "mean_tpm_quantile", prob = probs, value = quantile(sel$mean_tpm, probs)),
  data.table(statistic = "spearman_log10pi_log10tpm", prob = NA_real_, value = rho),
  data.table(statistic = "n_genes", prob = NA_real_, value = nrow(sel)))
fwrite(summ, file.path(OUT_DIR, "tables", "selection_summary.tsv"), sep = "\t")

message("selected ", nrow(sel), " genes (", sum(plus), " +, ", sum(!plus), " -)")
message("PI   0/5/25/50/75/95/100%: ",
        paste(signif(quantile(sel$pi, probs), 3), collapse = " / "))
message("TPM  0/5/25/50/75/95/100%: ",
        paste(signif(quantile(sel$mean_tpm, probs), 3), collapse = " / "))
message("Spearman log10(PI) vs log10(TPM): ", round(rho, 3))

cols <- c("region_id", "chr", "start", "end", "analysis_start", "analysis_end", "width",
          "pair", "pi_group", "gene_name", "ensg", "tx_id", "tss", "tss_col", "strand",
          "tes", "gene_len", "tpm1", "tpm2", "mean_tpm", "tssr_sum", "body_sum", "pi")
setorder(sel, chr, tss)
fwrite(sel[, ..cols], file.path(OUT_DIR, "tables", "tss_regions.tsv"), sep = "\t")
message("Wrote ", file.path(OUT_DIR, "tables", "tss_regions.tsv"))
