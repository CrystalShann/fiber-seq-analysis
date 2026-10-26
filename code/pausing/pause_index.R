#' RNA Pol II pausing index for GM12878 from GRO-seq, in hg38.
#'
#' Method (defined on Pol II ChIP-seq in the original, applied to GRO-seq here):
#'   TSS region (TSSR) = -50 .. +300 around the TSS
#'   gene body         = +300 past the TSS .. 3 kb past the TES
#'   density           = reads per million mapped reads per bp (rpm/bp) in each window
#'   PI                = TSSR density / gene-body density
#'   active TSS        = TSSR density >= 0.001 rpm/bp; genes below it get no PI and are
#'                       counted as non-paused
#'   paused            = PI > 2
#'   expressed         = FPKM > 1, FPKM = mean of the two RNA-seq replicates
#'
#' Genome builds:
#'   GRO-seq GSM1480326 (Core et al. 2014)       hg19 -> signal lifted to hg38 here
#'   TSS, GENCODE v46 Ensembl canonical          hg38
#'   TES, same canonical transcript, v46 GTF     hg38
#'   RNA-seq ENCSR820PHH gene_quan (RSEM)        GRCh38, GENCODE V29 (joined on version-stripped ENSG)
#'
#' The bigWigs hold one count per read at a single base, so the library size is the summed
#' signal over every hg19 chromosome on both strands.
#'
#' Output: <OUT_DIR>/pause_index.tsv, one row per canonical gene on chr1-22 and chrX.
#'         All coordinates are 1-based inclusive (tss = BED start + 11).
#'
#' Also writes dominant_tss/pause_index_dominant.tsv using CAGE-selected anchors.
#' Run through run_pause_index.sh, which sets the libcurl path rtracklayer needs.

suppressPackageStartupMessages({
  library(data.table); library(rtracklayer); library(GenomicRanges)
})

ANN     <- "/project/spott/cshan/annotations"
LCL_PUB <- file.path(ANN, "LCL_public_data")
TSS_BED <- file.path(ANN, "TSS_interval_gencodev46_Ensembl_canonical.bed")
GTF     <- file.path(ANN, "gencode.v46.annotation.gtf.gz")
CHAIN   <- file.path(ANN, "hg19ToHg38.over.chain")
GRO     <- c(`+` = file.path(LCL_PUB, "GRO-seq/GSM1480326_GM12878_GROseq_plus.bigWig"),
             `-` = file.path(LCL_PUB, "GRO-seq/GSM1480326_GM12878_GROseq_minus.bigWig"))
EXPR    <- file.path(LCL_PUB, "gene_quan_ENCSR820PHH",
                     c("ENCFF873VWU.gene_quan.tsv", "ENCFF345SHY.gene_quan.tsv"))  # rep 1, rep 2
OUT_DIR <- "/project/spott/cshan/fiber-seq/LCL_project/pausing"

# Both canonical and dominant anchors are quantified in the same GRO-seq pass.
opts <- list(cage = file.path(ANN, "fantom5_LCL/fantom5.hg38.LCL.consensus.TSS.bed.gz"),
             `tss-distance` = "250", `out-dir` = OUT_DIR)
argv <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% argv) {
  cat("Usage: run_pause_index.sh [--cage hg38_CTSS.bed.gz] [--tss-distance 250] [--out-dir DIR]\n")
  quit(status = 0L)
}
if (length(argv) %% 2L) stop("Expected --option value pairs")
for (i in seq_len(length(argv) %/% 2L) * 2L - 1L) {
  key <- sub("^--", "", argv[i])
  if (!startsWith(argv[i], "--") || !key %in% names(opts)) stop("Unknown option: ", argv[i])
  opts[[key]] <- argv[i + 1L]
}
OUT_DIR <- opts[["out-dir"]]
tss_distance <- as.numeric(opts[["tss-distance"]])
stopifnot(length(tss_distance) == 1L, is.finite(tss_distance),
          tss_distance >= 0, tss_distance %% 1 == 0, file.exists(opts$cage))
script_path <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1L])
source(file.path(dirname(normalizePath(script_path)), "dominant_tss/select_dominant_tss.R"))

CHROMS    <- paste0("chr", c(1:22, "X"))   # GRO-seq has no chrY
TSSR      <- c(-50L, 300L)
BODY_FROM <- 300L
BODY_PAST <- 3000L
MIN_DENS  <- 0.001                         # rpm/bp in the TSSR
PAUSED_PI <- 2
EXPR_FPKM <- 1

stopifnot(all(file.exists(c(TSS_BED, GTF, CHAIN, GRO, EXPR))))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- canonical TSS -----------------------------------------------------------
bed <- fread(TSS_BED, header = FALSE,
             col.names = c("chrom", "start", "end", "name", "score", "strand"))
stopifnot(all(bed$end - bed$start == 20L))          # [TSS-10, TSS+10), 0-based
bed[, c("gene_id", "tx_id", "gene_name", "gene_type") := tstrsplit(name, ";", keep = 1:4)]
genes <- bed[chrom %in% CHROMS,
             .(ensg = sub("\\..*", "", gene_id), gene_id, tx_id, gene_name, gene_type,
               chrom, tss = start + 11L, strand)]
message(nrow(genes), " canonical TSS on chr1-22, chrX")

# --- TES of the same canonical transcript ------------------------------------
tx <- fread(cmd = paste("zcat", GTF, "| awk -F'\t' '$3 == \"transcript\"' | cut -f1,4,5,7,9"),
            header = FALSE, col.names = c("chrom", "start", "end", "strand", "attr"),
            quote = "")
tx[, tx_id := sub('.*transcript_id "([^"]+)".*', "\\1", attr)]
all_tx <- copy(tx)
all_tx[, `:=`(gene_id = sub('.*gene_id "([^"]+)".*', "\\1", attr),
              gene_name = sub('.*gene_name "([^"]+)".*', "\\1", attr),
              gene_type = sub('.*gene_type "([^"]+)".*', "\\1", attr),
              canonical = grepl('tag "Ensembl_canonical"', attr, fixed = TRUE))]
all_tx[, `:=`(tss = fifelse(strand == "+", start, end),
              tes = fifelse(strand == "+", end, start))]
all_tx <- all_tx[chrom %in% CHROMS]
con <- gzfile(opts$cage, "rt")
cage <- fread(text = paste(readLines(con), collapse = "\n"), header = FALSE, select = 1:6,
              col.names = c("chrom", "start0", "end", "ctss_id", "score", "strand"))
close(con)
dominant <- select_dominant_tss(all_tx, cage, tss_distance)
if (!nrow(dominant)) stop("No CAGE-supported dominant TSSs within the matching distance")
dominant[, ensg := sub("\\..*", "", gene_id)]
message(nrow(dominant), " genes with a CAGE-supported dominant TSS; ",
        uniqueN(all_tx$gene_id) - nrow(dominant), " genes without a selected site")
tx <- tx[tx_id %in% genes$tx_id]
genes <- merge(genes, tx[, .(tx_id, tx_start = start, tx_end = end)], by = "tx_id")
stopifnot(genes[, all(tss == ifelse(strand == "+", tx_start, tx_end))])
genes[, tes := ifelse(strand == "+", tx_end, tx_start)]
genes[, c("tx_start", "tx_end") := NULL]

genes[, anchor_type := "canonical"]
dominant[, anchor_type := "dominant"]
genes <- rbindlist(list(genes, dominant), fill = TRUE)

# --- GM12878 expression ------------------------------------------------------
read_fpkm <- function(path) {
  d <- fread(path, select = c("gene_id", "FPKM"))
  d[grepl("^ENSG", gene_id), .(fpkm = max(FPKM)), by = .(ensg = sub("\\..*", "", gene_id))]
}
ex <- merge(read_fpkm(EXPR[1]), read_fpkm(EXPR[2]), by = "ensg",
            suffixes = c("_rep1", "_rep2"))
ex[, fpkm := (fpkm_rep1 + fpkm_rep2) / 2]
genes <- merge(genes, ex, by = "ensg", all.x = TRUE)
message(genes[, sum(!is.na(fpkm))], " of ", nrow(genes), " genes have a V29 FPKM")

# --- GRO-seq counts in each window -------------------------------------------
chain  <- import.chain(CHAIN)
bw_len <- seqlengths(BigWigFile(GRO[["+"]]))

#' Signal sum per window. `win` is a GRanges of 1-based inclusive windows, `sig` the
#' lifted hg38 signal for the matching strand.
window_sum <- function(win, sig) {
  out <- numeric(length(win))
  ov  <- findOverlaps(win, sig)
  if (!length(ov)) return(out)
  overlap_width <- width(pintersect(win[queryHits(ov)], sig[subjectHits(ov)]))
  agg <- tapply(sig$score[subjectHits(ov)] * overlap_width, queryHits(ov), sum)
  out[as.integer(names(agg))] <- as.numeric(agg)
  out
}

lib_reads <- 0
retained  <- list()
res <- rbindlist(lapply(names(bw_len), function(cc) {
  sig <- lapply(GRO, function(f) {
    gr <- import(f, which = GRanges(cc, IRanges(1L, bw_len[[cc]])))
    gr$score <- abs(gr$score)                        # minus-strand counts are stored negative
    gr
  })
  lib_reads <<- lib_reads + sum(sapply(sig, function(gr) sum(gr$score * width(gr))))
  g <- genes[chrom == cc]
  if (!nrow(g)) return(NULL)

  sig <- lapply(names(sig), function(s) {
    lo <- unlist(liftOver(sig[[s]], chain))
    lo <- lo[as.character(seqnames(lo)) == cc]
    retained[[paste(cc, s)]] <<- sum(lo$score * width(lo)) / sum(sig[[s]]$score * width(sig[[s]]))
    lo
  })
  names(sig) <- names(GRO)

  # windows oriented by transcription
  plus <- g$strand == "+"
  t_lo <- ifelse(plus, g$tss + TSSR[1],      g$tss - TSSR[2])
  t_hi <- ifelse(plus, g$tss + TSSR[2],      g$tss - TSSR[1])
  b_lo <- ifelse(plus, g$tss + BODY_FROM,    g$tes - BODY_PAST)
  b_hi <- ifelse(plus, g$tes + BODY_PAST,    g$tss - BODY_FROM)
  tssr <- GRanges(cc, IRanges(pmax(1L, t_lo), t_hi))
  body <- GRanges(cc, IRanges(pmax(1L, b_lo), b_hi))

  ts <- bs <- numeric(nrow(g))
  pause_peak <- rep(NA_integer_, nrow(g))
  for (s in c("+", "-")) {
    i <- which(g$strand == s)
    ts[i] <- window_sum(tssr[i], sig[[s]])
    bs[i] <- window_sum(body[i], sig[[s]])
    pause_peak[i] <- find_pause_peak(g$tss[i], s, sig[[s]])
  }
  data.table(g, pause_peak_rel = pause_peak,
             tssr_start = start(tssr), tssr_end = end(tssr), tssr_reads = ts,
             body_start = start(body), body_end = end(body), body_reads = bs)
}))

lib_m <- lib_reads / 1e6
res[, tssr_dens := tssr_reads / lib_m / (tssr_end - tssr_start + 1L)]
res[, body_dens := body_reads / lib_m / (body_end - body_start + 1L)]
res[, active    := tssr_dens >= MIN_DENS]
res[, pi        := fifelse(active, tssr_dens / body_dens, NA_real_)]
res[, paused    := active & pi > PAUSED_PI]          # no PI counts as non-paused
res[, expressed := fpkm > EXPR_FPKM]
setorder(res, chrom, tss)
dominant_res <- res[anchor_type == "dominant"]
dominant_res[, `:=`(anchor_tss = tss, anchor_tes = tes,
                    anchor_source = "CAGE", cage_file = normalizePath(opts$cage),
                    tss_match_max_distance = as.numeric(opts[["tss-distance"]]))]
dir.create(file.path(OUT_DIR, "dominant_tss"), recursive = TRUE, showWarnings = FALSE)
dominant_path <- file.path(OUT_DIR, "dominant_tss/pause_index_dominant.tsv")
fwrite(dominant_res, dominant_path, sep = "\t", na = "NA")
message("Wrote ", dominant_path, " (", nrow(dominant_res), " genes)")
# Keep the original canonical output schema and summaries.
res <- res[anchor_type == "canonical"]
res[, c("anchor_type", "pause_peak_rel", "cage_score", "ctss_id", "annotated_tss", "tss_distance") := NULL]
fwrite(res, file.path(OUT_DIR, "pause_index.tsv"), sep = "\t")

# --- summary -----------------------------------------------------------------
retained <- unlist(retained)
message("library: ", round(lib_m, 2), " M GRO-seq reads")
message("liftOver: ", round(100 * min(retained), 3), "% of signal retained on the worst chromosome (",
        names(which.min(retained)), ")")
message(res[, sum(active)], " genes with an active TSS (>= ", MIN_DENS, " rpm/bp); ",
        res[, sum(active & !is.finite(pi))], " of them have no gene-body signal (PI = Inf)")
report_paused <- function(e, label) {
  message(label, ": ", nrow(e), " expressed (FPKM > ", EXPR_FPKM, "), ", e[, sum(active)],
          " with a PI, ", e[, sum(paused)], " paused = ", round(100 * e[, mean(paused)], 1),
          "% of expressed, ", round(100 * e[active == TRUE, mean(paused)], 1), "% of those with a PI")
}
report_paused(res[expressed == TRUE], "all genes")
report_paused(res[expressed == TRUE & gene_type == "protein_coding"], "protein-coding")
e <- res[expressed == TRUE & is.finite(pi)]
message("PI 5/25/50/75/95% (expressed, finite): ",
        paste(signif(quantile(e$pi, c(.05, .25, .5, .75, .95)), 3), collapse = " / "))
message("Spearman PI vs FPKM (expressed, finite): ",
        round(cor(e$pi, e$fpkm, method = "spearman"), 3))
message("Wrote ", file.path(OUT_DIR, "pause_index.tsv"))
