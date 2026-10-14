#' RNA Pol II pausing index for GM12878, in hg38.
#'
#' The pausing index follows Day et al. 2016 (Genome Biology 17:120; PMID 27259512,
#' PMC4893286, doi:10.1186/s13059-016-0984-2):
#'
#'   PI = (TSSR signal / TSSR width) / (gene-body signal / gene-body width)
#'   TSSR      = -50 .. +300 around the TSS
#'   gene body = +300 .. TES + 3kb
#'
#' Day et al. compute this from RNAP2 ChIP-seq on hg19/RefSeq and publish no per-gene
#' table, so it is recomputed here. The nascent-transcription data are GM12878 GRO-seq,
#' GEO GSM1480326 (series GSE60454/GSE60456) from Core et al. 2014 (Nat Genet
#' 46:1311-20; PMID 25383968, doi:10.1038/ng.3142) -- a different study from Day et al.
#' Those bigWigs are hg19, so the signal is lifted to hg38 here and everything
#' downstream stays in hg38.
#'
#' TSS anchors are the Ensembl_canonical TSS (1 bp BED rows, so tss = start, 0-based)
#' and the TES is taken from the SAME canonical transcript in the v46 GTF, joined by
#' transcript ID. Gene-record coordinates must not be used: for ACTB the Ensembl gene
#' record puts the 5' end 33 kb from the real promoter and the PI comes out 20 instead
#' of 155.
#'
#' Output: <OUT_DIR>/tables/pause_index.tsv  (one row per protein-coding canonical gene)
#'         <OUT_DIR>/tables/liftover_qc.tsv  (signal retained per chromosome)
#'
#' rtracklayer will not load unless a newer libcurl precedes the system one on
#' LD_LIBRARY_PATH; run this through run_lcl_tss.sh, which sets it.

suppressPackageStartupMessages({
  library(data.table); library(rtracklayer); library(GenomicRanges)
})

ANN      <- "/project/spott/cshan/annotations"
LCL_PUB  <- file.path(ANN, "LCL_public_data")
TSS_BED  <- file.path(ANN, "gencodev46_Ensembl_canonical_TSS.bed")
GTF      <- file.path(ANN, "gencode.v46.annotation.gtf.gz")
CHAIN    <- file.path(ANN, "hg19ToHg38.over.chain")
GRO      <- c(`+` = file.path(LCL_PUB, "GRO-seq/GSM1480326_GM12878_GROseq_plus.bigWig"),
              `-` = file.path(LCL_PUB, "GRO-seq/GSM1480326_GM12878_GROseq_minus.bigWig"))
EXPR     <- file.path(LCL_PUB, "gene_quan_ENCSR820PHH",
                      c("ENCFF345SHY.gene_quan.tsv", "ENCFF873VWU.gene_quan.tsv"))
OUT_DIR  <- "/project/spott/cshan/fiber-seq/LCL_project/auto_correlation/lcl_tss"

CHROMS    <- paste0("chr", 1:22)
TSSR      <- c(-50L, 300L)    # Day et al. TSSR, TSS-relative, transcriptional
BODY_FROM <- 300L             # Day et al. body start, TSS-relative
BODY_PAST <- 3000L            # Day et al. body extension past the TES
ISO_FLANK <- 2000L            # a gene is "isolated" with no other canonical TSS within this

stopifnot(all(file.exists(c(TSS_BED, GTF, CHAIN, GRO, EXPR))))
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

# --- canonical TSS anchors ---------------------------------------------------
bed <- fread(TSS_BED, header = FALSE,
             col.names = c("chrom", "start", "end", "name", "score", "strand"))
bed[, c("gene_id", "tx_id", "gene_name", "gene_type", "tags") := tstrsplit(name, ";", keep = 1:5)]
stopifnot(all(bed$end - bed$start == 1L))          # 1 bp rows: tss = start (0-based)

# every canonical TSS, any gene type: the isolation reference
all_tss <- GRanges(bed$chrom, IRanges(bed$start + 1L, width = 1L))

pc <- bed[gene_type == "protein_coding" & grepl("Ensembl_canonical", tags) &
            chrom %in% CHROMS,
          .(ensg = sub("\\..*", "", gene_id), tx_id, gene_name, chrom,
            tss = start, strand)]
message(nrow(pc), " protein-coding canonical TSS on chr1-22")

# --- TES of the same canonical transcript ------------------------------------
tx <- fread(cmd = paste("zcat", GTF, "| awk -F'\t' '$3 == \"transcript\"' | cut -f1,4,5,7,9"),
            header = FALSE, col.names = c("chrom", "start", "end", "strand", "attr"),
            quote = "")
tx[, tx_id := sub('.*transcript_id "([^"]+)".*', "\\1", attr)]
tx <- tx[tx_id %in% pc$tx_id]
tx[, tes := ifelse(strand == "+", end, start)]     # 1-based TES of the canonical transcript
pc <- merge(pc, tx[, .(tx_id, tes)], by = "tx_id")
message(nrow(pc), " matched to a canonical-transcript TES")

# --- GM12878 expression, both replicates -------------------------------------
read_tpm <- function(path) {
  d <- fread(path, select = c(1, 6), col.names = c("gene_id", "tpm"))
  d <- d[grepl("^ENSG", gene_id)]
  d[, .(ensg = sub("\\..*", "", gene_id), tpm)][, .(tpm = max(tpm)), by = ensg]
}
ex <- merge(read_tpm(EXPR[1]), read_tpm(EXPR[2]), by = "ensg", suffixes = c("1", "2"))
ex[, mean_tpm := (tpm1 + tpm2) / 2]
pc <- merge(pc, ex, by = "ensg")
message(nrow(pc), " with GM12878 expression (version-stripped ENSG join)")

# --- annotations used by the downstream filters ------------------------------
pc[, n_near := countOverlaps(GRanges(chrom, IRanges(tss + 1L - ISO_FLANK,
                                                    tss + 1L + ISO_FLANK)), all_tss) - 1L]
pc[, gene_len := abs(tes - tss)]

# --- GRO-seq, lifted hg19 -> hg38, scored in the Day et al. windows ----------
chain   <- import.chain(CHAIN)
hg19len <- seqlengths(BigWigFile(GRO[["+"]]))

#' Signal sum per window. `win` is a GRanges of 1-based inclusive windows on `cc`;
#' `sig` is the lifted hg38 GRanges for the matching strand.
window_sum <- function(win, sig) {
  out <- numeric(length(win))
  if (!length(sig)) return(out)
  ov  <- findOverlaps(win, sig)
  if (!length(ov)) return(out)
  agg <- tapply(sig$score[subjectHits(ov)], queryHits(ov), sum)
  out[as.integer(names(agg))] <- as.numeric(agg)
  out
}

qc  <- list()
res <- rbindlist(lapply(CHROMS, function(cc) {
  g <- pc[chrom == cc]
  if (!nrow(g)) return(NULL)
  sig <- lapply(names(GRO), function(s) {
    gr  <- import(GRO[[s]], which = GRanges(cc, IRanges(1L, hg19len[[cc]])), as = "GRanges")
    tot <- sum(abs(gr$score) * width(gr))
    lo  <- unlist(liftOver(gr, chain))
    lo  <- lo[as.character(seqnames(lo)) == cc]
    lo$score <- abs(lo$score)
    qc[[length(qc) + 1L]] <<- data.table(
      chrom = cc, strand = s, n_hg19 = length(gr), n_hg38 = length(lo),
      signal_hg19 = tot, signal_hg38 = sum(lo$score * width(lo)),
      pct_retained = 100 * sum(lo$score * width(lo)) / tot)
    lo
  })
  names(sig) <- names(GRO)

  plus <- g$strand == "+"
  # Day et al. windows, oriented by transcription, as 1-based inclusive intervals
  t_lo <- ifelse(plus, g$tss + 1L + TSSR[1],  g$tss + 1L - TSSR[2])
  t_hi <- ifelse(plus, g$tss + 1L + TSSR[2],  g$tss + 1L - TSSR[1])
  b_lo <- ifelse(plus, g$tss + 1L + BODY_FROM, g$tes - BODY_PAST)
  b_hi <- ifelse(plus, g$tes + BODY_PAST,      g$tss + 1L - BODY_FROM)

  body_ok <- b_hi >= b_lo
  tssr <- GRanges(cc, IRanges(pmax(1L, t_lo), pmax(1L, t_hi)))
  body <- GRanges(cc, IRanges(pmax(1L, ifelse(body_ok, b_lo, 1L)),
                              pmax(1L, ifelse(body_ok, b_hi, 1L))))

  ts <- bs <- numeric(nrow(g))
  for (s in c("+", "-")) {
    i <- which(if (s == "+") plus else !plus)
    if (!length(i)) next
    ts[i] <- window_sum(tssr[i], sig[[s]])
    bs[i] <- window_sum(body[i], sig[[s]])
  }
  data.table(g,
             tssr_start = start(tssr), tssr_end = end(tssr), tssr_w = width(tssr),
             body_start = ifelse(body_ok, start(body), NA_integer_),
             body_end   = ifelse(body_ok, end(body),   NA_integer_),
             body_w     = ifelse(body_ok, width(body), NA_integer_),
             tssr_sum = ts, body_sum = ifelse(body_ok, bs, NA_real_))
}))

res[, pi := (tssr_sum / tssr_w) / (body_sum / body_w)]
setorder(res, chrom, tss)

qc <- rbindlist(qc)
fwrite(qc, file.path(OUT_DIR, "tables", "liftover_qc.tsv"), sep = "\t")
fwrite(res, file.path(OUT_DIR, "tables", "pause_index.tsv"), sep = "\t")

# --- checks ------------------------------------------------------------------
stopifnot(nrow(res) > 15000L,
          res[, all(tssr_w == TSSR[2] - TSSR[1] + 1L)],
          min(qc$pct_retained) > 99)

usable <- res[is.finite(pi) & tssr_sum >= 50 & body_sum >= 20 & mean_tpm >= 1 &
                tpm1 > 0 & tpm2 > 0]
rho <- cor(log10(usable$pi), log10(usable$mean_tpm), method = "spearman")

message("liftOver: ", round(min(qc$pct_retained), 3), "% signal retained on the worst chromosome (",
        qc[which.min(pct_retained), paste0(chrom, " ", strand)], "); the shortfall is hg19 sequence with no hg38 counterpart")
message(nrow(usable), " genes with a usable PI")
message("PI 5/25/50/75/95%: ",
        paste(round(quantile(usable$pi, c(.05, .25, .5, .75, .95)), 1), collapse = " / "))
message("Spearman log10(PI) vs log10(TPM): ", round(rho, 3))
if (rho > 0) warning("PI is positively correlated with expression; expected around -0.2")
message("Wrote ", file.path(OUT_DIR, "tables", "pause_index.tsv"))
