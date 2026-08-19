#!/usr/bin/env Rscript
# LiftOver a single public track from hg19 -> hg38 with rtracklayer.
#
# Input type is inferred from the file name:
#   *.narrowPeak[.gz]  interval liftOver, re-exported as (plain) narrowPeak
#   *.bw / *.bigWig    per-interval signal liftOver, re-exported as bigWig
#
# narrowPeak: only peaks that map to a single contiguous hg38 block are kept
#   (the standard conservative liftOver behaviour); peaks that fail to map or that
#   split across blocks are dropped and counted. The summit offset (col 10) is
#   clamped into the lifted peak so the output stays a valid narrowPeak.
#
# bigWig: signal intervals are lifted individually, then collapsed to a disjoint
#   per-base track with score-weighted coverage() so the result can be written as
#   bigWig. Within-block signal is preserved exactly; at the rare block boundaries
#   where two hg19 intervals land on overlapping hg38 bases the scores are summed.
#
# Usage:
#   liftover_tracks.R <chain> <hg38.chrom.sizes> <input> <output>
#
# Memory: a whole-genome bigWig is imported at once (a few GB RAM, ~2-3 min).
suppressMessages({
  library(rtracklayer)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L)
  stop("usage: liftover_tracks.R <chain> <chrom.sizes> <input> <output>")
chain_path  <- args[[1]]
sizes_path  <- args[[2]]
in_path     <- args[[3]]
out_path    <- args[[4]]

chain <- import.chain(chain_path)
cs    <- read.table(sizes_path, col.names = c("chrom", "size"))
si    <- Seqinfo(cs$chrom, cs$size, genome = "hg38")

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

# Restrict to hg38 seqlevels we have lengths for, attach those lengths, and trim
# any range that liftOver pushed past a chromosome end.
constrain_to_hg38 <- function(gr) {
  keep <- intersect(as.character(unique(seqnames(gr))), seqnames(si))
  gr   <- keepSeqlevels(gr, keep, pruning.mode = "coarse")
  seqlevels(gr) <- keep
  seqlengths(gr) <- seqlengths(si)[keep]
  gr <- trim(gr)
  gr[width(gr) > 0L]
}

is_narrowpeak <- grepl("\\.narrowPeak(\\.gz)?$", in_path, ignore.case = TRUE)
is_bigwig     <- grepl("\\.(bw|bigWig)$",        in_path, ignore.case = TRUE)

if (is_narrowpeak) {
  gr     <- import(in_path, format = "narrowPeak")
  lifted <- liftOver(gr, chain)
  nseg   <- elementNROWS(lifted)

  gr2 <- unlist(lifted[nseg == 1L])          # keep unique, contiguous mappings only
  gr2 <- constrain_to_hg38(gr2)

  # col 10 (peak) is the 0-based summit offset from chromStart; clamp into range.
  gr2$peak <- pmax(0L, pmin(gr2$peak, width(gr2) - 1L))
  gr2 <- sort(gr2)

  # Write the 10-column narrowPeak by hand: rtracklayer's narrowPeak *writer*
  # only emits BED6 and drops the signalValue/pValue/qValue/peak columns.
  str <- as.character(strand(gr2)); str[!str %in% c("+", "-")] <- "."
  df <- data.frame(chrom = as.character(seqnames(gr2)),
                   start = start(gr2) - 1L, end = end(gr2),
                   name = gr2$name, score = gr2$score, strand = str,
                   signalValue = gr2$signalValue, pValue = gr2$pValue,
                   qValue = gr2$qValue, peak = gr2$peak, check.names = FALSE)
  write.table(df, out_path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  message(sprintf("%s: %d peaks -> %d kept (%d unmapped, %d split-dropped)",
                  basename(in_path), length(gr), length(gr2),
                  sum(nseg == 0L), sum(nseg > 1L)))

} else if (is_bigwig) {
  gr     <- import.bw(in_path)
  lifted <- unlist(liftOver(gr, chain))
  lifted <- constrain_to_hg38(lifted)

  cov <- coverage(lifted, weight = "score")  # disjoint per-base signal (overlap-safe)
  export.bw(cov, out_path)
  message(sprintf("%s: %d intervals -> %d lifted -> bigWig",
                  basename(in_path), length(gr), length(lifted)))

} else {
  stop("unrecognized input type (expected *.narrowPeak[.gz], *.bw, or *.bigWig): ", in_path)
}
