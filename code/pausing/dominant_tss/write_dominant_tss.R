#' CAGE dominant TSS per gene -> annotations/fantom5_LCL/fantom5.hg38.LCL.dominant_TSS.bed.gz
#'
#' Sites: fantom5.hg38.LCL.consensus.TSS.bed.gz (build_ctss_consensus.R). Transcripts: GENCODE
#' v46 on chr1-22 and chrX. Selection rules: select_dominant_tss() below and README.md.
#' Output: BED6+6, sorted, bgzipped and tabix-indexed; header line starts with "#".
#'   #chrom start end gene_name cage_score strand gene_id gene_type tx_id tes annotated_tss
#'   tss_distance   (start 0-based, end = 1-based TSS; tes and annotated_tss 1-based)
#' Run: bash code/pausing/dominant_tss/run_write_dominant_tss.sh

suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })

ANN      <- "/project/spott/cshan/annotations"
GTF      <- file.path(ANN, "gencode.v46.annotation.gtf.gz")
CAGE     <- file.path(ANN, "fantom5_LCL/fantom5.hg38.LCL.consensus.TSS.bed.gz")
OUT      <- file.path(ANN, "fantom5_LCL/fantom5.hg38.LCL.dominant_TSS.bed")
BIN      <- "/project/spott/cshan/envs/dimelo/bin"
CHROMS   <- paste0("chr", c(1:22, "X"))
MAX_DIST <- 250L   # bp between a CAGE site and an annotated transcript TSS

# Dominant TSS per gene from CTSS sites (hg38 BED) and transcripts (1-based)
select_dominant_tss <- function(tx, cage, max_distance = 250L) {
  stopifnot(nrow(tx) > 0L, !anyDuplicated(tx$tx_id),
            !anyNA(cage[, .(chrom, start0, end, score, strand)]),
            all(cage$end == cage$start0 + 1), all(cage$start0 >= 0),
            all(cage$start0 %% 1 == 0), all(is.finite(cage$score)),
            all(cage$score >= 0), all(cage$strand %in% c("+", "-")))
  cage <- copy(cage[score > 0])
  if (anyDuplicated(cage[, .(chrom, start0, strand)])) stop("Duplicate CTSS coordinates")
  windows <- GRanges(tx$chrom, IRanges(pmax(1L, tx$tss - max_distance),
                                     tx$tss + max_distance), strand = tx$strand)
  sites <- GRanges(cage$chrom, IRanges(cage$start0 + 1L, cage$end), strand = cage$strand)
  hits <- findOverlaps(sites, windows, ignore.strand = FALSE)
  candidates <- data.table(site = queryHits(hits), transcript = subjectHits(hits))
  candidates[, `:=`(gene_id = tx$gene_id[transcript],
                    distance = abs(cage$start0[site] + 1L - tx$tss[transcript]))]
  candidates <- candidates[ifelse(tx$strand[transcript] == "+",
                                 cage$start0[site] + 1L <= tx$tes[transcript],
                                 cage$start0[site] + 1L >= tx$tes[transcript])]
  # Assign a CTSS to its nearest same-strand gene; drop ties between genes.
  if (nrow(candidates))
    candidates <- candidates[candidates[, .I[distance == min(distance)], by = site]$V1]
  candidates <- candidates[site %in% candidates[, .(n = uniqueN(gene_id)), by = site][n == 1L, site]]
  candidates[, `:=`(score = cage$score[site], canonical = tx$canonical[transcript],
                    position = cage$start0[site] + 1L, tx_id = tx$tx_id[transcript])]
  # Strongest base wins; deterministic ties prefer distance, canonical transcript,
  # then genomic coordinate and transcript ID. TES belongs to that transcript.
  setorderv(candidates, c("gene_id", "score", "distance", "canonical", "position", "tx_id"),
            c(1L, -1L, 1L, -1L, 1L, 1L))
  chosen <- candidates[!duplicated(gene_id)]
  out <- tx[chosen$transcript, .(gene_id, tx_id, gene_name, gene_type, chrom, strand, tes)]
  out[, `:=`(tss = chosen$position, annotated_tss = tx$tss[chosen$transcript],
             tss_distance = chosen$distance, cage_score = chosen$score,
             ctss_id = cage$ctss_id[chosen$site])]
  out <- out[(strand == "+" & tes >= tss) | (strand == "-" & tes <= tss)]
  out
}

tx <- fread(cmd = paste("zcat", GTF, "| awk -F'\t' '$3 == \"transcript\"' | cut -f1,4,5,7,9"),
            header = FALSE, col.names = c("chrom", "start", "end", "strand", "attr"), quote = "")
tx[, `:=`(tx_id = sub('.*transcript_id "([^"]+)".*', "\\1", attr),
          gene_id = sub('.*gene_id "([^"]+)".*', "\\1", attr),
          gene_name = sub('.*gene_name "([^"]+)".*', "\\1", attr),
          gene_type = sub('.*gene_type "([^"]+)".*', "\\1", attr),
          canonical = grepl('tag "Ensembl_canonical"', attr, fixed = TRUE))]
tx[, `:=`(tss = fifelse(strand == "+", start, end), tes = fifelse(strand == "+", end, start))]
tx <- tx[chrom %in% CHROMS]
cage <- fread(CAGE, header = FALSE, select = 1:6,
              col.names = c("chrom", "start0", "end", "ctss_id", "score", "strand"))

d <- select_dominant_tss(tx, cage, MAX_DIST)
message(nrow(d), " genes with a CAGE dominant TSS; ", uniqueN(tx$gene_id) - nrow(d), " without")

b <- d[, .(chrom, start = tss - 1L, end = tss, gene_name, cage_score, strand,
           gene_id, gene_type, tx_id, tes, annotated_tss, tss_distance)]
b <- b[order(chrom, start)]
setnames(b, "chrom", "#chrom")
fwrite(b, OUT, sep = "\t", quote = FALSE)
stopifnot(system2(file.path(BIN, "bgzip"), c("-f", OUT)) == 0L,
          system2(file.path(BIN, "tabix"), c("-f", "-p", "bed", paste0(OUT, ".gz"))) == 0L)
message("Wrote ", OUT, ".gz")
