# Input CTSS coordinates are hg38 BED; transcript coordinates are 1-based.
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

# Strongest GRO-seq base in transcription-relative +20..+80; no signal -> NA.
# Ties select the earliest relative base on either strand.
find_pause_peak <- function(tss, strand, signal) {
  result <- rep(NA_integer_, length(tss))
  if (!length(signal) || !length(tss)) return(result)
  offsets <- 20:80
  positions <- rep(tss, each = length(offsets)) +
    rep(if (strand == "+") offsets else -offsets, times = length(tss))
  query <- GRanges(as.character(seqnames(signal))[1L], IRanges(positions, positions))
  hits <- findOverlaps(query, signal, ignore.strand = TRUE)
  values <- numeric(length(positions))
  if (length(hits)) {
    sums <- tapply(signal$score[subjectHits(hits)], queryHits(hits), sum)
    values[as.integer(names(sums))] <- sums
  }
  values <- matrix(values, nrow = length(tss), byrow = TRUE)
  best <- max.col(values, ties.method = "first")
  detected <- values[cbind(seq_along(tss), best)] > 0
  result[detected] <- offsets[best[detected]]
  result
}
