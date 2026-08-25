#!/usr/bin/env Rscript
# Quantify THP-1 RNAPII (phospho-S5) ChIP-seq (GSE32324, hg19 BED6 aligned
# reads) at the early-response promoters: read coordinates are lifted
# hg19 -> hg38 (conservative liftOver as in
# code/enhancer_accessibility/liftover_tracks.R: reads that fail to map or
# split across blocks are dropped) and counted by overlap with the hg38
# TSS-centered windows. No genome-wide track is built.
#
# Per gene x sample: read counts and RPM densities in
#   promoter = TSS +/- 1000 bp        (the phasing window)
#   pause    = TSS -50..+300          (strand-aware)
#   body     = TSS +300..TES          (strand-aware)
# and pausing index PI = pause density / body density.
# RPM uses the sample's total mapped reads (before liftover) as library size.
#
# Usage:
#   /software/R-4.4.1-el8-x86_64/bin/Rscript quantify_polII_chip.R
# Output:
#   <output_dir>/PolII_chip_tss_counts.tsv
suppressMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
})

chain_path <- "/project/spott/cshan/annotations/hg19ToHg38.over.chain"
chip_dir   <- "/project/spott/cshan/fiber-seq/macrophage_project/annotations/PolII_ChIP"
output_dir <- "/project/spott/cshan/fiber-seq/macrophage_project/clustering_TSS_methods/hamming_hierarchical"

chip_samples <- c(control_rep1 = "GSM869208_THP1_control_1_RNAPII",
                  LPS_rep1     = "GSM869209_THP1_LPS_1_RNAPII",
                  control_rep2 = "GSM869210_THP1_control_2_RNAPII",
                  LPS_rep2     = "GSM869211_THP1_LPS_2_RNAPII")

# Promoter windows (TSS +/- 1000 bp, hg38) and GENCODE v49 gene spans, built
# by make_promoter_table.R (run it first).
promoters <- read.delim(file.path(output_dir, "promoter_regions_gencode_v49.tsv"))

region_grs <- list(
  promoter = GRanges(promoters$chr, IRanges(promoters$start, promoters$end)),
  pause = GRanges(promoters$chr, IRanges(
    ifelse(promoters$gene_strand == "+", promoters$tss - 50,  promoters$tss - 300),
    ifelse(promoters$gene_strand == "+", promoters$tss + 300, promoters$tss + 50))),
  body = GRanges(promoters$chr, IRanges(
    ifelse(promoters$gene_strand == "+", promoters$tss + 300, promoters$gene_start),
    ifelse(promoters$gene_strand == "+", promoters$gene_end,  promoters$tss - 300)))
)

chain      <- import.chain(chain_path)
keep_chrom <- unique(promoters$chr)
chunk_size <- 2e6L

out <- list()
for (s in names(chip_samples)) {
  bed <- file.path(chip_dir, paste0(chip_samples[[s]], ".bed.gz"))
  dt  <- fread(cmd = paste("zcat", shQuote(bed)), header = FALSE,
               select = 1:3, col.names = c("chrom", "start", "end"))
  n_total <- nrow(dt)
  dt <- dt[chrom %in% keep_chrom]
  gr_all <- GRanges(dt$chrom, IRanges(dt$start + 1L, dt$end))  # BED is 0-based
  rm(dt); invisible(gc())

  counts <- lapply(region_grs, function(gr) integer(length(gr)))
  n_kept <- 0L; n_lifted <- 0L
  for (cs in seq(1L, length(gr_all), by = chunk_size)) {
    ce  <- min(cs + chunk_size - 1L, length(gr_all))
    lifted <- liftOver(gr_all[cs:ce], chain)
    gr2 <- unlist(lifted[elementNROWS(lifted) == 1L])  # unique, contiguous only
    rm(lifted)
    n_kept   <- n_kept + (ce - cs + 1L)
    n_lifted <- n_lifted + length(gr2)
    for (rg in names(region_grs))
      counts[[rg]] <- counts[[rg]] + countOverlaps(region_grs[[rg]], gr2)
    rm(gr2); invisible(gc())
  }
  message(sprintf("%s: %d reads total, %d on target chromosomes, %d lifted",
                  chip_samples[[s]], n_total, n_kept, n_lifted))

  rpm_per_kb <- function(n, gr) n / n_total * 1e6 / (width(gr) / 1000)
  out[[s]] <- data.frame(
    gene           = promoters$gene,
    sample         = s,
    n_library      = n_total,
    promoter_count = counts$promoter,
    promoter_rpm_per_kb = rpm_per_kb(counts$promoter, region_grs$promoter),
    pause_count    = counts$pause,
    pause_rpm_per_kb    = rpm_per_kb(counts$pause, region_grs$pause),
    body_count     = counts$body,
    body_rpm_per_kb     = rpm_per_kb(counts$body, region_grs$body),
    row.names = NULL
  )
  out[[s]]$PI <- out[[s]]$pause_rpm_per_kb / out[[s]]$body_rpm_per_kb
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
res <- do.call(rbind, out)
rownames(res) <- NULL
write.table(res, file.path(output_dir, "PolII_chip_tss_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("written: ", file.path(output_dir, "PolII_chip_tss_counts.tsv"))
