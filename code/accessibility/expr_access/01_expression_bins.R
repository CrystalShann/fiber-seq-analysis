#' Expression bins + TSS anchors for the TSS m6A (SAM-seq accessibility)
#' metaprofiles.
#'
#' Builds one table with a strand-aware TSS anchor and an expression bin for
#' every protein-coding gene:
#'   * counts: STAR ReadsPerGene column 4 (reverse-stranded, the column chosen
#'     in code/RNA/LPS_timecourse_DESeq2_LRT.Rmd), UNFILTERED so silent genes
#'     stay in the table; only the 0/5/10/15 min samples (the timepoints with
#'     fiber-seq m6A) enter the mean;
#'   * TPM from union-exon gene lengths (gencode v46, the STAR annotation);
#'   * TSS: Ensembl_canonical transcript TSS (gencodev46_Ensembl_canonical_TSS
#'     .bed; the 20 bp interval is TSS +/- 10, so tss0 = start + 10, 0-based);
#'   * bins: mean TPM < 1 = not_expressed, else quartiles Q1(low)..Q4(high).
#'
#' Output: <OUT_DIR>/tss_expression_bins.tsv

suppressPackageStartupMessages(library(data.table))

ALIGNED_DIR <- "/project/spott/cshan/fiber-seq/code/RNA/aligned"
GTF_PATH    <- "/project/spott/cshan/annotations/gencode.v46.annotation.gtf.gz"
TSS_BED     <- "/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed"
OUT_DIR     <- "/project/spott/cshan/fiber-seq/macrophage_project/expr_access/tables"

USE_TIMES  <- c(0L, 5L, 10L, 15L)   # minutes of LPS; the fiber-seq timepoints
TPM_EXPR   <- 1                     # mean TPM below this = not_expressed
KEEP_CHROM <- paste0("chr", c(1:22, "X", "Y"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- counts: unfiltered gene counts for the 0-15 min samples -----------------
samples <- list.dirs(ALIGNED_DIR, recursive = FALSE, full.names = FALSE)
sm <- data.table(sample = samples)
sm[, c("rep", "time") := {
  core <- sub("^SP-AL-18s-R", "", sub("_S\\d+$", "", sample))
  .(as.integer(substr(core, 1, 1)), as.integer(substring(core, 2)))
}]
sm <- sm[time %in% USE_TIMES][order(time, rep)]
stopifnot(nrow(sm) == length(USE_TIMES) * 3L)
message("samples: ", paste(sm$sample, collapse = " "))

counts <- sapply(sm$sample, function(s) {
  f  <- file.path(ALIGNED_DIR, s, paste0(s, "_ReadsPerGene.out.tab"))
  dt <- fread(f, header = FALSE, select = c(1, 4),
              col.names = c("gene_id", "n"))
  dt <- dt[!grepl("^N_", gene_id)]
  setkey(dt, gene_id)
  dt$n
})
gene_ids <- fread(file.path(ALIGNED_DIR, sm$sample[1],
                            paste0(sm$sample[1], "_ReadsPerGene.out.tab")),
                  header = FALSE, select = 1)[!grepl("^N_", V1), sort(V1)]
rownames(counts) <- gene_ids
message(nrow(counts), " genes x ", ncol(counts), " samples")

# --- union-exon gene lengths from the gencode v46 GTF ------------------------
ex <- fread(cmd = paste("zcat", GTF_PATH, "| awk -F'\t' '$3 == \"exon\"'",
                        "| cut -f1,4,5,9"),
            header = FALSE, col.names = c("chrom", "start", "end", "attr"),
            quote = "")
ex[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", attr)][, attr := NULL]
setorder(ex, gene_id, start, end)
# union length: merge overlapping exons per gene (touching islands sum exactly)
ex[, prev_max := shift(cummax(end), 1L), by = gene_id]
ex[, island := cumsum(is.na(prev_max) | start > prev_max), by = gene_id]
glen <- ex[, .(s = min(start), e = max(end)), by = .(gene_id, island)][
  , .(len = sum(e - s + 1L)), by = gene_id]
stopifnot(all(gene_ids %in% glen$gene_id))
len <- glen[match(gene_ids, gene_id), len]

# --- TPM (denominator over all genes), mean across the 12 samples ------------
rpk <- counts / (len / 1e3)
tpm <- t(t(rpk) / colSums(rpk)) * 1e6
expr <- data.table(gene_id = gene_ids, mean_tpm = rowMeans(tpm))

# --- canonical TSS, protein-coding only --------------------------------------
bed <- fread(TSS_BED, header = FALSE,
             col.names = c("chrom", "start", "end", "name", "score", "strand"))
bed[, c("gene_id", "tx_id", "gene_name", "gene_type") :=
      tstrsplit(name, ";", keep = 1:4)]
bed[, tss := start + 10L]  # interval is TSS +/- 10 (verified against the GTF)
pc <- bed[gene_type == "protein_coding" & chrom %in% KEEP_CHROM,
          .(gene_id, gene_name, chrom, tss, strand)]
message(nrow(pc), " protein-coding genes on ", uniqueN(pc$chrom), " chromosomes")

out <- merge(pc, expr, by = "gene_id")
message(nrow(out), "/", nrow(pc), " protein-coding genes matched to counts")
stopifnot(nrow(out) / nrow(pc) > 0.99)

# --- expression bins ---------------------------------------------------------
qs <- out[mean_tpm >= TPM_EXPR, quantile(mean_tpm, c(0.25, 0.5, 0.75))]
out[, expr_bin := fifelse(
  mean_tpm < TPM_EXPR, "not_expressed",
  as.character(cut(mean_tpm, breaks = c(TPM_EXPR, qs, Inf),
                   labels = c("Q1_low", "Q2", "Q3", "Q4_high"),
                   include.lowest = TRUE, right = FALSE)))]
print(out[, .(n = .N, min_tpm = min(mean_tpm), max_tpm = max(mean_tpm)),
          by = expr_bin][order(min_tpm)])

setorder(out, chrom, tss)
fwrite(out[, .(gene_id, gene_name, chrom, tss, strand, mean_tpm, expr_bin)],
       file.path(OUT_DIR, "tss_expression_bins.tsv"), sep = "\t")
message("Wrote ", file.path(OUT_DIR, "tss_expression_bins.tsv"))
