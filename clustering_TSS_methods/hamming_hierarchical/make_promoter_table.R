#!/usr/bin/env Rscript
# Build the early-response promoter table from the GENCODE v49 annotation:
# one row per gene with the strand-aware GENCODE TSS (gene 5' end), the
# TSS +/- 1000 bp promoter window used for phasing (the paper's 2-kb
# Pol II-analysis geometry), and the gene span used for the ChIP pausing
# index. Restricted to primary chromosomes (drops the TNF alt-contig copies).
# IL8 is CXCL8 in GENCODE; the project's gene naming (IL8) is kept.
#
# Usage:
#   /software/R-4.4.1-el8-x86_64/bin/Rscript make_promoter_table.R
# Output:
#   <output_dir>/promoter_regions_gencode_v49.tsv
suppressMessages(library(data.table))

gtf_path   <- "/project/spott/cshan/annotations/gencode.v49.chr_patch_hapl_scaff.annotation.gtf"
output_dir <- "/project/spott/cshan/fiber-seq/macrophage_project/clustering_TSS_methods/hamming_hierarchical"
flank      <- 1000L

gene_symbols <- c(IL1B = "IL1B", TNF = "TNF", IL1A = "IL1A", PTGS2 = "PTGS2",
                  IL8 = "CXCL8", CXCL2 = "CXCL2", IL6 = "IL6", RELB = "RELB")

gtf <- fread(cmd = sprintf("awk -F'\t' '$3==\"gene\"' %s", shQuote(gtf_path)),
             header = FALSE,
             col.names = c("chr", "source", "feature", "start", "end",
                           "score", "strand", "frame", "attr"))
gtf <- gtf[chr %in% paste0("chr", c(1:22, "X", "Y"))]
gtf[, ensg      := sub('.*gene_id "([^"]+)".*',   "\\1", attr)]
gtf[, gene_name := sub('.*gene_name "([^"]+)".*', "\\1", attr)]
gtf <- gtf[gene_name %in% gene_symbols]

stopifnot(!anyDuplicated(gtf$gene_name),
          all(gene_symbols %in% gtf$gene_name))
gtf <- gtf[match(gene_symbols, gene_name)]

promoters <- data.table(
  gene        = names(gene_symbols),
  gencode_name = gtf$gene_name,
  ensg        = gtf$ensg,
  chr         = gtf$chr,
  tss         = ifelse(gtf$strand == "+", gtf$start, gtf$end),
  gene_start  = gtf$start,
  gene_end    = gtf$end,
  gene_strand = gtf$strand
)
promoters[, start := tss - flank]
promoters[, end   := tss + flank]
setcolorder(promoters, c("gene", "gencode_name", "ensg", "chr", "start", "end",
                         "tss", "gene_start", "gene_end", "gene_strand"))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
out <- file.path(output_dir, "promoter_regions_gencode_v49.tsv")
write.table(promoters, out, sep = "\t", quote = FALSE, row.names = FALSE)
message("written: ", out)
print(promoters)
