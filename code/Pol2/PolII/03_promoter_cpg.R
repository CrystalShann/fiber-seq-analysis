#!/usr/bin/env Rscript
# CpG methylation of every curated PolII example promoter, from the same full-span reads as
# 02_extract_read_features.R. The pb-CpG-tools extraction is pol2_cluster_cpg() of
# clustering_Pol2, sourced below and left unchanged: each gene is passed as its own "cluster",
# so its reads go into their own BAM and only its own window is summed.
#
# Run with the environment of run_gene_examples.sh, after 02_extract_read_features.R, e.g.
#   sbatch --account=pi-spott --partition=caslake --cpus-per-task=8 --mem=32G --time=04:00:00 \
#     --output=/project/spott/cshan/fiber-seq/results/logs/promoter_cpg_%j.log \
#     --wrap "bash code/Pol2/PolII/run_gene_examples.sh code/Pol2/PolII/03_promoter_cpg.R"
#
# Output: LCL_project/Pol2/PolII/read_features/tables/promoter_cpg.tsv, one row per gene
#   promoter window = TSS +/- FLANK (EXAMPLE_RANGE), the window of read_features
#   mod / cov = methylated / covering read calls summed over the window's CpGs (count mode)
#   cpg = mod / cov; n_cpg_sites = CpGs with >= 1 covering read; n_reads = reads in the BAM

dyn.load("/project/spott/cshan/envs/dimelo/lib/libcurl.so.4", local = FALSE)
source("/project/spott/cshan/fiber-seq/code/Pol2/PolII/pol2_pausing_functions.r")
UTILS <- file.path(PROJECT, "code/Pol2/clustering_Pol2/pol2_chip_cluster_utils")
source(file.path(UTILS, "pol2_clustering_functions.R"))   # pol2_gr(), pol2_required()
source(file.path(UTILS, "pol2_clustering_plots.R"))       # pol2_cluster_cpg()
setDTthreads(1L)

TABLES <- file.path(OUT_DIR, "read_features/tables")

main <- function() {
  cfg <- yaml::read_yaml(file.path(UTILS, "config.yaml"))
  genes <- fread(file.path(TABLES, "gene_features.tsv"))[status == "complete",
    .(group, gene_name, gene_id, chrom, tss, strand, tss_source, n_reads)]
  reads <- fread(file.path(TABLES, "read_features.tsv.gz"), select = c("gene", "RID", "sample_name"))
  stopifnot(nrow(genes) > 0L, setequal(genes$gene_name, reads$gene))

  # one read can sit in two nearby promoters, so the gene is part of the molecule ID
  x <- list(molecules = reads[, .(RID = paste(gene, RID, sep = "::"),
                                  original_RID = sub("^[^:]+::", "", RID),
                                  sample_name, tss_id = gene, cluster = gene)])
  prep <- list(tss = genes[, .(tss_id = gene_name, chrom,
                               window_start = tss + EXAMPLE_RANGE[1] - 1L, window_end = tss + EXAMPLE_RANGE[2])],
               samples = fread(SAMPLES)[, .(sample_name, merged_bam_file)])
  cpg <- pol2_cluster_cpg(cfg, x, prep)$by_tss

  out <- merge(genes, cpg[, .(gene_name = tss_id, mod, cov, cpg, n_cpg_sites, median_cpg_coverage)],
               by = "gene_name", all.x = TRUE)
  setcolorder(out, c("group", "gene_name"))
  setorder(out, group, -cpg)
  stopifnot(nrow(out) == nrow(genes))
  fwrite(out, file.path(TABLES, "promoter_cpg.tsv"), sep = "\t", na = "NA")
  message("Wrote ", nrow(out), " promoters; ", sum(is.na(out$cpg)), " without CpG calls")
}

main()
