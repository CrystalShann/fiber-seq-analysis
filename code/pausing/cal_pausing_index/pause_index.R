#' GM12878 RNA Pol II pausing index after Day et al. 2016 (Genome Biol 17:120, Methods
#' "Calculation of PI from RNAP2 ChIP-seq data"), on RefSeq isoforms.
#'
#' Inputs (all GRCh38)
#'   isoforms     annotations/hg38.refGene.gtf.gz (UCSC refGene) transcripts on chr1-22 and chrX;
#'                genes with isoforms on more than one chromosome or strand are dropped
#'   Pol II ChIP  ENCSR000EAD (Snyder lab, 8WG16; Day et al.'s GM12878 set, GSE31477) reps 1-7 and
#'                its input ENCSR000EAF reps 1-2, each pooled; ENCODE filtered alignments as BED
#'                (annotations/LCL_public_data/Pol2_chip/prepare_pol2_chip.sh)
#'   expression   ENCSR820PHH gene_quan FPKM, mean of the two replicates (GENCODE V29), joined on
#'                the GENCODE v46 gene_id of the RefSeq gene symbol
#'
#' Method (Day et al.)
#'   TSSR       TSS -50..+299 (L1 = 350 bp)
#'   gene body  TSS +300 .. TES +3000 (L2 = its length)
#'   reads      counted at their 5' end; reads with the same 5' end and strand in a replicate are
#'              collapsed into one (Day's duplicate removal; ENCODE already removed duplicates)
#'   density    reads / L / (library reads / 1e6) = rpm/bp; library = all collapsed reads
#'   signal     ChIP density - input density, in each region
#'   isoform    per gene, the annotated TSS with the strongest TSSR signal; among isoforms with
#'              that TSS, the longest. Day et al. fall back to H3K4me3 when no TSSR reaches
#'              0.001 rpm/bp; without GM12878 H3K4me3 ChIP/input here such genes get no PI, as do
#'              Day's genes that fail both
#'   PI         TSSR signal / gene-body signal; only when the TSSR signal is >= 0.001 rpm/bp and
#'              the gene-body signal is > 0
#'   paused     PI >= 2; no PI counts as not paused
#'   filters    near_other_gene: the gene (all isoforms) +/- 3 kb overlaps another RefSeq gene
#'              (either strand); short: chosen isoform TSS..TES < 1 kb
#'
#' analysis_set (PolII example genes): chosen isoform NM_ (protein-coding),
#'   chr1-22 (FiberHMM is autosomal), FPKM > 1, neither filter, has a PI. pi_bin = PI tertile
#'   within it.
#'
#' Output: LCL_project/pausing/pause_index.tsv, one row per RefSeq gene; gene_id = gene_name =
#'   RefSeq symbol; coordinates 1-based inclusive.
#' Run: sbatch code/pausing/cal_pausing_index/run_pause_index.sh

suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })

ANN      <- "/project/spott/cshan/annotations"
LCL_PUB  <- file.path(ANN, "LCL_public_data")
REFSEQ   <- file.path(ANN, "hg38.refGene.gtf.gz")
GTF      <- file.path(ANN, "gencode.v46.annotation.gtf.gz")
POL2     <- file.path(LCL_PUB, "Pol2_chip")
CHIP     <- file.path(POL2, paste0(c("ENCFF441GAS", "ENCFF306JCA", "ENCFF248WXH", "ENCFF590GEB",
                                     "ENCFF944QSC", "ENCFF811KVZ", "ENCFF603TJG"), ".bed.gz"))
INPUT    <- file.path(POL2, paste0(c("ENCFF863XTX", "ENCFF278NPZ"), ".bed.gz"))
EXPR     <- file.path(LCL_PUB, "gene_quan_ENCSR820PHH",
                      c("ENCFF873VWU.gene_quan.tsv", "ENCFF345SHY.gene_quan.tsv"))
OUT      <- "/project/spott/cshan/fiber-seq/LCL_project/pausing/pause_index.tsv"

CHROMS    <- paste0("chr", c(1:22, "X"))
AUTOSOMES <- paste0("chr", 1:22)
TSSR      <- c(-50L, 299L)
BODY_FROM <- 300L
BODY_PAST <- 3000L
PAD       <- 3000L                      # no other gene within this of the gene
MIN_LEN   <- 1000L                      # chosen isoform TSS..TES
MIN_DENS  <- 0.001                      # rpm/bp, background-subtracted TSSR
PAUSED_PI <- 2
EXPR_FPKM <- 1

stopifnot(all(file.exists(c(REFSEQ, GTF, CHIP, INPUT, EXPR))))

# --- RefSeq isoforms --------------------------------------------------------------
tx <- fread(cmd = paste("zcat", REFSEQ, "| awk -F'\t' '$3 == \"transcript\"' | cut -f1,4,5,7,9"),
            header = FALSE, col.names = c("chrom", "start", "end", "strand", "attr"), quote = "")
tx[, `:=`(gene = sub('.*gene_id "([^"]+)".*', "\\1", attr),
          tx_id = sub('.*transcript_id "([^"]+)".*', "\\1", attr), attr = NULL)]
tx <- tx[chrom %in% CHROMS]
all_genes <- tx[, .(start = min(start), end = max(end)), by = .(gene, chrom)]   # for the 3 kb filter
multi <- tx[, uniqueN(paste(chrom, strand)), by = gene][V1 > 1L, gene]
tx <- tx[!gene %in% multi]
tx[, `:=`(tss = fifelse(strand == "+", start, end), tes = fifelse(strand == "+", end, start),
          len = end - start + 1L)]
message(uniqueN(tx$gene), " RefSeq genes (", nrow(tx), " isoforms) on chr1-22, chrX; ",
        length(multi), " genes on several chromosomes / strands dropped")

# --- Pol II ChIP and input reads --------------------------------------------------
# 5' ends of one replicate, collapsed by position and strand
read_ends <- function(path) {
  r <- fread(cmd = paste("zcat", shQuote(path)), header = FALSE, select = c(1, 2, 3, 6),
             col.names = c("chrom", "start0", "end", "strand"))
  unique(r[, .(chrom, pos = fifelse(strand == "+", start0 + 1L, end), strand)])
}
# per-base 5'-end counts and library size of pooled replicates
read_library <- function(paths) {
  ends <- rbindlist(lapply(paths, read_ends))
  list(n = nrow(ends), cov = coverage(GRanges(ends$chrom, IRanges(ends$pos, width = 1L))))
}
# reads of a library in 1-based regions (coverage stops at the chromosome's last read, so regions
# beyond it have none)
region_reads <- function(lib, chrom, start, end) {
  out <- numeric(length(chrom))
  for (cc in intersect(unique(chrom), names(lib$cov))) {
    cov <- lib$cov[[cc]]
    i <- which(chrom == cc & start <= length(cov))
    out[i] <- viewSums(Views(cov, start[i], pmin(end[i], length(cov))))
  }
  out
}
rpm_bp <- function(reads, n, len) reads / len / (n / 1e6)
chip  <- read_library(CHIP)
input <- read_library(INPUT)
message("library reads (collapsed, pooled replicates): ChIP ", chip$n, ", input ", input$n)

# TSSR of an annotated TSS: list(start, end, chip, input, signal)
tssr <- function(chrom, tss, strand) {
  s <- pmax(1L, fifelse(strand == "+", tss + TSSR[1], tss - TSSR[2]))
  e <- fifelse(strand == "+", tss + TSSR[2], tss - TSSR[1])
  ch <- region_reads(chip, chrom, s, e); inp <- region_reads(input, chrom, s, e)
  list(start = s, end = e, chip = ch, input = inp,
       signal = rpm_bp(ch, chip$n, e - s + 1L) - rpm_bp(inp, input$n, e - s + 1L))
}

# --- optimal isoform: strongest TSSR, then the longest isoform with that TSS ------
cand <- unique(tx[, .(gene, chrom, strand, tss)])
cand[, tssr_dens := tssr(chrom, tss, strand)$signal]
setorder(cand, gene, -tssr_dens, tss)
best <- cand[!duplicated(gene)]
g <- merge(best[, .(gene, tss)], tx, by = c("gene", "tss"))
setorder(g, gene, -len, tx_id)
g <- g[!duplicated(gene)]
g <- merge(g, tx[, .(n_isoforms = .N, n_tss = uniqueN(tss)), by = gene], by = "gene")

# --- regions, densities, PI -------------------------------------------------------
ts <- tssr(g$chrom, g$tss, g$strand)
plus <- g$strand == "+"
g[, `:=`(tssr_start = ts$start, tssr_end = ts$end, chip_tssr = ts$chip, input_tssr = ts$input,
         tssr_dens = ts$signal,
         body_start = pmax(1L, fifelse(plus, tss + BODY_FROM, tes - BODY_PAST)),
         body_end   = fifelse(plus, tes + BODY_PAST, tss - BODY_FROM))]
g[, `:=`(chip_body = region_reads(chip, chrom, body_start, body_end),
         input_body = region_reads(input, chrom, body_start, body_end))]
g[, body_dens := rpm_bp(chip_body, chip$n, body_end - body_start + 1L) -
                 rpm_bp(input_body, input$n, body_end - body_start + 1L)]
g[, active := tssr_dens >= MIN_DENS]
g[, pi := fifelse(active & body_dens > 0, tssr_dens / body_dens, NA_real_)]
g[, paused := !is.na(pi) & pi >= PAUSED_PI]

# --- Day et al. gene filters ------------------------------------------------------
span <- tx[, .(chrom = chrom[1], start = min(start), end = max(end)), by = gene][match(g$gene, gene)]
hits <- findOverlaps(GRanges(span$chrom, IRanges(pmax(1L, span$start - PAD), span$end + PAD)),
                     GRanges(all_genes$chrom, IRanges(all_genes$start, all_genes$end)))
hits <- hits[all_genes$gene[subjectHits(hits)] != g$gene[queryHits(hits)]]
g[, near_other_gene := seq_len(.N) %in% queryHits(hits)]
g[, short := abs(tes - tss) + 1L < MIN_LEN]

# --- expression, analysis set, tertiles -------------------------------------------
gc <- fread(cmd = paste("zcat", GTF, "| awk -F'\t' '$3 == \"gene\"' | cut -f9"), header = FALSE,
            col.names = "attr", quote = "", sep = "\t")
gc <- unique(gc[, .(ensg = sub("\\..*", "", sub('.*gene_id "([^"]+)".*', "\\1", attr)),
                    gene = sub('.*gene_name "([^"]+)".*', "\\1", attr))])
gc <- gc[gene %in% gc[, .N, by = gene][N == 1L, gene]]          # symbols with one ENSG only
read_fpkm <- function(path) {
  d <- fread(path, select = c("gene_id", "FPKM"))
  d[grepl("^ENSG", gene_id), .(fpkm = max(FPKM)), by = .(ensg = sub("\\..*", "", gene_id))]
}
ex <- merge(read_fpkm(EXPR[1]), read_fpkm(EXPR[2]), by = "ensg", suffixes = c("_rep1", "_rep2"))
ex[, fpkm := (fpkm_rep1 + fpkm_rep2) / 2]
g <- merge(g, gc, by = "gene", all.x = TRUE)
g <- merge(g, ex, by = "ensg", all.x = TRUE)
g[, `:=`(gene_id = gene, gene_name = gene,
         gene_type = fifelse(startsWith(tx_id, "NM_"), "protein_coding", "non_coding"))]
g[, expressed := fpkm > EXPR_FPKM]
g[, analysis_set := gene_type == "protein_coding" & chrom %in% AUTOSOMES & (expressed %in% TRUE) &
      !near_other_gene & !short & !is.na(pi)]
cuts <- g[analysis_set == TRUE, quantile(pi, c(0, 1/3, 2/3, 1))]
g[analysis_set == TRUE, pi_bin := as.character(cut(pi, cuts, include.lowest = TRUE,
                                                   labels = c("low", "mid", "high")))]

out <- g[, .(gene_id, gene_name, ensg, gene_type, tx_id, n_isoforms, n_tss, chrom, tss, strand, tes,
             tssr_start, tssr_end, body_start, body_end, chip_tssr, chip_body, input_tssr, input_body,
             tssr_dens, body_dens, active, pi, paused, fpkm_rep1, fpkm_rep2, fpkm, expressed,
             near_other_gene, short, analysis_set, pi_bin)]
setorder(out, chrom, tss)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
fwrite(out, OUT, sep = "\t", na = "NA")

# --- summary ----------------------------------------------------------------------
message(nrow(out), " genes; ", out[n_tss > 1, .N], " with several annotated TSSs; ",
        out[, sum(is.na(fpkm))], " without an FPKM (symbol not matched to one GENCODE gene)")
message(out[, sum(near_other_gene)], " within 3 kb of another RefSeq gene, ", out[, sum(short)],
        " shorter than ", MIN_LEN, " bp, ", out[, sum(!near_other_gene & !short)], " pass both")
message(out[, sum(active)], " with TSSR >= ", MIN_DENS, " rpm/bp; ",
        out[, sum(active & body_dens <= 0)], " of them with gene-body signal <= 0 (no PI)")
pc <- out[gene_type == "protein_coding" & !near_other_gene & !short]
message("protein-coding, filters passed: ", nrow(pc), "; with a PI ", pc[, sum(!is.na(pi))], "; paused (PI >= ",
        PAUSED_PI, ") ", pc[, sum(paused)], " = ", round(100 * pc[, mean(paused)], 1), "% of all, ",
        round(100 * pc[!is.na(pi), mean(paused)], 1), "% of those with a PI")
a <- out[analysis_set == TRUE]
message("analysis set ", nrow(a), "; PI 5/25/50/75/95%: ",
        paste(signif(quantile(a$pi, c(.05, .25, .5, .75, .95)), 3), collapse = " / "),
        "; tertile cuts ", paste(signif(cuts[2:3], 3), collapse = ", "),
        "; Spearman PI vs FPKM ", round(cor(a$pi, a$fpkm, method = "spearman"), 3))
message("Wrote ", OUT)
