#!/usr/bin/env Rscript
# Run: bash code/pausing/PolII/run_gene_examples.sh code/pausing/PolII/04_gene_body_pol2.R
# Optional GENE_BODY_WORKERS (default 1). Outputs only TSV tables; no plots.
# Gene body: GENCODE v46 gene feature (introns included), 1-based inclusive.
# Cohort: all fibertools reads overlapping the gene, longest overlapping alignment
# per sample/read. Partial gene coverage is retained and explicitly reported.
# Pol2 candidate: unique FiberHMM footprint of 30 <= length <= 60 bp, fully
# contained in the gene AND selected read alignment. Size alone is not proof of Pol2.
# Peak fraction: candidates with >=1 bp overlap / all candidates; NA if none.
# A molecule overlapping two genes occurs once in each gene. Overall counts are
# therefore gene-molecule observations, not globally unique molecules/footprints.
# Existing functions are sourced unchanged. Only the input size-bin configuration
# is restricted to 30-60bp (inclusive endpoints; do not use assign_size_class()).

dyn.load("/project/spott/cshan/envs/dimelo/lib/libcurl.so.4", local = FALSE)
source("/project/spott/cshan/fiber-seq/code/pausing/PolII/pol2_pausing_functions.r")
setDTthreads(1L)
COMBINED_BINS <- "30-60bp"
BODY_OUT <- file.path(OUT_DIR, "gene_body_pol2/tables")
GTF <- file.path(ANN, "gencode.v46.annotation.gtf.gz")
PEAK_FILE <- file.path(LCL_PUB, "Pol2_chip/ENCFF886PSD.bed.gz")
safe_fraction <- function(n, d) ifelse(d > 0, n / d, NA_real_)

body_annotation <- function() {
  a <- fread(cmd = sprintf("gzip -dc %s | awk -F'\\t' '$3 == \"gene\"'", shQuote(GTF)),
             header = FALSE, sep = "\t", quote = "")
  a <- a[, .(gene_id = sub("\\.[0-9]+$", "", sub('.*gene_id "([^"]+)".*', '\\1', V9)),
             chrom = V1, gene_start = V4, gene_end = V5, strand = V7)]
  g <- merge(curated_example_genes(), a, by = c("gene_id", "chrom"), all.x = TRUE)
  stopifnot(nrow(g) == nrow(curated_example_genes()), !anyNA(g$gene_start),
            !anyDuplicated(g$gene_id), all(g$chrom %in% AUTOSOMES))
  g[, gene_length_bp := gene_end - gene_start + 1L]
  g
}

count_gene_body <- function(g, samples, peaks) {
  message("Counting ", g$gene_name)
  region <- GRanges(g$chrom, IRanges(g$gene_start, g$gene_end))
  reads <- rbindlist(lapply(samples$sample_name, function(s) {
    d <- as.data.table(read_ft_bed12(nuc_file(s, g$chrom), region, longest_alignment = TRUE))
    if (!nrow(d)) return(NULL)
    d[, .(sample_name = s, original_RID = RID, RID = paste(s, RID, sep = "::"),
          read_start = start + 1L, read_end = end)]
  }))
  if (!nrow(reads)) stop("No overlapping reads: ", g$gene_name)
  if (anyDuplicated(reads$original_RID)) stop("Ambiguous unqualified read IDs: ", g$gene_name)
  reads[, `:=`(gene_name = g$gene_name, group = g$group, gene_id = g$gene_id,
                chrom = g$chrom, gene_start = g$gene_start, gene_end = g$gene_end,
                strand = g$strand, gene_length_bp = g$gene_length_bp)]
  reads[, covered_gene_bp := pmin(read_end, gene_end) - pmax(read_start, gene_start) + 1L]
  reads[, `:=`(full_gene_span = read_start <= gene_start & read_end >= gene_end,
                n_pol2 = 0L, pol2_length_sum_bp = 0L, n_pol2_chip_overlap = 0L)]
  files <- file.path(COMBINED_ROOT, g$chrom,
                     sprintf("combined_%s_%s_fps.bed.gz", g$chrom, COMBINED_BINS))
  stopifnot(all(file.exists(files)))
  fp <- window_tf_footprints(g$chrom, g$gene_start, g$gene_end)
  if (!ncol(fp)) fp <- data.table(RID = character(), start0 = integer(), end0 = integer())
  fp <- fp[end0 - start0 >= 30L & end0 - start0 <= 60L &
             start0 >= g$gene_start - 1L & end0 <= g$gene_end]
  n_raw <- nrow(fp)
  fp <- merge(fp, reads[, .(original_RID, qualified_RID = RID, read_start, read_end)],
              by.x = "RID", by.y = "original_RID")
  fp <- fp[start0 >= read_start - 1L & end0 <= read_end]
  if (nrow(fp)) {
    fp[, chip_overlap := countOverlaps(GRanges(g$chrom, IRanges(start0 + 1L, end0)),
                                      peaks, ignore.strand = TRUE) > 0L]
    counts <- fp[, .(n_pol2 = .N, pol2_length_sum_bp = sum(end0 - start0),
                     n_pol2_chip_overlap = sum(chip_overlap)), by = .(RID = qualified_RID)]
    reads[counts, on = "RID", `:=`(n_pol2 = i.n_pol2,
          pol2_length_sum_bp = i.pol2_length_sum_bp, n_pol2_chip_overlap = i.n_pol2_chip_overlap)]
  }
  reads[, `:=`(fraction_pol2_chip_overlap = safe_fraction(n_pol2_chip_overlap, n_pol2),
                pol2_per_kb_covered = 1000 * n_pol2 / covered_gene_bp)]
  stopifnot(all(reads$covered_gene_bp > 0L), !anyDuplicated(reads$RID),
            all(reads$n_pol2_chip_overlap <= reads$n_pol2))
  summary <- cbind(g, summarize_body(reads))
  summary[, `:=`(n_body_footprints_before_read_filter = n_raw,
                  n_body_footprints_excluded_by_read_filter = n_raw - sum(reads$n_pol2))]
  fwrite(reads, file.path(BODY_OUT, "per_gene", paste0(g$gene_name, ".tsv.gz")), sep = "\t", na = "NA")
  list(reads = reads, summary = summary)
}

summarize_body <- function(d) {
  d[, .(n_gene_molecules = .N, n_full_gene_span = sum(full_gene_span),
         n_with_pol2 = sum(n_pol2 > 0L), fraction_with_pol2 = mean(n_pol2 > 0L),
         n_pol2 = sum(n_pol2), n_pol2_chip_overlap = sum(n_pol2_chip_overlap),
         fraction_pol2_chip_overlap = safe_fraction(sum(n_pol2_chip_overlap), sum(n_pol2)),
         mean_pol2_per_molecule = mean(n_pol2), median_pol2_per_molecule = median(n_pol2),
         covered_gene_bp = sum(as.numeric(covered_gene_bp)),
         pol2_per_kb_covered = 1000 * sum(n_pol2) / sum(as.numeric(covered_gene_bp)))]
}

main <- function() {
  dir.create(file.path(BODY_OUT, "per_gene"), recursive = TRUE, showWarnings = FALSE)
  genes <- body_annotation()
  samples <- fread(SAMPLES)
  # narrowPeak has floating-point signal fields; only its BED coordinates matter.
  peak_bed <- fread(PEAK_FILE, header = FALSE, select = 1:3)
  stopifnot(all(peak_bed$V2 >= 0L), all(peak_bed$V3 > peak_bed$V2))
  peaks <- GRanges(peak_bed$V1, IRanges(peak_bed$V2 + 1L, peak_bed$V3))
  workers <- as.integer(Sys.getenv("GENE_BODY_WORKERS", "1"))
  stopifnot(workers >= 1L, nrow(samples) == 31L, length(peaks) > 0L)
  fwrite(genes, file.path(BODY_OUT, "gene_annotation.tsv"), sep = "\t")
  result <- parallel::mclapply(seq_len(nrow(genes)), function(i)
    count_gene_body(genes[i], samples, peaks), mc.cores = workers, mc.preschedule = FALSE)
  if (any(vapply(result, inherits, logical(1), "try-error"))) stop("Gene worker failed; see log")
  reads <- rbindlist(lapply(result, `[[`, "reads"))
  summary <- rbindlist(lapply(result, `[[`, "summary"))
  groups <- rbindlist(lapply(unique(reads$group), function(k)
    cbind(data.table(group = k), summarize_body(reads[group == k]))))
  overall <- summarize_body(reads)
  for (nm in c("molecules", "genes", "groups", "overall")) {
    x <- switch(nm, molecules = reads, genes = summary, groups = groups, overall = overall)
    fwrite(x, file.path(BODY_OUT, paste0(nm, if (nm == "molecules") ".tsv.gz" else ".tsv")),
           sep = "\t", na = "NA")
  }
  fwrite(data.table(field = c("gene_annotation", "peak_file", "footprint_root", "footprint_bins",
                            "footprint_rule", "peak_overlap_rule", "molecule_cohort"),
     value = c(GTF, PEAK_FILE, COMBINED_ROOT, COMBINED_BINS,
               "30-60 bp inclusive; fully inside gene and chosen alignment; exact duplicates removed",
               ">=1 bp; each footprint counted at most once regardless of number of peaks",
               "All gene-overlapping reads; longest overlapping alignment per sample/read; zeros retained")),
     file.path(BODY_OUT, "provenance.tsv"), sep = "\t")
  print(groups)
  print(overall)
  message("Completed ", nrow(summary), " genes; tables: ", BODY_OUT)
}
if (sys.nframe() == 0L) main()
