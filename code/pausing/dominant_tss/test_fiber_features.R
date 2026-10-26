#!/usr/bin/env Rscript
# End-to-end R fixtures for 03_fiber_features.R. No new functions are defined.
# Uses Rsamtools to build indexed BAM/BED inputs; no reference FASTA exists.
# Run with the R library environment documented by run_fiber_features.sh.
suppressPackageStartupMessages({
  library(data.table)
  library(Rsamtools)
})
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
root <- tempfile("dominant-fibers-test-")
dir.create(root)
message("Fixture directory: ", root)
genes <- data.table(gene_id = c("plus", "minus"), chrom = "chr22", strand = c("+", "-"),
                    anchor_tss = c(4001L, 10001L), anchor_tes = c(14001L, 1L), pause_peak_rel = 40)
gene_path <- file.path(root, "genes.tsv")
fwrite(genes, gene_path, sep = "\t")
all_fp <- list(`10-30` = list(), `40-60` = list(), `60-80` = list())
sample_rows <- list()

for (sample in c("s1", "s2")) {
  fire_dir <- file.path(root, sample)
  extracted <- file.path(fire_dir, "extracted_results")
  dir.create(extracted, recursive = TRUE)
  sam_records <- m6a <- cpg <- segments <- list()
  for (gi in seq_len(nrow(genes))) {
    gene <- genes[gi]; tss <- gene$anchor_tss
    orient <- if (gene$strand == "+") 1L else -1L
    labels <- if (sample == "s1") c("good", "empty", "amb", "short", "secondary") else c("amb", "pool")
    for (label in labels) {
      rid <- paste(gene$gene_id, label, sep = "_")
      st <- if (orient == 1L) 2000L else 7000L
      en <- if (orient == 1L) 7000L else 12000L
      if (label == "short" || (sample == "s2" && label == "amb")) { st <- tss - 100L; en <- tss + 100L }
      spans <- data.table(start = st, end = en, flag = if (label == "secondary") 256L else 0L)
      if (label %in% c("good", "short") && sample == "s1")
        spans <- rbind(spans, data.table(start = tss - 2500L, end = tss + 2500L, flag = 2048L))
      for (ai in seq_len(nrow(spans))) {
        st <- spans$start[ai]; en <- spans$end[ai]; flag <- spans$flag[ai]
        sam_records[[length(sam_records) + 1L]] <- paste(c(rid, flag, "chr22", st + 1L, 60,
          paste0(en - st, "M"), "*", 0, 0, strrep("A", en - st), "*"), collapse = "\t")
        for (kind in c("m6a", "cpg")) {
          calls <- if (label == "empty") integer() else if (kind == "m6a") c(tss - 1L, tss + 1L) else c(tss + 99L, tss + 100L)
          if (label == "pool") calls <- if (kind == "m6a") c(tss + 3L, tss + 5L) else c(tss + 199L, tss + 200L)
          # Rejected alignments have distinct sites that must not enter the union.
          if (label %in% c("short", "secondary") || flag == 2048L || (sample == "s2" && label == "amb"))
            calls <- c(tss + 10L, tss + 12L)
          if (label != "empty") calls <- c(calls, tss + 2200L) # outside analysis windows
          positions <- sort(unique(c(st, en - 1L, calls[calls > st & calls < en - 1L])))
          record <- as.list(c("chr22", st, en, rid, 0, "+", st, en, "0,0,0", length(positions),
                              paste(c(0L, rep(1L, length(positions) - 1L)), collapse = ","), paste(positions - st, collapse = ",")))
          if (kind == "m6a") m6a[[length(m6a) + 1L]] <- record else cpg[[length(cpg) + 1L]] <- record
        }
      }
      if (sample == "s2" || !label %in% c("good", "amb")) next
      seg <- data.table(left = c(-1000, -200, -100, -10, 11, 30, 210, 390, 570),
        right = c(-850, -101, -11, 10, 29, 119, 309, 489, 669),
        rgb = c("169,169,169", "169,169,169", "147,112,219", "255,0,0", "147,112,219",
                rep("169,169,169", 4)))
      for (i in seq_len(nrow(seg))) {
        bounds <- sort(tss + orient * c(seg$left[i], seg$right[i]))
        segments[[length(segments) + 1L]] <- as.list(c("chr22", bounds[1] - 1L, bounds[2], rid,
          0, "+", bounds[1] - 1L, bounds[2], seg$rgb[i], 0.8, 1))
      }
      fps <- data.table(size = c("60-80", "40-60", "40-60", "10-30"),
                        left = c(-30, 20, 350, -100), right = c(29, 59, 389, -81))
      for (i in seq_len(nrow(fps))) {
        bounds <- sort(tss + orient * c(fps$left[i], fps$right[i]))
        record <- as.list(c("chr22", bounds[1] - 1L, bounds[2], rid))
        size <- fps$size[i]
        all_fp[[size]] <- c(all_fp[[size]], list(record, record)) # duplicate pooled call
      }
    }
  }
  sam_path <- file.path(fire_dir, "reads.sam")
  writeLines(c("@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:chr22\tLN:16000", unlist(sam_records)), sam_path)
  unsorted <- asBam(sam_path, destination = file.path(fire_dir, "unsorted"), overwrite = TRUE, indexDestination = FALSE)
  bam_path <- sortBam(unsorted, destination = file.path(fire_dir, "reads"))
  indexBam(bam_path)
  for (kind in c("m6a", "cpg", "fire")) {
    records <- switch(kind, m6a = m6a, cpg = cpg, fire = segments)
    path <- if (kind == "fire") file.path(extracted, paste0(sample, ".fire_all.bed")) else
      file.path(extracted, paste0(kind, "_by_chr"), paste0(sample, ".ft_extracted_", kind, ".chr22.bed"))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (length(records)) {
      tab <- rbindlist(records)
      tab[, `:=`(V2 = as.integer(V2), V3 = as.integer(V3))]
      setorder(tab, V1, V2, V3)
      fwrite(tab, path, sep = "\t", col.names = FALSE, quote = FALSE)
    } else file.create(path)
    bgzip(path, dest = paste0(path, ".gz"), overwrite = TRUE)
    indexTabix(paste0(path, ".gz"), format = "bed")
  }
  sample_rows[[length(sample_rows) + 1L]] <- data.table(sample_name = sample, fire_dir = fire_dir, bam_file = bam_path)
}
manifest <- file.path(root, "samples.tsv")
fwrite(rbindlist(sample_rows), manifest, sep = "\t")
fp_root <- file.path(root, "fp")
dir.create(file.path(fp_root, "chr22"), recursive = TRUE)
for (size in names(all_fp)) {
  path <- file.path(fp_root, "chr22", paste0("combined_chr22_", size, "bp_fps.bed"))
  tab <- rbindlist(all_fp[[size]])
  tab[, `:=`(V2 = as.integer(V2), V3 = as.integer(V3))]
  setorder(tab, V1, V2, V3)
  fwrite(tab, path, sep = "\t", col.names = FALSE, quote = FALSE)
  bgzip(path, dest = paste0(path, ".gz"), overwrite = TRUE)
  indexTabix(paste0(path, ".gz"), format = "bed")
}
args <- c(file.path(script_dir, "run_fiber_features.sh"), "--chrom", "chr22", "--genes", gene_path,
          "--samples", manifest, "--fp-root", fp_root, "--out-dir", file.path(root, "out"))
stopifnot(system2("bash", shQuote(args)) == 0L)
output <- file.path(root, "out", "fiber_features.chr22.tsv.gz")
con <- gzfile(output, "rt"); rows <- read.delim(con, check.names = FALSE); close(con)
stopifnot(nrow(rows) == 8L, !any(c("at_n", "cpg_n") %in% names(rows)),
          all(rows$m6a_site_n == 4L), all(rows$cpg_site_n == 4L))
good <- rows[endsWith(rows$rid, "good"), ]
stopifnot(nrow(good) == 2L, all(good$promoter_state == "open_fire"))
expected <- c(plus1_start = 30, plus1_size = 90, plus2_start = 210, minus1_end = -101,
  ndr_len = 130, tss_msp_len = 21, down_n_nuc = 3, down_dyad_spacing = 180,
  pic = 1, ppp = 1, tf_up_n = 1, elong_n = 1, m6a_n = 2, mcpg_n = 2,
  m6a_site_n = 4, cpg_site_n = 4)
for (field in names(expected)) stopifnot(all(good[[field]] == expected[[field]]))
stopifnot(all(abs(good$elong_per_kb - 1000 / 1701) < 1e-10),
          all(good$read_start == ifelse(good$gene_id == "plus", 2000L, 7000L)))
for (field in setdiff(names(good), c("gene_id", "rid", "read_start", "read_end")))
  stopifnot(isTRUE(all.equal(good[[field]][1], good[[field]][2])))
amb <- rows[endsWith(rows$rid, "amb"), ]
stopifnot(all(amb$fp_ambiguous == 1), all(is.na(amb$pic)), all(is.na(amb$elong_n)))
empty <- rows[endsWith(rows$rid, "empty"), ]
stopifnot(all(empty$m6a_n == 0), all(empty$mcpg_n == 0), all(empty$promoter_state == "other"),
          all(is.na(empty$frac_nuc)))
qc <- fread(file.path(root, "out", "fiber_qc.chr22.tsv"))
stopifnot(qc$reads_overlapping[1] == 8L, qc$reads_spanning_window[1] == 6L,
          qc$reads_lacking_fire_all[1] == 2L, qc$rows_written[1] == 6L,
          qc$reads_spanning_window[2] == 2L, qc$rows_written[2] == 2L)

# Unindexed pooled input must produce exactly the same features.
unlink(list.files(file.path(fp_root, "chr22"), pattern = "\\.tbi$", full.names = TRUE))
stopifnot(system2("bash", shQuote(args)) == 0L)
con <- gzfile(output, "rt"); streamed <- read.delim(con, check.names = FALSE); close(con)
stopifnot(identical(rows, streamed))

# Empty selection still produces a valid schema and per-sample zero QC.
stopifnot(system2("bash", shQuote(c(args, "--gene-ids", "absent"))) == 0L)
con <- gzfile(output, "rt"); no_rows <- read.delim(con, check.names = FALSE); close(con)
stopifnot(nrow(no_rows) == 0L, identical(names(no_rows), names(rows)))
qc <- fread(file.path(root, "out", "fiber_qc.chr22.tsv"))
stopifnot(all(qc$rows_written == 0L), all(qc$reads_overlapping == 0L))

# All eligible reads have zero calls: retain them with zero observed-site counts.
for (path in list.files(root, pattern = "ft_extracted_.*\\.bed$", recursive = TRUE, full.names = TRUE)) {
  tab <- fread(path, header = FALSE)
  tab[, `:=`(V10 = 2L, V11 = "0,1", V12 = paste0("0,", V3 - V2 - 1L))]
  fwrite(tab, path, sep = "\t", col.names = FALSE, quote = FALSE)
  bgzip(path, dest = paste0(path, ".gz"), overwrite = TRUE)
  unlink(paste0(path, ".gz.tbi"))
  indexTabix(paste0(path, ".gz"), format = "bed")
}
stopifnot(system2("bash", shQuote(c(args, "--gene-ids", "plus"))) == 0L)
con <- gzfile(output, "rt"); zero_calls <- read.delim(con, check.names = FALSE); close(con)
stopifnot(nrow(zero_calls) == 4L, all(zero_calls$m6a_n == 0L), all(zero_calls$mcpg_n == 0L),
          all(zero_calls$m6a_site_n == 0L), all(zero_calls$cpg_site_n == 0L))
unlink(root, recursive = TRUE)
message("PASS: no FASTA, pooled observed sites, primary-only coverage, mirrored strands, sentinels, local ambiguity, QC, streaming, empty outputs, zero-call genes")
