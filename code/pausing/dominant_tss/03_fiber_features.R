#!/usr/bin/env Rscript
# Per-fiber promoter features at dominant TSSs; one chromosome per invocation.
# Input genes: gene_id, chrom, strand, anchor_tss, anchor_tes, pause_peak_rel.
# TSS/TES are 1-based; relative intervals are inclusive and gene-strand oriented.
# BED read_start/read_end remain 0-based half-open. Read strand never orients features.
# Eligibility: mapped primary BAM alignment (neither secondary nor supplementary)
# spans the entire [-1000,+2000] window, and an m6A BED12 record matches its span.
# Coverage means reference alignment span, including CIGAR D/N; no MAPQ cutoff.
# No longest-alignment selection. All other features are restricted to that span.
# BED12 first/last sentinel blocks are removed by the existing shared parser;
# the first sentinel can have size zero. Validate 1-bp calls AFTER removal.
# Nucleosomes: grey segments of original length >= --nuc-min (90); sizes touching
# a window edge are NA. Occupancy/fractions use unions; frac_nuc includes all sizes.
# ACF: mean-centred binary qualifying-nucleosome occupancy, biased FFT estimate
# normalized at lag zero, over +150..+2000; strongest positive local max 140..250.
# Methylation: same observed-position definition as assemble_region_m6a() and
# the topic-model/Leiden matrices. m6a_site_n is the union of m6A call positions
# in the whole window across eligible fibers in all manifest samples; cpg_site_n
# is the union of CpG call positions within TSS +/-500 on those fibers. Counts
# are unique called coordinates per fiber; no reference-base filtering or C/G
# collapsing. These are observed-site counts, NOT all possible A/T or CpG sites.
# No FASTA is used. BAM header lengths provide chromosome-boundary checks.
# Pooled footprints: identical records deduplicated. fp_ambiguous is LOCAL: same
# RID in multiple manifest samples' m6A alignments overlapping this gene window,
# including non-primary/non-spanning records. No genome-wide audit is performed.
# A missing segmentation gives state 'other' and NA segmentation measurements.
# No new functions are defined: source the repository parsers and use inline code.
# Example: Rscript 03_fiber_features.R --chrom chr22 --gene-ids ENSG00000100320

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(Rsamtools)
  library(GenomicAlignments)
})

opts <- list(
  genes = "/project/spott/cshan/fiber-seq/LCL_project/pausing/dominant_tss/pause_index_dominant.tsv",
  samples = "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv",
  `fp-root` = "/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged/combined/joint_trained_tracks",
  `out-dir` = "/project/spott/cshan/fiber-seq/LCL_project/pausing/dominant_tss",
  parser = "/project/spott/cshan/fiber-seq/code/parsing_functions/parsing_footprints_functions.r",
  chrom = "", `window-start` = "-1000", `window-end` = "2000", `nuc-min` = "90",
  `gene-ids` = "")
argv <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% argv) {
  cat("Usage: Rscript 03_fiber_features.R --chrom chrN [--option value]\n",
      "Genes require gene_id, chrom, strand, anchor_tss, anchor_tes, pause_peak_rel.\n",
      "Samples require sample_name, fire_dir, bam_file (CSV or TSV).\n",
      "Options/defaults:\n", sep = "")
  print(opts)
  quit(status = 0L)
}
if (length(argv) %% 2L) stop("Options must be --name value pairs; see --help")
for (i in (seq_len(length(argv) %/% 2L) * 2L - 1L)) {
  key <- sub("^--", "", argv[i])
  if (!startsWith(argv[i], "--") || !key %in% names(opts)) stop("Unknown option: ", argv[i])
  opts[[key]] <- argv[i + 1L]
}
chrom <- opts$chrom
lo <- as.integer(opts[["window-start"]]); hi <- as.integer(opts[["window-end"]])
nuc_min <- as.integer(opts[["nuc-min"]])
stopifnot(chrom %in% paste0("chr", c(1:22, "X")), !anyNA(c(lo, hi, nuc_min)),
          lo <= -500L, hi >= 500L, nuc_min > 0L)
source(opts$parser)
genes <- fread(opts$genes)
samples <- fread(opts$samples)
stopifnot(all(c("gene_id", "chrom", "strand", "anchor_tss", "anchor_tes", "pause_peak_rel") %in% names(genes)),
          all(c("sample_name", "fire_dir", "bam_file") %in% names(samples)),
          nrow(samples) > 0L, !anyDuplicated(samples$sample_name))
target_chrom <- chrom
genes <- genes[chrom == target_chrom]
if (nzchar(opts[["gene-ids"]])) genes <- genes[gene_id %in% strsplit(opts[["gene-ids"]], ",", fixed = TRUE)[[1]]]
stopifnot(!anyDuplicated(genes$gene_id), !anyNA(genes[, .(gene_id, strand, anchor_tss, anchor_tes)]),
          all(genes$strand %in% c("+", "-")), all(genes$anchor_tss %% 1 == 0),
          all(genes$anchor_tes %% 1 == 0))
setorder(genes, anchor_tss, gene_id)
chr_len <- unname(scanBamHeader(samples$bam_file[1])[[1]]$targets[chrom])
stopifnot(length(chr_len) == 1L, !is.na(chr_len), chr_len > 0L)
genes[, window_start := anchor_tss + fifelse(strand == "+", lo, -hi)]
genes[, window_end := anchor_tss + fifelse(strand == "+", hi, -lo)]
message("Skipping ", genes[window_start < 1L | window_end > chr_len, .N], " chromosome-edge genes")
genes <- genes[window_start >= 1L & window_end <= chr_len]
samples[, m6a_path := file.path(fire_dir, "extracted_results", "m6a_by_chr",
                               paste0(sample_name, ".ft_extracted_m6a.", chrom, ".bed.gz"))]
samples[, cpg_path := file.path(fire_dir, "extracted_results", "cpg_by_chr",
                              paste0(sample_name, ".ft_extracted_cpg.", chrom, ".bed.gz"))]
samples[, fire_path := file.path(fire_dir, "extracted_results", paste0(sample_name, ".fire_all.bed.gz"))]
fp_paths <- file.path(opts[["fp-root"]], chrom, paste0("combined_", chrom, "_", c("10-30", "40-60", "60-80"), "bp_fps.bed.gz"))
for (p in c(samples$m6a_path, samples$cpg_path, samples$fire_path)) {
  if (!file.exists(p) || !file.exists(paste0(p, ".tbi"))) stop("Missing BED or tabix index: ", p)
}
for (si in seq_len(nrow(samples))) {
  for (kind in c("m6a", "cpg", "fire")) {
    present <- chrom %in% seqnamesTabix(samples[[paste0(kind, "_path")]][si])
    set(samples, i = si, j = paste0(kind, "_has_chrom"), value = present)
  }
}
if (!all(file.exists(fp_paths))) stop("Missing pooled footprint BED: ", paste(fp_paths[!file.exists(fp_paths)], collapse = ", "))
for (p in samples$bam_file) {
  if (!file.exists(p) || !any(file.exists(c(paste0(p, ".bai"), sub("\\.bam$", ".bai", p)))))
    stop("Missing BAM or BAI: ", p)
  sample_chr_len <- unname(scanBamHeader(p)[[1]]$targets[chrom])
  if (is.na(sample_chr_len) || sample_chr_len != chr_len) stop("Inconsistent chromosome length in BAM: ", p)
}

# Typed prototype also guarantees headers when no fibers qualify.
prototype <- data.table(gene_id = character(), sample_name = character(), rid = character(),
  read_start = integer(), read_end = integer(), hp = character(), ppp_window_source = character(),
  promoter_state = character(), tss_fire_score = double(), tss_msp_len = double(),
  tss_nuc_start = double(), tss_nuc_size = double(),
  plus1_start = double(), plus1_size = double(), plus2_start = double(), plus2_size = double(),
  plus3_start = double(), plus3_size = double(), minus1_start = double(), minus1_end = double(),
  minus1_size = double(), ndr_len = double(), down_n_nuc = integer(), down_dyad_spacing = double(),
  down_acf_lag = integer(), down_acf_peak = double(), pic = integer(), ppp = integer(),
  tf_up_n = integer(), elong_n = integer(), elong_per_kb = double(), fp_ambiguous = integer(),
  m6a_n = integer(), m6a_site_n = integer(), frac_nuc = double(), frac_msp = double(),
  mcpg_n = integer(), cpg_site_n = integer())
dir.create(opts[["out-dir"]], recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(opts[["out-dir"]], paste0("fiber_features.", chrom, ".tsv.gz"))
tmp_path <- paste0(out_path, ".", Sys.getpid(), ".tmp")
out <- gzfile(tmp_path, "wt")
write.table(prototype, out, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
qc <- samples[, .(sample_name, reads_overlapping = 0L, reads_spanning_window = 0L,
                  reads_lacking_fire_all = 0L, rows_written = 0L)]
seen_overlap <- seen_span <- seen_no_fire <- vector("list", nrow(samples))

# The supplied pooled tracks are ordinary gzip without tabix indexes. Stream
# each once and retain only footprints overlapping the union of gene windows.
# Indexed overrides can instead use the repository's tabix reader per window.
fp_cache <- vector("list", length(fp_paths))
windows <- reduce(GRanges(rep(chrom, nrow(genes)), IRanges(genes$window_start, genes$window_end)))
for (fi in seq_along(fp_paths)) {
  if (file.exists(paste0(fp_paths[fi], ".tbi")) || !nrow(genes)) next
  message("Streaming unindexed footprints: ", fp_paths[fi])
  con <- gzfile(fp_paths[fi], "rt")
  chunks <- list()
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    if (!length(lines)) next
    z <- fread(text = paste0(paste(lines, collapse = "\n"), "\n"), header = FALSE, sep = "\t", quote = "")
    if (ncol(z) != 4L) footprint_column_error(fp_paths[fi], "bed4_pooled", ncol(z))
    setnames(z, c("chr", "start", "end", "RID"))
    z <- z[chr == target_chrom]
    if (!nrow(z)) next
    keep <- overlapsAny(GRanges(z$chr, IRanges(z$start + 1L, z$end)), windows)
    if (any(keep)) chunks[[length(chunks) + 1L]] <- z[keep]
  }
  close(con)
  fp_cache[[fi]] <- unique(rbindlist(chunks))
}

for (gi in seq_len(nrow(genes))) {
  g <- genes[gi]
  tss <- g$anchor_tss; orient <- if (g$strand == "+") 1L else -1L
  ws <- g$window_start; we <- g$window_end
  region <- GRanges(chrom, IRanges(ws, we))
  observed_m6a <- observed_cpg <- integer()
  gene_rows <- vector("list", nrow(samples))
  pause <- g$pause_peak_rel
  ppp_lo <- if (is.finite(pause)) pause - 25 else 20
  ppp_hi <- if (is.finite(pause)) pause + 25 else 80
  elong_end <- min(2000, orient * (g$anchor_tes - tss))
  elong_len <- max(0, elong_end - 300 + 1)
  # Exclude the pause window from both elongation numerator and denominator.
  elong_len <- elong_len - max(0, floor(min(elong_end, ppp_hi)) - ceiling(max(300, ppp_lo)) + 1)

  beds <- vector("list", nrow(samples)); owners <- list()
  for (si in seq_len(nrow(samples))) {
    beds[[si]] <- if (samples$m6a_has_chrom[si]) unique(as.data.table(read_ft_bed12(samples$m6a_path[si], region,
                                                   longest_alignment = FALSE))) else data.table()
    if (nrow(beds[[si]])) owners[[si]] <- unique(beds[[si]][, .(RID, sample_name = samples$sample_name[si])])
  }
  owner_table <- rbindlist(owners)
  ambiguous <- if (nrow(owner_table)) owner_table[, .N, by = RID][N > 1L, RID] else character()
  footprints <- list()
  for (fi in seq_along(fp_paths)) {
    if (file.exists(paste0(fp_paths[fi], ".tbi"))) {
      z <- if (chrom %in% seqnamesTabix(fp_paths[fi])) as.data.table(read_tabix_region(fp_paths[fi], region)) else data.table()
    } else {
      cached <- fp_cache[[fi]]
      z <- if (nrow(cached)) copy(cached[start < we & end >= ws]) else data.table()
    }
    if (!nrow(z)) next
    if (ncol(z) != 4L) footprint_column_error(fp_paths[fi], "bed4_pooled", ncol(z))
    setnames(z, c("chr", "start", "end", "RID"))
    z[, size_class := c("10-30", "40-60", "60-80")[fi]]
    footprints[[fi]] <- z
  }
  fp <- unique(rbindlist(footprints))

  for (si in seq_len(nrow(samples))) {
    s <- samples[si]; bed <- beds[[si]]
    # Query primary records directly; supplementary m6A records cannot establish eligibility.
    bam <- scanBam(s$bam_file, param = ScanBamParam(which = region,
      flag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
                         isSupplementaryAlignment = FALSE), what = c("qname", "pos", "cigar")))[[1]]
    primary <- unique(data.table(RID = bam$qname, start = bam$pos - 1L,
      end = bam$pos - 1L + cigarWidthAlongReferenceSpace(bam$cigar)))
    if (anyDuplicated(primary$RID)) stop("Multiple primary alignments for a read in ", s$sample_name)
    seen_overlap[[si]] <- union(seen_overlap[[si]], primary$RID)
    primary <- primary[start <= ws - 1L & end >= we]
    seen_span[[si]] <- union(seen_span[[si]], primary$RID)
    if (!nrow(primary) || !nrow(bed)) next
    bed <- bed[primary, on = .(RID, start, end), nomatch = 0L]
    if (!nrow(bed)) next
    if (anyDuplicated(bed$RID)) stop("Conflicting m6A records at a primary span: ", s$sample_name)
    mods <- as.data.table(convert_ft_bed12_to_bed6(as.data.frame(bed),
      format = "bed12_fibertools", source = s$m6a_path))
    if (nrow(mods) && any(mods$end - mods$start != 1L)) stop("Non-1-bp m6A call after sentinel removal: ", s$m6a_path)
    if (nrow(mods)) {
      mods <- mods[end >= ws & end <= we]
      observed_m6a <- union(observed_m6a, mods$end)
    }
    cgbed <- if (s$cpg_has_chrom) unique(as.data.table(read_ft_bed12(s$cpg_path, region, longest_alignment = FALSE))) else data.table()
    cg <- data.table()
    if (nrow(cgbed)) {
      # CpG denominator uses the same emitted fibers as m6A, not extra CpG-only reads.
      cgbed <- cgbed[bed[, .(RID, start, end)], on = .(RID, start, end), nomatch = 0L]
      if (anyDuplicated(cgbed$RID)) stop("Conflicting CpG records at a primary span: ", s$sample_name)
      if (nrow(cgbed)) cg <- as.data.table(convert_ft_bed12_to_bed6(as.data.frame(cgbed),
        format = "bed12_fibertools", source = s$cpg_path))
      if (nrow(cg) && any(cg$end - cg$start != 1L)) stop("Non-1-bp CpG call after sentinel removal: ", s$cpg_path)
      if (nrow(cg)) {
        cg <- cg[end >= tss - 500L & end <= tss + 500L]
        observed_cpg <- union(observed_cpg, cg$end)
      }
    }
    fire <- if (s$fire_has_chrom) unique(as.data.table(read_tabix_region(s$fire_path, region))) else data.table()
    if (nrow(fire)) {
      if (ncol(fire) != 11L) stop("Expected 11 FIRE columns: ", s$fire_path)
      setnames(fire, FIRE_BED_COLS)
      fire[, class := fifelse(rgb == RGB_NUCLEOSOME, "nucleosome", fifelse(rgb == RGB_LINKER, "linker", "FIRE"))]
    }
    rows <- vector("list", nrow(bed))
    for (ri in seq_len(nrow(bed))) {
      b <- bed[ri]; rid_value <- b$RID
      row <- prototype[NA_integer_]
      row[, `:=`(gene_id = g$gene_id, sample_name = s$sample_name, rid = rid_value,
        read_start = b$start, read_end = b$end,
        ppp_window_source = if (is.finite(pause)) "gro_peak" else "fixed",
        promoter_state = "other",
        fp_ambiguous = as.integer(rid_value %in% ambiguous))]
      p <- if (nrow(mods)) unique(mods[RID == rid_value, start + 1L]) else integer()
      row[, m6a_n := length(p)]
      # Missing CpG record is unknown; a sentinel-only matching record means zero calls.
      if (nrow(cgbed) && rid_value %in% cgbed$RID) {
        p <- if (nrow(cg)) unique(cg[RID == rid_value, start + 1L]) else integer()
        row[, mcpg_n := length(p)]
      }
      seg <- if (nrow(fire)) copy(fire[RID == rid_value & start < b$end & end > b$start]) else data.table()
      if (!nrow(seg)) {
        seen_no_fire[[si]] <- union(seen_no_fire[[si]], rid_value)
      } else {
        tags <- unique(as.character(seg$tag[!is.na(seg$tag)]))
        row[, hp := if (length(tags) == 1L) tags else NA_character_]
        seg[, original_size := end - start]
        seg[, `:=`(left = pmax(start + 1L, ws, b$start + 1L), right = pmin(end, we, b$end))]
        seg <- seg[left <= right]
        seg[, `:=`(rel_start = if (orient == 1L) left - tss else tss - right,
                    rel_end = if (orient == 1L) right - tss else tss - left)]
        seg[, size := fifelse(left == ws | right == we, NA_real_, as.numeric(right - left + 1L))]
        setorder(seg, rel_start, rel_end)
        nuc <- seg[class == "nucleosome" & original_size >= nuc_min]
        at_tss <- seg[rel_start <= 10L & rel_end >= -10L]
        closed <- nuc[rel_start <= 10L & rel_end >= -10L]
        state <- if (nrow(closed)) "closed" else if (any(at_tss$class == "FIRE")) "open_fire" else if (any(at_tss$class == "linker")) "open_linker" else "other"
        row[, promoter_state := state]
        scores <- at_tss[class == "FIRE" & is.finite(fire_score), fire_score]
        if (length(scores)) row[, tss_fire_score := max(scores)]
        msp <- at_tss[class != "nucleosome"]
        row[, tss_msp_len := if (!nrow(msp)) 0 else if (anyNA(msp$size)) NA_real_ else max(msp$size)]
        covering <- nuc[rel_start <= 0 & rel_end >= 0]
        if (nrow(covering)) row[, `:=`(tss_nuc_start = covering$rel_start[1], tss_nuc_size = covering$size[1])]
        downstream <- nuc[rel_start > 10L]
        for (ni in seq_len(min(3L, nrow(downstream)))) {
          row[, (paste0("plus", ni, "_start")) := downstream$rel_start[ni]]
          row[, (paste0("plus", ni, "_size")) := downstream$size[ni]]
        }
        upstream <- nuc[rel_end < -10L]
        if (nrow(upstream)) {
          setorder(upstream, rel_end, rel_start)
          last <- upstream[.N]
          row[, `:=`(minus1_start = last$rel_start, minus1_end = last$rel_end, minus1_size = last$size)]
        }
        row[, ndr_len := plus1_start - minus1_end - 1L]
        dyads <- (nuc$rel_start + nuc$rel_end) / 2
        dyads <- sort(dyads[dyads >= 150 & dyads <= 2000])
        if (hi >= 2000L) {
          row[, down_n_nuc := length(dyads)]
          if (length(dyads) >= 3L) row[, down_dyad_spacing := median(diff(dyads))]
          occupancy <- numeric(1851L)
          for (ni in seq_len(nrow(nuc))) {
            a <- max(150L, nuc$rel_start[ni]); z <- min(2000L, nuc$rel_end[ni])
            if (a <= z) occupancy[seq.int(a, z) - 149L] <- 1
          }
          occupancy <- occupancy - mean(occupancy)
          if (sum(occupancy^2) > 0) {
            nfft <- 2L ^ ceiling(log2(2L * length(occupancy) - 1L))
            spectrum <- fft(c(occupancy, numeric(nfft - length(occupancy))))
            ac <- Re(fft(spectrum * Conj(spectrum), inverse = TRUE)) / nfft
            ac <- ac / ac[1]
            lags <- 140:250
            peaks <- lags[ac[lags + 1L] > 0 & ac[lags + 1L] > ac[lags] & ac[lags + 1L] >= ac[lags + 2L]]
            if (length(peaks)) {
              best <- peaks[which.max(ac[peaks + 1L])]
              row[, `:=`(down_acf_lag = best, down_acf_peak = ac[best + 1L])]
            }
          }
        }
        for (kind in c("nucleosome", "msp")) {
          intervals <- if (kind == "nucleosome") seg[class == "nucleosome"] else seg[class != "nucleosome"]
          fraction <- if (nrow(intervals)) sum(width(reduce(IRanges(intervals$rel_start, intervals$rel_end)))) / (hi - lo + 1L) else 0
          row[, (if (kind == "nucleosome") "frac_nuc" else "frac_msp") := fraction]
        }
      }
      if (!row$fp_ambiguous) {
        f <- if (nrow(fp)) copy(fp[RID == rid_value & start >= b$start & end <= b$end]) else data.table()
        row[, `:=`(pic = 0L, ppp = 0L, tf_up_n = 0L)]
        if (elong_len >= 300L && hi >= elong_end) row[, `:=`(elong_n = 0L, elong_per_kb = 0)]
        if (nrow(f)) {
          # BED centre in 0-based units is (start+end-1)/2; add one for TSS units.
          f[, rel := orient * ((start + end + 1) / 2 - tss)]
          row[, `:=`(pic = as.integer(any(f$size_class == "60-80" & abs(f$rel) <= 30)),
            ppp = as.integer(any(f$size_class == "40-60" & f$rel >= ppp_lo & f$rel <= ppp_hi)),
            tf_up_n = sum(f$size_class == "10-30" & f$rel >= -300 & f$rel <= -1))]
          if (elong_len >= 300L && hi >= elong_end) {
            count <- f[size_class == "40-60" & rel >= 300 & rel <= elong_end & !(rel >= ppp_lo & rel <= ppp_hi), .N]
            row[, `:=`(elong_n = count, elong_per_kb = count * 1000 / elong_len)]
          }
        }
      }
      rows[[ri]] <- row
    }
    block <- rbindlist(rows)
    gene_rows[[si]] <- block
    qc[si, rows_written := rows_written + nrow(block)]
  }
  # All samples share one observed-site universe for this gene. Buffer only this
  # gene's rows so a site found in the last sample also counts for the first.
  pooled <- rbindlist(gene_rows)
  if (nrow(pooled)) {
    pooled[, `:=`(m6a_site_n = length(observed_m6a), cpg_site_n = length(observed_cpg))]
    write.table(pooled, out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, na = "NA")
  }
  if (gi %% 25L == 0L || gi == nrow(genes)) message(chrom, ": ", gi, "/", nrow(genes), " genes; ", sum(qc$rows_written), " rows")
}
for (si in seq_len(nrow(samples))) {
  qc[si, `:=`(reads_overlapping = length(seen_overlap[[si]]),
              reads_spanning_window = length(seen_span[[si]]), reads_lacking_fire_all = length(seen_no_fire[[si]]))]
}
close(out)
if (!file.rename(tmp_path, out_path)) stop("Cannot finalize ", out_path)
fwrite(qc, file.path(opts[["out-dir"]], paste0("fiber_qc.", chrom, ".tsv")), sep = "\t", na = "NA")
message("Wrote ", out_path)
