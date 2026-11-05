# Paths, parameters, helpers and plots for the PolII example genes (run_gene_examples.R);
# pol2_functions.R and polII_pausing_plots.r merged. The shared parsing / plotting functions
# (code/parsing_functions) and the Leiden_Manhattan functions are sourced and used unchanged.
#
# Table coordinates are 1-based inclusive. Read-level positions are offsets from the gene's TSS
# in the direction of transcription (minus-strand genes are flipped), so offset 0 is the TSS base
# and positive offsets are downstream.
#
# plot_curated_gene() writes, per gene, plots/example_reads/<group>/<gene>/
#   unclustered/       <gene>_unclustered_m6a.pdf, <gene>_unclustered_footprints.pdf
#   clustered/         <gene>_clustered_m6a.pdf,   <gene>_clustered_footprints.pdf
#   read_level_plots/  <gene>_read_level.pdf: Tullius et al. Figure 1F, an R port of Data S1
#                      cell 54 (mmc3.zip)
# rtracklayer needs /project/spott/cshan/envs/dimelo/lib/libcurl.so.4 loaded first
# (dyn.load(..., local = FALSE), as in run_gene_examples.R).

suppressPackageStartupMessages({
  library(GenomicRanges); library(data.table)   # data.table last: its shift() masks IRanges'
})

PROJECT <- "/project/spott/cshan/fiber-seq"
ANN     <- "/project/spott/cshan/annotations"
LCL_PUB <- file.path(ANN, "LCL_public_data")
SHARED  <- "/project/spott/1_Shared_projects/LCL_Fiber_seq"

PLOTTING     <- file.path(PROJECT, "code/parsing_functions/plotting_functions.r")
PARSING      <- file.path(PROJECT, "code/parsing_functions/parsing_footprints_functions.r")
LEIDEN_DIR   <- file.path(PROJECT, "code/clustering_methods/Leiden_Manhattan")
PI_TABLE     <- file.path(PROJECT, "LCL_project/pausing/pause_index.tsv")   # code/pausing/cal_pausing_index/pause_index.R
OUT_DIR      <- file.path(PROJECT, "LCL_project/pausing/PolII")
DOMINANT_TSS <- file.path(ANN, "fantom5_LCL/fantom5.hg38.LCL.dominant_TSS.bed.gz")
CANONICAL_TSS <- file.path(ANN, "gencodev46_Ensembl_canonical_TSS.bed")
SAMPLES      <- file.path(SHARED, "Data/LCL_sample_metatable_merged_samples_31samples.csv")
# FiberHMM v1 footprints, all samples combined: closed size bins, together 10-80 bp (a footprint on
# a bin edge is in two files); plain gzip sorted by start, no index. The FiberHMM_v2 combined
# tracks are not usable here: one row per read (record span), not per footprint.
COMBINED_ROOT <- file.path(SHARED, "FiberHMM/merged/combined/joint_trained_tracks")
COMBINED_BINS <- c("10-30bp", "30-60bp", "60-80bp")
# per-base footprint coverage bigWigs of two of those bins, same folder (chr1-22)
COVERAGE_BINS <- c("30-60bp", "60-80bp")
coverage_file <- function(chrom, bin)
  file.path(COMBINED_ROOT, chrom, sprintf("combined_%s_%s_fp_cov.bw", chrom, bin))
# fibertools nucleosome calls per sample (ft extract --nuc, BED12 with 1-bp sentinel blocks, tabix)
nuc_file <- function(sample, chrom)
  file.path(SHARED, "FIRE/results", sample, "extracted_results/nuc_by_chr",
            sprintf("%s.ft_extracted_nuc.%s.bed.gz", sample, chrom))

CHAIN  <- file.path(ANN, "hg19ToHg38.over.chain")
GROSEQ <- c(`+` = file.path(LCL_PUB, "GRO-seq/GSM1480326_GM12878_GROseq_plus.bigWig"),     # hg19
            `-` = file.path(LCL_PUB, "GRO-seq/GSM1480326_GM12878_GROseq_minus.bigWig"))
DNASE  <- file.path(LCL_PUB, "DNase.ENCFF743ULW.bigWig")    # GM12878 ENCSR000EMT, GRCh38

AUTOSOMES    <- paste0("chr", 1:22)   # FiberHMM covers chr1-22 only
FLANK        <- 1000L                 # read window = TSS +/- FLANK
MIN_READS    <- 50L                   # example genes: fibers with a nucleosome record in the window
OTHER_TSS_BP <- 1000L                 # example genes: no other dominant TSS this close
N_EXAMPLES   <- 5L                    # example genes per PI bin
PI_BINS      <- c("low", "mid", "high")
CONTROL_GENES <- c("GADD45B", "HSPA1A", "HSPA1B", "HSPA8", "HSP90AA1", "DNAJB1",
                   "FOS", "JUNB", "EGR1", "IER2", "MYC")   # textbook paused genes
PLUS1_OFFSETS <- 0:299                # +1 nucleosome = first nucleosome-call base in here

# FiberHMM footprint size classes: assign_size_class() bins of parsing_footprints_functions.r
# matching COMBINED_BINS. A size on a shared break goes to the upper bin (assign_size_class()
# assigns bins in order, so 30 bp is 30-60 and 60 bp is 60-80, although its comment says lower).
TF_BREAKS  <- c(10L, 30L, 60L, 80L)
# fibertools nucleosomes kept, as load_region_results() / read_ft_nuc_region() (130-160 bp)
NUC_BREAKS <- c(130L, 160L)

EXAMPLE_RANGE <- c(-1000L, 1000L)    # inclusive offsets: 2,001 bases
READ_ROW_IN   <- 0.025               # PDF inches per read
X_LABEL       <- "bp from the TSS (upstream < 0; transcription ->)"

suppressPackageStartupMessages({ library(rtracklayer); library(ggplot2); library(patchwork) })
source(PLOTTING)   # theme_fiberseq(), save_figure(), plot_read_heatmap(), plot_group_profile(),
                   # plot_molecules(), stack_rows(), molecule_xlim(), cluster_id_palette(),
                   # MOLECULE_TRACK_COLORS
source(PARSING)    # assign_size_class(), read_tabix_region(), convert_ft_bed12_to_bed6(),
                   # assemble_region_m6a()
source(file.path(LEIDEN_DIR, "leiden_manhattan_functions.r"))   # leiden_manhattan_cluster()
source(file.path(LEIDEN_DIR, "leiden_manhattan_plots.r"))       # sample_palette(), plot_anchor(),
                   # plot_positions(), cluster_site_profiles() (plot_read_heatmap() helpers)

# footprint classes of the molecule panels, drawing order, named as in MOLECULE_TRACK_COLORS
FOOTPRINT_CLASSES <- c("nucleosome 130-160 bp",
                       paste("TF footprint", levels(assign_size_class(integer(), TF_BREAKS))))
stopifnot(all(FOOTPRINT_CLASSES %in% names(MOLECULE_TRACK_COLORS)))
COVERAGE_LABELS <- setNames(paste("TF footprint", sub("bp$", " bp", COVERAGE_BINS)), COVERAGE_BINS)
SIGNAL_TRACKS   <- c("DNase", "GRO-seq")


# TSS-relative first (lo) and last (hi) base of BED intervals [start0, end0)
to_relative <- function(start0, end0, tss, strand) {
  plus <- strand == "+"
  list(lo = as.integer(ifelse(plus, start0 + 1 - tss, tss - end0)),
       hi = as.integer(ifelse(plus, end0 - tss, tss - start0 - 1)))
}

# GM12878 GRO-seq of one chromosome lifted hg19 -> hg38:
# list(`+`, `-`) of GRanges with per-base read counts (minus stored negative)
lift_groseq <- function(chrom, chain) {
  lapply(GROSEQ, function(f) {
    len <- GenomeInfoDb::seqlengths(rtracklayer::BigWigFile(f))[[chrom]]
    gr <- rtracklayer::import(f, which = GRanges(chrom, IRanges(1L, len)))
    gr$score <- abs(gr$score)
    lo <- unlist(rtracklayer::liftOver(gr, chain))
    lo[as.character(seqnames(lo)) == chrom]
  })
}

# fibertools nucleosome calls of every LCL sample in chrom:start..end (1-based), pooled. One
# alignment per read (longest; palindromic pairs collapse), sentinel blocks dropped; a read in
# several samples is kept in the first.
# Output: list(reads = RID, sample, start0, end0 (read span); nucs = RID, start0, end0)
window_nucleosomes <- function(chrom, start, end) {
  region <- GRanges(chrom, IRanges(start, end))
  per_sample <- lapply(fread(SAMPLES)$sample_name, function(s) {
    f <- nuc_file(s, chrom)
    d <- read_tabix_region(f, region)
    if (!nrow(d)) return(NULL)
    names(d)[seq_along(FT_BED12_COLS)] <- FT_BED12_COLS
    d <- keep_longest_alignment(d)
    b <- convert_ft_bed12_to_bed6(d, format = "bed12_fibertools", source = f)
    list(reads = data.table(RID = d$RID, sample = s, start0 = d$start, end0 = d$end),
         nucs = data.table(RID = b$RID, sample = s, start0 = b$start, end0 = b$end))
  })
  reads <- rbindlist(lapply(per_sample, `[[`, "reads"))
  reads <- reads[!duplicated(RID)]
  nucs <- rbindlist(lapply(per_sample, `[[`, "nucs"))[reads[, .(RID, sample)], on = .(RID, sample), nomatch = 0L]
  list(reads = reads, nucs = nucs[, .(RID, start0, end0)])
}

# FiberHMM combined footprints (COMBINED_BINS, 10-80 bp) overlapping chrom:start..end (1-based).
# The files have no index, so each is streamed until its records start past `end`, in parallel.
# Output: RID, start0, end0 (unique)
window_tf_footprints <- function(chrom, start, end) {
  files <- file.path(COMBINED_ROOT, chrom, sprintf("combined_%s_%s_fps.bed.gz", chrom, COMBINED_BINS))
  cmds <- sprintf("zcat %s 2>/dev/null | awk -F'\\t' -v s=%d -v e=%d '$2 >= e {exit} $3 >= s'",
                  shQuote(files), as.integer(start), as.integer(end))
  out <- parallel::mclapply(cmds, function(x) {
    d <- fread(cmd = x, header = FALSE, sep = "\t")
    if (nrow(d)) setnames(d, c("chrom", "start0", "end0", "RID"))[, .(RID, start0, end0)]
  }, mc.cores = length(cmds))
  unique(rbindlist(out))
}

# Per-read matrix in the Tullius format: one row per read (the order of `reads`),
# one column per offset -FLANK..FLANK. A footprint base holds the footprint's
# length, any other base 0.
#   reads - RID; fps - RID, lo, hi, size (one gene)
footprint_matrix <- function(reads, fps) {
  m <- matrix(0L, nrow(reads), 2L * FLANK + 1L,
              dimnames = list(reads$RID, as.character(-FLANK:FLANK)))
  fill <- function(rid, lo, hi, value) {
    lo <- pmax(lo, -FLANK); hi <- pmin(hi, FLANK)
    keep <- lo <= hi & rid %in% reads$RID
    if (!any(keep)) return(invisible())
    len <- hi[keep] - lo[keep] + 1L
    cells <- cbind(rep(match(rid[keep], reads$RID), len), sequence(len, from = lo[keep] + FLANK + 1L))
    m[cells] <<- rep(value[keep], len)
  }
  fill(fps$RID, fps$lo, fps$hi, as.integer(fps$size))
  m
}

# +1 nucleosome: first offset in PLUS1_OFFSETS covered by a nucleosome call (m = footprint_matrix()
# of the nucleosome calls), NA when there is none
plus1_offset <- function(m) {
  sub <- (m[, as.character(PLUS1_OFFSETS), drop = FALSE] > 0) * 1L
  first <- max.col(sub, ties.method = "first")
  ifelse(rowSums(sub) > 0, PLUS1_OFFSETS[first], NA_integer_)
}


# --- example genes chosen by pausing index (not used by run_gene_examples.R) ----------

# Nearest other gene's CAGE dominant TSS (DOMINANT_TSS, any biotype, either strand) to each gene's
# TSS; the gene itself is matched by ENSG or symbol.
#   genes - gene_id, ensg, gene_name, chrom, tss, strand
# Output: gene_id, nearest_tss_bp, nearest_tss_gene, nearest_tss_orientation
nearest_other_tss <- function(genes) {
  dom <- fread(DOMINANT_TSS)[, .(dom_ensg = sub("\\..*", "", gene_id), dom_gene = gene_name,
                                 chrom = `#chrom`, pos = end, dom_strand = strand)]
  near <- merge(genes[, .(gene_id, ensg, gene_name, chrom, tss, strand)], dom, by = "chrom",
                allow.cartesian = TRUE)
  near <- near[(is.na(ensg) | dom_ensg != ensg) & dom_gene != gene_name]
  near[, dist := abs(pos - tss)]
  setorder(near, gene_id, dist)
  near[!duplicated(gene_id), .(gene_id, nearest_tss_bp = dist, nearest_tss_gene = dom_gene,
                               nearest_tss_orientation = fifelse(dom_strand == strand,
                                                                 "same strand", "opposite strand"))]
}

# Fibers (31 samples pooled, one alignment each) with an ft nucleosome record overlapping the
# gene's TSS +/- FLANK; g = one pause_index.tsv row
window_fiber_count <- function(g) {
  nrow(window_nucleosomes(g$chrom, g$tss - FLANK, g$tss + FLANK)$reads)
}

# Analysis-set genes with their nearest other TSS and matching scores:
#   candidate     no other gene's dominant TSS within OTHER_TSS_BP
#   in_fpkm_band  FPKM in the middle 50% of the analysis set
#   distance      |PI percentile - percentile of the tertile median| + |FPKM percentile - 0.5|
#   pi_tab - pause_index.tsv
example_candidates <- function(pi_tab) {
  an <- pi_tab[analysis_set == TRUE]
  an <- merge(an, nearest_other_tss(an), by = "gene_id", all.x = TRUE)
  an[, candidate := is.na(nearest_tss_bp) | nearest_tss_bp > OTHER_TSS_BP]
  band <- an[, quantile(fpkm, c(0.25, 0.75))]
  an[, in_fpkm_band := fpkm >= band[1] & fpkm <= band[2]]
  pi_ecdf <- ecdf(an$pi)
  an[, `:=`(pct_pi = frank(pi) / .N, pct_fpkm = frank(fpkm) / .N)]
  an[, target_pct_pi := pi_ecdf(median(pi)), by = pi_bin]
  an[, distance := abs(pct_pi - target_pct_pi) + abs(pct_fpkm - 0.5)]
  an[, pi_bin := factor(pi_bin, levels = PI_BINS)]
  setorder(an, pi_bin, distance)
  an[]
}

# Per PI tertile, walk the candidates in the FPKM band from the nearest to the matching target and
# keep the first N_EXAMPLES with >= MIN_READS fibers (window_fiber_count(), counted only for the
# genes looked at); adds n_reads
pick_examples <- function(an) {
  pool <- an[candidate & in_fpkm_band]
  picks <- rbindlist(lapply(PI_BINS, function(b) {
    kept <- list()
    for (i in which(pool$pi_bin == b)) {
      n <- window_fiber_count(pool[i])
      if (n >= MIN_READS) kept[[length(kept) + 1L]] <- pool[i][, n_reads := n]
      if (length(kept) == N_EXAMPLES) break
    }
    rbindlist(kept)
  }))
  picks[, group := as.character(pi_bin)][]
}

# Status of each CONTROL_GENES gene; the highest-PI one that "passes" (high tertile, no other TSS
# within OTHER_TSS_BP, >= MIN_READS fibers; FPKM band not required) is the positive control.
# Fibers are counted only for genes that pass the other rules.
control_status <- function(pi_tab, an) {
  controls <- merge(data.table(gene_name = CONTROL_GENES),
                    pi_tab[, .(gene_name, gene_id, chrom, tss, pi, fpkm, near_other_gene, short,
                               analysis_set, pi_bin)],
                    by = "gene_name", all.x = TRUE)
  controls <- merge(controls, an[, .(gene_id, nearest_tss_bp, candidate)], by = "gene_id", all.x = TRUE)
  controls[, status := fcase(is.na(gene_id), "not a RefSeq gene in pause_index.tsv",
                             !analysis_set & near_other_gene, "another gene within 3 kb",
                             !analysis_set & short, "shorter than 1 kb",
                             !analysis_set, "not in analysis set (expression, biotype, chrom or PI)",
                             pi_bin != "high", paste(pi_bin, "PI"),
                             !candidate, "other dominant TSS within 1 kb",
                             default = "passes")]
  controls[, n_reads := NA_integer_]
  for (i in which(controls$status == "passes")) controls[i, n_reads := window_fiber_count(controls[i])]
  controls[status == "passes" & n_reads < MIN_READS, status := paste("fewer than", MIN_READS, "fibers")]
  setorder(controls, -pi)
  controls[]
}

# picks plus the positive control (first "passes" row of control_status()), if any
add_control <- function(picks, an, controls) {
  control <- controls[status == "passes"]
  if (!nrow(control)) return(picks)
  row <- an[gene_id == control$gene_id[1]][, `:=`(n_reads = control$n_reads[1], group = "control")]
  rbind(picks, row, fill = TRUE)
}

# Write the chosen genes to tables/example_genes.tsv
write_example_table <- function(picks) {
  cols <- c("group", "gene_id", "gene_name", "chrom", "tss", "strand", "tes", "tx_id", "pi", "pi_bin",
            "fpkm", "tssr_dens", "body_dens", "n_reads",
            "nearest_tss_bp", "nearest_tss_gene", "nearest_tss_orientation",
            "candidate", "in_fpkm_band", "pct_pi", "pct_fpkm", "distance")
  path <- file.path(OUT_DIR, "tables/example_genes.tsv")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(picks[, ..cols], path, sep = "\t")
  invisible(path)
}


# --- curated housekeeping / cell type specific genes (run_gene_examples.R) ------------------

# Per-base summed score at `offsets` around gene g's TSS (direction of transcription) from
# intervals [start0, end0)
signal_at <- function(start0, end0, score, g, offsets) {
  out <- numeric(length(offsets))
  r <- to_relative(start0, end0, g$tss, rep(g$strand, length(start0)))
  lo <- pmax(r$lo, min(offsets)); hi <- pmin(r$hi, max(offsets))
  keep <- lo <= hi
  if (!any(keep)) return(out)
  len <- hi[keep] - lo[keep] + 1L
  s <- rowsum(rep(score[keep], len), sequence(len, from = lo[keep] - min(offsets) + 1L))
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

# DNase, GRO-seq (sense, lifted to hg38) and the FiberHMM combined COVERAGE_BINS footprint
# coverage (COVERAGE_LABELS columns; NA off chr1-22) at TSS -FLANK..FLANK
tss_tracks <- function(g) {
  offsets <- -FLANK:FLANK
  win <- GRanges(g$chrom, IRanges(g$tss - FLANK, g$tss + FLANK))
  gro <- subsetByOverlaps(lift_groseq(g$chrom, import.chain(CHAIN))[[g$strand]], win)
  dnase <- import(DNASE, which = win)
  out <- data.table(pos = offsets,
                    DNase = signal_at(start(dnase) - 1L, end(dnase), dnase$score, g, offsets),
                    `GRO-seq` = signal_at(start(gro) - 1L, end(gro), gro$score, g, offsets))
  for (b in COVERAGE_BINS) {
    value <- NA_real_
    if (g$chrom %in% AUTOSOMES) {
      cov <- import(coverage_file(g$chrom, b), which = win)
      value <- signal_at(start(cov) - 1L, end(cov), cov$score, g, offsets)
    }
    set(out, j = COVERAGE_LABELS[[b]], value = value)
  }
  out
}

# User-supplied groups; aliases and HLA component expansion remain explicit in the manifest.
curated_example_genes <- function() {
  hk <- strsplit("ACTB ACTG1 GAPDH B2M HPRT1 SDHA TFRC TBP EEF2 EIF2S3 GDI1 MYH9 SF3B1 SPCS3 RPS18 RPS27 POM121 AHSP HBS1L CYC1 PGK1 GUSB TUBB RPL13A ATP5B", " ")[[1]]
  ct <- strsplit("CD19 CD21 CD22 FCGR2B MS4A1 CD24 CD38 CD72 CD74 CD79A CD79B CD83 CD86 CD5 CD9 CD58 CD81 CD123 SYK LYN BTK BLNK BLK EBF PAX5 POU2AF1 SPIB SPI1 IRF4 IRF8 CEBPB BCL6 PRDM1 XBP1 MYC HLA-DM HLA-DO HLA-DP HLA-DQ HLA-DR IGHM IGHD IGHA1 IGKC IGLC2 LMO2 MYBL1 SERPINA9 MME HMCES ZFAND4 NLRP4 MUC4 SPP1 LAMB3 MERTK IKZF1", " ")[[1]]
  aliases <- list(ATP5B = "ATP5F1B", CD21 = "CR2", CD123 = "IL3RA", EBF = "EBF1",
    `HLA-DM` = c("HLA-DMA", "HLA-DMB"), `HLA-DO` = c("HLA-DOA", "HLA-DOB"),
    `HLA-DP` = c("HLA-DPA1", "HLA-DPB1"), `HLA-DQ` = c("HLA-DQA1", "HLA-DQB1"),
    `HLA-DR` = c("HLA-DRA", "HLA-DRB1", "HLA-DRB5"))
  rbindlist(lapply(c("housekeeping", "cell_type_specific"), function(group) {
    genes <- if (group == "housekeeping") hk else ct
    rbindlist(lapply(genes, function(g) data.table(group = group, requested_gene = g,
      gene_name = if (is.null(aliases[[g]])) g else aliases[[g]])))
  }))
}

# Prefer the output of code/pausing/dominant_tss/write_dominant_tss.R;
# fall back only when absent, using the canonical GENCODE v46 one-base BED.
resolve_example_tss <- function(gene, group, dominant = fread(DOMINANT_TSS)) {
  hit <- dominant[gene_name == gene]
  if (nrow(hit) > 1L) stop("multiple dominant tss found for ", gene, call. = FALSE)
  if (nrow(hit)) return(hit[, .(group = group, gene_name, gene_id, chrom = `#chrom`,
    tss = as.integer(end), strand, tss_source = "FANTOM5 dominant TSS")])
  canonical <- fread(CANONICAL_TSS, header = FALSE,
    col.names = c("chrom", "start0", "end", "annotation", "score", "strand"))
  canonical[, c("gene_id", "tx_id", "gene_name") := tstrsplit(annotation, ";", fixed = TRUE, keep = 1:3)]
  hit <- canonical[gene_name == gene]
  if (!nrow(hit)) stop("no dominant tss found for ", gene, "; no GENCODE canonical TSS found either", call. = FALSE)
  if (nrow(hit) != 1L) stop("multiple GENCODE canonical tss found for ", gene, call. = FALSE)
  stopifnot(hit$end - hit$start0 == 1L)
  hit[, .(group = group, gene_name, gene_id, chrom, tss = as.integer(end), strand,
           tss_source = "GENCODE v46 canonical TSS")]
}

gene_plot_label <- function(g) {
  if (!is.null(g$requested_gene) && g$requested_gene != g$gene_name)
    paste0(g$gene_name, " (", g$requested_gene, ")") else g$gene_name
}

# Title and subtitle of every figure: gene, group, view, TSS, read count and window
example_labels <- function(g, view, n) c(
  title = sprintf("%s | %s | %s", gene_plot_label(g), g$group, view),
  subtitle = sprintf("TSS %s:%d (%s) | %s\n%d full-span reads | TSS %d to +%d bp",
                     g$chrom, g$tss, g$strand, g$tss_source, n, EXAMPLE_RANGE[1], EXAMPLE_RANGE[2]))

# All panels use the same full-span m6A-read cohort, including reads without nucleosomes.
# Namespace IDs by sample to match assemble_region_m6a(), including footprint joins.
# Nucleosomes are the 130-160 bp (NUC_BREAKS) fibertools calls, as in load_region_results().
example_footprints <- function(g, dat) {
  region <- GRanges(g$chrom, IRanges(g$tss + EXAMPLE_RANGE[1], g$tss + EXAMPLE_RANGE[2]))
  meta <- as.data.table(dat$rids_df)
  if (anyDuplicated(meta$original_RID)) stop("Combined TF BED lacks sample IDs: ambiguous read names for ", g$gene_name)
  nucs <- rbindlist(lapply(unique(meta$sample_name), function(s) {
    d <- read_tabix_region(nuc_file(s, g$chrom), region)
    if (!nrow(d)) return(NULL)
    names(d)[seq_along(FT_BED12_COLS)] <- FT_BED12_COLS
    d <- keep_longest_alignment(d)
    b <- as.data.table(convert_ft_bed12_to_bed6(d, format = "bed12_fibertools", source = nuc_file(s, g$chrom)))
    b[, .(RID = paste(s, RID, sep = "::"), start0 = start, end0 = end)]
  }))
  if (!nrow(nucs)) nucs <- data.table(RID = character(), start0 = integer(), end0 = integer())
  nucs <- nucs[RID %in% meta$RID & !is.na(assign_size_class(end0 - start0, NUC_BREAKS))]
  tf <- if (g$chrom %in% AUTOSOMES) window_tf_footprints(g$chrom, start(region), end(region)) else
    data.table(RID = character(), start0 = integer(), end0 = integer())
  if (!nrow(tf)) tf <- data.table(RID = character(), start0 = integer(), end0 = integer())
  tf <- merge(tf, meta[, .(original_RID, qualified_RID = RID)], by.x = "RID", by.y = "original_RID", allow.cartesian = TRUE)
  tf <- tf[, .(RID = qualified_RID, start0, end0)]
  relative <- function(d) {
    r <- to_relative(d$start0, d$end0, g$tss, rep(g$strand, nrow(d)))
    d[, `:=`(lo = r$lo, hi = r$hi, size = end0 - start0)][]
  }
  list(reads = relative(meta[, .(RID, start0 = start - 1L, end0 = end)]),
       nucs = relative(nucs), fp = relative(tf), tracks = tss_tracks(g))
}

# Region row of plot_read_heatmap(), cluster_site_profiles() and plot_positions(): the
# EXAMPLE_RANGE window, TSS-anchored and flipped on the minus strand (plot_anchor())
example_region <- function(g) {
  data.frame(chr = g$chrom, analysis_start = g$tss + EXAMPLE_RANGE[1],
             analysis_end = g$tss + EXAMPLE_RANGE[2], width = diff(EXAMPLE_RANGE) + 1L,
             region_type = "promoter", strand = g$strand, tss = g$tss, annotation = "")
}

# plot_molecules() window in TSS-relative bp; start = end = 0 puts its dashed lines on the TSS
example_window <- function() {
  data.frame(chrom = "TSS", start = 0, end = 0,
             window_start = EXAMPLE_RANGE[1], window_end = EXAMPLE_RANGE[2])
}

# sample_palette() colours of the 31 LCL samples (the same across genes), for the samples present
example_sample_colors <- function(sample_names) {
  colors <- sample_palette(fread(SAMPLES)$sample_name)
  colors[names(colors) %in% sub("_.*$", "", sample_names)]
}

# Display rows in plot order: RID, cluster (factor; one "all reads" level when unclustered),
# start = row rank (the within-cluster sort key of plot_read_heatmap()), sample_name
example_rows <- function(read_order, meta, cluster = factor(rep("all reads", length(read_order)))) {
  data.frame(RID = read_order, cluster = cluster, start = seq_along(read_order),
             sample_name = as.character(meta$sample_name[match(read_order, meta$RID)]))
}

# Nucleosome and TF footprint intervals for plot_molecules(): TSS-relative, widened by 0.5 bp,
# classes named as in MOLECULE_TRACK_COLORS, nucleosomes first so TF footprints are drawn on top
footprint_tracks <- function(bundle) {
  tr <- rbind(bundle$nucs[, .(RID, lo, hi, class = rep(FOOTPRINT_CLASSES[1], .N))],
              bundle$fp[, .(RID, lo, hi, class = sprintf("TF footprint %s",
                                                       as.character(assign_size_class(size, TF_BREAKS))))])
  tr[, class := factor(class, levels = FOOTPRINT_CLASSES)]
  if (anyNA(tr$class)) stop("TF footprint sizes outside ", min(TF_BREAKS), "-", max(TF_BREAKS), " bp")
  setorder(tr, class)
  data.frame(RID = tr$RID, start = tr$lo - 0.5, end = tr$hi + 0.5, class = tr$class)
}

tss_vline <- function() geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3)

footprint_caption <- function(g) paste0(
  "Nucleosomes: fibertools ft nucleosome calls, 130-160 bp; TF footprints: FiberHMM combined 10-80 bp.",
  if (!g$chrom %in% AUTOSOMES) "\nFiberHMM TF footprints are unavailable on non-autosomal chromosomes." else "")

# m6A heatmap with plot_read_heatmap() (genomic layout): one slice per cluster with cluster and
# sample bars (sample bar only when unclustered), m6A fraction at observed sites per cluster on
# top, TSS dashed. rows = example_rows()
plot_example_m6a <- function(g, dat, rows, output, clustered) {
  clusters <- levels(rows$cluster)
  sample_colors <- example_sample_colors(rows$sample_name)
  labels <- example_labels(g, if (clustered) "Manhattan / Leiden" else "unclustered", nrow(rows))
  options <- list(column_title = paste(c(labels, "Top: m6A fraction at observed sites, no smoothing",
                                         X_LABEL), collapse = "\n"),
                  column_title_gp = grid::gpar(fontsize = 10))
  if (!clustered) options$left_annotation <- ComplexHeatmap::rowAnnotation(
    sample = factor(sub("_.*$", "", rows$sample_name), levels = names(sample_colors)),
    col = list(sample = sample_colors))
  ht <- plot_read_heatmap(list(assignments = rows, site_met_mat = dat$met_mat), example_region(g),
                          sample_colors, show_cluster_profiles = TRUE, heatmap_options = options)
  # offset 0 is column 1 - EXAMPLE_RANGE[1] in transcription order
  tss_x <- (0.5 - EXAMPLE_RANGE[1]) / (diff(EXAMPLE_RANGE) + 1)
  tss_line <- function() grid::grid.lines(x = c(tss_x, tss_x), y = c(0, 1),
                                          gp = grid::gpar(col = "#D55E00", lty = 2))
  save_figure(ht, output, width = 10, height = 6 + 0.55 * length(clusters), draw = function(x) {
    # room for the profile names plot_read_heatmap() draws right of the top annotation
    ComplexHeatmap::ht_opt(HEATMAP_LEGEND_PADDING = grid::unit(22, "mm"),
                           ANNOTATION_LEGEND_PADDING = grid::unit(22, "mm"))
    on.exit(ComplexHeatmap::ht_opt(RESET = TRUE))
    ComplexHeatmap::draw(x, merge_legend = TRUE)
    for (slice in seq_along(clusters)) ComplexHeatmap::decorate_heatmap_body("m6A", tss_line(), slice = slice)
    for (cl in clusters) ComplexHeatmap::decorate_annotation(paste0(cl, " m6A"), tss_line())
  })
}

# Footprints of the same rows with plot_molecules(): one facet per cluster with cluster and
# sample bars (sample bar only when unclustered), nucleosomes and TF footprints in
# MOLECULE_TRACK_COLORS; top: m6A fraction at observed sites per cluster (plot_group_profile()),
# as in plot_read_heatmap(show_cluster_profiles = TRUE)
plot_example_footprints <- function(g, dat, bundle, rows, output, clustered) {
  region <- example_region(g)
  rows$sample <- sub("_.*$", "", rows$sample_name)
  rows <- stack_rows(rows, "cluster")
  cluster_colors <- cluster_id_palette(levels(rows$cluster))
  sample_colors <- example_sample_colors(rows$sample_name)
  bars <- if (clustered) c("cluster", "sample") else "sample"
  ex <- example_window()
  lim <- molecule_xlim(ex, length(bars))
  breaks <- seq(EXAMPLE_RANGE[1], EXAMPLE_RANGE[2], 250)

  profile <- cluster_site_profiles(list(region = region, assignments = rows), dat$met_mat)
  observed <- as.numeric(colnames(dat$met_mat))[Matrix::colSums(dat$met_mat) > 0]
  profile <- profile[profile$pos %in% observed, , drop = FALSE]
  profile$cluster <- factor(profile$cluster, levels = levels(rows$cluster))
  top <- plot_group_profile(profile, group_col = "cluster", x_col = "pos", y_col = "met",
      style = "area", facet_col = "cluster", facet_layout = "grid",
      colors = cluster_colors, fill_colors = cluster_colors, coordinates = "relative", region = region,
      labels = list(x = NULL, y = "m6A fraction\n(observed sites)"),
      theme = theme_fiberseq("classic", base_size = 9)) +
    tss_vline() + scale_x_continuous(breaks = breaks) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    coord_cartesian(xlim = lim, expand = FALSE) +
    theme(legend.position = "none", strip.text.y = element_text(angle = 0),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())

  n <- table(rows$cluster)
  molecules <- plot_molecules(rows, footprint_tracks(bundle), ex, "cluster",
      setNames(sprintf("%s\n(n=%d)", names(n), as.integer(n)), names(n)),
      bar_cols = bars, bar_palette = c(cluster_colors, sample_colors),
      legend_extra = names(sample_colors), track_colors = MOLECULE_TRACK_COLORS[FOOTPRINT_CLASSES],
      breaks = breaks, theme_overrides = list(legend.position = "bottom")) +
    guides(fill = guide_legend(ncol = 9, byrow = TRUE)) + labs(x = X_LABEL)

  k <- length(n)
  heights <- c(0.65 * k, READ_ROW_IN * nrow(rows) + 0.15 * k + 1.2)
  labels <- example_labels(g, if (clustered) "Manhattan / Leiden" else "unclustered", nrow(rows))
  fig <- wrap_plots(top, molecules, ncol = 1, heights = heights) +
    plot_annotation(title = labels[["title"]], subtitle = labels[["subtitle"]],
                    caption = footprint_caption(g), theme = theme_fiberseq())
  save_figure(fig, output, width = 10, height = 1.6 + sum(heights), limitsize = FALSE)
}

# Tullius et al. Figure 1F (Data S1 cell 54) for one gene. Above: DNase and sense GRO-seq, each
# scaled to its window maximum. Middle: FiberHMM combined 30-60 / 60-80 bp footprint coverage
# bigWigs. Below: every full-span read with its nucleosomes and TF footprints in
# MOLECULE_TRACK_COLORS, ordered by +1 nucleosome (plus1_offset()); Tullius drops reads without
# one, here they are kept in a second facet, in read_order.
plot_read_level <- function(g, bundle, read_order, output) {
  rows <- data.table(RID = read_order)
  rows[, plus1 := plus1_offset(footprint_matrix(rows, bundle$nucs))]
  rows[, `:=`(group = factor(fifelse(is.na(plus1), "none", "plus1"), levels = c("plus1", "none")),
              rank = .I)]
  setorder(rows, group, plus1, rank)
  rows <- stack_rows(as.data.frame(rows), "group")
  n <- table(rows$group)
  facet_labels <- c(plus1 = sprintf("+1 nucleosome\n(n=%d)", n[["plus1"]]),
                    none = sprintf("no +1 nucleosome\n(n=%d)", n[["none"]]))
  ex <- example_window()
  breaks <- seq(EXAMPLE_RANGE[1], EXAMPLE_RANGE[2], 250)
  xs <- list(scale_x_continuous(breaks = breaks),
             coord_cartesian(xlim = molecule_xlim(ex, 0L), expand = FALSE))
  theme_rows <- theme_fiberseq(overrides = list(
    strip.background = element_blank(), strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.x = element_blank(), axis.line.x = element_blank()))

  tr <- bundle$tracks
  signal <- melt(tr[, c("pos", SIGNAL_TRACKS), with = FALSE], id.vars = "pos", variable.name = "track")
  signal[, value := if (max(value) > 0) value / max(value) else value, by = track]
  p_signal <- ggplot(signal, aes(pos, value)) + tss_vline() +
    geom_line(colour = "grey20", linewidth = 0.4) +
    facet_grid(track ~ ., switch = "y",
               labeller = as_labeller(c(DNase = "DNase", `GRO-seq` = "GRO-seq (sense)"))) +
    scale_y_continuous(breaks = c(0, 1), position = "right") + xs +
    labs(y = "scaled to max") + theme_rows

  p_cover <- if (g$chrom %in% AUTOSOMES) {
    cover <- melt(tr[, c("pos", unname(COVERAGE_LABELS)), with = FALSE], id.vars = "pos",
                  variable.name = "track")
    ggplot(cover, aes(pos, value, colour = track)) + tss_vline() + geom_line(linewidth = 0.4) +
      facet_grid(track ~ ., switch = "y", scales = "free_y") +
      scale_colour_manual(values = MOLECULE_TRACK_COLORS, guide = "none") +
      scale_y_continuous(n.breaks = 3, position = "right") + xs +
      labs(y = "footprint coverage") + theme_rows
  } else {
    ggplot() + annotate("text", x = 0, y = 0.5, size = 3,
                        label = "FiberHMM footprint coverage bigWigs cover chr1-22 only") +
      xs + labs(y = "footprint coverage") + theme_rows +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
  }

  molecules <- plot_molecules(rows, footprint_tracks(bundle), ex, "group", facet_labels,
      track_colors = MOLECULE_TRACK_COLORS[FOOTPRINT_CLASSES], breaks = breaks) +
    labs(x = X_LABEL) + theme(legend.position = "top")

  reads_height <- READ_ROW_IN * nrow(rows) + 1.1
  labels <- example_labels(g, "read level", nrow(rows))
  fig <- (p_signal / p_cover / molecules) + plot_layout(heights = c(1.1, 1.1, reads_height)) +
    plot_annotation(title = labels[["title"]], subtitle = labels[["subtitle"]],
      caption = paste0(
        "Above: DNase (ENCFF743ULW) and sense GRO-seq (GSM1480326, lifted hg19 -> hg38), each scaled to its window maximum.\n",
        "Middle: FiberHMM combined footprint coverage bigWigs (joint_trained_tracks).\n",
        "Below: reads ordered by +1 nucleosome (first nucleosome base in +", min(PLUS1_OFFSETS), "..+",
        max(PLUS1_OFFSETS), " bp), as in Tullius et al. Fig 1F; reads without one are shown below.\n",
        footprint_caption(g)),
      theme = theme_fiberseq())
  save_figure(fig, output, width = 10, height = 3.8 + reads_height, limitsize = FALSE)
}

plot_curated_gene <- function(g, gene_dir = file.path(OUT_DIR, "plots/example_reads", g$group, g$gene_name)) {
  region <- data.frame(region_id = g$gene_name, chr = g$chrom,
    start = g$tss + EXAMPLE_RANGE[1] - 1L, end = g$tss + EXAMPLE_RANGE[2],
    analysis_start = g$tss + EXAMPLE_RANGE[1], analysis_end = g$tss + EXAMPLE_RANGE[2],
    tss = g$tss, strand = g$strand, tss_source = g$tss_source)
  samples <- fread(SAMPLES)[, .(sample_name, fire_dir)]
  dat <- assemble_region_m6a(as.data.frame(samples), region, full_span = TRUE,
    positions = seq(region$analysis_start, region$analysis_end))
  meta <- dat$rids_df
  plain <- meta$RID[order(match(meta$sample_name, samples$sample_name), meta$start, meta$RID)]
  path <- function(view, suffix) file.path(gene_dir, view, paste0(g$gene_name, "_", suffix, ".pdf"))
  unclustered_rows <- example_rows(plain, meta)
  plot_example_m6a(g, dat, unclustered_rows, path("unclustered", "unclustered_m6a"), clustered = FALSE)
  result <- leiden_manhattan_cluster(dat$met_mat, dat$rids_df,
    region$analysis_start, region$analysis_end, window_size = 0L,
    k_neighbors = 10L, resolution = 1, seed = 1L)
  a <- result$assignments
  clustered <- a$RID[order(a$cluster, match(a$RID, plain))]
  clustered_rows <- example_rows(clustered, meta, a$cluster[match(clustered, a$RID)])
  plot_example_m6a(g, dat, clustered_rows, path("clustered", "clustered_m6a"), clustered = TRUE)
  tables <- file.path(OUT_DIR, "tables", g$group)
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)
  a$unclustered_row <- match(a$RID, plain)
  a$clustered_row <- match(a$RID, clustered)
  fwrite(a, file.path(tables, paste0(g$gene_name, "_read_clusters.tsv")), sep = "\t")
  fwrite(dat$qc, file.path(tables, paste0(g$gene_name, "_read_qc.tsv")), sep = "\t")
  saveRDS(list(region = region, params = result$params, assignments = a),
    file.path(tables, paste0(g$gene_name, "_clustering.rds")))
  bundle <- example_footprints(g, dat)
  plot_example_footprints(g, dat, bundle, unclustered_rows, path("unclustered", "unclustered_footprints"),
                          clustered = FALSE)
  plot_example_footprints(g, dat, bundle, clustered_rows, path("clustered", "clustered_footprints"),
                          clustered = TRUE)
  plot_read_level(g, bundle, plain, path("read_level_plots", "read_level"))
  data.table(n_reads = nrow(meta), n_clusters = result$n_clusters,
    tf_coverage = if (g$chrom %in% AUTOSOMES) "FiberHMM" else "unavailable (non-autosomal)")
}

# Record failures, process remaining genes, then raise collected errors.
run_curated_examples <- function(genes = curated_example_genes(), workers = 1L) {
  dir.create(file.path(OUT_DIR, "tables/status"), recursive = TRUE, showWarnings = FALSE)
  dominant <- fread(DOMINANT_TSS)
  manifest <- copy(genes)
  manifest[, `:=`(status = "pending", error = "", chrom = NA_character_, tss = NA_integer_,
                 strand = NA_character_, tss_source = NA_character_, n_reads = NA_integer_, n_clusters = NA_integer_)]
  path <- file.path(OUT_DIR, "tables/example_genes.tsv")
  fwrite(manifest, path, sep = "\t")
  run_one <- function(i) {
    row <- copy(manifest[i])
    message("Gene ", i, "/", nrow(manifest), ": ", row$gene_name)
    tryCatch({
      g <- resolve_example_tss(row$gene_name, row$group, dominant)
      g[, requested_gene := row$requested_gene]
      row[, `:=`(chrom = g$chrom, tss = g$tss, strand = g$strand, tss_source = g$tss_source)]
      stats <- plot_curated_gene(g)
      row[, `:=`(status = "complete", n_reads = stats$n_reads, n_clusters = stats$n_clusters)]
    }, error = function(e) {
      row[, `:=`(status = "error", error = conditionMessage(e))]
      message("ERROR: ", conditionMessage(e))
    })
    fwrite(row, file.path(OUT_DIR, "tables/status", paste0(row$gene_name, ".tsv")), sep = "\t")
    row
  }
  results <- parallel::mclapply(seq_len(nrow(manifest)), run_one,
    mc.cores = workers, mc.preschedule = FALSE)
  if (any(vapply(results, inherits, logical(1), "try-error"))) stop("A gene worker failed; see tables/status")
  manifest <- rbindlist(results)
  fwrite(manifest, path, sep = "\t")
  if (any(manifest$status == "error")) stop(paste(manifest[status == "error", paste(gene_name, error, sep = ": ")], collapse = "\n"), call. = FALSE)
  invisible(manifest)
}
