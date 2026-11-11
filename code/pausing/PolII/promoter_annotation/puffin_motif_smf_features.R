#!/usr/bin/env Rscript
# Features of TATA promoters on single fibers: m6A, +1 nucleosome and footprint occupancy, TATA
# promoters vs all other promoters with every curated PolII example gene pooled (no split by
# housekeeping / cell type specific / not expressed).
# Also tests sequence predictors of expression variability (Einarsson et al., 2022 adjusted
# log10 CV2): observed/expected CpG ratio (expected to predict low variability) and the Puffin
# TATA box / TBP site (expected to predict high variability), ranked against the other Puffin
# motif contributions.
# Reads only tables: Puffin class and TATA box position (04a_puffin_motif_scores.py), published
# variability (gene_annotations.tsv of expression_variability.R), the read features and 10 bp
# bins of 02_extract_read_features.R and the per gene 10 bp profiles of
# 03_cluster_gene_summary.R; CpG o/e from hg38. Nothing is re-extracted.
#
# Run (after 04a and expression_variability.R): bash code/pausing/PolII/run_gene_examples.sh \
#        code/pausing/PolII/promoter_annotation/puffin_motif_smf_features.R
#
# Outputs under LCL_project/pausing/PolII/promoter_annotation/
#   tables/puffin_motif_smf_genes.tsv       one row per gene: TATA box, m6A, +1, footprint features
#   tables/puffin_motif_smf_reads.tsv.gz    one row per read: TATA box state, +1 start, TATA -> +1
#   tables/puffin_motif_smf_tests.tsv       TATA vs other promoters, Wilcoxon per feature, BH
#   tables/puffin_motif_smf_profiles.tsv    mean gene profile per motif x (motif genes / all other genes)
#                                           x 10 bp bin x signal
#   tables/puffin_motif_smf_motif_tests.tsv each motif's genes vs all other genes, Wilcoxon per SMF
#                                           feature, BH over motif x feature
#   tables/puffin_motif_smf_motif_correlations.tsv
#                                           Spearman of each Puffin motif contribution vs core occupancy
#   tables/puffin_motif_smf_box_state.tsv   +1 fraction / position per TATA box state
#   tables/puffin_motif_smf_variability_predictors.tsv
#                                           Spearman of each sequence predictor vs adjusted log10 CV2
#   plots/puffin_motif_smf_features.pdf     one mean gene profile page per motif (motif genes vs all
#                                           other genes), motif contribution x occupancy heatmap, CpG o/e
#                                           and TATA box vs expression variability
#
# Motif class = Puffin group_excl_U1 of 04a (dominant motif, U1 snRNP excluded); each motif's
# genes are compared with all other genes.
# TATA box = tata_center of 04a (peak TATA+ motif activation within TSS -60..-10) for Puffin
# class TATA; other promoters have no TATA box and are read at CONTROL_POSITION instead.
# +1 start = plus1_start of 02: first nucleosome start within +0..+299 on reads whose TSS base is
# not covered by a nucleosome (NA otherwise). TATA -> +1 = +1 start - TATA box centre (upstream
# edge; the dyad is ~73 bp further). All genes with fiber data; every gene weighted equally.
# CpG o/e = CpG count x window length / (C count x G count) over TSS +/- CPG_FLANK bp
# (Gardiner-Garden & Frommer, 1987). Variability tests use genes matched to a published promoter.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rsamtools)
  library(Biostrings)
})
setDTthreads(1L)

BASE <- "/project/spott/cshan/fiber-seq/LCL_project/pausing/PolII"
PUFFIN_OUT <- file.path(BASE, "promoter_annotation")
GENOME <- "/project/spott/cshan/annotations/hg38.fa"
CPG_FLANK <- 200L
CONTROL_POSITION <- -30L
MIN_PLUS1_READS <- 10L
PROFILE_XLIM <- c(-500, 500)
PROMOTER_COLORS <- c(TATA = "#E41A1C", other = "#5b6b7a")  # TATA = Puffin's TATA colour
MOTIF_COLORS <- c(TATA = "#E41A1C", SP1 = "#377EB8", ETS = "#4DAF4A", NFY = "#984EA3", NRF1 = "#FF7F00",
                  CREB = "#A65628", YY1 = "#F781BF", ZNF143 = "#1B9E77")
OTHER_GENES <- "all other genes"
OTHER_COLOR <- "#5b6b7a"
MOTIF_FEATURES <- c("m6a_core", "m6a_upstream", "m6a_downstream", "nuc_core", "tf_10_30_core", "tf_30_60_core",
                    "tf_60_80_core", "occluded_fraction", "plus1_fraction", "plus1_start_median", "plus1_start_iqr")
CONTRIBUTIONS <- c(paste0(c("TATA", "SP", "ETS", "NFY", "NRF1", "CREB", "YY1", "ZNF143"), "_contribution"),
                   "tata_activation")
CORE_OCCUPANCY <- c(m6a_core = "m6A", nuc_core = "nucleosome", tf_10_30_core = "10-30 bp footprints",
                    tf_30_60_core = "30-60 bp footprints", tf_60_80_core = "60-80 bp footprints")
STATES <- c(nucleosome = "TATA box under a nucleosome", footprint = "TATA box footprinted (10-60 bp)",
            open = "TATA box accessible")
SIGNALS <- c(m6a = "m6A fraction", nuc = "nucleosome occupancy", tf_10_30 = "10-30 bp footprint occupancy",
             tf_30_60 = "30-60 bp footprint occupancy", tf_60_80 = "60-80 bp footprint occupancy")
FEATURES <- c(m6a_core = "m6A, core (-200..200)", m6a_upstream = "m6A, upstream (-500..-201)",
              m6a_downstream = "m6A, downstream (201..500)", m6a_tata_box = "m6A at TATA box bin",
              plus1_fraction = "+1 fraction", plus1_start_median = "+1 start, median (bp from TSS)",
              plus1_start_iqr = "+1 start, IQR across reads (bp)", tata_to_plus1 = "TATA box -> +1 edge (bp)",
              occluded_fraction = "TSS occluded fraction", nuc_core = "nucleosome, core",
              tf_10_30_core = "10-30 bp footprints, core", tf_30_60_core = "30-60 bp footprints, core",
              tf_60_80_core = "60-80 bp footprints, core",
              tata_footprint_fraction = "TATA box footprinted (read fraction)",
              tata_nucleosome_fraction = "TATA box under nucleosome (read fraction)")


############################################################
# genes and reads
############################################################
load_genes <- function() {
  pf <- fread(file.path(PUFFIN_OUT, "tables/04a_puffin_motif_scores.tsv"), check.names = FALSE)
  gf <- fread(file.path(BASE, "read_features/tables/gene_features.tsv"))[status == "complete"]
  g <- pf[gene_name %in% gf$gene_name, .(gene_name, group_excl_U1, motif = group_excl_U1, tata_score = `TATA +`, tata_center)]
  g <- merge(g, gf[, .(gene_name, group, n_reads, m6a_core, m6a_upstream, m6a_downstream, occluded_fraction,
                       nuc_core, tf_10_30_core, tf_30_60_core, tf_60_80_core)], by = "gene_name")
  g[, promoter := factor(fifelse(group_excl_U1 == "TATA", "TATA", "other"), names(PROMOTER_COLORS))]
  g[, tata_box := fifelse(promoter == "TATA", tata_center, CONTROL_POSITION)]
  g[]
}

# per read: nucleosome / footprint state of the 10 bp bin holding the TATA box, and the +1 start
load_reads <- function(genes) {
  reads <- fread(file.path(BASE, "read_features/tables/read_features.tsv.gz"),
                 select = c("gene", "RID", "nuc_at_tss", "plus1_start"))[gene %in% genes$gene_name]
  state <- rbindlist(lapply(genes$gene_name, function(gene) {
    b <- genes[gene_name == gene, (tata_box + 1000L) %/% 10L]
    x <- fread(file.path(BASE, "read_features/bins", paste0(gene, ".tsv.gz")))[bin == b]
    x[, .(gene, RID, state = fifelse(nuc_bases > 0, "nucleosome",
                                     fifelse(tf_10_30_bases + tf_30_60_bases > 0, "footprint", "open")))]
  }))
  r <- merge(reads, state, by = c("gene", "RID"))
  stopifnot(nrow(r) == nrow(reads))
  r <- merge(r, genes[, .(gene = gene_name, promoter, tata_box)], by = "gene")
  r[, state := factor(STATES[state], STATES)]
  r[, tata_to_plus1 := plus1_start - tata_box]
  r[]
}

# per gene 10 bp profiles of 03 (pooled count ratios), oriented by transcription
load_profiles <- function(genes) {
  p <- fread(file.path(BASE, "recurrent_patterns/tables/03_gene_bin_profiles.tsv.gz"))[gene %in% genes$gene_name]
  merge(p, genes[, .(gene = gene_name, promoter, motif, tata_box)], by = "gene")
}

gene_summary <- function(genes, reads, prof) {
  s <- reads[, .(plus1_fraction = mean(!is.na(plus1_start)), n_plus1_reads = sum(!is.na(plus1_start)),
                 plus1_start_median = as.numeric(median(plus1_start, na.rm = TRUE)),
                 plus1_start_iqr = IQR(plus1_start, na.rm = TRUE),
                 tata_nucleosome_fraction = mean(state == STATES[["nucleosome"]]),
                 tata_footprint_fraction = mean(state == STATES[["footprint"]]),
                 tata_open_fraction = mean(state == STATES[["open"]])), by = gene]
  box <- prof[bin == (tata_box + 1000L) %/% 10L, .(gene, m6a_tata_box = m6a)]
  g <- merge(merge(genes, s, by.x = "gene_name", by.y = "gene"), box, by.x = "gene_name", by.y = "gene")
  g[n_plus1_reads < MIN_PLUS1_READS, c("plus1_start_median", "plus1_start_iqr") := NA]
  g[, tata_to_plus1 := plus1_start_median - tata_box]
  g[order(promoter, gene_name)]
}

tata_tests <- function(g) {
  res <- rbindlist(lapply(names(FEATURES), function(f) {
    t <- g[promoter == "TATA", get(f)]; o <- g[promoter == "other", get(f)]
    t <- t[!is.na(t)]; o <- o[!is.na(o)]
    list(feature = f, n_tata = length(t), n_other = length(o), median_tata = median(t),
         median_other = median(o), p = suppressWarnings(wilcox.test(t, o))$p.value)
  }))
  res[, fdr := p.adjust(p, "BH")]
  res[]
}

# per motif: mean gene profile of that motif's genes and of all other genes
motif_profiles <- function(prof) {
  stopifnot(all(prof$motif %in% names(MOTIF_COLORS)))
  long <- melt(prof, id.vars = c("gene", "motif", "position"), measure.vars = names(SIGNALS),
               variable.name = "signal")[!is.na(value)]
  rbindlist(lapply(names(MOTIF_COLORS), function(m) {
    x <- copy(long)[, set := fifelse(motif == m, m, OTHER_GENES)]
    x[, .(mean = mean(value), se = sd(value) / sqrt(.N), n_genes = .N), by = .(set, signal, position)][
      , motif := m][]
  }))[, .(motif, set, signal, position, mean, se, n_genes)][order(motif, signal, set, position)]
}

motif_class_tests <- function(g) {
  res <- rbindlist(lapply(names(MOTIF_COLORS), function(m) rbindlist(lapply(MOTIF_FEATURES, function(f) {
    a <- g[motif == m, get(f)]; o <- g[motif != m, get(f)]
    a <- a[!is.na(a)]; o <- o[!is.na(o)]
    list(motif = m, feature = f, n_motif = length(a), n_other = length(o), median_motif = median(a),
         median_other = median(o),
         p = if (length(a) && length(o)) suppressWarnings(wilcox.test(a, o))$p.value else NA_real_)
  }))))
  res[, fdr := p.adjust(p, "BH")]
  res[]
}

# continuous Puffin motif contributions vs core occupancy, every gene with fiber data
motif_correlations <- function(g) {
  a <- fread(file.path(PUFFIN_OUT, "tables/gene_annotations.tsv"))[, c("gene_name", CONTRIBUTIONS), with = FALSE]
  x <- merge(g[, c("gene_name", names(CORE_OCCUPANCY)), with = FALSE], a, by = "gene_name")
  stopifnot(nrow(x) == nrow(g))
  res <- rbindlist(lapply(CONTRIBUTIONS, function(cc) rbindlist(lapply(names(CORE_OCCUPANCY), function(f) {
    ok <- is.finite(x[[cc]]) & is.finite(x[[f]])
    ct <- suppressWarnings(cor.test(x[[cc]][ok], x[[f]][ok], method = "spearman", exact = FALSE))
    list(contribution = cc, occupancy = f, n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
  }))))
  res[, fdr := p.adjust(p, "BH")]
  res[]
}

box_state_summary <- function(reads) {
  reads[, .(n_reads = .N, n_genes = uniqueN(gene), plus1_fraction = mean(!is.na(plus1_start)),
            tss_occluded_fraction = mean(nuc_at_tss),
            plus1_start_median = as.numeric(median(plus1_start, na.rm = TRUE))),
        by = .(promoter, state)][order(promoter, state)]
}


############################################################
# sequence predictors of expression variability
############################################################
cpg_oe <- function(chrom, tss) {
  seqs <- toupper(as.character(scanFa(FaFile(GENOME),
    GenomicRanges::GRanges(chrom, IRanges::IRanges(tss - CPG_FLANK, tss + CPG_FLANK)))))
  dna <- DNAStringSet(seqs)
  cg <- vcountPattern("CG", dna)
  base <- letterFrequency(dna, c("C", "G"))
  cg * width(dna) / (base[, "C"] * base[, "G"])
}

load_variability <- function() {
  a <- fread(file.path(PUFFIN_OUT, "tables/gene_annotations.tsv"))[match_status == "matched"]
  pf <- fread(file.path(PUFFIN_OUT, "tables/04a_puffin_motif_scores.tsv"), check.names = FALSE)
  v <- merge(a, pf[, .(gene_name, group_excl_U1, initiator_effect_tss, trinucleotide_effect_tss)], by = "gene_name")
  v[, cpg_oe := cpg_oe(chrom, tss)]
  v[, tata_class := as.integer(group_excl_U1 == "TATA")]
  v[, promoter := factor(fifelse(tata_class == 1L, "TATA", "other"), names(PROMOTER_COLORS))]
  v[]
}

PREDICTORS <- c(cpg_oe = "CpG o/e (TSS +/- 200 bp)", tata_activation = "TATA box activation (Puffin)",
                tata_class = "TATA promoter class (Puffin)", initiator_effect_tss = "Initiator effect (Puffin)",
                trinucleotide_effect_tss = "Trinucleotide effect (Puffin)",
                setNames(paste(c("TATA", "YY1", "ETS", "SP", "NFY", "CREB", "NRF1", "ZNF143"), "contribution (Puffin)"),
                         paste0(c("TATA", "YY1", "ETS", "SP", "NFY", "CREB", "NRF1", "ZNF143"), "_contribution")))

variability_predictors <- function(v) {
  res <- rbindlist(lapply(names(PREDICTORS), function(f) {
    ct <- suppressWarnings(cor.test(v[[f]], v$elife_adjusted_log10_CV2, method = "spearman", exact = FALSE))
    list(predictor = f, label = PREDICTORS[[f]], n = sum(is.finite(v[[f]])), rho = unname(ct$estimate), p = ct$p.value)
  }))
  res[, fdr := p.adjust(p, "BH")]
  res[order(-abs(rho))][, rank := .I][]
}

variability_plots <- function(v, pred) {
  lab <- function(f) with(pred[predictor == f], sprintf("Spearman rho = %.2f, p = %.2g, n = %d", rho, p, n))
  n <- v[, .N, by = promoter]
  nl <- setNames(sprintf("%s (n=%d)", n$promoter, n$N), n$promoter)[names(PROMOTER_COLORS)]
  cpg <- ggplot(v, aes(cpg_oe, elife_adjusted_log10_CV2, colour = promoter)) +
    geom_hline(yintercept = 0, linetype = 3, colour = "grey60") +
    geom_point(size = 2, alpha = 0.8) +
    scale_colour_manual(values = PROMOTER_COLORS, labels = nl, name = "Puffin class") +
    labs(x = PREDICTORS[["cpg_oe"]], y = "Adjusted log10(CV\u00b2)",
         title = "CpG o/e vs expression variability", subtitle = lab("cpg_oe")) +
    theme_bw(base_size = 11)
  tata <- ggplot(v, aes(promoter, elife_adjusted_log10_CV2, colour = promoter)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.7) +
    scale_colour_manual(values = PROMOTER_COLORS, guide = "none") +
    scale_x_discrete(labels = nl) +
    labs(x = NULL, y = "Adjusted log10(CV\u00b2)", title = "TATA box (TBP site) vs expression variability",
         subtitle = sprintf("%s; Puffin TATA activation: %s",
                            sub("Spearman", "class", lab("tata_class")), lab("tata_activation"))) +
    theme_bw(base_size = 11)
  rk <- copy(pred)[, label := factor(label, rev(label))]
  rank_plot <- ggplot(rk, aes(rho, label, fill = rho > 0)) +
    geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = "grey40") +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                      labels = c(`TRUE` = "higher variability", `FALSE` = "lower variability"), name = NULL) +
    labs(x = "Spearman rho with adjusted log10(CV\u00b2)", y = NULL,
         title = "Sequence predictors of expression variability, ranked by |rho|") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  list(cpg, tata, rank_plot)
}


############################################################
# figures
############################################################
motif_profile_plot <- function(mp, g, m) {
  d <- mp[motif == m & position >= PROFILE_XLIM[1] & position <= PROFILE_XLIM[2]]
  d[, signal := factor(SIGNALS[as.character(signal)], SIGNALS)]
  d[, set := factor(set, c(m, OTHER_GENES))]
  cols <- setNames(c(MOTIF_COLORS[[m]], OTHER_COLOR), c(m, OTHER_GENES))
  n_m <- g[motif == m, .N]
  labels <- setNames(c(sprintf("%s (n=%d)", m, n_m), sprintf("%s (n=%d)", OTHER_GENES, nrow(g) - n_m)), names(cols))
  p <- ggplot(d, aes(position, mean, colour = set, fill = set)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50")
  xlab <- "Position relative to TSS (bp, 10 bp bins)"
  if (m == "TATA") {
    p <- p + geom_vline(xintercept = g[motif == "TATA", median(tata_box)], linetype = 3, colour = cols[["TATA"]])
    xlab <- paste0(xlab, "; dotted red line = median TATA box")
  }
  p + geom_ribbon(aes(ymin = mean - se, ymax = mean + se), colour = NA, alpha = 0.2) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~ signal, scales = "free_y", ncol = 1, strip.position = "left") +
    scale_colour_manual(values = cols, labels = labels, name = NULL) +
    scale_fill_manual(values = cols, labels = labels, name = NULL) +
    labs(x = xlab, y = NULL, title = sprintf("Mean gene profiles (+/- se), %s promoters (n = %d genes) vs all other genes (n = %d)", m, n_m, nrow(g) - n_m)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", strip.background = element_blank(), strip.placement = "outside")
}

correlation_heatmap <- function(cr) {
  d <- copy(cr)
  d[, contribution := factor(sub("_contribution$", "", contribution),
                             rev(sub("_contribution$", "", CONTRIBUTIONS)))]
  d[, occupancy := factor(CORE_OCCUPANCY[occupancy], CORE_OCCUPANCY)]
  d[, cell := sprintf("%.2f%s", rho, fifelse(fdr < 0.05, "*", ""))]
  lim <- max(abs(d$rho), na.rm = TRUE)
  ggplot(d, aes(occupancy, contribution, fill = rho)) +
    geom_tile(colour = "white") + geom_text(aes(label = cell), size = 3.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-lim, lim), name = "Spearman rho") +
    labs(x = "Core occupancy (TSS +/- 200 bp)", y = "Puffin motif contribution",
         title = sprintf("Puffin motif contributions vs single-molecule occupancy (n = %d genes)", max(d$n)),
         subtitle = "* BH FDR < 0.05 over all contribution x occupancy pairs") +
    theme_bw(base_size = 11) + theme(panel.grid = element_blank())
}

main <- function() {
  genes <- load_genes()
  reads <- load_reads(genes)
  prof <- load_profiles(genes)
  g <- gene_summary(genes, reads, prof)
  tests <- tata_tests(g)
  mp <- motif_profiles(prof)
  motif_tests <- motif_class_tests(g)
  cr <- motif_correlations(g)
  states <- box_state_summary(reads)
  v <- load_variability()
  pred <- variability_predictors(v)
  dir.create(file.path(PUFFIN_OUT, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(PUFFIN_OUT, "plots"), recursive = TRUE, showWarnings = FALSE)
  fwrite(g, file.path(PUFFIN_OUT, "tables/puffin_motif_smf_genes.tsv"), sep = "\t", na = "NA")
  fwrite(reads[, .(gene, RID, promoter, tata_box, state, nuc_at_tss, plus1_start, tata_to_plus1)],
         file.path(PUFFIN_OUT, "tables/puffin_motif_smf_reads.tsv.gz"), sep = "\t", na = "NA")
  fwrite(tests, file.path(PUFFIN_OUT, "tables/puffin_motif_smf_tests.tsv"), sep = "\t", na = "NA")
  fwrite(mp, file.path(PUFFIN_OUT, "tables/puffin_motif_smf_profiles.tsv"), sep = "\t", na = "NA")
  fwrite(motif_tests, file.path(PUFFIN_OUT, "tables/puffin_motif_smf_motif_tests.tsv"), sep = "\t", na = "NA")
  fwrite(cr, file.path(PUFFIN_OUT, "tables/puffin_motif_smf_motif_correlations.tsv"), sep = "\t", na = "NA")
  fwrite(states, file.path(PUFFIN_OUT, "tables/puffin_motif_smf_box_state.tsv"), sep = "\t", na = "NA")
  fwrite(pred, file.path(PUFFIN_OUT, "tables/puffin_motif_smf_variability_predictors.tsv"), sep = "\t", na = "NA")
  grDevices::pdf(file.path(PUFFIN_OUT, "plots/puffin_motif_smf_features.pdf"), width = 13, height = 9, onefile = TRUE)
  on.exit(grDevices::dev.off())
  for (m in names(MOTIF_COLORS)) print(motif_profile_plot(mp, g, m))
  print(correlation_heatmap(cr))
  for (p in variability_plots(v, pred)) print(p)
}

if (sys.nframe() == 0L) main()
