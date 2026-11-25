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
# Run (after 04a and expression_variability.R): bash code/Pol2/PolII/run_gene_examples.sh \
#        code/Pol2/PolII/promoter_annotation/puffin_motif_smf_features.R
#
# Outputs under LCL_project/Pol2/PolII/promoter_annotation/
#   tables/puffin_motif_smf_genes.tsv       one row per gene: TATA box, m6A, +1, footprint features
#   tables/puffin_motif_smf_reads.tsv.gz    one row per read: TATA box state, +1 start, TATA -> +1
#   tables/puffin_motif_smf_tests.tsv       TATA vs other promoters, Wilcoxon per feature, BH
#   tables/puffin_motif_smf_profiles.tsv    mean gene profile per motif x (motif genes / all other genes)
#                                           x 10 bp bin x signal
#   tables/puffin_motif_smf_motif_tests.tsv each motif's genes vs all other genes, Wilcoxon per SMF
#                                           feature, BH over motif x feature
#   tables/puffin_motif_smf_box_state.tsv   +1 fraction / position per TATA box state
#   tables/puffin_motif_smf_variability_predictors.tsv
#                                           Spearman of each sequence predictor vs adjusted log10 CV2
#   plots/puffin_motif_smf_features.pdf     mean gene profiles of every motif (motif genes vs all other
#                                           genes) x signal in one figure
#   plots/tata_features.pdf                 CpG o/e and TATA box vs expression variability, sequence
#                                           predictors ranked by |rho|
#   tables/puffin_motif_smf_class_plot_tests.tsv
#                                           tests annotated on the motif class figures below: Kruskal-Wallis
#                                           over motif classes (BH over features) and each motif vs all
#                                           other genes (Wilcoxon, BH over motifs within a feature)
#   plots/puffin_motif_class_gene_dots.pdf  +1 start median and 10-30 bp core occupancy per gene by motif,
#                                           dots coloured by gene group
#   plots/puffin_motif_class_plus1_density.pdf
#                                           per read +1 start density per motif, every gene weighted equally
#   plots/puffin_motif_class_scatter.pdf    +1 start median vs 10-30 bp core occupancy per gene, motif medians
# Motif class figures: classes with >= MIN_CLASS_GENES genes, no adjustment for gene group or expression.
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

BASE <- "/project/spott/cshan/fiber-seq/LCL_project/Pol2/PolII"
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
MIN_CLASS_GENES <- 5L
GROUP_COLORS <- c(housekeeping = "#1B9E77", cell_type_specific = "#D95F02", not_expressed = "#7570B3")
CLASS_FEATURES <- c(plus1_start_median = "+1 start, median (bp)", plus1_fraction = "+1 fraction",
                    tf_10_30_core = "10-30 bp footprints", tf_30_60_core = "30-60 bp footprints",
                    tf_60_80_core = "60-80 bp footprints", tf_10_80_core = "10-80 bp footprints")
OTHER_COLOR <- "#5b6b7a"
MOTIF_FEATURES <- c("m6a_core", "m6a_upstream", "m6a_downstream", "nuc_core", "tf_10_30_core", "tf_30_60_core",
                    "tf_60_80_core", "occluded_fraction", "plus1_fraction", "plus1_start_median", "plus1_start_iqr")
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
# one row per gene with complete fiber features: Puffin dominant motif (motif), TATA score and
# centre from 04a, gene group and the gene-level fiber features of 02; promoter = TATA / other;
# tata_box = tata_center for TATA promoters, CONTROL_POSITION for all others
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

# per gene: +1 fraction, +1 start median / IQR (NA with < MIN_PLUS1_READS +1 reads), fraction of
# reads with the TATA box bin under a nucleosome / footprinted / open, m6A in the TATA box bin, and
# TATA box -> +1 distance, added to the load_genes() columns
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

# TATA vs other promoters for every feature in FEATURES: medians, Wilcoxon p, BH over features
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

# each motif's genes vs all other genes for every feature in MOTIF_FEATURES: medians, Wilcoxon
# p, BH over motif x feature
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

# per promoter type x TATA box state: reads, genes, +1 fraction, TSS occluded fraction and
# median +1 start, pooling reads across genes
box_state_summary <- function(reads) {
  reads[, .(n_reads = .N, n_genes = uniqueN(gene), plus1_fraction = mean(!is.na(plus1_start)),
            tss_occluded_fraction = mean(nuc_at_tss),
            plus1_start_median = as.numeric(median(plus1_start, na.rm = TRUE))),
        by = .(promoter, state)][order(promoter, state)]
}


############################################################
# sequence predictors of expression variability
############################################################
# observed / expected CpG ratio over TSS +/- CPG_FLANK bp of hg38 (vectorised over genes):
# CpG count x window length / (C count x G count)
cpg_oe <- function(chrom, tss) {
  seqs <- toupper(as.character(scanFa(FaFile(GENOME),
    GenomicRanges::GRanges(chrom, IRanges::IRanges(tss - CPG_FLANK, tss + CPG_FLANK)))))
  dna <- DNAStringSet(seqs)
  cg <- vcountPattern("CG", dna)
  base <- letterFrequency(dna, c("C", "G"))
  cg * width(dna) / (base[, "C"] * base[, "G"])
}

# genes matched to a published promoter (gene_annotations.tsv of expression_variability.R) with
# their published variability and Puffin contributions, plus Puffin class, initiator and
# trinucleotide effects, CpG o/e and a 0/1 TATA class
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

# Spearman rho of every sequence predictor in PREDICTORS with the published adjusted log10 CV2,
# BH over predictors, ranked by |rho|
variability_predictors <- function(v) {
  res <- rbindlist(lapply(names(PREDICTORS), function(f) {
    ct <- suppressWarnings(cor.test(v[[f]], v$elife_adjusted_log10_CV2, method = "spearman", exact = FALSE))
    list(predictor = f, label = PREDICTORS[[f]], n = sum(is.finite(v[[f]])), rho = unname(ct$estimate), p = ct$p.value)
  }))
  res[, fdr := p.adjust(p, "BH")]
  res[order(-abs(rho))][, rank := .I][]
}

# the three tata_features.pdf pages: CpG o/e vs variability (dots coloured by TATA / other),
# variability of TATA vs other promoters (box + dots), and every predictor's rho as a ranked bar
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
# motif class figures (no adjustment for gene group or expression)
############################################################
# genes of motif classes with >= MIN_CLASS_GENES genes, motifs ordered by median +1 start;
# tf_10_80_core = sum of the three size classes
class_genes <- function(g) {
  x <- copy(g)[, tf_10_80_core := tf_10_30_core + tf_30_60_core + tf_60_80_core]
  keep <- x[, .N, by = motif][N >= MIN_CLASS_GENES, motif]
  ord <- x[motif %in% keep, .(m = median(plus1_start_median, na.rm = TRUE)), by = motif][order(m), motif]
  list(all = x, genes = x[motif %in% keep][, motif := factor(motif, ord)][])
}

# Kruskal-Wallis over classes (BH over features); each motif vs all other genes of g
# (Wilcoxon, BH over motifs within a feature)
class_plot_tests <- function(cg) {
  x <- cg$genes; all <- cg$all
  kw <- rbindlist(lapply(names(CLASS_FEATURES), function(f) {
    y <- x[!is.na(get(f))]
    list(test = "kruskal_wallis", feature = f, motif = NA_character_, n_motif = nrow(y), n_other = NA_integer_,
         median_motif = NA_real_, median_other = NA_real_, p = kruskal.test(y[[f]], droplevels(y$motif))$p.value)
  }))[, fdr := p.adjust(p, "BH")]
  vs <- rbindlist(lapply(names(CLASS_FEATURES), function(f) rbindlist(lapply(levels(x$motif), function(m) {
    a <- all[motif == m & !is.na(get(f))][[f]]; o <- all[motif != m & !is.na(get(f))][[f]]
    list(test = "motif_vs_rest", feature = f, motif = m, n_motif = length(a), n_other = length(o),
         median_motif = median(a), median_other = median(o), p = suppressWarnings(wilcox.test(a, o))$p.value)
  }))))[, fdr := p.adjust(p, "BH"), by = feature]
  rbind(kw, vs)
}

# 1. one dot per gene (coloured by gene group) over a box per motif
class_gene_dots <- function(cg, tests) {
  f <- c("plus1_start_median", "tf_10_30_core")
  d <- melt(cg$genes, id.vars = c("gene_name", "motif", "group"), measure.vars = f, variable.name = "feature")[!is.na(value)]
  d[, feature := factor(CLASS_FEATURES[as.character(feature)], CLASS_FEATURES[f])]
  n <- cg$genes[, .N, by = motif]
  xl <- setNames(sprintf("%s
n=%d", n$motif, n$N), as.character(n$motif))
  top <- d[, .(y = max(value)), by = feature]
  lab <- merge(tests[test == "motif_vs_rest" & feature %in% f], data.table(feature = f), by = "feature")
  lab[, `:=`(motif = factor(motif, levels(cg$genes$motif)), feature = factor(CLASS_FEATURES[feature], CLASS_FEATURES[f]))]
  lab <- merge(lab, top, by = "feature")[, txt := sprintf("q=%.2g", fdr)]
  kw <- tests[test == "kruskal_wallis" & feature %in% f]
  ggplot(d, aes(motif, value)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, colour = "grey55", fill = "grey95") +
    geom_jitter(aes(colour = group), width = 0.18, height = 0, size = 1.6, alpha = 0.8) +
    geom_text(data = lab, aes(motif, y * 1.08, label = txt), size = 2.8) +
    facet_wrap(~ feature, scales = "free_y", nrow = 1) +
    scale_x_discrete(labels = xl) +
    scale_colour_manual(values = GROUP_COLORS, name = "gene group") +
    labs(x = "Puffin dominant motif (ordered by median +1 start)", y = NULL,
         title = "Promoter features by Puffin dominant motif (one dot per gene, not adjusted for gene group)",
         subtitle = sprintf("Kruskal-Wallis: +1 start FDR = %.2g, 10-30 bp footprints (TSS +/- 200) FDR = %.2g; q = each motif vs all other genes, Wilcoxon BH",
                            kw[feature == f[1], fdr], kw[feature == f[2], fdr])) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
}

# mean gene profile per motif (columns, in `motifs` order) vs all other genes (grey) for each of
# `signals` (rows), from motif_profiles(); dotted red line in the TATA column = median TATA box
motif_profile_grid <- function(mp, g, motifs, signals, title) {
  d <- mp[motif %in% motifs & signal %in% signals & position >= PROFILE_XLIM[1] & position <= PROFILE_XLIM[2]]
  d[, signal := factor(SIGNALS[as.character(signal)], SIGNALS[signals])]
  d[, motif := factor(motif, motifs)]
  d[, line := factor(fifelse(set == OTHER_GENES, OTHER_GENES, "motif genes"), c("motif genes", OTHER_GENES))]
  n <- g[motif %in% motifs, .N, by = motif]
  ml <- setNames(sprintf("%s (n=%d)", n$motif, n$N), as.character(n$motif))
  box <- if ("TATA" %in% motifs) data.table(motif = factor("TATA", motifs), x = g[motif == "TATA", median(tata_box)]) else NULL
  ggplot(d[line == OTHER_GENES], aes(position, mean)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    (if (!is.null(box)) geom_vline(data = box, aes(xintercept = x), linetype = 3, colour = MOTIF_COLORS[["TATA"]])) +
    geom_line(colour = "grey60", linewidth = 0.4) +
    geom_ribbon(data = d[line != OTHER_GENES], aes(ymin = mean - se, ymax = mean + se, fill = motif), alpha = 0.25) +
    geom_line(data = d[line != OTHER_GENES], aes(colour = motif), linewidth = 0.6) +
    facet_grid(signal ~ motif, scales = "free_y", switch = "y", labeller = labeller(motif = ml)) +
    scale_colour_manual(values = MOTIF_COLORS, guide = "none") + scale_fill_manual(values = MOTIF_COLORS, guide = "none") +
    scale_x_continuous(breaks = c(-500, 0, 500)) +
    labs(x = "Position relative to TSS (bp, 10 bp bins)", y = NULL, title = title) +
    theme_bw(base_size = 10) + theme(strip.background = element_blank(), strip.placement = "outside")
}

# 3. per read +1 start density per motif; genes with >= MIN_PLUS1_READS +1 reads, each gene weight 1
class_plus1_density <- function(reads, cg) {
  x <- merge(reads[!is.na(plus1_start), .(gene, plus1_start)], cg$genes[, .(gene = gene_name, motif)], by = "gene")
  x[, n := .N, by = gene]
  x <- x[n >= MIN_PLUS1_READS][, w := 1 / n][, w := w / sum(w), by = motif]
  med <- cg$genes[, .(m = median(plus1_start_median, na.rm = TRUE), N = sum(!is.na(plus1_start_median))), by = motif]
  ml <- setNames(sprintf("%s (n=%d genes)", med$motif, med$N), as.character(med$motif))
  ggplot(x, aes(plus1_start, weight = w, fill = motif, colour = motif)) +
    geom_density(alpha = 0.35, adjust = 0.6, linewidth = 0.4) +
    geom_vline(data = med, aes(xintercept = m, colour = motif), linetype = 2) +
    facet_wrap(~ motif, ncol = 1, strip.position = "left", labeller = labeller(motif = ml)) +
    scale_colour_manual(values = MOTIF_COLORS, guide = "none") + scale_fill_manual(values = MOTIF_COLORS, guide = "none") +
    labs(x = "+1 nucleosome start (bp from TSS; reads with an open TSS)", y = NULL,
         title = "Per read +1 start by Puffin dominant motif") +
    theme_bw(base_size = 10) +
    theme(strip.background = element_blank(), strip.placement = "outside", strip.text.y.left = element_text(angle = 0),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

# 5. +1 start median vs 10-30 bp core occupancy, one dot per gene, large point = motif median
class_scatter <- function(cg) {
  x <- cg$genes[!is.na(plus1_start_median)]
  med <- x[, .(plus1_start_median = median(plus1_start_median), tf_10_30_core = median(tf_10_30_core), N = .N), by = motif]
  ggplot(x, aes(plus1_start_median, tf_10_30_core, colour = motif)) +
    geom_point(size = 1.4, alpha = 0.45) +
    geom_point(data = med, size = 5, shape = 21, fill = "white", stroke = 1.4) +
    ggrepel::geom_text_repel(data = med, aes(label = sprintf("%s (n=%d)", motif, N)), size = 3.5,
                             box.padding = 1, point.padding = 0.8, min.segment.length = 0,
                             show.legend = FALSE, seed = 1) +
    scale_colour_manual(values = MOTIF_COLORS, name = NULL) +
    labs(x = CLASS_FEATURES[["plus1_start_median"]], y = "10-30 bp footprint occupancy (TSS +/- 200 bp)",
         title = "+1 start vs 10-30 bp footprint occupancy per gene; large circles = motif medians") +
    theme_bw(base_size = 11)
}

# 1. load genes, reads and gene profiles; summarise per gene; run the TATA and motif tests
# 2. load the published variability and test the sequence predictors
# 3. write the tables, then tata_features.pdf
# 4. motif classes with >= MIN_CLASS_GENES genes: class tests, the all-motif profile grid
#    (puffin_motif_smf_features.pdf), gene dots, +1 density and scatter PDFs
main <- function() {
  genes <- load_genes()
  reads <- load_reads(genes)
  prof <- load_profiles(genes)
  g <- gene_summary(genes, reads, prof)
  tests <- tata_tests(g)
  mp <- motif_profiles(prof)
  motif_tests <- motif_class_tests(g)
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
  fwrite(states, file.path(PUFFIN_OUT, "tables/puffin_motif_smf_box_state.tsv"), sep = "\t", na = "NA")
  fwrite(pred, file.path(PUFFIN_OUT, "tables/puffin_motif_smf_variability_predictors.tsv"), sep = "\t", na = "NA")
  grDevices::pdf(file.path(PUFFIN_OUT, "plots/tata_features.pdf"), width = 13, height = 9, onefile = TRUE)
  for (p in variability_plots(v, pred)) print(p)
  grDevices::dev.off()

  cg <- class_genes(g)
  ct <- class_plot_tests(cg)
  fwrite(ct, file.path(PUFFIN_OUT, "tables/puffin_motif_smf_class_plot_tests.tsv"), sep = "\t", na = "NA")
  nm <- nlevels(cg$genes$motif)
  all_motifs <- c(levels(cg$genes$motif), setdiff(names(MOTIF_COLORS), levels(cg$genes$motif)))
  ggsave(file.path(PUFFIN_OUT, "plots/puffin_motif_smf_features.pdf"),
         motif_profile_grid(mp, g, all_motifs, names(SIGNALS),
           "Mean gene profile per Puffin dominant motif (+/- se); grey = all other genes; dotted red = median TATA box"),
         width = 16, height = 10)
  ggsave(file.path(PUFFIN_OUT, "plots/puffin_motif_class_gene_dots.pdf"), class_gene_dots(cg, ct), width = 12, height = 5.5)
  ggsave(file.path(PUFFIN_OUT, "plots/puffin_motif_class_plus1_density.pdf"), class_plus1_density(reads, cg),
         width = 7, height = 1 + 1.1 * nm)
  ggsave(file.path(PUFFIN_OUT, "plots/puffin_motif_class_scatter.pdf"), class_scatter(cg), width = 8, height = 6.5)
}

if (sys.nframe() == 0L) main()
