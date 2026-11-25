#!/usr/bin/env Rscript
# Cluster level and gene level summary of the curated PolII example genes, in one PDF.
# Reads only tables: the per read x 10 bp bin counts of 02_extract_read_features.R and the
# family / cluster labels of 01_summarize_recurrent_patterns.R. Nothing is re-extracted and
# nothing is re-clustered.
#
# Run: bash code/Pol2/PolII/run_gene_examples.sh \
#        code/Pol2/PolII/recurrent_cluster_features/03_cluster_gene_summary.R
# Needs 02_extract_read_features.R to have been run first.
#
# Outputs under LCL_project/Pol2/PolII/recurrent_patterns/
#   tables/03_cluster_bin_profiles.tsv.gz    one row per original cluster x bin
#   tables/03_gene_bin_profiles.tsv.gz       one row per gene x bin
#   tables/03_archetype_bin_profiles.tsv.gz  one row per gene x archetype x bin
#   tables/03_gene_summary.tsv               one row per gene
#   plots/03_cluster_gene_summary.pdf        four pages
#
# Every rate is a pooled count ratio: m6A = calls / (reads x observed sites in the bin),
# occupancy = covered bases / (reads x bases in the bin)

# Archetypes are the per read occlusion / +1 classes of 02_extract_read_features.R

source("/project/spott/cshan/fiber-seq/code/Pol2/PolII/recurrent_cluster_features/01_summarize_recurrent_patterns.R")
source("/project/spott/cshan/fiber-seq/code/Pol2/PolII/02_extract_read_features.R")
setDTthreads(1L)

SIGNALS <- c(m6a = "m6A fraction (observed sites)", nuc = "nucleosome occupancy",
             setNames(paste("TF footprint", TF_CLASSES, "occupancy"), TF_COLS))
FAMILIES <- paste0("F", 1:6)
FAMILY_COLORS <- setNames(cluster_id_palette(paste0("cluster", 1:6)), FAMILIES)
ARCHETYPE_COLORS <- setNames(c("#4a9b3f", "#1f6fb4", "#e8871a", "#d0241f"), ARCHETYPES)
MIN_ARCHETYPE_READS <- 5L
ARCHETYPE_XLIM <- c(-400, 600)


############################################################
# one gene's 10 bp bin counts, pooled per cluster, per gene and per archetype
############################################################
# 1. read bins/<gene>.tsv.gz and the gene's read level group / archetype labels
# 2. attach the original cluster and family labels from 01, keeping reads of every cluster
# 3. sum calls, observed sites, covered bases and bases per bin, three times:
#    per cluster (with its family / eligible flag), over all reads of the gene, and per archetype
# 4. divide once to get the m6A fraction and the four occupancy rates
# 5. return the three tables: list(cluster, gene, archetype)
gene_profiles <- function(gene_id, reads, cluster_dir) {
  bins <- fread(file.path(ROOT, "bins", paste0(gene_id, ".tsv.gz")))
  labels <- fread(file.path(cluster_dir, gene_id, "read_features.tsv.gz"),
                  select = c("RID", "cluster", "family", "eligible"))
  r <- reads[gene == gene_id, .(RID, group, archetype)]
  stopifnot(setequal(bins$RID, r$RID), setequal(labels$RID, r$RID))
  d <- merge(merge(bins, r, by = "RID"), labels, by = "RID")
  rate <- function(by) {
    s <- d[, c(.(n_reads = uniqueN(RID), m6a_calls = sum(m6a_calls), m6a_sites = sum(m6a_sites),
                 bin_bases = sum(bin_bases)),
               lapply(.SD, sum)), by = c(by, "bin", "position"), .SDcols = paste0(COVER_COLS, "_bases")]
    s[, m6a := fifelse(m6a_sites > 0, m6a_calls / m6a_sites, NA_real_)]
    for (nm in COVER_COLS) set(s, j = nm, value = s[[paste0(nm, "_bases")]] / s$bin_bases)
    s[, c(by, "bin", "position", "n_reads", names(SIGNALS)), with = FALSE]
  }
  list(cluster = rate(c("gene", "group", "cluster", "family", "eligible")),
       gene = rate(c("gene", "group")),
       archetype = rate(c("gene", "group", "archetype")))
}


############################################################
# average bin profiles into one curve per family, gene group or archetype
############################################################

average_profiles <- function(d, key, min_reads = 1L)
  d[n_reads >= min_reads, lapply(.SD, mean, na.rm = TRUE), by = c(key, "bin", "position"),
    .SDcols = names(SIGNALS)]


############################################################
# line plot of one signal against bp from the TSS, one line per family / group / archetype
############################################################

profile_panel <- function(d, key, colors, signal, xlim = EXAMPLE_RANGE, ylim = NULL) {
  d <- d[position >= xlim[1] & position <= xlim[2]]
  if (is.null(ylim)) ylim <- range(d[[signal]], na.rm = TRUE)
  plot(NA, xlim = xlim, ylim = ylim, xlab = "bp from TSS", ylab = SIGNALS[[signal]],
       main = SIGNALS[[signal]], bty = "l", cex.main = 0.95)
  abline(v = 0, col = "grey60", lty = 2)
  for (k in names(colors)) {
    x <- d[get(key) == k][order(position)]
    if (nrow(x)) lines(x$position, x[[signal]], col = colors[[k]], lwd = 1.6)
  }
}


############################################################
# violin plot of one per gene value (e.g. occluded-TSS fraction), one violin per gene group
############################################################

violin_panel <- function(values, groups, colors, main, ylab) {
  levels_ <- names(colors)
  plot(NA, xlim = c(0.5, length(levels_) + 0.5), ylim = c(0, max(values, na.rm = TRUE) * 1.15),
       xaxt = "n", xlab = "", ylab = ylab, main = main, bty = "l")
  axis(1, at = seq_along(levels_), labels = levels_, cex.axis = 0.85)
  for (i in seq_along(levels_)) {
    v <- values[groups == levels_[i]]
    if (length(v) > 1L) {
      d <- density(v, from = max(0, min(v)), to = max(v))
      w <- 0.4 * d$y / max(d$y)
      polygon(c(i + w, rev(i - w)), c(d$x, rev(d$x)), col = adjustcolor(colors[[i]], 0.35), border = colors[[i]])
    }
    points(jitter(rep(i, length(v)), amount = 0.12), v, pch = 16, cex = 0.4,
           col = adjustcolor(colors[[i]], 0.6))
    segments(i - 0.4, median(v), i + 0.4, median(v), lwd = 2)
    text(i, max(values, na.rm = TRUE) * 1.1, sprintf("n=%d", length(v)), cex = 0.75)
  }
}


############################################################
# draw the four page PDF: family profiles, gene group profiles, occluded TSS, archetypes
############################################################
# 1. by family: m6A, nucleosome and three TF footprint profiles, eligible clusters averaged equally
# 2. by gene group: the same five signals, genes averaged equally
# 3. violins of the per gene fraction of fibers with a nucleosome over the TSS, by gene group
# 4. archetype profiles: open (A1 vs A2) and occluded (A3 vs A4) promoters, with vs without
#    a +1 nucleosome, m6A solid, nucleosome dotted
plot_summary_pages <- function(family, group, archetype, genes, path) {
  grDevices::pdf(path, width = 13, height = 7.5, onefile = TRUE)
  on.exit(grDevices::dev.off())

  for (page in list(list(d = family, key = "family", colors = FAMILY_COLORS,
                         head = "Profile families: eligible clusters averaged with equal weight"),
                    list(d = group, key = "group", colors = COLORS,
                         head = "Gene categories: genes averaged with equal weight"))) {
    par(mfrow = c(2, 3), mar = c(4, 4.2, 3, 1), oma = c(0, 0, 2.4, 0))
    for (s in names(SIGNALS)) profile_panel(page$d, page$key, page$colors, s)
    plot.new(); legend("center", names(page$colors), col = unlist(page$colors), lty = 1, lwd = 2, bty = "n")
    mtext(page$head, outer = TRUE, cex = 0.95)
  }

  par(mfrow = c(1, 1), mar = c(4, 5, 4, 2), oma = c(0, 0, 0, 0))
  violin_panel(genes$occluded_fraction, genes$group, COLORS,
               "Fibers with an occluded TSS", "fraction of fibers, nucleosome over the TSS base")

  par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1), oma = c(0, 0, 2.4, 0))
  for (pair in list(ARCHETYPES[1:2], ARCHETYPES[3:4])) {
    d <- archetype[archetype %in% pair & position >= ARCHETYPE_XLIM[1] & position <= ARCHETYPE_XLIM[2]]
    plot(NA, xlim = ARCHETYPE_XLIM, ylim = range(c(d$m6a, d$nuc), na.rm = TRUE), bty = "l",
         xlab = "bp from TSS", ylab = "solid m6A / dotted nucleosome",
         main = sprintf("%s promoters +/- +1", sub(".* ", "", sub(" /.*", "", pair[1]))))
    abline(v = 0, lty = 2)
    for (k in pair) {
      x <- d[archetype == k][order(position)]
      lines(x$position, x$m6a, col = ARCHETYPE_COLORS[[k]], lwd = 1.8)
      lines(x$position, x$nuc, col = ARCHETYPE_COLORS[[k]], lwd = 1.1, lty = 3)
    }
    legend("topright", pair, col = ARCHETYPE_COLORS[pair], lty = 1, lwd = 2, bty = "n", cex = 0.8)
  }
  mtext(sprintf("Per fiber classes; a gene joins a curve with at least %d reads in it", MIN_ARCHETYPE_READS),
        outer = TRUE, cex = 0.95)
}


############################################################
# run all genes: build the per cluster / gene / archetype profiles, write the tables, draw the PDF
############################################################
# 1. read the extracted read features and keep the genes 01 also reconstructed
# 2. gene_profiles() per gene, stacked into cluster, gene and archetype bin profile tables
# 3. per gene summary: read count, occluded-TSS fraction, +1 fraction, archetype fractions
# 4. write the four tables

main <- function() {
  cluster_dir <- file.path(OUT, "tables")
  reads <- fread(file.path(ROOT, "tables/read_features.tsv.gz"))
  done <- intersect(unique(reads$gene), list.dirs(cluster_dir, recursive = FALSE, full.names = FALSE))
  stopifnot(length(done) > 0L)
  message("Summarising ", length(done), " genes")
  parts <- lapply(done, gene_profiles, reads = reads, cluster_dir = cluster_dir)
  cluster_prof <- rbindlist(lapply(parts, `[[`, "cluster"))
  gene_prof <- rbindlist(lapply(parts, `[[`, "gene"))
  arch_prof <- rbindlist(lapply(parts, `[[`, "archetype"))

  genes <- reads[gene %in% done, .(n_reads = .N, occluded_fraction = mean(nuc_at_tss),
                                   plus1_fraction = mean(!is.na(plus1_start))), by = .(gene, group)]
  for (k in ARCHETYPES) {
    q <- reads[gene %in% done, .(f = mean(archetype == k)), by = gene]
    set(genes, j = paste0("fraction_", sub(" .*", "", k)), value = q$f[match(genes$gene, q$gene)])
  }
  genes[, group := factor(group, levels = GROUPS)]

  dir.create(file.path(OUT, "plots"), recursive = TRUE, showWarnings = FALSE)
  fwrite(cluster_prof, file.path(cluster_dir, "03_cluster_bin_profiles.tsv.gz"), sep = "\t", na = "NA")
  fwrite(gene_prof, file.path(cluster_dir, "03_gene_bin_profiles.tsv.gz"), sep = "\t", na = "NA")
  fwrite(arch_prof, file.path(cluster_dir, "03_archetype_bin_profiles.tsv.gz"), sep = "\t", na = "NA")
  fwrite(genes, file.path(cluster_dir, "03_gene_summary.tsv"), sep = "\t", na = "NA")

  plot_summary_pages(
    average_profiles(cluster_prof[eligible == TRUE & !is.na(family)], "family"),
    average_profiles(gene_prof, "group"),
    average_profiles(arch_prof[!is.na(archetype)], "archetype", MIN_ARCHETYPE_READS),
    genes, file.path(OUT, "plots/03_cluster_gene_summary.pdf"))
  message("Complete: ", nrow(cluster_prof), " cluster rows; ", nrow(gene_prof), " gene rows; ",
          nrow(arch_prof), " archetype rows; ", nrow(genes), " genes.")
}
if (sys.nframe() == 0L) main()
