# Plotting functions for code/pausing/PolII.
#
# plot_tss_example(gene): Tullius et al. Figure 1F for one gene at its FANTOM5 LCL dominant TSS, an R port of Data S1 cell 54
# (code/pausing/example_code/mmc3.zip).
#   top     GRO-seq (sense strand, standing in for PRO-seq), CAGE (sense) and DNase-seq, each
#           scaled to its maximum in the window. There is no GM12878 MNase-seq here; DNase marks
#           the NDR, i.e. the inverse of the MNase nucleosome signal. Ticks = called peaks.
#   middle  fraction of the plotted reads with a PPP, PIC or nucleosome call of the bottom panel
#           at each base, SMOOTH_BP running mean, counting peak-overlapping footprints in the size range
#   bottom  every read of the 31 LCL samples with a fibertools nucleosome record (ft extract --nuc,
#           window_nucleosomes()) in TSS +/- 1 kb that has a +1 nucleosome, sorted by its position
#           (reads without one are dropped, as in Tullius); no other read filtering. Nucleosomes
#           are the ft calls; TF footprints are the FiberHMM combined 10-80 bp footprints of those
#           reads (window_tf_footprints()), drawn on top: PPP = 30-60 bp, PIC = 60-80 bp,
#           only when overlapping a GRO-seq / CAGE peak, respectively (60 bp is 60-80);
#           other length-matched footprints are unknown (grey), every
#           other footprint its assign_size_class() TF class, coloured as in
#           MOLECULE_TRACK_COLORS. Grey line = the read's span without a footprint, as in
#           plot_molecules().
# Peaks follow Tullius call_peaks(): every sliding window of the search range whose maximum is
# above PEAK_MIN_COUNT and above 1/PEAK_REL of the range maximum contributes its first argmax.
# GRO-seq counts RNA fragment ends, so its peaks are coarser than PRO-seq's.
#
# Any chr1-22 gene of pause_index.tsv can be plotted. rtracklayer needs
# /project/spott/cshan/envs/dimelo/lib/libcurl.so.4 loaded first (dyn.load(..., local = FALSE), as
# in gene_pausing_example.rmd, which chooses the examples and plots them).
# Usage: source("/project/spott/cshan/fiber-seq/code/pausing/PolII/polII_pausing_plots.r")
#        plot_tss_example("GADD45B")   # -> plots/example_reads/GADD45B.pdf

source("/project/spott/cshan/fiber-seq/code/pausing/PolII/pol2_functions.R")
suppressPackageStartupMessages({ library(rtracklayer); library(ggplot2); library(patchwork) })
source(PLOTTING)   # theme_fiberseq(), save_figure(), MOLECULE_TRACK_COLORS
source(PARSING)    # assign_size_class(), read_tabix_region(), convert_ft_bed12_to_bed6()

EXAMPLE_RANGE  <- c(-400L, 400L)   # Tullius Fig 1F window
PEAK_SEARCH    <- list(groseq = c(from = 20L, to = 349L, window = 20L),   # Data S1 cell 53
                       cage   = c(from = -50L, to = 49L, window = 30L))
PEAK_MIN_COUNT <- 5
PEAK_REL       <- 20
READ_ROW_IN    <- 0.025            # PDF inches per read
TRACK_COLORS   <- c(`GRO-seq` = CLASS_COLORS[["PPP"]], CAGE = CLASS_COLORS[["PIC"]],
                    DNase = INK[["secondary"]])
TF_CLASSES     <- setdiff(levels(assign_size_class(integer(), TF_BREAKS)), unname(CALL_CLASS))
IDENTITY_COLORS <- c(CLASS_COLORS[c("PPP", "PIC")],
                     setNames(MOLECULE_TRACK_COLORS[paste("TF footprint", TF_CLASSES)], TF_CLASSES),
                     unknown = "grey65",
                     nucleosome = MOLECULE_TRACK_COLORS[["nucleosome 130-160 bp"]])
READ_LINE      <- "grey80"         # read backbone of plot_molecules()


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

# GRO-seq (sense, lifted to hg38), CAGE consensus (sense) and DNase at TSS -FLANK..FLANK
tss_tracks <- function(g) {
  offsets <- -FLANK:FLANK
  win <- GRanges(g$chrom, IRanges(g$tss - FLANK, g$tss + FLANK))
  gro <- subsetByOverlaps(lift_groseq(g$chrom, import.chain(CHAIN))[[g$strand]], win)
  cage <- fread(CAGE, header = FALSE, select = c(1:3, 5:6))[
    V1 == g$chrom & V6 == g$strand & V3 >= start(win) & V2 < end(win)]
  dnase <- import(DNASE, which = win)
  data.table(pos = offsets,
             `GRO-seq` = signal_at(start(gro) - 1L, end(gro), gro$score, g, offsets),
             CAGE = signal_at(cage$V2, cage$V3, cage$V5, g, offsets),
             DNase = signal_at(start(dnase) - 1L, end(dnase), dnase$score, g, offsets))
}

# Tullius call_peaks() on per-base signal x at offsets; sorted peak offsets
call_signal_peaks <- function(x, offsets, search) {
  keep <- offsets >= search[["from"]] & offsets <= search[["to"]]
  x <- x[keep]; pos <- offsets[keep]
  floor_rel <- max(x) / PEAK_REL
  peaks <- vapply(seq_along(x), function(i) {
    w <- x[i:min(length(x), i + search[["window"]] - 1L)]
    if (max(w) > PEAK_MIN_COUNT && max(w) > floor_rel) pos[i + which.max(w) - 1L] else NA_integer_
  }, integer(1))
  sort(unique(peaks[!is.na(peaks)]))
}

# Use the FANTOM5 dominant 1-bp TSS, matching stable ENSG first, then symbol.
# BED end is the 1-based TSS on either strand. Never silently fall back to RefSeq.
fantom_tss_gene <- function(g) {
  dom <- fread(DOMINANT_TSS)
  dom <- dom[`#chrom` == g$chrom & strand == g$strand]
  hit <- dom[sub("\\..*", "", gene_id) == g$ensg]
  if (!nrow(hit)) hit <- dom[gene_name == g$gene_name]
  if (nrow(hit) != 1L) stop(g$gene_name, ": expected one FANTOM5 dominant TSS; found ", nrow(hit))
  g <- copy(g)
  g[, `:=`(refseq_tss = tss, tss = as.integer(hit$end))]
  g
}

# A footprint must contain a called peak base and match the corresponding size class.
# Peaks are point coordinates from call_signal_peaks(), not broad peak intervals.
footprint_identity <- function(fps, peaks) {
  cls <- as.character(assign_size_class(fps$size, TF_BREAKS))
  out <- cls
  out[cls %in% unname(CALL_CLASS)] <- "unknown"
  overlaps_peak <- function(pos) {
    vapply(seq_len(nrow(fps)), function(i) any(pos >= fps$lo[i] & pos <= fps$hi[i]), logical(1))
  }
  out[cls == CALL_CLASS[["PPP"]] & overlaps_peak(peaks$groseq)] <- "PPP"
  out[cls == CALL_CLASS[["PIC"]] & overlaps_peak(peaks$cage)] <- "PIC"
  if (anyNA(out)) stop("footprint sizes outside ", min(TF_BREAKS), "-", max(TF_BREAKS), " bp: ",
                       paste(unique(fps$size[is.na(out)]), collapse = ", "))
  out
}

# Figure 1F for one chr1-22 gene of pause_index.tsv; writes plots/example_reads/<gene>.pdf,
# returns the plot
plot_tss_example <- function(gene) {
  g <- fread(PI_TABLE)[gene_name == gene]
  if (nrow(g) != 1L) stop(gene, ": ", nrow(g), " genes of that name in ", PI_TABLE)
  if (!g$chrom %in% AUTOSOMES) stop(gene, " is on ", g$chrom, "; FiberHMM covers chr1-22 only")
  g <- fantom_tss_gene(g)
  r1 <- EXAMPLE_RANGE[1]; r2 <- EXAMPLE_RANGE[2]

  # reads, their nucleosome calls and FiberHMM footprints at TSS +/- FLANK, TSS-relative
  relative <- function(d) {
    r <- to_relative(d$start0, d$end0, g$tss, rep(g$strand, nrow(d)))
    d[, `:=`(lo = r$lo, hi = r$hi, size = end0 - start0)][]
  }
  nuc   <- window_nucleosomes(g$chrom, g$tss - FLANK, g$tss + FLANK)
  reads <- relative(nuc$reads)
  nucs  <- relative(nuc$nucs)
  fp    <- relative(window_tf_footprints(g$chrom, g$tss - FLANK, g$tss + FLANK)[RID %in% reads$RID])

  # reads with a +1 nucleosome, nearest the TSS on top
  plus1 <- plus1_offset(footprint_matrix(reads, nucs))
  keep  <- !is.na(plus1)
  rids  <- reads$RID[keep][order(plus1[keep], reads$RID[keep])]
  n     <- length(rids)
  rows  <- data.table(RID = rids, y = rev(seq_len(n)))
  cover <- footprint_matrix(rows, reads[, .(RID, lo, hi, size = 1L)])

  tr <- tss_tracks(g)
  peaks <- list(groseq = call_signal_peaks(tr$`GRO-seq`, tr$pos, PEAK_SEARCH$groseq),
                cage   = call_signal_peaks(tr$CAGE, tr$pos, PEAK_SEARCH$cage))
  fp <- fp[RID %in% rids]
  fp[, identity := footprint_identity(.SD, peaks)]
  nucs <- nucs[RID %in% rids][, identity := "nucleosome"]

  # middle panel from the bottom panel's calls, over the reads covering each base
  calls <- c("PPP", "PIC", "nucleosome")
  enr <- rbindlist(lapply(calls, function(k) data.table(
    identity = k, pos = -FLANK:FLANK,
    frac = running_mean(call_fraction(
      footprint_matrix(rows, if (k == "nucleosome") nucs else fp[identity == k]), cover)))))

  # bottom panel: read spans, nucleosome calls, then FiberHMM footprints on top
  in_view <- function(d) merge(d[hi >= r1 & lo <= r2], rows, by = "RID")[
    , `:=`(x0 = pmax(lo, r1) - 0.5, x1 = pmin(hi, r2) + 0.5)][]
  spans <- in_view(reads)
  fp <- rbind(in_view(nucs), in_view(fp))

  # --- panels ---
  xs <- list(scale_x_continuous(expand = c(0, 0), breaks = seq(r1, r2, 100)),
             coord_cartesian(xlim = c(r1 - 0.5, r2 + 0.5)))
  tss_line <- geom_vline(xintercept = 0, colour = INK[["grid"]], linewidth = 0.3)
  theme_rows <- theme_fiberseq(overrides = list(
    strip.background = element_blank(), strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.x = element_blank(), axis.line.x = element_blank()))

  tl <- melt(tr[pos >= r1 & pos <= r2], id.vars = "pos", variable.name = "track")
  tl[, value := if (max(value) > 0) value / max(value) else value, by = track]
  tl[, track := factor(track, levels = names(TRACK_COLORS))]
  ticks <- data.table(track = factor(c(rep("GRO-seq", length(peaks$groseq)), rep("CAGE", length(peaks$cage))),
                                     levels = names(TRACK_COLORS)),
                      pos = c(peaks$groseq, peaks$cage))
  p_tracks <- ggplot(tl, aes(pos, value, colour = track)) + tss_line +
    geom_line(linewidth = 0.4) +
    geom_rug(data = ticks, aes(x = pos), inherit.aes = FALSE, sides = "b",
             colour = INK[["primary"]], linewidth = 0.3, length = unit(0.08, "npc")) +
    facet_grid(track ~ ., switch = "y",
               labeller = as_labeller(c(`GRO-seq` = "GRO-seq (sense)", CAGE = "CAGE (sense)",
                                        DNase = "DNase"))) +
    scale_colour_manual(values = TRACK_COLORS, guide = "none") +
    scale_y_continuous(breaks = c(0, 1), position = "right") + xs +
    labs(y = "scaled to max") + theme_rows

  el <- enr[pos >= r1 & pos <= r2]
  el[, identity := factor(identity, levels = calls)]
  p_enr <- ggplot(el, aes(pos, frac, colour = identity)) + tss_line +
    geom_line(linewidth = 0.4, na.rm = TRUE) +
    facet_grid(identity ~ ., switch = "y", scales = "free_y",
               labeller = as_labeller(c(PPP = "PPP calls", PIC = "PIC calls",
                                        nucleosome = "nucleosome calls"))) +
    scale_colour_manual(values = IDENTITY_COLORS, guide = "none") +
    scale_y_continuous(n.breaks = 3, position = "right") + xs +
    labs(y = "fraction of reads") + theme_rows

  identity_labels <- c(PPP = sprintf("PPP (%s)", CALL_CLASS[["PPP"]]),
                       PIC = sprintf("PIC (%s)", CALL_CLASS[["PIC"]]),
                       setNames(paste("TF footprint", TF_CLASSES), TF_CLASSES),
                       unknown = "unknown", nucleosome = "nucleosome")
  p_reads <- ggplot() +
    geom_segment(data = spans, aes(x = x0, xend = x1, y = y, yend = y),
                 colour = READ_LINE, linewidth = 0.15) +
    geom_rect(data = fp, aes(xmin = x0, xmax = x1, ymin = y - 0.4, ymax = y + 0.4, fill = identity)) +
    scale_fill_manual(values = IDENTITY_COLORS, limits = names(IDENTITY_COLORS),
                      labels = identity_labels, name = NULL, guide = guide_legend(nrow = 2)) +
    scale_y_continuous(expand = c(0, 0.5)) + xs +
    labs(x = "bp from the TSS (transcription ->)",
         y = sprintf("%d reads, by +1 nucleosome", n)) +
    theme_fiberseq(overrides = list(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                    axis.line.y = element_blank(), legend.position = "top"))

  tertile <- if (is.na(g$pi_bin)) "not in the analysis set" else paste(g$pi_bin, "tertile")
  p <- (p_tracks / p_enr / p_reads) +
    plot_layout(heights = c(1.2, 1.2, n * READ_ROW_IN)) +
    plot_annotation(
      title = sprintf("%s  %s:%d (%s strand), PI %.1f (%s), FPKM %.1f",
                      gene, g$chrom, g$tss, g$strand, g$pi, tertile, g$fpkm),
      subtitle = sprintf(paste0("%d of %d reads (31 LCLs pooled, ft nucleosome record in TSS +/- %d bp) have a ",
                                "+1 nucleosome in +%d..+%d and are shown\n%d GRO-seq and %d CAGE peaks (ticks)"),
                         n, nrow(reads), FLANK, min(PLUS1_OFFSETS), max(PLUS1_OFFSETS),
                         length(peaks$groseq), length(peaks$cage)),
      caption = paste0("Nucleosomes: fibertools; TF footprints: FiberHMM v1 combined tracks (10-80 bp). ",
                       "Grey line = read span without a footprint.\n",
                       "Tracks scaled to their window maximum; DNase stands in for MNase-seq (none for GM12878). ",
                       "TSS = FANTOM5 LCL dominant TSS. PPP/PIC require GRO-seq/CAGE peak overlap, respectively; ",
                       "grey length-matched footprints = unknown. PI retains its original RefSeq definition."),
      theme = theme_fiberseq())
  save_figure(p, file.path(OUT_DIR, "plots/example_reads", paste0(gene, ".pdf")),
              width = 7.5, height = 3.8 + n * READ_ROW_IN, limitsize = FALSE)
  invisible(p)
}
