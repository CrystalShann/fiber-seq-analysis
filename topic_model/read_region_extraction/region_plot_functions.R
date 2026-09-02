# Read-level region plots for the macrophage LPS series.
#
# Self-contained: depends only on region_data_utils.R (alongside this file) and CRAN /
# Bioconductor packages. Nothing here sources the fiberhub code.
#
# Three functions, deliberately decoupled:
#
#   load_region_results()  - read one region's extracted files across samples.
#   add_topic_clusters()   - OPTIONAL. Join the topic-model read clusters onto that
#                            result and regroup the pileup. Skip it and everything
#                            downstream still works, ungrouped.
#   plot_region_panels()   - build the stacked panels from whatever came out above.
#
# Input is whatever extract_region_result_macrophage.sh wrote to
#   <out_root>/<outname>/<sample>/parsed/

# Source region_data_utils.R before this file.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

if (!exists("read_ft_mod_region"))
  stop("source region_data_utils.R before region_plot_functions.R")


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

#' Load one region's extracted results across samples.
#'
#' @param region data.frame with columns chr, start, end (one row).
#' @param sample_names character vector of sample labels, e.g. c("LPS_0", "LPS_5", ...).
#' @param result_dirs the matching `parsed` directories, same order as sample_names.
#' @param include_cpg also load CpG calls. The methylation panel draws them as bars and
#'   the read panel as red tiles; FALSE keeps the plot to m6A only.
#' @param min_fiberHMM_fp_score minimum FiberHMM per-block score. The `tf` calls are
#'   already floored at 50; the `footprint` calls carry no score (all zero), so leave
#'   this at 0 if you ever point this at those.
#' @param fiberHMM_size_breaks footprint size bins. Calls outside the outer breaks get
#'   no class and are not drawn.
#' @param nucleosome_size_breaks size window for nucleosome calls, taken from the
#'   FiberHMM `footprint` bed. `ft fire --extract` (fibertools 0.8.2) does not emit
#'   nucleosome segments, so fire.bed cannot supply them.
#' @param nucleosome_label legend label for that class.
#' @return a list holding region, reads, rids_df, fire, fire_peaks, fps, fps_infire,
#'   nucs, pileup, size_levels, nuc_label, sample_names, group_col.
load_region_results <- function(region,
                                sample_names,
                                result_dirs,
                                include_cpg = FALSE,
                                fiberHMM_feature = "tf",
                                min_fiberHMM_fp_score = 50,
                                fiberHMM_size_breaks = c(10, 30, 60, 80),
                                nucleosome_feature = "footprint",
                                nucleosome_size_breaks = c(130, 160),
                                nucleosome_label = "130-160 bp nucleosome",
                                window_n = 10,
                                smooth_pileup = TRUE,
                                verbose = TRUE) {

  stopifnot(is.data.frame(region), nrow(region) == 1)
  stopifnot(length(sample_names) == length(result_dirs))

  region <- list(chr = as.character(region$chr),
                 start = as.numeric(region$start),
                 end = as.numeric(region$end))

  reads_l <- list(); fire_l <- list(); peaks_l <- list(); fps_l <- list(); nuc_l <- list()

  for (i in seq_along(sample_names)) {
    s <- sample_names[i]
    d <- result_dirs[i]
    if (!dir.exists(d)) stop("result dir not found: ", d)
    if (verbose) cat("Loading", s, "from", d, "\n")

    mods <- read_ft_mod_region(file.path(d, "extracted.m6a.bed.gz"),
                               region$start, region$end)
    if (nrow(mods) == 0) stop("no m6A calls in the region for ", s)
    mods$base <- "A"
    if (include_cpg) {
      cg <- read_ft_mod_region(file.path(d, "extracted.cpg.bed.gz"),
                               region$start, region$end)
      if (nrow(cg) > 0) { cg$base <- "CG"; mods <- rbind(mods, cg) }
    }
    mods$sample_name <- s
    reads_l[[i]] <- mods

    fr <- read_fire_region(file.path(d, "fire.bed"))
    fr$sample_name <- s
    fire_l[[i]] <- fr

    peaks_l[[i]] <- read_fire_peaks_region(file.path(d, "fire_peaks.bed"))

    fp <- read_fiberhmm_region(
      file.path(d, paste0("region.fiberhmm_", fiberHMM_feature, ".bed")),
      min_score = min_fiberHMM_fp_score,
      size_breaks = fiberHMM_size_breaks)
    if (nrow(fp) > 0) fp$sample_name <- s
    fps_l[[i]] <- fp

    # Nucleosomes. The `footprint` bed carries no per-block score (all zero), so no
    # score filter here.
    nuc <- read_fiberhmm_region(
      file.path(d, paste0("region.fiberhmm_", nucleosome_feature, ".bed")),
      min_score = 0,
      size_breaks = nucleosome_size_breaks)
    if (nrow(nuc) > 0) {
      nuc <- nuc[!is.na(nuc$class), , drop = FALSE]
      nuc$class <- factor(nucleosome_label, levels = nucleosome_label)
      nuc$sample_name <- s
    }
    nuc_l[[i]] <- nuc
  }

  reads <- do.call(rbind, reads_l)
  fire  <- do.call(rbind, fire_l)
  fps   <- do.call(rbind, Filter(function(x) nrow(x) > 0, fps_l))
  if (is.null(fps)) fps <- data.frame()
  nucs  <- do.call(rbind, Filter(function(x) nrow(x) > 0, nuc_l))
  if (is.null(nucs)) nucs <- data.frame()
  fire_peaks <- do.call(rbind, Filter(Negate(is.null), peaks_l))

  # A read name must not appear in two samples, or per-read joins become ambiguous.
  rid_sample <- unique(fire[, c("RID", "sample_name")])
  if (any(duplicated(rid_sample$RID)))
    stop(sum(duplicated(rid_sample$RID)), " read IDs appear in more than one sample")

  # keep reads that have both a modification call and a FIRE segmentation
  keep <- intersect(unique(reads$RID), unique(fire$RID))
  if (length(keep) == 0) stop("no reads shared between the m6A and FIRE data")
  reads <- reads[reads$RID %in% keep, , drop = FALSE]
  fire  <- fire[fire$RID %in% keep, , drop = FALSE]
  if (nrow(fps) > 0)  fps  <- fps[fps$RID %in% keep, , drop = FALSE]
  if (nrow(nucs) > 0) nucs <- nucs[nucs$RID %in% keep, , drop = FALSE]

  rids_df <- build_rids_df(fire)
  if (verbose)
    cat(nrow(rids_df), "reads in the region across", length(sample_names), "samples\n")

  size_levels <- levels(assign_size_class(numeric(0), fiberHMM_size_breaks))

  # Order reads by footprint similarity so they stay sorted within each facet. Reads
  # with no footprint in the region are not clustered and keep their initial
  # (alphabetical) order at the top. The pairwise distance is an O(n^2) loop, so this
  # is the slow step - roughly a minute at a few hundred reads.
  rid_levels <- order_reads_by_footprints(fps,
                                          sort(unique(as.character(rids_df$RID))),
                                          region$chr, region$start, region$end,
                                          size_levels)

  as_ordered <- function(df) {
    df$RID <- factor(as.character(df$RID), levels = rid_levels)
    df[order(df$RID), , drop = FALSE]
  }
  reads   <- as_ordered(reads)
  fire    <- as_ordered(fire)
  rids_df <- as_ordered(rids_df)
  if (nrow(fps) > 0)  fps  <- as_ordered(fps)
  if (nrow(nucs) > 0) nucs <- as_ordered(nucs)

  reads$base <- factor(reads$base, levels = c("A", "CG"))

  # Only the TF footprints are FIRE-polished. Nucleosomes are protected DNA and FIRE
  # elements are accessible patches, so intersecting the two would discard them all.
  fps_infire <- if (nrow(fps) > 0) subset_footprints_in_fire(fps, fire) else fps
  if (verbose) {
    cat(nrow(fps), "TF footprints,", nrow(fps_infire), "of them inside a FIRE element\n")
    cat(nrow(nucs), "nucleosome calls in",
        paste0(nucleosome_size_breaks, collapse = "-"), "bp\n")
  }

  pileup <- pileup_reads(reads, rids_df, region$start, region$end,
                         split_by = NULL, smooth = smooth_pileup, window_n = window_n)

  list(region = region,
       sample_names = sample_names,
       reads = reads,
       rids_df = rids_df,
       fire = fire,
       fire_peaks = fire_peaks,
       fps = fps,
       fps_infire = fps_infire,
       nucs = nucs,
       pileup = pileup,
       size_levels = size_levels,
       nuc_label = nucleosome_label,
       window_n = window_n,
       smooth_pileup = smooth_pileup,
       group_col = NULL)
}


# ---------------------------------------------------------------------------
# Attach topic-model clusters (optional)
# ---------------------------------------------------------------------------

#' Attach a per-read grouping to a region result and regroup the pileup.
#'
#' @param res result of load_region_results().
#' @param labels named character vector, read ID -> group label.
#' @param group_col name to give the attached column.
#' @param group_levels facet order; defaults to sorted unique labels.
#' @param drop_unassigned drop reads absent from `labels`.
add_read_groups <- function(res, labels, group_col, group_levels = NULL,
                            drop_unassigned = TRUE, verbose = TRUE) {

  if (is.null(names(labels)))
    stop("labels must be a named vector of read ID -> group")
  labels <- labels[!is.na(labels)]
  if (is.null(group_levels))
    group_levels <- sort(unique(unname(labels)))

  keep <- levels(res$rids_df$RID)            # preserves the footprint ordering
  n_before <- length(keep)
  if (drop_unassigned) keep <- keep[keep %in% names(labels)]
  if (length(keep) == 0)
    stop("no reads left after grouping by ", group_col)

  if (verbose)
    cat(n_before, "reads in region result;", length(keep), "grouped;",
        n_before - length(keep), "dropped.\n")

  attach_group <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)
    df <- df[as.character(df$RID) %in% keep, , drop = FALSE]
    df$RID <- factor(as.character(df$RID), levels = keep)
    df[[group_col]] <- factor(unname(labels[as.character(df$RID)]), levels = group_levels)
    df
  }
  for (nm in c("reads", "rids_df", "fire", "fps", "fps_infire", "nucs"))
    res[[nm]] <- attach_group(res[[nm]])

  # the pileup from load_region_results() is ungrouped; recompute it per group or the
  # methylation panel would show one pooled curve while the read panels are faceted
  res$pileup <- pileup_reads(res$reads, res$rids_df,
                             res$region$start, res$region$end,
                             split_by = group_col,
                             smooth = res$smooth_pileup,
                             window_n = res$window_n)
  names(res$pileup)[names(res$pileup) == "group"] <- group_col
  res$pileup[[group_col]] <- factor(res$pileup[[group_col]], levels = group_levels)

  if (verbose) print(table(res$rids_df[[group_col]], res$rids_df$sample_name))

  res$group_col <- group_col
  res
}


#' Join per-read topic-model clusters onto a region result.
#'
#' Reads a `read_topic_assignments_<outname>.tsv` written by the topic-model pipeline
#' and attaches its label to every read-level table.
#'
#' @param assignments_file path to read_topic_assignments_<outname>.tsv.
#' @param assignment_col column holding the label ("cluster" or "dominant_topic").
#' @param drop_unassigned drop reads absent from the assignments file. TRUE keeps the
#'   plotted read set identical to the topic model's.
add_topic_clusters <- function(res,
                               assignments_file,
                               assignment_col = "cluster",
                               group_col = "cluster",
                               drop_unassigned = TRUE,
                               verbose = TRUE) {

  if (!file.exists(assignments_file))
    stop("assignments file not found: ", assignments_file)

  assign_df <- data.table::fread(assignments_file, header = TRUE, data.table = FALSE)
  if (!all(c("RID", assignment_col) %in% colnames(assign_df)))
    stop("assignments file must have columns RID and ", assignment_col)
  if (any(duplicated(assign_df$RID)))
    stop("duplicated RIDs in ", assignments_file)

  if (verbose)
    cat(sum(!assign_df$RID %in% levels(res$rids_df$RID)),
        "assigned reads are absent from the region result.\n")

  add_read_groups(res,
                  labels = setNames(as.character(assign_df[[assignment_col]]), assign_df$RID),
                  group_col = group_col,
                  drop_unassigned = drop_unassigned,
                  verbose = verbose)
}


#' Split a region result by haplotype.
#'
#' The haplotype tag comes from fire.bed, i.e. the HP tag carried by the reads. Useful
#' at heterozygous sites, where the two alleles can differ in accessibility.
#'
#' @param drop_unknown drop reads whose haplotype is UNK (unphased).
add_haplotype_groups <- function(res, group_col = "haplotype",
                                 drop_unknown = TRUE, verbose = TRUE) {

  if (!"tag" %in% colnames(res$rids_df))
    stop("no haplotype tag in rids_df")

  tags <- setNames(as.character(res$rids_df$tag), as.character(res$rids_df$RID))
  if (verbose) {
    cat("haplotype tags:\n"); print(table(tags, useNA = "ifany"))
  }
  if (drop_unknown) tags <- tags[tags %in% c("H1", "H2")]
  if (length(tags) == 0)
    stop("no phased reads in this region")

  add_read_groups(res, labels = tags, group_col = group_col,
                  group_levels = sort(unique(unname(tags))),
                  drop_unassigned = TRUE, verbose = verbose)
}


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

#' Build the stacked region panels.
#'
#' @param res result of load_region_results(), optionally through add_topic_clusters().
#' @param group_col column to facet on; defaults to whatever add_topic_clusters() set,
#'   NULL for no faceting.
#' @param peaks_df optional intervals (chr, start, end, optionally group/score) to draw
#'   as a track under the methylation panel - e.g. motifs or footprints of interest.
#'   One row per interval, labelled on the left by its `group`. NULL (default) omits
#'   the panel entirely.
#' @param highlight optional c(start, end) drawn as a yellow band behind every panel,
#'   e.g. the region a test was run on.
#' @return list(panels = <named list of ggplots>, combined = <patchwork>).
plot_region_panels <- function(res,
                               plot_name = NULL,
                               subtitle = NULL,
                               group_col = res$group_col,
                               peaks_df = NULL,
                               highlight = NULL,
                               fiberHMM_cols = c("purple", "magenta", "#41c4e1"),
                               nucleosome_col = "darkblue",
                               linewidth = 0.5,
                               show_RIDs = FALSE,
                               expand = FALSE,
                               heights = c(pileup = 4, peaks = 3, fire_peaks = 1,
                                           reads = 8, fire_fiberHMM = 8)) {

  region <- res$region
  bp <- c("A" = "blue", "CG" = "red")
  xlims <- c(region$start, region$end)
  if (is.null(subtitle))
    subtitle <- paste0(region$chr, ":", region$start, "-", region$end)

  add_facet <- function(p) {
    if (is.null(group_col)) return(p)
    p + facet_grid(reformulate(".", group_col), scales = "free_y", space = "free_y")
  }
  read_arrow <- arrow(ends = res$rids_df$arrow_end, angle = 60,
                      length = unit(0.02, "inches"))

  gg <- list()

  ## 1. Methylation proportion -------------------------------------------------
  pil <- res$pileup
  gg$pileup <- ggplot(dplyr::filter(pil, base == "A"),
                      aes(x = pos, y = smooth_frac, colour = base)) +
    geom_col(data = dplyr::filter(pil, base == "CG"), aes(y = smooth_frac)) +
    geom_line() +
    theme_bw() +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_colour_manual(values = bp) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Met. prop.", title = plot_name, subtitle = subtitle)
  gg$pileup <- add_facet(gg$pileup) + coord_cartesian(xlim = xlims, expand = expand)

  ## 2. Optional interval track (motifs / footprints of interest) --------------
  if (!is.null(peaks_df) && nrow(peaks_df) > 0) {
    pk <- as.data.frame(peaks_df)
    colnames(pk)[1:3] <- c("chr", "start", "end")
    pk <- dplyr::filter(pk, as.character(chr) == region$chr,
                        end >= region$start, start <= region$end)
    if (nrow(pk) > 0) {
      # one row per interval, labelled on the left by its group (e.g. the cCRE class)
      if (is.null(pk$group)) pk$group <- "Peaks"
      pk$group <- factor(pk$group, levels = unique(pk$group))
      pk <- pk[order(pk$group, pk$start), , drop = FALSE]
      pk$y <- rev(seq_len(nrow(pk)))
      gg$peaks <- ggplot(pk, aes(x = start, xend = end, y = y, yend = y)) +
        geom_segment(linewidth = 1.5, colour = "black") +
        theme_bw() +
        labs(x = NULL, y = NULL) +
        scale_y_continuous(breaks = pk$y, labels = as.character(pk$group),
                           expand = expansion(mult = c(.15, .15))) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 7), panel.grid.minor = element_blank()) +
        coord_cartesian(xlim = xlims, expand = TRUE)
    }
  }

  ## 3. FIRE peaks -------------------------------------------------------------
  ## Peak calls are region-level, not per read, so this panel is never faceted; it
  ## stays a single row while the read panels below are split.
  if (!is.null(res$fire_peaks) && nrow(res$fire_peaks) > 0) {
    fp <- res$fire_peaks
    fp$y <- 1
    gg$fire_peaks <- ggplot(fp, aes(x = start, xend = end, y = y, yend = y,
                                    colour = logFDR)) +
      geom_segment(linewidth = 1) +
      theme_bw() +
      labs(x = NULL, y = NULL, colour = "-log10(FDR)") +
      scale_colour_gradient2(low = "black", mid = "gray", high = "red",
                             midpoint = 1.3, limits = c(0, 5.2),
                             oob = scales::squish) +
      scale_y_continuous(breaks = 1, labels = "FIRE peaks", expand = expansion(add = 1)) +
      theme(legend.key.size = unit(0.4, "cm"),
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 8)) +
      coord_cartesian(xlim = xlims, expand = TRUE)
  }

  ## 4. SMF reads --------------------------------------------------------------
  gg$reads <- ggplot(res$reads) +
    geom_segment(data = res$rids_df, aes(x = start, xend = end, y = RID, yend = RID),
                 colour = "grey50", linewidth = linewidth, arrow = read_arrow) +
    geom_tile(aes(x = pos, y = RID, fill = base)) +
    theme_bw() +
    scale_fill_manual(values = bp) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "SMF reads")
  if (!show_RIDs) gg$reads <- gg$reads + theme(axis.text.y = element_blank())
  gg$reads <- add_facet(gg$reads) + coord_cartesian(xlim = xlims, expand = expand)

  ## 5. FIRE-polished FiberHMM -------------------------------------------------
  ## Layers, bottom to top: read backbone, FIRE/linker segmentation from fire.bed,
  ## FiberHMM nucleosomes, then the FIRE-polished TF footprints.
  has_fps  <- !is.null(res$fps_infire) && nrow(res$fps_infire) > 0
  has_nucs <- !is.null(res$nucs) && nrow(res$nucs) > 0
  if (has_fps || has_nucs) {
    size_levels <- res$size_levels
    if (length(fiberHMM_cols) != length(size_levels))
      stop("fiberHMM_cols has ", length(fiberHMM_cols), " colours but there are ",
           length(size_levels), " size classes: ", paste(size_levels, collapse = ", "))
    names(fiberHMM_cols) <- size_levels
    # `nucleosome` covers fire.bed's own class, which fibertools 0.8.2 never emits but
    # older/newer versions do; it shares the colour of the FiberHMM nucleosome calls.
    class_cols <- c(FIRE = "#FF8C00", linker = "darkgray",
                    nucleosome = nucleosome_col, fiberHMM_cols)
    class_cols[res$nuc_label] <- nucleosome_col

    gg$fire_fiberHMM <- ggplot() +
      geom_segment(data = res$rids_df, aes(x = start, xend = end, y = RID, yend = RID),
                   colour = "grey70", linewidth = linewidth, arrow = read_arrow) +
      geom_segment(data = res$fire,
                   aes(x = start, xend = end, y = RID, yend = RID, colour = class),
                   linewidth = linewidth)
    if (has_nucs)
      gg$fire_fiberHMM <- gg$fire_fiberHMM +
        geom_segment(data = res$nucs,
                     aes(x = start, xend = end, y = RID, yend = RID, colour = class),
                     linewidth = linewidth)
    if (has_fps)
      gg$fire_fiberHMM <- gg$fire_fiberHMM +
        geom_segment(data = dplyr::filter(res$fps_infire, !is.na(class)),
                     aes(x = start, xend = end, y = RID, yend = RID, colour = class),
                     linewidth = linewidth)

    gg$fire_fiberHMM <- gg$fire_fiberHMM +
      theme_bw() +
      scale_colour_manual(values = class_cols) +
      labs(x = region$chr, y = "FIRE polished FiberHMM", colour = "class") +
      theme(legend.position = "right")
    if (!show_RIDs)
      gg$fire_fiberHMM <- gg$fire_fiberHMM + theme(axis.text.y = element_blank())
    gg$fire_fiberHMM <- add_facet(gg$fire_fiberHMM) +
      coord_cartesian(xlim = xlims, expand = expand)
  }

  ## Highlight band ------------------------------------------------------------
  ## Prepended to each panel's layers so it sits behind the data.
  if (!is.null(highlight)) {
    hl <- as.numeric(highlight)
    if (length(hl) != 2) stop("highlight must be c(start, end)")
    hl_layer <- annotate("rect", xmin = hl[1], xmax = hl[2], ymin = -Inf, ymax = Inf,
                         fill = "yellow", alpha = 0.3)
    gg <- lapply(gg, function(p) { p$layers <- c(list(hl_layer), p$layers); p })
  }

  ## Assemble, dropping any panel that is not available -----------------------
  panel_order <- c("pileup", "peaks", "fire_peaks", "reads", "fire_fiberHMM")
  present <- panel_order[panel_order %in% names(gg)]
  combined <- patchwork::wrap_plots(gg[present], ncol = 1,
                                    heights = unname(heights[present]))

  list(panels = gg, combined = combined)
}
