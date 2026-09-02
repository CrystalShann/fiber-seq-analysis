# Read-level co-accessibility figures for cCRE pairs.
#
# Independent implementation of the panels Kevin's get_region_combined_results()
# (fiberhub/code/plots.R) produces for a co-accessible CRE pair. Nothing of his is
# sourced. His version cannot be reused here in any case: rids_df requires a
# sample_name with exactly three underscore-delimited fields (LPS_0 has two),
# cluster = "fire_configs" hard-stops unless exactly two cluster_regions are passed,
# and the panels the published notebooks call (pileup_haps, dimelo_reads) no longer
# exist in the checked-out function.
#
# Everything here reads two tabix-indexed BEDs and nothing else - no per-region
# extraction, no ft, no BAM slicing:
#
#   <s>.read_spans.bed.gz                    chrom start end read_name mapq strand
#   <s>-v0.1-fire-elements.bed.gz            chrom start end read_name . strand
#                                            tstart tend rgb fdr HP     (11 cols)
#
# plus two display-only tracks from the existing `ft extract` run (BED12, blocks):
#   ft_result_dir/<s>/extracted_results/m6a_by_chr/...   methylated adenines
#   ft_result_dir/<s>/extracted_results/nuc_by_chr/...   nucleosomes
#
# Panels, top to bottom in the assembled figure:
#   m6a_prop      fraction of covering fibers methylated per A position, per
#                 timepoint (the topic model's empirical m6A profile)
#   cres          the cCRE track, the tested pair marked
#   m6a           one row per fiber: grey backbone + m6A marks
#   reads         one row per fiber: linker backbone + FiberHMM nucleosomes +
#                 FIRE elements + size-binned FiberHMM footprints, with a class
#                 legend (Kevin's fire_fiberHMM colour scheme)
#
# The tested pair is shaded on EVERY panel. (In Kevin's plots.R the highlight
# annotate('rect', ...) is live on 2 of 10 panels and commented out on all seven
# read-level ones.)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

TABIX_BIN <- "/project/spott/cshan/envs/dimelo/bin/tabix"

# per-fiber m6A and nucleosome calls, from the existing `ft extract` run. These are
# BED12, one row per fiber, blocks = the features. They are NOT the numerator: FIRE
# elements come from lizarraga_FIRE (see load_region). This tree also covers a
# superset of fibers, because it was extracted from the unfiltered BAM while the
# FIRE CRAM is -filtered - so it is used for display only.
FT_RESULT_DIR <- "/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/ft_result_dir"

# FiberHMM per-fiber calls (display only, like the ft tracks). Sample dirs drop
# the underscore (LPS_0 -> LPS0). firehmm_footprint holds the NUCLEOSOME calls
# (BED12 with NO sentinel blocks - every block is a real call);
# firehmm_tf/ft_by_size holds the size-binned footprints as plain BED6
# (chrom start end RID size strand), one row per footprint.
FIREHMM_DIR <- "/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract"
FP_SIZE_BINS <- c("size10-30", "size40-60", "size60-80")

# per-timepoint line colours for the aggregate panels, recycled over `samples`
# (the topic model's tp_colors in 1_m6a_promoter_topic_modelling.Rmd)
TP_COLS       <- c("#0072B2", "#E69F00", "#009E73", "#D55E00")

M6A_COL       <- "blue"
BACKBONE_COL  <- "grey70"
HIGHLIGHT_COL <- "#44d6d8"

# reads-panel classes in draw order (bottom to top; smallest footprints last so
# they stay visible). Colours are Kevin's fire_fiberHMM scheme:
# c(nucleosome = "darkblue", FIRE = "#FF8C00", linker = "darkgray") plus his
# fiberHMM footprint colours for the size bins.
SEG_COLS <- c("linker"       = "darkgray",
              "nucleosome"   = "darkblue",
              "FIRE element" = "#FF8C00",
              "size60-80"    = "purple",
              "size40-60"    = "#8B008B",
              "size10-30"    = "magenta")

# read configuration at the pair, in the order they are stacked (matches the 2x2)
CONFIG_LEVELS <- c("both accessible", "CRE1 only", "CRE2 only", "neither")
CONFIG_COLS   <- c("both accessible" = "#B2182B", "CRE1 only" = "#F4A582",
                   "CRE2 only" = "#92C5DE", "neither" = "grey70")


#' Tabix a BED region into a data.table. Returns an empty table when nothing is found.
tabix_region <- function(path, chrom, start, end, col_names) {
  if (!file.exists(path)) stop("file not found: ", path)
  q <- sprintf("%s:%d-%d", chrom, start + 1L, end)     # tabix is 1-based inclusive
  txt <- suppressWarnings(system2(TABIX_BIN, c(shQuote(path), shQuote(q)),
                                  stdout = TRUE, stderr = FALSE))
  if (length(txt) == 0)
    return(data.table(matrix(character(0), ncol = length(col_names),
                             dimnames = list(NULL, col_names))))
  dt <- data.table::fread(text = paste(txt, collapse = "\n"), header = FALSE,
                          sep = "\t", showProgress = FALSE)
  setnames(dt, seq_along(col_names), col_names)
  dt
}


#' Load fibers and their FIRE elements for one window, across timepoints.
#'
#' @param region list(chrom, start, end) in BED coordinates.
#' @param samples character vector of timepoint names.
#' @param root the co-accessibility output root (holds <s>/<s>.read_spans.bed.gz).
#' @param fire_root the FIRE tree (holds the per-read fire-elements beds).
#' @return list(region, samples, spans, elements).
load_region <- function(region, samples, root, fire_root, fire_ver = "v0.1") {
  span_cols <- c("chrom", "start", "end", "RID", "mapq", "strand")
  elem_cols <- c("chrom", "start", "end", "RID", "score", "strand",
                 "tstart", "tend", "rgb", "fdr", "HP")

  spans <- rbindlist(lapply(samples, function(s) {
    d <- tabix_region(file.path(root, s, paste0(s, ".read_spans.bed.gz")),
                      region$chrom, region$start, region$end, span_cols)
    if (nrow(d)) d[, sample_name := s]
    d
  }), fill = TRUE)

  elements <- rbindlist(lapply(samples, function(s) {
    d <- tabix_region(file.path(fire_root, s, paste0("additional-outputs-", fire_ver),
                                "fire-peaks",
                                paste0(s, "-", fire_ver, "-fire-elements.bed.gz")),
                      region$chrom, region$start, region$end, elem_cols)
    if (nrow(d)) d[, sample_name := s]
    d
  }), fill = TRUE)

  if (nrow(spans) == 0) stop("no reads in ", region$chrom, ":", region$start, "-", region$end)

  # a read name is unique within a timepoint (spans are primary alignments only), but
  # the same name can occur in two timepoints - key on both
  spans[, key := paste(sample_name, RID)]
  if (nrow(elements)) elements[, key := paste(sample_name, RID)]

  list(region = region, samples = samples, spans = spans, elements = elements)
}


#' Expand a BED12 `ft extract` slice into one row per block.
#'
#' In `ft extract` output the first and last block of every row are format
#' sentinels ft adds so the blocks span chromStart..chromEnd (verified on m6a and
#' nuc: first block size 0 at offset 0, last block size 1 ending exactly at
#' chromEnd). Both are dropped when `sentinels = TRUE`, matching
#' convert_ft_bed12_to_bed6() in topic_modelling_functions.r - kept, the trailing
#' one draws a fake 1-bp feature at every fiber's end. The FiberHMM extracts
#' (firehmm_tf / firehmm_footprint) have NO sentinels - every block is a real
#' footprint - so those are read with `sentinels = FALSE`.
#'
#' @param keys optional character vector of `paste(sample, RID)` to keep.
#' @return data.table(key, RID, start, end) in BED coordinates.
read_bed12_blocks <- function(path, chrom, start, end, sample_name, keys = NULL,
                              sentinels = TRUE) {
  # `key` is a reserved argument of data.table(), so a column of that name must be
  # assigned after construction, never passed to the constructor.
  empty <- function() {
    d <- data.table(RID = character(0), start = numeric(0), end = numeric(0))
    d[, key := character(0)][]
  }
  lo <- start; hi <- end        # region bounds, renamed so the column names below
                                # (also start/end) cannot shadow them

  cols <- c("chrom", "chromStart", "chromEnd", "RID", "score", "strand",
            "thickStart", "thickEnd", "rgb", "blockCount", "blockSizes", "blockStarts")
  dt <- tabix_region(path, chrom, lo, hi, cols)
  if (nrow(dt) == 0) return(empty())
  dt[, key := paste(sample_name, RID)]
  if (!is.null(keys)) dt <- dt[key %in% keys]
  if (nrow(dt) == 0) return(empty())

  out <- dt[, {
    sz <- as.numeric(strsplit(blockSizes, ",", fixed = TRUE)[[1]])
    st <- as.numeric(strsplit(blockStarts, ",", fixed = TRUE)[[1]])
    n <- min(length(sz), length(st))
    ok <- if (sentinels) setdiff(seq_len(n), c(1L, n)) else seq_len(n)
    ok <- ok[sz[ok] > 0]                  # drop any residual empty block
    .(start = chromStart + st[ok], end = chromStart + st[ok] + sz[ok])
  }, by = .(key, RID)]

  out[end > lo & start < hi]
}


#' Attach the per-fiber m6A track (`ft extract`), the FiberHMM nucleosome calls
#' (firehmm_footprint), and the size-binned FiberHMM footprints
#' (firehmm_tf/ft_by_size).
#'
#' Display only - the accessibility call used by the statistics stays the FIRE
#' elements loaded in load_region(). Call after load_region(); safe to skip, in
#' which case the corresponding panels/layers are simply omitted.
#'
#' @param keys restrict to these fibers (pass the ones actually being plotted; the
#'   m6A track is ~240 blocks per fiber and expands fast).
load_ft_tracks <- function(res, keys = NULL, ft_root = FT_RESULT_DIR,
                           hmm_root = FIREHMM_DIR) {
  region <- res$region
  grab <- function(kind) {
    rbindlist(lapply(res$samples, function(s) {
      p <- file.path(ft_root, s, "extracted_results", paste0(kind, "_by_chr"),
                     sprintf("%s.ft_extracted_%s.%s.bed.gz", s, kind, region$chrom))
      if (!file.exists(p)) {
        warning("missing ft extract track: ", p, call. = FALSE)
        return(NULL)
      }
      d <- read_bed12_blocks(p, region$chrom, region$start, region$end, s, keys)
      if (nrow(d)) d[, sample_name := s]
      d
    }), fill = TRUE)
  }
  # FiberHMM sample dirs drop the underscore (LPS_0 -> LPS0); keys/sample_name
  # keep the canonical name so they match the spans
  grab_hmm <- function(kind) {
    rbindlist(lapply(res$samples, function(s) {
      s2 <- gsub("_", "", s, fixed = TRUE)
      p <- file.path(hmm_root, paste0("firehmm_", kind), s2,
                     sprintf("%s_hmm_extracted_%s_%s.bed.gz", s2, kind, region$chrom))
      if (!file.exists(p)) {
        warning("missing FiberHMM track: ", p, call. = FALSE)
        return(NULL)
      }
      d <- read_bed12_blocks(p, region$chrom, region$start, region$end, s, keys,
                             sentinels = FALSE)
      if (nrow(d)) d[, sample_name := s]
      d
    }), fill = TRUE)
  }
  # size-binned FiberHMM footprints: plain BED6, one row per footprint
  grab_size <- function(bin) {
    rbindlist(lapply(res$samples, function(s) {
      s2 <- gsub("_", "", s, fixed = TRUE)
      p <- file.path(hmm_root, "firehmm_tf", "ft_by_size", bin, s2,
                     sprintf("%s_tf_%s_%s.bed.gz", s2, bin, region$chrom))
      if (!file.exists(p)) {
        warning("missing FiberHMM size track: ", p, call. = FALSE)
        return(NULL)
      }
      d <- tabix_region(p, region$chrom, region$start, region$end,
                        c("chrom", "start", "end", "RID", "size", "strand"))
      if (nrow(d) == 0) return(NULL)
      d[, key := paste(s, RID)]
      if (!is.null(keys)) d <- d[key %in% keys]
      if (nrow(d) == 0) return(NULL)
      d[, `:=`(sample_name = s, class = bin)]
      d[, .(key, RID, start, end, sample_name, class)]
    }), fill = TRUE)
  }
  res$m6a <- grab("m6a")
  res$nuc <- grab_hmm("footprint")   # the FiberHMM nucleosome calls
  res$size_fps <- rbindlist(lapply(FP_SIZE_BINS, grab_size), fill = TRUE)
  res
}


#' Per-fiber configuration at a cCRE pair.
#'
#' A fiber is accessible at a cCRE when one of its own FIRE elements overlaps that
#' cCRE by >= 1 bp - the same rule the statistics use. `read_rule` must match the
#' rule the table was computed under, or the figure and the 2x2 will disagree.
#'
#' @param cre1,cre2 list(start, end) in BED coordinates.
#' @param read_rule "any" (a fiber counts if it overlaps the cCRE at all, Kevin's
#'   rule) or "contain" (it must span the whole cCRE).
label_reads <- function(res, cre1, cre2, read_rule = c("any", "contain")) {
  read_rule <- match.arg(read_rule)
  sp <- res$spans
  el <- res$elements

  covers <- function(s, e) {
    if (read_rule == "contain") sp$start <= s & sp$end >= e
    else sp$start < e & sp$end > s
  }
  cov1 <- covers(cre1$start, cre1$end)
  cov2 <- covers(cre2$start, cre2$end)

  hit <- function(s, e) {
    if (nrow(el) == 0) return(character(0))
    unique(el$key[el$start < e & el$end > s])
  }
  a1 <- sp$key %in% hit(cre1$start, cre1$end)
  a2 <- sp$key %in% hit(cre2$start, cre2$end)

  shared <- cov1 & cov2
  cfg <- ifelse(!shared, NA_character_,
         ifelse( a1 &  a2, "both accessible",
         ifelse( a1 & !a2, "CRE1 only",
         ifelse(!a1 &  a2, "CRE2 only", "neither"))))

  # `key` is a reserved argument of data.table(), so assign that column afterwards
  out <- data.table(sample_name = sp$sample_name,
                    shared = shared, acc1 = a1, acc2 = a2,
                    config = factor(cfg, levels = CONFIG_LEVELS))
  out[, key := sp$key]
  setcolorder(out, "key")[]
}


#' Order fibers for the raster: by timepoint, then configuration, then position.
#'
#' The configuration ordering is the visual counterpart of the 2x2 table - it is what
#' Kevin's cluster = "fire_configs" achieves, without his config-string machinery.
order_reads <- function(res, labels, samples) {
  d <- merge(res$spans[, .(key, sample_name, start)], labels[, .(key, config)],
             by = "key", all.x = TRUE)
  d[, samp_rank := match(sample_name, samples)]
  d[, cfg_rank := ifelse(is.na(config), length(CONFIG_LEVELS) + 1L, as.integer(config))]
  setorder(d, samp_rank, cfg_rank, start)
  d$key
}


#' Build the stacked panels for one cCRE pair.
#'
#' @param res result of load_region().
#' @param cre1,cre2 list(start, end) - the tested pair, shaded on every panel.
#' @param cres optional data.table(start, end, CRE_ID) of other cCREs to draw.
#' @param drop_unshared drop fibers that do not count toward the 2x2 for this pair.
#' @return list(panels, combined).
plot_pair_panels <- function(res, cre1, cre2, labels, cres = NULL,
                             plot_name = NULL, subtitle = NULL,
                             drop_unshared = TRUE, linewidth = 1.1) {

  region <- res$region
  xlims <- c(region$start, region$end)
  samples <- res$samples
  if (is.null(subtitle))
    subtitle <- sprintf("%s:%d-%d", region$chrom, region$start, region$end)

  keep <- if (drop_unshared) labels[!is.na(config), key] else labels$key
  sp <- res$spans[key %in% keep]
  el <- if (nrow(res$elements)) res$elements[key %in% keep] else res$elements
  if (nrow(sp) == 0) stop("no fibers span both cCREs in this window")

  lev <- order_reads(res, labels, samples)
  lev <- lev[lev %in% keep]
  sp[, key := factor(key, levels = lev)]
  if (nrow(el)) el[, key := factor(key, levels = lev)]
  sp[, sample_name := factor(sample_name, levels = samples)]
  if (nrow(el)) el[, sample_name := factor(sample_name, levels = samples)]

  band <- function(p) {
    p +
      annotate("rect", xmin = cre1$start, xmax = cre1$end, ymin = -Inf, ymax = Inf,
               fill = HIGHLIGHT_COL, alpha = 0.25) +
      annotate("rect", xmin = cre2$start, xmax = cre2$end, ymin = -Inf, ymax = Inf,
               fill = HIGHLIGHT_COL, alpha = 0.25)
  }

  gg <- list()

  ## 2. cCRE track ------------------------------------------------------------
  if (!is.null(cres) && nrow(cres) > 0) {
    ct <- data.table::copy(cres)
    ct[, tested := (start == cre1$start & end == cre1$end) |
                   (start == cre2$start & end == cre2$end)]
    ct[, y := 1]
    gg$cres <- band(
      ggplot(ct, aes(x = start, xend = end, y = y, yend = y, colour = tested)) +
        geom_segment(linewidth = 3) +
        scale_colour_manual(values = c("FALSE" = "grey55", "TRUE" = "#B2182B"),
                            labels = c("FALSE" = "other cCRE", "TRUE" = "tested pair"),
                            name = NULL) +
        scale_y_continuous(breaks = 1, labels = "cCREs") +
        theme_bw(11) +
        theme(axis.text.y = element_text(size = 10), legend.position = "right") +
        labs(x = NULL, y = NULL)) +
      coord_cartesian(xlim = xlims)
  }

  raster_theme <- function(p, ylab) {
    band(p +
      facet_grid(sample_name ~ ., scales = "free_y", space = "free_y") +
      theme_bw(11) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank()) +
      labs(x = NULL, y = ylab)) +
      coord_cartesian(xlim = xlims)
  }

  ## 3. m6A fraction + raster -------------------------------------------------
  ## methylation fraction the way the topic model plots it (the empirical profile
  ## in 1_m6a_promoter_topic_modelling.Rmd): at each A position with an m6A call
  ## in >= 1 fiber, the fraction of covering fibers methylated - 1 = methylated,
  ## 0 = covered but unmethylated, non-covering fibers excluded; per bp, no binning
  if (!is.null(res$m6a) && nrow(res$m6a) > 0) {
    m <- res$m6a[key %in% keep]
    if (nrow(m) > 0) {
      m[, key := factor(key, levels = lev)]
      m[, sample_name := factor(sample_name, levels = samples)]

      prop <- rbindlist(lapply(samples, function(s) {
        ms <- m[sample_name == s]
        ss <- sp[sample_name == s]
        if (nrow(ms) == 0 || nrow(ss) == 0) return(NULL)
        d <- ms[, .(met_n = uniqueN(key)), by = .(pos = start)]
        d[, cov_n := vapply(pos, function(x) sum(ss$start <= x & ss$end > x),
                            numeric(1))]
        d[, .(pos, sample_name = s, frac = met_n / cov_n)]
      }))
      if (nrow(prop) > 0) {
        prop[, sample_name := factor(sample_name, levels = samples)]
        gg$m6a_prop <- band(
          ggplot(prop, aes(pos, frac, colour = sample_name)) +
            geom_line(linewidth = 0.4) +
            facet_grid(sample_name ~ .) +
            scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
            scale_colour_manual(values = rep(TP_COLS, length.out = length(samples))) +
            theme_bw(11) +
            theme(strip.text.y = element_text(angle = 0), legend.position = "none") +
            labs(x = NULL, y = "m6A fraction",
                 title = plot_name, subtitle = subtitle)) +
          coord_cartesian(xlim = xlims)
      }

      ## the raw single-molecule signal: one mark per methylated adenine
      gg$m6a <- raster_theme(
        ggplot() +
          geom_segment(data = sp, aes(x = start, xend = end, y = key, yend = key),
                       colour = BACKBONE_COL, linewidth = linewidth) +
          geom_tile(data = m, aes(x = start, y = key), width = 8, height = 0.85,
                    fill = M6A_COL),
        "m6A")
    }
  }

  ## 4. FiberHMM footprints ---------------------------------------------------
  ## one geom_segment over all classes so they share a legend (Kevin's
  ## fire_fiberHMM panel). Row order = draw order = SEG_COLS order: the linker
  ## backbone at the bottom, then the FiberHMM nucleosome calls, the FIRE
  ## elements (the accessibility call the 2x2 is built from), and the
  ## size-binned FiberHMM footprints on top, smallest last.
  track <- function(d, cls) {
    if (is.null(d) || nrow(d) == 0) return(NULL)
    d <- d[key %in% keep]
    if (nrow(d) == 0) return(NULL)
    d[, .(key = as.character(key), sample_name = as.character(sample_name),
          start, end, class = cls)]
  }
  sz <- res$size_fps
  seg <- rbindlist(c(list(track(sp, "linker"),
                          track(res$nuc, "nucleosome"),
                          track(el, "FIRE element")),
                     lapply(rev(FP_SIZE_BINS), function(b)
                       track(if (!is.null(sz) && nrow(sz)) sz[class == b] else NULL, b))))
  seg[, class := factor(class, levels = names(SEG_COLS))]
  seg[, key := factor(key, levels = lev)]
  seg[, sample_name := factor(sample_name, levels = samples)]
  setorder(seg, class)

  gg$reads <- raster_theme(
    ggplot(seg) +
      geom_segment(aes(x = start, xend = end, y = key, yend = key, colour = class),
                   linewidth = linewidth) +
      scale_colour_manual(values = SEG_COLS, name = NULL, drop = FALSE) +
      guides(colour = guide_legend(override.aes = list(linewidth = 3))),
    "FiberHMM footprints") +
    labs(x = region$chrom)

  order_ <- c("m6a_prop", "cres", "m6a", "reads")
  present <- order_[order_ %in% names(gg)]
  heights <- c(m6a_prop = 4, cres = 1.2, m6a = 11, reads = 13)[present]
  list(panels = gg,
       combined = patchwork::wrap_plots(gg[present], ncol = 1,
                                        heights = unname(heights)))
}


#' Configuration proportions per timepoint - the 2x2 read off the fibers.
plot_config_bars <- function(labels, samples, title = NULL) {
  d <- labels[!is.na(config)]
  d[, sample_name := factor(sample_name, levels = samples)]
  n_lab <- d[, .(n = .N), by = sample_name][, lab := paste0("n=", n)]

  ggplot(d, aes(x = sample_name, fill = config)) +
    geom_bar(position = "fill", width = 0.7) +
    geom_text(data = n_lab, aes(x = sample_name, y = 1.04, label = lab),
              inherit.aes = FALSE, size = 3) +
    scale_fill_manual(values = CONFIG_COLS, drop = FALSE) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.08),
                       breaks = c(0, .25, .5, .75, 1)) +
    theme_bw(11) +
    labs(x = NULL, y = "fibers spanning both cCREs", fill = "configuration",
         title = title)
}
