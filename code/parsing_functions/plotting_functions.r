# Shared plotting palettes
LEIDEN_CLUSTER_COLORS <- c(
  "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1",
  "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
  "green1", "yellow4", "yellow3", "darkorange4", "brown"
)

# timepoint colours (sequential, as in the spearman notebook)
LEIDEN_TIMEPOINT_COLORS <- c("LPS_0" = "#bdbdbd", "LPS_5" = "#6baed6",
                             "LPS_10" = "#2171b5", "LPS_15" = "#08306b")

# ---------------------------------------------------------------------------
# Colour each cluster level by its position: the first level gets the first
# LEIDEN_CLUSTER_COLORS colour, and so on. More than 25 levels are
# interpolated along the palette with colorRampPalette.
#
# Inputs:
#   levels - cluster labels, in display order
# Output:
#   named character vector of colours, names = levels
# ---------------------------------------------------------------------------
cluster_palette <- function(levels) {
  n <- length(levels)
  cols <- if (n <= length(LEIDEN_CLUSTER_COLORS)) {
    LEIDEN_CLUSTER_COLORS[seq_len(n)]
  } else {
    grDevices::colorRampPalette(LEIDEN_CLUSTER_COLORS)(n)
  }
  stats::setNames(cols, levels)
}

# ---------------------------------------------------------------------------
# Colour clusters by the number in their label, so "cluster3" always gets the
# third colour even when other clusters are absent. Falls back to position
# (as cluster_palette()) unless every label is "cluster<N>" with N >= 1.
#
# Inputs:
#   levels - cluster labels, e.g. c("cluster1", "cluster2")
# Output:
#   named character vector of colours, names = levels
# ---------------------------------------------------------------------------
cluster_id_palette <- function(levels) {
  idx <- suppressWarnings(as.integer(sub("^cluster", "", levels)))
  if (anyNA(idx) || any(!grepl("^cluster[0-9]+$", levels))) idx <- seq_along(levels)
  cols <- if (max(idx) <= length(LEIDEN_CLUSTER_COLORS)) LEIDEN_CLUSTER_COLORS[idx] else
    grDevices::colorRampPalette(LEIDEN_CLUSTER_COLORS)(max(idx))[idx]
  stats::setNames(cols, levels)
}


# ---------------------------------------------------------------------------
# Look up LEIDEN_TIMEPOINT_COLORS for timepoint labels written as "LPS_5",
# "5 min" or "5": the minutes are extracted and matched to "LPS_<minutes>".
#
# Inputs:
#   x - timepoint labels (character or factor; repeats allowed)
# Output:
#   named character vector of colours, one per unique label in order of
#   appearance, names = the labels as given; stops for a label with no colour
# ---------------------------------------------------------------------------
timepoint_palette <- function(x) {
  labels <- unique(as.character(x))
  minutes <- trimws(sub("min$", "", sub("^LPS_", "", labels)))
  colors <- unname(LEIDEN_TIMEPOINT_COLORS[paste0("LPS_", minutes)])
  unknown <- is.na(labels) | is.na(colors)
  if (any(unknown)) {
    shown <- ifelse(is.na(labels[unknown]), "<NA>", paste0('"', labels[unknown], '"'))
    stop("No LEIDEN_TIMEPOINT_COLORS entry for: ", paste(shown, collapse = ", "), call. = FALSE)
  }
  stats::setNames(colors, labels)
}


# ---------------------------------------------------------------------------
# Shared ggplot theme: a base theme plus optional element overrides.
#
# Inputs:
#   base        - "classic" (theme_classic), "bw" (theme_bw) or "cowplot"
#                 (cowplot::theme_cowplot)
#   base_size   - base font size
#   base_family - font family
#   overrides   - named list of ggplot2::theme() arguments added on top,
#                 e.g. list(legend.position = "bottom")
# Output:
#   ggplot2 theme object
# ---------------------------------------------------------------------------
theme_fiberseq <- function(base = c("classic", "bw", "cowplot"),
                           base_size = 9, base_family = "", overrides = list()) {
  base <- match.arg(base)
  value <- switch(base,
    classic = ggplot2::theme_classic(base_size = base_size, base_family = base_family),
    bw = ggplot2::theme_bw(base_size = base_size, base_family = base_family),
    cowplot = cowplot::theme_cowplot(font_size = base_size, font_family = base_family))
  if (length(overrides)) value <- value + do.call(ggplot2::theme, overrides)
  value
}


# ---------------------------------------------------------------------------
# Save a figure to `path`, creating its folder if needed.
#   draw = NULL : `plot` is a ggplot / patchwork / cowplot object written with
#                 ggplot2::ggsave
#   draw = fn   : for grid graphics (e.g. ComplexHeatmap) that ggsave cannot
#                 write. A PDF device is opened at `path` (must end in .pdf),
#                 draw(plot) and then after_draw(plot) are called, and the
#                 device is closed even if drawing fails
#
# Inputs:
#   plot          - the object to save (passed to draw() in callback mode)
#   path          - output file
#   width, height - figure size in `units` ("in", "cm", "mm" or "px")
#   device        - ggsave device, or in callback mode a PDF device function
#                   (default grDevices::pdf)
#   dpi, bg, limitsize - passed to ggsave (dpi / bg also used for previews)
#   draw          - optional drawing callback, see above
#   after_draw    - optional callback run after draw(), on the same device
#                   (e.g. to decorate a heatmap)
#   preview       - optional list(path, dpi, bg, embed): also render a PNG to
#                   preview$path (a temporary file when NULL) and pass its path
#                   to preview$embed()
#   on_saved      - optional callback called with `path` after saving
#   ...           - extra arguments for ggsave or the PDF device
# Output:
#   `path`, invisibly
# ---------------------------------------------------------------------------
save_figure <- function(plot, path, width, height,
                        device = NULL, units = "in", dpi = 300,
                        bg = NULL, limitsize = TRUE,
                        draw = NULL, after_draw = NULL,
                        preview = NULL, on_saved = NULL, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # writes one file: the main output, or the PNG preview when preview_image
  render <- function(filename, preview_image = FALSE, image_dpi = dpi, image_bg = bg) {
    if (is.null(draw)) {
      ggplot2::ggsave(filename, plot = plot, width = width, height = height,
        device = if (preview_image) "png" else device, units = units,
        dpi = image_dpi, bg = image_bg, limitsize = limitsize, ...)
    } else {
      stopifnot(is.function(draw))
      # base graphics devices take inches
      inches <- switch(units, "in" = 1, "cm" = 1 / 2.54, "mm" = 1 / 25.4,
                       "px" = 1 / image_dpi, stop("Unsupported units: ", units))
      background <- if (is.null(image_bg)) "white" else image_bg
      if (preview_image) {
        grDevices::png(filename, width = width * inches, height = height * inches,
          units = "in", res = image_dpi, type = "cairo", bg = background)
      } else {
        if (tolower(tools::file_ext(filename)) != "pdf")
          stop("A draw callback requires a PDF output")
        pdf_device <- if (is.null(device)) grDevices::pdf else device
        if (!is.function(pdf_device)) stop("For a draw callback, device must be a function")
        pdf_args <- list(filename, width = width * inches, height = height * inches, ...)
  
        if (!is.null(image_bg)) pdf_args$bg <- image_bg
        do.call(pdf_device, pdf_args)
      }
      opened_device <- grDevices::dev.cur()
      tryCatch({
        draw(plot)
        if (!is.null(after_draw)) after_draw(plot)
      }, finally = grDevices::dev.off(opened_device))
    }
  }
  render(path)
  if (!is.null(preview)) {
    temporary <- is.null(preview$path)
    preview_path <- if (temporary) tempfile(fileext = ".png") else preview$path
    if (temporary) on.exit(unlink(preview_path), add = TRUE)
    dir.create(dirname(preview_path), recursive = TRUE, showWarnings = FALSE)
    render(preview_path, preview_image = TRUE,
      image_dpi = if (is.null(preview$dpi)) dpi else preview$dpi,
      image_bg = if (is.null(preview$bg)) bg else preview$bg)
    if (!is.null(preview$embed)) preview$embed(preview_path)
  }
  if (!is.null(on_saved)) on_saved(path)
  invisible(path)
}


# ---------------------------------------------------------------------------
# Stacked / filled / dodged bar chart of `fill` categories within each `x`
# group. Without `y`, geom_bar counts rows (optionally weighted); with `y`,
# geom_col uses the values in that column.
#
# Inputs:
#   df       - data.frame to plot
#   x, fill  - column names for the bar groups and the stacked categories
#   colors   - fill colours, named by `fill` level
#   y        - optional column of bar heights (e.g. a proportion)
#   weight   - optional count weight column (only used when y is NULL)
#   position - "fill" (scale each bar to 1), "stack" or "dodge"
#   reverse  - reverse the stacking order
#   horizontal - flip to horizontal bars
#   width    - bar width
#   border   - optional list of extra bar arguments, e.g.
#              list(colour = "white", linewidth = 0.2)
#   labels   - list passed to ggplot2::labs(); an optional element
#              `inside = list(min_fraction, accuracy, vjust, color, size)`
#              writes percentage labels (of `y`) inside segments >= min_fraction
#   totals   - optional list(data, x, y, label, style): text drawn at height
#              `y` above each bar from the `label` column of `data`; `style`
#              is a list of extra geom_text arguments
#   scales   - list of ggplot components to add (scales, coords, ...)
#   legend   - extra scale_fill_manual() arguments (name, labels, ...)
#   theme    - optional theme added last
# Output:
#   ggplot object
# ---------------------------------------------------------------------------
plot_stacked_proportion <- function(df, x, fill, colors,
                                    y = NULL, weight = NULL,
                                    position = c("fill", "stack", "dodge"),
                                    reverse = FALSE, horizontal = FALSE,
                                    width = 0.75, border = NULL,
                                    labels = list(), totals = NULL,
                                    scales = list(), legend = list(), theme = NULL) {
  position <- match.arg(position)
  # build aes() from column names; y / weight only when given
  mapping <- ggplot2::aes(x = .data[[x]], fill = .data[[fill]])
  if (!is.null(y)) mapping$y <- ggplot2::aes(y = .data[[y]])$y
  if (!is.null(weight)) mapping$weight <- ggplot2::aes(weight = .data[[weight]])$weight
  placement <- switch(position,
    fill = ggplot2::position_fill(reverse = reverse),
    stack = ggplot2::position_stack(reverse = reverse),
    dodge = ggplot2::position_dodge())
  args <- c(list(width = width, position = placement), border)
  p <- ggplot2::ggplot(df, mapping) +
    do.call(if (is.null(y)) ggplot2::geom_bar else ggplot2::geom_col, args)
  # labels$inside is not a labs() argument; take it out before labs()
  inside <- labels$inside
  labels$inside <- NULL
  if (!is.null(inside)) {
    text <- as.data.frame(df)
    text$.proportion_label <- ifelse(text[[y]] >= inside$min_fraction,
      scales::percent(text[[y]], accuracy = inside$accuracy), "")
    p <- p + ggplot2::geom_text(data = text, ggplot2::aes(label = .proportion_label),
      position = ggplot2::position_stack(vjust = inside$vjust, reverse = reverse),
      color = inside$color, size = inside$size)
  }
  if (!is.null(totals)) {
    text <- as.data.frame(totals$data)
    text$.total_x <- text[[totals$x]]
    text$.total_y <- totals$y
    text$.total_label <- text[[totals$label]]
    p <- p + do.call(ggplot2::geom_text, c(list(data = text,
      mapping = ggplot2::aes(x = .total_x, y = .total_y, label = .total_label),
      inherit.aes = FALSE), totals$style))
  }
  if (horizontal) p <- p + ggplot2::coord_flip()
  p <- p + do.call(ggplot2::scale_fill_manual, c(list(values = colors), legend))
  p <- p + scales + do.call(ggplot2::labs, labels)
  if (!is.null(theme)) p <- p + theme
  p
}


# ---------------------------------------------------------------------------
# Genomic interval track (e.g. cCREs, peaks, tested regions): each interval is
# a rectangle or a horizontal segment on its own row.
#
# Inputs:
#   df             - one row per interval
#   window         - c(start, end) of the x axis
#   fill_col       - optional column coloured by (fill for rect, colour for
#                    segment)
#   highlight      - optional list(start, end, fill, alpha, position): a
#                    background band; position = "under" draws it beneath the
#                    intervals, otherwise on top. Also accepts c(start, end),
#                    a list of rectangles, or a data frame (start/end or xmin/xmax).
#   start_col, end_col - interval coordinate columns
#   row_col        - optional column giving each interval's y row (default 1)
#   geometry       - "rect" or "segment"
#   colors         - optional manual colours for fill_col
#   color_scale    - optional complete ggplot scale (used instead of colors)
#   height         - rectangle height in row units
#   linewidth, lineend - segment style
#   labels         - list passed to ggplot2::labs()
#   scales         - list of ggplot components to add; an element named
#                    `coord` replaces the default coord_cartesian(xlim = window)
#   theme          - optional theme added last
#   clip_intervals - TRUE trims intervals to the window
#   color          - optional fixed outline (rect) / line (segment) colour
# Output:
#   ggplot object
# ---------------------------------------------------------------------------
plot_interval_track <- function(df, window, fill_col = NULL, highlight = NULL,
                                start_col = "start", end_col = "end", row_col = NULL,
                                geometry = c("rect", "segment"),
                                colors = NULL, color_scale = NULL,
                                height = 0.7, linewidth = 1,
                                labels = list(), scales = list(), theme = NULL,
                                clip_intervals = FALSE, color = NULL, lineend = "butt") {
  geometry <- match.arg(geometry)
  # copy the configured columns to fixed internal names
  d <- as.data.frame(df)
  d$.interval_start <- d[[start_col]]
  d$.interval_end <- d[[end_col]]
  d$.interval_row <- if (is.null(row_col)) 1 else d[[row_col]]
  if (clip_intervals) {
    d$.interval_start <- pmax(d$.interval_start, window[1])
    d$.interval_end <- pmin(d$.interval_end, window[2])
  }
  if (geometry == "rect") {
    mapping <- ggplot2::aes(xmin = .interval_start, xmax = .interval_end,
      ymin = .interval_row - height / 2, ymax = .interval_row + height / 2)
    if (!is.null(fill_col)) mapping$fill <- ggplot2::aes(fill = .data[[fill_col]])$fill
    layer <- do.call(ggplot2::geom_rect, c(list(mapping = mapping),
      if (!is.null(color)) list(colour = color)))
  } else {
    mapping <- ggplot2::aes(x = .interval_start, xend = .interval_end,
      y = .interval_row, yend = .interval_row)
    if (!is.null(fill_col)) mapping$colour <- ggplot2::aes(colour = .data[[fill_col]])$colour
    layer <- do.call(ggplot2::geom_segment, c(list(mapping = mapping,
      linewidth = linewidth, lineend = lineend), if (!is.null(color)) list(colour = color)))
  }
  p <- ggplot2::ggplot(d) + layer
  if (!is.null(color_scale)) {
    p <- p + color_scale
  } else if (!is.null(colors)) {
    p <- p + if (geometry == "rect") ggplot2::scale_fill_manual(values = colors) else
      ggplot2::scale_colour_manual(values = colors)
  }
  p <- .fiberseq_add_highlights(p, highlight)
  coordinates <- scales$coord
  scales$coord <- NULL
  if (is.null(coordinates)) coordinates <- ggplot2::coord_cartesian(xlim = window)
  p <- p + scales
  if (!is.null(theme)) p <- p + theme
  p + do.call(ggplot2::labs, labels) + coordinates
}


# ---------------------------------------------------------------------------
# Profile of a value along the genome per group, e.g. the m6A fraction at
# each position for every timepoint or cluster, drawn as lines, areas,
# columns or ribbons.
#
# Inputs:
#   profile     - long table, one row per group x position
#   window      - optional c(start, end) x limits (coord_cartesian)
#   group_col   - optional column mapped to colour (and fill unless
#                 style = "line")
#   style       - default layers: "line", "area" (area + line), "column" or
#                 "ribbon" (ribbon from 0 + line)
#   x_col, y_col - position and value columns
#   facet_col   - optional column to facet by (one column of panels)
#   facet_layout - "wrap" (unchanged default) or "grid"; facet_scales and
#                  facet_space default to "fixed" (space applies to grid only)
#   series_col  - column defining each separate line (default group_col)
#   colors, fill_colors - optional manual colour / fill values
#   smooth_k    - > 1 replaces y by a centred rolling mean of smooth_k points
#                 within each facet/series; rows left NA at the edges are dropped
#   layers      - optional list of layer specs replacing the style defaults,
#                 each list(geom = "line"|"area"|"col"|"ribbon"|"vline", ...)
#                 where ... are arguments for that geom_*()
#   highlight   - optional list of annotate("rect") arguments (xmin, xmax,
#                 fill, alpha, ...) plus position = "under" to draw the band
#                 beneath the layers instead of on top. Also accepts start/end,
#                 a coordinate pair, a list of rectangles, or a data frame.
#   markers     - optional list of geom_vline() arguments
#   region      - region row, required when coordinates = "relative"
#   coordinates - "genomic", or "relative" to convert x to positions relative
#                 to the region anchor (TSS for promoters) with plot_positions()
#                 from leiden_manhattan_plots.r
#   labels      - list passed to ggplot2::labs()
#   scales      - list of ggplot components to add
#   legend      - extra arguments for the manual colour / fill scales
#   theme       - optional theme added last
#   mapping     - optional aes() replacing the default mapping
# Output:
#   ggplot object
# ---------------------------------------------------------------------------
plot_group_profile <- function(profile, window = NULL, group_col = NULL,
                               style = c("line", "area", "column", "ribbon"),
                               x_col = "pos", y_col = "fraction",
                               facet_col = NULL, series_col = group_col,
                               colors = NULL, fill_colors = NULL, smooth_k = 1L,
                               layers = NULL, highlight = NULL, markers = NULL,
                               region = NULL, coordinates = c("genomic", "relative"),
                               labels = list(), scales = list(), legend = list(),
                               theme = NULL, mapping = NULL,
                               facet_layout = c("wrap", "grid"),
                               facet_scales = "fixed", facet_space = "fixed") {
  style <- match.arg(style)
  coordinates <- match.arg(coordinates)
  facet_layout <- match.arg(facet_layout)
  d <- as.data.frame(profile)
  if (coordinates == "relative") {
    stopifnot(!is.null(region))
    d[[x_col]] <- plot_positions(d[[x_col]], region)
  }
  if (smooth_k > 1L) {
    # Keep input feature order: reversing promoter coordinates must not shift
    # the original rolling-mean bins or alter the treatment of edge NAs.
    groups <- unique(c(facet_col, series_col))
    indices <- if (!length(groups)) list(seq_len(nrow(d))) else
      split(seq_len(nrow(d)), interaction(d[groups], drop = TRUE))
    for (i in indices) d[[y_col]][i] <- as.numeric(stats::filter(
      d[[y_col]][i], rep(1 / smooth_k, smooth_k), sides = 2))
    d <- d[!is.na(d[[y_col]]), , drop = FALSE]
  }
  if (is.null(mapping)) {
    mapping <- ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])
    if (!is.null(group_col)) mapping$colour <- ggplot2::aes(colour = .data[[group_col]])$colour
    if (!is.null(group_col) && style != "line")
      mapping$fill <- ggplot2::aes(fill = .data[[group_col]])$fill
    if (!is.null(series_col)) mapping$group <- ggplot2::aes(group = .data[[series_col]])$group
  }
  p <- ggplot2::ggplot(d, mapping)
  if (is.null(layers)) {
    layers <- switch(style,
      line = list(list(geom = "line")),
      column = list(list(geom = "col")),
      area = list(list(geom = "area"), list(geom = "line")),
      ribbon = list(list(geom = "ribbon", mapping = ggplot2::aes(ymin = 0,
        ymax = .data[[y_col]])), list(geom = "line")))
  }
  # each spec becomes geom_<geom>(...) with the remaining entries as arguments
  for (spec in layers) {
    geom <- match.arg(spec$geom, c("line", "area", "col", "ribbon", "vline"))
    spec$geom <- NULL
    p <- p + do.call(getExportedValue("ggplot2", paste0("geom_", geom)), spec)
  }
  p <- .fiberseq_add_highlights(p, highlight)
  if (!is.null(markers)) p <- p + do.call(ggplot2::geom_vline, markers)
  if (!is.null(colors)) p <- p + do.call(ggplot2::scale_colour_manual,
    c(list(values = colors), legend))
  if (!is.null(fill_colors)) p <- p + do.call(ggplot2::scale_fill_manual,
    c(list(values = fill_colors), legend))
  if (!is.null(facet_col)) {
    if (facet_layout == "wrap") {
      p <- p + ggplot2::facet_wrap(stats::reformulate(facet_col), ncol = 1,
                                  scales = facet_scales)
    } else {
      p <- .fiberseq_facet(p, facet_col, facet_scales, facet_space)
    }
  }
  p <- p + scales
  if (!is.null(theme)) p <- p + theme
  if (length(labels)) p <- p + do.call(ggplot2::labs, labels)
  if (!is.null(window)) p <- p + ggplot2::coord_cartesian(xlim = window)
  p
}

# Normalize one rectangle, a list of rectangles, or a data frame of rectangles.
# Both start/end and annotate()'s xmin/xmax are accepted. Legacy helper calls
# still default to an overlay; the region builder supplies yellow/0.3/under.
.fiberseq_highlights <- function(highlight, defaults = list()) {
  if (is.null(highlight)) return(list())
  if (is.numeric(highlight)) {
    if (length(highlight) != 2L) stop("highlight must be c(start, end)")
    highlight <- list(start = highlight[1], end = highlight[2])
  }
  regions <- if (is.data.frame(highlight)) {
    lapply(seq_len(nrow(highlight)), function(i) as.list(highlight[i, , drop = FALSE]))
  } else if (is.list(highlight) && any(c("start", "xmin") %in% names(highlight))) {
    list(highlight)
  } else if (is.list(highlight)) {
    highlight
  } else stop("highlight must be a coordinate pair, rectangle list, or data frame")
  lapply(regions, function(h) {
    if (is.numeric(h) && length(h) == 2L) h <- list(start = h[1], end = h[2])
    if (!is.list(h)) stop("each highlight must describe a rectangle")
    for (nm in names(defaults)) if (is.null(h[[nm]])) h[nm] <- defaults[nm]
    if (!is.null(h$start)) { h$xmin <- h$start; h$start <- NULL }
    if (!is.null(h$end)) { h$xmax <- h$end; h$end <- NULL }
    if (is.null(h$xmin) || is.null(h$xmax)) stop("highlight needs start/end or xmin/xmax")
    if (is.null(h$ymin)) h$ymin <- -Inf
    if (is.null(h$ymax)) h$ymax <- Inf
    h
  })
}

.fiberseq_add_highlights <- function(p, highlight) {
  under <- list()
  over <- list()
  for (h in .fiberseq_highlights(highlight)) {
    behind <- identical(h$position, "under")
    h$position <- NULL
    layer <- do.call(ggplot2::annotate, c(list(geom = "rect"), h))
    if (behind) under <- c(under, list(layer)) else over <- c(over, list(layer))
  }
  p$layers <- c(under, p$layers, over)
  p
}

.fiberseq_facet <- function(p, group_col, scales = "free_y", space = "free_y") {
  if (is.null(group_col)) return(p)
  # Build the formula without requiring an attached stats package.
  p + ggplot2::facet_grid(stats::reformulate(".", group_col), scales = scales, space = space)
}

# Merge named display options while preserving explicit NULL overrides.
.fiberseq_options <- function(defaults, options) {
  if (!is.list(options) || (length(options) &&
      (is.null(names(options)) || any(!nzchar(names(options))) || anyDuplicated(names(options)))))
    stop("display options must be a uniquely named list")
  unknown <- setdiff(names(options), names(defaults))
  if (length(unknown)) stop("unknown display option: ", paste(unknown, collapse = ", "))
  defaults[names(options)] <- options
  defaults
}

#' Build stacked region figures from prepared read-level tables.
#'
#' res contains region (chr/start/end), reads (RID/pos/base), rids_df
#' (RID/start/end, arrow_end when arrows are enabled), fire/nucs/fps_infire
#' (RID/start/end/class), pileup (pos/base/smooth_frac), size_levels in ascending
#' size order, nuc_label and group_col. Read tables also carry the grouping column
#' and share RID factor levels. Optional fire_peaks contains start/end/logFDR.
#' Preparation and fraction calculations belong to the calling data helpers.
#'
#' Existing positional arguments and defaults retain the original region layout.
#' highlight accepts c(start, end), a list of rectangles, or a data frame, with
#' per-rectangle fill/alpha/position overrides. Defaults: yellow, 0.3, under.
#' peaks_df is chr/start/end/group for separate interval rows, or start/end/tested
#' for a single tested-pair cCRE row. It is supplied explicitly by the caller.
#'
#' Display option lists (defaults shown below in the implementation):
#' profile_options: colour_by (base/group), colors, y_col, linewidth,
#'   facet_scales, facet_space, require_data.
#' read_options: tile_width/tile_height (NULL preserves geom_tile sizing),
#'   backbone_col, arrows, base_colors, require_data.
#' peaks_options: variant (intervals/tested); tested is one row, linewidth 3,
#'   tested #B2182B / other grey55 with the original co-accessibility legend.
#' footprint_options: draw_order (NULL retains separate FIRE, nuc, TF layers;
#'   class names combine/sort them into one layer), backbone, always,
#'   legend_title, legend_drop, legend_linewidth.
#' theme_options: base, base_size and per-panel overrides named pileup, peaks,
#'   fire_peaks, reads, fire_fiberHMM. Defaults preserve the region themes.
#' heights and y_labels use canonical panel names. panel_names maps canonical
#' names to returned names, in returned-list order; it does not change stacking
#' order or height lookup.
#' @return list(panels = named ggplots, combined = patchwork).
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
                                           reads = 8, fire_fiberHMM = 8),
                               highlight_fill = "yellow",
                               highlight_alpha = 0.3,
                               highlight_position = c("under", "over"),
                               profile_options = list(),
                               read_options = list(),
                               peaks_options = list(),
                               footprint_options = list(),
                               theme_options = list(),
                               y_labels = c(pileup = "Met. prop.", reads = "SMF reads",
                                            fire_fiberHMM = "FIRE polished FiberHMM"),
                               x_label = res$region$chr,
                               panel_names = c(pileup = "pileup", peaks = "peaks",
                                 fire_peaks = "fire_peaks", reads = "reads",
                                 fire_fiberHMM = "fire_fiberHMM")) {
  region <- res$region
  xlims <- c(region$start, region$end)
  bp <- c(A = "blue", CG = "red")
  if (is.null(subtitle)) subtitle <- paste0(region$chr, ":", region$start, "-", region$end)
  highlight_position <- match.arg(highlight_position)
  highlights <- .fiberseq_highlights(highlight,
    list(fill = highlight_fill, alpha = highlight_alpha, position = highlight_position))
  po <- .fiberseq_options(list(colour_by = "base", colors = bp, y_col = "smooth_frac",
    linewidth = NULL, facet_scales = "free_y", facet_space = "free_y",
    require_data = FALSE), profile_options)
  po$colour_by <- match.arg(po$colour_by, c("base", "group"))
  ro <- .fiberseq_options(list(tile_width = NULL, tile_height = NULL,
    backbone_col = "grey50", arrows = TRUE, base_colors = bp, require_data = FALSE), read_options)
  pkopt <- .fiberseq_options(list(variant = "intervals"), peaks_options)
  pkopt$variant <- match.arg(pkopt$variant, c("intervals", "tested"))
  fo <- .fiberseq_options(list(draw_order = NULL, backbone = TRUE, always = FALSE,
    legend_title = "class", legend_drop = TRUE, legend_linewidth = NULL), footprint_options)
  th <- .fiberseq_options(list(base = "bw", base_size = 11, pileup = list(),
    peaks = list(), fire_peaks = list(), reads = list(), fire_fiberHMM = list()), theme_options)
  panel_theme <- function(panel, overrides = list()) {
    overrides[names(th[[panel]])] <- th[[panel]]
    theme_fiberseq(th$base, base_size = th$base_size, overrides = overrides)
  }
  read_arrow <- if (ro$arrows) grid::arrow(ends = res$rids_df$arrow_end, angle = 60,
    length = grid::unit(0.02, "inches")) else NULL
  has_rows <- function(x) !is.null(x) && nrow(x) > 0L
  gg <- list()

  # Fraction curve: no calculation or smoothing is done here.
  pil <- res$pileup
  if (!po$require_data || has_rows(pil)) {
    colour_col <- if (po$colour_by == "base") "base" else group_col
    if (is.null(colour_col)) stop("group-coloured profiles require group_col")
    line <- list(geom = "line")
    if (!is.null(po$linewidth)) line$linewidth <- po$linewidth
    layers <- list(line)
    profile_data <- pil
    if (po$colour_by == "base") {
      profile_data <- dplyr::filter(pil, base == "A")
      layers <- c(list(list(geom = "col", data = dplyr::filter(pil, base == "CG"),
        mapping = ggplot2::aes(y = .data[[po$y_col]]))), layers)
    }
    gg$pileup <- plot_group_profile(profile_data, group_col = colour_col,
      style = "line", y_col = po$y_col,
      mapping = ggplot2::aes(x = pos, y = .data[[po$y_col]], colour = .data[[colour_col]]),
      layers = layers, highlight = highlights, facet_col = group_col,
      facet_layout = "grid", facet_scales = po$facet_scales, facet_space = po$facet_space,
      colors = po$colors,
      scales = list(ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1))),
      theme = panel_theme("pileup", list(legend.position = "none")),
      labels = list(x = NULL, y = y_labels[["pileup"]], title = plot_name, subtitle = subtitle)) +
      ggplot2::coord_cartesian(xlim = xlims, expand = expand)
  }

  # cCRE/motif intervals: retain both callers' row layouts and scale settings.
  if (has_rows(peaks_df)) {
    pk <- as.data.frame(peaks_df)
    if (pkopt$variant == "tested") {
      pk$y <- 1
      gg$peaks <- plot_interval_track(pk, xlims, fill_col = "tested", row_col = "y",
        geometry = "segment", linewidth = 3, highlight = highlights,
        color_scale = ggplot2::scale_colour_manual(
          values = c("FALSE" = "grey55", "TRUE" = "#B2182B"),
          labels = c("FALSE" = "other cCRE", "TRUE" = "tested pair"), name = NULL),
        scales = list(y = ggplot2::scale_y_continuous(breaks = 1, labels = "cCREs")),
        labels = list(x = NULL, y = NULL),
        theme = panel_theme("peaks", list(axis.text.y = ggplot2::element_text(size = 10),
                                         legend.position = "right")))
    } else {
      colnames(pk)[1:3] <- c("chr", "start", "end")
      pk <- dplyr::filter(pk, as.character(chr) == region$chr,
                          end >= region$start, start <= region$end)
      if (nrow(pk) > 0) {
        if (is.null(pk$group)) pk$group <- "Peaks"
        pk$group <- factor(pk$group, levels = unique(pk$group))
        pk <- pk[order(pk$group, pk$start), , drop = FALSE]
        pk$y <- rev(seq_len(nrow(pk)))
        gg$peaks <- plot_interval_track(pk, xlims, row_col = "y", geometry = "segment",
          linewidth = 1.5, color = "black", highlight = highlights,
          labels = list(x = NULL, y = NULL),
          scales = list(y = ggplot2::scale_y_continuous(breaks = pk$y,
            labels = as.character(pk$group), expand = ggplot2::expansion(mult = c(.15, .15))),
            coord = ggplot2::coord_cartesian(xlim = xlims, expand = TRUE)),
          theme = panel_theme("peaks", list(axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 7),
            panel.grid.minor = ggplot2::element_blank())))
      }
    }
  }
  if (has_rows(res$fire_peaks)) {
    fp <- res$fire_peaks
    fp$y <- 1
    gg$fire_peaks <- plot_interval_track(fp, xlims, fill_col = "logFDR", row_col = "y",
      geometry = "segment", linewidth = 1, highlight = highlights,
      labels = list(x = NULL, y = NULL, colour = "-log10(FDR)"),
      color_scale = ggplot2::scale_colour_gradient2(low = "black", mid = "gray", high = "red",
        midpoint = 1.3, limits = c(0, 5.2), oob = scales::squish),
      scales = list(y = ggplot2::scale_y_continuous(breaks = 1, labels = "FIRE peaks",
        expand = ggplot2::expansion(add = 1)),
        coord = ggplot2::coord_cartesian(xlim = xlims, expand = TRUE)),
      theme = panel_theme("fire_peaks", list(legend.key.size = grid::unit(.4, "cm"),
        legend.text = ggplot2::element_text(size = 7), legend.title = ggplot2::element_text(size = 8))))
  }

  if (!ro$require_data || has_rows(res$reads)) {
    tiles <- list(mapping = ggplot2::aes(x = pos, y = RID, fill = base))
    if (!is.null(ro$tile_width)) tiles$width <- ro$tile_width
    if (!is.null(ro$tile_height)) tiles$height <- ro$tile_height
    p <- ggplot2::ggplot(res$reads) +
      ggplot2::geom_segment(data = res$rids_df,
        ggplot2::aes(x = start, xend = end, y = RID, yend = RID),
        colour = ro$backbone_col, linewidth = linewidth, arrow = read_arrow) +
      do.call(ggplot2::geom_tile, tiles) + ggplot2::scale_fill_manual(values = ro$base_colors)
    p <- .fiberseq_facet(.fiberseq_add_highlights(p, highlights), group_col)
    overrides <- list(legend.position = "none")
    if (!show_RIDs) overrides$axis.text.y <- ggplot2::element_blank()
    gg$reads <- p + panel_theme("reads", overrides) +
      ggplot2::labs(x = NULL, y = y_labels[["reads"]]) +
      ggplot2::coord_cartesian(xlim = xlims, expand = expand)
  }

  has_fps <- has_rows(res$fps_infire)
  has_nucs <- has_rows(res$nucs)
  if (fo$always || has_fps || has_nucs) {
    if (length(fiberHMM_cols) != length(res$size_levels))
      stop("fiberHMM_cols has ", length(fiberHMM_cols), " colours but there are ",
        length(res$size_levels), " size classes: ", paste(res$size_levels, collapse = ", "))
    names(fiberHMM_cols) <- res$size_levels
    class_cols <- c(FIRE = "#FF8C00", linker = "darkgray", nucleosome = nucleosome_col,
                    fiberHMM_cols)
    class_cols[res$nuc_label] <- nucleosome_col
    class_cols["FIRE element"] <- "#FF8C00"
    p <- ggplot2::ggplot()
    if (fo$backbone) p <- p + ggplot2::geom_segment(data = res$rids_df,
      ggplot2::aes(x = start, xend = end, y = RID, yend = RID),
      colour = "grey70", linewidth = linewidth, arrow = read_arrow)
    tracks <- list(res$fire)
    if (has_nucs) tracks <- c(tracks, list(res$nucs))
    if (has_fps) tracks <- c(tracks, list(dplyr::filter(res$fps_infire, !is.na(class))))
    if (!is.null(fo$draw_order)) {
      if (any(!fo$draw_order %in% names(class_cols))) stop("unknown footprint class in draw_order")
      # Keep only columns needed for drawing; input tables may have differing metadata.
      cols <- unique(c("RID", "start", "end", "class", group_col))
      seg <- dplyr::bind_rows(lapply(tracks, function(d) as.data.frame(d)[, cols, drop = FALSE]))
      if (any(!as.character(seg$class) %in% fo$draw_order)) stop("draw_order omits footprint classes")
      seg$class <- factor(seg$class, levels = fo$draw_order)
      seg <- seg[order(seg$class), , drop = FALSE]
      tracks <- list(seg)
      class_cols <- class_cols[fo$draw_order]
    }
    for (d in tracks) p <- p + ggplot2::geom_segment(data = d,
      ggplot2::aes(x = start, xend = end, y = RID, yend = RID, colour = class),
      linewidth = linewidth)
    p <- p + ggplot2::scale_colour_manual(values = class_cols,
      name = fo$legend_title, drop = fo$legend_drop)
    if (!is.null(fo$legend_linewidth)) p <- p + ggplot2::guides(
      colour = ggplot2::guide_legend(override.aes = list(linewidth = fo$legend_linewidth)))
    p <- .fiberseq_facet(.fiberseq_add_highlights(p, highlights), group_col)
    overrides <- list(legend.position = "right")
    if (!show_RIDs) overrides$axis.text.y <- ggplot2::element_blank()
    gg$fire_fiberHMM <- p + panel_theme("fire_fiberHMM", overrides) +
      ggplot2::labs(x = x_label, y = y_labels[["fire_fiberHMM"]], colour = fo$legend_title) +
      ggplot2::coord_cartesian(xlim = xlims, expand = expand)
  }
  panel_order <- c("pileup", "peaks", "fire_peaks", "reads", "fire_fiberHMM")
  present <- panel_order[panel_order %in% names(gg)]
  if (anyNA(heights[present]) || anyNA(panel_names[present]) || anyDuplicated(panel_names[present]))
    stop("heights and panel_names must cover the present canonical panels uniquely")
  combined <- patchwork::wrap_plots(gg[present], ncol = 1, heights = unname(heights[present]))
  gg <- gg[intersect(names(panel_names), names(gg))]
  names(gg) <- unname(panel_names[names(gg)])
  list(panels = gg, combined = combined)
}

# ---------------------------------------------------------------------------
# Shared read-by-position heatmap. Layouts preserve the original matrix encoding,
# row ordering, annotation styles and ComplexHeatmap defaults. Callers supply
# FIRE window labels/ticks so this definition has no notebook-specific globals.
# By default return the Heatmap object; draw=TRUE also applies FIRE boundaries.
#
# Rows are reads split into cluster slices; columns are positions; black = m6A.
#   layout = "genomic"  (default) every bp of region$analysis_start..end from
#                       res$site_met_mat; rows ordered by cluster, read start
#                       (or allele first with split_alleles). Annotations:
#                       cluster, sample, optional haplotype / allele; bottom
#                       coordinate ticks, TSS and SNP marks; minus-strand
#                       promoters are flipped so upstream is on the left
#   layout = "features" the clustering feature matrix res$feat_mat (NA = grey);
#                       rows ordered by cluster, read start; flipped like
#                       "genomic" for minus-strand promoters. Annotations:
#                       cluster, timepoint (sample_name)
#   layout = "fire"     every bp of window$window_start..window_end from
#                       res$site_met_mat, rows given by `rows` in their order.
#                       Annotations: cluster, time; a bottom bar marking the
#                       tested FIRE region and tick labels at `ticks`; with
#                       draw = TRUE, dashed lines at the tested region edges
# Needs plot_anchor(), plot_positions(), as_timepoint_factor(), window_snps(),
# fiberseq_category_palette() and LCL_HAPLOTYPE_COLORS from
# leiden_manhattan_plots.r.
#
# Inputs:
#   res                   - clustering result: assignments (RID, cluster,
#                           start, sample_name, ...), site_met_mat (reads x
#                           positions, rownames = assignments$RID) and, for
#                           "features", feat_mat
#   region                - region row ("genomic" / "features"): chr,
#                           analysis_start, analysis_end, width, and optionally
#                           region_type, strand, tss, annotation, focal_snp,
#                           focal_pos, ref, alt
#   sample_colors         - "genomic": colours named by sample label
#   include_haplotype     - "genomic": add a haplotype annotation, or the
#                           focal allele (assignments$allele_display) for
#                           region_type "top_asfire_het"
#   variants              - "genomic": SNP table (chr, pos, ref, alt,
#                           variant_id) marked below the heatmap when the
#                           region has no focal SNP
#   show_cluster_profiles - "genomic": add per-cluster m6A fraction lines on top
#   split_alleles         - "genomic": split rows by allele, then cluster
#   cluster_label         - "genomic": title text used with include_haplotype
#   sample_label_column   - "genomic": assignments column for sample labels
#                           (default: sample_name up to its first "_")
#   layout                - "genomic", "features" or "fire" (see above)
#   main                  - title ("features", "fire")
#   timepoint_cols        - timepoint colours ("features", "fire")
#   palette               - "features": function giving cluster colours
#   rows                  - "fire": one row per read in display order, with
#                           RID (= rownames of site_met_mat), cluster, time
#   window                - "fire": one-row table with chrom, window_start,
#                           window_end and the tested region start, end
#   cluster_colors        - "fire": cluster colours named by level
#   ticks                 - "fire": genomic positions to label on the x axis
#   use_raster            - optional override of ComplexHeatmap rasterising
#   heatmap_options       - named list overriding any ComplexHeatmap::Heatmap()
#                           argument
#   draw                  - FALSE returns the Heatmap; TRUE draws it on the
#                           current device
#   draw_options          - extra ComplexHeatmap::draw() arguments
# Output:
#   draw = FALSE: ComplexHeatmap Heatmap object; draw = TRUE: invisible NULL
#   (the heatmap is drawn as a side effect)
# ---------------------------------------------------------------------------
plot_read_heatmap <- function(res, region = res$region, sample_colors = NULL,
                              include_haplotype = FALSE, variants = NULL,
                              show_cluster_profiles = FALSE, split_alleles = FALSE,
                              cluster_label = "Saved m6A-defined clusters",
                              sample_label_column = NULL,
                              layout = c("genomic", "features", "fire"),
                              main = NULL, timepoint_cols = LEIDEN_TIMEPOINT_COLORS,
                              palette = cluster_palette, rows = NULL, window = NULL,
                              cluster_colors = NULL, ticks = NULL,
                              use_raster = NULL, heatmap_options = list(),
                              draw = FALSE, draw_options = list()) {
  layout <- match.arg(layout)
  boundary_positions <- NULL
  if (layout == "features") {
    # takes per read assignment from the clustering result
    df <- res$assignments
    df$sample_name <- as_timepoint_factor(df$sample_name, timepoint_cols)
    # sort by cluster
    # within each cluster, sort by read start position
    o  <- order(df$cluster, df$start)
    df <- df[o, ]
    # Takes the feature matrix used for  clustering
    # and puts its rows in exactly the same order as df
    mat <- res$feat_mat[df$RID, , drop = FALSE]
    if (plot_anchor(region)$direction == -1L)
      mat <- mat[, rev(seq_len(ncol(mat))), drop = FALSE]
    mat <- matrix(ifelse(is.na(mat), NA, ifelse(mat > 0, "m6A", "no m6A")),
                  nrow(mat), ncol(mat), dimnames = dimnames(mat))
  
    # create row annotations
    ha <- ComplexHeatmap::rowAnnotation(
      cluster   = df$cluster,
      timepoint = df$sample_name,
      col = list(cluster   = palette(levels(df$cluster)),
                 timepoint = timepoint_cols[levels(df$sample_name)]))
  
    # rows are reads, columns are m6a sites
    heatmap_args <- list(
      matrix = mat,
      name = "m6A call",
      col  = c("m6A" = "black", "no m6A" = "white"),
      na_col = "grey85",
      show_row_names = FALSE, show_column_names = FALSE,
      cluster_rows = FALSE, cluster_columns = FALSE,
      row_split = df$cluster, row_gap = grid::unit(0.6, "mm"),
      row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 8),
      width = grid::unit(11, "cm"), height = grid::unit(14, "cm"),
      use_raster = TRUE,
      column_title = main,
      column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
      left_annotation = ha)
  } else if (layout == "fire") {
    lr <- res
    ex <- window
    cl_cols <- cluster_colors
    time_cols <- timepoint_cols
    # one column per bp of the window; positions without a site stay 0
    pos  <- seq.int(ex$window_start, ex$window_end)
    site <- as.matrix(lr$site_met_mat[rows$RID, , drop = FALSE])
    m <- matrix(0L, nrow(rows), length(pos), dimnames = list(rows$RID, pos))
    m[, match(colnames(site), pos)] <- site
    n_cl  <- table(rows$cluster)
    in_region <- ifelse(pos > ex$start & pos <= ex$end, "tested FIRE region", "outside")
  
    heatmap_args <- list(
      matrix = m, name = "m6A", col = c("0" = "white", "1" = "black"),
      heatmap_legend_param = list(at = c(0, 1), labels = c("no m6A call", "m6A call"),
                                  title = "m6A"),
      cluster_rows = FALSE, cluster_columns = FALSE,
      row_split = rows$cluster,
      row_title = sprintf("%s\nn = %d", names(n_cl), as.integer(n_cl)),
      row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 8),
      row_gap = grid::unit(1, "mm"), border = TRUE,
      show_row_names = FALSE, show_column_names = FALSE,
      use_raster = TRUE, raster_quality = 4,
      left_annotation = ComplexHeatmap::rowAnnotation(
        cluster = rows$cluster, time = rows$time,
        col = list(cluster = cl_cols, time = time_cols),
        annotation_name_gp = grid::gpar(fontsize = 8),
        annotation_legend_param = list(cluster = list(title = "Leiden cluster"),
                                       time = list(title = "Time"))),
      bottom_annotation = ComplexHeatmap::HeatmapAnnotation(
        region = ComplexHeatmap::anno_simple(in_region, col = c("tested FIRE region" = "black", outside = "white"),
                             height = grid::unit(1.5, "mm")),
        position = ComplexHeatmap::anno_mark(at = match(ticks, pos), labels = scales::comma(ticks),
                             side = "bottom", labels_rot = 0,
                             labels_gp = grid::gpar(fontsize = 8)),
        show_annotation_name = c(region = TRUE, position = FALSE),
        annotation_name_gp = grid::gpar(fontsize = 8), annotation_label = c("tested region", "")),
      column_title = paste0(main, "\nLeiden + Manhattan clusters, reads pooled 0-15 min (",
                            ex$chrom, ", ", nrow(rows), " molecules)"),
      column_title_gp = grid::gpar(fontsize = 9))
  
    # tested region edges as fractions of the heatmap width, drawn after draw()
    boundary_positions <- match(c(ex$start, ex$end), pos) / length(pos)
    boundary_slices <- seq_along(n_cl)
  } else {
    assignments <- res$assignments
    if (split_alleles) stopifnot(include_haplotype, !is.null(assignments$allele_display))
    ordering <- if (split_alleles) order(assignments$allele_display, assignments$cluster, assignments$start, assignments$RID) else
      order(assignments$cluster, assignments$start, assignments$RID)
    assignments <- assignments[ordering, , drop = FALSE]
    sample_labels <- if (is.null(sample_label_column)) sub("_.*$", "", assignments$sample_name) else
      assignments[[sample_label_column]]
    if (!is.null(sample_label_column))
      stopifnot(length(sample_labels) == nrow(assignments), all(sample_labels %in% names(sample_colors)))
    sample <- factor(sample_labels, levels = names(sample_colors))
    annotation <- ComplexHeatmap::rowAnnotation(
      cluster = assignments$cluster, sample = sample,
      col = list(cluster = cluster_id_palette(levels(assignments$cluster)), sample = sample_colors))
    if (include_haplotype && identical(region$region_type, "top_asfire_het")) {
      annotation <- ComplexHeatmap::rowAnnotation(
        cluster = assignments$cluster, sample = sample, allele = assignments$allele_display,
        col = list(cluster = cluster_id_palette(levels(assignments$cluster)), sample = sample_colors,
          allele = fiberseq_category_palette(res$allele_display_levels)))
    } else if (include_haplotype) {
      annotation <- ComplexHeatmap::rowAnnotation(
        cluster = assignments$cluster, sample = sample, haplotype = assignments$haplotype,
        col = list(cluster = cluster_id_palette(levels(assignments$cluster)), sample = sample_colors,
          haplotype = LCL_HAPLOTYPE_COLORS))
    }
    # one column per bp of the analysis window; positions without a site stay 0
    met <- matrix(0L, nrow(assignments), region$width,
                    dimnames = list(assignments$RID, seq.int(region$analysis_start, region$analysis_end)))
    met[, match(colnames(res$site_met_mat), colnames(met))] <- as.matrix(res$site_met_mat[assignments$RID, , drop = FALSE])
    anchor <- plot_anchor(region)
    if (anchor$direction == -1L) met <- met[, rev(seq_len(ncol(met))), drop = FALSE]
    display <- met
    colors <- c("0" = "white", "1" = "black")
    legend <- list(at = c(0, 1), labels = c("no m6A call", "m6A"))
    legend_name <- "m6A"
    top <- NULL
    if (show_cluster_profiles) {
      cluster_colors <- cluster_id_palette(levels(assignments$cluster))
      profile_annotations <- lapply(levels(assignments$cluster), function(cluster) {
        selected <- assignments$cluster == cluster
        ComplexHeatmap::anno_lines(
          colMeans(met[selected, , drop = FALSE]),
          ylim = c(0, 1), gp = grid::gpar(col = cluster_colors[[cluster]], lwd = 0.7),
          axis_param = list(at = c(0, 0.5, 1), labels = c("0", ".5", "1")),
          height = grid::unit(12, "mm"))
      })
      names(profile_annotations) <- paste0(levels(assignments$cluster), " m6A")
      top <- do.call(ComplexHeatmap::HeatmapAnnotation, c(profile_annotations,
        list(annotation_name_gp = grid::gpar(fontsize = 8), gap = grid::unit(1.5, "mm"))))
    }
    # bottom annotations: 5 coordinate ticks, then TSS and SNP marks if present
    ticks <- unique(round(seq(1, ncol(display), length.out = 5L)))
    coordinates <- as.integer(colnames(display))
    tick_positions <- if (anchor$promoter) plot_positions(coordinates[ticks], region) else coordinates[ticks]
    bottom_parts <- list(coordinate = ComplexHeatmap::anno_mark(at = ticks,
      labels = format(tick_positions, scientific = FALSE, trim = TRUE),
      which = "column", side = "bottom", labels_gp = grid::gpar(fontsize = 8)))
    if (identical(region$region_type, "promoter") && !is.null(region$tss) &&
        !is.na(region$tss) && region$tss >= region$analysis_start && region$tss <= region$analysis_end) {
      bottom_parts$TSS <- ComplexHeatmap::anno_mark(
        at = match(region$tss, coordinates),
        labels = if (anchor$promoter) "TSS: 0" else paste0("TSS: ", region$chr, ":", region$tss),
        which = "column", side = "bottom", labels_gp = grid::gpar(fontsize = 8, col = "#D55E00"))
    }
    snps <- if (!is.null(region$focal_snp)) data.frame(pos = region$focal_pos,
      label = paste0(region$focal_snp, " ", region$ref, ">", region$alt)) else window_snps(variants, region)
    if (nrow(snps)) {
      bottom_parts$SNP <- ComplexHeatmap::anno_mark(
        at = match(snps$pos, coordinates),
        labels = if (anchor$promoter) paste0(snps$label, " (", plot_positions(snps$pos, region), " bp from TSS)") else snps$label,
        which = "column", side = "bottom", labels_gp = grid::gpar(fontsize = 7),
        link_gp = grid::gpar(col = "#984EA3"))
    }
    bottom <- do.call(ComplexHeatmap::HeatmapAnnotation, c(bottom_parts,
      list(annotation_name_gp = grid::gpar(fontsize = 8))))
    heatmap_args <- list(
      matrix = display, name = legend_name, col = colors, heatmap_legend_param = legend,
      cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE,
      row_split = if (split_alleles) data.frame(allele = assignments$allele_display,
        cluster = assignments$cluster) else assignments$cluster,
      row_gap = grid::unit(if (split_alleles) 3 else 0.6, "mm"),
      row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 8),
      show_row_names = FALSE, show_column_names = FALSE,
      use_raster = TRUE, raster_quality = 2, raster_resize_mat = FALSE,
      left_annotation = annotation,
      top_annotation = top, bottom_annotation = bottom,
      column_title = paste(c(region$annotation,
        if (include_haplotype) paste0(cluster_label, "; focal allele annotated per sample"),
        if (split_alleles) "Rows split by focal allele, then cluster",
        if (show_cluster_profiles) "Top: m6A call fraction per cluster (all full-span reads)",
        if (anchor$promoter) anchor$x_label,
        if (include_haplotype && !nrow(snps)) "No phased heterozygous SNP in this window",
        paste0(region$chr, ":", region$analysis_start, "-", region$analysis_end)), collapse = "\n"),
      column_title_gp = grid::gpar(fontsize = 11))
  }
  # caller overrides, then build (and optionally draw) the heatmap
  if (!is.null(use_raster)) heatmap_args$use_raster <- use_raster
  for (name in names(heatmap_options)) heatmap_args[name] <- heatmap_options[name]
  heatmap <- do.call(ComplexHeatmap::Heatmap, heatmap_args)
  if (!draw) return(heatmap)
  if (layout == "fire") draw_options <- utils::modifyList(list(
    merge_legend = TRUE, heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"), draw_options)
  do.call(ComplexHeatmap::draw, c(list(object = heatmap), draw_options))
  if (!is.null(boundary_positions)) {
    for (slice in boundary_slices) ComplexHeatmap::decorate_heatmap_body("m6A", slice = slice, {
      grid::grid.segments(x0 = boundary_positions, x1 = boundary_positions, y0 = 0, y1 = 1,
        gp = grid::gpar(lty = 2, lwd = 0.8, col = "#D55E00"))
    })
  }
  invisible(NULL)
}


# ===========================================================================
# cCRE-pair co-accessibility figures and one-call read-level region plots
#
# coaccess_region_result() and configuration_proportion_inputs() prepare a
# cCRE pair for plot_region_panels(). plot_coaccess_pair() and
# plot_region_example() load data with parsing_footprints_functions.r (source
# it first) and draw with plot_region_panels().
# ===========================================================================

# m6A tile and fiber backbone colours of the co-accessibility read panels
M6A_COL       <- "blue"
BACKBONE_COL  <- "grey70"

# configuration colours, named by CONFIG_LEVELS (parsing_footprints_functions.r)
CONFIG_COLS   <- c("both accessible" = "#B2182B", "CRE1 only" = "#F4A582",
                   "CRE2 only" = "#92C5DE", "neither" = "grey70")


# ---------------------------------------------------------------------------
# Adapt a cCRE pair's fibers to the plot_region_panels() input structure.
# Fibers are ordered by order_reads(); the canonical RID is the fiber key
# ("<sample> <read>", raw_RID keeps the read name) because the same read name
# may occur in several samples. No FIRE polishing or fraction smoothing.
#
# Inputs:
#   res           - load_region() + load_ft_tracks() result
#   cre1, cre2    - list(start, end) of the tested pair
#   labels        - label_reads() result (fibers with NA config are not drawn)
#   cres          - optional cCREs in the window (start, end, ...); a tested
#                   column is added for the cCRE track
#   drop_unshared - drop fibers that do not cover both cCREs
# Output:
#   list in the plot_region_panels() format (region, sample_names, reads,
#   rids_df, fire, fps, fps_infire, nucs, pileup, size_levels = FP_SIZE_BINS,
#   nuc_label, group_col = "sample_name", peaks_df, highlight = the pair)
# ---------------------------------------------------------------------------
coaccess_region_result <- function(res, cre1, cre2, labels, cres = NULL,
                                   drop_unshared = TRUE) {
  samples <- res$samples
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

  m <- if (!is.null(res$m6a) && nrow(res$m6a)) res$m6a[key %in% keep] else NULL
  if (!is.null(m) && nrow(m)) {
    m[, key := factor(key, levels = lev)]
    m[, sample_name := factor(sample_name, levels = samples)]
  }
  prop <- as.data.frame(coaccess_m6a_fraction(m, sp, samples))
  prop$base <- rep("A", nrow(prop))

  # Return consistently typed empty tables when an optional track is absent.
  as_region_rows <- function(d, cls = NULL) {
    if (is.null(d) || nrow(d) == 0) {
      out <- data.frame(RID = factor(character(), levels = lev), raw_RID = character(),
        start = numeric(), end = numeric(), sample_name = factor(character(), levels = samples))
      out$class <- character()
      return(out)
    }
    d <- as.data.frame(d)
    d <- d[as.character(d$key) %in% keep, , drop = FALSE]
    out <- data.frame(RID = factor(as.character(d$key), levels = lev),
      raw_RID = as.character(d$RID), start = d$start, end = d$end,
      sample_name = factor(as.character(d$sample_name), levels = samples))
    if (!is.null(cls)) out$class <- rep(cls, nrow(out))
    else if ("class" %in% names(d)) out$class <- as.character(d$class)
    out
  }
  reads <- as_region_rows(m)
  reads$pos <- reads$start
  reads$base <- factor(rep("A", nrow(reads)), levels = c("A", "CG"))
  fire <- rbind(as_region_rows(sp, "linker"), as_region_rows(el, "FIRE element"))
  nucs <- as_region_rows(res$nuc, "nucleosome")
  fps <- as_region_rows(res$size_fps)
  # The old drawing loop included exactly these bins, in descending size order.
  fps <- fps[fps$class %in% FP_SIZE_BINS, , drop = FALSE]
  fps$class <- factor(fps$class, levels = FP_SIZE_BINS)
  ct <- NULL
  if (!is.null(cres) && nrow(cres) > 0) {
    ct <- as.data.frame(cres)
    ct$tested <- (ct$start == cre1$start & ct$end == cre1$end) |
                 (ct$start == cre2$start & ct$end == cre2$end)
  }
  list(region = list(chr = res$region$chrom, start = res$region$start, end = res$region$end),
    sample_names = samples, reads = reads, rids_df = as_region_rows(sp),
    fire = fire, fire_peaks = NULL, fps = fps, fps_infire = fps, nucs = nucs,
    pileup = prop, size_levels = FP_SIZE_BINS, nuc_label = "nucleosome",
    group_col = "sample_name", peaks_df = ct,
    highlight = list(c(cre1$start, cre1$end), c(cre2$start, cre2$end)))
}


# ---------------------------------------------------------------------------
# Configuration proportions per sample (or facet) - the 2x2 read off the fibers.
#
# Inputs:
#   labels  - label_reads() result (sample_name holds the bar group)
#   samples - bar order
# Output:
#   list(rows = fibers covering both cCREs, totals = n per bar with label "n=<n>")
# ---------------------------------------------------------------------------
configuration_proportion_inputs <- function(labels, samples) {
  d <- data.table::copy(labels[!is.na(config)])
  d[, sample_name := factor(sample_name, levels = samples)]
  n_lab <- d[, .(n = .N), by = sample_name][, lab := paste0("n=", n)]
  list(rows = d, totals = n_lab)
}


# ---------------------------------------------------------------------------
# Read-level figure of one cCRE pair: m6A fraction per facet, cCRE track, m6A
# raster and footprint raster (linker, nucleosomes, FIRE elements, size-binned
# TF footprints) with the pair shaded on every panel, plus configuration bars.
# Macrophage: read spans, FiberHMM nucleosomes and size bins, one facet per
# timepoint. LCL: primary read spans, ft nucleosomes, FiberHMM v2 TF footprints,
# pooled / per-sample facets. Needs parsing_footprints_functions.r.
#
# Inputs:
#   pair        - one row with chrom, cre1_start, cre1_end, cre2_start, cre2_end
#                 (BED); CRE1/CRE2 keep the pair table's orientation
#   samples     - sample names
#   paths       - list(root, fire_root, ft_root, hmm_root, tf_root, tabix_bin,
#                 fire_ver, hmm_label), see load_region() and load_ft_tracks();
#                 entries a source does not use may be omitted
#   span_source - "read_spans" or "fire_all", see load_region()
#   nuc_source  - "fiberhmm" or "ft", see load_ft_tracks()
#   tf_source   - "by_size", "recalled_tf" or "none", see load_ft_tracks()
#   facet_by    - "sample" (one facet per sample), "pooled" (one facet; fiber
#                 keys keep their sample) or "hp" (H1 / H2 / UNK from fire_all)
#   colors      - facet colours of the profile lines; default cluster_palette()
#   cres        - optional cCREs (chrom, start, end, CRE_ID, CRE_label)
#   cre_track   - "tested": one row, tested pair red; "by_class": one row per class
#   cre_track_types - cCRE classes drawn (NULL = all); the tested pair is always kept
#   flank       - bp either side of the pair in the window
#   read_rule   - see label_reads()
#   plot_name   - title
#   subtitle    - text appended to the window coordinates ("chr:start-end  |  ...")
#   show_2x2    - add the 2x2 and Fisher OR of coaccess_2x2() to the subtitle
#   heights     - panel heights c(pileup, peaks, reads, fire_fiberHMM); NULL
#                 scales them with the facets and drawn fibers
#   max_fibers_per_facet - optional cap on the fibers drawn per facet, chosen at
#                 random after set.seed(seed); the 2x2 and bars use every fiber
#   seed        - seed of that downsampling
#   bar_by      - configuration bars per "facet" or per "sample"
#   pooled_bar  - add a "pooled" bar when there is more than one bar
#   bar_title   - bar chart title
#   y_labels    - panel y labels for plot_region_panels()
#   linewidth   - segment width of the read panels
#   fire_overlap_fraction - see label_reads(); must match the analysis table
# Output:
#   list(res, labels, counts = coaccess_2x2(), panels = plot_region_panels()
#   result, bars, region, height = suggested figure height in inches)
# ---------------------------------------------------------------------------
plot_coaccess_pair <- function(pair, samples, paths,
                               span_source = c("read_spans", "fire_all"),
                               nuc_source = c("fiberhmm", "ft"),
                               tf_source = c("by_size", "recalled_tf", "none"),
                               facet_by = c("sample", "pooled", "hp"),
                               colors = NULL, cres = NULL,
                               cre_track = c("tested", "by_class"), cre_track_types = NULL,
                               flank = 1000, read_rule = "any",
                               plot_name = NULL, subtitle = NULL, show_2x2 = TRUE,
                               heights = NULL, max_fibers_per_facet = NULL, seed = 1,
                               bar_by = c("facet", "sample"), pooled_bar = FALSE,
                               bar_title = NULL, y_labels = NULL, linewidth = 1.1,
                               fire_overlap_fraction = 0) {
  if (!exists("load_region", mode = "function"))
    stop("source parsing_footprints_functions.r before plot_coaccess_pair()")
  span_source <- match.arg(span_source)
  nuc_source <- match.arg(nuc_source)
  tf_source <- match.arg(tf_source)
  facet_by <- match.arg(facet_by)
  cre_track <- match.arg(cre_track)
  bar_by <- match.arg(bar_by)
  if (facet_by == "hp" && span_source != "fire_all")
    stop("facet_by = 'hp' needs span_source = 'fire_all' (read spans carry no haplotype)")

  region <- list(chrom = pair$chrom,
                 start = max(0, min(pair$cre1_start, pair$cre2_start) - flank),
                 end   = max(pair$cre1_end, pair$cre2_end) + flank)
  cre1 <- list(start = pair$cre1_start, end = pair$cre1_end)
  cre2 <- list(start = pair$cre2_start, end = pair$cre2_end)

  res <- load_region(region, samples, paths$root, paths$fire_root,
                     fire_ver = if (is.null(paths$fire_ver)) "v0.1" else paths$fire_ver,
                     tabix_bin = paths$tabix_bin, span_source = span_source)
  labs <- label_reads(res, cre1, cre2, read_rule = read_rule,
                      fire_overlap_fraction = fire_overlap_fraction)
  # display tracks only for the fibers that can be drawn
  track_args <- list(res, keys = labs[!is.na(config), key], ft_root = paths$ft_root,
                     hmm_root = paths$hmm_root, tabix_bin = paths$tabix_bin,
                     nuc_source = nuc_source, tf_source = tf_source, tf_root = paths$tf_root)
  if (!is.null(paths$hmm_label)) track_args$hmm_label <- paths$hmm_label
  res <- do.call(load_ft_tracks, track_args)
  counts <- coaccess_2x2(labs)

  # facet of every fiber: its sample, "pooled", or its haplotype
  facet_of <- switch(facet_by,
    sample = NULL,
    pooled = stats::setNames(rep("pooled", nrow(res$spans)), res$spans$key),
    hp     = stats::setNames(ifelse(is.na(res$spans$HP) | res$spans$HP == "", "UNK",
                                    res$spans$HP), res$spans$key))
  labs_facet <- labs
  if (!is.null(facet_of)) {
    present <- unique(unname(facet_of[labs$key[!is.na(labs$config)]]))
    hp_order <- c("H1", "H2", "UNK")
    facets <- if (facet_by == "pooled") "pooled" else
      c(intersect(hp_order, present), sort(setdiff(present, hp_order)))
    # replace the sample column of a fiber table by the fiber's facet
    relabel <- function(d) {
      if (!is.null(d) && nrow(d)) d[, sample_name := unname(facet_of[key])]
      d
    }
    for (nm in c("spans", "elements", "m6a", "nuc", "size_fps")) res[[nm]] <- relabel(res[[nm]])
    labs_facet <- data.table::copy(labs)[, sample_name := unname(facet_of[key])]
    res$samples <- facets
  }

  # optional per-facet downsampling of the drawn fibers only
  labs_plot <- labs_facet
  n_shared <- sum(!is.na(labs_facet$config))
  if (!is.null(max_fibers_per_facet)) {
    set.seed(seed)
    shared <- labs_facet[!is.na(config)]
    drop <- unlist(lapply(split(shared$key, shared$sample_name), function(k)
      if (length(k) > max_fibers_per_facet) setdiff(k, sample(k, max_fibers_per_facet)) else character(0)),
      use.names = FALSE)
    labs_plot <- data.table::copy(labs_facet)
    labs_plot[key %in% drop, config := NA]
  }
  n_drawn <- sum(!is.na(labs_plot$config))

  # cCREs in the window; cre_track_types filters them but keeps the tested pair
  cres_here <- NULL
  if (!is.null(cres)) {
    cres_here <- data.table::as.data.table(cres)[chrom == region$chrom & end > region$start &
                                                   start < region$end]
    if (!is.null(cre_track_types)) {
      tested <- (cres_here$start == cre1$start & cres_here$end == cre1$end) |
                (cres_here$start == cre2$start & cres_here$end == cre2$end)
      cres_here <- cres_here[CRE_label %in% cre_track_types | tested]
    }
  }
  region_res <- coaccess_region_result(res, cre1, cre2, labs_plot, cres = cres_here)

  peaks_df <- region_res$peaks_df
  peaks_variant <- "tested"
  if (cre_track == "by_class" && !is.null(cres_here) && nrow(cres_here)) {
    pk <- as.data.frame(cres_here)
    cls <- factor(pk$CRE_label, levels = c(intersect(CCRE_CLASSES, pk$CRE_label),
                                           setdiff(unique(pk$CRE_label), CCRE_CLASSES)))
    pk <- pk[order(cls, pk$start), , drop = FALSE]
    peaks_df <- data.frame(chr = pk$chrom, start = pk$start, end = pk$end,
                           group = as.character(pk$CRE_label))
    peaks_variant <- "intervals"
  }

  sub <- sprintf("%s:%d-%d", region$chrom, region$start, region$end)
  if (!is.null(subtitle)) sub <- paste0(sub, "  |  ", subtitle)
  if (show_2x2)
    sub <- paste0(sub, "\n", sprintf(paste0("%d fibers span both cCREs: both %d | CRE1 only %d | ",
                                            "CRE2 only %d | neither %d | Fisher OR (table + 1) = %.2f, p = %.2g"),
      counts$n_shared, counts$both, counts$cre1_only, counts$cre2_only, counts$neither,
      counts$fisher_or, counts$fisher_p))
  if (n_drawn < n_shared)
    sub <- paste0(sub, sprintf("\ndrawn: %d of %d fibers (at most %d per facet, seed %d)",
                               n_drawn, n_shared, max_fibers_per_facet, seed))

  n_facets <- length(res$samples)
  if (is.null(heights)) {
    raster <- max(2, 0.025 * n_drawn + 0.3 * n_facets)
    n_rows <- if (peaks_variant == "intervals") length(unique(peaks_df$group)) else 1
    heights <- c(pileup = max(1.5, 1.2 * n_facets), peaks = 0.3 + 0.3 * n_rows,
                 reads = raster, fire_fiberHMM = raster)
  }
  if (is.null(colors)) colors <- cluster_palette(res$samples)
  if (is.null(y_labels))
    y_labels <- c(pileup = "m6A fraction", reads = "m6A", fire_fiberHMM =
      if (nuc_source == "fiberhmm") "FiberHMM footprints" else if (tf_source == "none")
        "ft nucleosomes + FIRE" else "ft nucleosomes + FIRE + FiberHMM TF")
  draw_order <- c("linker", "nucleosome", "FIRE element",
                  if (tf_source != "none") rev(region_res$size_levels))

  read_theme <- list(axis.ticks.y = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank())
  p <- plot_region_panels(region_res, plot_name = plot_name, subtitle = sub,
    peaks_df = peaks_df, highlight = region_res$highlight,
    linewidth = linewidth, expand = TRUE, heights = heights,
    profile_options = list(colour_by = "group", colors = colors, y_col = "frac",
      linewidth = 0.4, facet_scales = "fixed", facet_space = "fixed", require_data = TRUE),
    read_options = list(tile_width = 8, tile_height = 0.85, backbone_col = BACKBONE_COL,
      arrows = FALSE, base_colors = c(A = M6A_COL), require_data = TRUE),
    peaks_options = list(variant = peaks_variant),
    footprint_options = list(draw_order = draw_order, backbone = FALSE, always = TRUE,
      legend_title = NULL, legend_drop = FALSE, legend_linewidth = 3),
    theme_options = list(base = "bw", base_size = 11,
      pileup = list(strip.text.y = ggplot2::element_text(angle = 0)),
      reads = read_theme, fire_fiberHMM = read_theme),
    y_labels = y_labels,
    panel_names = c(peaks = "cres", pileup = "m6a_prop", reads = "m6a", fire_fiberHMM = "reads"))

  # configuration bars from every fiber covering both cCREs
  bar_labels <- if (bar_by == "sample") labs else labs_facet
  bar_levels <- if (bar_by == "sample") samples else res$samples
  if (pooled_bar && length(bar_levels) > 1) {
    pooled <- data.table::copy(bar_labels[!is.na(config)])[, sample_name := "pooled"]
    bar_labels <- rbind(pooled, bar_labels)
    bar_levels <- c("pooled", bar_levels)
  }
  config_inputs <- configuration_proportion_inputs(bar_labels, bar_levels)
  if (is.null(bar_title))
    bar_title <- paste0(if (is.null(plot_name)) "cCRE pair" else plot_name,
                        ": configuration by ", if (bar_by == "sample") "sample" else facet_by)
  bar_theme <- if (length(bar_levels) > 6) theme_fiberseq("bw", base_size = 11,
    overrides = list(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))) else
    theme_fiberseq("bw", base_size = 11)
  bars <- plot_stacked_proportion(config_inputs$rows, "sample_name", "config", CONFIG_COLS,
    position = "fill", width = .7, legend = list(drop = FALSE),
    totals = list(data = config_inputs$totals, x = "sample_name", y = 1.04, label = "lab",
      style = list(size = 3)),
    scales = list(ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1.08),
      breaks = c(0, .25, .5, .75, 1))),
    labels = list(x = NULL, y = "fibers spanning both cCREs", fill = "configuration",
      title = bar_title),
    theme = bar_theme)

  list(res = res, labels = labs, counts = counts, panels = p, bars = bars, region = region,
       height = sum(heights) + 1.5)
}


# ---------------------------------------------------------------------------
# One-call read-level region plot: extract the region when a script is given
# (cached per sample, see extract_region_results()), load it with
# load_region_results(), optionally regroup the reads and draw it with
# plot_region_panels(). Needs parsing_footprints_functions.r.
#
# Inputs:
#   region            - one row with chr, start, end (1-based inclusive window)
#   sample_names      - sample labels
#   out_root, outname - results live in <out_root>/<outname>/<sample>/parsed
#   script            - extract_region_result_macrophage.sh /
#                       extract_region_result_lcl.sh; NULL uses existing folders
#   regenerate        - re-run the script even when the outputs exist
#   load_options      - further load_region_results() arguments
#   group             - optional function(res) -> res, e.g. calling
#                       add_read_groups() or add_topic_clusters()
#   panel_options     - plot_region_panels() arguments after res, as a list or
#                       a function(res) returning one (e.g. for a subtitle built
#                       from the grouped reads)
# Output:
#   list(res, plot = plot_region_panels() result, result_dirs)
# ---------------------------------------------------------------------------
plot_region_example <- function(region, sample_names, out_root, outname, script = NULL,
                                regenerate = FALSE, load_options = list(), group = NULL,
                                panel_options = list()) {
  if (!exists("load_region_results", mode = "function"))
    stop("source parsing_footprints_functions.r before plot_region_example()")
  result_dirs <- if (is.null(script)) file.path(out_root, outname, sample_names, "parsed") else
    extract_region_results(sample_names, region, outname, out_root, script, regenerate = regenerate)
  res <- do.call(load_region_results, c(list(region = region, sample_names = sample_names,
                                              result_dirs = result_dirs), load_options))
  if (!is.null(group)) res <- group(res)
  if (is.function(panel_options)) panel_options <- panel_options(res)
  list(res = res, plot = do.call(plot_region_panels, c(list(res), panel_options)),
       result_dirs = result_dirs)
}


# ---------------------------------------------------------------------------
# Per-cluster accessibility, molecule and heatmap figures for one set of TF
# footprint example regions (panel3C / panel3BC). For each example: an
# accessibility profile + molecule panel per Leiden cluster, saved alone,
# then combined side-by-side with the pooled read heatmap from
# plot_read_heatmap(layout = "fire"), saved combined. Uses cluster_palette(),
# plot_group_profile(), theme_fiberseq() and plot_read_heatmap() from this
# file; the remaining collaborators are supplied by the caller because they
# come from the calling pipeline (an extracted 05_vis_fire_freq.Rmd chunk and
# its own R-stage preamble), not from this file.
#
# Inputs:
#   examples              - one row per example region (selection_TF, start,
#                           end, ...), see reference_result()
#   reference_result      - function(ex) -> clustering result with
#                           $assignments (RID, cluster, time, start, sample_name)
#                           and site_met_mat, as built by the caller
#   show_reads            - function(result) printing a read/cluster summary
#   example_dir           - function(ex) -> folder with retained_molecule_tracks.rds
#   example_figure_path   - function(ex, panel, description) -> output PDF path
#   plot_footprint_panels - function(accessibility, molecules, ex, n_bars,
#                           n_groups, n_molecules, subtitle) combining an
#                           accessibility profile with a molecule panel
#   molecule_tracks, add_read_groups, stack_rows, plot_molecules -
#                           parsing_footprints_functions.r helpers for
#                           per-read track data and the molecule panel
#   time_cols, time_levels - timepoint colours and factor levels
#   example_theme_settings - theme_fiberseq() overrides shared across panels
#   window_breaks         - function(ex) -> x axis breaks for ex's window
#   example_title         - function(ex) -> figure title
# Output:
#   invisible NULL; for each example, saves the panel3C cluster_footprints
#   and panel3BC combined PDFs via example_figure_path()
# ---------------------------------------------------------------------------
plot_tf_cluster_profiles <- function(examples, reference_result, show_reads, example_dir,
                                     example_figure_path, plot_footprint_panels,
                                     molecule_tracks, add_read_groups, stack_rows, plot_molecules,
                                     time_cols, time_levels, example_theme_settings, window_breaks,
                                     example_title) {
  for (i in seq_len(nrow(examples))) {
    ex <- examples[i, , drop = FALSE]
    result <- reference_result(ex)
    show_reads(result)
    a <- result$assignments
    cl_cols <- cluster_palette(levels(a$cluster))
    res_ex <- readRDS(file.path(example_dir(ex), "retained_molecule_tracks.rds"))
    stopifnot(setequal(as.character(res_ex$rids_df$RID), a$RID))
    tracks <- molecule_tracks(res_ex)
    rows_C <- a[order(a$cluster, a$time, a$start, a$RID), ]
    res_cl <- add_read_groups(res_ex, labels = setNames(as.character(a$cluster), a$RID),
                              group_col = "cluster", group_levels = levels(a$cluster), verbose = FALSE)
    rows_C2 <- stack_rows(rows_C, "cluster")
    n_C2 <- table(rows_C2$cluster)
    lab_C2 <- setNames(sprintf("%s\nn = %d", names(n_C2), as.integer(n_C2)), names(n_C2))
    pC2_mol <- plot_molecules(rows_C2, tracks, ex, "cluster", lab_C2,
                              bar_cols = c("cluster", "time"),
                              bar_palette = c(cl_cols, time_cols), legend_extra = time_levels)
    pC2_accessibility <- local({
      .profile_args <- list(pileup = res_cl$pileup,
    ex = ex,
    group_col = "cluster",
    group_levels = levels(a$cluster),
    group_cols = cl_cols,
    n_reads = n_C2)
      pileup <- .profile_args$pileup
      ex <- .profile_args$ex
      group_col <- .profile_args$group_col
      group_levels <- .profile_args$group_levels
      group_cols <- .profile_args$group_cols
      n_reads <- .profile_args$n_reads
          stopifnot(all(group_levels %in% names(n_reads)))
          pileup$accessibility_group <- factor(as.character(pileup[[group_col]]),
                                                levels = group_levels)
          stopifnot(!anyNA(pileup$accessibility_group))
          strip <- setNames(sprintf("%s\nn = %d", group_levels,
                                    as.integer(n_reads[group_levels])), group_levels)
          plot_group_profile(pileup, group_col = "accessibility_group", style = "area",
            y_col = "smooth_frac", mapping = aes(pos, smooth_frac),
            highlight = list(xmin = ex$start, xmax = ex$end, fill = "grey88", position = "under"),
            layers = list(list(geom = "area", mapping = aes(fill = accessibility_group), alpha = 0.9),
              list(geom = "line", linewidth = 0.25, colour = "grey15"))) +
            facet_grid(accessibility_group ~ ., labeller = as_labeller(strip)) +
            scale_fill_manual(values = group_cols, guide = "none") +
            scale_x_continuous(breaks = window_breaks(ex), labels = scales::comma) +
            scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
            labs(x = NULL, y = "m6A fraction (+/- 5 bp mean)") +
            theme_fiberseq("classic", base_size = 9, overrides = example_theme_settings) +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                  axis.line.x = element_blank())

    })
    pC2 <- plot_footprint_panels(pC2_accessibility, pC2_mol, ex, n_bars = 2L,
      n_groups = nlevels(a$cluster), n_molecules = nrow(rows_C2),
      subtitle = sprintf(paste0("Accessibility by Leiden cluster; %s motif positions; tested FIRE; ",
                                "all %d full-span molecules by cluster\n",
                                "Left bars = cluster, time; molecules ordered by time then read start; ",
                                "dashed lines = tested region"), ex$selection_TF, nrow(rows_C2)))
    h_C2 <- pC2$height
    path_C2 <- example_figure_path(ex, "panel3C", "cluster_footprints")
    ggsave(path_C2, pC2$plot, width = 8.5, height = h_C2, limitsize = FALSE)
    print(head(res_cl$pileup))
    print(head(rows_C2[, c("RID", "cluster", "time", "start", "y")]))
    print(n_C2)

    h_C1 <- min(14, 3 + 0.02 * nrow(rows_C))
    h_C <- max(h_C1, h_C2)
    g_C1 <- grid::grid.grabExpr(plot_read_heatmap(result, layout = "fire", rows = rows_C, window = ex,
    cluster_colors = cl_cols, timepoint_cols = time_cols, ticks = window_breaks(ex),
    main = example_title(ex), draw = TRUE, use_raster = FALSE),
                                width = 9, height = h_C)
    pC <- cowplot::plot_grid(
      cowplot::plot_grid(NULL, NULL, nrow = 1, rel_widths = c(9, 8.5), labels = c("C1", "C2")),
      cowplot::plot_grid(g_C1, patchwork::patchworkGrob(pC2$plot), nrow = 1, rel_widths = c(9, 8.5)),
      ncol = 1, rel_heights = c(0.3, h_C))
    path_C <- example_figure_path(ex, "panel3BC", "combined")
    ggsave(path_C, pC, width = 17.5, height = h_C + 0.3, limitsize = FALSE)
    cat("Saved", path_C2, "and", path_C, "\n")
  }
  invisible(NULL)
}
