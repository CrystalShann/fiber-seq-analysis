source(file.path(snakemake@params[["scripts_dir"]], "common.R"))
cfg <- snakemake@config
dataset <- snakemake@wildcards[["dataset"]]
ds <- cfg$datasets[[dataset]]
load_shared(cfg, plots = TRUE)
source("/project/spott/cshan/fiber-seq/code/parsing_functions/plotting_functions.r", local = TRUE)
source(file.path(snakemake@params[["scripts_dir"]], "plot_helpers.R"))
source(file.path(snakemake@params[["scripts_dir"]], "plot_registry.R"))

result <- readRDS(snakemake@input[["clustering"]])
assembled <- readRDS(snakemake@input[["assembled"]])
plot_dir <- snakemake@params[["plot_dir"]]
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
# A Snakemake directory output belongs wholly to this rule. Remove stale PDFs
# after a panel list changes so the directory contains exactly the selected set.
stale <- list.files(plot_dir, pattern = "[.]pdf$", full.names = TRUE)
if (length(stale)) unlink(stale)
context <- workflow_plot_context(result, assembled, cfg, ds,
                                 dirname(snakemake@input[["clustering"]]))
panels <- unlist(context$ds$plots, use.names = FALSE)
manifest <- file.path(plot_dir, "panels.tsv")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
if (anyDuplicated(panels)) stop("Dataset plots contains duplicate panel names")
missing <- setdiff(panels, names(panel_registry))
if (length(missing)) stop("Unknown plot panels: ", paste(missing, collapse = ", "))
records <- list()
# Layout helpers may open a device; keep incidental Rplots.pdf out of outputs.
grDevices::pdf(NULL)
layout_device <- grDevices::dev.cur()
tryCatch({
for (name in panels) {
  entry <- panel_registry[[name]]
  filename <- entry$filename(context$result$region$region_id)
  if (!is.null(entry$enabled) && !entry$enabled(context)) {
    records[[name]] <- data.frame(panel = name, filename = filename, status = "outside_rank_limit")
    next
  }
  descriptor <- entry$render(context)
  path <- file.path(plot_dir, filename)
  save_figure(descriptor$plot, path, width = descriptor$width, height = descriptor$height,
    device = if (is.null(descriptor$draw)) grDevices::cairo_pdf else grDevices::pdf,
    draw = descriptor$draw, limitsize = FALSE,
    bg = if (is.null(descriptor$draw)) "white" else NULL)
  stopifnot(file.exists(path), file.info(path)$size > 0)
  records[[name]] <- data.frame(panel = name, filename = filename, status = descriptor$status)
}
}, finally = if (layout_device %in% grDevices::dev.list()) grDevices::dev.off(layout_device))
tab <- if (length(records)) do.call(rbind, records) else
  data.frame(panel = character(), filename = character(), status = character())
utils::write.table(tab, manifest, sep = "\t", quote = FALSE, row.names = FALSE)
invisible(tab)
