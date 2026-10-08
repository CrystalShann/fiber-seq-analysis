source(file.path(snakemake@params[["scripts_dir"]], "common.R"))
cfg <- snakemake@config
dataset <- snakemake@wildcards[["dataset"]]
ds <- cfg$datasets[[dataset]]
load_shared(cfg, plots = TRUE)
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
save_all_plots(context, plot_dir)
