#!/usr/bin/env Rscript
# Plot saved pooled assignments and signals without extraction or clustering.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(Rsamtools)
  library(grid)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(ComplexHeatmap)
})
project_root <- Sys.getenv("ENHANCER_PROJECT_ROOT", "/project/spott/cshan/fiber-seq")
table_dir <- file.path(project_root, "macrophage_project/enhancer/TF_co-occ/tables")
result_file <- Sys.getenv("ENHANCER_MANHATTAN_OUTPUT",
  file.path(table_dir, "enhancer_manhattan_results_pooled_capped.rds"))
plot_dir <- Sys.getenv("ENHANCER_MANHATTAN_PLOT_DIR",
  file.path(dirname(dirname(result_file)), "plots"))
source(file.path(project_root, "code/enhancer/manhattan/enhancer_manhattan_plots.R"))
source(file.path(project_root, "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r"))
saved <- readRDS(result_file)
stopifnot(!is.null(saved$results$pooled))
save_enhancer_manhattan_pdfs(saved$results$pooled, "pooled", plot_dir,
  dirname(result_file), expected_k = 100L, suffix = "_capped")
message("Replotted saved pooled Manhattan results: ", plot_dir)
