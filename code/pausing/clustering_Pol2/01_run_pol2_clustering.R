#!/usr/bin/env Rscript
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1L])
here <- dirname(normalizePath(script))
args <- commandArgs(trailingOnly = TRUE)
opt <- list(config = file.path(here, "pol2_chip_cluster_utils", "config.yaml"), stage = "all", chrom = NULL)
if ("--help" %in% args) {
  cat("00_run.sh [--config YAML] [--stage prepare|count|features|cluster|project|tracks|autocorrelation|nucleosome|expression|report|all] [--chrom chrN]\n")
  quit(status = 0)
}
if (length(args) %% 2) stop("Expected --option value pairs")
for (i in seq_len(length(args) %/% 2) * 2L - 1L) {
  key <- sub("^--", "", args[i])
  if (!key %in% names(opt)) stop("Unknown option: ", args[i])
  opt[[key]] <- args[i + 1L]
}
source(file.path(here, "pol2_chip_cluster_utils", "pol2_clustering_functions.R"))
source(file.path(here, "pol2_chip_cluster_utils", "pol2_molecule_tables.R"))
source(file.path(here, "pol2_chip_cluster_utils", "pol2_clustering_plots.R"))
cfg <- yaml::read_yaml(opt$config)
pol2_load_helpers(cfg$project)
if (!opt$stage %in% c("prepare", "count", "features", "cluster", "project", "tracks", "autocorrelation", "nucleosome", "expression", "report", "all")) stop("Invalid stage")
if (!is.null(opt$chrom) && opt$stage != "features") stop("--chrom is only for the features stage")
dir.create(cfg$output, recursive = TRUE, showWarnings = FALSE)
if (opt$stage == "count") {
  source(file.path(here, "03_count_pol2_eligible.R"), local = TRUE)
  quit(status = 0L)
}
if (opt$stage %in% c("prepare", "all")) source(file.path(here,"02_select_tss.R"),local=TRUE)
if (opt$stage %in% c("features", "all")) {
  chroms <- if (is.null(opt$chrom)) cfg$chromosomes else opt$chrom
  if (!all(chroms %in% cfg$chromosomes)) stop("Unknown chromosome")
  for (cc in chroms) pol2_features(cfg, cc)
}
if (opt$stage %in% c("cluster", "all")) pol2_cluster(cfg)
if (opt$stage %in% c("project", "all")) pol2_project(cfg)
if (opt$stage %in% c("tracks", "all")) pol2_tracks(cfg)

if (opt$stage %in% c("autocorrelation","nucleosome","expression")) {
  report <- readLines(file.path(here,"05_pol2_clustering_report.Rmd"))
  label <- paste0(opt$stage,"-analysis")
  start <- grep(paste0("^```\\{r ",label,","),report)
  stopifnot(length(start)==1L)
  end <- start + which(report[(start+1L):length(report)]=="```")[1L]
  recompute_autocorrelation <- opt$stage == "autocorrelation"
  recompute_nucleosome <- opt$stage == "nucleosome"
  eval(parse(text=report[(start+1L):(end-1L)]),envir=new.env())
}
if (opt$stage %in% c("report", "all")) pol2_report(cfg, here)
