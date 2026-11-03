#!/usr/bin/env Rscript

script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1L])
here <- dirname(normalizePath(script))
args <- commandArgs(trailingOnly = TRUE)
opt <- list(config = file.path(here, "pol2_chip_cluster_utils", "config.yaml"), stage = "all", chrom = NULL, gene = NULL)
if ("--help" %in% args) {
  cat("00_run.sh [--config YAML] [--stage prepare|count|features|cluster|project|tracks|autocorrelation|nucleosome|expression|observed-sites|motif-inputs|histone-inputs|render|promoter-motifs|motif-bars|report|all] [--chrom chrN] [--gene ARF4]\n")
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
source(file.path(here,"pol2_chip_cluster_utils/pol2_observed_profiles.R"))
if (!opt$stage %in% c("prepare", "count", "features", "cluster", "project", "tracks", "autocorrelation", "nucleosome", "expression", "observed-sites", "motif-inputs", "histone-inputs", "render", "promoter-motifs", "motif-bars", "report", "all")) stop("Invalid stage")
if (!is.null(opt$chrom) && opt$stage != "features") stop("--chrom is only for the features stage")
if (!is.null(opt$gene) && !opt$stage %in% c("promoter-motifs","observed-sites","motif-inputs")) stop("--gene is only for promoter-motifs, observed-sites or motif-inputs")
dir.create(cfg$output, recursive = TRUE, showWarnings = FALSE)
if (opt$stage == "observed-sites") {
  x <- readRDS(file.path(cfg$output,"tables/clustering.rds"))
  masks <- pol2_observed_masks(cfg,x,build=TRUE,genes=opt$gene)
  message("Prepared observed-site caches for ",nrow(masks)," promoters")
  quit(status=0L)
}
if (opt$stage %in% c("motif-inputs","histone-inputs")) {
  source(file.path(here,"08_prepare_report_inputs.R"),local=TRUE)
  quit(status=0L)
}
if (opt$stage == "render") {
  rmarkdown::render(file.path(here,"06_pol2_clustering_report.Rmd"),
    output_file="pol2_clustering_report.html",output_dir=cfg$output,
    params=list(output=cfg$output),envir=new.env(parent=globalenv()),quiet=FALSE)
  quit(status=0L)
}

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
  report <- readLines(file.path(here,"06_pol2_clustering_report.Rmd"))
  label <- paste0(opt$stage,"-analysis")
  start <- grep(paste0("^```\\{r ",label,","),report)
  stopifnot(length(start)==1L)
  end <- start + which(report[(start+1L):length(report)]=="```")[1L]
  recompute_autocorrelation <- opt$stage == "autocorrelation"
  recompute_nucleosome <- opt$stage == "nucleosome"
  eval(parse(text=report[(start+1L):(end-1L)]),envir=new.env())
  if (opt$stage == "autocorrelation") pol2_autocorrelation_plots(cfg, save = TRUE)
  if (opt$stage == "nucleosome") pol2_nucleosome_plots(cfg, save = TRUE)
}
if (opt$stage %in% c("report", "all")) pol2_report(cfg, here)

if (opt$stage == "motif-bars") source(file.path(here,"07_motif_cluster_bars.R"))
if (opt$stage == "promoter-motifs") {
  report <- readLines(file.path(here,"06_pol2_clustering_report.Rmd"))
  params <- list(output=cfg$output, example_gene=opt$gene)
  for (label in c("setup","promoter-motif-selection","promoter-motif-analysis")) {
    start <- grep(paste0("^```\\{r ",label,","),report)
    stopifnot(length(start)==1L)
    end <- start + which(report[(start+1L):length(report)]=="```")[1L]
    eval(parse(text=report[(start+1L):(end-1L)]))
  }
}
