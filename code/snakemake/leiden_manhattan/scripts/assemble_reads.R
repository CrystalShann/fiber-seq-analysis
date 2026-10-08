source(file.path(snakemake@params[["scripts_dir"]], "common.R"))
source(file.path(snakemake@params[["scripts_dir"]], "footprint_io.R"))
cfg <- snakemake@config
ds <- cfg$datasets[[snakemake@wildcards[["dataset"]]]]
load_shared(cfg, footprints = TRUE)
dat <- readRDS(snakemake@input[["raw"]])
if ("phased_heterozygous" %in% unlist(ds$read_filter)) {
  phase <- readRDS(snakemake@input[["phase"]][[1L]])
  dat <- lcl_filter_focal_heterozygotes(dat, dat$region, dat$selected_samples, phase)
}
tracks <- collect_footprints(ds, cfg, dat$selected_samples, dat$region, dat$rids_df)
dat$footprints <- tracks$records
dat$track_spec <- tracks[setdiff(names(tracks), "records")]
stopifnot(!anyNA(dat$met_mat), !anyDuplicated(dat$rids_df$RID),
  identical(rownames(dat$met_mat), as.character(dat$rids_df$RID)),
  all(dat$rids_df$start <= dat$region$analysis_start),
  all(dat$rids_df$end >= dat$region$analysis_end))
save_rds(dat, snakemake@output[["assembled"]])
