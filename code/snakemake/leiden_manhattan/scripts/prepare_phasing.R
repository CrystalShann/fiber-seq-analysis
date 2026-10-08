source(file.path(snakemake@params[["scripts_dir"]], "common.R"))
cfg <- snakemake@config
load_shared(cfg)
samples <- read_samples(snakemake@input[["samples"]])
if (!all(c("phasing_dir", "cell_line") %in% names(samples)))
  stop("Phased datasets require phasing_dir and cell_line in the sample sheet")
# The helper needs only read metadata here, not all the matrices simultaneously.
region_data <- lapply(snakemake@input[["raw"]], function(path) list(rids_df = readRDS(path)$rids_df))
cache <- cache_lcl_haplotags(samples, region_data, phasing_root = NULL,
  cache_dir = snakemake@params[["cache_dir"]], reuse = FALSE)
save_rds(cache, snakemake@output[["phase"]])
