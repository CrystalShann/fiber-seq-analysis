source(file.path(snakemake@params[["scripts_dir"]], "common.R"))
cfg <- snakemake@config
ds <- cfg$datasets[[snakemake@wildcards[["dataset"]]]]
load_shared(cfg)
samples <- read_samples(snakemake@input[["samples"]])
region <- read_region(snakemake@input[["regions"]], snakemake@wildcards[["region_id"]])
selected <- region_samples(region, samples)
paths <- vapply(seq_len(nrow(selected)), function(i) {
  values <- as.list(selected[i, , drop = FALSE])
  values$chr <- region$chr
  expand_template(ds$m6a, values)
}, character(1))
dat <- assemble_region_m6a(sample_table = selected, region = region,
  full_span = TRUE, matrix_dir = NULL, reuse = FALSE, m6a_paths = paths)
for (column in setdiff(names(selected), names(dat$rids_df))) {
  dat$rids_df[[column]] <- selected[[column]][match(dat$rids_df$sample_name, selected$sample_name)]
}
dat$region <- region
dat$sample_table <- samples
dat$selected_samples <- selected
dat$m6a_paths <- paths
save_rds(dat, snakemake@output[["raw"]])
