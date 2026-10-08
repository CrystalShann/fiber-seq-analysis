source(file.path(snakemake@params[["scripts_dir"]], "common.R"))
cfg <- snakemake@config
load_shared(cfg)
dat <- readRDS(snakemake@input[["assembled"]])
region <- dat$region
wc <- snakemake@wildcards
result <- leiden_manhattan_cluster(dat$met_mat, dat$rids_df,
  region_start = region$analysis_start, region_end = region$analysis_end,
  window_size = as.integer(wc[["window_size"]]), k_neighbors = as.integer(wc[["k"]]),
  resolution = as.numeric(wc[["resolution"]]), seed = as.integer(wc[["seed"]]),
  sigma = cfg$clustering$kernel_sigma)
meta <- dat$rids_df[match(result$assignments$RID, dat$rids_df$RID), , drop = FALSE]
for (column in setdiff(names(meta), names(result$assignments))) result$assignments[[column]] <- meta[[column]]
result$assignments$sample_name <- factor(result$assignments$sample_name, levels = dat$sample_table$sample_name)
result$assignments$region_id <- region$region_id
result$assignments$annotation <- region$annotation
result$region <- region
result$met_mat <- result$site_met_mat <- dat$met_mat
result$variants <- dat$variants
result$focal <- dat$focal
result$read_filter_audit <- dat$read_filter_audit
result$qc <- dat$qc
result$track_spec <- dat$track_spec
result$params$region_chr <- region$chr
result$params$region_id <- region$region_id
result$params$coordinate_system <- "1-based inclusive"
result$analysis_signature <- list(region = region, samples = dat$selected_samples,
  read_filter = cfg$datasets[[wc[["dataset"]]]]$read_filter,
  m6a_paths = dat$m6a_paths, params = result$params,
  R = R.version.string, igraph = as.character(utils::packageVersion("igraph")))
out <- snakemake@output
save_rds(result, out[["clustering"]])
write_tsv(result$assignments, out[["assignments"]])
edges <- igraph::as_data_frame(result$graph, what = "edges")
stopifnot(all(c("weight", "manhattan_distance") %in% names(edges)),
  all(is.finite(edges$manhattan_distance)), all(edges$manhattan_distance >= 0))
write_tsv(edges, out[["edges"]])
write_tsv(data.frame(cluster = rownames(result$profiles), result$profiles, check.names = FALSE), out[["profiles"]])
write_tsv(result$window_anno, out[["windows"]])
front <- intersect(c("RID", "sample_name", "timepoint", "cluster"), names(result$assignments))
write_tsv(data.frame(result$assignments[, front, drop = FALSE], result$feat_mat, check.names = FALSE), out[["features"]])
write_tsv(data.frame(parameter = names(result$params),
  value = vapply(result$params, as.character, character(1))), out[["parameters"]])
