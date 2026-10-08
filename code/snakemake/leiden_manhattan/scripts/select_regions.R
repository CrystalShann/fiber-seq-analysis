source(file.path(snakemake@params[["scripts_dir"]], "common.R"))
cfg <- snakemake@config
load_shared(cfg)
dataset <- snakemake@wildcards[["dataset"]]
ds <- cfg$datasets[[dataset]]
for (script in c("regions_tss.R", "regions_custom.R", "regions_bed.R", "regions_asfire_top.R"))
  source(file.path(snakemake@params[["scripts_dir"]], script))
selectors <- list(tss = select_tss_regions, custom = select_custom_regions,
                  bed = select_bed_regions,
                  asfire_top = select_asfire_top_regions)
unknown <- setdiff(names(ds$regions), names(selectors))
if (length(unknown)) stop("Unknown region type(s): ", paste(unknown, collapse = ", "))
regions <- dplyr::bind_rows(lapply(names(ds$regions), function(type) selectors[[type]](ds$regions[[type]])))
if (!nrow(regions)) stop("No regions selected for ", dataset)
if (anyNA(regions$region_id) || anyDuplicated(regions$region_id) ||
    any(!grepl("^[A-Za-z0-9_.-]+$", regions$region_id)))
  stop("Region IDs must be unique nonempty names containing letters, numbers, '_', '.', or '-'")
stopifnot(all(regions$analysis_start == regions$start + 1L),
          all(regions$analysis_end == regions$end), all(regions$width == regions$end - regions$start),
          all(regions$analysis_start > 0), all(regions$width > 0))
front <- c("region_id", "chr", "start", "end", "analysis_start", "analysis_end")
regions <- regions[, c(front, setdiff(names(regions), front)), drop = FALSE]
write_tsv(regions, snakemake@output[["regions"]])
