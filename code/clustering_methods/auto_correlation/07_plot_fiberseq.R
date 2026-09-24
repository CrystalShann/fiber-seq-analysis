#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(6L, 7L))
read_footprints_only <- length(args) == 7L
if (read_footprints_only) stopifnot(args[7] == "--read-footprints-only")
project <- args[1]
output <- normalizePath(args[2], mustWork = TRUE)
rid <- args[3]
seed <- as.integer(args[4])
resolution <- as.numeric(args[5])
k <- as.integer(args[6])
suppressPackageStartupMessages({
  library(dplyr); library(Matrix); library(GenomicRanges); library(Rsamtools)
  library(ggplot2); library(ComplexHeatmap); library(jsonlite)
})
source(file.path(project, "code/topic_model/topic_modelling_functions.r"))
source(file.path(project, "code/haplotype_phasing/LCL_phasing.r"))
source(file.path(project, "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r"))
data.table::setDTthreads(1L)
options(scipen = 999)
message(R.version.string)
message(paste(c("ComplexHeatmap", "ggplot2", "igraph"),
  vapply(c("ComplexHeatmap", "ggplot2", "igraph"), function(p) as.character(packageVersion(p)), ""), collapse = "; "))

# Shared export helper requires this path guard; scope it to this new run only.
lcl_output_path <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (!(identical(path, output) || startsWith(path, paste0(output, "/"))))
    stop("Plot output must stay inside the autocorrelation run")
  path
}

# Every input read is already phased and heterozygous before ACF computation.
# Use the original strict LCL display helper, including its allele assertions.
original_display_result <- fiberseq_display_result

payload <- jsonlite::fromJSON(paste(readLines(file("stdin"), warn = FALSE), collapse = "\n"), simplifyMatrix = TRUE)
region <- payload$region
stopifnot(is.data.frame(region), nrow(region) == 1L, region$region_id == rid)
region$annotation <- paste0(region$annotation, " | autocorrelation")
a <- payload$records
stopifnot(is.data.frame(a), all(a$region_id == rid))
a <- a[a$status == "clustered", , drop = FALSE]
a$original_RID <- a$RID
a$RID <- paste(a$sample_name, a$original_RID, sep = "::")
a$start <- a$read_start
a$end <- a$read_end
a$chr <- region$chr
a$sample_label <- sub("_.*$", "", a$sample_name)
a$cluster <- factor(a$cluster, levels = as.character(sort(unique(as.integer(a$cluster)))))
# Keep each ACF cluster's color stable when a display subset lacks a cluster.
all_cluster_levels <- levels(a$cluster)
original_cluster_palette <- cluster_id_palette
cluster_id_palette <- function(levels) original_cluster_palette(all_cluster_levels)[levels]
met <- graph <- NULL
if (!read_footprints_only) {
matched <- match(a$row_id, payload$row_ids)
stopifnot(!anyNA(matched), !anyDuplicated(a$RID), length(matched) == length(payload$row_ids))
offsets <- if (is.matrix(payload$call_offsets)) asplit(payload$call_offsets, 1L) else payload$call_offsets
offsets <- offsets[matched]
met <- Matrix::sparseMatrix(i = rep(seq_len(nrow(a)), lengths(offsets)),
  j = as.integer(unlist(offsets)) + 1L, x = 1,
  dims = c(nrow(a), region$width),
  dimnames = list(a$RID, as.character(seq.int(region$analysis_start, region$analysis_end))))
acfs <- payload$acf[matched, , drop = FALSE]
stopifnot(nrow(acfs) == nrow(a), ncol(acfs) <= region$width,
  all(is.finite(acfs)), all(abs(acfs[, 1] - 1) < 1e-8))
}

leiden_root <- file.path(project, "LCL_project/Leiden_manhattan")
focal <- data.frame(chr = region$chr, pos = region$focal_pos, variant_id = region$focal_snp,
  ref = region$ref, alt = region$alt)
stopifnot(all(a$haplotype %in% c("HP1", "HP2")),
  all(a$focal_genotype %in% c("0|1", "1|0")),
  all(a$allele_status == "phased_focal_genotype"))
a$heatmap_included <- TRUE

if (!read_footprints_only) {
edges <- as.data.frame(payload$edges)
edges$source <- a$RID[match(edges$source, a$row_id)]
edges$target <- a$RID[match(edges$target, a$row_id)]
stopifnot(!anyNA(edges), all(edges$weight > 0))
graph <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = a$RID)
}
result <- list(assignments = a, region = region, site_met_mat = met, graph = graph,
  n_clusters = nlevels(a$cluster), variants = focal, focal = focal,
  params = list(seed = seed, k_eff = k, resolution = resolution))
result <- fiberseq_display_result(result)
a <- result$assignments
short_colors <- sample_palette(a$sample_name)
full_colors <- setNames(unname(short_colors[sub("_.*$", "", unique(a$sample_name))]), unique(a$sample_name))
plot_dir <- file.path(output, "outputs", rid)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
save_gg <- function(plot, name, width, height) {
  ggplot2::ggsave(file.path(plot_dir, paste0(name, ".pdf")), plot,
    width = width, height = height, limitsize = FALSE, bg = "white")
}

# Use the same phased-heterozygote input samples for nucleosome extraction.
ft_root <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results"
nuc_samples <- data.frame(sample_name = unique(a$sample_name),
  fire_dir = file.path(ft_root, unique(a$sample_name)))
nucleosomes <- extract_nucleosomes(nuc_samples, region, a)

# Reuse read-only regional footprint caches only when they cover this exact
# requested window. Wider windows stream the source BED instead of truncating.
tracks_root <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged/combined/joint_trained_tracks"
footprint_sources <- list()
tf_records <- lapply(c("10-30", "40-60", "60-80"), function(size) {
  filename <- paste0("combined_", region$chr, "_", size, "bp_fps.bed.gz")
  cache <- file.path(leiden_root, "footprint_summary", paste0(filename, ".rds"))
  source_file <- file.path(tracks_root, region$chr, filename)
  cached <- if (file.exists(cache)) readRDS(cache) else NULL
  bounds <- if (!is.null(cached)) cached$signature$regions else NULL
  covered <- !is.null(bounds) && any(bounds$chr == region$chr & bounds$start <= region$start & bounds$end >= region$end)
  if (covered) {
    records <- cached$records
    used <- cache
  } else {
    if (!file.exists(source_file)) stop("Missing footprint source: ", source_file)
    program <- sprintf('BEGIN {FS=OFS="\t"} $1 == "%s" && $2 < %.0f && $3 > %.0f {print}',
      region$chr, region$end, region$start)
    command <- paste("gzip -cd", shQuote(source_file), "| awk", shQuote(program))
    lines <- system2("/bin/bash", c("-o", "pipefail", "-c", shQuote(command)), stdout = TRUE)
    if (!is.null(attr(lines, "status")) && attr(lines, "status") != 0L)
      stop("Footprint extraction failed: ", source_file)
    records <- if (length(lines)) data.table::fread(text = paste(lines, collapse = "\n"),
      header = FALSE, data.table = FALSE) else data.frame()
    if (!nrow(records)) records <- data.frame(chr = character(), start = integer(), end = integer(), original_RID = character())
    stopifnot(ncol(records) == 4L)
    names(records) <- c("chr", "start", "end", "original_RID")
    used <- source_file
  }
  records <- records[records$chr == region$chr & records$start < region$end & records$end > region$start, , drop = FALSE]
  # A physical PacBio RID can occur in both an original and a merged AL input.
  # Match its same footprint to every retained sample/RID row, preserving cohort.
  records <- merge(records, a[, c("original_RID", "RID")], by = "original_RID", sort = FALSE)
  records$size <- records$end - records$start
  records$track <- rep(paste0("FiberHMM_", size, "bp"), nrow(records))
  footprint_sources[[size]] <<- data.frame(track = paste0("FiberHMM_", size, "bp"), source = used, n_records = nrow(records))
  records
})
footprints <- dplyr::bind_rows(nucleosomes, dplyr::bind_rows(tf_records))
reads_plot <- plot_smf_reads(result, footprints, short_colors, include_haplotype = TRUE) +
  ggplot2::labs(title = paste(region$annotation, "Autocorrelation-defined Leiden clusters", sep = " | "),
    subtitle = paste0("One row per phased heterozygous molecule; k = ", k,
      "; resolution = ", resolution))
save_gg(reads_plot, "fig1_read_footprints", 12,
  max(7, .03 * nrow(a) + 3.5 + .3 * result$n_clusters))
if (read_footprints_only) {
  message("Read footprints: ", rid, "; ", nrow(a), " phased reads; ", nrow(footprints), " footprint records")
  quit(save = "no", status = 0L)
}
signals <- signal_profiles(result, footprints)
save_heatmap <- function(plot, name, width, height) {
  grDevices::pdf(file.path(plot_dir, paste0(name, ".pdf")), width = width, height = height)
  tryCatch(ComplexHeatmap::draw(plot, newpage = FALSE), finally = grDevices::dev.off())
}
heat_keep <- a$heatmap_included
if (any(heat_keep)) {
heat_result <- result
heat_result$assignments <- droplevels(a[heat_keep, , drop = FALSE])
heat_result$site_met_mat <- met[heat_keep, , drop = FALSE]
heat_result$n_clusters <- nlevels(heat_result$assignments$cluster)
heat_result$region$annotation <- paste0(region$annotation, "\nPhased heterozygous display: ",
  sum(heat_keep), " / ", nrow(a), " clustered reads")
# The original strict display helper verifies HP-linked focal alleles and
# removes Unknown from the legend, just as in the LCL Leiden notebook.
heat_result <- original_display_result(heat_result)
heat_a <- heat_result$assignments
heat_colors <- short_colors[names(short_colors) %in% sub("_.*$", "", heat_a$sample_name)]
met_heatmap <- function() plot_genomic_cluster_heatmap(heat_result, heat_result$region, heat_colors,
  include_haplotype = TRUE, variants = heat_result$variants, show_cluster_profiles = TRUE,
  cluster_label = "Autocorrelation-defined Leiden clusters")
save_heatmap(met_heatmap(), "heatmap_m6a_footprints", 14, 11 + .55 * heat_result$n_clusters)
ordering <- order(heat_a$cluster, heat_a$start, heat_a$RID)
acf_display <- acfs[heat_keep, -1, drop = FALSE][ordering, , drop = FALSE]
rownames(acf_display) <- heat_a$RID[ordering]
limit <- max(as.numeric(quantile(abs(acf_display), .99)), .01)
ticks <- unique(round(seq(1, ncol(acf_display), length.out = min(5L, ncol(acf_display)))))
acf_heatmap <- ComplexHeatmap::Heatmap(acf_display, name = "ACF",
  col = circlize::colorRamp2(c(-limit, 0, limit), c("#2166ac", "white", "#b2182b")),
  cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE,
  row_split = heat_a$cluster[ordering], show_row_names = FALSE, show_column_names = FALSE,
  use_raster = TRUE, raster_resize_mat = FALSE,
  column_title = "Autocorrelation: lag (bp)",
  bottom_annotation = ComplexHeatmap::HeatmapAnnotation(lag_bp = ComplexHeatmap::anno_mark(
    at = ticks, labels = as.character(ticks), which = "column", side = "bottom")),
  width = grid::unit(65, "mm"))
save_heatmap(acf_heatmap + met_heatmap(), "acf_and_m6a_heatmaps", 19, 11 + .55 * heat_result$n_clusters)
} else {
  empty <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
    label = "No resolved phased heterozygous reads for the heatmap") + ggplot2::theme_void()
  save_gg(empty, "heatmap_m6a_footprints", 10, 4)
  save_gg(empty, "acf_and_m6a_heatmaps", 10, 4)
}
save_gg(plot_knn_graph(result, main = paste0(region$annotation, "\nACF correlation-neighbor graph")),
  "knn_graph", 12, 9)
occupancy <- plot_signal_profile(signals$profiles, region, levels(a$cluster)) +
  ggplot2::labs(subtitle = "m6A calls and footprint occupancy; all clustered phased heterozygous reads")
save_gg(occupancy, "fig2_occupancy_by_cluster", 10, 1.4 * result$n_clusters + 2.5)
embedding <- data.frame(RID = a$RID, UMAP1 = a$umap1, UMAP2 = a$umap2)
save_fiberseq_plots(result, footprints, plot_dir, full_colors, embedding = embedding, plot_writer = save_gg)
message("Annotated ", nrow(a), " reads; known focal alleles: ", sum(a$allele_status == "phased_focal_genotype"),
  "; heatmap phased heterozygotes: ", sum(heat_keep), "; graph edges: ", igraph::ecount(graph))

# The shared helper creates an unused plots/ directory even with our flat writer.
# Remove only that empty directory; never delete unexpected contents.
unused_plot_dir <- file.path(plot_dir, "plots")
if (dir.exists(unused_plot_dir) &&
    !length(list.files(unused_plot_dir, all.files = TRUE, no.. = TRUE)))
  unlink(unused_plot_dir, recursive = TRUE)

cat("\nAUTOCOR_HEATMAP_READS=", sum(a$heatmap_included), "\n", sep = "")
