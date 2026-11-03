#!/usr/bin/env Rscript
# Native R orchestration: phased reads -> ACF -> Leiden -> report_data.rds.
acf_default_project <- "/project/spott/cshan/fiber-seq"

acf_load_functions <- function(project) {
  directory <- file.path(project, "code/clustering_methods/auto_correlation")
  env <- new.env(parent = globalenv())
  for (file in c("01_prepare_m6a_input.R", "02_compute_autocorrelations.R",
                 "03_cluster_autocorrelations.R", "04_summarize_autocorrelations.R"))
    sys.source(file.path(directory, file), envir = env)
  # Only load named function chunks; never execute setup, examples or plots.
  lines <- readLines(file.path(directory, "06_autocorrelation.Rmd"), warn = FALSE)
  starts <- grep("^```\\{r[[:space:]]+", lines)
  for (start in starts) {
    name <- sub("^```\\{r[[:space:]]+([^,}[:space:]]+).*", "\\1", lines[start])
    if (!startsWith(name, "function-") && !endsWith(name, "-function")) next
    end <- which(seq_along(lines) > start & lines == "```")[1L]
    if (is.na(end)) stop("Unclosed report chunk: ", name)
    eval(parse(text = lines[seq.int(start + 1L, end - 1L)]), envir = env)
  }
  stopifnot(is.function(env$prepare_report_data), is.function(env$autocor_payload_matrix))
  env
}

acf_save_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(".acf-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(value, temporary, compress = "gzip")
  if (!file.rename(temporary, path)) stop("Could not save ", path)
  invisible(path)
}

acf_selection <- function(project, output, requested, functions) {
  if (!is.null(requested)) {
    if (!file.exists(requested)) stop("Selection signature does not exist: ", requested)
    return(normalizePath(requested))
  }
  candidates <- c(file.path(output, "inputs/selection_signature.rds"),
    file.path(project, "LCL_project/Leiden_manhattan/summary tables/selection_signature.rds"),
    file.path(project, "LCL_project/auto_correlation/resolution04/inputs/selection_signature.rds"),
    file.path(project, "LCL_project/auto_correlation/resolution1/inputs/selection_signature.rds"))
  found <- candidates[file.exists(candidates)]
  if (length(found)) return(normalizePath(found[1L]))
  functions$recover_m6a_selection(file.path(project, "LCL_project/Leiden_manhattan"), candidates[1L])
  normalizePath(candidates[1L])
}

acf_input_signature <- function(region, selection, project, ft_root) {
  samples <- strsplit(region$contributing_samples, ",", fixed = TRUE)[[1L]]
  phase_root <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples"
  files <- unlist(lapply(samples, function(sample) {
    bed <- file.path(ft_root, sample, "extracted_results/m6a_by_chr",
      paste0(sample, ".ft_extracted_m6a.", region$chr, ".bed.gz"))
    vcf <- file.path(phase_root, sample, paste0(sample, ".5mC.6mA.aligned.phased.vcf.gz"))
    c(bed, paste0(bed, ".tbi"), vcf, paste0(vcf, ".tbi"),
      file.path(phase_root, sample, c("read-level-phasing.tsv", "blocks.tsv")),
      file.path(project,"LCL_project/Leiden_manhattan/phase_summary",paste0(sample,"_haplotags.rds")))
  }), use.names = FALSE)
  info <- file.info(files)
  sources <- file.path(project, c("code/clustering_methods/auto_correlation/01_prepare_m6a_input.R",
    "code/parsing_functions/parsing_footprints_functions.r",
    "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r",
    "code/haplotype_phasing/LCL_phasing.r", "code/haplotype_phasing/LCL_phased_m6a_input.r"))
  list(schema_version=1L, region=region, selection_md5=unname(tools::md5sum(selection)),
    ft_root=normalizePath(ft_root, mustWork=TRUE), source_md5=tools::md5sum(sources),
    files=data.frame(path=files, size=info$size, mtime=as.numeric(info$mtime), ctime=as.numeric(info$ctime)))
}

acf_validate_input <- function(input, region) {
  stopifnot(is.matrix(input$binary), is.data.frame(input$records),
    nrow(input$binary) == nrow(input$records), nrow(input$records) > 0L,
    ncol(input$binary) == region$width, all(is.finite(input$binary)),
    all(input$binary %in% c(0,1)), !anyNA(input$records$row_id),
    !anyDuplicated(input$records$row_id),
    identical(rownames(input$binary), input$records$row_id),
    all(input$records$region_id == region$region_id),
    all(input$records$haplotype %in% c("HP1","HP2")),
    all(input$records$focal_genotype %in% c("0|1","1|0")),
    all(input$records$allele_status == "phased_focal_genotype"))
  input$records$m6a_calls <- rowSums(input$binary)
  observed <- sum(colSums(input$binary) > 0)
  input$records$n_observed_sites <- observed
  input$records$m6a_call_fraction <- if (observed) input$records$m6a_calls/observed else NA_real_
  input
}

acf_prepare_input <- function(region, selection_path, project, ft_root, output, functions) {
  signature <- acf_input_signature(region, selection_path, project, ft_root)
  path <- file.path(output, "inputs", region$region_id, "m6a_input.rds")
  if (file.exists(path)) {
    saved <- readRDS(path)
    if (identical(saved$signature, signature)) {
      message(region$region_id, ": reusing R m6A input cache")
      return(acf_validate_input(saved, region))
    }
    message(region$region_id, ": input sources changed; rebuilding R input cache")
  }
  selection <- readRDS(selection_path)
  z <- functions$build_m6a_input(region, selection$sample_table, ft_root, project)
  records <- as.data.frame(z$read_meta)
  records$region_id <- region$region_id
  input <- acf_validate_input(list(binary=as.matrix(z$met_mat), records=records,
    signature=signature, phasing_sources=z$phasing_sources), region)
  if (!identical(signature, acf_input_signature(region, selection_path, project, ft_root)))
    stop("Input files changed during extraction; rerun before caching")
  acf_save_rds(input, path)
  input
}

acf_analyze_region <- function(input, region, sample_table, params, functions) {
  input <- acf_validate_input(input, region)
  result <- functions$autocorrelations(input$binary, params$n_features)
  records <- functions$annotate_alleles(input$records, region, sample_table)
  records$acf_valid <- result$valid
  peaks <- t(vapply(seq_len(nrow(records)), function(i)
    functions$repeat_peak(result$profiles[i, ]), numeric(2L)))
  records$positive_local_peak_140_250_bp <- peaks[,1L]
  records$peak_acf <- peaks[,2L]
  records$peak_band_complete <- params$n_features >= 252L
  clustering <- functions$cluster_profiles(result$profiles, result$valid,
    n_pcs=params$n_pcs, n_neighbors=params$n_neighbors,
    resolution=params$resolution, seed=params$seed)
  records$cluster <- clustering$labels
  records$status <- clustering$status
  records$umap1 <- clustering$embedding[,1L]
  records$umap2 <- clustering$embedding[,2L]
  summary <- functions$summarize_autocorrelations(result$profiles, records)
  alleles <- functions$summarize_alleles(result$profiles, records, region)
  edges <- clustering$edges
  edges$source <- records$row_id[edges$source]
  edges$target <- records$row_id[edges$target]
  params$k_eff <- if (is.null(clustering$info$actual_neighbors)) 0L else clustering$info$actual_neighbors
  list(schema_version=1L, region=region, records=records, row_ids=records$row_id,
    call_offsets=lapply(seq_len(nrow(input$binary)), function(i) which(input$binary[i,] == 1)-1L),
    acf=result$profiles, params=params,
    tables=c(list(cluster_acf_summary=summary$averages, cluster_summary=summary$summaries,
      sample_composition=summary$composition), alleles), edges=edges,
    clustering_info=clustering$info)
}

run_lcl_autocorrelation <- function(project=acf_default_project, output_dir=NULL,
    regions=NULL, window_size=2000L, n_features=window_size, n_pcs=50L,
    n_neighbors=10L, resolution=0.4, seed=0L, selection_signature=NULL,
    sample_table="/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv",
    ft_root="/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results",
    tracks_root="/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged/combined/joint_trained_tracks") {
  project <- normalizePath(project, mustWork=TRUE)
  if (is.null(output_dir)) output_dir <- file.path(project,"LCL_project/auto_correlation/resolution04")
  integers <- c(window_size,n_features,n_pcs,n_neighbors,seed)
  if (any(!is.finite(integers)) || any(integers != floor(integers)) || window_size < 3L ||
      n_features < 3L || n_features > window_size || n_pcs < 2L || n_neighbors < 2L ||
      seed < 0L || seed > .Machine$integer.max || !is.finite(resolution) || resolution <= 0)
    stop("Require integer window-size >=3; 3 <= n-features <= window-size; PCs/neighbors >=2; nonnegative seed; positive resolution")
  functions <- acf_load_functions(project)
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  output <- normalizePath(output_dir, mustWork=TRUE)
  selection <- acf_selection(project,output,selection_signature,functions)
  selected <- functions$load_m6a_regions(readRDS(selection), window_size=window_size)
  if (!is.null(regions) && length(setdiff(regions,selected$region_id)))
    stop("Unknown region IDs: ",paste(setdiff(regions,selected$region_id),collapse=", "))
  samples <- utils::read.csv(sample_table, colClasses="character",check.names=FALSE)
  stopifnot(all(c("sample_name","cell_line") %in% names(samples)), !anyDuplicated(samples$sample_name))
  params <- list(window_size=as.integer(window_size), n_features=as.integer(n_features),
    n_pcs=as.integer(n_pcs), n_neighbors=as.integer(n_neighbors),resolution=resolution,seed=as.integer(seed))
  directory <- file.path(project,"code/clustering_methods/auto_correlation")
  sources <- c(list.files(directory, pattern="\\.(R|Rmd)$", full.names=TRUE),
    file.path(project,c("code/parsing_functions/parsing_footprints_functions.r",
      "code/parsing_functions/plotting_functions.r",
      "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r",
      "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r",
      "code/haplotype_phasing/LCL_phasing.r", "code/haplotype_phasing/LCL_phased_m6a_input.r")))
  manifest <- list(schema_version=1L, implementation="native_R_acf_v1", params=params,
    region_ids=selected$region_id, selection_md5=unname(tools::md5sum(selection)),
    sample_table_md5=unname(tools::md5sum(sample_table)), source_md5=tools::md5sum(sources),
    ft_root=normalizePath(ft_root,mustWork=TRUE), tracks_root=normalizePath(tracks_root,mustWork=TRUE),
    R_version=as.character(getRversion()),
    packages=vapply(c("igraph","uwot","irlba","Matrix"),function(x) as.character(utils::packageVersion(x)),character(1)))
  path <- file.path(output,"inputs/report_run.rds")
  if (file.exists(path) && !identical(readRDS(path),manifest))
    message("Run settings or code changed; incompatible region outputs will be recomputed")
  acf_save_rds(manifest,path)
  selection_copy <- file.path(output,"inputs/selection_signature.rds")
  if (!file.exists(selection_copy) ||
      !identical(unname(tools::md5sum(selection)),unname(tools::md5sum(selection_copy)))) {
    if (!file.copy(selection,selection_copy,overwrite=TRUE)) stop("Could not copy selection signature")
  }
  if (!is.null(regions)) selected <- selected[selected$region_id %in% regions,,drop=FALSE]
  destinations <- character(nrow(selected))
  for (i in seq_len(nrow(selected))) {
    region <- selected[i,,drop=FALSE]
    destination <- file.path(output,"outputs",region$region_id,"report_data.rds")
    destinations[i] <- destination
    input <- acf_prepare_input(region,selection,project,ft_root,output,functions)
    if (file.exists(destination)) {
      saved <- readRDS(destination)
      if (identical(saved$run_manifest,manifest) && identical(saved$input_signature,input$signature)) {
        message(region$region_id,": report_data.rds already complete; skipping")
        next
      }
      message(region$region_id,": replacing previous report inputs with native R results")
    }
    message(region$region_id,": computing R ACF and clustering ",nrow(input$records)," phased reads")
    bundle <- acf_analyze_region(input,region,samples,params,functions)
    bundle$run_manifest <- manifest
    bundle$input_signature <- input$signature
    if (!identical(manifest$source_md5,tools::md5sum(sources))) stop("Code changed during run")
    functions$prepare_report_data(bundle,project,destination,ft_root=ft_root,tracks_root=tracks_root)
    message(region$region_id,": saved ",bundle$clustering_info$n_clusters," clusters; ",destination)
  }
  if (!identical(manifest$source_md5,tools::md5sum(sources))) stop("Code changed during run")
  message("Complete. Knit 06_autocorrelation.Rmd with output_dir=",output," to draw the figures.")
  invisible(setNames(destinations,selected$region_id))
}

acf_main <- function(args=commandArgs(trailingOnly=TRUE)) {
  if (any(args %in% c("--help","-h"))) {
    cat("Usage: Rscript 05_run_autocorrelation.R [--region ID (repeatable)] [--output-dir DIR]\n",
      "  [--window-size 2000] [--n-features N] [--n-pcs 50] [--n-neighbors 10]\n",
      "  [--resolution 0.4] [--seed 0] [--selection-signature FILE] [--sample-table CSV]\n",
      "  [--project DIR] [--ft-root DIR] [--tracks-root DIR]\n",
      "Default output: LCL_project/auto_correlation/resolution04; cached matrices are RDS.\n",sep="")
    return(invisible(NULL))
  }
  if (length(args) %% 2L) stop("Options require a value; use --help")
  values <- list(); regions <- character()
  numeric <- c("window_size","n_features","n_pcs","n_neighbors","resolution","seed")
  allowed <- c(numeric,"project","output_dir","selection_signature","sample_table","ft_root","tracks_root","region")
  for (i in seq_along(args)[seq_along(args) %% 2L == 1L]) {
    if (!startsWith(args[i],"--")) stop("Expected an option: ",args[i])
    name <- gsub("-","_",substring(args[i],3L),fixed=TRUE)
    if (!name %in% allowed) stop("Unknown option: ",args[i])
    value <- args[i+1L]
    if (name == "region") { regions <- c(regions,value); next }
    if (name %in% names(values)) stop("Repeated option: ",args[i])
    values[[name]] <- if (name %in% numeric) as.numeric(value) else value
  }
  if (length(regions)) values$regions <- unique(regions)
  do.call(run_lcl_autocorrelation,values)
}
if (sys.nframe() == 0L) acf_main()
