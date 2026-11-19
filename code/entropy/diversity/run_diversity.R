#!/usr/bin/env Rscript
# Run with the existing project R environment:
# bash code/pausing/PolII/run_gene_examples.sh code/entropy/diversity/run_diversity.R
suppressPackageStartupMessages({library(data.table); library(Matrix); library(ggplot2)})
main <- function() {
  setDTthreads(1L)
  project <- '/project/spott/cshan/fiber-seq'
  code <- file.path(project, 'code/entropy/diversity')
  out <- file.path(project, 'LCL_project/entropy/diversity')
  source(file.path(code, 'diversity_functions.R'))
  source(file.path(code, 'plotting_functions.R'))
  # GO-only mode reuses the completed diversity table; it never rebuilds matrices.
  # DIVERSITY_MODE=go bash code/pausing/PolII/run_gene_examples.sh code/entropy/diversity/run_diversity.R
  if (identical(Sys.getenv('DIVERSITY_MODE'), 'go')) {
    go <- rao_go_enrichment(fread(file.path(out, 'tables/gene_diversity.tsv')))
    go_out <- file.path(out, 'go_enrichment')
    for (subdir in c('tables', 'plots')) dir.create(file.path(go_out, subdir), recursive = TRUE, showWarnings = FALSE)
    for (name in names(go)) fwrite(go[[name]], file.path(go_out, 'tables', paste0('go_', name, '.tsv')),
      sep = '\t', na = 'NA')
    plot_rao_go(go$results, file.path(go_out, 'plots/go_biological_process.pdf'))
    print(go$summary)
    print(head(go$results[, c('rao_category','go_id','description','n_hits','fold_enrichment','p_adjusted_BH')], 15))
    message('GO enrichment complete: ', go_out)
    return(invisible(go))
  }
  source(file.path(project, 'code/parsing_functions/parsing_footprints_functions.r'))
  source(file.path(project, 'code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r'))
  validate_rao()
  genes <- as.data.frame(fread(file.path(project, 'code/pausing/PolII/example_gene_list.tsv')))
  groups <- c('housekeeping', 'cell_type_specific', 'not_expressed')
  stopifnot(nrow(genes) == 150L, !anyDuplicated(genes$gene_id),
    !anyDuplicated(genes$gene_name), all(genes$group %in% groups))
  samples <- as.data.frame(fread('/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv'))
  samples <- samples[, c('sample_name', 'fire_dir')]
  upstream <- as.data.frame(fread(file.path(project, 'LCL_project/pausing/PolII/read_features/tables/gene_features.tsv')))
  analyze_gene <- function(i) {
    gene <- genes[i, ]
    prefix <- file.path(project, 'LCL_project/pausing/PolII/recurrent_patterns/tables', gene$group, gene$gene_name)
    if (!file.exists(paste0(prefix, '_clustering.rds'))) {
      old <- upstream[upstream$gene_id == gene$gene_id, , drop = FALSE]
      # Preserve known upstream failures as missing, never infer zero fibers or clusters.
      stopifnot(nrow(old) == 1L, old$status == 'error',
        grepl('Insufficient full-span reads/m6A sites', old$error, fixed = TRUE),
        !file.exists(paste0(prefix, '_read_clusters.tsv')))
      message(gene$gene_name, ': NA; upstream failure: ', old$error)
      return(data.frame(gene_id = gene$gene_id, gene_name = gene$gene_name,
        gene_group = gene$group, n_fibers = NA_integer_, n_clusters = NA_integer_,
        rao_total = NA_real_, rao_between = NA_real_,
        between_cluster_share = NA_real_, q_over_g = NA_real_))
    }
    saved <- readRDS(paste0(prefix, '_clustering.rds'))
    region <- saved$region
    a <- as.data.frame(saved$assignments)
    tab <- as.data.frame(fread(paste0(prefix, '_read_clusters.tsv')))
    stopifnot(saved$params$window_size == 0L, region$chr == gene$chrom,
      region$analysis_start == region$tss - 1000L,
      region$analysis_end == region$tss + 1000L,
      saved$params$region_start == region$analysis_start,
      saved$params$region_end == region$analysis_end,
      !anyDuplicated(a$RID), !anyDuplicated(tab$RID), setequal(a$RID, tab$RID),
      identical(as.character(a$cluster), as.character(tab$cluster[match(a$RID, tab$RID)])))
    dat <- assemble_region_m6a(samples, region, full_span = TRUE,
      positions = seq.int(region$analysis_start, region$analysis_end),
      matrix_dir = NULL, reuse = FALSE)
    # Require the complete eligible cohort, rather than silently intersecting IDs.
    stopifnot(setequal(rownames(dat$met_mat), a$RID),
      nrow(dat$met_mat) == nrow(a), nrow(a) == saved$params$n_reads,
      all(dat$rids_df$start <= region$analysis_start),
      all(dat$rids_df$end >= region$analysis_end))
    mat <- bin_read_matrix(dat$met_mat, region$analysis_start, region$analysis_end, window_size = 0L)$mat
    stopifnot(ncol(mat) == 2001L, saved$params$n_features == ncol(mat))
    r <- rao_diversity(mat, a$cluster[match(rownames(mat), a$RID)])
    row <- cbind(data.frame(gene_id = gene$gene_id, gene_name = gene$gene_name,
      gene_group = gene$group), r)
    message(sprintf('[%d/150] %s: %d fibers, %d clusters, total=%.5f',
      i, gene$gene_name, r$n_fibers, r$n_clusters, r$rao_total))
    row
  }
  # Independent genes can be read concurrently; no sampling or random analysis steps.
  workers <- as.integer(Sys.getenv('DIVERSITY_WORKERS', '4'))
  stopifnot(length(workers) == 1L, !is.na(workers), workers >= 1L)
  rows <- parallel::mclapply(seq_len(nrow(genes)), analyze_gene,
    mc.cores = workers, mc.set.seed = FALSE)
  if (any(vapply(rows, inherits, logical(1), 'try-error'))) stop('Gene analysis failed; no outputs written')
  d <- do.call(rbind, rows)
  metrics <- c('rao_total', 'rao_between', 'between_cluster_share', 'q_over_g', 'n_clusters')
  tests <- compare_diversity_groups(d, metrics, groups)
  summary <- summarize_diversity(d, metrics, groups, tests$omnibus)
  for (subdir in c('tables', 'plots')) dir.create(file.path(out, subdir), recursive = TRUE, showWarnings = FALSE)
  fwrite(d, file.path(out, 'tables/gene_diversity.tsv'), sep = '\t', na = 'NA')
  fwrite(summary, file.path(out, 'tables/group_diversity_summary.tsv'), sep = '\t', na = 'NA')
  plot_diversity(d, file.path(out, 'plots'), groups)
  print(summary, row.names = FALSE)
  cat('\nExploratory gene-level tests: BH across five omnibus tests and across 15 pairwise tests.\n')
  print(tests$pairwise, row.names = FALSE)
  cat('\nExact cohort and decomposition checks passed for', sum(is.finite(d$rao_total)),
    'genes;', sum(!is.finite(d$rao_total)), 'upstream failures retained as NA.\n')
  invisible(d)
}
if (sys.nframe() == 0L) main()
