#!/usr/bin/env Rscript
# Molecule-state entropy from original promoter cluster proportions; no resampling.
suppressPackageStartupMessages({library(data.table)})
PROJECT <- '/project/spott/cshan/fiber-seq'
BASE <- file.path(PROJECT, 'LCL_project/pausing/PolII')
OUT <- file.path(PROJECT, 'LCL_project/entropy/molecule_cluster_entropy')

entropy_from_counts <- function(n) {
  stopifnot(all(is.finite(n)), all(n >= 0))
  n <- n[n > 0]; k <- length(n)
  if (!k) return(list(H = NA_real_, normalized = NA_real_, effective = NA_real_))
  p <- n / sum(n); stopifnot(abs(sum(p) - 1) < 1e-12)
  h <- -sum(p * log2(p))
  list(H = h, normalized = if (k > 1) h / log2(k) else NA_real_, effective = 2^h)
}

canonical_reads <- function(d) {
  d <- copy(d)
  if (!'sample_name' %in% names(d)) d[, sample_name := NA_character_]
  d[, RID := as.character(RID)]
  qualified <- !is.na(d$RID) & grepl('::', d$RID, fixed = TRUE)
  prefix <- ifelse(qualified, sub('::.*$', '', d$RID), NA_character_)
  missing_sample <- is.na(d$sample_name) | d$sample_name == ''
  d[missing_sample, sample_name := prefix[missing_sample]]
  d[, invalid_identity := is.na(RID) | RID == '' | is.na(sample_name) | sample_name == '' |
      (qualified & !is.na(prefix) & prefix != sample_name)]
  d[, fiber_id := ifelse(qualified, RID, paste(sample_name, RID, sep = '::'))]
  d
}

analyze_assignments <- function(a, cohort) {
  a <- canonical_reads(a); cohort <- canonical_reads(cohort)
  a[, cluster := as.character(cluster)]
  a[, valid_cluster := !is.na(cluster) & !tolower(trimws(cluster)) %in%
      c('', 'na', 'nan', 'none', 'unassigned', 'unclustered', '-1')]
  good <- a[invalid_identity == FALSE]
  assignments <- good[, {
    labels <- unique(cluster[valid_cluster])
    list(cluster = if(length(labels) == 1L) labels else NA_character_,
         assignment_rows = .N, invalid_label_rows = sum(!valid_cluster),
         ambiguous = length(labels) > 1L)
  }, by = .(fiber_id, sample_name)]
  universe <- unique(rbindlist(list(good[, .(fiber_id, sample_name)],
                                  cohort[invalid_identity == FALSE, .(fiber_id, sample_name)])))
  audit <- merge(universe, assignments, by = c('fiber_id', 'sample_name'), all.x = TRUE)
  audit[, assignment_status := fifelse(is.na(assignment_rows), 'missing_assignment',
    fifelse(ambiguous, 'ambiguous', fifelse(is.na(cluster), 'unassigned', 'assigned')))]
  counts <- audit[assignment_status == 'assigned', .(n_reads = .N), by = cluster]
  counts[, proportion := n_reads / sum(n_reads)]
  e <- entropy_from_counts(counts$n_reads)
  summary <- data.table(n_total_unique_reads = nrow(audit), n_assigned_reads = sum(counts$n_reads),
    n_clusters = nrow(counts), n_samples = uniqueN(audit$sample_name),
    shannon_entropy_bits = e$H, normalized_entropy = e$normalized, effective_states = e$effective)
  list(summary = summary, counts = counts, audit = audit)
}

main <- function() {
  setDTthreads(1L)
  td <- file.path(OUT, 'tables')
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  genes <- fread(file.path(PROJECT, 'code/pausing/PolII/example_gene_list.tsv'))
  stopifnot(!anyDuplicated(genes$gene_name))
  cohort <- fread(file.path(BASE, 'read_features/tables/read_features.tsv.gz'))
  results <- lapply(seq_len(nrow(genes)), function(i) {
    g <- genes[i]
    path <- file.path(BASE, 'recurrent_patterns/tables', g$group, paste0(g$gene_name, '_read_clusters.tsv'))
    exists <- file.exists(path)
    a <- if(exists) fread(path, colClasses = c(RID='character', cluster='character')) else
      data.table(RID=character(), sample_name=character(), cluster=character())
    stopifnot(all(c('RID','cluster') %in% names(a)))
    r <- analyze_assignments(a, cohort[gene == g$gene_name, .(RID, sample_name)])
    for (nm in names(r)) r[[nm]][, `:=`(gene_name = g$gene_name, gene_id = g$gene_id, group = g$group)]
    r$summary[, `:=`(assignment_file = path,
      status = if(n_assigned_reads > 0) 'valid' else 'no_valid_assignments')]
    r
  })
  write <- function(x,n) fwrite(x,file.path(td,n),sep='\t',na='NA')
  s <- rbindlist(lapply(results, `[[`, 'summary'), fill=TRUE)
  counts <- rbindlist(lapply(results, `[[`, 'counts'), fill=TRUE)
  stopifnot(all(abs(counts[, sum(proportion), by=gene_name]$V1 - 1) < 1e-12))
  write(s,'gene_molecule_cluster_entropy.tsv'); write(counts,'cluster_proportions.tsv')
  valid <- s[status == 'valid']
  med <- function(x) if(all(is.na(x))) NA_real_ else median(x,na.rm=TRUE)
  gs <- valid[, .(n_valid_genes=.N, median_shannon_entropy=median(shannon_entropy_bits),
    q25_shannon_entropy=quantile(shannon_entropy_bits,.25), q75_shannon_entropy=quantile(shannon_entropy_bits,.75),
    iqr_shannon_entropy=IQR(shannon_entropy_bits), median_normalized_entropy=med(normalized_entropy),
    median_effective_states=median(effective_states), median_cluster_number=as.numeric(median(n_clusters)),
    median_assigned_reads=as.numeric(median(n_assigned_reads))), by=group]
  write(gs,'gene_group_summary.tsv')
  tests <- rbindlist(lapply(c('shannon_entropy_bits','normalized_entropy','effective_states'),function(metric) {
    d <- valid[is.finite(get(metric))]; d[, value := get(metric)]
    kw <- kruskal.test(value ~ group, data=d)
    pairs <- combn(sort(unique(d$group)),2,simplify=FALSE)
    pair <- rbindlist(lapply(pairs,function(z) data.table(metric=metric,test='Wilcoxon rank-sum',
      comparison=paste(z,collapse=' vs '),p_value=wilcox.test(d[group==z[1],value],d[group==z[2],value],exact=FALSE)$p.value)))
    pair[, p_adjusted_BH := p.adjust(p_value,'BH')]
    rbind(data.table(metric=metric,test='Kruskal-Wallis',comparison='all groups',p_value=kw$p.value,p_adjusted_BH=NA_real_),pair)
  }))
  write(tests,'gene_group_tests.tsv')
  cat('Analyzed',nrow(s),'genes;',nrow(valid),'valid;',sum(valid$n_assigned_reads),'assigned gene/fiber observations.\n')
  print(gs)
}
if(sys.nframe()==0L) main()
