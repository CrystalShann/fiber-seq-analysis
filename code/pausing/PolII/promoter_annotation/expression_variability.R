#!/usr/bin/env Rscript
# Run: bash code/pausing/PolII/run_gene_examples.sh code/pausing/PolII/promoter_annotation/expression_variability.R
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
BASE <- '/project/spott/cshan/fiber-seq/LCL_project/pausing/PolII'
SOURCE_DIR <- '/project/spott/cshan/annotations/LCL_public_data/Einarsson_2022_eLife_80943'
SOURCE <- file.path(SOURCE_DIR, 'supplementary_file_2.tsv')
OUT <- file.path(BASE, 'promoter_annotation')
TABLES <- file.path(OUT, 'tables')
PLOTS <- file.path(OUT, 'plots')
FAMILIES <- c('ETS', 'SP', 'NFY', 'CREB', 'NRF1', 'ZNF143')

# TRUE where point lies in the closed interval [start, end] (vectorised over start/end)
contains <- function(start, end, point) start <= point & point <= end

# Figures for the genes matched to a published Einarsson promoter:
# 1. feature_boxplots.pdf: TATA activation, nucleosome occupancy and m6A fraction (TSS +/-200 bp)
#    per published adjusted-variability quartile, box + one jittered dot per gene
# 2. expression_m6a_nucleosome_variability.pdf: m6A fraction and nucleosome occupancy vs the
#    published adjusted log10(CV2), one dot per gene
comparison_figures <- function(d) {
  blue <- '#2166AC'
  d <- d[match_status == 'matched']

  qlevels <- c('0-25%','25-50%','50-75%','75-100%')
  stopifnot(all(d$elife_adjusted_quartile %in% qlevels))
  metrics <- c('tata_activation','nuc_core','m6a_core')
  labels <- c('TATA activation','Nucleosome occupancy (TSS ±200 bp)','m6A fraction (TSS ±200 bp)')
  pdf(file.path(PLOTS,'feature_boxplots.pdf'),width=12,height=4.5,useDingbats=FALSE)
  par(mfrow=c(1,3),mar=c(5,9,3,1),oma=c(3,0,2,0))
  set.seed(80943)
  for(k in seq_along(metrics)) {
    vals <- d[[metrics[k]]]
    groups <- rev(qlevels)
    samples <- lapply(groups,function(q) vals[d$elife_adjusted_quartile==q & is.finite(vals)])
    names(samples) <- paste0(groups,' (n=',lengths(samples),')')
    boxplot(samples,horizontal=TRUE,las=1,outline=FALSE,ylim=range(vals[is.finite(vals)]),col='#F0F2F2',border='#444444',
      xlab=labels[k],main=labels[k],cex.axis=.75,cex.main=.85,boxwex=.5)
    title(ylab='Adjusted log10 squared coefficient of variation',line=6.5,cex.lab=.85)
    for(j in seq_along(samples)) {
      points(samples[[j]],j+runif(length(samples[[j]]),-.14,.14),pch=16,
        col=adjustcolor(blue,.4),cex=.65)
      if(length(samples[[j]])) segments(median(samples[[j]]),j-.25,median(samples[[j]]),j+.25,lwd=2)
    }
  }
  mtext('Promoter features across published variability quartiles',outer=TRUE,side=3,cex=1.15,font=2)
  mtext('Expression-adjusted variability metric from Einarsson et al., 2022',outer=TRUE,side=1,line=1,cex=.8)
  dev.off()

  pdf(file.path(PLOTS,'expression_m6a_nucleosome_variability.pdf'),width=13,height=6,useDingbats=FALSE)
  par(mfrow=c(1,2),mar=c(5,5,3,1),oma=c(2,0,2,0))
  for(key in c('m6a_core','nuc_core')) {
    z <- d[is.finite(get(key)) & is.finite(elife_adjusted_log10_CV2)]
    plot(z[[key]],z$elife_adjusted_log10_CV2,pch=21,bg=adjustcolor(blue,.5),col='#444444',cex=1.1,
      xlab=if(key=='m6a_core')'m6A fraction (TSS ±200 bp)' else 'Nucleosome occupancy (TSS ±200 bp)',
      ylab='Adjusted log10(CV²)',
      main=paste(if(key=='m6a_core')'m6A fraction' else 'Nucleosome occupancy',sprintf('(n=%d)',nrow(z))))
    abline(h=0,lty=3,col='gray60')
  }
  mtext('Promoter chromatin vs expression variability',outer=TRUE,side=3,font=2)
  mtext('Expression-adjusted variability metric from Einarsson et al., 2022',outer=TRUE,side=1,cex=.8)
  dev.off()
}

# 1. read the fiber gene features and the Puffin 04a scores, check they describe the same genes
#    at the same TSS (fiber tss = puffin_tss0 + 1)
# 2. per gene: keep the fiber features and Puffin contributions, summing + and - orientations
#    into one contribution per motif family
# 3. read the published table (Einarsson et al., 2022 supplementary file 2) and parse its promoter
#    IDs (chr:start-end;strand)
# 4. match each gene to the published promoter on the same chromosome and strand whose interval
#    holds the fiber TSS and whose Ensembl gene ID agrees; exactly one hit = matched (variability
#    columns copied as elife_*), more = ambiguous, none = unmatched
# 5. write gene_annotations.tsv (one row per gene) and gene_variability_qc.tsv (the candidate
#    promoters behind every match), then the figures unless ANNOTATION_TABLES_ONLY=1
main <- function() {
  stopifnot(file.exists(SOURCE))
  dir.create(TABLES, recursive = TRUE, showWarnings = FALSE)
  dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)
  gf_path <- file.path(BASE, 'read_features/tables/gene_features.tsv')
  pf_path <- file.path(TABLES, '04a_puffin_motif_scores.tsv')
  gf <- fread(gf_path)
  pf <- fread(pf_path, check.names = FALSE)
  stopifnot(!anyDuplicated(pf$gene_name), all(gf$gene_name %in% pf$gene_name))
  pf <- pf[match(gf$gene_name, gene_name)]
  stopifnot(nrow(gf) == nrow(pf), !anyDuplicated(gf$gene_name),
    !anyDuplicated(pf$gene_name), identical(gf$gene_name, pf$gene_name),
    all(gf$chrom == pf$chrom), all(gf$strand == pf$strand),
    all(gf$tss == pf$puffin_tss0 + 1), contains(10,20,10), contains(10,20,20),
    !contains(10,20,9), !contains(10,20,21))
  d <- gf[, .(gene_name, gene_id, chrom, strand, tss, fiber_status = status,
    fiber_error = error, n_reads, n_samples, nuc_core, m6a_core)]
  d[, puffin_tss0 := pf$puffin_tss0]
  orientations <- c('TATA +', 'YY1 +', as.vector(t(outer(FAMILIES, c('+', '-'), paste))))
  for (col in c(orientations, 'tata_activation')) set(d, j = col, value = pf[[col]])
  d[, `:=`(TATA_contribution = pf[['TATA +']], YY1_contribution = pf[['YY1 +']])]
  for (f in FAMILIES) {
    value <- pf[[paste(f, '+')]] + pf[[paste(f, '-')]]
    set(d, j = paste0(f, '_contribution'), value = value)
    stopifnot(isTRUE(all.equal(value, d[[paste(f, '+')]] + d[[paste(f, '-')]])))
  }
  features <- c(paste0(c('TATA','ETS','YY1','SP','NFY','CREB','NRF1','ZNF143'), '_contribution'), 'tata_activation')
  stopifnot(all(is.finite(as.matrix(d[, features, with = FALSE]))),
    all(d$nuc_core >= 0 & d$nuc_core <= 1, na.rm = TRUE),
    all(d$m6a_core >= 0 & d$m6a_core <= 1, na.rm = TRUE))
  # R write.table's unnamed row-name field is absent from the TSV header.
  pub <- as.data.table(read.delim(SOURCE, row.names = 1, check.names = FALSE), keep.rownames = 'published_promoter')
  stopifnot(!anyDuplicated(pub$published_promoter),
    all(grepl('^chr[^:]+:[0-9]+-[0-9]+;[+-]$', pub$published_promoter)))
  parts <- tstrsplit(pub$published_promoter, '[:;]')
  interval <- tstrsplit(parts[[2]], '-', fixed = TRUE)
  pub[, `:=`(chrom = parts[[1]], start1 = as.integer(interval[[1]]),
    end1 = as.integer(interval[[2]]), strand = parts[[3]], stable_gene_id = sub('\\.[0-9]+$', '', geneID))]
  cols <- c('geneID','mean','median','CV','log10_CV2','adjusted_log10_CV2','adjusted_quartile')
  for (col in cols) set(d, j = paste0('elife_', col), value = if (is.character(pub[[col]])) NA_character_ else NA_real_)
  d[, `:=`(published_promoter = NA_character_, published_start1 = NA_integer_, published_end1 = NA_integer_)]
  audit <- vector('list', nrow(d))
  for (i in seq_len(nrow(d))) {
    g <- d[i]
    same <- pub[chrom == g$chrom & strand == g$strand]
    positional <- same[contains(start1, end1, g$tss)]
    candidates <- positional[stable_gene_id == sub('\\.[0-9]+$', '', g$gene_id)]
    alternative <- same[start1 < g$tss & g$tss <= end1 & stable_gene_id == sub('\\.[0-9]+$', '', g$gene_id)]
    status <- if (nrow(candidates) == 1) 'matched' else if (nrow(candidates) > 1) 'ambiguous' else 'unmatched'
    d[i, `:=`(match_status = status, candidate_count = nrow(candidates))]
    audit[[i]] <- data.table(gene_name = g$gene_name, tss1 = g$tss,
      positional_candidates = paste(positional$published_promoter, collapse = '|'),
      identity_matched_candidates = paste(candidates$published_promoter, collapse = '|'),
      coordinate_sensitivity_same = setequal(candidates$published_promoter, alternative$published_promoter))
    if (status == 'matched') {
      d[i, `:=`(published_promoter = candidates$published_promoter,
        published_start1 = candidates$start1, published_end1 = candidates$end1)]
      for (col in cols) set(d, i = i, j = paste0('elife_', col), value = candidates[[col]])
    }
  }
  fwrite(d, file.path(TABLES, 'gene_annotations.tsv'), sep = '\t', na = 'NA')
  fwrite(rbindlist(audit), file.path(TABLES, 'gene_variability_qc.tsv'), sep = '\t')
  if (Sys.getenv('ANNOTATION_TABLES_ONLY') != '1') comparison_figures(d)
}
main()
