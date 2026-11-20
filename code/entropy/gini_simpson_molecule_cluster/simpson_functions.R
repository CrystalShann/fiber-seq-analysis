#!/usr/bin/env Rscript
# Gini-Simpson diversity of original, gene-specific Fiber-seq clusters.
# Run in the project's R environment:
# bash code/pausing/PolII/run_gene_examples.sh code/entropy/gini_simpson_molecule_cluster/simpson_functions.R
# GINI_MODE=diversity-only validates assignments and exports diversity without feature extraction.

# 1. Packages and paths. All outputs are isolated from existing entropy results.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
PROJECT <- '/project/spott/cshan/fiber-seq'
CODE <- file.path(PROJECT,'code/entropy/gini_simpson_molecule_cluster')
OUT <- file.path(PROJECT,'LCL_project/entropy/gini_simpson_molecule_cluster')
POLII <- file.path(PROJECT,'LCL_project/pausing/PolII')
SHANNON <- file.path(PROJECT,'LCL_project/entropy/molecule_cluster_entropy/tables')
GROUPS <- c('housekeeping','cell_type_specific','not_expressed')

# Feature settings confirmed by the user.
# Fixed per the request: all assigned fibers, original cluster labels, no resampling.
FOOTPRINT_WINDOWS <- list(tss=c(-100L,100L),downstream=c(0L,300L))
FOOTPRINT_SIZES <- list(fp40_60=c(40L,60L),fp60_80=c(60L,80L),fp10_80=c(10L,80L))
FOOTPRINT_MIN_OVERLAP <- .5
# 60 bp belongs to fp60_80, not fp40_60; 80 bp is included. At least 50% of the footprint length must overlap.
# Existing continuous m6A features: upstream [-500,-201], core [-200,200], downstream [201,500].
# No new binary accessibility rule is defined.


# Calculation functions. No clustering or resampling.
simpson_indices <- function(counts) {
  stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts == floor(counts)))
  counts <- counts[counts > 0]; k <- length(counts)
  if (!k) return(list(gini_simpson=NA_real_, normalized_gini_simpson=NA_real_, inverse_simpson=NA_real_))
  p <- counts/sum(counts); concentration <- sum(p^2); g <- 1-concentration
  stopifnot(abs(sum(p)-1)<1e-12)
  list(gini_simpson=g, normalized_gini_simpson=if(k>1) g/(1-1/k) else NA_real_, inverse_simpson=1/concentration)
}

finite_mean <- function(x) if(any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
finite_median <- function(x) if(any(is.finite(x))) as.numeric(median(x[is.finite(x)])) else NA_real_

# Same FiberHMM v1 calls as existing entropy helpers, with new exact size ranges.
# Coordinates supplied here are 1-based inclusive relative to the oriented TSS.
footprint_presence <- function(reads, fp, tss, strand) {
  rr <- rel_interval(reads$start, reads$end, tss, strand)
  ff <- copy(fp)
  if(nrow(ff)) ff[, c('lo','hi') := rel_interval(start+1L,end,tss,strand)]
  else ff[, `:=`(lo=integer(),hi=integer())]
  ff[, size := end-start]
  windows <- FOOTPRINT_WINDOWS
  sizes <- FOOTPRINT_SIZES
  rbindlist(lapply(names(windows),function(w) {
    win <- windows[[w]]
    eligible <- reads$RID[rr$rel_start<=win[1] & rr$rel_end>=win[2]]
    rbindlist(lapply(names(sizes),function(s) {
      bounds <- sizes[[s]]
      keep <- ff$size>=bounds[1] & if(s=='fp40_60') ff$size<bounds[2] else ff$size<=bounds[2]
      overlap <- pmax(0L, pmin(ff$hi,win[2])-pmax(ff$lo,win[1])+1L)
      n <- uniqueN(ff[keep & RID %in% eligible & overlap >= FOOTPRINT_MIN_OVERLAP * size,RID])
      data.table(feature=paste(s,w,sep='_'),size_class=s,window=w,window_start=win[1],window_end=win[2],
        minimum_footprint_overlap_fraction=FOOTPRINT_MIN_OVERLAP,
        n_eligible_fibers=length(eligible),n_fibers_with_footprint=n,
        fraction=if(length(eligible)) n/length(eligible) else NA_real_)
    }))
  }))
}

feature_associations <- function(d, features) {
  r <- rbindlist(lapply(features,function(f) rbindlist(lapply(c('all',GROUPS),function(gr) {
    z <- d[is.finite(gini_simpson) & is.finite(get(f))]
    if(gr!='all') z <- z[gene_group==gr]
    rho <- if(nrow(z)>=3 && uniqueN(z$gini_simpson)>1 && uniqueN(z[[f]])>1)
      cor(z$gini_simpson,z[[f]],method='spearman') else NA_real_
    data.table(feature=f,stratum=gr,n_promoters=nrow(z),spearman_rho=rho)
  }))))
  # Correlation only; no p-values.
  r
}

group_comparisons <- function(d) {
  rbindlist(lapply(c('gini_simpson','normalized_gini_simpson'),function(f) {
    z <- d[is.finite(get(f))]; z[, value := get(f)]
    kw <- kruskal.test(value~gene_group,data=z)
    pw <- rbindlist(lapply(combn(GROUPS,2,simplify=FALSE),function(g) {
      data.table(metric=f,test='Wilcoxon rank-sum',comparison=paste(g,collapse=' vs '),
        p_value=suppressWarnings(wilcox.test(z[gene_group==g[1],value],z[gene_group==g[2],value],exact=FALSE)$p.value))
    }))
    pw[, p_adjusted_BH := p.adjust(p_value,'BH')]
    rbind(data.table(metric=f,test='Kruskal-Wallis',comparison='all groups',p_value=kw$p.value,p_adjusted_BH=NA_real_),pw)
  }))
}

save_analysis_table <- function(d,name) fwrite(d,file.path(OUT,'tables',name),sep='\t',na='NA')

# Lowest, median and highest Gini-Simpson promoters.
representative_promoters <- function(d) {
  sorted <- d[is.finite(gini_simpson)][order(gini_simpson,gene_name)]
  representatives <- sorted[c(1L,ceiling(.N/2),.N)]
  representatives[,representative:=c('Low','Intermediate','High')]
  representatives
}

main <- function() {
  dir.create(file.path(OUT,'tables'),recursive=TRUE,showWarnings=FALSE)
  mode <- Sys.getenv('GINI_MODE','full')
  stopifnot(mode %in% c('full','diversity-only'))
  # 2. Example gene list: exactly one existing promoter per gene in this dataset.
  genes <- fread(file.path(PROJECT,'code/pausing/PolII/example_gene_list.tsv'))
  stopifnot(!anyDuplicated(genes$gene_name),!anyDuplicated(genes$gene_id),all(genes$group %in% GROUPS))
  locations <- fread(file.path(POLII,'read_features/tables/gene_features.tsv'))
  stopifnot(!anyDuplicated(locations$gene_name))
  genes <- merge(genes[,.(gene_id,gene_name,gene_group=group)],
    locations[,.(gene_name,chrom,tss,strand,tss_source)],by='gene_name',all.x=TRUE)
  genes[, promoter_id:=ifelse(!is.na(tss),paste(gene_id,chrom,tss,strand,sep=':'),paste0(gene_id,':unresolved'))]
  stopifnot(!anyDuplicated(genes$promoter_id))
  # 3. The requested old paths have moved; use current tables without modifying them.
  proportions <- fread(file.path(SHANNON,'cluster_proportions.tsv'))
  summary <- fread(file.path(SHANNON,'gene_molecule_cluster_entropy.tsv'))
  saved_reads <- fread(file.path(POLII,'read_features/tables/read_features.tsv.gz'))
  stopifnot(!anyDuplicated(saved_reads[,.(gene,RID)]))
  original <- new.env(parent=globalenv())
  sys.source(file.path(PROJECT,'code/entropy/molecule_cluster/molecule_cluster_entropy.R'),envir=original)
  # 4. Validate existing proportions, gene annotations, and original per-read counts.
  proportions <- proportions[gene_name %in% genes$gene_name]
  stopifnot(!anyDuplicated(proportions[,.(gene_name,cluster)]),
    all(proportions$n_reads>0), all(proportions$n_reads==floor(proportions$n_reads)),
    all(proportions$gene_id==genes$gene_id[match(proportions$gene_name,genes$gene_name)]),
    all(proportions$group==genes$gene_group[match(proportions$gene_name,genes$gene_name)]))
  proportions[, recalculated:=n_reads/sum(n_reads),by=gene_name]
  stopifnot(all(abs(proportions$proportion-proportions$recalculated)<1e-12))
  proportions[,proportion:=recalculated];proportions[,recalculated:=NULL]
  assigned <- vector('list',nrow(genes)); rows <- vector('list',nrow(genes))
  for(i in seq_len(nrow(genes))) {
    g <- genes[i]
    path <- file.path(POLII,'recurrent_patterns/tables',g$gene_group,paste0(g$gene_name,'_read_clusters.tsv'))
    a <- if(file.exists(path)) fread(path,colClasses=c(RID='character',cluster='character')) else
      data.table(RID=character(),sample_name=character(),cluster=character())
    checked <- original$analyze_assignments(a,saved_reads[gene==g$gene_name,.(RID,sample_name)])
    old <- proportions[gene_name==g$gene_name][order(cluster)]
    new <- checked$counts[order(cluster)]
    stopifnot(identical(as.character(old$cluster),as.character(new$cluster)),
      identical(as.integer(old$n_reads),as.integer(new$n_reads)))
    assigned[[i]] <- checked$audit[assignment_status=='assigned',fiber_id]
    # 5/6. Raw, normalized and inverse Simpson; preserve original cluster labels.
    n <- sum(old$n_reads); k <- nrow(old)
    dominant <- if(k) old[order(-n_reads,cluster)][1] else NULL
    sh <- summary[gene_name==g$gene_name]
    # Genes failing the molecule-cluster read/cluster filter (status != valid) get no diversity
    indices <- if(sh$status=='valid') simpson_indices(old$n_reads) else simpson_indices(integer())
    stopifnot(nrow(sh)==1L,n==sh$n_assigned_reads,k==sh$n_clusters)
    if(k) stopifnot(abs(-sum(old$proportion*log2(old$proportion))-sh$shannon_entropy_bits)<1e-10)
    rows[[i]] <- cbind(g,data.table(n_fibers_assigned=n,n_clusters=k,
      dominant_cluster_id=if(k) dominant$cluster else NA_character_,
      dominant_cluster_fraction=if(k) dominant$proportion else NA_real_,
      n_dominant_ties=if(k) sum(old$n_reads==max(old$n_reads)) else 0L,
      shannon_entropy=sh$shannon_entropy_bits),as.data.table(indices))
  }
  d <- rbindlist(rows)
  counts <- merge(proportions,genes[,.(gene_name,promoter_id)],by='gene_name')
  setorder(counts,gene_name,-n_reads,cluster)
  counts[,abundance_rank:=seq_len(.N),by=promoter_id]
  setnames(counts,'group','gene_group')
  save_analysis_table(d,'gini_simpson_molecule_cluster.tsv')
  save_analysis_table(counts,'cluster_proportions.tsv')
  # 11. Promoter-level comparisons; every promoter contributes once, no depth matching.
  group_summary <- d[is.finite(gini_simpson),.(n_valid_genes=.N,
    median_gini_simpson=median(gini_simpson),q25=quantile(gini_simpson,.25),q75=quantile(gini_simpson,.75),
    iqr=IQR(gini_simpson),median_normalized_gini_simpson=finite_median(normalized_gini_simpson),
    median_inverse_simpson=median(inverse_simpson),median_cluster_count=finite_median(n_clusters),
    median_assigned_fibers=finite_median(n_fibers_assigned)),by=gene_group]
  save_analysis_table(group_summary,'simpson_group_summary.tsv')
  save_analysis_table(group_comparisons(d),'simpson_group_tests.tsv')
  message('Diversity: ',sum(is.finite(d$gini_simpson)),' valid genes; ',sum(d$n_fibers_assigned),' assigned gene/fiber observations.')
  if(mode=='diversity-only') return(invisible(d))

  # 7/8. Existing +1 and m6A features, restricted to the exact assigned read IDs.
  # +1 = first 130-160 bp nucleosome start at +1..+299 excluding TSS-overlapping calls.
  # NA plus1_start means no qualifying +1, not a missing positional value to impute.
  feature_rows <- lapply(seq_len(nrow(genes)),function(i) {
    g <- genes[i]; ids <- assigned[[i]]
    r <- saved_reads[gene==g$gene_name & RID %in% ids]
    if(length(ids) && !setequal(r$RID,ids)) stop('Missing saved per-read features for ',g$gene_name)
    pos <- r$plus1_start[is.finite(r$plus1_start)]
    stopifnot(all(pos>=1 & pos<=299))
    if(nrow(r)) stopifnot(all(r$chrom==g$chrom),all(r$tss==g$tss),all(r$strand==g$strand))
    data.table(promoter_id=g$promoter_id,n_feature_fibers=nrow(r),n_plus1_fibers=length(pos),
      plus1_fraction=if(nrow(r)) length(pos)/nrow(r) else NA_real_,
      plus1_start_mean=if(length(pos)) mean(pos) else NA_real_,
      tss_occlusion_fraction=finite_mean(as.numeric(r$nuc_at_tss)),
      m6a_promoter_proximal=finite_mean(r$m6a_upstream),m6a_tss_centered=finite_mean(r$m6a_core),
      m6a_downstream=finite_mean(r$m6a_downstream),n_m6a_proximal=sum(is.finite(r$m6a_upstream)),
      n_m6a_tss=sum(is.finite(r$m6a_core)),n_m6a_downstream=sum(is.finite(r$m6a_downstream)))
  })
  features <- rbindlist(feature_rows)
  # 9. New footprint fractions from original sample-specific calls. No subsampling.
  nuc <- new.env(parent=globalenv())
  sys.source(file.path(PROJECT,'code/entropy/nucleosome/nucleosome_entropy.R'),envir=nuc)
  # Reuse the original relative coordinate transformation unchanged.
  rel_interval <<- nuc$rel_interval
  workers <- as.integer(Sys.getenv('GINI_WORKERS','4'))
  stopifnot(is.finite(workers),workers>=1)
  fp_rows <- parallel::mclapply(seq_len(nrow(genes)),function(i) {
    g <- genes[i]; r <- saved_reads[gene==g$gene_name & RID %in% assigned[[i]]]
    message('Footprints: ',g$gene_name)
    tryCatch({
      fp <- if(nrow(r)) nuc$load_short_footprints(g,r,window=c(-100L,300L)) else
        data.table(RID=character(),start=integer(),end=integer())
      a <- footprint_presence(r,fp,g$tss,g$strand)
      a[,`:=`(promoter_id=g$promoter_id,gene_name=g$gene_name,
        status=if(nrow(r)) 'complete' else 'no_assigned_reads',error='')]
      a
    },error=function(e) data.table(promoter_id=g$promoter_id,gene_name=g$gene_name,status='error',error=conditionMessage(e)))
  },mc.cores=workers,mc.preschedule=FALSE)
  stopifnot(!any(vapply(fp_rows,inherits,logical(1),'try-error')))
  fp <- rbindlist(fp_rows,fill=TRUE)
  save_analysis_table(fp,'short_footprint_features.tsv')
  wide <- dcast(fp[status=='complete'],promoter_id~feature,value.var='fraction')
  # 10. All feature joins use unique promoter IDs, not duplicated gene values.
  d <- merge(d,features,by='promoter_id',all.x=TRUE)
  d <- merge(d,wide,by='promoter_id',all.x=TRUE)
  stopifnot(nrow(d)==nrow(genes),!anyDuplicated(d$promoter_id))
  feature_names <- c('plus1_fraction','plus1_start_mean','tss_occlusion_fraction',
    'm6a_promoter_proximal','m6a_tss_centered','m6a_downstream',
    'fp40_60_tss','fp60_80_tss','fp10_80_tss','fp40_60_downstream','fp60_80_downstream','fp10_80_downstream')
  for(f in setdiff(feature_names,names(d))) d[,(f):=NA_real_]
  save_analysis_table(d,'chromatin_features.tsv')
  associations <- feature_associations(d,feature_names)
  save_analysis_table(associations,'simpson_feature_associations.tsv')
  # 12. Representative promoters (figures are drawn in entropy_features.rmd).
  save_analysis_table(representative_promoters(d),'representative_promoters.tsv')
  if(any(fp$status=='error')) warning('Some footprint inputs failed; see short_footprint_features.tsv; missing values were not zero-filled.')
  message('Completed Gini-Simpson analysis in ',OUT)
}
if(sys.nframe()==0L) main()
