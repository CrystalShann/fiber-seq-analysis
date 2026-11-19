#!/usr/bin/env Rscript
# Same selection and count-mode CpG calculation as 03_promoter_cpg.R.
# Shared functions stay unchanged. This wrapper batches BAM queries, parallelizes
# independent promoters, checkpoints intermediate files, and records empty CpG
# output as NA rather than aborting all genes. No plots.
# Run with run_gene_examples.sh; CPG_COMPARE_WORKERS defaults to 4.
dyn.load('/project/spott/cshan/envs/dimelo/lib/libcurl.so.4',local=FALSE)
source('/project/spott/cshan/fiber-seq/code/pausing/PolII/pol2_pausing_functions.r')
utils <- file.path(PROJECT,'code/pausing/clustering_Pol2/pol2_chip_cluster_utils')
source(file.path(utils,'pol2_clustering_functions.R'))
setDTthreads(1L)
prepare_comparison_cpg <- function() {
  cfg <- yaml::read_yaml(file.path(utils,'config.yaml'))
  workers <- as.integer(Sys.getenv('CPG_COMPARE_WORKERS','4'))
  tables <- file.path(OUT_DIR,'read_features/tables')
  genes <- fread(file.path(tables,'gene_features.tsv'))[status=='complete',
    .(group,gene_name,gene_id,chrom,tss,strand,tss_source,n_reads)]
  reads <- fread(file.path(tables,'read_features.tsv.gz'),select=c('gene','RID','sample_name'))
  stopifnot(setequal(reads$gene,genes$gene_name),workers>=1L)
  manifest <- merge(reads[,.(RID=paste(gene,RID,sep='::'),original_RID=sub('^[^:]+::','',RID),
                              sample_name,tss_id=gene,cluster=gene)],
    genes[,.(tss_id=gene_name,chrom,window_start=tss+EXAMPLE_RANGE[1]-1L,window_end=tss+EXAMPLE_RANGE[2])],by='tss_id')
  manifest <- merge(manifest,fread(SAMPLES)[,.(sample_name,merged_bam_file)],by='sample_name')
  for(s in names(cfg$merged_bam_overrides)) manifest[sample_name==s,merged_bam_file:=cfg$merged_bam_overrides[[s]]]
  pol2_required(unique(manifest$merged_bam_file))
  stopifnot(nrow(manifest)==nrow(reads),!anyDuplicated(manifest$RID))
  fingerprint <- tempfile();fwrite(manifest,fingerprint,sep='\t')
  cache <- file.path('/scratch/midway3/cshan/tmp',paste0('promoter_cpg_compare_',tools::md5sum(fingerprint)))
  dir.create(cache,recursive=TRUE,showWarnings=FALSE);file.copy(fingerprint,file.path(cache,'manifest.tsv'),overwrite=TRUE);unlink(fingerprint)
  message('Checkpoint directory: ',cache)
  # Each group owns different gene BAMs; each sample/window is queried once.
  chunks <- split(genes$gene_name,rep(seq_len(workers),length.out=nrow(genes)))
  extract <- parallel::mclapply(seq_along(chunks),function(i) {
    selected <- chunks[[i]]
    ready <- file.exists(file.path(cache,paste0(selected,'.bam.bai')))
    selected <- selected[!ready]
    if(!length(selected)) return(0L)
    path <- file.path(cache,paste0('manifest_',i,'.tsv'))
    fwrite(manifest[cluster %in% selected],path,sep='\t')
    status <- system2(cfg$python,shQuote(c(file.path(PROJECT,'code/pausing/PolII/promoter_annotation/extract_promoter_bams_by_window.py'),
       '--manifest',path,'--output',cache)),stdout=file.path(cache,paste0('extraction_',i,'.log')),
       stderr=file.path(cache,paste0('extraction_',i,'.err')))
    if(status!=0L) stop('BAM extraction failed: chunk ',i)
    message('Finished BAM chunk ',i)
    status
  },mc.cores=workers)
  if(any(vapply(extract,inherits,logical(1),'try-error'))) stop('BAM extraction failed; see checkpoints')
  status <- system2(cfg$python,shQuote(c(file.path(PROJECT,'code/pausing/PolII/promoter_annotation/prepare_cpg_chromosome_inputs.py'),
    '--manifest',file.path(cache,'manifest.tsv'),'--reference',cfg$reference,'--output',cache)))
  if(status!=0L) stop('Chromosome input preparation failed')
  one <- function(i) {
    g <- genes[i];cl <- g$gene_name;prefix <- file.path(cache,cl)
    result_file <- paste0(prefix,'.summary.tsv')
    if(file.exists(result_file)) return(fread(result_file))
    status <- system2(cfg$pb_cpg_tools,shQuote(c('--bam',paste0(prefix,'.chrom.bam'),
      '--output-prefix',prefix,'--pileup-mode','count','--modsites-mode','reference','--ref',file.path(cache,'references',paste0(g$chrom,'.fa')),
      '--min-coverage','1','--threads','1')),stdout=paste0(prefix,'.log'),stderr=paste0(prefix,'.log'))
    if(status!=0L) stop('CpG extraction failed: ',cl)
    path <- paste0(prefix,'.combined.bed.gz')
    con <- gzfile(path,'rt');header <- readLines(con,n=50L);close(con)
    header <- header[startsWith(header,'#chrom\t')]
    expected <- c('#chrom','begin','end','mod_score','type','cov','mod_count','unmod_count','avg_mod_score','avg_unmod_score')
    stopifnot(length(header)==1L,identical(strsplit(header,'\t',fixed=TRUE)[[1]],expected))
    bed <- suppressWarnings(fread(cmd=paste('zcat',shQuote(path),"| sed '/^#/d'"),header=FALSE))
    out <- data.table(gene_name=cl,mod=0,cov=0,cpg=NA_real_,n_cpg_sites=0L,
                      median_cpg_coverage=NA_real_,cpg_status='no_covered_sites')
    if(nrow(bed)) {
      stopifnot(ncol(bed)==length(expected))
      setnames(bed,names(bed)[1:8],c('chrom','start','end','mod_score','type','cov','mod','unmod'))
      stopifnot(all(bed$cov>=1),all(bed$mod>=0),all(bed$mod<=bed$cov),all(abs(bed$mod+bed$unmod-bed$cov)<1e-5))
      h <- findOverlaps(pol2_gr(bed$chrom,bed$start,bed$end),
        pol2_gr(g$chrom,g$tss+EXAMPLE_RANGE[1]-1L,g$tss+EXAMPLE_RANGE[2]),ignore.strand=TRUE)
      z <- bed[unique(queryHits(h))]
      if(nrow(z)) out <- z[,.(gene_name=cl,mod=sum(mod),cov=sum(cov),cpg=sum(mod)/sum(cov),
         n_cpg_sites=.N,median_cpg_coverage=median(as.numeric(cov)),cpg_status='complete')]
    }
    fwrite(out,result_file,sep='\t',na='NA');message('CpG ',cl,': ',out$cpg_status)
    out
  }
  results <- parallel::mclapply(seq_len(nrow(genes)),one,mc.cores=workers,mc.preschedule=FALSE)
  if(any(vapply(results,inherits,logical(1),'try-error'))) stop('CpG extraction failed; completed genes checkpointed')
  out <- merge(genes,rbindlist(results),by='gene_name',all.x=TRUE)
  stopifnot(nrow(out)==nrow(genes));setorder(out,group,-cpg)
  fwrite(out,file.path(tables,'promoter_cpg.tsv'),sep='\t',na='NA')
  message('Wrote ',nrow(out),' promoters; ',sum(is.na(out$cpg)),' without CpG calls')
}
prepare_comparison_cpg()
