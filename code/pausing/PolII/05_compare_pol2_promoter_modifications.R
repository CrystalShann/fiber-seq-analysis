#!/usr/bin/env Rscript
# Tables only. Existing functions and upstream files are sourced/read unchanged.
# Run with run_gene_examples.sh. PROMOTER_COMPARE_STAGE=extract prepares Pol2
# counts; default summarize requires these plus 03_promoter_cpg.R output.
# Promoter: TSS +/-1000, 2001 bp inclusive; same full-span molecule cohort as
# read_features and promoter_cpg. Pol2 candidates: 30-60 bp inclusive, fully
# inside promoter; ChIP overlap >=1 bp. CpG is pooled per gene, NOT per molecule.
# m6a_total is the existing mean over sites observed in at least one cohort read.
# Correlations use genes as observations. BH adjustment spans all reported tests.
source('/project/spott/cshan/fiber-seq/code/pausing/PolII/04_gene_body_pol2.R')
CMP_OUT <- file.path(OUT_DIR, 'promoter_pol2_modifications/tables')
FEATURES <- file.path(OUT_DIR, 'read_features/tables')

gene_promoter_counts <- function(g, reads, peaks) {
  message('Promoter ', g$gene_name)
  d <- copy(reads[gene == g$gene_name])
  lo <- g$tss + EXAMPLE_RANGE[1]; hi <- g$tss + EXAMPLE_RANGE[2]
  stopifnot(nrow(d) == g$n_reads, all(d$start <= lo & d$end >= hi))
  d[, original_RID := sub('^[^:]+::', '', RID)]
  stopifnot(!anyDuplicated(d$original_RID))
  fp <- window_tf_footprints(g$chrom, lo, hi)
  if (!ncol(fp)) fp <- data.table(RID=character(), start0=integer(), end0=integer())
  fp <- fp[end0-start0 >= 30 & end0-start0 <= 60 & start0 >= lo-1L & end0 <= hi]
  fp <- fp[RID %in% d$original_RID]
  d[, `:=`(promoter_start=lo, promoter_end=hi, promoter_pol2=0L, promoter_pol2_chip=0L)]
  if (nrow(fp)) {
    fp[, chip := countOverlaps(GRanges(g$chrom, IRanges(start0+1L,end0)), peaks, ignore.strand=TRUE)>0L]
    z <- fp[, .(promoter_pol2=.N, promoter_pol2_chip=sum(chip)), by=.(original_RID=RID)]
    d[z, on='original_RID', `:=`(promoter_pol2=i.promoter_pol2, promoter_pol2_chip=i.promoter_pol2_chip)]
  }
  d[, promoter_chip_fraction := safe_fraction(promoter_pol2_chip,promoter_pol2)]
  fwrite(d, file.path(CMP_OUT,'per_gene',paste0(g$gene_name,'.tsv.gz')), sep='\t',na='NA')
  invisible(NULL)
}

compare_main <- function() {
  dir.create(file.path(CMP_OUT,'per_gene'),recursive=TRUE,showWarnings=FALSE)
  genes <- fread(file.path(FEATURES,'gene_features.tsv'))
  fwrite(genes[status!='complete',.(gene_name,group,status,error)],file.path(CMP_OUT,'excluded_genes.tsv'),sep='\t')
  genes <- genes[status=='complete']
  if (Sys.getenv('PROMOTER_COMPARE_STAGE','summarize')=='extract') {
    reads <- fread(file.path(FEATURES,'read_features.tsv.gz'),
      select=c('gene','RID','sample_name','start','end','m6a_total','m6a_core'))
    bed <- fread(PEAK_FILE,header=FALSE,select=1:3)
    peaks <- GRanges(bed$V1,IRanges(bed$V2+1L,bed$V3))
    result <- parallel::mclapply(seq_len(nrow(genes)),function(i)
      gene_promoter_counts(genes[i],reads,peaks),mc.cores=as.integer(Sys.getenv('PROMOTER_COMPARE_WORKERS','4')))
    if(any(vapply(result,inherits,logical(1),'try-error'))) stop('Promoter extraction failed')
    return(invisible(NULL))
  }
  molecules <- rbindlist(lapply(genes$gene_name,function(g)
    fread(file.path(CMP_OUT,'per_gene',paste0(g,'.tsv.gz')))))
  stopifnot(!anyDuplicated(molecules[,.(gene,RID)]))
  body_reads <- fread(file.path(BODY_OUT,'molecules.tsv.gz'),
    select=c('gene_name','RID','read_start','read_end','n_pol2','n_pol2_chip_overlap','covered_gene_bp','pol2_per_kb_covered'))
  setnames(body_reads,c('gene_name','n_pol2','n_pol2_chip_overlap','covered_gene_bp','pol2_per_kb_covered'),
    c('gene','body_pol2','body_pol2_chip','body_covered_bp','body_pol2_per_kb'))
  molecules <- merge(molecules,body_reads,by=c('gene','RID'),all.x=TRUE)
  molecules[,same_body_alignment:=!is.na(read_start) & start==read_start & end==read_end]
  molecules[same_body_alignment==FALSE, `:=`(body_pol2=NA_integer_,body_pol2_chip=NA_integer_,
                                      body_covered_bp=NA_integer_,body_pol2_per_kb=NA_real_)]
  s <- molecules[,.(n_promoter_molecules=.N,promoter_pol2=sum(promoter_pol2),
    promoter_pol2_chip=sum(promoter_pol2_chip),
    promoter_chip_fraction=safe_fraction(sum(promoter_pol2_chip),sum(promoter_pol2)),
    promoter_pol2_per_kb=sum(promoter_pol2)/(.N*2.001),
    promoter_fraction_with_pol2=mean(promoter_pol2>0),
    promoter_m6a=mean(m6a_total),promoter_m6a_core=mean(m6a_core),
    n_matched_body_molecules=sum(!is.na(body_pol2)),
    matched_body_pol2_per_kb=1000*sum(body_pol2,na.rm=TRUE)/sum(as.numeric(body_covered_bp),na.rm=TRUE)),by=.(gene_name=gene)]
  cpg <- fread(file.path(FEATURES,'promoter_cpg.tsv'))
  stopifnot(setequal(cpg$gene_name,genes$gene_name))
  check <- merge(genes[,.(gene_name,chrom,tss,strand,n_reads)],cpg,by='gene_name',suffixes=c('_m6a','_cpg'))
  for(k in c('chrom','tss','strand','n_reads')) stopifnot(identical(check[[paste0(k,'_m6a')]],check[[paste0(k,'_cpg')]]))
  s <- merge(genes[,.(gene_name,gene_id,group,chrom,tss,strand,tss_source,tpm_mean)],s,by='gene_name')
  s <- merge(s,cpg[,.(gene_name,cpg_mod=mod,cpg_cov=cov,promoter_cpg=cpg,n_cpg_sites,cpg_status)],by='gene_name')
  body <- fread(file.path(BODY_OUT,'genes.tsv'))[,.(gene_name,body_pol2_per_kb=pol2_per_kb_covered,
    body_chip_fraction=fraction_pol2_chip_overlap,body_mean_pol2=mean_pol2_per_molecule)]
  s <- merge(s,body,by='gene_name',all.x=TRUE)
  s[, `:=`(promoter_start=tss+EXAMPLE_RANGE[1],promoter_end=tss+EXAMPLE_RANGE[2])]
  stopifnot(nrow(s)==nrow(genes), all(s$n_promoter_molecules==genes$n_reads[match(s$gene_name,genes$gene_name)]))
  stopifnot(max(abs(s$promoter_m6a-genes$m6a_total[match(s$gene_name,genes$gene_name)]))<1e-12,
            all(s$promoter_cpg>=0 & s$promoter_cpg<=1,na.rm=TRUE),
            all(s$cpg_mod<=s$cpg_cov,na.rm=TRUE))
  groups <- s[,.(n_genes=.N,n_promoter_molecules=sum(n_promoter_molecules),
    promoter_pol2=sum(promoter_pol2),promoter_chip_fraction=sum(promoter_pol2_chip)/sum(promoter_pol2),
    mean_gene_promoter_pol2_per_kb=mean(promoter_pol2_per_kb),
    mean_gene_body_pol2_per_kb=mean(body_pol2_per_kb),
    mean_gene_promoter_m6a=mean(promoter_m6a),
    median_gene_promoter_m6a=median(promoter_m6a),
    n_genes_with_cpg=sum(is.finite(promoter_cpg)),
    mean_gene_promoter_cpg=mean(promoter_cpg,na.rm=TRUE),
    median_gene_promoter_cpg=median(promoter_cpg,na.rm=TRUE),
    pooled_promoter_cpg=sum(cpg_mod,na.rm=TRUE)/sum(cpg_cov,na.rm=TRUE)),by=group]
  pairs <- CJ(pol2_metric=c('body_pol2_per_kb','body_chip_fraction','matched_body_pol2_per_kb',
                            'promoter_pol2_per_kb','promoter_chip_fraction'),
              modification=c('promoter_m6a','promoter_cpg'))
  corrs <- rbindlist(lapply(c('all',unique(s$group)),function(grp) {
    z <- if(grp=='all') s else s[group==grp]
    rbindlist(lapply(seq_len(nrow(pairs)),function(i) {
      x <- z[[pairs$pol2_metric[i]]]; y <- z[[pairs$modification[i]]]; keep <- is.finite(x)&is.finite(y)
      ct <- if(sum(keep)>=3 && uniqueN(x[keep])>1 && uniqueN(y[keep])>1)
        suppressWarnings(cor.test(x[keep],y[keep],method='spearman',exact=FALSE)) else NULL
      data.table(group=grp,pol2_metric=pairs$pol2_metric[i],modification=pairs$modification[i],
                 n_genes=sum(keep),spearman_rho=if(is.null(ct)) NA_real_ else unname(ct$estimate),
                 p_value=if(is.null(ct)) NA_real_ else ct$p.value)
    }))
  }))
  corrs[,p_adj_BH:=p.adjust(p_value,method='BH')]
  for(nm in c('genes','groups','correlations','molecules'))
    fwrite(switch(nm,genes=s,groups=groups,correlations=corrs,molecules=molecules),
      file.path(CMP_OUT,paste0(nm,if(nm=='molecules') '.tsv.gz' else '.tsv')),sep='\t',na='NA')
  fwrite(data.table(field=c('window','m6a','cpg','pol2','groups','correlations'),
    definition=c('TSS -1000 through +1000 inclusive; same full-span read cohort as read_features',
      'Existing m6a_total: per-read fraction over sites observed in at least one cohort molecule; averaged per gene',
      'pb-CpG-tools count-mode methylated calls / covered calls; same input promoter molecules, default MAPQ >=1 filtering; pooled per gene; no covered calls is NA',
      'Unique 30-60 bp inclusive footprints fully contained in promoter; peak overlap >=1 bp',
      'mean_gene columns weight each gene equally; pooled CpG weights covered calls; ChIP fraction weights footprints',
      'Spearman across genes; all and within each group; BH correction across all reported tests')),
    file.path(CMP_OUT,'definitions.tsv'),sep='\t')
  print(groups);print(corrs[group=='all']);message('Wrote ',nrow(s),' matched promoter comparisons')
}
if(sys.nframe()==0L) compare_main()
