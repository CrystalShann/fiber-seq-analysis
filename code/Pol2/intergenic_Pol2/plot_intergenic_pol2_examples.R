# Intergenic Pol II figures (PDF only): 500-bp consensus-cluster examples, then summaries by
# FIRE / footprint rank.
#   Rscript --vanilla plot_intergenic_pol2_examples.R              # examples + summary bars
#   Rscript --vanilla plot_intergenic_pol2_examples.R --bars-only  # summary figures, no clustering
# Reuse shared helpers without modifying them.
# Output policy: PDF figures only; no previews, tables, notes, or matrix caches.
suppressPackageStartupMessages({library(data.table); library(Matrix); library(ggplot2); library(ComplexHeatmap)})
root <- '/project/spott/cshan/fiber-seq'
source(file.path(root, 'code/parsing_functions/parsing_footprints_functions.r'))
# leiden_manhattan_plots.r also sources code/parsing_functions/plotting_functions.r
source(file.path(root, 'code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r'))
source(file.path(root, 'code/clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r'))
tables <- file.path(root, 'LCL_project/Pol2/intergenic_Pol2/tables')
out <- file.path(root, 'LCL_project/Pol2/intergenic_Pol2/plots')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
regions <- fread(file.path(tables, 'chromatin_regions.tsv.gz'), na.strings='')

#################################################
# 500-bp consensus examples (centre +/- 250 bp, the same windows as the rank metrics)
#################################################
# Each locus independently: unbinned m6A, Manhattan k=15, resolution=0.5,
# consensus seeds 1:100, tau=0.5. All full-span fibers; no footprint selection.
plot_consensus_examples <- function() {
  samples <- as.data.frame(fread('/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv')[,.(sample_name,fire_dir)])
  sample_colors <- sample_palette(samples$sample_name)
  ids <- c('promoter_chr8:94719547-94719834','pol2_10778','pol2_40246','pol2_79242')
  titles <- c('Promoter: PLS / active TSS','Intergenic: dELS / active enhancer 1','Intergenic: dELS / upstream TSS flank','Closed atlas: dELS / quiescent')
  # Resolve the actual center-overlapping ENCODE cCRE, not the whole analysis window.
  ccres <- fread('/project/spott/cshan/annotations/GRCh38-cCREs.bed',header=FALSE,
    col.names=c('chrom','start','end','dhs','ccre_id','class'))
  fp_root <- '/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/FiberHMM_v2/results/results_v1_trained_model'
  for (i in seq_along(ids)) {
   x <- regions[region_id==ids[i]]
   r <- data.frame(region_id=ids[i],chr=x$chrom,start=x$center-250L,end=x$center+250L,
                   analysis_start=x$center-249L,analysis_end=x$center+250L,width=500L,
                   region_type='intergenic',annotation=titles[i])
   cre <- ccres[chrom==r$chr & start<=x$center & end>x$center]
   stopifnot(nrow(cre)==1L,cre$class==x$ccre_center)
   message('CRE ',cre$ccre_id,' ',cre$class,': ',cre$start,'-',cre$end)
   message('Assembling ',ids[i],': ',r$chr,':',r$start,'-',r$end)
   a <- assemble_region_m6a(samples,r,full_span=TRUE,matrix_dir=NULL)
   m <- a$met_mat
   feat <- bin_read_matrix(m,r$analysis_start,r$analysis_end,window_size=0)$mat
   knn <- manhattan_knn_graph(feat,k_neighbors=15,verbose=FALSE)
   cons <- leiden_consensus_partition(knn$graph,resolution=0.5,seed=1,consensus=TRUE,n_partitions=100,tau=.5)
   if (!cons$converged) stop('Consensus did not converge for ',ids[i])
   cl <- relabel_clusters_by_size(cons$membership[rownames(m)])
   assignments <- a$rids_df
   assignments$cluster <- cl
   result <- list(assignments=assignments,site_met_mat=m,region=r)
   nuc <- extract_nucleosomes(samples,r,assignments,min_size=130L,max_size=160L)
   nuc_profiles <- footprint_profiles(nuc,'ft_nuc_130-160bp',assignments,r)
   cluster_cols <- cluster_id_palette(levels(cl))
   title <- sprintf('%s | %s\n%s:%s-%s (500 bp; %d fibers)\nConsensus Leiden: k=15, resolution=0.5, 100 seeds, tau=0.5; %d clusters',
      titles[i],ids[i],r$chr,format(r$analysis_start,big.mark=',',scientific=FALSE),
      format(r$analysis_end,big.mark=',',scientific=FALSE),nrow(m),nlevels(cl))
   # The genomic layout supplies cluster/sample annotations and colored cluster profiles.
   # Its row order matches plot_smf_reads(): cluster, read start, RID.
   ht <- plot_read_heatmap(result,layout='genomic',sample_colors=sample_colors,
      show_cluster_profiles=TRUE,cluster_colors=cluster_id_palette(levels(cl)),
      heatmap_options=list(column_title=title,width=grid::unit(170,'mm'),height=grid::unit(145,'mm'),
        row_title=sprintf('%s (n=%d)',levels(cl),as.integer(table(cl))),
        raster_by_magick=FALSE))
   # Retain the pooled aggregate above the colored cluster-specific aggregate tracks.
   pos <- seq.int(r$analysis_start,r$analysis_end)
   profile <- rep(NA_real_,length(pos))
   profile[match(as.integer(colnames(m)),pos)] <- Matrix::colMeans(m)
   pooled <- HeatmapAnnotation(pooled=methylation_area_annotation(profile,'#444444',grid::unit(16,'mm')),
      annotation_label='Pooled m6A',annotation_name_gp=grid::gpar(fontsize=8))
   cre_values <- ifelse(pos>cre$start & pos<=cre$end,cre$class,'outside')
   cre_top <- HeatmapAnnotation(CRE=anno_simple(cre_values,
      col=c(outside='white',setNames(CCRE_COLORS[[cre$class]],cre$class)),height=grid::unit(4,'mm')),
      annotation_label=paste(cre$ccre_id,cre$class),annotation_name_gp=grid::gpar(fontsize=8))
   top <- c(pooled,ht@top_annotation,cre_top)
   profile_names <- c('pooled',paste0(levels(cl),' m6A'))
   # Highlight exact cCRE BED bounds on each genomic-position panel.
   cre_left <- max(0,(cre$start-r$start)/r$width)
   cre_right <- min(1,(cre$end-r$start)/r$width)
   highlight_grid <- function() {
     grid::grid.rect(x=grid::unit((cre_left+cre_right)/2,'npc'),width=grid::unit(cre_right-cre_left,'npc'),
       gp=grid::gpar(fill='#FFD54F',alpha=.15,col=NA))
     grid::grid.segments(x0=grid::unit(c(cre_left,cre_right),'npc'),x1=grid::unit(c(cre_left,cre_right),'npc'),
       y0=0,y1=1,gp=grid::gpar(col='#B8860B',lty=2,lwd=.6))
   }
   o <- order(assignments$cluster,assignments$start,assignments$RID)
   annotation <- rowAnnotation(cluster=cl[o],
      sample=factor(sub('_.*$','',assignments$sample_name[o]),levels=names(sample_colors)),
      col=list(cluster=cluster_id_palette(levels(cl)),sample=sample_colors),
      annotation_legend_param=list(cluster=list(nrow=1),
        sample=list(ncol=6,labels_gp=grid::gpar(fontsize=8),grid_width=grid::unit(3,'mm'))))
   ht <- plot_read_heatmap(result,layout='genomic',sample_colors=sample_colors,
      show_cluster_profiles=FALSE,
      heatmap_options=list(column_title=title,top_annotation=top,left_annotation=annotation,
        width=grid::unit(170,'mm'),height=grid::unit(145,'mm'),
        row_title=sprintf('%s (n=%d)',levels(cl),as.integer(table(cl))),raster_by_magick=FALSE))
   # e.g. promoter_Pol2_chr8:94,719,702 or intergenic_FIRE-D10_fp-D8_chr1:161,738,901
   # (group, or FIRE / footprint deciles for accessible peaks, + window centre)
   prefix <- if (x$group=='intergenic accessible') sprintf('intergenic_FIRE-D%d_fp-D%d',x$fire_decile,x$fp_decile) else
     gsub(' ','_',sub('Pol II','Pol2',x$group))
   stem <- paste0(prefix,'_',r$chr,':',format(x$center,big.mark=',',scientific=FALSE))
   region_out <- file.path(out,stem)
   dir.create(region_out,recursive=TRUE,showWarnings=FALSE)
   save_figure(ht,file.path(region_out,paste0(stem,'.pdf')),width=12,height=18,
      draw=function(h) {
        ComplexHeatmap::draw(h,heatmap_legend_side='bottom',annotation_legend_side='bottom',merge_legend=TRUE)
        for (k in seq_len(nlevels(cl))) ComplexHeatmap::decorate_heatmap_body('m6A',highlight_grid(),slice=k)
        for (nm in profile_names) ComplexHeatmap::decorate_annotation(nm,highlight_grid())
      })
   message('Extracting footprints for ',ids[i])
   tf <- lapply(seq_len(nrow(samples)),function(j) {
     sam <- samples$sample_name[j]
     reads <- assignments[assignments$sample_name==sam,,drop=FALSE]
     path <- file.path(fp_root,sam,paste0(sam,'.recalled_tf.bed.gz'))
     raw <- read_footprint_region(path,'bed15_fiberhmm_tf',r,reads$original_RID)
     if (!nrow(raw)) return(NULL)
     b <- as.data.table(convert_ft_bed12_to_bed6(raw,format='bed15_fiberhmm_tf',source=path))
     b[,size:=end-start]
     b <- b[start<r$end & end>r$start & size>=10 & size<=80]
     b[,original_RID:=as.character(RID)]
     b[,RID:=paste(sam,original_RID,sep='::')]
     b[,track:=fcase(size<30,'TF_10_29',size<=60,'TF_30_60',default='TF_61_80')]
     b[,.(RID,original_RID,chr,start,end,size,track)]
   })
   records <- as.data.frame(rbindlist(c(list(nuc),tf),fill=TRUE))
   stopifnot(all(records$RID %in% assignments$RID))
   tracks <- c('ft_nuc_130-160bp','TF_10_29','TF_30_60','TF_61_80')
   colors <- setNames(unname(LCL_TRACK_COLORS[c('ft_nuc_130-160bp','FiberHMM_10-30bp','FiberHMM_40-60bp','FiberHMM_60-80bp')]),tracks)
   labels <- setNames(c('Nucleosome 130-160 bp','Footprint 10-29 bp','Footprint 30-60 bp','Footprint 61-80 bp'),tracks)
   p <- plot_smf_reads(result,records,sample_colors=sample_colors,
       tracks=tracks,track_colors=colors,track_labels=labels)+
       labs(title=paste(titles[i],ids[i],sep=' | '),
         subtitle=sprintf('500-bp window; %d fibers; same consensus clusters and row order as m6A heatmap',nrow(m)),
         caption='Fibertools nucleosomes (130-160 bp); sample-matched FiberHMM footprints. Footprint size does not identify Pol II.')
   anchor <- plot_anchor(r)
   bounds <- c(anchor$left-.085*(anchor$right-anchor$left),anchor$right+.5)
   # cCRE bounds clipped to the window
   cleft <- max(cre$start+.5-anchor$anchor,anchor$left)
   cright <- min(cre$end+.5-anchor$anchor,anchor$right)
   band <- list(xmin=cleft,xmax=cright,fill='#FFD54F',alpha=.15,position='under')
   p <- p + annotate('rect',xmin=cleft,xmax=cright,ymin=-Inf,ymax=Inf,fill='#FFD54F',alpha=.15)+
        geom_vline(xintercept=c(cleft,cright),linetype='dashed',color='#B8860B',linewidth=.3)
   cp <- data.frame(start=cleft,end=cright,class=cre$class)
   cre_panel <- plot_interval_track(cp,window=bounds,fill_col='class',colors=CCRE_COLORS,
      labels=list(x=NULL,y=NULL,title=paste('ENCODE cCRE:',cre$ccre_id,cre$class)),
      scales=list(scale_x_continuous(limits=bounds,expand=c(0,0))),theme=theme_fiberseq())+
      theme(legend.position='none',axis.text.y=element_blank(),axis.ticks.y=element_blank())
   np <- plot_group_profile(nuc_profiles,window=bounds,group_col='cluster',facet_col='cluster',
      style='area',colors=cluster_cols,fill_colors=cluster_cols,region=r,coordinates='relative',
      highlight=band,scales=list(scale_x_continuous(limits=bounds,expand=c(0,0)),
        scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=c(0,0))),
      labels=list(x=NULL,y='Nucleosome occupancy',title='Per-cluster nucleosome profiles (130-160 bp)'),
      theme=theme_fiberseq(base_size=9)+theme(legend.position='none'))
   combined <- cowplot::plot_grid(cre_panel,np,p,ncol=1,align='v',axis='lr',rel_heights=c(.7,3.5,10))
   save_figure(combined,file.path(region_out,paste0(stem,'_footprints.pdf')),width=12,height=19)
   # Remove superseded flat-folder copies only after both new PDFs are saved.
   old <- file.path(out,paste0(stem,c('.pdf','_footprints.pdf')))
   unlink(old[file.exists(old)])
   message(sprintf('DONE %s: fibers=%d; clusters=%s; consensus rounds=%d; footprint blocks=%d',
     ids[i],nrow(m),paste(as.integer(table(cl)),collapse=','),cons$n_rounds,nrow(records)))
  }
}
if (!'--bars-only' %in% commandArgs(TRUE)) plot_consensus_examples()

#################################################
# summary figures by FIRE / footprint rank
#################################################
# x = promoter, D1 (lowest) ... D10 (highest decile among accessible peaks), closed;
# one facet per ranking axis
bins <- c('promoter',paste0('D',1:10),'closed')
bin_cols <- setNames(c('#386CB0',colorRampPalette(c('#FDD49E','#B30000'))(10),'#999999'),bins)
axis_labels <- c(FIRE='Ranked by FIRE frequency',footprint='Ranked by 30-60 bp footprint frequency')
rank_axis <- function(x) x[,`:=`(bin=factor(bin,bins),axis=factor(axis,names(axis_labels)))][]
rank_facet <- facet_wrap(~axis,ncol=1,labeller=labeller(axis=axis_labels))
rank_caption <- paste('Deciles of accessible intergenic peaks (>=50% inside the LCL FIRE union);',
  '500-bp windows around the centroid; pooled 31 LCL samples.')
state_labels <- c(
 '1_TssA'='Active TSS', '2_TssFlnk'='Flanking active TSS',
 '3_TssFlnkU'='Upstream TSS flank', '4_TssFlnkD'='Downstream TSS flank',
 '5_Tx'='Strong transcription', '6_TxWk'='Weak transcription',
 '7_EnhG1'='Genic enhancer 1', '8_EnhG2'='Genic enhancer 2',
 '9_EnhA1'='Active enhancer 1', '10_EnhA2'='Active enhancer 2',
 '11_EnhWk'='Weak enhancer', '12_ZNF/Rpts'='Zinc-finger genes / repeats',
 '13_Het'='Heterochromatin', '14_TssBiv'='Bivalent / poised TSS',
 '15_EnhBiv'='Bivalent enhancer', '16_ReprPC'='Strong Polycomb repression',
 '17_ReprPCWk'='Weak Polycomb repression', '18_Quies'='Quiescent / low signal',
 'no_state'='No lifted state')

g <- rank_axis(fread(file.path(tables,'m6a_summary.tsv')))
p <- ggplot(g,aes(bin,median_m6a_prop,fill=bin))+geom_col(width=.75)+
 geom_errorbar(aes(ymin=q25,ymax=q75),width=.25,linewidth=.3)+rank_facet+scale_fill_manual(values=bin_cols)+
 labs(x=NULL,y='Median per-region m6A metric (IQR)',title='m6A by FIRE and footprint rank',
 caption=paste0(rank_caption,'\nClosed: 5,000 sampled peaks. Denominator: fibers x observed-call positions, not all callable adenines.'))+
 theme_fiberseq(base_size=11)+theme(legend.position='none')
save_figure(p,file.path(out,'m6a_by_rank.pdf'),width=9,height=7)

r <- rank_axis(fread(file.path(tables,'rank_bin_summary.tsv')))
p <- ggplot(r[bin!='promoter'],aes(bin,frac_pol2_chip,fill=bin))+geom_col(width=.75)+
 geom_text(aes(label=format(n,big.mark=',')),vjust=-.3,size=2.6)+rank_facet+scale_fill_manual(values=bin_cols)+
 scale_y_continuous(expand=expansion(mult=c(0,.1)))+
 labs(x=NULL,y='Fraction with a Pol II ChIP summit within 250 bp',title='GM12878 Pol II ChIP by FIRE and footprint rank',
 caption=paste0(rank_caption,'\nLabels: peaks per bin. Promoters omitted (Pol II ChIP summits within 250 bp of a TSS by definition).'))+
 theme_fiberseq(base_size=11)+theme(legend.position='none')
save_figure(p,file.path(out,'pol2_chip_by_rank.pdf'),width=9,height=7)

cc <- rank_axis(fread(file.path(tables,'ccre_summary.tsv'))[rule=='center'])
cc <- melt(cc,id.vars=c('rule','axis','bin','n'),variable.name='class',value.name='pct')
ccre_levels <- c('PLS','pELS','dELS','CA-H3K4me3','CA-CTCF','CA-TF','CA','TF','none')
cc <- cc[!is.na(pct)][,class:=factor(class,ccre_levels)]
p <- ggplot(cc,aes(bin,pct,fill=class))+geom_col(width=.8)+rank_facet+
 scale_fill_manual(values=c(CCRE_COLORS[setdiff(ccre_levels,'none')],none='grey85'),name='cCRE at centre')+
 scale_y_continuous(expand=c(0,0))+
 labs(x=NULL,y='% of peaks',title='ENCODE cCRE class at the centre by FIRE and footprint rank',caption=rank_caption)+
 theme_fiberseq(base_size=11)
save_figure(p,file.path(out,'ccre_by_rank.pdf'),width=10,height=7)

hm <- rank_axis(fread(file.path(tables,'chromhmm_center_summary.tsv')))
hm <- melt(hm,id.vars=c('axis','bin','n'),variable.name='state',value.name='pct')
hm[,state:=factor(state,rev(names(state_labels)))]
stopifnot(!anyNA(hm$state))
p <- ggplot(hm,aes(bin,state,fill=pct))+geom_tile(color='white')+
 geom_text(aes(label=ifelse(pct>=1,sprintf('%.0f',pct),'')),size=2.3)+
 facet_wrap(~axis,ncol=2,labeller=labeller(axis=axis_labels))+
 scale_y_discrete(labels=state_labels)+scale_fill_gradient(low='white',high='#B30000',name='% of peaks')+
 labs(x=NULL,y=NULL,title='E116 ChromHMM state at the centre by FIRE and footprint rank',
 caption=paste0(rank_caption,'\nNumbers: % of peaks in the bin (shown when >= 1%).'))+
 theme_fiberseq(base_size=10)
save_figure(p,file.path(out,'chromhmm_by_rank.pdf'),width=14,height=7)

gr <- fread(file.path(tables,'fire_fp_grid_summary.tsv'))
gr <- melt(gr,id.vars=c('fire_decile','fp_decile','n'),variable.name='metric',value.name='value')
grid_plot <- function(m,title) {
  ggplot(gr[metric==m],aes(factor(fire_decile),factor(fp_decile),fill=value))+geom_tile(color='white')+
   geom_text(aes(label=n),size=2.2)+scale_fill_viridis_c(name=NULL)+coord_equal()+
   labs(x='FIRE frequency decile',y='Footprint frequency decile',title=title)+theme_fiberseq(base_size=10)
}
p <- cowplot::plot_grid(grid_plot('frac_pol2_chip','Fraction with a Pol II ChIP summit within 250 bp'),
  grid_plot('median_m6a_prop','Median m6A metric'),nrow=1)
p <- cowplot::plot_grid(p,cowplot::ggdraw()+cowplot::draw_label(
  paste('Accessible intergenic peaks; numbers = peaks per tile.',rank_caption),size=8),ncol=1,rel_heights=c(1,.06))
save_figure(p,file.path(out,'fire_fp_grid.pdf'),width=13,height=6.5)
