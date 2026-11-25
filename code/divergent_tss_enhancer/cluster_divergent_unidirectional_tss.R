# Joint, center-aligned TSS chromatin clustering. Source without running;
# Rscript .../cluster_divergent_unidirectional_tss.R [--pilot] executes main().
ROOT <- '/project/spott/cshan/fiber-seq'
source(file.path(ROOT, 'code/parsing_functions/parsing_footprints_functions.r'))
source(file.path(ROOT, 'code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r'))
source(file.path(ROOT, 'code/clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r'))
suppressPackageStartupMessages({library(ggplot2); library(data.table)})

########
# Settings: input paths, output folder and clustering parameters
########

# 1. Point to the PINTS TSS calls, the 31-sample LCL table, chromosome sizes
#    and the GRO-cap bigWig list (inputs.tsv).
# 2. Set the output folder (clustering_results).
# 3. Set the parameters: 5,000 reads per class, seed, k = 15 neighbours,
#    Leiden resolution 0.5, 100 consensus partitions, tau 0.5, at most 10 rounds.
# 4. orient_strand: promoters are pooled, so minus-strand unidirectional TSSs
#    are flipped (transcription to the right). Divergent pairs keep genome
#    orientation: minus TSS on the left, plus TSS on the right in every pair.
tss_config <- function() list(
  pints = file.path(ROOT, 'LCL_project/divergent_tss_enhancer/GM12878_GROcap_PINTS_outputs'),
  output = file.path(ROOT, 'LCL_project/divergent_tss_enhancer/clustering_results'),
  samples = '/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv',
  sizes = '/project/spott/cshan/annotations/hg38.chrom.sizes',
  inputs = file.path(ROOT, 'code/divergent_tss_enhancer/inputs.tsv'),
  n = 5000L, seed = 1L, k = 15L, resolution = 0.5, partitions = 100L,
  tau = 0.5, max_rounds = 10L, orient_strand = TRUE)

########
# Build the 1-kb windows for divergent and unidirectional TSSs
########

# 1. Read the PINTS divergent pairs and keep the high-confidence ones; the
#    centre is the midpoint between the plus and minus TSS (rounded down).
# 2. Read the PINTS unidirectional peaks; the centre is the peak summit
#    (column 8) and the peak strand says whether it is a plus or minus TSS.
# 3. TSS strand: "+" for divergent pairs, the peak strand for unidirectional.
# 4. Make the window centre - 500 to centre + 500 (BED, 0-based start); when
#    orienting, minus-strand windows move 1 bp right (centre -499..+500) so
#    that, flipped, they cover TSS -500..+499 like the others.
# 5. Drop windows that run off a chromosome end.
# 6. Add the 1-based analysis start/end and check every window is 1,000 bp.
tss_regions <- function(cfg) {
  d <- fread(file.path(cfg$pints, 'divergent_regions.tsv'), data.table = FALSE)
  d <- d[d$high_confidence == 1, ]
  u <- fread(file.path(cfg$pints, 'GM12878_GROcap_GRCh38_1_unidirectional_peaks.bed'), header = FALSE, data.table = FALSE)
  stopifnot(ncol(u) >= 9, all(u$V6 %in% c('+', '-')))
  r <- rbind(data.frame(region_id = d$region_id, chr = d$chrom, class = 'divergent',
    center = floor((d$plus_tss_0based + d$minus_tss_0based)/2),
    plus = d$plus_tss_0based, minus = d$minus_tss_0based, tss_strand = '+'),
    data.frame(region_id = paste0('unidirectional_', seq_len(nrow(u))), chr = u$V1,
    class = 'unidirectional', center = u$V8,
    plus = ifelse(u$V6 == '+', u$V8, NA), minus = ifelse(u$V6 == '-', u$V8, NA), tss_strand = u$V6))
  sizes <- fread(cfg$sizes, header = FALSE)
  shift <- if (cfg$orient_strand) as.integer(r$tss_strand == '-') else 0L
  r$start <- r$center - 500L + shift; r$end <- r$center + 500L + shift
  keep <- !is.na(match(r$chr, sizes$V1)) & r$start >= 0 & r$end <= sizes$V2[match(r$chr, sizes$V1)]
  message('Boundary-excluded windows: ', sum(!keep))
  r <- r[keep, ]; r$analysis_start <- r$start + 1L; r$analysis_end <- r$end
  r$width <- 1000L
  stopifnot(!anyDuplicated(r$region_id), all(r$end-r$start == 1000))
  r
}

########
# Keep a uniform random sample of at most n reads, added in batches
########

# 1. Give each new read a random number (its priority).
# 2. Stack the new reads under the reads already kept.
# 3. Keep the n reads with the highest priorities, so the sample is uniform
#    whatever the order of the loci or how many reads each locus has.
# 4. Count every read ever offered (seen).
reservoir_add <- function(old, feature, meta, n) {
  keys <- runif(nrow(meta))
  if (is.null(old)) old <- list(feature = feature[FALSE,,drop=FALSE], meta = meta[FALSE,,drop=FALSE], keys = numeric(), seen = 0)
  feature <- rbind(old$feature, feature); meta <- rbind(old$meta, meta)
  keys <- c(old$keys, keys); take <- head(order(keys, decreasing = TRUE), n)
  list(feature = feature[take,,drop=FALSE], meta = meta[take,,drop=FALSE],
       keys = keys[take], seen = old$seen + length(keys)-length(old$keys))
}

########
# Read m6A at every window and sample 5,000 full-span reads per class
########

# 1. For each window, read the m6A calls of every read that spans the whole
#    1 kb, from all 31 samples (assemble_region_m6a()).
# 2. If the shared reader refuses the window (fewer than 3 reads or no m6A
#    call), read it again with the same parsers so those reads can still be
#    sampled; skip windows with no spanning read at all.
# 3. Bin the per-bp calls into 40 bins of 25 bp (fraction of bases with a call).
# 4. Label each read with its window, class and TSS strand; read ID =
#    window::sample::read. Features stay in genome orientation here; the
#    clustering flips the minus-strand reads.
# 5. Add the reads to that class's random sample (reservoir_add()).
# 6. Stop if either class ends with fewer than 5,000 reads; return the binned
#    features, the read table and the number of reads seen per class.
tss_sample <- function(regions, samples, cfg) {
  set.seed(cfg$seed); pools <- list(); skipped <- 0L
  for (i in seq_len(nrow(regions))) {
    r <- regions[i,,drop=FALSE]
    dat <- tryCatch(assemble_region_m6a(samples, r, full_span = TRUE,
      positions = seq.int(r$analysis_start, r$analysis_end), matrix_dir = NULL),
      error = function(e) {
        if (grepl('^Insufficient full-span reads/', conditionMessage(e))) return(NULL)
        stop(e)
      })
    # The shared assembler requires >=3 reads and >=1 call. Recover small/zero-call
    # windows with the same underlying parsers to retain their sampling eligibility.
    if (is.null(dat)) {
      query <- GenomicRanges::GRanges(r$chr, IRanges::IRanges(r$analysis_start, r$analysis_end))
      mm <- rr <- list()
      for (j in seq_len(nrow(samples))) {
        path <- extracted_path(samples[j,,drop=FALSE], r$chr, 'm6a')
        b <- read_ft_bed12(path, query, longest_alignment = TRUE)
        b <- b[b$start <= r$start & b$end >= r$end,,drop=FALSE]
        if (!nrow(b)) next
        original <- as.character(b$RID); b$RID <- paste(samples$sample_name[j],original,sep='::')
        md <- data.frame(RID=b$RID, original_RID=original, chr=b$chr,
          start=b$start+1L,end=b$end,strand=b$strand,sample_name=samples$sample_name[j],score=b$score)
        blocks <- convert_ft_bed12_to_bed6(b,format='bed12_fibertools',source=path)
        calls <- data.frame(RID=blocks$RID,pos=blocks$end)
        mm[[j]] <- build_met_mat(calls,md,seq.int(r$analysis_start,r$analysis_end)); rr[[j]] <- md
      }
      if (!length(rr)) {skipped <- skipped+1L; next}
      dat <- list(met_mat=do.call(rbind,mm),rids_df=dplyr::bind_rows(rr))
    }
    stopifnot(all(dat$rids_df$start <= r$analysis_start), all(dat$rids_df$end >= r$analysis_end))
    colnames(dat$met_mat) <- as.character(-500:499)
    feat <- bin_read_matrix(dat$met_mat, -500L, 499L, 25L)$mat
    md <- dat$rids_df; md$molecule_id <- md$RID
    md$region_id <- r$region_id; md$class <- r$class; md$tss_strand <- r$tss_strand
    md$RID <- paste(r$region_id,md$molecule_id,sep='::')
    rownames(feat) <- md$RID
    stopifnot(ncol(feat)==40L, !anyNA(feat))
    pools[[r$class]] <- reservoir_add(pools[[r$class]],feat,md,cfg$n)
    if (i %% 100L == 0L) message('Extracted ',i,'/',nrow(regions),' windows')
  }
  message('No-coverage windows: ',skipped)
  for (cl in c('divergent','unidirectional')) {
    message(cl, ': eligible observations = ', pools[[cl]]$seen)
    if (is.null(pools[[cl]]) || nrow(pools[[cl]]$meta) < cfg$n) stop('Fewer than requested reads in ',cl)
  }
  list(feature=do.call(rbind,lapply(pools,`[[`,'feature')),
       meta=dplyr::bind_rows(lapply(pools,`[[`,'meta')),
       eligible=vapply(pools,`[[`,numeric(1),'seen'))
}

########
# Per-bp m6A of the sampled reads, for the heatmap and profiles
########

# 1. For each window with sampled reads, re-read those reads' m6A calls from
#    each sample's fibertools file.
# 2. Check every sampled read is found exactly once.
# 3. Keep the calls inside the window and record (read, position); flip the
#    positions of minus-strand reads when the clustering was oriented.
# 4. Build a reads x 1,000 bp 0/1 matrix (columns = bp from the TSS / centre).
# 5. Check that binning this matrix gives back exactly the clustering features.
tss_site_m6a <- function(regions, samples, result) {
  a <- result$assignments; hits <- list()
  used <- regions[regions$region_id %in% a$region_id, ]
  for (i in seq_len(nrow(used))) {
    r <- used[i,,drop=FALSE]
    query <- GenomicRanges::GRanges(r$chr,IRanges::IRanges(r$analysis_start,r$analysis_end))
    ix <- which(a$region_id==r$region_id)
    for (s in unique(a$sample_name[ix])) {
      path <- extracted_path(samples[samples$sample_name==s,,drop=FALSE],r$chr,'m6a')
      rows <- ix[a$sample_name[ix]==s]
      b <- read_ft_bed12(path,query,longest_alignment=TRUE)
      b <- b[as.character(b$RID) %in% a$original_RID[rows],,drop=FALSE]
      stopifnot(nrow(b)==length(rows),!anyDuplicated(b$RID))
      blocks <- convert_ft_bed12_to_bed6(b,format='bed12_fibertools',source=path)
      blocks <- blocks[blocks$end>=r$analysis_start & blocks$end<=r$analysis_end,,drop=FALSE]
      j <- blocks$end-r$analysis_start+1L
      if (isTRUE(result$params$orient_strand) && r$tss_strand=='-') j <- 1001L-j
      hits[[length(hits)+1L]] <- data.frame(i=rows[match(as.character(blocks$RID),a$original_RID[rows])],j=j)
    }
    if (i%%500==0) message('Per-bp m6A: ',i,'/',nrow(used),' contributing loci')
  }
  h <- unique(do.call(rbind,hits))
  site <- Matrix::sparseMatrix(i=h$i,j=h$j,x=1,dims=c(nrow(a),1000L),
    dimnames=list(a$RID,as.character(-500:499)))
  stopifnot(isTRUE(all.equal(bin_read_matrix(as.matrix(site),-500L,499L,25L)$mat,
    result$feat_mat[a$RID,,drop=FALSE],check.attributes=FALSE)))
  site
}

########
# GRO-cap signal, TSS marks and nucleosome occupancy for the sampled windows
########

# 1. For each window with sampled reads, read the plus- and minus-strand
#    GRO-cap signal and add it to that class's running total.
# 2. Record the window's own major TSS(s) and every PINTS bidirectional-peak
#    major TSS that falls inside it.
# 3. Count, per bp, how many bidirectional peaks cover the window.
# 4. Read the fibertools nucleosome calls of the sampled reads (all sizes)
#    and mark the bases each nucleosome covers, one row per read.
# 5. When orienting, flip minus-strand windows: every track is reversed,
#    offsets change sign and the GRO-cap strands swap (+ = sense).
# 6. Return mean GRO-cap / peak profiles per class, the TSS marks, the
#    reads x bp nucleosome matrix and the number of windows per class.
tss_annotations <- function(regions, samples, assignments, cfg) {
  used <- regions[regions$region_id %in% assignments$region_id, ]
  inputs <- fread(cfg$inputs, data.table=FALSE)
  bw <- setNames(lapply(c('+','-'),function(s) rtracklayer::BigWigFile(inputs$source_path[inputs$strand==s])),c('+','-'))
  bid <- fread(file.path(cfg$pints,'GM12878_GROcap_GRCh38_1_bidirectional_peaks.bed'),header=FALSE,data.table=FALSE)
  bg <- GenomicRanges::GRanges(bid$V1,IRanges::IRanges(bid$V2+1L,bid$V3))
  sums <- list(); markers <- list()
  nuc <- matrix(0,nrow(assignments),1000,dimnames=list(assignments$RID,as.character(-500:499)))
  for (i in seq_len(nrow(used))) {
    r <- used[i,,drop=FALSE]; cl <- r$class
    flip <- isTRUE(cfg$orient_strand) && r$tss_strand=='-'
    o <- function(v) if (flip) rev(v) else v
    sign <- if (flip) -1L else 1L
    if (is.null(sums[[cl]])) sums[[cl]] <- list(n=0,plus=numeric(1000),minus=numeric(1000),peaks=numeric(1000))
    sums[[cl]]$n <- sums[[cl]]$n+1
    query <- GenomicRanges::GRanges(r$chr,IRanges::IRanges(r$analysis_start,r$analysis_end))
    for (strand in c('+','-')) {
      signal <- rtracklayer::import(bw[[strand]],which=query,as='GRanges')
      v <- numeric(1000)
      if (length(signal)) for (j in seq_along(signal)) {
        lo <- max(1L, BiocGenerics::start(signal)[j]-r$start)
        hi <- min(1000L,BiocGenerics::end(signal)[j]-r$start)
        v[lo:hi] <- abs(signal$score[j])
      }
      key <- if (strand=='+') 'plus' else 'minus'
      shown <- if (flip) setdiff(c('+','-'),strand) else strand
      shown_key <- if (shown=='+') 'plus' else 'minus'
      sums[[cl]][[shown_key]] <- sums[[cl]][[shown_key]]+o(v)
      major <- r[[key]]
      if (!is.na(major)) markers[[length(markers)+1L]] <- data.frame(class=cl,strand=shown,pos=sign*(major-r$center),kind='Focal major TSS',region_id=r$region_id)
    }
    hits <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(query,bg))
    v <- numeric(1000)
    for (j in hits) {
      lo <- max(1L,bid$V2[j]-r$start+1L); hi <- min(1000L,bid$V3[j]-r$start)
      v[lo:hi] <- v[lo:hi]+1
      for (strand in c('+','-')) {
        pos <- sign*(bid[j,if(strand=='+') 5 else 6]-r$center)
        shown <- if (flip) setdiff(c('+','-'),strand) else strand
        if (pos>= -500 && pos<500) markers[[length(markers)+1L]] <- data.frame(class=cl,strand=shown,pos=pos,kind='All bidirectional major TSSs',region_id=r$region_id)
      }
    }
    sums[[cl]]$peaks <- sums[[cl]]$peaks+o(v)
    ix <- which(assignments$region_id==r$region_id)
    a <- assignments[ix,,drop=FALSE]; a$RID <- a$molecule_id
    ss <- samples[samples$sample_name %in% a$sample_name,,drop=FALSE]
    fp <- extract_nucleosomes(ss,r,a,min_size=NULL,max_size=NULL)
    if (nrow(fp)) for (j in seq_len(nrow(fp))) {
      row <- ix[match(fp$RID[j],a$RID)]
      lo <- max(1L,fp$start[j]-r$start+1L); hi <- min(1000L,fp$end[j]-r$start)
      if(lo<=hi) nuc[row,if (flip) 1001L-(lo:hi) else lo:hi] <- 1
    }
    if (i%%100==0) message('Annotated ',i,'/',nrow(used),' contributing loci')
  }
  profile <- dplyr::bind_rows(lapply(names(sums),function(cl) {
    x <- sums[[cl]]
    dplyr::bind_rows(lapply(c('plus','minus','peaks'),function(track)
      data.frame(class=cl,pos=-500:499,track=track,value=x[[track]]/x$n)))
  }))
  list(profiles=profile,markers=dplyr::bind_rows(markers),nucleosomes=nuc,
       n_loci=vapply(sums,`[[`,numeric(1),'n'))
}

########
# Draw and save all figures as PDFs
########

# 1. pooled_tss_grocap.pdf: TSS marks per 25 bp, bidirectional-peak coverage
#    and mean GRO-cap signal, divergent vs unidirectional.
# 2. cluster_profiles_composition.pdf: one row per cluster with the per-bp m6A
#    (cluster colour) over nucleosome occupancy (grey), plus bars of the
#    sample and TSS-class mix in that cluster (Leiden profile layout).
# 3. cluster_composition.pdf: fraction of each class's reads in each cluster.
# 4. sampled_read_heatmaps.pdf: Leiden per-bp m6A heatmap of all sampled
#    reads, both classes together, split and labelled by cluster, with
#    per-cluster m6A profiles on top.
# 5. Return the per-cluster profile table.
tss_plots <- function(result, ann, composition, output) {
  colors <- cluster_id_palette(levels(result$assignments$cluster))
  x_label <- if (isTRUE(result$params$orient_strand))
    'Position relative to TSS / pair centre (bp); minus-strand TSSs flipped' else
    'Position relative to window centre (bp)'
  strand_colors <- c('+'='#0072B2','-'='#D55E00')
  theme_set(theme_bw(base_size=11))
  guide <- geom_vline(xintercept=0,linetype=3,color='grey50')
  marker <- ann$markers
  marker$bin <- floor((marker$pos+500)/25)*25-487.5
  m <- as.data.table(marker)[,.(value=.N),by=.(class,strand,kind,bin)]
  m$value <- m$value/ann$n_loci[m$class]
  p1 <- ggplot(m,aes(bin,value,color=strand))+geom_line()+guide+
    facet_grid(kind~class)+scale_color_manual(values=strand_colors)+
    labs(x=NULL,y='Major TSSs / locus / 25 bp',color='Strand')+coord_cartesian(xlim=c(-500,500))
  p2 <- ggplot(subset(ann$profiles,track=='peaks'),aes(pos,value))+geom_area(fill='grey65')+
    guide+facet_wrap(~class,nrow=1)+labs(x=NULL,y='Bidirectional intervals / locus')
  signals <- subset(ann$profiles,track!='peaks'); signals$strand <- ifelse(signals$track=='plus','+','-')
  p3 <- ggplot(signals,aes(pos,value,color=strand))+geom_line()+guide+
    facet_grid(strand~class)+scale_color_manual(values=strand_colors)+
    labs(x=x_label,y='Mean GRO-cap signal',color='Strand')
  ggplot2::ggsave(file.path(output,'pooled_tss_grocap.pdf'),cowplot::plot_grid(p1,p2,p3,ncol=1,rel_heights=c(2,1,2)),width=12,height=12)
  a <- result$assignments
  # Leiden plot_fiberseq_profiles() layout, both classes pooled: per-bp m6A
  # (cluster colour) over nucleosome occupancy, then sample / TSS-class bars.
  met <- cluster_site_profiles(list(assignments=a,region=data.frame(analysis_start=-500L,analysis_end=499L)),
    result$site_met_mat)
  sample_colors <- sample_palette(a$sample_name)
  class_colors <- fiberseq_category_palette(a$class)
  prof <- dplyr::bind_rows(lapply(levels(a$cluster),function(cl) rbind(
    data.frame(cluster=cl,pos=-500:499,track='m6A',fraction=met$met[met$cluster==cl]),
    data.frame(cluster=cl,pos=-500:499,track='nucleosome',fraction=colMeans(ann$nucleosomes[a$RID[a$cluster==cl],,drop=FALSE])))))
  rows <- lapply(levels(a$cluster),function(cl) {
    ix <- which(a$cluster==cl); d <- prof[prof$cluster==cl,]
    nuc <- d[d$track=='nucleosome',]
    p <- plot_group_profile(d,style='ribbon',mapping=aes(pos,fraction),
      layers=list(list(geom='ribbon',data=nuc,mapping=aes(ymin=0,ymax=fraction),fill='#808080',alpha=.18,color=NA),
        list(geom='line',data=nuc,color='#666666',linewidth=.55),
        list(geom='line',data=d[d$track=='m6A',],color=colors[[cl]],linewidth=.45)))+
      geom_vline(xintercept=0,color='grey40',linetype='dashed',linewidth=.35)+
      scale_x_continuous(limits=c(-500,499),expand=expansion(mult=0))+
      scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+
      labs(title=paste0(cl,' (n = ',length(ix),')'),x=x_label,y='Read fraction')+
      theme_bw(base_size=9)+theme(legend.position='none',plot.title=element_text(color=colors[[cl]],face='bold'))
    ########
    # One horizontal stacked bar: category mix within this cluster
    ########

    # 1. Count the reads of each category (sample or TSS class).
    # 2. Turn the counts into fractions of the cluster.
    # 3. Draw them as one 0-100% bar in the Leiden composition style.
    bar <- function(values,cols,heading) {
      tab <- as.data.frame(table(category=factor(values,levels=names(cols))))
      tab$fraction <- tab$Freq/sum(tab$Freq); tab$bar <- 1
      plot_stacked_proportion(tab,'bar','category',cols,y='fraction',position='stack',reverse=TRUE,
        horizontal=TRUE,width=.6,legend=list(drop=FALSE),
        scales=list(scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),labels=scales::percent,expand=expansion(mult=0)),
          scale_x_continuous(breaks=NULL)),
        labels=list(title=heading,subtitle=paste0('n = ',length(ix)),x=NULL,y='Within cluster'),
        theme=theme_bw(base_size=9)+theme(legend.position='none',panel.grid=element_blank(),
          plot.margin=margin(5.5,16,5.5,5.5)))
    }
    cowplot::plot_grid(p,bar(sub('_.*$','',a$sample_name[ix]),sample_colors,'LCL sample'),
      bar(a$class[ix],class_colors,'TSS class'),nrow=1,rel_widths=c(2.7,1,1.5),align='h',axis='tb')
  })
  header <- cowplot::ggdraw()+cowplot::draw_label(paste0('Divergent + unidirectional TSS windows pooled (',
    nrow(a),' sampled reads)\nCluster color: m6A calls; grey line / light grey fill: nucleosome occupancy (all ft nuc calls)'),size=10)
  legends <- list(fiberseq_legend(sample_colors,'LCL sample',6),fiberseq_legend(class_colors,'TSS class',3))
  heights <- c(.85,rep(1.35,length(rows)),vapply(legends,fiberseq_legend_height,numeric(1)))
  ggplot2::ggsave(file.path(output,'cluster_profiles_composition.pdf'),
    cowplot::plot_grid(plotlist=c(list(header),rows,legends),ncol=1,rel_heights=heights),width=16,height=sum(heights))
  comp <- as.data.frame(composition); comp$cluster <- factor(comp$cluster,levels(a$cluster))
  p5 <- plot_stacked_proportion(comp,'class','cluster',colors,y='fraction',
    border=list(colour='white',linewidth=.2),
    labels=list(x=NULL,y='Fraction of sampled reads',fill='Cluster',
      inside=list(min_fraction=.03,accuracy=1,vjust=.5,color='white',size=3)),
    scales=list(scale_y_continuous(labels=scales::percent,expand=expansion(mult=c(0,.02)))))
  ggplot2::ggsave(file.path(output,'cluster_composition.pdf'),p5,width=5,height=5)
  # Leiden per-bp m6A heatmap, both classes pooled; rows follow the
  # heatmap's own order (cluster, read start, RID), annotated by cluster only.
  shown <- a[order(a$cluster,a$start,a$RID),,drop=FALSE]
  heatmap <- plot_read_heatmap(list(assignments=a,site_met_mat=result$site_met_mat),
    region=data.frame(chr='pooled',analysis_start=-500L,analysis_end=499L,width=1000L),
    sample_colors=sample_palette(a$sample_name),show_cluster_profiles=TRUE,cluster_colors=colors,
    heatmap_options=list(
      left_annotation=ComplexHeatmap::rowAnnotation(cluster=shown$cluster,col=list(cluster=colors)),
      column_title=paste0('Divergent + unidirectional TSS windows pooled (',nrow(a),' sampled reads)\n',
        'Top: m6A at observed sites per cluster; filled, no smoothing\n',x_label)))
  stopifnot(identical(rownames(heatmap@matrix),shown$RID))
  save_figure(heatmap,file.path(output,'sampled_read_heatmaps.pdf'),width=14,
    height=11+0.55*nlevels(a$cluster),draw=function(x) ComplexHeatmap::draw(x,newpage=FALSE))
  invisible(prof)
}

########
# Quick checks of binning and sampling before the real run
########

# 1. Make a toy 3-read matrix and check the 25-bp binning gives the expected
#    values; check the divergent centre rounds down.
# 2. Run the random sampler twice with the same seed and check it returns the
#    same reads, no duplicates, and counts every read seen.
# 3. Check that flipping a minus-strand read reverses it, leaves plus-strand
#    reads alone, and gives the same features whether flipped before or after
#    binning.
tss_self_test <- function() {
  m <- matrix(0,3,1000,dimnames=list(paste0('r',1:3),as.character(-500:499)))
  m[1,1:25]<-1; m[2,1000]<-1
  f <- bin_read_matrix(m,-500,499,25)$mat
  stopifnot(ncol(f)==40,f[1,1]==1,f[2,40]==1/25,all(f[3,]==0),floor((10+13)/2)==11)
  md <- data.frame(RID=rownames(m))
  set.seed(1); x <- reservoir_add(NULL,f,md,2)
  set.seed(1); y <- reservoir_add(NULL,f,md,2)
  stopifnot(identical(x,y),nrow(x$meta)==2,!anyDuplicated(x$meta$RID),x$seen==3)
  s <- c('+','-','+'); o <- orient_rows_by_strand(m,s)
  stopifnot(o[2,1]==1,o[2,1000]==0,identical(o[-2,],m[-2,]),
    isTRUE(all.equal(bin_read_matrix(o,-500,499,25)$mat,orient_rows_by_strand(f,s))))
  message('Coordinate, binning, orientation and reproducible reservoir checks passed')
}

########
# Run the whole pipeline: sample, cluster, annotate, save, plot
########

# 1. Load libcurl for rtracklayer and run the self-test.
# 2. Build the windows and read the sample table; --pilot keeps 6 windows
#    per class and 20 reads per class.
# 3. Sample the reads (tss_sample()) and check 5,000 per class, unique IDs.
# 4. Cluster both classes together: consensus Leiden on the 25-bp features,
#    minus-strand TSS reads flipped first (orient_strand).
# 5. Copy the read details into the assignments; if consensus did not
#    converge, save clustering.rds and stop.
# 6. Add per-bp m6A (tss_site_m6a()) and the GRO-cap / nucleosome
#    annotations (tss_annotations()); count reads per class x cluster.
# 7. Pilot: draw the figures into a temporary folder, then delete it.
#    Full run: save clustering.rds, read_assignments.tsv and
#    cluster_composition.tsv in tables/, figures in plots/.
tss_main <- function(pilot=FALSE,cfg=tss_config()) {
  # Match the existing PolII workflow's libcurl setup for rtracklayer.
  dyn.load('/project/spott/cshan/envs/dimelo/lib/libcurl.so.4', local=FALSE)
  stopifnot(requireNamespace('rtracklayer',quietly=TRUE))
  data.table::setDTthreads(1L); tss_self_test()
  regions <- tss_regions(cfg); samples <- read.csv(cfg$samples,stringsAsFactors=FALSE)
  if(pilot) {cfg$n<-20L; regions<-dplyr::bind_rows(lapply(split(regions,regions$class),head,6L))}
  dat <- tss_sample(regions,samples,cfg)
  stopifnot(all(table(dat$meta$class)==cfg$n),!anyDuplicated(dat$meta$RID),
    identical(rownames(dat$feature),dat$meta$RID))
  # Already-binned features: window_size=0 prevents a second binning step.
  result <- leiden_manhattan_consensus_cluster(dat$feature,dat$meta,-500L,499L,
    window_size=0,k_neighbors=cfg$k,resolution=cfg$resolution,seed=cfg$seed,
    consensus=TRUE,n_partitions=cfg$partitions,tau=cfg$tau,max_rounds=cfg$max_rounds,
    orient_strand=cfg$orient_strand,strand_col='tss_strand')
  message('Consensus converged: ',result$params$converged,'; rounds: ',result$params$n_rounds)
  for (key in setdiff(names(dat$meta),names(result$assignments))) result$assignments[[key]]<-dat$meta[[key]]
  result$config<-cfg; result$regions<-regions[regions$region_id %in% dat$meta$region_id,]
  result$eligible_observations<-dat$eligible; result$params$feature_bin_bp<-25L
  tables <- file.path(cfg$output,'tables'); plots <- file.path(cfg$output,'plots')
  if (!isTRUE(result$params$converged)) {
    if(!pilot) {dir.create(tables,recursive=TRUE,showWarnings=FALSE); saveRDS(result,file.path(tables,'clustering.rds'))}
    stop('Consensus did not converge; final figures not generated')
  }
  result$site_met_mat <- tss_site_m6a(regions,samples,result)
  ann <- tss_annotations(regions,samples,result$assignments,cfg)
  composition <- as.data.table(result$assignments)[,.N,by=.(class,cluster)]
  composition[,fraction:=N/sum(N),by=class]
  result$annotations<-ann; result$composition<-composition
  if(pilot) {
    tmp<-tempfile('tss-pilot-');dir.create(tmp);on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
    tss_plots(result,ann,composition,tmp)
    message('PILOT PASSED: extraction, consensus, annotations and all figures; temporary outputs removed')
  } else {
    dir.create(tables,recursive=TRUE,showWarnings=FALSE); dir.create(plots,showWarnings=FALSE)
    saveRDS(result,file.path(tables,'clustering.rds'))
    fwrite(result$assignments,file.path(tables,'read_assignments.tsv'),sep='\t')
    fwrite(composition,file.path(tables,'cluster_composition.tsv'),sep='\t')
    tss_plots(result,ann,composition,plots)
    message('COMPLETE: ',cfg$output)
  }
  invisible(result)
}
if (sys.nframe()==0L) tss_main(pilot='--pilot' %in% commandArgs(trailingOnly=TRUE))
