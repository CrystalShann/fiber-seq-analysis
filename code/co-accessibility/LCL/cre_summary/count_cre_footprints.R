# Per-CRE footprint counts read by LCL_coaccess_cre_summary.Rmd 

#   covered     the fiber's aligned blocks cover every base of the CRE
#   open        one FIRE element of the fiber still covers >= 50% of the CRE after deletions
#   informative open, with a FiberHMM record on the chromosome
#   footprint   10-80 bp FiberHMM block with block score >= 50, inside a FIRE element of the same fiber



suppressPackageStartupMessages({library(data.table); library(GenomicRanges)})
setDTthreads(1L)

source('/project/spott/cshan/fiber-seq/code/parsing_functions/parsing_footprints_functions.r')
root <- Sys.getenv('LCL_COACCESS_ROOT','/project/spott/cshan/fiber-seq/LCL_project/co-accessibility')
staging <- file.path(root,'footprints')
tabix <- '/project/spott/cshan/envs/dimelo/bin/tabix'

# recurrent footprint position: >= 15% of open fibers, >= 5 fibers at the peak
min_frac <- .15; min_fibers <- 5L

##############################
# load CREs and assign coordinates to each CRE by laying them end to end --> virtual axis
##############################
read_cres <- function() {
  cre <- fread(file.path(root,'cres.tsv.gz'),select=c('CRE_ID','chrom','start','end'))
  stopifnot(nrow(cre)>0,!anyDuplicated(cre$CRE_ID),all(cre$end>cre$start))
  # offset creates axis for laying CREs
  cre[,`:=`(width=end-start,offset=cumsum(c(0,head(end-start,-1))))]
}

##############################
# read per chr from bgzipped, tabix-indexed BED file
##############################
read_chrom <- function(path,chrom,contigs,cols) {
  if(!chrom %in% contigs) return(data.table())
  x <- fread(cmd=paste(tabix,shQuote(path),chrom),header=FALSE,sep='\t',colClasses=cols,showProgress=FALSE)
  if(!nrow(x)) data.table() else x
}

##############################
# get aligned width of the query intervals
##############################
aligned_width <- function(q,blocks) {
  m <- merge(q,blocks[,.(rid=RID,bs=start,be=end)],by='rid',allow.cartesian=TRUE)
  w <- m[,.(w=sum(pmax(0,pmin(e,be)-pmax(s,bs)))),by=id]
  ans <- numeric(nrow(q));ans[w$id] <- w$w;ans
}

##############################
# process one sample and one chr to get per CRE chr counts with the # of align footprint bases on virtual axis
##############################

# output a list with 1) table with one row per CRE on this chr 2) IRanges with 
# fiber’s merged footprint run inside one CRE, placed on the axis

# 1. table with one row per CRE on this chr:
# CRE_ID
# n_informative: fibers that cover the CRE, are open there, and have a FiberHMM footprint
# n_footprints: footprint × fiber hits touching the CRE (≥ 1 shared base)
# footprint_bp: footprinted bases inside the CRE, summed over fibers (if one fiber has two footprints covering the same base, 
# that base is counted once)
# fibers_with_fp: informative fibers with at least one footprint in the CRE



count_chrom <- function(s,chrom,cg,contigs) {
  # output table
  rows <- data.table(CRE_ID=cg$CRE_ID,n_informative=0L,n_footprints=0L,footprint_bp=0L,fibers_with_fp=0L)
  empty <- list(rows=rows,global=IRanges())
  
  # reads every alignment record and expand each read into one row per block
  # columns: chr, start, end, RID, score, strand
  bl <- read_chrom(s$blocks_path,chrom,contigs$blocks,list(character=c(4,11,12)))
  if(!nrow(bl)) return(empty)
  bl <- bl[!duplicated(V4)]   
  blocks <- as.data.table(convert_ft_bed12_to_bed6(as.data.frame(bl),drop_sentinels=FALSE))
  
  # BED --> GRanges
  cre_gr <- GRanges(chrom,IRanges(cg$start+1L,cg$end))
  
  # covered: aligned blocks cover every base of the CRE
  h <- findOverlaps(cre_gr,GRanges(chrom,IRanges(blocks$start+1L,blocks$end)))
  q <- queryHits(h);b <- subjectHits(h)
  
  # find every CRE-aligned block pair that overlaps
  cov <- data.table(cre=q,rid=blocks$RID[b],w=pmin(cg$end[q],blocks$end[b])-pmax(cg$start[q],blocks$start[b]))
  cov <- cov[,.(w=sum(w)),by=.(cre,rid)][w==cg$width[cre],.(cre,rid)]
  if(!nrow(cov)) return(empty)
  
  # open: a FIRE element whose aligned bases inside the CRE cover >= 50% of it
  fi <- read_chrom(s$elements_path,chrom,contigs$fire,list(character=4))
  if(!nrow(fi)) return(empty)
  fi <- fi[V4 %chin% cov$rid,.(rid=V4,start=V2,end=V3)]
  
  # reads the FIRE elements, only keep the ones that cover at least one CRE
  h <- findOverlaps(cre_gr,GRanges(chrom,IRanges(fi$start+1L,fi$end)))
  q <- queryHits(h);f <- subjectHits(h)
  
  # find every CRE-FIRE overlap
  cand <- data.table(cre=q,rid=fi$rid[f],s=pmax(cg$start[q],fi$start[f]),e=pmin(cg$end[q],fi$end[f]))
  cand <- cand[e-s>=.5*cg$width[cre]] 
  
  # raw overlap first (bedtools -f 0.5), only keeps a CRE/FIRE pair 
  cand[,id:=.I];cand[,aligned:=aligned_width(cand,blocks)]
  acc <- merge(unique(cand[aligned>0 & aligned>=.5*cg$width[cre],.(cre,rid)]),cov,by=c('cre','rid'))
  
  # informative: open with a FiberHMM called footprint on this chr
  tf <- read_chrom(s$tf_path,chrom,contigs$tf,list(character=c(4,11:15)))
  if(!nrow(tf)) return(empty)
  tf <- tf[V4 %chin% acc$rid]
  inf <- acc[rid %chin% tf$V4]
  if(!nrow(inf)) return(empty)
  n_inf <- inf[,.N,by=cre];rows[n_inf$cre,n_informative:=n_inf$N]
  
  # filter footprints to keep them if they are 10-80bp, and with block score >=50
    # expand FiberHMM BED file into one row per footprint block
  fp <- as.data.table(convert_ft_bed12_to_bed6(as.data.frame(tf),keep_block_scores=TRUE,format='bed15_fiberhmm_tf'))
  fp <- fp[end-start>=10 & end-start<=80 & block_score>=50]
  if(!nrow(fp)) return(list(rows=rows,global=IRanges()))
  fp <- as.data.table(subset_footprints_in_fire(as.data.frame(fp),
    data.frame(RID=fi$rid,start=fi$start,end=fi$end,class='FIRE')))
  if(!nrow(fp)) return(list(rows=rows,global=IRanges()))
  
  # footprints overlapping the CRE on each informative fiber
    # footprints that share >= 1 base with the CRE
  inf[,`:=`(pair=.I,cs=cg$start[cre],ce=cg$end[cre])]
  m <- merge(inf,fp[,.(rid=RID,fs=start,fe=end)],by='rid',allow.cartesian=TRUE)[fe>cs & fs<ce]
  if(!nrow(m)) return(list(rows=rows,global=IRanges()))
  n_fp <- m[,.N,by=cre];rows[n_fp$cre,n_footprints:=n_fp$N]
  
  # merge each fiber's footprints inside the CRE: every pair gets its own segment of a virtual axis
  # (width + 1 apart, so pairs never touch), then the merged runs go to cres.tsv.gz coordinates
  pair_off <- cumsum(c(0L,head(cg$width[inf$cre]+1L,-1L)))
  
  # each footprint is clipped to its CRE, converted to a position relative to the CRE start
  merged <- reduce(IRanges(pair_off[m$pair]+pmax(m$fs,m$cs)-m$cs+1L,pair_off[m$pair]+pmin(m$fe,m$ce)-m$cs))
  
  # map each merged interval to the origin coordinates
  mp <- findInterval(start(merged)-1L,pair_off);mc <- inf$cre[mp]
  rel <- start(merged)-1L-pair_off[mp]
  bp <- data.table(cre=mc,w=width(merged))[,.(w=sum(w)),by=cre];rows[bp$cre,footprint_bp:=as.integer(bp$w)]
  nf <- data.table(cre=inf$cre[unique(mp)])[,.N,by=cre];rows[nf$cre,fibers_with_fp:=nf$N]
  global <- IRanges(cg$offset[mc]+rel+1L,width=width(merged))
  list(rows=rows,global=global)
}


##############################
# loop over each sample and chr and calls count_chrom()
##############################
count_sample <- function(i) {
  samples <- fread(file.path(root,'samples.tsv'))
  stopifnot(i>=1L,i<=nrow(samples))
  s <- samples[i]
  cre <- read_cres()
  contigs <- lapply(c(blocks=s$blocks_path,fire=s$elements_path,tf=s$tf_path),
    function(p) system2(tabix,c('-l',shQuote(p)),stdout=TRUE))
  total <- sum(cre$width);base <- integer(total);rows <- list()
  for(ch in unique(cre$chrom)) {
    t0 <- Sys.time()
    ans <- count_chrom(s,ch,cre[chrom==ch],contigs)
    rows[[ch]] <- ans$rows
    if(length(ans$global)) base <- base+as.integer(coverage(ans$global,width=total))
    message(s$sample_name,' ',ch,': ',sum(ans$rows$n_informative),' open informative fiber-CRE pairs, ',
      round(as.numeric(difftime(Sys.time(),t0,units='mins')),1),' min')
  }
  rows <- rbindlist(rows)
  stopifnot(identical(rows$CRE_ID,cre$CRE_ID),all(rows$fibers_with_fp<=rows$n_informative),
    all(rows$footprint_bp<=rows$n_informative*cre$width))
  dir.create(staging,showWarnings=FALSE)
  out <- file.path(staging,paste0(s$sample_name,'.rds'))
  saveRDS(list(sample_name=s$sample_name,total=total,base=base,rows=rows),paste0(out,'.partial'))
  stopifnot(file.rename(paste0(out,'.partial'),out))
}

write_tsv <- function(x,name) {
  path <- file.path(root,name);tmp <- paste0(path,'.partial')
  fwrite(x,tmp,sep='\t',na='NA',compress=if(grepl('\\.gz$',name)) 'gzip' else 'none')
  stopifnot(file.rename(tmp,path))
}

##############################
# combine the results for all samples and find recurrent footprint positions
##############################


# pools all samples’ per-CRE footprint counts, finds recurrent footprint positions, 
# and writes a per-CRE summary table and a positions table


# find recurrent footprint positions
# 1. At each base of the CRE, count how many fibers have a footprint covering that base
# 2. A base is marked if at least 15% of the CRE’s informative fibers are footprinted there
# 3. Consecutive marked bases form one candidate position. The stretch ends wherever the count drops below 15%.
# Two stretches separated by even a single unmarked base are treated as separate positions
# A stretch is only kept if its peak has at least 5 fibers
# 4. Identify the peak by taking the single base with the most footprinted fibers


combine <- function() {
  
  # Load and pool the per-sample results
  samples <- fread(file.path(root,'samples.tsv'))
  files <- file.path(staging,paste0(samples$sample_name,'.rds'))
  if(any(!file.exists(files))) stop('Missing staging files: ',paste(basename(files[!file.exists(files)]),collapse=', '))
  cre <- read_cres();total <- sum(cre$width)
  base <- integer(total);counts <- NULL
  for(f in files) {
    x <- readRDS(f)
    stopifnot(x$total==total,identical(x$rows$CRE_ID,cre$CRE_ID))
    base <- base+x$base
    m <- as.matrix(x$rows[,-1L]);counts <- if(is.null(counts)) m else counts+m
  }
  pooled <- data.table(CRE_ID=cre$CRE_ID,counts)
  # recurrent positions: runs of bases with >= min_frac of open fibers, >= min_fibers at the peak
  # peak is highest fiber count at a single base
  d <- data.table(i=rep(seq_len(nrow(cre)),cre$width),pos=sequence(cre$width)-1L,n=base,
    open=rep(pooled$n_informative,cre$width))
  d[,hot:=open>0 & n>=min_frac*open]
  d[,run:=rleid(i,hot)]
  pos <- d[hot==TRUE,.(from=min(pos),to=max(pos)+1L,peak_fibers=max(n),open=open[1]),by=.(i,run)]
  pos <- pos[peak_fibers>=min_fibers]
  peak <- d[,.(peak_fibers=max(n)),by=i]
  rm(d);invisible(gc())
  positions <- data.table(CRE_ID=cre$CRE_ID[pos$i],chrom=cre$chrom[pos$i],start=cre$start[pos$i]+pos$from,
    end=cre$start[pos$i]+pos$to,peak_fibers=pos$peak_fibers,n_informative=pos$open,
    peak_pct=100*pos$peak_fibers/pos$open)
  summary <- cbind(cre[,.(CRE_ID,chrom,start,end,width)],pooled[,-1L])
  summary[,`:=`(footprints_per_fiber=n_footprints/n_informative,
    pct_bases_footprinted=100*footprint_bp/(n_informative*width),
    peak_pct=100*peak$peak_fibers[match(seq_len(.N),peak$i)]/n_informative,
    n_positions=tabulate(pos$i,nbins=.N))]
  summary[n_informative==0,`:=`(footprints_per_fiber=NA_real_,pct_bases_footprinted=NA_real_,peak_pct=NA_real_)]
  stopifnot(all(summary$pct_bases_footprinted<=100,na.rm=TRUE),all(summary$peak_pct<=100,na.rm=TRUE))
  write_tsv(summary,'cre_footprint_summary.tsv.gz')
  write_tsv(positions,'cre_footprint_positions.tsv.gz')
  unlink(staging,recursive=TRUE)
  message('Combined ',length(files),' samples: ',nrow(summary),' CREs, ',nrow(positions),' recurrent positions')
}

args <- commandArgs(TRUE)
if(length(args)==2L && args[1]=='count') count_sample(as.integer(args[2])) else
if(length(args)==1L && args[1]=='combine') combine() else
  stop('Usage: Rscript count_cre_footprints.R count <sample_index> | combine')
