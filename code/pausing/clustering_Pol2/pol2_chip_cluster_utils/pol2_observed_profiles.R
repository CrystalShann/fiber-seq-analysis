# Display adapters only: preserve clustering matrices and shared plotting functions.
# Observed positions are pooled over ALL full-window reads at each TSS, before
# footprint eligibility or clustering. Keep a reusable, independently checked mask.
pol2_observed_masks <- function(cfg, x, build = FALSE, genes = NULL) {
  tables <- file.path(cfg$output,"tables")
  prep <- readRDS(file.path(tables,"prepared.rds"))
  stopifnot(identical(prep$signature,x$signature))
  example_genes <- c(cfg$locus_genes,pol2_select_example_promoters(cfg)$gene_name)
  sites <- if (is.null(genes)) prep$tss[tss_id %in% x$molecules$tss_id | gene_name %in% example_genes] else prep$tss[gene_name %in% genes]
  if (!nrow(sites) || (!is.null(genes) && !all(genes %in% sites$gene_name))) stop("Unknown example promoter")
  directory <- file.path(tables,"observed_site_profiles")
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  worker <- function(cc) {
    ss <- sites[chrom==cc]
    paths <- pol2_ft_path(prep$samples,cc)
    key <- list(sites=ss[,.(tss_id,tss,strand,window_start,window_end)],
      inputs=pol2_stamp(c(paths,paste0(paths,".tbi"))),version=1L)
    path <- file.path(directory,paste0(cc,"_mask.rds"))
    if (file.exists(path)) {
      old <- readRDS(path)
      if (isTRUE(all.equal(old$key,key,tolerance=0))) return(old$mask)
      if (!is.null(genes) && identical(old$key$inputs,key$inputs) &&
          all(ss$tss_id %in% rownames(old$mask)) &&
          isTRUE(all.equal(old$key$sites[match(ss$tss_id,old$key$sites$tss_id)],key$sites,tolerance=0)))
        return(old$mask[ss$tss_id,,drop=FALSE])
    }
    if (!is.null(genes)) {
      path <- file.path(directory,paste0(cc,"_",paste(sort(genes),collapse="_"),"_mask.rds"))
      if (file.exists(path)) {
        old <- readRDS(path)
        if (isTRUE(all.equal(old$key,key,tolerance=0))) return(old$mask)
      }
    }
    if (!build) stop("Missing or stale observed-site cache for ",cc,
      ". Run: sbatch --cpus-per-task=4 code/pausing/clustering_Pol2/00_run.sh --stage observed-sites",call.=FALSE)
    mask <- matrix(FALSE,nrow(ss),2000,dimnames=list(ss$tss_id,as.character(-1000:999)))
    # Query reads in the same 100-TSS batches as pol2_features(): each read's longest
    # alignment is chosen within one batch query, so other batches change n_positions.
    all_sites <- prep$tss[chrom==cc]
    for (idx in split(seq_len(nrow(all_sites)),ceiling(seq_len(nrow(all_sites))/100L))) {
      block <- all_sites[idx]
      rows <- which(ss$tss_id %in% block$tss_id)
      if (!length(rows)) next
      regions <- GenomicRanges::reduce(pol2_gr(block$chrom,block$window_start,block$window_end))
      for (path_in in paths) {
        bed <- as.data.table(read_ft_bed12(path_in,regions,longest_alignment=TRUE))
        if (!nrow(bed)) next
        for (j in rows) {
          site <- ss[j]
          pool <- bed[start <= site$window_start & end >= site$window_end]
          if (!nrow(pool)) next
          calls <- as.data.table(convert_ft_bed12_to_bed6(as.data.frame(pool),format="bed12_fibertools"))
          if (!nrow(calls)) next
          stopifnot(all(calls$end-calls$start==1L))
          pos <- unique((calls[end > site$window_start & end <= site$window_end,end]-site$tss)*ifelse(site$strand=="+",1L,-1L))
          stopifnot(all(pos %in% -1000:999))
          mask[j,pos+1001L] <- TRUE
        }
      }
      message(cc,": recovered observed sites for ",max(idx),"/",nrow(all_sites)," promoters")
    }
    expected <- unique(x$molecules[,.(tss_id,n_positions)])
    ii <- match(rownames(mask),expected$tss_id)
    stopifnot(all(rowSums(mask)[!is.na(ii)]==expected$n_positions[ii[!is.na(ii)]]))
    temporary <- tempfile(tmpdir=directory)
    saveRDS(list(key=key,mask=mask),temporary)
    if (!file.rename(temporary,path)) stop("Could not save mask cache")
    mask
  }
  cores <- max(1L,min(4L,as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","4"))))
  chroms <- unique(sites$chrom)
  result <- if (build) parallel::mclapply(chroms,worker,mc.cores=cores,mc.preschedule=FALSE) else
    lapply(chroms,worker)
  failed <- vapply(result,inherits,logical(1),"try-error")
  if (any(failed)) stop("Observed-site extraction failed for ",paste(chroms[failed],collapse=", "),": ",
    paste(unique(vapply(result[failed],as.character,character(1))),collapse="; "))
  do.call(rbind,result)
}

# Fractions at each bp: sum(calls) / reads whose TSS has this observed position.
# Zero-denominator positions are missing, never artificial zero-accessibility sites.
pol2_observed_profile <- function(calls, assignments, masks, tss_ids) {
  calls <- as.matrix(calls[assignments$RID,,drop=FALSE])
  opportunities <- masks[tss_ids,,drop=FALSE]
  stopifnot(nrow(calls)==nrow(opportunities),identical(colnames(calls),colnames(opportunities)),
    all(calls[!opportunities]==0))
  rbindlist(lapply(levels(assignments$cluster),function(cl) {
    ii <- which(assignments$cluster==cl)
    denominator <- colSums(opportunities[ii,,drop=FALSE])
    numerator <- colSums(calls[ii,,drop=FALSE])
    data.table(cluster=cl,pos=as.integer(colnames(calls)),n_calls=numerator,n_observed_reads=denominator,
      met=ifelse(denominator>0,numerator/denominator,NA_real_))
  }))
}

pol2_observed_heatmap <- function(res, sample_colors, profile) {
  groups <- levels(res$assignments$cluster)
  colors <- cluster_id_palette(groups)
  top <- lapply(groups,function(cl) methylation_area_annotation(
    profile[cluster==cl][match(-1000:999,pos),met],colors[[cl]]))
  names(top) <- paste0(groups," m6A")
  top <- do.call(ComplexHeatmap::HeatmapAnnotation,c(top,list(annotation_name_gp=grid::gpar(fontsize=8),gap=grid::unit(1.5,"mm"))))
  plot_read_heatmap(res,pol2_plot_region(),sample_colors,show_cluster_profiles=FALSE,
    heatmap_options=list(top_annotation=top,column_title=paste0(
      "Pol II clustering: ",nrow(res$assignments)," molecules\nTop: m6A at observed sites; no smoothing\n",POL2_X_LABEL)))
}

pol2_observed_loci <- function(cfg,x,annotation,masks) {
  prep <- readRDS(file.path(cfg$output,"tables/prepared.rds"))
  stopifnot(identical(prep$signature,x$signature))
  models <- load_canonical_gene_models(cfg$gene_models)
  eqtl <- if (is.null(cfg$eqtl_pairs)) NULL else fread(cfg$eqtl_pairs,select=c("variant_id","gene_id","pval_nominal"))
  # Private adapters reuse original layout and saved-input reader without changing
  # any shared function. This is a display-only render of the saved analysis.
  reader <- pol2_locus_result
  reader_env <- new.env(parent=environment(reader))
  reader_env$pol2_read_prepared <- function(config) prep
  environment(reader) <- reader_env
  setNames(lapply(cfg$locus_genes,function(gene) {
    site <- prep$tss[gene_name==gene]
    renderer <- plot_locus_molecules
    renderer_env <- new.env(parent=environment(renderer))
    renderer_env$cluster_site_profiles <- function(res,met_mat) {
      out <- as.data.table(cluster_site_profiles(res,met_mat))
      relative <- (out$pos-site$tss)*ifelse(site$strand=="+",1L,-1L)
      keep <- masks[site$tss_id,match(relative,as.integer(colnames(masks)))]
      out[which(keep)]
    }
    environment(renderer) <- renderer_env
    figure <- pol2_locus_figure
    figure_env <- new.env(parent=environment(figure))
    figure_env$pol2_locus_result <- function(...) {
      res <- reader(...)
      relative <- as.integer(colnames(masks))[masks[site$tss_id,]]
      res$observed_positions <- site$tss + ifelse(site$strand=="+",1L,-1L)*relative
      res
    }
    figure_env$plot_locus_molecules <- renderer
    environment(figure) <- figure_env
    figure(cfg,x,annotation,gene,models,eqtl)
  }),cfg$locus_genes)
}
