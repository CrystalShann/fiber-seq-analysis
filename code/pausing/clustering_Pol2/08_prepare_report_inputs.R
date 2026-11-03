# Preparation stages only. Sourced by 01_run_pol2_clustering.R.
if (opt$stage == "motif-inputs") {
  x <- readRDS(file.path(cfg$output,"tables/clustering.rds"))
  example_gene <- opt$gene
  tables <- file.path(cfg$output, "tables")
  tabdir <- file.path(tables, "promoter_motif_cluster_examples")
  shared_read_cache <- file.path(tabdir,"promoter_read_data.rds")
  if (!is.null(example_gene)) tabdir <- file.path(tabdir,example_gene)
  plotdir <- file.path(cfg$output, "plots", "promoter_motif_cluster_examples")
  dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)
  selected <- pol2_select_example_promoters(cfg)
  if (!is.null(example_gene)) {
    selected <- selected[gene_name==example_gene]
    if (nrow(selected)!=1L) stop("--gene must be one of the seven selected example promoters")
  }
  prep <- readRDS(file.path(tables, "prepared.rds"))
  promoters <- prep$tss[match(selected$tss_id, tss_id)]
  stopifnot(identical(promoters$tss_id, selected$tss_id))
  families <- c("TATA","YY1","NFY","SP","ETS","ZNF143","NRF1","CREB")
  # Scan explicit genomic windows with the standalone scanner; no pooled coordinates.
  promoter_path <- file.path(tabdir, "selected_promoter_windows.tsv")
  fwrite(promoters, promoter_path, sep="\t")
  status <- system2(cfg$python, c(shQuote(file.path(here,"05_motif_scan.py")),
    "--promoters",shQuote(promoter_path),"--output",shQuote(tabdir),
    "--reference",shQuote(cfg$reference)))
  stopifnot(status == 0L)
  motif_path <- file.path(tabdir, "matrix_motif_hits.tsv.gz")
  hits <- fread(motif_path)
  stopifnot(all(hits$tss_id %in% selected$tss_id), all(hits$category %in% families))
  hits[, `:=`(width = end-start, xmin = relative_start-.5, xmax = relative_end+.5)]
  fwrite(selected, file.path(tabdir, "selected_promoters.tsv"), sep = "\t")
  fwrite(data.table(category = families, status = "available"),
         file.path(tabdir,"motif_category_status.tsv"),sep="\t")
  # Saved local data avoids repeatedly querying 31 sample files on report renders.
  cache_path <- file.path(tabdir,"promoter_read_data.rds")
  input_paths <- c(file.path(tables,c("prepared.rds","clustering.rds","tracks/molecule_projection.tsv.gz")),
    file.path(tables,"features",paste0(unique(promoters$chrom),".rds")))
  cache_key <- list(tss_id=selected$tss_id, inputs=pol2_stamp(input_paths))
  if (!is.null(example_gene) && !file.exists(cache_path) && file.exists(shared_read_cache)) {
    shared <- readRDS(shared_read_cache)
    stopifnot(all(selected$tss_id %in% names(shared$data)),
      identical(shared$key$inputs,pol2_stamp(shared$key$inputs$path)))
    saveRDS(list(key=cache_key,data=shared$data[selected$tss_id]),cache_path)
  }
  if (file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    stopifnot(identical(cached$key,cache_key))
    local_data <- cached$data
  } else {
    projection <- fread(file.path(tables,"tracks/molecule_projection.tsv.gz"))[tss_id %in% selected$tss_id]
    local_data <- list()
    for (cc in unique(promoters$chrom)) {
      feature <- readRDS(file.path(tables,"features",paste0(cc,".rds")))
      stopifnot(identical(feature$signature$prepared,prep$signature))
      for (tid in promoters[chrom == cc,tss_id]) {
        meta <- copy(feature$molecules[tss_id == tid])
        meta[, cluster := projection$cluster[match(observation_id,projection$observation_id)]]
        meta[, clustered := projection$clustered[match(observation_id,projection$observation_id)]]
        stopifnot(!anyNA(meta$cluster),!anyDuplicated(meta$RID))
        ii <- match(meta$observation_id,x$molecules$observation_id)
        stopifnot(all(meta$cluster[!is.na(ii)] == as.character(x$molecules$cluster[ii[!is.na(ii)]])))
        calls <- feature$sites[meta$observation_id,,drop=FALSE]
        stopifnot(identical(colnames(calls),as.character(-1000:999)))
        local_data[[tid]] <- list(meta=meta,calls=calls)
      }
    }
    all_meta <- rbindlist(lapply(local_data,`[[`,"meta"))
    all_meta <- merge(all_meta,promoters[,.(tss_id,window_start,window_end)],by="tss_id")
    blocks <- list()
    for (sid in unique(all_meta$sample_name)) {
      sm <- prep$samples[sample_name == sid]
      mm <- all_meta[sample_name == sid]
      pp <- promoters[tss_id %in% mm$tss_id]
      regions <- GenomicRanges::reduce(pol2_gr(pp$chrom,pp$window_start,pp$window_end))
      tf <- pol2_read_footprints(pol2_tf_path(sm,cfg$footprint_dir),regions)
      if (nrow(tf)) tf[, kind := "tf"]
      nuc <- rbindlist(lapply(unique(pp$chrom),function(cc) {
        path <- file.path(sm$fire_dir,"extracted_results/nuc_by_chr",paste0(sid,".ft_extracted_nuc.",cc,".bed.gz"))
        pol2_required(c(path,paste0(path,".tbi")))
        raw <- read_tabix_region(path,regions[as.character(GenomicRanges::seqnames(regions)) == cc])
        if (!nrow(raw)) return(NULL)
        raw <- raw[raw[[4]] %in% mm$original_RID,,drop=FALSE]
        if (!nrow(raw)) return(NULL)
        b <- as.data.table(convert_ft_bed12_to_bed6(raw,format="bed12_fibertools",longest_alignment=TRUE,source=path))
        unique(b[,.(chrom=chr,start,end,original_RID=RID,kind="nuc")])
      }))
      b <- rbindlist(list(tf,nuc),fill=TRUE)
      b <- merge(b,mm[,.(original_RID,observation_id,tss_id,tss,strand,mol_chrom=chrom,window_start,window_end)],by="original_RID",allow.cartesian=TRUE)
      b <- b[chrom == mol_chrom & end > window_start & start < window_end]
      b[, length := end-start]
      b[, class := fcase(kind=="nuc" & length>=130 & length<=160,"nucleosome 130-160 bp",
        kind=="tf" & length>=10 & length<30,"TF footprint 10-30 bp",
        kind=="tf" & length>=30 & length<=60,"TF footprint 30-60 bp",
        kind=="tf" & length>60 & length<=80,"TF footprint 60-80 bp")]
      b <- b[!is.na(class)]
      orientation <- ifelse(b$strand=="+",1,-1)
      b[, `:=`(rel_start=pmin(orientation*(start+.5-tss),orientation*(end+.5-tss)),
                rel_end=pmax(orientation*(start+.5-tss),orientation*(end+.5-tss)))]
      blocks[[sid]] <- unique(b[,.(observation_id,tss_id,chrom,start,end,length,class,kind,rel_start,rel_end)])
      message("Loaded promoter footprints: ",sid)
    }
    blocks <- rbindlist(blocks)
    for (tid in names(local_data)) local_data[[tid]]$blocks <- blocks[tss_id == tid]
    saveRDS(list(key=cache_key,data=local_data),cache_path)
  }

  message("Prepared motif scan and read cache for ",paste(selected$gene_name,collapse=", "))
} else if (opt$stage == "histone-inputs") {
  tables <- file.path(cfg$output,"tables")
  assignments <- fread(file.path(tables,"molecule_assignments.tsv.gz"))
  profiles <- pol2_histone_profiles(cfg,assignments)
  temporary <- tempfile(tmpdir=tables,fileext=".tsv.gz")
  fwrite(profiles,temporary,sep="\t")
  stopifnot(file.rename(temporary,file.path(tables,"histone_cluster_profiles.tsv.gz")))
  saveRDS(list(inputs=pol2_stamp(c(file.path(tables,"molecule_assignments.tsv.gz"),unlist(cfg$histone_marks))),
    flank=2000L,bin=50L),file.path(tables,"histone_cluster_profiles_provenance.rds"))
  message("Prepared histone profiles: ",nrow(profiles)," rows")
}
