#!/usr/bin/env Rscript
# Run the full active/inactive comparison from the saved footprint input snapshot.
run_enhancer_site_pairs <- function() {
  suppressPackageStartupMessages({library(data.table);library(GenomicRanges);library(Rsamtools)})
  setDTthreads(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","4")))
  project_root <- Sys.getenv("ENHANCER_PROJECT_ROOT","/project/spott/cshan/fiber-seq")
  output_dir <- Sys.getenv("ENHANCER_COOCC_OUTPUT_DIR",file.path(project_root,"macrophage_project/enhancer/TF_co-occ"))
  table_dir <- file.path(output_dir,"tables"); plot_dir <- file.path(output_dir,"plots")
  dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
  dir.create(plot_dir,recursive=TRUE,showWarnings=FALSE)
  input <- Sys.getenv("ENHANCER_FOOTPRINT_RESULTS",file.path(table_dir,"enhancer_footprint_results.rds"))
  input_md5 <- unname(tools::md5sum(input))
  footprint_result <- readRDS(input)
  analysis_enhancers <- as.data.table(footprint_result$inputs$enhancers)
  sample_table <- as.data.table(footprint_result$inputs$samples)
  footprint_counts <- as.data.table(footprint_result$footprint_counts)
  extract_root <- footprint_result$inputs$extract_root
  ft_result_dir <- footprint_result$inputs$ft_result_dir
  minimum_reads <- as.integer(Sys.getenv("ENHANCER_PAIR_MIN_READS","10"))
  min_pair_reads <- minimum_reads
  site_window_bp <- as.integer(Sys.getenv("ENHANCER_PAIR_WINDOW_BP","50"))
  min_site_observations <- as.integer(Sys.getenv("ENHANCER_PAIR_MIN_SITE_OBSERVATIONS","5"))
  stopifnot(minimum_reads>=1L,site_window_bp>=1L,min_site_observations>=1L)
  eligible_enhancers <- footprint_counts[,.(keep=.N==nrow(sample_table) &&
    all(n_reads>=minimum_reads)),by=enhancer_id][keep==TRUE,enhancer_id]
  hmm_file <- function(sample,chr,kind) {
    compact <- gsub("_","",sample,fixed=TRUE)
    file.path(extract_root,paste0("firehmm_",kind),compact,
              paste0(compact,"_hmm_extracted_",kind,"_",chr,".bed.gz"))
  }
  source(file.path(project_root,"code/enhancer/co_occ_functions/run_enhancer_footprints.R"),local=TRUE)
  source(file.path(project_root,"code/enhancer/co_occ_functions/enhancer_site_pair_functions.R"),local=TRUE)
  source(file.path(project_root,"code/enhancer/co_occ_functions/enhancer_site_pair_plots.R"),local=TRUE)
  message("Comparing ",length(eligible_enhancers)," eligible enhancers across both classes")
  result <- calculate_enhancer_site_pairs(analysis_enhancers,sample_table,footprint_counts,
    footprint_result,eligible_enhancers,site_window_bp,min_site_observations,min_pair_reads,
    cobinding_alpha=0.01,cobinding_significance="q_value")
  result$comparisons <- compare_enhancer_pair_counts(result$enhancers)
  result$inputs <- footprint_result$inputs
  result$input_rds <- normalizePath(input); result$input_md5 <- input_md5
  result$completed_at <- Sys.time(); result$job_id <- Sys.getenv("SLURM_JOB_ID","interactive")
  stopifnot(identical(input_md5,unname(tools::md5sum(input))))
  write_table <- function(x,name) {
    temp <- tempfile(".pair_table_",tmpdir=table_dir)
    fwrite(x,temp,sep="\t",na="NA")
    if(!file.rename(temp,file.path(table_dir,name))) stop("Could not publish ",name)
  }
  write_table(result$pairs,"05_site_pair_cooccupancy.tsv")
  write_table(result$enhancers,"05_cobinding_per_enhancer.tsv")
  write_table(result$sites,"05_cobinding_site_inventory.tsv")
  write_table(result$comparisons,"05_active_inactive_cobinding_tests.tsv")
  save_enhancer_cobinding_plots(result,result$comparisons,plot_dir)
  temp <- tempfile(".pair_results_",tmpdir=table_dir)
  saveRDS(result,temp)
  stopifnot(file.rename(temp,file.path(table_dir,"enhancer_site_pair_results.rds")))
  print(result$comparisons[metric=="standardized_pairs"])
  message("Completed active/inactive site-pair analysis at ",Sys.time())
  invisible(result)
}
if (sys.nframe()==0L) run_enhancer_site_pairs()
