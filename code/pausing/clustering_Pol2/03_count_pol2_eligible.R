# Count footprint-positive full-window molecules independently of m6A denominators.
project <- "/project/spott/cshan/fiber-seq"
here <- file.path(project, "code/pausing/clustering_Pol2")
source(file.path(here, "pol2_chip_cluster_utils", "pol2_clustering_functions.R"))
pol2_load_helpers(project)
if (!exists("cfg", inherits = FALSE)) cfg <- yaml::read_yaml(file.path(here, "pol2_chip_cluster_utils", "config.yaml"))
prep <- pol2_read_prepared(cfg)
tf_paths <- pol2_tf_path(prep$samples, cfg$footprint_dir)
count_chrom <- function(chromosome) {
  sites <- prep$tss[chrom == chromosome]
  if (!nrow(sites)) return(data.table())
  windows <- pol2_gr(sites$chrom, sites$window_start, sites$window_end)
  query <- GenomicRanges::reduce(windows)
  paths <- pol2_ft_path(prep$samples, chromosome)
  pol2_required(c(paths, paste0(paths, ".tbi")))
  matches <- lapply(seq_len(nrow(prep$samples)), function(i) {
    fp <- pol2_read_footprints(tf_paths[i], query)
    if (!nrow(fp)) return(NULL)
    fp[, sample_name := prep$samples$sample_name[i]]
    candidate <- unique(pol2_candidates(sites, fp)[, .(tss_id, original_RID)])
    if (!nrow(candidate)) return(NULL)
    raw <- unlist(Rsamtools::scanTabix(Rsamtools::TabixFile(paths[i]), param = query), use.names = FALSE)
    if (!length(raw)) return(NULL)
    # Only alignment spans/IDs are needed for eligibility, not the long call lists.
    first_four <- unique(sub("^([^\t]*\t[^\t]*\t[^\t]*\t[^\t]*).*", "\\1", raw, perl = TRUE))
    rm(raw)
    bed <- fread(text = paste(first_four, collapse = "\n"), header = FALSE,
                 col.names = c("chr", "start", "end", "RID"))
    bed <- bed[RID %in% candidate$original_RID]
    if (!nrow(bed)) return(NULL)
    bed <- as.data.table(keep_longest_alignment(as.data.frame(bed)))
    hit <- findOverlaps(windows, pol2_gr(bed$chr, bed$start, bed$end), type = "within", ignore.strand = TRUE)
    if (!length(hit)) return(NULL)
    out <- data.table(tss_id = sites$tss_id[queryHits(hit)], original_RID = bed$RID[subjectHits(hit)])
    out <- merge(out, candidate, by = c("tss_id", "original_RID"))
    out[, `:=`(sample_name = prep$samples$sample_name[i], cell_line = prep$samples$cell_line[i])]
    out
  })
  out <- unique(rbindlist(matches))
  if (!nrow(out)) return(data.table())
  out[, RID := paste(sample_name, original_RID, sep = "::")]
  message(chromosome, ": ", uniqueN(out$RID), " unique eligible molecules")
  out
}
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
parts <- parallel::mclapply(cfg$chromosomes, count_chrom, mc.cores = cores, mc.preschedule = FALSE)
if (any(vapply(parts, inherits, logical(1), "try-error"))) stop("A chromosome count failed")
eligible <- unique(rbindlist(parts))
pol2_save(list(signature = prep$signature, molecules = eligible,
              n_molecules = if (nrow(eligible)) uniqueN(eligible$RID) else 0L),
          file.path(cfg$output, "tables", "eligible_molecules.rds"))
cat("TOTAL UNIQUE ELIGIBLE MOLECULES:", uniqueN(eligible$RID), "\n",
    "ELIGIBLE MOLECULE-TSS PAIRS:", nrow(eligible), "\n",
    "CONTRIBUTING SAMPLES:", uniqueN(eligible$sample_name), "of", nrow(prep$samples), "\n",
    "TSSs WITH ELIGIBLE MOLECULES:", uniqueN(eligible$tss_id), "\n")
