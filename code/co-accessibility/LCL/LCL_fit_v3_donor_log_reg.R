# Fits the donor-adjusted logistic regression for every CRE pair and writes the
# tables that LCL_co-access_v3_donor_log_reg.rmd reads (knitting does no fitting).
# Submitted by LCL_run_v3_donor_log_reg.sh:  Rscript LCL_fit_v3_donor_log_reg.R <cores> [chrom]
args <- commandArgs(trailingOnly = TRUE)
params <- list(
  root = "/project/spott/cshan/fiber-seq/LCL_project/co-accessibility",
  code_root = "/project/spott/cshan/fiber-seq/code",
  min_reads_per_donor = 2L,
  cores = as.integer(args[1]),
  chrom = if (length(args) >= 2L) args[2] else NULL)

suppressPackageStartupMessages(library(data.table))
source(file.path(params$code_root, "parsing_functions", "stats_testing_functions.r"), local = TRUE)
data.table::setDTthreads(1L)

## Load inputs
chrom_scope <- params$chrom
cores <- as.integer(params$cores)
input_paths <- c(
  counts = file.path(params$root, "coaccess", "LCL_pair_sample_counts.tsv.gz"),
  pooled = file.path(params$root, "coaccess", "LCL_coaccess_pairs.tsv.gz"),
  sex = file.path(params$root, "coaccess", "noXCI", "LCL_sample_sex_noXCI.tsv"))
input_md5 <- tools::md5sum(input_paths)
sex <- fread(input_paths[["sex"]])
pooled <- fread(input_paths[["pooled"]])
annotation_columns <- c("CRE_pair", "CRE1", "CRE2", "chrom", "cre1_start", "cre1_end",
  "cre2_start", "cre2_end", "dist", "CRE_pair_labels", "gene_ids", "gene_names", "n_genes")

counts <- fread(input_paths[["counts"]])

counts <- counts[CRE_pair %in% pooled$CRE_pair]
ann <- pooled[, .(CRE_pair, OR_pooled = OR, fdr_pooled = fdr, n_shared_reads_pooled = n_shared_reads,
                   co_access_pooled = co_access, expected_co_access_pooled = expected_co_access)]
annotations <- merge(pooled[, ..annotation_columns], ann, by = "CRE_pair", sort = FALSE)
rm(pooled)
if (!is.null(chrom_scope)) {
  counts <- counts[chrom == chrom_scope]
  annotations <- annotations[chrom == chrom_scope]
}
setorder(annotations, chrom, CRE_pair)
fdr_scope <- if (is.null(chrom_scope)) "genome-wide" else paste0("chromosome:", chrom_scope)
output_base <- file.path(params$root, "coaccess", "donor_adjusted")
output_dir <- if (is.null(chrom_scope)) output_base else
  file.path(output_base, "smoke_tests", chrom_scope)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat(nrow(annotations), "pairs; FDR scope:", fdr_scope, "\n")

## Aggregate by donor
# donor_all uses every sample on every chr; donor_male_x keeps only chrX rows from male samples
donor_all <- aggregate_donors(counts)
male_samples <- sex[sex == "male (XY)", sample_name]
stopifnot(length(male_samples) > 0L)
donor_male_x <- aggregate_donors(counts[chrom == "chrX" & sample_name %in% male_samples])
rm(counts)

print(head(donor_all))
print(head(donor_male_x))

## Fit the regression model
# turns any R object into a unique code that changes whenever the object’s contents change
# check if saved results still match the current settings
hash_object <- function(x) {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(x, tmp, version = 3)
  unname(tools::md5sum(tmp))
}

##########################################
# fits the donor-adjusted logistic regression model on every CRE pair
##########################################

# d = donor aggregated table for this chr
# male only = exclude female chr x

# Fits two logistic regressions on those rows:
# Null: CRE2 ~ factor(donor), where CRE2 depends only on donor.
# Full: CRE2 ~ factor(donor) + x, where x is CRE1 status (1 = open, 0 = closed)

fit_chromosome <- function(d, chromosome, male_only = FALSE) {
  # remove duplicate CRE ID
  expected_ids <- sort(unique(d$CRE_pair), method = "radix")
  # build the sample path
  fit_summary_file <- file.path(fit_summary_dir, paste0(chromosome,
    if (male_only) "_male_only" else "_all_donors", ".rds"))
  if (file.exists(fit_summary_file)) {
    saved <- readRDS(fit_summary_file)
    stopifnot(identical(saved$signature, summary_signature),
              !anyDuplicated(saved$fits$CRE_pair),
              identical(sort(saved$fits$CRE_pair, method = "radix"), expected_ids))
    message("Reusing summary: ", basename(fit_summary_file))
    return(saved$fits)
  }
  # Sorts the table by pair, then donor --> all rows of one pair are together
  setorder(d, CRE_pair, donor)
  groups <- d[, .(first = .I[1], last = .I[.N]), by = CRE_pair]
  message("Fitting ", nrow(groups), " pairs on ", chromosome,
          if (male_only) " (male only)" else " (all donors)", " using ", cores, " cores")

  # fit logistic model
  ####################
  # CRE2 accessible is the outcome, CRE1 accessible as predictor, with separate baseline for each door

  # loop over # of pairs on each chr, and takes one row per donor with four counts
  # (co_closed, CRE1_access, CRE2_access, co_access)
  fits <- parallel::mclapply(seq_len(nrow(groups)), function(i) {
    fit_pair(d[seq.int(groups$first[i], groups$last[i])], params$min_reads_per_donor)
  }, mc.cores = cores, mc.preschedule = TRUE, mc.set.seed = FALSE)
  if (any(vapply(fits, inherits, logical(1), what = "try-error")) ||
      any(vapply(fits, is.null, logical(1))))
    stop("A fitting worker failed on ", chromosome, "; no pairs were discarded")
  fits <- rbindlist(fits)
  stopifnot(nrow(fits) == length(expected_ids), !anyDuplicated(fits$CRE_pair),
            setequal(fits$CRE_pair, expected_ids))
  saveRDS(list(signature = summary_signature, fits = fits), fit_summary_file)
  fits
}

##########################################
# run logistic model for all chr or male only
##########################################
run_analysis <- function(analysis_name) {
  chromosomes <- unique(annotations$chrom)
  # stact fitted result into one tab;e
  fits <- rbindlist(lapply(chromosomes, function(chr) {
    male_only <- analysis_name == "noXCI" && chr == "chrX"
    d <- if (male_only) copy(donor_male_x) else copy(donor_all[chrom == chr])
    fit_chromosome(d, chr, male_only)
  }))
  # add annotation
  out <- merge(annotations, fits, by = "CRE_pair", all.x = TRUE, sort = FALSE)
  stopifnot(nrow(out) == nrow(annotations), !anyNA(out$fit_status),
            !anyDuplicated(out$CRE_pair), setequal(out$CRE_pair, annotations$CRE_pair))
  out[, chrX_rule := fifelse(analysis_name == "noXCI" & chrom == "chrX",
                             "male samples only", "all samples")]
  if (analysis_name == "noXCI")
    out[chrom == "chrX" & fit_status == "no_coverage" &
          !is.na(n_shared_reads_pooled) & n_shared_reads_pooled > 0,
        fit_status := "no_coverage_noXCI"]
  # remove pairs without reads
  excluded_zero <- sum(out$fit_status %in% c("no_coverage", "no_coverage_noXCI"))
  message(analysis_name, ": excluding ", excluded_zero, " zero-coverage pairs after cohort filtering")
  out <- out[!fit_status %in% c("no_coverage", "no_coverage_noXCI")]
  # count the number of pairs
  family_size <- sum(out$n_donors_informative >= 1L)
  valid <- which(is.finite(out$p_adj_one_sided))
  stopifnot(all(out$p_adj_one_sided[valid] >= 0 & out$p_adj_one_sided[valid] <= 1),
            all(out$n_donors_informative[valid] >= 1L))
  out[, fdr_adj := NA_real_]
  if (length(valid)) out[valid, fdr_adj := adjust_bh(p_adj_one_sided, n = family_size)]
  # add metadata + full and null model formulas for each pair
  out[, `:=`(analysis = analysis_name, min_reads_per_donor = as.integer(params$min_reads_per_donor),
    model = fifelse(n_donors_informative == 1L, "CRE2 ~ 1 + CRE1", "CRE2 ~ factor(donor) + CRE1"),
    null_model = fifelse(n_donors_informative == 1L, "CRE2 ~ 1", "CRE2 ~ factor(donor)"),
    fdr_scope = fdr_scope, bh_family_size = family_size, run_signature = summary_signature,
    ci_level = 0.95, glm_maxit = 100L,
    read_rule = "contain", fire_overlap_fraction = 0.5,
    input_coverage_mode = "aligned_blocks_v1")]
  # sort pairs by FDR ascending
  ranks <- out[is.finite(p_adj_one_sided) & !is.na(OR_adj), .(CRE_pair, fdr_adj, OR_adj)]
  setorderv(ranks, c("fdr_adj", "OR_adj", "CRE_pair"), c(1L, -1L, 1L))
  out[, adjusted_rank := match(CRE_pair, ranks$CRE_pair)]
  setorder(out, chrom, CRE_pair)
  out
}

# save the result
summary_signature <- hash_object(list(inputs = input_md5,
  min_reads_per_donor = params$min_reads_per_donor,
  statistics_md5 = tools::md5sum(file.path(params$code_root, "parsing_functions", "stats_testing_functions.r")),
  functions = lapply(list(aggregate_donors, capture_fit, fit_pair, fit_chromosome, run_analysis), deparse),
  R = R.version.string, logistf = as.character(packageVersion("logistf"))))
fit_summary_dir <- file.path(output_base, "summary", summary_signature)
dir.create(fit_summary_dir, recursive = TRUE, showWarnings = FALSE)

## Run A: all chromosomes and all donors
all_chr_results <- run_analysis("all_chr")

## Run B: exclude female chrX samples
noxci_results <- run_analysis("noXCI")

## Outputs
stopifnot(identical(unname(tools::md5sum(input_paths)), unname(input_md5)))
result_sets <- list(all_chr = all_chr_results, noXCI = noxci_results)
for (a in names(result_sets)) {
  path <- file.path(output_dir, a, paste0("LCL_coaccess_pairs_donor_adjusted_", a, ".tsv.gz"))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(result_sets[[a]], path, sep = "\t", na = "NA")
  cat(a, ":", nrow(result_sets[[a]]), "rows ->", path, "\n")
}
