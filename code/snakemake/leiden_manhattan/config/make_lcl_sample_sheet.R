#!/usr/bin/env Rscript
# Build config/samples_lcl.tsv from the LCL sample metatable used by leiden_LCL.Rmd.
#
#   Rscript make_lcl_sample_sheet.R [metatable.csv] [phasing_root] [output.tsv]
#
# Columns written (same layout as samples_macrophage.tsv where they overlap):
#   sample_name   e.g. AL10_bc2178_19130 (unchanged; used in read IDs "<sample_name>::<read>")
#   sample_label  prefix before the first "_" (e.g. AL10), as used by sample_palette()
#   cell_line     VCF sample column used by LCL_phasing.r
#   fire_dir      root holding extracted_results/{m6a,nuc}_by_chr/
#   phasing_dir   <phasing_root>/<sample_name>, holding blocks.tsv,
#                 read-level-phasing.tsv and <sample_name>.5mC.6mA.aligned.phased.vcf.gz
# Any other metatable columns are kept after these.

args <- commandArgs(trailingOnly = TRUE)
metatable <- if (length(args) >= 1) args[[1]] else
  "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv"
phasing_root <- if (length(args) >= 2) args[[2]] else
  "/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples"
output <- if (length(args) >= 3) args[[3]] else "samples_lcl.tsv"

x <- read.csv(metatable, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("sample_name", "fire_dir", "cell_line")
missing <- setdiff(required, names(x))
if (length(missing)) stop("Metatable lacks columns: ", paste(missing, collapse = ", "))
if (anyDuplicated(x$sample_name)) stop("Duplicate sample_name values in ", metatable)

x$sample_label <- sub("_.*$", "", x$sample_name)
if (anyDuplicated(x$sample_label)) stop("sample_label (prefix before '_') is not unique")
x$phasing_dir <- file.path(phasing_root, x$sample_name)

front <- c("sample_name", "sample_label", "cell_line", "fire_dir", "phasing_dir")
x <- x[, c(front, setdiff(names(x), front)), drop = FALSE]

# Report, but do not fail on, inputs that are missing on this machine.
expected <- c(
  file.path(x$fire_dir, "extracted_results", "m6a_by_chr"),
  file.path(x$fire_dir, "extracted_results", "nuc_by_chr"),
  file.path(x$phasing_dir, "blocks.tsv"),
  file.path(x$phasing_dir, "read-level-phasing.tsv"),
  file.path(x$phasing_dir, paste0(x$sample_name, ".5mC.6mA.aligned.phased.vcf.gz")))
absent <- expected[!file.exists(expected)]
if (length(absent))
  warning(length(absent), " expected input path(s) not found, e.g. ", absent[[1]], call. = FALSE)

write.table(x, output, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", nrow(x), " samples to ", normalizePath(output))
