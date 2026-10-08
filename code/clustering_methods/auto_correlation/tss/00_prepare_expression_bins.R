#!/usr/bin/env Rscript

###############################################
# build expression bins 
###############################################

# assigns every protein coding gene an expression bin from RNA-seq, using the mean TPM over the four LPS timepoints

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[1] != "--out-dir")
  stop("Usage: 00_prepare_expression_bins.R --out-dir OUTPUT/inputs")
out <- normalizePath(args[2], mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", script_arg)))
project <- normalizePath(file.path(here, "../../../.."))
source_file <- file.path(project, "code/accessibility/expr_access/01_expression_bins.R")
aligned <- file.path(project, "code/RNA/aligned")
inputs <- c(source_file, file.path(here, "00_prepare_expression_bins.R"),
  "/project/spott/cshan/annotations/gencode.v46.annotation.gtf.gz",
  "/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed",
  list.files(aligned, "ReadsPerGene.out.tab$", recursive = TRUE, full.names = TRUE))
if (!length(inputs) || any(!file.exists(inputs))) stop("Missing expression prerequisites")
signature <- list(md5 = tools::md5sum(inputs), R = R.version.string,
                  data_table = as.character(packageVersion("data.table")))
target <- file.path(out, "tss_expression_bins.tsv")
manifest <- file.path(out, "expression_bins.manifest.rds")
if (file.exists(manifest) && file.exists(target)) {
  previous <- tryCatch(readRDS(manifest), error = function(e) NULL)
  if (!is.null(previous) && identical(previous$signature, signature) &&
      identical(previous$output_md5, tools::md5sum(target))) {
    message("Expression preparation: reuse validated original-code output ", target)
    quit(status = 0L)
  }
}
expressions <- parse(source_file)
env <- new.env(parent = globalenv())
n_redirected <- 0L
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("OUT_DIR"))) {
    env$OUT_DIR <- out
    n_redirected <- n_redirected + 1L
  } else eval(expr, envir = env)
}
stopifnot(n_redirected == 1L, file.exists(target))
temporary <- paste0(manifest, ".tmp")
saveRDS(list(signature = signature, output_md5 = tools::md5sum(target)), temporary)
stopifnot(file.rename(temporary, manifest))
message("Prepared expression table using unchanged source: ", target)
