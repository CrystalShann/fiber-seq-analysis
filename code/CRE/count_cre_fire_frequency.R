#!/usr/bin/env Rscript
# Explicit counting stage; never run by the summary notebook.
# Rscript count_cre_fire_frequency.R [--test]
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% "--test")) stop("Usage: count_cre_fire_frequency.R [--test]")
params <- list(test = "--test" %in% args)

suppressPackageStartupMessages({library(data.table)})

out_root <- Sys.getenv("CRE_ACCESS_ROOT", "/project/spott/cshan/fiber-seq/LCL_project/CRE_access")
sample_sheet <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv"
bedtools <- "/project/spott/cshan/envs/bedtools/bin/bedtools"
region_dir <- file.path(out_root, "regions")
per_read <- file.path(out_root, "per_read")
run_dir <- out_root
sets <- c("active_cre", "low_dnase_cre")
chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
samples <- fread(sample_sheet)$sample_name
if (params$test) {
  run_dir <- file.path(out_root, "test_chr21"); chroms <- "chr21"; samples <- samples[1]
}

inputs <- CJ(s = samples, chrom = chroms, kind = c("reads", "fire"), sorted = FALSE)[
  , file.path(per_read, s, sprintf("%s.%s.%s.bed.gz", s, chrom, kind))]
if (!all(file.exists(inputs)))
  stop("Missing per-read BEDs (run cre_fire_freq_dnase.sh chrom first), e.g. ", inputs[!file.exists(inputs)][1])

# unique (region id, read name) pairs; -a is BED4, so the read name is column 8
pairs <- function(regions, b, opt = "") {
  x <- fread(cmd = sprintf("%s intersect -sorted %s -wa -wb -a %s -b %s | cut -f4,8", bedtools, opt, regions, b),
             header = FALSE, sep = "\t", colClasses = "character")
  if (!nrow(x)) return(data.table(id = character(), read = character()))
  unique(setnames(x, c("id", "read")))
}

by_chrom_dir <- file.path(run_dir, "by_chrom")
dir.create(by_chrom_dir, recursive = TRUE, showWarnings = FALSE)
# Reset region settings when this chunk is rerun in an interactive R session.
sets <- c("active_cre", "low_dnase_cre")
region_dir <- file.path(Sys.getenv("CRE_ACCESS_ROOT",
  "/project/spott/cshan/fiber-seq/LCL_project/CRE_access"), "regions")
region_beds <- file.path(region_dir, paste0(sets, ".bed"))
missing_regions <- region_beds[file.access(region_beds, 4L) != 0L]
if (length(missing_regions))
  stop("Missing or unreadable region BEDs: ", paste(missing_regions, collapse = ", "),
       ". Generate them with: bash code/CRE/cre_fire_freq_dnase.sh regions")
regions <- NULL
tmp <- tempfile("cre-access."); dir.create(tmp)

for (chrom in chroms) {
  tables <- file.path(by_chrom_dir, sprintf("%s.%s.tsv", sets, chrom))
  chrom_inputs <- c(region_beds, grep(sprintf("\\.%s\\.(reads|fire)\\.bed\\.gz$", chrom), inputs, value = TRUE))
  if (all(file.exists(tables)) && min(file.mtime(tables)) > max(file.mtime(chrom_inputs))) next
  if (is.null(regions)) regions <- setNames(lapply(region_beds, fread, header = FALSE, sep = "\t"), sets)

  for (set in sets) {
    bed4 <- file.path(tmp, paste0(set, ".bed4"))
    reg <- regions[[set]][V1 == chrom]
    fwrite(reg[, 1:4], bed4, sep = "\t", col.names = FALSE)
    out <- data.table(chrom = reg$V1, start = reg$V2, end = reg$V3, id = reg$V4)
    if (ncol(reg) >= 10) out[, class := reg$V10]
    for (s in samples) {
      f <- function(kind) file.path(per_read, s, sprintf("%s.%s.%s.bed.gz", s, chrom, kind))
      cover <- pairs(bed4, f("reads"), "-f 1.0")
      fire <- cover[pairs(bed4, f("fire")), on = .(id, read), nomatch = NULL]
      nr <- paste0("n_reads_", s); nf <- paste0("n_fire_", s)
      out[, (nr) := 0L][, (nf) := 0L]
      out[cover[, .N, by = id], on = "id", (nr) := i.N]
      out[fire[, .N, by = id], on = "id", (nf) := i.N]
    }
    nr <- paste0("n_reads_", samples); nf <- paste0("n_fire_", samples)
    out[, n_reads := as.integer(rowSums(.SD)), .SDcols = nr]
    out[, n_fire := as.integer(rowSums(.SD)), .SDcols = nf]
    out[, fire_freq := fifelse(n_reads > 0, n_fire / n_reads, NA_real_)]
    setcolorder(out, c(intersect(c("chrom", "start", "end", "id", "class"), names(out)),
                       "n_reads", "n_fire", "fire_freq", as.vector(rbind(nr, nf))))
    bad <- sum(out$n_fire > out$n_reads) + sum(mapply(function(r, f) sum(out[[f]] > out[[r]]), nr, nf))
    cat(sprintf("[%s] %-20s regions %7d | n_reads > 0: %7d | median n_reads %5.0f | n_fire > n_reads: %d\n",
                chrom, set, nrow(out), sum(out$n_reads > 0), median(out$n_reads), bad))
    fwrite(out, file.path(by_chrom_dir, sprintf("%s.%s.tsv", set, chrom)), sep = "\t", na = "NA", quote = FALSE)
  }
}
unlink(tmp, recursive = TRUE)
rm(regions)

for (set in sets) {
  tabs <- lapply(file.path(by_chrom_dir, sprintf("%s.%s.tsv", set, chroms)), fread)
  if (length(unique(lapply(tabs, names))) != 1) stop(set, ": chromosome tables have different sample columns")
  fwrite(rbindlist(tabs), file.path(run_dir, sprintf("fire_freq.%s.tsv", set)), sep = "\t", na = "NA", quote = FALSE)
}
rm(tabs)
cat(sprintf("Combined %d chromosome(s) into %s/fire_freq.<set>.tsv\n", length(chroms), run_dir))
