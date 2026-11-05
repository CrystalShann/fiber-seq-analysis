source("/project/spott/cshan/fiber-seq/code/compartments/compartments_functions.r")
bins <- make_bins(); samples <- read_samples(); chunks <- make_chunks(bins)
files <- file.path(OUT, "per_read", paste0(chunks$chunk, ".tsv.gz"))
nuc_files <- file.path(OUT, "tables/nucleosome_fraction", paste0(samples$sample_name, ".tsv.gz"))
if (!all(file.exists(c(files, nuc_files)))) stop("Incomplete arrays: missing chunk or sample outputs")
# Aggregate one chunk at a time; never load genome-wide per-read records together.
rates <- distributions <- vector("list", length(files))
pooled_rates <- pooled_distributions <- rates
for (i in seq_along(files)) {
  message(files[i]); x <- fread(files[i])
  stopifnot(!anyDuplicated(x[, .(read, sample, bin)]), all(x$n_calls <= x$n_positions))
  rates[[i]] <- bin_rates(x); pooled_rates[[i]] <- bin_rates(x, TRUE)
  distributions[[i]] <- read_rate_distributions(x); pooled_distributions[[i]] <- read_rate_distributions(x, TRUE)
}
complete <- function(x, pooled, counts) {
  grid <- if (pooled) bins[, .(bin)] else CJ(bin = bins$bin, sample = samples$sample_name)
  x <- merge(grid, x, by = if (pooled) "bin" else c("bin", "sample"), all.x = TRUE, sort = FALSE)
  for (col in counts) set(x, which(is.na(x[[col]])), col, 0)
  merge(bins, x, by = "bin", sort = FALSE)
}
for (pooled in c(FALSE, TRUE)) {
  suffix <- if (pooled) "pooled" else "per_sample"
  x <- complete(rbindlist(if (pooled) pooled_rates else rates), pooled, c("n_reads", "n_calls", "n_positions"))
  write_table(x, file.path(OUT, "tables", paste0("bin_m6a_rate_", suffix, ".tsv.gz")))
  x <- complete(rbindlist(if (pooled) pooled_distributions else distributions), pooled, "n_reads")
  write_table(x, file.path(OUT, "tables", paste0("read_m6a_rate_distributions_", suffix, ".tsv.gz")))
}
x <- rbindlist(lapply(nuc_files, fread))
stopifnot(nrow(x) == nrow(bins) * nrow(samples), !anyDuplicated(x[, .(bin, sample)]),
  setequal(x$bin, bins$bin), setequal(x$sample, samples$sample_name))
write_table(merge(bins, x, by = "bin", sort = FALSE), file.path(OUT, "tables/bin_nucleosome_fraction_per_sample.tsv.gz"))
p <- x[, .(nuc_sum = sum(nuc_sum), total_sum = sum(total_sum)), by = bin]
p <- merge(bins, p, by = "bin", sort = FALSE)
p[, `:=`(fraction = fifelse(total_sum > 0, nuc_sum / total_sum, NA_real_),
  mean_total_depth = total_sum / ((end - start) * nrow(samples)))]
write_table(p, file.path(OUT, "tables/bin_nucleosome_fraction_pooled.tsv.gz"))
