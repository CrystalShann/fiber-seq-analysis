source("/project/spott/cshan/fiber-seq/code/compartments/compartments_functions.r")
load_bigwig()
samples <- read_samples(); bins <- make_bins(); chunks <- make_chunks(bins)
args <- commandArgs(TRUE)
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", if (length(args)) args[1] else "1"))
stopifnot(!is.na(task), task >= 1, task <= nrow(samples))
sample <- samples[task]
rows <- lapply(seq_len(nrow(chunks)), function(i) {
  c <- chunks[i]; message(sample$sample_name, " ", c$chunk)
  nucleosome_fraction(sample, bins[chr == c$chr & start >= c$start & end <= c$end])
})
write_table(rbindlist(rows), file.path(OUT, "tables/nucleosome_fraction", paste0(sample$sample_name, ".tsv.gz")))
