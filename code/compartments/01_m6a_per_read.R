source("/project/spott/cshan/fiber-seq/code/compartments/compartments_functions.r")
bins <- make_bins(); chunks <- make_chunks(bins); samples <- read_samples()
args <- commandArgs(TRUE)
if (length(args) && args[1] == "--list") { print(chunks); cat("Tasks:", nrow(chunks), "\n"); quit() }
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", if (length(args)) args[1] else "1"))
stopifnot(!is.na(task), task >= 1, task <= nrow(chunks))
c <- chunks[task]
windows <- seq(c$start, c$end - 1L, by = 100000L)
rows <- lapply(windows, function(s) {
  message(c$chr, ":", s)
  b <- bins[chr == c$chr & start >= s & end <= min(s + 100000L, c$end)]
  window_per_read(samples, region_row(c$chr, s, max(b$end)), b)
})
write_table(rbindlist(rows), file.path(OUT, "per_read", paste0(c$chunk, ".tsv.gz")))
