#!/usr/bin/env Rscript
# Run from the repository root; fixtures and plot data remain in memory.
suppressPackageStartupMessages(library(data.table))
lines <- readLines('code/enhancer/enhancer_co-occ.Rmd')
for (label in c('function-load-enhancer-example', 'function-plot-enhancer-example')) {
  start <- grep(paste0('^```\\{r ', label, '[,}]'), lines)
  stopifnot(length(start) == 1L)
  end <- start + which(lines[(start + 1L):length(lines)] == '```')[1L]
  eval(parse(text = lines[seq.int(start + 1L, end - 1L)]))
}
options(warn = 2)
# Same read, two alignments with different block counts: take only the selected
# alignment's real blocks, excluding the first/last extraction sentinels.
raw <- data.table(V1 = 'chr1', V2 = c(0L, 100L), V3 = c(100L, 200L),
  V4 = 'read_a', V5 = 0, V6 = '+', V7 = c(0L, 100L), V8 = c(100L, 200L),
  V9 = 0, V10 = c(3L, 4L), V11 = c('1,1,1', '1,1,1,1'),
  V12 = c('0,10,99', '0,20,40,99'))
reads <- data.table(read_id = 'read_a', start0 = 100, end0 = 200)
blocks <- enhancer_example_m6a_blocks(raw, reads)
stopifnot(identical(blocks$start0, c(120, 140)),
          identical(blocks$end0, c(121, 141)))
# Duplicate records do not double-count calls; sentinel-only reads have no marks.
stopifnot(identical(enhancer_example_m6a_blocks(rbind(raw, raw), reads), blocks))
raw[2L, `:=`(V10 = 2L, V11 = '1,1', V12 = '0,99')]
stopifnot(nrow(enhancer_example_m6a_blocks(raw, reads)) == 0L)

# Fractions use covering displayed fibers, deduplicate calls, and leave
# uncovered positions/timepoints missing instead of counting them as zero.
example <- list(window = c(0L, 6L), samples = data.table(sample = c('LPS_0', 'LPS_5')),
  spans = data.table(sample = 'LPS_0', read_id = c('a', 'b'),
    start0 = c(0L, 2L), end0 = c(4L, 5L)),
  m6a = data.table(sample = 'LPS_0', read_id = c('a', 'a', 'b'),
    start0 = c(2L, 2L, 4L), end0 = c(3L, 3L, 5L)))
pileup <- enhancer_example_pileup(example)
stopifnot(pileup[sample == 'LPS_0' & position == 2L, fraction] == 0.5,
          pileup[sample == 'LPS_0' & position == 4L, fraction] == 1,
          is.na(pileup[sample == 'LPS_0' & position == 5L, fraction]),
          all(is.na(pileup[sample == 'LPS_5', fraction])),
          all(pileup$modified <= pileup$coverage))
cat('PASS: selected-alignment m6A blocks, sentinel removal, duplicate calls and coverage-aware pileup.\n')
