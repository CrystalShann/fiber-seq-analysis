# Coordinates in output: BED [start,end); observed call positions: 1-based.
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
PROJECT <- Sys.getenv("COMPARTMENTS_PROJECT", "/project/spott/cshan/fiber-seq")
source(file.path(PROJECT, "code/parsing_functions/parsing_footprints_functions.r"))
source(file.path(PROJECT, "code/parsing_functions/plotting_functions.r"))
setDTthreads(1L)
CHROMS <- strsplit(Sys.getenv("COMPARTMENTS_CHROMS", paste0("chr", 1:22, collapse = ",")), ",")[[1]]
stopifnot(length(CHROMS) > 0, !anyDuplicated(CHROMS), all(CHROMS %in% paste0("chr", 1:22)))
OUT <- Sys.getenv("COMPARTMENTS_OUT", file.path(PROJECT, "LCL_project/compartments",
  if (!identical(CHROMS, paste0("chr", 1:22))) paste(CHROMS, collapse = "_" ) else ""))
SAMPLES_FILE <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv"
RAO_DIR <- "/project/spott/cshan/annotations/LCL_public_data/rao2014_GM12878/hg38"
read_samples <- function() {
  x <- fread(SAMPLES_FILE)
  stopifnot(nrow(x) == 31L, !anyDuplicated(x$sample_name))
  x
}
write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp", if (grepl("gz$", path)) ".gz" else "")
  fwrite(x, tmp, sep = "\t", na = "NA")
  if (!file.rename(tmp, path)) stop("Cannot rename output: ", path)
}
make_bins <- function(sizes = fread("/project/spott/cshan/annotations/hg38.chrom.sizes", header = FALSE),
                      chroms = CHROMS, width = 20000L) {
  rbindlist(lapply(chroms, function(ch) {
    len <- sizes[[2]][match(ch, sizes[[1]])]
    stopifnot(length(len) == 1L, !is.na(len), len >= width)
    s <- seq.int(0L, (len %/% width - 1L) * width, by = width)
    data.table(chr = ch, start = s, end = s + width,
      bin = paste0(ch, ":", s, "-", s + width), midpoint = s + width / 2)
  }))
}
make_chunks <- function(bins, width = 5000000L) {
  bins[, .(start = min(start), end = max(end)), by = .(chr, chunk_index = start %/% width)][,
    chunk := paste0(chr, "_", start, "_", end)]
}
region_row <- function(chr, start, end) data.frame(chr = chr, start = start, end = end,
  analysis_start = start + 1L, analysis_end = end, region_id = paste(chr, start, end, sep = "_"))
assemble_safe <- function(sample, region, positions = NULL) {
  tryCatch(assemble_region_m6a(sample, region, "m6a", full_span = FALSE, positions = positions),
    error = function(e) { if (startsWith(conditionMessage(e), "Insufficient ")) NULL else stop(e) })
}
empty_reads <- function() data.table(read = character(), sample = character(), bin = character(),
  n_calls = numeric(), n_positions = integer(), overlap_bp = integer())
window_per_read <- function(samples, window, bins, assemble = assemble_safe) {
  # Retain only calls and spans after each sample, avoiding 31 dense matrices in memory.
  parts <- lapply(seq_len(nrow(samples)), function(j) {
    a <- assemble(samples[j, ], window)
    if (is.null(a)) return(NULL)
    trip <- summary(a$met_mat)
    trip <- trip[!is.na(trip$x) & trip$x > 0, ]
    list(meta = a$rids_df, calls = data.table(i = trip$i,
      pos = as.integer(colnames(a$met_mat))[trip$j]))
  })
  sites <- sort(unique(as.integer(unlist(lapply(parts, function(a) a$calls$pos)))))
  rows <- list()
  for (a in parts) {
    if (is.null(a)) next
    m <- a$meta
    for (k in seq_len(nrow(bins))) {
      b <- bins[k]
      lo <- pmax(m$start, b$start + 1L); hi <- pmin(m$end, b$end)
      keep <- which(lo <= hi)
      if (!length(keep)) next
      p <- sites[sites > b$start & sites <= b$end]
      counts <- tabulate(a$calls$i[a$calls$pos > b$start & a$calls$pos <= b$end], nrow(m))
      denom <- findInterval(hi[keep], p) - findInterval(lo[keep] - 1L, p)
      rows[[length(rows) + 1L]] <- data.table(read = m$RID[keep], sample = m$sample_name[keep],
        bin = b$bin, n_calls = counts[keep], n_positions = denom, overlap_bp = hi[keep] - lo[keep] + 1L)
    }
  }
  x <- if (length(rows)) rbindlist(rows) else empty_reads()
  stopifnot(all(x$n_calls <= x$n_positions), all(x$n_calls >= 0), !anyDuplicated(x[, .(read, sample, bin)]))
  x
}
bin_rates <- function(reads, pooled = FALSE) {
  groups <- if (pooled) "bin" else c("bin", "sample")
  reads[, .(n_calls = sum(n_calls), n_positions = sum(n_positions), n_reads = .N), by = groups][,
    rate := fifelse(n_positions > 0, n_calls / n_positions, NA_real_)]
}
read_rate_distributions <- function(reads, pooled = FALSE) {
  groups <- if (pooled) "bin" else c("bin", "sample")
  reads[overlap_bp >= 5000 & n_positions > 0, {
    r <- n_calls / n_positions
    q <- quantile(r, c(.05, .25, .5, .75, .95), names = FALSE)
    .(n_reads = .N, mean = mean(r), sd = sd(r), q05 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q95 = q[5])
  }, by = groups]
}
load_bigwig <- function() {
  curl <- "/project/spott/cshan/envs/dimelo/lib/libcurl.so.4"
  if (file.exists(curl)) dyn.load(curl, local = FALSE)
  loadNamespace("rtracklayer")
}
coverage_sums <- function(path, bins) {
  if (!file.exists(path)) stop("Missing bigWig: ", path)
  gr <- GenomicRanges::GRanges(bins$chr, IRanges::IRanges(bins$start + 1L, bins$end))
  chunk <- GenomicRanges::GRanges(bins$chr[1], IRanges::IRanges(min(bins$start) + 1L, max(bins$end)))
  v <- rtracklayer::import(path, which = chunk)
  h <- GenomicRanges::findOverlaps(gr, v)
  i <- S4Vectors::queryHits(h); j <- S4Vectors::subjectHits(h)
  sums <- numeric(nrow(bins))
  if (length(i)) {
    w <- pmin(bins$end[i], BiocGenerics::end(v)[j]) - pmax(bins$start[i] + 1L, BiocGenerics::start(v)[j]) + 1L
    z <- data.table(i = i, value = w * v$score[j])[, .(value = sum(value)), by = i]
    sums[z$i] <- z$value
  }
  sums
}
nucleosome_fraction <- function(sample, bins) {
  sums <- lapply(c("nucleosome", "linker", "fire"), function(track)
    coverage_sums(file.path(sample$fire_dir, "trackHub-v0.1/bw", paste0("all.", track, ".coverage.bw")), bins))
  total <- Reduce(`+`, sums)
  data.table(bin = bins$bin, sample = sample$sample_name, nuc_sum = sums[[1]], total_sum = total,
    mean_total_depth = total / (bins$end - bins$start), fraction = fifelse(total > 0, sums[[1]] / total, NA_real_))
}
normalize_bins <- function(x, bins, samples, signal) {
  value <- if (signal == "m6a") "rate" else "fraction"
  qc <- if (signal == "m6a") "n_reads" else "mean_total_depth"
  m <- matrix(NA_real_, nrow(bins), length(samples), dimnames = list(bins$bin, samples))
  stopifnot(!anyDuplicated(x[, .(bin, sample)]))
  ij <- cbind(match(x$bin, bins$bin), match(x$sample, samples))
  stopifnot(!anyNA(ij))
  v <- x[[value]]; v[!is.finite(v) | v <= 0 | is.na(x[[qc]]) | x[[qc]] < if (signal == "m6a") 10 else 5] <- NA_real_
  m[ij] <- log(v)
  if (any(colSums(is.finite(m)) < 2L)) stop("Too few usable bins in a sample for normalization")
  limma::normalizeQuantiles(m)
}
pair_domain_class <- function(left, right, domains, edges) {
  # Edges strictly between midpoints; midpoint membership is BED half-open.
  e <- sort(unique(edges))
  across <- findInterval(right - .5, e) > findInterval(left, e)
  common <- rep(FALSE, length(left))
  # Each left midpoint needs only the greatest end of domains containing it.
  max_end <- rep(-Inf, length(left))
  for (k in seq_len(nrow(domains))) {
    inside <- left >= domains$start[k] & left < domains$end[k]
    max_end[inside] <- pmax(max_end[inside], domains$end[k])
  }
  common <- right < max_end
  fifelse(across, "across", fifelse(common, "within", NA_character_))
}
lag_correlations <- function(m, bins, domains, edges, max_lag = 100L, min_samples = 20L) {
  if (nrow(m) < 2L) return(data.table())
  # Cache domain membership once, then classify each shifted pair in O(n).
  max_end <- rep(-Inf, nrow(bins))
  for (k in seq_len(nrow(domains))) {
    i <- which(bins$midpoint >= domains$start[k] & bins$midpoint < domains$end[k])
    max_end[i] <- pmax(max_end[i], domains$end[k])
  }
  e <- sort(unique(edges))
  rbindlist(lapply(seq_len(min(max_lag, nrow(m) - 1L)), function(lag) {
    i <- seq_len(nrow(m) - lag); j <- i + lag
    a <- m[i, , drop = FALSE]; b <- m[j, , drop = FALSE]
    ok <- is.finite(a) & is.finite(b); n <- rowSums(ok)
    a[!ok] <- 0; b[!ok] <- 0
    a <- a - rowSums(a) / pmax(n, 1); b <- b - rowSums(b) / pmax(n, 1)
    a[!ok] <- 0; b[!ok] <- 0
    r <- rowSums(a * b) / sqrt(rowSums(a * a) * rowSums(b * b))
    r[n < min_samples | !is.finite(r)] <- NA_real_
    r <- pmax(-1, pmin(1, r))
    across <- findInterval(bins$midpoint[j] - .5, e) > findInterval(bins$midpoint[i], e)
    cls <- fifelse(across, "across", fifelse(bins$midpoint[j] < max_end[i], "within", NA_character_))
    data.table(bin = bins$bin[i], lag = lag, r = r, n_samples = n, class = cls)
  }))
}
boundary_summary <- function(pairs) {
  pairs[!is.na(class) & is.finite(r), .(n_pairs = .N, mean = mean(r), median = median(r),
    q25 = quantile(r, .25), q75 = quantile(r, .75), iqr = IQR(r)), by = .(lag, class)]
}
