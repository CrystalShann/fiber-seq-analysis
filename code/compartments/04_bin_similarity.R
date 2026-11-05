source("/project/spott/cshan/fiber-seq/code/compartments/compartments_functions.r")
bins <- make_bins(); samples <- read_samples(); rao <- load_rao_hg38(RAO_DIR)
write_table(rao$domains, file.path(OUT, "tables/rao_domains.tsv"))
write_table(rao$edges, file.path(OUT, "tables/rao_edges.tsv"))
all_summary <- list()
for (signal in c("m6a", "nucleosome")) {
  input <- if (signal == "m6a") "bin_m6a_rate_per_sample.tsv.gz" else "bin_nucleosome_fraction_per_sample.tsv.gz"
  m <- normalize_bins(fread(file.path(OUT, "tables", input)), bins, samples$sample_name, signal)
  write_table(cbind(bins, as.data.table(m)), file.path(OUT, "tables", paste0("normalized_", signal, ".tsv.gz")))
  q <- rbindlist(lapply(seq_len(ncol(m)), function(j) data.table(sample = colnames(m)[j],
    probability = seq(0, 1, .01), value = quantile(m[, j], seq(0, 1, .01), na.rm = TRUE))))
  write_table(q, file.path(OUT, "tables", paste0("normalization_quantiles_", signal, ".tsv")))
  # Retain only r/lag/class for exact genome-wide medians and IQRs.
  parts <- lapply(CHROMS, function(ch) {
    i <- which(bins$chr == ch)
    p <- lag_correlations(m[i, , drop = FALSE], bins[i], rao$domains[chrom == ch], rao$edges[chrom == ch, pos])
    write_table(p, file.path(OUT, "tables/pair_correlations", signal, paste0(ch, ".tsv.gz")))
    p[!is.na(class) & is.finite(r), .(lag, class, r)]
  })
  s <- boundary_summary(rbindlist(parts))
  grid <- CJ(lag = 1:100, class = c("within", "across"))
  s <- merge(grid, s, all.x = TRUE); s[is.na(n_pairs), n_pairs := 0L]
  delta <- dcast(s, lag ~ class, value.var = "mean")
  delta[, within_minus_across := within - across]
  s <- merge(s, delta[, .(lag, within_minus_across)], by = "lag")
  s[, `:=`(signal = signal, distance_bp = lag * 20000L)]
  all_summary[[signal]] <- s
}
write_table(rbindlist(all_summary), file.path(OUT, "tables/boundary_summary.tsv"))
