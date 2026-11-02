# Small fixtures exercise selection/denominator semantics, not a second clustering implementation.
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1L])
here <- dirname(normalizePath(script))
source(file.path(here, "pol2_clustering_functions.R"))
project <- normalizePath(file.path(here, "../../.."))
pol2_load_helpers(project)
site <- pol2_windows(data.table(chrom = "chr1", tss = 2001L, strand = "+", tss_id = "plus",
                                gene_id = "ENSGTEST.1", gene_name = "TEST", pol2_tss_signal = 4))
minus <- pol2_windows(data.table(chrom = "chr1", tss = 2001L, strand = "-", tss_id = "minus",
                                 gene_id = "ENSGTEST.1", gene_name = "TEST", pol2_tss_signal = 4))
stopifnot(site$window_start == 1000, site$window_end == 3000,
          minus$window_start == 1001, minus$window_end == 3001)
# TSS region -50..+299 (1-based 1951..2300 on "+", 1702..2051 on "-"), as BED.
tssr <- pol2_tssr(rbindlist(list(site, minus)))
stopifnot(identical(tssr$start, c(1950L, 1701L)), identical(tssr$end, c(2300L, 2051L)))

# 40 sites in an eligible molecule, zero in another, and 40 different sites in a
# nonqualifying full-window molecule. All 80 sites define the denominator;
# only the two qualifying molecules can be clustering rows.
make_bed <- function(id, positions, sample = "s1", start0 = 500L, end0 = 4000L) {
  positions <- sort(positions)
  data.table(chr = "chr1", start = start0, end = end0, RID = id, score = 60, strand = "+",
             x1 = start0, x2 = end0, rgb = "0", blockCount = length(positions) + 2L,
             blockSizes = paste(rep(1L, length(positions) + 2L), collapse = ","),
             blockStarts = paste(c(0L, positions - start0 - 1L, end0 - start0 - 1L), collapse = ","),
             original_RID = id, sample_name = sample, cell_line = sample)
}
relative <- seq(-975L, 975L, by = 50L)
bed <- rbindlist(list(make_bed("positive", site$tss + relative), make_bed("zero", integer()),
                     make_bed("excluded", site$tss + relative + 1L)))
candidates <- data.table(original_RID = c("positive", "zero"), sample_name = "s1", tss_id = "plus",
                         chrom = "chr1", start = 2020L, end = 2060L)
x <- pol2_locus_features(bed, candidates, site)
stopifnot(nrow(x$features) == 2, ncol(x$features) == 40,
          all(x$denominators$n_sites == 2), !"excluded" %in% x$molecules$original_RID,
          identical(x$molecules$m6a_proportion, c(0.5, 0)),
          all(x$features[1, ] == 0.5), all(x$features[2, ] == 0), all(x$molecules$pol2_tss_signal == 4))
# Per-bp calls keep the molecule's own calls on the -1000..999 axis.
stopifnot(nrow(x$sites) == 40L, all(x$sites$i == 1L), identical(sort(x$sites$j), relative + 1001L))

# Mirrored positions at a minus-strand TSS produce the same features.
bed_minus <- rbindlist(list(make_bed("positive", minus$tss - relative), make_bed("zero", integer()),
                          make_bed("excluded", minus$tss - relative - 1L)))
y <- pol2_locus_features(bed_minus, candidates, minus)
stopifnot(identical(unname(x$features), unname(y$features)), identical(sort(y$sites$j), sort(x$sites$j)))
# An empty bin in the entire full-window pool does not become a fabricated zero.
empty_bin <- make_bed("positive", site$tss + relative[-1L])
stopifnot(is.null(pol2_locus_features(empty_bin, candidates, site)))
# Other full-window molecules support a zero-call eligible molecule.
zero <- pol2_locus_features(bed, candidates[original_RID == "zero"], site)
stopifnot(nrow(zero$features) == 1, all(zero$features == 0),
          all(zero$denominators$n_sites == 2))
# Partial-span molecule cannot affect features or denominators.
partial <- make_bed("positive", site$tss + relative[-1L], start0 = 1100L)
stopifnot(is.null(pol2_locus_features(partial, candidates, site)))
partial_extra <- make_bed("partial", site$tss + relative[-1L] + 2L, start0 = 1100L)
with_partial <- pol2_locus_features(rbindlist(list(bed, partial_extra)), candidates, site)
stopifnot(identical(x$features, with_partial$features),
          identical(x$denominators, with_partial$denominators))
# A footprint qualifies only the read of its own sample, even if a name repeats.
two <- rbindlist(list(bed, make_bed("positive", site$tss + relative, sample = "s2")))
in_s2 <- copy(candidates)[, sample_name := "s2"]
own <- pol2_locus_features(two, in_s2[original_RID == "positive"], site)
stopifnot(nrow(own$molecules) == 1L, own$molecules$sample_name == "s2")
stopifnot(is.null(pol2_locus_features(bed, in_s2, site)))

# Size rule on actual intervals; peaks are no longer a footprint requirement.
fp <- data.table(chrom = "chr1", start = c(2020L, 1940L, 2300L, 2020L, 2020L, 2020L, 2020L),
                 end = c(2060L, 1980L, 2340L, 2049L, 2081L, 2101L, 2080L), original_RID = letters[1:7],
                 sample_name = "s1")
fp <- rbind(fp, fp[original_RID == "g"])
qual <- pol2_candidates(rbindlist(list(site, minus)), fp)
stopifnot(nrow(qual) == 8L, setequal(qual$original_RID, c("a", "b", "c", "g")),
          qual[tss_id == "plus" & original_RID == "b", midpoint_relative] < 0,
          qual[tss_id == "minus" & original_RID == "a", midpoint_relative] < 0,
          nrow(qual[original_RID == "g"]) == 2L,
          all(qual[original_RID == "g", footprint_class] == "30-60 bp"))
stopifnot(!nrow(pol2_candidates(site, fp[, `:=`(start = 5000L, end = 5040L)])))

# Signal averages over covered bases and retains the covered fraction.
tmp <- tempfile(fileext = ".bw")
gr <- GRanges("chr1", IRanges(c(2021L, 2031L), c(2030L, 2040L)), score = c(2, 4))
GenomeInfoDb::seqlengths(gr) <- c(chr1 = 5000L)
rtracklayer::export(gr, tmp, format = "BigWig")
sig <- pol2_signal(data.table(chrom = "chr1", start = c(2020L, 3000L), end = c(2060L, 3040L)), tmp)
stopifnot(sig$pol2_mean[1] == 3, sig$pol2_covered_fraction[1] == 0.5,
          is.na(sig$pol2_mean[2]), sig$pol2_covered_fraction[2] == 0)
unlink(tmp)

# PI joins by Ensembl gene ID (versions ignored); anchors are not compared.
path <- tempfile(fileext = ".tsv")
fwrite(data.table(ensg = c("ENSGTEST", "ENSGOTHER"), tss = 9999L, pi = c(3, NA)), path, sep = "\t")
pi <- pol2_attach_pi(x$molecules, path)
stopifnot(all(pi$pi == 3), all(pi$pi_status == "matched"))
other <- copy(x$molecules)[, gene_id := "ENSGOTHER.1"]
stopifnot(all(pol2_attach_pi(other, path)$pi_status == "no_pi"))
unlink(path)
stopifnot(all(pol2_attach_pi(x$molecules, NULL)$pi_status == "not_supplied"))

# TSS signal ratios point the right way and skip missing coverage; clusters are
# ordered by median PI over their TSSs (each TSS once), clusters without PI last.
m <- data.table(cluster = rep(c("cluster1", "cluster2", "cluster3"), c(40L, 40L, 20L)),
                tss_id = c(rep(c("t1", "t2"), 20L), paste0("u", 1:40), paste0("v", 1:20)),
                pol2_tss_signal = c(seq(20, 30, length.out = 40L), seq(1, 3, length.out = 40L),
                                    rep(c(10, NA), 10L)),
                pi = c(rep(c(50, 1), 20L), rep(5, 40L), rep(NA, 20L)))
s <- pol2_cluster_summary(m)
overall <- mean(m$pol2_tss_signal, na.rm = TRUE)
stopifnot(identical(s$cluster, c("cluster2", "cluster1", "cluster3")),
          s[cluster == "cluster1", log2_ratio] > 0, s[cluster == "cluster2", log2_ratio] < 0,
          all.equal(s[cluster == "cluster3", log2_ratio], log2(10 / overall)),
          s[cluster == "cluster3", n_with_signal] == 10L, s[cluster == "cluster3", n] == 20L,
          s[cluster == "cluster1", median_pi] == 25.5, s[cluster == "cluster1", n_tss_with_pi] == 2L,
          s[cluster == "cluster3", n_tss_with_pi] == 0L, s[cluster == "cluster1", signif] == "***")

# Global sampling selects unique physical molecules, not read/TSS pairs.
pool <- data.table(RID = rep(paste0("sample::r", seq_len(6000L)), each = 2L))
selected <- pol2_sample_rows(pool, 5000L, 1L)
stopifnot(length(selected) == 5000L, uniqueN(pool$RID[selected]) == 5000L,
          identical(selected, pol2_sample_rows(pool, 5000L, 1L)),
          length(pol2_sample_rows(pool[1:20], 5000L, 1L)) == 10L)

message("All footprint-selection, orientation, denominator, signal, PI, TSS-signal and random-sampling fixtures passed")
