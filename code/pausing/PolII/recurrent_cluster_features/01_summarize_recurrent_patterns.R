#!/usr/bin/env Rscript

# Reconstruct existing read-cluster profiles 

# bash code/pausing/PolII/run_gene_examples.sh code/pausing/PolII/recurrent_cluster_features/01_summarize_recurrent_patterns.R


# Outputs:
#   recurrent_patterns/tables/gene_features.tsv       one row per input gene
#   recurrent_patterns/tables/<gene>/read_features.tsv.gz
#   recurrent_patterns/tables/<gene>/cluster_features.tsv
#   recurrent_patterns/tables/<gene>/cluster_profiles.tsv.gz
#   recurrent_patterns/plots/recurrent_patterns_summary.pdf

# Per-gene tables retain all original reads/clusters, including excluded clusters

# Criteria for clusters to be computed for cluster level feature clustering:
# a gene must have >= 50 reads, each cluster needs to have 20 reads,  <= 4	of the 40 × 50 bp bins, at most 4 may have no measurable m6A


# Family is NA for excluded clusters; their original cluster label is retained.

# Methods: fixed +/-1 kb TSS window, strand oriented, 40 bins (last includes +1000)
# m6A opportunities are sites observed among all full-span reads of that gene
# Nucleosomes retain the 130-160 bp class 

# Ward.D2 on m6A profiles Euclidean Ward linkage
# Families are labelled by increasing mean m6A
# Gene >=50 reads; cluster >=20 reads and <=4 unobserved bins; presence >=10%

# Family  fractions exclude ineligible clusters, so family fractions can sum to less than one. k=6 is descriptive;
# silhouettes for k=2..8 show whether discrete families are well separated.
#
# Gene table: original metadata/status plus pooled read-level signal means,
# family read fractions, retained cluster counts, dominant family,
# and gene-average 50 bp m6A/nucleosome profiles (bin01=-1000..-951, ..., bin40=
# +950..+1000). Gene-average profiles include all clusters and weight by reads

suppressPackageStartupMessages({ library(Matrix); library(data.table) })
setDTthreads(1L)
PROJECT <- "/project/spott/cshan/fiber-seq"
BASE <- file.path(PROJECT, "LCL_project/pausing/PolII")
OUT <- file.path(BASE, "recurrent_patterns")
source(file.path(PROJECT, "code/parsing_functions/parsing_footprints_functions.r"))
GROUPS <- c("housekeeping", "cell_type_specific", "not_expressed")
COLORS <- setNames(c("#2878B5", "#D77A21", "#777777"), GROUPS)
POS <- seq(-975, 975, 50)
METRICS <- c("m6a_total", "m6a_core", "m6a_upstream", "m6a_downstream", "m6a_far_up",
             "m6a_far_down", "nuc_at_tss", "nuc_core", "nuc_total")

############################################################
# one gene's saved read clusters rebuilt into profiles and features
############################################################
# 1. read the saved read-cluster table and re-assemble the gene's full-span m6A reads
# 2. orient by transcription and mark the sites observed in at least one read
# 3. paint the 130-160 bp fibertools nucleosome calls onto each read
# 4. per cluster x 50 bp bin: m6A at observed sites and nucleosome coverage
# 5. per read: regional m6A rates, nucleosome coverage and the +1 nucleosome
# 6. average those to cluster and gene features

extract_gene <- function(g, samples) {
  gene <- g$gene_name
  message("Reconstructing ", gene)
  a <- fread(file.path(OUT, "tables", g$group, paste0(gene, "_read_clusters.tsv")))
  # each row is a read, with read ID, sample name and cluster
  region <- data.frame(chr = g$chrom, start = g$tss - 1001L, end = g$tss + 1000L,
    analysis_start = g$tss - 1000L, analysis_end = g$tss + 1000L, region_id = gene)
  dat <- assemble_region_m6a(as.data.frame(samples), region, "m6a", full_span = TRUE,
    positions = seq(region$analysis_start, region$analysis_end))
  stopifnot(!anyDuplicated(a$RID), setequal(a$RID, dat$rids_df$RID),
            nrow(a) == g$n_reads, uniqueN(a$cluster) == g$n_clusters)
  mat <- as.matrix(dat$met_mat[a$RID, , drop = FALSE])
  
  # orient strand
  # +, offset = pos - tss
  # -, offset = tss - pos
  offset <- (as.integer(colnames(mat)) - g$tss) * if (g$strand == "+") 1L else -1L
  ord <- order(offset); offset <- offset[ord]; mat <- mat[, ord, drop = FALSE]
  observed <- colSums(mat) > 0
  nuc <- matrix(FALSE, nrow(a), length(offset))
  gr <- GenomicRanges::GRanges(g$chrom, IRanges::IRanges(region$analysis_start, region$analysis_end))
  for (sn in unique(a$sample_name)) {
    path <- file.path(samples[sample_name == sn, fire_dir], "extracted_results/nuc_by_chr",
                      paste0(sn, ".ft_extracted_nuc.", g$chrom, ".bed.gz"))
    bed <- read_ft_bed12(path, gr, longest_alignment = TRUE)
    if (!nrow(bed)) next
    b <- as.data.table(convert_ft_bed12_to_bed6(bed, format = "bed12_fibertools", source = path))
    b <- b[end - start >= 130 & end - start <= 160]
    b[, rid := match(paste(sn, RID, sep = "::"), a$RID)]
    b <- b[!is.na(rid)]
    if (g$strand == "+") b[, `:=`(lo = start + 1L - g$tss, hi = end - g$tss)] else
      b[, `:=`(lo = g$tss - end, hi = g$tss - start - 1L)]
    
    # convert each nucleosome into an offset interval
    b[, `:=`(lo = pmax(lo, -1000L), hi = pmin(hi, 1000L))]
    b <- b[lo <= hi]
    if (nrow(b)) {
      lens <- b$hi - b$lo + 1L
      nuc[cbind(rep(b$rid, lens), sequence(lens, from = b$lo + 1001L))] <- TRUE
    }
  }
  bin <- pmin((offset + 1000L) %/% 50L, 39L)
  profiles <- rbindlist(lapply(unique(a$cluster), function(cl) {
    ix <- which(a$cluster == cl)
    rbindlist(lapply(0:39, function(k) {
      sites <- bin == k & observed
      data.table(gene = gene, cluster = cl, bin = k, position = POS[k + 1L],
        n_reads = length(ix), n_observed_sites = sum(sites),
        m6a = if (any(sites)) mean(mat[ix, sites, drop = FALSE]) else NA_real_,
        nucleosome = mean(nuc[ix, bin == k, drop = FALSE]))
    }))
  }))
  rate <- function(lo, hi) {
    ix <- offset >= lo & offset <= hi & observed
    if (any(ix)) rowMeans(mat[, ix, drop = FALSE]) else rep(NA_real_, nrow(a))
  }
  plus <- nuc[, offset >= 0 & offset <= 299, drop = FALSE]
  first <- max.col(plus * 1L, ties.method = "first") - 1L
  first[rowSums(plus) == 0] <- NA_integer_
  a[, `:=`(group = g$group, gene = gene, m6a_total = rate(-1000, 1000),
    m6a_core = rate(-200, 200), m6a_upstream = rate(-500, -201),
    m6a_downstream = rate(201, 500), m6a_far_up = rate(-1000, -501),
    m6a_far_down = rate(501, 1000), nuc_at_tss = nuc[, offset == 0],
    nuc_core = rowMeans(nuc[, abs(offset) <= 200, drop = FALSE]),
    nuc_total = rowMeans(nuc), plus1 = first)]
  clusters <- a[, c(list(n_reads = .N, n_samples = uniqueN(sample_name)),
                     lapply(.SD, mean)), by = .(group, gene, cluster), .SDcols = METRICS]
  gene_features <- a[, c(list(gene = unique(gene), n_samples = uniqueN(sample_name),
    plus1_fraction = mean(!is.na(plus1)), plus1_median = median(plus1, na.rm = TRUE)),
    lapply(.SD, mean)), .SDcols = METRICS]
  list(clusters = clusters, profiles = profiles, genes = gene_features, reads = a)
}

############################################################
# cross-gene profile families and their recurrence
############################################################
# 1. pool every gene's clusters and set missing bins to 0
# 2. keep eligible clusters: gene >= 50 reads, cluster >= 20 reads, <= 4 unobserved bins
# 3. Ward.D2 on the 40-bin m6A profiles, cut at k = 6, relabel by increasing mean m6A
# 4. silhouettes for k = 2..8 as a separation check
# 5. flag low-signal, promoter-enriched and broad clusters
# 6. per gene: read fraction and presence per family, plus the dominant family
# 7. all of a gene's clusters are combined into one profile per gene, then those gene profiles are averaged into one profile per gene group
# 8. merge everything onto the manifest
summarize_patterns <- function(clusters, profiles, manifest, gene_features) {
  a <- copy(clusters)
  setorder(a, group, gene, cluster)
  a[, gene_reads := sum(n_reads), by = gene]
  
  # Each gene gets one 40-bin m6A profile and one 40-bin nucleosome profile. 
  # Each is the average of that gene's cluster profiles, with each cluster weighted by its share of the gene's reads 
  a[, read_fraction := n_reads / gene_reads]
  a[, core_enrichment := m6a_core - (m6a_far_up + m6a_far_down) / 2]
  profile_matrix <- function(value) {
    w <- dcast(profiles, gene + cluster ~ bin, value.var = value)
    key <- paste(w$gene, w$cluster, sep = "::")
    m <- as.matrix(w[match(paste(a$gene, a$cluster, sep = "::"), key), as.character(0:39), with = FALSE])
    stopifnot(nrow(m) == nrow(a), ncol(m) == 40L)
    m
  }
  stopifnot(!anyDuplicated(a[, .(gene, cluster)]),
            !anyDuplicated(profiles[, .(gene, cluster, bin)]), nrow(profiles) == 40L * nrow(a))
  x <- profile_matrix("m6a"); nuc <- profile_matrix("nucleosome")
  a[, missing_profile_bins := rowSums(!is.finite(x))]
  x[!is.finite(x)] <- 0
  a[, eligible := gene_reads >= 50 & n_reads >= 20 & missing_profile_bins <= 4]
  fit <- which(a$eligible)
  stopifnot(length(fit) > 8L, all(is.finite(x[fit, ])),
            all(x[fit, ] >= 0 & x[fit, ] <= 1), all(is.finite(nuc)), all(nuc >= 0 & nuc <= 1))
  distance <- dist(x[fit, ])
  tree <- hclust(distance, method = "ward.D2")
  sensitivity <- rbindlist(lapply(2:8, function(k) {
    label <- cutree(tree, k)
    data.table(k = k, mean_silhouette = mean(cluster::silhouette(label, distance)[, "sil_width"]))
  }))
  label <- cutree(tree, 6)
  centers <- t(vapply(1:6, function(k) colMeans(x[fit[label == k], , drop = FALSE]), numeric(40)))
  order_mean <- order(rowMeans(centers))
  centers <- centers[order_mean, , drop = FALSE]
  
  # Ineligible clusters are not assigned a family and never enter recurrence counts
  a[, family := NA_integer_]
  a[fit, family := match(label, order_mean)]
  a[, `:=`(low_signal = apply(x, 1, quantile, .9, na.rm = TRUE) < .25 &
                         apply(x, 1, max, na.rm = TRUE) < .40,
            promoter_enriched = core_enrichment >= .15 & m6a_core >= .30,
            broad_signal = rowMeans(x >= .30) >= .40)]
  completed <- manifest[status == "complete", .(group, gene = gene_name, n_reads)]
  fractions <- copy(completed)
  for (k in 1:6) {
    q <- a[eligible & family == k, .(fraction = sum(read_fraction)), by = gene]
    v <- q$fraction[match(fractions$gene, q$gene)]; v[is.na(v)] <- 0
    set(fractions, j = paste0("F", k, "_read_fraction"), value = v)
    set(fractions, j = paste0("F", k, "_present"), value = fractions$n_reads >= 50 & v >= .10)
  }
  for (pattern in c("low_signal", "promoter_enriched", "broad_signal")) {
    q <- a[eligible & get(pattern) == TRUE, .(fraction = sum(read_fraction)), by = gene]
    v <- q$fraction[match(fractions$gene, q$gene)]; v[is.na(v)] <- 0
    set(fractions, j = paste0(pattern, "_read_fraction"), value = v)
  }
  fractions[, gene_eligible := n_reads >= 50]
  fraction_cols <- paste0("F", 1:6, "_read_fraction")
  fm <- as.matrix(fractions[, ..fraction_cols])
  fractions[, `:=`(retained_read_fraction = rowSums(fm),
    dominant_family = fifelse(rowSums(fm) > 0, paste0("F", max.col(fm, ties.method = "first")), NA_character_))]
  count <- a[, .(n_eligible_clusters = sum(eligible)), by = gene]
  fractions <- merge(fractions, count, by = "gene", sort = FALSE)
  recurrence <- rbindlist(lapply(1:6, function(k) {
    fractions[gene_eligible == TRUE, .(gene_prevalence = mean(get(paste0("F", k, "_present"))),
      mean_gene_read_fraction = mean(get(paste0("F", k, "_read_fraction")))), by = group][, family := k]
  }))
  # Include every cluster when forming read-weighted gene profiles 
  gene_profiles <- rbindlist(lapply(completed$gene, function(gene_id) {
    ix <- which(a$gene == gene_id); w <- a$read_fraction[ix]
    data.table(gene = gene_id, position = POS,
      m6a = as.numeric(crossprod(w, x[ix, , drop = FALSE])),
      nuc = as.numeric(crossprod(w, nuc[ix, , drop = FALSE])))
  }))
  category <- merge(gene_profiles, completed[n_reads >= 50, .(gene, group)], by = "gene")
  category <- category[, .(m6a = mean(m6a), nuc = mean(nuc)), by = .(group, position)]
  features <- merge(gene_features, fractions[, !c("group", "n_reads"), with = FALSE], by = "gene", sort = FALSE)
  for (signal in c("m6a", "nuc")) {
    wide <- dcast(gene_profiles, gene ~ position, value.var = signal)
    setnames(wide, as.character(POS), sprintf("%s_bin%02d", signal, 1:40))
    features <- merge(features, wide, by = "gene", sort = FALSE)
  }
  features <- merge(manifest, features, by.x = "gene_name", by.y = "gene", all.x = TRUE, sort = FALSE)
  features <- features[match(manifest$gene_name, gene_name)]
  features[, reconstruction_status := fifelse(status == "complete", "complete", "not_run_input_failed")]
  stopifnot(nrow(features) == nrow(manifest), !anyDuplicated(features$gene_name),
    all(fractions$retained_read_fraction <= 1 + 1e-12))
  list(features = features, clusters = a, x = x, nuc = nuc, centers = centers,
       recurrence = recurrence, category = category, sensitivity = sensitivity)
}

############################################################
# the four-page summary figure
############################################################
# 1. the six family profiles with their nucleosome coverage
# 2. family prevalence and mean read fraction per gene category
# 3. category mean m6A and nucleosome profiles
# 4. silhouette curve and core m6A against core nucleosome per cluster
plot_summary <- function(s, path) {
  grDevices::pdf(path, width = 12, height = 7, onefile = TRUE)
  on.exit(grDevices::dev.off())
  frame <- function(title, ylab = "Fraction") {
    plot(NA, xlim = c(-1000, 1000), ylim = c(0, 1), main = title,
         xlab = "bp from TSS", ylab = ylab, bty = "l")
    abline(v = 0, col = "grey60", lty = 2)
  }
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (k in 1:6) {
    ix <- which(s$clusters$eligible & s$clusters$family == k)
    frame(sprintf("F%d: %d clusters", k, length(ix)))
    matlines(POS, t(s$x[ix, , drop = FALSE]), col = adjustcolor("grey50", .045), lty = 1, lwd = .6)
    lines(POS, s$centers[k, ], lwd = 2)
    lines(POS, colMeans(s$nuc[ix, , drop = FALSE]), col = "#2356a0", lwd = 2)
    if (k == 1) legend("topright", c("m6A", "130-160 bp nucleosome coverage"),
      col = c("black", "#2356a0"), lty = 1, lwd = 2, cex = .65, bty = "n")
  }
  mtext("Exploratory profile families: absolute m6A patterns, Ward clustering (k=6)", outer = TRUE)
  par(mfrow = c(1, 2), mar = c(4, 4, 4, 1), oma = c(0, 0, 3, 0))
  for (metric in c("gene_prevalence", "mean_gene_read_fraction")) {
    m <- t(vapply(GROUPS, function(g) {
      q <- s$recurrence[group == g]; q[[metric]][match(1:6, q$family)]
    }, numeric(6)))
    barplot(m, beside = TRUE, col = COLORS, names.arg = paste0("F", 1:6),
      ylim = c(0, 1), border = NA, ylab = "Fraction",
      main = if (metric == "gene_prevalence") "Genes with >=10% reads in family" else "Mean read fraction per gene")
    if (metric == "mean_gene_read_fraction") legend("topright", GROUPS, fill = COLORS, bty = "n", cex = .8)
  }
  mtext("Genes >=50 reads; clusters >=20 reads; genes weighted equally", outer = TRUE)
  par(mfrow = c(1, 2), mar = c(4, 4, 4, 1), oma = c(0, 0, 3, 0))
  for (signal in c("m6a", "nuc")) {
    frame(if (signal == "m6a") "m6A rate at observed sites" else "130-160 bp nucleosome coverage")
    for (g in GROUPS) {
      d <- s$category[group == g][order(position)]
      lines(d$position, d[[signal]], col = COLORS[[g]], lwd = 2)
    }
    if (signal == "m6a") legend("topright", GROUPS, col = COLORS, lty = 1, lwd = 2, bty = "n", cex = .8)
  }
  mtext("Category means: equal weight per gene", outer = TRUE)
  par(mfrow = c(1, 2), mar = c(4, 4, 4, 1), oma = c(0, 0, 2, 0))
  plot(s$sensitivity$k, s$sensitivity$mean_silhouette, type = "b", pch = 16, bty = "l",
       xlab = "Number of profile families", ylab = "Mean silhouette", main = "Are families well separated?")
  d <- s$clusters[eligible == TRUE]
  plot(d$m6a_core, d$nuc_core, col = adjustcolor(COLORS[d$group], .35), pch = 16, cex = .6,
       xlim = c(0, 1), ylim = c(0, 1), bty = "l", xlab = "Core m6A rate (-200..+200)",
       ylab = "Core nucleosome coverage", main = "Corresponding nucleosome footprints")
  legend("topright", GROUPS, col = COLORS, pch = 16, bty = "n", cex = .8)
}

############################################################
# run every completed gene, then write the tables and the PDF
############################################################
# 1. read the 31-sample metatable and the gene manifest
# 2. reconstruct the completed genes in parallel, raising any worker error
# 3. summarize_patterns() over all of them
# 4. write gene_features.tsv and the per-gene read, cluster and profile tables
# 5. draw recurrent_patterns_summary.pdf
main <- function() {
  samples <- fread("/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv")
  manifest <- fread(file.path(OUT, "tables/example_genes.tsv"))
  stopifnot(!anyDuplicated(manifest$gene_name), !anyDuplicated(samples$sample_name),
            nrow(samples) == 31L, all(manifest$group %in% GROUPS))
  workers <- as.integer(Sys.getenv("RECURRENT_WORKERS", Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
  stopifnot(length(workers) == 1L, !is.na(workers), workers >= 1L)
  idx <- which(manifest$status == "complete")
  result <- parallel::mclapply(idx, function(i) {
    tryCatch(extract_gene(manifest[i], samples), error = function(e)
      list(error = paste(manifest$gene_name[i], conditionMessage(e), sep = ": ")))
  }, mc.cores = workers, mc.preschedule = FALSE)
  failed <- vapply(result, function(x) inherits(x, "try-error") || !is.null(x$error), logical(1))
  if (any(failed)) stop(paste(vapply(result[failed], function(x)
    if (inherits(x, "try-error")) as.character(x) else x$error, character(1)), collapse = "\n"))
  s <- summarize_patterns(rbindlist(lapply(result, `[[`, "clusters")),
    rbindlist(lapply(result, `[[`, "profiles")), manifest, rbindlist(lapply(result, `[[`, "genes")))
  dir.create(file.path(OUT, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(OUT, "plots"), recursive = TRUE, showWarnings = FALSE)
  fwrite(s$features, file.path(OUT, "tables/gene_features.tsv"), sep = "\t", na = "NA")
  for (r in result) {
    gene_id <- unique(r$profiles$gene)
    stopifnot(length(gene_id) == 1L, !is.na(gene_id))
    folder <- file.path(OUT, "tables", gene_id)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    cl <- copy(s$clusters[gene == gene_id])
    cl[, family := fifelse(is.na(family), NA_character_, paste0("F", family))]
    reads <- merge(r$reads, cl[, .(cluster, family, eligible)], by = "cluster", all.x = TRUE, sort = FALSE)
    setorder(reads, clustered_row)
    pp <- merge(r$profiles, cl[, .(cluster, family, eligible)], by = "cluster", all.x = TRUE, sort = FALSE)
    setorder(pp, cluster, bin)
    fwrite(reads, file.path(folder, "read_features.tsv.gz"), sep = "\t", na = "NA")
    fwrite(cl, file.path(folder, "cluster_features.tsv"), sep = "\t", na = "NA")
    fwrite(pp, file.path(folder, "cluster_profiles.tsv.gz"), sep = "\t", na = "NA")
  }
  plot_summary(s, file.path(OUT, "plots/recurrent_patterns_summary.pdf"))
  message("Complete: ", length(idx), " reconstructed genes; ", nrow(s$clusters), " original clusters; ",
          sum(s$clusters$eligible), " eligible clusters; ", nrow(s$features), " gene table rows.")
}
if (sys.nframe() == 0L) main()
