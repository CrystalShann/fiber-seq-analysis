#!/usr/bin/env Rscript
# Li et al. (2026), doi:10.1093/nar/gkaf1385, adapted to Fiber-seq occupancy.
# Run: bash code/entropy/nucleosome/run_entropy.sh [--by-cluster | --gene-only] [gene names]
# Refresh footprint fractions without recalculating entropy: --footprints-only
# Uses the existing PolII full-span promoter cohort and resolved TSSs.
# Entropy is calculated from 500 nucleotide occupancy fractions, never from bins.

suppressPackageStartupMessages(library(data.table))
PROJECT <- "/project/spott/cshan/fiber-seq"
BASE <- file.path(PROJECT, "LCL_project/Pol2/PolII")
SAMPLES <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv"
# This also sources the shared parsing functions; existing files stay unchanged.
source(file.path(PROJECT, "code/nucleosome_positioning/process_nuc_pos.R"), local = TRUE)
setDTthreads(1L)

SHORT_FP_ROOT <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged"
SHORT_FP_BREAKS <- c(10L, 30L, 60L, 80L)
SHORT_FP_WINDOW <- c(-100L, 100L) # inclusive, 201 bases; entropy stays at +1..+500
# Genes need > MIN_READS saved reads and > MIN_CLUSTERS Leiden clusters for entropy
MIN_READS <- 20L
MIN_CLUSTERS <- 1L

#########
# Convert FiberHMM v1 BED10 rows to one row per footprint
#######
# 1. Rebuild BED12 columns with neutral score/strand
# 2. Expand blocks with the shared converter, keeping first/last blocks
# Legacy FiberHMM v1 is BED10: chr,start,end,RID,thickStart,thickEnd,rgb,
# blockCount,blockStarts,blockSizes. Supply neutral score/strand fields to the
# shared BED12 converter. These are real blocks, including the first and last;
# unlike fibertools BED12, there are no sentinel blocks to remove.
decode_short_footprints <- function(raw) {
  stopifnot(ncol(raw) == 10L)
  bed <- data.frame(chr = raw[[1]], start = raw[[2]], end = raw[[3]], RID = raw[[4]],
    score = 0, strand = ".", x1 = raw[[5]], x2 = raw[[6]], rgb = raw[[7]],
    blockCount = raw[[8]], blockSizes = raw[[10]], blockStarts = raw[[9]])
  as.data.table(convert_ft_bed12_to_bed6(bed, drop_sentinels = FALSE, validate_blocks = TRUE))
}

#########
# Load FiberHMM footprints of cohort fibers around a TSS
#######
# 1. Build the genomic region of the TSS window
# 2. Per sample, read the tabix footprint file in that region
# 3. Keep exactly the saved cohort fibers (error if missing or duplicated)
# 4. Expand to footprints and assign size classes
load_short_footprints <- function(g, reads, window = SHORT_FP_WINDOW, root = SHORT_FP_ROOT) {
  stopifnot(length(window) == 2L, window[1] <= window[2], g$strand %in% c("+", "-"))
  endpoints <- g$tss + if (g$strand == "+") window else -window
  region <- GenomicRanges::GRanges(g$chrom, IRanges::IRanges(min(endpoints), max(endpoints)))
  parts <- lapply(unique(reads$sample_name), function(sn) {
    path <- file.path(root, sn, "joint_trained_results", paste0(sn, "_", g$chrom, "_fp.bed.gz"))
    if (!all(file.exists(c(path, paste0(path, ".tbi"))))) stop("Missing footprint file/index: ", path)
    raw <- read_tabix_region(path, region)
    wanted <- reads[sample_name == sn, RID]
    if (!nrow(raw)) stop("Missing footprint records for cohort: ", sn)
    stopifnot(ncol(raw) == 10L)
    raw <- unique(raw[paste(sn, raw[[4]], sep = "::") %in% wanted, , drop = FALSE])
    ids <- paste(sn, raw[[4]], sep = "::")
    if (!setequal(ids, wanted) || anyDuplicated(ids))
      stop("Missing or ambiguous footprint records for cohort: ", sn)
    b <- decode_short_footprints(raw)
    b[, RID := paste(sn, RID, sep = "::")]
    b[, footprint_class := assign_size_class(end - start, SHORT_FP_BREAKS)]
    b[!is.na(footprint_class), .(RID, start, end, footprint_class)]
  })
  if (!length(parts)) return(data.table(RID = character(), start = integer(), end = integer(),
                                       footprint_class = character()))
  unique(rbindlist(parts))
}

#########
# Fraction of fibers with a footprint in the TSS window, per size class
#######
# 1. Keep fibers spanning the whole window
# 2. Keep footprints fully inside the window and assign size classes
# 3. Count each fiber once per class and divide by spanning fibers
# Each retained fiber contributes at most one success per size class. A footprint
# qualifies only if all its bases lie within the inclusive window. Denominator:
# fibers spanning the ENTIRE footprint window, whether or not they have a call.
short_footprint_fraction <- function(reads, footprints, tss, strand, window = SHORT_FP_WINDOW) {
  stopifnot(!anyDuplicated(reads$RID), strand %in% c("+", "-"), length(window) == 2L,
            window[1] <= window[2])
  r <- rel_interval(reads$start, reads$end, tss, strand)
  eligible <- reads$RID[r$rel_start <= window[1] & r$rel_end >= window[2]]
  f <- copy(footprints)
  f[, c("lo", "hi") := rel_interval(start + 1L, end, tss, strand)]
  f[, footprint_class := assign_size_class(end - start, SHORT_FP_BREAKS)]
  hits <- unique(f[RID %in% eligible & lo >= window[1] & hi <= window[2] & !is.na(footprint_class),
                   .(RID, footprint_class)])
  classes <- levels(assign_size_class(integer(), SHORT_FP_BREAKS))
  data.table(footprint_class = classes, n_fibers = nrow(reads), n_evaluable_fibers = length(eligible),
    n_fibers_with_footprint = vapply(classes, function(cl) sum(hits$footprint_class == cl), integer(1)))[
      , fraction_with_footprint := if (length(eligible)) n_fibers_with_footprint / length(eligible) else NA_real_][]
}

#########
# Footprint fractions for all genes, written to a table
#######
# 1. Per gene with complete entropy, load footprints and compute fractions
# 2. Record status and errors for the other genes
# 3. Write short_footprint_fractions.tsv; stop if any gene failed
write_short_footprint_fractions <- function(entropy_genes, reads, table_dir, workers = 1L) {
  classes <- levels(assign_size_class(integer(), SHORT_FP_BREAKS))
  run_one <- function(i) {
    g <- entropy_genes[i]
    message("Short footprints: ", g$gene_name)
    a <- copy(reads[gene == g$gene_name])
    result <- data.table(footprint_class = classes, n_fibers = nrow(a),
      n_evaluable_fibers = NA_integer_, n_fibers_with_footprint = NA_integer_, fraction_with_footprint = NA_real_)
    status <- "entropy_unavailable"; error <- "Gene has no valid entropy result"
    if (g$status == "complete") tryCatch({
      stopifnot(nrow(a) == g$n_fibers, all(a$tss == g$tss), all(a$chrom == g$chrom), all(a$strand == g$strand))
      fp <- load_short_footprints(g, a)
      result <- short_footprint_fraction(a, fp, g$tss, g$strand)
      status <- if (all(result$n_evaluable_fibers > 0)) "complete" else "no_full_window_fibers"
      error <- ""
    }, error = function(e) { status <<- "error"; error <<- conditionMessage(e) })
    cbind(data.table(gene_name = g$gene_name, gene_id = g$gene_id, group = g$group,
      chrom = g$chrom, tss = g$tss, strand = g$strand, entropy_status = g$status,
      nucleosome_entropy_bits = g$nucleosome_entropy_bits,
      footprint_window_start_offset = SHORT_FP_WINDOW[1], footprint_window_end_offset = SHORT_FP_WINDOW[2],
      footprint_window_bp = diff(SHORT_FP_WINDOW) + 1L, footprint_inclusion = "fully_contained", footprint_status = status, error), result)
  }
  result <- parallel::mclapply(seq_len(nrow(entropy_genes)), run_one, mc.cores = workers)
  if (any(vapply(result, inherits, logical(1), "try-error"))) stop("Footprint worker failed")
  result <- rbindlist(result)
  stopifnot(!anyDuplicated(result[, .(gene_name, footprint_class)]))
  fwrite(result, file.path(table_dir, "short_footprint_fractions.tsv"), sep = "\t", na = "NA")
  failures <- unique(result[footprint_status == "error", .(gene_name, error)])
  if (nrow(failures)) {
    print(failures)
    stop("Footprint extraction failed for ", nrow(failures), " genes; see short_footprint_fractions.tsv")
  }
  invisible(result)
}

#########
# Load saved Leiden cluster labels for one gene
#######
# 1. Find the per-gene label file
# 2. Check the labels match the saved fiber cohort
# 3. Return RID, cluster and eligible in cohort order
# Use the saved per-gene Leiden labels.
entropy_cluster_labels <- function(g, reads, cluster_root) {
  path <- file.path(cluster_root, g$gene_name, "read_features.tsv.gz")
  if (!file.exists(path)) path <- file.path(cluster_root, g$group, paste0(g$gene_name, "_read_clusters.tsv"))
  if (!file.exists(path)) stop("No saved Leiden cluster labels for ", g$gene_name)
  labels <- fread(path)
  if (!all(c("RID", "cluster") %in% names(labels))) stop("Cluster labels lack RID/cluster")
  if (!"eligible" %in% names(labels)) labels[, eligible := NA]
  if (anyDuplicated(labels$RID) || !setequal(labels$RID, reads$RID) || anyNA(labels$cluster))
    stop("Cluster labels do not match saved fiber cohort")
  labels[match(reads$RID, RID), .(RID, cluster, eligible)]
}

#########
# Nucleosome entropy from per-base counts over +1..+500
#######
# 1. Occupancy = nucleosome fibers / covering fibers at each base
# 2. Set status: no fibers, incomplete coverage, zero occupancy or complete
# 3. If complete, normalize occupancy to probabilities and compute H in bits
# 4. Return total, entropy, status and the per-base profile
# Missing coverage is unknown, not zero occupancy. Require all 500 positions to
# be observed; otherwise a smaller effective window would change the estimand.
entropy_profile <- function(nucleosome_fibers, covering_fibers) {
  stopifnot(length(nucleosome_fibers) == 500L, length(covering_fibers) == 500L,
            all(is.finite(nucleosome_fibers)), all(is.finite(covering_fibers)),
            all(nucleosome_fibers >= 0), all(covering_fibers >= nucleosome_fibers))
  occupancy <- rep(NA_real_, 500L)
  observed <- covering_fibers > 0
  occupancy[observed] <- nucleosome_fibers[observed] / covering_fibers[observed]
  total <- if (all(observed)) sum(occupancy) else NA_real_
  status <- if (!any(observed)) "no_fibers" else if (!all(observed)) "incomplete_coverage" else
    if (total == 0) "zero_nucleosome_occupancy" else "complete"
  p <- rep(NA_real_, 500L)
  h <- NA_real_
  if (status == "complete") {
    p <- occupancy / total
    h <- -sum(p[p > 0] * log2(p[p > 0]))
  }
  list(total = total, entropy = h, status = status,
       profile = data.table(offset = 1:500, covering_fibers, nucleosome_fibers,
                            occupancy, probability = p))
}

#########
# Per-base nucleosome occupancy and entropy for a fiber set
#######
# 1. Count covering and occupied fibers per base with occupancy_profile
# 2. Fill missing positions with zero depth
# 3. Compute entropy with entropy_profile
# Reuse the existing occupancy engine at 1 bp resolution. It counts each fiber
# once per base and divides by fibers whose alignments cover that base. Explicitly
# restore absent positions to depth zero so gaps cannot silently shrink the window.
nucleotide_profile <- function(footprints, reads, anchors) {
  counts <- data.table(pos = 1:500, n_reads = 0L, n_occupied = 0L)
  if (nrow(reads)) {
    p <- occupancy_profile(footprints, reads, anchors, by = "region_id",
                           bin_size = 1L, metric = "occupancy")
    counts[p, on = "pos", `:=`(n_reads = i.n_reads, n_occupied = i.n_occupied)]
  }
  entropy_profile(counts$n_occupied, counts$n_reads)
}

#########
# Load 130-160 bp nucleosome calls of cohort fibers
#######
# 1. Per sample, read ft nuc_by_chr records in TSS +/- 1000 bp
# 2. Keep records matching the saved cohort reads and alignments
# 3. Expand blocks and keep 130-160 bp nucleosomes
# 4. Return nucleosomes with sample::RID
# Query original calls, retaining the exact alignment used for the saved cohort.
# A saved read without a BED record (fibertools writes none for reads with almost
# no m6A) counts as a covering fiber with no nucleosomes, like a present record
# with no retained blocks. Sentinel blocks are removed by the shared parser
# before the existing inclusive 130-160 bp filter is applied.
load_nucleosomes <- function(g, reads, samples) {
  region <- GenomicRanges::GRanges(g$chrom,
    IRanges::IRanges(g$tss - 1000L, g$tss + 1000L))
  out <- lapply(unique(reads$sample_name), function(sn) {
    root <- samples[sample_name == sn, fire_dir]
    if (length(root) != 1L) stop("Missing or duplicate sample metadata: ", sn)
    path <- file.path(root, "extracted_results/nuc_by_chr",
                      paste0(sn, ".ft_extracted_nuc.", g$chrom, ".bed.gz"))
    b <- as.data.table(read_ft_bed12(path, region, longest_alignment = FALSE))
    wanted <- reads[sample_name == sn]
    if (!nrow(b)) return(NULL)
    key <- paste(sn, b$RID, sep = "::")
    ix <- match(key, wanted$RID)
    b <- unique(b[which(!is.na(ix) & b$start + 1L == wanted$start[ix] & b$end == wanted$end[ix])])
    if (anyDuplicated(b$RID))
      stop("Ambiguous nucleosome alignment for saved cohort: ", sn)
    if (!nrow(b)) return(NULL)
    blocks <- as.data.table(convert_ft_bed12_to_bed6(b, format = "bed12_fibertools", source = path))
    blocks <- blocks[!is.na(assign_size_class(end - start, c(130L, 160L)))]
    if (nrow(blocks) && any(blocks$start < b$start[match(blocks$RID, b$RID)] |
                          blocks$end > b$end[match(blocks$RID, b$RID)]))
      stop("Nucleosome call outside its alignment: ", path)
    blocks[, .(RID = paste(sn, RID, sep = "::"), start, end)]
  })
  out <- rbindlist(out)
  if (!nrow(out)) return(data.table(RID = character(), start = integer(), end = integer()))
  out
}

#########
# Mean of finite values
#######
# 1. Drop non-finite values
# 2. Return their mean, or NA if none are left
finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
#########
# Standard deviation of finite values
#######
# 1. Drop non-finite values
# 2. Return their sd, or NA if fewer than 2 are left
finite_sd <- function(x) if (sum(is.finite(x)) >= 2L) sd(x[is.finite(x)]) else NA_real_

#########
# Entropy and nucleosome summary for a subset of fibers
#######
# 1. Compute the per-base occupancy profile and entropy of the subset
# 2. Add coverage, +1 nucleosome, m6A and TSS-occupancy summaries
# 3. Return one summary row and the per-base profile
summarize_subset <- function(g, reads, footprints, anchors, ix, level, cluster = NA_character_,
                             eligible = NA) {
  a <- reads[ix]
  result <- nucleotide_profile(footprints[RID %in% a$RID], a, anchors)
  depth <- result$profile$covering_fibers
  positions <- g$tss + if (g$strand == "+") 1:500 else -(1:500)
  id <- data.table(gene_name = g$gene_name, gene_id = g$gene_id, group = g$group,
                   level, cluster, eligible)
  row <- cbind(id, data.table(chrom = g$chrom, tss = g$tss, strand = g$strand,
    tss_source = g$tss_source, window_start = min(positions),
    window_end = max(positions), window_bp = 500L,
    n_fibers = sum(a$rel_read_start <= 500L & a$rel_read_end >= 1L), n_samples = uniqueN(a$sample_name),
    n_full_window_fibers = sum(a$rel_read_start <= 1L & a$rel_read_end >= 500L),
    min_coverage = min(depth), mean_coverage = mean(depth),
    max_coverage = max(depth), n_uncovered_positions = sum(depth == 0),
    total_nucleosome_occupancy_signal = result$total,
    nucleosome_entropy_bits = result$entropy,
    normalized_nucleosome_entropy = result$entropy / log2(500), # H / maximum entropy over 500 bases
    status = result$status, error = "",
    n_plus1_fibers = sum(is.finite(a$plus1_start)),
    plus1_start_mean_bp = finite_mean(a$plus1_start),
    plus1_start_sd_bp = finite_sd(a$plus1_start),
    n_plus1_nearest_midpoint_fibers = sum(is.finite(a$plus1_mid)),
    plus1_nearest_midpoint_sd_bp = finite_sd(a$plus1_mid),
    promoter_m6a_core = finite_mean(a$m6a_core),
    tss_nucleosome_free_fraction = finite_mean(!a$nuc_at_tss),
    tpm_rep1 = g$tpm_rep1, tpm_rep2 = g$tpm_rep2, tpm_mean = g$tpm_mean))
  profile <- cbind(id[rep(1L, 500L)], result$profile,
                   genomic_position = positions)
  list(summary = row, profile = profile)
}

#########
# Gene- and cluster-level nucleosome entropy for one gene
#######
# 1. Check TSS and saved read metadata; return NA if below the read/cluster filter
# 2. Load nucleosomes and convert them to TSS-relative coordinates
# 3. Check re-extracted +1 and TSS occupancy against saved features
# 4. Summarize all fibers, then each Leiden cluster
calculate_gene <- function(g, reads, samples, cluster_root, by_cluster = TRUE) {
  if (is.na(g$tss) || !g$strand %in% c("+", "-")) stop("No resolved TSS")
  if (nrow(reads) && (any(reads$chrom != g$chrom | reads$tss != g$tss | reads$strand != g$strand) ||
                     anyDuplicated(reads$RID))) stop("Inconsistent saved read metadata")
  if (nrow(reads) && (nrow(reads) <= MIN_READS ||
                      uniqueN(entropy_cluster_labels(g, reads, cluster_root)$cluster) <= MIN_CLUSTERS))
    return(list(summary = cbind(g[, .(gene_name, gene_id, group, chrom, tss, strand, tss_source)],
      data.table(level = "gene", n_fibers = nrow(reads), nucleosome_entropy_bits = NA_real_,
                 status = "below_min_reads_or_clusters", error = "")), profile = NULL))
  nucs <- load_nucleosomes(g, reads, samples)
  anchors <- build_region_anchors(NULL, NULL, NULL, relative_window = c(1L, 500L),
    anchor_table = data.table(region_id = g$gene_name, chr = g$chrom, anchor = g$tss,
                              strand = g$strand, region_type = "promoter", gene = g$gene_name, ensg = g$gene_id))
  assignments <- reads[, .(region_id = gene, RID, sample_name, read_start = start, read_end = end)]
  relative_reads <- build_read_table(assignments, anchors)
  nucs[, c("rel_start", "rel_end") := rel_interval(start + 1L, end, g$tss, g$strand)]
  nucs[, `:=`(region_id = g$gene_name, fp_class = "nucleosome", rel_mid = (rel_start + rel_end) / 2)]
  # Match the existing promoter positioning window for supporting metrics only.
  nucs <- nucs[rel_end >= -1000L & rel_start <= 1000L]
  plus <- nucs[rel_start >= 1L & rel_start <= 299L][order(rel_start), .(start = rel_start[1]), by = RID]
  starts <- plus$start[match(reads$RID, plus$RID)]
  if (!isTRUE(all.equal(as.numeric(starts), as.numeric(reads$plus1_start))))
    stop("Re-extracted +1 calls disagree with saved features")
  # Shared nearest-downstream-midpoint definition (can include a TSS-crossing
  # nucleosome) is distinct from the saved +1-start-in-1..299 definition above.
  positions <- read_position_summary(nucs, relative_reads, include_gap = FALSE)
  positions <- positions[match(reads$RID, RID)]
  if (!identical(positions$anchor_occupied, reads$nuc_at_tss))
    stop("Re-extracted TSS occupancy disagrees with saved features")
  reads[, `:=`(region_id = gene, rel_read_start = positions$rel_read_start,
                rel_read_end = positions$rel_read_end, plus1_mid = positions$plus1_mid)]
  parts <- list(summarize_subset(g, reads, nucs, anchors, seq_len(nrow(reads)), "gene"))
  label_path <- file.path(cluster_root, g$gene_name, "read_features.tsv.gz")
  labels_available <- file.exists(label_path) ||
    file.exists(file.path(cluster_root, g$group, paste0(g$gene_name, "_read_clusters.tsv")))
  if (by_cluster && nrow(reads)) {
    labels <- entropy_cluster_labels(g, reads, cluster_root)
    for (cl in sort(unique(labels$cluster))) {
      ix <- which(labels$cluster == cl)
      eligible <- unique(labels$eligible[ix])
      stopifnot(length(eligible) == 1L)
      parts[[length(parts) + 1L]] <- summarize_subset(g, reads, nucs, anchors, ix, "cluster", cl, eligible)
    }
  }
  summary <- rbindlist(lapply(parts, `[[`, "summary"))
  summary[, `:=`(cluster_labels_available = labels_available,
                 cohort = "existing_PolII_TSS_plus_minus_1000_full_span")]
  if (!nrow(reads)) summary[, error := g$source_error]
  list(summary = summary, profile = rbindlist(lapply(parts, `[[`, "profile")))
}

# PDF comparisons are generated separately by entropy_features.rmd.

#########
# Nucleosome entropy for the gene panel
#######
# 1. Parse options and select genes
# 2. Load gene TSSs, saved reads and sample sheet
# 3. With --footprints-only, recompute footprint fractions only
# 4. Otherwise compute entropy per gene in parallel and write tables
# 5. Compute and write short footprint fractions
main <- function() {
  genes <- fread(file.path(PROJECT, "code/Pol2/PolII/example_gene_list.tsv"))
  args <- commandArgs(trailingOnly = TRUE)
  flags <- args[startsWith(args, "--")]
  if (any(!flags %in% c("--by-cluster", "--gene-only", "--footprints-only"))) stop("Unknown option")
  if (all(c("--by-cluster", "--gene-only") %in% flags)) stop("Choose --by-cluster or --gene-only")
  by_cluster <- !"--gene-only" %in% flags # preserves the previous cluster outputs by default
  selected <- args[!startsWith(args, "--")]
  if (length(selected)) {
    if (any(!selected %in% genes$gene_name)) stop("Unknown gene(s): ", paste(setdiff(selected, genes$gene_name), collapse = ", "))
    genes <- genes[gene_name %in% selected]
  }
  stopifnot(!anyDuplicated(genes$gene_name))
  metadata <- fread(file.path(BASE, "read_features/tables/gene_features.tsv"))
  stopifnot(!anyDuplicated(metadata$gene_name), all(genes$gene_name %in% metadata$gene_name))
  metadata <- metadata[match(genes$gene_name, gene_name)]
  stopifnot(identical(genes$gene_id, metadata$gene_id), identical(genes$chrom, metadata$chrom))
  genes[, `:=`(tss = metadata$tss, strand = metadata$strand,
                tss_source = metadata$tss_source, source_error = metadata$error)]
  reads <- fread(file.path(BASE, "read_features/tables/read_features.tsv.gz"))
  samples <- fread(SAMPLES)
  out <- Sys.getenv("ENTROPY_OUT", file.path(PROJECT, "LCL_project/entropy/nucleosome_entropy"))
  workers <- as.integer(Sys.getenv("ENTROPY_WORKERS", Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
  stopifnot(!is.na(workers), workers >= 1L)
  table_dir <- file.path(out, "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  if ("--footprints-only" %in% flags) {
    saved <- fread(file.path(table_dir, "gene_entropy.tsv"))
    if (length(selected)) saved <- saved[gene_name %in% selected]
    if (!nrow(saved)) stop("No matching saved gene entropies for footprint extraction")
    write_short_footprint_fractions(saved, reads, table_dir, workers)
    return(invisible(NULL))
  }
  run_one <- function(i) {
    g <- genes[i]
    message("Entropy: ", g$gene_name)
    tryCatch(calculate_gene(g, copy(reads[gene == g$gene_name]), samples,
                           file.path(BASE, "recurrent_patterns/tables"), by_cluster = by_cluster),
      error = function(e) list(summary = cbind(g[, .(gene_name, gene_id, group, chrom, tss, strand, tss_source)],
        data.table(level = "gene", n_fibers = NA_integer_, total_nucleosome_occupancy_signal = NA_real_,
                   nucleosome_entropy_bits = NA_real_, status = "error", error = conditionMessage(e))),
        profile = NULL))
  }
  result <- parallel::mclapply(seq_len(nrow(genes)), run_one, mc.cores = workers)
  if (any(vapply(result, inherits, logical(1), "try-error"))) stop("Entropy worker failed")
  summary <- rbindlist(lapply(result, `[[`, "summary"), fill = TRUE)
  profiles <- rbindlist(lapply(result, `[[`, "profile"), fill = TRUE)
  stopifnot(nrow(summary[level == "gene"]) == nrow(genes))
  fwrite(summary, file.path(table_dir, "nucleosome_entropy.tsv"), sep = "\t", na = "NA")
  fwrite(summary[level == "gene"], file.path(table_dir, "gene_entropy.tsv"), sep = "\t", na = "NA")
  fwrite(profiles, file.path(table_dir, "occupancy_per_base.tsv.gz"), sep = "\t", na = "NA")
  failures <- summary[status == "error"]
  if (nrow(failures)) {
    print(failures[, .(gene_name, error)])
    stop("Input errors in ", nrow(failures), " genes; inspect ", file.path(table_dir, "gene_entropy.tsv"))
  }
  write_short_footprint_fractions(summary[level == "gene"], reads, table_dir, workers)
  message("Wrote ", nrow(summary[level == "gene"]), " genes, ", nrow(summary[level == "cluster"]),
          " clusters to ", out)
  print(summary[level == "gene", .N, by = status])
}

if (sys.nframe() == 0L) main()
