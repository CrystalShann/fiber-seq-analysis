# All full-window molecules define observed sites at each TSS; only molecules with a
# 30-60 bp FiberHMM footprint in the window define clustering rows.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Matrix)
  library(data.table)
})

###########################
# Load shared tools
############################
# 1. Input: Project folder.
# 2. Output: Shared tools available for the workflow.
#
# 3. Steps:
#    - Find the shared parsing, statistics and clustering files
#    - Load them into the working session.
pol2_load_helpers <- function(project) {
  source(file.path(project, "code/parsing_functions/parsing_footprints_functions.r"), local = .GlobalEnv)
  source(file.path(project, "code/parsing_functions/stats_testing_functions.r"), local = .GlobalEnv)
  source(file.path(project, "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r"), local = .GlobalEnv)
  source(file.path(project, "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r"), local = .GlobalEnv)
}

###########################
# Convert genomic coordinates
############################
# 1. Input: Chromosomes and zero-based interval boundaries.
# 2. Output: Genomic ranges with one-based starts.
#
# 3. Steps:
#    - Shift each start by one base
#    - Keep the chromosome and end position.
pol2_gr <- function(chrom, start0, end0) GRanges(chrom, IRanges(start0 + 1L, end0))
###########################
# Check required files
############################
# 1. Input: Paths to required files.
# 2. Output: An error if any file is absent.
#
# 3. Steps:
#    - Check each path
#    - List missing files and stop if needed.
pol2_required <- function(paths) {
  absent <- paths[!file.exists(paths)]
  if (length(absent)) stop("Missing input: ", paste(absent, collapse = ", "))
}
###########################
# Record file details
############################
# 1. Input: Paths to input files.
# 2. Output: A table of full paths, sizes and modification times.
#
# 3. Steps:
#    - Check that the files exist
#    - Collect their details.
pol2_stamp <- function(paths) {
  pol2_required(paths)
  x <- file.info(paths)
  data.frame(path = normalizePath(paths), size = x$size, mtime = as.numeric(x$mtime))
}
###########################
# Identify workflow inputs
############################
# 1. Input: Workflow settings.
# 2. Output: A record of settings, input files and source-file checksums.
#
# 3. Steps:
#    - Select settings that affect the features
#    - Record input details and checksums.
pol2_signature <- function(cfg) {
  keys <- c("samples", "expected_samples", "tss", "chrom_sizes", "footprint_dir", "peaks", "signal", "chromosomes", "bin_bp")
  list(version = 2L, config = cfg[keys], inputs = pol2_stamp(unlist(cfg[c("samples", "tss", "chrom_sizes", "peaks", "signal")])),
       source = tools::md5sum(file.path(cfg$project, c("code/pausing/clustering_Pol2/pol2_clustering_functions.R",
         "code/parsing_functions/parsing_footprints_functions.r", "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r"))))
}
###########################
# Save a result safely
############################
# 1. Input: An object and its destination.
# 2. Output: A saved R data file.
#
# 3. Steps:
#    - Create the destination folder
#    - Write a temporary file
#    - Rename it to the final destination.
pol2_save <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Cannot save ", path)
}

###########################
# Set promoter windows
############################
# 1. Input: TSS locations and strands.
# 2. Output: The TSS table with promoter boundaries.
#
# 3. Steps:
#    - Place the window from 1,000 bases upstream to 999 bases downstream
#    - Orient the boundaries to each strand.
pol2_windows <- function(tss) {
  tss <- copy(tss)
  tss[, `:=`(window_start = ifelse(strand == "+", tss - 1001L, tss - 1000L),
             window_end = ifelse(strand == "+", tss + 999L, tss + 1000L))]
  tss
}

###########################
# Set Pol II signal windows
############################
# 1. Input: TSS locations and strands.
# 2. Output: Genomic intervals for measuring promoter signal.
#
# 3. Steps:
#    - Use 50 bases upstream through 299 bases downstream
#    - Convert the window to genomic coordinates.
pol2_tssr <- function(tss) {
  data.table(chrom = tss$chrom,
             start = ifelse(tss$strand == "+", tss$tss - 51L, tss$tss - 300L),
             end = ifelse(tss$strand == "+", tss$tss + 299L, tss$tss + 50L))
}


###########################
# Find qualifying footprints
############################
# 1. Input: Promoter windows and footprint calls.
# 2. Output: Matching footprints linked to their TSSs.
#
# 3. Steps:
#    - Keep footprints 30–60 bases long
#    - Find overlaps with promoter windows
#    - Record strand-relative midpoints.
pol2_candidates <- function(tss, fp) {
  fp <- unique(copy(fp))
  if (!nrow(fp) || !nrow(tss)) return(data.table())
  fp <- fp[end - start >= 30L & end - start <= 60L]
  if (!nrow(fp)) return(data.table())
  fp[, footprint_class := "30-60 bp"]
  fg <- pol2_gr(fp$chrom, fp$start, fp$end)
  hits <- findOverlaps(fg, pol2_gr(tss$chrom, tss$window_start, tss$window_end), ignore.strand = TRUE)
  if (!length(hits)) return(data.table())
  out <- copy(fp[queryHits(hits)])
  ti <- subjectHits(hits)
  out[, tss_id := tss$tss_id[ti]]
  out[, midpoint_relative := ((start + 1 + end) / 2 - tss$tss[ti]) * ifelse(tss$strand[ti] == "+", 1, -1)]
  out
}

###########################
# Measure signal over intervals
############################
# 1. Input: Genomic intervals and a bigWig signal file.
# 2. Output: Intervals with mean signal and covered fractions.
#
# 3. Steps:
#    - Read overlapping signal
#    - Weight each value by its overlap length
#    - Divide by covered bases.
pol2_signal <- function(fp, path) {
  fp <- copy(fp)
  if (!nrow(fp)) return(fp)
  gr <- pol2_gr(fp$chrom, fp$start, fp$end)
  signal <- rtracklayer::import(path, which = GenomicRanges::reduce(gr))
  h <- findOverlaps(gr, signal, ignore.strand = TRUE)
  weighted <- numeric(nrow(fp)); covered <- numeric(nrow(fp))
  if (length(h)) {
    q <- queryHits(h); s <- subjectHits(h)
    w <- pmin(end(gr)[q], end(signal)[s]) - pmax(start(gr)[q], start(signal)[s]) + 1
    score <- as.numeric(signal$score[s]); good <- is.finite(score)
    z <- data.table(i = q[good], w = w[good], signal = w[good] * score[good])[, .(w = sum(w), signal = sum(signal)), by = i]
    covered[z$i] <- z$w; weighted[z$i] <- z$signal
  }
  fp[, `:=`(pol2_mean = ifelse(covered > 0, weighted / covered, NA_real_),
             pol2_covered_fraction = covered / (end - start))]
  fp
}

###########################
# Prepare samples and TSSs
############################
# 1. Input: Workflow settings and configured input files.
# 2. Output: Saved sample and promoter information.
#
# 3. Steps:
#    - Check samples and TSS coordinates
#    - Keep promoters within chromosome bounds and Pol II peaks
#    - Add promoter signal and save.
pol2_prepare <- function(cfg) {
  sig <- pol2_signature(cfg)
  samples <- fread(cfg$samples)
  stopifnot(nrow(samples) == cfg$expected_samples, !anyDuplicated(samples$sample_name),
            all(c("sample_name", "fire_dir", "cell_line") %in% names(samples)), cfg$bin_bp == 50)
  tf_paths <- pol2_tf_path(samples, cfg$footprint_dir)
  pol2_required(c(tf_paths, paste0(tf_paths, ".tbi")))
  tss <- fread(cfg$tss)
  setnames(tss, names(tss)[1L], "chrom")
  stopifnot(all(c("start", "end", "strand", "gene_id", "gene_name") %in% names(tss)),
            all(tss$end == tss$start + 1), all(tss$strand %in% c("+", "-")))
  tss <- tss[chrom %in% cfg$chromosomes]
  tss[, `:=`(tss = end, tss_id = paste(gene_id, chrom, end, strand, sep = ":"))]
  stopifnot(!anyDuplicated(tss$tss_id))
  tss <- pol2_windows(tss)
  sizes <- fread(cfg$chrom_sizes, header = FALSE, select = 1:2, col.names = c("chrom", "length"))
  tss <- merge(tss, sizes, by = "chrom")
  tss <- tss[window_start >= 0 & window_end <= length]
  peaks <- fread(cfg$peaks, header = FALSE, select = 1:3, col.names = c("chrom", "start", "end"))
  tss <- tss[countOverlaps(pol2_gr(chrom, start, end), pol2_gr(peaks$chrom, peaks$start, peaks$end), ignore.strand = TRUE) > 0]
  signal <- pol2_signal(pol2_tssr(tss), cfg$signal)
  tss[, pol2_tss_signal := signal$pol2_mean]
  setorder(tss, chrom, tss, strand, gene_id)
  prepared <- list(signature = sig, samples = samples, tss = tss)
  pol2_save(prepared, file.path(cfg$output, "tables", "prepared.rds"))
  message("Prepared ", nrow(tss), " TSSs in a Pol II peak; all ", nrow(samples), " samples will be pooled")
  invisible(prepared)
}

###########################
# Load prepared inputs
############################
# 1. Input: Workflow settings.
# 2. Output: Previously prepared samples and promoters.
#
# 3. Steps:
#    - Read the saved inputs
#    - Check that their signature still matches.
pol2_read_prepared <- function(cfg) {
  path <- file.path(cfg$output, "tables", "prepared.rds"); pol2_required(path)
  x <- readRDS(path)
  if (!identical(x$signature, pol2_signature(cfg))) stop("Inputs or feature code changed; run prepare again")
  x
}

###########################
# Locate m6A files
############################
# 1. Input: Sample information and a chromosome.
# 2. Output: One m6A file path per sample.
#
# 3. Steps:
#    - Use each sample's extraction folder
#    - Add the sample and chromosome names.
pol2_ft_path <- function(samples, chrom) file.path(samples$fire_dir, "extracted_results/m6a_by_chr",
  paste0(samples$sample_name, ".ft_extracted_m6a.", chrom, ".bed.gz"))
###########################
# Locate TF footprint files
############################
# 1. Input: Sample information and the footprint folder.
# 2. Output: One TF footprint file path per sample.
#
# 3. Steps:
#    - Find each sample's folder
#    - Add its recalled-footprint filename.
pol2_tf_path <- function(samples, dir) file.path(dir, samples$sample_name,
  paste0(samples$sample_name, ".recalled_tf.bed.gz"))


###########################
# Read footprint intervals
############################
# 1. Input: An indexed footprint file and genomic regions.
# 2. Output: Unique footprint intervals with read identities.
#
# 3. Steps:
#    - Read the requested regions
#    - Expand footprint blocks
#    - Remove duplicate intervals.
pol2_read_footprints <- function(path, regions) {
  raw <- read_tabix_region(path, regions)
  if (!nrow(raw)) return(data.table())
  fp <- as.data.table(convert_ft_bed12_to_bed6(raw, format = "bed15_fiberhmm_tf", source = path))
  unique(fp[, .(chrom = chr, start, end, original_RID = RID)])
}


###########################
# Select eligible molecules
############################
# 1. Input: Reads, qualifying footprints and one promoter.
# 2. Output: Reads with a qualifying footprint that span the whole window.
#
# 3. Steps:
#    - Match sample and read identities
#    - Require coverage of both window boundaries.
pol2_eligible_reads <- function(bed, candidate, site) {
  keys <- paste(candidate$sample_name, candidate$original_RID, sep = "::")
  bed[paste(sample_name, original_RID, sep = "::") %in% keys &
        start <= site$window_start & end >= site$window_end]
}


###########################
# Build one promoter's features
############################
# 1. Input: Reads, qualifying footprints and one promoter.
# 2. Output: Accessibility features, base-level calls, molecule details and footprints.
#
# 3. Steps:
#    - Find observed m6A positions using full-window reads
#    - Select eligible molecules
#    - Build strand-oriented calls and 50-base bins
#    - Attach identities and denominators.
pol2_locus_features <- function(bed, candidate, site) {
  pool <- copy(bed[start <= site$window_start & end >= site$window_end])
  if (!nrow(pool)) return(NULL)
  pool[, RID := paste(sample_name, original_RID, sep = "::")]
  calls <- as.data.table(convert_ft_bed12_to_bed6(as.data.frame(pool[, 1:12]), format = "bed12_fibertools"))
  if (!nrow(calls)) return(NULL)
  stopifnot(all(calls$end - calls$start == 1L))
  calls <- calls[end > site$window_start & end <= site$window_end,
                 .(RID, relative = (end - site$tss) * ifelse(site$strand == "+", 1L, -1L))]
  if (!nrow(calls)) return(NULL)
  positions <- sort(unique(calls$relative))
  rows <- pol2_eligible_reads(pool, candidate, site)
  if (!nrow(rows)) return(NULL)
  pairs <- unique(calls[RID %in% rows$RID, .(i = match(RID, rows$RID), j = match(relative, positions))])
  stopifnot(!anyNA(pairs))
  mat <- sparseMatrix(i = pairs$i, j = pairs$j, x = 1, dims = c(nrow(rows), length(positions)),
                      dimnames = list(rows$RID, as.character(positions)))
  binned <- bin_read_matrix(mat, -1000L, 999L, window_size = 50L)
  if (any(binned$window_anno$n_sites == 0L)) return(NULL)
  stopifnot(ncol(binned$mat) == 40L, !anyNA(binned$mat))
  meta <- rows[, .(sample_name, original_RID, RID, read_start = start + 1L, read_end = end, cell_line)]
  meta[, `:=`(tss_id = site$tss_id, gene_id = site$gene_id, gene_name = site$gene_name,
              chrom = site$chrom, tss = site$tss, strand = site$strand, pol2_tss_signal = site$pol2_tss_signal,
              m6a_calls = as.numeric(Matrix::rowSums(mat)), n_positions = length(positions),
              m6a_proportion = as.numeric(Matrix::rowMeans(mat)))]
  meta[, observation_id := paste(tss_id, RID, sep = "::")]
  rownames(binned$mat) <- meta$observation_id
  fp <- merge(candidate, meta[, .(sample_name, original_RID, observation_id)],
              by = c("sample_name", "original_RID"))
  fp[, footprint_length := end - start]
  denominators <- as.data.table(binned$window_anno)
  denominators[, tss_id := site$tss_id]
  list(features = binned$mat, sites = data.table(i = pairs$i, j = positions[pairs$j] + 1001L),
       molecules = meta, footprints = fp, denominators = denominators)
}

###########################
# Save chromosome features
############################
# 1. Input: Workflow settings and a chromosome.
# 2. Output: A saved chromosome feature cache.
#
# 3. Steps:
#    - Reuse a matching cache when available
#    - Read promoters in batches
#    - Build each promoter's features
#    - Combine and save the results.
pol2_features <- function(cfg, chrom) {
  prep <- pol2_read_prepared(cfg)
  chromosome <- chrom
  sites <- prep$tss[chrom == chromosome]
  paths <- pol2_ft_path(prep$samples, chrom)
  tf_paths <- pol2_tf_path(prep$samples, cfg$footprint_dir)
  stamp <- if (nrow(sites)) pol2_stamp(c(paths, tf_paths, paste0(c(paths, tf_paths), ".tbi"))) else NULL
  signature <- list(prepared = prep$signature, ft = stamp)
  path <- file.path(cfg$output, "tables", "features", paste0(chrom, ".rds"))
  if (file.exists(path) && identical(readRDS(path)$signature, signature)) {
    message("Reusing ", chrom, " features"); return(invisible(NULL))
  }
  chunks <- list()
  eligible_RIDs <- list()
  for (idx in split(seq_len(nrow(sites)), ceiling(seq_len(nrow(sites)) / 100L))) {
    block <- sites[idx]
    regions <- GenomicRanges::reduce(pol2_gr(block$chrom, block$window_start, block$window_end))
    reads <- lapply(seq_len(nrow(prep$samples)), function(i) {
      b <- as.data.table(read_ft_bed12(paths[i], regions, longest_alignment = TRUE))
      if (nrow(b)) b[, `:=`(original_RID = RID, sample_name = prep$samples$sample_name[i], cell_line = prep$samples$cell_line[i])]
      fp <- pol2_read_footprints(tf_paths[i], regions)
      if (nrow(fp)) fp[, sample_name := prep$samples$sample_name[i]]
      list(bed = b, footprints = fp)
    })
    bed <- rbindlist(lapply(reads, `[[`, "bed"))
    candidates <- pol2_candidates(block, rbindlist(lapply(reads, `[[`, "footprints")))
    rm(reads)
    if (!nrow(bed) || !nrow(candidates)) next
    for (i in seq_len(nrow(block))) {
      site <- block[i]
      candidate <- candidates[tss_id == site$tss_id]
      eligible <- pol2_eligible_reads(bed, candidate, site)
      if (!nrow(eligible)) next
      eligible_RIDs[[length(eligible_RIDs) + 1L]] <- paste(eligible$sample_name, eligible$original_RID, sep = "::")
      chunk <- pol2_locus_features(bed, candidate, site)
      if (!is.null(chunk)) chunks[[length(chunks) + 1L]] <- chunk
    }
  }
  features <- if (length(chunks)) do.call(rbind, lapply(chunks, `[[`, "features")) else NULL
  sites_mat <- NULL
  if (length(chunks)) {
    n_rows <- vapply(chunks, function(x) nrow(x$features), integer(1))
    offsets <- cumsum(c(0L, head(n_rows, -1L)))
    ij <- rbindlist(Map(function(x, o) x$sites[, .(i = i + o, j)], chunks, offsets))
    sites_mat <- sparseMatrix(i = ij$i, j = ij$j, dims = c(sum(n_rows), 2000L),
                              dimnames = list(rownames(features), as.character(-1000:999)))
  }
  out <- list(signature = signature, eligible_RIDs = unique(unlist(eligible_RIDs)),
              features = features, sites = sites_mat,
              molecules = rbindlist(lapply(chunks, `[[`, "molecules")),
              footprints = rbindlist(lapply(chunks, `[[`, "footprints")),
              denominators = rbindlist(lapply(chunks, `[[`, "denominators")))
  pol2_save(out, path)
  message(chrom, ": ", length(out$eligible_RIDs), " footprint-positive full-span molecules; saved ",
          nrow(out$molecules), " molecule-TSS observations with complete bin features")
  invisible(out)
}

###########################
# Estimate clustering memory
############################
# 1. Input: Number of molecules and feature columns.
# 2. Output: Estimated memory use in GiB.
#
# 3. Steps:
#    - Allow for pairwise distances and feature storage
#    - Add working space.
pol2_memory_gb <- function(n, p = 40) (32 * as.double(n)^2 + 8 * as.double(n) * p) / 1024^3 + 2


###########################
# Sample distinct molecules
############################
# 1. Input: Molecule metadata, sample size and random seed.
# 2. Output: Selected row numbers.
#
# 3. Steps:
#    - Sample physical molecules without replacement
#    - Choose one TSS observation per selected molecule.
pol2_sample_rows <- function(meta, n = 5000L, seed = 1L) {
  stopifnot(length(n) == 1L, is.finite(n), n >= 3L, n == as.integer(n), !anyNA(meta$RID))
  set.seed(seed)
  molecules <- unique(meta$RID)
  selected <- molecules[sample.int(length(molecules), min(n, length(molecules)), replace = FALSE)]
  rows <- split(seq_len(nrow(meta)), meta$RID)
  vapply(selected, function(id) {
    possible <- rows[[id]]
    possible[sample.int(length(possible), 1L)]
  }, integer(1), USE.NAMES = FALSE)
}

###########################
# Cluster pooled molecules
############################
# 1. Input: Workflow settings and saved chromosome features.
# 2. Output: Saved clusters, molecule identities, calls and profiles.
#
# 3. Steps:
#    - Check the feature caches
#    - Sample distinct molecules and check memory
#    - Build the neighbour graph and run Leiden clustering
#    - Add footprint signal and save.
pol2_cluster <- function(cfg) {
  prep <- pol2_read_prepared(cfg)
  files <- file.path(cfg$output, "tables", "features", paste0(cfg$chromosomes, ".rds"))
  pol2_required(files)
  chunks <- lapply(files, readRDS)
  footprint_eligible_n <- length(unique(unlist(lapply(chunks, `[[`, "eligible_RIDs"))))
  message("Across all samples: ", footprint_eligible_n, " unique footprint-positive full-span molecules")
  for (i in seq_along(chunks)) {
    x <- chunks[[i]]
    if (!identical(x$signature$prepared, prep$signature)) stop("Stale chromosome features: ", files[i])
    if (!is.null(x$signature$ft) && !identical(x$signature$ft, pol2_stamp(x$signature$ft$path)))
      stop("m6A or footprint inputs changed: ", files[i])
  }
  mat <- do.call(rbind, lapply(chunks, `[[`, "features"))
  if (is.null(mat) || nrow(mat) < 3L) stop("Fewer than three eligible observations")
  sites <- do.call(rbind, Filter(Negate(is.null), lapply(chunks, `[[`, "sites")))
  meta <- rbindlist(lapply(chunks, `[[`, "molecules"))
  stopifnot(identical(rownames(mat), meta$observation_id), identical(rownames(sites), rownames(mat)),
            !anyDuplicated(meta$observation_id), !anyNA(mat), all(mat >= 0 & mat <= 1))
  eligible_n <- uniqueN(meta$RID)
  if (eligible_n < 3L) stop("Fewer than three eligible physical molecules")
  selected <- pol2_sample_rows(meta, n = cfg$sample_n, seed = cfg$seed)
  mat <- mat[selected, , drop = FALSE]
  sites <- sites[selected, , drop = FALSE]
  meta <- meta[selected]
  stopifnot(!anyDuplicated(meta$RID), identical(rownames(mat), meta$observation_id))
  estimate <- pol2_memory_gb(nrow(mat), ncol(mat))
  limit <- cfg$max_memory_gb
  slurm_mb <- suppressWarnings(as.numeric(Sys.getenv("SLURM_MEM_PER_NODE", "")))
  if (is.finite(slurm_mb)) limit <- min(limit, 0.8 * slurm_mb / 1024)
  if (estimate > limit) stop(sprintf(
    "Exact global clustering of %s sampled molecules needs an estimated %.1f GiB (budget %.1f). Features are saved; no approximate graph was substituted.", nrow(mat), estimate, limit))
  message("One pooled clustering: ", nrow(mat), " molecules randomly selected from ", eligible_n,
          "; estimated working memory ", round(estimate, 1), " GiB")
  knn <- manhattan_knn_graph(mat, k_neighbors = cfg$k_neighbors, sigma = NULL)
  part <- leiden_partition(knn$graph, resolution = cfg$resolution, seed = cfg$seed)
  membership <- part$membership[rownames(mat)]
  meta[, cluster := as.character(relabel_clusters_by_size(membership))]
  profiles <- cluster_feature_profiles(mat, relabel_clusters_by_size(membership))
  footprints <- rbindlist(lapply(chunks, `[[`, "footprints"))[observation_id %in% meta$observation_id]
  # Bulk Pol II signal is only measured over the footprints of sampled molecules.
  signal <- pol2_signal(unique(footprints[, .(chrom, start, end)]), cfg$signal)
  footprints <- merge(footprints, signal, by = c("chrom", "start", "end"), sort = FALSE)
  result <- list(signature = prep$signature, feature_stamps = pol2_stamp(files),
    config = cfg, eligible_molecules = eligible_n, footprint_eligible_molecules = footprint_eligible_n,
    molecules = meta, features = mat, sites = sites, profiles = profiles,
    graph = knn$graph, sigma = knn$sigma, quality = part$quality,
    footprints = footprints,
    denominators = rbindlist(lapply(chunks, `[[`, "denominators"))[tss_id %in% meta$tss_id])
  pol2_save(result, file.path(cfg$output, "tables", "clustering.rds"))
  message("Saved ", length(unique(meta$cluster)), " shared clusters")
  invisible(result)
}


###########################
# Attach pausing indices
############################
# 1. Input: Molecule or TSS metadata and an optional pausing-index table.
# 2. Output: The metadata with pausing indices and match status.
#
# 3. Steps:
#    - Remove gene-version suffixes
#    - Match gene identities
#    - Mark matched and unmatched entries.
pol2_attach_pi <- function(meta, path) {
  out <- copy(meta)
  out[, `:=`(pi = NA_real_, pi_status = "not_supplied")]
  if (is.null(path)) return(out)
  pol2_required(path)
  pi_data <- fread(path)
  if (!all(c("ensg", "pi") %in% names(pi_data))) stop("PI table needs: ensg, pi")
  pi_data <- pi_data[is.finite(pi) & grepl("^ENSG", ensg)]
  pi_data[, gene_key := sub("\\..*$", "", ensg)]
  if (anyDuplicated(pi_data$gene_key)) stop("Duplicate Ensembl IDs with a PI")
  hit <- match(sub("\\..*$", "", out$gene_id), pi_data$gene_key)
  out[, pi := pi_data$pi[hit]]
  out[, pi_status := fifelse(is.na(hit), "no_pi", "matched")]
  out
}


###########################
# Compare promoter signal
############################
# 1. Input: Molecule metadata with clusters and promoter signal.
# 2. Output: Cluster signal averages, relative enrichment and adjusted significance.
#
# 3. Steps:
#    - Compare each cluster mean with the overall mean
#    - Test each cluster against the remaining molecules
#    - Adjust for multiple comparisons.
pol2_signal_enrichment <- function(m) {
  m <- m[is.finite(pol2_tss_signal)]
  overall <- mean(m$pol2_tss_signal)
  out <- m[, .(n_with_signal = .N, mean_tss_signal = mean(pol2_tss_signal)), by = cluster]
  out[, log2_ratio := log2(mean_tss_signal / overall)]
  out[, p := vapply(cluster, function(cl) stats::wilcox.test(
    m[cluster == cl, pol2_tss_signal], m[cluster != cl, pol2_tss_signal])$p.value, numeric(1))]
  out[, fdr := adjust_bh(p)]
  out[, signif := fcase(fdr < 0.001, "***", fdr < 0.01, "**", fdr < 0.05, "*", default = "NS")]
  out
}


###########################
# Summarize clusters
############################
# 1. Input: Molecule metadata with pausing indices and signal.
# 2. Output: Cluster counts, enrichment and display order.
#
# 3. Steps:
#    - Count molecules and summarize signal
#    - Calculate median pausing index across distinct TSSs
#    - Order clusters by pausing index.
pol2_cluster_summary <- function(m) {
  pi <- unique(m[is.finite(pi), .(cluster, tss_id, pi)])[, .(median_pi = median(pi), n_tss_with_pi = .N), by = cluster]
  out <- merge(m[, .(n = .N), by = cluster], pol2_signal_enrichment(m), by = "cluster", all.x = TRUE)
  out <- merge(out, pi, by = "cluster", all.x = TRUE)
  out[is.na(n_tss_with_pi), n_tss_with_pi := 0L]
  out[, size_rank := as.integer(sub("^cluster", "", cluster))]
  setorder(out, median_pi, size_rank, na.last = TRUE)
  out[, size_rank := NULL]
  out[, display_order := .I]
  out
}

###########################
# Build the clustering report
############################
# 1. Input: Workflow settings and the script folder.
# 2. Output: Saved report tables, figures and an HTML report.
#
# 3. Steps:
#    - Check saved clustering inputs
#    - Add annotations and summaries
#    - Save figures and run details
#    - Render the report.
pol2_report <- function(cfg, here) {
  tables <- file.path(cfg$output, "tables")
  pol2_required(file.path(tables, "clustering.rds"))
  x <- readRDS(file.path(tables, "clustering.rds"))
  if (!identical(x$signature, pol2_signature(cfg))) stop("Clustering uses different inputs or features")
  cluster_keys <- c("sample_n", "seed", "k_neighbors", "resolution")
  if (!identical(cfg[cluster_keys], x$config[cluster_keys])) stop("Clustering parameters changed; rerun cluster")
  m <- pol2_attach_pi(x$molecules, cfg$pi_table)
  summary <- pol2_cluster_summary(m)
  m[, cluster := factor(cluster, levels = summary$cluster)]
  fp <- merge(x$footprints, m[, .(observation_id, cluster)], by = "observation_id", sort = FALSE)
  fwrite(m, file.path(tables, "molecule_assignments.tsv.gz"), sep = "\t")
  fwrite(fp, file.path(tables, "candidate_footprints.tsv.gz"), sep = "\t")
  fwrite(x$denominators, file.path(tables, "bin_denominators.tsv.gz"), sep = "\t")
  fwrite(summary, file.path(tables, "cluster_summary.tsv"), sep = "\t")
  fractions <- m[, .(n = .N), by = .(sample_name, cell_line, cluster)]
  fractions[, fraction := n / sum(n), by = sample_name]
  fwrite(fractions, file.path(tables, "sample_cluster_fractions.tsv"), sep = "\t")
  tss_counts <- m[, .(n = .N, pi = pi[1L], pi_status = pi_status[1L], pol2_tss_signal = pol2_tss_signal[1L]),
                  by = .(tss_id, gene_id, gene_name, cluster)]
  tss_counts[, fraction := n / sum(n), by = tss_id]
  fwrite(tss_counts, file.path(tables, "tss_cluster_fractions.tsv.gz"), sep = "\t")
  means <- data.table(cluster = rep(rownames(x$profiles), each = ncol(x$profiles)),
                      relative_bp = rep(as.numeric(colnames(x$profiles)), times = nrow(x$profiles)),
                      m6a_proportion = as.vector(t(x$profiles)))
  fwrite(means, file.path(tables, "cluster_profiles.tsv"), sep = "\t")
  annotation <- list(molecules = m, footprints = fp, profiles = means, summary = summary)
  pol2_save(annotation, file.path(tables, "annotations.rds"))
  pol2_save_figures(x, annotation, cfg)
  yaml::write_yaml(cfg, file.path(cfg$output, "run_config.yaml"))
  capture.output(sessionInfo(), file = file.path(cfg$output, "session_info.txt"))
  rmarkdown::render(file.path(here, "pol2_clustering_report.Rmd"), output_file = "pol2_clustering_report.html", output_dir = cfg$output,
                    params = list(output = cfg$output), envir = new.env(parent = globalenv()), quiet = TRUE)
  message("Saved pooled report to ", file.path(cfg$output, "pol2_clustering_report.html"))
}
