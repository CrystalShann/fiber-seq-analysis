# Per-molecule helpers used after clustering, without touching the cached
# features or clustering: pol2_knn_clusters() (cluster of a non-clustered
# molecule, used by the locus figures) and the `tracks` stage (nucleosome and
# TF footprints of the clustered molecules).
# Load pol2_clustering_functions.R and pol2_load_helpers() first.

# Cluster of each row of `features` (molecules x 40 bins): the majority cluster
# of its k nearest reference molecules by Manhattan distance, as in the Leiden
# graph; ties go to the nearest tied cluster.
###########################
# Assign nearby clusters
############################
# 1. Input: New features, reference features, reference labels and neighbour count.
# 2. Output: One cluster label for each new molecule.
#
# 3. Steps:
#    - Measure Manhattan distances
#    - Take the nearest reference molecules
#    - Use their majority label and resolve ties by proximity.
pol2_knn_clusters <- function(features, reference, reference_cluster, k) {
  stopifnot(ncol(features) == ncol(reference), length(reference_cluster) == nrow(reference))
  d <- matrix(0, nrow(features), nrow(reference))
  for (j in seq_len(ncol(reference))) d <- d + abs(outer(features[, j], reference[, j], "-"))
  apply(d, 1L, function(r) {
    nearest <- reference_cluster[order(r)[seq_len(k)]]
    votes <- table(nearest)
    nearest[nearest %in% names(votes)[votes == max(votes)]][1L]
  })
}

###########################
# Load saved clusters
############################
# 1. Input: Workflow settings and prepared inputs.
# 2. Output: The saved clustering result.
#
# 3. Steps:
#    - Read the clustering file
#    - Check the input signature and molecule order.
pol2_read_clustering <- function(cfg, prep = pol2_read_prepared(cfg)) {
  x <- readRDS(file.path(cfg$output, "tables", "clustering.rds"))
  if (!identical(x$signature, prep$signature)) stop("Clustering uses different inputs or features")
  stopifnot(identical(rownames(x$features), x$molecules$observation_id))
  x
}

# tables/molecule_footprints.tsv.gz: the footprints of each clustered molecule
# inside its TSS window, from its own sample: FiberHMM v2 recalled_tf blocks
# binned as TF footprints 10-30 (< 30), 30-60 (inclusive, the clustering rule)
# and 60-80 bp, and fibertools nucleosome calls of 130-160 bp (nuc_by_chr,
# longest alignment, sentinel blocks dropped). Other sizes are left out. Class
# names are MOLECULE_TRACK_COLORS names; rel_start/rel_end are strand-oriented
# TSS-relative bounds.
###########################
# Extract molecule footprints
############################
# 1. Input: Workflow settings and saved clusters.
# 2. Output: A saved table of strand-relative footprint tracks.
#
# 3. Steps:
#    - Read each molecule's own sample
#    - Keep overlapping TF and nucleosome calls in the selected size ranges
#    - Orient intervals around the TSS and save.
pol2_tracks <- function(cfg) {
  prep <- pol2_read_prepared(cfg)
  x <- pol2_read_clustering(cfg, prep)
  m <- merge(x$molecules[, .(observation_id, tss_id, sample_name, original_RID)],
             prep$tss[, .(tss_id, chrom, tss, strand, window_start, window_end)], by = "tss_id")
  samples <- prep$samples
  tf_paths <- pol2_tf_path(samples, cfg$footprint_dir)
  cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  parts <- parallel::mclapply(seq_len(nrow(samples)), function(i) {
    ms <- m[sample_name == samples$sample_name[i]]
    if (!nrow(ms)) return(NULL)
    regions <- GenomicRanges::reduce(pol2_gr(ms$chrom, ms$window_start, ms$window_end))
    tf <- pol2_read_footprints(tf_paths[i], regions)
    nuc <- rbindlist(lapply(unique(ms$chrom), function(chromosome) {
      path <- file.path(samples$fire_dir[i], "extracted_results/nuc_by_chr",
                        paste0(samples$sample_name[i], ".ft_extracted_nuc.", chromosome, ".bed.gz"))
      raw <- read_tabix_region(path, regions[as.character(GenomicRanges::seqnames(regions)) == chromosome])
      if (nrow(raw)) raw <- raw[raw[[4]] %in% ms$original_RID, , drop = FALSE]
      if (!nrow(raw)) return(NULL)
      b <- as.data.table(convert_ft_bed12_to_bed6(raw, format = "bed12_fibertools", longest_alignment = TRUE,
                                                  source = path))
      unique(b[, .(chrom = chr, start, end, original_RID = RID)])
    }))
    fp <- rbindlist(list(if (nrow(tf)) tf[, kind := "tf"], if (nrow(nuc)) nuc[, kind := "nuc"]), fill = TRUE)
    if (!nrow(fp)) return(NULL)
    fp <- merge(fp, ms[, .(original_RID, observation_id, mol_chrom = chrom, tss, strand, window_start, window_end)],
                by = "original_RID")
    fp <- fp[chrom == mol_chrom & end > window_start & start < window_end]
    fp[, length := end - start]
    fp[, class := fcase(kind == "nuc" & length >= 130L & length <= 160L, "nucleosome 130-160 bp",
                        kind == "tf" & length >= 10L & length < 30L, "TF footprint 10-30 bp",
                        kind == "tf" & length >= 30L & length <= 60L, "TF footprint 30-60 bp",
                        kind == "tf" & length > 60L & length <= 80L, "TF footprint 60-80 bp")]
    orient <- ifelse(fp$strand == "+", 1, -1)
    a <- orient * (fp$start + 0.5 - fp$tss)
    b <- orient * (fp$end + 0.5 - fp$tss)
    fp[, `:=`(rel_start = pmin(a, b), rel_end = pmax(a, b))]
    fp[!is.na(class), .(observation_id, chrom, start, end, length, class, rel_start, rel_end)]
  }, mc.cores = cores, mc.preschedule = FALSE)
  if (any(vapply(parts, inherits, logical(1), "try-error"))) stop("A sample's footprint extraction failed")
  out <- rbindlist(parts)
  dir.create(file.path(cfg$output,"tables","tracks"),recursive=TRUE,showWarnings=FALSE)
  fwrite(out, file.path(cfg$output, "tables", "tracks", "molecule_footprints.tsv.gz"), sep = "\t")
  message("Saved ", nrow(out), " footprints of ", uniqueN(out$observation_id), " clustered molecules")
  invisible(out)
}

# Project in bounded batches to avoid allocating chromosome-wide distance matrices.
###########################
# Assign remaining molecules
############################
# 1. Input: Workflow settings and saved feature caches.
# 2. Output: A saved table of molecule cluster assignments.
#
# 3. Steps:
#    - Keep existing cluster labels
#    - Assign remaining molecules from nearby references in small batches
#    - Check identities and save.
pol2_project <- function(cfg) {
  prep <- pol2_read_prepared(cfg)
  x <- pol2_read_clustering(cfg, prep)
  cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  parts <- parallel::mclapply(cfg$chromosomes, function(cc) {
    f <- readRDS(file.path(cfg$output, "tables", "features", paste0(cc, ".rds")))
    stopifnot(identical(f$signature$prepared, prep$signature))
    m <- copy(f$molecules)
    message("Projecting ", cc, ": ", nrow(m), " observations")
    hit <- match(m$observation_id, x$molecules$observation_id)
    m[, `:=`(clustered = !is.na(hit), cluster = as.character(x$molecules$cluster[hit]))]
    todo <- which(is.na(hit))
    for (batch in split(todo, ceiling(seq_along(todo) / 200L))) {
      m$cluster[batch] <- pol2_knn_clusters(f$features[m$observation_id[batch], , drop = FALSE],
        x$features, x$molecules$cluster, cfg$k_neighbors)
    }
    message("Finished ", cc)
    m[, .(observation_id, tss_id, sample_name, RID, clustered, cluster)]
  }, mc.cores = cores, mc.preschedule = FALSE)
  if (any(vapply(parts, inherits, logical(1), "try-error"))) stop("Molecule projection failed")
  out <- rbindlist(parts)
  stopifnot(!anyNA(out$cluster), !anyDuplicated(out$observation_id))
  dir.create(file.path(cfg$output,"tables","tracks"),recursive=TRUE,showWarnings=FALSE)
  fwrite(out, file.path(cfg$output, "tables", "tracks", "molecule_projection.tsv.gz"), sep = "\t")
  message("Projected ", nrow(out), " observations of ", uniqueN(out$RID), " physical molecules")
  invisible(out)
}

###########################
# Select example promoters
############################
# 1. Input: Workflow settings.
# 2. Output: One promoter per cluster with its read counts.
#
# 3. Steps:
#    - Count all eligible reads at each TSS, using the saved cluster projections
#    - Keep protein-coding TSSs with at least 100 reads present in every sample
#    - Per cluster, choose the TSS with the highest fraction of its reads in that cluster.
pol2_select_example_promoters <- function(cfg) {
  prep <- readRDS(file.path(cfg$output, "tables", "prepared.rds"))
  p <- fread(file.path(cfg$output, "tables", "tracks", "molecule_projection.tsv.gz"))
  totals <- p[, .(n_total = .N, n_samples = uniqueN(sample_name)), by = tss_id]
  counts <- p[, .(n_cluster = .N), by = .(tss_id, cluster)][totals, on = "tss_id"]
  counts[, fraction := n_cluster / n_total]
  counts <- merge(counts, prep$tss[, .(tss_id, gene_name, gene_type)], by = "tss_id")
  counts <- counts[gene_type == "protein_coding" & n_total >= 100L & n_samples == nrow(prep$samples)]
  setorder(counts, -fraction, -n_total, tss_id)
  selected <- counts[, .SD[1L], by = cluster]
  selected[order(as.integer(sub("^cluster", "", cluster)))]
}
