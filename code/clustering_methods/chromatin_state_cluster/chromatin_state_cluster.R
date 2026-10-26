#!/usr/bin/env Rscript
# Shared, TSS-aligned chromatin states across all genes and manifest samples.
# Three independently invoked stages: prepare (one chromosome), cluster (one or
# all broad states), report (join assignments to original promoter tables).
# Existing source files are loaded, never changed; no helper functions defined.
#
# Authoritative rows: Step 3 gene_id/sample_name/rid/read_start/read_end and
# promoter_state. read_start/end are BED coordinates of a primary alignment.
# Never reselect the longest alignment or downsample the discovery population.
# Features: 25-bp bins over inclusive -1000..+2000, oriented by GENE strand;
# four occupancy channels: >=90-bp nucleosomes, shorter nucleosomes, FIRE,
# linker. Distance weights = bin width / (2 * full window width), per channel.
# No FASTA, m6A, CpG, expression, or PI enters the clustering distance.
# Profiles with uncovered bases, conflicting classes, or a broad-state mismatch
# are retained as QC failures with no subcluster, never imputed to zero.
#
# Within each broad state: approximate Manhattan KNN -> exponential affinities
# -> existing Leiden implementation -> merge/split to EXACTLY three groups.
# Merges use closest weighted mean profiles; splits use deterministic distant
# seeds in the group with greatest total L1 dispersion. These final groups are
# adjusted Leiden partitions, not unconstrained Leiden communities. Fewer than
# three distinct usable profiles is an explicit insufficient-variation result.
# RcppAnnoy validation uses a bounded, uniformly sampled diagnostic subset only;
# the actual discovery graph always includes every usable gene/fiber row.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(Rsamtools)
  library(igraph)
})
repo <- "/project/spott/cshan/fiber-seq"
opts <- list(stage = "", chrom = "", state = "all",
  chroms = paste(paste0("chr", c(1:22, "X", "Y")), collapse = ","),
  `feature-dir` = file.path(repo, "LCL_project/pausing/dominant_tss"),
  genes = file.path(repo, "LCL_project/pausing/dominant_tss/pause_index_dominant.tsv"),
  samples = "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv",
  `out-dir` = file.path(repo, "LCL_project/pausing/dominant_tss/chromatin_state_cluster"),
  `window-start` = "-1000", `window-end` = "2000", `bin-size` = "25", `nuc-min` = "90",
  k = "10", trees = "50", `search-k` = "5000", resolution = "1", seed = "1",
  `validate-size` = "2000", `plot-per-cluster` = "100")
argv <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% argv) {
  cat("Usage: Rscript chromatin_state_cluster.R --stage prepare|cluster|report [--option value]\n",
      "prepare requires --chrom; cluster/report use --chroms (or a single --chrom).\n",
      "cluster accepts --state closed|open_fire|open_linker|other|all.\n",
      "Exactly three final subclusters per sufficiently variable broad state.\n", sep = "")
  print(opts)
  quit(status = 0L)
}
if (length(argv) %% 2L) stop("Options must be --name value pairs")
for (i in seq_len(length(argv) %/% 2L) * 2L - 1L) {
  name <- sub("^--", "", argv[i])
  if (!startsWith(argv[i], "--") || !name %in% names(opts)) stop("Unknown option: ", argv[i])
  opts[[name]] <- argv[i + 1L]
}
stopifnot(opts$stage %in% c("prepare", "cluster", "report"))
broad_states <- c("closed", "open_fire", "open_linker", "other")
stopifnot(opts$state %in% c(broad_states, "all"))
chroms <- if (nzchar(opts$chrom)) opts$chrom else strsplit(opts$chroms, ",", fixed = TRUE)[[1]]
stopifnot(length(chroms) > 0L, !anyDuplicated(chroms), all(chroms %in% paste0("chr", c(1:22, "X", "Y"))))
lo <- as.integer(opts[["window-start"]]); hi <- as.integer(opts[["window-end"]])
bin_size <- as.integer(opts[["bin-size"]]); nuc_min <- as.integer(opts[["nuc-min"]])
k <- as.integer(opts$k); trees <- as.integer(opts$trees); search_k <- as.integer(opts[["search-k"]])
seed <- as.integer(opts$seed); resolution <- as.numeric(opts$resolution)
validate_size <- as.integer(opts[["validate-size"]]); plot_n <- as.integer(opts[["plot-per-cluster"]])
stopifnot(!anyNA(c(lo, hi, bin_size, nuc_min, k, trees, search_k, seed, validate_size, plot_n)),
  lo <= -10L, hi >= 10L, bin_size > 0L, nuc_min > 0L, k >= 2L,
  trees > 0L, search_k >= k + 1L, validate_size >= 0L, plot_n > 0L,
  is.finite(resolution), resolution > 0)
channels <- c("nuc_ge_min", "nuc_short", "fire", "linker")
bins <- data.table(rel_start = seq.int(lo, hi, by = bin_size))
bins[, rel_end := pmin(rel_start + bin_size - 1L, hi)]
bins[, width := rel_end - rel_start + 1L]
feature_anno <- rbindlist(rep(list(bins), length(channels)))
feature_anno[, channel := rep(channels, each = nrow(bins))]
feature_anno[, feature := paste(channel, rel_start, rel_end, sep = ":")]
distance_weights <- feature_anno$width / (2 * (hi - lo + 1L))
profile_params <- list(window_start = lo, window_end = hi, bin_size = bin_size,
                       nuc_min = nuc_min, channels = channels)
source(file.path(repo, "code/parsing_functions/parsing_footprints_functions.r"))
source(file.path(repo, "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r"))
source(file.path(repo, "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r"))
out_dir <- opts[["out-dir"]]
for (directory in c("profiles", "states", "tables", "plots"))
  dir.create(file.path(out_dir, directory), recursive = TRUE, showWarnings = FALSE)

if (opts$stage == "prepare") {
  if (!nzchar(opts$chrom)) stop("prepare requires --chrom, one task per chromosome")
  if (!file.exists(opts$genes) || file.access(opts$genes, 4L) != 0L)
    stop("Missing or unreadable --genes table: ", opts$genes,
      "\nThe default is an upstream Step 2 output; this workflow does not create it.",
      "\nPass --genes with the same anchor table used to generate the input fiber_features table.",
      "\nRequired columns: gene_id, chrom, strand, anchor_tss (1-based).")
  chrom <- opts$chrom
  input <- file.path(opts[["feature-dir"]], paste0("fiber_features.", chrom, ".tsv.gz"))
  fibers <- fread(cmd = paste("gzip -dc --", shQuote(input)), na.strings = "NA")
  input_columns <- names(fibers)
  required <- c("gene_id", "sample_name", "rid", "read_start", "read_end", "promoter_state")
  stopifnot(all(required %in% input_columns), !anyNA(fibers[, ..required]),
            !anyDuplicated(fibers[, .(gene_id, sample_name, rid)]),
            all(fibers$promoter_state %in% broad_states))
  reserved <- c(".key", ".row", ".chrom", "profile_status", "leiden_community", "promoter_subcluster", "clustering_status")
  if (any(reserved %in% input_columns)) stop("Expected original Step 3 input, without clustering/internal columns")
  samples <- fread(opts$samples)
  genes <- fread(opts$genes)
  stopifnot(all(c("sample_name", "fire_dir") %in% names(samples)), !anyDuplicated(samples$sample_name),
            all(c("gene_id", "chrom", "strand", "anchor_tss") %in% names(genes)))
  target_chrom <- chrom
  genes <- genes[chrom == target_chrom]
  stopifnot(!anyDuplicated(genes$gene_id), all(fibers$gene_id %in% genes$gene_id),
            all(fibers$sample_name %in% samples$sample_name), all(genes$strand %in% c("+", "-")))
  fibers[, `:=`(.key = paste(gene_id, sample_name, rid, sep = "::"),
                 .row = seq_len(.N), .chrom = chrom, profile_status = "segmentation_missing")]
  if (anyDuplicated(fibers$.key)) stop("Composite key collision")
  profiles <- matrix(NA_real_, nrow(fibers), nrow(feature_anno),
                     dimnames = list(fibers$.key, feature_anno$feature))
  for (sample_value in unique(fibers$sample_name)) {
    sample_row <- samples[sample_name == sample_value]
    path <- file.path(sample_row$fire_dir, "extracted_results", paste0(sample_value, ".fire_all.bed.gz"))
    if (!file.exists(path) || !file.exists(paste0(path, ".tbi"))) stop("Missing FIRE input/index: ", path)
    if (!chrom %in% seqnamesTabix(path)) next
    for (gene in unique(fibers[sample_name == sample_value, gene_id])) {
      ids <- fibers[sample_name == sample_value & gene_id == gene, .row]
      anchor <- genes[gene_id == gene]
      tss <- anchor$anchor_tss
      if (!is.finite(tss) || tss %% 1 != 0L) stop("Invalid anchor for ", gene)
      ws <- tss + if (anchor$strand == "+") lo else -hi
      we <- tss + if (anchor$strand == "+") hi else -lo
      if (ws < 1L || any(fibers$read_start[ids] > ws - 1L | fibers$read_end[ids] < we))
        stop("Input primary spans do not cover the configured window: ", gene, "/", sample_value)
      region <- data.frame(start = ws - 1L, end = we, analysis_start = ws,
                           analysis_end = we, width = we - ws + 1L)
      segments <- unique(as.data.table(read_tabix_region(path, GRanges(chrom, IRanges(ws, we)))))
      if (!nrow(segments)) next
      if (ncol(segments) != 11L) stop("Expected 11 FIRE columns: ", path)
      setnames(segments, FIRE_BED_COLS)
      segments <- segments[RID %in% fibers$rid[ids]]
      if (!nrow(segments)) next
      if (anyNA(segments[, .(start, end, rgb)]) || any(segments$end <= segments$start))
        stop("Invalid FIRE interval: ", path)
      # Length classification precedes clipping, exactly as in Step 3.
      segments[, channel := fifelse(rgb == RGB_NUCLEOSOME,
        fifelse(end - start >= nuc_min, "nuc_ge_min", "nuc_short"),
        fifelse(rgb == RGB_LINKER, "linker", "fire"))]
      assignments <- data.frame(RID = fibers$rid[ids])
      coverage <- matrix(0L, length(ids), hi - lo + 1L)
      local <- matrix(0, length(ids), nrow(feature_anno))
      broad_hits <- matrix(FALSE, length(ids), length(channels), dimnames = list(NULL, channels))
      for (ci in seq_along(channels)) {
        # Existing helper unions intervals, so same-class duplicates count once.
        occupancy <- occupancy_matrix(segments[channel == channels[ci]], assignments, region)
        if (anchor$strand == "-") occupancy <- occupancy[, ncol(occupancy):1L, drop = FALSE]
        coverage <- coverage + occupancy
        broad_hits[, ci] <- rowSums(occupancy[, seq.int(-10L, 10L) - lo + 1L, drop = FALSE]) > 0L
        for (bi in seq_len(nrow(bins))) {
          columns <- seq.int(bins$rel_start[bi], bins$rel_end[bi]) - lo + 1L
          local[, (ci - 1L) * nrow(bins) + bi] <- rowMeans(occupancy[, columns, drop = FALSE])
        }
      }
      status <- rep("ready", length(ids))
      status[rowSums(coverage == 0L) > 0L] <- "segmentation_incomplete"
      status[rowSums(coverage) == 0L] <- "segmentation_missing"
      status[rowSums(coverage > 1L) > 0L] <- "segmentation_conflict"
      calculated <- ifelse(broad_hits[, "nuc_ge_min"], "closed",
        ifelse(broad_hits[, "fire"], "open_fire", ifelse(broad_hits[, "linker"], "open_linker", "other")))
      status[status == "ready" & calculated != fibers$promoter_state[ids]] <- "broad_state_mismatch"
      valid <- which(status == "ready")
      if (length(valid)) profiles[ids[valid], ] <- local[valid, , drop = FALSE]
      set(fibers, i = ids, j = "profile_status", value = status)
    }
    message(chrom, " ", sample_value, ": profiles prepared")
  }
  metadata_path <- file.path(out_dir, "profiles", paste0("rows.", chrom, ".tsv.gz"))
  profile_path <- file.path(out_dir, "profiles", paste0("profiles.", chrom, ".rds"))
  saveRDS(list(profiles = profiles, keys = fibers$.key, annotation = feature_anno,
               params = profile_params, input_columns = input_columns, source = normalizePath(input)),
          paste0(profile_path, ".tmp"), compress = FALSE)
  fwrite(fibers, paste0(metadata_path, ".tmp.gz"), sep = "\t", na = "NA")
  stopifnot(file.rename(paste0(profile_path, ".tmp"), profile_path),
            file.rename(paste0(metadata_path, ".tmp.gz"), metadata_path))
  fwrite(fibers[, .N, by = .(sample_name, promoter_state, profile_status)],
         file.path(out_dir, "profiles", paste0("qc.", chrom, ".tsv")), sep = "\t")
  message("Prepared ", nrow(fibers), " rows; ", sum(fibers$profile_status == "ready"), " usable profiles")
  quit(status = 0L)
}

# Both clustering and reporting require the explicitly selected chromosome set.
metadata_paths <- file.path(out_dir, "profiles", paste0("rows.", chroms, ".tsv.gz"))
profile_paths <- file.path(out_dir, "profiles", paste0("profiles.", chroms, ".rds"))
if (!all(file.exists(c(metadata_paths, profile_paths)))) stop("Prepare every selected chromosome first")
signature_paths <- normalizePath(c(metadata_paths, profile_paths))
signature_info <- file.info(signature_paths)
signature <- data.frame(path = signature_paths, size = signature_info$size,
                         mtime = as.numeric(signature_info$mtime))

if (opts$stage == "cluster") {
  if (!requireNamespace("RcppAnnoy", quietly = TRUE)) stop("RcppAnnoy is required")
  states <- if (opts$state == "all") broad_states else opts$state
  for (state in states) {
    message("Building shared subclusters for ", state)
    rows <- list()
    for (ci in seq_along(chroms)) {
      tab <- fread(cmd = paste("gzip -dc --", shQuote(metadata_paths[ci])), na.strings = "NA")
      rows[[ci]] <- tab[promoter_state == state & profile_status == "ready", .(.key, .row, .chrom)]
    }
    rows <- rbindlist(rows)
    setorder(rows, .key)
    stopifnot(!anyDuplicated(rows$.key))
    mat <- matrix(NA_real_, nrow(rows), nrow(feature_anno), dimnames = list(rows$.key, feature_anno$feature))
    for (ci in seq_along(chroms)) {
      loaded <- readRDS(profile_paths[ci])
      if (!identical(loaded$params, profile_params) ||
          !identical(as.data.frame(loaded$annotation), as.data.frame(feature_anno)))
        stop("Profile settings differ; use the same window/bin/nuc options as prepare")
      dest <- which(rows$.chrom == chroms[ci])
      if (length(dest)) {
        stopifnot(identical(loaded$keys[rows$.row[dest]], rows$.key[dest]))
        mat[dest, ] <- loaded$profiles[rows$.row[dest], , drop = FALSE]
      }
      rm(loaded)
    }
    stopifnot(!anyNA(mat), all(is.finite(mat)))
    for (j in seq_len(ncol(mat))) mat[, j] <- mat[, j] * distance_weights[j]
    # Probe distinct rows in bounded blocks instead of stringifying a huge matrix.
    probes <- mat[integer(), , drop = FALSE]
    for (first in seq.int(1L, max(1L, nrow(mat)), by = 1000L)) {
      if (!nrow(mat)) break
      last <- min(nrow(mat), first + 999L)
      probes <- unique(rbind(probes, mat[first:last, , drop = FALSE]))
      if (nrow(probes) >= 3L) break
    }
    assigned <- rows[, .(.key, leiden_community = NA_character_, promoter_subcluster = NA_character_,
                         clustering_status = "insufficient_variation")]
    state_path <- file.path(out_dir, "states", state)
    info <- list(state = state, signature = signature, profile_params = profile_params,
                 options = opts, n_reads = nrow(mat), status = "insufficient_variation")
    if (nrow(probes) < 3L) {
      if (nrow(assigned)) {
        fwrite(assigned, paste0(state_path, ".assignments.tsv.gz"), sep = "\t", na = "NA")
      } else {
        connection <- gzfile(paste0(state_path, ".assignments.tsv.gz"), "wt")
        writeLines(".key\tleiden_community\tpromoter_subcluster\tclustering_status", connection)
        close(connection)
      }
      saveRDS(info, paste0(state_path, ".model.rds"))
      message(state, ": fewer than three distinct profiles; no artificial groups assigned")
      next
    }
    # Approximate graph construction stays inline here; existing graph helpers
    # are unchanged. A separate small validation graph is never used for fitting.
    scopes <- c(if (validate_size >= 3L) "validation", "main")
    validation <- NULL
    for (scope in scopes) {
      if (scope == "validation") {
        set.seed(seed)
        select <- sort(sample.int(nrow(mat), min(validate_size, nrow(mat))))
        x <- mat[select, , drop = FALSE]
      } else x <- mat
      n <- nrow(x); k_eff <- min(k, n - 1L)
      index <- new(RcppAnnoy::AnnoyManhattan, as.integer(ncol(x)))
      index$setSeed(seed)
      for (i in seq_len(n)) index$addItem(as.integer(i - 1L), as.numeric(x[i, ]))
      index$build(trees)
      neighbors <- matrix(NA_integer_, n, k_eff)
      distances <- matrix(NA_real_, n, k_eff)
      for (i in seq_len(n)) {
        hit <- index$getNNsByItemList(as.integer(i - 1L), as.integer(k_eff + 1L), search_k, TRUE)
        candidates <- hit$item + 1L
        keep <- which(candidates != i & !duplicated(candidates))
        keep <- keep[order(hit$distance[keep], candidates[keep])]
        if (length(keep) < k_eff) stop("Too few Annoy neighbors; increase --search-k")
        keep <- keep[seq_len(k_eff)]
        neighbors[i, ] <- candidates[keep]
        distances[i, ] <- hit$distance[keep]
        if (i %% 100000L == 0L) message(state, " ", scope, ": ", i, "/", n, " neighbor queries")
      }
      index$unload()
      sigma <- mean(distances)
      if (!is.finite(sigma) || sigma <= 0) sigma <- 1
      from <- rep(seq_len(n), each = k_eff); to <- as.vector(t(neighbors))
      a <- pmin(from, to); b <- pmax(from, to)
      keep <- !duplicated(a * as.double(n) + b)
      edge_distance <- as.vector(t(distances))[keep]
      graph <- igraph::graph_from_data_frame(data.frame(from = a[keep], to = b[keep],
        weight = exp(-edge_distance / sigma), manhattan_distance = edge_distance), directed = FALSE,
        vertices = data.frame(name = as.character(seq_len(n)), key = rownames(x)))
      part <- leiden_partition(graph, resolution = resolution, seed = seed)
      initial <- as.integer(part$membership[as.character(seq_len(n))])
      if (scope == "validation") {
        exact <- manhattan_knn_graph(x, k_neighbors = k, verbose = FALSE)
        recall <- radius_fraction <- numeric(n)
        for (i in seq_len(n)) {
          exact_ids <- order(exact$dist[i, ])[seq_len(k_eff)]
          recall[i] <- length(intersect(exact_ids, neighbors[i, ])) / k_eff
          radius_fraction[i] <- mean(exact$dist[i, neighbors[i, ]] <= max(exact$dist[i, exact_ids]) + 1e-10)
        }
        exact_part <- leiden_partition(exact$graph, resolution = resolution, seed = seed)
        validation <- data.table(promoter_state = state, n_validation = n,
          mean_neighbor_recall = mean(recall), fraction_within_exact_k_radius = mean(radius_fraction),
          adjusted_rand = igraph::compare(initial, as.integer(exact_part$membership[rownames(x)]), method = "adjusted.rand"))
        fwrite(validation, paste0(state_path, ".neighbor_validation.tsv"), sep = "\t")
        rm(exact, exact_part)
      }
      rm(index, from, to, a, b, neighbors, distances, edge_distance, x)
    }
    # Save graph separately so report does not need to deserialize its edges.
    saveRDS(graph, paste0(state_path, ".graph.rds"), compress = FALSE)
    graph_diag <- list(vertices = igraph::vcount(graph), edges = igraph::ecount(graph),
      components = igraph::components(graph)$no, sigma = sigma, k_eff = k_eff,
      leiden_quality = part$quality, initial_communities = length(unique(initial)))
    rm(graph)
    membership <- initial
    adjustments <- list()
    if (length(unique(membership)) > 3L) {
      levels_initial <- levels(factor(membership))
      centers <- cluster_feature_profiles(mat, factor(membership, levels = levels_initial))
      rownames(centers) <- levels_initial
      counts <- setNames(as.numeric(table(factor(membership, levels = levels_initial))), levels_initial)
      mapping <- setNames(levels_initial, levels_initial)
      while (nrow(centers) > 3L) {
        d <- as.matrix(dist(centers, method = "manhattan"))
        d[lower.tri(d, diag = TRUE)] <- Inf
        pair <- arrayInd(which.min(d), dim(d))
        keep_id <- rownames(centers)[pair[1]]; drop_id <- rownames(centers)[pair[2]]
        adjustments[[length(adjustments) + 1L]] <- data.table(action = "merge", from = drop_id, into = keep_id)
        centers[keep_id, ] <- (centers[keep_id, ] * counts[keep_id] + centers[drop_id, ] * counts[drop_id]) /
          (counts[keep_id] + counts[drop_id])
        counts[keep_id] <- counts[keep_id] + counts[drop_id]
        centers <- centers[rownames(centers) != drop_id, , drop = FALSE]
        counts <- counts[names(counts) != drop_id]
        mapping[mapping == drop_id] <- keep_id
      }
      membership <- as.integer(mapping[as.character(initial)])
    }
    while (length(unique(membership)) < 3L) {
      ids <- sort(unique(membership))
      dispersion <- numeric(length(ids))
      for (ci in seq_along(ids)) {
        positions <- which(membership == ids[ci])
        center <- colMeans(mat[positions, , drop = FALSE])
        for (first in seq.int(1L, length(positions), by = 2048L)) {
          chunk <- positions[first:min(length(positions), first + 2047L)]
          dispersion[ci] <- dispersion[ci] + sum(abs(sweep(mat[chunk, , drop = FALSE], 2L, center, "-")))
        }
      }
      selected <- ids[which.max(dispersion)]
      if (max(dispersion) <= 0) stop("Cannot split identical profiles into meaningful groups")
      positions <- which(membership == selected)
      center <- colMeans(mat[positions, , drop = FALSE])
      delta <- numeric(length(positions))
      for (first in seq.int(1L, length(positions), by = 2048L)) {
        chunk <- first:min(length(positions), first + 2047L)
        delta[chunk] <- rowSums(abs(sweep(mat[positions[chunk], , drop = FALSE], 2L, center, "-")))
      }
      seed1 <- positions[which.max(delta)]
      for (first in seq.int(1L, length(positions), by = 2048L)) {
        chunk <- first:min(length(positions), first + 2047L)
        delta[chunk] <- rowSums(abs(sweep(mat[positions[chunk], , drop = FALSE], 2L, mat[seed1, ], "-")))
      }
      seed2 <- positions[which.max(delta)]
      new_id <- max(membership) + 1L
      for (first in seq.int(1L, length(positions), by = 2048L)) {
        chunk <- first:min(length(positions), first + 2047L)
        second <- rowSums(abs(sweep(mat[positions[chunk], , drop = FALSE], 2L, mat[seed2, ], "-")))
        membership[positions[chunk[second < delta[chunk]]]] <- new_id
      }
      adjustments[[length(adjustments) + 1L]] <- data.table(action = "split", from = as.character(selected), into = as.character(new_id))
    }
    ordering <- data.table(group = membership, read_key = rows$.key)[, .(n = .N, first_key = min(read_key)), by = group]
    setorder(ordering, -n, first_key)
    ordering[, label := paste0(state, "_S", seq_len(.N))]
    labels <- ordering$label[match(membership, ordering$group)]
    stopifnot(length(unique(labels)) == 3L, !anyNA(labels))
    assigned[, `:=`(leiden_community = paste0(state, "_L", initial),
                    promoter_subcluster = labels, clustering_status = "assigned")]
    means <- cluster_feature_profiles(mat, factor(labels, levels = ordering$label))
    means <- sweep(means, 2L, distance_weights, "/")
    rownames(means) <- ordering$label
    mean_table <- rbindlist(rep(list(feature_anno), 3L))
    mean_table[, `:=`(promoter_state = state, promoter_subcluster = rep(ordering$label, each = ncol(mat)),
                     n_reads = rep(ordering$n, each = ncol(mat)), mean_occupancy = as.vector(t(means)))]
    # Representative = closest fiber to its mean profile, not a claimed medoid.
    representatives <- display <- list()
    set.seed(seed)
    for (ci in seq_len(3L)) {
      positions <- which(labels == ordering$label[ci])
      delta <- numeric(length(positions))
      for (first in seq.int(1L, length(positions), by = 2048L)) {
        chunk <- first:min(length(positions), first + 2047L)
        center <- means[ci, ] * distance_weights
        delta[chunk] <- rowSums(abs(sweep(mat[positions[chunk], , drop = FALSE], 2L, center, "-")))
      }
      chosen <- positions[order(delta, rows$.key[positions])[seq_len(min(3L, length(positions)))]]
      representatives[[ci]] <- list(keys = rows$.key[chosen], profiles = sweep(mat[chosen, , drop = FALSE], 2L, distance_weights, "/"))
      display[[ci]] <- sort(positions[sample.int(length(positions), min(plot_n, length(positions)))])
    }
    display_rows <- unlist(display, use.names = FALSE)
    info$status <- "assigned"
    info$graph <- graph_diag
    info$validation <- validation
    info$centers <- means
    info$mean_profiles <- mean_table
    info$counts <- ordering
    info$adjustments <- if (length(adjustments)) rbindlist(adjustments) else data.table(action = character(), from = character(), into = character())
    info$representatives <- representatives
    info$display <- list(keys = rows$.key[display_rows], labels = labels[display_rows],
      profiles = sweep(mat[display_rows, , drop = FALSE], 2L, distance_weights, "/"))
    fwrite(assigned, paste0(state_path, ".assignments.tsv.gz"), sep = "\t", na = "NA")
    fwrite(mean_table, paste0(state_path, ".mean_profiles.tsv"), sep = "\t")
    fwrite(info$adjustments, paste0(state_path, ".adjustments.tsv"), sep = "\t")
    saveRDS(info, paste0(state_path, ".model.rds"))
    message(state, ": ", nrow(rows), " fibers -> ", paste(ordering$label, ordering$n, collapse = "; "))
    rm(mat, rows, probes)
  }
  quit(status = 0L)
}

if (opts$stage == "report") {
  models <- assignments <- list()
  for (state in broad_states) {
    state_path <- file.path(out_dir, "states", state)
    model <- readRDS(paste0(state_path, ".model.rds"))
    if (!identical(model$signature, signature) || !identical(model$profile_params, profile_params))
      stop("Stale/different chromosome or profile inputs for ", state, "; rerun clustering")
    models[[state]] <- model
    assignments[[state]] <- fread(cmd = paste("gzip -dc --", shQuote(paste0(state_path, ".assignments.tsv.gz"))), na.strings = "NA")
  }
  assignments <- rbindlist(assignments)
  stopifnot(!anyDuplicated(assignments$.key))
  summaries <- sample_composition <- gene_composition <- joint_composition <- qc <- list()
  metrics <- c("tss_fire_score", "tss_msp_len", "tss_nuc_start", "tss_nuc_size",
    "plus1_start", "plus1_size", "plus2_start", "plus2_size", "plus3_start", "plus3_size",
    "minus1_start", "minus1_end", "minus1_size", "ndr_len", "down_n_nuc", "down_dyad_spacing",
    "down_acf_lag", "down_acf_peak", "pic", "ppp", "tf_up_n", "elong_n", "elong_per_kb",
    "frac_nuc", "frac_msp", "m6a_n", "m6a_site_n", "mcpg_n", "cpg_site_n", "m6a_observed_fraction", "cpg_observed_fraction")
  for (ci in seq_along(chroms)) {
    tab <- fread(cmd = paste("gzip -dc --", shQuote(metadata_paths[ci])), na.strings = "NA")
    original_columns <- setdiff(names(tab), c(".key", ".row", ".chrom", "profile_status"))
    hit <- match(tab$.key, assignments$.key)
    if (any(tab$profile_status == "ready" & is.na(hit))) stop("Missing assignments for prepared rows")
    tab[, `:=`(leiden_community = assignments$leiden_community[hit],
                promoter_subcluster = assignments$promoter_subcluster[hit], clustering_status = profile_status)]
    tab[!is.na(hit), clustering_status := assignments$clustering_status[hit[!is.na(hit)]]]
    output_columns <- c(original_columns, "leiden_community", "promoter_subcluster", "clustering_status")
    fwrite(tab[, ..output_columns], file.path(out_dir, "tables", paste0("fiber_features.clustered.", chroms[ci], ".tsv.gz")), sep = "\t", na = "NA")
    qc[[ci]] <- tab[, .(n_reads = .N), by = .(promoter_state, clustering_status)]
    sample_composition[[ci]] <- tab[, .(n_reads = .N), by = .(sample_name, promoter_state, promoter_subcluster, clustering_status)]
    gene_composition[[ci]] <- tab[, .(n_reads = .N), by = .(gene_id, promoter_state, promoter_subcluster, clustering_status)]
    joint_composition[[ci]] <- tab[, .(n_reads = .N), by = .(gene_id, sample_name, promoter_state, promoter_subcluster, clustering_status)]
    if (all(c("m6a_n", "m6a_site_n") %in% names(tab)))
      tab[, m6a_observed_fraction := fifelse(m6a_site_n > 0, m6a_n / m6a_site_n, NA_real_)]
    if (all(c("mcpg_n", "cpg_site_n") %in% names(tab)))
      tab[, cpg_observed_fraction := fifelse(cpg_site_n > 0, mcpg_n / cpg_site_n, NA_real_)]
    present <- intersect(metrics, names(tab))
    summaries[[ci]] <- tab[clustering_status == "assigned", c("gene_id", "sample_name", "promoter_state", "promoter_subcluster", present), with = FALSE]
  }
  for (kind in c("sample", "gene", "gene_sample")) {
    pieces <- switch(kind, sample = sample_composition, gene = gene_composition, gene_sample = joint_composition)
    group <- switch(kind, sample = "sample_name", gene = "gene_id", gene_sample = c("gene_id", "sample_name"))
    totals <- rbindlist(pieces)[, .(n_reads = sum(n_reads)), by = c(group, "promoter_state", "promoter_subcluster", "clustering_status")]
    totals[, fraction_all_rows := n_reads / sum(n_reads), by = group]
    totals[, fraction_within_broad := n_reads / sum(n_reads), by = c(group, "promoter_state")]
    fwrite(totals, file.path(out_dir, "tables", paste0("composition_by_", kind, ".tsv")), sep = "\t", na = "NA")
  }
  qc <- rbindlist(qc)[, .(n_reads = sum(n_reads)), by = .(promoter_state, clustering_status)]
  fwrite(qc, file.path(out_dir, "tables", "clustering_qc.tsv"), sep = "\t")
  summaries <- rbindlist(summaries, fill = TRUE)
  counts <- summaries[, .(n_reads = .N, n_genes = uniqueN(gene_id), n_samples = uniqueN(sample_name)),
                       by = .(promoter_state, promoter_subcluster)]
  fwrite(counts, file.path(out_dir, "tables", "subcluster_counts.tsv"), sep = "\t")
  stats <- list()
  for (metric in intersect(metrics, names(summaries))) {
    values <- summaries[, .(promoter_state, promoter_subcluster, value = as.numeric(get(metric)))]
    stats[[metric]] <- values[, .(metric = metric, n_observed = sum(is.finite(value)),
      mean = if (any(is.finite(value))) mean(value[is.finite(value)]) else NA_real_,
      median = if (any(is.finite(value))) median(value[is.finite(value)]) else NA_real_,
      q25 = if (any(is.finite(value))) unname(quantile(value[is.finite(value)], .25)) else NA_real_,
      q75 = if (any(is.finite(value))) unname(quantile(value[is.finite(value)], .75)) else NA_real_),
      by = .(promoter_state, promoter_subcluster)]
  }
  fwrite(rbindlist(stats), file.path(out_dir, "tables", "subcluster_feature_summary.tsv"), sep = "\t", na = "NA")
  # Small plotting subsets are for display only; all rows participated in fitting.
  colors <- c(nuc_ge_min = "#666666", nuc_short = "#BBBBBB", fire = "#E69F00", linker = "#9370DB")
  cluster_colors <- c("#0072B2", "#D55E00", "#009E73")
  pdf(file.path(out_dir, "plots", "state_profiles.pdf"), width = 12, height = 10)
  for (state in broad_states) {
    model <- models[[state]]
    if (model$status != "assigned") next
    par(mfrow = c(4, 1), mar = c(3, 4, 3, 1))
    for (channel_value in channels) {
      data <- model$mean_profiles[channel == channel_value]
      plot(range(c(lo, hi)), c(0, 1), type = "n", xlab = "Position relative to TSS (bp)",
           ylab = "Mean occupancy", main = paste(state, channel_value))
      abline(v = 0, lty = 2, col = "grey50")
      for (j in seq_len(3L)) {
        sub <- data[promoter_subcluster == model$counts$label[j]]
        lines((sub$rel_start + sub$rel_end) / 2, sub$mean_occupancy, col = cluster_colors[j], lwd = 2)
      }
      legend("topright", model$counts$label, col = cluster_colors, lty = 1, bty = "n")
    }
  }
  dev.off()
  pdf(file.path(out_dir, "plots", "fiber_profile_heatmaps.pdf"), width = 14, height = 11)
  for (state in broad_states) {
    model <- models[[state]]
    if (model$status != "assigned") next
    par(mfrow = c(4, 1), mar = c(3, 5, 3, 1))
    for (channel in channels) {
      columns <- which(feature_anno$channel == channel)
      z <- model$display$profiles[, columns, drop = FALSE]
      image(x = seq_len(ncol(z)), y = seq_len(nrow(z)), z = t(z), zlim = c(0, 1),
        col = colorRampPalette(c("white", colors[channel]))(100), axes = FALSE,
        xlab = "TSS-relative position (bp)", ylab = "Displayed fibers", main = paste(state, channel))
      ticks <- unique(c(1L, which.min(abs(bins$rel_start)), seq.int(1L, nrow(bins), by = max(1L, floor(nrow(bins) / 6))), nrow(bins)))
      axis(1, at = ticks, labels = bins$rel_start[ticks])
      lengths <- rle(model$display$labels)$lengths
      boundaries <- cumsum(lengths)
      abline(h = head(boundaries, -1L) + .5, col = "black")
      axis(2, at = boundaries - lengths / 2 + .5, labels = rle(model$display$labels)$values, las = 1, cex.axis = .7)
      box()
    }
  }
  dev.off()
  pdf(file.path(out_dir, "plots", "representative_fibers.pdf"), width = 12, height = 10)
  for (state in broad_states) {
    model <- models[[state]]
    if (model$status != "assigned") next
    for (ci in seq_len(3L)) {
      selected <- model$representatives[[ci]]
      for (ri in seq_len(nrow(selected$profiles))) {
        par(mfrow = c(4, 1), mar = c(3, 4, 3, 1))
        for (channel in channels) {
          columns <- which(feature_anno$channel == channel)
          plot((bins$rel_start + bins$rel_end) / 2, selected$profiles[ri, columns], type = "s", ylim = c(0, 1),
            col = colors[channel], xlab = "Position relative to TSS (bp)", ylab = channel,
            main = paste(model$counts$label[ci], selected$keys[ri]))
          abline(v = 0, lty = 2, col = "grey50")
        }
      }
    }
  }
  dev.off()
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "tables", "session_info.txt"))
  message("Wrote enriched promoter tables, summaries, composition tables, and profile plots")
}
