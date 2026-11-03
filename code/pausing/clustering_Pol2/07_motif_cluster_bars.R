#!/usr/bin/env Rscript
# Preserve the original 5000-read cluster bars; reuse exported cohort/blocks.
project <- "/project/spott/cshan/fiber-seq"
base <- file.path(project,"code/pausing/clustering_Pol2")
cfg <- yaml::read_yaml(file.path(base,"pol2_chip_cluster_utils/config.yaml"))
old_tables <- file.path(cfg$output,"tables/motif_puffin_scan")
if (!all(file.exists(file.path(old_tables,c("cohort.tsv.gz","same_molecule_footprints.tsv.gz","assigned_promoters.tsv"))))) {
# Export only the original saved cohort; recover available footprint quality.
project <- "/project/spott/cshan/fiber-seq"
base <- file.path(project, "code/pausing/clustering_Pol2")
source(file.path(base, "pol2_chip_cluster_utils/pol2_clustering_functions.R"))
source(file.path(base, "pol2_chip_cluster_utils/pol2_molecule_tables.R"))
source(file.path(project, "code/parsing_functions/parsing_footprints_functions.r"))
cfg <- yaml::read_yaml(file.path(base, "pol2_chip_cluster_utils/config.yaml"))
tables <- file.path(cfg$output, "tables")
out <- file.path(tables, "motif_puffin_scan")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
# Deliberately use saved preparation, not current configuration to regenerate it.
prep <- readRDS(file.path(tables, "prepared.rds"))
x <- pol2_read_clustering(cfg, prep)
a <- fread(file.path(tables, "molecule_assignments.tsv.gz"))
a <- a[match(x$molecules$observation_id, observation_id)]
stopifnot(!anyNA(a$observation_id), !anyDuplicated(a$observation_id),
          identical(a$cluster, as.character(x$molecules$cluster)),
          identical(a$sample_name, x$molecules$sample_name),
          identical(a$original_RID, x$molecules$original_RID),
          identical(a$tss_id, x$molecules$tss_id))
a[, saved_order := .I]
p <- prep$tss[tss_id %in% a$tss_id]
stopifnot(all(a$tss_id %in% p$tss_id))
fwrite(p, file.path(out, "assigned_promoters.tsv"), sep = "\t")
fwrite(prep$samples, file.path(out, "samples.tsv"), sep = "\t")
tracks_path <- file.path(tables, "tracks/molecule_footprints.tsv.gz")
pol2_required(tracks_path)
tracks <- fread(tracks_path)[observation_id %in% a$observation_id]
fp <- unique(tracks[startsWith(class, "TF footprint") & length >= 10 & length <= 80])
fp <- merge(fp, a[, .(observation_id, sample_name, original_RID, tss_id)], by = "observation_id")
fp[, footprint_id := paste(sample_name, original_RID, chrom, start, end, sep = "::")]
fp[, `:=`(block_score = NA_real_, block_edge_left = NA_real_, block_edge_right = NA_real_,
          read_score = NA_real_, quality_status = "unavailable",
          inside_same_read_fire = FALSE)]
a[, `:=`(tf_footprint_evaluable = FALSE, tf_footprint_input_status = "unresolved",
          fire_evaluable = FALSE, fire_input_status = "unresolved")]
source_info <- list()
for (sample_id in unique(a$sample_name)) {
  sm <- prep$samples[sample_name == sample_id]
  path <- pol2_tf_path(sm, cfg$footprint_dir)
  fire_path <- file.path(sm$fire_dir, "additional-outputs-v0.1/fire-peaks",
                         paste0(sample_id, "-v0.1-fire-elements.bed.gz"))
  ids <- a[sample_name == sample_id]
  source_info[[paste0(sample_id, "::tf")]] <- data.table(
    sample_name = sample_id, source_type = "fiberhmm_tf", path = path,
    exists = file.exists(path), index_exists = file.exists(paste0(path, ".tbi")))
  source_info[[paste0(sample_id, "::fire")]] <- data.table(
    sample_name = sample_id, source_type = "fire_elements", path = fire_path,
    exists = file.exists(fire_path), index_exists = file.exists(paste0(fire_path, ".tbi")))
  pp <- p[tss_id %in% ids$tss_id]
  query <- GenomicRanges::reduce(pol2_gr(pp$chrom, pp$window_start, pp$window_end))
  if (!file.exists(path) || !file.exists(paste0(path, ".tbi"))) {
    a[sample_name == sample_id, tf_footprint_input_status := "missing_source_or_index"]
    message("Missing footprint source or index: ", path)
  } else {
    raw <- read_tabix_region(path, query)
    if (!nrow(raw)) {
      a[sample_name == sample_id, tf_footprint_input_status := "source_read_record_absent"]
    } else {
      raw <- raw[raw[[4]] %in% ids$original_RID, , drop = FALSE]
      found <- unique(raw[[4]])
      # Missing read records remain unresolved, not zero-footprint observations.
      a[sample_name == sample_id & original_RID %in% found,
        `:=`(tf_footprint_evaluable = TRUE, tf_footprint_input_status = "source_read_record_found")]
      a[sample_name == sample_id & !original_RID %in% found,
        tf_footprint_input_status := "source_read_record_absent"]
      if (nrow(raw)) {
        q <- as.data.table(convert_ft_bed12_to_bed6(raw, format = "bed15_fiberhmm_tf",
                          keep_block_scores = TRUE, source = path))
        # Explicit BED15 parsing preserves input alignment/block order and all blocks.
        # Carry the additional per-block source fields exactly; filter below.
        edge_left <- strsplit(sub(",$", "", as.character(raw[[14]])), ",", fixed = TRUE)
        edge_right <- strsplit(sub(",$", "", as.character(raw[[15]])), ",", fixed = TRUE)
        stopifnot(all(lengths(edge_left) == raw[[10]]), all(lengths(edge_right) == raw[[10]]))
        q[, `:=`(block_edge_left = as.numeric(unlist(edge_left, use.names = FALSE)),
                 block_edge_right = as.numeric(unlist(edge_right, use.names = FALSE)))]
        q <- unique(q[, .(original_RID = RID, chrom = chr, start, end, block_score,
                         block_edge_left, block_edge_right, read_score = as.numeric(score))])
        q[, sample_name := sample_id]
        # Preserve cached full boundaries. Match quality to that exact individual block.
        fp[q, on = .(sample_name, original_RID, chrom, start, end),
           `:=`(block_score = i.block_score, block_edge_left = i.block_edge_left,
                block_edge_right = i.block_edge_right, read_score = i.read_score,
                quality_status = "source_block_score_and_edges")]
      }
      message(sample_id, ": recovered TF source records for ", length(found), "/", nrow(ids), " sampled molecules")
    }
  }

  if (!file.exists(fire_path) || !file.exists(paste0(fire_path, ".tbi"))) {
    a[sample_name == sample_id, fire_input_status := "missing_source_or_index"]
    message("Missing FIRE source or index: ", fire_path)
  } else {
    # The fire-elements file contains positive FIRE calls. With an available
    # file, absence of a read/row is an evaluable zero rather than missing data.
    a[sample_name == sample_id,
      `:=`(fire_evaluable = TRUE, fire_input_status = "source_available_no_fire_element")]
    fire_raw <- read_tabix_region(fire_path, query)
    if (nrow(fire_raw)) {
      fire <- as.data.table(fire_raw)
      stopifnot(ncol(fire) >= length(FIRE_BED_COLS))
      setnames(fire, seq_along(FIRE_BED_COLS), FIRE_BED_COLS)
      fire <- unique(fire[RID %in% ids$original_RID])
      fire[, `:=`(RID = as.character(RID),
                  class = factor(fifelse(rgb == RGB_NUCLEOSOME, "nucleosome",
                    fifelse(rgb == RGB_LINKER, "linker", "FIRE")),
                    levels = c("nucleosome", "linker", "FIRE")))]
      found_fire <- unique(fire$RID)
      a[sample_name == sample_id & original_RID %in% found_fire,
        fire_input_status := "source_available_fire_element_found"]
      candidate_rows <- which(fp$sample_name == sample_id)
      if (length(candidate_rows) && nrow(fire)) {
        candidate_fp <- data.frame(row_id = candidate_rows,
          RID = fp$original_RID[candidate_rows], start = fp$start[candidate_rows],
          end = fp$end[candidate_rows])
        inside <- subset_footprints_in_fire(candidate_fp, as.data.frame(fire))
        fp[inside$row_id, inside_same_read_fire := TRUE]
      }
      message(sample_id, ": FIRE elements found for ", length(found_fire), "/", nrow(ids),
              " sampled molecules; absent calls retained as evaluable zero")
    }
  }
}
a[, `:=`(footprint_evaluable = tf_footprint_evaluable,
          footprint_input_status = tf_footprint_input_status)]
fwrite(a, file.path(out, "cohort.tsv.gz"), sep = "\t")
fwrite(fp, file.path(out, "same_molecule_footprints.tsv.gz"), sep = "\t")
fwrite(rbindlist(source_info), file.path(out, "footprint_sources.tsv"), sep = "\t")
paths <- c(file.path(tables, c("clustering.rds", "prepared.rds", "molecule_assignments.tsv.gz")),
           tracks_path, file.path(base, "pol2_chip_cluster_utils/config.yaml"))
stamp <- pol2_stamp(paths)
stamp$md5 <- unname(tools::md5sum(paths))
fwrite(stamp, file.path(out, "cohort_input_provenance.tsv"), sep = "\t")
cat("Original cohort:", nrow(a), "molecules;", uniqueN(a$tss_id), "TSSs;",
    uniqueN(a$sample_name), "samples\n")

}
if (!file.exists(file.path(old_tables,"matrix_motif_hits.tsv.gz"))) {
  status <- system2(cfg$python,shQuote(c(file.path(base,"05_motif_scan.py"),
    "--promoters",file.path(old_tables,"assigned_promoters.tsv"),"--output",old_tables,
    "--reference",cfg$reference)))
  if (status != 0L) stop("Promoter motif scan failed")
}
# Descriptive same-molecule overlap, supporting tables and cluster bars only.
project <- "/project/spott/cshan/fiber-seq"
base <- file.path(project, "code/pausing/clustering_Pol2")
source(file.path(base, "pol2_chip_cluster_utils/pol2_clustering_functions.R"))
source(file.path(project, "code/parsing_functions/plotting_functions.r"))
source(file.path(base, "pol2_chip_cluster_utils/pol2_clustering_plots.R"))
cfg <- yaml::read_yaml(file.path(base, "pol2_chip_cluster_utils/config.yaml"))
out <- file.path(cfg$output, "tables/motif_puffin_scan")
plots <- file.path(cfg$output, "plots/motif_puffin_scan")
dir.create(plots, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(out, name), sep = "\t", na = "NA")
save_plot <- function(p, name, w = 10, h = 6) {
  save_figure(p, file.path(plots, paste0(name, ".pdf")), width = w, height = h)
}
m <- fread(file.path(out, "cohort.tsv.gz"))
lookup <- fread(file.path(out, "promoter_sequence_lookup.tsv"))
manifest <- fread(file.path(out, "motif_manifest.tsv"))
sites <- fread(file.path(out, "motif_sites.tsv.gz"))
fp_all <- fread(file.path(out, "same_molecule_footprints.tsv.gz"))
# Support uses full individual 10-30-bp blocks, without score/FIRE filters.
fp <- fp_all[length >= 10L & length <= 30L]
m[, footprint_evaluable := tf_footprint_evaluable]
m[, footprint_input_status := tf_footprint_input_status]
stopifnot(nrow(m) == 5000L, !anyDuplicated(m[, .(sample_name, original_RID)]),
          all(m$footprint_evaluable))
categories <- c("SP", "ETS", "ZNF143", "NRF1", "CREB", "TATA", "YY1", "NFY", "U1")
available <- unique(manifest[status == "available", category])
clusters <- pol2_numeric_clusters(m$cluster)
colors <- cluster_id_palette(clusters)
m <- merge(m, lookup[, .(tss_id, sequence_id, usable_sequence, sequence_status)], by = "tss_id", sort = FALSE)
m[, cluster := factor(cluster, levels = clusters)]
m[, n_detected_tf_footprints := 0L]
nf <- fp[, .(n_detected_tf_footprints = uniqueN(footprint_id)), by = observation_id]
m[nf, on = "observation_id", n_detected_tf_footprints := i.n_detected_tf_footprints]
m[footprint_evaluable == FALSE, n_detected_tf_footprints := NA_integer_]
setorder(m, saved_order)
stopifnot(!anyDuplicated(m$observation_id))
site_cols <- c("tss_id", "site_id", "category", "start", "end", "relative_start", "relative_end", "relative_center",
               "genomic_strand", "orientation_relative_to_transcription", "orientations_relative_to_transcription",
               "representative_matrix_id", "matrix_ids", "n_matrix_hits")
ms <- merge(m, sites[, ..site_cols], by = "tss_id", allow.cartesian = TRUE)
ms[, molecule_site_id := paste(observation_id, site_id, sep = "::")]
ms[, `:=`(evaluable = usable_sequence & footprint_evaluable,
          supported = FALSE, max_overlap_fraction = 0, supporting_footprint_ids = "")]
# Use observation IDs as independent sequence names: overlaps cannot cross fibers.
# Each observation encodes sample, original read and assigned TSS, checked at export.
hit <- findOverlaps(pol2_gr(ms$observation_id, ms$start, ms$end),
                    pol2_gr(fp$observation_id, fp$start, fp$end), ignore.strand = TRUE)
ov <- cbind(ms[queryHits(hit), .(molecule_site_id, observation_id, sample_name, original_RID,
             tss_id, cluster, category, site_id, motif_start = start, motif_end = end, evaluable)],
            fp[subjectHits(hit), .(footprint_id, footprint_start = start, footprint_end = end,
                                  footprint_length = length, block_score, block_edge_left, block_edge_right,
                                  read_score, quality_status, inside_same_read_fire,
                                  passes_notebook_footprint_filter)])
ov[, overlap_bp := pmin(motif_end, footprint_end) - pmax(motif_start, footprint_start)]
ov[, overlap_fraction := overlap_bp / (motif_end - motif_start)]
ov[, supported := overlap_fraction >= .5]
support <- ov[, .(supported = any(supported),
                 max_overlap_fraction = max(overlap_fraction),
                 supporting_footprint_ids = paste(unique(footprint_id[supported]), collapse = ";")), by = molecule_site_id]
ms[support, on = "molecule_site_id", `:=`(supported = i.supported,
   max_overlap_fraction = i.max_overlap_fraction, supporting_footprint_ids = i.supporting_footprint_ids)]
ms[evaluable == FALSE, `:=`(supported = NA, max_overlap_fraction = NA_real_)]
save_table(ov, "molecule_motif_footprint_overlaps.tsv.gz")
save_table(ms, "molecule_motif_sites.tsv.gz")

# Complete molecule x category universe, including absent motifs and unavailable U1.
status <- m[rep(seq_len(.N), each = length(categories))]
status[, category := rep(categories, times = nrow(m))]
status[, `:=`(category_available = category %in% available, n_sites = 0L,
              n_supported_sites = 0L)]
totals <- ms[, .(n_sites = .N, n_supported_sites = sum(supported, na.rm = TRUE)), by = .(observation_id, category)]
status[totals, on = .(observation_id, category), `:=`(n_sites = i.n_sites,
      n_supported_sites = i.n_supported_sites)]
status[, occurrence_evaluable := usable_sequence & category_available]
status[, motif_present := n_sites > 0]
status[occurrence_evaluable == FALSE, `:=`(motif_present = NA, n_sites = NA_integer_)]
status[, footprint_occurrence_evaluable := occurrence_evaluable & footprint_evaluable]
status[, support_evaluable := footprint_occurrence_evaluable & !is.na(motif_present) & motif_present]
status[, supported := n_supported_sites > 0]
# Motif-absent, evaluable molecules contribute zero to footprint occurrence.
status[footprint_occurrence_evaluable == FALSE, `:=`(supported = NA, n_supported_sites = NA_integer_)]
status[, status := fcase(!category_available, "unavailable_category", !usable_sequence, "unusable_sequence",
                        !footprint_evaluable, "unresolved_footprint_input", !motif_present, "motif_absent",
                        supported, "footprint_supported_ge50pct",
                        default = "no_qualifying_footprint_detected")]
save_table(status, "molecule_motif_status.tsv.gz")
save_table(status[status %in% c("unavailable_category", "unusable_sequence", "unresolved_footprint_input")],
           "unresolved_observations.tsv.gz")

summary <- status[, .(category_available = all(category_available),
  n_original_molecules = .N, n_promoters = uniqueN(sequence_id), n_tss_records = uniqueN(tss_id), n_samples = uniqueN(sample_name),
  footprint_occurrence_denominator = if (all(category_available)) .N else 0L,
  footprint_occurrence_numerator = sum(supported, na.rm = TRUE),
  n_unresolved_sequence = sum(!usable_sequence), n_unresolved_footprint = sum(!footprint_evaluable),
  n_motif_positive_promoters = uniqueN(sequence_id[!is.na(motif_present) & motif_present]),
  n_motif_positive_samples = uniqueN(sample_name[!is.na(motif_present) & motif_present])), by = .(cluster, category)]
ratio <- function(n, d) ifelse(d > 0, n/d, NA_real_)
summary[, footprint_occurrence_fraction := ratio(footprint_occurrence_numerator, footprint_occurrence_denominator)]
save_table(summary, "cluster_motif_summary.tsv")

bars <- summary[, .(cluster, category, numerator = footprint_occurrence_numerator,
                    denominator = footprint_occurrence_denominator, fraction = footprint_occurrence_fraction)]
bars[, label := ifelse(is.na(fraction), "NA", sprintf("%.1f%%", 100*fraction))]
bars[, cluster := factor(cluster, levels = clusters)]
save_table(bars, "cluster_footprint_proportions.tsv")
for (cat in available) {
  d <- bars[category == cat]
  p <- plot_stacked_proportion(d, x = "cluster", fill = "cluster", colors = colors,
        y = "fraction", position = "stack", labels = list(x = "Cluster", y = "Proportion of molecules")) +
    geom_text(aes(label = label), vjust = -.4, size = 3) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0,.15))) +
    labs(title = paste(cat, "motif-associated footprints"),
         subtitle = "At least one >=50%-supported motif / all original reads; 10-30 bp; no score or FIRE filter") +
    theme_fiberseq("classic", base_size = 10) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, paste0(cat, "_cluster_bars"), 7, 4.5)
}

capture.output(sessionInfo(), file = file.path(out, "R_session_info.txt"))
provenance <- jsonlite::read_json(file.path(out, "provenance.json"), simplifyVector = FALSE)
provenance$footprint_selection <- "Individual same-read FiberHMM TF blocks of 10-30 bp inclusive; no block-score or FIRE-containment filter. Full boundaries retained; blocks never merged."
provenance$overlap <- "One individual footprint covers >=50% of the full representative motif interval on the same read; motif length is the denominator. Each distinct read counts once per category."
provenance$promoter_balanced_occurrence <- NULL
provenance$summary <- list(n_original_reads = nrow(m),
  n_molecule_site_observations = nrow(ms), n_individual_overlaps = nrow(ov),
  n_footprint_evaluable_molecules = sum(m$footprint_evaluable),
  n_usable_sequence_molecules = sum(m$usable_sequence),
  n_length_eligible_footprints = nrow(fp), support_threshold = .5,
  cluster_bar_denominator = "All original reads in the cluster, including motif-absent, unsupported and sequence-unresolved reads. Sequence-unresolved reads remain unknown, not negative; their counts are recorded in cluster_motif_summary.tsv and figure labels. All 5000 have footprint input (asserted).",
  display = "Cluster bars only: one bar per original cluster, one PDF per available motif family.")
source_paths <- fread(file.path(out, "footprint_sources.tsv"))$path
provenance$footprint_sources <- pol2_stamp(source_paths[file.exists(source_paths)])
code_paths <- file.path(base,c("05_motif_scan.py","07_motif_cluster_bars.R"))
provenance$analysis_code <- data.frame(path = code_paths, md5 = unname(tools::md5sum(code_paths)))
shared_paths <- file.path(project, c("code/parsing_functions/plotting_functions.r",
                                   "code/parsing_functions/parsing_footprints_functions.r"))
provenance$shared_functions <- data.frame(path = shared_paths, md5 = unname(tools::md5sum(shared_paths)))
provenance$figure_format <- "Cluster-bar PDFs only; supported distinct reads divided by all original reads in each cluster."
jsonlite::write_json(provenance, file.path(out, "provenance.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
cat("Summarized", nrow(m), "original molecules,", uniqueN(m$tss_id), "TSSs,", uniqueN(m$sample_name),
    "samples;", nrow(ms), "molecule-site observations;", nrow(ov), "individual overlaps\n")

