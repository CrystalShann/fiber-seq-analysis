library(data.table)


accessibility_file <- "/project/spott/cshan/fiber-seq/results/PolII/m6a_pausing_quartiles/AL10_bc2178_19130_1kb_bin10_modthresh0.9_all/tables/AL10_bc2178_19130_tss_read_accessibility.tsv.gz"
footprint_dir <- "/project/spott/cshan/fiber-seq/results/PolII/footprint_summary_beds/FIRE_nucleosome"
out_dir <- "/project/spott/cshan/fiber-seq/results/PolII/m6a_pausing_quartiles/AL10_bc2178_19130_1kb_bin10_modthresh0.9_all/tss_accessible_FIRE_nucleosome_1kb"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

normalize_read_id <- function(x) {
  sub("^.*?(m[0-9]+_[0-9]+_[0-9]+_s[0-9]+/.+)$", "\\1", x, perl = TRUE)
}

message("Reading accessibility table: ", accessibility_file)
acc <- fread(
  accessibility_file,
  select = c(
    "read_id", "tss_uid", "gene_id", "gene_id_base", "gene_name",
    "chromosome", "tss_coordinate", "tss_strand", "PI",
    "pausing_group", "tss_accessible"
  )
)

acc[, read_core := normalize_read_id(read_id)]
acc[, tss_coordinate := as.integer(tss_coordinate)]
acc[, tss_accessible := as.logical(tss_accessible)]

accessible_tss_reads <- acc[
  tss_accessible == TRUE,
  .(
    read_id, read_core, tss_uid, gene_id, gene_id_base, gene_name,
    chrom = chromosome, tss_coordinate, tss_strand, PI, pausing_group,
    tss_window_start = tss_coordinate - 1000L,
    tss_window_end = tss_coordinate + 1000L
  )
]

rm(acc)
gc()

message("Accessible TSS-read rows: ", nrow(accessible_tss_reads))

chroms <- sort(unique(accessible_tss_reads$chrom))
footprint_counts_by_chrom <- vector("list", length(chroms))
names(footprint_counts_by_chrom) <- chroms
total_footprints_on_accessible_reads <- 0L
total_hits <- 0L

for (chrom_name in chroms) {
  chrom_accessible <- accessible_tss_reads[chrom == chrom_name]
  if (!nrow(chrom_accessible)) next

  bed_file <- file.path(footprint_dir, sprintf("%s_footprints.bed.gz", chrom_name))
  if (!file.exists(bed_file)) {
    message("  [WARN] Missing footprint BED for ", chrom_name, ": ", bed_file)
    next
  }

  read_id_file <- tempfile(pattern = paste0(chrom_name, "_accessible_reads_"), tmpdir = out_dir)
  fwrite(
    data.table(read_core = unique(chrom_accessible$read_core)),
    read_id_file,
    col.names = FALSE
  )
  on.exit(unlink(read_id_file), add = TRUE)

  message("Reading ", chrom_name, " footprints for ", uniqueN(chrom_accessible$read_core), " accessible reads")
  filter_cmd <- sprintf(
    "zcat %s 2>/dev/null | awk 'NR==FNR {ids[$1]=1; next} ($4 in ids) {print $1 \"\\t\" $4 \"\\t\" $5 \"\\t\" $6}' %s -",
    shQuote(bed_file),
    shQuote(read_id_file)
  )

  footprints <- suppressWarnings(fread(
    cmd = filter_cmd,
    col.names = c("chrom", "read_core", "footprint_class", "midpoint"),
    colClasses = c("character", "character", "character", "numeric"),
    showProgress = FALSE
  ))

  unlink(read_id_file)
  total_footprints_on_accessible_reads <- total_footprints_on_accessible_reads + nrow(footprints)
  if (!nrow(footprints)) next

  footprint_points <- footprints[
    !is.na(midpoint),
    .(
      read_core,
      chrom,
      footprint_start = midpoint,
      footprint_end = midpoint,
      footprint_class
    )
  ]
  rm(footprints)
  if (!nrow(footprint_points)) next

  tss_windows <- chrom_accessible[, .(
    read_core,
    chrom,
    tss_window_start,
    tss_window_end,
    tss_uid
  )]

  setkey(footprint_points, read_core, chrom, footprint_start, footprint_end)
  setkey(tss_windows, read_core, chrom, tss_window_start, tss_window_end)

  hits <- foverlaps(
    footprint_points,
    tss_windows,
    by.x = c("read_core", "chrom", "footprint_start", "footprint_end"),
    by.y = c("read_core", "chrom", "tss_window_start", "tss_window_end"),
    type = "within",
    nomatch = 0L
  )

  total_hits <- total_hits + nrow(hits)
  if (nrow(hits)) {
    footprint_counts_by_chrom[[chrom_name]] <- hits[, .(
      n_footprints_1kb = .N,
      n_FIRE_nucleosome_1kb = sum(footprint_class == "FIRE_nucleosome", na.rm = TRUE),
      footprint_classes_1kb = paste(sort(unique(footprint_class)), collapse = ",")
    ), by = .(read_core, tss_uid)]
  }

  rm(footprint_points, tss_windows, hits, chrom_accessible)
  gc()
}

footprint_counts <- rbindlist(footprint_counts_by_chrom, use.names = TRUE, fill = TRUE)
if (!nrow(footprint_counts)) {
  footprint_counts <- data.table(
    read_core = character(),
    tss_uid = character(),
    n_footprints_1kb = integer(),
    n_FIRE_nucleosome_1kb = integer(),
    footprint_classes_1kb = character()
  )
}

message("Footprints on accessible reads: ", total_footprints_on_accessible_reads)
message("Footprint hits within TSS +/- 1 kb: ", total_hits)
rm(footprint_counts_by_chrom)
gc()

read_tss_counts <- merge(
  accessible_tss_reads,
  footprint_counts,
  by = c("read_core", "tss_uid"),
  all.x = TRUE
)

read_tss_counts[is.na(n_footprints_1kb), n_footprints_1kb := 0L]
read_tss_counts[is.na(n_FIRE_nucleosome_1kb), n_FIRE_nucleosome_1kb := 0L]
read_tss_counts[is.na(footprint_classes_1kb), footprint_classes_1kb := "none"]

tss_summary <- read_tss_counts[, .(
  n_accessible_reads = .N,
  n_accessible_reads_with_footprint_1kb = sum(n_footprints_1kb > 0),
  total_footprints_1kb = sum(n_footprints_1kb),
  total_FIRE_nucleosomes_1kb = sum(n_FIRE_nucleosome_1kb),
  mean_footprints_per_accessible_read = mean(n_footprints_1kb),
  median_footprints_per_accessible_read = as.numeric(median(n_footprints_1kb))
), by = .(
  tss_uid, gene_id, gene_id_base, gene_name,
  chrom, tss_coordinate, tss_strand, PI, pausing_group
)]

pausing_summary <- read_tss_counts[, .(
  n_accessible_tss_read_rows = .N,
  n_accessible_tss_read_rows_with_footprint_1kb = sum(n_footprints_1kb > 0),
  total_footprints_1kb = sum(n_footprints_1kb),
  total_FIRE_nucleosomes_1kb = sum(n_FIRE_nucleosome_1kb),
  mean_footprints_per_accessible_read = mean(n_footprints_1kb),
  median_footprints_per_accessible_read = as.numeric(median(n_footprints_1kb))
), by = pausing_group]

setorder(read_tss_counts, pausing_group, chrom, tss_coordinate, read_core)
setorder(tss_summary, pausing_group, chrom, tss_coordinate)
setorder(pausing_summary, pausing_group)

fwrite(
  read_tss_counts,
  file.path(out_dir, "accessible_tss_read_FIRE_nucleosome_counts_1kb.tsv.gz"),
  sep = "\t"
)

fwrite(
  tss_summary,
  file.path(out_dir, "accessible_read_FIRE_nucleosome_counts_by_tss.tsv"),
  sep = "\t"
)

fwrite(
  pausing_summary,
  file.path(out_dir, "accessible_read_FIRE_nucleosome_counts_by_pausing_group.tsv"),
  sep = "\t"
)

message("Wrote tables to ", out_dir)
print(pausing_summary)
