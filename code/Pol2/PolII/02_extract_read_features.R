#!/usr/bin/env Rscript
# Per-read and per-gene features of the curated PolII example genes, extracted straight from
# the data: no read clustering and no dependence on recurrent_cluster_features outputs.
# All the work is done by pol2_pausing_functions.r and code/parsing_functions, sourced below
# and left unchanged.
#
# Run with the environment of run_gene_examples.sh, e.g.
#   sbatch --account=pi-spott --partition=caslake --cpus-per-task=8 --mem=48G --time=08:00:00 \
#     --output=/project/spott/cshan/fiber-seq/results/logs/extract_read_features_%j.log \
#     --wrap "bash code/Pol2/PolII/run_gene_examples.sh code/Pol2/PolII/02_extract_read_features.R"
# Worker count: READ_FEATURE_WORKERS, else SLURM_CPUS_PER_TASK, else 1.
#
# Outputs under LCL_project/Pol2/PolII/read_features/
#   tables/read_features.tsv.gz   one row per read of every gene
#   tables/gene_features.tsv      one row per input gene, failures kept with their error
#   tables/status/<gene>.tsv      that gene's row, written as soon as it finishes
#   bins/<gene>.tsv.gz            one row per read x 10 bp bin, raw counts
#
# Definitions follow 01_summarize_recurrent_patterns.R so the two agree read for read:
# window = TSS +/- FLANK oriented by transcription, full-span reads only, m6A rates are means
# over the sites observed in at least one read of that gene, nucleosomes are the fibertools
# 130-160 bp calls and TF footprints the FiberHMM combined 10-80 bp calls, split by
# assign_size_class(TF_BREAKS). The bin files hold counts, not rates, so any later grouping
# (cluster, family, gene, archetype) is a sum of numerators over a sum of denominators.
#
# Occlusion and the +1 nucleosome are per read: a read is occluded when a nucleosome covers the
# TSS base, and +1-enriched when a nucleosome *starts* within PLUS1_OFFSETS without covering the
# TSS, so the two are independent. plus1 is also kept as plus1_offset() defines it.

dyn.load("/project/spott/cshan/envs/dimelo/lib/libcurl.so.4", local = FALSE)
source("/project/spott/cshan/fiber-seq/code/Pol2/PolII/pol2_pausing_functions.r")
setDTthreads(1L)

ROOT    <- file.path(OUT_DIR, "read_features")
BIN     <- 10L                                   # bin width of the profile tables
NBINS   <- (2L * FLANK) %/% BIN                  # 200 bins; the last one also holds +FLANK
BIN_POS <- seq(EXAMPLE_RANGE[1] + BIN %/% 2L, EXAMPLE_RANGE[2] - BIN %/% 2L, BIN)
TF_CLASSES  <- levels(assign_size_class(integer(), TF_BREAKS))
TF_COLS     <- paste0("tf_", gsub("[- ]", "_", sub(" bp$", "", TF_CLASSES)))
COVER_COLS  <- c("nuc", TF_COLS)
ARCHETYPES  <- c("A1 open / +1-depleted",     "A2 open / +1-enriched",
                 "A3 occluded / +1-depleted", "A4 occluded / +1-enriched")
# regions of the per-read m6A rates, as in 01_summarize_recurrent_patterns.R
M6A_REGIONS <- list(total = c(-1000L, 1000L), core = c(-200L, 200L), upstream = c(-500L, -201L),
                    downstream = c(201L, 500L), far_up = c(-1000L, -501L), far_down = c(501L, 1000L))


############################################################
# nucleosome and TF footprint intervals of one gene
############################################################
# 1. read each sample's fibertools nucleosome BED12 over the window, keep the longest alignment
#    per read, expand to BED6 and keep the NUC_BREAKS (130-160 bp) calls
# 2. stream the FiberHMM combined 10-80 bp TF footprints of the window (window_tf_footprints())
#    and map their read names onto the sample-namespaced RIDs of the m6A cohort
# 3. put reads, nucleosomes and TF footprints on the TSS-relative axis with to_relative()
read_bundle <- function(g, dat) {
  region <- GRanges(g$chrom, IRanges(g$tss + EXAMPLE_RANGE[1], g$tss + EXAMPLE_RANGE[2]))
  meta <- as.data.table(dat$rids_df)
  if (anyDuplicated(meta$original_RID)) stop("ambiguous read names for ", g$gene_name)
  nucs <- rbindlist(lapply(unique(meta$sample_name), function(s) {
    d <- read_ft_bed12(nuc_file(s, g$chrom), region, longest_alignment = TRUE)
    if (!nrow(d)) return(NULL)
    b <- as.data.table(convert_ft_bed12_to_bed6(d, format = "bed12_fibertools", source = nuc_file(s, g$chrom)))
    b[, .(RID = paste(s, RID, sep = "::"), start0 = start, end0 = end)]
  }))
  if (!nrow(nucs)) nucs <- data.table(RID = character(), start0 = integer(), end0 = integer())
  nucs <- nucs[RID %in% meta$RID & !is.na(assign_size_class(end0 - start0, NUC_BREAKS))]
  tf <- if (g$chrom %in% AUTOSOMES) window_tf_footprints(g$chrom, start(region), end(region)) else
    data.table(RID = character(), start0 = integer(), end0 = integer())
  if (!nrow(tf)) tf <- data.table(RID = character(), start0 = integer(), end0 = integer())
  tf <- merge(tf, meta[, .(original_RID, qualified_RID = RID)], by.x = "RID", by.y = "original_RID",
              allow.cartesian = TRUE)[, .(RID = qualified_RID, start0, end0)]
  relative <- function(d) {
    r <- to_relative(d$start0, d$end0, g$tss, rep(g$strand, nrow(d)))
    d[, `:=`(lo = r$lo, hi = r$hi, size = end0 - start0)][]
  }
  list(reads = relative(meta[, .(RID, start0 = start - 1L, end0 = end)]),
       nucs = relative(nucs), fp = relative(tf))
}


############################################################
# per read features and bin counts of one gene
############################################################
# 1. open the TSS +/- EXAMPLE_RANGE window and assemble its full-span m6A reads
# 2. orient the m6A matrix by transcription and mark the sites observed in at least one read
# 3. read_bundle() for the nucleosome and TF footprint intervals of the same reads
# 4. footprint_matrix() per cover class, reduced to a logical per-base coverage matrix
# 5. per read: regional m6A rates over observed sites, mean coverage per class over the whole
#    window and the core, plus1_offset(), the +1 nucleosome start and the occlusion archetype
# 6. per read x 10 bp bin: m6A calls, observed sites, and covered bases per class
# 7. return the read rows, the bin rows and the gene's averages

extract_gene <- function(g, samples) {
  message("Extracting ", g$gene_name)
  region <- data.frame(region_id = g$gene_name, chr = g$chrom,
    start = g$tss + EXAMPLE_RANGE[1] - 1L, end = g$tss + EXAMPLE_RANGE[2],
    analysis_start = g$tss + EXAMPLE_RANGE[1], analysis_end = g$tss + EXAMPLE_RANGE[2])
  dat <- assemble_region_m6a(as.data.frame(samples), region, "m6a", full_span = TRUE,
    positions = seq(region$analysis_start, region$analysis_end))
  bundle <- read_bundle(g, dat)
  reads <- bundle$reads
  mat <- as.matrix(dat$met_mat[reads$RID, , drop = FALSE])
  offset <- (as.integer(colnames(mat)) - g$tss) * if (g$strand == "+") 1L else -1L
  ord <- order(offset); offset <- offset[ord]; mat <- mat[, ord, drop = FALSE]
  stopifnot(identical(offset, EXAMPLE_RANGE[1]:EXAMPLE_RANGE[2]))
  observed <- colSums(mat) > 0

  cover <- c(list(nuc = bundle$nucs),
             setNames(lapply(TF_CLASSES, function(k) bundle$fp[assign_size_class(size, TF_BREAKS) == k]), TF_COLS))
  nuc_m <- footprint_matrix(reads, bundle$nucs)
  cover_m <- lapply(cover, function(d) footprint_matrix(reads, d) > 0)

  a <- data.table(gene = g$gene_name, group = g$group, chrom = g$chrom, tss = g$tss,
                  strand = g$strand, RID = reads$RID,
                  sample_name = dat$rids_df$sample_name[match(reads$RID, dat$rids_df$RID)],
                  start = reads$start0 + 1L, end = reads$end0)
  for (nm in names(M6A_REGIONS)) {
    ix <- offset >= M6A_REGIONS[[nm]][1] & offset <= M6A_REGIONS[[nm]][2] & observed
    set(a, j = paste0("m6a_", nm),
        value = if (any(ix)) rowMeans(mat[, ix, drop = FALSE]) else NA_real_)
  }
  for (nm in names(cover_m)) {
    set(a, j = paste0(nm, "_total"), value = rowMeans(cover_m[[nm]]))
    set(a, j = paste0(nm, "_core"), value = rowMeans(cover_m[[nm]][, abs(offset) <= 200, drop = FALSE]))
  }
  a[, nuc_at_tss := cover_m$nuc[, offset == 0L]]
  a[, plus1 := plus1_offset(nuc_m)]
  # +1 nucleosome start, ignoring a nucleosome that already covers the TSS
  starts <- bundle$nucs[lo >= min(PLUS1_OFFSETS) & lo <= max(PLUS1_OFFSETS) & !(lo <= 0L & hi >= 0L),
                        .(plus1_start = min(lo)), by = RID]
  a[, plus1_start := starts$plus1_start[match(RID, starts$RID)]]
  a[, archetype := factor(ARCHETYPES[1L + 2L * nuc_at_tss + !is.na(plus1_start)], levels = ARCHETYPES)]

  bin <- pmin((offset - EXAMPLE_RANGE[1]) %/% BIN, NBINS - 1L)
  by_bin <- function(m) t(rowsum(t(m * 1L), bin))
  bins <- data.table(gene = g$gene_name, RID = rep(reads$RID, NBINS),
                     bin = rep(0:(NBINS - 1L), each = nrow(reads)),
                     position = rep(BIN_POS, each = nrow(reads)),
                     m6a_sites = rep(as.vector(rowsum(as.integer(observed), bin)), each = nrow(reads)),
                     bin_bases = rep(as.vector(rowsum(rep(1L, length(bin)), bin)), each = nrow(reads)),
                     m6a_calls = as.vector(by_bin(mat)))
  for (nm in names(cover_m)) set(bins, j = paste0(nm, "_bases"), value = as.vector(by_bin(cover_m[[nm]])))

  metrics <- setdiff(names(a), c("gene", "group", "chrom", "tss", "strand", "RID", "sample_name",
                                 "start", "end", "archetype"))
  gene_row <- a[, c(list(n_reads = .N, n_samples = uniqueN(sample_name),
                         occluded_fraction = mean(nuc_at_tss),
                         plus1_fraction = mean(!is.na(plus1_start))),
                    lapply(.SD, mean, na.rm = TRUE)), .SDcols = metrics]
  for (k in ARCHETYPES) set(gene_row, j = paste0("fraction_", sub(" .*", "", k)),
                            value = mean(a$archetype == k))
  list(reads = a, bins = bins, gene = gene_row)
}


############################################################
# main
############################################################
# 1. take the curated gene list and the 31 sample metatable
# 2. per gene: resolve its TSS (FANTOM5 dominant, else GENCODE canonical), run extract_gene(),
#    write bins/<gene>.tsv.gz and status/<gene>.tsv, and record the error if it fails
# 3. run the genes in parallel, then write the pooled read table and the gene table
# 4. raise the collected errors once every gene has been attempted
main <- function() {
  genes <- curated_example_genes()
  samples <- fread(SAMPLES)[, .(sample_name, fire_dir)]
  dominant <- fread(DOMINANT_TSS)
  for (d in c("tables/status", "bins")) dir.create(file.path(ROOT, d), recursive = TRUE, showWarnings = FALSE)
  workers <- as.integer(Sys.getenv("READ_FEATURE_WORKERS", Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
  stopifnot(nrow(genes) > 0L, nrow(samples) == 31L, !is.na(workers), workers >= 1L)

  run_one <- function(i) {
    row <- copy(as.data.table(genes)[i])
    row[, `:=`(status = "pending", error = "", chrom = NA_character_, tss = NA_integer_,
               strand = NA_character_, tss_source = NA_character_)]
    out <- NULL
    tryCatch({
      g <- resolve_example_tss(row$gene_name, row$group, dominant)
      row[, `:=`(chrom = g$chrom, tss = g$tss, strand = g$strand, tss_source = g$tss_source)]
      out <- extract_gene(g, samples)
      fwrite(out$bins, file.path(ROOT, "bins", paste0(g$gene_name, ".tsv.gz")), sep = "\t", na = "NA")
      row <- cbind(row[, status := "complete"], out$gene)
    }, error = function(e) {
      row[, `:=`(status = "error", error = conditionMessage(e))]
      message("ERROR ", row$gene_name, ": ", conditionMessage(e))
    })
    fwrite(row, file.path(ROOT, "tables/status", paste0(row$gene_name, ".tsv")), sep = "\t", na = "NA")
    list(reads = out$reads, gene = row)
  }
  result <- parallel::mclapply(seq_len(nrow(genes)), run_one, mc.cores = workers, mc.preschedule = FALSE)
  if (any(vapply(result, inherits, logical(1), "try-error")))
    stop("a gene worker failed; see read_features/tables/status")

  gene_features <- rbindlist(lapply(result, `[[`, "gene"), fill = TRUE)
  reads <- rbindlist(lapply(result, `[[`, "reads"))
  stopifnot(nrow(gene_features) == nrow(genes), !anyDuplicated(reads[, .(gene, RID)]))
  fwrite(reads, file.path(ROOT, "tables/read_features.tsv.gz"), sep = "\t", na = "NA")
  fwrite(gene_features, file.path(ROOT, "tables/gene_features.tsv"), sep = "\t", na = "NA")
  message("Complete: ", sum(gene_features$status == "complete"), "/", nrow(genes), " genes; ",
          nrow(reads), " reads.")
  # Too few full-span reads is a limit of the input, not a failure: the gene stays in the gene
  # table with its error and the run still finishes.
  if (any(gene_features$status == "error"))
    message("Skipped:\n", paste(gene_features[status == "error", paste(gene_name, error, sep = ": ")],
                                collapse = "\n"))
}
if (sys.nframe() == 0L) main()
