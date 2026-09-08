# LCL windows, m6A matrices, pooled Leiden analysis, haplotype annotation and result tables.
# The shared Manhattan/Leiden algorithm remains in the non-LCL helper.
source("/project/spott/cshan/fiber-seq/code/clustering_TSS_methods/Leiden_Manhattan/leiden_manhattan_functions.r")

# Binary calls for full-span LCL reads; zero-call molecules remain in the matrix.
lcl_sparse_m6a_matrix <- function(reads, rids_df, window_start, window_end) {
  stopifnot(!anyDuplicated(rids_df$RID), all(rids_df$start <= window_start),
            all(rids_df$end >= window_end))
  sites <- sort(unique(reads$pos[reads$pos >= window_start & reads$pos <= window_end]))
  pairs <- unique(data.frame(i = match(reads$RID, rids_df$RID), j = match(reads$pos, sites)))
  pairs <- pairs[!is.na(pairs$i) & !is.na(pairs$j), , drop = FALSE]
  Matrix::sparseMatrix(i = pairs$i, j = pairs$j, x = rep(1, nrow(pairs)),
    dims = c(nrow(rids_df), length(sites)), dimnames = list(as.character(rids_df$RID), as.character(sites)))
}

construct_lcl_windows <- function(coordinate_regions, promoter_genes, tss_bed,
                                  tss_flank = 1000L) {
  stopifnot(length(tss_flank) == 1L, tss_flank > 0, tss_flank == as.integer(tss_flank))
  required <- c("region_id", "chr", "start", "end", "annotation")
  stopifnot(all(required %in% names(coordinate_regions)))
  tss <- data.table::fread(tss_bed, header = FALSE, data.table = FALSE)
  stopifnot(ncol(tss) == 6L, all(tss[[3]] - tss[[2]] == 20L))
  names(tss) <- c("chr", "bed_start", "bed_end", "annotation", "score", "strand")
  fields <- strsplit(tss$annotation, ";", fixed = TRUE)
  tss$gene <- vapply(fields, `[`, character(1), 3L)
  tss$ensg <- vapply(fields, `[`, character(1), 1L)
  tss$tss0 <- tss$bed_start + 10L
  promoters <- lapply(promoter_genes, function(gene) {
    selected <- tss[tss$gene == gene, ]
    if (nrow(selected) != 1L) stop("Expected one canonical TSS for ", gene)
    data.frame(region_id = paste0(gene, "_promoter"), chr = selected$chr,
               start = selected$tss0 - tss_flank, end = selected$tss0 + tss_flank,
               annotation = paste(gene, "promoter"), gene = gene,
               ensg = selected$ensg, strand = selected$strand,
               tss = selected$tss0 + 1L, region_type = "promoter")
  })
  coordinates <- as.data.frame(coordinate_regions)
  coordinates$gene <- NA_character_
  coordinates$ensg <- NA_character_
  coordinates$strand <- "*"
  coordinates$tss <- NA_integer_
  coordinates$region_type <- ifelse(grepl("dsQTL", coordinates$annotation), "dsQTL", "region")
  regions <- dplyr::bind_rows(coordinates, dplyr::bind_rows(promoters))
  stopifnot(nrow(regions) > 0L, !anyNA(regions$start), !anyNA(regions$end),
            all(regions$start >= 0), all(regions$end > regions$start),
            all(regions$start == as.integer(regions$start)),
            all(regions$end == as.integer(regions$end)),
            !anyDuplicated(regions$region_id),
            all(grepl("^[A-Za-z0-9_-]+$", regions$region_id)))
  regions$width <- regions$end - regions$start
  regions$analysis_start <- regions$start + 1L
  regions$analysis_end <- regions$end
  regions$coordinate_system <- "BED: 0-based, half-open"
  regions
}

lcl_source_signature <- function(paths) {
  info <- file.info(paths)
  if (anyNA(info$size)) stop("Missing input: ", paste(paths[is.na(info$size)], collapse = ", "))
  data.frame(path = normalizePath(paths), size = info$size, mtime = as.numeric(info$mtime))
}

lcl_extracted_path <- function(sample_table, chromosome, feature) {
  file.path(sample_table$fire_dir, "extracted_results", paste0(feature, "_by_chr"),
            paste0(sample_table$sample_name, ".ft_extracted_", feature, ".", chromosome, ".bed.gz"))
}

read_lcl_bed12 <- function(path, region) {
  if (!file.exists(paste0(path, ".tbi"))) stop("Missing tabix index for ", path)
  query <- GenomicRanges::GRanges(region$chr,
    IRanges::IRanges(region$analysis_start, region$analysis_end))
  records <- Rsamtools::scanTabix(Rsamtools::TabixFile(path), param = query)[[1]]
  if (!length(records)) return(data.frame())
  bed <- data.table::fread(text = paste(records, collapse = "\n"), header = FALSE,
                          data.table = FALSE)
  if (ncol(bed) != 12L) stop("Expected fibertools BED12 in ", path)
  names(bed) <- c("chr", "start", "end", "original_RID", "score", "strand",
                  "thick_start", "thick_end", "rgb", "blockCount", "blockSizes", "blockStarts")
  bed <- bed[order(bed$original_RID, -(bed$end - bed$start)), ]
  bed[!duplicated(bed$original_RID), , drop = FALSE]
}

lcl_bed12_blocks <- function(bed) {
  empty <- data.frame(RID = character(), original_RID = character(),
                       start = integer(), end = integer(), size = integer())
  if (!nrow(bed)) return(empty)
  blocks <- lapply(seq_len(nrow(bed)), function(row_index) {
    record <- bed[row_index, ]
    sizes <- as.integer(strsplit(sub(",$", "", record$blockSizes), ",", fixed = TRUE)[[1]])
    offsets <- as.integer(strsplit(sub(",$", "", record$blockStarts), ",", fixed = TRUE)[[1]])
    if (length(sizes) != record$blockCount || length(offsets) != record$blockCount ||
        anyNA(sizes) || anyNA(offsets)) stop("Invalid BED12 blocks: ", record$original_RID)
    if (length(sizes) <= 2L) return(empty)
    selected <- seq.int(2L, length(sizes) - 1L)
    data.frame(RID = record$RID, original_RID = record$original_RID,
               start = record$start + offsets[selected],
               end = record$start + offsets[selected] + sizes[selected], size = sizes[selected])
  })
  dplyr::bind_rows(blocks)
}

assemble_lcl_region_m6a <- function(sample_table, region, matrix_dir, reuse = TRUE) {
  region_matrix_dir <- file.path(matrix_dir, region$region_id)
  dir.create(region_matrix_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- lcl_extracted_path(sample_table, region$chr, "m6a")
  signature <- list(version = 1L, region = region, samples = sample_table,
                    sources = lcl_source_signature(c(paths, paste0(paths, ".tbi"))))
  cache <- file.path(region_matrix_dir, paste0(region$region_id, "_m6a_matrix.rds"))
  if (reuse && file.exists(cache)) {
    previous <- readRDS(cache)
    if (identical(previous$signature, signature)) return(previous)
  }
  reads_list <- metadata_list <- qc_list <- vector("list", nrow(sample_table))
  for (sample_index in seq_len(nrow(sample_table))) {
    sample_name <- sample_table$sample_name[sample_index]
    bed <- read_lcl_bed12(paths[sample_index], region)
    overlapping <- nrow(bed)
    if (overlapping) bed <- bed[bed$start <= region$start & bed$end >= region$end, , drop = FALSE]
    qc_list[[sample_index]] <- data.frame(sample_name = sample_name,
      overlapping_reads = overlapping, full_span_reads = nrow(bed))
    if (!nrow(bed)) next
    bed$RID <- paste(sample_name, bed$original_RID, sep = "::")
    metadata_list[[sample_index]] <- data.frame(
      RID = bed$RID, original_RID = bed$original_RID, chr = bed$chr,
      start = bed$start + 1L, end = bed$end, strand = bed$strand,
      sample_name = sample_name, score = bed$score)
    blocks <- lcl_bed12_blocks(bed)
    if (any(blocks$size != 1L)) stop("Non-single-base m6A block in ", paths[sample_index])
    blocks <- blocks[blocks$end >= region$analysis_start & blocks$end <= region$analysis_end, ]
    reads_list[[sample_index]] <- data.frame(RID = blocks$RID, pos = blocks$end)
  }
  rids_df <- dplyr::bind_rows(metadata_list)
  reads <- dplyr::bind_rows(reads_list)
  if (nrow(rids_df) < 3L || !nrow(reads)) stop("Insufficient full-span reads/m6A sites: ", region$region_id)
  met_mat <- lcl_sparse_m6a_matrix(reads, rids_df, region$analysis_start, region$analysis_end)
  stopifnot(!anyNA(met_mat), !anyDuplicated(rownames(met_mat)))
  result <- list(met_mat = met_mat, rids_df = rids_df, qc = dplyr::bind_rows(qc_list),
                  signature = signature)
  saveRDS(result, cache)
  data.table::fwrite(result$qc, file.path(region_matrix_dir, paste0(region$region_id, "_read_qc.tsv")), sep = "\t")
  result
}


# HiPhase and VCF annotation; haplotypes never define clustering groups.
lcl_phasing_paths <- function(sample_name, phasing_root) {
  directory <- file.path(phasing_root, sample_name)
  c(vcf = file.path(directory, paste0(sample_name, ".5mC.6mA.aligned.phased.vcf.gz")),
    tags = file.path(directory, "read-level-phasing.tsv"),
    blocks = file.path(directory, "blocks.tsv"))
}

lcl_vcf_sample_column <- function(path, cell_line) {
  connection <- gzfile(path, "rt")
  on.exit(close(connection))
  repeat {
    header <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(header)) stop("VCF has no #CHROM header: ", path)
    if (startsWith(header, "#CHROM\t")) break
  }
  columns <- strsplit(header, "\t", fixed = TRUE)[[1]]
  selected <- match(cell_line, columns)
  if (is.na(selected) || selected < 10L) return(NA_integer_)
  selected
}

cache_lcl_haplotags <- function(sample_table, region_data, phasing_root, cache_dir,
                                reuse = TRUE) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  read_info <- dplyr::bind_rows(lapply(region_data, function(dataset) dataset$rids_df))
  results <- setNames(vector("list", nrow(sample_table)), sample_table$sample_name)
  for (sample_index in seq_len(nrow(sample_table))) {
    sample <- sample_table[sample_index, , drop = FALSE]
    paths <- lcl_phasing_paths(sample$sample_name, phasing_root)
    required <- c(unname(paths), paste0(paths[["vcf"]], ".tbi"))
    original_ids <- sort(unique(read_info$original_RID[read_info$sample_name == sample$sample_name]))
    signature <- list(version = 1L, cell_line = sample$cell_line, reads = original_ids,
      required = required, sources = lcl_source_signature(required[file.exists(required)]))
    cache <- file.path(cache_dir, paste0(sample$sample_name, "_haplotags.rds"))
    if (reuse && file.exists(cache)) {
      previous <- tryCatch(readRDS(cache), error = function(e) {
        message("Rebuilding unreadable haplotag cache: ", cache)
        NULL
      })
      if (!is.null(previous) && identical(previous$signature, signature)) {
        results[[sample$sample_name]] <- previous
        next
      }
    }
    message("Reading haplotags for ", sample$sample_name)
    result <- list(available = FALSE, reason = "missing_phasing_files", paths = paths,
      missing_files = required[!file.exists(required)], signature = signature)
    if (all(file.exists(required))) {
      result$sample_column <- lcl_vcf_sample_column(paths[["vcf"]], sample$cell_line)
      result$reason <- "cell_line_absent_from_vcf"
      if (!is.na(result$sample_column)) {
        blocks <- data.table::fread(paths[["blocks"]], data.table = FALSE,
                                    colClasses = "character")
        stopifnot(all(c("sample_name", "chrom", "phase_block_id", "start", "end") %in% names(blocks)))
        blocks <- blocks[blocks$sample_name == sample$cell_line, , drop = FALSE]
        blocks$start <- as.integer(blocks$start)
        blocks$end <- as.integer(blocks$end)
        stopifnot(!anyNA(blocks$start), !anyNA(blocks$end))
        tags <- data.frame(sample_name = character(), chrom = character(),
          phase_block_id = character(), read_name = character(), haplotag = integer())
        if (length(original_ids)) {
          ids_file <- tempfile(fileext = ".txt")
          writeLines(original_ids, ids_file)
          program <- "BEGIN { FS = OFS = \"\\t\" } NR == FNR { wanted[$1] = 1; next } FNR == 1 || ($5 in wanted) { print }"
          tags <- tryCatch(data.table::fread(cmd = paste("awk", shQuote(program),
            shQuote(ids_file), shQuote(paths[["tags"]])), colClasses = "character",
            data.table = FALSE), finally = unlink(ids_file))
          stopifnot(all(c("sample_name", "chrom", "phase_block_id", "read_name", "haplotag") %in% names(tags)))
          tags <- tags[tags$sample_name == sample$cell_line, , drop = FALSE]
          tags$haplotag <- as.integer(tags$haplotag)
          stopifnot(!anyNA(tags$haplotag), all(tags$haplotag %in% c(1L, 2L)))
        }
        result$available <- TRUE
        result$reason <- "available"
        result$blocks <- blocks
        result$tags <- unique(tags)
      }
    }
    temporary_cache <- tempfile(pattern = "haplotags_", tmpdir = cache_dir)
    saveRDS(result, temporary_cache)
    if (!file.rename(temporary_cache, cache)) stop("Could not replace haplotag cache: ", cache)
    results[[sample$sample_name]] <- result
  }
  results
}

lcl_phase_variants <- function(phase, region, phase_set, query_start, query_end) {
  empty <- data.frame(chr = character(), pos = integer(), variant_id = character(),
    ref = character(), alt = character(), genotype = character(), phase_set = character(),
    HP1_allele = character(), HP2_allele = character())
  query <- GenomicRanges::GRanges(region$chr, IRanges::IRanges(query_start, query_end))
  raw <- Rsamtools::scanTabix(Rsamtools::TabixFile(phase$paths[["vcf"]]), param = query)[[1]]
  variants <- lapply(raw, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < phase$sample_column) stop("Truncated VCF record")
    keys <- strsplit(fields[9], ":", fixed = TRUE)[[1]]
    values <- strsplit(fields[phase$sample_column], ":", fixed = TRUE)[[1]]
    genotype <- values[match("GT", keys)]
    block <- values[match("PS", keys)]
    if (is.na(genotype) || is.na(block) || block != phase_set || !grepl("|", genotype, fixed = TRUE)) return(NULL)
    allele_indices <- suppressWarnings(as.integer(strsplit(genotype, "|", fixed = TRUE)[[1]]))
    alleles <- c(fields[4], strsplit(fields[5], ",", fixed = TRUE)[[1]])
    if (length(allele_indices) != 2L || anyNA(allele_indices) ||
        any(allele_indices < 0L | allele_indices >= length(alleles)) ||
        allele_indices[1] == allele_indices[2]) return(NULL)
    data.frame(chr = fields[1], pos = as.integer(fields[2]), variant_id = fields[3],
      ref = fields[4], alt = fields[5], genotype = genotype, phase_set = block,
      HP1_allele = alleles[allele_indices[1] + 1L], HP2_allele = alleles[allele_indices[2] + 1L])
  })
  dplyr::bind_rows(empty, dplyr::bind_rows(variants))
}

annotate_lcl_haplotypes <- function(dataset, region, sample_table, phase_cache,
                                    variant_flank = 50000L) {
  stopifnot(variant_flank >= 0L)
  reads <- dataset$rids_df
  reads$haplotype <- "unphased"
  reads$phase_set <- NA_character_
  reads$phase_status <- "unassigned_read"
  reads$phase_source <- "HiPhase read-level-phasing.tsv + sample VCF"
  variants <- list()
  anchor <- if (is.na(region$tss)) floor(mean(c(region$analysis_start, region$analysis_end))) else region$tss
  for (sample_name in sample_table$sample_name) {
    read_indices <- which(reads$sample_name == sample_name)
    if (!length(read_indices)) next
    phase <- phase_cache[[sample_name]]
    if (!phase$available) {
      reads$phase_status[read_indices] <- phase$reason
      next
    }
    blocks <- phase$blocks
    blocks <- blocks[blocks$chrom == region$chr & blocks$start <= anchor & blocks$end >= anchor, ]
    blocks <- unique(blocks[, c("phase_block_id", "start", "end"), drop = FALSE])
    if (nrow(blocks) != 1L) {
      reads$phase_status[read_indices] <- if (nrow(blocks)) "ambiguous_anchor_phase_block" else "no_anchor_phase_block"
      next
    }
    phase_set <- blocks$phase_block_id
    local_variants <- lcl_phase_variants(phase, region, phase_set,
      max(1L, blocks$start, region$analysis_start - variant_flank),
      min(blocks$end, region$analysis_end + variant_flank))
    if (!nrow(local_variants)) {
      reads$phase_status[read_indices] <- "no_local_phased_heterozygous_variant"
      next
    }
    local_variants$sample_name <- sample_name
    local_variants$cell_line <- sample_table$cell_line[match(sample_name, sample_table$sample_name)]
    local_variants$region_id <- region$region_id
    variants[[sample_name]] <- local_variants
    tags <- phase$tags
    tags <- tags[tags$chrom == region$chr & tags$phase_block_id == phase_set, , drop = FALSE]
    reads$phase_set[read_indices] <- phase_set
    for (read_index in read_indices) {
      haplotags <- unique(tags$haplotag[tags$read_name == reads$original_RID[read_index]])
      if (length(haplotags) == 1L) {
        reads$haplotype[read_index] <- paste0("HP", haplotags)
        reads$phase_status[read_index] <- "assigned"
      } else if (length(haplotags) > 1L) {
        reads$phase_status[read_index] <- "conflicting_haplotags"
      }
    }
  }
  list(reads = reads, variants = dplyr::bind_rows(variants))
}

# Collapse sample-level variants to exact in-window SNP positions and REF/ALT bases.
lcl_window_snps <- function(variants, region) {
  empty <- data.frame(pos = integer(), label = character())
  if (is.null(variants) || !nrow(variants)) return(empty)
  snps <- variants[variants$chr == region$chr & variants$pos >= region$analysis_start &
    variants$pos <= region$analysis_end & nchar(variants$ref) == 1L &
    grepl("^[ACGT](,[ACGT])*$", variants$alt), , drop = FALSE]
  if (!nrow(snps)) return(empty)
  snps$id <- ifelse(is.na(snps$variant_id) | snps$variant_id == ".", "SNP", snps$variant_id)
  snps$label <- paste0(snps$id, " ", snps$pos, " ", snps$ref, ">", snps$alt)
  labels <- tapply(snps$label, snps$pos, function(x) paste(unique(x), collapse = "; "))
  data.frame(pos = as.integer(names(labels)), label = unname(labels))
}

# Persist the pooled m6A analysis, independently of any footprint extraction.
run_lcl_clustering <- function(regions, sample_table, matrix_dir, output_dir,
                                window_size = 0L, k_neighbors = 10L,
                                leiden_resolution = 1, leiden_seed = 1L,
                                kernel_sigma = NULL, reuse_cache = TRUE) {
  dir.create(file.path(output_dir, "summary tables"), recursive = TRUE, showWarnings = FALSE)
  run_summary <- list()
  result_paths <- setNames(character(nrow(regions)), regions$region_id)
  for (region_index in seq_len(nrow(regions))) {
    region <- regions[region_index, , drop = FALSE]
    message("Clustering ", region$region_id)
    dat <- assemble_lcl_region_m6a(sample_table, region, matrix_dir, reuse_cache)
    result <- leiden_manhattan_cluster(
      dat$met_mat, dat$rids_df, region$analysis_start, region$analysis_end,
      window_size = window_size, k_neighbors = k_neighbors,
      resolution = leiden_resolution, sigma = kernel_sigma, seed = leiden_seed,
      impute_missing = FALSE)
    result$site_met_mat <- dat$met_mat
    result$region <- region
    result$assignments$original_RID <- dat$rids_df$original_RID[
      match(result$assignments$RID, dat$rids_df$RID)]
    result$assignments$sample_label <- sub("_.*$", "", result$assignments$sample_name)
    result$assignments$region_id <- region$region_id
    result$assignments$annotation <- region$annotation
    result$params$region_chr <- region$chr
    result$params$coordinate_system <- "1-based inclusive"
    region_dir <- file.path(output_dir, region$region_id, paste0("bin", window_size),
                            paste0("k", k_neighbors), paste0("resolution", leiden_resolution), "tables")
    dir.create(region_dir, recursive = TRUE, showWarnings = FALSE)
    result_paths[region$region_id] <- file.path(region_dir, "clustering.rds")
    saveRDS(result, result_paths[region$region_id])
    data.table::fwrite(result$assignments, file.path(region_dir, "read_cluster_assignments.tsv"), sep = "\t")
    data.table::fwrite(result$window_anno, file.path(region_dir, "window_annotation.tsv"), sep = "\t")
    data.table::fwrite(cluster_site_profiles(result, dat$met_mat),
                       file.path(region_dir, "aggregate_accessibility.tsv"), sep = "\t")
    counts <- result$assignments %>% count(cluster, sample_name, sample_label)
    data.table::fwrite(counts, file.path(region_dir, "cluster_sample_counts.tsv"), sep = "\t")
    data.table::fwrite(data.frame(parameter = names(result$params),
      value = vapply(result$params, as.character, character(1))),
      file.path(region_dir, "run_parameters.tsv"), sep = "\t")
    run_summary[[region_index]] <- data.frame(region, n_reads = nrow(dat$met_mat),
      n_zero_call_reads = sum(Matrix::rowSums(dat$met_mat) == 0),
      n_features = ncol(result$feat_mat), n_clusters = result$n_clusters,
      bin = window_size, k_neighbors = k_neighbors, k_eff = result$params$k_eff,
      resolution = leiden_resolution, seed = leiden_seed)
    data.table::fwrite(dplyr::bind_rows(run_summary), file.path(output_dir, "summary tables", "run_summary.tsv"), sep = "\t")
  }
  saveRDS(result_paths, file.path(output_dir, "summary tables", "result_paths.rds"))
  result_paths
}

lcl_methylation_by_cluster <- function(result) {
  dplyr::bind_rows(lapply(levels(result$assignments$cluster), function(cluster) {
    ids <- result$assignments$RID[result$assignments$cluster == cluster]
    m <- result$site_met_mat[ids, , drop = FALSE]
    n_calls <- sum(m > 0)
    denominator <- nrow(m) * ncol(m)
    data.frame(cluster = cluster, n_reads = nrow(m), n_sites = ncol(m),
      n_m6a_calls = n_calls, denominator = denominator,
      methylation_proportion = if (denominator) n_calls / denominator else NA_real_)
  }))
}

annotate_lcl_results <- function(result_paths, sample_table, phasing_root,
                                  haplotype_output_dir, reuse_cache = TRUE) {
  dir.create(file.path(haplotype_output_dir, "summary tables"), recursive = TRUE, showWarnings = FALSE)
  datasets <- lapply(result_paths, function(path) list(rids_df = readRDS(path)$assignments))
  phase_cache <- cache_lcl_haplotags(sample_table, datasets, phasing_root,
    file.path(haplotype_output_dir, "phase_cache"), reuse = reuse_cache)
  paths <- result_paths
  rows <- list()
  for (rid in names(result_paths)) {
    message("Annotating pooled reads: ", rid)
    result <- readRDS(result_paths[[rid]])
    phased <- annotate_lcl_haplotypes(list(rids_df = result$assignments),
      result$region, sample_table, phase_cache)
    result$assignments <- phased$reads
    result$group_id <- "pooled"
    result$variants <- phased$variants
    out <- file.path(haplotype_output_dir, rid, "pooled",
      paste0("bin", result$params$window_size), paste0("k", result$params$k_neighbors),
      paste0("resolution", result$params$resolution), "tables")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    paths[[rid]] <- file.path(out, "clustering.rds")
    saveRDS(result, paths[[rid]])
    data.table::fwrite(phased$reads, file.path(out, "read_cluster_assignments.tsv"), sep = "\t")
    data.table::fwrite(phased$reads, file.path(out, paste0(rid, "_read_haplotypes.tsv.gz")), sep = "\t")
    if (ncol(phased$variants)) data.table::fwrite(phased$variants,
      file.path(out, paste0(rid, "_phased_variants.tsv")), sep = "\t")
    snps <- lcl_window_snps(phased$variants, result$region)
    data.table::fwrite(snps, file.path(out, "plotted_snps.tsv"), sep = "\t")
    rows[[rid]] <- data.frame(region_id = rid, grouping = "pooled", n_reads = nrow(phased$reads),
      n_clusters = result$n_clusters, n_HP1 = sum(phased$reads$haplotype == "HP1"),
      n_HP2 = sum(phased$reads$haplotype == "HP2"), n_unphased = sum(phased$reads$haplotype == "unphased"),
      n_plotted_snps = nrow(snps))
  }
  data.table::fwrite(dplyr::bind_rows(rows), file.path(haplotype_output_dir, "summary tables", "haplotype_run_summary.tsv"), sep = "\t")
  saveRDS(paths, file.path(haplotype_output_dir, "summary tables", "result_paths.rds"))
  paths
}

write_lcl_result_index <- function(output_dir, include_haplotype = FALSE) {
  summary <- data.table::fread(file.path(output_dir, "summary tables", "run_summary.tsv"), data.table = FALSE)
  lines <- c("# LCL Leiden results", "",
    "m6A-only clustering and heatmaps. Fig 1 shows only nucleosomes (130-160 bp) and TF footprints (10-30, 40-60, 60-80 bp), in TNF colors.", "")
  for (hap in c(FALSE, if (include_haplotype) TRUE)) {
    lines <- c(lines, if (hap) "## Haplotype annotation of pooled clusters" else "## All samples pooled", "",
      paste0("| Region | Reads | Clusters | m6A heatmap | Read footprints | ",
        if (hap) "Occupancy lines | " else "", "Aggregate profile | Methylation bars |"),
      paste0("| --- | ---: | ---: | --- | --- | ", if (hap) "--- | " else "", "--- | --- |"))
    for (i in seq_len(nrow(summary))) {
      r <- summary[i, ]
      directory <- file.path(r$region_id, if (hap) "pooled" else "",
        paste0("bin", r$bin), paste0("k", r$k_neighbors), paste0("resolution", r$resolution), "plots")
      if (hap) directory <- file.path("..", "haplotype", directory)
      artifacts <- file.path(directory, c(if (hap) "heatmap_m6a_footprints.pdf" else "heatmap_m6a.pdf",
        "fig1_read_footprints.pdf", if (hap) "fig2_occupancy_by_cluster.pdf", "aggregate_profile.pdf", "methylation_proportion_by_cluster.pdf"))
      links <- ifelse(file.exists(file.path(output_dir, artifacts)), paste0("[PDF](", artifacts, ")"), "Pending")
      lines <- c(lines, paste0("| ", paste(c(r$annotation, r$n_reads, r$n_clusters, links), collapse = " | "), " |"))
    }
    lines <- c(lines, "")
  }
  lines <- c(lines, "m6A profiles connect observed-site proportions without smoothing or inserted zero-valued sites.",
    "Aggregate profiles use all lines for pooled results; haplotype aggregate profiles use footprint bars and an m6A line.",
    "Methylation bars average the observed m6A call fraction across sites in each cluster; they are not A/T-opportunity normalized.",
    "HP1 is light blue; HP2 is light yellow. Haplotype annotations never split or recluster the pooled reads.",
    "Notebook: `code/clustering_TSS_methods/Leiden_Manhattan/LCL/leiden_LCL.Rmd`.")
  writeLines(lines, file.path(output_dir, "index.md"))
  invisible(lines)
}
