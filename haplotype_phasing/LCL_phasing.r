# ---------------------------------------------------------------------------------
# filter for reads that can be phased at the lead AS-FIRE SNP
# 1. Find each sample's phase block containing the lead SNP
# 2. Assign each read HP1/HP2 from read level pahsing
# 3. Use the phased VCF genotype at the lead SNP to determine whether HP1/HP2 carries reference or alternative allele
# 	a. Only keep reads that have phased heterozygous focal SNPs
# 	b. Reads with HP1/HP2 
# 4. Subset m6a matrix to these reads
# ---------------------------------------------------------------------------------

# inputs
# 1. dat: read table and m6a matrix for the rergion
# 2. region: selected AS-FIRE region with the lead SNP
# 3. sample table: LCL samples containing the SNP
# 4. phase cache: phase blocks and read level haplotags

lcl_filter_focal_heterozygotes <- function(dat, region, sample_table, phase_cache) {
  # define lead SNP: chr, pos, variant id, ref, alt
  focal <- data.frame(chr = region$chr, pos = region$focal_pos,
    variant_id = region$focal_snp, ref = region$ref, alt = region$alt)
  # assign haplotags (HP1/HP2) to reads 
  phased <- annotate_lcl_haplotypes(dat, region, sample_table, phase_cache)
  # convert haplotags to genotype
  reads <- lcl_allele_labels(phased$reads, phased$variants, focal)
  reads$focal_genotype <- NA_character_
  reads$retained_for_analysis <- FALSE
  # restrcit variants to the lead SNP
  v <- phased$variants
  if (nrow(v)) {
    v <- v[v$chr == region$chr & v$pos == region$focal_pos &
      v$ref == region$ref & v$alt == region$alt & v$genotype %in% c("0|1", "1|0"), , drop = FALSE]
    # sample and phase block key
    keys <- paste(v$sample_name, v$phase_set, sep = "::")
    for (i in seq_len(nrow(reads))) {
      gt <- unique(v$genotype[keys == paste(reads$sample_name[i], reads$phase_set[i], sep = "::")])
      if (length(gt) == 1L) reads$focal_genotype[i] <- gt
    }
    reads$retained_for_analysis <- !is.na(reads$focal_genotype) &
      reads$haplotype %in% c("HP1", "HP2") & reads$allele_status == "phased_focal_genotype"
  }
  kept <- reads[reads$retained_for_analysis, , drop = FALSE]
  if (nrow(kept) < 3L || length(unique(kept$allele_label)) != 2L)
    stop("Selected peak has insufficient resolved full-span reads or lacks one allele: ", region$region_id)
  dat$rids_df <- kept
  dat$met_mat <- dat$met_mat[kept$RID, , drop = FALSE]
  dat$met_mat <- dat$met_mat[, Matrix::colSums(dat$met_mat) > 0, drop = FALSE]
  if (ncol(dat$met_mat) < 2L) stop("Fewer than two observed m6A sites after filtering: ", region$region_id)
  dat$read_filter_audit <- reads
  dat$variants <- phased$variants
  dat$focal <- focal
  dat
}


# ---------------------------------------------------------------------------------
# locate VCF file for each sample
# ---------------------------------------------------------------------------------

lcl_phasing_paths <- function(sample_name, phasing_root) {
  directory <- file.path(phasing_root, sample_name)
  c(vcf = file.path(directory, paste0(sample_name, ".5mC.6mA.aligned.phased.vcf.gz")),
    tags = file.path(directory, "read-level-phasing.tsv"),
    blocks = file.path(directory, "blocks.tsv"))
}


# ---------------------------------------------------------------------------------
# find the genotype column for the requested cell line
# ---------------------------------------------------------------------------------

# different LCL samples may come from different cell lines

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


# ---------------------------------------------------------------------------------
# collects each sample’s read IDs, loads its phase-block table and matching read haplotags
# ---------------------------------------------------------------------------------

# the output is saved as a RDS file - a list with blocks and read level haplotypes

cache_lcl_haplotags <- function(sample_table, region_data, phasing_root, cache_dir,
                                reuse = TRUE) {
  cache_dir <- lcl_output_path(cache_dir)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  read_info <- dplyr::bind_rows(lapply(region_data, function(dataset) dataset$rids_df))
  results <- setNames(vector("list", nrow(sample_table)), sample_table$sample_name)
  for (sample_index in seq_len(nrow(sample_table))) {
    sample <- sample_table[sample_index, , drop = FALSE]
    paths <- lcl_phasing_paths(sample$sample_name, phasing_root)
    required <- c(unname(paths), paste0(paths[["vcf"]], ".tbi"))
    original_ids <- sort(unique(read_info$original_RID[read_info$sample_name == sample$sample_name]))
    signature <- list(version = 1L, cell_line = sample$cell_line, reads = original_ids,
      required = required)
    cache <- file.path(cache_dir, paste0(sample$sample_name, "_haplotags.rds"))
    if (reuse && file.exists(cache)) {
      previous <- tryCatch(readRDS(cache), error = function(e) {
        message("Rebuilding unreadable haplotag cache: ", cache)
        NULL
      })
      if (!is.null(previous) && identical(previous$signature[names(signature)], signature) &&
          identical(previous$missing_files, required[!file.exists(required)])) {
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
          ids_file <- tempfile(tmpdir = cache_dir, fileext = ".txt")
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

# ---------------------------------------------------------------------------------
# extract variants' genotype from VCF within the selected window and phase block
# ---------------------------------------------------------------------------------


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
    
    # match genotype carried by HP1 and HP2 using the VCF file
    genotype <- values[match("GT", keys)]
    block <- values[match("PS", keys)]
    if (is.na(genotype) || is.na(block) || block != phase_set || !grepl("|", genotype, fixed = TRUE)) return(NULL)
    allele_indices <- suppressWarnings(as.integer(strsplit(genotype, "|", fixed = TRUE)[[1]]))
    alleles <- c(fields[4], strsplit(fields[5], ",", fixed = TRUE)[[1]])
    if (length(allele_indices) != 2L || anyNA(allele_indices) ||
        any(allele_indices < 0L | allele_indices >= length(alleles))) return(NULL)
    data.frame(chr = fields[1], pos = as.integer(fields[2]), variant_id = fields[3],
      ref = fields[4], alt = fields[5], genotype = genotype, phase_set = block,
      HP1_allele = alleles[allele_indices[1] + 1L], HP2_allele = alleles[allele_indices[2] + 1L])
  })
  dplyr::bind_rows(empty, dplyr::bind_rows(variants))
}


# ---------------------------------------------------------------------------------
# match each read to its sample's phase block and read level haplotag
# ---------------------------------------------------------------------------------

# return annotated reads and associated variants

annotate_lcl_haplotypes <- function(dataset, region, sample_table, phase_cache,
                                    variant_flank = 50000L) {
  stopifnot(variant_flank >= 0L)
  reads <- dataset$rids_df
  reads$haplotype <- "unphased"
  reads$phase_set <- NA_character_
  reads$phase_status <- "unassigned_read"
  reads$phase_source <- "HiPhase read-level-phasing.tsv + sample VCF"
  variants <- list()
  anchor <- if (!is.null(region$focal_pos) && !is.na(region$focal_pos)) region$focal_pos else
    if (is.na(region$tss)) floor(mean(c(region$analysis_start, region$analysis_end))) else region$tss
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
      reads$phase_status[read_indices] <- "no_local_phased_variant"
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

# ---------------------------------------------------------------------------------
# Converts each read’s sample specific HP1/HP2 assignment into lead SNP allele
# ---------------------------------------------------------------------------------
lcl_allele_labels <- function(reads, variants, focal) {
  reads$allele_label <- "unphased"
  reads$allele_status <- "unphased"
  phased <- reads$haplotype %in% c("HP1", "HP2") & !is.na(reads$phase_set)
  reads$allele_label[phased] <- paste0(sub("_.*$", "", reads$sample_name[phased]),
    ":", reads$haplotype[phased], " (PS ", reads$phase_set[phased], ")")
  reads$allele_status[phased] <- "sample_local_haplotype"
  if (!nrow(focal) || !nrow(variants)) return(reads)
  v <- variants[variants$chr == focal$chr[1] & variants$pos == focal$pos[1] &
                  variants$ref == focal$ref[1] & variants$alt == focal$alt[1], , drop = FALSE]
  for (i in which(phased)) {
    hit <- v[v$sample_name == reads$sample_name[i] & v$phase_set == reads$phase_set[i], , drop = FALSE]
    bases <- unique(hit[[paste0(reads$haplotype[i], "_allele")]])
    if (length(bases) == 1L && !is.na(bases) && bases %in% c("A", "C", "G", "T")) {
      reads$allele_label[i] <- paste0(focal$variant_id[1], ": ", bases)
      reads$allele_status[i] <- "phased_focal_genotype"
    } else {
      reads$allele_label[i] <- paste0("unresolved at SNP | ", reads$allele_label[i])
      reads$allele_status[i] <- "focal_genotype_unresolved"
    }
  }
  reads
}
