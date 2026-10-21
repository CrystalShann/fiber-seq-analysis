# Shared functions for consensus construction and peak-TPM-filtered divergent TSS pairing.
suppressPackageStartupMessages(library(data.table))

validate_gzip <- function(path) {
    if (file.info(path)$size <= 0) 
        stop("Empty input: ", path)
    con <- gzfile(path, "rb")
    on.exit(close(con))
    withCallingHandlers(while (length(readBin(con, "raw", 1048576L))) {
    }, warning = function(w) stop(conditionMessage(w)))
}

read_gzip_table <- function(path, ...) {
    con <- gzfile(path, "rt")
    on.exit(close(con))
    fread(text = paste(readLines(con, warn = FALSE), collapse = "\n"), sep = "\t", quote = "", na.strings = "", 
        ...)
}

read_gencode_genes <- function(path, chunk_size = 50000L) {
    con <- if (grepl("\\.gz$", path)) 
        gzfile(path, "rt")
    else file(path, "rt")
    on.exit(close(con))
    parts <- list()
    repeat {
        lines <- readLines(con, n = chunk_size, warn = FALSE)
        if (!length(lines)) 
            break
        lines <- lines[grepl("^[^#][^\t]*\t[^\t]*\tgene\t", lines)]
        if (!length(lines)) 
            next
        x <- fread(text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE, quote = "")
        attr <- function(key) sub(paste0(".*(?:^|; )", key, " \"([^\"]+)\".*"), "\\1", x$V9, perl = TRUE)
        parts[[length(parts) + 1L]] <- data.table(gene_id = sub("\\.[0-9]+$", "", attr("gene_id")), gene_id_version = attr("gene_id"), 
            gene_name = attr("gene_name"), gene_type = attr("gene_type"), gene_chrom = x$V1, gene_strand = x$V7)
    }
    genes <- rbindlist(parts)
    if (!nrow(genes) || anyDuplicated(genes$gene_id)) 
        stop("Empty or duplicate GENCODE genes")
    genes
}

collapse_values <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = "|")

extract_ids <- function(x, pattern) {
    if (is.na(x) || !nzchar(x)) 
        return(character())
    unique(regmatches(x, gregexpr(pattern, x, perl = TRUE))[[1]])
}

make_id_maps <- function(hgnc) {
    stopifnot(all(c("hgnc_id", "entrez_id", "ensembl_gene_id") %in% names(hgnc)))
    hgnc <- hgnc[!is.na(ensembl_gene_id) & nzchar(ensembl_gene_id)]
    list(hgnc = split(hgnc$ensembl_gene_id, hgnc$hgnc_id), entrez = split(hgnc$ensembl_gene_id, as.character(hgnc$entrez_id)))
}

resolve_gene <- function(hgnc_ids, entrez_ids, maps) {
    h <- extract_ids(hgnc_ids, "HGNC:[0-9]+")
    e <- extract_ids(entrez_ids, "[0-9]+")
    hg <- unique(unlist(maps$hgnc[h], use.names = FALSE))
    eg <- unique(unlist(maps$entrez[e], use.names = FALSE))
    hg <- unique(unlist(lapply(hg, extract_ids, pattern = "ENSG[0-9]+")))
    eg <- unique(unlist(lapply(eg, extract_ids, pattern = "ENSG[0-9]+")))
    reason <- if (length(h) > 1L || length(e) > 1L) 
        "multiple_gene_ids"
    else if (length(hg) > 1L || length(eg) > 1L) 
        "ambiguous_mapping"
    else if (length(hg) && length(eg) && !identical(hg, eg)) 
        "conflicting_mapping"
    else if (!length(hg) && !length(eg)) 
        "unmapped_gene"
    else "mapped"
    list(gene_id = if (reason == "mapped") {
        if (length(hg)) hg else eg
    } else NA_character_, mapping_method = if (reason != "mapped") NA_character_ else if (length(hg)) "HGNC" else "Entrez", 
        mapping_reason = reason)
}

map_promoter_peaks <- function(bed, peak_names, annot, genes, hgnc) {
    stopifnot(!anyDuplicated(bed$peak_id), !anyDuplicated(peak_names$peak_id))
    a <- annot[, .(hgnc_ids = collapse_values(hgnc_ids), entrez_ids = collapse_values(entrez_ids), annotated_symbols = collapse_values(annotated_symbols)), 
        by = peak_id]
    x <- merge(bed, peak_names[, .(peak_id, peak_name)], by = "peak_id", all.x = TRUE, sort = FALSE)
    x <- merge(x, a, by = "peak_id", all.x = TRUE, sort = FALSE)
    maps <- make_id_maps(hgnc)
    mapped <- rbindlist(Map(function(h, e) resolve_gene(h, e, maps), x$hgnc_ids, x$entrez_ids))
    x <- cbind(x, mapped)
    x <- merge(x, genes, by = "gene_id", all.x = TRUE, sort = FALSE)
    x[, `:=`(exclusion_reason, mapping_reason)]
    x[mapping_reason == "mapped", `:=`(exclusion_reason, "retained")]
    x[mapping_reason == "mapped" & is.na(gene_chrom), `:=`(exclusion_reason, "absent_from_gencode")]
    x[mapping_reason == "mapped" & !is.na(gene_chrom) & (chrom != gene_chrom | strand != gene_strand), `:=`(exclusion_reason, 
        "chromosome_or_strand_mismatch")]
    multi_symbol <- vapply(x$annotated_symbols, function(s) {
        if (is.na(s) || !nzchar(s)) 
            return(FALSE)
        length(unique(strsplit(s, "[,|]")[[1]])) > 1L
    }, logical(1))
    x[(!is.na(peak_name) & grepl(",", peak_name, fixed = TRUE)) | multi_symbol, `:=`(exclusion_reason, "multiple_gene_names")]
    x[is.na(peak_name) | !nzchar(peak_name) | grepl("^hg_[0-9]+[.][0-9]+$", peak_name) | grepl("^(?:p[0-9]*@)?(?:hg[0-9]+::)?chr[[:alnum:]_]+:", 
        peak_name, perl = TRUE), `:=`(exclusion_reason, "unannotated_peak")]
    x
}

annotate_pairs <- function(pairs) {
    x <- copy(pairs)
    x[, `:=`(anchor_start0 = as.integer(floor((minus_tss0 + plus_tss0)/2)), directionality = abs(plus_signal - 
        minus_signal)/(plus_signal + minus_signal), dominant_partner = fifelse(plus_signal > minus_signal, "plus", 
        fifelse(minus_signal > plus_signal, "minus", "tie")), dominant_peak_id = fifelse(plus_signal > minus_signal, 
        plus_peak_id, fifelse(minus_signal > plus_signal, minus_peak_id, NA_character_)), dominant_gene_id = fifelse(plus_signal > 
        minus_signal, plus_gene_id, fifelse(minus_signal > plus_signal, minus_gene_id, NA_character_)), anchor_strand = fifelse(plus_signal > 
        minus_signal, "+", fifelse(minus_signal > plus_signal, "-", ".")), pair_name = paste(fifelse(!is.na(minus_gene_name), 
        minus_gene_name, minus_peak_id), fifelse(!is.na(plus_gene_name), plus_gene_name, plus_peak_id), sep = "–"))]
    x[, `:=`(anchor_end0, anchor_start0 + 1L)]
    x
}

validate_main_pairs <- function(x, threshold, min_distance, max_distance) {
    stopifnot(!anyDuplicated(x$pair_id), !anyDuplicated(x[, .(chrom, minus_tss0)]), all(x$distance_bp >= min_distance & 
        x$distance_bp <= max_distance), all(x$minus_signal >= threshold & x$plus_signal >= threshold), all(x$minus_peak_id != 
        x$plus_peak_id), all(x$anchor_start0 >= 0L & x$anchor_end0 == x$anchor_start0 + 1L), all(x$directionality >= 
        0 & x$directionality <= 1))
    invisible(TRUE)
}

read_ctss <- function(path, threshold = 1) {
    x <- read_gzip_table(path, header = FALSE)
    if (ncol(x) < 6L) 
        stop("CTSS input requires at least BED6")
    x <- x[, 1:6]
    setnames(x, c("chrom", "tss0", "end0", "peak_id", "signal", "strand"))
    if (!is.numeric(x$signal) || any(!is.finite(x$signal)) || any(x$signal < 0) || any(x$signal%%1 != 0)) 
        stop("CTSS scores must be finite nonnegative integer counts")
    stopifnot(is.finite(threshold), threshold > 0, all(x$end0 == x$tss0 + 1), all(is.finite(x$tss0) & x$tss0 >= 
        0 & x$tss0%%1 == 0), all(x$strand %in% c("+", "-")), !anyNA(x$peak_id), !anyDuplicated(x$peak_id), !anyDuplicated(x[, 
        .(chrom, tss0, strand)]))
    x[, `:=`(gene_id = NA_character_, gene_name = NA_character_, gene_type = NA_character_, annotation_status = "no_CAGE_peak_overlap")]
    list(all = x, expressed = x[signal >= threshold])
}

read_replicate_counts <- function(path, chunk_size = 200000L) {
    con <- gzfile(path, "rt")
    on.exit(close(con))
    parts <- list()
    withCallingHandlers(repeat {
        lines <- readLines(con, n = chunk_size, warn = FALSE)
        if (!length(lines)) 
            break
        tab <- fread(text = paste0(paste(lines, collapse = "\n"), "\n"), header = FALSE, sep = "\t", quote = "", 
            select = c(1L, 2L, 3L, 5L, 6L))
        setnames(tab, c("chrom", "start0", "end0", "count", "strand"))
        if (!is.numeric(tab$count) || anyNA(tab) || any(!is.finite(tab$count)) || any(tab$count <= 0 | tab$count%%1 != 
            0) || any(tab$end0 != tab$start0 + 1) || any(tab$start0 < 0 | tab$start0%%1 != 0) || any(!tab$strand %in% 
            c("+", "-"))) 
            stop("Invalid single-base CTSS record: ", path)
        parts[[length(parts) + 1L]] <- tab
    }, warning = function(w) stop(conditionMessage(w)))
    x <- rbindlist(parts)
    if (!nrow(x)) 
        stop("Empty CTSS file: ", path)
    if (anyDuplicated(x[, .(chrom, start0, strand)])) 
        stop("Duplicate CTSS within replicate: ", path)
    x
}

consensus_bed <- function(common, minimum_count = 1) {
    stopifnot(length(minimum_count) == 1L, is.finite(minimum_count), minimum_count >= 1, minimum_count%%1 == 0)
    common[minimum_replicate_count >= minimum_count, .(chrom, start0, end0, name = paste0(chrom, ":", start0, "..", 
        end0, ",", strand), score = summed_count, strand, thickStart = start0, thickEnd = end0, itemRgb = "0,0,0")]
}

# Stream only the three GM12878 expression columns, averaging zeros as well.
read_gm12878_peak_tpm <- function(path, chunk_size = 500L) {
  ids <- c("CNhs12331", "CNhs12332", "CNhs12333")
  con <- gzfile(path, "rt"); on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("Missing TPM matrix header")
    if (!startsWith(line, "#")) break
  }
  header <- strsplit(line, "\t", fixed = TRUE)[[1]]
  columns <- vapply(ids, function(id) {
    hit <- which(startsWith(header, "tpm.") & grepl(id, header, fixed = TRUE) & grepl("GM12878", header, fixed = TRUE))
    if (length(hit) != 1L) stop("Expected exactly one GM12878 TPM column: ", id)
    hit
  }, integer(1))
  message("Averaging peak TPM across ", paste(ids, collapse = ", "), "; zeros included")
  parts <- list()
  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    x <- fread(text = paste0(paste(lines, collapse = "\n"), "\n"), header = FALSE,
      select = unname(c(1L, columns)), sep = "\t", quote = "", colClasses = "character", fill = FALSE)
    setnames(x, c("peak_id", ids))
    x <- x[!grepl("^[0-9]+STAT:", peak_id)]
    for (id in ids) {
      values <- suppressWarnings(as.numeric(x[[id]]))
      if (any(!is.finite(values) | values < 0)) stop("Invalid expression value in ", id)
      set(x, j = id, value = values)
    }
    x[, mean_peak_tpm := rowMeans(.SD), .SDcols = ids]
    parts[[length(parts) + 1L]] <- x[, .(peak_id, mean_peak_tpm)]
  }
  result <- rbindlist(parts)
  if (!nrow(result) || anyDuplicated(result$peak_id)) stop("Empty or duplicate TPM peak IDs")
  result
}

# Gene mapping is descriptive only and remains entirely in memory.
prepare_gene_annotations <- function(input_dir, gencode, bed) {
  nm <- read_gzip_table(file.path(input_dir, "human_phase1and2_CAGE_Peak_name.txt.gz"), header = TRUE, fill = TRUE)
  setnames(nm, names(nm)[1:2], c("peak_id", "peak_name"))
  ann <- read_gzip_table(file.path(input_dir, "hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt.gz"),
    header = TRUE, fill = TRUE, colClasses = "character")
  setnames(ann, c("#CAGE_Peak_ID", "HGNC/MGI_ID", "GeneID", "Gene_symbol"),
    c("peak_id", "hgnc_ids", "entrez_ids", "annotated_symbols"))
  hgnc <- fread(file.path(input_dir, "hgnc_complete_set.txt"), sep = "\t", quote = "", colClasses = "character", na.strings = "")
  genes <- read_gencode_genes(gencode)
  map_promoter_peaks(bed, nm[peak_id %in% bed$peak_id], ann[peak_id %in% bed$peak_id], genes, hgnc)
}

# Retain a CTSS only when its base overlaps an expressed peak on its own strand.
# Unknown/ambiguous gene assignments never exclude a qualifying CTSS.
filter_and_annotate_ctss <- function(ctss, reference) {
  ref <- copy(reference)
  ref[, `:=`(start = peak_start0, end = peak_end0 - 1L)]
  setkey(ref, chrom, strand, start, end)
  query <- ctss[, .(peak_id, chrom, strand, start = tss0, end = tss0)]
  hits <- foverlaps(query, ref, by.x = c("chrom", "strand", "start", "end"), nomatch = 0L)
  x <- copy(ctss[peak_id %in% hits$peak_id])
  x[, `:=`(cage_peak_ids = NA_character_, mean_peak_tpm = NA_real_, annotation_status = "CAGE_peak_unresolved")]
  if (!nrow(hits)) return(x)
  support <- hits[, .(cage_peak_ids = paste(sort(unique(cage_peak_id)), collapse = "|"),
                      mean_peak_tpm = max(mean_peak_tpm)), by = peak_id]
  x[support, on = "peak_id", `:=`(cage_peak_ids = i.cage_peak_ids, mean_peak_tpm = i.mean_peak_tpm)]
  confident <- hits[exclusion_reason == "retained" & !is.na(gene_id)]
  assignments <- confident[, .(n_genes = uniqueN(gene_id), gene_id = gene_id[1L],
    gene_name = gene_name[1L], gene_type = gene_type[1L]), by = peak_id]
  x[peak_id %in% assignments[n_genes > 1L, peak_id], annotation_status := "ambiguous_gene"]
  x[assignments[n_genes == 1L], on = "peak_id", `:=`(gene_id = i.gene_id, gene_name = i.gene_name,
    gene_type = i.gene_type, annotation_status = "retained")]
  x
}

# Independent nearest-plus pairing: no rejected-candidate or decision tables.
find_divergent_pairs <- function(peaks, min_distance = 1L, max_distance = 500L) {
  parts <- list()
  for (ch in sort(unique(peaks$chrom))) {
    minus <- peaks[chrom == ch & strand == "-"][order(tss0, peak_id)]
    plus <- peaks[chrom == ch & strand == "+"][order(tss0, peak_id)]
    if (!nrow(minus) || !nrow(plus)) next
    j <- findInterval(minus$tss0 + min_distance - 1L, plus$tss0) + 1L
    keep <- which(j <= nrow(plus) & !is.na(plus$tss0[j]) & plus$tss0[j] - minus$tss0 <= max_distance)
    if (length(keep)) parts[[length(parts) + 1L]] <- data.table(minus_peak_id = minus$peak_id[keep], plus_peak_id = plus$peak_id[j[keep]])
  }
  pairs <- if (length(parts)) rbindlist(parts) else data.table(minus_peak_id = character(), plus_peak_id = character())
  fields <- c("peak_id", "gene_id", "gene_name", "gene_type", "tss0", "signal", "annotation_status", "cage_peak_ids", "mean_peak_tpm")
  for (side in c("minus", "plus")) {
    tab <- copy(peaks[, ..fields]); setnames(tab, paste0(side, "_", fields))
    pairs <- merge(pairs, tab, by = paste0(side, "_peak_id"), all.x = TRUE, sort = FALSE)
  }
  pairs[, `:=`(chrom = peaks$chrom[match(minus_peak_id, peaks$peak_id)],
    distance_bp = plus_tss0 - minus_tss0, pair_id = paste(minus_peak_id, plus_peak_id, sep = "__"))]
  setorder(pairs, chrom, minus_tss0, plus_tss0, minus_peak_id, plus_peak_id)
  pairs
}

write_anchor_bed <- function(pairs, path) {
  first <- pairs[, .(chrom, start = anchor_start0, end = anchor_end0,
    name = pair_name, score = 0L, strand = anchor_strand)]
  extra <- c("pair_id", "minus_gene_id", "plus_gene_id", "minus_gene_type", "plus_gene_type",
    "minus_peak_id", "plus_peak_id", "minus_tss0", "plus_tss0", "minus_signal", "plus_signal",
    "distance_bp", "directionality", "dominant_partner", "dominant_peak_id", "dominant_gene_id",
    "signal_unit", "minus_annotation_status", "plus_annotation_status",
    "minus_cage_peak_ids", "plus_cage_peak_ids", "minus_mean_peak_tpm", "plus_mean_peak_tpm")
  bed <- cbind(first, pairs[, ..extra])
  part <- paste0(path, ".part")
  fwrite(bed, part, sep = "\t", col.names = FALSE, quote = FALSE, na = ".")
  if (!file.rename(part, path)) stop("Cannot finalize BED: ", path)
}

run_ctss_workflow <- function(cfg) {
  ctss <- read_ctss(file.path(cfg$input_dir, cfg$ctss_file), cfg$score_threshold)$expressed
  reference <- read_gzip_table(file.path(cfg$input_dir, "fantom5.hg38.CAGE_peaks.bed.gz"), header = FALSE)
  if (ncol(reference) != 9L) stop("Expected reference CAGE BED9")
  setnames(reference, c("chrom", "peak_start0", "peak_end0", "peak_id", "bed_score", "strand", "tss0", "tss_end0", "rgb"))
  expression <- read_gm12878_peak_tpm(file.path(cfg$input_dir, cfg$tpm_file), cfg$chunk_size)
  index <- match(reference$peak_id, expression$peak_id)
  if (anyNA(index)) stop("Reference CAGE peaks missing from TPM matrix")
  reference[, mean_peak_tpm := expression$mean_peak_tpm[index]]
  reference <- reference[mean_peak_tpm >= cfg$tpm_threshold]
  if (!nrow(reference)) stop("No reference CAGE peaks meet the TPM threshold")
  message(nrow(reference), " CAGE peaks pass mean TPM >= ", cfg$tpm_threshold)
  message("Preparing descriptive gene labels from local references")
  mapped <- prepare_gene_annotations(cfg$input_dir, cfg$gencode, reference)
  reference <- mapped[, .(cage_peak_id = peak_id, chrom, peak_start0, peak_end0, strand,
    mean_peak_tpm, gene_id, gene_name, gene_type, exclusion_reason)]
  pairing <- filter_and_annotate_ctss(ctss, reference)
  message(nrow(pairing), " of ", nrow(ctss), " consensus CTSS positions overlap qualifying peaks on the same strand")
  pairs <- annotate_pairs(find_divergent_pairs(pairing, cfg$min_distance, cfg$max_distance))
  pairs[, signal_unit := "CAGE_tag_count"]
  validate_main_pairs(pairs, cfg$score_threshold, cfg$min_distance, cfg$max_distance)
  stopifnot(all(pairs$minus_mean_peak_tpm >= cfg$tpm_threshold & pairs$plus_mean_peak_tpm >= cfg$tpm_threshold))
  path <- file.path(cfg$output_dir, "divergent_promoter_anchors.bed")
  tables <- list(divergent_pairs.tsv = pairs, pairing_peaks.tsv = pairing,
    ctss_annotation_counts.tsv = pairing[, .(count = .N), by = annotation_status][order(annotation_status)])
  for (name in names(tables)) {
    destination <- file.path(cfg$output_dir, name)
    fwrite(tables[[name]], paste0(destination, ".part"), sep = "\t", quote = FALSE, na = ".")
    if (!file.rename(paste0(destination, ".part"), destination)) stop("Cannot finalize table: ", destination)
  }
  write_anchor_bed(pairs, path)
  message("Wrote ", nrow(pairs), " pairs to ", path)
  invisible(list(pairs = pairs, pairing = pairing, expressed_peaks = reference, input_CTSS = nrow(ctss)))
}

intersect_replicate_counts <- function(paths) {
  stopifnot(length(paths) >= 2L, !anyDuplicated(paths))
  common <- NULL
  for (i in seq_along(paths)) {
    message("Reading replicate ", i)
    x <- read_replicate_counts(paths[i])
    setnames(x, "count", paste0("rep", i, "_count"))
    common <- if (is.null(common)) x else merge(common, x, by = c("chrom", "start0", "end0", "strand"), sort = FALSE)
    rm(x); gc(verbose = FALSE)
  }
  count_cols <- paste0("rep", seq_along(paths), "_count")
  common[, summed_count := rowSums(.SD), .SDcols = count_cols]
  common[, minimum_replicate_count := do.call(pmin, .SD), .SDcols = count_cols]
  setorder(common, chrom, start0, strand)
  common
}
