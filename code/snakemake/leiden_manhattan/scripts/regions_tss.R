# Canonical GENCODE promoters. Stored start/end use BED coordinates.
select_tss_regions <- function(spec) {
  required <- c("tss_bed", "flank_bp", "genes")
  if (!all(required %in% names(spec))) stop("TSS configuration needs: ", paste(required, collapse = ", "))
  flank <- as.numeric(spec$flank_bp)
  if (length(flank) != 1L || !is.finite(flank) || flank < 1 || flank != as.integer(flank))
    stop("regions.tss.flank_bp must be a positive integer")
  convention <- if (is.null(spec$width_convention)) "half_open_2000bp" else spec$width_convention
  if (!convention %in% c("half_open_2000bp", "inclusive"))
    stop("Unknown TSS width_convention: ", convention)
  genes <- unlist(spec$genes, use.names = TRUE)
  if (!length(genes) || is.null(names(genes)) || any(!nzchar(names(genes))))
    stop("regions.tss.genes must map region IDs to GENCODE gene names")
  bed <- read.delim(spec$tss_bed, header = FALSE, stringsAsFactors = FALSE,
                    quote = "", comment.char = "", check.names = FALSE)
  if (ncol(bed) != 6L) stop("Expected six columns in canonical TSS BED: ", spec$tss_bed)
  names(bed) <- c("chr", "start", "end", "annotation", "score", "strand")
  annotations <- strsplit(as.character(bed$annotation), ";", fixed = TRUE)
  if (any(lengths(annotations) < 3L)) stop("Invalid canonical TSS annotation in ", spec$tss_bed)
  bed$gencode_name <- vapply(annotations, `[[`, character(1), 3L)
  bed$ensg <- sub("\\..*$", "", vapply(annotations, `[[`, character(1), 1L))
  bed <- bed[bed$chr %in% paste0("chr", c(1:22, "X", "Y")) & bed$gencode_name %in% genes, , drop = FALSE]
  if (anyDuplicated(bed$gencode_name)) stop("Multiple canonical TSS records for configured gene(s)")
  idx <- match(as.character(genes), bed$gencode_name)
  if (anyNA(idx)) stop("Canonical TSS not found for: ", paste(genes[is.na(idx)], collapse = ", "))
  bed <- bed[idx, , drop = FALSE]
  if (any(as.integer(bed$end) - as.integer(bed$start) != 1L))
    stop("Canonical TSS BED rows must be 1 bp (start = 0-based TSS): ", spec$tss_bed)
  tss <- as.integer(bed$start) + 1L
  first <- tss - as.integer(flank)
  last <- tss + as.integer(flank) - if (convention == "half_open_2000bp") 1L else 0L
  if (any(first < 1L)) stop("Promoter window extends before chromosome position 1")
  sample_mode <- if (is.null(spec$samples)) "all" else spec$samples
  if (!identical(sample_mode, "all")) stop("TSS regions require samples: all")
  data.frame(region_id = names(genes), chr = as.character(bed$chr), start = first - 1L,
             end = last, analysis_start = first, analysis_end = last, width = last - first + 1L,
             region_type = "promoter", strand = as.character(bed$strand), tss = tss,
             annotation = paste(names(genes), "canonical promoter"), gene = names(genes),
             gencode_name = as.character(genes), ensg = bed$ensg,
             samples_mode = sample_mode,
             coordinate_system = "BED: 0-based, half-open; analysis: 1-based inclusive",
             stringsAsFactors = FALSE, row.names = NULL)
}
