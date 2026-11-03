# Match scipy.signal.find_peaks: midpoint (rounded down) for a flat local peak;
# endpoints excluded; equal-height peaks resolve to the first/lower lag.
repeat_peak <- function(profile, lower = 140L, upper = 250L) {
  missing <- c(lag = NA_real_, height = NA_real_)
  if (length(profile) <= upper + 1L || any(!is.finite(profile))) return(missing)
  runs <- rle(profile)
  if (length(runs$values) < 3L) return(missing)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  middle <- seq.int(2L, length(runs$values) - 1L)
  candidates <- middle[runs$values[middle] > runs$values[middle - 1L] &
                         runs$values[middle] > runs$values[middle + 1L]]
  indices <- floor((starts[candidates] + ends[candidates]) / 2)
  lags <- indices - 1L
  indices <- indices[lags >= lower & lags <= upper & profile[indices] > 0]
  if (!length(indices)) return(missing)
  peak <- indices[which.max(profile[indices])]
  c(lag = peak - 1L, height = profile[peak])
}

summarize_autocorrelations <- function(profiles, records) {
  records <- as.data.frame(records)
  stopifnot(nrow(profiles) == nrow(records))
  clustered <- records[records$status == "clustered", , drop = FALSE]
  clusters <- unique(as.character(clustered$cluster))
  clusters <- clusters[order(as.integer(clusters))]
  samples <- sort(unique(clustered$sample_name))
  averages <- data.frame(region_id=character(), cluster=character(), lag_bp=integer(),
    mean_acf=numeric(), sd_acf=numeric(), sem_acf=numeric(), n_reads=integer())
  summaries <- data.frame(region_id=character(), cluster=character(), n_reads=integer(),
    fraction_clustered=numeric(), mean_m6a_call_fraction=numeric(),
    positive_local_peak_140_250_bp=numeric(), peak_acf=numeric())
  composition <- data.frame(region_id=character(), cluster=character(), sample_name=character(),
    n_reads=integer(), fraction_within_cluster=numeric(), fraction_within_sample=numeric())
  for (cluster in clusters) {
    keep <- records$cluster == cluster & records$status == "clustered"
    subset <- profiles[keep, , drop = FALSE]
    n <- nrow(subset)
    average <- colMeans(subset)
    sd <- if (n > 1L) apply(subset, 2L, stats::sd) else rep(NA_real_, ncol(profiles))
    peak <- repeat_peak(average)
    region <- records$region_id[1L]
    summaries <- rbind(summaries, data.frame(region_id=region, cluster=cluster,
      n_reads=n, fraction_clustered=n/nrow(clustered),
      mean_m6a_call_fraction=mean(records$m6a_call_fraction[keep], na.rm=TRUE),
      positive_local_peak_140_250_bp=unname(peak[1L]), peak_acf=unname(peak[2L])))
    averages <- rbind(averages, data.frame(region_id=region, cluster=cluster,
      lag_bp=seq_len(ncol(profiles))-1L, mean_acf=average, sd_acf=sd,
      sem_acf=sd/sqrt(n), n_reads=n))
    for (sample in samples) {
      count <- sum(clustered$cluster == cluster & clustered$sample_name == sample)
      composition <- rbind(composition, data.frame(region_id=region, cluster=cluster,
        sample_name=sample, n_reads=count, fraction_within_cluster=count/n,
        fraction_within_sample=count/sum(clustered$sample_name == sample)))
    }
  }
  list(averages=averages, summaries=summaries, composition=composition)
}

annotate_alleles <- function(records, region, sample_table) {
  records <- as.data.frame(records)
  stopifnot(nrow(region) == 1L, region$ref %in% c("A", "C", "G", "T"),
    region$alt %in% c("A", "C", "G", "T"), region$ref != region$alt,
    all(records$haplotype %in% c("HP1", "HP2")),
    all(records$focal_genotype %in% c("0|1", "1|0")),
    all(records$allele_status == "phased_focal_genotype"),
    all(c("sample_name", "cell_line") %in% names(sample_table)),
    !anyDuplicated(sample_table$sample_name))
  records$cell_line <- sample_table$cell_line[match(records$sample_name, sample_table$sample_name)]
  if (anyNA(records$cell_line) || any(!nzchar(records$cell_line)))
    stop("Every input sample needs an explicit biological cell_line in sample metadata")
  alt <- (records$haplotype == "HP2" & records$focal_genotype == "0|1") |
         (records$haplotype == "HP1" & records$focal_genotype == "1|0")
  records$allele_group <- ifelse(alt, "ALT", "REF")
  records$allele_base <- ifelse(alt, region$alt, region$ref)
  if (!all(records$allele_label == paste0(region$focal_snp, ": ", records$allele_base)))
    stop("HP/genotype mapping disagrees with the saved focal allele label")
  records
}

summarize_alleles <- function(profiles, records, region) {
  records <- as.data.frame(records)
  stopifnot(is.matrix(profiles), ncol(profiles) > 0L, nrow(records) > 0L,
    nrow(profiles) == nrow(records), !anyDuplicated(records$row_id),
    all(records$allele_group %in% c("REF", "ALT")))
  valid <- records$acf_valid
  clustered <- records$status == "clustered"
  stopifnot(!anyNA(valid), !any(clustered & !valid),
    all(is.finite(profiles[valid, , drop = FALSE])))
  records$allele_analysis_included <- TRUE
  # One physical molecule can appear in multiple sample libraries. Count once
  # only after verifying its biological annotation and ACF agree.
  duplicates <- split(seq_len(nrow(records)), records$RID)
  for (ix in duplicates[lengths(duplicates) > 1L]) {
    for (field in c("cell_line", "allele_group", "status", "cluster", "acf_valid"))
      if (length(unique(records[[field]][ix])) != 1L)
        stop("Conflicting duplicate molecule ", records$RID[ix[1L]], ": ", field)
    for (i in ix[-1L]) {
      x <- profiles[ix[1L], ]; y <- profiles[i, ]
      equal <- (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & abs(x-y) <= 1e-8 + 1e-5*abs(y))
      if (!all(equal)) stop("Conflicting ACFs for duplicate molecule ", records$RID[i])
    }
    records$allele_analysis_included[ix[-1L]] <- FALSE
  }
  peaks <- t(vapply(seq_len(nrow(profiles)), function(i) repeat_peak(profiles[i, ]), numeric(2L)))
  records$positive_local_peak_140_250_bp <- peaks[, 1L]
  records$peak_acf <- peaks[, 2L]
  complete <- ncol(profiles) >= 252L
  records$peak_band_complete <- complete
  clusters <- sort(unique(as.integer(records$cluster[clustered])))
  coverage <- curves <- features <- list()
  proportions <- data.frame(allele_group=character(), cluster=character(),
    n_in_cluster=integer(), n_clustered=integer(), fraction=numeric())
  for (allele in c("REF", "ALT")) {
    original <- records$allele_group == allele
    included <- original & records$allele_analysis_included
    vm <- included & valid
    cm <- included & clustered
    n_unique <- sum(included); n_valid <- sum(vm); n_clustered <- sum(cm)
    coverage[[allele]] <- data.frame(allele_group=allele, n_input=sum(original),
      n_unique=n_unique, n_duplicates_removed=sum(original)-n_unique,
      n_valid=n_valid, n_invalid=n_unique-n_valid, n_clustered=n_clustered,
      n_valid_not_clustered=n_valid-n_clustered, peak_band_complete=complete)
    average <- if (n_valid) colMeans(profiles[vm, , drop=FALSE]) else rep(NA_real_, ncol(profiles))
    curves[[allele]] <- data.frame(allele_group=allele,
      lag_bp=seq_len(ncol(profiles))-1L, mean_acf=average, n_reads=n_valid)
    for (cluster in clusters) {
      count <- sum(cm & records$cluster == as.character(cluster))
      proportions <- rbind(proportions, data.frame(allele_group=allele,
        cluster=as.character(cluster), n_in_cluster=count, n_clustered=n_clustered,
        fraction=if (n_clustered) count/n_clustered else NA_real_))
    }
    n_peak <- sum(!is.na(records$positive_local_peak_140_250_bp[vm]))
    features[[allele]] <- data.frame(allele_group=allele, n_valid=n_valid,
      n_with_peak=if (complete) n_peak else NA_real_,
      n_without_peak=if (complete) n_valid-n_peak else NA_real_,
      peak_fraction=if (complete && n_valid) n_peak/n_valid else NA_real_,
      median_peak_lag=if (n_peak) stats::median(records$positive_local_peak_140_250_bp[vm],na.rm=TRUE) else NA_real_,
      median_peak_height=if (n_peak) stats::median(records$peak_acf[vm],na.rm=TRUE) else NA_real_,
      peak_band_complete=complete)
  }
  tables <- list(allele_coverage=do.call(rbind,coverage),
    allele_cluster_proportions=proportions, allele_acf_summary=do.call(rbind,curves),
    allele_peak_features=do.call(rbind,features),
    allele_read_audit=records[,c("row_id", "RID", "sample_name", "cell_line",
      "allele_group", "allele_base", "acf_valid", "status", "cluster",
      "allele_analysis_included", "peak_band_complete", "positive_local_peak_140_250_bp", "peak_acf")])
  lapply(tables, function(table) {
    rownames(table) <- NULL
    data.frame(region_id=rep(region$region_id,nrow(table)),table,check.names=FALSE)
  })
}
