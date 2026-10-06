# Common site discovery and co-binding analysis for active/inactive enhancers.
# Coordinates are 0-based half-open; discovery and tests are exploratory.
load_site_pair_fibers <- function(loci, samples) {
  stopifnot(nrow(loci) > 0L, data.table::uniqueN(loci$chr) == 1L,
            !anyDuplicated(loci$enhancer_id), !anyDuplicated(samples$sample))
  locus_ranges <- bed_ranges(loci)
  query_ranges <- GenomicRanges::reduce(locus_ranges)
  chromosome_name <- loci$chr[[1L]]
  coverage_parts <- midpoint_parts <- vector("list", nrow(samples))
  empty_coverage <- data.table(enhancer_id = character(), sample = character(),
                               read_id = character())
  empty_midpoints <- data.table(enhancer_id = character(), sample = character(),
                                read_id = character(), midpoint0 = numeric())
  for (sample_index in seq_len(nrow(samples))) {
    sample_name <- samples$sample[[sample_index]]
    methylation_file <- file.path(ft_result_dir, sample_name, "extracted_results",
      "m6a_by_chr", paste0(sample_name, ".ft_extracted_m6a.", chromosome_name, ".bed.gz"))
    raw_fibers <- tabix_bed(methylation_file, query_ranges)
    if (!nrow(raw_fibers)) next
    fiber_spans <- unique(data.table(
      chr = as.character(raw_fibers[[1L]]), start0 = as.integer(raw_fibers[[2L]]),
      end0 = as.integer(raw_fibers[[3L]]), read_id = as.character(raw_fibers[[4L]])))
    fiber_spans[, alignment_width := end0 - start0]
    data.table::setorderv(fiber_spans, c("read_id", "alignment_width"), c(1L, -1L))
    fiber_spans <- unique(fiber_spans, by = "read_id")
    span_hits <- GenomicRanges::findOverlaps(locus_ranges, bed_ranges(fiber_spans),
      type = "within", ignore.strand = TRUE)
    spanning_fibers <- unique(data.table(
      enhancer_id = loci$enhancer_id[S4Vectors::queryHits(span_hits)],
      read_id = fiber_spans$read_id[S4Vectors::subjectHits(span_hits)]))
    if (!nrow(spanning_fibers)) next
    spanning_fibers[, sample := sample_name]
    coverage_parts[[sample_index]] <- spanning_fibers

    short_calls <- expand_hmm_blocks(tabix_bed(
      hmm_file(sample_name, chromosome_name, "tf"), query_ranges))
    short_calls <- unique(short_calls[length > 0L & length < 60L],
      by = c("chr", "read_id", "start0", "end0"))
    if (!nrow(short_calls)) next
    short_calls[, midpoint0 := floor((as.numeric(start0) + end0) / 2)]
    midpoint_ranges <- GenomicRanges::GRanges(short_calls$chr,
      IRanges::IRanges(short_calls$midpoint0 + 1L, short_calls$midpoint0 + 1L))
    midpoint_hits <- GenomicRanges::findOverlaps(midpoint_ranges, locus_ranges,
      ignore.strand = TRUE)
    assigned_calls <- data.table(
      enhancer_id = loci$enhancer_id[S4Vectors::subjectHits(midpoint_hits)],
      read_id = short_calls$read_id[S4Vectors::queryHits(midpoint_hits)],
      midpoint0 = short_calls$midpoint0[S4Vectors::queryHits(midpoint_hits)])
    # Preserve distinct footprint observations, but exclude partial-span fibers.
    midpoint_parts[[sample_index]] <- merge(assigned_calls, spanning_fibers,
      by = c("enhancer_id", "read_id"), all = FALSE, sort = FALSE)
  }
  list(
    coverage = rbindlist(c(list(empty_coverage), coverage_parts), use.names = TRUE),
    midpoints = rbindlist(c(list(empty_midpoints), midpoint_parts), use.names = TRUE))
}

discover_footprint_sites <- function(midpoints, window_bp, minimum_observations) {
  ordered_midpoints <- sort(midpoints)
  if (!length(ordered_midpoints))
    return(data.table(site_index = integer(), site_center = numeric(),
                      n_observations = integer()))
  cluster_index <- integer(length(ordered_midpoints))
  current_cluster <- 1L
  first_midpoint <- ordered_midpoints[[1L]]
  for (position_index in seq_along(ordered_midpoints)) {
    if (ordered_midpoints[[position_index]] - first_midpoint > window_bp) {
      current_cluster <- current_cluster + 1L
      first_midpoint <- ordered_midpoints[[position_index]]
    }
    cluster_index[[position_index]] <- current_cluster
  }
  # Compare with the FIRST midpoint, not the previous midpoint or moving median.
  site_members <- data.table(site_index = cluster_index, midpoint0 = ordered_midpoints)
  site_members[, .(site_center = stats::median(midpoint0), n_observations = .N),
    by = site_index][n_observations >= minimum_observations]
}

empty_site_pair_results <- function() {
  data.table(enhancer_id = character(), site_1_center = numeric(),
    site_2_center = numeric(), site_distance_bp = numeric(), sample = character(),
    time_min = integer(), n_tf_tf = integer(), n_tf_naked = integer(),
    n_naked_tf = integer(), n_naked_naked = integer(), total = integer(),
    obs_exp = numeric(), p_value = numeric(), q_value = numeric(), status = character())
}

count_site_pair_states <- function(enhancer_key, sites, coverage, midpoints,
                                   samples, window_bp, minimum_pair_reads) {
  if (nrow(sites) < 2L) return(empty_site_pair_results())
  site_combinations <- utils::combn(seq_len(nrow(sites)), 2L)
  first_sites <- site_combinations[1L, ]
  second_sites <- site_combinations[2L, ]
  sample_results <- vector("list", nrow(samples))
  for (sample_index in seq_len(nrow(samples))) {
    sample_name <- samples$sample[[sample_index]]
    spanning_ids <- unique(coverage[sample == sample_name, read_id])
    sample_calls <- midpoints[sample == sample_name]
    bound_by_site <- matrix(FALSE, nrow = length(spanning_ids), ncol = nrow(sites))
    for (site_index in seq_len(nrow(sites))) {
      center_position <- sites$site_center[[site_index]]
      bound_ids <- sample_calls[abs(midpoint0 - center_position) <= window_bp / 2,
                               read_id]
      bound_by_site[, site_index] <- spanning_ids %in% bound_ids
    }
    # Boolean membership enforces the TF-any priority and prevents multiple
    # footprints on one fiber from inflating a site's occupancy count.
    site_totals <- colSums(bound_by_site)
    joint_totals <- crossprod(1L * bound_by_site)
    both_bound <- as.integer(joint_totals[cbind(first_sites, second_sites)])
    first_only <- as.integer(site_totals[first_sites] - both_bound)
    second_only <- as.integer(site_totals[second_sites] - both_bound)
    pair_counts <- data.table(
      enhancer_id = enhancer_key,
      site_1_center = sites$site_center[first_sites],
      site_2_center = sites$site_center[second_sites],
      site_distance_bp = sites$site_center[second_sites] - sites$site_center[first_sites],
      sample = sample_name, time_min = samples$time_min[[sample_index]],
      n_tf_tf = both_bound, n_tf_naked = first_only, n_naked_tf = second_only,
      n_naked_naked = length(spanning_ids) - both_bound - first_only - second_only,
      total = length(spanning_ids), obs_exp = NA_real_, p_value = NA_real_,
      q_value = NA_real_, status = "not_testable")
    stopifnot(all(pair_counts$n_naked_naked >= 0L),
      all(pair_counts[, n_tf_tf + n_tf_naked + n_naked_tf + n_naked_naked == total]))
    pair_counts[total > 0L & n_tf_tf + n_tf_naked > 0L & n_tf_tf + n_naked_tf > 0L,
      obs_exp := (n_tf_tf / total) /
        (((n_tf_tf + n_tf_naked) / total) * ((n_tf_tf + n_naked_tf) / total))]
    for (pair_index in seq_len(nrow(pair_counts))) {
      pair_row <- pair_counts[pair_index]
      contingency <- matrix(c(pair_row$n_tf_tf, pair_row$n_tf_naked,
                               pair_row$n_naked_tf, pair_row$n_naked_naked), nrow = 2L)
      if (pair_row$total >= minimum_pair_reads &&
          all(rowSums(contingency) > 0L) && all(colSums(contingency) > 0L)) {
        pair_p <- cobinding_fisher_p(contingency)
        data.table::set(pair_counts, pair_index, "p_value", pair_p)
        data.table::set(pair_counts, pair_index, "status", "tested")
      }
    }
    sample_results[[sample_index]] <- pair_counts
  }
  rbindlist(c(list(empty_site_pair_results()), sample_results), use.names = TRUE)
}

calculate_enhancer_site_pairs <- function(analysis_enhancers, sample_table,
    footprint_counts, footprint_result, eligible_enhancers,
    site_window_bp = 50L, min_site_observations = 5L, min_pair_reads = 10L,
    cobinding_alpha = 0.01, cobinding_significance = "q_value") {
  stopifnot(site_window_bp>0, min_site_observations>=1L, min_pair_reads>=1L,
            cobinding_alpha > 0, cobinding_alpha < 1,
            cobinding_significance %in% c("p_value", "q_value"))
pair_enhancers <- copy(analysis_enhancers[
  enhancer_class %in% c("active", "inactive") & enhancer_id %in% eligible_enhancers])
stopifnot(!anyDuplicated(pair_enhancers$enhancer_id),
          !anyDuplicated(sample_table$sample), !anyDuplicated(sample_table$time_min))
pair_reference_counts <- copy(footprint_counts)[, .(enhancer_id, sample, n_reads)]
stopifnot(!anyDuplicated(pair_reference_counts[, .(enhancer_id, sample)]))
# Analysis 3 chose read alignments in these exact chromosome batches, before
# either the active-class or minimum-coverage filter. Preserve their membership
# and order so a competing alignment cannot appear/disappear through rebatching.
pair_input_enhancers <- as.data.table(footprint_result$inputs$enhancers)
stopifnot(nrow(pair_input_enhancers) > 0L,
          !anyDuplicated(pair_input_enhancers$enhancer_id),
          all(pair_enhancers$enhancer_id %in% pair_input_enhancers$enhancer_id))
site_pair_parts <- list()
site_inventory <- list()
for (chromosome_name in unique(pair_enhancers$chr)) {
  chromosome_loci <- pair_input_enhancers[chr == chromosome_name]
  for (batch_start in seq.int(1L, nrow(chromosome_loci), by = 100L)) {
    batch_loci <- chromosome_loci[batch_start:min(nrow(chromosome_loci), batch_start + 99L)]
    batch_targets <- batch_loci[enhancer_id %in% pair_enhancers$enhancer_id]
    if (!nrow(batch_targets)) next
    message("Site-pair batch: ", chromosome_name, " / ", batch_start, " at ", Sys.time())
    batch_fibers <- load_site_pair_fibers(batch_loci, sample_table)
    batch_fibers$coverage <- batch_fibers$coverage[
      enhancer_id %in% batch_targets$enhancer_id]
    batch_fibers$midpoints <- batch_fibers$midpoints[
      enhancer_id %in% batch_targets$enhancer_id]
    span_totals <- batch_fibers$coverage[, .(loaded_reads = .N),
                                        by = .(enhancer_id, sample)]
    denominator_audit <- merge(CJ(enhancer_id = batch_targets$enhancer_id,
                                  sample = sample_table$sample),
      pair_reference_counts, by = c("enhancer_id", "sample"), all.x = TRUE, sort = FALSE)
    denominator_audit <- merge(denominator_audit, span_totals,
      by = c("enhancer_id", "sample"), all.x = TRUE, sort = FALSE)
    denominator_audit[is.na(loaded_reads), loaded_reads := 0L]
    mismatched_denominators <- denominator_audit[
      is.na(n_reads) | n_reads != loaded_reads]
    if (nrow(mismatched_denominators))
      stop("Site-pair spanning-read counts differ from Analysis 3. ",
           "The original enhancer batches were reused; check that the saved ",
           "footprint results and current m6A inputs correspond. ",
           "First mismatches (n_reads = Analysis 3; loaded_reads = site-pair):\n",
           paste(capture.output(print(utils::head(mismatched_denominators, 6L))),
                 collapse = "\n"))
    for (enhancer_key in batch_targets$enhancer_id) {
      locus_coverage <- batch_fibers$coverage[enhancer_id == enhancer_key]
      locus_midpoints <- batch_fibers$midpoints[enhancer_id == enhancer_key]
      retained_sites <- discover_footprint_sites(locus_midpoints$midpoint0,
        site_window_bp, min_site_observations)
      site_inventory[[length(site_inventory) + 1L]] <- data.table(
        enhancer_id = enhancer_key, n_sites = nrow(retained_sites),
        n_possible_pairs = choose(nrow(retained_sites), 2L))
      if (nrow(retained_sites) < 2L) next
      site_pair_parts[[length(site_pair_parts) + 1L]] <- count_site_pair_states(
        enhancer_key, retained_sites, locus_coverage, locus_midpoints, sample_table,
        site_window_bp, min_pair_reads)
    }
  }
}
site_pair_cooccupancy <- rbindlist(c(list(empty_site_pair_results()), site_pair_parts),
                                  use.names = TRUE)

  site_pair_cooccupancy[, enhancer_class := pair_enhancers$enhancer_class[
    match(enhancer_id, pair_enhancers$enhancer_id)]]
  pairs <- annotate_cobinding_pairs(site_pair_cooccupancy, cobinding_alpha,
                                    cobinding_significance)
  inventory <- rbindlist(c(list(data.table(enhancer_id=character(), n_sites=integer(),
    n_possible_pairs=numeric())), site_inventory))
  enhancers <- summarize_enhancer_pairs(pairs, inventory, pair_enhancers,
                                       footprint_counts, min_pair_reads)
  list(pairs=pairs, enhancers=enhancers, sites=inventory,
       parameters=list(site_window_bp=site_window_bp,
         min_site_observations=min_site_observations, min_pair_reads=min_pair_reads,
         cobinding_alpha=cobinding_alpha, cobinding_significance=cobinding_significance))
}

annotate_cobinding_pairs <- function(pairs, alpha=0.01, significance="q_value") {
  x <- copy(as.data.table(pairs))
  stopifnot(significance %in% c("p_value","q_value"), alpha > 0, alpha < 1)
  x[, q_value := NA_real_]
  x[status == "tested", q_value := p.adjust(p_value, method="BH"), by=time_min]
  x[, cobinding_class := "untestable"]
  x[status == "tested", cobinding_class := "no significant association"]
  x[status == "tested" & get(significance) < alpha & obs_exp > 1,
    cobinding_class := "enriched"]
  x[status == "tested" & get(significance) < alpha & obs_exp < 1,
    cobinding_class := "depleted"]
  x[, `:=`(occupancy_both=n_tf_tf/total, occupancy_site1_only=n_tf_naked/total,
    occupancy_site2_only=n_naked_tf/total, occupancy_neither=n_naked_naked/total)]
  setorder(x, enhancer_id, site_1_center, site_2_center, time_min)
  x
}

# Expected number of distinct pairs observed in a fixed-size, uniformly sampled
# subset of spanning fibers. Pair dependence does not affect this expectation:
# linearity of expectation applies without an independence assumption.
rarefied_pair_probability <- function(n_both, n_reads, depth) {
  stopifnot(length(depth)==1L, depth>=1L, depth%%1==0,
    all(n_reads>=depth), all(n_both>=0 & n_both<=n_reads))
  pmin(1, pmax(0, -expm1(lchoose(n_reads-n_both,depth)-lchoose(n_reads,depth))))
}

summarize_enhancer_pairs <- function(pairs, inventory, enhancers, counts, depth=10L) {
  metadata <- enhancers[,.(enhancer_id, enhancer_class, length_bp=end0-start0)]
  stopifnot(!anyDuplicated(metadata$enhancer_id), all(metadata$length_bp>0))
  base <- merge(counts[enhancer_id %in% metadata$enhancer_id,
    .(enhancer_id,sample,time_min,n_reads)], metadata, by="enhancer_id")
  base <- merge(base, inventory, by="enhancer_id", all.x=TRUE)
  stopifnot(!anyNA(base$n_sites), !anyNA(base$n_possible_pairs), all(base$n_reads>=depth))
  summary <- pairs[,.(n_cobound_pairs=sum(n_tf_tf>0),
    n_tested_pairs=sum(status=="tested"), n_enriched_pairs=sum(cobinding_class=="enriched"),
    n_depleted_pairs=sum(cobinding_class=="depleted"),
    n_unclassified_pairs=sum(cobinding_class=="no significant association"),
    cobinding_events=sum(as.numeric(n_tf_tf)),
    standardized_pairs=sum(rarefied_pair_probability(n_tf_tf,total,depth))),
    by=.(enhancer_id,sample)]
  base <- merge(base,summary,by=c("enhancer_id","sample"),all.x=TRUE)
  zero_cols <- c("n_cobound_pairs","n_tested_pairs","n_enriched_pairs",
    "n_depleted_pairs","n_unclassified_pairs","cobinding_events","standardized_pairs")
  for (col in zero_cols) set(base,which(is.na(base[[col]])),col,0)
  base[, `:=`(standardized_depth=depth, mean_pairs_per_fiber=cobinding_events/n_reads,
    fraction_pairs_cobound=fifelse(n_possible_pairs>0,n_cobound_pairs/n_possible_pairs,NA_real_))]
  stopifnot(all(base$n_cobound_pairs<=base$n_possible_pairs),
    all(base$standardized_pairs<=base$n_cobound_pairs+1e-8))
  setorder(base, enhancer_id,time_min)
  base
}

compare_enhancer_pair_counts <- function(summary) {
  metrics <- c("standardized_pairs","n_cobound_pairs","mean_pairs_per_fiber",
               "fraction_pairs_cobound")
  empty <- data.table(time_min=integer(), metric=character(), method=character(),
    n_active=integer(), n_inactive=integer(), median_active=numeric(),
    median_inactive=numeric(), estimate=numeric(), statistic=numeric(), p_value=numeric(),
    q_value=numeric(), status=character())
  rows <- list()
  for (tm in sort(unique(summary$time_min))) for (metric in metrics) {
    d <- copy(summary[time_min==tm & is.finite(get(metric))])
    a <- d[enhancer_class=="active",get(metric)]
    b <- d[enhancer_class=="inactive",get(metric)]
    for (method in c("wilcoxon_active_greater","adjusted_HC3_active_greater")) {
      row <- data.table(time_min=as.integer(tm),metric=metric,method=method,
        n_active=length(a),n_inactive=length(b),
        median_active=if(length(a)) median(a) else NA_real_,
        median_inactive=if(length(b)) median(b) else NA_real_,
        estimate=NA_real_,statistic=NA_real_,p_value=NA_real_,q_value=NA_real_,
        status="insufficient_data")
      if (length(a)>=2L && length(b)>=2L) {
        if (uniqueN(c(a,b))==1L) {
          row[,`:=`(estimate=0,statistic=0,p_value=1,status="constant_outcome")]
        } else if (method=="wilcoxon_active_greater") {
          fit <- wilcox.test(a,b,alternative="greater",exact=FALSE)
          row[,`:=`(estimate=median(a)-median(b),statistic=unname(fit$statistic),
                     p_value=fit$p.value,status="tested")]
        } else {
          d[,`:=`(response=log1p(get(metric)),
            active=as.integer(enhancer_class=="active"),
            log_reads=log1p(n_reads),log_length=log(length_bp),
            log_possible_pairs=log1p(n_possible_pairs))]
          covariates <- c("log_reads","log_length","log_possible_pairs")
          covariates <- covariates[vapply(covariates,function(col) uniqueN(d[[col]])>1L,logical(1))]
          fit <- lm(reformulate(c("active",covariates),response="response"),data=d)
          # A singular design does not identify a unique adjusted class effect.
          if (fit$rank < length(coef(fit))) {
            row[,status:="singular_design"]
          } else if (df.residual(fit)>1L && all(hatvalues(fit)<1-1e-8)) {
            variance <- sandwich::vcovHC(fit,type="HC3")
            se <- sqrt(variance["active","active"])
            effect <- unname(coef(fit)["active"])
            if (is.finite(se) && se>0) row[,`:=`(estimate=effect,
              statistic=effect/se,p_value=pt(effect/se,df=df.residual(fit),lower.tail=FALSE),
              status="tested")]
          }
        }
      }
      rows[[length(rows)+1L]] <- row
    }
  }
  out <- rbindlist(c(list(empty),rows))
  out[is.finite(p_value),q_value:=p.adjust(p_value,"BH"),by=.(metric,method)]
  out
}

# Repeated small contingency tables are common. Cache exact Fisher p-values;
# reset after 200,000 unique tables to bound memory during large full-genome runs.
cobinding_fisher_p <- local({
  cache <- new.env(hash=TRUE, parent=emptyenv())
  n_cached <- 0L
  function(contingency) {
    key <- paste(as.integer(contingency),collapse=":")
    if (exists(key,envir=cache,inherits=FALSE)) return(get(key,envir=cache))
    p <- stats::fisher.test(contingency,alternative="two.sided")$p.value
    if (n_cached>=200000L) {
      cache <<- new.env(hash=TRUE,parent=emptyenv()); n_cached <<- 0L
    }
    assign(key,p,envir=cache); n_cached <<- n_cached+1L
    p
  }
})
