# Sampling and contribution diagnostics; these do not transform m6A features.

validate_sampling_integer <- function(value, name, minimum = 1L) {
  if (length(value) != 1L || !is.finite(value) || value < minimum ||
      value > .Machine$integer.max || value != floor(value))
    stop(name, " must be an integer >= ", minimum)
  as.integer(value)
}

enhancer_sampling_diagnostics <- function(metadata, class_id, stage, samples) {
  x <- data.table::as.data.table(metadata)
  counts <- x[, .N, by = enhancer_id]$N
  out <- data.table::data.table(enhancer_class = class_id, stage = stage,
    n_fibers = nrow(x), n_enhancers = length(counts),
    fibers_per_enhancer_min = if (length(counts)) min(counts) else 0L,
    fibers_per_enhancer_median = if (length(counts)) stats::median(counts) else 0,
    fibers_per_enhancer_max = if (length(counts)) max(counts) else 0L)
  # Distribution is over contributing enhancers, not zero-coverage loci.
  for (sample_id in samples)
    out[, (paste0("n_fibers_", sample_id)) := sum(x$sample_name == sample_id)]
  out
}

sample_capped_enhancer_fibers <- function(metadata, cap_per_timepoint = 3L,
                                        cap_per_enhancer = 10L,
                                        n_fibers_per_class = 10000L, seed = 1L,
                                        samples = c("LPS_0", "LPS_5", "LPS_10", "LPS_15"),
                                        class_id = NULL) {
  cap_per_timepoint <- validate_sampling_integer(cap_per_timepoint, "cap_per_timepoint")
  cap_per_enhancer <- validate_sampling_integer(cap_per_enhancer, "cap_per_enhancer")
  n_fibers_per_class <- validate_sampling_integer(n_fibers_per_class, "n_fibers_per_class")
  seed <- validate_sampling_integer(seed, "seed", 0L)
  x <- data.table::copy(data.table::as.data.table(metadata))
  stopifnot(all(c("RID", "enhancer_id", "enhancer_class", "sample_name") %in% names(x)),
    !anyNA(x[, .(RID, enhancer_id, enhancer_class, sample_name)]), !anyDuplicated(x$RID),
    !anyDuplicated(samples), all(x$sample_name %in% samples))
  if (is.null(class_id)) {
    stopifnot(data.table::uniqueN(x$enhancer_class) == 1L)
    class_id <- as.character(x$enhancer_class[1L])
  }
  stopifnot(length(class_id) == 1L, !is.na(class_id), all(x$enhancer_class == class_id))
  x[, input_row := seq_len(.N)]
  # Canonical ordering makes the same seed select the same RIDs even if the
  # extractor changes row order. Return indices in the original matrix order.
  data.table::setorder(x, enhancer_id, sample_name, RID)
  set.seed(seed)
  time_capped <- x[, .SD[if (.N > cap_per_timepoint)
    sample.int(.N, cap_per_timepoint) else seq_len(.N)],
    by = .(enhancer_id, sample_name)]
  # One random fiber from each available timepoint per round. Randomize ties
  # every round so a partially filled last round does not favor early times.
  time_capped[, sampling_round := sample.int(.N), by = .(enhancer_id, sample_name)]
  capped <- time_capped[, .SD[if (.N > cap_per_enhancer)
    order(sampling_round, stats::runif(.N))[seq_len(cap_per_enhancer)] else seq_len(.N)],
    by = enhancer_id]
  if (nrow(capped) < n_fibers_per_class)
    warning(class_id, ": capped pool has ", nrow(capped), " fibers; target is ",
      n_fibers_per_class, ". Using all available fibers without relaxing caps.", call. = FALSE)
  chosen <- if (nrow(capped) > n_fibers_per_class)
    sample.int(nrow(capped), n_fibers_per_class) else seq_len(nrow(capped))
  selected <- sort(capped$input_row[chosen])
  stopifnot(!anyDuplicated(selected), length(selected) <= n_fibers_per_class,
    all(capped[, .N, by = .(enhancer_id, sample_name)]$N <= cap_per_timepoint),
    all(capped[, .N, by = enhancer_id]$N <= cap_per_enhancer))
  diagnostics <- data.table::rbindlist(list(
    enhancer_sampling_diagnostics(x, class_id, "before_capping", samples),
    enhancer_sampling_diagnostics(time_capped, class_id, "after_timepoint_cap", samples),
    enhancer_sampling_diagnostics(capped, class_id, "after_capping", samples),
    enhancer_sampling_diagnostics(capped[chosen], class_id, "selected", samples)))
  list(selected = selected, diagnostics = diagnostics,
    info = list(n_eligible = nrow(x), n_after_timepoint_cap = nrow(time_capped),
      n_capped = nrow(capped), n_selected = length(selected),
      cap_per_timepoint = cap_per_timepoint, cap_per_enhancer = cap_per_enhancer,
      n_fibers_per_class = n_fibers_per_class, seed = seed,
      method = "capped_enhancer_timepoint_round_robin"))
}

enhancer_cluster_diagnostics <- function(assignments, class_id, samples) {
  a <- data.table::copy(data.table::as.data.table(assignments))
  a[, cluster := as.character(cluster)]
  a[is.na(cluster), cluster := "unclustered"]
  contributions <- a[, .N, by = .(cluster, enhancer_id)]
  out <- contributions[, .(n_fibers = sum(N), n_enhancers = .N,
    largest_enhancer_fraction = if (.N) max(N) / sum(N) else NA_real_), by = cluster]
  out[, enhancer_class := class_id]
  out[, `:=`(flag_dominant_enhancer = largest_enhancer_fraction > 0.10,
    flag_few_enhancers = n_enhancers < 50L)]
  out[, flagged := flag_dominant_enhancer | flag_few_enhancers]
  for (sample_id in samples) {
    counts <- a[sample_name == sample_id, .N, by = cluster]
    values <- counts$N[match(out$cluster, counts$cluster)]
    values[is.na(values)] <- 0L
    out[, (paste0("n_fibers_", sample_id)) := values]
  }
  if ("enhancer_class" %in% names(a)) {
    for (enhancer_group in c("active", "inactive")) {
      counts <- a[enhancer_class == enhancer_group, .N, by = cluster]
      values <- counts$N[match(out$cluster, counts$cluster)]
      values[is.na(values)] <- 0L
      out[, (paste0("n_fibers_", enhancer_group)) := values]
    }
  }
  data.table::setcolorder(out, c("enhancer_class", "cluster", setdiff(names(out), c("enhancer_class", "cluster"))))
  out[order(suppressWarnings(as.integer(sub("^cluster", "", cluster))), na.last = TRUE)]
}
