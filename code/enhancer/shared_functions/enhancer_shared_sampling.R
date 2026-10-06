# Shared immutable input contract for pooled Manhattan and ACF analyses.

atomic_save_enhancer_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- tempfile(".enhancer_shared_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(partial), add = TRUE)
  saveRDS(object, partial)
  if (!file.rename(partial, path)) stop("Could not publish results: ", path)
  invisible(path)
}

validate_shared_enhancer_fibers <- function(shared) {
  stopifnot(all(c("mat", "metadata", "footprints", "qc", "fiber_sampling",
    "sampling_diagnostics", "inputs", "parameters", "sample_id") %in% names(shared)))
  a <- data.table::as.data.table(shared$metadata)
  p <- shared$parameters
  stopifnot(inherits(shared$mat, "sparseMatrix"), ncol(shared$mat) == 1000L,
    identical(colnames(shared$mat), as.character(-500:499)),
    identical(as.character(rownames(shared$mat)), as.character(a$RID)),
    nrow(shared$mat) == nrow(a), !anyNA(shared$mat), all(shared$mat@x == 1),
    !anyDuplicated(a$RID), !anyNA(a[, .(RID, enhancer_id, enhancer_class, sample_name)]),
    all(a$enhancer_class %in% c("active", "inactive")),
    all(a$sample_name %in% shared$inputs$samples),
    all(a$enhancer_id %in% shared$inputs$regions$enhancer_id),
    all(a[, data.table::uniqueN(enhancer_class), by = enhancer_id]$V1 == 1L),
    identical(p$position_bp, -500:499), identical(p$signal, "raw_binary_m6a"),
    identical(p$sampling_method, "capped_enhancer_timepoint_round_robin"),
    length(shared$sample_id) == 1L, !is.na(shared$sample_id), nzchar(shared$sample_id),
    all(a[, .N, by = .(enhancer_id, sample_name)]$N <= p$cap_per_timepoint),
    all(a[, .N, by = enhancer_id]$N <= p$cap_per_enhancer))
  counts <- a[, .N, by = enhancer_class]
  info <- data.table::as.data.table(shared$fiber_sampling)
  actual <- counts$N[match(info$enhancer_class, counts$enhancer_class)]
  actual[is.na(actual)] <- 0L
  stopifnot(nrow(info) == 2L, setequal(info$enhancer_class, c("active", "inactive")),
    all(info$n_selected == actual),
    all(info$n_selected == pmin(info$n_capped, p$n_fibers_per_class)),
    all(info$cap_per_timepoint == p$cap_per_timepoint),
    all(info$cap_per_enhancer == p$cap_per_enhancer),
    all(info$n_fibers_per_class == p$n_fibers_per_class),
    all(info$seed == p$sampling_seed), all(info$method == p$sampling_method))
  stopifnot(all(shared$footprints$RID %in% a$RID))
  invisible(shared)
}

load_shared_enhancer_fibers <- function(path) {
  if (!file.exists(path)) stop("Shared fiber sample is missing; run the preparation job first: ", path)
  before <- unname(tools::md5sum(path))
  shared <- readRDS(path)
  after <- unname(tools::md5sum(path))
  if (is.na(before) || !identical(before, after)) stop("Shared fiber sample changed while loading: ", path)
  validate_shared_enhancer_fibers(shared)
  shared$md5 <- before
  shared$path <- normalizePath(path)
  shared
}
