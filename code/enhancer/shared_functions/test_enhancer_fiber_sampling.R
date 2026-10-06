#!/usr/bin/env Rscript
# Run from the project root. Synthetic checks; no extraction or Slurm jobs.
suppressPackageStartupMessages(library(data.table))
source("code/enhancer/shared_functions/enhancer_fiber_sampling.R")
samples <- paste0("LPS_", c(0, 5, 10, 15))
make_fibers <- function(enhancer, counts, class_id = "inactive") {
  rbindlist(lapply(seq_along(samples), function(i) data.table(
    RID = paste(enhancer, samples[i], seq_len(counts[i]), sep = "::"),
    enhancer_id = rep(enhancer, counts[i]), enhancer_class = rep(class_id, counts[i]),
    sample_name = rep(samples[i], counts[i]))))
}
pool <- rbindlist(list(make_fibers("huge", rep(1000L, 4)),
  make_fibers("uneven", c(100L, 1L, 2L, 0L)),
  make_fibers("tiny", c(1L, 0L, 0L, 1L))))
warnings <- character()
fit <- withCallingHandlers(sample_capped_enhancer_fibers(pool), warning = function(w) {
  warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
})
selected <- pool[fit$selected]
stopifnot(length(warnings) == 1L, grepl("18 fibers", warnings),
  nrow(selected) == 18L, !anyDuplicated(selected$RID),
  max(selected[, .N, by = .(enhancer_id, sample_name)]$N) <= 3L,
  max(selected[, .N, by = enhancer_id]$N) <= 10L,
  identical(sort(selected[enhancer_id == "huge", .N, by = sample_name]$N), c(2L, 2L, 3L, 3L)),
  all(pool[enhancer_id == "tiny", RID] %in% selected$RID),
  identical(fit$diagnostics$n_fibers, c(4105L, 20L, 18L, 18L)),
  identical(fit$diagnostics$n_enhancers, rep(3L, 4)),
  fit$diagnostics[stage == "before_capping", fibers_per_enhancer_max] == 4000L,
  fit$diagnostics[stage == "after_capping", fibers_per_enhancer_max] == 10L)

# RIDs selected are reproducible even if source rows arrive in another order.
again <- suppressWarnings(sample_capped_enhancer_fibers(pool))
reversed <- pool[nrow(pool):1L]
reordered <- suppressWarnings(sample_capped_enhancer_fibers(reversed))
stopifnot(identical(fit$selected, again$selected),
  setequal(selected$RID, reversed$RID[reordered$selected]))
other_seed <- suppressWarnings(sample_capped_enhancer_fibers(pool, seed = 2L))
stopifnot(!setequal(selected$RID, pool$RID[other_seed$selected]))

# A partial final round is unbiased; exhausted timepoints are skipped.
balanced <- make_fibers("balanced", rep(8L, 4))
short_time <- make_fibers("short_time", c(1L, 3L, 3L, 3L))
for (s in 1:12) {
  even <- sample_capped_enhancer_fibers(balanced, cap_per_enhancer = 7L,
    n_fibers_per_class = 7L, seed = s)
  stopifnot(identical(sort(balanced[even$selected, .N, by = sample_name]$N), c(1L, 2L, 2L, 2L)))
  scarce <- sample_capped_enhancer_fibers(short_time, cap_per_enhancer = 8L,
    n_fibers_per_class = 8L, seed = s)
  stopifnot(identical(sort(short_time[scarce$selected, .N, by = sample_name]$N), c(1L, 2L, 2L, 3L)))
}

# Final class target is sampled after the caps and can only reduce contributions.
large_pool <- rbindlist(lapply(1:120, function(i) make_fibers(paste0("e", i), rep(5L, 4))))
small_target <- sample_capped_enhancer_fibers(large_pool, n_fibers_per_class = 100L)
stopifnot(small_target$info$n_capped == 1200L, length(small_target$selected) == 100L,
  max(large_pool[small_target$selected, .N, by = enhancer_id]$N) <= 10L,
  max(large_pool[small_target$selected, .N, by = .(enhancer_id, sample_name)]$N) <= 3L)

# Empty/single-fiber pools retain their actual sizes without relaxing caps.
for (n in 0:1) {
  few <- pool[seq_len(n)]
  result <- suppressWarnings(sample_capped_enhancer_fibers(few, class_id = "inactive"))
  stopifnot(length(result$selected) == n, all(result$diagnostics$n_fibers == n))
}
invalid <- try(sample_capped_enhancer_fibers(pool, cap_per_timepoint = 1.5), silent = TRUE)
stopifnot(inherits(invalid, "try-error"))

# Exact flag boundaries: 10% and 50 enhancers are allowed; >10% or <50 flag.
assignments <- rbindlist(list(
  data.table(cluster = "cluster1", enhancer_id = c(rep("dominant", 10), paste0("a", 1:90))),
  data.table(cluster = "cluster2", enhancer_id = c(rep("dominant", 11), paste0("b", 1:89))),
  data.table(cluster = "cluster3", enhancer_id = paste0("c", 1:49)),
  data.table(cluster = "cluster4", enhancer_id = paste0("d", 1:50))))
assignments[, sample_name := rep(samples, length.out = .N)]
diagnostics <- enhancer_cluster_diagnostics(assignments, "inactive", samples)
stopifnot(identical(diagnostics$flagged, c(FALSE, TRUE, TRUE, FALSE)),
  identical(diagnostics$n_enhancers, c(91L, 90L, 49L, 50L)),
  all(rowSums(as.matrix(diagnostics[, paste0("n_fibers_", samples), with = FALSE])) == diagnostics$n_fibers))
empty_diagnostics <- enhancer_cluster_diagnostics(assignments[0L], "inactive", samples)
stopifnot(nrow(empty_diagnostics) == 0L, "flagged" %in% names(empty_diagnostics))
assignments[, enhancer_class := rep(c("active", "inactive"), length.out = .N)]
pooled_diagnostics <- enhancer_cluster_diagnostics(assignments, "pooled", samples)
stopifnot(all(pooled_diagnostics$n_fibers_active + pooled_diagnostics$n_fibers_inactive ==
  pooled_diagnostics$n_fibers))

# The two methods must load the same ordered rows and preserve class labels.
source("code/enhancer/shared_functions/enhancer_shared_sampling.R")
active_pool <- copy(pool)
active_pool[, `:=`(RID = paste0("active_", RID), enhancer_id = paste0("active_", enhancer_id),
  enhancer_class = "active")]
active_fit <- suppressWarnings(sample_capped_enhancer_fibers(active_pool))
metadata <- rbind(active_pool[active_fit$selected], selected)
sample_info <- rbind(as.data.table(active_fit$info)[, enhancer_class := "active"],
  as.data.table(fit$info)[, enhancer_class := "inactive"])
shared <- list(mat = Matrix::sparseMatrix(i = seq_len(nrow(metadata)), j = rep(501L, nrow(metadata)),
    x = rep(1, nrow(metadata)), dims = c(nrow(metadata), 1000L),
    dimnames = list(metadata$RID, as.character(-500:499))), metadata = metadata,
  footprints = data.table(RID = character()), qc = data.table(), fiber_sampling = sample_info,
  sampling_diagnostics = rbind(active_fit$diagnostics, fit$diagnostics),
  inputs = list(samples = samples, regions = unique(metadata[, .(enhancer_id, enhancer_class)])),
  parameters = list(cap_per_timepoint = 3L, cap_per_enhancer = 10L, n_fibers_per_class = 10000L,
    sampling_seed = 1L, sampling_method = "capped_enhancer_timepoint_round_robin",
    position_bp = -500:499, signal = "raw_binary_m6a"), sample_id = "synthetic-shared")
shared_path <- tempfile(fileext = ".rds")
atomic_save_enhancer_rds(shared, shared_path)
manhattan_input <- load_shared_enhancer_fibers(shared_path)
acf_input <- load_shared_enhancer_fibers(shared_path)
stopifnot(identical(manhattan_input$md5, acf_input$md5),
  identical(manhattan_input$mat, acf_input$mat),
  identical(manhattan_input$metadata, acf_input$metadata),
  setequal(acf_input$metadata$enhancer_class, c("active", "inactive")))
reordered_shared <- shared
reordered_shared$metadata <- shared$metadata[nrow(metadata):1L]
stopifnot(inherits(try(validate_shared_enhancer_fibers(reordered_shared), silent = TRUE), "try-error"))
invalid_caps <- shared
invalid_caps$parameters$cap_per_enhancer <- 1L
stopifnot(inherits(try(validate_shared_enhancer_fibers(invalid_caps), silent = TRUE), "try-error"))
unlink(shared_path)
cat("PASS: capped sampling, balance, reproducibility, diagnostics and shared input identity\n")
