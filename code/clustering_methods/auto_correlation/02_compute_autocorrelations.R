# Mean-centered, variance-normalized linear ACF; columns are consecutive bases.
# Lag k uses sum(x_centered[t] * x_centered[t+k]) / sum(x_centered^2).
autocorrelations <- function(binary, n_features = NULL, positions = NULL, window = NULL) {
  x <- as.matrix(binary)
  if (length(dim(binary)) != 2L || ncol(x) < 2L)
    stop("Expected reads x consecutive genomic bases")
  if (!is.null(window)) {
    if (is.null(positions)) stop("Window selection requires a position for each input column")
    if (length(positions) != ncol(x) || length(window) != 2L ||
        anyNA(c(positions, window)) || !all(diff(positions) == 1) ||
        window[1L] > window[2L] || !all(window %in% positions))
      stop("Expected consecutive 1-bp positions and covered inclusive window bounds")
    x <- x[, positions >= window[1L] & positions <= window[2L], drop = FALSE]
    if (ncol(x) < 2L) stop("Autocorrelation requires at least two bases")
  }
  if (is.null(n_features)) n_features <- ncol(x)
  if (length(n_features) != 1L || !is.numeric(n_features) || !is.finite(n_features) ||
      n_features != floor(n_features) || n_features < 1L || n_features > ncol(x))
    stop("n_features must be an integer in [1, read width]")
  if (!is.numeric(x) && !is.logical(x)) stop("Input must contain only fully covered binary 0/1 calls")
  if (any(!is.finite(x)) || !all(x %in% c(0, 1)))
    stop("Input must contain only fully covered binary 0/1 calls")
  valid <- unname(rowSums(x) > 0 & rowSums(x) < ncol(x))
  profiles <- matrix(NA_real_, nrow(x), n_features)
  for (i in which(valid)) {
    # stats::acf uses the fixed lag-zero denominator, without overlap correction.
    profiles[i, ] <- as.numeric(stats::acf(x[i, ], lag.max = n_features - 1L,
      type = "correlation", plot = FALSE, demean = TRUE)$acf)
  }
  list(profiles = profiles, valid = valid)
}
