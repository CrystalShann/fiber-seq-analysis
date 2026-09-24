"""mean-centered, variance-normalized linear autocorrelation."""

import numpy as np


def autocorrelations(binary, n_features=None):
    """Return all nonnegative lags of the linear autocorrelation function for each read
    ACF(k) = sum((x[t]-mean(x))*(x[t+k]-mean(x))) / (N * var(x))
    n_features: number of nonnegative lags to return, default is all features (read width)
    """
    x = np.asarray(binary)
    # define input matrix to be 2D, with shape (n_reads, read_width)
    if x.ndim != 2 or x.shape[1] < 2:
        raise ValueError("Expected reads x consecutive genomic bases")
    # define the number of features to return, default is all features (read width)
    n_features = x.shape[1] if n_features is None else n_features
    # check for binary matrix and valid n_features
    if not isinstance(n_features, (int, np.integer)) or not 1 <= n_features <= x.shape[1]:
        raise ValueError("n_features must be an integer in [1, read width]")
    if not np.isfinite(x).all() or not np.isin(x, [0, 1]).all():
        raise ValueError("Input must contain only fully covered binary 0/1 calls")
    # subtract the mean methylation level for every each separately
    centered = x.astype(np.float64) - x.mean(axis=1, keepdims=True)
    # calculate the denominator for normalization, which is N * var(x) for each read
    # square the centered values and sum across the read width 
    denominator = np.einsum("ij,ij->i", centered, centered)
    valid = denominator > 0
    profiles = np.full((len(x), n_features), np.nan)
    # compute the autocorrelation for each read with nonzero variance
    for i in np.flatnonzero(valid):
        full = np.correlate(centered[i], centered[i], mode="full")
        profiles[i] = full[x.shape[1] - 1:x.shape[1] - 1 + n_features] / denominator[i]
    return profiles, valid
