"""mean-centered, variance-normalized linear autocorrelation."""

import numpy as np


def autocorrelations(binary, n_features=None, positions=None, window=None):
    """Return all nonnegative lags of the linear autocorrelation function for each read
    ACF(k) = sum((x[t]-mean(x))*(x[t+k]-mean(x))) / (N * var(x))
    positions/window: optional 1-bp column coordinates and inclusive analysis bounds.
    n_features: number of nonnegative lags to return, default is all features (read width)
    """
    x = np.asarray(binary)
    # define input matrix to be 2D, with shape (n_reads, read_width)
    if x.ndim != 2 or x.shape[1] < 2:
        raise ValueError("Expected reads x consecutive genomic bases")
    # Optional closed genomic/strand-relative bounds; defaults retain the full input.
    if window is not None:
        if positions is None:
            raise ValueError("Window selection requires a position for each input column")
        positions = np.asarray(positions)
        bounds = np.asarray(window)
        if (positions.shape != (x.shape[1],) or bounds.shape != (2,) or
                not np.all(np.diff(positions) == 1) or bounds[0] > bounds[1] or
                not np.isin(bounds, positions).all()):
            raise ValueError("Expected consecutive 1-bp positions and covered inclusive window bounds")
        x = x[:, (positions >= bounds[0]) & (positions <= bounds[1])]
        if x.shape[1] < 2:
            raise ValueError("Autocorrelation requires at least two bases")
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
