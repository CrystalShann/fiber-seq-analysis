"""Memory-bounded Manhattan neighbors for pooled, binary per-base m6A rows.

The exact option matches R order(distance), including ties by input row order.
The approximate option uses NNDescent candidates; every retained distance is
recomputed as exact L1 on the original binary features. Neither method drops
or merges identical rows or samples observations.
"""

import gc
import numpy as np
from scipy import sparse


def _binary_csr(binary):
    mat = sparse.csr_matrix(binary, dtype=np.float32, copy=True)
    mat.sum_duplicates()
    mat.eliminate_zeros()
    mat.sort_indices()
    if mat.ndim != 2 or mat.shape[1] < 1:
        raise ValueError("binary must be a two-dimensional matrix with features")
    if not np.all(np.isfinite(mat.data)) or np.any(mat.data != 1):
        raise ValueError("binary must contain only finite 0/1 values")
    return mat


def _exact_rows(mat, rows, k_eff, row_sums, working_memory_mb, batch_size):
    # Sparse multiplication can briefly hold both CSR and dense output. Reserve
    # 16 bytes per pair, rather than budgeting only the final float32 array.
    n = mat.shape[0]
    budget = int(working_memory_mb * 1024**2)
    block_size = min(batch_size, budget // (16 * n))
    if block_size < 1:
        raise ValueError("working_memory_mb is too small for one distance row")
    neighbors = np.empty((len(rows), k_eff), dtype=np.int64)
    distances = np.empty((len(rows), k_eff), dtype=np.float32)
    # Convert once: CSR @ CSC would reconvert the full right-hand matrix
    # inside every sparse multiplication.
    transpose = mat.transpose().tocsr()
    for left in range(0, len(rows), block_size):
        block_rows = rows[left:left + block_size]
        overlap = (mat[block_rows] @ transpose).toarray()
        overlap *= -2
        overlap += row_sums[block_rows, None]
        overlap += row_sums[None, :]
        overlap[np.arange(len(block_rows)), block_rows] = np.inf
        for local, row_dist in enumerate(overlap):
            cutoff = np.partition(row_dist, k_eff - 1)[k_eff - 1]
            smaller = np.flatnonzero(row_dist < cutoff)
            equal = np.flatnonzero(row_dist == cutoff)[:k_eff - len(smaller)]
            chosen = np.concatenate((smaller, equal))
            chosen = chosen[np.lexsort((chosen, row_dist[chosen]))]
            neighbors[left + local] = chosen
            distances[left + local] = row_dist[chosen]
        del overlap
    return neighbors, distances, block_size


def build_manhattan_knn(binary, k_neighbors=10, method="approximate", seed=1,
                        n_jobs=4, working_memory_mb=256, batch_size=256):
    """Return 0-based non-self neighbor indices and exact binary L1 distances.

    ``method='exact'`` evaluates all pairs in bounded blocks (quadratic time).
    ``method='approximate'`` uses sparse Manhattan NNDescent (bounded candidate
    lists) and recomputes distances on its candidates. Incomplete candidate
    rows are repaired with an exact search and counted in the metadata.
    ``working_memory_mb`` bounds exact-distance temporary blocks; it is not
    a total-memory limit for the input, neighbor graph, or NNDescent index.
    """
    mat = _binary_csr(binary)
    n, n_features = mat.shape
    if n < 3:
        raise ValueError("At least three observations are required")
    if int(k_neighbors) != k_neighbors or k_neighbors < 1:
        raise ValueError("k_neighbors must be a positive integer")
    if int(batch_size) != batch_size or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")
    if int(n_jobs) != n_jobs or n_jobs < 1:
        raise ValueError("n_jobs must be a positive integer")
    if not np.isfinite(working_memory_mb) or working_memory_mb <= 0:
        raise ValueError("working_memory_mb must be finite and positive")
    if method not in ("exact", "approximate"):
        raise ValueError("method must be 'exact' or 'approximate'")
    k_eff = min(int(k_neighbors), n - 1)
    batch_size = int(batch_size)
    row_sums = np.asarray(mat.sum(axis=1)).ravel()
    metadata = dict(n_observations=n, n_features=n_features, seed=int(seed),
                    requested_k=int(k_neighbors), working_memory_mb=float(working_memory_mb),
                    batch_size=batch_size, n_jobs=int(n_jobs), repaired_rows=0)
    if method == "exact":
        neighbors, distances, exact_block_size = _exact_rows(
            mat, np.arange(n), k_eff, row_sums, working_memory_mb, batch_size)
        metadata["exact_block_size"] = exact_block_size
    else:
        import pynndescent
        candidate_k = min(n, max(30, k_eff + 1))
        metadata.update(candidate_k=candidate_k, pynndescent_version=pynndescent.__version__)
        index = pynndescent.NNDescent(
            mat, metric="manhattan", n_neighbors=candidate_k,
            random_state=int(seed), n_jobs=int(n_jobs), low_memory=True,
            compressed=False)
        candidates, _ = index.neighbor_graph
        del index
        gc.collect()
        neighbors = np.empty((n, k_eff), dtype=np.int64)
        distances = np.empty((n, k_eff), dtype=np.float32)
        repair = []
        # Recalculate all candidate distances in small sparse batches. This also
        # avoids trusting any approximate or uninitialized distance entries.
        for left in range(0, n, batch_size):
            right = min(n, left + batch_size)
            candidate_block = candidates[left:right]
            valid = (candidate_block >= 0) & (candidate_block < n)
            safe_candidates = np.where(valid, candidate_block, 0)
            source_rows = np.repeat(np.arange(left, right), candidate_k)
            target_rows = safe_candidates.ravel()
            overlaps = np.asarray(mat[source_rows].multiply(mat[target_rows]).sum(axis=1)).ravel()
            candidate_dist = (row_sums[source_rows] + row_sums[target_rows] - 2 * overlaps)
            candidate_dist = candidate_dist.reshape(candidate_block.shape)
            for local, observation in enumerate(range(left, right)):
                keep = valid[local] & (candidate_block[local] != observation)
                ids, unique_indices = np.unique(candidate_block[local, keep], return_index=True)
                d = candidate_dist[local, keep][unique_indices]
                if len(ids) < k_eff:
                    repair.append(observation)
                    continue
                order = np.lexsort((ids, d))[:k_eff]
                neighbors[observation] = ids[order]
                distances[observation] = d[order]
        del candidates
        if repair:
            repaired_neighbors, repaired_distances, exact_block_size = _exact_rows(
                mat, np.asarray(repair), k_eff, row_sums, working_memory_mb, batch_size)
            neighbors[repair] = repaired_neighbors
            distances[repair] = repaired_distances
            metadata.update(repaired_rows=len(repair), exact_block_size=exact_block_size)
    if (np.any(neighbors < 0) or np.any(neighbors >= n)
            or np.any(neighbors == np.arange(n)[:, None])
            or not np.all(np.isfinite(distances)) or np.any(distances < 0)):
        raise RuntimeError("Invalid neighbor indices or Manhattan distances")
    if k_eff > 1 and np.any(np.diff(np.sort(neighbors, axis=1), axis=1) == 0):
        raise RuntimeError("Duplicate neighbor IDs remain")
    metadata["mean_knn_distance"] = float(distances.mean(dtype=np.float64))
    return dict(indices=neighbors, distances=distances, method=method,
                k_eff=k_eff, metadata=metadata)
