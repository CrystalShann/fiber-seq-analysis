"""Cluster ACF shapes """

import numpy as np

# Parameters: number of principal components, number of neighbors, Leiden resolution
def cluster_profiles(profiles, valid, *, n_pcs=50, n_neighbors=10, resolution=0.4, seed=0):
    import scanpy as sc

    if n_pcs < 2 or n_neighbors < 2 or resolution <= 0:
        raise ValueError("Require >=2 PCs, >=2 neighbors and positive Leiden resolution")
    profiles = np.asarray(profiles)
    valid = np.asarray(valid, dtype=bool)
    if profiles.ndim != 2 or valid.shape != (len(profiles),):
        raise ValueError("Profiles and validity mask are not aligned")
    if not np.isfinite(profiles[valid]).all():
        raise ValueError("Valid profiles must be finite")

    # prepare cluster labels, valid = True for profiles with sufficient reads and non-zero variance
    labels = np.full(len(profiles), "", dtype=object)
    status = np.where(valid, "insufficient_reads", "zero_variance").astype(object)
    # prepare UMAP coordinates and graph edges
    # graph edges: read_i, read_j, weight
    embedding = np.full((len(profiles), 2), np.nan)
    edges = np.empty((0, 3))
    info = {"method": "PCA -> correlation neighbors -> Leiden; UMAP for display only",
            "requested_pcs": n_pcs, "requested_neighbors": n_neighbors,
            "resolution": resolution, "seed": seed, "n_clustered": 0}
    indices = np.flatnonzero(valid)
    if len(indices) < 4:
        return labels, status, embedding, info, edges
    values = profiles[indices]
    if np.allclose(values, values[0], rtol=0, atol=1e-12):
        status[indices] = "identical_profiles"
        return labels, status, embedding, info, edges
    # Lag 0 is retained to match the reference feature set. It becomes a zero
    # column when PCA centers the features and thus does not drive clustering

    # decide how many PCs to use, based on the number of valid reads and the number of features
    # require at least 3 PCs to cluster
    actual_pcs = min(n_pcs, len(indices) - 1, profiles.shape[1] - 1)
    if actual_pcs < 2:
        raise ValueError("At least three lag features are required for clustering")
    adata = sc.AnnData(values.copy())

    # perform PCA, neighbors, Leiden clustering, and UMAP embedding 
    # zero_center=True is required for correlation distance to be equivalent to Euclidean distance on PCA scores
    sc.pp.pca(adata, n_comps=actual_pcs, zero_center=True, svd_solver="arpack", random_state=seed)
    # extract PCA cooridinates and filter out reads with non-finite PCA scores or zero variance
    pca = adata.obsm["X_pca"]
    usable = np.isfinite(pca).all(axis=1) & (np.std(pca, axis=1) > 1e-12)
    status[indices[~usable]] = "undefined_pca_correlation"
    indices = indices[usable]
    adata = adata[usable].copy()
    if len(indices) < 4:
        return labels, status, embedding, info, edges
    actual_neighbors = min(n_neighbors, len(indices) - 1)
    # construct a neighbor graph using correlation distance on PCA scores, then cluster with Leiden
    sc.pp.neighbors(adata, n_neighbors=actual_neighbors, metric="correlation",
                    use_rep="X_pca", n_pcs=actual_pcs, random_state=seed,
                    method="umap", transformer="sklearn")
    if not np.isfinite(adata.obsp["connectivities"].data).all():
        raise ValueError("Nonfinite neighbor graph")
    # run leiden clustering on the directed graph, using the actual number of neighbors and PCs
    sc.tl.leiden(adata, resolution=resolution, random_state=seed,
                 flavor="leidenalg", directed=True, use_weights=True, n_iterations=-1)
    # Random initialization also works for small/disconnected graphs.
    sc.tl.umap(adata, random_state=seed, init_pos="random")
    labels[indices] = adata.obs["leiden"].astype(str).to_numpy()
    status[indices] = "clustered"
    embedding[indices] = adata.obsm["X_umap"]
    # Carry the actual Scanpy connectivity graph into the R plotting step.
    # Only one copy of each symmetric edge is sent; no graph matrix is saved.
    from scipy.sparse import triu
    graph = triu(adata.obsp["connectivities"], k=1).tocoo()
    graph.eliminate_zeros()
    edges = np.column_stack([indices[graph.row], indices[graph.col], graph.data])
    info.update(actual_pcs=actual_pcs, actual_neighbors=actual_neighbors,
                n_clustered=len(indices), n_clusters=len(set(labels[indices])),
                pca_variance_explained=float(adata.uns["pca"]["variance_ratio"].sum()),
                neighbor_metric="correlation", neighbor_representation="X_pca",
                leiden_flavor="leidenalg", leiden_directed=True)
    return labels, status, embedding, info, edges
