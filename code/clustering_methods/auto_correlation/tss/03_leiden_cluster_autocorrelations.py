#!/usr/bin/env python3
"""Jointly cluster all sampled expression bins with the unchanged parent workflow."""

import argparse
from importlib.metadata import version
from pathlib import Path

import numpy as np

import tss_common
from tss_common import (
    assert_alignment, fingerprint, finish_stage, load_parent, log,
    prepare_dirs, read_metadata, stage_valid, write_json,
)


DEFAULT_OUT = Path("/project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss")
PARENT = Path(__file__).resolve().parent.parent / "04_cluster_autocorrelations.py"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-neighbors", type=int, default=10)
    parser.add_argument("--resolution", type=float, default=0.4)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--force", action="store_true", help="Recompute this stage even if its cache is valid.")
    args = parser.parse_args()
    if (args.n_pcs < 2 or args.n_neighbors < 2 or not np.isfinite(args.resolution)
            or args.resolution <= 0 or args.seed < 0):
        parser.error("Require >=2 PCs, >=2 neighbors, positive resolution and nonnegative seed")
    out = args.out_dir.resolve()
    prepare_dirs(out)
    binary_path = out / "intermediate/binary_m6a.npy"
    metadata_path = out / "tables/sampled_molecules.tsv"
    rows_path = out / "intermediate/row_ids.npy"
    acf_path = out / "intermediate/acf.npy"
    valid_path = out / "intermediate/acf_valid.npy"
    parameters = {"n_pcs": args.n_pcs, "n_neighbors": args.n_neighbors,
                  "resolution": args.resolution, "seed": args.seed, "metric": "correlation",
                  "joint_expression_bins": True,
                  "package_versions": {p: version(p) for p in
                      ("numpy", "scipy", "scanpy", "anndata", "scikit-learn", "umap-learn", "leidenalg", "igraph")}}
    signature = fingerprint(
        [binary_path, metadata_path, rows_path, acf_path, valid_path,
         Path(__file__), Path(tss_common.__file__), PARENT], parameters,
    )
    if not args.force and stage_valid(out, "03_clustering", signature):
        log("03 Leiden: validated outputs already exist; skipping")
        return

    binary = np.load(binary_path, mmap_mode="r", allow_pickle=False)
    metadata = read_metadata(metadata_path)
    row_ids = np.load(rows_path, allow_pickle=False)
    assert_alignment(binary, metadata, row_ids)
    profiles = np.load(acf_path, mmap_mode="r", allow_pickle=False)
    valid = np.load(valid_path, allow_pickle=False)
    expected_valid = (binary.sum(axis=1) > 0) & (binary.sum(axis=1) < binary.shape[1])
    if profiles.ndim != 2 or len(profiles) != len(metadata) or profiles.shape[1] < 3:
        raise ValueError("Aligned ACF profiles with at least three lag columns are required for Leiden")
    if valid.dtype != np.bool_ or valid.shape != (len(metadata),) or not np.array_equal(valid, expected_valid):
        raise ValueError("ACF validity mask is not aligned to the sampled nonzero-variance molecules")
    if not np.isfinite(profiles[valid]).all() or not np.isnan(profiles[~valid]).all():
        raise ValueError("Unexpected nonfinite valid ACFs or finite zero-variance ACFs")
    if not np.allclose(profiles[valid, 0], 1, rtol=0, atol=1e-10):
        raise ValueError("Valid ACF lag 0 must equal one")
    log(f"03 Leiden: jointly clustering {int(valid.sum()):,} valid profiles from {len(metadata):,} "
        f"sampled reads; {args.n_pcs} PCs, {args.n_neighbors} neighbors, resolution {args.resolution}")
    parent = load_parent("04_cluster_autocorrelations.py")
    labels, status, embedding, info, edges = parent.cluster_profiles(
        profiles, valid, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors,
        resolution=args.resolution, seed=args.seed,
    )
    n = len(metadata)
    if labels.shape != (n,) or status.shape != (n,) or embedding.shape != (n, 2):
        raise ValueError("Parent clustering changed the sampled row count")
    clustered = status == "clustered"
    if not clustered.any() or info.get("n_clusters", 0) < 1:
        raise RuntimeError("Leiden produced no clusters: " + str(dict(zip(*np.unique(status, return_counts=True)))))
    if (clustered & ~valid).any() or np.any(status[~valid] != "zero_variance"):
        raise ValueError("Zero-variance molecules were not excluded from Leiden clustering")
    if not np.isfinite(embedding[clustered]).all() or np.any(labels[clustered] == ""):
        raise ValueError("Clustered reads must have labels and finite UMAP coordinates")
    if not np.isnan(embedding[~clustered]).all() or np.any(labels[~clustered] != ""):
        raise ValueError("Unclustered reads must retain missing UMAP coordinates and empty parent labels")
    if edges.ndim != 2 or edges.shape[1] != 3 or not np.isfinite(edges).all():
        raise ValueError("Invalid parent connectivity graph")
    if len(edges):
        indices = edges[:, :2]
        if (np.any(indices != np.floor(indices)) or np.any(indices < 0) or np.any(indices >= n)
                or np.any(edges[:, 2] <= 0) or not clustered[indices.astype(int)].all()):
            raise ValueError("Connectivity edges do not index clustered rows in the original sample order")
    records = metadata.copy()
    records["cluster"] = np.where(clustered, labels, "Unclustered")
    records["status"] = status
    records["acf_valid"] = valid
    records["umap1"] = embedding[:, 0]
    records["umap2"] = embedding[:, 1]
    assert_alignment(binary, records, row_ids)
    outputs = [out / "tables/clustered_molecules.tsv", out / "intermediate/umap.npy",
               out / "intermediate/graph_edges.npy", out / "validation/clustering_info.json"]
    records.to_csv(outputs[0], sep="\t", index=False, na_rep="NA")
    np.save(outputs[1], embedding, allow_pickle=False)
    np.save(outputs[2], edges, allow_pickle=False)
    info.update(
        n_input=n, n_acf_valid=int(valid.sum()), n_zero_variance=int((~valid).sum()),
        n_unclustered=int((~clustered).sum()), n_graph_edges=len(edges),
        joint_expression_bins=True, expression_bins=sorted(records.expr_bin.unique().tolist()),
        status_counts={str(k): int(v) for k, v in zip(*np.unique(status, return_counts=True))},
        zero_variance_excluded=True, metadata_and_matrix_row_order_preserved=True,
        parent_function=str(PARENT) + ":cluster_profiles", package_versions=parameters["package_versions"],
    )
    write_json(outputs[3], info)
    finish_stage(out, "03_clustering", signature, outputs, details=info)
    log(f"03 Leiden: finished; {info['n_clusters']} joint clusters, {int(clustered.sum()):,} assigned reads")


if __name__ == "__main__":
    main()
