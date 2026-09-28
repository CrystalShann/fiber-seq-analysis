#!/usr/bin/env python3
"""Prepare full-resolution ACF tables, global heatmap ranks and expression comparisons."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import tss_common
from tss_common import (
    assert_alignment, fingerprint, finish_stage, load_parent, log,
    prepare_dirs, read_metadata, stage_valid, write_json,
)


DEFAULT_OUT = Path("/project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss")
PARENT = Path(__file__).resolve().parent.parent / "05_summarize_autocorrelations.py"
BINS = ("Q1_low", "Q2", "Q3", "Q4_high")


def aggregate(profiles, mask, group_column, group):
    """Summaries preserve every retained lag; mask includes only valid ACF rows."""
    subset = profiles[mask]
    width = profiles.shape[1]
    return pd.DataFrame({
        group_column: group, "lag_bp": np.arange(width),
        "mean_acf": subset.mean(axis=0) if len(subset) else np.full(width, np.nan),
        "median_acf": np.median(subset, axis=0) if len(subset) else np.full(width, np.nan),
        "n_reads": len(subset),
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--force", action="store_true", help="Recompute this stage even if its cache is valid.")
    args = parser.parse_args()
    out = args.out_dir.resolve()
    prepare_dirs(out)
    binary_path = out / "intermediate/binary_m6a.npy"
    sampled_path = out / "tables/sampled_molecules.tsv"
    metadata_path = out / "tables/clustered_molecules.tsv"
    rows_path = out / "intermediate/row_ids.npy"
    acf_path = out / "intermediate/acf.npy"
    valid_path = out / "intermediate/acf_valid.npy"
    embedding_path = out / "intermediate/umap.npy"
    signature = fingerprint(
        [binary_path, sampled_path, metadata_path, rows_path, acf_path, valid_path, embedding_path,
         Path(__file__), Path(tss_common.__file__), PARENT],
        {"expression_bins": list(BINS), "heatmap_rank": "numeric cluster then descending m6a_count then row_index",
         "heatmap_matrix_order": "original sampled row order", "numpy_version": np.__version__,
         "pandas_version": pd.__version__},
    )
    if not args.force and stage_valid(out, "04_plot_tables", signature):
        log("04 Tables: validated outputs already exist; skipping")
        return

    binary = np.load(binary_path, mmap_mode="r", allow_pickle=False)
    row_ids = np.load(rows_path, allow_pickle=False)
    profiles = np.load(acf_path, mmap_mode="r", allow_pickle=False)
    valid = np.load(valid_path, allow_pickle=False)
    embedding = np.load(embedding_path, allow_pickle=False)
    sampled = read_metadata(sampled_path)
    records = read_metadata(metadata_path)
    # Metadata loading preserves literal gene/sample names; decode only numeric missing values.
    for column in ("umap1", "umap2"):
        records[column] = pd.to_numeric(records[column].replace("NA", np.nan), errors="raise")
    assert_alignment(binary, sampled, row_ids)
    assert_alignment(binary, records, row_ids)
    pd.testing.assert_frame_equal(records[sampled.columns], sampled,
                                  check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12)
    if profiles.ndim != 2 or len(profiles) != len(records) or valid.shape != (len(records),):
        raise ValueError("ACF rows are not aligned with the sampled metadata")
    expected_valid = (binary.sum(axis=1) > 0) & (binary.sum(axis=1) < binary.shape[1])
    if (valid.dtype != np.bool_ or not np.array_equal(valid, records.acf_valid.to_numpy())
            or not np.array_equal(valid, expected_valid)):
        raise ValueError("Metadata ACF validity flags differ from the saved validity mask")
    if not np.isfinite(profiles[valid]).all() or not np.isnan(profiles[~valid]).all():
        raise ValueError("Valid ACFs must be finite and zero-variance ACFs must be NaN")
    if not np.allclose(profiles[valid, 0], 1, rtol=0, atol=1e-10):
        raise ValueError("Lag 0 must equal one in valid ACF rows")
    if embedding.shape != (len(records), 2) or not np.allclose(
            records[["umap1", "umap2"]].to_numpy(), embedding, equal_nan=True, rtol=1e-12, atol=1e-12):
        raise ValueError("UMAP coordinates are not aligned to the sampled metadata")
    if set(records.expr_bin) != set(BINS):
        raise ValueError("Plotting requires exactly Q1_low, Q2, Q3, Q4_high expression bins")
    clustered = records.status.eq("clustered").to_numpy()
    if not clustered.any() or (clustered & ~valid).any():
        raise ValueError("Plotting requires at least one cluster and no clustered zero-variance reads")
    if not records.loc[~clustered, "cluster"].eq("Unclustered").all():
        raise ValueError("All excluded molecules must retain the explicit Unclustered label")
    cluster_ids = sorted(records.loc[clustered, "cluster"].unique().tolist(), key=int)
    # Always retain an Unclustered composition category so its zero count is explicit.
    composition_groups = cluster_ids + ["Unclustered"]
    cluster_order = {cluster: i for i, cluster in enumerate(composition_groups)}
    rank_order = records.assign(_cluster_order=records.cluster.map(cluster_order)).sort_values(
        ["_cluster_order", "m6a_count", "row_index"], ascending=[True, False, True], kind="stable",
    ).index.to_numpy()
    ranks = np.empty(len(records), dtype=np.int64)
    ranks[rank_order] = np.arange(1, len(records) + 1)
    records["heatmap_rank"] = ranks
    if not np.array_equal(np.sort(ranks), np.arange(1, len(records) + 1)):
        raise ValueError("Heatmap ranks must be a unique one-based permutation")
    log(f"04 Tables: preparing {len(records):,} original-order rows and {len(cluster_ids)} global clusters")
    parent = load_parent("05_summarize_autocorrelations.py")
    cluster_curves, expression_curves, summaries, composition = [], [], [], []
    for cluster in composition_groups:
        mask = records.cluster.eq(cluster).to_numpy()
        curve = aggregate(profiles, mask & valid, "cluster", cluster)
        if cluster != "Unclustered":
            cluster_curves.append(curve)
        peak_lag, peak_height = parent.repeat_peak(curve.mean_acf.to_numpy())
        n = int(mask.sum())
        counts = records.loc[mask, "m6a_count"]
        summaries.append({
            "cluster": cluster, "n_reads": n, "n_acf_valid": int((mask & valid).sum()),
            "fraction_of_sample": n / len(records),
            "mean_m6a_count": counts.mean() if n else np.nan,
            "median_m6a_count": counts.median() if n else np.nan,
            "positive_local_peak_140_250_bp": peak_lag, "peak_acf": peak_height,
            "peak_band_complete": profiles.shape[1] >= 252,
        })
        for expr_bin in BINS:
            in_bin = records.expr_bin.eq(expr_bin).to_numpy()
            count = int((mask & in_bin).sum())
            total_bin = int(in_bin.sum())
            composition.append({
                "cluster": cluster, "expr_bin": expr_bin, "n_reads": count,
                "n_bin_total": total_bin, "fraction_within_bin": count / total_bin,
                "fraction_within_cluster": count / n if n else np.nan,
            })
    for expr_bin in BINS:
        expression_curves.append(aggregate(
            profiles, valid & records.expr_bin.eq(expr_bin).to_numpy(), "expr_bin", expr_bin,
        ))
    composition_table = pd.DataFrame(composition)
    fractions = composition_table.groupby("expr_bin").fraction_within_bin.sum()
    if not np.allclose(fractions.to_numpy(), 1, rtol=0, atol=1e-12):
        raise ValueError("Cluster composition must account for all sampled molecules in each bin")
    outputs = [out / "tables/plot_metadata.tsv", out / "tables/acf_heatmap.tsv.gz",
               out / "tables/cluster_acf.tsv", out / "tables/expression_acf.tsv",
               out / "tables/cluster_composition.tsv", out / "tables/cluster_summary.tsv",
               out / "validation/plot_table_validation.json"]
    records.to_csv(outputs[0], sep="\t", index=False, na_rep="NA")
    heatmap = pd.DataFrame(profiles, columns=[f"lag_{i}" for i in range(profiles.shape[1])])
    heatmap.insert(0, "row_index", records.row_index.to_numpy())
    # Full precision, all lags and every original row; display sorting happens only in R.
    heatmap.to_csv(outputs[1], sep="\t", index=False, na_rep="NA",
                   compression={"method": "gzip", "compresslevel": 6, "mtime": 0}, chunksize=128)
    pd.concat(cluster_curves, ignore_index=True).to_csv(outputs[2], sep="\t", index=False, na_rep="NA")
    pd.concat(expression_curves, ignore_index=True).to_csv(outputs[3], sep="\t", index=False, na_rep="NA")
    composition_table.to_csv(outputs[4], sep="\t", index=False, na_rep="NA")
    pd.DataFrame(summaries).to_csv(outputs[5], sep="\t", index=False, na_rep="NA")
    report = {
        "n_reads": len(records), "n_acf_valid": int(valid.sum()),
        "n_clustered": int(clustered.sum()), "n_unclustered": int((~clustered).sum()),
        "n_clusters": len(cluster_ids), "n_lags": profiles.shape[1],
        "metadata_and_matrix_row_order_preserved": True, "lag0_equals_one": True,
        "heatmap_rank_unique_one_based": True,
        "heatmap_sort": "numeric global Leiden cluster; descending m6a count; original row_index breaks ties",
        "heatmap_matrix_original_row_order": True, "cluster_labels_shared_across_bins": True,
        "composition_denominator": "all sampled molecules in each bin, including Unclustered",
        "composition_fractions_sum_to_one": True,
        "expression_acf_population": "all nonzero-variance sampled reads, including any valid unclustered reads",
        "repeat_peak_parent_function": str(PARENT) + ":repeat_peak",
    }
    write_json(outputs[6], report)
    finish_stage(out, "04_plot_tables", signature, outputs, details=report)
    log("04 Tables: finished; full-resolution matrix, curves, compositions and metadata are ready for R")


if __name__ == "__main__":
    main()
