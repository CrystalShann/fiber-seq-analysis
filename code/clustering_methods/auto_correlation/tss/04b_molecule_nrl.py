#!/usr/bin/env python3
"""Per-molecule NRL and regularity metrics on the fixed-window ACFs (tss_nrl.py).

Uses the saved parent ACF rows (every 1 bp lag of the window) and the joint
Leiden labels. The raw ACF is unchanged; only the decay metrics use the
unbiased N/(N-k) rescaling restricted to lags <= decay_max_lag.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import tss_common
import tss_nrl
from tss_common import (BINS, DEFAULT_OUT, assert_alignment, fingerprint, finish_stage, load_parent,
                        log, prepare_dirs, read_metadata, stage_valid, window_width, write_json)

PARENT = Path(__file__).resolve().parent.parent / "05_summarize_autocorrelations.py"
META = ["read_id", "sample", "timepoint", "lps_minutes", "gene_id", "gene_name", "chrom", "tss", "strand",
        "expr_bin", "mean_tpm", "m6a_count", "window", "window_offset_start", "window_offset_end", "cluster"]
METRICS = ["nrl_bp", "nrl_status", "nrl_peak_height", "peak2_lag_bp", "peak2_height", "peak_ratio",
           "repeat_peak_lag", "repeat_peak_value", "decay_length_bp", "fit_period_bp", "fit_r2",
           "fit_converged", "damping_lag_bp"]


def summarize(table, group_column, groups):
    rows = []
    for group in groups:
        sub = table[table[group_column].eq(group)]
        peak = sub.nrl_bp.notna()
        row = {group_column: group, "n": len(sub), "n_acf_valid": int(sub.nrl_status.ne("zero_variance").sum()),
               "n_with_nrl_peak": int(peak.sum()), "fraction_with_nrl_peak": peak.mean() if len(sub) else np.nan}
        for column in ("nrl_bp", "decay_length_bp", "damping_lag_bp", "peak_ratio", "fit_r2"):
            values = sub[column].dropna()
            row[f"median_{column}"] = values.median() if len(values) else np.nan
            row[f"iqr_low_{column}"] = values.quantile(0.25) if len(values) else np.nan
            row[f"iqr_high_{column}"] = values.quantile(0.75) if len(values) else np.nan
        row["fraction_fit_converged"] = sub.fit_converged.mean() if len(sub) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--window-name", default=None, help="Window label for the table; default = folder name.")
    tss_nrl.add_arguments(parser)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    tss_nrl.validate_arguments(parser, args)
    out = args.out_dir.resolve()
    prepare_dirs(out)
    binary_path = out / "intermediate/binary_m6a.npy"
    rows_path = out / "intermediate/row_ids.npy"
    acf_path = out / "intermediate/acf.npy"
    valid_path = out / "intermediate/acf_valid.npy"
    metadata_path = out / "tables/clustered_molecules.tsv"
    params = tss_nrl.metric_parameters(args)
    signature = fingerprint([binary_path, rows_path, acf_path, valid_path, metadata_path, Path(__file__),
                             Path(tss_common.__file__), Path(tss_nrl.__file__), PARENT],
                            {**params, "window_name": args.window_name, "min_prominence": args.min_prominence,
                             "prominence_quantile": args.prominence_quantile, "n_null": args.n_null})
    if not args.force and stage_valid(out, "04b_nrl", signature):
        log("04b NRL: validated outputs already exist; skipping")
        return
    binary = np.load(binary_path, mmap_mode="r", allow_pickle=False)
    row_ids = np.load(rows_path, allow_pickle=False)
    metadata = read_metadata(metadata_path)
    assert_alignment(binary, metadata, row_ids)
    width = window_width(metadata)
    profiles = np.load(acf_path, mmap_mode="r", allow_pickle=False)
    valid = np.load(valid_path, allow_pickle=False)
    if profiles.shape[0] != len(metadata) or valid.shape != (len(metadata),):
        raise ValueError("ACF rows are not aligned with the sampled metadata")
    if profiles.shape[1] < width:
        log(f"04b NRL: ACF has {profiles.shape[1]} lags of a {width}-bp window (--max-lag was set)")
    decay_max_lag = args.decay_max_lag if args.decay_max_lag is not None else tss_nrl.default_decay_max_lag(width)
    if decay_max_lag >= profiles.shape[1]:
        raise ValueError("decay_max_lag exceeds the saved lags")
    window_name = args.window_name or out.name
    # min_prominence: calibrated on position-shuffled copies of this window's own
    # molecules at this window length (see tss_nrl.calibrate_prominence), unless overridden.
    acf_parent = load_parent("03_compute_autocorrelations.py")
    min_prominence, calibration = tss_nrl.resolve_prominence(
        args, binary, lambda rows: acf_parent.autocorrelations(rows, n_features=profiles.shape[1]))
    log(f"04b NRL: min_prominence = {min_prominence:.4f} ({calibration['method']})")
    log(f"04b NRL: {len(metadata):,} molecules, window {width} bp, lags 0-{profiles.shape[1] - 1} at 1 bp; "
        f"NRL band {args.nrl_min}-{args.nrl_max}, min_lag {args.min_lag}, prominence {min_prominence:.4f}, "
        f"decay lags <= {decay_max_lag}, flat threshold {args.flat_threshold}")
    parent = load_parent("05_summarize_autocorrelations.py")
    records = []
    for i in range(len(metadata)):
        metrics = tss_nrl.molecule_metrics(
            np.asarray(profiles[i]), width, min_prominence, repeat_peak=parent.repeat_peak,
            **{**params, "decay_max_lag": decay_max_lag})
        records.append({column: metrics[column] for column in METRICS})
        if (i + 1) % 2000 == 0:
            log(f"04b NRL: {i + 1:,} molecules done")
    metrics_table = pd.DataFrame(records)
    table = metadata[[c for c in META if c != "window"]].copy()
    table.insert(META.index("window"), "window", window_name)
    table = pd.concat([table.reset_index(drop=True), metrics_table], axis=1)[META + METRICS]
    if not table.nrl_status.isin(tss_nrl.STATUSES).all():
        raise ValueError("Unknown NRL status")
    if not np.array_equal(table.nrl_status.eq("zero_variance").to_numpy(), ~valid):
        raise ValueError("zero_variance status must coincide with the invalid ACF rows")
    if not (table.nrl_bp.notna() == table.nrl_status.eq("ok")).all():
        raise ValueError("nrl_bp must be present exactly for status ok")
    if not table.nrl_bp.dropna().between(args.nrl_min, args.nrl_max).all():
        raise ValueError("NRL outside the configured band")
    clusters = sorted(table.cluster.unique(), key=lambda c: (c == "Unclustered", int(c) if c != "Unclustered" else 0))
    outputs = [out / "tables/nrl_per_molecule.tsv", out / "tables/nrl_summary_by_bin.tsv",
               out / "tables/nrl_summary_by_cluster.tsv", out / "tables/nrl_prominence_calibration.tsv",
               out / "validation/nrl_validation.json"]
    table.to_csv(outputs[0], sep="\t", index=False, na_rep="NA", float_format="%.6g")
    for path, column, groups in ((outputs[1], "expr_bin", list(BINS)), (outputs[2], "cluster", clusters)):
        summary = summarize(table, column, groups)
        summary["min_prominence"] = min_prominence
        summary.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.6g")
    pd.DataFrame([{"window": window_name, "acf_length_bp": int(profiles.shape[1]), **calibration}]).to_csv(
        outputs[3], sep="\t", index=False, na_rep="NA", float_format="%.6g")
    report = {"n_molecules": len(table), "window": window_name, "window_width_bp": width,
              "lags": [0, int(profiles.shape[1] - 1)], "lag_resolution_bp": 1,
              "status_counts": table.nrl_status.value_counts().to_dict(),
              "fraction_with_nrl_peak": float(table.nrl_bp.notna().mean()),
              "median_nrl_bp": float(table.nrl_bp.median()) if table.nrl_bp.notna().any() else None,
              "min_prominence": min_prominence, "prominence_calibration": calibration,
              "parameters": {**params, "decay_max_lag": decay_max_lag, "min_prominence": min_prominence},
              "decay_rescaling": "acf[k] * N / (N - k) for k <= decay_max_lag, decay metrics only; raw ACF unchanged",
              "nrl_reference": "Abdulhay et al. 2020 eLife (secondary-peak scan)",
              "decay_reference": "Baldi et al. 2018 Mol Cell (ACF damping and array regularity)"}
    write_json(outputs[4], report)
    finish_stage(out, "04b_nrl", signature, outputs, details=report)
    log(f"04b NRL: finished; {report['status_counts']}")


if __name__ == "__main__":
    main()
