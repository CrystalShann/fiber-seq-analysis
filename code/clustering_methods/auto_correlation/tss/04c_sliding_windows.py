#!/usr/bin/env python3
"""Sliding strand-oriented windows over this task's own sampled molecules.

Never resamples: every molecule of the region contributes to every window of
the region's span (at most 2,000 bp). For each window the oriented binary row
is sliced, the parent ``autocorrelations()`` gives lags 0..win_width-1 at 1 bp,
and ``tss_nrl`` metrics are applied (decay metrics on lags <= floor(N/2)).
Windows that overlap the NDR core [-100, +100) are flagged
``overlaps_ndr_core``. A sequence control repeats the analysis on the
strand-oriented reference A/T indicator of the same windows.
"""
import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

import tss_common
import tss_nrl
from tss_common import (BINS, DEFAULT_OUT, assert_alignment, fingerprint, finish_stage, load_parent,
                        log, prepare_dirs, read_metadata, stage_valid, window_width, write_json)

REF_FA = "/project/spott/reference/human/GRCh38/hg38.fa"  # as in code/accessibility/expr_access/02_tss_m6a_profiles.py
PARENT = Path(__file__).resolve().parent.parent / "03_compute_autocorrelations.py"
META = ["read_id", "sample", "timepoint", "lps_minutes", "gene_id", "gene_name", "chrom", "tss", "strand",
        "expr_bin", "mean_tpm", "cluster"]
LONG_METRICS = ["nrl_bp", "nrl_status", "nrl_peak_height", "peak2_height", "peak_ratio",
                "decay_length_bp", "damping_lag_bp", "fit_r2"]
NDR_CORE = (-100, 100)


def median_or_nan(series):
    values = series.dropna()
    return float(values.median()) if len(values) else np.nan


def mean_or_nan(series):
    return float(series.mean()) if len(series) else np.nan


def region_label(offset_start, offset_end):
    if offset_start < 0 < offset_end:
        return "span_2kb" if offset_end - offset_start == 2000 else "span"
    return "upstream" if offset_end <= 0 else "downstream"


def window_starts(offset_start, offset_end, win_width, win_step):
    span = offset_end - offset_start
    if win_width > span:
        raise ValueError(f"win_width {win_width} exceeds the {span}-bp span")
    if span > 2000:
        raise ValueError("The sliding span is the task's own window and never exceeds 2,000 bp")
    starts = list(range(offset_start, offset_end - win_width + 1, win_step))
    if not starts or any(s < offset_start or s + win_width > offset_end for s in starts):
        raise ValueError("Every sliding window must lie entirely inside the span")
    return starts


def at_indicator(fasta, metadata, width):
    """Strand-oriented A/T indicator (1 = reference base A or T) for every molecule's window."""
    indicator = np.zeros((len(metadata), width), dtype=np.uint8)
    for row in metadata.itertuples(index=False):
        sequence = fasta.fetch(row.chrom, int(row.window_start), int(row.window_end)).upper()
        if len(sequence) != width:
            raise ValueError(f"Reference window length mismatch for {row.read_id}")
        values = np.frombuffer(sequence.encode(), dtype=np.uint8)
        at = ((values == ord("A")) | (values == ord("T"))).astype(np.uint8)
        indicator[row.row_index] = at[::-1] if row.strand == "-" else at
    return indicator


def calibration_rows(binary, starts, offset_start, win_width, n_null, seed):
    """Up to n_null (read, window) slices drawn across all windows for the prominence calibration."""
    rng = np.random.default_rng(seed)
    total = len(binary) * len(starts)
    picks = np.sort(rng.choice(total, min(n_null, total), replace=False))
    rows = np.empty((len(picks), win_width), dtype=np.uint8)
    for k, pick in enumerate(picks):
        read, window = divmod(int(pick), len(starts))
        j0 = starts[window] - offset_start
        rows[k] = binary[read, j0: j0 + win_width]
    return rows


def window_metrics(matrix, parent, win_width, min_prominence, params, chunk=1000):
    """Parent ACF + tss_nrl metrics for one window slice, in read chunks."""
    records = []
    for start in range(0, len(matrix), chunk):
        block = np.ascontiguousarray(matrix[start: start + chunk])
        acf, valid = parent.autocorrelations(block, n_features=win_width)
        if acf.shape != (len(block), win_width):
            raise ValueError("Parent ACF did not return exactly win_width lags")
        if valid.any() and not np.allclose(acf[valid, 0], 1, rtol=0, atol=1e-10):
            raise ValueError("Lag 0 must equal one for valid windows")
        for i in range(len(block)):
            metrics = tss_nrl.molecule_metrics(acf[i], win_width, min_prominence, **params)
            metrics["m6a_count_win"] = int(block[i].sum())
            metrics["acf_valid"] = bool(valid[i])
            records.append(metrics)
        yield start, acf, valid, records
        records = []


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--win-width", type=int, default=500)
    parser.add_argument("--win-step", type=int, default=100, help="Window step (100 or 200 are typical).")
    parser.add_argument("--ref", type=Path, default=Path(REF_FA), help="Reference FASTA for the A/T control.")
    tss_nrl.add_arguments(parser)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    tss_nrl.validate_arguments(parser, args)
    if args.win_width < 10 or args.win_step < 1:
        parser.error("Require win_width >= 10 and win_step >= 1")
    out = args.out_dir.resolve()
    prepare_dirs(out)
    binary_path = out / "intermediate/binary_m6a.npy"
    rows_path = out / "intermediate/row_ids.npy"
    metadata_path = out / "tables/clustered_molecules.tsv"
    if not args.ref.is_file() or not Path(str(args.ref) + ".fai").is_file():
        raise SystemExit(f"Reference FASTA and .fai are required for the A/T control: {args.ref}")
    params = tss_nrl.metric_parameters(args)
    decay_max_lag = args.decay_max_lag if args.decay_max_lag is not None else tss_nrl.default_decay_max_lag(args.win_width)
    decay_max_lag = min(decay_max_lag, args.win_width - 1)
    params["decay_max_lag"] = decay_max_lag
    signature = fingerprint([binary_path, rows_path, metadata_path, Path(__file__), Path(tss_common.__file__),
                             Path(tss_nrl.__file__), PARENT, str(args.ref) + ".fai"],
                            {**params, "win_width": args.win_width, "win_step": args.win_step,
                             "ref": str(args.ref.resolve()), "min_prominence": args.min_prominence,
                             "prominence_quantile": args.prominence_quantile, "n_null": args.n_null})
    if not args.force and stage_valid(out, "04c_sliding", signature):
        log("04c sliding: validated outputs already exist; skipping")
        return
    binary = np.load(binary_path, mmap_mode="r", allow_pickle=False)
    row_ids = np.load(rows_path, allow_pickle=False)
    metadata = read_metadata(metadata_path)
    assert_alignment(binary, metadata, row_ids)
    width = window_width(metadata)
    if metadata.orientation.ne("transcriptional").any():
        raise ValueError("Sliding windows require the strand-oriented (transcriptional) binary matrix")
    offset_start = int(metadata.window_offset_start.unique()[0])
    offset_end = int(metadata.window_offset_end.unique()[0])
    if metadata.window_offset_start.nunique() != 1 or metadata.window_offset_end.nunique() != 1 \
            or offset_end - offset_start != width:
        raise ValueError("All molecules must share one TSS-relative window")
    region = region_label(offset_start, offset_end)
    starts = window_starts(offset_start, offset_end, args.win_width, args.win_step)
    log(f"04c sliding: region {region} span [{offset_start}, {offset_end}), {len(starts)} windows of "
        f"{args.win_width} bp (step {args.win_step}); lags 0-{args.win_width - 1} at 1 bp; "
        f"decay lags <= {decay_max_lag}; {len(metadata):,} molecules")
    parent = load_parent("03_compute_autocorrelations.py")
    # min_prominence for the win_width-bp ACFs: calibrated on position-shuffled
    # (read, window) slices of this region's own molecules (tss_nrl.calibrate_prominence),
    # and shared by the m6A windows and the A/T control; --min-prominence overrides.
    null_rows = calibration_rows(binary, starts, offset_start, args.win_width, args.n_null, seed=0)
    min_prominence, calibration = tss_nrl.resolve_prominence(
        args, null_rows, lambda rows: parent.autocorrelations(rows, n_features=args.win_width))
    log(f"04c sliding: min_prominence = {min_prominence:.4f} ({calibration['method']})")
    with pysam.FastaFile(str(args.ref)) as fasta:
        at_matrix = at_indicator(fasta, metadata, width)
    meta = metadata[META].reset_index(drop=True)
    long_path = out / "tables/sliding_nrl_long.tsv.gz"
    mean_rows, control_rows, long_counts = [], [], 0
    read_ids_seen = set()
    with gzip.open(long_path, "wt", compresslevel=6) as handle:
        first = True
        for win_start in starts:
            win_end = win_start + args.win_width
            j0 = win_start - offset_start
            overlaps = win_start < NDR_CORE[1] and win_end > NDR_CORE[0]
            label = dict(region=region, win_start=win_start, win_end=win_end, win_center=(win_start + win_end) / 2,
                         overlaps_ndr_core=overlaps)
            acf_sum = {b: np.zeros(args.win_width) for b in BINS}
            acf_n = {b: 0 for b in BINS}
            at_sum = {b: np.zeros(args.win_width) for b in BINS}
            at_n = {b: 0 for b in BINS}
            control = []
            for start, acf, valid, records in window_metrics(binary[:, j0:j0 + args.win_width], parent,
                                                              args.win_width, min_prominence, params):
                table = pd.DataFrame(records)
                block_meta = meta.iloc[start: start + len(table)].reset_index(drop=True)
                for b in BINS:
                    mask = block_meta.expr_bin.eq(b).to_numpy() & valid
                    acf_sum[b] += acf[mask].sum(axis=0)
                    acf_n[b] += int(mask.sum())
                table = pd.concat([block_meta, pd.DataFrame([label] * len(table)),
                                   table[["m6a_count_win"] + LONG_METRICS]], axis=1)
                table.to_csv(handle, sep="\t", index=False, header=first, na_rep="NA", float_format="%.6g")
                first = False
                long_counts += len(table)
                read_ids_seen.update(block_meta.read_id.tolist())
            for start, acf, valid, records in window_metrics(at_matrix[:, j0:j0 + args.win_width], parent,
                                                              args.win_width, min_prominence, params):
                block_meta = meta.iloc[start: start + len(records)]
                for b in BINS:
                    mask = block_meta.expr_bin.eq(b).to_numpy() & valid
                    at_sum[b] += acf[mask].sum(axis=0)
                    at_n[b] += int(mask.sum())
                control.append(pd.concat([block_meta[["expr_bin"]].reset_index(drop=True),
                                          pd.DataFrame(records)[["nrl_bp", "nrl_status", "decay_length_bp",
                                                                 "damping_lag_bp"]]], axis=1))
            control = pd.concat(control, ignore_index=True)
            for b in list(BINS) + ["all"]:
                sub = control if b == "all" else control[control.expr_bin.eq(b)]
                control_rows.append({"expr_bin": b, **label, "n": len(sub),
                                     "fraction_with_nrl_peak": mean_or_nan(sub.nrl_bp.notna()),
                                     "median_nrl_bp": median_or_nan(sub.nrl_bp),
                                     "median_decay_length_bp": median_or_nan(sub.decay_length_bp),
                                     "median_damping_lag_bp": median_or_nan(sub.damping_lag_bp)})
            for b in BINS:
                for signal, sums, counts in (("m6a", acf_sum, acf_n), ("at_control", at_sum, at_n)):
                    mean = sums[b] / counts[b] if counts[b] else np.full(args.win_width, np.nan)
                    mean_rows.append(pd.DataFrame({"signal": signal, "expr_bin": b, **label,
                                                   "lag_bp": np.arange(args.win_width), "mean_acf": mean,
                                                   "n_reads": counts[b]}))
            log(f"04c sliding: window [{win_start}, {win_end}) done" + (" (overlaps NDR core)" if overlaps else ""))
    if read_ids_seen != set(metadata.read_id):
        raise ValueError("Sliding long table does not contain exactly the sampled molecules")
    if long_counts != len(metadata) * len(starts):
        raise ValueError("Expected one long row per (read, window)")
    long = pd.read_csv(long_path, sep="\t", keep_default_na=False, na_values=["NA"],
                       dtype={"cluster": str, "read_id": str, "gene_id": str})
    long["has_peak"] = long.nrl_bp.notna()
    long["m6a_density"] = long.m6a_count_win / args.win_width
    summary = []
    keys = ["expr_bin", "timepoint", "region", "win_start", "win_end", "win_center", "overlaps_ndr_core"]
    pooled = long.assign(timepoint="all")
    for _, sub in pd.concat([long, pooled], ignore_index=True).groupby(keys, sort=False):
        nrl = sub.nrl_bp.dropna()
        summary.append({**dict(zip(keys, _)), "n": len(sub), "n_with_nrl_peak": int(sub.has_peak.sum()),
                        "fraction_with_nrl_peak": mean_or_nan(sub.has_peak), "median_nrl_bp": median_or_nan(nrl),
                        "iqr_low_nrl_bp": nrl.quantile(0.25) if len(nrl) else np.nan,
                        "iqr_high_nrl_bp": nrl.quantile(0.75) if len(nrl) else np.nan,
                        "median_decay_length_bp": median_or_nan(sub.decay_length_bp),
                        "median_damping_lag_bp": median_or_nan(sub.damping_lag_bp),
                        "mean_m6a_density": mean_or_nan(sub.m6a_density), "n_genes": sub.gene_id.nunique()})
    summary = pd.DataFrame(summary)
    summary["min_prominence"] = min_prominence
    control_table = pd.DataFrame(control_rows)
    control_table["min_prominence"] = min_prominence
    outputs = [long_path, out / "tables/sliding_summary.tsv", out / "tables/sliding_mean_acf.tsv.gz",
               out / "tables/sliding_at_control.tsv", out / "tables/sliding_prominence_calibration.tsv",
               out / "validation/sliding_validation.json"]
    summary.to_csv(outputs[1], sep="\t", index=False, na_rep="NA", float_format="%.6g")
    pd.concat(mean_rows, ignore_index=True).to_csv(outputs[2], sep="\t", index=False, na_rep="NA",
                                                   float_format="%.8g", compression={"method": "gzip", "mtime": 0})
    control_table.to_csv(outputs[3], sep="\t", index=False, na_rep="NA", float_format="%.6g")
    pd.DataFrame([{"region": region, "acf_length_bp": args.win_width, **calibration}]).to_csv(
        outputs[4], sep="\t", index=False, na_rep="NA", float_format="%.6g")
    report = {"region": region, "span": [offset_start, offset_end], "n_molecules": len(metadata),
              "n_windows": len(starts), "win_width": args.win_width, "win_step": args.win_step,
              "window_starts": starts, "all_windows_inside_span": True,
              "n_windows_overlapping_ndr_core": int(sum(s < NDR_CORE[1] and s + args.win_width > NDR_CORE[0] for s in starts)),
              "lags": [0, args.win_width - 1], "lag_resolution_bp": 1,
              "long_rows": long_counts, "read_ids_match_sampled_molecules": True,
              "status_counts": long.nrl_status.value_counts().to_dict(),
              "min_prominence": min_prominence, "prominence_calibration": calibration,
              "parameters": {**params, "min_prominence": min_prominence}, "reference_fasta": str(args.ref.resolve()),
              "orientation": "transcriptional; A/T control uses the same reversal for - strand genes",
              "limitations": ["500 bp holds only about 2-3 nucleosome periods: single-molecule NRL and decay "
                              "estimates are noisy; compare medians between groups",
                              "at lag k only N - k base pairs contribute to the ACF",
                              "windows overlapping [-100, +100) partly reflect the NDR edge rather than spacing "
                              "(Clarkson et al. 2019 NAR)"]}
    write_json(outputs[5], report)
    finish_stage(out, "04c_sliding", signature, outputs, details=report)
    log(f"04c sliding: finished; {len(starts)} windows, {long_counts:,} (read, window) rows")


if __name__ == "__main__":
    main()
