#!/usr/bin/env python3
"""Compute the unchanged parent ACF for the balanced, previously sampled reads."""

import argparse
from pathlib import Path

import numpy as np

import tss_common
from tss_common import (
    assert_alignment, fingerprint, finish_stage, load_parent, log,
    prepare_dirs, read_metadata, stage_valid, write_json,
)


DEFAULT_OUT = Path("/project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss")
PARENT = Path(__file__).resolve().parent.parent / "03_compute_autocorrelations.py"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--max-lag", type=int, default=1999,
                        help="Inclusive maximum lag; default returns all 2,000 nonnegative lags.")
    parser.add_argument("--force", action="store_true", help="Recompute this stage even if its cache is valid.")
    args = parser.parse_args()
    if not 0 <= args.max_lag <= 1999:
        parser.error("--max-lag must be in [0, 1999] for the 2-kb input window")
    out = args.out_dir.resolve()
    prepare_dirs(out)
    binary_path = out / "intermediate/binary_m6a.npy"
    metadata_path = out / "tables/sampled_molecules.tsv"
    rows_path = out / "intermediate/row_ids.npy"
    signature = fingerprint(
        [binary_path, metadata_path, rows_path, Path(__file__), Path(tss_common.__file__), PARENT],
        {"max_lag": args.max_lag, "numpy_version": np.__version__},
    )
    if not args.force and stage_valid(out, "02_acf", signature):
        log("02 ACF: validated outputs already exist; skipping")
        return

    binary = np.load(binary_path, mmap_mode="r", allow_pickle=False)
    metadata = read_metadata(metadata_path)
    row_ids = np.load(rows_path, allow_pickle=False)
    assert_alignment(binary, metadata, row_ids)
    if binary.shape[1] != 2000 or not len(binary):
        raise ValueError("ACF input must contain sampled molecules in exactly 2,000 genomic bases")
    log(f"02 ACF: computing parent autocorrelations() for {len(binary):,} sampled reads, "
        f"lags 0–{args.max_lag}; no signal transformations")
    parent = load_parent("03_compute_autocorrelations.py")
    profiles, valid = parent.autocorrelations(binary, n_features=args.max_lag + 1)

    expected_valid = (binary.sum(axis=1) > 0) & (binary.sum(axis=1) < binary.shape[1])
    if profiles.shape != (len(metadata), args.max_lag + 1) or valid.shape != (len(metadata),):
        raise ValueError("Parent ACF changed the sampled row count or lag count")
    if valid.dtype != np.bool_ or not np.array_equal(valid, expected_valid):
        raise ValueError("ACF validity flags do not match the binary rows with nonzero variance")
    if not np.isfinite(profiles[valid]).all() or not np.isnan(profiles[~valid]).all():
        raise ValueError("Valid ACF rows must be finite; zero-variance ACF rows must be entirely NaN")
    if not np.allclose(profiles[valid, 0], 1.0, rtol=0, atol=1e-10):
        raise ValueError("Lag 0 differs from one for nonzero-variance reads")
    assert_alignment(binary, metadata, row_ids)
    outputs = [out / "intermediate/acf.npy", out / "intermediate/acf_valid.npy",
               out / "validation/acf_validation.json"]
    np.save(outputs[0], profiles, allow_pickle=False)
    np.save(outputs[1], valid, allow_pickle=False)
    report = {
        "n_reads": len(metadata), "n_valid": int(valid.sum()),
        "n_zero_variance": int((~valid).sum()), "matrix_shape": list(profiles.shape),
        "max_lag_bp": args.max_lag, "lag0_equals_one": True,
        "lag0_max_absolute_error": float(np.max(np.abs(profiles[valid, 0] - 1))) if valid.any() else None,
        "zero_variance_flagged": True, "invalid_rows_all_nan": True,
        "metadata_and_matrix_row_order_preserved": True,
        "parent_function": str(PARENT) + ":autocorrelations",
        "acf_definition": "sum((x[t]-mean(x))*(x[t+k]-mean(x))) / (N*var(x))",
        "binary_orientation": "increasing genomic coordinate; no strand reversal",
    }
    write_json(outputs[2], report)
    finish_stage(out, "02_acf", signature, outputs, details=report)
    log(f"02 ACF: finished; {int(valid.sum()):,} valid and {int((~valid).sum()):,} zero-variance reads")


if __name__ == "__main__":
    main()
