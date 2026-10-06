#!/usr/bin/env python3
"""Synthetic-signal tests for tss_nrl.py; CPU only, run with python or pytest.

Signals are 0/1 m6A rows built from 147 bp footprints (low call rate) and
40 bp linkers (high call rate), so the true NRL is 187 bp. ACFs come from the
parent ``autocorrelations()`` at every 1 bp lag, and ``min_prominence`` is
calibrated from position-shuffled copies of the tested rows exactly as the
pipeline scripts do.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tss_nrl  # noqa: E402
from tss_common import load_parent  # noqa: E402

PARENT = load_parent("03_compute_autocorrelations.py")
ACF = PARENT.autocorrelations
NUC, LINK = 147, 40
TRUE_NRL = NUC + LINK


def regular_array(rng, width, jitter_sd=0.0, p_nuc=0.03, p_link=0.5):
    """One binary row: random phase, footprint/linker alternation, optional linker jitter."""
    row = np.zeros(width, dtype=np.uint8)
    position = -int(rng.integers(0, TRUE_NRL))
    while position < width:
        linker = LINK if jitter_sd == 0 else int(max(5, round(rng.normal(LINK, jitter_sd))))
        for length, p in ((NUC, p_nuc), (linker, p_link)):
            lo, hi = max(position, 0), min(position + length, width)
            if hi > lo:
                row[lo:hi] = rng.random(hi - lo) < p
            position += length
    return row


def metrics_for(rows, min_prominence=None, **params):
    rows = np.asarray(rows)
    if min_prominence is None:
        min_prominence, _ = tss_nrl.calibrate_prominence(rows, ACF)
    acf, valid = ACF(rows)
    assert acf.shape == (len(rows), rows.shape[1])
    return [tss_nrl.molecule_metrics(acf[i], rows.shape[1], min_prominence, **params) for i in range(len(rows))]


def test_calibration_tracks_noise_floor():
    rng = np.random.default_rng(0)
    thresholds = {}
    for width in (2000, 500):
        rows = (rng.random((300, width)) < 0.3).astype(np.uint8)
        thresholds[width], details = tss_nrl.calibrate_prominence(rows, ACF, seed=0)
        assert details["n_null_rows"] == 300 and details["window_width"] == width
        assert 0 < thresholds[width] < 1
    # Noise scales as 1/sqrt(N): the shorter window needs the larger prominence.
    assert thresholds[500] > thresholds[2000], thresholds


def test_regular_array_recovers_nrl():
    rng = np.random.default_rng(1)
    for width in (2000, 500):
        rows = np.stack([regular_array(rng, width) for _ in range(60)])
        found = [m["nrl_bp"] for m in metrics_for(rows) if m["nrl_status"] == "ok"]
        assert len(found) >= 0.9 * len(rows), f"width {width}: only {len(found)} peaks"
        median = float(np.median(found))
        assert abs(median - TRUE_NRL) <= 5, f"width {width}: median NRL {median}"
        assert np.mean(np.abs(np.array(found) - TRUE_NRL) <= 5) >= 0.5, f"width {width}: NRL spread"


def test_decay_length_decreases_with_jitter():
    rng = np.random.default_rng(2)
    reference = np.stack([regular_array(rng, 2000) for _ in range(150)])
    threshold, _ = tss_nrl.calibrate_prominence(reference, ACF)
    medians = []
    for jitter in (0, 10, 20, 40):
        rows = np.stack([regular_array(rng, 2000, jitter_sd=jitter) for _ in range(150)])
        decay = [m["decay_length_bp"] for m in metrics_for(rows, threshold) if np.isfinite(m["decay_length_bp"])]
        assert len(decay) >= 100, f"jitter {jitter}: too few converged fits ({len(decay)})"
        medians.append(float(np.median(decay)))
    assert all(a > b for a, b in zip(medians, medians[1:])), f"decay medians not decreasing: {medians}"


def test_random_signal_mostly_no_peak():
    rng = np.random.default_rng(3)
    rows = (rng.random((200, 2000)) < 0.3).astype(np.uint8)
    statuses = [m["nrl_status"] for m in metrics_for(rows)]
    fraction = statuses.count("no_peak") / len(statuses)
    assert fraction > 0.5, f"random rows: only {fraction:.2f} no_peak"


def test_constant_row_zero_variance():
    for value in (0, 1):
        rows = np.full((2, 2000), value, dtype=np.uint8)
        statuses = [m["nrl_status"] for m in metrics_for(rows, 0.1)]
        assert statuses == ["zero_variance", "zero_variance"]
        assert all(np.isnan(m["decay_length_bp"]) for m in metrics_for(rows, 0.1))


def test_rescaling_and_lag_count():
    acf = np.linspace(1, 0, 11)
    rescaled = tss_nrl.rescale_unbiased(acf, 10, 5)
    assert len(rescaled) == 6 and np.isclose(rescaled[5], acf[5] * 10 / 5)
    rng = np.random.default_rng(4)
    rows = np.stack([regular_array(rng, 500) for _ in range(3)])
    acf, _ = ACF(rows, n_features=500)
    assert acf.shape[1] == 500 and np.allclose(acf[:, 0], 1)
    assert tss_nrl.find_nrl(acf[0][:100], 60, 0.1)[-1] == "window_too_short"


if __name__ == "__main__":
    for name, test in sorted(globals().items()):
        if name.startswith("test_") and callable(test):
            test()
            print(f"PASS {name}", flush=True)
    print("All tss_nrl tests passed")
