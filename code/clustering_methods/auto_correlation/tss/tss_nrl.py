"""Per-molecule nucleosome repeat length (NRL) and regularity from one ACF row.

Pure functions, no I/O. Every function takes one autocorrelation row ``acf``
(lags 0..max_lag at 1 bp, as returned by the parent ``autocorrelations()``)
and the window width ``n`` the row was computed from.

NRL (SAMOSA-style secondary-peak scan; Abdulhay et al. 2020, eLife):
    lags below ``min_lag`` are ignored to skip the lag-0 shoulder;
    ``scipy.signal.find_peaks`` with ``prominence=min_prominence`` finds local
    maxima; the NRL is the first positive local maximum in [nrl_min, nrl_max].
    The second peak is the first positive local maximum within +/-50 bp of
    2 x NRL that lies inside the available lags.

    ``min_prominence`` is calibrated from the data, not fixed: a per-base
    single-molecule ACF has sampling noise of about 1/sqrt(N) (0.02 for
    N = 2000, 0.034 for N = 500), so a small fixed value such as 0.01 makes
    every random 0/1 row yield a "peak" and puts the first local maximum of
    real reads on a noise bump at the lower band edge. ``calibrate_prominence``
    shuffles the 0/1 positions within each molecule's own row (same m6A
    count, no spatial structure), computes the parent ACF of the shuffled
    rows, takes the largest prominence of any positive local maximum in
    [nrl_min, nrl_max] (lags >= min_lag) per row, and returns the
    ``quantile`` (default 0.95) of these null maxima. By construction at most
    about 5% of structureless rows of that window length pass. Each script
    calibrates at the ACF length it analyses (fixed window in 04b, sliding
    window width in 04c) and records the value; ``--min-prominence`` is an
    optional fixed override.

Regularity from the ACF decay (Baldi et al. 2018, Mol Cell): the parent ACF
divides by N, so even a perfectly periodic signal decays as (N - k) / N. The
decay metrics only use the unbiased rescaling acf[k] * N / (N - k) restricted
to lags <= ``decay_max_lag`` (default floor(N / 2)); the raw ACF is returned
unchanged everywhere else.

Status codes: ``ok``, ``no_peak``, ``zero_variance``, ``window_too_short``,
``peak_negative``.
"""
import math

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

STATUSES = ("ok", "no_peak", "zero_variance", "window_too_short", "peak_negative")
FALLBACK_PERIOD = 190
DEFAULT_PROMINENCE_QUANTILE = 0.95
DEFAULT_N_NULL = 2000
DECAY_BOUND_FACTOR = 10  # decay_length_bp is censored at 10 x decay_max_lag
METRIC_COLUMNS = (
    "nrl_bp", "nrl_status", "nrl_peak_height", "peak2_lag_bp", "peak2_height", "peak_ratio",
    "repeat_peak_lag", "repeat_peak_value", "peak1_height", "decay_length_bp", "fit_period_bp",
    "fit_r2", "fit_converged", "damping_lag_bp",
)


def default_decay_max_lag(n):
    return n // 2


def rescale_unbiased(acf, n, decay_max_lag):
    """acf[k] * N / (N - k) for k <= decay_max_lag (unbiased ACF estimate)."""
    acf = np.asarray(acf, dtype=np.float64)
    if decay_max_lag >= n or decay_max_lag >= len(acf):
        raise ValueError("decay_max_lag must be below the window width and the available lags")
    k = np.arange(decay_max_lag + 1)
    return acf[: decay_max_lag + 1] * (n / (n - k))


def null_peak_prominences(rows, acf_fn, min_lag, nrl_min, nrl_max, seed=0):
    """Largest prominence of a positive in-band local maximum per position-shuffled row.

    Each row's 0/1 values are permuted (same count, no spatial structure); the
    ACF comes from ``acf_fn`` (the parent ``autocorrelations``). Rows with no
    positive in-band peak contribute 0; zero-variance rows are NaN.
    """
    rng = np.random.default_rng(seed)
    rows = np.asarray(rows)
    shuffled = np.stack([rng.permutation(row) for row in rows])
    acf, valid = acf_fn(shuffled)
    maxima = np.full(len(rows), np.nan)
    for i in np.flatnonzero(valid):
        peaks, properties = find_peaks(acf[i, min_lag:], prominence=0)
        peaks = peaks + min_lag
        keep = (peaks >= nrl_min) & (peaks <= nrl_max) & (acf[i, peaks] > 0)
        maxima[i] = float(properties["prominences"][keep].max()) if keep.any() else 0.0
    return maxima


def calibrate_prominence(rows, acf_fn, min_lag=60, nrl_min=120, nrl_max=300,
                         quantile=DEFAULT_PROMINENCE_QUANTILE, n_null=DEFAULT_N_NULL, seed=0):
    """Data-driven min_prominence: the ``quantile`` of the shuffled-row null maxima.

    Up to ``n_null`` rows are drawn (seeded) from ``rows``; the threshold is the
    quantile of ``null_peak_prominences`` over them, so at most about
    (1 - quantile) of structureless rows of this window length yield a peak.
    Returns (threshold, details).
    """
    rows = np.asarray(rows)
    if rows.ndim != 2 or not len(rows):
        raise ValueError("Calibration needs a 2-D 0/1 matrix with at least one row")
    rng = np.random.default_rng(seed)
    if len(rows) > n_null:
        rows = rows[np.sort(rng.choice(len(rows), n_null, replace=False))]
    null = null_peak_prominences(rows, acf_fn, min_lag, nrl_min, nrl_max, seed)
    null = null[np.isfinite(null)]
    if not len(null):
        raise ValueError("Calibration needs at least one nonzero-variance row")
    threshold = float(np.quantile(null, quantile))
    details = {"method": ("quantile of the maximum positive in-band peak prominence of the parent ACF of "
                          "position-shuffled copies of the sampled rows (same m6A count per row)"),
               "window_width": int(rows.shape[1]), "n_null_rows": int(len(null)), "quantile": float(quantile),
               "seed": int(seed), "min_lag": int(min_lag), "band": [int(nrl_min), int(nrl_max)],
               "null_median": float(np.median(null)), "null_q90": float(np.quantile(null, 0.9)),
               "null_q99": float(np.quantile(null, 0.99)), "min_prominence": threshold,
               "expected_noise_sd_1_over_sqrt_n": float(1 / np.sqrt(rows.shape[1]))}
    return threshold, details


def find_nrl(acf, min_lag, min_prominence, nrl_min=120, nrl_max=300):
    """Return (nrl_lag, nrl_height, peak2_lag, peak2_height, status) on the raw ACF.

    ``min_prominence`` must be supplied: use ``calibrate_prominence`` or an explicit override.
    """
    acf = np.asarray(acf, dtype=np.float64)
    if min_prominence is None or not np.isfinite(min_prominence) or min_prominence < 0:
        raise ValueError("min_prominence must be a finite nonnegative number (calibrated or overridden)")
    max_lag = len(acf) - 1
    if not np.isfinite(acf).all():
        return np.nan, np.nan, np.nan, np.nan, "zero_variance"
    if max_lag < nrl_min or max_lag <= min_lag:
        return np.nan, np.nan, np.nan, np.nan, "window_too_short"
    segment = acf[min_lag:]
    peaks, _ = find_peaks(segment, prominence=min_prominence)
    peaks = peaks + min_lag
    in_band = peaks[(peaks >= nrl_min) & (peaks <= nrl_max)]
    if not len(in_band):
        return np.nan, np.nan, np.nan, np.nan, "no_peak"
    positive = in_band[acf[in_band] > 0]
    if not len(positive):
        return np.nan, np.nan, np.nan, np.nan, "peak_negative"
    nrl = int(positive[0])
    second = peaks[(peaks >= 2 * nrl - 50) & (peaks <= 2 * nrl + 50) & (peaks <= max_lag) & (acf[peaks] > 0)]
    if len(second):
        peak2 = int(second[0])
        return nrl, float(acf[nrl]), peak2, float(acf[peak2]), "ok"
    return nrl, float(acf[nrl]), np.nan, np.nan, "ok"


def _damped_cosine(k, amplitude, decay, period, phase, offset):
    return amplitude * np.exp(-k / decay) * np.cos(2 * np.pi * k / period + phase) + offset


def fit_damped_cosine(rescaled, min_lag, decay_max_lag, nrl_min, nrl_max, nrl=np.nan):
    """Fit A*exp(-k/lambda)*cos(2*pi*k/P + phi) + c over [min_lag, decay_max_lag].

    Bounds: P in [nrl_min, nrl_max], lambda in (0, DECAY_BOUND_FACTOR * decay_max_lag].
    P starts at the NRL when found. A decay length at the upper bound means no
    decay is resolvable within the fitted lags (the window is too short to
    distinguish a long decay from none), so it is a censored value, not an
    estimate. Returns (decay_length_bp, fit_period_bp, fit_r2, fit_converged);
    failures give NaN and False, never an exception.
    """
    failed = (np.nan, np.nan, np.nan, False)
    k = np.arange(min_lag, decay_max_lag + 1, dtype=np.float64)
    y = np.asarray(rescaled, dtype=np.float64)[min_lag: decay_max_lag + 1]
    if len(k) < 6 or not np.isfinite(y).all() or np.ptp(y) == 0:
        return failed
    period0 = float(nrl) if np.isfinite(nrl) and nrl_min <= nrl <= nrl_max else float(np.clip(FALLBACK_PERIOD, nrl_min, nrl_max))
    amplitude0 = max(float(np.max(np.abs(y))), 1e-3)
    p0 = [amplitude0, max(decay_max_lag / 2.0, 1.0), period0, 0.0, float(np.mean(y))]
    lower = [0.0, 1e-3, float(nrl_min), -np.pi, -2.0]
    upper = [10.0, float(DECAY_BOUND_FACTOR * decay_max_lag), float(nrl_max), np.pi, 2.0]
    try:
        params, _ = curve_fit(_damped_cosine, k, y, p0=p0, bounds=(lower, upper), maxfev=5000)
    except (RuntimeError, ValueError, TypeError):
        return failed
    if not np.isfinite(params).all():
        return failed
    fitted = _damped_cosine(k, *params)
    residual = float(np.sum((y - fitted) ** 2))
    total = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - residual / total if total > 0 else np.nan
    return float(params[1]), float(params[2]), r2, True


def damping_lag(rescaled, decay_max_lag, period, flat_threshold=0.05):
    """Smallest lag after which the |ACF| envelope stays below flat_threshold.

    The envelope is a centred rolling maximum of |rescaled ACF| over one period,
    evaluated on lags 1..decay_max_lag (lag 0 is identically one and excluded).
    NaN when the envelope is still above the threshold at decay_max_lag.
    """
    values = np.abs(np.asarray(rescaled, dtype=np.float64)[1: decay_max_lag + 1])
    if not len(values) or not np.isfinite(values).all():
        return np.nan
    half = max(int(round(period)) // 2, 0)
    envelope = np.empty(len(values))
    for i in range(len(values)):
        envelope[i] = values[max(0, i - half): i + half + 1].max()
    below = envelope < flat_threshold
    if not below[-1]:
        return np.nan
    # First index from which every later envelope value is below the threshold.
    above = np.flatnonzero(~below)
    first = int(above[-1]) + 1 if len(above) else 0
    return float(first + 1)  # values[0] is lag 1


def molecule_metrics(acf, n, min_prominence, min_lag=60, nrl_min=120, nrl_max=300,
                     decay_max_lag=None, flat_threshold=0.05, repeat_peak=None):
    """All per-molecule metrics for one raw ACF row of a width-n window.

    min_prominence: the calibrated (or overridden) find_peaks prominence.
    repeat_peak: optional parent ``repeat_peak()`` (140-250 bp) for comparison.
    """
    acf = np.asarray(acf, dtype=np.float64)
    if decay_max_lag is None:
        decay_max_lag = default_decay_max_lag(n)
    decay_max_lag = int(min(decay_max_lag, len(acf) - 1, n - 1))
    nrl, nrl_height, peak2, peak2_raw, status = find_nrl(acf, min_lag, min_prominence, nrl_min, nrl_max)
    result = {column: np.nan for column in METRIC_COLUMNS}
    result.update(nrl_bp=nrl, nrl_status=status, nrl_peak_height=nrl_height, peak2_lag_bp=peak2,
                  fit_converged=False)
    if repeat_peak is not None and status != "zero_variance":
        lag, value = repeat_peak(acf)
        result.update(repeat_peak_lag=lag, repeat_peak_value=value)
    if status == "zero_variance":
        return result
    # Decay metrics: unbiased rescaling, lags <= decay_max_lag only.
    rescaled = rescale_unbiased(acf, n, decay_max_lag)
    if status == "ok":
        result["peak1_height"] = float(acf[nrl] * n / (n - nrl))
        if np.isfinite(peak2):
            result["peak2_height"] = float(acf[peak2] * n / (n - peak2))
            result["peak_ratio"] = result["peak2_height"] / result["peak1_height"]
    decay_length, period, r2, converged = fit_damped_cosine(
        rescaled, min_lag, decay_max_lag, nrl_min, nrl_max, nrl)
    result.update(decay_length_bp=decay_length, fit_period_bp=period, fit_r2=r2, fit_converged=converged)
    result["damping_lag_bp"] = damping_lag(
        rescaled, decay_max_lag, nrl if status == "ok" else FALLBACK_PERIOD, flat_threshold)
    return result


def metric_parameters(args):
    """Pick the tss_nrl keyword arguments out of an argparse namespace (min_prominence set later)."""
    return dict(min_lag=args.min_lag, nrl_min=args.nrl_min, nrl_max=args.nrl_max,
                decay_max_lag=args.decay_max_lag, flat_threshold=args.flat_threshold)


def resolve_prominence(args, rows, acf_fn, seed=0):
    """The calibrated threshold for ``rows`` or the explicit ``--min-prominence`` override."""
    if args.min_prominence is not None:
        return float(args.min_prominence), {"method": "fixed override (--min-prominence)",
                                            "min_prominence": float(args.min_prominence),
                                            "window_width": int(np.asarray(rows).shape[1])}
    return calibrate_prominence(rows, acf_fn, args.min_lag, args.nrl_min, args.nrl_max,
                                args.prominence_quantile, args.n_null, seed)


def add_arguments(parser):
    parser.add_argument("--min-lag", type=int, default=60, help="Ignore lags below this (lag-0 shoulder).")
    parser.add_argument("--min-prominence", type=float, default=None,
                        help="Fixed find_peaks prominence override; default: calibrate from shuffled rows.")
    parser.add_argument("--prominence-quantile", type=float, default=DEFAULT_PROMINENCE_QUANTILE,
                        help="Quantile of the shuffled-row null maxima used as min_prominence.")
    parser.add_argument("--n-null", type=int, default=DEFAULT_N_NULL,
                        help="Maximum number of rows shuffled for the calibration.")
    parser.add_argument("--nrl-min", type=int, default=120, help="Lowest NRL lag considered (bp).")
    parser.add_argument("--nrl-max", type=int, default=300, help="Highest NRL lag considered (bp).")
    parser.add_argument("--decay-max-lag", type=int, default=None,
                        help="Largest lag for the decay metrics; default floor(N/2).")
    parser.add_argument("--flat-threshold", type=float, default=0.05,
                        help="|ACF| envelope level that counts as damped.")


def validate_arguments(parser, args):
    if (args.min_lag < 1 or (args.min_prominence is not None and args.min_prominence < 0)
            or not 0.5 <= args.prominence_quantile < 1 or args.n_null < 10
            or args.nrl_min < 1 or args.nrl_max <= args.nrl_min
            or (args.decay_max_lag is not None and args.decay_max_lag < args.min_lag + 5)
            or not math.isfinite(args.flat_threshold) or args.flat_threshold <= 0):
        parser.error("Require min_lag >= 1, prominence >= 0 if given, 0.5 <= prominence_quantile < 1, "
                     "n_null >= 10, 1 <= nrl_min < nrl_max, decay_max_lag >= min_lag + 5 and a positive flat threshold")
