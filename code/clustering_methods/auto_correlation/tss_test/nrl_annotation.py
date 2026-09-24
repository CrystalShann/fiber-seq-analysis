"""Annotation only. Nothing in this module is an input to clustering."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.stats import linregress


@dataclass(frozen=True)
class AnnotationParameters:
    rolling_bp: int = 33
    prominence: float = 0.02
    amplitude: float = 0.02
    same_sign_spacing_bp: int = 100
    min_lag_bp: int = 50
    min_overlap_bp: int = 500
    nucleosome_min_bp: int = 140
    nucleosome_max_bp: int = 250
    min_positive_peaks: int = 3
    min_r_squared: float = 0.90
    max_gap_cv: float = 0.25
    max_residual_period_fraction: float = 0.20


def centered_mean33(signal):
    """Full 33-base windows only: centers 16..1983 for a 2,000-base input.

    No padding, partial windows, or imputation. The 16 edge bases on each
    side have undefined means and are omitted from the annotation ACF.
    """
    x = np.asarray(signal, dtype=float)
    if x.shape != (2000,) or not np.isfinite(x).all():
        raise ValueError("Annotation requires one finite 2-kb signal")
    return np.convolve(x, np.ones(33) / 33, mode="valid")


def annotation_acf(smoothed):
    """Biased linear ACF, sum centered products / sum centered squares.

    Equivalent to statsmodels acf(adjusted=False, fft=False); deliberately
    separate from the imported binary-only clustering ACF implementation.
    """
    centered = smoothed - np.mean(smoothed)
    denominator = centered @ centered
    if np.ptp(smoothed) < 1e-12:
        return np.full(len(smoothed), np.nan)
    return np.correlate(centered, centered, mode="full")[len(centered)-1:] / denominator


def features_from_acf(acf, params=AnnotationParameters()):
    """All qualified extrema in the supported lag range; no count cap.

    Fit consecutive accepted positive peaks to ordinal numbers 1..K with
    a free intercept. Missing cycles are NOT silently imputed. Gap variation,
    residuals and a biological slope band separately flag unreliable fits.
    """
    p = params
    result = dict(nucleosome_peak_lag_bp=np.nan, n_positive_peaks=0,
                  n_extrema=0, nrl_bp=np.nan, regression_intercept_bp=np.nan,
                  regression_r_squared=np.nan, peak_gap_cv=np.nan,
                  max_residual_period_fraction=np.nan, extrema_abs_sum=np.nan,
                  regularity=np.nan, reliable_nrl=False,
                  annotation_status="zero_variance", positive_peak_lags=[],
                  negative_peak_lags=[], accepted_extrema_lags=[])
    if not np.isfinite(acf).all():
        return result
    last = len(acf) - p.min_overlap_bp
    # Search the full ACF before lag filtering: eligible boundary extrema
    # still have their actual neighbors. Global prominence (wlen=None).
    positive, _ = find_peaks(acf, prominence=p.prominence,
                            distance=p.same_sign_spacing_bp, height=p.amplitude)
    negative, _ = find_peaks(-acf, prominence=p.prominence,
                            distance=p.same_sign_spacing_bp, height=p.amplitude)
    positive = positive[(positive >= p.min_lag_bp) & (positive <= last)]
    negative = negative[(negative >= p.min_lag_bp) & (negative <= last)]
    extrema = np.sort(np.r_[positive, negative])
    band = positive[(positive >= p.nucleosome_min_bp) & (positive <= p.nucleosome_max_bp)]
    result.update(annotation_status="too_few_peaks", n_positive_peaks=len(positive),
                  n_extrema=len(extrema), positive_peak_lags=positive.tolist(),
                  negative_peak_lags=negative.tolist(), accepted_extrema_lags=extrema.tolist(),
                  extrema_abs_sum=float(np.abs(acf[extrema]).sum()),
                  regularity=float(np.abs(acf[extrema]).mean()) if len(extrema) else np.nan)
    if len(band):
        result["nucleosome_peak_lag_bp"] = int(band[np.argmax(acf[band])])
    if len(positive) >= 2:
        numbers = np.arange(1, len(positive) + 1)
        fit = linregress(numbers, positive)
        gaps = np.diff(positive)
        residual = np.max(np.abs(positive - (fit.intercept + fit.slope * numbers))) / fit.slope
        gap_cv = float(np.std(gaps) / np.mean(gaps))
        result.update(nrl_bp=float(fit.slope), regression_intercept_bp=float(fit.intercept),
                      regression_r_squared=float(fit.rvalue**2), peak_gap_cv=gap_cv,
                      max_residual_period_fraction=float(residual))
        reasons = []
        if len(positive) < p.min_positive_peaks:
            reasons.append("too_few_peaks")
        if not len(band):
            reasons.append("no_nucleosome_peak")
        if not p.nucleosome_min_bp <= fit.slope <= p.nucleosome_max_bp:
            reasons.append("period_out_of_band")
        if fit.rvalue**2 < p.min_r_squared:
            reasons.append("low_r_squared")
        if gap_cv > p.max_gap_cv:
            reasons.append("irregular_peak_spacing")
        if residual > p.max_residual_period_fraction:
            reasons.append("large_residual")
        result.update(reliable_nrl=not reasons, annotation_status=";".join(reasons) or "reliable")
    return result


def annotate(binary, params=AnnotationParameters()):
    if params.rolling_bp != 33:
        raise ValueError("This workflow fixes annotation smoothing at 33 bp")
    return pd.DataFrame([features_from_acf(annotation_acf(centered_mean33(row)), params)
                         for row in binary])
