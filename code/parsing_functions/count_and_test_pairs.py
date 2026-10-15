"""Shared co-accessibility statistical tests; no analysis runs on import."""

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio as conditional_odds_ratio

def count_pairs(pairs, coverage, accessible):
    """Return integer counts: neither, CRE1-only, CRE2-only, both.

    Coverage and accessibility map CRE identifiers to sets of fiber IDs.
    """
    counts = np.zeros((len(pairs), 4), dtype=np.int64)
    empty = set()
    # loop over the pairs, and pair.CRE1 and pair.CRE2 are CRE ids
    for i, pair in enumerate(pairs.itertuples(index=False)):
        # shared reads = reads covering both CREs
        shared = coverage.get(pair.CRE1, empty) & coverage.get(pair.CRE2, empty)
        # a = shared reads where CRE1 is acessible, b = shared reads where CRE2 is accessible
        a = accessible.get(pair.CRE1, empty) & shared
        b = accessible.get(pair.CRE2, empty) & shared
        both = len(a & b)
        # co_closed (neither accessible) is shared − |a| − |b| + both
          # Adding both back is needed because reads in both a and b were subtracted twice
        # CRE1_access (CRE1 only) is |a| − both
        # CRE2_access (CRE2 only) is |b| − both
        # co_access (both) is both
        
        counts[i] = len(shared) - len(a) - len(b) + both, len(a) - both, len(b) - both, both
    return counts

def bh(p):
    """BH adjustment of one complete family of finite p-values."""
    p = np.asarray(p, dtype=float)
    if not len(p):
        return p.copy()
    # sort p values from smallest to largest
    order = np.argsort(p)
    # takes each sorted p-value and multiplies it by n / rank. The smallest p gets rank 1 and the largest gets rank n
    # [::-1] largest p values come first
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    # Running minimum from the largest ranked p-value.
    out = np.empty(len(p)); out[order] = np.minimum(adjusted, 1)
    return out

def test_counts(counts, pseudocount):
    """Two-sided Fisher p-values and conditional ORs on pooled count tables.

    Input columns: neither, CRE1-only, CRE2-only, both. Add the explicit
    pseudocount after pooling the caller's cohort; also return raw p-values.
    Identical count tables are tested once.
    """
    # only keep distinct read counts to run fisher on
    unique, inverse = np.unique(counts, axis=0, return_inverse=True)
    rows = []
    for closed, first, second, both in unique:
        # raw count
        raw = np.array([[closed, second], [first, both]], dtype=np.int64)
        # +1 to every count to avoid division by -
        corrected = raw + pseudocount
        # pval is the two-sided Fisher p-value on the +1 table
        # pval_raw is the two-sided Fisher p-value on the raw table,
        # fisher_estimate is the conditional maximum-likelihood odds ratio of the +1 table
        rows.append((fisher_exact(corrected, alternative="two-sided").pvalue,
                     conditional_odds_ratio(corrected, kind="conditional").statistic,
                     fisher_exact(raw, alternative="two-sided").pvalue))
    stats = np.array(rows, dtype=float)[inverse] if rows else np.empty((0, 3))
    return pd.DataFrame(stats, columns=["pval", "fisher_estimate", "pval_raw"])

def summarize_coaccess_counts(result, pseudocount=1):
    """Add counts, effect sizes, expectations and BH columns to a pair table.

    The caller supplies the complete BH family. Tests must already be present.
    """
    result = result.copy()
    pc = pseudocount
    counts = result[["co_closed", "CRE1_access", "CRE2_access", "co_access"]].to_numpy()
    result["n_shared_reads"] = counts.sum(axis=1)
    result["co_inaccess"] = result.co_closed
    with np.errstate(divide="ignore", invalid="ignore"):
        result["OR"] = ((result.co_access.astype(float) + pc) * (result.co_closed + pc)) / ((result.CRE1_access + pc) * (result.CRE2_access + pc))
        result["or_haldane"] = ((result.co_access + 0.5) * (result.co_closed + 0.5)) / ((result.CRE1_access + 0.5) * (result.CRE2_access + 0.5))
    a = (result.co_access + result.CRE1_access).astype(float)
    b = result.co_access + result.CRE2_access
    result["expected_co_access"] = a * b / result.n_shared_reads
    result["observed_expected_ratio"] = result.co_access / result.expected_co_access.replace(0, np.nan)
    result["excess_fraction"] = (result.co_access - result.expected_co_access) / result.n_shared_reads
    result["zero_access_cre1"] = a == 0
    result["zero_access_cre2"] = b == 0
    result["fdr"] = bh(result.pval)
    result["fdr_raw"] = bh(result.pval_raw)
    return result
