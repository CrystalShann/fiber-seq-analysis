#!/usr/bin/env python3
"""Same-molecule co-actuation (codependency) of FIRE peak pairs across LPS timepoints.

Global analysis: every FIRE peak pair co-spanned by the same fiber, genome-wide,
scored separately per timepoint against a fixed peak universe (the LPS_0
filtered peaks) so pair identities match across the timecourse. Each timepoint
is tested by itself; timepoints are not tested against a baseline.

Per timepoint, a fiber's state at a peak is:
    covered  : the read fully contains the peak      (FIRE_read_intersect, -F 1)
    actuated : covered AND an MSP >=50%-overlaps it  (MSP_150_FIRE_intersect)
For every peak pair co-spanned by a fiber, counts [neither, first, second, both]
are accumulated over fibers; pairs where each peak is actuated on at least one
of the co-spanning fibers (no coverage minimum by default; see --min-cov) get
    score      = (obs - p1*p2) * 4      p1/p2/obs from the same co-spanning fibers
    fs_score   = (obs - act1*act2) * 4  expectation from overall peak actuations
    odds_ratio = ((both+1)*(neither+1)) / ((first+1)*(second+1))
                 Haldane-Anscombe; robust to the marginal actuation shifts
                 between timepoints, so pair strength stays comparable
    fisher_p   = two-sided Fisher exact test of the 2x2 table with a +1
                 pseudocount in every cell: the null is that the two peaks'
                 open/closed states are independent across the co-spanning
                 fibers, given each peak's marginal rate
Distances are midpoint distances, binned by --bin-len with the <2*bin_len bins
merged into bin 1.

Inputs (from intersect_msp_fire.sh), per sample s:
    <root>/LPS_0/LPS_0.FIRE_peaks_pass_coverage.autosomes.bed.gz   peak universe
    <root>/<s>/<s>.FIRE_read_intersect.vsLPS0peaks.bed.gz          cols 1-4 = peak + read
    <root>/<s>/<s>.MSP_<LEN>_FIRE_intersect.vsLPS0peaks.bed.gz     cols 1-4 = peak + read
(all samples including LPS_0: the .vsLPS0peaks read intersects exclude
secondary/supplementary alignments, which otherwise fabricate distant pairs)

Outputs in --out-dir:
    peak_actuation_by_timepoint.tsv     per-peak nCov/nMSP/actuation x timepoints
    codep_<s>.csv.gz                    per-pair counts + scores + fisher_p
    bin_codep_<s>_avg.csv / _median.csv distance-binned score summaries
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

#############################
# calculate co-dependency stats
#############################
# compares observed joint opening with expected joint opening if peaks are independent
def codep_scores(counts):
    """score: (observed both - p1*p2) * 4, per pair (counts = n_pairs x 4)."""
    total = counts.sum(axis=1)
    # prob that peak 1 is actuated: first + both
    p1 = (counts[:, 1] + counts[:, 3]) / total
    # prob that peak 2 is actuated: second + both
    p2 = (counts[:, 2] + counts[:, 3]) / total
    # observed fraction when both peaks are open
    obs = counts[:, 3] / total
    # score > 0 -> peaks open together more than expected
    return (obs - p1 * p2) * 4


#############################
# calculate co-actuation score using expected co-actuation rate
#############################
# Among molecules that span both FIRE peaks, is the observed frequency of both peaks
# being open higher or lower than expected from their overall accessibility levels at that timepoint?
# counts = n_pairs x [neither, first, second, both]
# exp = expected probability that both peaks are actuated, usually act[i] * act[j]

def codep_scores_fsExp(counts, exp):
    """expectation supplied from overall peak actuations."""
    total = counts.sum(axis=1)
    # observed fraction of co-spanning fibers where both peaks are actuated
    obs = counts[:, 3] / total
    return (obs - exp) * 4


#############################
# Find read associated with each FIRE peak
#############################

# for FIRE-read intersect file, maps each FIRE peak to a fixed integer ID,
# checks that every peak belongs to the common LPS_0 peak universe

def load_pairs(path, pid_to_idx, chrom=None):
    """Read an intersect bed (cols 1-4 = peak chrom/start/end + read name) into
    a deduplicated DataFrame of (pidx, read)."""
    # use the LPS_0.FIRE_read_intersect.vsLPS0peaks.bed.gz
      ## Reads that fully cover the FIRE peak = denominator
    # directory mapping peak coordinates to integer index
    # option to map one selected chr
    df = pd.read_csv(
        path, sep="\t", header=None, usecols=[0, 1, 2, 3],
        names=["chrom", "start", "end", "read"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64, "read": str},
    )
    if chrom:
        df = df[df["chrom"] == chrom]

    # create a unique string ID for each peak from its coordinates = peak ID
      ## chr1:1000-1200
    # convert peak ID into its integer index
    pid = df["chrom"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
    df["pidx"] = pid.map(pid_to_idx)
    if df["pidx"].isna().any():
        n = int(df["pidx"].isna().sum())
        raise ValueError(f"{path}: {n} rows reference peaks absent from the universe")
    df["pidx"] = df["pidx"].astype(np.int64)
    # only keep peak ID and read name
    return df[["pidx", "read"]].drop_duplicates()


#############################
# Vectorized co-spanned pair machinery
#############################
# All co-spanned peak pairs are enumerated as flat index arrays; counting the
# four molecule states then reduces to boolean indexing + one bincount.

def co_spanned_index_pairs(rr, pp, chrom_codes):
    """All (a, b) index pairs into rr/pp where one fiber covers two peaks on
    the same chromosome. rr must be grouped by fiber with pp ascending within
    each fiber, so pp[a] < pp[b] for every pair (matches the old sorted(pl)).

    Pairs on different chromosomes are skipped: they can only arise from
    residual alignment artifacts, never from one contiguous molecule."""
    n = len(rr)
    if n == 0:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    # fiber group boundaries
    starts = np.flatnonzero(np.r_[True, rr[1:] != rr[:-1]])
    counts = np.diff(np.r_[starts, n])
    # for each group size k, expand a precomputed C(k,2) combination template
    a_parts, b_parts = [], []
    for k in np.unique(counts):
        if k < 2:
            continue
        ia, ib = np.triu_indices(int(k), 1)
        s = starts[counts == k]
        a_parts.append((s[:, None] + ia[None, :]).ravel())
        b_parts.append((s[:, None] + ib[None, :]).ravel())
    if not a_parts:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    idx_a = np.concatenate(a_parts)
    idx_b = np.concatenate(b_parts)
    same_chrom = chrom_codes[pp[idx_a]] == chrom_codes[pp[idx_b]]
    return idx_a[same_chrom], idx_b[same_chrom]


def pair_state_counts(pair_idx, n_pairs, state_a, state_b):
    """[neither, first, second, both] per pair from per-incidence boolean states.

    cat = 0 neither, 1 first only, 2 second only, 3 both"""
    cat = state_a.astype(np.int64) + 2 * state_b.astype(np.int64)
    return np.bincount(pair_idx * 4 + cat, minlength=n_pairs * 4).reshape(n_pairs, 4)


def write_binned(df, bin_len, out_avg, out_median):
    with_scores = df[df["score"].notna()]
    for stat, out in (("mean", out_avg), ("median", out_median)):
        binned = with_scores.groupby("dist_bin")["score"].agg(stat).to_frame("score")
        binned["dist_bin"] = binned.index
        binned["dist"] = binned["dist_bin"] * bin_len
        binned["count"] = with_scores["dist_bin"].value_counts().sort_index()
        binned.to_csv(out, index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fire-msp-root", default="/project/spott/cshan/fiber-seq/macrophage_project/FIRE_MSP")
    ap.add_argument("--out-dir", default=None, help="default: <root>/codependency")
    ap.add_argument("--timepoints", nargs="+", default=["LPS_0", "LPS_5", "LPS_10", "LPS_15"])
    ap.add_argument("--len", type=int, default=150, dest="msp_len")
    ap.add_argument("--min-act", type=float, default=0.0)
    ap.add_argument("--min-cov", type=int, default=0,
                    help="minimum co-spanning fibers per pair (0 = no minimum)")
    ap.add_argument("--bin-len", type=int, default=100)
    ap.add_argument("--chrom", default=None, help="restrict to one chromosome (testing)")
    args = ap.parse_args()

    root = Path(args.fire_msp_root)
    out_dir = Path(args.out_dir) if args.out_dir else root / "codependency"
    out_dir.mkdir(parents=True, exist_ok=True)

    ## Peak universe: the LPS_0 filtered peaks.
    peaks_file = root / "LPS_0" / "LPS_0.FIRE_peaks_pass_coverage.autosomes.bed.gz"
    peaks = pd.read_csv(
        peaks_file, sep="\t", header=None,
        names=["chrom", "start", "end", "name", "prop_acc", "strand"],
    )
    if args.chrom:
        peaks = peaks[peaks["chrom"] == args.chrom].reset_index(drop=True)
    peaks["pid"] = peaks["chrom"] + ":" + peaks["start"].astype(str) + "-" + peaks["end"].astype(str)
    peaks["mid"] = peaks["start"] + ((peaks["end"] - peaks["start"]) / 2).round().astype(int)
    n_peaks = len(peaks)
    pid_to_idx = dict(zip(peaks["pid"], range(n_peaks)))
    chrom_codes = pd.factorize(peaks["chrom"])[0]
    print(f"peak universe: {n_peaks} LPS_0 peaks" + (f" on {args.chrom}" if args.chrom else ""), flush=True)

    mids = peaks["mid"].to_numpy()
    chroms_arr = peaks["chrom"].to_numpy()
    starts_arr = peaks["start"].to_numpy()
    ends_arr = peaks["end"].to_numpy()
    props_arr = peaks["prop_acc"].to_numpy()
    ids_arr = peaks["pid"].to_numpy()

    ## Load per-timepoint fiber states.
    tp = {}
    for s in args.timepoints:
        cov_file = root / s / f"{s}.FIRE_read_intersect.vsLPS0peaks.bed.gz"
        msp_file = root / s / f"{s}.MSP_{args.msp_len}_FIRE_intersect.vsLPS0peaks.bed.gz"
        for f in (cov_file, msp_file):
            if not f.is_file():
                sys.exit(f"ERROR: missing input {f}")

        cov = load_pairs(cov_file, pid_to_idx, args.chrom)
        reads, rid = np.unique(cov["read"].to_numpy(), return_inverse=True)
        pidx = cov["pidx"].to_numpy()
        cov_keys = rid.astype(np.int64) * n_peaks + pidx

        msp = load_pairs(msp_file, pid_to_idx, args.chrom)
        read_to_rid = dict(zip(reads, range(len(reads))))
        m_rid = msp["read"].map(read_to_rid)
        msp = msp[m_rid.notna()]
        m_keys = m_rid.dropna().astype(np.int64).to_numpy() * n_peaks + msp["pidx"].to_numpy()
        # boolean per covered (fiber,peak) cell: actuated requires covered
        act_bool = np.isin(cov_keys, m_keys)

        n_cov = np.bincount(pidx, minlength=n_peaks)
        n_msp = np.bincount(pidx[act_bool], minlength=n_peaks)
        with np.errstate(divide="ignore", invalid="ignore"):
            act = np.where(n_cov > 0, n_msp / np.maximum(n_cov, 1), np.nan)

        tp[s] = dict(rid=rid.astype(np.int64), pidx=pidx, act_bool=act_bool,
                     n_cov=n_cov, n_msp=n_msp, act=act)
        print(f"{s}: {len(cov)} covered (fiber,peak) pairs, "
              f"{int(act_bool.sum())} actuated, {len(reads)} fibers", flush=True)

    ## Per-peak actuation table + min_act mask
    act_tbl = peaks[["pid", "chrom", "start", "end", "prop_acc"]].copy()
    for s in args.timepoints:
        act_tbl[f"{s}_nCov"] = tp[s]["n_cov"]
        act_tbl[f"{s}_nMSP"] = tp[s]["n_msp"]
        act_tbl[f"{s}_act"] = tp[s]["act"]
    act_mat = act_tbl[[f"{s}_act" for s in args.timepoints]].to_numpy()
    with np.errstate(invalid="ignore"):
        max_act = np.nanmax(act_mat, axis=1)
    act_tbl["max_act"] = max_act
    pass_mask = np.nan_to_num(max_act, nan=0.0) >= args.min_act
    act_tbl["pass_min_act"] = pass_mask
    act_tbl.to_csv(out_dir / "peak_actuation_by_timepoint.tsv", sep="\t", index=False)
    print(f"peaks passing min_act >= {args.min_act} (max over timepoints): "
          f"{int(pass_mask.sum())}/{n_peaks}", flush=True)

    ## Pair counting + scoring + Fisher testing, each timepoint by itself.
    for s in args.timepoints:
        d = tp[s]
        keep_cell = pass_mask[d["pidx"]]
        rr, pp, aa = d["rid"][keep_cell], d["pidx"][keep_cell], d["act_bool"][keep_cell]
        # group cells by fiber, peaks ascending within each fiber
        order = np.lexsort((pp, rr))
        rr, pp, aa = rr[order], pp[order], aa[order]

        idx_a, idx_b = co_spanned_index_pairs(rr, pp, chrom_codes)
        pair_key = pp[idx_a].astype(np.int64) * n_peaks + pp[idx_b]
        uniq_keys, pair_idx = np.unique(pair_key, return_inverse=True)
        counts_all = pair_state_counts(pair_idx, len(uniq_keys), aa[idx_a], aa[idx_b])

        # a pair is scored only with >= min_cov co-spanning fibers AND each peak
        # actuated on at least one of them; pairs where either peak is never
        # actuated among the shared fibers are dropped
        totals = counts_all.sum(axis=1)
        act1 = counts_all[:, 1] + counts_all[:, 3]
        act2 = counts_all[:, 2] + counts_all[:, 3]
        kept = (totals >= args.min_cov) & (act1 > 0) & (act2 > 0)
        n_kept = int(kept.sum())

        counts = counts_all[kept]
        pi = (uniq_keys[kept] // n_peaks).astype(np.int64)
        pj = (uniq_keys[kept] % n_peaks).astype(np.int64)

        score = codep_scores(counts)
        act_arr = d["act"]
        fs = codep_scores_fsExp(counts, act_arr[pi] * act_arr[pj])
        # Haldane-Anscombe odds ratio: margin-robust effect size
        odds = ((counts[:, 3] + 1.0) * (counts[:, 0] + 1.0)) \
             / ((counts[:, 1] + 1.0) * (counts[:, 2] + 1.0))
        dist = np.abs(mids[pi] - mids[pj]).astype(np.int64)
        dist_bin = np.maximum(dist // args.bin_len, 1).astype(int)

        df = pd.DataFrame({
            "motif_1": ids_arr[pi], "chrom_1": chroms_arr[pi],
            "start_1": starts_arr[pi], "end_1": ends_arr[pi],
            "prop_1": props_arr[pi], "act_1": act_arr[pi],
            "motif_2": ids_arr[pj], "chrom_2": chroms_arr[pj],
            "start_2": starts_arr[pj], "end_2": ends_arr[pj],
            "prop_2": props_arr[pj], "act_2": act_arr[pj],
            "neither": counts[:, 0], "first": counts[:, 1],
            "second": counts[:, 2], "both": counts[:, 3],
            "score": score, "fs_score": fs, "odds_ratio": odds,
            "dist": dist, "dist_bin": dist_bin,
        })

        ## Per-pair significance: two-sided Fisher exact test on the 2x2 table
        ## [[neither, second], [first, both]] with a +1 pseudocount in every
        ## cell. Identical tables share one test, so only unique tables are
        ## computed.
        uniq_tab, tab_inv = np.unique(counts, axis=0, return_inverse=True)
        pvals = np.empty(len(uniq_tab))
        for t, (ne, fi, se, bo) in enumerate(uniq_tab):
            pvals[t] = fisher_exact([[ne + 1, se + 1], [fi + 1, bo + 1]])[1]
        df["fisher_p"] = pvals[tab_inv] if n_kept > 0 else np.nan

        print(f"{s}: {len(uniq_keys)} co-spanned pairs, {n_kept} testable "
              f"(both peaks actuated on shared fibers), "
              f"mean score {df['score'].mean():.4f}, "
              f"{int((df['fisher_p'] < 0.05).sum())} pairs with fisher_p < 0.05",
              flush=True)

        df.to_csv(out_dir / f"codep_{s}.csv.gz", index=False)
        write_binned(df, args.bin_len,
                     out_dir / f"bin_codep_{s}_avg.csv",
                     out_dir / f"bin_codep_{s}_median.csv")


if __name__ == "__main__":
    main()
