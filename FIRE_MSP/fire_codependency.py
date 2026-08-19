#!/usr/bin/env python3
"""Same-molecule co-actuation (codependency) of FIRE peak pairs across LPS timepoints.

per-timepoint scoring against a fixed peak universe (the LPS_0 filtered peaks) so pair identities match
across the timecourse

Per timepoint, a fiber's state at a peak is:
    covered  : the read fully contains the peak      (FIRE_read_intersect, -F 1)
    actuated : covered AND an MSP >=50%-overlaps it  (MSP_150_FIRE_intersect)
For every peak pair co-spanned by a fiber, counts [neither, first, second, both]
are accumulated over fibers; pairs with >= --min-cov fibers get the example's
score = (obs - p1*p2) * 4, where p1/p2/obs come from the same co-spanning fibers,
plus fs_score, which uses the peaks' overall per-timepoint actuations as the
expectation. Distances are midpoint distances, binned by --bin-len with the
<2*bin_len bins merged into bin 1 (as in the example).

Inputs (from intersect_msp_fire.sh), per sample s:
    <root>/LPS_0/LPS_0.FIRE_peaks_pass_coverage.autosomes.bed.gz   peak universe
    <root>/<s>/<s>.FIRE_read_intersect.vsLPS0peaks.bed.gz          cols 1-4 = peak + read
    <root>/<s>/<s>.MSP_<LEN>_FIRE_intersect.vsLPS0peaks.bed.gz     cols 1-4 = peak + read
(all samples including LPS_0: the .vsLPS0peaks read intersects exclude
secondary/supplementary alignments, which otherwise fabricate distant pairs)

Outputs in --out-dir:
    peak_actuation_by_timepoint.tsv     per-peak nCov/nMSP/actuation x timepoints
    codep_<s>.csv.gz                    per-pair counts + scores, one file per timepoint
    bin_codep_<s>_avg.csv / _median.csv distance-binned score summaries
    codep_paired_across_time.csv.gz     pairs measurable at ALL timepoints, wide scores
    bin_codep_<s>_shuffle_avg.csv       (--shuffle-control) same-marginals independence
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

#############################
# calculate co-dependency stats
#############################
# compares observed joint opening with expected joinit opening if peaks are independent
def codep_score(counts):
    """score: (observed both - p1*p2) * 4, from co-spanning fiber counts."""
    total = sum(counts)
    # prob that peak 1 is actuated: first + both
    p1 = (counts[1] + counts[3]) / total
    # prob that peak 2 is actuated: second + both
    p2 = (counts[2] + counts[3]) / total
    # observed fraction when both peaks are open
    obs = counts[3] / total
    # score > 0 → peaks open together more than expected
    return (obs - p1 * p2) * 4


#############################
# calculate co-actuation score using expected co-actuation rate 
#############################
# Among molecules that span both FIRE peaks, is the observed frequency of both peaks 
# being open higher or lower than expected from their overall accessibility levels at that timepoint?
# counts = [neither, first, second, both]
# exp = expected probability that both peaks are actuated, usually act[i] * act[j]

def codep_score_fsExp(counts, exp):
    """expectation supplied from overall peak actuations."""
    # total number of fibers that span both peaks
    total = sum(counts)
    # observed fraction of co-spanning fibers where both peaks are actuated
    obs = counts[3] / total
    # expectation uses each peak's overall actuation 
    # across all fibers covering that peak at that timepoint
    return (obs - exp) * 4


#############################
# Find read associated with each FIRE peak 
#############################

# for FIRE–read intersect file, maps each FIRE peak to a fixed integer ID, 
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
# Builds co-actuation counts for every FIRE peak pair 
#############################

# For every pair of peaks that occurs on the same spanning fiber, 
# count how many molecules show each of the four possible states

# A closed, B closed -> neither
# A open,   B closed -> first
# A closed, B open   -> second
# A open,   B open   -> both

# 
# {
#     (10, 20): [35, 10, 15, 40],
#     (10, 30): [20, 12, 8, 50]
# }
# 
# 35 fibers = neither
# 10 fibers = peak 10 only
# 15 fibers = peak 20 only
# 40 fibers = both

def count_pairs(read_peaks, actuated, pass_mask, n_peaks, chrom_codes):
    """Accumulate [neither, first, second, both] per co-spanned peak pair.

    Pairs on different chromosomes are skipped: they can only arise from
    residual alignment artifacts, never from one contiguous molecule."""
    
    pair_counts = {}
    
    # loop through every read, only keep peaks that pass the analysis filter
      ## pass_mask = true
    # if a read does not cover 2 peaks, the read is skipped
    for rid, plist in read_peaks.items():
        pl = sorted(p for p in plist if pass_mask[p])
        if len(pl) < 2:
            continue
        
        ##
        # determine whether each covered peak is open
        ##
        
        # generates one T/F value for each peak
        # states = encode (read, peak) combination into one integer
          # rid = 42
          # n_peaks = 100000
          # p = 10
          # states = read 42 + peak 10
        # for the FIRE peaks that are actuated, look for if there is MSP
          # Yes = T
          # No = F
          
        states = [(rid * n_peaks + p) in actuated for p in pl]
        
        # generate every possible pair of peaks on that fiber
        for a in range(len(pl) - 1):
            for b in range(a + 1, len(pl)):
              # generate every possible pair of peaks on that fiber
                if chrom_codes[pl[a]] != chrom_codes[pl[b]]:
                    continue
                # conver the two open/closed states into a category
                # states[a] = False
                # states[b] = False
                
                # cat = 0   neither = 0 + 0
                # cat = 1   first only = 1 + 0
                # cat = 2   second only = 0 + 2
                # cat = 3   both
                
                cat = (1 if states[a] else 0) + (2 if states[b] else 0)
                
                # define the peak pair
                key = (pl[a], pl[b])
                arr = pair_counts.get(key)
                if arr is None:
                    arr = [0, 0, 0, 0]
                    pair_counts[key] = arr
                # cat: 0 neither, 1 first only, 2 second only, 3 both
                arr[cat] += 1
    return pair_counts


#############################
# Convert state counts into co-dependency score
#############################

# takes the raw four-state counts for each FIRE peak pair,
# and turns them into a scored table containing the co-dependency score,
# overall-actuation-based score, genomic distance, and distance bin

def score_pairs(pair_counts, peaks, act, min_cov, bin_len):
    """Score pairs with >= min_cov fibers; returns a DataFrame."""
    
    # take midpoint of every FIRE peak, chr, start, end, fire_coverage / coverage in LPS_0,
    # and peak ID into array
    mids = peaks["mid"].to_numpy()
    chroms = peaks["chrom"].to_numpy()
    starts = peaks["start"].to_numpy()
    ends = peaks["end"].to_numpy()
    # actuation/accessibility value stored in the fixed LPS_0 FIRE peak
    props = peaks["prop_acc"].to_numpy()
    ids = peaks["pid"].to_numpy()

    rows = []
    # loop through every FIRE peak pair in pair counts
    # calculate total number of fibers covering both peaks --> calculate co-dependency score
    for (i, j), c in pair_counts.items():
        total = c[0] + c[1] + c[2] + c[3]
        if total < min_cov:
            continue
        score = codep_score(c)
        # calculate fs_score = expected using each peak's overall actuation
        fs = codep_score_fsExp(c, act[i] * act[j])
        # calculate distance between two peak midpoints
        dist = abs(int(mids[i]) - int(mids[j]))
        rows.append((
            ids[i], chroms[i], starts[i], ends[i], props[i], act[i],
            ids[j], chroms[j], starts[j], ends[j], props[j], act[j],
            c[0], c[1], c[2], c[3], score, fs, dist,
        ))
    # convert rows into tables
      ## motif = peak ID
    df = pd.DataFrame(rows, columns=[
        "motif_1", "chrom_1", "start_1", "end_1", "prop_1", "act_1",
        "motif_2", "chrom_2", "start_2", "end_2", "prop_2", "act_2",
        "neither", "first", "second", "both", "score", "fs_score", "dist",
    ])
    # divide the pair's genomic distance by the bin size
    df["dist_bin"] = np.floor(df["dist"] / bin_len).astype(int)
    df.loc[df["dist_bin"] < 1, "dist_bin"] = 1
    return df


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
    ap.add_argument("--min-act", type=float, default=0.2)
    ap.add_argument("--min-cov", type=int, default=8)
    ap.add_argument("--bin-len", type=int, default=100)
    ap.add_argument("--chrom", default=None, help="restrict to one chromosome (testing)")
    ap.add_argument("--shuffle-control", action="store_true")
    ap.add_argument("--seed", type=int, default=42)
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
        cov_keys = rid.astype(np.int64) * n_peaks + cov["pidx"].to_numpy()

        msp = load_pairs(msp_file, pid_to_idx, args.chrom)
        read_to_rid = dict(zip(reads, range(len(reads))))
        m_rid = msp["read"].map(read_to_rid)
        msp = msp[m_rid.notna()]
        m_keys = m_rid.dropna().astype(np.int64).to_numpy() * n_peaks + msp["pidx"].to_numpy()
        actuated_keys = np.intersect1d(m_keys, cov_keys)  # actuated requires covered

        n_cov = np.bincount(cov["pidx"].to_numpy(), minlength=n_peaks)
        n_msp = np.bincount(actuated_keys % n_peaks, minlength=n_peaks)
        with np.errstate(divide="ignore", invalid="ignore"):
            act = np.where(n_cov > 0, n_msp / np.maximum(n_cov, 1), np.nan)

        read_peaks = {}
        for r, p in zip(rid, cov["pidx"].to_numpy()):
            read_peaks.setdefault(int(r), []).append(int(p))

        tp[s] = dict(read_peaks=read_peaks, actuated=set(actuated_keys.tolist()),
                     n_cov=n_cov, n_msp=n_msp, act=act, n_reads=len(reads))
        print(f"{s}: {len(cov)} covered (fiber,peak) pairs, "
              f"{len(actuated_keys)} actuated, {len(reads)} fibers", flush=True)

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

    ## Pair counting + scoring per timepoint.
    per_t = {}
    for s in args.timepoints:
        pair_counts = count_pairs(tp[s]["read_peaks"], tp[s]["actuated"], pass_mask, n_peaks, chrom_codes)
        df = score_pairs(pair_counts, peaks, tp[s]["act"], args.min_cov, args.bin_len)
        per_t[s] = df
        df.to_csv(out_dir / f"codep_{s}.csv.gz", index=False)
        write_binned(df, args.bin_len,
                     out_dir / f"bin_codep_{s}_avg.csv",
                     out_dir / f"bin_codep_{s}_median.csv")
        print(f"{s}: {len(pair_counts)} co-spanned pairs, "
              f"{len(df)} with >= {args.min_cov} fibers, "
              f"mean score {df['score'].mean():.4f}", flush=True)

        if args.shuffle_control:
            rng = np.random.default_rng(args.seed)
            covering = {}
            for r, plist in tp[s]["read_peaks"].items():
                for p in plist:
                    covering.setdefault(p, []).append(r)
            shuf = set()
            for p, k in enumerate(tp[s]["n_msp"]):
                if k > 0 and pass_mask[p]:
                    for r in rng.choice(covering[p], size=int(k), replace=False):
                        shuf.add(int(r) * n_peaks + p)
            s_counts = count_pairs(tp[s]["read_peaks"], shuf, pass_mask, n_peaks, chrom_codes)
            s_df = score_pairs(s_counts, peaks, tp[s]["act"], args.min_cov, args.bin_len)
            write_binned(s_df, args.bin_len,
                         out_dir / f"bin_codep_{s}_shuffle_avg.csv",
                         out_dir / f"bin_codep_{s}_shuffle_median.csv")
            print(f"{s}: shuffle control mean score {s_df['score'].mean():.4f}", flush=True)

    ## Pairs measurable at every timepoint: wide scores + delta
    first = args.timepoints[0]
    last = args.timepoints[-1]
    wide = None
    for s in args.timepoints:
        d = per_t[s][["motif_1", "motif_2", "dist", "dist_bin",
                      "neither", "first", "second", "both", "score"]].copy()
        d[f"n_{s}"] = d[["neither", "first", "second", "both"]].sum(axis=1)
        d = d.rename(columns={"score": f"score_{s}"})
        d = d[["motif_1", "motif_2", "dist", "dist_bin", f"n_{s}", f"score_{s}"]]
        wide = d if wide is None else wide.merge(
            d.drop(columns=["dist", "dist_bin"]), on=["motif_1", "motif_2"], how="inner")
    if wide is not None and len(args.timepoints) > 1:
        wide["delta_score"] = wide[f"score_{last}"] - wide[f"score_{first}"]
    wide.to_csv(out_dir / "codep_paired_across_time.csv.gz", index=False)
    print(f"pairs measurable at all {len(args.timepoints)} timepoints: {len(wide)}", flush=True)


if __name__ == "__main__":
    main()
