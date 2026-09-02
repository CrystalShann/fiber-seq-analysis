#!/usr/bin/env python3
"""Per-read m6A fractions over the 1000 bp entropy windows.

For every window in entropy_windows_tss.tsv (TSS +/- 500) and
entropy_windows_fire.tsv (FIRE-region midpoint +/- 500, both from
01_entropy_windows.R), and every FULLY-SPANNING read (aligned span covering
the whole window), per timepoint:

    n_sites   reference A/T positions in the window (m6A-callable sites)
    n_mod     m6A calls landing on those sites
    frac      n_mod / n_sites

This is the per-read "average m6A density" of the SAM-seq accessibility
heterogeneity analysis (Leduque et al. 2024); 03_entropy_violin.R pools the
rows across timepoints and computes the 4-bin Shannon entropy per window.

Parsing follows 02_tss_m6a_profiles.py / 04_m6a_binned.py: the first and last
BED12 blocks are ft extract sentinels (dropped), split reads are deduped
keeping the longest alignment per read name, and calls in the window but not
on a reference A/T are counted off-site and excluded. Because reads must span
the whole window, covered A/T sites == all window A/T sites, so every read of
a window shares the same denominator. Reads with zero calls in the window
still get a row (frac 0) -- they are the closed chromatin the entropy needs.

Outputs in --out-dir:
    read_fractions_tss.tsv.gz    one row per (gene window, timepoint, read)
    read_fractions_fire.tsv.gz   one row per (FIRE window, timepoint, read)
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

warnings.filterwarnings("ignore", category=DeprecationWarning)

SAMPLES = ["LPS_0", "LPS_5", "LPS_10", "LPS_15"]
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

FT_ROOT = "/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/ft_result_dir"
REF_FA = "/project/spott/reference/human/GRCh38/hg38.fa"
TAB_DIR = "/project/spott/cshan/fiber-seq/macrophage_project/clustering_TSS_methods/entropy/tables"

WIN = 1000   # every window is exactly this wide (asserted); the containment
             # searchsorted in fractions_for_chrom relies on the constant width


def m6a_path(ft_root, sample, chrom):
    return Path(ft_root) / sample / "extracted_results" / "m6a_by_chr" / \
        f"{sample}.ft_extracted_m6a.{chrom}.bed.gz"


def merge_intervals(starts, ends):
    """Merge sorted, possibly overlapping [start, end) windows for fetching."""
    out = []
    for s, e in zip(starts, ends):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


###############################
# Per-chromosome A/T-site arrays (shared by all timepoints)
###############################
def build_sites(sub, chrom_seq):
    """Per-window sorted reference A/T positions, concatenated with offsets.
    Windows may overlap (bidirectional promoters), so sites are duplicated per
    window; within one window they stay sorted by genomic position."""
    ws = sub["win_start"].to_numpy()
    we = sub["win_end"].to_numpy()
    sites_l = []
    for r in range(len(sub)):
        seg = np.frombuffer(chrom_seq[ws[r]:we[r]].encode(), dtype=np.uint8)
        sites_l.append(ws[r] + np.nonzero((seg == ord("A")) | (seg == ord("T")))[0])
    nsite = np.array([a.size for a in sites_l], dtype=np.int64)
    return {
        "ws": ws, "we": we,
        "sites": np.concatenate(sites_l) if sites_l else np.empty(0, np.int64),
        "offsets": np.concatenate([[0], np.cumsum(nsite)]),
        "n_sites": nsite,
    }


###############################
# Pass 1: longest alignment per read name over the merged windows
###############################
def survey_reads(tb, chrom, merged):
    best = {}
    for s, e in merged:
        for line in tb.fetch(chrom, s, e):
            f = line.split("\t", 4)
            st = int(f[1])
            ln = int(f[2]) - st
            rid = f[3]
            cur = best.get(rid)
            if cur is None or ln > cur[0]:
                best[rid] = (ln, st)
    return best


###############################
# Pass 2: (window, read, n_mod) for one (chromosome, timepoint)
###############################
def fractions_for_chrom(tb, chrom, merged, best, b):
    ws, we = b["ws"], b["we"]
    sites, offsets = b["sites"], b["offsets"]
    rows = []
    off_site = 0
    n_kept = 0

    for ms, me in merged:
        for line in tb.fetch(chrom, ms, me):
            f = line.split("\t")
            st = int(f[1])
            en = int(f[2])
            rid = f[3]
            if best.get(rid) != (en - st, st):
                continue        # shorter split alignment, or already counted
            best[rid] = None
            n_kept += 1

            # windows fully inside [st, en): ws >= st and ws <= en - WIN
            r0 = np.searchsorted(ws, st, side="left")
            r1 = np.searchsorted(ws, en - WIN, side="right")
            if r1 <= r0:
                continue

            # m6A calls: drop the sentinel first/last blocks; 1 bp blocks, so
            # the modified base (0-based) is span start + blockStart (sorted)
            starts = np.fromstring(f[11], dtype=np.int64, sep=",")[1:-1]
            pos = st + starts

            for r in range(r0, r1):
                o0, o1 = offsets[r], offsets[r + 1]
                if o0 == o1:
                    continue    # window with no A/T sites (warned in main)
                wsl = sites[o0:o1]
                a = np.searchsorted(pos, ws[r], side="left")
                c = np.searchsorted(pos, we[r], side="left")
                p = pos[a:c]
                n_mod = 0
                if p.size:
                    idx = np.searchsorted(wsl, p)
                    on = (idx < wsl.size) & (wsl[np.minimum(idx, wsl.size - 1)] == p)
                    n_mod = int(on.sum())
                    off_site += int((~on).sum())   # in window, not a ref A/T
                rows.append((r, rid, n_mod))
    return rows, off_site, n_kept


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tss-windows", default=f"{TAB_DIR}/entropy_windows_tss.tsv")
    ap.add_argument("--fire-windows", default=f"{TAB_DIR}/entropy_windows_fire.tsv")
    ap.add_argument("--ft-root", default=FT_ROOT)
    ap.add_argument("--ref", default=REF_FA)
    ap.add_argument("--out-dir", default=TAB_DIR)
    ap.add_argument("--timepoints", nargs="+", default=SAMPLES)
    ap.add_argument("--chrom", default=None, help="restrict to one chromosome (testing)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tss = pd.read_csv(args.tss_windows, sep="\t")
    tss = tss.rename(columns={"gene_id": "window_id"})
    tss = tss[["window_id", "chrom", "win_start", "win_end"]]
    tss["source"] = "tss"
    fire = pd.read_csv(args.fire_windows, sep="\t")
    fire = fire.rename(columns={"region_id": "window_id"})
    fire = fire[["window_id", "chrom", "win_start", "win_end"]]
    fire["source"] = "fire"
    win = pd.concat([tss, fire], ignore_index=True)

    chroms = [args.chrom] if args.chrom else CHROMS
    win = win[win["chrom"].isin(chroms)].reset_index(drop=True)
    assert (win["win_end"] - win["win_start"] == WIN).all(), \
        f"every window must be exactly {WIN} bp"

    fa = pysam.FastaFile(args.ref)
    chrom_len = {c: fa.get_reference_length(c) for c in chroms}
    edge = (win["win_start"] < 0) | (win["win_end"] > win["chrom"].map(chrom_len))
    if edge.any():
        print(f"dropping {edge.sum()} windows leaving the chromosome", flush=True)
        win = win[~edge].reset_index(drop=True)
    print(f"{len(win)} windows ({(win['source'] == 'tss').sum()} tss, "
          f"{(win['source'] == 'fire').sum()} fire) on "
          f"{win['chrom'].nunique()} chromosomes; {WIN} bp, spanning reads only",
          flush=True)

    for s in args.timepoints:
        for chrom in chroms:
            if not m6a_path(args.ft_root, s, chrom).is_file():
                sys.exit(f"ERROR: missing input {m6a_path(args.ft_root, s, chrom)}")

    frames = []
    stats = {s: {"off_site": 0, "n_kept": 0} for s in args.timepoints}
    for chrom in chroms:
        sub = win[win["chrom"] == chrom].sort_values("win_start").reset_index(drop=True)
        if sub.empty:
            continue
        chrom_seq = fa.fetch(chrom).upper()
        b = build_sites(sub, chrom_seq)
        del chrom_seq
        if (b["n_sites"] == 0).any():
            print(f"WARNING: {(b['n_sites'] == 0).sum()} windows on {chrom} "
                  f"have no A/T sites; they get no rows", flush=True)
        merged = merge_intervals(b["ws"], b["we"])

        for s in args.timepoints:
            with pysam.TabixFile(str(m6a_path(args.ft_root, s, chrom))) as tb:
                best = survey_reads(tb, chrom, merged)
                rows, off, kept = fractions_for_chrom(tb, chrom, merged, best, b)
            stats[s]["off_site"] += off
            stats[s]["n_kept"] += kept
            print(f"{chrom} {s}: kept {kept} reads, {len(rows)} spanning "
                  f"(read, window) pairs; off-site calls {off}", flush=True)
            if not rows:
                continue
            df = pd.DataFrame(rows, columns=["widx", "RID", "n_mod"])
            df["window_id"] = sub["window_id"].to_numpy()[df["widx"]]
            df["source"] = sub["source"].to_numpy()[df["widx"]]
            df["n_sites"] = b["n_sites"][df["widx"]]
            df["timepoint"] = s
            frames.append(df.drop(columns="widx"))

    cols = ["window_id", "timepoint", "RID", "n_sites", "n_mod", "frac"]
    out = (pd.concat(frames, ignore_index=True) if frames
           else pd.DataFrame(columns=["RID", "n_mod", "window_id", "source",
                                      "n_sites", "timepoint"]))
    out["frac"] = out["n_mod"] / out["n_sites"]
    assert (out["n_mod"] <= out["n_sites"]).all()
    assert out["frac"].between(0, 1).all()

    for src, fname in (("tss", "read_fractions_tss.tsv.gz"),
                       ("fire", "read_fractions_fire.tsv.gz")):
        sl = out.loc[out["source"] == src, cols]
        path = out_dir / fname
        sl.to_csv(path, sep="\t", index=False, float_format="%.6g")
        print(f"wrote {path} ({len(sl)} rows, {sl['window_id'].nunique()} windows)",
              flush=True)

    for s in args.timepoints:
        print(f"{s}: kept reads {stats[s]['n_kept']}; off-site calls "
              f"{stats[s]['off_site']}", flush=True)


if __name__ == "__main__":
    main()
