#!/usr/bin/env python3
"""Average m6A level in 10 bp bins across the shared FIRE peak universe.

For every region of the union peak universe (same fire_peaks_union.bed as
01_fire_frequency.py), tiled into 10 bp bins, and every timepoint:

    n_sites   reference A/T positions in the bin (m6A-callable sites)
    n_cov     read x site events: sum over the bin's A/T sites of read depth
              (aligned-span depth) at that site
    n_mod     called m6A read x site events (calls landing on a reference A/T)
    m6a       n_mod / n_cov   (NA when n_cov == 0)

Region-level values are the same counts summed over the region's bins, plus
n_reads = distinct reads whose aligned span overlaps the region (the
coverage-filter quantity for 05, analogous to n_reads in 01_fire_frequency.py).
Per (read, region) overlap, the read's own m6A fraction over the region
(n_mod / covered A/T sites) is also emitted -- the test unit for the per-read
Wilcoxon in 05 (reads, not read x site events, are the independent unit).

Reads come from the per-chromosome ft extract BED12 (FiberHMM-recalled BAMs),
streamed twice per (chromosome, timepoint); nothing intermediate is written.
Parsing follows read_ft_mod_region (region_data_utils.R): the first and last
blocks are ft extract sentinels (0 bp at read start, 1 bp at read end) and are
dropped; split reads are deduped keeping the longest alignment per read name.
Calls on a non-A/T reference base (read-vs-reference mismatches, ~1%) are
counted as off-site and excluded from n_mod.

Outputs in --out-dir:
    m6a_bins_long.tsv.gz     one row per (bin, timepoint)
    m6a_region_long.tsv.gz   one row per (region, timepoint)
    m6a_reads_long.tsv.gz    one row per (read, region, timepoint) overlap
"""

import argparse
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

warnings.filterwarnings("ignore", category=DeprecationWarning)

SAMPLES = ["LPS_0", "LPS_5", "LPS_10", "LPS_15"]
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

# Per-read m6A calls from the FiberHMM-recalled BAMs (NOT the FIRE CRAM: the
# extracts and their read spans are one universe, so numerator and denominator
# stay consistent).
FT_ROOT = "/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/ft_result_dir"
REF_FA = "/project/spott/reference/human/GRCh38/hg38.fa"

N_UNION_REGIONS = 143070  # rows in fire_peaks_union.bed
BIN_SIZE = 10


def open_cut(path, fields):
    """Stream `zcat path | cut -f<fields>`; returns (procs, byte-line iterator)."""
    p1 = subprocess.Popen(["zcat", str(path)], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["cut", "-f", fields], stdin=p1.stdout,
                          stdout=subprocess.PIPE)
    p1.stdout.close()
    return (p1, p2), p2.stdout


def close_cut(procs, path):
    for p in procs:
        if p.wait() != 0:
            sys.exit(f"ERROR: streaming {path} failed (exit {p.returncode})")


###############################
# Per-chromosome bin / A/T-site arrays from the union regions
###############################
def build_bins_and_sites(sub, chrom_seq):
    """Tile each region into BIN_SIZE bins (last bin truncated at the region
    end) and locate every reference A/T position inside them."""
    bin_start_l, bin_end_l, bin_region_l, bin_idx_l = [], [], [], []
    sites_l, site_bin_l = [], []
    first_bin = 0
    for ri, (s, e) in enumerate(zip(sub["start"], sub["end"])):
        bs = np.arange(s, e, BIN_SIZE, dtype=np.int64)
        k = len(bs)
        bin_start_l.append(bs)
        bin_end_l.append(np.minimum(bs + BIN_SIZE, e))
        bin_region_l.append(np.full(k, ri, dtype=np.int64))
        bin_idx_l.append(np.arange(k, dtype=np.int64))
        seg = np.frombuffer(chrom_seq[s:e].encode(), dtype=np.uint8)
        at = s + np.nonzero((seg == ord("A")) | (seg == ord("T")))[0]
        sites_l.append(at)
        site_bin_l.append(first_bin + (at - s) // BIN_SIZE)
        first_bin += k
    b = {
        "bin_start": np.concatenate(bin_start_l),
        "bin_end": np.concatenate(bin_end_l),
        "bin_region": np.concatenate(bin_region_l),
        "bin_idx": np.concatenate(bin_idx_l),
        "sites": np.concatenate(sites_l),      # sorted: regions are sorted + merged
        "site_bin": np.concatenate(site_bin_l),
        "n_bins": first_bin,
    }
    b["n_sites"] = np.bincount(b["site_bin"], minlength=b["n_bins"]).astype(np.int64)
    assert (b["bin_end"] - b["bin_start"]).sum() == (sub["end"] - sub["start"]).sum()
    return b


###############################
# Pass 1: longest alignment per read name (dedup rule of read_ft_mod_region)
###############################
def survey_reads(path):
    """{read name: (span_length, start)} of the longest alignment per read."""
    best = {}
    procs, stream = open_cut(path, "2,3,4")
    for line in stream:
        f = line.split(b"\t")
        st = int(f[0])
        ln = int(f[1]) - st
        rid = f[2].rstrip()
        cur = best.get(rid)
        if cur is None or ln > cur[0]:
            best[rid] = (ln, st)
    close_cut(procs, path)
    return best


###############################
# Pass 2: per-bin n_cov / n_mod and per-(read, region) counts for one chromosome
###############################
def counts_for_chrom(path, best, b, region_starts, region_ends):
    """Coverage uses a difference array over the sorted A/T-site index (a read
    covers a contiguous run of sites), so it is O(1) per read; only reads
    overlapping >= 1 region pay for parsing their blockStarts."""
    cov_diff = np.zeros(b["sites"].size + 1, dtype=np.int64)
    n_mod = np.zeros(b["n_bins"], dtype=np.int64)
    rr_region, rr_sites, rr_mod = [], [], []  # one entry per (read, region) overlap
    off_site = 0
    n_kept = 0

    sites = b["sites"]
    procs, stream = open_cut(path, "2-4,12")
    for line in stream:
        f = line.split(b"\t")
        st = int(f[0])
        en = int(f[1])
        rid = f[2]
        cur = best.get(rid)
        if cur is None or cur != (en - st, st):
            continue                      # shorter split alignment, or already counted
        best[rid] = None
        n_kept += 1

        # regions overlapped by the aligned span (any overlap >= 1 bp)
        r0 = np.searchsorted(region_ends, st, side="right")
        r1 = np.searchsorted(region_starts, en, side="left")
        if r1 <= r0:
            continue

        # span-depth over the in-bin A/T sites: contiguous site-index run
        i0 = np.searchsorted(sites, st, side="left")
        i1 = np.searchsorted(sites, en, side="left")
        cov_diff[i0] += 1
        cov_diff[i1] -= 1

        # m6A calls: drop the sentinel first/last blocks; 1 bp blocks, so the
        # modified base (0-based) is span start + blockStart
        starts = np.fromstring(f[3].rstrip().decode(), sep=",", dtype=np.int64)[1:-1]
        pos = st + starts
        j = np.searchsorted(region_starts, pos, side="right") - 1
        in_reg = (j >= 0) & (pos < region_ends[np.maximum(j, 0)])
        pin = pos[in_reg]
        idx = np.searchsorted(sites, pin)
        on = sites[np.minimum(idx, sites.size - 1)] == pin
        np.add.at(n_mod, b["site_bin"][idx[on]], 1)
        off_site += int((~on).sum())      # in a region but not on a reference A/T
        jon = j[in_reg][on]               # region index of each counted call

        # per-read counts for each overlapped region (usually 1)
        for r in range(r0, r1):
            a = max(st, region_starts[r])
            c = min(en, region_ends[r])
            k0 = np.searchsorted(sites, a, side="left")
            k1 = np.searchsorted(sites, c, side="left")
            rr_region.append(r)
            rr_sites.append(k1 - k0)
            rr_mod.append(int((jon == r).sum()))
    close_cut(procs, path)

    site_depth = np.cumsum(cov_diff[:-1])
    n_cov = np.bincount(b["site_bin"], weights=site_depth,
                        minlength=b["n_bins"]).astype(np.int64)
    assert (n_mod <= n_cov).all()
    rr = (np.array(rr_region, dtype=np.int64),
          np.array(rr_sites, dtype=np.int64),
          np.array(rr_mod, dtype=np.int64))
    return n_cov, n_mod, rr, off_site, n_kept


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--union-bed",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/co-accessibility/universe/fire_peaks_union.bed")
    ap.add_argument("--ft-root", default=FT_ROOT,
                    help="m6A extracts at <ft-root>/<s>/extracted_results/m6a_by_chr/")
    ap.add_argument("--ref", default=REF_FA)
    ap.add_argument("--out-dir",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/methylation")
    ap.add_argument("--timepoints", nargs="+", default=SAMPLES)
    ap.add_argument("--chrom", default=None, help="restrict to one chromosome (testing)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    regions = pd.read_csv(args.union_bed, sep="\t", header=None,
                          names=["chrom", "start", "end"])
    if args.chrom is None and len(regions) != N_UNION_REGIONS:
        sys.exit(f"ERROR: {args.union_bed} has {len(regions)} rows, "
                 f"expected {N_UNION_REGIONS} (rerun 01_make_cre_universe.sh?)")
    chroms = [args.chrom] if args.chrom else CHROMS
    regions = regions[regions["chrom"].isin(chroms)].reset_index(drop=True)
    regions["region_id"] = (regions["chrom"] + ":" + regions["start"].astype(str)
                            + "-" + regions["end"].astype(str))
    assert regions["region_id"].is_unique
    print(f"union regions: {len(regions)} on {regions['chrom'].nunique()} chromosomes; "
          f"bin size {BIN_SIZE} bp", flush=True)

    for s in args.timepoints:
        for chrom in chroms:
            f = Path(args.ft_root) / s / "extracted_results" / "m6a_by_chr" / \
                f"{s}.ft_extracted_m6a.{chrom}.bed.gz"
            if not f.is_file():
                sys.exit(f"ERROR: missing input {f}")

    fa = Fasta(args.ref)
    bins_by_tp = {s: [] for s in args.timepoints}
    regs_by_tp = {s: [] for s in args.timepoints}
    reads_by_tp = {s: [] for s in args.timepoints}
    stats = {s: {"off_site": 0, "n_kept": 0} for s in args.timepoints}
    n_bins_total = 0

    for chrom in chroms:
        sub = regions[regions["chrom"] == chrom].reset_index(drop=True)
        if sub.empty:
            continue
        chrom_seq = str(fa[chrom][:]).upper()
        b = build_bins_and_sites(sub, chrom_seq)
        del chrom_seq
        n_bins_total += b["n_bins"]
        region_starts = sub["start"].to_numpy()
        region_ends = sub["end"].to_numpy()

        bins_base = pd.DataFrame({
            "chrom": chrom,
            "bin_start": b["bin_start"],
            "bin_end": b["bin_end"],
            "region_id": sub["region_id"].to_numpy()[b["bin_region"]],
            "bin_idx": b["bin_idx"],
            "n_sites": b["n_sites"],
        })
        reg_n_sites = np.bincount(b["bin_region"], weights=b["n_sites"],
                                  minlength=len(sub)).astype(np.int64)

        for s in args.timepoints:
            path = Path(args.ft_root) / s / "extracted_results" / "m6a_by_chr" / \
                f"{s}.ft_extracted_m6a.{chrom}.bed.gz"
            best = survey_reads(path)
            n_cov, n_mod, rr, off, kept = counts_for_chrom(
                path, best, b, region_starts, region_ends)
            stats[s]["off_site"] += off
            stats[s]["n_kept"] += kept
            print(f"{chrom} {s}: alignments {len(best)}, kept {kept}; "
                  f"off-site calls {off}", flush=True)

            part = bins_base.copy()
            part["timepoint"] = s
            part["n_cov"] = n_cov
            part["n_mod"] = n_mod
            bins_by_tp[s].append(part)

            rpart = sub.copy()
            rpart["timepoint"] = s
            rpart["n_sites"] = reg_n_sites
            rpart["n_cov"] = np.bincount(b["bin_region"], weights=n_cov,
                                         minlength=len(sub)).astype(np.int64)
            rpart["n_mod"] = np.bincount(b["bin_region"], weights=n_mod,
                                         minlength=len(sub)).astype(np.int64)
            rpart["n_reads"] = np.bincount(rr[0], minlength=len(sub)).astype(np.int64)
            # per-read sums must reproduce the pooled bin-level counts
            assert np.array_equal(
                np.bincount(rr[0], weights=rr[1], minlength=len(sub)).astype(np.int64),
                rpart["n_cov"].to_numpy())
            assert np.array_equal(
                np.bincount(rr[0], weights=rr[2], minlength=len(sub)).astype(np.int64),
                rpart["n_mod"].to_numpy())
            regs_by_tp[s].append(rpart)

            reads_by_tp[s].append(pd.DataFrame({
                "chrom": chrom,
                "region_id": sub["region_id"].to_numpy()[rr[0]],
                "timepoint": s,
                "n_sites_cov": rr[1],
                "n_mod": rr[2],
            }))

    # long tables, timepoint-major with genomic order within (as in 01)
    bins_long = pd.concat([pd.concat(bins_by_tp[s], ignore_index=True)
                           for s in args.timepoints], ignore_index=True)
    bins_long = bins_long[["chrom", "bin_start", "bin_end", "region_id", "bin_idx",
                           "timepoint", "n_sites", "n_cov", "n_mod"]]
    bins_long["m6a"] = bins_long["n_mod"] / bins_long["n_cov"].where(bins_long["n_cov"] > 0)
    assert len(bins_long) == n_bins_total * len(args.timepoints)
    assert bins_long["m6a"].dropna().between(0, 1).all()
    assert (bins_long["m6a"].isna() == (bins_long["n_cov"] == 0)).all()

    reg_long = pd.concat([pd.concat(regs_by_tp[s], ignore_index=True)
                          for s in args.timepoints], ignore_index=True)
    reg_long["m6a"] = reg_long["n_mod"] / reg_long["n_cov"].where(reg_long["n_cov"] > 0)
    assert len(reg_long) == len(regions) * len(args.timepoints)
    assert (reg_long["n_mod"] <= reg_long["n_cov"]).all()

    for s in args.timepoints:
        sl = reg_long[reg_long["timepoint"] == s]
        q = np.percentile(sl["n_reads"], [25, 50, 75])
        print(f"{s}: kept reads {stats[s]['n_kept']}; off-site calls {stats[s]['off_site']}; "
              f"n_reads quartiles {q[0]:.0f}/{q[1]:.0f}/{q[2]:.0f}; "
              f"median region m6a {sl['m6a'].median():.4f}; "
              f"zero-coverage regions {(sl['n_cov'] == 0).sum()}", flush=True)

    bins_out = out_dir / "m6a_bins_long.tsv.gz"
    bins_long.to_csv(bins_out, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {bins_out} ({len(bins_long)} rows)", flush=True)

    reg_out = out_dir / "m6a_region_long.tsv.gz"
    reg_long.to_csv(reg_out, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {reg_out} ({len(reg_long)} rows)", flush=True)

    reads_long = pd.concat([pd.concat(reads_by_tp[s], ignore_index=True)
                            for s in args.timepoints], ignore_index=True)
    reads_long["frac"] = reads_long["n_mod"] / \
        reads_long["n_sites_cov"].where(reads_long["n_sites_cov"] > 0)
    assert len(reads_long) == reg_long["n_reads"].sum()
    assert (reads_long["n_mod"] <= reads_long["n_sites_cov"]).all()
    assert (reads_long["frac"].isna() == (reads_long["n_sites_cov"] == 0)).all()

    reads_out = out_dir / "m6a_reads_long.tsv.gz"
    reads_long.to_csv(reads_out, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {reads_out} ({len(reads_long)} rows)", flush=True)


if __name__ == "__main__":
    main()
