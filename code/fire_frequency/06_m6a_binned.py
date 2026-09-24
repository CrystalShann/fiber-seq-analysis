#!/usr/bin/env python3
"""Average m6A level in 10 bp bins across the shared FIRE peak universe.

For every region of the union peak universe (same fire_peaks_union.bed as
03_fire_frequency.py), binned into 10 bp bins, and every timepoint:

    n_sites   reference A/T positions in the bin 
    n_cov     read x site events: sum over the bin's A/T sites of read depth
              (aligned-span depth) at that site
    n_mod     called m6A read x site events (calls landing on a reference A/T)
    m6a       n_mod / n_cov   (NA when n_cov == 0)


n_reads = distinct reads whose aligned span overlaps the region 


Outputs in --out-dir:
    m6a_bins_long.tsv.gz     one row per (bin, timepoint)
    m6a_region_long.tsv.gz   one row per (region, timepoint)
    m6a_reads_long.tsv.gz    one row per (read, region, timepoint) overlap,
                             read_id = original Fiber-seq read name
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

# Per-read m6A calls from the FiberHMM-recalled BAMs 
FT_ROOT = "/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/ft_result_dir"
REF_FA = "/project/spott/reference/human/GRCh38/hg38.fa"

BIN_SIZE = 10

# read BED file
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
# Per-chromosome bin and A/T-site from the union regions
###############################
def build_bins_and_sites(sub, chrom_seq):
  # For every bin, get the start and end coordinate, region window, bin number within region
  # For A/T site, get the genomic position and which bin it belongs to
    """Tile each region into BIN_SIZE bins and locate every reference A/T position inside them."""
    bin_start_l, bin_end_l, bin_region_l, bin_idx_l = [], [], [], []
    sites_l, site_bin_l = [], []
    first_bin = 0
    # loop over every union FIRE region on given chr
    for ri, (s, e) in enumerate(zip(sub["start"], sub["end"])):
        # create starts every 10 bp
        bs = np.arange(s, e, BIN_SIZE, dtype=np.int64)
        # number of bins
        k = len(bs)
        bin_start_l.append(bs)
        # bin end = bin start + 10
        # truncate the final bin at FIRE region boundary
        bin_end_l.append(np.minimum(bs + BIN_SIZE, e))
        # region which FIRE region each bin belongs to (one region have multiple bins)
        bin_region_l.append(np.full(k, ri, dtype=np.int64))
        # record the bin's index within the FIRE region
        bin_idx_l.append(np.arange(k, dtype=np.int64))
        
        # retrieve the ref seq for the FIRE region
        seg = np.frombuffer(chrom_seq[s:e].encode(), dtype=np.uint8)
        # find all positions where ref seq is A or T
        at = s + np.nonzero((seg == ord("A")) | (seg == ord("T")))[0]
        sites_l.append(at)
        # assign A/T position to a specific bin
          # ex: A/T at 1025 (1025-1000)/10 = 2 -> belong to bin 2
        site_bin_l.append(first_bin + (at - s) // BIN_SIZE)
        
        # process next FIRE region
        first_bin += k
        
    b = {
        # corrdinates of every 10 bp bin
        "bin_start": np.concatenate(bin_start_l),
        "bin_end": np.concatenate(bin_end_l),
        # which FIRE region each bin belongs to
        "bin_region": np.concatenate(bin_region_l),
        # bin number within a FIRE region
        "bin_idx": np.concatenate(bin_idx_l),
        # Every callable A/T position in all FIRE regions on this chromosome
        "sites": np.concatenate(sites_l),   
        # For each A/T position, which 10-bp bin contains it.
        "site_bin": np.concatenate(site_bin_l),
        # total number of bins on this chromosome
        "n_bins": first_bin,
    }
    # how many A/T sites belong to each bin
    b["n_sites"] = np.bincount(b["site_bin"], minlength=b["n_bins"]).astype(np.int64)
    assert (b["bin_end"] - b["bin_start"]).sum() == (sub["end"] - sub["start"]).sum()
    return b


###############################
# keep longest alignment per read from ft extract 
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
# per-bin n_cov / n_mod and per-(read, region) counts for one chromosome
###############################


# path: the chromosome-specific ft extract BED12 file
# best: output from survey_reads(), telling you which alignment to keep for each read
# b: information about the 10-bp FIRE bins and all reference A/T positions
# region_starts: start coordinates of FIRE union regions
# region_ends: end coordinates of FIRE union regions


def counts_for_chrom(path, best, b, region_starts, region_ends):
    """counts how many reference A/T sites are covered by reads (n_cov), and how many valid m6a calls fall on A/T sites (n_mod)."""
    
    # difference array used to calculate coverage at every reference A/T site
    cov_diff = np.zeros(b["sites"].size + 1, dtype=np.int64)
    
    # counter for every 10 bp bins to store number of m6a calls in each bin
    n_mod = np.zeros(b["n_bins"], dtype=np.int64)
    
    
    # rr_region = which FIRE region
    # rr_sites  = how many reference A/T sites this read covers in that region
    # rr_mod    = how many m6A calls this read has in that region
    
    # rr_read    = original Fiber-seq read name
    rr_region, rr_sites, rr_mod, rr_read = [], [], [], []
    
    # m6A calls that occur inside a FIRE region but do not correspond to a reference A/T base -> exclued from n_mod
    off_site = 0
    # how many unique selected read alignments are kept.
    n_kept = 0


    # array of every reference A/T position inside the FIRE regions on this chromosome
    sites = b["sites"]
    
    # read ft extracted m6a BED file for alignment start, end, read name and blockStarts
    procs, stream = open_cut(path, "2-4,12")
    for line in stream:
        f = line.split(b"\t")
        st = int(f[0])
        en = int(f[1])
        rid = f[2]
        cur = best.get(rid)
        if cur is None or cur != (en - st, st):
            continue                    
        best[rid] = None
        n_kept += 1

        # find the first FIRE region whose end occurs after the read starts -> all FIRE regions overlapping the read by at least 1 bo
        
        # FIRE regions upstram of the read
        r0 = np.searchsorted(region_ends, st, side="right")
        # FIRE regions downstream of the read
        r1 = np.searchsorted(region_starts, en, side="left")
        # stop processing the read if there is no overlapping FIRE regions
        if r1 <= r0:
            continue

        # Count read coverage of all A/T site
        
        # first callable A/T site at or after the read start
        i0 = np.searchsorted(sites, st, side="left")
        i1 = np.searchsorted(sites, en, side="left")
        
        # coverage increases by 1
        cov_diff[i0] += 1
        cov_diff[i1] -= 1

        # extract m6a positions
        starts = np.fromstring(f[3].rstrip().decode(), sep=",", dtype=np.int64)[1:-1]
        # convert each m6a position from relative start pos to genomic coordinate
        pos = st + starts
        
        # for every m6a coordinate, find the most recent FIRE region start preceding the coordinate
        
          # find the corresponding A/T site
          # find which 10-bp bin that site belongs to
          # add one m6A call to that bin
          
        # Verify that each m6A position is actually inside that region
        j = np.searchsorted(region_starts, pos, side="right") - 1
        in_reg = (j >= 0) & (pos < region_ends[np.maximum(j, 0)])
        pin = pos[in_reg]
        # For each m6A coordinate, find its potential location in the sorted reference A/T
        idx = np.searchsorted(sites, pin)
        on = sites[np.minimum(idx, sites.size - 1)] == pin
        np.add.at(n_mod, b["site_bin"][idx[on]], 1)
        # how many m6A calls were inside FIRE regions but not located on a reference A/T
        off_site += int((~on).sum())   
        jon = j[in_reg][on]             

        # loop through every FIRE region span by the read by at least 1np
        for r in range(r0, r1):
            # beginning of the actual read-region overlap
            a = max(st, region_starts[r])
            # end of the actual overlap
            c = min(en, region_ends[r])
            # first callable A/T position inside the overlap
            k0 = np.searchsorted(sites, a, side="left")
            # position after the last callable A/T site inside the overlap
            k1 = np.searchsorted(sites, c, side="left")
            # record FIRE region index 
            rr_region.append(r)
            # how many reference A/T positions this read covers inside this FIRE region
            rr_sites.append(k1 - k0)
            # how many valid m6A calls from this read belong to this region
            rr_mod.append(int((jon == r).sum()))
            rr_read.append(rid.decode())
    close_cut(procs, path)

    site_depth = np.cumsum(cov_diff[:-1])
    
    # Sum A/T-site coverage within each 10-bp bin
    n_cov = np.bincount(b["site_bin"], weights=site_depth,
                        minlength=b["n_bins"]).astype(np.int64)
    assert (n_mod <= n_cov).all()
    rr = (np.array(rr_region, dtype=np.int64),
          np.array(rr_sites, dtype=np.int64),
          np.array(rr_mod, dtype=np.int64),
          np.array(rr_read, dtype=object))
    return n_cov, n_mod, rr, off_site, n_kept


def main():
  # set up CLI
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--union-bed",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/universe/fire_peaks_union.bed")
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
    n_union_regions = len(regions)  # rows in fire_peaks_union.bed
    chroms = [args.chrom] if args.chrom else CHROMS
    regions = regions[regions["chrom"].isin(chroms)].reset_index(drop=True)
    regions["region_id"] = (regions["chrom"] + ":" + regions["start"].astype(str)
                            + "-" + regions["end"].astype(str))
    assert regions["region_id"].is_unique
    print(f"union regions: {len(regions)} of {n_union_regions} in {args.union_bed} "
          f"on {regions['chrom'].nunique()} chromosomes; bin size {BIN_SIZE} bp", flush=True)

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

            # Build the bin-level result table
            part = bins_base.copy()
            part["timepoint"] = s
            part["n_cov"] = n_cov
            part["n_mod"] = n_mod
            bins_by_tp[s].append(part)

            # Build whole-region results
            rpart = sub.copy()
            rpart["timepoint"] = s
            rpart["n_sites"] = reg_n_sites
            
            # Aggregate bin coverage into whole FIRE region
            rpart["n_cov"] = np.bincount(b["bin_region"], weights=n_cov,
                                         minlength=len(sub)).astype(np.int64)
                                         
            # Aggregate m6A calls into whole regions
            rpart["n_mod"] = np.bincount(b["bin_region"], weights=n_mod,
                                         minlength=len(sub)).astype(np.int64)
                                         
            # Calculate number of reads overlapping each FIRE region
            rpart["n_reads"] = np.bincount(rr[0], minlength=len(sub)).astype(np.int64)
          
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
                "read_id": rr[3],
            }))

    # Combine all chromosomes into the bin-level table
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
    reads_long = reads_long[["chrom", "region_id", "timepoint", "n_sites_cov",
                             "n_mod", "frac", "read_id"]]
    assert len(reads_long) == reg_long["n_reads"].sum()
    assert (reads_long["n_mod"] <= reads_long["n_sites_cov"]).all()
    assert (reads_long["frac"].isna() == (reads_long["n_sites_cov"] == 0)).all()

    reads_out = out_dir / "m6a_reads_long.tsv.gz"
    reads_long.to_csv(reads_out, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {reads_out} ({len(reads_long)} rows)", flush=True)


if __name__ == "__main__":
    main()
