#!/usr/bin/env python3
"""Per-region FIRE frequency across the LPS timecourse.

For every region of the shared FIRE peak universe (union of the four timepoints'
peaks, from 01_make_cre_universe.sh) and every timepoint:

    n_reads   distinct reads whose aligned span overlaps the region (>= 1 bp)
    n_fire    those reads carrying >= 1 FIRE element overlapping the region (>= 1 bp)
    freq      n_fire / n_reads   (NA when n_reads == 0)

Among moelcules that were observed at this shared FIRE region, what fraction contained
a FIRE element somewhere in that region?

Inputs (from co-accessibility 01/02):
    <spans-root>/universe/fire_peaks_union.bed      chrom start end   (BED3, merged)
    <spans-root>/<s>/<s>.read_spans.bed.gz          aligned read spans (tabix)
    <FIRE>/<s>/additional-outputs-v0.1/fire-peaks/<s>-v0.1-fire-elements.bed.gz

Outputs in --out-dir:
    fire_frequency_long.tsv.gz   one row per (region, timepoint)
    fire_frequency_wide.tsv.gz   one row per region, per-timepoint columns
"""

import argparse
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

SAMPLES = ["LPS_0", "LPS_5", "LPS_10", "LPS_15"]
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

BEDTOOLS = "/project/spott/cshan/envs/bedtools/bin/bedtools"
TABIX = "/project/spott/cshan/envs/dimelo/bin/tabix"
# Per-read FIRE element calls (FDR <= 0.05) from the FIRE pipeline, same numerator
# source as 03_coaccess_cres.py; see the note there on why NOT FiberHMM/ft_result_dir.
FIRE_ROOT = "/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE"

N_UNION_REGIONS = 143070  # rows in fire_peaks_union.bed

###############################
# function for running bedtools intersect

# BEDTOOLS, "intersect",
#     "-a", "regions.bed",
#     "-b", "reads.bed",
#     "-wa", "-wb"


###############################
def run_intersect(args_list, desc):
    p = subprocess.run(args_list, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"ERROR: {desc} failed:\n{p.stderr[:2000]}")
    return p.stdout


###############################
# Count n_reads and n_fire for one chromosome, returning a dictionary of region_id -> (n_reads, n_fire)
###############################
# for one chr at one timepoint, it asks for every shared FIRE region
    # 1. how many unique reads overlap this region (>= 1 bp)
    # 2. how many of those reads have a FIRE element overlapping this region (>= 1 bp)

def counts_for_chrom(chrom, regions, spans_path, fe_path, tmpdir):
    """{region_id: (n_reads, n_fire)} for one chromosome
    """
    # subet to a single chr
    sub = regions[regions["chrom"] == chrom]
    if sub.empty:
        return {}, 0
    reg_chr = tmpdir / f"regions.{chrom}.bed"
    # write the regions bed for this chromosome, 4 columns (chrom, start, end, region_id)
    sub[["chrom", "start", "end", "region_id"]].to_csv(
        reg_chr, sep="\t", header=False, index=False)

    # extract the read spans and read levelfire elements for this chr
    spans_chr = tmpdir / f"spans.{chrom}.bed"
    with open(spans_chr, "w") as fh:
        # tabix fire-elements.bed.gz chr1
        p = subprocess.run([TABIX, str(spans_path), chrom], stdout=fh, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        sys.exit(f"ERROR: tabix on spans for {chrom} failed:\n{p.stderr[:2000]}")

    fe_chr = tmpdir / f"fe.{chrom}.bed"
    with open(fe_chr, "w") as fh:
        p = subprocess.run([TABIX, str(fe_path), chrom], stdout=fh, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        sys.exit(f"ERROR: tabix on fire elements for {chrom} failed:\n{p.stderr[:2000]}")

    # Reads covering each region: any overlap >= 1 bp
    cov = defaultdict(set)
    # run bedtools intersect -a regions.bed -b spans.bed -wa -wb
    for line in run_intersect(
            [BEDTOOLS, "intersect", "-a", str(reg_chr), "-b", str(spans_chr), "-wa", "-wb"],
            f"read intersect {chrom}").splitlines():
        f = line.split("\t")
        cov[f[3]].add(f[7])          # region_id, read name

    # FIRE elements covering each region: any overlap >= 1 bp
    # count which FIRE positive reads overlap each region
    acc = defaultdict(set)
    for line in run_intersect(
            [BEDTOOLS, "intersect", "-a", str(reg_chr), "-b", str(fe_chr), "-wa", "-wb"],
            f"fire element intersect {chrom}").splitlines():
        f = line.split("\t")
        acc[f[3]].add(f[7])          # region_id, read name

    # restrict FIRE reads to reads that exist in the read spans (cov)
    n_orphan = 0 # count reads that are FIRE positive but not in the read spans (cov)
    counts = {}

    # loop over each region that has least one covering read
    for r, reads in cov.items():
        # get the set of FIRE positive reads for this region, if any
        fire = acc.get(r, set())
        before = len(fire)
        fire &= reads
        n_orphan += before - len(fire)
        counts[r] = (len(reads), len(fire))
    for r in acc:
        if r not in cov:
            n_orphan += len(acc[r])

    for f in (reg_chr, spans_chr, fe_chr):
        f.unlink(missing_ok=True)
    return counts, n_orphan


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--union-bed",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/co-accessibility/universe/fire_peaks_union.bed")
    ap.add_argument("--spans-root",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/co-accessibility",
                    help="read spans at <spans-root>/<s>/<s>.read_spans.bed.gz")
    ap.add_argument("--out-dir",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency")
    ap.add_argument("--timepoints", nargs="+", default=SAMPLES)
    ap.add_argument("--chrom", default=None, help="restrict to one chromosome (testing)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmpdir = out_dir / ".tmp"
    tmpdir.mkdir(exist_ok=True)

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
    print(f"union regions: {len(regions)} on {regions['chrom'].nunique()} chromosomes",
          flush=True)

    long_parts = []
    for s in args.timepoints:
        spans = Path(args.spans_root) / s / f"{s}.read_spans.bed.gz"
        fe = Path(FIRE_ROOT) / s / "additional-outputs-v0.1" / "fire-peaks" / \
            f"{s}-v0.1-fire-elements.bed.gz"
        for f in (spans, fe):
            if not f.is_file():
                sys.exit(f"ERROR: missing input {f}")

        n_reads = np.zeros(len(regions), dtype=np.int64)
        n_fire = np.zeros(len(regions), dtype=np.int64)
        row_of = {rid: i for i, rid in enumerate(regions["region_id"])}
        orphans = 0
        for chrom in chroms:
            counts, n_orph = counts_for_chrom(chrom, regions, spans, fe, tmpdir)
            orphans += n_orph
            for rid, (nr, nf) in counts.items():
                i = row_of[rid]
                n_reads[i] = nr
                n_fire[i] = nf

        assert (n_fire <= n_reads).all()
        part = regions.copy()
        part["timepoint"] = s
        part["n_reads"] = n_reads
        part["n_fire"] = n_fire
        long_parts.append(part)
        q = np.percentile(n_reads, [25, 50, 75])
        print(f"{s}: orphan element hits {orphans}; "
              f"n_reads quartiles {q[0]:.0f}/{q[1]:.0f}/{q[2]:.0f}; "
              f"zero-coverage regions {(n_reads == 0).sum()}", flush=True)

    long_df = pd.concat(long_parts, ignore_index=True)
    # NaN where n_reads == 0 (Series.where masks the zero denominators)
    long_df["freq"] = long_df["n_fire"] / long_df["n_reads"].where(long_df["n_reads"] > 0)
    assert len(long_df) == len(regions) * len(args.timepoints)
    assert long_df["freq"].dropna().between(0, 1).all()
    assert (long_df["freq"].isna() == (long_df["n_reads"] == 0)).all()

    long_out = out_dir / "fire_frequency_long.tsv.gz"
    long_df.to_csv(long_out, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {long_out} ({len(long_df)} rows)", flush=True)

    # wide pivot, keeping the genomic order of the union bed
    wide = regions.copy()
    for s, part in zip(args.timepoints, long_parts):
        nr = part["n_reads"].to_numpy()
        nf = part["n_fire"].to_numpy()
        wide[f"{s}_n_reads"] = nr
        wide[f"{s}_n_fire"] = nf
        wide[f"{s}_freq"] = np.where(nr > 0, nf / np.maximum(nr, 1), np.nan)
    wide_out = out_dir / "fire_frequency_wide.tsv.gz"
    wide.to_csv(wide_out, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {wide_out} ({len(wide)} rows)", flush=True)


if __name__ == "__main__":
    main()
