#!/usr/bin/env python3
"""Per-region FIRE frequency across the LPS timecourse.

For every region of the shared FIRE peak universe (union of the four timepoints'
peaks, from 01_make_fire_universe.sh) and every timepoint:

    n_reads   distinct reads whose aligned span fully covers the region
    n_fire    fully covering reads carrying >= 1 FIRE element that covers >= 50%
              of the region
    freq      n_fire / n_reads   (NA when n_reads == 0)

Among moelcules that were observed at this shared FIRE region, what fraction contained
a FIRE element somewhere in that region?

Inputs (from 01_make_fire_universe.sh and 02_read_spans.sh):
    <spans-root>/universe/fire_peaks_union.bed      chrom start end   (BED3, merged)
    <spans-root>/<s>/<s>.read_spans.bed.gz          aligned read spans (tabix)
    <FIRE>/<s>/additional-outputs-v0.1/fire-peaks/<s>-v0.1-fire-elements.bed.gz

Outputs in --out-dir:
    fire_frequency_long.tsv.gz   one row per (region, timepoint)
    fire_frequency_wide.tsv.gz   one row per region, per-timepoint columns
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "parsing_functions"))
import fiber_regions as fr

SAMPLES = ["LPS_0", "LPS_5", "LPS_10", "LPS_15"]
# Per-read FIRE element calls (FDR <= 0.05) from the FIRE pipeline, same numerator
# source as ../co-accessibility/macrophage/03_coaccess_cres.py; see the note there on why NOT FiberHMM/ft_result_dir.
FIRE_ROOT = fr.MACROPHAGE_FIRE_ROOT


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--union-bed",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/universe/fire_peaks_union.bed")
    ap.add_argument("--spans-root",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency",
                    help="read spans at <spans-root>/<s>/<s>.read_spans.bed.gz")
    ap.add_argument("--out-dir",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency")
    ap.add_argument("--timepoints", nargs="+", default=SAMPLES)
    ap.add_argument("--read-rule", choices=fr.READ_RULES, default="span",
                    help="span: read span contains the region; aligned: aligned blocks "
                         "cover every base (needs <s>.aligned_blocks.bed.gz)")
    ap.add_argument("--fire-overlap-fraction", type=float, default=0.5,
                    help="minimum fraction of the region covered by ONE FIRE element")
    ap.add_argument("--chrom", default=None, help="restrict to one chromosome (testing)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmpdir = out_dir / ".tmp"
    tmpdir.mkdir(exist_ok=True)

    regions = pd.read_csv(args.union_bed, sep="\t", header=None,
                          names=["chrom", "start", "end"])
    n_union_regions = len(regions)  # rows in fire_peaks_union.bed
    chroms = [args.chrom] if args.chrom else fr.CHROMS
    regions = regions[regions["chrom"].isin(chroms)].reset_index(drop=True)
    regions["region_id"] = (regions["chrom"] + ":" + regions["start"].astype(str)
                            + "-" + regions["end"].astype(str))
    assert regions["region_id"].is_unique
    print(f"union regions: {len(regions)} of {n_union_regions} in {args.union_bed} "
          f"on {regions['chrom'].nunique()} chromosomes", flush=True)

    regions_bed = tmpdir / "regions.bed"
    regions[["chrom", "start", "end", "region_id"]].to_csv(
        regions_bed, sep="\t", header=False, index=False)

    long_parts = []
    for s in args.timepoints:
        spans = fr.reads_path(args.spans_root, s, args.read_rule)
        fe = fr.macrophage_elements_path(s, FIRE_ROOT)
        for f in (spans, fe):
            if not f.is_file():
                sys.exit(f"ERROR: missing input {f}")

        # n_reads = reads covering the whole region; n_fire = those with ONE FIRE
        # element covering >= --fire-overlap-fraction of the region
        n_reads = np.zeros(len(regions), dtype=np.int64)
        n_fire = np.zeros(len(regions), dtype=np.int64)
        row_of = {rid: i for i, rid in enumerate(regions["region_id"])}
        orphans = 0
        for chrom in chroms:
            cov, acc, n_orph = fr.region_reads(chrom, regions_bed, spans, fe, tmpdir,
                                               args.read_rule, args.fire_overlap_fraction)
            # same orphan tally as before: FIRE hits on reads that do not cover the region
            orphans += n_orph
            for rid, reads in cov.items():
                n_reads[row_of[rid]] = len(reads)
                n_fire[row_of[rid]] = len(acc.get(rid, ()))

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

    regions_bed.unlink()
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
