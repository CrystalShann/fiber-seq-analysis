#!/usr/bin/env python3
"""Per-region FIRE element and FiberHMM footprint counts across the LPS timecourse.

Same regions and fibers as 03_fire_frequency.py, so the two long tables join
one-to-one on region_id + timepoint. For every union region and timepoint:

    n_reads              distinct reads whose aligned span fully covers the region
    n_fire               spanning reads with >= 1 FIRE element covering >= 50% of the
                         region (= fire_n_reads; equals 03's n_fire)
    <k>_n_elements       feature-k elements on spanning reads
    <k>_n_reads          spanning reads with >= 1 feature-k element
    <k>_n_elements_gated / <k>_n_reads_gated
                         the same, restricted to the region's FIRE-positive reads
                         (only for features with a gate)

Features (FEATURES registry; a new feature only needs a new entry):
    fire   FIRE elements (same file as 03), element covers >= 50% of the region
    fp     FiberHMM TF footprints, score >= 50 and 10-80 bp, footprint lies 100%
           inside the region; not required to sit inside a FIRE element; gate = fire

Rates (NA when the denominator is 0):
    fire_freq                fire_n_reads / n_reads        (= 03's freq)
    fire_elements_per_read   fire_n_elements / n_reads
    fp_per_read              fp_n_elements / n_reads
    fp_per_fire_read         fp_n_elements_gated / n_fire  (main footprint frequency)
    fp_freq_given_fire       fp_n_reads_gated / n_fire
    fp_per_fire_read_per_kb  fp_per_fire_read / (width / 1000)

Outputs in --out-dir:
    feature_count_long.tsv.gz   one row per (region, timepoint)
    feature_count_wide.tsv.gz   one row per region, per-timepoint columns
"""

import argparse
import os
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

SAMPLES = ["LPS_0", "LPS_5", "LPS_10", "LPS_15"]
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

BEDTOOLS = os.environ.get("BEDTOOLS", "/project/spott/cshan/envs/bedtools/bin/bedtools")
TABIX = os.environ.get("TABIX", "/project/spott/cshan/envs/dimelo/bin/tabix")
READ_MIN_REGION_FRACTION = 1.0
FIRE_ROOT = "/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE"


def run_intersect(args_list, desc):
    p = subprocess.run(args_list, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"ERROR: {desc} failed:\n{p.stderr[:2000]}")
    return p.stdout


def tabix_awk(path, chrom, out_bed, program, awk_vars=None):
    """tabix <path> <chrom> | awk <program> > out_bed, failing if either step fails."""
    awk_cmd = ["awk"]
    for k, v in (awk_vars or {}).items():
        awk_cmd += ["-v", f"{k}={v}"]
    with open(out_bed, "w") as fh:
        tb = subprocess.Popen([TABIX, str(path), chrom], stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        aw = subprocess.Popen(awk_cmd + [program], stdin=tb.stdout, stdout=fh,
                              stderr=subprocess.PIPE, text=True)
        tb.stdout.close()
        aw_err = aw.communicate()[1]
        tb_err = tb.stderr.read().decode()
        tb.wait()
    if tb.returncode != 0 or aw.returncode != 0:
        sys.exit(f"ERROR: loading {path} {chrom} failed:\n{tb_err[:1000]}{aw_err[:1000]}")


###############################
# loaders: write a BED4 (chrom start end read_name) for one chromosome
###############################
def load_bed4(path, chrom, out_bed, args):
    """Intervals as-is, first four columns."""
    tabix_awk(path, chrom, out_bed, 'BEGIN { FS = OFS = "\\t" } { print $1, $2, $3, $4 }')


# FiberHMM BED12+3, one line per read: $2 read start, $4 read name, $10 blockCount,
# $11 sizes, $12 offsets, $13 scores (comma-separated). Every block is a real call
# (no sentinel blocks, see split_footprints_by_size.sh). $13 is the per-block tq
# (round(LLR*10), 0-255); $14/$15 are edge scores. fiberhmm-extract already drops
# tq < 50, so --min-score only bites above 50.
EXPLODE_BLOCKS = r'''
BEGIN { FS = OFS = "\t" }
{
    n = split($11, sz, ","); m = split($12, off, ","); k = split($13, sc, ",")
    if (sz[n] == "") n--; if (off[m] == "") m--; if (sc[k] == "") k--
    if (n != $10 || m != n || k != n) { print "block count mismatch: " $4 > "/dev/stderr"; exit 1 }
    for (i = 1; i <= n; i++) {
        s = sz[i] + 0
        if (s < min_size || s > max_size || sc[i] + 0 < min_score) continue
        start = $2 + off[i]
        print $1, start, start + s, $4
    }
}
'''


def load_fiberhmm_blocks(path, chrom, out_bed, args):
    """One BED4 line per block passing score >= min_score and min_size <= size <= max_size."""
    tabix_awk(path, chrom, out_bed, EXPLODE_BLOCKS,
              {"min_score": args.min_score, "min_size": args.min_size,
               "max_size": args.max_size})


def fiberhmm_label(s):
    return s.replace("_", "")       # LPS_5 -> LPS5


###############################
# feature registry. Rules are bedtools intersect flags with regions as -a:
#   -f x   element covers >= x of the region     -F x   >= x of the element is inside
# gate: also count this feature on the region's reads positive for the gate feature
# (the gate must come earlier in the dict).
###############################
FEATURES = {
    "fire": {
        "path": lambda s, chrom, args: Path(FIRE_ROOT) / s / "additional-outputs-v0.1" /
                "fire-peaks" / f"{s}-v0.1-fire-elements.bed.gz",
        "load": load_bed4,
        "flag": "-f", "min_frac": 0.5,          # identical to 03
        "gate": None,
    },
    "fp": {
        "path": lambda s, chrom, args: Path(args.tf_root) / fiberhmm_label(s) /
                f"{fiberhmm_label(s)}_hmm_extracted_tf_{chrom}.bed.gz",
        "load": load_fiberhmm_blocks,
        "flag": "-F", "min_frac": 1.0,
        "gate": "fire",
    },
}


def count_columns():
    cols = ["n_reads"]
    for k, spec in FEATURES.items():
        cols += [f"{k}_n_elements", f"{k}_n_reads"]
        if spec["gate"]:
            cols += [f"{k}_n_elements_gated", f"{k}_n_reads_gated"]
    return cols


def feature_hits(regions_bed, feature_bed, flag, min_frac):
    """{region_id: [read name per element]} for elements passing the bedtools rule.

    regions_bed is BED4 with region_id in column 4; feature_bed has the read name in
    column 4, so it is field 7 of the -wa -wb output.
    """
    hits = defaultdict(list)
    for line in run_intersect(
            [BEDTOOLS, "intersect", "-a", str(regions_bed), "-b", str(feature_bed),
             flag, str(min_frac), "-wa", "-wb"],
            f"intersect {feature_bed.name}").splitlines():
        f = line.split("\t")
        hits[f[3]].append(f[7])
    return hits


def counts_for_chrom(chrom, regions, s, spans_path, args, tmpdir):
    """({region_id: {count column: value}}, {feature: dropped hits}) for one chromosome.

    Only regions with >= 1 spanning read are returned; missing columns are 0.
    """
    sub = regions[regions["chrom"] == chrom]
    if sub.empty:
        return {}, {k: 0 for k in FEATURES}
    reg_chr = tmpdir / f"regions.{chrom}.bed"
    sub[["chrom", "start", "end", "region_id"]].to_csv(
        reg_chr, sep="\t", header=False, index=False)

    # Denominator: reads whose span covers 100% of the region (identical to 03).
    spans_chr = tmpdir / f"spans.{chrom}.bed"
    load_bed4(spans_path, chrom, spans_chr, args)
    cov = {r: set(names) for r, names in
           feature_hits(reg_chr, spans_chr, "-f", READ_MIN_REGION_FRACTION).items()}
    spans_chr.unlink()

    counts = {r: {"n_reads": len(reads)} for r, reads in cov.items()}
    positive = {}       # feature -> {region_id: spanning reads with >= 1 element}
    dropped = {}
    for k, spec in FEATURES.items():
        feat_chr = tmpdir / f"{k}.{chrom}.bed"
        spec["load"](spec["path"](s, chrom, args), chrom, feat_chr, args)
        hits = feature_hits(reg_chr, feat_chr, spec["flag"], spec["min_frac"])
        feat_chr.unlink()

        gate = positive[spec["gate"]] if spec["gate"] else None
        positive[k] = {}
        dropped[k] = 0
        for r, names in hits.items():
            span = cov.get(r, set())
            kept = [n for n in names if n in span]
            dropped[k] += len(names) - len(kept)
            if not kept:
                continue
            reads = set(kept)
            positive[k][r] = reads
            c = counts[r]
            c[f"{k}_n_elements"] = len(kept)
            c[f"{k}_n_reads"] = len(reads)
            if gate is not None:
                g = gate.get(r, set())
                c[f"{k}_n_elements_gated"] = sum(n in g for n in kept)
                c[f"{k}_n_reads_gated"] = len(reads & g)

    reg_chr.unlink()
    return counts, dropped


def add_rates(df):
    """Per-read rates; NA (never 0) where the denominator is 0."""
    def ratio(num, den):
        return df[num] / df[den].where(df[den] > 0)
    df["fire_freq"] = ratio("fire_n_reads", "n_reads")
    df["fire_elements_per_read"] = ratio("fire_n_elements", "n_reads")
    df["fp_per_read"] = ratio("fp_n_elements", "n_reads")
    df["fp_per_fire_read"] = ratio("fp_n_elements_gated", "n_fire")
    df["fp_freq_given_fire"] = ratio("fp_n_reads_gated", "n_fire")
    df["fp_per_fire_read_per_kb"] = df["fp_per_fire_read"] / (df["width"] / 1000)
    return df


def check(df):
    for k, spec in FEATURES.items():
        assert (df[f"{k}_n_reads"] <= df["n_reads"]).all(), k
        assert (df[f"{k}_n_reads"] <= df[f"{k}_n_elements"]).all(), k
        if spec["gate"]:
            assert (df[f"{k}_n_elements_gated"] <= df[f"{k}_n_elements"]).all(), k
            assert (df[f"{k}_n_reads_gated"] <= df[f"{k}_n_reads"]).all(), k
    for col, den in [("fire_freq", "n_reads"), ("fp_freq_given_fire", "n_fire")]:
        assert df[col].dropna().between(0, 1).all(), col
        assert (df[col].isna() == (df[den] == 0)).all(), col


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--union-bed",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/universe/fire_peaks_union.bed")
    ap.add_argument("--spans-root",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency",
                    help="read spans at <spans-root>/<s>/<s>.read_spans.bed.gz")
    ap.add_argument("--tf-root",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_tf",
                    help="FiberHMM calls at <tf-root>/LPS<t>/LPS<t>_hmm_extracted_tf_<chrom>.bed.gz")
    ap.add_argument("--out-dir",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency")
    ap.add_argument("--timepoints", nargs="+", default=SAMPLES)
    ap.add_argument("--chrom", default=None, help="restrict to one chromosome (testing)")
    ap.add_argument("--min-score", type=float, default=50, help="FiberHMM block score >= this")
    ap.add_argument("--min-size", type=int, default=10, help="footprint size >= this (bp)")
    ap.add_argument("--max-size", type=int, default=80, help="footprint size <= this (bp)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    regions = pd.read_csv(args.union_bed, sep="\t", header=None,
                          names=["chrom", "start", "end"])
    n_union_regions = len(regions)
    chroms = [args.chrom] if args.chrom else CHROMS
    regions = regions[regions["chrom"].isin(chroms)].reset_index(drop=True)
    regions["region_id"] = (regions["chrom"] + ":" + regions["start"].astype(str)
                            + "-" + regions["end"].astype(str))
    assert regions["region_id"].is_unique
    regions["width"] = regions["end"] - regions["start"]
    regions = regions[["region_id", "chrom", "start", "end", "width"]]
    region_chroms = [c for c in chroms if (regions["chrom"] == c).any()]
    print(f"union regions: {len(regions)} of {n_union_regions} in {args.union_bed} "
          f"on {len(region_chroms)} chromosomes", flush=True)

    # fail before any counting if an input is missing
    for s in args.timepoints:
        need = [Path(args.spans_root) / s / f"{s}.read_spans.bed.gz"]
        need += [spec["path"](s, c, args) for spec in FEATURES.values() for c in region_chroms]
        for f in need:
            if not f.is_file():
                sys.exit(f"ERROR: missing input {f}")

    cols = count_columns()
    row_of = {rid: i for i, rid in enumerate(regions["region_id"])}
    long_parts = []
    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix=".tmp_feature_count_", dir=out_dir) as tmp:
        tmpdir = Path(tmp)
        for s in args.timepoints:
            spans = Path(args.spans_root) / s / f"{s}.read_spans.bed.gz"
            arr = {c: np.zeros(len(regions), dtype=np.int64) for c in cols}
            dropped = dict.fromkeys(FEATURES, 0)
            for chrom in region_chroms:
                counts, drop = counts_for_chrom(chrom, regions, s, spans, args, tmpdir)
                for k in FEATURES:
                    dropped[k] += drop[k]
                for rid, c in counts.items():
                    i = row_of[rid]
                    for col, v in c.items():
                        arr[col][i] = v
                print(f"  {s} {chrom} done ({time.time() - t0:.0f}s)", flush=True)

            part = regions.copy()
            part["timepoint"] = s
            part["n_reads"] = arr["n_reads"]
            part["n_fire"] = arr["fire_n_reads"]
            for col in cols[1:]:
                part[col] = arr[col]
            part = add_rates(part)
            check(part)
            long_parts.append(part)
            covered = part["n_reads"] > 0
            print(f"{s}: dropped non-spanning hits "
                  + ", ".join(f"{k} {v}" for k, v in dropped.items())
                  + f"; zero-coverage regions {(~covered).sum()}"
                  f"; zero-FIRE regions (covered) {(covered & (part['n_fire'] == 0)).sum()}",
                  flush=True)

    long_df = pd.concat(long_parts, ignore_index=True)
    assert len(long_df) == len(regions) * len(args.timepoints)
    long_out = out_dir / "feature_count_long.tsv.gz"
    long_df.to_csv(long_out, sep="\t", index=False, float_format="%.6g", na_rep="NA")
    print(f"wrote {long_out} ({len(long_df)} rows)", flush=True)

    # wide pivot, keeping the genomic order of the union bed
    value_cols = [c for c in long_df.columns if c not in regions.columns and c != "timepoint"]
    wide = regions.copy()
    for s, part in zip(args.timepoints, long_parts):
        for c in value_cols:
            wide[f"{s}_{c}"] = part[c].to_numpy()
    wide_out = out_dir / "feature_count_wide.tsv.gz"
    wide.to_csv(wide_out, sep="\t", index=False, float_format="%.6g", na_rep="NA")
    print(f"wrote {wide_out} ({len(wide)} rows)", flush=True)


if __name__ == "__main__":
    main()
