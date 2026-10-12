#!/usr/bin/env python3
"""Per-read FiberHMM TF footprint counts and per-region footprint-count distributions.

Same regions, spanning reads, FIRE rule and footprint filters as 03b_feature_count.py,
whose loaders and intersect rules are imported, not copied. For every union region,
timepoint and spanning read:

    fire_positive   1 when a FIRE element covers >= 50% of the region (03b's rule), else 0
    fp_count_read   FiberHMM TF footprints (score >= 50, 10-80 bp) lying 100% inside the
                    region on this read; 0 when there are none. Blocks are never merged,
                    and a record repeated on the same read (same start and end) counts once

Per region and timepoint, over all spanning reads (no suffix), the FIRE-positive reads
(_fire) and the FIRE-negative reads (_nonfire):

    n_reads                 reads in the group (n_reads_fire = n_fire)
    n_fp0 n_fp1 n_fp2 n_fp3plus   reads with 0, 1, 2 and >= 3 footprints
    p_fp0 ... p_fp3plus     n_fp* / n_reads; NA when n_reads == 0
    mean_fp_per_read        mean fp_count_read; NA when n_reads == 0
    var_fp_per_read         sample variance (n - 1) of fp_count_read; NA when n_reads < 2
    fp_count_all            sum of fp_count_read (all reads: = 03b's fp_n_elements)

Every region gets a row, also with no spanning reads. The counts are checked against
--feature-count (03b's feature_count_long.tsv.gz): n_reads, n_fire, fp_n_elements,
fp_n_reads, fp_n_elements_gated and fp_n_reads_gated. Disagreements are written out, one
row per region, timepoint and column, and never dropped.

If regions overlap, a read or footprint counts toward every region it qualifies for.

Outputs in --out-dir (default: the project fire_frequency folder, or its
ft_dist_test_<chrom> subfolder with --chrom):
    footprint_count_per_read.tsv.gz           one row per (region, timepoint, read);
                                              not written with --no-per-read
    footprint_distribution_long.tsv.gz        one row per (region, timepoint)
    footprint_distribution_mismatches.tsv.gz  disagreements with --feature-count
    footprint_distribution_qc.tsv             per timepoint and chromosome: intersect
                                              records, what was removed and why
"""

import argparse
import csv
import gzip
import importlib.util
import io
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

# 03b's loaders, input paths and intersect rules, so both scripts count the same things
# (no __pycache__ in the code folder)
sys.dont_write_bytecode = True
_spec = importlib.util.spec_from_file_location(
    "feature_count", Path(__file__).resolve().parent / "03b_feature_count.py")
fc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fc)

PROJECT_DIR = "/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency"
REGION_COLS = ["region_id", "chrom", "start", "end", "width"]
PER_READ_COLS = REGION_COLS + ["timepoint", "read_id", "fire_positive", "fp_count_read"]
CATEGORIES = ["n_fp0", "n_fp1", "n_fp2", "n_fp3plus"]
GROUPS = {"": None, "_fire": 1, "_nonfire": 0}      # column suffix -> fire_positive value
KEY = ["region_id", "read_id"]

# our column -> the 03b column it must equal
REFERENCE = {"n_reads": "n_reads", "n_fire": "n_fire", "fp_count_all": "fp_n_elements",
             "fp_n_reads": "fp_n_reads", "fp_count_all_fire": "fp_n_elements_gated",
             "fp_n_reads_fire": "fp_n_reads_gated"}
REMOVED_COLS = ["span_duplicate_records", "fire_nonspanning_records",
                "fp_duplicate_records", "fp_nonspanning_records"]
MISMATCH_COLS = ["region_id", "timepoint", "column", "ft_dist", "feature_count"] + REMOVED_COLS


def distribution_columns():
    cols = ["n_reads", "n_fire"]
    for sfx in GROUPS:
        if sfx:
            cols.append("n_reads" + sfx)
        cols += [c + sfx for c in CATEGORIES]
        cols += ["p_" + c[2:] + sfx for c in CATEGORIES]
        cols += [c + sfx for c in ["mean_fp_per_read", "var_fp_per_read", "fp_count_all"]]
    return cols


def spans_path(args, s):
    return Path(args.spans_root) / s / f"{s}.read_spans.bed.gz"


def intersect(regions_bed, feature_bed, flag, min_frac, cols):
    """bedtools intersect with the regions as -a (-wa -wb) as a DataFrame.

    cols maps a column name to its output field; the feature beds are BED4 with the
    read name in column 4, so the feature's fields are 4-7.
    """
    out = fc.run_intersect(
        [fc.BEDTOOLS, "intersect", "-a", str(regions_bed), "-b", str(feature_bed),
         flag, str(min_frac), "-wa", "-wb"], f"intersect {feature_bed.name}")
    if not out:
        return pd.DataFrame({c: pd.Series(dtype=object if c in KEY else "int64")
                             for c in cols})
    df = pd.read_csv(io.StringIO(out), sep="\t", header=None, usecols=list(cols.values()),
                     quoting=csv.QUOTE_NONE,
                     dtype={f: str for c, f in cols.items() if c in KEY})
    return df.rename(columns={f: c for c, f in cols.items()})[list(cols)]


def n_overlapping(sub):
    """Regions that overlap an earlier region of the chromosome, sorted by start."""
    s = sub.sort_values(["start", "end"])
    return int((s["start"] < s["end"].cummax().shift()).sum())


###############################
# one row per spanning read for one chromosome and timepoint
###############################
def per_read_for_chrom(chrom, sub, s, args, tmpdir):
    """(per-read table, qc counts, per-region removed records) for one chromosome."""
    reg = tmpdir / f"regions.{chrom}.bed"
    sub[["chrom", "start", "end", "region_id"]].to_csv(reg, sep="\t", header=False,
                                                       index=False)
    fire_spec, fp_spec = fc.FEATURES["fire"], fc.FEATURES["fp"]
    spans_bed = tmpdir / f"spans.{chrom}.bed"
    fire_bed = tmpdir / f"fire.{chrom}.bed"
    fp_bed = tmpdir / f"fp.{chrom}.bed"
    fc.load_bed4(spans_path(args, s), chrom, spans_bed, args)
    fire_spec["load"](fire_spec["path"](s, chrom, args), chrom, fire_bed, args)
    fp_spec["load"](fp_spec["path"](s, chrom, args), chrom, fp_bed, args)

    ids = {"region_id": 3, "read_id": 7}
    span = intersect(reg, spans_bed, "-f", fc.READ_MIN_REGION_FRACTION, ids)
    fire = intersect(reg, fire_bed, fire_spec["flag"], fire_spec["min_frac"], ids)
    fp = intersect(reg, fp_bed, fp_spec["flag"], fp_spec["min_frac"],
                   {"region_id": 3, "fp_start": 5, "fp_end": 6, "read_id": 7})
    for f in (reg, spans_bed, fire_bed, fp_bed):
        f.unlink()

    # a read listed twice for a region, or a footprint record repeated on a read,
    # counts once; FIRE elements and footprints on non-spanning reads are excluded
    span_dup = span.duplicated()
    spanning = pd.MultiIndex.from_frame(span.loc[~span_dup, KEY])
    fire_on_span = pd.MultiIndex.from_frame(fire[KEY]).isin(spanning)
    fp_dup = fp.duplicated()
    fp_on_span = pd.MultiIndex.from_frame(fp[KEY]).isin(spanning)
    fp_kept = fp[~fp_dup & fp_on_span]

    fire_reads = fire.loc[fire_on_span, KEY].drop_duplicates().assign(fire_positive=1)
    fp_counts = fp_kept.groupby(KEY).size().rename("fp_count_read").reset_index()
    per_read = (span[~span_dup]
                .merge(fire_reads, on=KEY, how="left")
                .merge(fp_counts, on=KEY, how="left")
                .fillna({"fire_positive": 0, "fp_count_read": 0})
                .astype({"fire_positive": "int64", "fp_count_read": "int64"}))
    assert not per_read.duplicated(KEY).any()
    # footprint multiplicity summed over reads = footprints counted
    assert per_read["fp_count_read"].sum() == len(fp_kept)

    removed = pd.concat({
        "span_duplicate_records": span.loc[span_dup, "region_id"].value_counts(),
        "fire_nonspanning_records": fire.loc[~fire_on_span, "region_id"].value_counts(),
        "fp_duplicate_records": fp.loc[fp_dup, "region_id"].value_counts(),
        "fp_nonspanning_records": fp.loc[~fp_dup & ~fp_on_span, "region_id"].value_counts(),
    }, axis=1).reindex(columns=REMOVED_COLS).fillna(0).astype("int64")

    qc = {"timepoint": s, "chrom": chrom, "n_regions": len(sub),
          "n_overlapping_regions": n_overlapping(sub),
          "spanning_records": len(span), "span_duplicate_records": int(span_dup.sum()),
          "fire_records": len(fire),
          "fire_nonspanning_records": int((~fire_on_span).sum()),
          "fp_records": len(fp), "fp_duplicate_records": int(fp_dup.sum()),
          "fp_nonspanning_records": int((~fp_dup & ~fp_on_span).sum()),
          "fp_counted": len(fp_kept)}
    return per_read, qc, removed


###############################
# per-region distributions and checks
###############################
def distribution(per_read, region_ids):
    """Distribution columns over all, FIRE-positive and FIRE-negative reads, one row per
    region_id; regions without reads get zero counts and NA rates."""
    x = per_read["fp_count_read"]
    flags = pd.DataFrame({"region_id": per_read["region_id"],
                          "fire_positive": per_read["fire_positive"], "fp": x,
                          "n_fp0": x == 0, "n_fp1": x == 1, "n_fp2": x == 2,
                          "n_fp3plus": x >= 3})
    parts = []
    for sfx, fire in GROUPS.items():
        d = flags if fire is None else flags[flags["fire_positive"] == fire]
        g = d.groupby("region_id")
        out = g[CATEGORIES].sum()
        out.insert(0, "n_reads", g.size())
        out["fp_count_all"] = g["fp"].sum()
        out = out.reindex(region_ids).fillna(0).astype("int64")
        n = out["n_reads"].where(out["n_reads"] > 0)
        for c in CATEGORIES:
            out["p_" + c[2:]] = out[c] / n
        # from the per-read counts, not the collapsed categories
        out["mean_fp_per_read"] = g["fp"].mean().reindex(region_ids)
        out["var_fp_per_read"] = g["fp"].var(ddof=1).reindex(region_ids)
        parts.append(out.add_suffix(sfx))
    dist = pd.concat(parts, axis=1)
    dist["n_fire"] = dist["n_reads_fire"]
    return dist[distribution_columns()].rename_axis("region_id").reset_index()


def check(dist):
    """Internal consistency of one distribution chunk; a failure stops the run."""
    assert (dist["n_fire"] <= dist["n_reads"]).all()
    assert (dist["n_reads_fire"] + dist["n_reads_nonfire"] == dist["n_reads"]).all()
    assert (dist["fp_count_all_fire"] + dist["fp_count_all_nonfire"]
            == dist["fp_count_all"]).all()
    for sfx in GROUPS:
        n = dist["n_reads" + sfx]
        assert (dist[[c + sfx for c in CATEGORIES]].sum(axis=1) == n).all(), sfx
        p = dist[["p_" + c[2:] + sfx for c in CATEGORIES]]
        assert (p.isna() | ((p >= 0) & (p <= 1))).all().all(), sfx
        assert (p.isna().all(axis=1) == (n == 0)).all(), sfx
        assert (p.notna().all(axis=1) == (n > 0)).all(), sfx
        assert np.allclose(p[n > 0].sum(axis=1), 1), sfx
        mean = dist["mean_fp_per_read" + sfx]
        assert (mean.isna() == (n == 0)).all(), sfx
        assert (dist["var_fp_per_read" + sfx].isna() == (n < 2)).all(), sfx
        assert np.allclose((mean * n)[n > 0], dist["fp_count_all" + sfx][n > 0]), sfx


def compare(dist, ref, removed, s):
    """One row per region and 03b column where dist and 03b disagree; regions missing
    from either table are mismatches too."""
    ours = dist.set_index("region_id")
    ours = ours.assign(fp_n_reads=ours["n_reads"] - ours["n_fp0"],
                       fp_n_reads_fire=ours["n_reads_fire"] - ours["n_fp0_fire"])
    index = ours.index.union(ref.index, sort=False)
    ours, ref = ours.reindex(index), ref.reindex(index)
    rows = []
    for col, ref_col in REFERENCE.items():
        bad = (~(ours[col] == ref[ref_col])).to_numpy()   # NA on either side mismatches
        if bad.any():
            rows.append(pd.DataFrame({"region_id": index[bad], "timepoint": s,
                                      "column": ref_col,
                                      "ft_dist": ours[col].to_numpy()[bad],
                                      "feature_count": ref[ref_col].to_numpy()[bad]}))
    if not rows:
        return pd.DataFrame(columns=MISMATCH_COLS)
    out = pd.concat(rows, ignore_index=True)
    return out.join(removed.reindex(out["region_id"]).fillna(0).astype("int64")
                    .reset_index(drop=True))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--union-bed", default=f"{PROJECT_DIR}/universe/fire_peaks_union.bed")
    ap.add_argument("--spans-root", default=PROJECT_DIR,
                    help="read spans at <spans-root>/<s>/<s>.read_spans.bed.gz")
    ap.add_argument("--tf-root",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_tf",
                    help="FiberHMM calls at <tf-root>/LPS<t>/LPS<t>_hmm_extracted_tf_<chrom>.bed.gz")
    ap.add_argument("--feature-count", default=f"{PROJECT_DIR}/feature_count_long.tsv.gz",
                    help="03b_feature_count.py long table to validate against")
    ap.add_argument("--out-dir", default=None,
                    help=f"default {PROJECT_DIR}, or its ft_dist_test_<chrom> with --chrom")
    ap.add_argument("--timepoints", nargs="+", default=fc.SAMPLES)
    ap.add_argument("--chrom", default=None, help="restrict to one chromosome (testing)")
    ap.add_argument("--min-score", type=float, default=50, help="FiberHMM block score >= this")
    ap.add_argument("--min-size", type=int, default=10, help="footprint size >= this (bp)")
    ap.add_argument("--max-size", type=int, default=80, help="footprint size <= this (bp)")
    ap.add_argument("--no-per-read", action="store_true",
                    help="write only the distributions, not footprint_count_per_read.tsv.gz")
    args = ap.parse_args()

    if args.out_dir:
        out_dir = Path(args.out_dir)
    elif args.chrom:
        out_dir = Path(PROJECT_DIR) / f"ft_dist_test_{args.chrom}"
    else:
        out_dir = Path(PROJECT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"output: {out_dir}", flush=True)

    # the regions exactly as 03b builds them
    try:
        regions = pd.read_csv(args.union_bed, sep="\t", header=None,
                              names=["chrom", "start", "end"])
    except pd.errors.EmptyDataError:
        regions = pd.DataFrame({"chrom": pd.Series(dtype=str),
                                "start": pd.Series(dtype="int64"),
                                "end": pd.Series(dtype="int64")})
    n_union_regions = len(regions)
    chroms = [args.chrom] if args.chrom else fc.CHROMS
    regions = regions[regions["chrom"].isin(chroms)].reset_index(drop=True)
    regions["region_id"] = (regions["chrom"] + ":" + regions["start"].astype(str)
                            + "-" + regions["end"].astype(str))
    if not regions["region_id"].is_unique:
        sys.exit(f"ERROR: duplicate intervals in {args.union_bed}; region_id is not unique")
    regions["width"] = regions["end"] - regions["start"]
    regions = regions[REGION_COLS]
    region_chroms = [c for c in chroms if (regions["chrom"] == c).any()]
    print(f"union regions: {len(regions)} of {n_union_regions} in {args.union_bed} "
          f"on {len(region_chroms)} chromosomes", flush=True)
    if not region_chroms:
        print(f"NOTE: no regions on {', '.join(chroms)}; writing empty tables", flush=True)

    # fail before any counting if an input is missing
    need = [Path(args.feature_count)]
    for s in args.timepoints:
        need.append(spans_path(args, s))
        need += [fc.FEATURES[k]["path"](s, c, args) for k in ("fire", "fp")
                 for c in region_chroms]
    for f in need:
        if not f.is_file():
            sys.exit(f"ERROR: missing input {f}")

    ref = pd.read_csv(args.feature_count, sep="\t",
                      usecols=["region_id", "chrom", "timepoint", *REFERENCE.values()])
    ref = ref[ref["chrom"].isin(region_chroms) & ref["timepoint"].isin(args.timepoints)]
    if ref.duplicated(["region_id", "timepoint"]).any():
        sys.exit(f"ERROR: duplicate region_id x timepoint rows in {args.feature_count}")
    ref_of = {k: g.set_index("region_id") for k, g in ref.groupby(["timepoint", "chrom"])}
    empty_ref = ref.iloc[:0].set_index("region_id")

    row_of = pd.Series(np.arange(len(regions)), index=regions["region_id"])
    per_read_out = out_dir / "footprint_count_per_read.tsv.gz"
    per_read_fh = None
    if not args.no_per_read:
        per_read_fh = gzip.open(per_read_out, "wt", compresslevel=6)
        per_read_fh.write("\t".join(PER_READ_COLS) + "\n")
    dist_parts, mismatch_parts, qc_rows = [], [], []
    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix=".tmp_ft_dist_", dir=out_dir) as tmp:
        tmpdir = Path(tmp)
        for s in args.timepoints:
            for chrom in region_chroms:
                sub = regions[regions["chrom"] == chrom]
                per_read, qc, removed = per_read_for_chrom(chrom, sub, s, args, tmpdir)
                if qc["n_overlapping_regions"]:
                    print(f"  NOTE: {qc['n_overlapping_regions']} {chrom} regions overlap "
                          "another; their reads and footprints count in each", flush=True)

                dist = distribution(per_read, sub["region_id"])
                check(dist)
                mism = compare(dist, ref_of.get((s, chrom), empty_ref), removed, s)
                qc["mismatch_regions"] = mism["region_id"].nunique()
                qc_rows.append(qc)
                if len(mism):
                    mismatch_parts.append(mism)
                dist_parts.append(sub.merge(dist, on="region_id").assign(timepoint=s))

                if per_read_fh is not None:
                    per_read = per_read.merge(sub, on="region_id").assign(timepoint=s)
                    per_read["_row"] = per_read["region_id"].map(row_of)
                    per_read.sort_values(["_row", "read_id"])[PER_READ_COLS].to_csv(
                        per_read_fh, sep="\t", header=False, index=False)
                print(f"  {s} {chrom}: {len(per_read)} spanning read records, "
                      f"{qc['fp_counted']} footprints, {qc['mismatch_regions']} mismatching "
                      f"regions ({time.time() - t0:.0f}s)", flush=True)
    if per_read_fh is not None:
        per_read_fh.close()
        print(f"wrote {per_read_out}", flush=True)

    # genomic order of the union bed within each timepoint, as feature_count_long
    cols = REGION_COLS + ["timepoint"] + distribution_columns()
    if dist_parts:
        long_df = pd.concat(dist_parts, ignore_index=True)
        long_df["_t"] = long_df["timepoint"].map({s: i for i, s in enumerate(args.timepoints)})
        long_df["_row"] = long_df["region_id"].map(row_of)
        long_df = long_df.sort_values(["_t", "_row"])[cols]
    else:
        long_df = pd.DataFrame(columns=cols)
    assert len(long_df) == len(regions) * len(args.timepoints)
    assert not long_df.duplicated(["region_id", "timepoint"]).any()
    long_out = out_dir / "footprint_distribution_long.tsv.gz"
    long_df.to_csv(long_out, sep="\t", index=False, float_format="%.6g", na_rep="NA")
    print(f"wrote {long_out} ({len(long_df)} rows)", flush=True)

    mismatches = (pd.concat(mismatch_parts, ignore_index=True) if mismatch_parts else
                  pd.DataFrame(columns=MISMATCH_COLS))
    mism_out = out_dir / "footprint_distribution_mismatches.tsv.gz"
    mismatches.to_csv(mism_out, sep="\t", index=False, na_rep="NA")
    qc_df = pd.DataFrame(qc_rows, columns=["timepoint", "chrom", "n_regions",
                                           "n_overlapping_regions", "spanning_records",
                                           "span_duplicate_records", "fire_records",
                                           "fire_nonspanning_records", "fp_records",
                                           "fp_duplicate_records", "fp_nonspanning_records",
                                           "fp_counted", "mismatch_regions"])
    qc_out = out_dir / "footprint_distribution_qc.tsv"
    qc_df.to_csv(qc_out, sep="\t", index=False)
    print(f"wrote {mism_out} ({len(mismatches)} rows) and {qc_out}", flush=True)

    for s, g in qc_df.groupby("timepoint", sort=False):
        print(f"{s}: removed span duplicates {g['span_duplicate_records'].sum()}, "
              f"FIRE hits on non-spanning reads {g['fire_nonspanning_records'].sum()}, "
              f"footprint duplicates {g['fp_duplicate_records'].sum()}, "
              f"footprints on non-spanning reads {g['fp_nonspanning_records'].sum()}; "
              f"mismatching regions vs feature_count {g['mismatch_regions'].sum()}",
              flush=True)
    if len(mismatches):
        print(f"WARNING: {len(mismatches)} disagreements with {args.feature_count}; "
              f"see {mism_out}", flush=True)


if __name__ == "__main__":
    main()
