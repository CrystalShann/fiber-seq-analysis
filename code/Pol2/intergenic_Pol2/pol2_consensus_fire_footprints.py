#!/usr/bin/env python3
"""FIRE accessibility and 30-60 bp footprints at Langen et al. 2023 intergenic
RNAPII consensus peaks, one LCL sample at a time.

Counts use a fixed window of centroid +/- HALF_WINDOW bp (the same 500-bp windows as
intergenic_vs_tss_pol2_chromatin_states.r), not the peak interval, so footprint frequency
does not depend on peak width. For every consensus peak (chr1-22, X, Y):
    n_reads      reads covering the window (--read-rule span: read span contains it;
                 aligned: aligned blocks cover every base)
    n_fire       of those, reads with ONE FIRE element covering >= --fire-overlap-fraction
    n_fp         of those, reads with >= 1 FiberHMM v2 footprint block of --fp-sizes bp
                 with >= --fp-fraction of the footprint inside the window (1.0 = entirely inside)
    n_fp_blocks  number of such footprint blocks (on n_reads reads)

Usage:
    python pol2_consensus_fire_footprints.py --sample AL6_bc2162_19093 [--out-dir <dir>]
    python pol2_consensus_fire_footprints.py --summarize [--out-dir <dir>]

Summary (pol2_consensus_summary.tsv), per group (all, annotated, each category):
    frac_accessible                 peaks (peak interval) >= --fire-overlap-fraction inside an
                                    LCL FIRE peak
    median_fire_freq, median_fp_freq  among accessible peaks (window counts)
    frac_accessible_in_pol2_chip    accessible peaks with a GM12878 Pol II ChIP IDR peak summit
                                    (ENCFF886PSD, strongest summit per peak) inside the window,
                                    i.e. within HALF_WINDOW bp of the centroid
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "parsing_functions"))
import fiber_regions as fr

ATLAS = Path("/project/spott/cshan/annotations/LCL_public_data/intergenic_pol2_Langen_2023")
MANIFEST = Path("/project/spott/cshan/fiber-seq/LCL_project/co-accessibility/"
                "CRE_co-accessibility_FIRE/universe/sample_manifest.tsv")
FIRE_UNION = MANIFEST.parent / "fire_peaks_union.bed"
FP_ROOT = Path("/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/FiberHMM_v2/"
               "results/results_v1_trained_model")
# GM12878 RNAPII ChIP-seq IDR thresholded peaks, ENCSR000EAD (8WG16; Day et al. 2016 reference)
POL2_CHIP = Path("/project/spott/cshan/annotations/LCL_public_data/Pol2_chip/ENCFF886PSD.bed.gz")
HALF_WINDOW = 250


def load_consensus():
    """consensuses.bed (IDs) joined row-wise to the annotated bed (category;biotype)."""
    cons = pd.read_csv(ATLAS / "pol2_consensuses" / "consensuses.bed", sep="\t", header=None,
                       usecols=[0, 1, 2, 3, 4, 6], names=["chrom", "start", "end", "cid", "n_datasets", "center"])
    ann = pd.read_csv(ATLAS / "RNAP2_consensuses_annotated.bed", sep="\t", header=None,
                      usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "label"])
    assert (cons[["chrom", "start", "end"]].values == ann[["chrom", "start", "end"]].values).all()
    cons["category"] = ann["label"].str.split(";").str[0]
    cons["biotype"] = ann["label"].str.split(";", n=1).str[1]
    cons = cons[cons["chrom"].isin(fr.CHROMS)].reset_index(drop=True)
    cons["region_id"] = "pol2_" + cons["cid"].astype(str)
    cons["win_start"] = (cons["center"] - HALF_WINDOW).clip(lower=0)
    cons["win_end"] = cons["center"] + HALF_WINDOW
    return cons


def run_sample(sample, out_dir, args):
    man = pd.read_csv(MANIFEST, sep="\t").set_index("sample_name")
    reads = man.at[sample, "blocks_path" if args.read_rule == "aligned" else "spans_path"]
    fe = man.at[sample, "elements_path"]
    fp = FP_ROOT / sample / f"{sample}.recalled_tf.bed.gz"
    regions = load_consensus()
    tmpdir = out_dir / f".tmp_{sample}"
    tmpdir.mkdir(parents=True, exist_ok=True)
    regions_bed = tmpdir / "regions.bed"
    regions[["chrom", "win_start", "win_end", "region_id"]].to_csv(regions_bed, sep="\t", header=False, index=False)

    cols = {k: np.zeros(len(regions), dtype=np.int64) for k in ("n_reads", "n_fire", "n_fp", "n_fp_blocks")}
    row_of = {rid: i for i, rid in enumerate(regions["region_id"])}
    for chrom in fr.CHROMS:
        if not (regions["chrom"] == chrom).any():
            continue
        cov, acc, _ = fr.region_reads(chrom, regions_bed, reads, fe, tmpdir,
                                      args.read_rule, args.fire_overlap_fraction)
        fps = fr.footprints_in_regions(chrom, regions_bed, fp, tmpdir, sizes=tuple(args.fp_sizes),
                                       fp_fraction=args.fp_fraction)
        for rid, rd in cov.items():
            i = row_of[rid]
            cols["n_reads"][i] = len(rd)
            cols["n_fire"][i] = len(acc.get(rid, ()))
            hits = {r: n for r, n in fps.get(rid, {}).items() if r in rd}
            cols["n_fp"][i] = len(hits)
            cols["n_fp_blocks"][i] = sum(hits.values())
        print(f"{sample} {chrom} done", flush=True)
    regions_bed.unlink()
    tmpdir.rmdir()
    out = regions[["region_id"]].assign(sample=sample, **cols)
    assert (out["n_fire"] <= out["n_reads"]).all() and (out["n_fp"] <= out["n_reads"]).all()
    (out_dir / "per_sample").mkdir(exist_ok=True)
    out.to_csv(out_dir / "per_sample" / f"{sample}.tsv.gz", sep="\t", index=False)


def overlapping(regions, bed, fraction=None):
    """region_ids with >= fraction of the region inside `bed` (None = >= 1 bp)."""
    cmd = [fr.BEDTOOLS, "intersect", "-u", "-a", "stdin", "-b", str(bed)]
    if fraction is not None:
        cmd[3:3] = ["-f", str(fraction)]
    p = subprocess.run(cmd, input=regions[["chrom", "start", "end", "region_id"]].to_csv(
                           sep="\t", header=False, index=False),
                       capture_output=True, text=True, check=True)
    return regions["region_id"].isin({l.split("\t")[3] for l in p.stdout.splitlines()})


def chip_summits():
    """GM12878 Pol II ChIP summits; narrowPeak repeats a peak once per summit: keep the strongest."""
    pk = pd.read_csv(POL2_CHIP, sep="\t", header=None, usecols=[0, 1, 2, 6, 9],
                     names=["chrom", "start", "end", "signal", "summit"])
    pk = pk.sort_values("signal", ascending=False).drop_duplicates(["chrom", "start", "end"])
    return pk.assign(pos=pk["start"] + pk["summit"])[["chrom", "pos"]]


def summit_in_window(regions, summits):
    """True when a summit lies in [win_start, win_end), i.e. within HALF_WINDOW bp of the centroid."""
    hit = np.zeros(len(regions), dtype=bool)
    for chrom, idx in regions.groupby("chrom").groups.items():
        pos = np.sort(summits.loc[summits["chrom"] == chrom, "pos"].to_numpy())
        r = regions.loc[idx]
        hit[idx.to_numpy()] = (np.searchsorted(pos, r["win_end"].to_numpy())
                               > np.searchsorted(pos, r["win_start"].to_numpy()))
    return hit


def summarize(out_dir, fire_fraction):
    regions = load_consensus()
    # accessible: >= fire_fraction of the consensus peak inside the 31-sample LCL FIRE peak union
    regions["accessible"] = overlapping(regions, FIRE_UNION, fire_fraction)
    regions["in_pol2_chip"] = summit_in_window(regions, chip_summits())

    files = sorted((out_dir / "per_sample").glob("*.tsv.gz"))
    pooled = (pd.concat([pd.read_csv(f, sep="\t") for f in files])
              .groupby("region_id")[["n_reads", "n_fire", "n_fp", "n_fp_blocks"]].sum())
    d = regions.join(pooled, on="region_id")
    d["fire_freq"] = d["n_fire"] / d["n_reads"].where(d["n_reads"] > 0)
    d["fp_freq"] = d["n_fp"] / d["n_reads"].where(d["n_reads"] > 0)
    print(f"samples pooled: {len(files)}; consensus peaks chr1-22,X,Y: {len(d)}")
    d.to_csv(out_dir / "pooled_per_peak.tsv.gz", sep="\t", index=False)

    def tab(g):
        acc = g[g["accessible"]]
        return pd.Series({
            "n_peaks": len(g),
            "frac_accessible": g["accessible"].mean(),
            "n_accessible": len(acc),
            "median_fire_freq": acc["fire_freq"].median(),
            "median_fp_freq": acc["fp_freq"].median(),
            "frac_accessible_in_pol2_chip": acc["in_pol2_chip"].mean(),
        })
    groups = [tab(d).rename("All")]
    groups.append(tab(d[d["category"] != "Unannotated"]).rename("Annotated"))
    cat = d.assign(category=d["category"].str.split("/")).explode("category")
    groups += [tab(g).rename(c) for c, g in cat.groupby("category")]
    table = pd.DataFrame(groups).rename_axis("group").reset_index()
    for c in [c for c in table if c.startswith("n_")]:
        table[c] = table[c].astype(int)
    table.to_csv(out_dir / "pol2_consensus_summary.tsv", sep="\t", index=False, float_format="%.4f")
    pd.set_option("display.width", 250, "display.max_columns", 20)
    print(table.round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample")
    ap.add_argument("--summarize", action="store_true")
    ap.add_argument("--out-dir", default="/project/spott/cshan/fiber-seq/LCL_project/Pol2/intergenic_Pol2/tables")
    ap.add_argument("--read-rule", choices=fr.READ_RULES, default="span")
    ap.add_argument("--fire-overlap-fraction", type=float, default=0.5)
    ap.add_argument("--fp-sizes", type=int, nargs=2, default=[30, 60], help="inclusive bp range")
    ap.add_argument("--fp-fraction", type=float, default=1.0,
                    help="minimum fraction of the footprint inside the window (1.0 = entirely inside)")
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.summarize:
        summarize(out_dir, args.fire_overlap_fraction)
    else:
        run_sample(args.sample, out_dir, args)


if __name__ == "__main__":
    main()
