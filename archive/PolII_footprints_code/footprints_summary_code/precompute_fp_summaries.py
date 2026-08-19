#!/usr/bin/env python3
"""Precompute per-read footprint summaries for fast plotting.

Reads all FiberHMM footprint tabix files and the FIRE nucleosome parsed bed
files, classifies each footprint, computes its midpoint, and writes one
gzipped TSV per (fp_class, chromosome) combination.

Rows are streamed directly to disk — no full-chromosome in-memory accumulation
— so the script stays memory-efficient on large genome-wide tabix files.

Output layout
-------------
  <out_dir>/
  ├── PPP/                        25-55 bp footprints
  │   └── {chrom}_footprints.tsv.gz
  ├── PIC/                        61-80 bp footprints
  │   └── {chrom}_footprints.tsv.gz
  ├── nucleosome_140-160/         140-160 bp footprints
  │   └── {chrom}_footprints.tsv.gz
  ├── FIRE_nucleosome/            FIRE parsed nucleosome calls
  │   └── {chrom}_footprints.tsv.gz
  └── unknown/                    everything else (<25, 56-60, 81-139 bp)
      └── {chrom}_footprints.tsv.gz

Columns (all files)
-------------------
chrom       : chromosome
read_name   : CCS read identifier
fp_start    : footprint start (0-based)
fp_end      : footprint end
fp_mid      : (fp_start + fp_end) / 2  — midpoint
fp_class    : class label matching the subfolder

At plot time, load only the class file(s) you need for the chosen sort mode:
  nuc_mid  : PPP or FIRE_nucleosome; min |fp_mid - tss|
  nuc_edge : PPP or FIRE_nucleosome; min(|fp_start - tss|, |fp_end - tss|)
  ppp_mid  : PPP; min |fp_mid - tss|
  ppp_edge : PPP; min(|fp_start - tss|, |fp_end - tss|)

Usage
-----
  python3 precompute_fp_summaries.py                       # all chromosomes
  python3 precompute_fp_summaries.py --chroms chr1 chr22   # specific chroms
  python3 precompute_fp_summaries.py --out-dir /custom/dir
"""

import argparse
import gzip
import os

import pysam


# ── Input paths ───────────────────────────────────────────────────────────────

PEAKS_DIR = (
    "/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged/combined/"
    "joint_trained_peaks/"
)

FP_FILES = {
    "10-30":   PEAKS_DIR + "combined_all_chrs_10-30bp_fps.bed.gz",
    "30-45":   PEAKS_DIR + "combined_all_chrs_30-45bp_fps.bed.gz",
    "45-60":   PEAKS_DIR + "combined_all_chrs_45-60bp_fps.bed.gz",
    "60-80":   PEAKS_DIR + "combined_all_chrs_60-80bp_fps.bed.gz",
    "140-160": PEAKS_DIR + "combined_all_chrs_140-160bp_fps.bed.gz",
}

FIRE_NUC_BED_DIR = (
    "/project/spott/cshan/fiber-seq/results/PolII/FIRE_nuc_by_chr_combined_sample/"
    "parsed_bed_files/sorted_indexed"
)

DEFAULT_OUT_DIR = (
    "/project/spott/cshan/fiber-seq/results/PolII/footprint_summaries"
)

ALL_CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

# fp_class label → output subfolder name
CLASS_SUBDIRS = {
    "PPP":                  "PPP",
    "PIC":                  "PIC",
    "nucleosome (140-160)": "nucleosome_140-160",
    "FIRE_nucleosome":      "FIRE_nucleosome",
    "unknown":              "unknown",
}

HEADER = "chrom\tread_name\tfp_start\tfp_end\tfp_mid\tfp_class\n"


# ── Classification (mirrors classify_by_size in plot_polii_footprints.py) ─────

def classify_by_size(size: int) -> str:
    if 25 <= size <= 55:
        return "PPP"
    if 60 < size <= 80:
        return "PIC"
    if 140 <= size <= 160:
        return "nucleosome (140-160)"
    return "unknown"


# ── Per-chromosome streaming writer ──────────────────────────────────────────

def collect_chrom(chrom: str, out_dir: str) -> dict:
    """Stream footprints for *chrom* into per-class gzipped TSV files.

    Opens one file handle per fp_class, streams all sources in a single pass,
    and closes handles on exit.  Returns a {class: count} dict.
    """
    # Ensure all class subdirectories exist
    for subdir in CLASS_SUBDIRS.values():
        os.makedirs(os.path.join(out_dir, subdir), exist_ok=True)

    # Open one gzip handle per class
    handles = {
        cls: gzip.open(
            os.path.join(out_dir, subdir, f"{chrom}_footprints.tsv.gz"), "wt"
        )
        for cls, subdir in CLASS_SUBDIRS.items()
    }
    for fh in handles.values():
        fh.write(HEADER)

    counts = {cls: 0 for cls in CLASS_SUBDIRS}

    try:
        # ── FiberHMM size-range footprint files ───────────────────────────
        for label, fp_file in FP_FILES.items():
            if not os.path.exists(fp_file):
                print(f"  [WARN] missing FP file: {fp_file}", flush=True)
                continue
            tbx = pysam.TabixFile(fp_file)
            n = 0
            try:
                for row in tbx.fetch(chrom):
                    f = row.split("\t")
                    fp_start, fp_end, read_name = int(f[1]), int(f[2]), f[3]
                    fp_mid = (fp_start + fp_end) / 2.0
                    cls = classify_by_size(fp_end - fp_start)
                    handles[cls].write(
                        f"{chrom}\t{read_name}\t{fp_start}\t{fp_end}"
                        f"\t{fp_mid}\t{cls}\n"
                    )
                    counts[cls] += 1
                    n += 1
            except (ValueError, KeyError):
                pass
            finally:
                tbx.close()
            print(f"  {label:8s}: {n:>10,} rows", flush=True)

        # ── FIRE nucleosome per-chromosome bed ────────────────────────────
        fire_nuc_path = os.path.join(
            FIRE_NUC_BED_DIR, f"all_samples_nuc_{chrom}_parsed_sorted.bed.gz"
        )
        if os.path.exists(fire_nuc_path):
            tbx = pysam.TabixFile(fire_nuc_path)
            n = 0
            cls = "FIRE_nucleosome"
            try:
                for row in tbx.fetch(chrom):
                    f = row.split("\t")
                    if len(f) < 4:
                        continue
                    nuc_start, nuc_end, read_name = int(f[1]), int(f[2]), f[3]
                    nuc_mid = (nuc_start + nuc_end) / 2.0
                    handles[cls].write(
                        f"{chrom}\t{read_name}\t{nuc_start}\t{nuc_end}"
                        f"\t{nuc_mid}\t{cls}\n"
                    )
                    counts[cls] += 1
                    n += 1
            except (ValueError, KeyError):
                pass
            finally:
                tbx.close()
            print(f"  {'FIRE_nuc':8s}: {n:>10,} rows", flush=True)
        else:
            print(f"  [WARN] missing FIRE nuc bed: {fire_nuc_path}", flush=True)

    finally:
        for fh in handles.values():
            fh.close()

    return counts


# ── Main ──────────────────────────────────────────────────────────────────────

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Precompute per-read footprint summaries (one file per class per chrom)."
    )
    p.add_argument(
        "--chroms",
        nargs="+",
        default=ALL_CHROMS,
        metavar="CHROM",
        help="Chromosomes to process (default: all autosomes + chrX chrY).",
    )
    p.add_argument(
        "--out-dir",
        default=DEFAULT_OUT_DIR,
        help=f"Root output directory (default: {DEFAULT_OUT_DIR}).",
    )
    args = p.parse_args(argv)

    print(f"Output directory : {args.out_dir}")
    print(f"Subfolders       : {', '.join(CLASS_SUBDIRS.values())}")
    print(f"Chromosomes      : {' '.join(args.chroms)}\n", flush=True)

    total_rows = 0
    for chrom in args.chroms:
        print(f"[{chrom}] collecting footprints ...", flush=True)
        counts = collect_chrom(chrom, args.out_dir)
        n = sum(counts.values())
        cls_summary = "  ".join(
            f"{CLASS_SUBDIRS[k]}: {v:,}" for k, v in sorted(counts.items()) if v
        )
        print(f"[{chrom}] total {n:,}  |  {cls_summary}\n", flush=True)
        total_rows += n

    print(f"Done. Total rows written: {total_rows:,}")


if __name__ == "__main__":
    main()
