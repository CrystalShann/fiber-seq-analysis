#!/usr/bin/env python3
"""Haplotype-selective status of every paired cCRE in every sample, from FIRE's own peak test.

1. cCREs = universe/cre_universe.bed entries that are CRE1 or CRE2 of a pooled pair.
2. Per sample, <fire_dir>/<sample>-fire-v0.1-hap-differences.bed.gz (missing/empty = unphased,
   no rows): FIRE peaks that passed its per-haplotype coverage filter, with a two-sided Fisher
   test of hap1 vs hap2 accessibility (p_value, BH p_adjust) and diff = hap1 - hap2.
3. bedtools intersect, >=1 bp, on peak_start/peak_end; one cCRE can overlap several peaks.
4. A peak is haplotype-selective at p_adjust <= 0.05, favouring hap1 if diff > 0, hap2 if diff < 0.
5. status: conflict (significant peaks in both directions), hap1 / hap2, or ns.
6. Representative peak = smallest p_adjust (ties: smallest peak_start); its values are copied.
7. One row per cCRE per sample; a cCRE overlapping no tested peak in a sample has no row (untested).

hap1/hap2 are that sample's HP tags: never compare them across samples, and treat them as
consistent only within a phase block.
"""

import argparse
from pathlib import Path

import pandas as pd

from hsca_functions import ROOT, cre_status, read_hap_differences, sample_manifest
from LCL_coaccess_cres import METATABLE

COLUMNS = ["sample_name", "donor", "CRE_ID", "status", "n_hap_peaks", "n_hsca_peaks", "p_value", "p_adjust",
           "diff", "hap1_frac_acc", "hap2_frac_acc", "coverage_H1", "coverage_H2"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--sample-metatable", default=METATABLE)
    ap.add_argument("--out-dir", type=Path)
    args = ap.parse_args()
    out_dir = args.out_dir or args.root / "hsca" / "outputs"

    # 1. cCREs tested in at least one pooled pair
    pairs = pd.read_csv(args.root / "coaccess" / "LCL_coaccess_pairs.tsv.gz", sep="\t", usecols=["CRE1", "CRE2"])
    cres = pd.read_csv(args.root / "universe" / "cre_universe.bed", sep="\t", header=None,
                       names=["chrom", "start", "end", "CRE_ID", "CRE_label"])
    cres = cres[cres.CRE_ID.isin(set(pairs.CRE1) | set(pairs.CRE2))]
    print(f"{len(cres)} paired cCREs", flush=True)

    # 2-7. per sample
    samples = sample_manifest(args.sample_metatable, args.root)
    tables = []
    for s in samples.itertuples(index=False):
        peaks = read_hap_differences(Path(s.fire_dir) / f"{s.sample_name}-fire-v0.1-hap-differences.bed.gz")
        status = cre_status(cres, peaks)
        print(f"{s.sample_name}: {len(peaks)} tested peaks -> {len(status)} cCREs, "
              f"{int(status.status.ne('ns').sum()) if len(status) else 0} haplotype-selective", flush=True)
        if len(status):
            tables.append(status.assign(sample_name=s.sample_name, donor=s.donor))
    result = pd.concat(tables, ignore_index=True)[COLUMNS]

    out_dir.mkdir(parents=True, exist_ok=True)
    destination = out_dir / "LCL_hsca_cre_status.tsv.gz"
    partial = destination.with_name(destination.name + ".partial")
    result.to_csv(partial, sep="\t", index=False, compression="gzip")
    partial.replace(destination)
    print(f"Completed: {len(result)} rows, {result.sample_name.nunique()} samples -> {destination}", flush=True)


if __name__ == "__main__":
    main()
