#!/usr/bin/env python3
"""Re-pool existing co-accessibility counts for the 27-sample HSCA cohort."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from hsca_functions import CELLS, PHASED_METATABLE, ROOT, phased_samples
from LCL_coaccess_cres import bh_min_shared
from LCL_coaccess_noXCI import FEMALE, GENES, KEYS, PARAMS, STATS, call_sex, write_atomic
from count_and_test_pairs import summarize_coaccess_counts, test_counts


# inputs: 
#   
# 1. The 31-sample pooled pair table. Only its pair coordinates, gene annotations and settings (for example the pseudocount) 
#   are used (no test stats)
# 2. The per-sample count table. For every cCRE pair and every sample it has four numbers: fibers with both cCREs closed, 
#   only the first open, only the second open, and both open.

#################################################
# re-pool one cohort and recompute Fisher tests and BH
#################################################
# Sum the four cells only after selecting the samples (and the chrX sex rule).
# Reuse upstream coordinates/annotations, but no upstream statistics.
# Keep the upstream pseudocount and >=20 shared-read BH threshold.
def pool_counts(counts, annotations, pseudocount):
    counts = counts.assign(contributing=counts.n_shared_reads > 0)
    grouped = counts.groupby("CRE_pair", observed=True)
    summed = grouped[CELLS + ["contributing"]].sum()
    summed["n_samples"] = grouped.size()
    summed = summed[summed[CELLS].sum(axis=1) > 0].rename(columns={"contributing": "n_contributing_samples"})
    result = annotations.merge(summed, left_on="CRE_pair", right_index=True, validate="one_to_one").reset_index(drop=True)
    result = pd.concat([result, test_counts(result[CELLS].to_numpy(dtype=np.int64), pseudocount)], axis=1)
    return bh_min_shared(summarize_coaccess_counts(result, pseudocount))


#################################################
# read upstream per-sample counts and write HSCA-only pooled tables
#################################################
# The upstream 31-sample files stay intact. Both outputs use only the phased sheet;
# noXCI additionally excludes female chrX rows before pooling and BH correction.
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--sample-metatable", default=PHASED_METATABLE)
    ap.add_argument("--out-dir", type=Path)
    args = ap.parse_args()
    samples = phased_samples(args.sample_metatable, args.root)
    out_dir = args.out_dir or args.root.parent / "hsca" / "outputs"
    source = pd.read_csv(args.root / "coaccess/LCL_coaccess_pairs.tsv.gz", sep="\t")
    if source.pseudocount.nunique() != 1 or source.analysis_scope.ne("genome-wide").any():
        raise ValueError("Expected genome-wide upstream counts with one pseudocount")
    pc = int(source.pseudocount.iloc[0])
    annotations = source[KEYS + PARAMS + GENES]
    counts = pd.read_csv(args.root / "coaccess/LCL_pair_sample_counts.tsv.gz", sep="\t",
                         dtype={"CRE_pair": "category", "sample_name": "category", "chrom": "category", "donor": "category"})
    counts = counts[counts.sample_name.isin(samples.sample_name)].copy()
    if set(counts.sample_name) != set(samples.sample_name):
        raise ValueError("Upstream counts must include every phased sample")
    if counts.duplicated(["CRE_pair", "sample_name"]).any():
        raise ValueError("Duplicate pair/sample rows in upstream counts")
    sex = call_sex(counts)
    female = set(sex.loc[sex.sex == FEMALE, "sample_name"])
    pooled = pool_counts(counts, annotations, pc)
    if not pooled.n_samples.eq(27).all():
        raise ValueError("Each pooled pair must include all 27 sample rows")
    pooled["analysis_scope"] = "genome-wide, 27 phased samples"
    print(f"Pooled {len(pooled)} pairs from 27 phased samples", flush=True)
    noxci = pool_counts(counts[~((counts.chrom == "chrX") & counts.sample_name.isin(female))], annotations, pc)
    expected = np.where(noxci.chrom == "chrX", 27 - len(female), 27)
    if not np.array_equal(noxci.n_samples, expected):
        raise ValueError("Missing sample rows in noXCI pooling")
    noxci["analysis_scope"] = "genome-wide, 27 phased samples, chrX male-only"
    noxci["chrX_rule"] = np.where(noxci.chrom == "chrX", "male phased samples only", "all phased samples")
    noxci = noxci.drop(columns="co_inaccess").rename(columns={c: c + "_noXCI" for c in STATS})
    out_dir.mkdir(parents=True, exist_ok=True)
    write_atomic(pooled, out_dir / "LCL_hsca_coaccess_pairs.tsv.gz")
    write_atomic(noxci, out_dir / "LCL_hsca_coaccess_pairs_noXCI.tsv.gz")
    write_atomic(sex, out_dir / "LCL_hsca_sample_sex.tsv")
    write_atomic(samples, out_dir / "LCL_hsca_sample_manifest.tsv")
    print(f"Saved pooled HSCA comparisons to {out_dir}", flush=True)


if __name__ == "__main__":
    main()
