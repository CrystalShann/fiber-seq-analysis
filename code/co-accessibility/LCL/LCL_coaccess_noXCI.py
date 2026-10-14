#!/usr/bin/env python3
"""Rebuild the pooled LCL co-accessibility test with X-inactivation excluded.

Female chrX fibers are a mixture of closed Xi and open Xa molecules, which
inflates a pooled 2x2 odds ratio without any co-accessibility. chrX pairs
are therefore re-pooled from male samples only (one active X); autosomal
and chrY pairs keep all 31 samples

"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from LCL_coaccess_cres import CELLS, ROOT, bh, test_counts

SUFFIX = "_noXCI"
SEX_DEPTH_CUT = 0.9  # two X copies sit near 1, one copy near 0.5
SCOPE = "genome-wide, chrX male-only"
FEMALE, MALE = "female (XX)", "male (XY)"
KEYS = ["CRE_pair", "CRE1", "CRE2", "dist", "cre1_start", "cre1_end", "cre2_start", "cre2_end",
        "chrom", "CRE_pair_labels"]
PARAMS = ["dataset", "coverage_mode", "read_rule", "fire_overlap_fraction", "min_dist", "max_dist",
          "pseudocount"]
GENES = ["n_genes", "gene_ids", "gene_names"]
STATS = CELLS + ["pval", "fisher_estimate", "pval_raw", "n_shared_reads", "n_contributing_samples",
                 "n_samples", "OR", "or_haldane", "expected_co_access", "observed_expected_ratio",
                 "excess_fraction", "zero_access_cre1", "zero_access_cre2", "fdr", "fdr_raw"]


#################################################
# donor sex from sequencing coverage
#################################################

# For each sample, it compares the average number of shared fibers per pair on chrX with the average on autosomes. 
# A ratio near 1 means two X copies (female), and a ratio near 0.5 means one X copy (male). 
# It returns a table with one row per sample

def call_sex(counts):
    # select chrX, chrY, and everything else is autosomes
    cls = np.select([counts.chrom == "chrX", counts.chrom == "chrY"], ["chrX", "chrY"], "autosome")
    # calculate the total spanning fibers and the number of pairs
    depth = counts.assign(cls=cls).groupby(["sample_name", "cls"], observed=True).n_shared_reads.agg(["sum", "size"])
    # divides total fibers by number of pairs to get the average shared fibers per pair
    depth = (depth["sum"] / depth["size"]).unstack()
  
    sex = pd.DataFrame({"x_depth": depth.chrX / depth.autosome}).reset_index()
    
    # samples at or above 0.9 are female, and all others are male
    if sex.x_depth.between(0.8, 1.0, inclusive="neither").any():
        raise ValueError("A sample has chrX depth between 0.8 and 1.0; sex cannot be called from coverage")
    sex["sex"] = np.where(sex.x_depth >= SEX_DEPTH_CUT, FEMALE, MALE)
    sex["sample_name"] = sex.sample_name.astype(str)
    donors = counts[["sample_name", "donor"]].drop_duplicates().astype(str)
    sex = sex.merge(donors, on="sample_name", validate="one_to_one")
    return sex[["sample_name", "donor", "x_depth", "sex"]]

#################################################
# save completed tables
#################################################
def write_atomic(table, destination):
    temporary = destination.with_name(destination.name + ".partial")
    table.to_csv(temporary, sep="\t", index=False,
                 compression="gzip" if destination.suffix == ".gz" else None)
    temporary.replace(destination)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=ROOT)
    args = ap.parse_args()
    coaccess = args.root / "coaccess"
    out_dir = coaccess / "noXCI"

    pooled = pd.read_csv(coaccess / "LCL_coaccess_pairs.tsv.gz", sep="\t", keep_default_na=False,
                         float_precision="round_trip")
    if pooled.analysis_scope.ne("genome-wide").any() or pooled.pseudocount.nunique() != 1:
        raise ValueError("Expected the genome-wide pooled table with one pseudocount")
    pc = int(pooled.pseudocount.iloc[0])
    counts = pd.read_csv(coaccess / "LCL_pair_sample_counts.tsv.gz", sep="\t",
                         dtype={"CRE_pair": "category", "chrom": "category", "sample_name": "category",
                                "donor": "category"})
    if counts.sample_name.nunique() != 31:
        raise ValueError("Per-sample counts must cover all 31 samples")

    sex = call_sex(counts)
    female = set(sex.sample_name[sex.sex == FEMALE])
    print(f"{len(female)} female and {len(sex) - len(female)} male samples", flush=True)

    # drop the female chrX rows, then pool the remaining per-sample cells
    counts = counts[~((counts.chrom == "chrX") & counts.sample_name.isin(female))]
    counts = counts.assign(contributing=counts.n_shared_reads > 0)
    grouped = counts.groupby("CRE_pair", observed=True)
    rebuilt = grouped[CELLS + ["contributing"]].sum()
    rebuilt["n_samples"] = grouped.size()
    rebuilt = rebuilt[rebuilt[CELLS].sum(axis=1) > 0].rename(columns={"contributing": "n_contributing_samples"})
    rebuilt.index = rebuilt.index.astype(str)
    missing = rebuilt.index.difference(pooled.CRE_pair)
    if len(missing):
        raise ValueError(f"{len(missing)} re-pooled pairs are absent from the pooled table")

    # pooled order and annotations; only pairs that still have shared fibers are tested
    result = pooled[KEYS + PARAMS + GENES].merge(rebuilt, left_on="CRE_pair", right_index=True, how="inner")
    result = result.reset_index(drop=True)
    expected_n = np.where(result.chrom == "chrX", len(sex) - len(female), len(sex))
    if not np.array_equal(result.n_samples.to_numpy(), expected_n):
        raise ValueError("Each pair must be pooled over 15 male samples on chrX and 31 samples elsewhere")

    # statistics exactly as in LCL_coaccess_cres.main()
    cells = result[CELLS].to_numpy(dtype=np.int64)
    result = pd.concat([result, test_counts(cells, pc)], axis=1)
    result["n_shared_reads"] = cells.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        result["OR"] = ((result.co_access.astype(float) + pc) * (result.co_closed + pc)) / ((result.CRE1_access + pc) * (result.CRE2_access + pc))
        result["or_haldane"] = ((result.co_access + 0.5) * (result.co_closed + 0.5)) / ((result.CRE1_access + 0.5) * (result.CRE2_access + 0.5))
    a = (result.co_access + result.CRE1_access).astype(float)
    b = result.co_access + result.CRE2_access
    result["expected_co_access"] = a * b / result.n_shared_reads
    result["observed_expected_ratio"] = result.co_access / result.expected_co_access.replace(0, np.nan)
    result["excess_fraction"] = (result.co_access - result.expected_co_access) / result.n_shared_reads
    result["zero_access_cre1"] = a == 0
    result["zero_access_cre2"] = b == 0
    result["fdr"] = bh(result.pval)
    result["fdr_raw"] = bh(result.pval_raw)

    # outside chrX nothing was re-pooled, so everything except the BH family must match
    reference = pooled.set_index("CRE_pair").loc[result.CRE_pair]
    other = (result.chrom != "chrX").to_numpy()
    for column in CELLS + ["n_shared_reads", "n_contributing_samples", "pval", "pval_raw", "fisher_estimate"]:
        if not np.array_equal(result[column].to_numpy()[other], reference[column].to_numpy()[other]):
            raise ValueError(f"Non-chrX {column} differs from the pooled table")
    if (result.n_shared_reads.to_numpy()[~other] > reference.n_shared_reads.to_numpy()[~other]).any():
        raise ValueError("A chrX pair gained fibers after dropping female samples")

    result["chrX_rule"] = np.where(result.chrom == "chrX", "male samples only", "all samples")
    result["analysis_scope"] = SCOPE
    result = result.rename(columns={c: c + SUFFIX for c in STATS})
    result = result[KEYS + ["chrX_rule"] + [c + SUFFIX for c in STATS] + ["dataset", "analysis_scope"]
                    + PARAMS[1:] + GENES]

    out_dir.mkdir(parents=True, exist_ok=True)
    write_atomic(sex, out_dir / "LCL_sample_sex_noXCI.tsv")
    write_atomic(result, out_dir / "LCL_coaccess_pairs_noXCI.tsv.gz")
    positive = (result.fdr_noXCI < .05) & (result.OR_noXCI > 1) & \
        (result.co_access_noXCI > result.expected_co_access_noXCI)
    print(f"Completed: {len(result)} pairs ({int((result.chrom == 'chrX').sum())} chrX); "
          f"{int(positive.sum())} positive at FDR <0.05, "
          f"{int((positive & (result.chrom == 'chrX')).sum())} on chrX -> {out_dir}", flush=True)


if __name__ == "__main__":
    main()
