#!/usr/bin/env python3
"""REF/ALT accessibility skew of LCL cCREs and pairs in the 27 phased samples.

For a cCRE, open means accessible; for a pair, open means both cCREs are accessible.
Read bases at the selected heterozygous SNP define REF and ALT groups.
Within donors and before pooling, require >= MIN_READS reads for EACH allele per unit.
Run two-sided Fisher tests with BH separately per donor or across pooled units.
Write four result tables: cCRE/pair allele tests within donors and pooled across donors.
"""

import argparse
import multiprocessing
import os
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from hsca_functions import (AUTOSOMES, CELLS, MIN_READS, ROOT, allele_counts,
                            call_status, choose_cre_variants, choose_pair_variants, count_pairs, fisher_bh,
                            het_snvs, incidence_for_chrom, phased_samples,
                            read_alleles)
from hsca_functions import PHASED_METATABLE

SHARED = {}  # pairs and cCREs, inherited by the forked workers
ALLELE_COUNTS = ["ref_open", "ref_closed", "alt_open", "alt_closed"]

#################################################
# per-sample workers
#################################################

def sample_het_snvs(sample):
    """Pass 1: this sample's het SNVs inside paired cCREs."""
    return het_snvs(sample.phased_vcf_file, sample.cell_line, SHARED["cres"]).assign(sample_name=sample.sample_name)


def process_sample(sample, reference, variants, pair_variants):
    """Pass 2: REF/ALT counts of one sample, chromosome by chromosome.

    reference = this sample's rows of coaccess/LCL_pair_sample_counts.tsv.gz; the pooled 2x2
    cells recomputed here must equal them, or the read rule has drifted.
    """
    pairs, cres = SHARED["pairs"], SHARED["cres"]
    cre_allele, pair_allele, empty = [], [], set()
    with tempfile.TemporaryDirectory(prefix="hsca-", dir=os.environ.get("SLURM_TMPDIR")) as tmp:
        tmp = Path(tmp)
        for chrom, chrom_pairs in pairs.groupby("chrom", sort=False):
            # the same cCRE BED5 and read rule as LCL_coaccess_cres.main()
            cre_bed = tmp / "cre.bed"
            cres[cres.CRE_ID.isin(set(chrom_pairs.CRE1) | set(chrom_pairs.CRE2))].to_csv(
                cre_bed, sep="\t", index=False, header=False)
            cov, acc, _ = incidence_for_chrom(chrom, cre_bed, sample.blocks_path, sample.elements_path, 0.5, tmp)
            cells = count_pairs(chrom_pairs, cov, acc)
            if not np.array_equal(cells, reference.loc[chrom_pairs.CRE_pair, CELLS].to_numpy()):
                raise ValueError(f"{sample.sample_name} {chrom}: counts differ from LCL_pair_sample_counts")

            # REF vs ALT at the chosen SNVs this sample is het for
            v = variants[variants.chrom == chrom]
            calls = read_alleles(sample.cram_path, chrom, v)
            cre_allele.append(allele_counts(
                ((c, var, cov.get(c, empty), acc.get(c, empty)) for c, var in zip(v.CRE_ID, v.variant_id)), calls))
            pv = pair_variants[pair_variants.chrom == chrom]
            pair_allele.append(allele_counts(
                ((p, var, cov.get(c1, empty) & cov.get(c2, empty), acc.get(c1, empty) & acc.get(c2, empty))
                 for p, var, c1, c2 in zip(pv.CRE_pair, pv.variant_id, pv.CRE1, pv.CRE2)), calls))
            print(f"{sample.sample_name} {chrom}: {len(pv)} allele pairs, {len(v)} SNVs", flush=True)
            del cov, acc
    tag = lambda frames: pd.concat(frames, ignore_index=True).assign(sample_name=sample.sample_name, donor=sample.donor)
    return tag(cre_allele), tag(pair_allele)

#################################################
# tests
#################################################

def within_sample(counts, rows, cols, positive, negative):
    """>= MIN_READS per row of the 2x2, two-sided Fisher, BH per sample, effect = row1 - row2."""
    (a1, b1), (a2, b2) = rows
    counts = counts[(counts[a1] + counts[b1] >= MIN_READS) & (counts[a2] + counts[b2] >= MIN_READS)]
    counts = fisher_bh(counts, [a1, b1, a2, b2], group="sample_name")
    counts[cols[0]] = counts[a1] / (counts[a1] + counts[b1])
    counts[cols[1]] = counts[a2] / (counts[a2] + counts[b2])
    counts["diff"] = counts[cols[0]] - counts[cols[1]]
    counts["status"] = call_status(counts["diff"], counts.p_adjust, positive, negative)
    return counts


def pooled_alleles(by_sample, variants, key):
    """Pool only donor/unit rows with >= MIN_READS reads of EACH allele.

    Apply the same coverage rule as within_sample before summing counts.
    n_samples_both_alleles counts donors passing this coverage rule.
    """
    ref_n = by_sample.ref_open + by_sample.ref_closed
    alt_n = by_sample.alt_open + by_sample.alt_closed
    both = by_sample[(ref_n >= MIN_READS) & (alt_n >= MIN_READS)]
    alt_minus_ref = both.alt_open / (both.alt_open + both.alt_closed) - both.ref_open / (both.ref_open + both.ref_closed)
    consistency = both.assign(alt_higher=alt_minus_ref > 0, ref_higher=alt_minus_ref < 0).groupby("unit").agg(
        n_samples_alt_higher=("alt_higher", "sum"), n_samples_ref_higher=("ref_higher", "sum"),
        n_samples_both_alleles=("unit", "size"))
    pooled = both.groupby("unit")[ALLELE_COUNTS].sum().join(consistency).reset_index().rename(columns={"unit": key})
    pooled = variants.merge(pooled, on=key, how="left", validate="one_to_one")
    filled = ALLELE_COUNTS + ["n_samples_alt_higher", "n_samples_ref_higher", "n_samples_both_alleles"]
    pooled[filled] = pooled[filled].fillna(0).astype(np.int64)
    ref_n, alt_n = pooled.ref_open + pooled.ref_closed, pooled.alt_open + pooled.alt_closed
    pooled = fisher_bh(pooled, ["alt_open", "alt_closed", "ref_open", "ref_closed"],
                       testable=(ref_n >= MIN_READS) & (alt_n >= MIN_READS)).rename(columns={"p_adjust": "fdr"})
    with np.errstate(divide="ignore", invalid="ignore"):
        pooled["ref_frac_open"] = pooled.ref_open / ref_n
        pooled["alt_frac_open"] = pooled.alt_open / alt_n
    pooled["diff"] = pooled.alt_frac_open - pooled.ref_frac_open
    return pooled


def write(table, destination):
    partial = destination.with_name(destination.name + ".partial")
    table.to_csv(partial, sep="\t", index=False, compression="gzip")
    partial.replace(destination)
    print(f"{len(table)} rows -> {destination}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--sample-metatable", default=PHASED_METATABLE)
    ap.add_argument("--out-dir", type=Path)
    ap.add_argument("--chrom", choices=AUTOSOMES, help="Smoke test: this chromosome only (BH within it)")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    args = ap.parse_args()
    out_dir = args.out_dir or args.root.parent / "hsca" / "outputs"
    chroms = [args.chrom] if args.chrom else AUTOSOMES

    samples = phased_samples(args.sample_metatable, args.root)
    pairs = pd.read_csv(args.root / "coaccess" / "LCL_coaccess_pairs.tsv.gz", sep="\t",
                        usecols=["CRE_pair", "CRE1", "CRE2", "chrom", "cre1_start", "cre1_end", "cre2_start", "cre2_end"])
    pairs = pairs[pairs.chrom.isin(chroms)].reset_index(drop=True)
    cres = pd.read_csv(args.root / "universe" / "cre_universe.bed", sep="\t", header=None,
                       names=["chrom", "start", "end", "CRE_ID", "CRE_label"])
    cres = cres[cres.CRE_ID.isin(set(pairs.CRE1) | set(pairs.CRE2))].reset_index(drop=True)
    SHARED.update(pairs=pairs, cres=cres)
    print(f"{len(samples)} phased samples, {len(pairs)} pairs, {len(cres)} cCREs on {len(chroms)} autosomes", flush=True)

    reference = pd.read_csv(args.root / "coaccess" / "LCL_pair_sample_counts.tsv.gz", sep="\t",
                            usecols=["CRE_pair", "sample_name"] + CELLS, dtype={"sample_name": "category"})
    reference = reference[reference.sample_name.isin(samples.sample_name) & reference.CRE_pair.isin(pairs.CRE_pair)]
    reference = {s: g.set_index("CRE_pair") for s, g in reference.groupby("sample_name", observed=True)}

    rows = [row for _, row in samples.iterrows()]  # Series pickle to the workers; itertuples rows do not
    with ProcessPoolExecutor(args.workers, mp_context=multiprocessing.get_context("fork")) as pool:
        # pass 1: het SNVs and the variant used for each cCRE and pair
        het = pd.concat(pool.map(sample_het_snvs, rows), ignore_index=True)
        cre_variants = choose_cre_variants(het, cres)
        pair_variants = choose_pair_variants(pairs, cre_variants)
        print(f"{het.variant_id.nunique()} het SNVs; {len(cre_variants)} cCREs and {len(pair_variants)} pairs "
              f"with a variant", flush=True)
        is_het = lambda table, s: table[table.het_samples.str.split(";").apply(lambda x: s in x)]

        # pass 2: per-sample counts
        jobs = [pool.submit(process_sample, s, reference[s.sample_name], is_het(cre_variants, s.sample_name),
                            is_het(pair_variants, s.sample_name)) for s in rows]
        cre_allele, pair_allele = (pd.concat(parts, ignore_index=True) for parts in zip(*(j.result() for j in jobs)))

    pair_keys = pairs[["CRE_pair", "CRE1", "CRE2", "chrom"]]
    first = ["sample_name", "donor"]

    # within sample, ALT vs REF
    allele_rows = (("alt_open", "alt_closed"), ("ref_open", "ref_closed"))
    allele_out = ALLELE_COUNTS + ["ref_frac_open", "alt_frac_open", "diff", "p_value", "p_adjust", "status"]
    cre_by_sample = within_sample(cre_allele, allele_rows, ("alt_frac_open", "ref_frac_open"), "alt", "ref")
    cre_by_sample = cre_by_sample.rename(columns={"unit": "CRE_ID"})[first + ["CRE_ID", "variant_id"] + allele_out]
    pair_by_sample = within_sample(pair_allele, allele_rows, ("alt_frac_open", "ref_frac_open"), "alt", "ref")
    pair_by_sample = pair_by_sample.rename(columns={"unit": "CRE_pair"}).merge(
        pair_variants[["CRE_pair", "variant_cre"]], on="CRE_pair").merge(pair_keys, on="CRE_pair")
    pair_by_sample = pair_by_sample[first + ["CRE_pair", "CRE1", "CRE2", "variant_cre", "variant_id"] + allele_out]

    # across samples, ALT vs REF
    pooled_out = ["variant_id", "n_het_samples"] + ALLELE_COUNTS + ["ref_frac_open", "alt_frac_open", "diff",
                  "p_value", "fdr", "n_samples_alt_higher", "n_samples_ref_higher", "n_samples_both_alleles"]
    cre_pooled = pooled_alleles(cre_allele, cre_variants[["CRE_ID", "variant_id", "n_het_samples"]], "CRE_ID")
    cre_pooled = cre_pooled[["CRE_ID"] + pooled_out]
    pair_pooled = pooled_alleles(pair_allele, pair_variants[["CRE_pair", "CRE1", "CRE2", "variant_cre",
                                                             "variant_id", "n_het_samples"]], "CRE_pair")
    pair_pooled = pair_pooled[["CRE_pair", "CRE1", "CRE2", "variant_cre"] + pooled_out]

    out_dir.mkdir(parents=True, exist_ok=True)
    write(cre_by_sample, out_dir / "LCL_hsca_cre_allele_by_sample.tsv.gz")
    write(pair_by_sample, out_dir / "LCL_hsca_pair_allele_by_sample.tsv.gz")
    write(cre_pooled, out_dir / "LCL_hsca_cre_allele.tsv.gz")
    write(pair_pooled, out_dir / "LCL_hsca_pair_allele.tsv.gz")
    print(f"Completed: {int((cre_pooled.fdr < 0.05).sum())} cCREs and {int((pair_pooled.fdr < 0.05).sum())} pairs allele-skewed "
          f"across samples", flush=True)


if __name__ == "__main__":
    main()
