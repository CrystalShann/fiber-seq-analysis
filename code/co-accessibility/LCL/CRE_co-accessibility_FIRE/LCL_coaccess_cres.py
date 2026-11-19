#!/usr/bin/env python3
"""Genome-wide pooled co-accessibility across all 31 LCL samples.

Adapted from macrophage/03_coaccess_cres.py. Candidate cCREs have >=50% of
their length inside ONE merged LCL FIRE peak and share a canonical TSS +/-10 kb window. The interval gap
must satisfy 500 < gap < 20000. Shared reads have aligned bases covering the
FULL length of both cCREs: every base of each cCRE must
fall inside the read's M/= /X blocks, so a read whose alignment gap or end lands
inside either cCRE is not counted for that pair.
A cCRE is accessible when ONE FIRE element on that same read covers >=50% of
the cCRE length after intersecting with CIGAR M/= /X blocks. D/N gaps never
contribute coverage or accessibility; separate FIRE elements are not combined.

Sum raw 2x2 cells across samples before testing. Preserve macrophage two-sided
Fisher tests on both table+1 and raw counts; apply BH once across unique pooled
pairs with >= 20 fibers spanning both cCREs (others keep p-values, FDR is NA).
Donor information is metadata only; there are no LCL timepoints.
"""

import argparse
import gzip
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "parsing_functions"))
from count_and_test_pairs import bh, count_pairs, test_counts, summarize_coaccess_counts
import fiber_regions as fr

ROOT = Path(os.environ.get("LCL_COACCESS_ROOT", "/project/spott/cshan/fiber-seq/LCL_project/co-accessibility/CRE_co-accessibility_FIRE"))
METATABLE = os.environ.get("LCL_SAMPLE_METATABLE", "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv")
CHROMS = fr.CHROMS
CELLS = ["co_closed", "CRE1_access", "CRE2_access", "co_access"]
COVERAGE_MODE = "aligned_blocks_v1"
# coverage filter: only pairs with >= 20 fibers spanning both cCREs enter the BH family
MIN_SHARED_READS = 20


def bh_min_shared(result):
    # pairs below MIN_SHARED_READS keep their p-values, but fdr / fdr_raw are NA
    tested = result.n_shared_reads >= MIN_SHARED_READS
    for p, q in (("pval", "fdr"), ("pval_raw", "fdr_raw")):
        result[q] = np.nan
        result.loc[tested, q] = bh(result.loc[tested, p])
    result["min_shared_reads"] = MIN_SHARED_READS
    return result

#################################################
# build sample meta table
#################################################
# read the sample metatable and FIRE element and FIRE peaks file path for all 31 samples
def sample_manifest(path, root):
    samples = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"sample_name", "cell_line", "fire_dir"}
    if not required.issubset(samples.columns):
        raise ValueError(f"Sample CSV requires {sorted(required)}")
    if len(samples) != 31 or samples.sample_name.nunique() != 31:
        raise ValueError("LCL analysis requires exactly 31 unique samples; no samples are silently skipped")
    if samples[list(required)].eq("").any().any():
        raise ValueError("Empty sample_name, cell_line, or fire_dir in sample CSV")
    for name in samples.sample_name:
        if any(c in name for c in "/\\\t\n\r"):
            raise ValueError(f"Invalid sample identifier: {name!r}")
    return fr.add_lcl_fire_paths(samples, root)

#################################################
# check every file exists in sample meta table
#################################################
def require_files(paths):
    missing = [str(p) for p in paths if not Path(p).is_file() or Path(p).stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing/empty inputs:\n" + "\n".join(missing))

#################################################
# builds the candidate CRE pairs test
#################################################
# a pair must share a gene's TSS window and should be 500-20,000 bp distance

# /project/spott/cshan/fiber-seq/LCL_project/co-accessibility/CRE_co-accessibility_FIRE/universe/cre_universe.bed -> 
  # this the the overlap ENCODE cCREs /project/spott/cshan/annotations/GRCh38-cCREs.bedand merged FIRE peaks across all 31 samples 
  # /project/spott/cshan/fiber-seq/LCL_project/co-accessibility/CRE_co-accessibility_FIRE/universe/fire_peaks_union.bed
  
  
# TSS: /project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed


def build_pairs(cre, membership, min_dist, max_dist):
    """Gene memberships plus one stable coordinate-oriented record per pair."""
    # look up key by cCRE ID
    pos = cre.set_index("CRE_ID")
    records = []
    # loop over genes, each group hold rows for the cCREs that fall in that gene's TSS window, 10kb on each side
    for gene_id, group in membership.groupby("gene_id", sort=False):
        # sort the gene's CRE ID by chr, start, end
        ids = sorted(group.CRE_ID.unique(), key=lambda c: (pos.at[c, "chrom"], pos.at[c, "start"], pos.at[c, "end"], c))
        # select the first row and get every pair once
        gene = group.iloc[0]
        for i, c1 in enumerate(ids):
            for c2 in ids[i + 1:]:
                if pos.at[c1, "chrom"] != pos.at[c2, "chrom"]:
                    continue
                  # distance from the end the c1 (first CRE) to the start of c2, overlapping CREs get 0
                gap = max(0, int(pos.at[c2, "start"]) - int(pos.at[c1, "end"]))
                if gap >= max_dist:
                    break
                if gap > min_dist:
                    records.append((c1, c2, gap, gene_id, gene.gene_name, gene.transcript_type))
    # pairwise cCRE id table, and add CRE coordinates
    gm = pd.DataFrame(records, columns=["CRE1", "CRE2", "dist", "gene_id", "gene_name", "transcript_type"])
    gm["CRE_pair"] = gm.CRE1 + "-" + gm.CRE2
    pairs = gm.drop_duplicates("CRE_pair")[["CRE_pair", "CRE1", "CRE2", "dist"]].reset_index(drop=True)
    if not pairs.empty:
        for side in (1, 2):
            entries = pos.loc[pairs[f"CRE{side}"]]
            pairs[f"cre{side}_start"] = entries.start.to_numpy()
            pairs[f"cre{side}_end"] = entries.end.to_numpy()
        pairs["chrom"] = pos.loc[pairs.CRE1, "chrom"].to_numpy()
        pairs["CRE_pair_labels"] = pos.loc[pairs.CRE1, "CRE_label"].to_numpy() + "-" + pos.loc[pairs.CRE2, "CRE_label"].to_numpy()
    return pairs, gm

#################################################
# spanning reads and accessible reads for every CRE
#################################################
# fr.region_reads (parsing_functions/fiber_regions.py). With rule="aligned" the read's
# aligned blocks must cover every base of the cCRE, and ONE FIRE element must cover
# >= fraction of the cCRE after clipping to aligned blocks.

# Example. The cCRE is 200 bp, so the threshold is 100 bp. One FIRE element covers 160 bp of it,
# but a deletion removes 70 bp of that stretch, leaving 90 aligned bp. That is below 100,
# so the read is not accessible. If the aligned part were 110 bp, it would be

#################################################
# build 2x2 table for every CRE pair, for one sample and on one chr
#################################################

# table1: one pair, one sample
# co_closed	CRE1_access	CRE2_access	co_access
# 60	20	10	10
# 
# table2: 2x2 table
# 	CRE2 closed	CRE2 accessible
# CRE1 closed	60	10
# CRE1 accessible	20	10
# 
# table3: pooling across samples
# sample	co_closed	CRE1_access	CRE2_access	co_access	shared
# sample 1	60	20	10	10	100
# sample 2	45	15	8	12	80
# sample 3	30	10	5	5	50
# pooled	135	45	23	27	230
# 
# Expected co-access if independent = (27 + 45) × (27 + 23) / 230 = 72 × 50 / 230 ≈ 15.7
# Observed co-access = 27, so observed/expected ≈ 1.7. The two cCREs are accessible together more often than chance.
# OR = (27 + 1)(135 + 1) / ((45 + 1)(23 + 1)) = 28 × 136 / (46 × 24) ≈ 3.45

# Shared implementation: parsing_functions/count_and_test_pairs.py

#################################################
# BH correction
#################################################
# Shared implementation: parsing_functions/count_and_test_pairs.py

#################################################
# Fisher test on 2x2 table
#################################################
# Shared implementation: parsing_functions/count_and_test_pairs.py

#################################################
# Parse options and validates inputs
#################################################

# --fire-overlap-fraction: Minimum fraction of a cCRE's length that must be covered by ONE FIRE,
# element on the same read, counting only aligned bases, for the cCRE to be called accessible on that read

# --min-dist and --max-dist default to 500 and 20000: distance between CRE pairs

# --pseudocount defaults to 1

def main():
    # build command line parser
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--sample-metatable", default=METATABLE)
    ap.add_argument("--out-dir", type=Path)
    ap.add_argument("--prepare-manifest", action="store_true", help="Validate all31 source inputs and write universe/sample_manifest.tsv only")
    ap.add_argument("--fire-overlap-fraction", type=float, default=0.5, help="Fraction of cCRE covered by ONE same-read FIRE element")
    ap.add_argument("--read-rule", choices=fr.READ_RULES, default="aligned",
                    help="aligned (default): aligned blocks cover every cCRE base; span: read span contains the cCRE")
    ap.add_argument("--min-dist", type=int, default=500)
    ap.add_argument("--max-dist", type=int, default=20000)
    ap.add_argument("--pseudocount", type=int, default=1)
    ap.add_argument("--chrom", choices=CHROMS, help="Smoke test only: FDR scope becomes this chromosome")
    args = ap.parse_args()
    if not (0 < args.fire_overlap_fraction <= 1) or not (0 <= args.min_dist < args.max_dist) or args.pseudocount < 0:
        ap.error("Require 0 < FIRE fraction <= 1, 0 <= min-dist < max-dist, and pseudocount >= 0")
    samples = sample_manifest(args.sample_metatable, args.root)
    uni = args.root / "universe"
    manifest_path = uni / "sample_manifest.tsv"
    if args.prepare_manifest:
        uni.mkdir(parents=True, exist_ok=True)
        # creates folder
        (uni / "universe.complete").unlink(missing_ok=True)
        source_files = [p for col in ("cram_path", "peaks_path", "elements_path") for p in samples[col]]
        source_files += [p + ".crai" for p in samples.cram_path]
        source_files += [p + ".tbi" for col in ("peaks_path", "elements_path") for p in samples[col]]
        require_files(source_files)
        uni.mkdir(parents=True, exist_ok=True)
        samples.to_csv(manifest_path, sep="\t", index=False)
        samples[["sample_name", "cram_path", "peaks_path", "elements_path", "spans_path"]].to_csv(
            uni / "sample_inputs.tsv", sep="\t", index=False, header=False)
        print(f"Validated all {len(samples)} samples -> {manifest_path}")
        return
    require_files([manifest_path, uni / "universe.complete", uni / "cre_universe.bed", uni / "cre_gene_map.tsv.gz"])
    saved_manifest = pd.read_csv(manifest_path, sep="\t", dtype=str, keep_default_na=False)
    compare = ["sample_name", "fire_dir", "donor", "spans_path", "elements_path"]
    if not saved_manifest[compare].equals(samples[compare]):
        raise ValueError("Manifest differs from current sample CSV/root; rebuild LCL_fire_universe.sh")
    reads_col = "blocks_path" if args.read_rule == "aligned" else "spans_path"
    coverage_mode = COVERAGE_MODE if args.read_rule == "aligned" else "read_spans"
    require_files([p for col in (reads_col, "elements_path") for p in samples[col]] +
                  [p + ".tbi" for col in (reads_col, "elements_path") for p in samples[col]])
    for sample in (samples.itertuples(index=False) if args.read_rule == "aligned" else ()):
        marker = Path(sample.blocks_path + ".source.tsv")
        if not marker.exists() or marker.read_text().splitlines()[:2] != [sample.cram_path, COVERAGE_MODE]:
            raise ValueError(f"Rebuild aligned blocks for {sample.sample_name}: missing/current provenance required")
        if Path(sample.blocks_path).stat().st_mtime < Path(sample.cram_path).stat().st_mtime:
            raise ValueError(f"Aligned blocks predate CRAM for {sample.sample_name}")
    cre = pd.read_csv(uni / "cre_universe.bed", sep="\t", header=None,
                     names=["chrom", "start", "end", "CRE_ID", "CRE_label"])
    memberships = pd.read_csv(uni / "cre_gene_map.tsv.gz", sep="\t", keep_default_na=False)
    if not cre.CRE_ID.is_unique or (cre.end <= cre.start).any():
        raise ValueError("cCRE IDs must be unique and intervals nonempty")
    chroms = [args.chrom] if args.chrom else CHROMS
    cre = cre[cre.chrom.isin(chroms)].reset_index(drop=True)
    memberships = memberships[memberships.CRE_ID.isin(cre.CRE_ID)]
    pairs, gene_map = build_pairs(cre, memberships, args.min_dist, args.max_dist)
    if pairs.empty:
        raise ValueError("No candidate pairs after gene-window and distance filters")
    out_dir = args.out_dir or (args.root / "coaccess" if not args.chrom else args.root / "smoke" / args.chrom)
    out_dir.mkdir(parents=True, exist_ok=True)
    pooled = np.zeros((len(pairs), 4), dtype=np.int64)
    contributing_samples = np.zeros(len(pairs), dtype=np.int64)
    qc = []
    sample_counts_tmp = out_dir / "LCL_pair_sample_counts.tsv.gz.partial"
    print(f"All31 pooled: {len(pairs)} distinct pairs on {pairs.chrom.nunique()} chromosomes", flush=True)
    with tempfile.TemporaryDirectory(prefix="lcl-coaccess-", dir=os.environ.get("SLURM_TMPDIR")) as temp, \
            gzip.open(sample_counts_tmp, "wt") as sample_output:
        temp = Path(temp)
        first_output = True
        for chrom in chroms:
            idx = pairs.index[pairs.chrom == chrom].to_numpy()
            if not len(idx):
                continue
            chrom_pairs = pairs.loc[idx]
            cre_bed = temp / "cre.bed"
            involved = set(chrom_pairs.CRE1) | set(chrom_pairs.CRE2)
            cre[cre.CRE_ID.isin(involved)].to_csv(cre_bed, sep="\t", index=False, header=False)
            for sample in samples.itertuples(index=False):
                cov, acc, orphans = fr.region_reads(chrom, cre_bed, getattr(sample, reads_col), sample.elements_path,
                                                    temp, args.read_rule, args.fire_overlap_fraction)
                cells = count_pairs(chrom_pairs, cov, acc)
                pooled[idx] += cells
                shared = cells.sum(axis=1)
                contributing_samples[idx] += shared > 0
                per_sample = pd.DataFrame(cells, columns=CELLS)
                per_sample.insert(0, "CRE_pair", chrom_pairs.CRE_pair.to_numpy())
                per_sample["chrom"] = chrom
                per_sample["sample_name"] = sample.sample_name
                per_sample["donor"] = sample.donor
                per_sample["n_shared_reads"] = shared
                per_sample.to_csv(sample_output, sep="\t", index=False, header=first_output)
                first_output = False
                qc.append((chrom, sample.sample_name, sample.donor, len(cov), orphans, int(shared.sum())))
                print(f"{chrom} {sample.sample_name}: {int((shared > 0).sum())} covered pairs", flush=True)
                del cov, acc

    keep = pooled.sum(axis=1) > 0
    counts = pooled[keep]
    result = pairs.loc[keep].reset_index(drop=True)
    result = pd.concat([result, pd.DataFrame(counts, columns=CELLS), test_counts(counts, args.pseudocount)], axis=1)
    result = bh_min_shared(summarize_coaccess_counts(result, args.pseudocount))
    result.insert(result.columns.get_loc("n_shared_reads") + 1,
                  "n_contributing_samples", contributing_samples[keep])
    result["dataset"] = "LCL"
    result["analysis_scope"] = args.chrom or "genome-wide"
    result["n_samples"] = len(samples)
    result["coverage_mode"] = coverage_mode
    result["read_rule"] = "contain"
    for field in ("fire_overlap_fraction", "min_dist", "max_dist", "pseudocount"):
        result[field] = getattr(args, field)
    annotations = gene_map.groupby("CRE_pair", sort=False).agg(
        n_genes=("gene_id", "nunique"),
        gene_ids=("gene_id", lambda x: ";".join(sorted(set(x)))),
        gene_names=("gene_name", lambda x: ";".join(sorted(set(x)))))
    result = result.merge(annotations, on="CRE_pair", how="left", validate="one_to_one")
    gene_result = gene_map[["CRE_pair", "gene_id", "gene_name", "transcript_type"]].drop_duplicates().merge(
        result, on="CRE_pair", validate="many_to_one")
    outputs = {"LCL_coaccess_pairs.tsv.gz": result, "LCL_coaccess_stat.tsv.gz": gene_result,
               "LCL_count_qc.tsv.gz": pd.DataFrame(qc, columns=["chrom", "sample_name", "donor", "n_cres_covered", "orphan_hits", "fiber_pair_observations"])}
    for name, table in outputs.items():
        destination = out_dir / name
        temporary = out_dir / (name + ".partial")
        table.to_csv(temporary, sep="\t", index=False, compression="gzip")
        temporary.replace(destination)
    sample_counts_tmp.replace(out_dir / "LCL_pair_sample_counts.tsv.gz")
    samples.to_csv(out_dir / "sample_manifest.tsv", sep="\t", index=False)
    (out_dir / "RERUN_PENDING.txt").unlink(missing_ok=True)
    print(f"Completed: {len(result)} pooled pairs; {int((result.fdr < .05).sum())} with FDR <0.05 -> {out_dir}", flush=True)


if __name__ == "__main__":
    main()
