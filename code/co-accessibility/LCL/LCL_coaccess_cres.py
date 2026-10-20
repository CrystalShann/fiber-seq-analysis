#!/usr/bin/env python3
"""Genome-wide pooled co-accessibility across all 31 LCL samples.

Adapted from macrophage/03_coaccess_cres.py. Candidate cCREs overlap the union
of LCL FIRE peaks and share a canonical TSS +/-10 kb window. The interval gap
must satisfy 500 < gap < 20000. Shared reads have aligned bases covering the
FULL length of both cCREs: every base of each cCRE must
fall inside the read's M/= /X blocks, so a read whose alignment gap or end lands
inside either cCRE is not counted for that pair.
A cCRE is accessible when ONE FIRE element on that same read covers >=50% of
the cCRE length after intersecting with CIGAR M/= /X blocks. D/N gaps never
contribute coverage or accessibility; separate FIRE elements are not combined.

Sum raw 2x2 cells across samples before testing. Preserve macrophage two-sided
Fisher tests on both table+1 and raw counts; apply BH once across unique pooled
pairs. Donor information is metadata only; there are no LCL timepoints.
"""

import argparse
import gzip
import os
import sys
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "parsing_functions"))
from count_and_test_pairs import count_pairs, test_counts, summarize_coaccess_counts

ROOT = Path(os.environ.get("LCL_COACCESS_ROOT", "/project/spott/cshan/fiber-seq/LCL_project/co-accessibility"))
METATABLE = os.environ.get("LCL_SAMPLE_METATABLE", "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv")
BEDTOOLS = "/project/spott/cshan/envs/bedtools/bin/bedtools"
TABIX = "/project/spott/cshan/envs/dimelo/bin/tabix"
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
CELLS = ["co_closed", "CRE1_access", "CRE2_access", "co_access"]
COVERAGE_MODE = "aligned_blocks_v1"

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
    samples["donor"] = samples.cell_line
    samples["cram_path"] = [str(Path(r.fire_dir) / f"{r.sample_name}-fire-v0.1-filtered.cram") for r in samples.itertuples()]
    samples["peaks_path"] = [str(Path(r.fire_dir) / f"{r.sample_name}-fire-v0.1-peaks.bed.gz") for r in samples.itertuples()]
    samples["elements_path"] = [str(Path(r.fire_dir) / "additional-outputs-v0.1/fire-peaks" / f"{r.sample_name}-v0.1-fire-elements.bed.gz") for r in samples.itertuples()]
    samples["spans_path"] = [str(root / s / f"{s}.read_spans.bed.gz") for s in samples.sample_name]
    samples["blocks_path"] = [str(root / s / f"{s}.aligned_blocks.bed.gz") for s in samples.sample_name]
    return samples

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

# /project/spott/cshan/fiber-seq/LCL_project/co-accessibility/universe/cre_universe.bed -> 
  # this the the overlap ENCODE cCREs /project/spott/cshan/annotations/GRCh38-cCREs.bedand merged FIRE peaks across all 31 samples 
  # /project/spott/cshan/fiber-seq/LCL_project/co-accessibility/universe/fire_peaks_union.bed
  
  
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
# run bedtools and outputs error message
#################################################
def command_lines(command):
    """Stream bedtools output; do not hold every overlap line in memory."""
    with tempfile.TemporaryFile(mode="w+t") as err:
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=err, text=True)
        try:
            for line in proc.stdout:
                yield line.rstrip("\n").split("\t")
        finally:
            proc.stdout.close()
            code = proc.wait()
        if code:
            err.seek(0)
            raise RuntimeError(f"{command[0]} failed: {err.read()[:4000]}")
          
#################################################
# Convert bed12 to bed3
#################################################
# /project/spott/cshan/fiber-seq/LCL_project/co-accessibility/<sample>/<sample>.aligned_blocks.bed.gz
def aligned_blocks(fields):
    """Decode BED12 (without zeros at the end and beginning) from bamtobed -bed12 -splitD to a list of genomic coordinates"""
    # read the read's outer start and end, and the number of blocks from bed12
    start, end, count = int(fields[1]), int(fields[2]), int(fields[9])
    # size = size of the block
    # offsets = distance to the start 
    sizes = [int(x) for x in fields[10].rstrip(",").split(",")]
    offsets = [int(x) for x in fields[11].rstrip(",").split(",")]
    # block start = read start + offset
    # block end - block start + size
    blocks = [(start + offset, start + offset + size) for offset, size in zip(offsets, sizes)]
    if len(sizes) != count or len(offsets) != count or not blocks or any(
            s < start or e > end or e < s or (i and s < blocks[i-1][1])
            for i, (s, e) in enumerate(blocks)):
        raise ValueError(f"Invalid aligned BED12 blocks for {fields[3]}")
    return [(s, e) for s, e in blocks if e > s]

#################################################
# counts how many bases of an interval [start, end) are covered by a read’s aligned blocks
#################################################

# Example: A read has blocks (1000, 1200) and (1300, 1500), and the cCRE is (1150, 1350).
# Block 1: min(1200, 1350) − max(1000, 1150) = 1200 − 1150 = 50
# Block 2: min(1500, 1350) − max(1300, 1150) = 1350 − 1300 = 50
# Total = 100 bp

def aligned_overlap(blocks, start, end):
    # overlap length between the block and the interval
    return sum(max(0, min(e, end) - max(s, start)) for s, e in blocks)

#################################################
# get spanning reads, and accessible reads for every CRE
#################################################
def incidence_for_chrom(chrom, cre_bed, blocks_path, elements_path, fire_fraction, tmpdir):
    """Intersect each FIRE with aligned blocks."""
    # pull only chrom out of the sample's aligned blocks file and its FIRE elements file
    spans = tmpdir / "spans.bed"
    elements = tmpdir / "elements.bed"
    for source, destination in ((blocks_path, spans), (elements_path, elements)):
        with destination.open("w") as output:
            subprocess.run([TABIX, str(source), chrom], stdout=output, check=True)
    # maps cCRE ID to read ID
    coverage = defaultdict(set)
    blocks_by_read = {}
    # which read covers the CRE by at least 1bp
    #################################################
    # find every read that has outer span overlaps a cCRE, holds CRE and read side by side
    # -wa -wb prints both the cCRE record and the read record on one line
    command = [BEDTOOLS, "intersect", "-a", str(cre_bed), "-b", str(spans), "-wa", "-wb"]
    # pull out the read ID, and take the CRE end and start, 
    # count how many bases of CRE fall inside the read's aligned block
    for fields in command_lines(command):
        rid = fields[8]  # BED5 cCRE followed by BED12 primary alignment
        if rid not in blocks_by_read:
            blocks_by_read[rid] = aligned_blocks(fields[5:17])
        start, end = int(fields[1]), int(fields[2])
        overlap = aligned_overlap(blocks_by_read[rid], start, end)
        # the read's aligned blocks must cover every base of the CRE
        if overlap == end - start:
          # if the read pass, write the read ID in CRE coverage
            coverage[fields[3]].add(rid)
            
    # For the reads that cover CRE by at least 1bp, on which read is each CRE accessible?
    #################################################
    # accessible = FIRE element overlap with 50% of the CRE length
    
    # Example. The cCRE is 200 bp, so the threshold is 100 bp. One FIRE element covers 160 bp of it, 
    # but a deletion removes 70 bp of that stretch, leaving 90 aligned bp. That is below 100, 
    # so the read is not accessible. If the aligned part were 110 bp, it would be
    
    accessible = defaultdict(set)
    
    for fields in command_lines([BEDTOOLS, "intersect", "-a", str(cre_bed), "-b", str(elements),
                                  "-f", str(fire_fraction), "-wa", "-wb"]):
        rid = fields[8]
        # take only FIRE element that sits inside the CRE
        start, end = max(int(fields[1]), int(fields[6])), min(int(fields[2]), int(fields[7]))
        # count the # of aligned blocks, using the read's aligned blocks from coverage step
        overlap = aligned_overlap(blocks_by_read.get(rid, ()), start, end)
        if overlap > 0 and overlap >= fire_fraction * (int(fields[2]) - int(fields[1])):
            accessible[fields[3]].add(rid)
    # For each cCRE, only keep accessible regions that apear in CRE coverage read     
    orphan_hits = 0
    for cre_id in accessible:
        before = len(accessible[cre_id])
        accessible[cre_id] &= coverage.get(cre_id, set())
        orphan_hits += before - len(accessible[cre_id])
    spans.unlink(); elements.unlink()
    return coverage, accessible, orphan_hits

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
    require_files([p for col in ("blocks_path", "elements_path") for p in samples[col]] +
                  [p + ".tbi" for col in ("blocks_path", "elements_path") for p in samples[col]])
    for sample in samples.itertuples(index=False):
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
                cov, acc, orphans = incidence_for_chrom(chrom, cre_bed, sample.blocks_path, sample.elements_path,
                                                       args.fire_overlap_fraction, temp)
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
    result = summarize_coaccess_counts(result, args.pseudocount)
    result.insert(result.columns.get_loc("n_shared_reads") + 1,
                  "n_contributing_samples", contributing_samples[keep])
    result["dataset"] = "LCL"
    result["analysis_scope"] = args.chrom or "genome-wide"
    result["n_samples"] = len(samples)
    result["coverage_mode"] = COVERAGE_MODE
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
