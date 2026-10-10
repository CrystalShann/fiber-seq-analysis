#!/usr/bin/env python3
"""Genome-wide pooled co-accessibility across all 31 LCL samples.

Adapted from macrophage/03_coaccess_cres.py. Candidate cCREs overlap the union
of LCL FIRE peaks and share a canonical TSS +/-10 kb window. The interval gap
must satisfy 500 < gap < 20000. Shared reads have aligned bases at both cCREs.
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
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio as conditional_odds_ratio

ROOT = Path(os.environ.get("LCL_COACCESS_ROOT", "/project/spott/cshan/fiber-seq/LCL_project/co-accessibility"))
METATABLE = os.environ.get("LCL_SAMPLE_METATABLE", "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv")
BEDTOOLS = "/project/spott/cshan/envs/bedtools/bin/bedtools"
TABIX = "/project/spott/cshan/envs/dimelo/bin/tabix"
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
CELLS = ["co_closed", "CRE1_access", "CRE2_access", "co_access"]
COVERAGE_MODE = "aligned_blocks_v1"


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


def require_files(paths):
    missing = [str(p) for p in paths if not Path(p).is_file() or Path(p).stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing/empty inputs:\n" + "\n".join(missing))


def build_pairs(cre, membership, min_dist, max_dist):
    """Gene memberships plus one stable coordinate-oriented record per pair."""
    pos = cre.set_index("CRE_ID")
    records = []
    for gene_id, group in membership.groupby("gene_id", sort=False):
        ids = sorted(group.CRE_ID.unique(), key=lambda c: (pos.at[c, "chrom"], pos.at[c, "start"], pos.at[c, "end"], c))
        gene = group.iloc[0]
        for i, c1 in enumerate(ids):
            for c2 in ids[i + 1:]:
                if pos.at[c1, "chrom"] != pos.at[c2, "chrom"]:
                    continue
                gap = max(0, int(pos.at[c2, "start"]) - int(pos.at[c1, "end"]))
                if gap >= max_dist:
                    break
                if gap > min_dist:
                    records.append((c1, c2, gap, gene_id, gene.gene_name, gene.transcript_type))
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


def aligned_blocks(fields):
    """Decode sentinel-free BED12 from bamtobed -bed12 -splitD."""
    start, end, count = int(fields[1]), int(fields[2]), int(fields[9])
    sizes = [int(x) for x in fields[10].rstrip(",").split(",")]
    offsets = [int(x) for x in fields[11].rstrip(",").split(",")]
    blocks = [(start + offset, start + offset + size) for offset, size in zip(offsets, sizes)]
    if len(sizes) != count or len(offsets) != count or not blocks or any(
            s < start or e > end or e < s or (i and s < blocks[i-1][1])
            for i, (s, e) in enumerate(blocks)):
        raise ValueError(f"Invalid aligned BED12 blocks for {fields[3]}")
    # bedtools emits zero-size blocks at terminal/adjacent D/N operators.
    # They contain no aligned bases and are not an error or a QC exclusion.
    return [(s, e) for s, e in blocks if e > s]


def aligned_overlap(blocks, start, end):
    return sum(max(0, min(e, end) - max(s, start)) for s, e in blocks)


def incidence_for_chrom(chrom, cre_bed, blocks_path, elements_path, read_rule, fire_fraction, tmpdir):
    """Intersect each original FIRE with its own primary aligned blocks."""
    spans = tmpdir / "spans.bed"
    elements = tmpdir / "elements.bed"
    for source, destination in ((blocks_path, spans), (elements_path, elements)):
        with destination.open("w") as output:
            subprocess.run([TABIX, str(source), chrom], stdout=output, check=True)
    coverage = defaultdict(set)
    blocks_by_read = {}
    command = [BEDTOOLS, "intersect", "-a", str(cre_bed), "-b", str(spans), "-wa", "-wb"]
    for fields in command_lines(command):
        rid = fields[8]  # BED5 cCRE followed by BED12 primary alignment.
        if rid not in blocks_by_read:
            blocks_by_read[rid] = aligned_blocks(fields[5:17])
        start, end = int(fields[1]), int(fields[2])
        overlap = aligned_overlap(blocks_by_read[rid], start, end)
        if overlap > 0 and (read_rule == "any" or overlap == end - start):
            coverage[fields[3]].add(rid)
    accessible = defaultdict(set)
    for fields in command_lines([BEDTOOLS, "intersect", "-a", str(cre_bed), "-b", str(elements),
                                  "-f", str(fire_fraction), "-wa", "-wb"]):
        # The BED interval fraction is only a cheap upper-bound prefilter.
        # Sum aligned pieces of THIS element, never pieces of different elements.
        rid = fields[8]
        start, end = max(int(fields[1]), int(fields[6])), min(int(fields[2]), int(fields[7]))
        overlap = aligned_overlap(blocks_by_read.get(rid, ()), start, end)
        if overlap > 0 and overlap >= fire_fraction * (int(fields[2]) - int(fields[1])):
            accessible[fields[3]].add(rid)
    orphan_hits = 0
    for cre_id in accessible:
        before = len(accessible[cre_id])
        accessible[cre_id] &= coverage.get(cre_id, set())
        orphan_hits += before - len(accessible[cre_id])
    spans.unlink(); elements.unlink()
    return coverage, accessible, orphan_hits


def count_pairs(pairs, coverage, accessible):
    counts = np.zeros((len(pairs), 4), dtype=np.int64)
    empty = set()
    for i, pair in enumerate(pairs.itertuples(index=False)):
        shared = coverage.get(pair.CRE1, empty) & coverage.get(pair.CRE2, empty)
        a = accessible.get(pair.CRE1, empty) & shared
        b = accessible.get(pair.CRE2, empty) & shared
        both = len(a & b)
        counts[i] = len(shared) - len(a) - len(b) + both, len(a) - both, len(b) - both, both
    return counts


def bh(p):
    p = np.asarray(p, dtype=float)
    if not len(p):
        return p.copy()
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p)); out[order] = np.minimum(adjusted, 1)
    return out


def test_counts(counts, pseudocount):
    """Add pseudocount only once, AFTER pooling all 31 samples."""
    unique, inverse = np.unique(counts, axis=0, return_inverse=True)
    rows = []
    for closed, first, second, both in unique:
        raw = np.array([[closed, second], [first, both]], dtype=np.int64)
        corrected = raw + pseudocount
        rows.append((fisher_exact(corrected, alternative="two-sided").pvalue,
                     conditional_odds_ratio(corrected, kind="conditional").statistic,
                     fisher_exact(raw, alternative="two-sided").pvalue))
    stats = np.array(rows, dtype=float)[inverse] if rows else np.empty((0, 3))
    return pd.DataFrame(stats, columns=["pval", "fisher_estimate", "pval_raw"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--sample-metatable", default=METATABLE)
    ap.add_argument("--out-dir", type=Path)
    ap.add_argument("--prepare-manifest", action="store_true", help="Validate all31 source inputs and write universe/sample_manifest.tsv only")
    ap.add_argument("--read-rule", choices=["any", "contain"], default="any")
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
                                                       args.read_rule, args.fire_overlap_fraction, temp)
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
    # No per-donor filtering: any shared fiber contributes to the pooled table.
    keep = pooled.sum(axis=1) > 0
    counts = pooled[keep]
    result = pairs.loc[keep].reset_index(drop=True)
    result = pd.concat([result, pd.DataFrame(counts, columns=CELLS), test_counts(counts, args.pseudocount)], axis=1)
    result["n_shared_reads"] = counts.sum(axis=1)
    result["n_contributing_samples"] = contributing_samples[keep]
    result["co_inaccess"] = result.co_closed
    pc = args.pseudocount
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
    result["dataset"] = "LCL"
    result["analysis_scope"] = args.chrom or "genome-wide"
    result["n_samples"] = len(samples)
    result["coverage_mode"] = COVERAGE_MODE
    for field in ("read_rule", "fire_overlap_fraction", "min_dist", "max_dist", "pseudocount"):
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
