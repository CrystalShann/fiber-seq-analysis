#!/usr/bin/env python3
"""Extract FiberHMM calls for the saved sample; no sampling or clustering."""
import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

import tss_common as common
from tss_common import (PROJECT, DEFAULT_OUT, assert_alignment, fingerprint, finish_stage,
                        log, prepare_dirs, read_metadata, stage_valid)

SAMPLER = Path(__file__).with_name("01_sample_tss_molecules.py")
spec = importlib.util.spec_from_file_location("tss_sampler_coordinates", SAMPLER)
sampler = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sampler)
TRACKS = {"footprint": ("Nucleosome >90 bp", 1), "tf": ("TF <60 bp", 2)}


def footprint_path(root, sample, chrom, kind):
    compact = sample.replace("_", "")
    return root / f"firehmm_{kind}" / compact / f"{compact}_hmm_extracted_{kind}_{chrom}.bed.gz"


def source_states(paths):
    return {str(p): {"size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns,
                     "index_sha256": common.sha256(str(p) + ".tbi")} for p in paths}


def parse_blocks(line):
    f = line.split("\t")
    if len(f) < 12:
        raise ValueError("Expected FiberHMM BED12")
    start, end, count = int(f[1]), int(f[2]), int(f[9])
    if start < 0 or end <= start or count < 0:
        raise ValueError("Invalid FiberHMM span/block count")
    if count == 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    sizes = np.array([int(x) for x in f[10].rstrip(",").split(",")], dtype=np.int64)
    offsets = np.array([int(x) for x in f[11].rstrip(",").split(",")], dtype=np.int64)
    if (len(sizes) != count or len(offsets) != count or (sizes <= 0).any()
            or (offsets < 0).any() or (start + offsets + sizes > end).any()):
        raise ValueError(f"Malformed FiberHMM blocks: {f[3]}")
    # All FiberHMM blocks are real calls, including first and last.
    return start + offsets, sizes


def oriented_intervals(row, starts, ends):
    # Map the two INCLUDED endpoint bases with the sampler's existing helper.
    # Restore a half-open column interval after ordering the mapped endpoints.
    first = sampler.oriented_column(row.strand, row.tss, row.window_offset_start, starts)
    last = sampler.oriented_column(row.strand, row.tss, row.window_offset_start, ends - 1)
    lo, hi = np.minimum(first, last), np.maximum(first, last) + 1
    if not np.array_equal(hi - lo, ends - starts):
        raise ValueError("Footprint width changed under orientation")
    for positions, columns in ((starts, first), (ends - 1, last)):
        back = sampler.column_genomic_position(row.strand, row.window_start, row.window_end, columns)
        if not np.array_equal(back, positions):
            raise ValueError("Footprint endpoint failed genomic/oriented round trip")
    return lo.astype(np.int64), hi.astype(np.int64)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--extract-root", type=Path, default=PROJECT / "macrophage_project/FiberHMM/extract")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    out = args.out_dir.resolve()
    prepare_dirs(out)
    metadata_path = out / "tables/plot_metadata.tsv"
    binary_path = out / "intermediate/binary_m6a.npy"
    ids_path = out / "intermediate/row_ids.npy"
    metadata = read_metadata(metadata_path)
    binary = np.load(binary_path, mmap_mode="r", allow_pickle=False)
    assert_alignment(binary, metadata, np.load(ids_path, allow_pickle=False))
    if not metadata.orientation.eq("transcriptional").all() or not metadata.strand.isin(["+", "-"]).all():
        raise ValueError("Expected strand-oriented sample metadata")
    sources = [(sample, chrom, kind, footprint_path(args.extract_root, sample, chrom, kind))
               for sample, chrom in metadata.groupby(["sample", "chrom"]).groups for kind in TRACKS]
    paths = [p for _, _, _, p in sources]
    states = source_states(paths)
    signature = fingerprint([metadata_path, binary_path, ids_path, __file__, common.__file__, SAMPLER],
                            {"sources": states, "orientation": "transcriptional", "pysam": pysam.__version__})
    if not args.force and stage_valid(out, "04b_footprints", signature):
        return
    categories = np.zeros(binary.shape, dtype=np.uint8)
    width = binary.shape[1]
    records, match_rows, mirrored_checks = [], [], []
    has_footprint = np.zeros(len(metadata), dtype=bool)
    helpers = sampler.extraction_helpers()
    # Process nucleosomes before TF; m6A overrides both after all sources.
    for sample, chrom, kind, path in sources:
        # merge_intervals requires ascending genomic starts; matrix writes still use row_index.
        selected = metadata[metadata["sample"].eq(sample) & metadata.chrom.eq(chrom)].sort_values("window_start")
        wanted = set(selected.read_id)
        candidates = {rid: set() for rid in wanted}
        with pysam.TabixFile(str(path)) as tb:
            for left, right in helpers.merge_intervals(selected.window_start, selected.window_end):
                for line in tb.fetch(chrom, int(left), int(right)):
                    f = line.split("\t", 4)
                    if f[3] in wanted:
                        candidates[f[3]].add(line)
        matched = 0
        for row in selected.itertuples(index=False):
            lines = candidates[row.read_id]
            exact = [line for line in sorted(lines)
                     if tuple(map(int, line.split("\t", 3)[1:3])) == (row.read_start, row.read_end)]
            # FiberHMM records can be trimmed to the outermost calls (especially TF).
            # A unique read/sample/chrom record is sufficient; use the read span
            # to resolve multiple records, never arbitrarily choose among them.
            chosen = next(iter(lines)) if len(lines) == 1 else exact[0] if len(exact) == 1 else None
            status = "matched" if chosen is not None else "ambiguous" if lines else "no_matching_record"
            match_rows.append((row.read_id, row.row_index, sample, chrom, kind,
                               len(lines), len(exact), status))
            if status != "matched":
                continue
            matched += 1
            starts, sizes = parse_blocks(chosen)
            keep = ((sizes < 60) if kind == "tf" else (sizes > 90))
            keep &= (starts < row.window_end) & (starts + sizes > row.window_start)
            sizes, starts = sizes[keep], starts[keep]
            ends = np.minimum(starts + sizes, row.window_end)
            starts = np.maximum(starts, row.window_start)
            if not len(starts):
                continue
            lo, hi = oriented_intervals(row, starts, ends)
            if (lo < 0).any() or (hi > width).any() or (hi <= lo).any():
                raise ValueError("Oriented footprint outside matrix")
            # Interval union via a difference array: vectorized, no per-base loop.
            delta = np.zeros(width + 1, dtype=np.int32)
            np.add.at(delta, lo, 1)
            np.add.at(delta, hi, -1)
            covered = np.cumsum(delta[:-1]) > 0
            categories[row.row_index, covered] = TRACKS[kind][1]
            has_footprint[row.row_index] = True
            records.extend((row.read_id, row.row_index, TRACKS[kind][0], int(s), int(e), int(z))
                           for s, e, z in zip(starts, ends, sizes))
            if row.strand == "-" and len(mirrored_checks) < 10 and row.read_id not in {x[0] for x in mirrored_checks}:
                # Independent full-mask reversal checks half-open endpoints and interior bases.
                genomic_delta = np.zeros(width + 1, dtype=np.int32)
                np.add.at(genomic_delta, starts - row.window_start, 1)
                np.add.at(genomic_delta, ends - row.window_start, -1)
                if not np.array_equal(covered, (np.cumsum(genomic_delta[:-1]) > 0)[::-1]):
                    raise ValueError("Minus-strand footprint differs from mirrored genomic mask")
                mirrored_checks.append((row.read_id, int(starts[0]), int(ends[0]), int(lo[0]), int(hi[0])))
        log(f"footprints: {sample} {chrom} {kind}: {matched}/{len(selected)} record matches")
    categories[binary == 1] = 3
    checks = {
        "shape_matches_binary": categories.shape == binary.shape,
        "codes_0_to_3": bool(np.isin(categories, [0, 1, 2, 3]).all()),
        "m6a_count_per_row": bool(np.array_equal((categories == 3).sum(axis=1), metadata.m6a_count.to_numpy())),
        "m6a_position_per_base": bool(np.array_equal(categories == 3, binary == 1)),
        "minus_strand_mirror": len(mirrored_checks) > 0 or not metadata.strand.eq("-").any(),
        "sources_unchanged": source_states(paths) == states,
    }
    if not all(checks.values()):
        raise ValueError(f"Footprint validation failed: {checks}")
    matrix_path = out / "intermediate/footprint_categories.bin"
    categories.tofile(matrix_path)  # uint8, C (row-major) order
    shape_path = out / "intermediate/footprint_categories_shape.tsv"
    pd.DataFrame([{"n_rows": len(metadata), "n_cols": width, "dtype": "uint8", "order": "C",
                   "orientation": "transcriptional"}]).to_csv(shape_path, sep="\t", index=False)
    records_path = out / "tables/footprint_records.tsv.gz"
    pd.DataFrame(records, columns=["read_id", "row_index", "track", "start", "end", "size"]).to_csv(
        records_path, sep="\t", index=False, compression="gzip")
    matches = pd.DataFrame(match_rows, columns=["read_id", "row_index", "sample", "chrom", "kind",
                                               "n_records", "n_exact_records", "status"])
    matches_path = out / "tables/footprint_matches.tsv"
    matches.to_csv(matches_path, sep="\t", index=False)
    summary = {"n_molecules": len(metadata), "n_with_footprint": int(has_footprint.sum()),
               "fraction_with_footprint": float(has_footprint.mean()),
               "n_no_matching_record": int(matches.loc[matches.status.eq("no_matching_record"), "read_id"].nunique()),
               "n_multiple_records": int(matches.loc[matches.n_records.gt(1), "read_id"].nunique()),
               "n_ambiguous_matches": int(matches.loc[matches.status.eq("ambiguous"), "read_id"].nunique()),
               "minus_strand_reads_checked": len(mirrored_checks), **checks}
    summary_path = out / "tables/footprint_validation.tsv"
    pd.DataFrame([summary]).to_csv(summary_path, sep="\t", index=False)
    mirrored_path = out / "tables/footprint_minus_strand_checks.tsv"
    pd.DataFrame(mirrored_checks, columns=["read_id", "genomic_start", "genomic_end",
                                          "oriented_column_start", "oriented_column_end"]).to_csv(
        mirrored_path, sep="\t", index=False)
    log(f"Footprint validation PASSED: {summary}")
    finish_stage(out, "04b_footprints", signature,
                 [matrix_path, shape_path, records_path, matches_path, summary_path, mirrored_path], summary)


if __name__ == "__main__":
    main()
