#!/usr/bin/env python3
"""Survey fully spanning molecules, assign unique reads, then sample before ACF.

Only the selected molecules have their m6A blocks parsed or a binary row built.
The disk-backed pool is resumable at chromosome/sample transaction boundaries.

Windows are TSS-relative and strand-oriented: offsets are transcriptional,
negative = upstream, positive = downstream. For the interval [s, e), end
exclusive, with tss the zero-based TSS base (BED start + 10):
    + strand: genomic [tss + s, tss + e)
    - strand: genomic [tss - e + 1, tss - s + 1)
Tabix scans, interval merging and full-span eligibility run on these genomic
intervals. Binary rows are built in genomic order and reversed for - genes, so
column j always equals TSS-relative offset s + j.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sqlite3

import numpy as np
import pandas as pd
import pysam

import tss_common as common
from tss_common import (BINS, DEFAULT_OUT, PROJECT, assert_alignment, fingerprint,
                        finish_stage, log, prepare_dirs, stage_valid, write_json)

EXTRACT_SCRIPT = PROJECT / "code/accessibility/expr_access/02_tss_m6a_profiles.py"
CANONICAL = "/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed"
EXPRESSION = PROJECT / "macrophage_project/expr_access/tables/tss_expression_bins.tsv"
FT_ROOT = PROJECT / "macrophage_project/FiberHMM/extract/ft_result_dir"
CHROMS = [f"chr{x}" for x in list(range(1, 23)) + ["X", "Y"]]


def extraction_helpers():
    spec = importlib.util.spec_from_file_location("tss_existing_extraction", EXTRACT_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def stable_key(seed, domain, *values):
    payload = json.dumps([seed, domain, *values], separators=(",", ":"))
    return hashlib.blake2b(payload.encode(), digest_size=16).hexdigest()


def genomic_window(strand, tss, offset_start, offset_end):
    """Genomic [start, end) of the TSS-relative interval [offset_start, offset_end)."""
    plus = np.asarray(strand) == "+"
    tss = np.asarray(tss)
    start = np.where(plus, tss + offset_start, tss - offset_end + 1)
    end = np.where(plus, tss + offset_end, tss - offset_start + 1)
    return start, end


def oriented_column(strand, tss, offset_start, position):
    """Column of genomic base ``position`` in the oriented row (offset - offset_start)."""
    offset = np.where(np.asarray(strand) == "+", np.asarray(position) - tss, tss - np.asarray(position))
    return offset - offset_start


def column_genomic_position(strand, window_start, window_end, column):
    """Genomic base shown in oriented column ``column``."""
    return np.where(np.asarray(strand) == "+", np.asarray(window_start) + column,
                    np.asarray(window_end) - 1 - column)


def eligible_gene_range(window_starts, width, read_start, read_end):
    """[lo, hi) indices of genes whose genomic window the read fully spans.

    window_starts must be sorted; all windows share ``width``. A read spans
    [ws, ws + width) iff read_start <= ws and ws <= read_end - width.
    """
    lo = np.searchsorted(window_starts, read_start, side="left")
    hi = np.searchsorted(window_starts, read_end - width, side="right")
    return lo, hi


def canonical_genes(bed_path, bins_path, chromosomes, window_start, window_end):
    bed = pd.read_csv(bed_path, sep="\t", header=None,
                      names=["chrom", "start", "end", "name", "score", "strand"])
    fields = bed.name.str.split(";", expand=True)
    if fields.shape[1] < 5:
        raise ValueError("Canonical BED must encode gene ID, transcript, name, type and tags")
    bed["gene_id"], bed["gene_name"], bed["gene_type"] = fields[0], fields[2], fields[3]
    bed["canonical"] = fields[4].str.split(",").apply(lambda tags: "Ensembl_canonical" in tags)
    bed = bed[bed.gene_type.eq("protein_coding") & bed.canonical & bed.chrom.isin(chromosomes)].copy()
    if not (bed.end - bed.start).eq(20).all():
        raise ValueError("Canonical BED intervals must be 20 bp with TSS = start + 10")
    bed["tss"] = bed.start + 10
    keys = ["gene_id", "gene_name", "chrom", "tss", "strand"]
    if bed.gene_id.duplicated().any():
        raise ValueError("Canonical BED has duplicate protein-coding gene IDs")
    bins = pd.read_csv(bins_path, sep="\t")
    bins = bins[bins.expr_bin.isin(BINS) & bins.chrom.isin(chromosomes)].copy()
    if bins.empty or bins.gene_id.duplicated().any():
        raise ValueError("Expression table is empty or duplicates expressed gene IDs")
    joined = bins.merge(bed[keys + ["gene_type", "canonical"]], on=keys,
                        how="left", validate="one_to_one", indicator=True)
    if not joined._merge.eq("both").all():
        raise ValueError("Expression genes/coordinates do not match protein-coding canonical BED: "
                         + str(joined.loc[joined._merge.ne("both"), "gene_id"].tolist()[:10]))
    if (not np.isfinite(joined.mean_tpm).all() or (joined.mean_tpm < 1).any()
            or not joined.strand.isin(["+", "-"]).all()):
        raise ValueError("Invalid expression/strand metadata")
    joined = joined.drop(columns="_merge")
    # TSS-relative, strand-oriented offsets [s, e): negative upstream, positive downstream.
    joined["window_offset_start"], joined["window_offset_end"] = window_start, window_end
    joined["window_start"], joined["window_end"] = genomic_window(
        joined.strand, joined.tss, window_start, window_end)
    joined = joined[joined.window_start >= 0]
    joined = joined.sort_values(["chrom", "tss", "gene_id"]).reset_index(drop=True)
    joined["gene_index"] = np.arange(len(joined))
    return joined


def source_states(sources):
    # Hashing tens of GB of BED data would defeat indexed extraction. Use their
    # size/mtime plus a full digest of every tabix index; small inputs are hashed.
    return {str(path): {"size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns,
                       "index_sha256": common.sha256(str(path) + ".tbi")}
            for _, _, path in sources}


UPSERT = """INSERT INTO pool VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
ON CONFLICT(read_id) DO UPDATE SET
sample=excluded.sample,chrom=excluded.chrom,read_start=excluded.read_start,
read_end=excluded.read_end,span=excluded.span,gene_index=excluded.gene_index,
expr_bin=excluded.expr_bin,assignment_key=excluded.assignment_key,
selection_key=excluded.selection_key,line_sha=excluded.line_sha,read_strand=excluded.read_strand
WHERE excluded.span>pool.span OR
(excluded.span=pool.span AND excluded.assignment_key<pool.assignment_key)"""


def build_pool(out, genes, sources, signature, seed, helpers, window_start, window_end, force=False):
    target = out / "intermediate/eligible_pool.sqlite"
    if not force and stage_valid(out, "01_pool", signature):
        return target
    pending = out / "intermediate/eligible_pool.building.sqlite"
    progress = out / "validation/01_pool.inprogress.json"
    try:
        resume = not force and json.loads(progress.read_text()) == signature and pending.is_file()
    except (OSError, ValueError):
        resume = False
    if not resume:
        # Only this stage's incomplete, regenerable file is replaced.
        pending.unlink(missing_ok=True)
        write_json(progress, signature)
    db = sqlite3.connect(pending)
    db.execute("PRAGMA cache_size=-65536")
    db.execute("""CREATE TABLE IF NOT EXISTS pool (
        read_id TEXT PRIMARY KEY, sample TEXT, chrom TEXT, read_start INTEGER,
        read_end INTEGER, span INTEGER, gene_index INTEGER, expr_bin TEXT,
        assignment_key TEXT, selection_key TEXT, line_sha TEXT, read_strand TEXT)""")
    db.execute("CREATE TABLE IF NOT EXISTS completed (source TEXT PRIMARY KEY)")
    db.commit()
    completed = {r[0] for r in db.execute("SELECT source FROM completed")}
    for sample, chrom, path in sources:
        if str(path) in completed:
            log(f"pool: resume completed {sample} {chrom}")
            continue
        # Anchor on the genomic window start: both strands share one width, so
        # full-span eligibility is the same test on genomic coordinates.
        sub = genes[genes.chrom.eq(chrom)].sort_values(["window_start", "gene_id"])
        anchors = sub.window_start.to_numpy()
        width = window_end - window_start
        if not (sub.window_end.to_numpy() - anchors == width).all():
            raise ValueError("All genomic windows must share the TSS-relative width")
        indices = sub.gene_index.to_numpy()
        gene_ids = sub.gene_id.to_numpy()
        bins = sub.expr_bin.to_numpy()
        windows = helpers.merge_intervals(sub.window_start, sub.window_end)
        batch, n_records = [], 0
        log(f"pool: scan {sample} {chrom}, {len(sub)} TSSs / {len(windows)} merged windows")
        with db, pysam.TabixFile(str(path)) as tb:
            for ws, we in windows:
                for line in tb.fetch(chrom, int(ws), int(we)):
                    # Deliberately do not parse the potentially large m6A block arrays.
                    f = line.split("\t", 6)
                    st, en, rid = int(f[1]), int(f[2]), f[3]
                    if not rid or rid == "." or st < 0 or en <= st or f[5] not in ("+", "-", "."):
                        raise ValueError(f"Malformed read ID/span/strand in {path}: {rid}")
                    # Read spans genomic [ws, ws + width) <=> st <= ws <= en - width
                    lo, hi = eligible_gene_range(anchors, width, st, en)
                    if hi <= lo:
                        continue
                    line_sha = hashlib.sha256(line.encode()).hexdigest()
                    key, j = min((stable_key(seed, "assignment", rid, sample, chrom,
                                            st, en, gene_ids[j], line_sha), j) for j in range(lo, hi))
                    batch.append((rid, sample, chrom, st, en, en-st, int(indices[j]), str(bins[j]),
                                  key, stable_key(seed, "sampling", rid),
                                  line_sha, f[5]))
                    n_records += 1
                    if len(batch) >= 10000:
                        db.executemany(UPSERT, batch)
                        batch.clear()
            if batch:
                db.executemany(UPSERT, batch)
            db.execute("INSERT INTO completed VALUES (?)", (str(path),))
        log(f"pool: committed {sample} {chrom}, {n_records:,} eligible fetched records (before dedup)")
    db.execute("CREATE INDEX IF NOT EXISTS bin_rank ON pool(expr_bin,selection_key,read_id)")
    db.commit()
    if db.execute("PRAGMA quick_check").fetchone()[0] != "ok":
        raise RuntimeError("Eligible pool SQLite integrity check failed")
    counts = dict(db.execute("SELECT expr_bin,COUNT(*) FROM pool GROUP BY expr_bin"))
    db.close()
    target.parent.mkdir(parents=True, exist_ok=True)
    pending.replace(target)
    count_path = out / "tables/eligible_pool_counts.tsv"
    pd.DataFrame({"expr_bin": BINS, "n_unique_assigned_eligible": [counts.get(b, 0) for b in BINS]}).to_csv(
        count_path, sep="\t", index=False)
    finish_stage(out, "01_pool", signature, [target, count_path],
                 {"pool_counts": counts, "assignment": "longest eligible alignment; seeded gene/tie assignment"})
    return target


def sampled_binary(selected, sources, helpers, width):
    binary = np.zeros((len(selected), width), dtype=np.uint8)
    found = np.zeros(len(selected), dtype=bool)
    lookup = {row.read_id: row for row in selected.itertuples(index=False)}
    for sample, chrom, path in sources:
        sub = selected[selected["sample"].eq(sample) & selected.chrom.eq(chrom)].sort_values("window_start")
        if sub.empty:
            continue
        windows = helpers.merge_intervals(sub.window_start, sub.window_end)
        with pysam.TabixFile(str(path)) as tb:
            for ws, we in windows:
                for line in tb.fetch(chrom, int(ws), int(we)):
                    prefix = line.split("\t", 4)
                    row = lookup.get(prefix[3])
                    if row is None or found[row.row_index] or row.sample != sample or row.chrom != chrom:
                        continue
                    if hashlib.sha256(line.encode()).hexdigest() != row.line_sha:
                        continue
                    f = line.split("\t")
                    if len(f) != 12:
                        raise ValueError(f"Expected BED12 for {row.read_id}")
                    starts = np.fromstring(f[11].rstrip(","), sep=",", dtype=np.int64)
                    sizes = np.fromstring(f[10].rstrip(","), sep=",", dtype=np.int64)
                    if len(starts) != int(f[9]) or len(sizes) != len(starts) or len(starts) < 2:
                        raise ValueError(f"Malformed BED12 blocks: {row.read_id}")
                    calls = starts[1:-1] + row.read_start
                    if (not np.all(sizes[1:-1] == 1) or np.any(calls < row.read_start)
                            or np.any(calls >= row.read_end)):
                        raise ValueError(f"Non-base/interior m6A call outside read span: {row.read_id}")
                    offsets = calls[(calls >= row.window_start) & (calls < row.window_end)] - row.window_start
                    binary[row.row_index, offsets] = 1
                    if row.strand == "-":
                        # Genomic order -> transcriptional order: column j = offset s + j.
                        binary[row.row_index] = binary[row.row_index][::-1]
                    found[row.row_index] = True
        log(f"binary: extracted {sample} {chrom}, {int(found[sub.row_index].sum())}/{len(sub)} selected molecules")
    if not found.all():
        raise RuntimeError(f"Could not recover {int((~found).sum())} selected exact BED records")
    return binary


def check_offset_zero_maps_to_tss(selected):
    """Pure-arithmetic check: on both strands, the column of offset 0 is genomic base tss."""
    s = selected.window_offset_start.to_numpy()
    column = oriented_column(selected.strand, selected.tss.to_numpy(), s, selected.tss.to_numpy())
    if not np.array_equal(column, -s):
        raise RuntimeError("Offset 0 does not map to column -window_offset_start")
    back = column_genomic_position(selected.strand, selected.window_start.to_numpy(),
                                   selected.window_end.to_numpy(), column)
    if not np.array_equal(back, selected.tss.to_numpy()):
        raise RuntimeError("The column for offset 0 does not map back to the genomic TSS base")
    for strand in ("+", "-"):
        if not selected.strand.eq(strand).any():
            raise RuntimeError(f"No sampled {strand}-strand molecules; orientation check needs both strands")
    return {"offset0_column": "-window_offset_start", "checked_rows": int(len(selected)),
            "strands": sorted(selected.strand.unique().tolist())}


def verify_orientation_on_reads(binary, selected, sources, helpers, per_strand=5):
    """Re-fetch several real reads per strand and check every m6A call's oriented column."""
    picks = pd.concat([selected[selected.strand.eq(strand)].head(per_strand) for strand in ("+", "-")])
    if picks.strand.nunique() != 2:
        raise RuntimeError("Orientation verification needs sampled reads on both strands")
    verified = []
    for row in picks.itertuples(index=False):
        path = next(p for smp, chrom, p in sources if smp == row.sample and chrom == row.chrom)
        with pysam.TabixFile(str(path)) as tb:
            record = next(line for line in tb.fetch(row.chrom, int(row.window_start), int(row.window_end))
                          if hashlib.sha256(line.encode()).hexdigest() == row.line_sha)
        f = record.split("\t")
        calls = np.fromstring(f[11].rstrip(","), sep=",", dtype=np.int64)[1:-1] + row.read_start
        inside = calls[(calls >= row.window_start) & (calls < row.window_end)]
        columns = oriented_column(row.strand, row.tss, row.window_offset_start, inside)
        expected = np.zeros(binary.shape[1], dtype=np.uint8)
        expected[columns] = 1
        if not np.array_equal(expected, binary[row.row_index]):
            raise RuntimeError(f"Oriented row differs from re-fetched calls: {row.read_id}")
        zero = -row.window_offset_start
        tss_is_call = bool(np.isin(row.tss, inside))
        if 0 <= zero < binary.shape[1] and bool(binary[row.row_index, zero]) != tss_is_call:
            raise RuntimeError(f"Offset-0 column disagrees with the genomic TSS base: {row.read_id}")
        if not np.array_equal(column_genomic_position(row.strand, row.window_start, row.window_end, columns),
                              inside):
            raise RuntimeError(f"Column -> genomic mapping is not invertible: {row.read_id}")
        verified.append({"read_id": row.read_id, "strand": row.strand, "tss": int(row.tss),
                         "n_calls_in_window": int(len(inside)), "tss_base_called": tss_is_call})
    log(f"orientation: verified {len(verified)} real reads ({per_strand} per strand) against re-fetched calls")
    return verified


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--bins-tsv", type=Path, default=EXPRESSION)
    ap.add_argument("--canonical-bed", type=Path, default=Path(CANONICAL))
    ap.add_argument("--ft-root", type=Path, default=FT_ROOT)
    ap.add_argument("--per-bin", type=int, default=2500)
    ap.add_argument("--window-start", type=int, default=-1000,
                    help="Window start as a TSS-relative transcriptional offset (inclusive; negative = upstream).")
    ap.add_argument("--window-end", type=int, default=1000,
                    help="Window end as a TSS-relative transcriptional offset (exclusive; positive = downstream).")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--timepoints", nargs="+", default=["LPS_0", "LPS_5", "LPS_10", "LPS_15"])
    ap.add_argument("--chrom", nargs="+", default=CHROMS)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--test", action="store_true",
                    help="Also re-fetch several real reads per strand and verify their oriented columns.")
    args = ap.parse_args()
    if (args.per_bin < 1 or args.seed < 0 or len(set(args.timepoints)) != len(args.timepoints)
            or len(set(args.chrom)) != len(args.chrom) or not set(args.chrom).issubset(CHROMS)
            or not set(args.timepoints).issubset({"LPS_0", "LPS_5", "LPS_10", "LPS_15"})):
        ap.error("Require positive per-bin count, nonnegative seed, unique chr1-22/X/Y and LPS_0/5/10/15")
    width = args.window_end - args.window_start
    if width < 3:
        ap.error("--window-end must exceed --window-start by at least three bases")
    out = args.out_dir.resolve()
    prepare_dirs(out)
    helpers = extraction_helpers()
    genes = canonical_genes(args.canonical_bed, args.bins_tsv, args.chrom, args.window_start, args.window_end)
    sources = [(sample, chrom, helpers.m6a_path(args.ft_root, sample, chrom).resolve())
               for sample in args.timepoints for chrom in args.chrom if chrom in set(genes.chrom)]
    states = source_states(sources)  # fails clearly if any required BED/index is missing
    params = {"seed": args.seed, "chromosomes": args.chrom, "samples": args.timepoints,
              "window_start": args.window_start, "window_end": args.window_end, "width": width,
              "source_states": states, "pysam": pysam.__version__}
    pool_signature = fingerprint([args.canonical_bed, args.bins_tsv, __file__, common.__file__, EXTRACT_SCRIPT], params)
    gene_path = out / "tables/eligible_genes.tsv"
    genes.to_csv(gene_path, sep="\t", index=False)
    pool = build_pool(out, genes, sources, pool_signature, args.seed, helpers,
                      args.window_start, args.window_end, args.force)
    signature = fingerprint([pool, gene_path, __file__, common.__file__], {"per_bin": args.per_bin, **params})
    if not args.force and stage_valid(out, "01_sampling", signature):
        return
    with sqlite3.connect(f"file:{pool}?mode=ro", uri=True) as db:
        counts = dict(db.execute("SELECT expr_bin, COUNT(*) FROM pool GROUP BY expr_bin"))
        log(f"Unique, assigned fully spanning pool counts: {counts}")
        if any(counts.get(b, 0) < args.per_bin for b in BINS):
            report = {"required_per_bin": args.per_bin, "available": {b: counts.get(b, 0) for b in BINS}}
            write_json(out / "validation/insufficient_pool.json", report)
            raise SystemExit("Insufficient unique eligible molecules; no oversampling. " + json.dumps(report))
        parts = [pd.read_sql_query("SELECT * FROM pool WHERE expr_bin=? ORDER BY selection_key,read_id LIMIT ?",
                                   db, params=(b, args.per_bin)) for b in BINS]
    selected = pd.concat(parts, ignore_index=True).merge(
        genes.drop(columns=["expr_bin", "chrom"]), on="gene_index", how="left", validate="many_to_one", sort=False)
    selected.insert(0, "row_index", np.arange(len(selected)))
    selected["timepoint"] = selected["sample"]
    selected["lps_minutes"] = selected["sample"].str.extract(r"LPS_(\d+)$", expand=False).astype(int)
    selected["orientation"] = "transcriptional"
    if not ((selected.read_start <= selected.window_start) & (selected.read_end >= selected.window_end)).all():
        raise RuntimeError("Sampling admitted a partially spanning molecule")
    binary = sampled_binary(selected, sources, helpers, width)
    selected["m6a_count"] = binary.sum(axis=1)
    row_ids = selected.read_id.to_numpy(dtype=str)
    assert_alignment(binary, selected, row_ids)
    orientation_check = check_offset_zero_maps_to_tss(selected)
    if args.test:
        orientation_check["verified_reads"] = verify_orientation_on_reads(binary, selected, sources, helpers)
    if source_states(sources) != states:
        raise RuntimeError("Source BEDs/indexes changed during extraction")
    binary_path, ids_path = out / "intermediate/binary_m6a.npy", out / "intermediate/row_ids.npy"
    metadata_path = out / "tables/sampled_molecules.tsv"
    np.save(binary_path, binary, allow_pickle=False)
    np.save(ids_path, row_ids, allow_pickle=False)
    selected.to_csv(metadata_path, sep="\t", index=False)
    validation = {"all_protein_coding_canonical": bool(selected.gene_type.eq("protein_coding").all() and selected.canonical.all()),
                  "not_expressed_excluded": bool(selected.expr_bin.isin(BINS).all()),
                  "all_fully_spanning": True, "unique_physical_read_ids": True,
                  "balanced_counts": selected.expr_bin.value_counts().to_dict(),
                  "n_molecules": len(selected), "binary_shape": list(binary.shape),
                  "row_order_verified": True, "seed": args.seed,
                  "window_offsets_from_tss": [args.window_start, args.window_end],
                  "orientation": "transcriptional",
                  "orientation_definition": ("column j = TSS-relative offset window_offset_start + j; negative = "
                                             "upstream, positive = downstream; + strand genomic [tss+s, tss+e), "
                                             "- strand genomic [tss-e+1, tss-s+1) with the row reversed"),
                  "acf_orientation_note": ("the ACF is invariant to sequence reversal, so per-read ACF values "
                                           "change only through which reads and bases fall in each window; "
                                           "heatmaps, sliding windows and other position-based outputs depend "
                                           "on orientation"),
                  "reference_base_filtering": "none",
                  "offset_zero_check": orientation_check,
                  "pool_counts": counts}
    validation_path = out / "validation/sampling_validation.json"
    write_json(validation_path, validation)
    finish_stage(out, "01_sampling", signature, [binary_path, ids_path, metadata_path, validation_path], validation)


if __name__ == "__main__":
    main()
