#!/usr/bin/env python3
"""Survey fully spanning molecules, assign unique reads, then sample before ACF.

Only the selected molecules have their m6A blocks parsed or a binary row built.
The disk-backed pool is resumable at chromosome/sample transaction boundaries.
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


def canonical_genes(bed_path, bins_path, chromosomes):
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
    joined = joined[joined.tss >= 1000].drop(columns="_merge")
    joined = joined.sort_values(["chrom", "tss", "gene_id"]).reset_index(drop=True)
    joined["gene_index"] = np.arange(len(joined))
    joined["window_start"], joined["window_end"] = joined.tss - 1000, joined.tss + 1000
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


def build_pool(out, genes, sources, signature, seed, helpers, force=False):
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
        sub = genes[genes.chrom.eq(chrom)].sort_values("tss")
        anchors = sub.tss.to_numpy()
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
                    lo = np.searchsorted(anchors, st + 1000, side="left")
                    hi = np.searchsorted(anchors, en - 1000, side="right")
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


def sampled_binary(selected, sources, helpers):
    binary = np.zeros((len(selected), 2000), dtype=np.uint8)
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
                    found[row.row_index] = True
        log(f"binary: extracted {sample} {chrom}, {int(found[sub.row_index].sum())}/{len(sub)} selected molecules")
    if not found.all():
        raise RuntimeError(f"Could not recover {int((~found).sum())} selected exact BED records")
    return binary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--bins-tsv", type=Path, default=EXPRESSION)
    ap.add_argument("--canonical-bed", type=Path, default=Path(CANONICAL))
    ap.add_argument("--ft-root", type=Path, default=FT_ROOT)
    ap.add_argument("--per-bin", type=int, default=2500)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--timepoints", nargs="+", default=["LPS_0", "LPS_5", "LPS_10", "LPS_15"])
    ap.add_argument("--chrom", nargs="+", default=CHROMS)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if (args.per_bin < 1 or args.seed < 0 or len(set(args.timepoints)) != len(args.timepoints)
            or len(set(args.chrom)) != len(args.chrom) or not set(args.chrom).issubset(CHROMS)
            or not set(args.timepoints).issubset({"LPS_0", "LPS_5", "LPS_10", "LPS_15"})):
        ap.error("Require positive per-bin count, nonnegative seed, unique chr1-22/X/Y and LPS_0/5/10/15")
    out = args.out_dir.resolve()
    prepare_dirs(out)
    helpers = extraction_helpers()
    genes = canonical_genes(args.canonical_bed, args.bins_tsv, args.chrom)
    sources = [(sample, chrom, helpers.m6a_path(args.ft_root, sample, chrom).resolve())
               for sample in args.timepoints for chrom in args.chrom if chrom in set(genes.chrom)]
    states = source_states(sources)  # fails clearly if any required BED/index is missing
    params = {"seed": args.seed, "chromosomes": args.chrom, "samples": args.timepoints,
              "width": 2000, "source_states": states, "pysam": pysam.__version__}
    pool_signature = fingerprint([args.canonical_bed, args.bins_tsv, __file__, common.__file__, EXTRACT_SCRIPT], params)
    gene_path = out / "tables/eligible_genes.tsv"
    genes.to_csv(gene_path, sep="\t", index=False)
    pool = build_pool(out, genes, sources, pool_signature, args.seed, helpers, args.force)
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
    if not ((selected.read_start <= selected.window_start) & (selected.read_end >= selected.window_end)).all():
        raise RuntimeError("Sampling admitted a partially spanning molecule")
    binary = sampled_binary(selected, sources, helpers)
    selected["m6a_count"] = binary.sum(axis=1)
    row_ids = selected.read_id.to_numpy(dtype=str)
    assert_alignment(binary, selected, row_ids)
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
                  "orientation": "ascending genomic positions; no strand reversal; no reference-base filtering",
                  "pool_counts": counts}
    validation_path = out / "validation/sampling_validation.json"
    write_json(validation_path, validation)
    finish_stage(out, "01_sampling", signature, [binary_path, ids_path, metadata_path, validation_path], validation)


if __name__ == "__main__":
    main()
