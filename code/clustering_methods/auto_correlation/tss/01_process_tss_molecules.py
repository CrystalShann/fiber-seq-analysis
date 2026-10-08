#!/usr/bin/env python3
"""Macrophage canonical-TSS autocorrelation: every Python processing stage in one script.

Stages, run in order on one TSS region (all share the --out-dir working folder):
    sample     fully spanning molecules per canonical protein-coding TSS, unique
               read assignment, balanced seeded sample, strand-oriented 0/1 matrix
    acf        unchanged parent autocorrelations() at every 1 bp lag of the window
    cluster    unchanged parent cluster_profiles(): PCA -> correlation kNN -> Leiden -> UMAP
    tables     mean/median ACF per cluster and expression bin, composition, heatmap order
    footprints FiberHMM nucleosome/TF calls mapped into the oriented matrix frame
    nrl        per-molecule NRL with a prominence threshold calibrated on A/T-shuffled
               molecules, and group mean-ACF peak tests against group-sized shuffled means

Windows are TSS-relative and strand-oriented (transcriptional): negative = upstream,
positive = downstream, end exclusive. For [s, e) and the zero-based TSS base tss
(the 1 bp BED start): + strand genomic [tss + s, tss + e); - strand genomic
[tss - e + 1, tss - s + 1) with the row reversed, so column j = offset s + j.

Parent scripts in .. are imported unchanged: 03_compute_autocorrelations.py,
04_cluster_autocorrelations.py and 05_summarize_autocorrelations.py.
"""
import argparse
from datetime import datetime
import hashlib
import importlib.util
from importlib.metadata import version
import json
import os
from pathlib import Path
import sqlite3
import sys

import numpy as np
import pandas as pd
import pysam
from scipy.signal import find_peaks
from scipy.stats import median_abs_deviation

HERE = Path(__file__).resolve().parent
PARENT_DIR = HERE.parent
PROJECT = HERE.parents[3]
DEFAULT_OUT = PROJECT / "macrophage_project/auto_correlation/tss"
PARENT_ACF = PARENT_DIR / "03_compute_autocorrelations.py"
PARENT_CLUSTER = PARENT_DIR / "04_cluster_autocorrelations.py"
PARENT_SUMMARY = PARENT_DIR / "05_summarize_autocorrelations.py"
EXTRACT_SCRIPT = PROJECT / "code/accessibility/expr_access/02_tss_m6a_profiles.py"
CANONICAL = "/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed"
EXPRESSION = PROJECT / "macrophage_project/expr_access/tables/tss_expression_bins.tsv"
FT_ROOT = PROJECT / "macrophage_project/FiberHMM/extract/ft_result_dir"
EXTRACT_ROOT = PROJECT / "macrophage_project/FiberHMM/extract"
REF_FA = "/project/spott/reference/human/GRCh38/hg38.fa"  # as in code/accessibility/expr_access/02_tss_m6a_profiles.py
BINS = ("Q1_low", "Q2", "Q3", "Q4_high")
CHROMS = [f"chr{x}" for x in list(range(1, 23)) + ["X", "Y"]]
TRACKS = {"footprint": ("Nucleosome >90 bp", 1), "tf": ("TF <60 bp", 2)}
SHUFFLE_MODES = ("at", "all")
NRL_STATUSES = ("ok", "no_peak", "zero_variance", "window_too_short", "peak_negative")
NRL_METRICS = ["nrl_bp", "nrl_status", "nrl_peak_height", "peak2_lag_bp", "peak2_height",
               "peak_ratio", "repeat_peak_lag", "repeat_peak_value"]
NRL_META = ["read_id", "sample", "timepoint", "lps_minutes", "gene_id", "gene_name", "chrom", "tss", "strand",
            "expr_bin", "mean_tpm", "m6a_count", "window", "window_offset_start", "window_offset_end", "cluster"]
GROUP_COLUMNS = ["group_peak_lag_bp", "group_peak_height", "group_peak_prominence", "group_peak_pvalue"]


# ---------------------------------------------------------------------------
# Shared I/O and content-checked stage records (no ACF/clustering algorithms)
# ---------------------------------------------------------------------------

def log(message):
    print(f"[{datetime.now().isoformat(timespec='seconds')}] {message}", flush=True)


def load_parent(path):
    name = "tss_parent_" + Path(path).stem
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


def extraction_helpers():
    spec = importlib.util.spec_from_file_location("tss_existing_extraction", EXTRACT_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(paths, parameters):
    return {"inputs": {str(Path(p).resolve()): sha256(p) for p in paths},
            "parameters": json.loads(json.dumps(parameters, default=str))}


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True, default=str) + "\n")
    temporary.replace(path)


def prepare_dirs(outdir):
    for name in ("intermediate", "tables", "validation", "logs", "plots"):
        (Path(outdir) / name).mkdir(parents=True, exist_ok=True)


def finish_stage(outdir, stage, signature, outputs, details=None):
    """Record the stage's input fingerprint and output digests; stages are never skipped."""
    outdir = Path(outdir).resolve()
    digests = {}
    for item in outputs:
        path = Path(item)
        if not path.is_absolute():
            path = outdir / path
        digests[str(path.resolve().relative_to(outdir))] = sha256(path)
    write_json(outdir / "validation" / (stage + ".manifest.json"),
               {"signature": signature, "outputs": digests, "details": details or {}})
    log(f"{stage}: complete ({len(outputs)} verified outputs)")


def read_metadata(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False,
                       dtype={"cluster": str, "read_id": str, "gene_id": str})


def window_width(metadata):
    widths = (metadata.window_end - metadata.window_start).unique()
    if len(widths) != 1 or widths[0] < 3:
        raise ValueError("All sampled molecules must share one window of at least three bases")
    return int(widths[0])


def assert_alignment(binary, metadata, row_ids):
    if binary.shape != (len(metadata), window_width(metadata)):
        raise ValueError("Binary matrix must have one window-width row per molecule")
    if not np.array_equal(metadata.row_index.to_numpy(), np.arange(len(metadata))):
        raise ValueError("Metadata row_index must match matrix order exactly")
    if metadata.read_id.duplicated().any() or not np.array_equal(
            metadata.read_id.astype(str).to_numpy(), row_ids.astype(str)):
        raise ValueError("Duplicate physical read IDs or changed row order")
    if not np.isin(binary, [0, 1]).all():
        raise ValueError("Binary matrix contains values other than zero and one")
    if not np.array_equal(binary.sum(axis=1), metadata.m6a_count.to_numpy()):
        raise ValueError("m6A counts and binary matrix differ")


def load_matrix(out):
    """Oriented binary matrix, row IDs and (metadata-free) alignment inputs of the working folder."""
    binary = np.load(out / "intermediate/binary_m6a.npy", mmap_mode="r", allow_pickle=False)
    row_ids = np.load(out / "intermediate/row_ids.npy", allow_pickle=False)
    return binary, row_ids


def load_acf(out, metadata):
    profiles = np.load(out / "intermediate/acf.npy", mmap_mode="r", allow_pickle=False)
    valid = np.load(out / "intermediate/acf_valid.npy", allow_pickle=False)
    if profiles.ndim != 2 or len(profiles) != len(metadata) or valid.shape != (len(metadata),):
        raise ValueError("ACF rows are not aligned with the sampled metadata")
    if not np.isfinite(profiles[valid]).all() or not np.isnan(profiles[~valid]).all():
        raise ValueError("Valid ACFs must be finite and zero-variance ACFs must be NaN")
    if valid.any() and not np.allclose(profiles[valid, 0], 1, rtol=0, atol=1e-10):
        raise ValueError("Lag 0 must equal one in valid ACF rows")
    return profiles, valid


# ---------------------------------------------------------------------------
# Stage: sample
# ---------------------------------------------------------------------------

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
    if not (bed.end - bed.start).eq(1).all():
        raise ValueError("Canonical BED intervals must be 1 bp with TSS = start")
    bed["tss"] = bed.start
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


def source_states(paths):
    # Hashing tens of GB of BED data would defeat indexed extraction. Use their
    # size/mtime plus a full digest of every tabix index; small inputs are hashed.
    return {str(path): {"size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns,
                        "index_sha256": sha256(str(path) + ".tbi")} for path in paths}


UPSERT = """INSERT INTO pool VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
ON CONFLICT(read_id) DO UPDATE SET
sample=excluded.sample,chrom=excluded.chrom,read_start=excluded.read_start,
read_end=excluded.read_end,span=excluded.span,gene_index=excluded.gene_index,
expr_bin=excluded.expr_bin,assignment_key=excluded.assignment_key,
selection_key=excluded.selection_key,line_sha=excluded.line_sha,read_strand=excluded.read_strand
WHERE excluded.span>pool.span OR
(excluded.span=pool.span AND excluded.assignment_key<pool.assignment_key)"""


def build_pool(out, genes, sources, signature, seed, helpers, window_start, window_end):
    """Disk-backed pool of fully spanning reads; resumable at sample/chromosome boundaries."""
    target = out / "intermediate/eligible_pool.sqlite"
    pending = out / "intermediate/eligible_pool.building.sqlite"
    progress = out / "validation/sample_pool.inprogress.json"
    try:
        resume = json.loads(progress.read_text()) == signature and pending.is_file()
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
                    batch.append((rid, sample, chrom, st, en, en - st, int(indices[j]), str(bins[j]),
                                  key, stable_key(seed, "sampling", rid), line_sha, f[5]))
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
    finish_stage(out, "sample_pool", signature, [target, count_path],
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


def stage_sample(args, out):
    width = args.window_end - args.window_start
    helpers = extraction_helpers()
    genes = canonical_genes(args.canonical_bed, args.bins_tsv, args.chrom, args.window_start, args.window_end)
    sources = [(sample, chrom, helpers.m6a_path(args.ft_root, sample, chrom).resolve())
               for sample in args.timepoints for chrom in args.chrom if chrom in set(genes.chrom)]
    states = source_states([p for _, _, p in sources])  # fails clearly if any required BED/index is missing
    params = {"seed": args.seed, "chromosomes": args.chrom, "samples": args.timepoints,
              "window_start": args.window_start, "window_end": args.window_end, "width": width,
              "source_states": states, "pysam": pysam.__version__}
    pool_signature = fingerprint([args.canonical_bed, args.bins_tsv, __file__, EXTRACT_SCRIPT], params)
    gene_path = out / "tables/eligible_genes.tsv"
    genes.to_csv(gene_path, sep="\t", index=False)
    pool = build_pool(out, genes, sources, pool_signature, args.seed, helpers, args.window_start, args.window_end)
    signature = fingerprint([pool, gene_path, __file__], {"per_bin": args.per_bin, **params})
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
    if source_states([p for _, _, p in sources]) != states:
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
                                           "heatmaps and other position-based outputs depend on orientation"),
                  "reference_base_filtering": "none",
                  "offset_zero_check": orientation_check,
                  "pool_counts": counts}
    validation_path = out / "validation/sampling_validation.json"
    write_json(validation_path, validation)
    finish_stage(out, "sample", signature, [binary_path, ids_path, metadata_path, validation_path], validation)


# ---------------------------------------------------------------------------
# Stage: acf (unchanged parent autocorrelations)
# ---------------------------------------------------------------------------

def stage_acf(args, out):
    binary_path, metadata_path, rows_path = (out / "intermediate/binary_m6a.npy",
                                             out / "tables/sampled_molecules.tsv", out / "intermediate/row_ids.npy")
    metadata = read_metadata(metadata_path)
    width = window_width(metadata)
    max_lag = width - 1 if args.max_lag is None else args.max_lag
    if not 0 <= max_lag < width:
        raise SystemExit(f"--max-lag must be in [0, {width - 1}] for the {width}-bp input window")
    signature = fingerprint([binary_path, metadata_path, rows_path, Path(__file__), PARENT_ACF],
                            {"max_lag": max_lag, "numpy_version": np.__version__})
    binary, row_ids = load_matrix(out)
    assert_alignment(binary, metadata, row_ids)
    if not len(binary):
        raise ValueError("ACF input must contain sampled molecules")
    log(f"acf: computing parent autocorrelations() for {len(binary):,} sampled reads, "
        f"lags 0–{max_lag}; no signal transformations")
    parent = load_parent(PARENT_ACF)
    profiles, valid = parent.autocorrelations(binary, n_features=max_lag + 1)
    expected_valid = (binary.sum(axis=1) > 0) & (binary.sum(axis=1) < binary.shape[1])
    if profiles.shape != (len(metadata), max_lag + 1) or valid.shape != (len(metadata),):
        raise ValueError("Parent ACF changed the sampled row count or lag count")
    if valid.dtype != np.bool_ or not np.array_equal(valid, expected_valid):
        raise ValueError("ACF validity flags do not match the binary rows with nonzero variance")
    if not np.isfinite(profiles[valid]).all() or not np.isnan(profiles[~valid]).all():
        raise ValueError("Valid ACF rows must be finite; zero-variance ACF rows must be entirely NaN")
    if not np.allclose(profiles[valid, 0], 1.0, rtol=0, atol=1e-10):
        raise ValueError("Lag 0 differs from one for nonzero-variance reads")
    outputs = [out / "intermediate/acf.npy", out / "intermediate/acf_valid.npy", out / "validation/acf_validation.json"]
    np.save(outputs[0], profiles, allow_pickle=False)
    np.save(outputs[1], valid, allow_pickle=False)
    saved = np.load(outputs[0], mmap_mode="r", allow_pickle=False)
    if saved.shape != (len(metadata), max_lag + 1):
        raise ValueError(f"Saved ACF rows must have exactly max_lag + 1 = {max_lag + 1} columns")
    report = {
        "n_reads": len(metadata), "n_valid": int(valid.sum()),
        "n_zero_variance": int((~valid).sum()), "matrix_shape": list(profiles.shape),
        "window_width_bp": width, "max_lag_bp": max_lag, "lag0_equals_one": True,
        "saved_columns_equal_max_lag_plus_one": True, "lag_resolution_bp": 1,
        "lag0_max_absolute_error": float(np.max(np.abs(profiles[valid, 0] - 1))) if valid.any() else None,
        "zero_variance_flagged": True, "invalid_rows_all_nan": True,
        "metadata_and_matrix_row_order_preserved": True,
        "parent_function": str(PARENT_ACF) + ":autocorrelations",
        "acf_definition": "sum((x[t]-mean(x))*(x[t+k]-mean(x))) / (N*var(x))",
        "binary_orientation": "transcriptional (column j = TSS-relative offset window_offset_start + j)",
        "acf_orientation_note": ("the ACF is invariant to sequence reversal, so per-read ACF values change "
                                 "only through which reads and bases fall in each strand-oriented window"),
    }
    write_json(outputs[2], report)
    finish_stage(out, "acf", signature, outputs, details=report)
    log(f"acf: finished; {int(valid.sum()):,} valid and {int((~valid).sum()):,} zero-variance reads")


# ---------------------------------------------------------------------------
# Stage: cluster (unchanged parent cluster_profiles, all expression bins jointly)
# ---------------------------------------------------------------------------

def stage_cluster(args, out):
    binary_path, metadata_path, rows_path = (out / "intermediate/binary_m6a.npy",
                                             out / "tables/sampled_molecules.tsv", out / "intermediate/row_ids.npy")
    parameters = {"n_pcs": args.n_pcs, "n_neighbors": args.n_neighbors,
                  "resolution": args.resolution, "seed": args.seed, "metric": "correlation",
                  "joint_expression_bins": True,
                  "package_versions": {p: version(p) for p in
                      ("numpy", "scipy", "scanpy", "anndata", "scikit-learn", "umap-learn", "leidenalg", "igraph")}}
    signature = fingerprint([binary_path, metadata_path, rows_path, out / "intermediate/acf.npy",
                             out / "intermediate/acf_valid.npy", Path(__file__), PARENT_CLUSTER], parameters)
    binary, row_ids = load_matrix(out)
    metadata = read_metadata(metadata_path)
    assert_alignment(binary, metadata, row_ids)
    profiles, valid = load_acf(out, metadata)
    expected_valid = (binary.sum(axis=1) > 0) & (binary.sum(axis=1) < binary.shape[1])
    if profiles.shape[1] < 3:
        raise ValueError("Aligned ACF profiles with at least three lag columns are required for Leiden")
    if valid.dtype != np.bool_ or not np.array_equal(valid, expected_valid):
        raise ValueError("ACF validity mask is not aligned to the sampled nonzero-variance molecules")
    log(f"cluster: jointly clustering {int(valid.sum()):,} valid profiles from {len(metadata):,} "
        f"sampled reads; {args.n_pcs} PCs, {args.n_neighbors} neighbors, resolution {args.resolution}")
    parent = load_parent(PARENT_CLUSTER)
    labels, status, embedding, info, edges = parent.cluster_profiles(
        profiles, valid, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors,
        resolution=args.resolution, seed=args.seed)
    n = len(metadata)
    if labels.shape != (n,) or status.shape != (n,) or embedding.shape != (n, 2):
        raise ValueError("Parent clustering changed the sampled row count")
    clustered = status == "clustered"
    if not clustered.any() or info.get("n_clusters", 0) < 1:
        raise RuntimeError("Leiden produced no clusters: " + str(dict(zip(*np.unique(status, return_counts=True)))))
    if (clustered & ~valid).any() or np.any(status[~valid] != "zero_variance"):
        raise ValueError("Zero-variance molecules were not excluded from Leiden clustering")
    if not np.isfinite(embedding[clustered]).all() or np.any(labels[clustered] == ""):
        raise ValueError("Clustered reads must have labels and finite UMAP coordinates")
    if not np.isnan(embedding[~clustered]).all() or np.any(labels[~clustered] != ""):
        raise ValueError("Unclustered reads must retain missing UMAP coordinates and empty parent labels")
    if edges.ndim != 2 or edges.shape[1] != 3 or not np.isfinite(edges).all():
        raise ValueError("Invalid parent connectivity graph")
    if len(edges):
        indices = edges[:, :2]
        if (np.any(indices != np.floor(indices)) or np.any(indices < 0) or np.any(indices >= n)
                or np.any(edges[:, 2] <= 0) or not clustered[indices.astype(int)].all()):
            raise ValueError("Connectivity edges do not index clustered rows in the original sample order")
    records = metadata.copy()
    records["cluster"] = np.where(clustered, labels, "Unclustered")
    records["status"] = status
    records["acf_valid"] = valid
    records["umap1"] = embedding[:, 0]
    records["umap2"] = embedding[:, 1]
    assert_alignment(binary, records, row_ids)
    outputs = [out / "tables/clustered_molecules.tsv", out / "intermediate/umap.npy",
               out / "intermediate/graph_edges.npy", out / "validation/clustering_info.json"]
    records.to_csv(outputs[0], sep="\t", index=False, na_rep="NA")
    np.save(outputs[1], embedding, allow_pickle=False)
    np.save(outputs[2], edges, allow_pickle=False)
    info.update(
        n_input=n, n_acf_valid=int(valid.sum()), n_zero_variance=int((~valid).sum()),
        n_unclustered=int((~clustered).sum()), n_graph_edges=len(edges),
        joint_expression_bins=True, expression_bins=sorted(records.expr_bin.unique().tolist()),
        status_counts={str(k): int(v) for k, v in zip(*np.unique(status, return_counts=True))},
        zero_variance_excluded=True, metadata_and_matrix_row_order_preserved=True,
        parent_function=str(PARENT_CLUSTER) + ":cluster_profiles", package_versions=parameters["package_versions"])
    write_json(outputs[3], info)
    finish_stage(out, "cluster", signature, outputs, details=info)
    log(f"cluster: finished; {info['n_clusters']} joint clusters, {int(clustered.sum()):,} assigned reads")


# ---------------------------------------------------------------------------
# Stage: tables (plot tables, heatmap ranks, composition, cluster summaries)
# ---------------------------------------------------------------------------

def aggregate(profiles, mask, group_column, group):
    """Summaries preserve every retained lag; mask includes only valid ACF rows."""
    subset = profiles[mask]
    width = profiles.shape[1]
    return pd.DataFrame({
        group_column: group, "lag_bp": np.arange(width),
        "mean_acf": subset.mean(axis=0) if len(subset) else np.full(width, np.nan),
        "median_acf": np.median(subset, axis=0) if len(subset) else np.full(width, np.nan),
        "n_reads": len(subset)})


def stage_tables(args, out):
    binary_path, sampled_path, rows_path = (out / "intermediate/binary_m6a.npy",
                                            out / "tables/sampled_molecules.tsv", out / "intermediate/row_ids.npy")
    metadata_path, embedding_path = out / "tables/clustered_molecules.tsv", out / "intermediate/umap.npy"
    signature = fingerprint(
        [binary_path, sampled_path, metadata_path, rows_path, out / "intermediate/acf.npy",
         out / "intermediate/acf_valid.npy", embedding_path, Path(__file__), PARENT_SUMMARY],
        {"expression_bins": list(BINS), "heatmap_rank": "numeric cluster then descending m6a_count then row_index",
         "heatmap_matrix_order": "original sampled row order", "numpy_version": np.__version__,
         "pandas_version": pd.__version__})
    binary, row_ids = load_matrix(out)
    embedding = np.load(embedding_path, allow_pickle=False)
    sampled = read_metadata(sampled_path)
    records = read_metadata(metadata_path)
    # Metadata loading preserves literal gene/sample names; decode only numeric missing values.
    for column in ("umap1", "umap2"):
        records[column] = pd.to_numeric(records[column].replace("NA", np.nan), errors="raise")
    assert_alignment(binary, sampled, row_ids)
    assert_alignment(binary, records, row_ids)
    pd.testing.assert_frame_equal(records[sampled.columns], sampled,
                                  check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12)
    profiles, valid = load_acf(out, records)
    expected_valid = (binary.sum(axis=1) > 0) & (binary.sum(axis=1) < binary.shape[1])
    if (valid.dtype != np.bool_ or not np.array_equal(valid, records.acf_valid.to_numpy())
            or not np.array_equal(valid, expected_valid)):
        raise ValueError("Metadata ACF validity flags differ from the saved validity mask")
    if embedding.shape != (len(records), 2) or not np.allclose(
            records[["umap1", "umap2"]].to_numpy(), embedding, equal_nan=True, rtol=1e-12, atol=1e-12):
        raise ValueError("UMAP coordinates are not aligned to the sampled metadata")
    if set(records.expr_bin) != set(BINS):
        raise ValueError("Plotting requires exactly Q1_low, Q2, Q3, Q4_high expression bins")
    clustered = records.status.eq("clustered").to_numpy()
    if not clustered.any() or (clustered & ~valid).any():
        raise ValueError("Plotting requires at least one cluster and no clustered zero-variance reads")
    if not records.loc[~clustered, "cluster"].eq("Unclustered").all():
        raise ValueError("All excluded molecules must retain the explicit Unclustered label")
    cluster_ids = sorted(records.loc[clustered, "cluster"].unique().tolist(), key=int)
    # Always retain an Unclustered composition category so its zero count is explicit.
    composition_groups = cluster_ids + ["Unclustered"]
    cluster_order = {cluster: i for i, cluster in enumerate(composition_groups)}
    rank_order = records.assign(_cluster_order=records.cluster.map(cluster_order)).sort_values(
        ["_cluster_order", "m6a_count", "row_index"], ascending=[True, False, True], kind="stable").index.to_numpy()
    ranks = np.empty(len(records), dtype=np.int64)
    ranks[rank_order] = np.arange(1, len(records) + 1)
    records["heatmap_rank"] = ranks
    if not np.array_equal(np.sort(ranks), np.arange(1, len(records) + 1)):
        raise ValueError("Heatmap ranks must be a unique one-based permutation")
    log(f"tables: preparing {len(records):,} original-order rows and {len(cluster_ids)} global clusters")
    parent = load_parent(PARENT_SUMMARY)
    width = binary.shape[1]
    cluster_curves, expression_curves, summaries, composition = [], [], [], []
    for cluster in composition_groups:
        mask = records.cluster.eq(cluster).to_numpy()
        curve = aggregate(profiles, mask & valid, "cluster", cluster)
        if cluster != "Unclustered":
            cluster_curves.append(curve)
        peak_lag, peak_height = parent.repeat_peak(curve.mean_acf.to_numpy())
        n = int(mask.sum())
        counts = records.loc[mask, "m6a_count"]
        summaries.append({
            "cluster": cluster, "n_reads": n, "n_acf_valid": int((mask & valid).sum()),
            "n_zero_variance": int((mask & ~valid).sum()),
            "n_unclustered": int((mask & ~clustered).sum()),
            "fraction_of_sample": n / len(records),
            "fraction_clustered": n / int(clustered.sum()) if cluster != "Unclustered" else np.nan,
            "mean_m6a_count": counts.mean() if n else np.nan,
            "median_m6a_count": counts.median() if n else np.nan,
            "mean_m6a_call_fraction": counts.mean() / width if n else np.nan,
            "positive_local_peak_140_250_bp": peak_lag, "peak_acf": peak_height,
            "peak_band_complete": profiles.shape[1] >= 252})
        for expr_bin in BINS:
            in_bin = records.expr_bin.eq(expr_bin).to_numpy()
            count = int((mask & in_bin).sum())
            total_bin = int(in_bin.sum())
            composition.append({
                "cluster": cluster, "expr_bin": expr_bin, "n_reads": count,
                "n_bin_total": total_bin, "fraction_within_bin": count / total_bin,
                "fraction_within_cluster": count / n if n else np.nan,
                "n_bin_zero_variance": int((in_bin & ~valid).sum()),
                "n_bin_unclustered": int((in_bin & ~clustered).sum())})
    # Zero-variance and otherwise unclustered molecules per expression bin and timepoint.
    unclustered_rows = []
    for expr_bin in BINS:
        for timepoint in sorted(records.timepoint.unique()):
            cell = records.expr_bin.eq(expr_bin).to_numpy() & records.timepoint.eq(timepoint).to_numpy()
            unclustered_rows.append({
                "expr_bin": expr_bin, "timepoint": timepoint, "n_reads": int(cell.sum()),
                "n_clustered": int((cell & clustered).sum()), "n_zero_variance": int((cell & ~valid).sum()),
                "n_unclustered_other": int((cell & ~clustered & valid).sum()),
                "n_unclustered": int((cell & ~clustered).sum()),
                "fraction_unclustered": float((cell & ~clustered).sum() / cell.sum()) if cell.any() else np.nan})
    unclustered_table = pd.DataFrame(unclustered_rows)
    if unclustered_table.n_reads.sum() != len(records) or unclustered_table.n_unclustered.sum() != int((~clustered).sum()):
        raise ValueError("Unclustered counts by bin and timepoint must account for every sampled molecule")
    for expr_bin in BINS:
        expression_curves.append(aggregate(profiles, valid & records.expr_bin.eq(expr_bin).to_numpy(), "expr_bin", expr_bin))
    composition_table = pd.DataFrame(composition)
    fractions = composition_table.groupby("expr_bin").fraction_within_bin.sum()
    if not np.allclose(fractions.to_numpy(), 1, rtol=0, atol=1e-12):
        raise ValueError("Cluster composition must account for all sampled molecules in each bin")
    outputs = [out / "tables/plot_metadata.tsv", out / "tables/acf_heatmap.tsv.gz",
               out / "tables/cluster_acf.tsv", out / "tables/expression_acf.tsv",
               out / "tables/cluster_composition.tsv", out / "tables/cluster_summary.tsv",
               out / "tables/unclustered_counts.tsv", out / "validation/plot_table_validation.json"]
    records.to_csv(outputs[0], sep="\t", index=False, na_rep="NA")
    heatmap = pd.DataFrame(profiles, columns=[f"lag_{i}" for i in range(profiles.shape[1])])
    heatmap.insert(0, "row_index", records.row_index.to_numpy())
    # Full precision, all lags and every original row; display sorting happens only in R.
    heatmap.to_csv(outputs[1], sep="\t", index=False, na_rep="NA",
                   compression={"method": "gzip", "compresslevel": 6, "mtime": 0}, chunksize=128)
    pd.concat(cluster_curves, ignore_index=True).to_csv(outputs[2], sep="\t", index=False, na_rep="NA")
    pd.concat(expression_curves, ignore_index=True).to_csv(outputs[3], sep="\t", index=False, na_rep="NA")
    composition_table.to_csv(outputs[4], sep="\t", index=False, na_rep="NA")
    pd.DataFrame(summaries).to_csv(outputs[5], sep="\t", index=False, na_rep="NA")
    unclustered_table.to_csv(outputs[6], sep="\t", index=False, na_rep="NA")
    report = {
        "n_reads": len(records), "n_acf_valid": int(valid.sum()),
        "n_clustered": int(clustered.sum()), "n_unclustered": int((~clustered).sum()),
        "n_clusters": len(cluster_ids), "n_lags": profiles.shape[1],
        "metadata_and_matrix_row_order_preserved": True, "lag0_equals_one": True,
        "heatmap_rank_unique_one_based": True,
        "heatmap_sort": "numeric global Leiden cluster; descending m6a count; original row_index breaks ties",
        "heatmap_matrix_original_row_order": True, "cluster_labels_shared_across_bins": True,
        "composition_denominator": "all sampled molecules in each bin, including Unclustered",
        "composition_fractions_sum_to_one": True,
        "expression_acf_population": "all nonzero-variance sampled reads, including any valid unclustered reads",
        "repeat_peak_parent_function": str(PARENT_SUMMARY) + ":repeat_peak",
        "unclustered_counts": "unclustered_counts.tsv: per expression bin and timepoint"}
    write_json(outputs[7], report)
    finish_stage(out, "tables", signature, outputs, details=report)
    log("tables: finished; full-resolution matrix, curves, compositions and metadata are ready for R")


# ---------------------------------------------------------------------------
# Stage: footprints (FiberHMM calls in the oriented matrix frame)
# ---------------------------------------------------------------------------

def footprint_path(root, sample, chrom, kind):
    compact = sample.replace("_", "")
    return root / f"firehmm_{kind}" / compact / f"{compact}_hmm_extracted_{kind}_{chrom}.bed.gz"


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
    # Map the two INCLUDED endpoint bases with the sampler's coordinate helper.
    # Restore a half-open column interval after ordering the mapped endpoints.
    first = oriented_column(row.strand, row.tss, row.window_offset_start, starts)
    last = oriented_column(row.strand, row.tss, row.window_offset_start, ends - 1)
    lo, hi = np.minimum(first, last), np.maximum(first, last) + 1
    if not np.array_equal(hi - lo, ends - starts):
        raise ValueError("Footprint width changed under orientation")
    for positions, columns in ((starts, first), (ends - 1, last)):
        back = column_genomic_position(row.strand, row.window_start, row.window_end, columns)
        if not np.array_equal(back, positions):
            raise ValueError("Footprint endpoint failed genomic/oriented round trip")
    return lo.astype(np.int64), hi.astype(np.int64)


def stage_footprints(args, out):
    metadata_path = out / "tables/plot_metadata.tsv"
    binary_path, ids_path = out / "intermediate/binary_m6a.npy", out / "intermediate/row_ids.npy"
    metadata = read_metadata(metadata_path)
    binary, row_ids = load_matrix(out)
    assert_alignment(binary, metadata, row_ids)
    if not metadata.orientation.eq("transcriptional").all() or not metadata.strand.isin(["+", "-"]).all():
        raise ValueError("Expected strand-oriented sample metadata")
    sources = [(sample, chrom, kind, footprint_path(args.extract_root, sample, chrom, kind))
               for sample, chrom in metadata.groupby(["sample", "chrom"]).groups for kind in TRACKS]
    paths = [p for _, _, _, p in sources]
    states = source_states(paths)
    signature = fingerprint([metadata_path, binary_path, ids_path, __file__],
                            {"sources": states, "orientation": "transcriptional", "pysam": pysam.__version__})
    categories = np.zeros(binary.shape, dtype=np.uint8)
    width = binary.shape[1]
    records, match_rows, mirrored_checks = [], [], []
    has_footprint = np.zeros(len(metadata), dtype=bool)
    helpers = extraction_helpers()
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
            match_rows.append((row.read_id, row.row_index, sample, chrom, kind, len(lines), len(exact), status))
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
        "sources_unchanged": source_states(paths) == states}
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
    finish_stage(out, "footprints", signature,
                 [matrix_path, shape_path, records_path, matches_path, summary_path, mirrored_path], summary)


# ---------------------------------------------------------------------------
# Per-molecule NRL (SAMOSA secondary-peak scan) and shuffled nulls
# ---------------------------------------------------------------------------

def find_nrl(acf, min_lag, min_prominence, nrl_min=120, nrl_max=300):
    """Return (nrl_lag, nrl_height, peak2_lag, peak2_height, status) on the raw parent ACF.

    Lags below min_lag are ignored (lag-0 shoulder); find_peaks with the calibrated
    prominence finds local maxima; the NRL is the first positive local maximum in
    [nrl_min, nrl_max]; the second peak is the first positive local maximum within
    +/-50 bp of 2 x NRL inside the available lags (Abdulhay et al. 2020, eLife).
    """
    acf = np.asarray(acf, dtype=np.float64)
    if min_prominence is None or not np.isfinite(min_prominence) or min_prominence < 0:
        raise ValueError("min_prominence must be a finite nonnegative number (calibrated or overridden)")
    max_lag = len(acf) - 1
    if not np.isfinite(acf).all():
        return np.nan, np.nan, np.nan, np.nan, "zero_variance"
    if max_lag < nrl_min or max_lag <= min_lag:
        return np.nan, np.nan, np.nan, np.nan, "window_too_short"
    peaks, _ = find_peaks(acf[min_lag:], prominence=min_prominence)
    peaks = peaks + min_lag
    in_band = peaks[(peaks >= nrl_min) & (peaks <= nrl_max)]
    if not len(in_band):
        return np.nan, np.nan, np.nan, np.nan, "no_peak"
    positive = in_band[acf[in_band] > 0]
    if not len(positive):
        return np.nan, np.nan, np.nan, np.nan, "peak_negative"
    nrl = int(positive[0])
    second = peaks[(peaks >= 2 * nrl - 50) & (peaks <= 2 * nrl + 50) & (peaks <= max_lag) & (acf[peaks] > 0)]
    if len(second):
        peak2 = int(second[0])
        return nrl, float(acf[nrl]), peak2, float(acf[peak2]), "ok"
    return nrl, float(acf[nrl]), np.nan, np.nan, "ok"


def molecule_metrics(acf, min_prominence, min_lag, nrl_min, nrl_max, repeat_peak=None):
    """Per-molecule NRL metrics; peak_ratio = raw second-peak height / raw NRL-peak height."""
    acf = np.asarray(acf, dtype=np.float64)
    nrl, nrl_height, peak2, peak2_height, status = find_nrl(acf, min_lag, min_prominence, nrl_min, nrl_max)
    result = {column: np.nan for column in NRL_METRICS}
    result.update(nrl_bp=nrl, nrl_status=status, nrl_peak_height=nrl_height,
                  peak2_lag_bp=peak2, peak2_height=peak2_height)
    if status == "zero_variance":
        return result
    if repeat_peak is not None:
        lag, value = repeat_peak(acf)
        result.update(repeat_peak_lag=lag, repeat_peak_value=value)
    if status == "ok" and np.isfinite(peak2) and nrl_height > 0:
        result["peak_ratio"] = peak2_height / nrl_height
    return result


def at_mask(fasta, metadata, width):
    """Strand-oriented A/T indicator (True = reference base A or T) per molecule window."""
    mask = np.zeros((len(metadata), width), dtype=bool)
    for row in metadata.itertuples(index=False):
        sequence = fasta.fetch(row.chrom, int(row.window_start), int(row.window_end)).upper()
        if len(sequence) != width:
            raise ValueError(f"Reference window length mismatch for {row.read_id}")
        values = np.frombuffer(sequence.encode(), dtype=np.uint8)
        at = (values == ord("A")) | (values == ord("T"))
        mask[row.row_index] = at[::-1] if row.strand == "-" else at
    return mask


def allowed_positions(binary, mode, at=None):
    """Boolean matrix of positions a shuffled call may occupy, and a count check."""
    allowed = np.asarray(at, dtype=bool) if mode == "at" else np.ones(binary.shape, dtype=bool)
    if allowed.shape != binary.shape:
        raise ValueError("Allowed-position mask must match the binary matrix")
    if (binary.sum(axis=1) > allowed.sum(axis=1)).any():
        raise ValueError("A molecule has more m6A calls than allowed shuffle positions")
    return allowed


def shuffle_rows(binary, allowed, rng):
    """Same m6A count per row, calls placed uniformly among that row's allowed positions."""
    out = np.zeros(binary.shape, dtype=np.uint8)
    for i in range(len(binary)):
        count = int(binary[i].sum())
        if count:
            out[i, rng.choice(np.flatnonzero(allowed[i]), count, replace=False)] = 1
    return out


def check_shuffle(original, shuffled, allowed, mode):
    """Shuffled rows keep the m6A count and, in at mode, call only allowed positions."""
    if not np.array_equal(np.asarray(original).sum(axis=1), shuffled.sum(axis=1)):
        raise ValueError("Shuffled rows changed the m6A count")
    if mode == "at" and (shuffled.astype(bool) & ~allowed).any():
        raise ValueError("Shuffled rows placed calls outside A/T positions")


def band_peak(acf, min_lag, nrl_min, nrl_max):
    """Most prominent positive local maximum of one ACF row inside the NRL band.

    Returns (lag, height, prominence); (nan, nan, 0.0) when there is no such
    peak and (nan, nan, nan) when the row is not finite.
    """
    acf = np.asarray(acf, dtype=np.float64)
    if not np.isfinite(acf).all():
        return np.nan, np.nan, np.nan
    peaks, properties = find_peaks(acf[min_lag:], prominence=0)
    peaks = peaks + min_lag
    keep = (peaks >= nrl_min) & (peaks <= nrl_max) & (acf[peaks] > 0)
    if not keep.any():
        return np.nan, np.nan, 0.0
    prominences = properties["prominences"][keep]
    best = int(np.argmax(prominences))
    lag = int(peaks[keep][best])
    return lag, float(acf[lag]), float(prominences[best])


def acf_in_chunks(acf_fn, rows, chunk=2000):
    profiles, valid = [], []
    for start in range(0, len(rows), chunk):
        p, v = acf_fn(np.ascontiguousarray(rows[start: start + chunk]))
        profiles.append(p)
        valid.append(v)
    return np.concatenate(profiles), np.concatenate(valid)


def calibrate_prominence(binary, allowed, mode, acf_fn, n_null, quantile, min_lag, nrl_min, nrl_max, seed):
    """Calibrated min_prominence: the quantile of the largest positive in-band peak
    prominence of parent ACFs of one shuffle of up to n_null seeded molecules.

    Returns (threshold, details, table) with one table row per null molecule.
    """
    rng = np.random.default_rng(seed)
    picked = np.arange(len(binary))
    if len(binary) > n_null:
        picked = np.sort(rng.choice(len(binary), n_null, replace=False))
    shuffled = shuffle_rows(binary[picked], allowed[picked], rng)
    check_shuffle(binary[picked], shuffled, allowed[picked], mode)
    acf, valid = acf_in_chunks(acf_fn, shuffled)
    if acf.shape != (len(picked), binary.shape[1]):
        raise ValueError("Null ACF shape differs from the molecule window")
    null = np.full(len(picked), np.nan)
    for i in np.flatnonzero(valid):
        null[i] = band_peak(acf[i], min_lag, nrl_min, nrl_max)[2]
    table = pd.DataFrame({"row_index": picked, "acf_valid": valid, "null_max_prominence": null})
    finite = null[np.isfinite(null)]
    if not len(finite):
        raise ValueError("Calibration needs at least one nonzero-variance shuffled molecule")
    threshold = float(np.quantile(finite, quantile))
    details = {"method": (f"quantile of the largest positive in-band peak prominence of the parent ACF of "
                          f"{mode}-shuffled copies of sampled molecules (same m6A count per molecule)"),
               "shuffle_mode": mode, "window_width": int(binary.shape[1]), "n_null_rows": int(len(finite)),
               "quantile": float(quantile), "seed": int(seed), "min_lag": int(min_lag),
               "band": [int(nrl_min), int(nrl_max)],
               "null_median": float(np.median(finite)), "null_q90": float(np.quantile(finite, 0.9)),
               "null_q95": float(np.quantile(finite, 0.95)), "null_q99": float(np.quantile(finite, 0.99)),
               "min_prominence": threshold}
    return threshold, details, table


def shuffle_pool(binary, allowed, mode, acf_fn, n_shuffles, seed, chunk=2000):
    """Parent ACFs of n_shuffles seeded shuffles of every molecule: (n_molecules, n_shuffles, n_lags)."""
    rng = np.random.default_rng(seed)
    pool = np.empty((len(binary), n_shuffles, binary.shape[1]), dtype=np.float64)
    pool_valid = np.empty((len(binary), n_shuffles), dtype=bool)
    for s in range(n_shuffles):
        for start in range(0, len(binary), chunk):
            block = binary[start: start + chunk]
            shuffled = shuffle_rows(block, allowed[start: start + chunk], rng)
            check_shuffle(block, shuffled, allowed[start: start + chunk], mode)
            acf, valid = acf_fn(shuffled)
            pool[start: start + len(block), s] = acf
            pool_valid[start: start + len(block), s] = valid
    return pool, pool_valid


def group_null(pool, pool_valid, members, n_repeats, min_lag, nrl_min, nrl_max, rng):
    """Null peak prominences of n_repeats group-sized means drawn (with replacement) from members' shuffles."""
    members = np.asarray(members)
    if not len(members):
        return np.full(n_repeats, np.nan)
    if not pool_valid[members].all():
        raise ValueError("Group null members must have valid shuffled ACFs")
    null = np.empty(n_repeats)
    for r in range(n_repeats):
        drawn = members[rng.integers(0, len(members), len(members))]
        if not np.isin(drawn, members).all():
            raise ValueError("Group null drew a molecule outside the group")
        shuffles = rng.integers(0, pool.shape[1], len(members))
        null[r] = band_peak(pool[drawn, shuffles].mean(axis=0), min_lag, nrl_min, nrl_max)[2]
    return null


def empirical_p(observed, null):
    """(1 + #null >= observed) / (1 + n_null); NaN when either side is undefined."""
    null = np.asarray(null)
    null = null[np.isfinite(null)]
    if not np.isfinite(observed) or not len(null):
        return np.nan
    return float((1 + np.sum(null >= observed)) / (1 + len(null)))


def summarize_nrl(table, group_column, groups, tests):
    rows = []
    for group in groups:
        sub = table[table[group_column].eq(group)]
        peak = sub.nrl_bp.notna()
        nrl = sub.nrl_bp.dropna()
        row = {group_column: group, "n": len(sub), "n_acf_valid": int(sub.nrl_status.ne("zero_variance").sum()),
               "n_zero_variance": int(sub.nrl_status.eq("zero_variance").sum()),
               "n_unclustered": int(sub.cluster.eq("Unclustered").sum()),
               "n_with_nrl_peak": int(peak.sum()), "fraction_with_nrl_peak": peak.mean() if len(sub) else np.nan,
               "median_nrl_bp": nrl.median() if len(nrl) else np.nan,
               "iqr_low_nrl_bp": nrl.quantile(0.25) if len(nrl) else np.nan,
               "iqr_high_nrl_bp": nrl.quantile(0.75) if len(nrl) else np.nan,
               "mad_nrl_bp": float(median_abs_deviation(sub.nrl_bp.to_numpy(dtype=float), nan_policy="omit"))
               if len(nrl) else np.nan,
               "median_peak_ratio": sub.peak_ratio.median() if sub.peak_ratio.notna().any() else np.nan}
        test = tests[(tests.group_type == group_column) & (tests.group == group)]
        for column in GROUP_COLUMNS:
            row[column] = float(test[column].iloc[0]) if len(test) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def group_tests(profiles, valid, metadata, groupings, pool, pool_valid, n_repeats, params, seed):
    """Observed in-band peak of each group's mean ACF and its empirical p-value against group-sized null means."""
    rng = np.random.default_rng(seed)
    rows = []
    for group_type, groups in groupings:
        for group in groups:
            members = np.flatnonzero(metadata[group_type].eq(group).to_numpy() & valid)
            if len(members):
                lag, height, prominence = band_peak(profiles[members].mean(axis=0), **params)
                null = group_null(pool, pool_valid, members, n_repeats, rng=rng, **params)
            else:
                lag = height = prominence = np.nan
                null = np.full(n_repeats, np.nan)
            finite = null[np.isfinite(null)]
            rows.append({"group_type": group_type, "group": group, "n_valid": len(members),
                         "group_peak_lag_bp": lag, "group_peak_height": height,
                         "group_peak_prominence": prominence, "n_null_groups": int(len(finite)),
                         "n_null_ge_observed": int(np.sum(finite >= prominence)) if np.isfinite(prominence) else np.nan,
                         "group_peak_pvalue": empirical_p(prominence, null),
                         "null_median_prominence": float(np.median(finite)) if len(finite) else np.nan,
                         "null_q95_prominence": float(np.quantile(finite, 0.95)) if len(finite) else np.nan})
            log(f"nrl: group {group_type}={group}: n_valid={len(members)}, peak lag={lag}, "
                f"prominence={prominence}, p={rows[-1]['group_peak_pvalue']}")
    return pd.DataFrame(rows)


def stage_nrl(args, out):
    binary_path, rows_path = out / "intermediate/binary_m6a.npy", out / "intermediate/row_ids.npy"
    metadata_path = out / "tables/clustered_molecules.tsv"
    params = dict(min_lag=args.min_lag, nrl_min=args.nrl_min, nrl_max=args.nrl_max)
    signature = fingerprint([binary_path, rows_path, out / "intermediate/acf.npy", out / "intermediate/acf_valid.npy",
                             metadata_path, Path(__file__), PARENT_SUMMARY, PARENT_ACF]
                            + ([str(args.ref) + ".fai"] if args.null_shuffle == "at" else []),
                            {**params, "window_name": args.window_name, "min_prominence": args.min_prominence,
                             "prominence_quantile": args.prominence_quantile, "n_null": args.n_null,
                             "null_shuffle": args.null_shuffle, "ref": str(args.ref.resolve()),
                             "n_shuffles_per_molecule": args.n_shuffles_per_molecule,
                             "n_null_groups": args.n_null_groups, "seed": args.seed})
    binary, row_ids = load_matrix(out)
    metadata = read_metadata(metadata_path)
    assert_alignment(binary, metadata, row_ids)
    width = window_width(metadata)
    profiles, valid = load_acf(out, metadata)
    if profiles.shape[1] != width:
        raise ValueError(f"Expected all {width} lags of the window; the ACF has {profiles.shape[1]}")
    if metadata.orientation.ne("transcriptional").any():
        raise ValueError("The null requires the strand-oriented (transcriptional) binary matrix")
    window_name = args.window_name or out.name
    binary = np.ascontiguousarray(binary)
    acf_parent = load_parent(PARENT_ACF)
    acf_fn = lambda rows: acf_parent.autocorrelations(rows, n_features=width)

    # Null positions: A/T bases of each molecule's own oriented reference window (default) or all.
    at = None
    calls_off_at = np.nan
    if args.null_shuffle == "at":
        with pysam.FastaFile(str(args.ref)) as fasta:
            at = at_mask(fasta, metadata, width)
        calls_off_at = float((binary.astype(bool) & ~at).sum() / max(1, binary.sum()))
        log(f"nrl: A/T mask fetched; {calls_off_at:.4%} of observed m6A calls lie on reference G/C")
    allowed = allowed_positions(binary, args.null_shuffle, at)

    # Threshold: calibrated on shuffled copies of this region's own molecules unless overridden.
    if args.min_prominence is not None:
        min_prominence = float(args.min_prominence)
        calibration = {"method": "fixed override (--min-prominence)", "shuffle_mode": args.null_shuffle,
                       "window_width": width, "min_prominence": min_prominence}
        null_table = pd.DataFrame(columns=["row_index", "acf_valid", "null_max_prominence"])
    else:
        min_prominence, calibration, null_table = calibrate_prominence(
            binary, allowed, args.null_shuffle, acf_fn, args.n_null, args.prominence_quantile,
            seed=args.seed, **params)
    log(f"nrl: min_prominence = {min_prominence:.4f} ({calibration['method']})")
    log(f"nrl: {len(metadata):,} molecules, window {width} bp, lags 0-{width - 1} at 1 bp; "
        f"NRL band {args.nrl_min}-{args.nrl_max}, min_lag {args.min_lag}, prominence {min_prominence:.4f}")

    parent = load_parent(PARENT_SUMMARY)
    records = []
    for i in range(len(metadata)):
        metrics = molecule_metrics(np.asarray(profiles[i]), min_prominence, repeat_peak=parent.repeat_peak, **params)
        records.append({column: metrics[column] for column in NRL_METRICS})
        if (i + 1) % 2000 == 0:
            log(f"nrl: {i + 1:,} molecules done")
    table = metadata[[c for c in NRL_META if c != "window"]].copy()
    table.insert(NRL_META.index("window"), "window", window_name)
    table = pd.concat([table.reset_index(drop=True), pd.DataFrame(records)], axis=1)[NRL_META + NRL_METRICS]
    if not table.nrl_status.isin(NRL_STATUSES).all():
        raise ValueError("Unknown NRL status")
    if not np.array_equal(table.nrl_status.eq("zero_variance").to_numpy(), ~valid):
        raise ValueError("zero_variance status must coincide with the invalid ACF rows")
    if not (table.nrl_bp.notna() == table.nrl_status.eq("ok")).all():
        raise ValueError("nrl_bp must be present exactly for status ok")
    if not table.nrl_bp.dropna().between(args.nrl_min, args.nrl_max).all():
        raise ValueError("NRL outside the configured band")

    # Group-level null: every molecule shuffled n times; group-sized means drawn from the group's own shuffles.
    log(f"nrl: shuffling every molecule {args.n_shuffles_per_molecule} times for the group null")
    pool, pool_valid = shuffle_pool(binary, allowed, args.null_shuffle, acf_fn, args.n_shuffles_per_molecule, args.seed + 1)
    if pool.shape != (len(metadata), args.n_shuffles_per_molecule, width):
        raise ValueError("Shuffle pool shape differs from (molecules, shuffles, lags)")
    if not np.array_equal(pool_valid.all(axis=1), valid) or not np.array_equal(pool_valid.any(axis=1), valid):
        raise ValueError("Shuffled validity must coincide with the original nonzero-variance molecules")
    clusters = sorted(table.cluster.unique(), key=lambda c: (c == "Unclustered", int(c) if c != "Unclustered" else 0))
    tests = group_tests(profiles, valid, metadata, [("expr_bin", list(BINS)), ("cluster", clusters)],
                        pool, pool_valid, args.n_null_groups, params, args.seed + 2)
    del pool

    outputs = [out / "tables/nrl_per_molecule.tsv", out / "tables/nrl_summary_by_bin.tsv",
               out / "tables/nrl_summary_by_cluster.tsv", out / "tables/nrl_prominence_calibration.tsv",
               out / "tables/nrl_null_prominences.tsv", out / "tables/group_peak_tests.tsv",
               out / "validation/nrl_validation.json"]
    table.to_csv(outputs[0], sep="\t", index=False, na_rep="NA", float_format="%.6g")
    for path, column, groups in ((outputs[1], "expr_bin", list(BINS)), (outputs[2], "cluster", clusters)):
        summary = summarize_nrl(table, column, groups, tests)
        summary["min_prominence"] = min_prominence
        summary.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.6g")
    pd.DataFrame([{"window": window_name, "acf_length_bp": int(profiles.shape[1]), **calibration}]).to_csv(
        outputs[3], sep="\t", index=False, na_rep="NA", float_format="%.6g")
    null_table = null_table.copy()
    picked = null_table.row_index.to_numpy(dtype=int)
    null_table.insert(1, "read_id", metadata.read_id.to_numpy()[picked] if len(null_table) else [])
    null_table.insert(2, "expr_bin", metadata.expr_bin.to_numpy()[picked] if len(null_table) else [])
    null_table["min_prominence"] = min_prominence
    null_table.to_csv(outputs[4], sep="\t", index=False, na_rep="NA", float_format="%.6g")
    tests.to_csv(outputs[5], sep="\t", index=False, na_rep="NA", float_format="%.6g")
    report = {"n_molecules": len(table), "window": window_name, "window_width_bp": width,
              "lags": [0, int(profiles.shape[1] - 1)], "lag_resolution_bp": 1,
              "status_counts": table.nrl_status.value_counts().to_dict(),
              "fraction_with_nrl_peak": float(table.nrl_bp.notna().mean()),
              "median_nrl_bp": float(table.nrl_bp.median()) if table.nrl_bp.notna().any() else None,
              "min_prominence": min_prominence, "prominence_calibration": calibration,
              "null_shuffle": args.null_shuffle, "reference_fasta": str(args.ref.resolve()),
              "fraction_observed_calls_on_reference_gc": calls_off_at,
              "n_shuffles_per_molecule": args.n_shuffles_per_molecule, "n_null_groups": args.n_null_groups,
              "shuffle_checks": {"m6a_count_preserved": True, "calls_only_at_allowed_positions": True,
                                 "group_null_draws_only_group_molecules": True},
              "parameters": {**params, "min_prominence": min_prominence, "seed": args.seed},
              "acf": "unchanged parent autocorrelations(); all lags of the window at 1 bp; no smoothing",
              "nrl_reference": "Abdulhay et al. 2020 eLife (secondary-peak scan)"}
    write_json(outputs[6], report)
    finish_stage(out, "nrl", signature, outputs, details=report)
    log(f"nrl: finished; {report['status_counts']}")


# ---------------------------------------------------------------------------
# Command line: all stages in order on one region
# ---------------------------------------------------------------------------

def parse_arguments():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    g = ap.add_argument_group("sample")
    g.add_argument("--bins-tsv", type=Path, default=EXPRESSION)
    g.add_argument("--canonical-bed", type=Path, default=Path(CANONICAL))
    g.add_argument("--ft-root", type=Path, default=FT_ROOT)
    g.add_argument("--per-bin", type=int, default=2500)
    g.add_argument("--window-start", type=int, default=-1000,
                   help="Window start as a TSS-relative transcriptional offset (inclusive; negative = upstream).")
    g.add_argument("--window-end", type=int, default=1000,
                   help="Window end as a TSS-relative transcriptional offset (exclusive; positive = downstream).")
    g.add_argument("--seed", type=int, default=0, help="Sampling, clustering and null seed.")
    g.add_argument("--timepoints", nargs="+", default=["LPS_0", "LPS_5", "LPS_10", "LPS_15"])
    g.add_argument("--chrom", nargs="+", default=CHROMS)
    g = ap.add_argument_group("acf")
    g.add_argument("--max-lag", type=int, default=None,
                   help="Inclusive maximum lag; default returns all nonnegative lags (window width - 1).")
    g = ap.add_argument_group("cluster")
    g.add_argument("--n-pcs", type=int, default=50)
    g.add_argument("--n-neighbors", type=int, default=10)
    g.add_argument("--resolution", type=float, default=0.4)
    g = ap.add_argument_group("footprints")
    g.add_argument("--extract-root", type=Path, default=EXTRACT_ROOT)
    g = ap.add_argument_group("nrl")
    g.add_argument("--window-name", default=None, help="Window label for the NRL table; default = folder name.")
    g.add_argument("--min-lag", type=int, default=60, help="Ignore lags below this (lag-0 shoulder).")
    g.add_argument("--nrl-min", type=int, default=120, help="Lowest NRL lag considered (bp).")
    g.add_argument("--nrl-max", type=int, default=300, help="Highest NRL lag considered (bp).")
    g.add_argument("--min-prominence", type=float, default=None,
                   help="Fixed find_peaks prominence override; default: calibrate from shuffled molecules.")
    g.add_argument("--prominence-quantile", type=float, default=0.95,
                   help="Quantile of the shuffled-molecule null maxima used as min_prominence.")
    g.add_argument("--n-null", type=int, default=2000, help="Maximum number of molecules shuffled for the calibration.")
    g.add_argument("--ref", type=Path, default=Path(REF_FA), help="Reference FASTA for the A/T-restricted shuffle.")
    g.add_argument("--null-shuffle", choices=SHUFFLE_MODES, default="at",
                   help="Shuffle calls among A/T positions of the oriented reference window (at) or all positions.")
    g.add_argument("--n-shuffles-per-molecule", type=int, default=5,
                   help="Shuffles of every sampled molecule for the group-level null.")
    g.add_argument("--n-null-groups", type=int, default=200,
                   help="Group-sized null means per group for the mean-ACF peak p-value.")
    args = ap.parse_args()
    if (args.per_bin < 1 or args.seed < 0 or len(set(args.timepoints)) != len(args.timepoints)
            or len(set(args.chrom)) != len(args.chrom) or not set(args.chrom).issubset(CHROMS)
            or not set(args.timepoints).issubset({"LPS_0", "LPS_5", "LPS_10", "LPS_15"})):
        ap.error("Require positive per-bin count, nonnegative seed, unique chr1-22/X/Y and LPS_0/5/10/15")
    if args.window_end - args.window_start < 3:
        ap.error("--window-end must exceed --window-start by at least three bases")
    if args.n_pcs < 2 or args.n_neighbors < 2 or not np.isfinite(args.resolution) or args.resolution <= 0:
        ap.error("Require >=2 PCs, >=2 neighbors and a finite positive resolution")
    if (args.min_lag < 1 or (args.min_prominence is not None and args.min_prominence < 0)
            or not 0.5 <= args.prominence_quantile < 1 or args.n_null < 10
            or args.nrl_min < 1 or args.nrl_max <= args.nrl_min
            or args.n_shuffles_per_molecule < 1 or args.n_null_groups < 10):
        ap.error("Require min_lag >= 1, prominence >= 0 if given, 0.5 <= prominence_quantile < 1, n_null >= 10, "
                 "1 <= nrl_min < nrl_max, n_shuffles_per_molecule >= 1 and n_null_groups >= 10")
    if args.null_shuffle == "at" and (not args.ref.is_file() or not Path(str(args.ref) + ".fai").is_file()):
        ap.error(f"Reference FASTA and .fai are required for the A/T-restricted shuffle: {args.ref}")
    return args


def main():
    args = parse_arguments()
    out = args.out_dir.resolve()
    prepare_dirs(out)
    for name, stage in (("sample", stage_sample), ("acf", stage_acf), ("cluster", stage_cluster),
                        ("tables", stage_tables), ("footprints", stage_footprints), ("nrl", stage_nrl)):
        log(f"==== stage {name} ====")
        stage(args, out)
    log("All processing stages finished: " + str(out))


if __name__ == "__main__":
    main()
