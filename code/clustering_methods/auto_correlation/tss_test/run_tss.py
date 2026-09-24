#!/usr/bin/env python3
"""Matched protein-coding canonical TSS workflow. All eligible TSSs by default; --smoke for a small test."""

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[3]
OUTPUT_ROOT = PROJECT / "LCL_project/auto_correlation/tss_test"
TSS_BED = Path("/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed")
FAI = Path("/project/spott/reference/human/GRCh38/hg38.fa.fai")
FT_ROOT = PROJECT / "macrophage_project/FiberHMM/extract/ft_result_dir"
SAMPLES = ["LPS_0", "LPS_5", "LPS_10", "LPS_15"]


def existing_module(filename):
    path = HERE.parent / filename
    spec = importlib.util.spec_from_file_location("tss_shared_" + path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextmanager
def runtime():
    """Only disposable library caches; no disk-backed analysis objects."""
    with tempfile.TemporaryDirectory(prefix="tss_autocorrelation_") as temporary:
        overrides = {key: "1" for key in (
            "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS")}
        overrides.update(NUMBA_CACHE_DIR=temporary + "/numba", MPLCONFIGDIR=temporary + "/mpl",
                         TMPDIR=temporary, PYTHONDONTWRITEBYTECODE="1")
        previous = {key: os.environ.get(key) for key in overrides}
        os.environ.update(overrides)
        try:
            yield
        finally:
            for key, value in previous.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value


class RBridge:
    def __enter__(self):
        env = os.environ.copy()
        rscript = env.get("AUTOCOR_RSCRIPT") or shutil.which("Rscript")
        if rscript is None:
            rscript = "/software/R-4.4.1-el8-x86_64/bin/Rscript"
            env.setdefault("R_LIBS_USER", str(Path.home() / "R/x86_64-pc-linux-gnu-library/4.4"))
            env["LD_LIBRARY_PATH"] = os.pathsep.join([
                "/software/openblas-0.3.29-el8-x86_64/lib",
                "/software/glpk-5.0-el8-x86_64/lib", env.get("LD_LIBRARY_PATH", "")])
        self.process = subprocess.Popen([rscript, "--vanilla", str(HERE / "tss_input.R"), str(PROJECT)],
                                        text=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, env=env)
        return self

    def request(self, **payload):
        self.process.stdin.write(json.dumps(payload) + "\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError("R input bridge failed; see its stderr above")
        return json.loads(line)

    def __exit__(self, kind, value, traceback):
        if self.process.poll() is None:
            self.process.stdin.close()
        try:
            code = self.process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            self.process.terminate()
            code = self.process.wait(timeout=10)
        self.process.stdout.close()
        if kind is None and code:
            raise RuntimeError(f"R bridge exited with status {code}")


def define_tss(bed_path=TSS_BED, fai_path=FAI):
    import numpy as np
    import pandas as pd
    bed = pd.read_csv(bed_path, sep="\t", header=None,
                      names=["chrom", "start", "end", "name", "score", "strand"])
    if not bed.end.sub(bed.start).eq(20).all() or not bed.strand.isin(["+", "-"]).all():
        raise ValueError("Expected canonical TSS +/-10-bp BED intervals with known strands")
    fields = bed.name.str.split(";", expand=True)
    bed[["gene_id", "tx_id", "gene_name", "gene_type"]] = fields.iloc[:, :4]
    if fields.shape[1] < 5 or not fields[4].str.contains("Ensembl_canonical").all():
        raise ValueError("Input contains noncanonical TSS records")
    bed["tss0"] = bed.start + 10
    bed["window_start0"] = bed.tss0 - 1000
    bed["window_end0"] = bed.tss0 + 1000
    bed["analysis_start"] = bed.window_start0 + 1
    bed["analysis_end"] = bed.window_end0
    bed["width"] = 2000
    # Collapse duplicate annotations of the same gene/TSS/strand, retaining
    # distinct TSSs of a gene and coincident TSSs belonging to different genes.
    bed = bed.drop_duplicates(["gene_id", "chrom", "tss0", "strand"]).copy()
    bed["region_id"] = bed.gene_id + ":" + bed.chrom + ":" + bed.tss0.astype(str) + ":" + bed.strand
    lengths = pd.read_csv(fai_path, sep="\t", header=None, usecols=[0, 1], index_col=0)[1]
    end = bed.chrom.map(lengths)
    bed["window_eligible"] = end.notna() & bed.window_start0.ge(0) & bed.window_end0.le(end)
    bed["window_status"] = np.where(bed.window_eligible, "eligible", "out_of_reference_bounds")
    return bed.reset_index(drop=True)


def expression_annotation(regions, mapping):
    import pandas as pd
    expression = pd.DataFrame(mapping).rename(columns={"tss": "tss0"})
    keys = ["gene_id", "chrom", "tss0", "strand"]
    expression = expression[keys + ["mean_tpm", "expr_bin"]].drop_duplicates()
    regions = regions.merge(expression, on=keys, how="left", validate="one_to_one")
    regions["expr_bin"] = regions.expr_bin.fillna("expression_unmapped")
    return regions


def select_analysis_regions(regions):
    """Filter BEFORE extraction: protein coding, matched expression, full window.

    Silent protein-coding genes are retained. The global RNA expression/bin
    definitions are not recomputed after filtering the canonical TSS cohort.
    """
    import numpy as np
    import pandas as pd
    from expression_summary import EXPRESSION_BINS

    coding = regions.gene_type.eq("protein_coding")
    tpm = pd.to_numeric(regions.mean_tpm, errors="coerce")
    matched = regions.expr_bin.isin(EXPRESSION_BINS) & np.isfinite(tpm) & tpm.ge(0)
    eligible = regions.window_eligible
    keep = coding & matched & eligible
    audit = dict(input_tss=len(regions), excluded_noncoding=int((~coding).sum()),
                 excluded_coding_unmatched=int((coding & ~matched).sum()),
                 excluded_coding_matched_out_of_bounds=int((coding & matched & ~eligible).sum()),
                 retained_tss=int(keep.sum()))
    return regions.loc[keep].copy().reset_index(drop=True), audit


def run_analysis(args):
    import numpy as np
    import pandas as pd
    from nrl_annotation import annotate, AnnotationParameters
    from expression_summary import summarize_expression
    from plots import make_figures

    acf_module = existing_module("03_compute_autocorrelations.py")
    cluster_module = existing_module("04_cluster_autocorrelations.py")
    shared_summary = existing_module("05_summarize_autocorrelations.py")
    regions = define_tss(args.tss_bed, args.fai)
    if not args.all_tss:
        wanted = args.smoke_genes
        if not 1 <= len(wanted) <= 3 or len(set(wanted)) != len(wanted):
            raise ValueError("Smoke test requires 1-3 distinct gene names")
        regions = regions[regions.gene_name.isin(wanted)].copy()
        if len(regions) != len(wanted) or set(regions.gene_name) != set(wanted):
            raise ValueError("Each smoke gene must identify exactly one canonical TSS")
        if not regions.window_eligible.all():
            raise ValueError("Smoke TSS outside reference bounds")
        print(f"SMOKE: exactly {len(regions)} TSSs; at most {args.smoke_reads} reads/sample/TSS", flush=True)
    chunks, audits = [], []
    with RBridge() as bridge:
        print("Loading the unchanged global expression-bin definitions in memory", flush=True)
        expr = bridge.request(mode="expression")
        regions = expression_annotation(regions, expr["mapping"])
        regions, tss_filter_audit = select_analysis_regions(regions)
        print("Pre-extraction TSS filter: " + str(tss_filter_audit), flush=True)
        if regions.empty:
            raise ValueError("No eligible protein-coding TSSs with matched expression annotations")
        if not args.all_tss and set(regions.gene_name) != set(args.smoke_genes):
            raise ValueError("Every smoke gene must be protein coding and have matched expression annotations")
        print(f"Expression quartiles: {expr['quartiles']}; {len(expr['rna_samples'])} RNA samples", flush=True)
        for region in regions.to_dict("records"):
            # Pass only extraction coordinates, avoiding NaN in the JSON request.
            request_region = {key: region[key] for key in ("region_id", "chrom", "analysis_start", "analysis_end")}
            result = bridge.request(mode="extract", region=request_region, samples=SAMPLES,
                                    ft_root=str(args.ft_root),
                                    max_reads=None if args.all_tss else args.smoke_reads, seed=args.seed)
            chunks.extend(result["records"])
            audits.extend(result["audit"])
            print(f"{region['gene_name']} {region['chrom']}:{region['window_start0']}-{region['window_end0']} "
                  f"({region['strand']}; BED): {len(result['records'])} retained", flush=True)
    if not chunks:
        raise ValueError("No fully spanning molecules at the requested TSSs")
    records = pd.DataFrame(chunks).rename(columns={"strand": "read_strand"})
    offsets = records.pop("call_offsets")
    records = records.merge(regions, on="region_id", how="left", validate="many_to_one", sort=False)
    if records.row_id.duplicated().any():
        raise AssertionError("Duplicate TSS/sample/read identifier")
    if not records.gene_type.eq("protein_coding").all() or records.mean_tpm.isna().any() \
            or records.expr_bin.eq("expression_unmapped").any():
        raise AssertionError("Noncoding or unmatched TSS reached the ACF input")
    binary = np.zeros((len(records), 2000), dtype=np.uint8)
    for i, positions in enumerate(offsets):
        positions = np.atleast_1d(positions).astype(int)
        if len(positions) and (positions.min() < 0 or positions.max() >= 2000):
            raise ValueError("m6A offset outside the exact 2-kb window")
        binary[i, positions] = 1
    assert (records.read_start <= records.analysis_start).all()
    assert (records.read_end >= records.analysis_end).all()
    records["m6a_calls"] = binary.sum(axis=1)
    records["m6a_call_fraction"] = records.m6a_calls / 2000
    print("Computing unchanged raw binary ACF and pooled PCA/correlation-kNN/weighted Leiden", flush=True)
    profiles, valid = acf_module.autocorrelations(binary)  # all lags 0..1999
    binary.flags.writeable = False
    profiles.flags.writeable = False
    records["acf_valid"] = valid
    labels, status, embedding, info, edges = cluster_module.cluster_profiles(
        profiles, valid, n_pcs=50, n_neighbors=10, resolution=0.4, seed=args.seed)
    records["cluster"], records["status"] = labels, status
    records["umap1"], records["umap2"] = embedding[:, 0], embedding[:, 1]
    records["raw_peak_140_250_bp"] = [shared_summary.repeat_peak(p)[0] for p in profiles]
    # A single pooled cluster vocabulary is necessary for expression-bin
    # composition. The original per-region summary is given that pool's ID.
    cluster_summary = shared_summary.summarize(profiles, records.assign(region_id="all_tss_pool"))
    fingerprint = hashlib.sha256(profiles).hexdigest()
    print("Annotating each molecule using the separate 33-bp branch", flush=True)
    annotation = annotate(binary)
    records = pd.concat([records, annotation], axis=1)
    assert hashlib.sha256(profiles).hexdigest() == fingerprint
    summaries = summarize_expression(records, regions, profiles)
    print(f"Clustered {info['n_clustered']}/{len(records)}; "
          f"reliable NRL {int(records.reliable_nrl.sum())}/{len(records)}", flush=True)
    if not args.all_tss:
        assert records.region_id.nunique() == len(regions), "Smoke TSS lacks spanning reads"
        assert (records.groupby("region_id").acf_valid.sum() >= 4).all()
        assert info["n_clustered"] >= 4 and len(edges) > 0, "Clustering did not execute"
        assert profiles.shape[1] == 2000 and np.allclose(profiles[valid, 0], 1)
        assert records.mean_tpm.notna().all(), "Smoke genes must exercise expression mapping"
        assert records.n_extrema.gt(0).any(), "Smoke reads did not exercise peak detection"
        assert records.nrl_bp.notna().any(), "Smoke reads did not exercise period regression"
    # Both modes run the identical plotting code. Smoke PDFs remain in BytesIO.
    pages = make_figures(records, regions, binary, profiles, summaries, cluster_summary,
                         None if not args.all_tss else args.pdf_dir)
    print(f"PASS: TSS definition, extraction, raw ACF, PCA/kNN/Leiden, smoothing, extrema, "
          f"NRL/regularity, expression summaries, and {pages} PDF figures", flush=True)
    print("NRL status counts: " + str(records.annotation_status.value_counts().to_dict()), flush=True)
    return dict(regions=regions, records=records, binary=binary, raw_acf=profiles,
                summaries=summaries, cluster_summary=cluster_summary, cluster_info=info,
                edges=edges, audit=pd.DataFrame(audits), expression=expr,
                tss_filter_audit=tss_filter_audit,
                annotation_parameters=AnnotationParameters(), plot_count=pages)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--smoke", dest="all_tss", action="store_false",
                      help="Test only 1-3 TSSs; render PDFs in memory")
    mode.add_argument("--all-tss", dest="all_tss", action="store_true",
                      help="Run all expression-matched protein-coding canonical TSSs (default; potentially large RAM usage)")
    parser.set_defaults(all_tss=True)
    parser.add_argument("--smoke-genes", nargs="+", default=["ISG15", "ACTB", "HBB"])
    parser.add_argument("--smoke-reads", type=int, default=8, help="Per TSS/sample, hard maximum 16")
    parser.add_argument("--tss-bed", type=Path, default=TSS_BED)
    parser.add_argument("--fai", type=Path, default=FAI)
    parser.add_argument("--ft-root", type=Path, default=FT_ROOT)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--pdf-dir", type=Path,
                        help=f"New PDF directory under {OUTPUT_ROOT}; relative paths use this root "
                             "(default: figures/run_<UTC timestamp>_<PID>)")
    args = parser.parse_args(argv)
    if not 1 <= args.smoke_reads <= 16:
        parser.error("--smoke-reads must be 1..16")
    if args.all_tss:
        if args.pdf_dir is None:
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
            args.pdf_dir = OUTPUT_ROOT / "figures" / f"run_{stamp}_{os.getpid()}"
        elif not args.pdf_dir.is_absolute():
            args.pdf_dir = OUTPUT_ROOT / args.pdf_dir
        args.pdf_dir = args.pdf_dir.resolve()
        if not args.pdf_dir.is_relative_to(OUTPUT_ROOT.resolve()):
            parser.error(f"--pdf-dir must be inside {OUTPUT_ROOT}")
        if args.pdf_dir.exists():
            parser.error("--pdf-dir must not already exist (protect existing outputs)")
    elif args.pdf_dir is not None:
        parser.error("Smoke plots remain in memory; omit --pdf-dir")
    return args


if __name__ == "__main__":
    arguments = parse_args()
    if arguments.all_tss:
        print(f"ALL TSS: matched protein-coding cohort; final PDF directory = {arguments.pdf_dir}", flush=True)
    with runtime():
        run_analysis(arguments)
