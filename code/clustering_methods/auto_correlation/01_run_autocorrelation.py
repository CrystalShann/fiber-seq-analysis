#!/usr/bin/env python3
"""Run LCL autocorrelation with configurable SNP-centered windows and features."""

import argparse
import hashlib
from importlib import import_module
import io
import json
import os
from pathlib import Path
import shutil
import subprocess


for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[variable] = "1"
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/fiberseq-autocorrelation-numba")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/fiberseq-autocorrelation-matplotlib")

import numpy as np
import pandas as pd

acf = import_module("03_compute_autocorrelations")
clustering = import_module("04_cluster_autocorrelations")
summaries = import_module("05_summarize_autocorrelations")
repeat_plots = import_module("09_haplotype_repeat_plots")
PROJECT = Path(__file__).resolve().parents[3]
FUNCTIONS = PROJECT / "code/topic_model/topic_modelling_functions.r"
PHASED_INPUT = PROJECT / "code/haplotype_phasing/LCL_phased_m6a_input.r"
R_BUILDER = Path(__file__).with_name("02_build_m6a_input.r")
SIGNATURE = PROJECT / "LCL_project/Leiden_manhattan/summary tables/selection_signature.rds"
FT_ROOT = Path("/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results")
OUTPUT_ROOT = PROJECT / "LCL_project/auto_correlation"
DEFAULT_OUTPUT = OUTPUT_ROOT / "resolution04"
RECOVERED_SIGNATURE = OUTPUT_ROOT / "resolution1/inputs/selection_signature.rds"
SAMPLE_TABLE = Path("/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv")
R_PLOTS = Path(__file__).with_name("07_plot_fiberseq.R")


def input_directory(output, region_id=None):
    path = output / "inputs"
    return path if region_id is None else path / region_id


def result_directory(output, region_id=None):
    path = output / "outputs"
    return path if region_id is None else path / region_id


def r_environment():
    """Use an activated R environment, or the existing RCC R 4.4 installation."""
    env = os.environ.copy()
    rscript = env.get("AUTOCOR_RSCRIPT") or shutil.which("Rscript")
    if rscript is None:
        rscript = "/software/R-4.4.1-el8-x86_64/bin/Rscript"
        user_lib = Path.home() / "R/x86_64-pc-linux-gnu-library/4.4"
        if user_lib.is_dir():
            env.setdefault("R_LIBS_USER", str(user_lib))
        libraries = ["/software/openblas-0.3.29-el8-x86_64/lib",
                     "/software/glpk-5.0-el8-x86_64/lib"]
        env["LD_LIBRARY_PATH"] = os.pathsep.join(libraries + [env.get("LD_LIBRARY_PATH", "")])
    if not Path(rscript).is_file():
        raise FileNotFoundError("Activate R with dplyr, Matrix, GenomicRanges and Rsamtools, or set AUTOCOR_RSCRIPT")
    return rscript, env


def run_m6a_builder(*args):
    """Receive the R builder's output through a pipe, without temporary matrices."""
    rscript, env = r_environment()
    result = subprocess.run(
        [rscript, "--vanilla", str(R_BUILDER), *map(str, args)],
        env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if result.returncode:
        raise RuntimeError(f"m6A input construction failed:\n{result.stderr}")
    return result


def load_regions(window_size=2000, selection_signature=SIGNATURE):
    """Let the R input builder validate and center the selected LCL regions."""
    result = run_m6a_builder("--regions", selection_signature, window_size)
    return pd.read_csv(io.StringIO(result.stdout), sep="\t")


def load_m6a_binary(region, selection_signature=SIGNATURE):
    command = [FUNCTIONS, selection_signature, region.region_id, FT_ROOT, region.width]
    result = run_m6a_builder(*command)
    records = pd.read_csv(io.StringIO(result.stdout), sep="\t", keep_default_na=False,
                          dtype={"row_id": str, "sample_name": str, "RID": str,
                                 "strand": str, "call_offsets": str})
    matrix = np.zeros((len(records), region.width), dtype=np.uint8)
    for i, offsets in enumerate(records.pop("call_offsets")):
        if offsets:
            indices = np.array([int(x) for x in offsets.split(",")])
            if (indices < 0).any() or (indices >= region.width).any():
                raise ValueError("m6A input contains a position outside the window")
            matrix[i, indices] = 1
    if records.row_id.duplicated().any():
        raise ValueError("Duplicate sample/read identifiers")
    if not (records.haplotype.isin(["HP1", "HP2"]).all()
            and records.focal_genotype.isin(["0|1", "1|0"]).all()
            and records.allele_status.eq("phased_focal_genotype").all()):
        raise ValueError("LCL autocorrelation requires resolved phased heterozygous reads")
    records.insert(0, "region_id", region.region_id)
    records["m6a_calls"] = matrix.sum(axis=1)
    records["m6a_call_fraction"] = matrix.mean(axis=1)
    return matrix, records, result.stderr


def plot_fiberseq(output, region, binary, profiles, records, edges, args, info):
    chosen = records.status.eq("clustered").to_numpy()
    if not chosen.any():
        return 0
    ids = records.row_id.to_numpy()
    payload = dict(region=json.loads(pd.DataFrame([region._asdict()]).to_json(orient="records")),
                   records=json.loads(records.to_json(orient="records", double_precision=15)),
                   row_ids=ids[chosen].tolist(),
                   call_offsets=[np.flatnonzero(row).tolist() for row in binary[chosen]],
                   acf=profiles[chosen].tolist(),
                   edges=dict(source=ids[edges[:, 0].astype(int)].tolist(),
                              target=ids[edges[:, 1].astype(int)].tolist(),
                              weight=edges[:, 2].tolist()))
    rscript, env = r_environment()
    result = subprocess.run([rscript, "--vanilla", str(R_PLOTS), str(PROJECT),
        str(output), region.region_id, str(args.seed), str(args.resolution),
        str(info["actual_neighbors"])], input=json.dumps(payload, allow_nan=False),
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    (result_directory(output, region.region_id) / "plotting.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"Fiber-seq plotting failed for {region.region_id}:\n{result.stderr}")
    prefix = "AUTOCOR_HEATMAP_READS="
    counts = [line[len(prefix):] for line in result.stdout.splitlines() if line.startswith(prefix)]
    if len(counts) != 1:
        raise RuntimeError(f"Missing plotting completion count for {region.region_id}")
    return int(counts[0])


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_allele_results(output, region, profiles, records):
    """Compute allele tables in memory and save only their final PDF plots."""
    tables = summaries.summarize_alleles(profiles, records, region)
    summaries.plot_alleles(result_directory(output, region.region_id), region, tables)
    return tables


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region", action="append", help="Exact saved region_id; repeat to select several (default: all selected regions)")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT,
                        help="Run root containing inputs/ and outputs/ (default: auto_correlation/resolution04)")
    parser.add_argument("--window-size", type=int, default=2000, help="SNP-centered genomic width in bp (default: 2000)")
    parser.add_argument("--n-features", type=int, default=None,
                        help="Retain the first N ACF lags, including lag 0 (default: all window-size lags)")
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-neighbors", type=int, default=10)
    parser.add_argument("--resolution", type=float, default=0.4,
                        help="Leiden resolution (default: 0.4); set explicitly to override")
    parser.add_argument("--selection-signature", type=Path, default=SIGNATURE if SIGNATURE.is_file() else RECOVERED_SIGNATURE,
                        help="Saved top-10 selection RDS containing regions and sample metadata")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--sample-table", type=Path, default=SAMPLE_TABLE, help="CSV mapping sample_name to biological cell_line")
    args = parser.parse_args(argv)
    if args.window_size < 3 or args.n_pcs < 2 or args.n_neighbors < 2 or not np.isfinite(args.resolution) or args.resolution <= 0:
        parser.error("Require window-size >=3, n-pcs/n-neighbors >=2 and finite positive resolution")
    n_features = args.window_size if args.n_features is None else args.n_features
    if not 3 <= n_features <= args.window_size:
        parser.error("n-features must be between 3 and window-size, inclusive")
    sample_table = pd.read_csv(args.sample_table, dtype=str, keep_default_na=False)
    output = args.output_dir.resolve()
    if output != OUTPUT_ROOT.resolve() and OUTPUT_ROOT.resolve() not in output.parents:
        parser.error(f"Output must be inside {OUTPUT_ROOT}")
    if output.exists() and any(output.iterdir()):
        parser.error("Output directory is not empty; choose a new --output-dir subdirectory")
    regions = load_regions(args.window_size, args.selection_signature)
    if args.region:
        unknown = set(args.region) - set(regions.region_id)
        if unknown:
            parser.error(f"Unknown region IDs: {sorted(unknown)}")
        regions = regions[regions.region_id.isin(args.region)]
    for subdir in ("inputs", "outputs"):
        (output / subdir).mkdir(parents=True, exist_ok=True)
    # Keep region input folders for extraction logs; source data stays at its original path.
    for region in regions.itertuples(index=False):
        input_directory(output, region.region_id).mkdir()
        result_directory(output, region.region_id).mkdir()
    sources = [FUNCTIONS, PHASED_INPUT, R_BUILDER, args.selection_signature, args.sample_table, R_PLOTS,
               PROJECT / "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r",
               PROJECT / "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r",
               PROJECT / "code/haplotype_phasing/LCL_phasing.r",
               *sorted(Path(__file__).parent.glob("*.py"))]
    source_hashes = {str(p): sha256(p) for p in sources}
    print(f"Resolution={args.resolution}; window={args.window_size}; lags=0–{n_features - 1}; "
          f"PCs={args.n_pcs}; neighbors={args.n_neighbors}; seed={args.seed}. "
          "Tables and matrices remain in memory; saving PDFs and diagnostic logs only.", flush=True)
    locus_rows = []
    for region in regions.itertuples(index=False):
        rid = region.region_id
        print(f"{rid}: building m6A input from original BEDs in memory", flush=True)
        binary, records, extraction_log = load_m6a_binary(region, args.selection_signature)
        (input_directory(output, rid) / "extraction.log").write_text(extraction_log)
        profiles, valid = acf.autocorrelations(binary, args.n_features)
        records = summaries.annotate_alleles(records, region, sample_table)
        records["acf_valid"] = valid
        peaks = [summaries.repeat_peak(p) if v else (np.nan, np.nan) for p, v in zip(profiles, valid)]
        records["positive_local_peak_140_250_bp"] = [p[0] for p in peaks]
        records["peak_acf"] = [p[1] for p in peaks]
        records["peak_band_complete"] = n_features >= 252
        print(f"{rid}: {len(records)} fully spanning reads; clustering {valid.sum()} valid ACFs", flush=True)
        labels, status, embedding, info, edges = clustering.cluster_profiles(profiles, valid, n_pcs=args.n_pcs,
            n_neighbors=args.n_neighbors, resolution=args.resolution, seed=args.seed)
        records["cluster"], records["status"] = labels, status
        records["umap1"], records["umap2"] = embedding[:, 0], embedding[:, 1]
        avg, stats, counts = summaries.summarize(profiles, records)
        summaries.plot_region(result_directory(output, rid), region, binary, profiles, records, avg, counts)
        allele_tables = write_allele_results(output, region, profiles, records)
        repeat_rows = repeat_plots.summarize_repeats(profiles, records, region, allele_tables)
        repeat_plots.plot_haplotype_repeats(result_directory(output, rid), region, profiles, repeat_rows)
        locus_rows.extend(repeat_plots.compact_locus_rows(repeat_rows))
        n_heatmap = plot_fiberseq(output, region, binary, profiles, records, edges, args, info)
        print(f"{rid}: wrote {len(stats)} clusters' figures; {n_heatmap} heatmap reads", flush=True)
    repeat_plots.plot_locus_heatmap(result_directory(output), locus_rows)
    if any(sha256(p) != digest for p, digest in source_hashes.items()):
        raise RuntimeError("An input/code source changed during the run")
    print(f"Complete: {output}", flush=True)


if __name__ == "__main__":
    main()
