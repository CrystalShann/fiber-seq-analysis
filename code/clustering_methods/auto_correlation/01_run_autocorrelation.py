#!/usr/bin/env python3
"""Run LCL autocorrelation with configurable SNP-centered windows and features."""

import argparse
from functools import lru_cache
import hashlib
from importlib import import_module
import inspect
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[variable] = "1"
os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(tempfile.gettempdir()) / "fiberseq-autocorrelation-numba"))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "fiberseq-autocorrelation-matplotlib"))

import numpy as np
import pandas as pd

acf = import_module("03_compute_autocorrelations")
clustering = import_module("04_cluster_autocorrelations")
summaries = import_module("05_summarize_autocorrelations")
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
R_REPORT = Path(__file__).with_name("06_autocorrelation.Rmd")


def input_directory(output, region_id=None):
    path = output / "inputs"
    return path if region_id is None else path / region_id


def result_directory(output, region_id=None):
    path = output / "outputs"
    return path if region_id is None else path / region_id


@lru_cache(maxsize=1)
def r_environment():
    """Choose and check one R environment for input and report preparation."""
    rcc_rscript = Path("/software/R-4.4.1-el8-x86_64/bin/Rscript")
    requested = os.environ.get("AUTOCOR_RSCRIPT")
    preferred = requested or shutil.which("Rscript")
    candidates = [preferred] if preferred else []
    if not requested:
        candidates.append(str(rcc_rscript))
    probe = ('packages <- c("dplyr", "Matrix", "GenomicRanges", "Rsamtools", '
             '"IRanges", "jsonlite", "data.table", "igraph"); '
             'failures <- vapply(packages, function(package) tryCatch({ '
             'loadNamespace(package); "" }, error=function(e) '
             'paste0(package, ": ", conditionMessage(e))), character(1)); '
             'if (any(nzchar(failures))) stop(paste(failures[nzchar(failures)], '
             'collapse="\\n"), call.=FALSE); cat(R.version.string)')
    failures, checked = [], set()
    for candidate in candidates:
        executable = shutil.which(candidate)
        if executable is None:
            failures.append(f"{candidate}: executable not found")
            continue
        rscript = str(Path(executable).resolve())
        if rscript in checked:
            continue
        checked.add(rscript)
        env = os.environ.copy()
        env.setdefault("R_LIBS_USER", "/project/spott/cshan/Rlibs/%p-library/%v")
        if Path(rscript) == rcc_rscript.resolve():
            env.pop("R_HOME", None)
            user_lib = Path("/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4")
            if user_lib.is_dir():
                env["R_LIBS_USER"] = os.pathsep.join(dict.fromkeys(
                    [str(user_lib), env["R_LIBS_USER"]]))
            libraries = ["/software/gcc-12.2.0-el8-x86_64/lib64",
                         "/software/openblas-0.3.29-el8-x86_64/lib",
                         "/software/glpk-5.0-el8-x86_64/lib"]
            env["LD_LIBRARY_PATH"] = os.pathsep.join(
                libraries + [env.get("LD_LIBRARY_PATH", "")])
        try:
            result = subprocess.run([rscript, "--vanilla", "-e", probe], env=env,
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
        except (OSError, subprocess.TimeoutExpired) as error:
            failures.append(f"{rscript}: {error}")
            continue
        if result.returncode:
            failures.append(f"{rscript}:\n{result.stderr.strip()}")
            continue
        if failures:
            print("R on PATH could not load the required packages; using the RCC installation.", flush=True)
        print(f"Using {rscript} ({result.stdout.strip()})", flush=True)
        return rscript, env
    raise RuntimeError("No usable R environment for autocorrelation preparation.\n"
                       + "\n".join(failures)
                       + "\nSet AUTOCOR_RSCRIPT to an Rscript with the required packages "
                         "and R_LIBS_USER to its matching package library.")


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


def m6a_input_provenance(region, selection_signature):
    """Track extraction inputs without hashing large BED or phasing files."""
    sources = [R_BUILDER, FUNCTIONS, PHASED_INPUT,
               PROJECT / "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r",
               PROJECT / "code/haplotype_phasing/LCL_phasing.r"]
    phase_root = Path("/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples")
    phase_cache = PROJECT / "LCL_project/Leiden_manhattan/phase_summary"
    inputs = []
    for sample in str(region.contributing_samples).split(","):
        bed = FT_ROOT / sample / "extracted_results/m6a_by_chr" / f"{sample}.ft_extracted_m6a.{region.chr}.bed.gz"
        vcf = phase_root / sample / f"{sample}.5mC.6mA.aligned.phased.vcf.gz"
        inputs.extend([bed, Path(str(bed) + ".tbi"), vcf, Path(str(vcf) + ".tbi"),
                       phase_root / sample / "read-level-phasing.tsv",
                       phase_root / sample / "blocks.tsv", phase_cache / f"{sample}_haplotags.rds"])
    file_stats = {}
    for path in inputs:
        try:
            stat = path.stat()
            file_stats[str(path)] = dict(size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                                         ctime_ns=stat.st_ctime_ns)
        except FileNotFoundError:
            file_stats[str(path)] = None
    return dict(schema_version=1, selection_sha256=sha256(selection_signature),
                region=json.loads(pd.DataFrame([region._asdict()]).to_json(orient="records"))[0],
                ft_root=str(FT_ROOT.resolve()),
                extraction_sha256=hashlib.sha256(inspect.getsource(load_m6a_binary).encode()).hexdigest(),
                source_sha256={str(path): sha256(path) for path in sources},
                input_files=file_stats)


def load_or_build_m6a_binary(output, region, selection_signature):
    """Cache extraction before clustering so failed report preparation can resume."""
    destination = input_directory(output, region.region_id) / "m6a_input.npz"
    provenance = m6a_input_provenance(region, selection_signature)
    if destination.exists():
        try:
            with np.load(destination, allow_pickle=False) as saved:
                previous = json.loads(str(saved["provenance"].item()))
                if previous == provenance:
                    matrix = saved["matrix"]
                    records = pd.read_json(io.StringIO(str(saved["records"].item())), orient="table")
                    dtypes = json.loads(str(saved["record_dtypes"].item()))
                    records = records.astype(dtypes)
                    log = str(saved["extraction_log"].item())
                    if (matrix.dtype != np.uint8 or matrix.shape != (len(records), region.width)
                            or not np.isin(matrix, [0, 1]).all()
                            or records.row_id.isna().any() or records.row_id.duplicated().any()
                            or not records.region_id.eq(region.region_id).all()):
                        raise ValueError("matrix dimensions, values, or read identifiers are invalid")
                    # Recompute these exact binary-matrix summaries after the JSON round trip.
                    records["m6a_calls"] = matrix.sum(axis=1)
                    records["m6a_call_fraction"] = matrix.mean(axis=1)
                    print(f"{region.region_id}: reusing saved m6A matrix ({len(records)} reads)", flush=True)
                    return matrix, records, log
        except Exception as error:
            raise RuntimeError(f"Cannot read cached m6A input {destination}: {error}. "
                               "Remove this cache file and rerun to rebuild it from the original inputs.") from error
        print(f"{region.region_id}: m6A cache inputs changed; rebuilding", flush=True)
    else:
        print(f"{region.region_id}: building m6A input from original BEDs", flush=True)
    matrix, records, log = load_m6a_binary(region, selection_signature)
    if m6a_input_provenance(region, selection_signature) != provenance:
        raise RuntimeError(f"m6A inputs changed during extraction for {region.region_id}; rerun before caching")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(prefix=".m6a-input-", suffix=".npz",
                                         dir=destination.parent, delete=False) as handle:
            temporary = Path(handle.name)
            np.savez_compressed(handle, matrix=matrix,
                records=records.to_json(orient="table", index=False, double_precision=15),
                record_dtypes=json.dumps({column: str(dtype) for column, dtype in records.dtypes.items()}),
                provenance=json.dumps(provenance, sort_keys=True), extraction_log=log)
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    print(f"{region.region_id}: saved m6A matrix to {destination}", flush=True)
    return matrix, records, log


def save_report_data(output, region, binary, profiles, records, edges, args, info,
                     averages, stats, counts, allele_tables):
    """Save portable R inputs; plotting and repeat-length analysis live in the Rmd."""
    ids = records.row_id.to_numpy()
    def table_data(table):
        return json.loads(table.to_json(orient="records", double_precision=15))

    region_data = dict(schema_version=1, region=table_data(pd.DataFrame([region._asdict()])),
                   records=json.loads(records.to_json(orient="records", double_precision=15)),
                   row_ids=ids.tolist(),
                   call_offsets=[np.flatnonzero(row).tolist() for row in binary],
                   acf=json.loads(pd.DataFrame(profiles).to_json(orient="values", double_precision=15)),
                   params=dict(seed=args.seed, k_eff=info.get("actual_neighbors", 0),
                               resolution=args.resolution, window_size=args.window_size,
                               n_features=profiles.shape[1], n_pcs=args.n_pcs),
                   tables=dict(cluster_acf_summary=table_data(averages),
                               cluster_summary=table_data(stats),
                               sample_composition=table_data(counts),
                               **{name: table_data(table) for name, table in allele_tables.items()}),
                   edges=dict(source=ids[edges[:, 0].astype(int)].tolist(),
                              target=ids[edges[:, 1].astype(int)].tolist(),
                              weight=edges[:, 2].tolist()))
    rscript, env = r_environment()
    destination = result_directory(output, region.region_id) / "report_data.rds"
    # Select definition chunks explicitly: changing include/purl in the notebook
    # must never run setup, data loading, figures, or knitting during preparation.
    region_data["report_functions"] = report_function_code()
    prepare = ('args <- commandArgs(trailingOnly=TRUE); stopifnot(length(args) == 2L); '
               'region_data <- jsonlite::fromJSON(paste(readLines(file("stdin"), warn=FALSE), '
               'collapse="\\n"), simplifyMatrix=TRUE); '
               'definitions <- region_data$report_functions; region_data$report_functions <- NULL; '
               'eval(parse(text=definitions), envir=.GlobalEnv); '
               'prepare_report_data(region_data, project=args[[1L]], destination=args[[2L]])')
    # Rscript inserts --args itself; passing another --args shifts the values
    # returned by commandArgs() and can turn the project directory into the target.
    result = subprocess.run([rscript, "--vanilla", "-e", prepare,
        str(PROJECT), str(destination)], input=json.dumps(region_data, allow_nan=False),
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    (result_directory(output, region.region_id) / "report_preparation.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"Report data preparation failed for {region.region_id}:\n{result.stderr}")
    if not destination.is_file() or not destination.stat().st_size:
        raise RuntimeError(f"Missing report data for {region.region_id}")
    return destination


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def report_function_code():
    """Extract named function/constant chunks without evaluating the notebook."""
    chunks = re.findall(r"^```\{r\s+([^,}\s]+)[^\n]*\}\n(.*?)^```\s*$",
                        R_REPORT.read_text(), re.MULTILINE | re.DOTALL)
    definitions = "\n\n".join(code for name, code in chunks
                                if name.startswith("function-") or name.endswith("-function"))
    for name in ("prepare_report_data", "autocor_payload_matrix", "empty_footprints"):
        if not re.search(rf"\b{name}\s*<-\s*function\b", definitions):
            raise ValueError(f"Missing report function {name} in {R_REPORT}")
    return definitions


def selection_for_run(output, requested):
    """Reuse a selection, or recover its original metadata from saved results."""
    if requested is not None:
        if not requested.is_file():
            raise FileNotFoundError(f"Selection signature does not exist: {requested}")
        return requested.resolve()
    saved = input_directory(output) / "selection_signature.rds"
    for candidate in (saved, SIGNATURE, RECOVERED_SIGNATURE):
        if candidate.is_file():
            return candidate.resolve()
    print("Recovering original region/sample metadata from saved Leiden results", flush=True)
    run_m6a_builder("--recover-selection", PROJECT / "LCL_project/Leiden_manhattan", saved)
    if not saved.is_file():
        raise FileNotFoundError("Selection recovery produced no signature; supply --selection-signature")
    return saved.resolve()


def prepare_run(output, selection, sample_table_path, args, n_features, regions):
    """Resume missing report bundles without overwriting figures or mixing parameters."""
    manifest = input_directory(output) / "report_run.json"
    parameters = dict(window_size=args.window_size, n_features=n_features,
                      n_pcs=args.n_pcs, n_neighbors=args.n_neighbors,
                      resolution=args.resolution, seed=args.seed,
                      selection_sha256=sha256(selection), sample_table_sha256=sha256(sample_table_path),
                      region_ids=regions.region_id.tolist())
    if manifest.exists():
        previous = json.loads(manifest.read_text())
        if previous != parameters:
            changed = sorted(key for key in set(previous) | set(parameters)
                             if previous.get(key) != parameters.get(key))
            raise ValueError("Existing report run has different settings/inputs: " + ", ".join(changed)
                             + ". Choose a new --output-dir to keep the runs separate.")
    else:
        untracked = list(result_directory(output).glob("*/report_data.rds"))
        if untracked:
            raise ValueError("Existing report_data.rds files have no report_run.json provenance; "
                             "choose a new --output-dir rather than replacing them.")
        input_directory(output).mkdir(parents=True, exist_ok=True)
        result_directory(output).mkdir(parents=True, exist_ok=True)
        manifest.write_text(json.dumps(parameters, indent=2) + "\n")
    for rid in regions.region_id:
        input_directory(output, rid).mkdir(parents=True, exist_ok=True)
        result_directory(output, rid).mkdir(parents=True, exist_ok=True)


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
    parser.add_argument("--selection-signature", type=Path, default=None,
                        help="Selection RDS; default reuses saved inputs or recovers original metadata from Leiden results")
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
    try:
        args.selection_signature = selection_for_run(output, args.selection_signature)
        regions = load_regions(args.window_size, args.selection_signature)
        available_regions = regions.copy()
        if args.region:
            unknown = set(args.region) - set(regions.region_id)
            if unknown:
                parser.error(f"Unknown region IDs: {sorted(unknown)}")
        prepare_run(output, args.selection_signature, args.sample_table, args, n_features, available_regions)
    except (ValueError, FileNotFoundError, RuntimeError) as error:
        parser.error(str(error))
    if args.region:
        regions = regions[regions.region_id.isin(args.region)]
    sources = [FUNCTIONS, PHASED_INPUT, R_BUILDER, args.selection_signature, args.sample_table, R_REPORT,
               PROJECT / "code/clustering_methods/Leiden_Manhattan/leiden_manhattan_functions.r",
               PROJECT / "code/haplotype_phasing/LCL_phasing.r",
               *sorted(Path(__file__).parent.glob("*.py"))]
    source_hashes = {str(p): sha256(p) for p in sources}
    print(f"Resolution={args.resolution}; window={args.window_size}; lags=0–{n_features - 1}; "
          f"PCs={args.n_pcs}; neighbors={args.n_neighbors}; seed={args.seed}. "
          "Saving m6A matrix caches, summary PDFs, portable report_data.rds files, and diagnostic logs.", flush=True)
    for region in regions.itertuples(index=False):
        rid = region.region_id
        report_path = result_directory(output, rid) / "report_data.rds"
        if report_path.is_file() and report_path.stat().st_size:
            print(f"{rid}: report_data.rds already complete; skipping", flush=True)
            continue
        preserve_figures = any(result_directory(output, rid).iterdir())
        binary, records, extraction_log = load_or_build_m6a_binary(output, region, args.selection_signature)
        log_name = "report_extraction.log" if preserve_figures else "extraction.log"
        (input_directory(output, rid) / log_name).write_text(extraction_log)
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
        if preserve_figures:
            allele_tables = summaries.summarize_alleles(profiles, records, region)
        else:
            summaries.plot_region(result_directory(output, rid), region, binary, profiles, records, avg, counts)
            allele_tables = write_allele_results(output, region, profiles, records)
        report_data = save_report_data(output, region, binary, profiles, records, edges, args, info,
                                      avg, stats, counts, allele_tables)
        print(f"{rid}: wrote {len(stats)} clusters' summaries and {report_data.name}", flush=True)
    if any(sha256(p) != digest for p, digest in source_hashes.items()):
        raise RuntimeError("An input/code source changed during the run")
    print(f"Complete: {output}. Knit {R_REPORT.name} to draw the R figures and print tables.", flush=True)


if __name__ == "__main__":
    main()
