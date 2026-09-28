"""Shared I/O and content-checked stage manifests; no ACF/clustering algorithms."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
from datetime import datetime

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
PROJECT = HERE.parents[3]
DEFAULT_OUT = PROJECT / "macrophage_project/auto_correlation/tss"
BINS = ("Q1_low", "Q2", "Q3", "Q4_high")


def log(message):
    print(f"[{datetime.now().isoformat(timespec='seconds')}] {message}", flush=True)


def load_parent(filename):
    path = PARENT / filename
    name = "tss_parent_" + path.stem
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


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


def stage_valid(outdir, stage, signature):
    manifest = Path(outdir) / "validation" / (stage + ".manifest.json")
    try:
        previous = json.loads(manifest.read_text())
        if previous["signature"] != signature:
            return False
        for name, digest in previous["outputs"].items():
            path = Path(outdir) / name
            if not path.is_file() or sha256(path) != digest:
                return False
        log(f"{stage}: reuse verified completed stage")
        return True
    except (OSError, ValueError, KeyError):
        return False


def finish_stage(outdir, stage, signature, outputs, details=None):
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


def assert_alignment(binary, metadata, row_ids):
    if binary.shape != (len(metadata), 2000):
        raise ValueError("Binary matrix must have one 2000-base row per molecule")
    if not np.array_equal(metadata.row_index.to_numpy(), np.arange(len(metadata))):
        raise ValueError("Metadata row_index must match matrix order exactly")
    if metadata.read_id.duplicated().any() or not np.array_equal(
            metadata.read_id.astype(str).to_numpy(), row_ids.astype(str)):
        raise ValueError("Duplicate physical read IDs or changed row order")
    if not np.isin(binary, [0, 1]).all():
        raise ValueError("Binary matrix contains values other than zero and one")
    if not np.array_equal(binary.sum(axis=1), metadata.m6a_count.to_numpy()):
        raise ValueError("m6A counts and binary matrix differ")
