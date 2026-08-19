#!/usr/bin/env python3
"""Convert precomputed footprint summary TSVs to sorted, tabix-indexed BED files.

The plotting summaries from ``precompute_fp_summaries.py`` are gzipped TSV files
partitioned by footprint class and chromosome.  They are useful, but not
coordinate-indexed, so querying one small locus requires scanning a whole
chromosome file.

This script converts each summary file into a BED-like file:

    chrom  fp_start  fp_end  read_name  fp_class  fp_mid

The output is sorted, bgzip-compressed, and indexed with tabix.  A plotting
script can then use ``pysam.TabixFile(...).fetch(chrom, start, end)`` to read
only the requested window.
"""

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys


DEFAULT_INPUT_DIR = "/project/spott/cshan/fiber-seq/results/PolII/footprint_summaries"
DEFAULT_OUTPUT_DIR = "/project/spott/cshan/fiber-seq/results/PolII/footprint_summary_beds"
ALL_CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

CLASS_SUBDIRS = {
    "PPP": "PPP",
    "PIC": "PIC",
    "nucleosome (140-160)": "nucleosome_140-160",
    "FIRE_nucleosome": "FIRE_nucleosome",
    "unknown": "unknown",
}

CLASS_ALIASES = {
    "ppp": "PPP",
    "pic": "PIC",
    "nucleosome": "nucleosome (140-160)",
    "nucleosome_140-160": "nucleosome (140-160)",
    "nucleosome (140-160)": "nucleosome (140-160)",
    "fire_nucleosome": "FIRE_nucleosome",
    "fire-nucleosome": "FIRE_nucleosome",
    "fire": "FIRE_nucleosome",
    "unknown": "unknown",
    "all": "all",
}

REQUIRED_COLUMNS = {"chrom", "read_name", "fp_start", "fp_end", "fp_mid", "fp_class"}


def normalize_classes(values):
    classes = []
    for value in values:
        key = str(value).strip().lower()
        cls = CLASS_ALIASES.get(key)
        if cls is None:
            choices = ", ".join(["all"] + list(CLASS_SUBDIRS))
            raise ValueError(f"Unknown class '{value}'. Choose one of: {choices}")
        if cls == "all":
            return list(CLASS_SUBDIRS)
        classes.append(cls)
    return list(dict.fromkeys(classes))


def require_tool(name):
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(
            f"Required command '{name}' was not found on PATH. "
            "On this cluster, try: module load htslib"
        )
    return path


def convert_one(input_path, output_path, tmp_dir=None, force=False):
    """Convert one gzipped summary TSV to sorted BED.gz + tabix index."""
    index_path = output_path + ".tbi"
    if os.path.exists(output_path) and os.path.exists(index_path) and not force:
        return {"status": "skipped", "rows": None, "bad_rows": None}

    if not os.path.exists(input_path):
        return {"status": "missing", "rows": 0, "bad_rows": 0}

    require_tool("bgzip")
    require_tool("tabix")
    require_tool("sort")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    tmp_root = tmp_dir or os.path.dirname(output_path)
    os.makedirs(tmp_root, exist_ok=True)

    base = os.path.basename(output_path)
    # Keep the final compressed temp beside the destination so os.replace()
    # stays on one filesystem even when --tmp-dir points at /tmp.
    tmp_output_path = os.path.join(os.path.dirname(output_path), f".{base}.tmp")

    rows = 0
    bad_rows = 0
    try:
        with gzip.open(input_path, "rt", newline="") as in_fh:
            reader = csv.DictReader(in_fh, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"{input_path} has no header")
            missing = REQUIRED_COLUMNS.difference(reader.fieldnames)
            if missing:
                raise ValueError(
                    f"{input_path} is missing required columns: {sorted(missing)}"
                )

            with open(tmp_output_path, "wb") as compressed_fh:
                sort_cmd = ["sort", "-T", tmp_root, "-k1,1", "-k2,2n", "-k3,3n"]
                sort_proc = subprocess.Popen(
                    sort_cmd,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                )
                bgzip_proc = subprocess.Popen(
                    ["bgzip", "-c"],
                    stdin=sort_proc.stdout,
                    stdout=compressed_fh,
                )
                sort_proc.stdout.close()

                try:
                    for row in reader:
                        try:
                            chrom = row["chrom"]
                            start = int(row["fp_start"])
                            end = int(row["fp_end"])
                            read_name = row["read_name"]
                            fp_class = row["fp_class"]
                            fp_mid = row["fp_mid"]
                        except (TypeError, ValueError, KeyError):
                            bad_rows += 1
                            continue

                        if not chrom or not read_name or end < start:
                            bad_rows += 1
                            continue

                        sort_proc.stdin.write(
                            (
                                f"{chrom}\t{start}\t{end}\t{read_name}"
                                f"\t{fp_class}\t{fp_mid}\n"
                            ).encode()
                        )
                        rows += 1
                finally:
                    if sort_proc.stdin:
                        sort_proc.stdin.close()

                bgzip_return = bgzip_proc.wait()
                sort_return = sort_proc.wait()

        if sort_return != 0:
            raise RuntimeError(f"sort failed for {input_path}")
        if bgzip_return != 0:
            raise RuntimeError(f"bgzip failed for {input_path}")

        os.replace(tmp_output_path, output_path)
        subprocess.run(["tabix", "-f", "-p", "bed", output_path], check=True)
    finally:
        if os.path.exists(tmp_output_path):
            os.remove(tmp_output_path)

    return {"status": "created", "rows": rows, "bad_rows": bad_rows}


def build_argparser():
    parser = argparse.ArgumentParser(
        description="Convert footprint summary TSV.gz files to sorted BED.gz + tabix indexes."
    )
    parser.add_argument(
        "--input-dir",
        default=DEFAULT_INPUT_DIR,
        help=f"Root directory from precompute_fp_summaries.py (default: {DEFAULT_INPUT_DIR}).",
    )
    parser.add_argument(
        "--out-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output root directory for BED.gz files (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--classes",
        nargs="+",
        default=["all"],
        help="Classes to process: PPP, PIC, nucleosome_140-160, FIRE_nucleosome, unknown, all.",
    )
    parser.add_argument(
        "--chroms",
        nargs="+",
        default=ALL_CHROMS,
        help="Chromosomes to process (default: chr1-chr22 chrX chrY).",
    )
    parser.add_argument(
        "--tmp-dir",
        default=None,
        help="Directory for temporary unsorted BED files (default: output class directory).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing BED.gz and .tbi outputs.",
    )
    return parser


def main(argv=None):
    args = build_argparser().parse_args(argv)
    classes = normalize_classes(args.classes)

    print(f"Input directory : {args.input_dir}")
    print(f"Output directory: {args.out_dir}")
    print(f"Classes         : {', '.join(classes)}")
    print(f"Chromosomes     : {' '.join(args.chroms)}")
    if args.tmp_dir:
        print(f"Temporary dir   : {args.tmp_dir}")
    print("")

    total_created = 0
    total_skipped = 0
    total_rows = 0
    total_bad = 0

    for cls in classes:
        subdir = CLASS_SUBDIRS[cls]
        for chrom in args.chroms:
            input_path = os.path.join(args.input_dir, subdir, f"{chrom}_footprints.tsv.gz")
            output_path = os.path.join(args.out_dir, subdir, f"{chrom}_footprints.bed.gz")

            result = convert_one(
                input_path,
                output_path,
                tmp_dir=args.tmp_dir,
                force=args.force,
            )
            status = result["status"]
            rows = result["rows"]
            bad_rows = result["bad_rows"]

            if status == "created":
                total_created += 1
                total_rows += rows
                total_bad += bad_rows
                print(
                    f"[created] {subdir}/{chrom}: {rows:,} rows"
                    + (f" ({bad_rows:,} skipped malformed)" if bad_rows else "")
                )
            elif status == "skipped":
                total_skipped += 1
                print(f"[skipped] {subdir}/{chrom}: output exists")
            else:
                print(f"[missing] {subdir}/{chrom}: {input_path}", file=sys.stderr)

    print("")
    print(f"Created files: {total_created}")
    print(f"Skipped files: {total_skipped}")
    print(f"Rows written : {total_rows:,}")
    print(f"Bad rows     : {total_bad:,}")


if __name__ == "__main__":
    main()
