#!/usr/bin/env python3
"""Validate hg38 GRO-cap inputs and export divergent TSSs from PINTS calls.

Use the validate subcommand before PINTS, then summarize after peak calling.
Keep outputs flat: do not write comparison files, metadata, checksums,
LATEST_RUN.txt, symlinks, or nested run directories.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import numpy as np
import pyBigWig


def prepare(manifest):
    records = list(csv.DictReader(manifest.open(), delimiter="\t"))
    assert len(records) == 2 and {r["strand"] for r in records} == {"+", "-"}
    chrom_sizes = None
    reference = dict(line.split() for line in
                     Path("/project/spott/cshan/annotations/hg38.chrom.sizes").read_text().splitlines())
    for row in records:
        source = Path(row["source_path"])
        assert source.is_file() and not source.is_symlink(), source
        assert source.parent == Path("/project/spott/cshan/annotations/LCL_public_data/GRO-cap")
        assert row["assembly"] == "GRCh38"
        positions = reads = 0
        with pyBigWig.open(str(source)) as bw:
            assert bw.isBigWig()
            sizes = bw.chroms()
            assert all(sizes[c] == int(reference[c]) for c in sizes), "Chromosome lengths do not match hg38"
            assert set(sizes) == {"chr" + str(c) for c in list(range(1, 23)) + ["X", "Y"]}
            if chrom_sizes is not None:
                assert sizes == chrom_sizes, "Plus/minus chromosome dictionaries differ"
            chrom_sizes = sizes
            for chrom in sizes:
                entries = bw.intervals(chrom)
                if entries is None:
                    continue
                values = np.asarray(entries, dtype=np.float64)
                assert np.all(values[:, 1] - values[:, 0] == 1), "Input is not single-base signal"
                signal = values[:, 2]
                assert np.isfinite(signal).all() and np.all(signal == np.rint(signal))
                assert np.all(signal >= 0) if row["strand"] == "+" else np.all(signal <= 0)
                positions += int(np.count_nonzero(signal))
                reads += int(np.abs(signal).sum())
        print(f"Validated {row['accession']} ({row['strand']}): hg38, "
              f"{positions} signal positions, {reads} read starts", flush=True)


def summarize(native, output, plus_bw, minus_bw):
    output.mkdir(parents=True, exist_ok=True)
    columns = ["region_id", "chrom", "start_0based", "end_exclusive", "pints_confidence",
               "plus_tss_0based", "minus_tss_0based", "plus_minus_distance_bp",
               "major_tss_divergent", "high_confidence", "plus_tss_reads", "minus_tss_reads"]
    confidence_counts, chrom_counts = Counter(), Counter()
    total = high = nondivergent_major = 0
    with native.open() as src, pyBigWig.open(str(plus_bw)) as pl, pyBigWig.open(str(minus_bw)) as mn, \
            (output / "divergent_regions.tsv").open("w") as table, \
            (output / "divergent_regions.bed").open("w") as regions, \
            (output / "high_confidence_divergent_regions.bed").open("w") as strict, \
            (output / "major_tss.bed").open("w") as tss, \
            (output / "high_confidence_major_tss.bed").open("w") as strict_tss:
        writer = csv.DictWriter(table, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for line in src:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            assert len(fields) == 6, "Unexpected PINTS output schema"
            chrom, start, end, confidence, plus, minus = fields
            start, end, plus, minus = map(int, (start, end, plus, minus))
            assert 0 <= start <= plus < end <= pl.chroms(chrom)
            assert start <= minus < end
            # PINTS chooses each strand's strongest summit after merging overlapping pairs.
            # These two representative summits can occasionally be convergent; flag them.
            divergent = plus > minus
            keep = confidence == "Stringent(qval)" and divergent
            plus_reads = pl.values(chrom, plus, plus + 1)[0]
            minus_reads = abs(mn.values(chrom, minus, minus + 1)[0])
            assert plus_reads > 0 and minus_reads > 0
            total += 1
            high += int(keep)
            nondivergent_major += int(not divergent)
            confidence_counts[confidence] += 1
            chrom_counts[chrom] += 1
            region_id = f"GM12878_divergent_{total:06d}"
            writer.writerow(dict(zip(columns, [region_id, chrom, start, end, confidence,
                plus, minus, plus - minus, int(divergent), int(keep), plus_reads, minus_reads])))
            # Native PINTS's six columns are NOT standard BED6. Export valid BED6 separately.
            bed = f"{chrom}\t{start}\t{end}\t{region_id}\t0\t.\n"
            regions.write(bed)
            if keep:
                strict.write(bed)
            for position, strand in sorted(((plus, "+"), (minus, "-"))):
                bed = f"{chrom}\t{position}\t{position + 1}\t{region_id}:{strand}\t0\t{strand}\n"
                tss.write(bed)
                if keep:
                    strict_tss.write(bed)
    summary = dict(assembly="GRCh38", experiment="ENCSR719WNF", biological_replicates=1,
                   pints_divergent_regions=total, high_confidence_divergent_regions=high,
                   major_tss_records=2 * total, high_confidence_major_tss_records=2 * high,
                   regions_with_nondivergent_representative_tss=nondivergent_major,
                   confidence_counts=dict(confidence_counts), regions_by_chromosome=dict(chrom_counts),
                   high_confidence_rule="confidence exactly Stringent(qval) and plus major TSS > minus major TSS",
                   coordinates="0-based, half-open; TSS positions are 0-based",
                   note="Divergent transcription alone does not distinguish promoters from enhancers.")
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("validate", help="Check hg38 input tracks without writing files")
    validate.add_argument("--manifest", type=Path, default=Path(__file__).with_name("inputs.tsv"))
    export = commands.add_parser("summarize", help="Export TSS tables directly into the output folder")
    export.add_argument("--native", type=Path, required=True)
    export.add_argument("--output", type=Path, required=True)
    export.add_argument("--plus-bw", type=Path, required=True)
    export.add_argument("--minus-bw", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "validate":
        prepare(args.manifest)
    else:
        summarize(args.native, args.output, args.plus_bw, args.minus_bw)
