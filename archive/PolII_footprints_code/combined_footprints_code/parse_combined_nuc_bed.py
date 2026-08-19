import os
import pysam
from collections import defaultdict

FIRE_NUC_SORTED = (
    "/project/spott/cshan/fiber-seq/results/PolII/FIRE_nuc_by_chr_combined_sample/"
    "combined_nuc_beds/sorted_indexed"
)
PARSED_OUT_DIR = (
    "/project/spott/cshan/fiber-seq/results/PolII/FIRE_nuc_by_chr_combined_sample/parsed_bed_files"
)


def parse_nuc_blocks(row):
    """Return (chrom, start, end, read_name) for each real nucleosome block in a bed12 row.

    Skips the first and last sentinel blocks that mark fiber start/end.
    """
    f = row.split("\t")
    if len(f) < 12:
        return []
    chrom = f[0]
    fiber_start = int(f[1])
    read_name = f[3]
    sizes  = [int(x) for x in f[10].rstrip(",").split(",") if x]
    starts = [int(x) for x in f[11].rstrip(",").split(",") if x]
    return [
        (chrom, fiber_start + s, fiber_start + s + bl, read_name)
        for s, bl in zip(starts[1:-1], sizes[1:-1])
    ]


def collect_reads_parsed(locus_chrom, view_start, view_end):
    """Collect footprints per read with correctly parsed bed12 FIRE nucleosomes.

    Drop-in replacement for collect_reads() in plot_polii_footprints.py once the
    parsed files are available; use with the sorted_indexed parsed bed4 files.
    """
    from plot_polii_footprints import FP_FILES, FP_COLORS, classify_by_size

    reads = defaultdict(list)

    for label, fp_file in FP_FILES.items():
        tbx = pysam.TabixFile(fp_file)
        try:
            for row in tbx.fetch(locus_chrom, view_start, view_end):
                f = row.split("\t")
                fp_start, fp_end, read_name = int(f[1]), int(f[2]), f[3]
                fp_size = fp_end - fp_start
                cls = classify_by_size(fp_size)
                reads[read_name].append((fp_start, fp_end, cls, FP_COLORS[cls]))
        except (ValueError, KeyError):
            pass
        tbx.close()

    parsed_nuc_path = os.path.join(
        PARSED_OUT_DIR, "sorted_indexed",
        f"all_samples_nuc_{locus_chrom}_parsed_sorted.bed.gz"
    )
    if os.path.exists(parsed_nuc_path):
        tbx = pysam.TabixFile(parsed_nuc_path)
        try:
            for row in tbx.fetch(locus_chrom, view_start, view_end):
                f = row.split("\t")
                if len(f) < 4:
                    continue
                nuc_start, nuc_end, read_name = int(f[1]), int(f[2]), f[3]
                reads[read_name].append(
                    (nuc_start, nuc_end, "FIRE_nucleosome", FP_COLORS["FIRE_nucleosome"])
                )
        except (ValueError, KeyError):
            pass
        tbx.close()
    return reads


# ── parse all chromosomes and write flat bed4 files ───────────────────────────
os.makedirs(PARSED_OUT_DIR, exist_ok=True)

for fname in sorted(os.listdir(FIRE_NUC_SORTED)):
    if not fname.endswith("_sorted.bed.gz") or fname.endswith(".tbi"):
        continue
    in_path  = os.path.join(FIRE_NUC_SORTED, fname)
    base     = fname.replace("_sorted.bed.gz", "")
    out_path = os.path.join(PARSED_OUT_DIR, f"{base}_parsed.bed")

    print(f"Parsing {fname} -> {os.path.basename(out_path)}", end=" ... ", flush=True)
    count = 0
    with open(out_path, "w") as out_fh:
        tbx = pysam.TabixFile(in_path)
        for row in tbx.fetch():
            for chrom, start, end, read_name in parse_nuc_blocks(row):
                out_fh.write(f"{chrom}\t{start}\t{end}\t{read_name}\n")
                count += 1
        tbx.close()
    print(f"{count:,} nucleosomes written")