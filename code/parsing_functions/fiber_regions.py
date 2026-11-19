"""Shared per-read region functions for Fiber-seq analyses (FIRE accessibility, footprints).

Python counterpart of parsing_footprints_functions.r. Analysis scripts give the input
paths and the rules; the functions here do the tabix / bedtools work.

Read rules for "a read covers a region" (rule=):
    "span"     the read's outer aligned span contains the whole region
               (bedtools intersect -f 1.0 against <s>.read_spans.bed.gz, BED6)
    "aligned"  every base of the region lies in the read's aligned M/=/X blocks, so a
               deletion or N gap inside the region breaks coverage
               (against <s>.aligned_blocks.bed.gz, BED12 from bamtobed -bed12 -splitD)

FIRE accessibility: ONE FIRE element on the same read covers >= fire_fraction of the
region length (bedtools -f on the region). Under "aligned", the element/region overlap
is clipped to the read's aligned blocks before the fraction is applied. Separate FIRE
elements are never summed. Accessible reads are restricted to covering reads; hits on
non-covering reads are counted as orphans.

Import with:
    sys.path.insert(0, str(<path to code/parsing_functions>))
    import fiber_regions as fr
"""

import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

BEDTOOLS = "/project/spott/cshan/envs/bedtools/bin/bedtools"
TABIX = "/project/spott/cshan/envs/dimelo/bin/tabix"
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
READ_RULES = ("span", "aligned")
MACROPHAGE_FIRE_ROOT = Path("/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE")


#################################################
# command helpers
#################################################
def command_lines(command):
    """Stream a command's stdout as split TSV fields; raise on a non-zero exit."""
    with tempfile.TemporaryFile(mode="w+t") as err:
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=err, text=True)
        try:
            for line in proc.stdout:
                yield line.rstrip("\n").split("\t")
        finally:
            proc.stdout.close()
            code = proc.wait()
        if code:
            err.seek(0)
            raise RuntimeError(f"{command[0]} failed: {err.read()[:4000]}")


def tabix_to_file(source, chrom, destination):
    """tabix <source> <chrom> > <destination>"""
    with open(destination, "w") as out:
        subprocess.run([TABIX, str(source), chrom], stdout=out, check=True)
    return destination


def regions_for_chrom(regions_bed, chrom, destination):
    """Rows of a BED file on one chromosome -> destination; returns (path, n_rows, n_cols)."""
    n, ncol = 0, 0
    with open(regions_bed) as src, open(destination, "w") as out:
        for line in src:
            if line.split("\t", 1)[0] == chrom:
                out.write(line)
                n += 1
                ncol = ncol or len(line.rstrip("\n").split("\t"))
    return destination, n, ncol


#################################################
# BED12 blocks
#################################################
def bed12_blocks(fields, validate=True):
    """Genomic (start, end) blocks of one BED12 record (fields 0-11); zero-width blocks dropped."""
    start, end, count = int(fields[1]), int(fields[2]), int(fields[9])
    sizes = [int(x) for x in fields[10].rstrip(",").split(",")]
    offsets = [int(x) for x in fields[11].rstrip(",").split(",")]
    blocks = [(start + o, start + o + s) for o, s in zip(offsets, sizes)]
    if validate and (len(sizes) != count or len(offsets) != count or not blocks or any(
            s < start or e > end or e < s or (i and s < blocks[i - 1][1])
            for i, (s, e) in enumerate(blocks))):
        raise ValueError(f"Invalid BED12 blocks for {fields[3]}")
    return [(s, e) for s, e in blocks if e > s]


def expand_blocks(fields, sizes=None, min_score=None, drop_sentinels=False):
    """Footprint blocks of one BED12+ record as (start, end, score) tuples.

    sizes          (lo, hi) inclusive block-size filter, or None
    min_score      keep blocks with column-13 per-block score >= min_score (FiberHMM tf), or None
    drop_sentinels drop the first and last block (ft extract BED12 sentinels)
    """
    base, count = int(fields[1]), int(fields[9])
    bsizes = fields[10].rstrip(",").split(",")
    offsets = fields[11].rstrip(",").split(",")
    scores = fields[12].rstrip(",").split(",") if len(fields) > 12 else [None] * count
    if not len(bsizes) == len(offsets) == len(scores) == count:
        raise ValueError(f"blockCount mismatch for {fields[3]}")
    out = []
    for i, (size, offset, score) in enumerate(zip(bsizes, offsets, scores)):
        if drop_sentinels and (i == 0 or i == count - 1):
            continue
        size = int(size)
        if sizes is not None and not sizes[0] <= size <= sizes[1]:
            continue
        if min_score is not None and float(score) < min_score:
            continue
        start = base + int(offset)
        out.append((start, start + size, score))
    return out


def aligned_overlap(blocks, start, end):
    """Bases of [start, end) inside a read's aligned blocks."""
    return sum(max(0, min(e, end) - max(s, start)) for s, e in blocks)


#################################################
# reads covering / accessible at each region
#################################################
def region_reads(chrom, regions_bed, reads_path, elements_path, tmpdir,
                 rule="span", fire_fraction=0.5):
    """({region_id: covering reads}, {region_id: FIRE-accessible reads}, n_orphan) on one chromosome.

    regions_bed    BED with region_id in column 4 (any chromosomes; only `chrom` is used)
    reads_path     tabix'd read spans (rule="span") or aligned blocks BED12 (rule="aligned")
    elements_path  tabix'd FIRE elements, read name in column 4; None skips accessibility
    """
    if rule not in READ_RULES:
        raise ValueError(f"rule must be one of {READ_RULES}")
    tmpdir = Path(tmpdir)
    reg, n, ncol = regions_for_chrom(regions_bed, chrom, tmpdir / f"regions.{chrom}.bed")
    if n == 0:
        reg.unlink()
        return {}, {}, 0
    rid_col = ncol + 3                      # read name follows the region columns
    reads = tabix_to_file(reads_path, chrom, tmpdir / f"reads.{chrom}.bed")

    cov = defaultdict(set)
    blocks_by_read = {}
    if rule == "span":
        for f in command_lines([BEDTOOLS, "intersect", "-a", str(reg), "-b", str(reads),
                                "-f", "1.0", "-wa", "-wb"]):
            cov[f[3]].add(f[rid_col])
    else:
        for f in command_lines([BEDTOOLS, "intersect", "-a", str(reg), "-b", str(reads), "-wa", "-wb"]):
            rid = f[rid_col]
            if rid not in blocks_by_read:
                blocks_by_read[rid] = bed12_blocks(f[ncol:ncol + 12])
            start, end = int(f[1]), int(f[2])
            if aligned_overlap(blocks_by_read[rid], start, end) == end - start:
                cov[f[3]].add(rid)

    acc = defaultdict(set)
    n_orphan = 0
    if elements_path is not None:
        elements = tabix_to_file(elements_path, chrom, tmpdir / f"elements.{chrom}.bed")
        for f in command_lines([BEDTOOLS, "intersect", "-a", str(reg), "-b", str(elements),
                                "-f", str(fire_fraction), "-wa", "-wb"]):
            rid = f[rid_col]
            if rule == "aligned":
                s, e = max(int(f[1]), int(f[ncol + 1])), min(int(f[2]), int(f[ncol + 2]))
                overlap = aligned_overlap(blocks_by_read.get(rid, ()), s, e)
                if not (overlap > 0 and overlap >= fire_fraction * (int(f[2]) - int(f[1]))):
                    continue
            acc[f[3]].add(rid)
        for r in acc:
            before = len(acc[r])
            acc[r] &= cov.get(r, set())
            n_orphan += before - len(acc[r])
        elements.unlink()
    reg.unlink()
    reads.unlink()
    return cov, acc, n_orphan


def footprints_in_regions(chrom, regions_bed, footprint_path, tmpdir, sizes=None,
                          min_score=None, drop_sentinels=False, fp_fraction=1.0):
    """{region_id: {read: n_blocks}} for footprint blocks overlapping each region.

    fp_fraction    minimum fraction of the FOOTPRINT inside the region (bedtools -F);
                   1.0 = footprint entirely inside, None = any overlap
    """
    tmpdir = Path(tmpdir)
    reg, n, _ = regions_for_chrom(regions_bed, chrom, tmpdir / f"fp_regions.{chrom}.bed")
    hits = defaultdict(lambda: defaultdict(int))
    if n == 0:
        reg.unlink()
        return hits
    fp = tmpdir / f"fp.{chrom}.bed"
    with open(fp, "w") as out:
        for f in command_lines([TABIX, str(footprint_path), chrom]):
            for s, e, _ in expand_blocks(f, sizes, min_score, drop_sentinels):
                out.write(f"{chrom}\t{s}\t{e}\t{f[3]}\n")
    cmd = [BEDTOOLS, "intersect", "-a", str(reg), "-b", str(fp), "-wa", "-wb"]
    if fp_fraction is not None:
        cmd[-2:-2] = ["-F", str(fp_fraction)]
    for f in command_lines(cmd):
        hits[f[3]][f[-1]] += 1
    reg.unlink()
    fp.unlink()
    return hits


#################################################
# per-sample input paths
#################################################
def add_lcl_fire_paths(samples, root):
    """Add donor and per-sample FIRE / read paths to an LCL metatable (sample_name, cell_line, fire_dir)."""
    root = Path(root)
    samples["donor"] = samples.cell_line
    samples["cram_path"] = [str(Path(r.fire_dir) / f"{r.sample_name}-fire-v0.1-filtered.cram") for r in samples.itertuples()]
    samples["peaks_path"] = [str(Path(r.fire_dir) / f"{r.sample_name}-fire-v0.1-peaks.bed.gz") for r in samples.itertuples()]
    samples["elements_path"] = [str(Path(r.fire_dir) / "additional-outputs-v0.1/fire-peaks" / f"{r.sample_name}-v0.1-fire-elements.bed.gz") for r in samples.itertuples()]
    samples["spans_path"] = [str(root / s / f"{s}.read_spans.bed.gz") for s in samples.sample_name]
    samples["blocks_path"] = [str(root / s / f"{s}.aligned_blocks.bed.gz") for s in samples.sample_name]
    return samples


def macrophage_elements_path(sample, fire_root=MACROPHAGE_FIRE_ROOT):
    return Path(fire_root) / sample / "additional-outputs-v0.1" / "fire-peaks" / f"{sample}-v0.1-fire-elements.bed.gz"


def reads_path(root, sample, rule):
    """<root>/<s>/<s>.read_spans.bed.gz (span) or <s>.aligned_blocks.bed.gz (aligned)."""
    name = "read_spans" if rule == "span" else "aligned_blocks"
    return Path(root) / sample / f"{sample}.{name}.bed.gz"
