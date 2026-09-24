#!/usr/bin/env python3
"""Per-molecule footprint matrices for THP-1 Fiber-seq, in the Tullius et al. layout.

The two notebooks in this directory both need the same thing: for a genomic
anchor, a reads x positions matrix whose value at each base is the size of the
footprint covering it. That is the representation the Tullius Data S1/S3
notebooks operate on, so their plotting and quantification functions port over
directly.

Matrix encoding (identical to the paper's preprocessed `genes_df`):
    size (>=10)  base is covered by a footprint of that size
    0            read covers this base, no footprint called
    2            read does not cover this base (paper's sentinel)

Footprint sources, on the FiberHMM output for this project:
    10-80 bp   firehmm_tf/ft_by_size/size{10-30,40-60,60-80} - BED6, one row per
               footprint. 40-60 = promoter-proximal paused Pol II (PPP),
               60-80 = pre-initiation complex (PIC).
    >=90 bp    firehmm_footprint BED12 blocks - nucleosomes. 99% of these blocks
               are >=90 bp, so the cut cleanly separates the two callers and
               nothing is double-counted. Sizes 30-40 and 80-90 have no caller
               and are absent from the matrix by construction.

The same nucleosomes are also available exploded as BED6 in
firehmm_footprint/ft_by_size/size90plus (split_footprints_by_size.sh footprint),
but this module keeps reading the BED12: it needs the per-read alignment span
from those rows to know which bases a read covers, which the BED6 drops.
"""

import colorsys
import glob
import os

import matplotlib
import numpy as np
import pandas as pd
import pyBigWig
import pysam

EXTRACT_ROOT = "/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract"
ANCHOR_DIR = "/project/spott/cshan/fiber-seq/macrophage_project/annotations/enhancer_anchors"
PUBLIC_DIR = "/project/spott/cshan/fiber-seq/macrophage_project/annotations/enhancer_public_data"

ENHANCER_ANCHORS = f"{ANCHOR_DIR}/thp1_enhancer_anchors.bed"
PROMOTER_ANCHORS = f"{ANCHOR_DIR}/thp1_promoter_anchors.bed"

TF_SIZE_BINS = ["10-30", "40-60", "60-80"]
NUC_MIN_SIZE = 90       # firehmm_footprint blocks at or above this are nucleosomes
NOT_COVERED = 2         # paper's sentinel for "read does not reach this base"

PPP_RANGE = (40, 60)    # promoter-proximal paused Pol II
PIC_RANGE = (60, 80)    # pre-initiation complex
NUC_RANGE = (90, 10_000)
TF_RANGE = (10, 30)


# --------------------------------------------------------------------------- #
# Colour maps - ports of the paper's make_cp(), which interpolates in HSL via
# the `colour` package. That package is not installed here, so the same HSL
# interpolation is done with colorsys to keep the palette identical.
# --------------------------------------------------------------------------- #
def _hsl_range(hex_start, hex_end, n):
    """n colours from hex_start to hex_end inclusive, interpolated in HSL."""
    r1, g1, b1 = matplotlib.colors.to_rgb(hex_start)
    r2, g2, b2 = matplotlib.colors.to_rgb(hex_end)
    h1, l1, s1 = colorsys.rgb_to_hls(r1, g1, b1)
    h2, l2, s2 = colorsys.rgb_to_hls(r2, g2, b2)
    if n == 1:
        return [(r1, g1, b1, 1.0)]
    out = []
    for i in range(n):
        f = i / (n - 1)
        r, g, b = colorsys.hls_to_rgb(
            h1 + (h2 - h1) * f, l1 + (l2 - l1) * f, s1 + (s2 - s1) * f
        )
        out.append((r, g, b, 1.0))
    return out


def make_cp(colors, sizes):
    """Port of the paper's make_cp: transparent below sizes[0], then gradients."""
    cp = [(0.0, 0.0, 0.0, 0.0)] * sizes[0]
    for i in range(len(colors) - 1):
        cp += _hsl_range(colors[i], colors[i + 1], sizes[i + 1] - sizes[i])
    return matplotlib.colors.ListedColormap(cp)


# The paper's footprint-size palette. Used with vmin=0, vmax=100, so the colour
# index is the footprint size in bp and anything >=100 bp is nucleosome green.
POL_CP = make_cp(
    ["#ed306f", "#e42cf6", "#f3aafa", "#265de6", "#33a040"],
    [20, 45, 55, 70, 100],
)
# Grey overlay used to block out footprints with no transcription interpretation.
MASK_CP = make_cp(["#949494", "#949494"], [1, 1000])

# Line colours for the metaprofiles, matching the paper's figure panels.
LINE_COLORS = {"PPP": "orchid", "PIC": "cornflowerblue", "nucleosome": "mediumseagreen"}

# Fill colours for the public bulk tracks, kept away from the footprint palette.
TRACK_COLORS = {"ATAC": "#555555", "H3K27ac": "#c0504d"}


# --------------------------------------------------------------------------- #
# Anchors
# --------------------------------------------------------------------------- #
def available_chroms(timepoint="LPS0"):
    """Chromosomes FiberHMM produced output for.

    The GENCODE-derived promoter anchors include unplaced scaffolds, which have no
    per-chromosome BED, so anchor sweeps must be filtered to this set first.
    """
    prefix = f"{timepoint}_hmm_extracted_footprint_"
    return {
        os.path.basename(p)[len(prefix):-len(".bed.gz")]
        for p in glob.glob(f"{EXTRACT_ROOT}/firehmm_footprint/{timepoint}/{prefix}*.bed.gz")
    }


def load_enhancer_anchors(path=ENHANCER_ANCHORS):
    """FANTOM5 consensus TSSs inside an ENCODE cCRE ELS (build_enhancer_annotation.sh)."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "chrom", "start", "end", "enhancer_id", "f5_score", "strand",
            "encode_class", "ccre_id", "enh_start", "enh_end", "dnase",
        ],
    )
    df["anchor"] = df["start"]
    return df


def load_promoter_anchors(path=PROMOTER_ANCHORS):
    """Protein-coding TSSs, THP-1 CAGE where available (build_promoter_anchors.sh)."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "gene_id", "gene_name", "strand", "tss_source"],
    )
    df["anchor"] = df["start"]
    return df


# --------------------------------------------------------------------------- #
# Public bulk tracks - LIMA THP-1 LPS series, lifted to hg38 by
# run_liftover_public_data.sh (ATAC GSE201351, H3K27ac ChIP GSE201352).
#
# The public series was sampled at 0/30/60/90/120/240/360/1440 min and this
# project's Fiber-seq at LPS0/5/10/15, so no timepoint is actually shared. The
# two are paired by rank - nth Fiber-seq timepoint to nth public timepoint - and
# the tracks keep the *public* label (H3K27ac_30min, not H3K27ac_LPS5) so the
# mismatch stays visible on the figure.
# --------------------------------------------------------------------------- #
PUBLIC_TIMEPOINT = {"LPS0": "0000", "LPS5": "0030", "LPS10": "0060", "LPS15": "0090"}

BIGWIG_GLOBS = {
    "ATAC": f"{PUBLIC_DIR}/atac_hg38/*_LIMA_ATAC_THP1_WT_LPIF_S_{{tp}}_*.bw",
    "H3K27ac": f"{PUBLIC_DIR}/chip_hg38/*_LIMA_ChIP_h3k27ac_THP1_WT_LPIF_S_{{tp}}_*.bw",
}


def public_track_label(assay, timepoint):
    """('H3K27ac', 'LPS5') -> 'H3K27ac_30min'."""
    return f"{assay}_{int(PUBLIC_TIMEPOINT[timepoint])}min"


_BIGWIG_CACHE = {}


def _bigwig(path):
    bw = _BIGWIG_CACHE.get(path)
    if bw is None:
        bw = _BIGWIG_CACHE[path] = pyBigWig.open(path)
    return bw


def bigwig_profiles(assay, chrom, anchor, window, strand="."):
    """Per-base public signal around `anchor`, one profile per Fiber-seq timepoint.

    Returns {timepoint: array of length 2*window} over [anchor-window,
    anchor+window), replicate bigWigs averaged and missing values set to 0. As
    with build_read_matrix the profile is reversed for '-' strand anchors, so
    positive offsets are downstream. All four timepoints come back together so
    the caller can put them on a shared y-axis.
    """
    width = 2 * window
    out = {}
    for tp, public_tp in PUBLIC_TIMEPOINT.items():
        reps = sorted(glob.glob(BIGWIG_GLOBS[assay].format(tp=public_tp)))
        if not reps:
            raise FileNotFoundError(
                f"no {assay} bigWig for public timepoint {public_tp} "
                f"({BIGWIG_GLOBS[assay].format(tp=public_tp)})"
            )
        vals = []
        for path in reps:
            bw = _bigwig(path)
            v = np.zeros(width)
            if chrom in bw.chroms():
                start = max(0, anchor - window)
                end = min(bw.chroms()[chrom], anchor + window)
                if end > start:
                    v[start - (anchor - window):end - (anchor - window)] = np.nan_to_num(
                        np.array(bw.values(chrom, start, end), dtype=float)
                    )
            vals.append(v)
        prof = np.mean(vals, axis=0)
        out[tp] = prof[::-1] if strand == "-" else prof
    return out


# --------------------------------------------------------------------------- #
# Footprint fetching
# --------------------------------------------------------------------------- #
def _tf_path(size_bin, chrom, timepoint):
    return (
        f"{EXTRACT_ROOT}/firehmm_tf/ft_by_size/size{size_bin}/{timepoint}/"
        f"{timepoint}_tf_size{size_bin}_{chrom}.bed.gz"
    )


def _nuc_path(chrom, timepoint):
    return (
        f"{EXTRACT_ROOT}/firehmm_footprint/{timepoint}/"
        f"{timepoint}_hmm_extracted_footprint_{chrom}.bed.gz"
    )


# Anchor-by-anchor sweeps reopen the same handful of files tens of thousands of
# times, so keep the handles around (a few per chromosome).
_TABIX_CACHE = {}


def _tabix(path):
    tbx = _TABIX_CACHE.get(path)
    if tbx is None:
        tbx = _TABIX_CACHE[path] = pysam.TabixFile(path)
    return tbx


def fetch_footprints(chrom, start, end, timepoint="LPS0"):
    """Collect per-read footprints overlapping [start, end).

    Returns (spans, footprints):
        spans      read_name -> (read_start, read_end) alignment span
        footprints read_name -> list of (fp_start, fp_end, size)
    """
    spans = {}
    footprints = {}

    tbx = _tabix(_nuc_path(chrom, timepoint))
    try:
        for row in tbx.fetch(chrom, start, end):
            f = row.split("\t")
            read = f[3]
            r_start, r_end = int(f[1]), int(f[2])
            spans[read] = (r_start, r_end)
            block_sizes = f[10].rstrip(",").split(",")
            block_starts = f[11].rstrip(",").split(",")
            hits = footprints.setdefault(read, [])
            for bs, st in zip(block_sizes, block_starts):
                size = int(bs)
                if size < NUC_MIN_SIZE:
                    continue
                fp_start = r_start + int(st)
                if fp_start < end and fp_start + size > start:
                    hits.append((fp_start, fp_start + size, size))
    except (ValueError, KeyError):
        pass

    for size_bin in TF_SIZE_BINS:
        tbx = _tabix(_tf_path(size_bin, chrom, timepoint))
        try:
            for row in tbx.fetch(chrom, start, end):
                f = row.split("\t")
                read = f[3]
                if read not in spans:
                    continue
                footprints.setdefault(read, []).append(
                    (int(f[1]), int(f[2]), int(f[4]))
                )
        except (ValueError, KeyError):
            pass

    return spans, footprints


def build_read_matrix(chrom, anchor, window, strand=".", timepoint="LPS0",
                      require_anchor_covered=True):
    """Reads x positions footprint-size matrix centred on `anchor`.

    Columns are stringified offsets from the anchor over [-window, window), in
    transcription-unit orientation: for a '-' strand anchor the window is
    reversed so that positive offsets are always downstream. Enhancer anchors
    are unstranded ('.') and are treated as '+'.
    """
    view_start = max(0, anchor - window)
    view_end = anchor + window
    spans, footprints = fetch_footprints(chrom, view_start, view_end, timepoint)

    width = 2 * window
    rows, names = [], []
    for read, (r_start, r_end) in spans.items():
        if require_anchor_covered and not (r_start <= anchor < r_end):
            continue
        arr = np.full(width, NOT_COVERED, dtype=np.int32)
        # bases of this window the read actually covers
        lo = max(r_start, view_start) - view_start
        hi = min(r_end, view_end) - view_start
        if hi <= lo:
            continue
        arr[lo:hi] = 0
        # nucleosomes first, then the shorter TF-scale calls drawn over them
        for fp_start, fp_end, size in sorted(footprints.get(read, []), key=lambda x: -x[2]):
            a = max(fp_start, view_start) - view_start
            b = min(fp_end, view_end) - view_start
            if b > a:
                arr[a:b] = size
        rows.append(arr)
        names.append(read)

    if not rows:
        return pd.DataFrame(columns=[str(i) for i in range(-window, window)])

    mat = np.vstack(rows)
    if strand == "-":
        mat = mat[:, ::-1]
    return pd.DataFrame(mat, index=names, columns=[str(i) for i in range(-window, window)])


# --------------------------------------------------------------------------- #
# Ports of the Tullius Data S1/S3 analysis functions
# --------------------------------------------------------------------------- #
def slice_positions(mat, start, end):
    """Columns of `mat` for offsets [start, end)."""
    return mat.loc[:, [str(i) for i in range(start, end)]]


def meta_line_plot(mat, start, end, fp_size_min, fp_size_max):
    """Fraction of covering reads with a footprint in [fp_size_min, fp_size_max).

    Port of the paper's meta_line_plot. Bounds are inclusive-low/exclusive-high
    to match how this project's size bins are defined, where the paper used
    strict inequalities on both sides.
    """
    covered = slice_positions(mat, start, end).ne(NOT_COVERED).sum(axis=0).to_numpy()
    sub = slice_positions(mat, start, end)
    hit = ((sub >= fp_size_min) & (sub < fp_size_max)).sum(axis=0).to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(covered > 0, hit / np.maximum(covered, 1), np.nan)


def size_position_density(mat, start, end, min_size, max_size, bin_size):
    """Footprint size (rows) x position (columns) count matrix.

    Port of the paper's plot_meta_heatmap: for each size bin, count how many
    reads carry a footprint of that size at each position.
    """
    sub = slice_positions(mat, start, end)
    out, index = [], []
    for s in range(min_size // bin_size, max_size // bin_size):
        lo, hi = s * bin_size, (s + 1) * bin_size
        out.append(((sub >= lo) & (sub < hi)).sum(axis=0).to_numpy())
        index.append(lo)
    return pd.DataFrame(out, index=index, columns=sub.columns)


def dist_to_fp(mat, start, end, min_size, max_size):
    """Offset of the first base covered by a footprint in [min_size, max_size).

    Port of the paper's dist_to_fp, used to order reads in the browser plots.
    Reads with no such footprint in the range are dropped, as in the paper.
    """
    sub = slice_positions(mat, start, end)
    hit = (sub >= min_size) & (sub < max_size)
    hit = hit.loc[hit.sum(axis=1) > 0]
    if hit.empty:
        return pd.Series(dtype=float)
    return hit.idxmax(axis=1).astype(int)


def order_reads_by_nucleosome(mat, start=0, end=300, min_size=NUC_MIN_SIZE):
    """Reorder `mat` the way the paper's example-locus figures do.

    Reads are sorted by how far downstream of the anchor the first nucleosome
    begins; reads with no downstream nucleosome in the range are dropped.
    """
    d = dist_to_fp(mat, start, end, min_size, NUC_RANGE[1])
    return mat.loc[d.sort_values().index]


def smooth(arr, width):
    """Boxcar smooth and trim the padding, as the paper does around each profile."""
    box = np.ones(width) / width
    return np.convolve(np.nan_to_num(arr), box, mode="same")[width:len(arr) - width]


def mask_layer(mat, max_size=PPP_RANGE[0]):
    """Grey overlay for footprints too short to be a Pol II complex (<40 bp).

    The paper greyed out small footprints that did not overlap a PRO-seq or
    CAGE peak. There is no matched transcription data for THP-1, so the
    equivalent here is to grey the TF-scale (10-30 bp) calls and leave the
    PPP/PIC/nucleosome size classes coloured.
    """
    return ((mat >= TF_RANGE[0]) & (mat < max_size)).astype(int) * 500
