#!/usr/bin/env python3
"""Plot macrophage nucleosome footprints per read, faceted by topic-model cluster.

Combines two existing pieces:
  * per-read cluster labels from the m6A topic model
    (macrophage_project/topic_model/<outname>/read_topic_assignments_<outname>.tsv),
  * mono-nucleosome (130-160 bp) FiberHMM footprints, extracted the same way as
    process_nuc_pos.R (block sizes in [FP_SIZE_MIN, FP_SIZE_MAX), midpoint in window).

For each promoter gene, reads are drawn as rows of nucleosome rectangles in a
browser-style panel, one facet per topic cluster. Footprints are placed by their
absolute genomic coordinate; reads within each cluster are ordered exactly as the
m6A cluster heatmap orders them (metA_mat / topic-table row order) - no sorting.
All reads shown in the m6A heatmap are drawn (no span / max-reads filter); a read
with no nucleosome footprint in the window shows just its alignment baseline.

Anchor:
  * if a genomic interval is given (default for the 9 promoters), use its midpoint;
  * otherwise try the CAGE dominant TSS for the gene, then the gencode TSS.

Example:
  python3 plot_nuc_by_topic_cluster.py --gene IL1B --window-size 1000
"""

import argparse
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pysam

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
TOPIC_DIR = "/project/spott/cshan/fiber-seq/macrophage_project/topic_model"
FP_ROOT = "/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_footprint"
CAGE_PATH = "/project/spott/cshan/fiber-seq/macrophage_project/annotations/fantom5_thp1/THP1.consensus.CAGE_peaks.dominantPerGene.bed.gz"
GENCODE_TSS_PATH = "/project/spott/cshan/annotations/TSS.gencode.v49.bed"
OUT_DIR = "/project/spott/cshan/fiber-seq/macrophage_project/nucleosome_positioning/results/nuc_by_topic_cluster"

# Mono-nucleosome size filter, identical to process_nuc_pos.R.
FP_SIZE_MIN = 130
FP_SIZE_MAX = 160

# Selected promoter genes (same intervals as process_nuc_pos.R / the topic-model Rmd).
SELECTED_PROMOTERS = [
    ("IL1B", "chr2", 112836323, 112837385),
    ("TNF", "chr6", 31575265, 31575824),
    ("IL1A", "chr2", 112793233, 112794969),
    ("PTGS2", "chr1", 186679985, 186681163),
    ("IL8", "chr4", 73740006, 73741332),
    ("CXCL2", "chr4", 74098863, 74099651),
    ("IL6", "chr7", 22726376, 22727094),
    ("RELB", "chr19", 45001011, 45002417),
    ("rs2836882", "chr21", 39094468, 39094712),
]

# Symbol -> gencode ENSG (IL8 = CXCL8). rs2836882 is a SNP with no gencode gene.
SELECTED_ENSG = {
    "IL1B": "ENSG00000125538", "TNF": "ENSG00000232810", "IL1A": "ENSG00000115008",
    "PTGS2": "ENSG00000073756", "IL8": "ENSG00000169429", "CXCL2": "ENSG00000081041",
    "IL6": "ENSG00000136244", "RELB": "ENSG00000104856",
}

# Cluster colours: colors_25 from the topic-model Rmd (Kevin Luo's plots.R), the same
# palette used for RowSideColors in the m6A cluster heatmap. Indexed by cluster number - 1.
# R color names resolved to hex via col2rgb (dodgerblue2, green4, black).
CLUSTER_COLORS = [
    "#1C86EE", "#E31A1C", "#008B00", "#6A3D9A", "#FF7F00",
    "#000000", "#FFD700", "#7EC0EE", "#FB9A99", "#90EE90",
]


# ---------------------------------------------------------------------------
# TSS / anchor resolution
# ---------------------------------------------------------------------------

def load_gencode_tss(path, ensg_map):
    """gencode TSS point + strand for the selected ENSGs (version-stripped)."""
    want = {v: k for k, v in ensg_map.items()}
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            ensg = f[3].split(".")[0]
            gene = want.get(ensg)
            if gene and gene not in out:
                out[gene] = (f[0], int(f[1]), f[5])  # chrom, tss(0-based), strand
    return out


def cage_tss_for_gene(ensg, cage_path):
    """First CAGE dominant peak whose gene_id base matches ensg -> (chrom, pos, strand)."""
    tbx = pysam.TabixFile(cage_path)
    try:
        for chrom in tbx.contigs:
            for row in tbx.fetch(chrom):
                f = row.split("\t")
                if f[3].split(".")[0] == ensg:
                    tbx.close()
                    return f[0], int(f[1]), f[5]
    finally:
        tbx.close()
    return None


def resolve_anchor(gene, interval, ensg, gencode_tss, use_interval=True):
    """Resolve (center, strand, source).

    interval given + use_interval -> promoter-interval midpoint (strand from gencode).
    otherwise -> CAGE dominant TSS, else gencode TSS.
    """
    strand = gencode_tss.get(gene, (None, None, "+"))[2]
    if interval is not None and use_interval:
        _chrom, start, end = interval
        return (start + end) // 2, strand, "interval_midpoint"

    if ensg is not None:
        cage = cage_tss_for_gene(ensg, CAGE_PATH)
        if cage is not None:
            return cage[1], cage[2], "cage_tss"
    if gene in gencode_tss:
        _c, tss, strand = gencode_tss[gene]
        return tss, strand, "gencode_tss"
    raise ValueError(f"No interval, CAGE, or gencode TSS available for {gene}.")


# ---------------------------------------------------------------------------
# Footprint extraction (mirrors process_nuc_pos.R)
# ---------------------------------------------------------------------------

def _sample_to_dir(sample_name):
    """topic-table sample_name 'LPS_0' -> FiberHMM dir 'LPS0'."""
    return sample_name.replace("_", "")


def collect_nuc_footprints(read_sample, locus_chrom, view_start, view_end):
    """Nucleosome footprints per read within the window, from each read's sample BED.

    read_sample: dict RID -> sample_name (topic table). Only these reads are kept.
    Parses BED12 blocks, keeps 130 <= size < 160, keeps footprints whose midpoint
    falls inside [view_start, view_end]. Returns dict RID -> [(fp_start, fp_end)].
    """
    reads = defaultdict(list)
    wanted_by_dir = defaultdict(set)
    for rid, sample in read_sample.items():
        wanted_by_dir[_sample_to_dir(sample)].add(rid)

    for sample_dir, wanted in wanted_by_dir.items():
        bed = os.path.join(
            FP_ROOT, sample_dir,
            f"{sample_dir}_hmm_extracted_footprint_{locus_chrom}.bed.gz",
        )
        if not os.path.exists(bed):
            print(f"  [WARN] missing footprint BED: {bed}")
            continue
        tbx = pysam.TabixFile(bed)
        try:
            for row in tbx.fetch(locus_chrom, view_start, view_end):
                f = row.split("\t")
                rid = f[3]
                if rid not in wanted:
                    continue
                cstart = int(f[1])
                sizes = f[10].split(",")
                starts = f[11].split(",")
                for sz, st in zip(sizes, starts):
                    if not sz or not st:
                        continue
                    size = int(sz)
                    if size < FP_SIZE_MIN or size >= FP_SIZE_MAX:
                        continue
                    fp_start = cstart + int(st)
                    fp_end = fp_start + size
                    fp_mid = fp_start + size // 2
                    if view_start <= fp_mid <= view_end:
                        reads[rid].append((fp_start, fp_end))
        except (ValueError, KeyError):
            pass
        tbx.close()
    return reads


# ---------------------------------------------------------------------------
# Topic table
# ---------------------------------------------------------------------------

def load_topic_table(outname):
    """RID -> (sample_name, cluster, read_start, read_end, order) from the topic TSV.

    `order` is the row index in the TSV, which is the metA_mat / fit$L row order the
    m6A cluster heatmap uses within each cluster (heatmap orders by order(clusters),
    stable within a cluster). Reusing it makes the footprint rows match the heatmap.
    """
    path = os.path.join(TOPIC_DIR, outname, f"read_topic_assignments_{outname}.tsv")
    header = None
    rows = {}
    idx = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = {name: i for i, name in enumerate(f)}
                continue
            rid = f[header["RID"]]
            rows[rid] = (
                f[header["sample_name"]],
                f[header["cluster"]],
                int(f[header["start"]]),
                int(f[header["end"]]),
                idx,
            )
            idx += 1
    return rows


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_gene(gene, chrom, start, end, window_size,
              gencode_tss, use_interval, out_dir):
    outname = f"{gene}_{chrom}_{start}_{end}"
    topic = load_topic_table(outname)

    ensg = SELECTED_ENSG.get(gene)
    interval = (chrom, start, end)
    center, strand, tss_source = resolve_anchor(
        gene, interval, ensg, gencode_tss, use_interval=use_interval)
    view_start = max(0, center - window_size)
    view_end = center + window_size

    print(f"Gene   : {gene}")
    print(f"Locus  : {chrom}:{view_start:,}-{view_end:,}  strand={strand}")
    print(f"Anchor : {center:,}  ({tss_source})")
    print(f"Window : {view_end - view_start:,} bp")

    read_sample = {rid: v[0] for rid, v in topic.items()}
    reads = collect_nuc_footprints(read_sample, chrom, view_start, view_end)

    # Every read in the topic table is a heatmap read (a metA_mat row). Include ALL of
    # them (no span / max-reads filter), grouped by cluster and ordered by the
    # topic-table / metA_mat row order, exactly as the m6A heatmap orders reads.
    by_cluster = defaultdict(list)
    for rid, (_sample, cluster, r_start, r_end, order) in topic.items():
        by_cluster[cluster].append((order, rid, reads.get(rid, []), r_start, r_end))

    clusters = sorted(by_cluster)
    for cl in clusters:
        by_cluster[cl].sort(key=lambda t: t[0])

    n_by_cluster = {cl: len(by_cluster[cl]) for cl in clusters}
    total = sum(n_by_cluster.values())
    with_nuc = sum(1 for cl in clusters for t in by_cluster[cl] if t[2])
    print(f"Reads (all heatmap reads): {total:,}  ({with_nuc:,} with a nucleosome in window)  "
          + ", ".join(f"{cl}={n}" for cl, n in n_by_cluster.items()))
    if total == 0:
        print(f"  [SKIP] {gene}: no reads to plot")
        return None

    # one facet per cluster, heights proportional to read counts
    heights = [max(1, n_by_cluster[cl]) for cl in clusters]
    fig = plt.figure(figsize=(14, 0.6 + 0.05 * total + 0.4 * len(clusters)))
    gs = plt.GridSpec(len(clusters), 1, height_ratios=heights, hspace=0.12)
    axes = [fig.add_subplot(gs[i]) for i in range(len(clusters))]

    for ax, cl in zip(axes, clusters):
        color = CLUSTER_COLORS[(int(cl.replace("cluster", "")) - 1) % len(CLUSTER_COLORS)]
        read_list = by_cluster[cl]
        n = len(read_list)
        for idx, (_order, rid, footprints, r_start, r_end) in enumerate(read_list):
            y = n - 1 - idx
            # read baseline = alignment span (topic table), clipped to the view
            ax.plot([max(r_start, view_start), min(r_end, view_end)], [y, y],
                    color="black", lw=0.3, zorder=1)
            for fp_start, fp_end in footprints:
                ax.add_patch(plt.Rectangle(
                    (fp_start, y - 0.4), fp_end - fp_start, 0.8,
                    linewidth=0, facecolor=color, alpha=0.9, zorder=2))
        ax.axvline(center, color="black", lw=1, ls="--", alpha=0.5)
        ax.set_xlim(view_start, view_end)
        ax.set_ylim(-1, max(1, n))
        ax.set_yticks([])
        ax.set_ylabel(f"{cl}\n(n={n})", fontsize=8, rotation=0, ha="right", va="center")
        # absolute genomic coordinates (no offset), matching the m6A heatmap/met-profile axis
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x):,}"))

    axes[0].set_title(
        f"{gene}  |  {chrom}:{view_start:,}-{view_end:,}  ({strand})  |  "
        f"nucleosome {FP_SIZE_MIN}-{FP_SIZE_MAX} bp  |  reads in heatmap order",
        fontsize=11)
    for ax in axes[:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(axes[-1].get_xticklabels(), rotation=30, ha="right")
    axes[-1].set_xlabel(f"Genomic position  ({chrom})", fontsize=9)

    os.makedirs(out_dir, exist_ok=True)
    out_base = os.path.join(out_dir, f"{gene}_nuc_by_cluster")
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_base}.pdf\n")
    return out_base + ".pdf"


def build_argparser():
    p = argparse.ArgumentParser(
        description="Nucleosome footprints per read, faceted by topic cluster.")
    p.add_argument("--gene", default=None,
                   help="Gene to plot (default: all 9 promoter genes).")
    p.add_argument("--window-size", type=int, default=1000,
                   help="Bases on each side of the anchor (default: 1000).")
    p.add_argument("--no-interval", action="store_true",
                   help="Ignore the promoter interval; anchor at CAGE TSS then gencode TSS.")
    p.add_argument("--out-dir", default=OUT_DIR, help="Directory to write PDF output.")
    return p


def main(argv=None):
    args = build_argparser().parse_args(argv)
    gencode_tss = load_gencode_tss(GENCODE_TSS_PATH, SELECTED_ENSG)

    genes = SELECTED_PROMOTERS
    if args.gene is not None:
        genes = [g for g in SELECTED_PROMOTERS if g[0] == args.gene]
        if not genes:
            raise ValueError(f"Gene '{args.gene}' is not in SELECTED_PROMOTERS.")

    for gene, chrom, start, end in genes:
        plot_gene(
            gene, chrom, start, end,
            window_size=args.window_size,
            gencode_tss=gencode_tss,
            use_interval=not args.no_interval,
            out_dir=args.out_dir,
        )


if __name__ == "__main__":
    main()
