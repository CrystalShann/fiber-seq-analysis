#!/usr/bin/env python3
"""Plot Pol II footprint tracks around a single gene.

This is the notebook code from PolII_FP.ipynb turned into a reusable module.
You can:
  1. run it from the command line with --gene ACTB
  2. import it from another script or notebook and call plot_gene_footprints()

Example:
  python3 plot_polii_footprints_precomputed_sort.py --gene ACTB --window-size 2000 \
      --fp-class FIRE_nucleosome --plot-mode midpoint
"""

import argparse
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import pyBigWig


PEAKS_DIR = "/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged/combined/joint_trained_peaks/"
TRACKS_DIR = "/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged/combined/joint_trained_tracks/"
FIRE_NUC_DIR = (
    "/project/spott/cshan/fiber-seq/results/PolII/FIRE_nuc_by_chr_combined_sample/"
    "parsed_bed_files"
)

# Non-overlapping footprint size range files
FP_FILES = {
    "10-30": PEAKS_DIR + "combined_all_chrs_10-30bp_fps.bed.gz",
    "30-45": PEAKS_DIR + "combined_all_chrs_30-45bp_fps.bed.gz",
    "45-60": PEAKS_DIR + "combined_all_chrs_45-60bp_fps.bed.gz",
    "60-80": PEAKS_DIR + "combined_all_chrs_60-80bp_fps.bed.gz",
    "140-160": PEAKS_DIR + "combined_all_chrs_140-160bp_fps.bed.gz",
}

GROCAP_CTSS = "/project/spott/cshan/annotations/GRO-cap/grocap.hg38.GM12878.combined.CTSS.bed.gz"
CAGE_PEAKS = "/project/spott/cshan/annotations/fantom5/fantom5.hg38.LCL.consensus.CAGE_peaks.withGene.bed.gz"
GENCODE_GTF = "/project/spott/cshan/annotations/gencode.v49.primary_assembly.annotation.gtf"

GENES_DF_CACHE = "/project/spott/cshan/fiber-seq/results/PolII/annotations/gencode_protein_coding_genes_with_TSS_TES.csv"
ANNOTATION_DIR = "/project/spott/cshan/fiber-seq/results/PolII/annotations"
PAUSING_INDEX_PRINCIPAL_CSV = f"{ANNOTATION_DIR}/pausing_index_principal.csv"
CLOSEST_CAGE_TSS_CSV = f"{ANNOTATION_DIR}/closest_CAGE_TSS_to_genes_in_pausing_index_principal.csv"
DEFAULT_PRINCIPAL_WITH_CAGE_CSV = f"{ANNOTATION_DIR}/pausing_index_principal_with_CAGE_TSS_all_genes.tsv"
DEFAULT_OUT_DIR = "/project/spott/cshan/fiber-seq/results/PolII/plots"
DEFAULT_FP_SUMMARY_DIR = "/project/spott/cshan/fiber-seq/results/PolII/footprint_summaries"
DEFAULT_FP_BED_DIR = "/project/spott/cshan/fiber-seq/results/PolII/footprint_summary_beds"

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

FP_COLORS = {
    "PPP": "#FF69B4",
    "PIC": "#4169E1",
    "nucleosome (140-160)": "#2E8B57",
    "FIRE_nucleosome": "#FFA500",
    "FIRE 80-100": "#08311CD2",
    "FIRE 100-120": "#00BCD4",
    "FIRE 120-140": "#856CDF",
    "unknown": "#AAAAAA",
}

TRACK_SIZES = ["10-30", "30-60", "60-80", "140-160"]
TRACK_COLORS = {
    "10-30": "#AAAAAA",
    "30-60": "#FF69B4",
    "60-80": "#4169E1",
    "140-160": "#2E8B57",
}
TRACK_LABELS = {
    "10-30": "10-30 bp",
    "30-60": "40-60 bp (PPP)",
    "60-80": "60-80 bp (PIC)",
    "140-160": "140-160 bp (nucleosome)",
}
N_BINS = 500
NUCLEOSOME_SORT_CLASSES = {"FIRE_nucleosome", "nucleosome (140-160)"}
DEFAULT_FIRE_NUC_LENGTH_RANGE = (140, 160)

FIRE_NUC_BW_DIR = "/project/spott/cshan/fiber-seq/results/PolII/FIRE_nuc_by_chr_combined_sample/bigwig_files"

FIRE_FP_PEAKS_DIR = "/project/spott/cshan/fiber-seq/results/PolII/FIRE_combined_footprints/joint_trained_peaks/"
FIRE_FP_TRACKS_DIR = "/project/spott/cshan/fiber-seq/results/PolII/FIRE_combined_footprints/joint_trained_tracks/"

FIRE_FP_FILES = {
    "80-100": FIRE_FP_PEAKS_DIR + "combined_all_chrs_80-100bp_fps.bed.gz",
    "100-120": FIRE_FP_PEAKS_DIR + "combined_all_chrs_100-120bp_fps.bed.gz",
    "120-140": FIRE_FP_PEAKS_DIR + "combined_all_chrs_120-140bp_fps.bed.gz",
}
FIRE_FP_SIZES = ["80-100", "100-120", "120-140"]
FIRE_FP_TRACK_COLORS = {
    "80-100": "#E67E00",
    "100-120": "#00BCD4",
    "120-140": "#FFBD00",
}
FIRE_FP_TRACK_LABELS = {
    "80-100": "80-100 bp (FIRE)",
    "100-120": "100-120 bp (FIRE)",
    "120-140": "120-140 bp (FIRE)",
}


def _read_table_auto(path):
    """Read comma/tab-separated tables with delimiter auto-detection."""
    df = pd.read_csv(path, sep=None, engine="python")
    if len(df.columns) == 1:
        only_col = str(df.columns[0])
        if "\t" in only_col:
            return pd.read_csv(path, sep="\t")
        if "," in only_col:
            return pd.read_csv(path, sep=",")
    return df


def build_principal_with_cage_fallback_df(
    pausing_csv=PAUSING_INDEX_PRINCIPAL_CSV,
    closest_cage_csv=CLOSEST_CAGE_TSS_CSV,
    out_csv=DEFAULT_PRINCIPAL_WITH_CAGE_CSV,
):
    """
    Build/save a principal table with all genes preserved.

    Rules:
    - Keep all rows/columns from pausing_index_principal.
    - Use closest CAGE TSS where available.
    - Fall back to GENCODE TSS for genes without CAGE.
    """
    principal_df = _read_table_auto(pausing_csv)
    cage_df = _read_table_auto(closest_cage_csv)

    cage_df = cage_df.rename(columns={"nearest_gene_name": "gene_name"})
    if "start" in cage_df.columns and "cage_tss" not in cage_df.columns:
        cage_df["cage_tss"] = cage_df["start"]

    if "distance_to_gene_tss" in cage_df.columns:
        cage_df["distance_to_gene_tss"] = pd.to_numeric(
            cage_df["distance_to_gene_tss"], errors="coerce"
        )
        cage_df = cage_df.sort_values("distance_to_gene_tss", kind="stable")
    cage_df = cage_df.drop_duplicates(subset=["gene_name"], keep="first")

    keep_cols = ["gene_name", "cage_tss", "distance_to_gene_tss", "chrom", "strand"]
    keep_cols = [c for c in keep_cols if c in cage_df.columns]
    cage_per_gene = cage_df[keep_cols].copy()
    cage_per_gene = cage_per_gene.rename(
        columns={"chrom": "cage_chrom", "strand": "cage_strand"}
    )

    merged_df = principal_df.merge(cage_per_gene, on="gene_name", how="left")
    if "tss_gencode" not in merged_df.columns and "tss" in merged_df.columns:
        merged_df["tss_gencode"] = merged_df["tss"]

    merged_df["cage_tss"] = pd.to_numeric(merged_df.get("cage_tss"), errors="coerce")
    merged_df["tss_gencode"] = pd.to_numeric(merged_df["tss_gencode"], errors="coerce")
    merged_df["tss"] = merged_df["cage_tss"].fillna(merged_df["tss_gencode"]).astype("Int64")

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    merged_df.to_csv(out_csv, sep="\t", index=False)
    return merged_df


def parse_gencode_genes(gtf_file, gene_type="protein_coding", min_length=1000):
    """Parse gene records from a GENCODE GTF."""
    genes = []
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            if f'gene_type "{gene_type}"' not in attrs:
                continue
            chrom = fields[0]
            gstart = int(fields[3]) - 1
            gend = int(fields[4])
            strand = fields[6]
            if gend - gstart < min_length:
                continue
            gene_id = ""
            gene_name = ""
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_id"):
                    gene_id = attr.split('"')[1]
                if attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1]
            tss = gstart if strand == "+" else gend
            tes = gend if strand == "+" else gstart
            genes.append(
                dict(
                    gene_id=gene_id,
                    gene_name=gene_name,
                    chrom=chrom,
                    strand=strand,
                    gstart=gstart,
                    gend=gend,
                    tss=tss,
                    tes=tes,
                    gene_length=gend - gstart,
                )
            )
    return pd.DataFrame(genes)


def load_genes_df(gtf_file=GENCODE_GTF, cache_path=GENES_DF_CACHE, genes_csv=None):
    """Load gene annotations from user CSV or cached/parsing GTF fallback."""
    if genes_csv is not None:
        print(f"Using user genes table: {genes_csv}")
        return _read_table_auto(genes_csv)

    # Default behavior: always use full-gene principal+CAGE fallback table.
    if os.path.exists(DEFAULT_PRINCIPAL_WITH_CAGE_CSV):
        print(f"Using default genes table: {DEFAULT_PRINCIPAL_WITH_CAGE_CSV}")
        return _read_table_auto(DEFAULT_PRINCIPAL_WITH_CAGE_CSV)
    if os.path.exists(PAUSING_INDEX_PRINCIPAL_CSV) and os.path.exists(CLOSEST_CAGE_TSS_CSV):
        print(f"Building default genes table: {DEFAULT_PRINCIPAL_WITH_CAGE_CSV}")
        return build_principal_with_cage_fallback_df()

    if os.path.exists(cache_path):
        return pd.read_csv(cache_path)
    genes_df = parse_gencode_genes(gtf_file)
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    genes_df.to_csv(cache_path, index=False)
    return genes_df


def _validate_gene_table(genes_df, tss_column="tss", gene_column="gene_name"):
    """Validate required columns and selected TSS/gene columns."""
    required = {"chrom", "strand", gene_column, tss_column}
    missing = [c for c in required if c not in genes_df.columns]
    if missing:
        raise ValueError(
            f"Missing required column(s): {missing}. "
            f"Available columns: {list(genes_df.columns)}"
        )


def find_gene(genes_df, gene_name, gene_column="gene_name"):
    """Find a gene by exact or case-insensitive match in gene_column."""
    hit = genes_df[genes_df[gene_column] == gene_name]
    if len(hit) == 0:
        hit = genes_df[
            genes_df[gene_column].astype(str).str.upper() == str(gene_name).upper()
        ]
    if len(hit) == 0:
        raise ValueError(
            f"Gene '{gene_name}' not found in column '{gene_column}'. "
            f"Check spelling or browse that column in your gene table."
        )
    primary = hit[hit.chrom.str.match(r"^chr[0-9XY]+$")]
    return primary.iloc[0] if len(primary) else hit.iloc[0]


def _resolve_tss(gene, requested_tss_column):
    """
    Resolve TSS position and report source.

    Behavior:
    - If requested column is cage_tss and value is missing, fall back to tss_gencode.
    - Otherwise use requested column directly.
    """
    if requested_tss_column == "cage_tss":
        cage_val = pd.to_numeric(pd.Series([gene.get("cage_tss")]), errors="coerce").iloc[0]
        if pd.notna(cage_val):
            return int(cage_val), "cage_tss"
        if "tss_gencode" in gene:
            gencode_val = pd.to_numeric(pd.Series([gene.get("tss_gencode")]), errors="coerce").iloc[0]
            if pd.notna(gencode_val):
                return int(gencode_val), "tss_gencode (fallback)"
        raise ValueError(
            "Requested tss_column='cage_tss', but both cage_tss and tss_gencode are missing for this gene."
        )

    tss_val = pd.to_numeric(pd.Series([gene.get(requested_tss_column)]), errors="coerce").iloc[0]
    if pd.isna(tss_val):
        raise ValueError(
            f"TSS value is missing for column '{requested_tss_column}' in selected gene row."
        )
    return int(tss_val), requested_tss_column


def classify_by_size(size):
    if 40 <= size <= 60:
        return "PPP"
    if 60 < size <= 80:
        return "PIC"
    if 140 <= size <= 160:
        return "nucleosome (140-160)"
    return "unknown"


def collect_reads(locus_chrom, view_start, view_end, show_fire_nuc=True):
    """Collect footprints per read in the requested window.

    show_fire_nuc=True  → include orange FIRE nucleosome marks (parsed nuc files);
                          FIRE_combined_footprints per-read marks are NOT added.
    show_fire_nuc=False → include FIRE_combined_footprints per-read marks (80-140 bp);
                          orange FIRE nucleosome marks are NOT added.
    """
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

    if show_fire_nuc:
        fire_nuc_path = os.path.join(
            FIRE_NUC_DIR, "sorted_indexed",
            f"all_samples_nuc_{locus_chrom}_parsed_sorted.bed.gz",
        )
        if os.path.exists(fire_nuc_path):
            tbx = pysam.TabixFile(fire_nuc_path)
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
    else:
        for label, fp_file in FIRE_FP_FILES.items():
            cls = f"FIRE {label}"
            tbx = pysam.TabixFile(fp_file)
            try:
                for row in tbx.fetch(locus_chrom, view_start, view_end):
                    f = row.split("\t")
                    fp_start, fp_end, read_name = int(f[1]), int(f[2]), f[3]
                    reads[read_name].append(
                        (fp_start, fp_end, cls, FP_COLORS[cls])
                    )
            except (ValueError, KeyError):
                pass
            tbx.close()
    return reads


def _normalize_fp_class(fp_class):
    key = str(fp_class).strip()
    normalized = CLASS_ALIASES.get(key.lower())
    if normalized is None:
        choices = ", ".join(["all"] + list(CLASS_SUBDIRS))
        raise ValueError(f"Unknown fp_class '{fp_class}'. Choose one of: {choices}")
    return normalized


def _summary_classes(fp_class):
    fp_class = _normalize_fp_class(fp_class)
    if fp_class == "all":
        return list(CLASS_SUBDIRS)
    return [fp_class]


def collect_reads_from_summaries(
    locus_chrom,
    view_start,
    view_end,
    fp_class="FIRE_nucleosome",
    summaries_dir=DEFAULT_FP_SUMMARY_DIR,
):
    """Collect precomputed summary footprints per read for a class or all classes."""
    reads = defaultdict(list)
    for cls in _summary_classes(fp_class):
        path = os.path.join(
            summaries_dir,
            CLASS_SUBDIRS[cls],
            f"{locus_chrom}_footprints.tsv.gz",
        )
        if not os.path.exists(path):
            print(f"  [WARN] missing summary file for {cls}: {path}")
            continue

        chunks = pd.read_csv(
            path,
            sep="\t",
            usecols=["read_name", "fp_start", "fp_end", "fp_class"],
            chunksize=500_000,
        )
        color = FP_COLORS[cls]
        for df in chunks:
            df = df[(df["fp_end"] >= view_start) & (df["fp_start"] <= view_end)]
            for row in df.itertuples(index=False):
                reads[row.read_name].append(
                    (int(row.fp_start), int(row.fp_end), cls, color)
                )
    return reads


def collect_reads_from_indexed_beds(
    locus_chrom,
    view_start,
    view_end,
    fp_class="FIRE_nucleosome",
    beds_dir=DEFAULT_FP_BED_DIR,
):
    """Collect precomputed summary footprints per read using tabix-indexed BED files.

    BED columns: chrom, fp_start, fp_end, read_name, fp_class, fp_mid
    """
    reads = defaultdict(list)
    for cls in _summary_classes(fp_class):
        path = os.path.join(
            beds_dir,
            CLASS_SUBDIRS[cls],
            f"{locus_chrom}_footprints.bed.gz",
        )
        if not os.path.exists(path):
            print(f"  [WARN] missing indexed BED for {cls}: {path}")
            continue
        color = FP_COLORS[cls]
        tbx = pysam.TabixFile(path)
        try:
            for row in tbx.fetch(locus_chrom, view_start, view_end):
                f = row.split("\t")
                fp_start, fp_end, read_name = int(f[1]), int(f[2]), f[3]
                reads[read_name].append((fp_start, fp_end, cls, color))
        except (ValueError, KeyError):
            pass
        tbx.close()
    return reads


def fetch_bw(bw_path, chrom, start, end, n_bins=N_BINS):
    try:
        bw = pyBigWig.open(bw_path)
        vals = bw.stats(chrom, start, end, type="mean", nBins=n_bins)
        bw.close()
        vals = np.array([v if v is not None else 0.0 for v in vals])
        return np.linspace(start, end, n_bins), vals
    except Exception:
        return np.linspace(start, end, n_bins), np.zeros(n_bins)


def _parse_length_range(length_range):
    """Normalize length range input to (min_len, max_len), or None."""
    if length_range is None:
        return None
    if isinstance(length_range, str):
        value = length_range.strip().lower()
        if value in {"", "none", "off", "false"}:
            return None
        if "-" not in value:
            raise ValueError("Length range must be formatted as 'MIN-MAX', e.g. '140-160'.")
        min_len, max_len = value.split("-", 1)
        length_range = (int(min_len), int(max_len))
    if len(length_range) != 2:
        raise ValueError("Length range must contain exactly two values: min and max.")
    min_len, max_len = int(length_range[0]), int(length_range[1])
    if min_len < 0 or max_len < min_len:
        raise ValueError("Length range must satisfy 0 <= min <= max.")
    return min_len, max_len


def _filter_fire_nucleosome_by_length(reads, length_range):
    """Drop FIRE_nucleosome footprints outside the requested length range."""
    length_range = _parse_length_range(length_range)
    if length_range is None:
        return reads

    min_len, max_len = length_range
    filtered_reads = defaultdict(list)
    for read_name, footprints in reads.items():
        for fp_start, fp_end, cls, color in footprints:
            if cls == "FIRE_nucleosome":
                fp_len = fp_end - fp_start
                if fp_len < min_len or fp_len > max_len:
                    continue
            filtered_reads[read_name].append((fp_start, fp_end, cls, color))
    return filtered_reads


def _relative_to_tss(pos, tss_pos, strand):
    """Return strand-aware TSS-relative position."""
    if strand == "+":
        return pos - tss_pos
    if strand == "-":
        return tss_pos - pos
    raise ValueError("strand must be '+' or '-'.")


def _relative_interval_to_tss(fp_start, fp_end, tss_pos, strand):
    """Return a footprint interval in strand-aware TSS-relative coordinates."""
    if strand == "+":
        return fp_start - tss_pos, fp_end - tss_pos
    if strand == "-":
        return tss_pos - fp_end, tss_pos - fp_start
    raise ValueError("strand must be '+' or '-'.")


def _nucleosome_distance_order_key(footprints, sort_cls, tss_pos, strand):
    """Sort by the TSS-spanning gap between -1 and +1 nucleosome edges."""
    if sort_cls in NUCLEOSOME_SORT_CLASSES:
        sort_classes = {sort_cls}
    else:
        sort_classes = NUCLEOSOME_SORT_CLASSES

    upstream = []
    downstream = []
    for fp_start, fp_end, cls, _ in footprints:
        if cls not in sort_classes:
            continue
        rel_start, rel_end = _relative_interval_to_tss(
            fp_start,
            fp_end,
            tss_pos,
            strand,
        )
        if rel_end <= 0:
            upstream.append(rel_end)
        elif rel_start >= 0:
            downstream.append(rel_start)

    if not upstream or not downstream:
        return (float("inf"), float("inf"), float("inf"))

    minus1_tss_edge = max(upstream)
    plus1_tss_edge = min(downstream)
    gap = plus1_tss_edge - minus1_tss_edge
    return (gap, plus1_tss_edge, abs(minus1_tss_edge))


def _footprint_order_key(
    footprints,
    sort_cls,
    tss_pos,
    plot_mode="midpoint",
    strand="+",
    sort_mode="nearest",
):
    """Sort reads by selected class midpoint/edge relative to the TSS.

    Uses footprints of sort_cls only unless sort_cls='all'; reads with no matching
    footprint sort last.
    """
    if sort_mode == "nucleosome_distance":
        return _nucleosome_distance_order_key(footprints, sort_cls, tss_pos, strand)

    candidates = []
    for fp_start, fp_end, cls, _ in footprints:
        if sort_cls != "all" and cls != sort_cls:
            continue
        if plot_mode == "midpoint":
            anchor = (fp_start + fp_end) / 2
            dist = abs(anchor - tss_pos)
        elif plot_mode == "edge":
            anchor = fp_end if strand == "+" else fp_start
            dist = abs(anchor - tss_pos)
        else:
            raise ValueError("plot_mode must be 'midpoint' or 'edge'.")

        if sort_mode == "nearest":
            candidates.append((dist, anchor))
        elif sort_mode == "plus1":
            anchor_rel = _relative_to_tss(anchor, tss_pos, strand)
            if anchor_rel > 0:
                candidates.append((anchor_rel, anchor))
        else:
            raise ValueError(
                "sort_mode must be 'nearest', 'plus1', or 'nucleosome_distance'."
            )
    if candidates:
        return min(candidates)
    return (float("inf"), min(fp[0] for fp in footprints))


def plot_gene_footprints(
    gene_name="ACTB",
    window_size=2000,
    max_reads=100,
    min_window_span_fraction=0.8,
    out_dir=DEFAULT_OUT_DIR,
    genes_df=None,
    genes_csv=None,
    tss_column="tss",
    gene_column="gene_name",
    show=True,
    show_fire_nuc=True,
    use_precomputed_summaries=True,
    fp_class="FIRE_nucleosome",
    plot_mode="midpoint",
    sort_mode="nearest",
    fire_nuc_length_range=DEFAULT_FIRE_NUC_LENGTH_RANGE,
    fp_bed_dir=DEFAULT_FP_BED_DIR,
    save_pdf=False,
):
    """Create the gene-centric Pol II footprint plot and save it as PDF."""
    if genes_df is None:
        genes_df = load_genes_df(genes_csv=genes_csv)

    _validate_gene_table(genes_df, tss_column=tss_column, gene_column=gene_column)
    gene = find_gene(genes_df, gene_name, gene_column=gene_column)
    locus_chrom = gene.chrom
    locus_center, tss_source = _resolve_tss(gene, tss_column)
    locus_strand = gene.strand
    view_start = max(0, locus_center - window_size)
    view_end = locus_center + window_size

    gene_id = gene["gene_id"] if "gene_id" in gene else "NA"
    print(f"Gene   : {gene_name}  ({gene_id})")
    print(f"Locus  : {locus_chrom}:{view_start:,}-{view_end:,}  strand={locus_strand}")
    print(f"TSS    : {locus_center:,}")
    print(f"TSS col: {tss_column}")
    print(f"TSS src: {tss_source}")
    print(f"Window : {view_end - view_start:,} bp")
    print(f"FP data: {'indexed BED summaries' if use_precomputed_summaries else 'tabix footprint beds'}")
    if use_precomputed_summaries:
        fp_class = _normalize_fp_class(fp_class)
        print(f"Class  : {fp_class}")
        print(f"Mode   : {plot_mode}")
        print(f"Sort   : {sort_mode}")
        print(f"FIRE nuc length: {_parse_length_range(fire_nuc_length_range)}")

    if not (0 < min_window_span_fraction <= 1):
        raise ValueError("min_window_span_fraction must be in the interval (0, 1].")

    if plot_mode not in {"midpoint", "edge"}:
        raise ValueError("plot_mode must be 'midpoint' or 'edge'.")
    if sort_mode not in {"nearest", "plus1", "nucleosome_distance"}:
        raise ValueError(
            "sort_mode must be 'nearest', 'plus1', or 'nucleosome_distance'."
        )

    if use_precomputed_summaries:
        read_spans = collect_reads_from_indexed_beds(
            locus_chrom,
            view_start,
            view_end,
            fp_class="all",
            beds_dir=fp_bed_dir,
        )
        reads = collect_reads_from_indexed_beds(
            locus_chrom,
            view_start,
            view_end,
            fp_class=fp_class,
            beds_dir=fp_bed_dir,
        )
    else:
        reads = collect_reads(locus_chrom, view_start, view_end, show_fire_nuc=show_fire_nuc)
        read_spans = reads
    read_spans = _filter_fire_nucleosome_by_length(read_spans, fire_nuc_length_range)
    if reads is read_spans:
        reads = read_spans
    else:
        reads = _filter_fire_nucleosome_by_length(reads, fire_nuc_length_range)
    print(list(reads.items())[:5])

    span_margin = int(round(window_size * (1 - min_window_span_fraction)))
    min_left = view_start + span_margin
    max_right = view_end - span_margin

    span_filtered_reads = {}
    for read_name, footprints in reads.items():
        span_footprints = read_spans.get(read_name, footprints)
        r_start = min(fp[0] for fp in span_footprints)
        r_end = max(fp[1] for fp in span_footprints)
        if r_start <= min_left and r_end >= max_right:
            # Store all-class footprints for display; sorting uses sort_cls filter internally.
            span_filtered_reads[read_name] = read_spans.get(read_name, footprints)

    sort_cls = fp_class if use_precomputed_summaries else (
        "FIRE_nucleosome" if show_fire_nuc else "nucleosome (140-160)"
    )
    read_list = sorted(
        span_filtered_reads.items(),
        key=lambda kv: _footprint_order_key(
            kv[1],
            sort_cls,
            locus_center,
            plot_mode,
            strand=locus_strand,
            sort_mode=sort_mode,
        ),
    )
    read_list = read_list[:max_reads]

    print(f"Reads with {fp_class} footprints in window: {len(reads):,}")
    print(
        f"Reads spanning >={min_window_span_fraction:.0%} of window "
        f"(start<={min_left:,}, end>={max_right:,}): {len(span_filtered_reads):,}"
    )
    print(f"Displaying: {len(read_list)}")

    by_cls = {}
    for _, fps in reads.items():
        for _, _, cls, _ in fps:
            by_cls[cls] = by_cls.get(cls, 0) + 1
    for k, v in sorted(by_cls.items()):
        print(f"  {k:12s}: {v:,} footprints")

    chr_track_dir = f"{TRACKS_DIR}/{locus_chrom}/"

    # (slab, color, display_label, bw_path) — order determines y-position in top panel
    track_defs = [
        (s, TRACK_COLORS[s], TRACK_LABELS[s],
         f"{chr_track_dir}/combined_{locus_chrom}_{s}bp_fp_cov.bw")
        for s in TRACK_SIZES
    ]
    if show_fire_nuc:
        track_defs.append((
            "FIRE_nuc",
            FP_COLORS["FIRE_nucleosome"],
            "FIRE nucleosome",
            f"{FIRE_NUC_BW_DIR}/all_samples_nuc_{locus_chrom}_parsed_sorted.bw",
        ))
    else:
        fire_fp_chr_dir = f"{FIRE_FP_TRACKS_DIR}/{locus_chrom}"
        track_defs += [
            (s, FIRE_FP_TRACK_COLORS[s], FIRE_FP_TRACK_LABELS[s],
             f"{fire_fp_chr_dir}/combined_{locus_chrom}_{s}bp_fp_cov.bw")
            for s in FIRE_FP_SIZES
        ]

    n_reads = len(read_list)
    n_top = len(track_defs)
    fig = plt.figure(figsize=(14, 3.5 + 0.10 * n_reads))
    gs = plt.GridSpec(
        2,
        1,
        height_ratios=[0.45 * n_top, max(2.5, 0.10 * n_reads)],
        hspace=0.04,
    )
    ax_top = fig.add_subplot(gs[0])
    ax_bot = fig.add_subplot(gs[1], sharex=ax_top)

    for i, (_, color, label, bw_path) in enumerate(track_defs):
        pos, vals = fetch_bw(bw_path, locus_chrom, view_start, view_end)
        vmax = vals.max() or 1
        offset = i * 1.1
        ax_top.fill_between(
            pos,
            offset,
            offset + vals / vmax,
            color=color,
            alpha=0.85,
            label=label,
        )

    ax_top.axvline(locus_center, color="black", lw=1, ls="--", alpha=0.5)
    ax_top.set_yticks(np.arange(n_top) * 1.1 + 0.5)
    ax_top.set_yticklabels([td[2] for td in track_defs], fontsize=8)
    ax_top.set_ylabel("Footprint size")
    ax_top.set_title(
        f"{gene_name}  |  {locus_chrom}:{view_start:,}-{view_end:,}  ({locus_strand})",
        fontsize=11,
    )
    ax_top.legend(loc="lower right", fontsize=7, ncol=2)

    # Explicit y-mapping so the first sorted read is always drawn on top.
    for idx, (read_name, footprints) in enumerate(read_list):
        y = n_reads - 1 - idx
        span_footprints = read_spans.get(read_name, footprints)
        r_start = min(fp[0] for fp in span_footprints)
        r_end = max(fp[1] for fp in span_footprints)
        ax_bot.plot([r_start, r_end], [y, y], color="black", lw=0.4, zorder=1)
        for fp_start, fp_end, cls, color in footprints:
            ax_bot.add_patch(
                plt.Rectangle(
                    (fp_start, y - 0.35),
                    fp_end - fp_start,
                    0.70,
                    linewidth=0,
                    facecolor=color,
                    alpha=0.85,
                    zorder=2,
                )
            )

    ax_bot.axvline(locus_center, color="black", lw=1, ls="--", alpha=0.5)
    ax_bot.set_xlim(view_start, view_end)
    ax_bot.set_ylim(-1, n_reads)
    ax_bot.set_yticks([])
    ax_bot.set_ylabel(f"Reads  (n={n_reads})", fontsize=9)
    ax_bot.set_xlabel(f"Genomic position  ({locus_chrom})", fontsize=9)

    used_cls = {cls for fps in reads.values() for _, _, cls, _ in fps}
    legend_patches = [
        plt.matplotlib.patches.Patch(color=FP_COLORS[cls], label=cls)
        for cls in FP_COLORS
        if cls in used_cls
    ]
    ax_bot.legend(handles=legend_patches, fontsize=7, loc="lower right", ncol=2)

    plt.setp(ax_top.get_xticklabels(), visible=False)
    plt.tight_layout()

    if save_pdf:
        os.makedirs(out_dir, exist_ok=True)
        if use_precomputed_summaries:
            class_tag = CLASS_SUBDIRS.get(fp_class, "all")
            out_base = f"{out_dir}/{gene_name}_{class_tag}_{plot_mode}_{sort_mode}_footprints"
        else:
            out_base = f"{out_dir}/{gene_name}_footprints"
        plt.savefig(out_base + ".pdf", bbox_inches="tight")
        print(f"Saved: {out_base}.pdf")
    else:
        out_base = None

    if show:
        plt.show()
    else:
        plt.close(fig)

    return {
        "gene": gene,
        "out_base": out_base,
        "reads": reads,
        "read_list": read_list,
        "view_start": view_start,
        "view_end": view_end,
        "locus_chrom": locus_chrom,
        "locus_center": locus_center,
        "read_spans": read_spans,
    }


def build_argparser():
    p = argparse.ArgumentParser(description="Plot Pol II footprints around a gene.")
    p.add_argument("--gene", default="ACTB", help="Gene name to plot.")
    p.add_argument(
        "--genes-csv",
        default=None,
        help=(
            "Optional user-provided gene table. "
            "If omitted, script uses/saves principal_with_CAGE_fallback table."
        ),
    )
    p.add_argument(
        "--tss-column",
        default="tss",
        help="Column name in the selected gene table to use as TSS (default: tss).",
    )
    p.add_argument(
        "--gene-column",
        default="gene_name",
        help="Column name used to match --gene (default: gene_name).",
    )
    p.add_argument("--window-size", type=int, default=2000, help="Bases to show on each side of TSS.")
    p.add_argument("--max-reads", type=int, default=100, help="Maximum reads to display.")
    p.add_argument(
        "--min-window-span-fraction",
        type=float,
        default=0.8,
        help=(
            "Minimum fraction of the plotted window a read must span to be shown "
            "(default: 0.8)."
        ),
    )
    p.add_argument(
        "--out-dir",
        default=DEFAULT_OUT_DIR,
        help="Directory to write PDF output.",
    )
    p.add_argument(
        "--no-show",
        action="store_true",
        help="Save the plot without opening an interactive window.",
    )
    p.add_argument(
        "--no-fire-nuc",
        action="store_true",
        help=(
            "Hide orange FIRE nucleosome marks and show FIRE_combined_footprints "
            "(80-140 bp) bigwig tracks and per-read marks instead."
        ),
    )
    p.add_argument(
        "--fp-bed-dir",
        default=DEFAULT_FP_BED_DIR,
        help=(
            "Root directory of tabix-indexed BED summary files from "
            f"index_fp_summaries_as_bed.py (default: {DEFAULT_FP_BED_DIR})."
        ),
    )
    p.add_argument(
        "--fp-class",
        default="FIRE_nucleosome",
        help=(
            "Footprint class to plot/sort from precomputed summaries. "
            "Accepted values: PPP, PIC, nucleosome_140-160, FIRE_nucleosome, "
            "unknown, all. Default: FIRE_nucleosome."
        ),
    )
    p.add_argument(
        "--plot-mode",
        choices=["midpoint", "edge"],
        default="midpoint",
        help=(
            "Sort selected reads by the selected class footprint midpoint or "
            "strand-specific edge relative to the TSS. For edge mode, use the right "
            "edge on + strand genes and the left edge on - strand genes "
            "(default: midpoint)."
        ),
    )
    p.add_argument(
        "--sort-mode",
        choices=["nearest", "plus1", "nucleosome_distance"],
        default="nearest",
        help=(
            "Sort selected reads by the nearest selected-class footprint on either side "
            "of the TSS, or by the nearest downstream/+1 selected-class footprint only "
            "using strand-aware TSS-relative coordinates, or by the distance between "
            "the selected-class -1 and +1 nucleosome edges across the TSS "
            "(default: nearest)."
        ),
    )
    p.add_argument(
        "--fire-nuc-length-range",
        default="140-160",
        help=(
            "Keep FIRE_nucleosome footprints only within this length range, formatted "
            "as MIN-MAX bp. Use 'none' to disable this filter (default: 140-160)."
        ),
    )
    p.add_argument(
        "--save-pdf",
        action="store_true",
        help="Save the plot as a PDF file (default: do not save).",
    )
    p.add_argument(
        "--use-tabix-beds",
        action="store_true",
        help=(
            "Use the original tabix BED collection path instead of indexed BED summaries. "
            "--fp-class and --fp-bed-dir are ignored in this mode."
        ),
    )
    return p


def main(argv=None):
    args = build_argparser().parse_args(argv)
    plot_gene_footprints(
        gene_name=args.gene,
        window_size=args.window_size,
        max_reads=args.max_reads,
        min_window_span_fraction=args.min_window_span_fraction,
        out_dir=args.out_dir,
        genes_df=None,
        genes_csv=args.genes_csv,
        tss_column=args.tss_column,
        gene_column=args.gene_column,
        show=not args.no_show,
        show_fire_nuc=not args.no_fire_nuc,
        use_precomputed_summaries=not args.use_tabix_beds,
        fp_class=args.fp_class,
        plot_mode=args.plot_mode,
        sort_mode=args.sort_mode,
        fire_nuc_length_range=args.fire_nuc_length_range,
        fp_bed_dir=args.fp_bed_dir,
        save_pdf=args.save_pdf,
    )


if __name__ == "__main__":
    main()
