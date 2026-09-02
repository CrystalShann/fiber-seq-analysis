#!/usr/bin/env python3
"""m6A (SAM-seq accessibility) metaprofiles around protein-coding TSS,
stratified by expression bin, per LPS timepoint.

For every gene in tss_expression_bins.tsv (from 01_expression_bins.R) a
TSS +/- WINDOW window is tiled into BIN_SIZE bins on the STRAND-AWARE axis
(rel = pos - tss on '+', tss - pos on '-'; so upstream is always negative and
the sense direction always points right). Per (window, bin) and timepoint:

    n_sites   reference A/T positions in the bin (m6A-callable sites)
    n_cov     read x site events: sum over the bin's A/T sites of aligned-span
              depth at that site
    n_mod     called m6A read x site events (calls landing on a reference A/T)
    m6a       n_mod / n_cov

i.e. the same accessibility definition as 04_m6a_binned.py; parsing follows it
too: the first and last BED12 blocks are ft extract sentinels (dropped), split
reads are deduped keeping the longest alignment per read name, and calls not on
a reference A/T are counted off-site and excluded. Reads are pulled through the
tabix index over merged TSS windows (two fetch passes, mirroring the two zcat
passes of 04), so the dedup universe is the promoter-overlapping reads only.

Outputs in --out-dir:
    tss_m6a_profile_by_expr_bin.tsv.gz  one row per (timepoint, expr_bin, bin)
    tss_m6a_gene_totals.tsv.gz          one row per (gene, timepoint), whole
                                        window, for QC / per-gene scatter
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

warnings.filterwarnings("ignore", category=DeprecationWarning)

SAMPLES = ["LPS_0", "LPS_5", "LPS_10", "LPS_15"]
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

FT_ROOT = "/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/ft_result_dir"
REF_FA = "/project/spott/reference/human/GRCh38/hg38.fa"

WINDOW = 1000          # bp each side of the TSS
BIN_SIZE = 10
NBINS = 2 * WINDOW // BIN_SIZE   # rel in [-WINDOW, WINDOW)


def m6a_path(ft_root, sample, chrom):
    return Path(ft_root) / sample / "extracted_results" / "m6a_by_chr" / \
        f"{sample}.ft_extracted_m6a.{chrom}.bed.gz"


def merge_intervals(starts, ends):
    """Merge sorted, possibly overlapping [start, end) windows for fetching."""
    out = []
    for s, e in zip(starts, ends):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


###############################
# Per-chromosome window / A/T-site arrays (shared by all timepoints)
###############################
def build_windows(sub, chrom_seq):
    """Per-window A/T sites with their strand-aware bin. Windows may overlap
    (bidirectional promoters), so sites are duplicated per window; within one
    window they stay sorted by genomic position."""
    ws = (sub["tss"] - WINDOW).to_numpy()
    we = (sub["tss"] + WINDOW + 1).to_numpy()   # fetch extent; bins trim edges
    tss = sub["tss"].to_numpy()
    minus = (sub["strand"] == "-").to_numpy()

    sites_l, k_l = [], []
    for r in range(len(sub)):
        seg = np.frombuffer(chrom_seq[ws[r]:we[r]].encode(), dtype=np.uint8)
        at = ws[r] + np.nonzero((seg == ord("A")) | (seg == ord("T")))[0]
        rel = (tss[r] - at) if minus[r] else (at - tss[r])
        k = (rel + WINDOW) // BIN_SIZE
        keep = (k >= 0) & (k < NBINS)
        sites_l.append(at[keep])
        k_l.append(k[keep])

    nsite = np.array([a.size for a in sites_l], dtype=np.int64)
    b = {
        "ws": ws, "we": we,
        "sites": np.concatenate(sites_l) if len(sites_l) else np.empty(0, np.int64),
        "site_k": np.concatenate(k_l) if len(k_l) else np.empty(0, np.int64),
        "offsets": np.concatenate([[0], np.cumsum(nsite)]),
    }
    b["site_win"] = np.repeat(np.arange(len(sub), dtype=np.int64), nsite)
    b["site_glob"] = b["site_win"] * NBINS + b["site_k"]
    b["n_sites"] = np.bincount(b["site_glob"],
                               minlength=len(sub) * NBINS).astype(np.int64)
    return b


###############################
# Pass 1: longest alignment per read name over the merged windows
###############################
def survey_reads(tb, chrom, merged):
    best = {}
    for s, e in merged:
        for line in tb.fetch(chrom, s, e):
            f = line.split("\t", 4)
            st = int(f[1])
            ln = int(f[2]) - st
            rid = f[3]
            cur = best.get(rid)
            if cur is None or ln > cur[0]:
                best[rid] = (ln, st)
    return best


###############################
# Pass 2: per-(window, bin) n_cov / n_mod for one (chromosome, timepoint)
###############################
def counts_for_chrom(tb, chrom, merged, best, b):
    nwin = len(b["ws"])
    cov_diff = np.zeros(b["sites"].size + 1, dtype=np.int64)
    n_mod = np.zeros(nwin * NBINS, dtype=np.int64)
    n_reads = np.zeros(nwin, dtype=np.int64)
    off_site = 0
    n_kept = 0

    ws, we, sites, offsets = b["ws"], b["we"], b["sites"], b["offsets"]
    site_k = b["site_k"]
    for ms, me in merged:
        for line in tb.fetch(chrom, ms, me):
            f = line.split("\t")
            st = int(f[1])
            en = int(f[2])
            rid = f[3]
            if best.get(rid) != (en - st, st):
                continue        # shorter split alignment, or already counted
            best[rid] = None
            n_kept += 1

            r0 = np.searchsorted(we, st, side="right")
            r1 = np.searchsorted(ws, en, side="left")
            if r1 <= r0:
                continue

            # m6A calls: drop the sentinel first/last blocks; 1 bp blocks, so
            # the modified base (0-based) is span start + blockStart (sorted)
            starts = np.fromstring(f[11], dtype=np.int64, sep=",")[1:-1]
            pos = st + starts

            for r in range(r0, r1):
                o0, o1 = offsets[r], offsets[r + 1]
                wsl = sites[o0:o1]
                # span-depth over the window's A/T sites: contiguous site run
                i0 = np.searchsorted(wsl, st, side="left")
                i1 = np.searchsorted(wsl, en, side="left")
                cov_diff[o0 + i0] += 1
                cov_diff[o0 + i1] -= 1
                n_reads[r] += 1

                a = np.searchsorted(pos, ws[r], side="left")
                c = np.searchsorted(pos, we[r], side="left")
                p = pos[a:c]
                if p.size:
                    idx = np.searchsorted(wsl, p)
                    on = (idx < wsl.size) & (wsl[np.minimum(idx, wsl.size - 1)] == p)
                    np.add.at(n_mod, r * NBINS + site_k[o0 + idx[on]], 1)
                    # in the window but not on a reference A/T, or trimmed off
                    # the strand-aware bin range at the window edge
                    off_site += int((~on).sum())

    site_depth = np.cumsum(cov_diff[:-1])
    n_cov = np.bincount(b["site_glob"], weights=site_depth,
                        minlength=nwin * NBINS).astype(np.int64)
    assert (n_mod <= n_cov).all()
    return n_cov, n_mod, n_reads, off_site, n_kept


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bins-tsv",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/expr_access/tables/tss_expression_bins.tsv")
    ap.add_argument("--ft-root", default=FT_ROOT)
    ap.add_argument("--ref", default=REF_FA)
    ap.add_argument("--out-dir",
                    default="/project/spott/cshan/fiber-seq/macrophage_project/expr_access/tables")
    ap.add_argument("--timepoints", nargs="+", default=SAMPLES)
    ap.add_argument("--chrom", default=None, help="restrict to one chromosome (testing)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(args.bins_tsv, sep="\t")
    chroms = [args.chrom] if args.chrom else CHROMS
    genes = genes[genes["chrom"].isin(chroms)].reset_index(drop=True)

    fa = pysam.FastaFile(args.ref)
    chrom_len = {c: fa.get_reference_length(c) for c in chroms}
    edge = (genes["tss"] < WINDOW) | \
           (genes["tss"] + WINDOW + 1 > genes["chrom"].map(chrom_len))
    if edge.any():
        print(f"dropping {edge.sum()} genes whose window leaves the chromosome",
              flush=True)
        genes = genes[~edge].reset_index(drop=True)
    print(f"{len(genes)} genes on {genes['chrom'].nunique()} chromosomes; "
          f"window +/-{WINDOW} bp, {NBINS} x {BIN_SIZE} bp bins", flush=True)

    for s in args.timepoints:
        for chrom in chroms:
            if not m6a_path(args.ft_root, s, chrom).is_file():
                sys.exit(f"ERROR: missing input {m6a_path(args.ft_root, s, chrom)}")

    expr_bins = sorted(genes["expr_bin"].unique())
    prof_cov = {s: np.zeros((len(expr_bins), NBINS), dtype=np.int64) for s in args.timepoints}
    prof_mod = {s: np.zeros((len(expr_bins), NBINS), dtype=np.int64) for s in args.timepoints}
    prof_sites = np.zeros((len(expr_bins), NBINS), dtype=np.int64)
    gene_rows = []
    stats = {s: {"off_site": 0, "n_kept": 0} for s in args.timepoints}

    for chrom in chroms:
        sub = genes[genes["chrom"] == chrom].sort_values("tss").reset_index(drop=True)
        if sub.empty:
            continue
        chrom_seq = fa.fetch(chrom).upper()
        b = build_windows(sub, chrom_seq)
        del chrom_seq
        merged = merge_intervals(b["ws"], b["we"])
        grp = np.searchsorted(expr_bins, sub["expr_bin"].to_numpy())
        sites_mat = b["n_sites"].reshape(len(sub), NBINS)
        np.add.at(prof_sites, grp, sites_mat)

        for s in args.timepoints:
            with pysam.TabixFile(str(m6a_path(args.ft_root, s, chrom))) as tb:
                best = survey_reads(tb, chrom, merged)
                n_cov, n_mod, n_reads, off, kept = counts_for_chrom(
                    tb, chrom, merged, best, b)
            stats[s]["off_site"] += off
            stats[s]["n_kept"] += kept
            print(f"{chrom} {s}: kept {kept} reads over {len(sub)} windows; "
                  f"off-site calls {off}", flush=True)

            cov_mat = n_cov.reshape(len(sub), NBINS)
            mod_mat = n_mod.reshape(len(sub), NBINS)
            np.add.at(prof_cov[s], grp, cov_mat)
            np.add.at(prof_mod[s], grp, mod_mat)

            gene_rows.append(pd.DataFrame({
                "gene_id": sub["gene_id"], "gene_name": sub["gene_name"],
                "chrom": chrom, "tss": sub["tss"], "strand": sub["strand"],
                "expr_bin": sub["expr_bin"], "mean_tpm": sub["mean_tpm"],
                "timepoint": s,
                "n_sites": sites_mat.sum(1), "n_cov": cov_mat.sum(1),
                "n_mod": mod_mat.sum(1), "n_reads": n_reads,
            }))

    n_genes = genes.groupby("expr_bin").size()
    bin_start = np.arange(NBINS) * BIN_SIZE - WINDOW
    prof = pd.concat([
        pd.DataFrame({
            "timepoint": s, "expr_bin": g,
            "bin_start": bin_start, "bin_end": bin_start + BIN_SIZE,
            "n_genes": n_genes[g],
            "n_sites": prof_sites[gi],
            "n_cov": prof_cov[s][gi], "n_mod": prof_mod[s][gi],
        })
        for s in args.timepoints for gi, g in enumerate(expr_bins)
    ], ignore_index=True)
    prof["m6a"] = prof["n_mod"] / prof["n_cov"].where(prof["n_cov"] > 0)
    assert prof["m6a"].dropna().between(0, 1).all()

    for s in args.timepoints:
        sl = prof[prof["timepoint"] == s]
        print(f"{s}: kept reads {stats[s]['n_kept']}; off-site calls "
              f"{stats[s]['off_site']}; mean m6a {sl['m6a'].mean():.4f}", flush=True)

    prof_out = out_dir / "tss_m6a_profile_by_expr_bin.tsv.gz"
    prof.to_csv(prof_out, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {prof_out} ({len(prof)} rows)", flush=True)

    gene_tot = pd.concat(gene_rows, ignore_index=True)
    gene_tot["m6a"] = gene_tot["n_mod"] / gene_tot["n_cov"].where(gene_tot["n_cov"] > 0)
    tot_out = out_dir / "tss_m6a_gene_totals.tsv.gz"
    gene_tot.to_csv(tot_out, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {tot_out} ({len(gene_tot)} rows)", flush=True)


if __name__ == "__main__":
    main()
