#!/usr/bin/env python3
"""Same-molecule co-accessibility of ENCODE cCRE pairs around gene TSSs, per LPS timepoint.

Reimplements coaccess_fire_CREs_combined_samples_around_genes.R and its helper
test_coaccess_fire_elements() from process_fiberseq_data.R. Nothing of Kevin's is
sourced; this is an independent implementation of the same definitions.

    elements        ENCODE cCREs overlapping a FIRE peak, cCRE coordinates kept
                    CRE_ID = accession1.accession2.CRE_label
    gene window     TSS +/- 10 kb, pairs enumerated within a window
    pair distance   GenomicRanges gap distance, MIN_DIST < d < MAX_DIST (500 / 20000)
    shared reads    reads overlapping BOTH cCREs (Kevin: two subsetByOverlaps + intersect
                    on read id, i.e. ANY overlap >= 1 bp) -- see --read-rule
    accessible      the read carries >= 1 FIRE element overlapping that cCRE (>= 1 bp)
    test            two-sided Fisher exact on the 2x2 + a pseudocount of 1 in every cell

The one structural difference from Kevin is that each timepoint is counted and tested
BY ITSELF, where he pools 17 LCL samples into a single test. Timepoints are never
tested against each other here.

2x2 layout (rows = element 1 state, cols = element 2 state, FALSE first), matching
his table(fire_region1, fire_region2) with levels forced to c(FALSE, TRUE):

                 elem2 FALSE      elem2 TRUE
    elem1 FALSE  co_closed        CRE2_access
    elem1 TRUE   CRE1_access      co_access

Inputs (from 01_make_cre_universe.sh and 02_read_spans.sh):
    <root>/universe/cre_universe.bed              chrom start end CRE_ID CRE_label
    <root>/universe/cre_gene_map.tsv.gz           CRE_ID -> gene windows
    <root>/universe/cre_in_timepoint_peaks.tsv.gz per-timepoint peak membership
    <root>/<s>/<s>.read_spans.bed.gz              aligned read spans (tabix)
    <FIRE>/<s>/additional-outputs-v0.1/fire-peaks/<s>-v0.1-fire-elements.bed.gz

Outputs in --out-dir:
    <s>_coaccess_stat.tsv.gz    one row per (gene, cCRE pair) - Kevin's coaccess_stat_df
    <s>_coaccess_pairs.tsv.gz   one row per distinct cCRE pair, with BH FDR
"""

import argparse
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio as conditional_odds_ratio

SAMPLES = ["LPS_0", "LPS_5", "LPS_10", "LPS_15"]
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

BEDTOOLS = "/project/spott/cshan/envs/bedtools/bin/bedtools"
TABIX = "/project/spott/cshan/envs/dimelo/bin/tabix"
# The accessibility numerator. This is the FIRE pipeline's own per-read element call
# (FDR <= 0.05), genome-wide and tabix-indexed - the same information Kevin reads out
# of fire.bed's FIRE-class segments, already computed. Symlinked as `lizarraga_FIRE`
# in the code directory.
#
# Deliberately NOT the numerator: FiberHMM/extract/ft_result_dir. That tree holds raw
# m6A / CpG / nucleosome calls with no FIRE scoring, and was extracted from the
# unfiltered BAM, so it covers ~8% more fibers than the -filtered FIRE CRAM. It is
# used for display only, by coaccess_plot_functions.R.
FIRE_ROOT = "/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE"


def gr_distance(s1, e1, s2, e2):
    """GenomicRanges::distance() on 0-based half-open intervals.

    The gap between the two ranges: 0 when they overlap or are bookended, otherwise
    the number of bases strictly between them. This is NOT a midpoint distance.
    """
    if s2 >= e1:
        return s2 - e1
    if s1 >= e2:
        return s1 - e2
    return 0


def build_pairs(cre, cre_gene, min_dist, max_dist):
    """All within-gene-window cCRE pairs passing the distance filter.

    Kevin enumerates expand.grid() per gene, drops self-pairs, and de-duplicates on a
    sorted key. His surviving orientation is an artifact of that construction, which
    makes CRE1/CRE2 - and therefore CRE1_access/CRE2_access - arbitrary. Here CRE1 is
    always the leftmost by (start, end, CRE_ID), so the columns are interpretable and
    the pair key is stable across timepoints.

    Returns a DataFrame with one row per (gene, pair).
    """
    pos = {r.CRE_ID: (r.chrom, r.start, r.end) for r in cre.itertuples()}
    rows = []
    for gene_id, grp in cre_gene.groupby("gene_id", sort=False):
        ids = grp["CRE_ID"].unique()
        if len(ids) < 2:
            continue
        # sort by coordinate so pair orientation is deterministic
        ids = sorted(ids, key=lambda c: (pos[c][1], pos[c][2], c))
        g0 = grp.iloc[0]
        for i in range(len(ids)):
            c1 = ids[i]
            _, s1, e1 = pos[c1]
            for j in range(i + 1, len(ids)):
                c2 = ids[j]
                _, s2, e2 = pos[c2]
                d = gr_distance(s1, e1, s2, e2)
                if min_dist < d < max_dist:
                    rows.append((g0.gene_id, g0.gene_name, g0.transcript_type,
                                 c1, c2, d))
    return pd.DataFrame(rows, columns=["gene_id", "gene_name", "transcript_type",
                                       "CRE1", "CRE2", "dist"])


def run_intersect(args_list, desc):
    p = subprocess.run(args_list, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"ERROR: {desc} failed:\n{p.stderr[:2000]}")
    return p.stdout


def incidence_for_chrom(chrom, cre_bed_path, spans_path, fe_path, read_rule, tmpdir):
    """{CRE_ID: set(read_id)} for covered and for accessible, on one chromosome."""
    cre_chr = tmpdir / f"cre.{chrom}.bed"
    with open(cre_chr, "w") as fh:
        n = 0
        for line in open(cre_bed_path):
            if line.split("\t", 1)[0] == chrom:
                fh.write(line)
                n += 1
    if n == 0:
        return {}, {}, 0

    spans_chr = tmpdir / f"spans.{chrom}.bed"
    with open(spans_chr, "w") as fh:
        p = subprocess.run([TABIX, str(spans_path), chrom], stdout=fh, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        sys.exit(f"ERROR: tabix on spans for {chrom} failed:\n{p.stderr[:2000]}")

    fe_chr = tmpdir / f"fe.{chrom}.bed"
    with open(fe_chr, "w") as fh:
        p = subprocess.run([TABIX, str(fe_path), chrom], stdout=fh, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        sys.exit(f"ERROR: tabix on fire elements for {chrom} failed:\n{p.stderr[:2000]}")

    # reads covering each cCRE. Kevin's rule is any overlap (subsetByOverlaps default
    # minoverlap = 1L); --read-rule contain requires the read to span the whole cCRE.
    cmd = [BEDTOOLS, "intersect", "-a", str(cre_chr), "-b", str(spans_chr), "-wa", "-wb"]
    if read_rule == "contain":
        cmd += ["-f", "1.0"]
    cov = defaultdict(set)
    for line in run_intersect(cmd, f"read intersect {chrom}").splitlines():
        f = line.split("\t")
        cov[f[3]].add(f[8])          # CRE_ID, read name

    # a read is accessible at a cCRE if one of ITS FIRE elements overlaps it
    acc = defaultdict(set)
    for line in run_intersect(
            [BEDTOOLS, "intersect", "-a", str(cre_chr), "-b", str(fe_chr), "-wa", "-wb"],
            f"fire element intersect {chrom}").splitlines():
        f = line.split("\t")
        acc[f[3]].add(f[8])          # CRE_ID, read name

    # Under --read-rule contain a FIRE element can overlap the cCRE on a read that does
    # not span all of it, so accessibility is not automatically a subset of coverage.
    n_orphan = 0
    for c in list(acc):
        before = len(acc[c])
        acc[c] &= cov.get(c, set())
        n_orphan += before - len(acc[c])

    for f in (cre_chr, spans_chr, fe_chr):
        f.unlink(missing_ok=True)
    return cov, acc, n_orphan


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="/project/spott/cshan/fiber-seq/macrophage_project/co-accessibility")
    ap.add_argument("--out-dir", default=None, help="default: <root>/coaccess")
    ap.add_argument("--timepoints", nargs="+", default=SAMPLES)
    ap.add_argument("--min-dist", type=int, default=500,
                    help="Kevin: min_dist = 500 (strict)")
    ap.add_argument("--max-dist", type=int, default=20000,
                    help="Kevin's combined-samples script: 20000 (strict). His "
                         "single-sample script uses 10000")
    ap.add_argument("--pseudocount", type=int, default=1,
                    help="added to every cell before the Fisher test, as Kevin does")
    ap.add_argument("--read-rule", choices=["any", "contain"], default="any",
                    help="'any' (default) reproduces Kevin's subsetByOverlaps; "
                         "'contain' requires the read to span the whole cCRE")
    ap.add_argument("--kevin-compat", action="store_true",
                    help="drop pairs where either cCRE is never accessible among the "
                         "shared reads, as Kevin's NA path silently does")
    ap.add_argument("--chrom", default=None, help="restrict to one chromosome (testing)")
    args = ap.parse_args()

    root = Path(args.root)
    out_dir = Path(args.out_dir) if args.out_dir else root / "coaccess"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmpdir = out_dir / ".tmp"
    tmpdir.mkdir(exist_ok=True)

    uni = root / "universe"
    cre_bed = uni / "cre_universe.bed"
    for f in (cre_bed, uni / "cre_gene_map.tsv.gz", uni / "cre_in_timepoint_peaks.tsv.gz"):
        if not f.is_file():
            sys.exit(f"ERROR: missing input {f} (run 01_make_cre_universe.sh first)")

    cre = pd.read_csv(cre_bed, sep="\t", header=None,
                      names=["chrom", "start", "end", "CRE_ID", "CRE_label"])
    cre_gene = pd.read_csv(uni / "cre_gene_map.tsv.gz", sep="\t")
    flags = pd.read_csv(uni / "cre_in_timepoint_peaks.tsv.gz", sep="\t").set_index("CRE_ID")

    chroms = [args.chrom] if args.chrom else CHROMS
    cre = cre[cre["chrom"].isin(chroms)].reset_index(drop=True)
    cre_gene = cre_gene[cre_gene["chrom"].isin(chroms)]
    print(f"cCRE universe: {len(cre)} elements on {cre['chrom'].nunique()} chromosomes",
          flush=True)

    pairs = build_pairs(cre, cre_gene, args.min_dist, args.max_dist)
    if pairs.empty:
        sys.exit("ERROR: no pairs after the gene-window and distance filters")
    meta = cre.set_index("CRE_ID")
    pairs["chrom"] = meta.loc[pairs["CRE1"], "chrom"].to_numpy()
    n_uniq_pairs = pairs.groupby(["CRE1", "CRE2"]).ngroups
    print(f"gene-anchored pairs: {len(pairs)}  ({n_uniq_pairs} distinct cCRE pairs) "
          f"with {args.min_dist} < gap < {args.max_dist}", flush=True)

    for s in args.timepoints:
        spans = root / s / f"{s}.read_spans.bed.gz"
        fe = Path(FIRE_ROOT) / s / "additional-outputs-v0.1" / "fire-peaks" / \
            f"{s}-v0.1-fire-elements.bed.gz"
        for f in (spans, fe):
            if not f.is_file():
                sys.exit(f"ERROR: missing input {f}")

        # Score one chromosome at a time and keep only the counts. Holding the
        # incidence for every chromosome at once would mean ~9M (cCRE, read name)
        # string entries genome-wide; per chromosome it stays a few hundred thousand.
        empty = set()
        counts = np.zeros((len(pairs), 4), dtype=np.int64)   # co_closed, CRE1, CRE2, co_access
        rows_by_chrom = pairs.groupby("chrom", sort=False).indices
        orphans = 0
        n_cre_seen = 0
        for chrom in chroms:
            idx = rows_by_chrom.get(chrom)
            if idx is None or len(idx) == 0:
                continue
            cov, acc, n_orph = incidence_for_chrom(chrom, cre_bed, spans, fe,
                                                   args.read_rule, tmpdir)
            orphans += n_orph
            n_cre_seen += len(cov)
            c1s = pairs["CRE1"].to_numpy()[idx]
            c2s = pairs["CRE2"].to_numpy()[idx]
            for k, i in enumerate(idx):
                shared = cov.get(c1s[k], empty) & cov.get(c2s[k], empty)
                if not shared:
                    continue
                a1 = acc.get(c1s[k], empty) & shared
                a2 = acc.get(c2s[k], empty) & shared
                both = len(a1 & a2)
                counts[i] = (len(shared) - len(a1) - len(a2) + both,   # neither
                             len(a1) - both,                            # first only
                             len(a2) - both,                            # second only
                             both)
            del cov, acc
        print(f"{s}: incidence built for {n_cre_seen} cCREs "
              f"({orphans} element hits dropped as not covered under --read-rule {args.read_rule})",
              flush=True)

        co_closed, cre1_access, cre2_access, co_access = counts.T
        n_shared = counts.sum(axis=1)

        # Fisher exact and the conditional MLE odds ratio on the +pseudocount table,
        # matching fisher.test(contingency_table + pseudocount). Identical tables share
        # one test, so only unique tables are computed.
        pc = args.pseudocount
        uniq, inv = np.unique(counts, axis=0, return_inverse=True)
        pvals = np.empty(len(uniq))
        fest = np.empty(len(uniq))
        praw = np.empty(len(uniq))
        for t, (ne, fi, se, bo) in enumerate(uniq):
            tab = [[int(ne) + pc, int(se) + pc], [int(fi) + pc, int(bo) + pc]]
            pvals[t] = fisher_exact(tab)[1]
            fest[t] = conditional_odds_ratio(tab, kind="conditional").statistic
            raw = [[int(ne), int(se)], [int(fi), int(bo)]]
            praw[t] = fisher_exact(raw)[1]

        df = pd.DataFrame({
            # --- Kevin's coaccess_stat_df columns, his names and order ---
            "gene_name": pairs["gene_name"].to_numpy(),
            "CRE1": pairs["CRE1"].to_numpy(),
            "CRE2": pairs["CRE2"].to_numpy(),
            "CRE_pair": pairs["CRE1"].to_numpy() + "-" + pairs["CRE2"].to_numpy(),
            "CRE_pair_labels": (meta.loc[pairs["CRE1"], "CRE_label"].to_numpy() + "-" +
                                meta.loc[pairs["CRE2"], "CRE_label"].to_numpy()),
            "dist": pairs["dist"].to_numpy(),
            "n_shared_reads": n_shared,
            "co_access": co_access,
            "co_closed": co_closed,
            "CRE1_access": cre1_access,
            "CRE2_access": cre2_access,
            "fisher_estimate": fest[inv],
            "pval": pvals[inv],
            "OR": ((co_access + pc) * (co_closed + pc)) /
                  ((cre1_access + pc) * (cre2_access + pc)),
            # --- added ---
            # Kevin names this cell co_closed when he builds coaccess_stat_df
            # (coaccess_...R:251) but filters on co_inaccess in example_figures.Rmd.
            # Both names are emitted so either of his snippets runs unchanged.
            "co_inaccess": co_closed,
            "timepoint": s,
            "gene_id": pairs["gene_id"].to_numpy(),
            "transcript_type": pairs["transcript_type"].to_numpy(),
            "chrom": pairs["chrom"].to_numpy(),
            "cre1_start": meta.loc[pairs["CRE1"], "start"].to_numpy(),
            "cre1_end": meta.loc[pairs["CRE1"], "end"].to_numpy(),
            "cre2_start": meta.loc[pairs["CRE2"], "start"].to_numpy(),
            "cre2_end": meta.loc[pairs["CRE2"], "end"].to_numpy(),
            "pval_raw": praw[inv],
            "or_haldane": ((co_access + 0.5) * (co_closed + 0.5)) /
                          ((cre1_access + 0.5) * (cre2_access + 0.5)),
            "zero_access_cre1": (co_access + cre1_access) == 0,
            "zero_access_cre2": (co_access + cre2_access) == 0,
            "cre1_in_tp_peaks": flags.loc[pairs["CRE1"], f"in_{s}_peaks"].to_numpy(),
            "cre2_in_tp_peaks": flags.loc[pairs["CRE2"], f"in_{s}_peaks"].to_numpy(),
            "read_rule": args.read_rule,
            "min_dist": args.min_dist,
            "max_dist": args.max_dist,
            "pseudocount": pc,
        })

        # Kevin's fisher_res is NA when either region has no FIRE element among the
        # shared reads, and those rows are silently dropped downstream. That discards
        # exactly the constitutively-closed-partner cases, so they are kept here and
        # flagged; --kevin-compat restores his row set.
        n_all = len(df)
        df = df[df["n_shared_reads"] > 0]
        if args.kevin_compat:
            df = df[~(df["zero_access_cre1"] | df["zero_access_cre2"])]
        df = df.reset_index(drop=True)

        stat_path = out_dir / f"{s}_coaccess_stat.tsv.gz"
        df.to_csv(stat_path, sep="\t", index=False)

        # One row per distinct cCRE pair. The 2x2 depends only on the two cCREs, so
        # this is a dedup plus the gene list. BH is computed HERE and only here: the
        # gene-anchored table repeats a pair once per containing gene, so correcting
        # on it would over-count the tests.
        g = (df.groupby("CRE_pair", sort=False)
               .agg(n_genes=("gene_id", "nunique"),
                    gene_ids=("gene_id", lambda v: ";".join(sorted(set(v)))),
                    gene_names=("gene_name", lambda v: ";".join(sorted(set(v)))))
               .reset_index())
        uniq_rows = df.drop_duplicates("CRE_pair").drop(
            columns=["gene_id", "gene_name", "transcript_type"])
        pairs_df = uniq_rows.merge(g, on="CRE_pair", how="left")
        for col, src in (("fdr", "pval"), ("fdr_raw", "pval_raw")):
            p = pairs_df[src].to_numpy()
            order = np.argsort(p)
            q = np.minimum.accumulate(
                (p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
            out = np.empty_like(q)
            out[order] = np.minimum(q, 1)
            pairs_df[col] = out
        pairs_path = out_dir / f"{s}_coaccess_pairs.tsv.gz"
        pairs_df.to_csv(pairs_path, sep="\t", index=False)

        print(f"{s}: {len(df)}/{n_all} gene-anchored rows with shared reads, "
              f"{len(pairs_df)} distinct pairs, "
              f"{int((pairs_df['pval'] < 0.05).sum())} with pval < 0.05, "
              f"{int((pairs_df['fdr'] < 0.05).sum())} at FDR < 0.05", flush=True)

    try:
        tmpdir.rmdir()
    except OSError:
        pass


if __name__ == "__main__":
    main()
