"""Gene/TSS summaries in memory; molecules are never expression replicates."""

import numpy as np
import pandas as pd

EXPRESSION_BINS = ["not_expressed", "Q1_low", "Q2", "Q3", "Q4_high"]


def summarize_expression(records, regions, raw_acf):
    r = records.copy()
    r["reliable_nrl_bp"] = r.nrl_bp.where(r.reliable_nrl)
    r["composition_cluster"] = r.cluster.where(r.status.eq("clustered"), "unclustered")
    keys = ["region_id", "sample_name"]
    per_sample = r.groupby(keys, sort=False).agg(
        n_reads=("row_id", "size"), n_reliable=("reliable_nrl", "sum"),
        nrl_bp=("reliable_nrl_bp", "median"), regularity=("regularity", "median"),
        reliable_fraction=("reliable_nrl", "mean"))
    metrics = ["nrl_bp", "regularity", "reliable_fraction"]
    per_tss = per_sample.groupby("region_id")[metrics].mean()
    totals = per_sample.groupby("region_id")[["n_reads", "n_reliable"]].sum()
    per_tss = regions.merge(per_tss.join(totals), on="region_id", how="left", validate="one_to_one")
    per_tss[["n_reads", "n_reliable"]] = per_tss[["n_reads", "n_reliable"]].fillna(0).astype(int)
    per_tss["covered"] = per_tss.n_reads.gt(0)
    per_tss["has_reliable_nrl"] = per_tss.n_reliable.gt(0)
    # Equal sample weight within a TSS, then equal TSS weight within a gene.
    # Missing sample coverage stays missing, never a fabricated zero signal.
    per_gene = per_tss.groupby("gene_id", sort=False).agg(
        gene_name=("gene_name", "first"), mean_tpm=("mean_tpm", "first"),
        expr_bin=("expr_bin", "first"), n_tss=("region_id", "size"),
        n_covered_tss=("covered", "sum"), n_reads=("n_reads", "sum"),
        n_reliable=("n_reliable", "sum"), nrl_bp=("nrl_bp", "mean"),
        regularity=("regularity", "mean"), reliable_fraction=("reliable_fraction", "mean"),
        fraction_tss_with_reliable=("has_reliable_nrl", "mean"))
    per_gene["has_reliable_nrl"] = per_gene.n_reliable.gt(0)

    # Average curves in the same hierarchy. Never weight expression curves by
    # molecule count. Each gene contributes at most one curve to its bin.
    tss_curves = {}
    for (rid, _), ix in r.groupby(keys, sort=False).indices.items():
        ix = ix[r.iloc[ix].acf_valid.to_numpy()]
        if len(ix):
            tss_curves.setdefault(rid, []).append(raw_acf[ix].mean(axis=0))
    gene_curves = {}
    gene_by_region = regions.set_index("region_id").gene_id
    for rid, curves in tss_curves.items():
        gene_curves.setdefault(gene_by_region[rid], []).append(np.mean(curves, axis=0))
    gene_curves = {gene: np.mean(curves, axis=0) for gene, curves in gene_curves.items()}

    # Include explicit zero membership for every cluster at every observed
    # TSS/sample, plus unclustered reads; fractions sum to one per gene.
    fractions = r.groupby(keys + ["composition_cluster"]).size().unstack(fill_value=0)
    fractions = fractions.div(fractions.sum(axis=1), axis=0)
    fractions = fractions.groupby("region_id").mean()
    fractions["gene_id"] = gene_by_region.reindex(fractions.index)
    gene_composition = fractions.groupby("gene_id").mean()
    if not np.allclose(gene_composition.sum(axis=1), 1):
        raise AssertionError("Cluster composition must sum to one within each gene")
    return dict(tss_sample=per_sample, tss=per_tss, gene=per_gene,
                gene_acf=gene_curves, gene_composition=gene_composition)
