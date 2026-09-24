"""Final PDF plots only. Smoke mode renders every figure into memory."""

import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from expression_summary import EXPRESSION_BINS
from nrl_annotation import centered_mean33, annotation_acf


def make_figures(records, regions, binary, raw_acf, summaries, cluster_summary, pdf_dir=None):
    if pdf_dir is not None:
        pdf_dir.mkdir(parents=True, exist_ok=False)
    written = []

    def save(fig, name):
        fig.tight_layout()
        if pdf_dir is None:
            with io.BytesIO() as buffer:
                fig.savefig(buffer, format="pdf", bbox_inches="tight")
                data = buffer.getvalue()
                assert data.startswith(b"%PDF-") and b"%%EOF" in data[-1024:]
        else:
            fig.savefig(pdf_dir / (name + ".pdf"), bbox_inches="tight")
        plt.close(fig)
        written.append(name)

    def empty(ax, label):
        ax.text(0.5, 0.5, label, ha="center", va="center", transform=ax.transAxes)

    genes = summaries["gene"]
    mapped = genes[genes.expr_bin.isin(EXPRESSION_BINS)]
    colors = dict(zip(EXPRESSION_BINS, ["#888888", "#4477aa", "#66ccee", "#228833", "#cc6677"]))
    fig, ax = plt.subplots(figsize=(9, 4))
    for group in EXPRESSION_BINS:
        curves = [summaries["gene_acf"][gene] for gene in mapped.index[mapped.expr_bin.eq(group)]
                  if gene in summaries["gene_acf"]]
        if curves:
            mean = np.mean(curves, axis=0)
            ax.plot(np.arange(1, len(mean)), mean[1:], label=f"{group}: {len(curves)} genes",
                    color=colors[group], lw=0.9)
    if ax.lines:
        ax.legend(fontsize=8)
    else:
        empty(ax, "No genes with mapped expression and valid ACF")
    ax.axhline(0, color="black", lw=0.4)
    ax.set(xlabel="Lag (bp)", ylabel="Raw binary ACF", title="Equal gene weight; equal TSS and sample means")
    save(fig, "aggregate_raw_acf_by_expression")

    for metric, ylabel in [("nrl_bp", "Reliable NRL (bp)"), ("regularity", "Mean absolute accepted-extremum ACF")]:
        fig, ax = plt.subplots(figsize=(9, 4))
        # Empirical CDFs remain interpretable for a one-gene smoke-test bin.
        for group in EXPRESSION_BINS:
            values = np.sort(mapped.loc[mapped.expr_bin.eq(group), metric].dropna())
            if len(values):
                ax.step(values, np.arange(1, len(values) + 1) / len(values), where="post",
                        marker=".", label=f"{group}: {len(values)} genes", color=colors[group])
        if ax.lines:
            ax.legend(fontsize=8)
        else:
            empty(ax, "No eligible gene summaries")
        ax.set(xlabel=ylabel, ylabel="Cumulative fraction of genes", ylim=(0, 1.05),
               title="Gene summaries by expression bin")
        save(fig, metric + "_distribution_by_expression")

        fig, ax = plt.subplots(figsize=(7, 5))
        for group in EXPRESSION_BINS:
            sub = mapped[mapped.expr_bin.eq(group)].dropna(subset=[metric, "mean_tpm"])
            if len(sub):
                ax.scatter(sub.mean_tpm, sub[metric], s=18, alpha=0.6,
                           label=f"{group}: {len(sub)} genes", color=colors[group], rasterized=True)
        ax.set_xscale("symlog", linthresh=1)
        ax.set(xlabel="Mean TPM across 0/5/10/15 min RNA samples (symlog)", ylabel=ylabel,
               title="One observation per gene")
        if ax.collections:
            ax.legend(fontsize=8)
        else:
            empty(ax, "No eligible gene summaries")
        save(fig, metric + "_versus_expression")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    bins = [b for b in EXPRESSION_BINS if mapped.expr_bin.eq(b).any()]
    for j, group in enumerate(bins):
        sub = mapped[mapped.expr_bin.eq(group)]
        fractions = sub.reliable_fraction.dropna()
        if len(fractions):
            axes[0].bar(j, fractions.mean(), color=colors[group])
            axes[0].text(j, fractions.mean() + 0.025, f"{len(fractions)} genes", ha="center", fontsize=7)
        tss = summaries["tss"]
        tss = tss[tss.expr_bin.eq(group) & tss.window_eligible]
        for ax, subset in [(axes[1], tss), (axes[2], tss[tss.covered])]:
            if len(subset):
                count = int(subset.has_reliable_nrl.sum())
                ax.bar(j, count / len(subset), color=colors[group])
                ax.text(j, count / len(subset) + .025, f"{count}/{len(subset)}", ha="center", fontsize=7)
    for ax, title in zip(axes, ["Reliable read fraction: mean across genes",
                                "TSS with >=1 reliable read / all eligible TSS",
                                "TSS with >=1 reliable read / covered TSS"]):
        ax.set(ylim=(0, 1.18), ylabel="Fraction", title=title)
        ax.set_xticks(range(len(bins)), bins, rotation=30, ha="right")
    save(fig, "reliable_nrl_fractions")

    clustered = records[records.status.eq("clustered")]
    clusters = sorted(clustered.cluster.unique(), key=int)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for cluster in clusters:
        sub = clustered[clustered.cluster.eq(cluster)]
        axes[0].scatter(sub.umap1, sub.umap2, label=f"C{cluster}", s=8, alpha=.6, rasterized=True)
    if clusters:
        axes[0].legend(fontsize=7, ncol=2)
    for ax, field in [(axes[1], "nrl_bp"), (axes[2], "regularity")]:
        sub = clustered[clustered.reliable_nrl] if field == "nrl_bp" else clustered
        sub = sub.dropna(subset=[field, "umap1", "umap2"])
        if len(sub):
            points = ax.scatter(sub.umap1, sub.umap2, c=sub[field], s=9, cmap="viridis", rasterized=True)
            fig.colorbar(points, ax=ax, label=field)
        else:
            empty(ax, "No reliable NRL" if field == "nrl_bp" else "No accepted extrema")
    for ax, title in zip(axes, ["Unchanged raw-ACF Leiden clusters", "Reliable NRL annotation", "Regularity annotation"]):
        ax.set(title=title, xlabel="UMAP 1", ylabel="UMAP 2")
    save(fig, "leiden_annotations")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, field in [(axes[0], "nrl_bp"), (axes[1], "regularity")]:
        for j, cluster in enumerate(clusters):
            sub = clustered[clustered.cluster.eq(cluster)]
            if field == "nrl_bp":
                sub = sub[sub.reliable_nrl]
            values = sub[field].dropna().to_numpy()
            if len(values):
                q25, median, q75 = np.quantile(values, [.25, .5, .75])
                ax.errorbar(j, median, yerr=[[median - q25], [q75 - median]], fmt="o", capsize=3)
                ax.annotate(f"n={len(values)}", (j, q75), fontsize=7)
        ax.set_xticks(range(len(clusters)), [f"C{c}" for c in clusters])
        ax.set(xlabel="Existing Leiden cluster", ylabel=field, title="Molecule median and IQR (descriptive)")
    save(fig, "cluster_nrl_regularity")

    fig, ax = plt.subplots(figsize=(9, 4))
    averages = cluster_summary[0]
    for cluster in clusters:
        sub = averages[averages.cluster.eq(cluster) & averages.lag_bp.gt(0)]
        ax.plot(sub.lag_bp, sub.mean_acf, label=f"C{cluster}", lw=.8)
    if clusters:
        ax.legend(fontsize=7)
    ax.set(xlabel="Lag (bp)", ylabel="Raw ACF", title="Existing cluster summary function; molecule means")
    save(fig, "cluster_raw_acf")

    composition = summaries["gene_composition"].join(mapped[["expr_bin"]], how="inner")
    by_bin = composition.groupby("expr_bin").mean().reindex(bins)
    fig, ax = plt.subplots(figsize=(9, 4))
    bottom = np.zeros(len(bins))
    for cluster in by_bin.columns:
        values = by_bin[cluster].to_numpy()
        ax.bar(np.arange(len(bins)), values, bottom=bottom, label=cluster)
        bottom += values
    ax.set_xticks(range(len(bins)), bins)
    ax.set(ylim=(0, 1), ylabel="Mean within-gene read fraction", xlabel="Expression bin",
           title="Cluster composition; each covered gene has equal weight")
    if len(by_bin.columns):
        ax.legend(title="Cluster", fontsize=7, ncol=3)
    save(fig, "cluster_composition_by_expression")

    # At most three representative TSSs, independent of production cohort size.
    # The first covered TSSs in BED order are shown; no enrichment for regularity.
    shown_regions = regions[regions.region_id.isin(records.region_id.unique())].head(3)
    for region in shown_regions.itertuples():
        ix = np.flatnonzero(records.region_id.eq(region.region_id))
        # Limit display only; every retained row was clustered and annotated.
        ix = sorted(ix, key=lambda i: (records.iloc[i].cluster or "~", records.iloc[i].row_id))[:80]
        oriented = binary[ix] if region.strand == "+" else binary[ix, ::-1]
        # Exact even-width interval: plus -1000..999; minus -999..1000.
        rel = np.arange(-1000, 1000) if region.strand == "+" else np.arange(-999, 1001)
        fig, axes = plt.subplots(3, 2, figsize=(13, 10))
        axes[0, 0].imshow(oriented, aspect="auto", interpolation="nearest", cmap="Greys",
                          extent=(rel[0]-.5, rel[-1]+.5, len(ix)-.5, -.5), rasterized=True)
        axes[0, 0].axvline(0, color="red", lw=.6)
        axes[0, 0].set(xlabel="Strand-oriented distance from TSS (bp)", ylabel="Molecules (cluster order)",
                       title="Original binary m6A; display subset")
        limit = max(float(np.nanquantile(np.abs(raw_acf[ix, 1:]), .99)), .01)
        axes[0, 1].imshow(raw_acf[ix, 1:], aspect="auto", interpolation="nearest", cmap="RdBu_r",
                          vmin=-limit, vmax=limit, extent=(.5, 1999.5, len(ix)-.5, -.5), rasterized=True)
        axes[0, 1].set(xlabel="Lag (bp)", title="Raw ACF; identical molecule order")
        axes[1, 0].plot(rel, oriented.mean(axis=0), lw=.6, label="Raw mean")
        axes[1, 0].plot(rel[16:-16], centered_mean33(oriented.mean(axis=0)), label="33-bp mean (display)")
        axes[1, 0].axvline(0, color="black", lw=.5)
        axes[1, 0].legend(fontsize=7)
        axes[1, 0].set(xlabel="Strand-oriented distance from TSS (bp)", ylabel="m6A call fraction")
        # First displayed molecule with most accepted peaks, for a diagnostic
        # view of extrema and regression; this is explicitly a selected example.
        representative = max(ix, key=lambda i: (records.iloc[i].reliable_nrl, records.iloc[i].n_positive_peaks))
        row = records.iloc[representative]
        acf = annotation_acf(centered_mean33(binary[representative]))
        axes[1, 1].plot(np.arange(1, len(acf)), acf[1:], lw=.8)
        for field, marker, color in [("positive_peak_lags", "^", "red"), ("negative_peak_lags", "v", "blue")]:
            positions = np.asarray(row[field], dtype=int)
            axes[1, 1].scatter(positions, acf[positions], marker=marker, c=color, s=18)
        axes[1, 1].axvspan(1468, 1967, color="gray", alpha=.15, label="<500 paired positions")
        axes[1, 1].set(xlabel="Lag (bp)", ylabel="33-bp annotation ACF", title="Selected example: most supported peaks")
        positions = np.asarray(row.positive_peak_lags)
        numbers = np.arange(1, len(positions) + 1)
        axes[2, 1].scatter(numbers, positions)
        if np.isfinite(row.nrl_bp):
            axes[2, 1].plot(numbers, row.regression_intercept_bp + row.nrl_bp * numbers, color="red")
        axes[2, 1].set(xlabel="Positive peak number", ylabel="Peak lag (bp)",
                       title=f"NRL={row.nrl_bp:.1f}; R²={row.regression_r_squared:.3f}; {row.annotation_status}")
        example = binary[representative] if region.strand == "+" else binary[representative, ::-1]
        axes[2, 0].plot(rel[16:-16], centered_mean33(example))
        axes[2, 0].axvline(0, color="black", lw=.5)
        axes[2, 0].set(xlabel="Strand-oriented distance from TSS (bp)", ylabel="33-bp mean m6A",
                       title=f"Same example molecule; {row.sample_name}; C{row.cluster}")
        fig.suptitle(f"{region.gene_name} | {region.chrom}:{region.window_start0}-{region.window_end0} "
                     f"({region.strand}, BED) | {region.expr_bin}", fontsize=11)
        save(fig, "representative_" + region.gene_id + "_" + str(region.tss0))
    return len(written)
