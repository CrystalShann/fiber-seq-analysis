"""Cluster averages, descriptive repeat peaks, pooled focal-allele summaries and figures."""

import numpy as np
import pandas as pd
from scipy.signal import find_peaks


def repeat_peak(profile, lower=140, upper=250):
    """Strongest positive *local* peak in a declared band; no forced period."""
    if len(profile) <= upper + 1 or not np.isfinite(profile).all():
        return np.nan, np.nan
    peaks, _ = find_peaks(profile)
    peaks = peaks[(peaks >= lower) & (peaks <= upper) & (profile[peaks] > 0)]
    if not len(peaks):
        return np.nan, np.nan
    peak = peaks[np.argmax(profile[peaks])]
    return int(peak), float(profile[peak])


def summarize(profiles, records):
    clustered = records[records.status == "clustered"]
    averages, summaries, composition = [], [], []
    cluster_ids = sorted(clustered.cluster.unique(), key=int)
    samples = sorted(clustered.sample_name.unique())
    for cluster in cluster_ids:
        in_cluster = records.cluster.eq(cluster).to_numpy() & records.status.eq("clustered").to_numpy()
        subset = profiles[in_cluster]
        mean = subset.mean(axis=0)
        sd = subset.std(axis=0, ddof=1) if len(subset) > 1 else np.full(profiles.shape[1], np.nan)
        lag, height = repeat_peak(mean)
        summaries.append(dict(region_id=records.region_id.iloc[0], cluster=cluster,
                              n_reads=len(subset), fraction_clustered=len(subset) / len(clustered),
                              mean_m6a_call_fraction=records.loc[in_cluster, "m6a_call_fraction"].mean(),
                              positive_local_peak_140_250_bp=lag, peak_acf=height))
        averages.append(pd.DataFrame(dict(region_id=records.region_id.iloc[0], cluster=cluster,
                                          lag_bp=np.arange(profiles.shape[1]), mean_acf=mean,
                                          sd_acf=sd, sem_acf=sd / np.sqrt(len(subset)), n_reads=len(subset))))
        for sample in samples:
            c = clustered.cluster.eq(cluster)
            s = clustered.sample_name.eq(sample)
            a = int((c & s).sum())
            composition.append(dict(region_id=records.region_id.iloc[0], cluster=cluster, sample_name=sample,
                                    n_reads=a, fraction_within_cluster=a / len(subset),
                                    fraction_within_sample=a / int(s.sum())))
    avg = pd.concat(averages, ignore_index=True) if averages else pd.DataFrame(columns=[
        "region_id", "cluster", "lag_bp", "mean_acf", "sd_acf", "sem_acf", "n_reads"])
    stats = pd.DataFrame(summaries, columns=["region_id", "cluster", "n_reads", "fraction_clustered",
                                           "mean_m6a_call_fraction", "positive_local_peak_140_250_bp", "peak_acf"])
    counts = pd.DataFrame(composition, columns=["region_id", "cluster", "sample_name", "n_reads",
                                               "fraction_within_cluster", "fraction_within_sample"])
    return avg, stats, counts


def plot_region(directory, region, binary, profiles, records, averages, composition):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    directory.mkdir(parents=True, exist_ok=True)
    title = f"{region.region_id}\n{region.chr}:{region.analysis_start:,}-{region.analysis_end:,} (1-based)"
    clustered = records[records.status == "clustered"]
    cluster_ids = sorted(clustered.cluster.unique(), key=int)
    palette = ["#1c86ee", "#E31A1C", "#008b00", "#6A3D9A", "#FF7F00", "black", "#ffd700"]
    colors = {c: palette[i % len(palette)] for i, c in enumerate(cluster_ids)}

    def save(fig, name):
        fig.suptitle(title, fontsize=10)
        fig.tight_layout()
        fig.savefig(directory / (name + ".pdf"), bbox_inches="tight")
        plt.close(fig)

    if clustered.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No clusterable profiles\n" + records.status.value_counts().to_string(),
                ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        save(fig, "qc")
        return

    fig, axes = plt.subplots(2, 1, figsize=(10, 7))
    for c in cluster_ids:
        a = averages[(averages.cluster == c) & (averages.lag_bp > 0)]
        axes[0].plot(a.lag_bp, a.mean_acf, color=colors[c], label=f"C{c} (n={int(a.n_reads.iloc[0])})", lw=0.9)
    axes[0].axhline(0, color="gray", lw=0.5)
    axes[0].set(xlabel="Lag (bp)", ylabel="Mean ACF", title="Cluster averages; no smoothing")
    axes[0].legend(fontsize=8, ncol=4)
    avg_matrix = averages.pivot(index="cluster", columns="lag_bp", values="mean_acf").loc[cluster_ids]
    avg_matrix = avg_matrix.loc[:, avg_matrix.columns > 0]
    bound = max(float(np.abs(avg_matrix.to_numpy()).max()), 0.01)
    im = axes[1].imshow(avg_matrix, aspect="auto", interpolation="nearest", cmap="RdBu_r",
                        vmin=-bound, vmax=bound, extent=(0.5, profiles.shape[1]-0.5, len(cluster_ids)-0.5, -0.5))
    axes[1].set_yticks(range(len(cluster_ids)), [f"C{c}" for c in cluster_ids])
    axes[1].set(xlabel="Lag (bp)", ylabel="Cluster")
    fig.colorbar(im, ax=axes[1], label="Mean ACF", shrink=0.6)
    save(fig, "cluster_average_acf")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for c in cluster_ids:
        m = clustered[clustered.cluster == c]
        axes[0].scatter(m.umap1, m.umap2, s=7, color=colors[c], label=f"C{c}", rasterized=True)
    axes[0].legend(fontsize=8, markerscale=2)
    axes[0].set(title="ACF clusters", xlabel="UMAP 1", ylabel="UMAP 2")
    im = axes[1].scatter(clustered.umap1, clustered.umap2, s=7,
                         c=clustered.m6a_call_fraction, cmap="viridis", rasterized=True)
    axes[1].set(title=f"m6A calls / {region.width:,} bases", xlabel="UMAP 1", ylabel="UMAP 2")
    fig.colorbar(im, ax=axes[1], label="m6A call fraction")
    save(fig, "umap")

    fractions = composition.pivot(index="sample_name", columns="cluster", values="fraction_within_sample")[cluster_ids]
    fig, ax = plt.subplots(figsize=(12, 5))
    fractions.plot.bar(stacked=True, ax=ax, color=[colors[c] for c in cluster_ids], width=0.85)
    ax.set(xlabel="Sample", ylabel="Fraction of clustered molecules", ylim=(0, 1))
    ax.tick_params(axis="x", labelsize=7)
    ax.legend(title="Cluster", fontsize=7, ncol=3)
    save(fig, "sample_composition")


def annotate_alleles(records, region, sample_table):
    """Map sample-local HP tags to the focal REF/ALT and biological cell line."""
    result = records.copy().reset_index(drop=True)
    if (region.ref not in ('A','C','G','T') or region.alt not in ('A','C','G','T') or region.ref == region.alt
            or not result.haplotype.isin(['HP1', 'HP2']).all()
            or not result.focal_genotype.isin(['0|1', '1|0']).all()
            or not result.allele_status.eq('phased_focal_genotype').all()):
        raise ValueError('Allele comparisons require resolved phased heterozygous SNP reads')
    if sample_table.sample_name.duplicated().any():
        raise ValueError('Duplicate samples in sample metadata')
    result['cell_line'] = result.sample_name.map(sample_table.set_index('sample_name').cell_line)
    if result.cell_line.isna().any() or result.cell_line.eq('').any():
        raise ValueError('Every input sample needs an explicit biological cell_line in sample metadata')
    alt = ((result.haplotype.eq('HP2') & result.focal_genotype.eq('0|1'))
           | (result.haplotype.eq('HP1') & result.focal_genotype.eq('1|0')))
    result['allele_group'] = np.where(alt, 'ALT', 'REF')
    result['allele_base'] = np.where(alt, region.alt, region.ref)
    expected = region.focal_snp + ': ' + result.allele_base
    if not result.allele_label.eq(expected).all():
        raise ValueError('HP/genotype mapping disagrees with the saved focal allele label')
    return result


def summarize_alleles(profiles, records, region):
    """Pool unique molecules by focal allele, with no sample-level stratification."""
    profiles = np.asarray(profiles)
    records = records.copy().reset_index(drop=True)
    if (profiles.ndim != 2 or profiles.shape[1] < 1 or records.empty
            or len(records) != len(profiles) or records.row_id.duplicated().any()):
        raise ValueError('Profiles and unique read records must be aligned')
    if not records.allele_group.isin(['REF', 'ALT']).all():
        raise ValueError('Every read must have a resolved REF or ALT allele')
    valid = records.acf_valid.to_numpy(dtype=bool)
    clustered = records.status.eq('clustered').to_numpy()
    if (clustered & ~valid).any():
        raise ValueError('Clustered molecules must have valid ACFs')
    if not np.isfinite(profiles[valid]).all():
        raise ValueError('Valid ACF profiles must be finite')
    # Sample metadata is used only to validate duplicate physical molecules.
    records['allele_analysis_included'] = True
    for _, group in records.groupby('RID', sort=False):
        if len(group) < 2:
            continue
        for col in ('cell_line', 'allele_group', 'status', 'cluster', 'acf_valid'):
            if group[col].nunique(dropna=False) != 1:
                raise ValueError(f'Conflicting duplicate molecule {group.RID.iloc[0]}: {col}')
        ix = group.index.to_numpy()
        if not all(np.allclose(profiles[ix[0]], profiles[i], equal_nan=True) for i in ix[1:]):
            raise ValueError(f'Conflicting ACFs for duplicate molecule {group.RID.iloc[0]}')
        records.loc[ix[1:], 'allele_analysis_included'] = False
    included = records.allele_analysis_included.to_numpy()
    band_complete = profiles.shape[1] >= 252
    peaks = [repeat_peak(p) if v else (np.nan, np.nan) for p, v in zip(profiles, valid)]
    records['positive_local_peak_140_250_bp'] = [p[0] for p in peaks]
    records['peak_acf'] = [p[1] for p in peaks]
    records['peak_band_complete'] = band_complete
    clusters = sorted(records.loc[clustered, 'cluster'].astype(str).unique(), key=int)
    coverage, proportions, curves, peak_rows = [], [], [], []
    for allele in ('REF', 'ALT'):
        original = records.allele_group.eq(allele).to_numpy()
        mask = included & original
        vm, cm = mask & valid, mask & clustered
        n_unique, n_valid, n_clustered = int(mask.sum()), int(vm.sum()), int(cm.sum())
        coverage.append(dict(allele_group=allele, n_input=int(original.sum()),
            n_unique=n_unique, n_duplicates_removed=int(original.sum())-n_unique,
            n_valid=n_valid, n_invalid=n_unique-n_valid, n_clustered=n_clustered,
            n_valid_not_clustered=n_valid-n_clustered, peak_band_complete=band_complete))
        mean = profiles[vm].mean(axis=0) if n_valid else np.full(profiles.shape[1], np.nan)
        curves.append(pd.DataFrame(dict(allele_group=allele,
            lag_bp=np.arange(profiles.shape[1]), mean_acf=mean, n_reads=n_valid)))
        for cluster in clusters:
            count = int((cm & records.cluster.astype(str).eq(cluster).to_numpy()).sum())
            proportions.append(dict(allele_group=allele, cluster=cluster,
                n_in_cluster=count, n_clustered=n_clustered,
                fraction=count/n_clustered if n_clustered else np.nan))
        n_peak = int(records.loc[vm, 'positive_local_peak_140_250_bp'].notna().sum())
        peak_rows.append(dict(allele_group=allele, n_valid=n_valid,
            n_with_peak=n_peak if band_complete else np.nan,
            n_without_peak=n_valid-n_peak if band_complete else np.nan,
            peak_fraction=n_peak/n_valid if band_complete and n_valid else np.nan,
            median_peak_lag=records.loc[vm, 'positive_local_peak_140_250_bp'].median() if n_peak else np.nan,
            median_peak_height=records.loc[vm, 'peak_acf'].median() if n_peak else np.nan,
            peak_band_complete=band_complete))
    tables = dict(
        allele_coverage=pd.DataFrame(coverage),
        allele_cluster_proportions=pd.DataFrame(proportions, columns=[
            'allele_group', 'cluster', 'n_in_cluster', 'n_clustered', 'fraction']),
        allele_acf_summary=pd.concat(curves, ignore_index=True),
        allele_peak_features=pd.DataFrame(peak_rows),
        allele_read_audit=records[['row_id', 'RID', 'sample_name', 'cell_line',
            'allele_group', 'allele_base', 'acf_valid', 'status', 'cluster',
            'allele_analysis_included', 'peak_band_complete',
            'positive_local_peak_140_250_bp', 'peak_acf']].copy())
    for table in tables.values():
        table.insert(0, 'region_id', region.region_id)
    return tables


def plot_alleles(directory, region, tables):
    """Four REF/ALT-only figures; all valid unique molecules receive equal weight."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    directory.mkdir(parents=True, exist_ok=True)
    alleles = ('REF', 'ALT')
    colors = {'REF': '#2878B5', 'ALT': '#D55E00'}
    labels = [f'REF ({region.ref})', f'ALT ({region.alt})']
    title = f'{region.region_id}\n{region.focal_snp}: {region.ref} / {region.alt}'
    coverage = tables['allele_coverage'].set_index('allele_group').reindex(alleles)

    def finish(fig, name):
        fig.suptitle(title, fontsize=11)
        fig.tight_layout(rect=(0, .05, 1, .92))
        fig.savefig(directory / f'{name}.pdf', bbox_inches='tight')
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, column, heading in zip(axes, ['n_valid', 'n_clustered'],
                                   ['Valid ACF reads', 'Clustered reads']):
        bars = ax.bar(labels, coverage[column], color=[colors[a] for a in alleles])
        ax.bar_label(bars, padding=3)
        ax.set(ylabel='Unique molecules', title=heading)
        ax.set_ylim(0, max(1, coverage[column].max()) * 1.18)
    finish(fig, 'allele_coverage')

    proportions = tables['allele_cluster_proportions']
    clusters = sorted(proportions.cluster.unique(), key=int)
    fig, ax = plt.subplots(figsize=(max(7, len(clusters)*1.2), 4.5))
    for i, allele in enumerate(alleles):
        d = proportions[proportions.allele_group.eq(allele)].set_index('cluster').reindex(clusters)
        ax.bar(np.arange(len(clusters))+(i-.5)*.36, d.fraction, width=.34,
               color=colors[allele], label=labels[i])
    ax.set(xticks=np.arange(len(clusters)), xticklabels=[f'C{c}' for c in clusters],
           xlabel='ACF cluster', ylabel='Fraction of clustered molecules within allele', ylim=(0, 1.05))
    ax.legend()
    if not clusters:
        ax.text(.5, .5, 'No clustered molecules', ha='center', transform=ax.transAxes)
    missing = [a for a in alleles if not coverage.loc[a, 'n_clustered']]
    if missing and clusters:
        ax.text(.5, .96, 'No clustered reads: '+', '.join(missing), ha='center', va='top', transform=ax.transAxes)
    finish(fig, 'allele_cluster_proportions')

    # Keep the full-lag REF/ALT view and add two zooms of the same mean curves.
    curves = tables['allele_acf_summary']
    fig = plt.figure(figsize=(12, 8))
    grid = fig.add_gridspec(2, 2)
    axes = [fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
    maximum_lag = int(curves.lag_bp.max())
    missing = [a for a in alleles if not coverage.loc[a, 'n_valid']]
    views = [('All retained lags', 0, maximum_lag), ('Lags 0–500 bp', 0, 500),
             ('Lags 500–2000 bp', 500, 2000)]
    for ax, (heading, lower, upper) in zip(axes, views):
        for i, allele in enumerate(alleles):
            d = curves[curves.allele_group.eq(allele) & curves.lag_bp.between(lower, upper)].sort_values('lag_bp')
            ax.plot(d.lag_bp, d.mean_acf, color=colors[allele], label=labels[i], lw=1)
        ax.set(title=heading, xlabel='Lag (bp)', ylabel='Mean ACF', xlim=(lower, max(lower + 1, upper)))
        ax.legend()
        if maximum_lag < lower:
            ax.text(.5, .5, 'Lag range not available', ha='center', transform=ax.transAxes)
        elif missing:
            ax.text(.5, .5, 'No valid ACF reads: '+', '.join(missing), ha='center', transform=ax.transAxes)
    fig.text(.5, .01, f'Available lags: 0–{maximum_lag} bp. Each panel uses the same pooled REF/ALT means; '
             'vertical scales are independent.', ha='center', fontsize=9)
    finish(fig, 'allele_average_acf')

    features = tables['allele_peak_features'].set_index('allele_group').reindex(alleles)
    audit = tables['allele_read_audit']
    audit = audit[audit.allele_analysis_included & audit.acf_valid]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))
    bars = axes[0].bar(labels, features.peak_fraction, color=[colors[a] for a in alleles])
    axes[0].set(ylabel='Fraction with a qualifying peak', ylim=(0, 1.15))
    for bar, value in zip(bars, features.peak_fraction):
        if np.isfinite(value):
            axes[0].text(bar.get_x()+bar.get_width()/2, value+.025, f'{value:.2f}', ha='center')
    for ax, column, ylabel in zip(axes[1:],
            ['positive_local_peak_140_250_bp', 'peak_acf'],
            ['Molecule peak lag (bp)', 'Molecule peak ACF height']):
        for i, allele in enumerate(alleles):
            values = audit.loc[audit.allele_group.eq(allele), column].dropna().to_numpy()
            if len(values):
                box = ax.boxplot([values], positions=[i], widths=.45, patch_artist=True,
                    medianprops={'color': 'black'}, manage_ticks=False)
                box['boxes'][0].set_facecolor(colors[allele])
                box['boxes'][0].set_alpha(.65)
            else:
                ax.text(i, .5, 'No peaks', ha='center', transform=ax.get_xaxis_transform())
        ax.set(xticks=[0, 1], xticklabels=labels, xlim=(-.6, 1.6), ylabel=ylabel)
    axes[1].set_ylim(135, 255)
    if not features.peak_band_complete.all():
        for ax in axes:
            ax.clear()
            ax.text(.5, .5, '140–250 bp peak range not fully available', ha='center',
                    transform=ax.transAxes, fontsize=9)
            ax.set_axis_off()
    else:
        for i, allele in enumerate(alleles):
            if not features.loc[allele, 'n_valid']:
                axes[0].text(i, .5, 'No valid reads', ha='center', transform=axes[0].get_xaxis_transform())
    fig.text(.5, .01, 'Peak fraction: all valid reads. Distributions: peak-bearing reads only; boxes show median and quartiles.',
             ha='center', fontsize=9)
    finish(fig, 'allele_peak_features')
