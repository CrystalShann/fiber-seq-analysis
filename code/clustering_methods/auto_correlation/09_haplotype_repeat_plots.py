"""ACF-based repeat-length 

strongest positive local maximum in 120–250 bp, prominence
>= 0.01 ACF units. Candidate maxima are >= 20 bp apart. Higher-order markers
require an observed qualifying peak within +/-20 bp of an integer multiple.

"""

import textwrap

import numpy as np
from scipy.signal import find_peaks

ALLELES = ("REF", "ALT")
COLORS = {"REF": "#2878B5", "ALT": "#D55E00"}
LOWER, UPPER = 120, 250
PROMINENCE, DISTANCE, MULTIPLE_TOLERANCE = 0.01, 20, 20


def repeat_candidates(profile):
    """Return primary and observed higher-order lags; undefined is never zero."""
    profile = np.asarray(profile)
    if len(profile) <= UPPER + 1 or not np.isfinite(profile).all():
        return np.nan, np.array([], dtype=int)
    peaks, _ = find_peaks(profile, prominence=PROMINENCE, distance=DISTANCE)
    peaks = peaks[profile[peaks] > 0]
    band = peaks[(peaks >= LOWER) & (peaks <= UPPER)]
    if not len(band):
        return np.nan, np.array([], dtype=int)
    primary = int(band[np.argmax(profile[band])])
    higher = []
    for multiple in range(2, (len(profile) - 1) // primary + 1):
        nearby = peaks[np.abs(peaks - multiple * primary) <= MULTIPLE_TOLERANCE]
        if len(nearby):
            higher.append(int(nearby[np.argmax(profile[nearby])]))
    return primary, np.array(higher, dtype=int)


def summarize_repeats(profiles, records, region, allele_tables):
    """Reuse the allele layer's validated physical-molecule deduplication mask."""
    profiles = np.asarray(profiles)
    audit = allele_tables["allele_read_audit"]
    if not np.array_equal(records.row_id.to_numpy(), audit.row_id.to_numpy()):
        raise ValueError("ACFs and allele audit must be in the same row order")
    included = audit.allele_analysis_included.to_numpy(dtype=bool)
    valid = audit.acf_valid.to_numpy(dtype=bool)
    results = []
    for allele in ALLELES:
        selected = included & audit.allele_group.eq(allele).to_numpy()
        indices = np.flatnonzero(selected & valid)
        values = profiles[indices]
        median = np.median(values, axis=0) if len(indices) else np.full(profiles.shape[1], np.nan)
        primary, higher = repeat_candidates(median)
        per_fiber = np.array([repeat_candidates(p)[0] for p in values], dtype=float)
        # Stable sort: defined peak lags ascending, then undefined peaks last.
        order = np.argsort(np.nan_to_num(per_fiber, nan=np.inf), kind="stable")
        results.append(dict(region_id=region.region_id, allele=allele,
            base=region.ref if allele == "REF" else region.alt,
            median=median, primary=primary, higher=higher, peak_lags=per_fiber,
            indices=indices[order], sorted_peaks=per_fiber[order],
            n_valid=len(indices), n_invalid=int((selected & ~valid).sum()),
            n_unique=int(selected.sum()), n_duplicates=int((~included & audit.allele_group.eq(allele).to_numpy()).sum())))
    return results


def pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def save_figure(fig, path, title, note):
    plt = pyplot()
    fig.suptitle(title, fontsize=11)
    lines = [part for line in note.splitlines() for part in textwrap.wrap(line, width=145)]
    fig.text(.5, .015, "\n".join(lines), ha="center", va="bottom", fontsize=8)
    bottom = (len(lines) * 10 + 14) / (72 * fig.get_figheight())
    fig.tight_layout(rect=(0, bottom, 1, .92))
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_locus_heatmap(directory, rows, filename="locus_haplotype_median_acf.pdf"):
    """Rows are locus x phased allele; columns are candidate lags, not estimates."""
    plt = pyplot()
    lags = np.arange(LOWER, UPPER + 1)
    matrix = np.full((len(rows), len(lags)), np.nan)
    for i, row in enumerate(rows):
        available = lags < len(row["median"])
        matrix[i, available] = row["median"][lags[available]]
    finite = np.abs(matrix[np.isfinite(matrix)])
    limit = max(float(finite.max()), .01) if len(finite) else .01
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#dddddd")
    fig, ax = plt.subplots(figsize=(11, max(3.5, .38 * len(rows) + 2)))
    im = ax.imshow(np.ma.masked_invalid(matrix), aspect="auto", interpolation="nearest",
        extent=(LOWER-.5, UPPER+.5, len(rows)-.5, -.5), cmap=cmap, vmin=-limit, vmax=limit)
    ax.set_yticks(range(len(rows)), [f'{r["region_id"]} | {r["allele"]} ({r["base"]}), n={r["n_valid"]}' for r in rows], fontsize=7)
    ax.set(xlabel="Candidate repeat length / ACF lag (bp)", ylabel="Locus | phased allele")
    for i, row in enumerate(rows):
        if np.isfinite(row["primary"]):
            ax.scatter(row["primary"], i, marker="*", s=75, color="#ffd34e", edgecolor="black", linewidth=.5)
    fig.colorbar(im, ax=ax, label="Median raw ACF")
    save_figure(fig, directory / filename, "Locus × haplotype median autocorrelation",
        "120–250 bp; common color scale. Stars: qualifying primary peaks. Gray: unavailable.\n"
        "Every unique retained fiber with a valid ACF contributes, including unclustered fibers.")


def plot_haplotype_repeats(directory, region, profiles, rows):
    plt = pyplot()
    directory.mkdir(parents=True, exist_ok=True)
    max_lag = profiles.shape[1] - 1
    lags = np.arange(max_lag + 1)
    title = f"{region.region_id}\nPhased focal alleles: REF {region.ref} / ALT {region.alt}"
    counts = "; ".join(f'{r["allele"]}: {r["n_valid"]} valid, {r["n_invalid"]} undefined ACF, {r["n_duplicates"]} duplicates removed' for r in rows)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharey="col")
    for i, row in enumerate(rows):
        for j, (lower, upper) in enumerate(((0, max_lag), (100, 500))):
            ax = axes[i, j]
            mask = (lags >= lower) & (lags <= upper)
            ax.plot(lags[mask], row["median"][mask], color=COLORS[row["allele"]], lw=1,
                label=f'{row["allele"]} median (n={row["n_valid"]})')
            ax.axvspan(LOWER, UPPER, color="#61a788", alpha=.13, label="Primary search: 120–250 bp")
            if np.isfinite(row["primary"]):
                primary = int(row["primary"])
                ax.scatter(primary, row["median"][primary], marker="*", s=130,
                    color="#ffd34e", edgecolor="black", zorder=4,
                    label=f"Candidate repeat length: {primary} bp")
                higher = row["higher"][(row["higher"] >= lower) & (row["higher"] <= upper)]
                if len(higher):
                    ax.scatter(higher, row["median"][higher], marker="o", s=22,
                        facecolors="none", edgecolors=COLORS[row["allele"]], label="Observed peaks near multiples")
            else:
                message = "No valid ACF fibers" if not row["n_valid"] else (
                    "Search band incomplete" if max_lag < UPPER + 1 else "No qualifying primary peak")
                ax.text(.98, .5, message, ha="right", transform=ax.transAxes, fontsize=8)
            ax.axhline(0, color="gray", lw=.5)
            ax.set(xlim=(lower, max(lower+1, upper)), xlabel="Lag (bp)", ylabel="Median raw ACF",
                title=f'{row["allele"]}: ' + ("all retained lags" if j == 0 else "nucleosome-scale view"))
            ax.legend(fontsize=7, loc="upper right")
    save_figure(fig, directory / "haplotype_median_acf.pdf", title,
        "Descriptive ACF-based repeat-length candidates; no smoothing. Positive maxima: prominence ≥0.01, spacing ≥20 bp.\n"
        "Higher-order markers require an observed peak within ±20 bp of an integer multiple. " + counts)

    ordered = np.concatenate([r["indices"] for r in rows])
    ordered_peaks = np.concatenate([r["sorted_peaks"] for r in rows])
    fig, axes = plt.subplots(1, 2, figsize=(12, max(5, min(10, len(ordered)/60))))
    if len(ordered):
        values = profiles[ordered, 1:]
        limit = max(float(np.quantile(np.abs(values), .99)), .01)
        for ax, upper in zip(axes, (max_lag, min(500, max_lag))):
            im = ax.imshow(values[:, :upper], aspect="auto", interpolation="nearest", cmap="RdBu_r",
                vmin=-limit, vmax=limit, extent=(.5, upper+.5, len(ordered)-.5, -.5), rasterized=True)
            visible = np.isfinite(ordered_peaks) & (ordered_peaks <= upper)
            ax.scatter(ordered_peaks[visible], np.flatnonzero(visible), s=2, color="black", rasterized=True)
            centers, labels, offset = [], [], 0
            for row in rows:
                n = row["n_valid"]
                if n:
                    centers.append(offset + (n-1)/2)
                    labels.append(f'{row["allele"]} (n={n})')
                    if offset:
                        ax.axhline(offset-.5, color="black", lw=.7)
                    offset += n
            ax.set_yticks(centers, labels)
            ax.set(xlabel="Lag (bp); lag 0 omitted", ylabel="Fibers: allele, then primary peak lag")
            fig.colorbar(im, ax=ax, label="Raw ACF (99th-percentile color limit)", shrink=.7)
    else:
        for ax in axes:
            ax.text(.5, .5, "No valid ACF fibers", ha="center", transform=ax.transAxes)
            ax.set_axis_off()
    save_figure(fig, directory / "haplotype_fiber_acf_heatmap.pdf", title,
        "All valid unique fibers, independent of Leiden assignment. No-peak fibers sort last within each allele.\n"
        "Dots: qualifying primary peaks; color limits clip the largest 1% of absolute ACF values for display. " + counts)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    labels = []
    for i, row in enumerate(rows):
        values = row["peak_lags"][np.isfinite(row["peak_lags"])]
        labels.append(f'{row["allele"]}\n{len(values)}/{row["n_valid"]} with peak')
        if len(values):
            box = axes[0].boxplot([values], positions=[i], widths=.5, patch_artist=True,
                manage_ticks=False, medianprops={"color": "black"})
            box["boxes"][0].set_facecolor(COLORS[row["allele"]])
            axes[1].hist(values, bins=np.arange(LOWER-.5, UPPER+5, 5), histtype="step", linewidth=1.5,
                weights=np.full(len(values), 1/row["n_valid"]), color=COLORS[row["allele"]], label=row["allele"])
        else:
            axes[0].text(i, (LOWER+UPPER)/2, "No estimates", ha="center", fontsize=8)
    axes[0].set(xticks=[0, 1], xticklabels=labels, xlim=(-.6, 1.6), ylim=(LOWER-5, UPPER+5),
        ylabel="Per-fiber candidate repeat length (bp)")
    axes[1].set(xlim=(LOWER-1, UPPER+1), xlabel="Per-fiber primary peak lag (bp)",
        ylabel="Fraction of all valid fibers per 5-bp bin")
    if axes[1].patches:
        axes[1].legend()
    save_figure(fig, directory / "haplotype_repeat_length_distribution.pdf", title,
        "Strongest qualifying positive local maximum in 120–250 bp; no smoothing. No peak means no estimate.\n"
        "Boxes include peak-bearing fibers only; histogram denominators include all valid fibers. " + counts)
    plot_locus_heatmap(directory, rows)


def compact_locus_rows(rows):
    """Retain only aggregate curves for the cross-locus plot, never fiber matrices."""
    return [{k: row[k] for k in ("region_id", "allele", "base", "median", "primary", "n_valid")} for row in rows]
