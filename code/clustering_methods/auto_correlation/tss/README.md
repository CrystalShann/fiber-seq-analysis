# Macrophage canonical-TSS autocorrelation

This pipeline samples molecules **before** constructing binary signals or
calculating autocorrelations. The full run selects exactly 10,000 distinct raw
read IDs: 2,500 each from `Q1_low`, `Q2`, `Q3`, and `Q4_high`. All expression bins
are clustered together. Parent scripts are imported and remain unchanged.

Three TSS regions are analysed, one per SLURM array task, each in its own output
folder. Offsets are **TSS-relative and strand-oriented** (transcriptional):
negative = upstream, positive = downstream, end exclusive.

| Array task | Folder | Window (TSS-relative) | Width |
|---|---|---|---|
| 0 | `2000_tss` | `[-1000, +1000)` | 2,000 bp |
| 1 | `upstream_1000_tss` | `[-1000, -100)` | 900 bp |
| 2 | `downstream_1000_tss` | `[+100, +1000)` | 900 bp |

`tss_autocorrelation.Rmd` describes the method and embeds every region's
published PDFs and summary tables in one report, like `../06_autocorrelation.Rmd`.
Knitting it never runs the pipeline; submit the runner first.

## Submit

Run from any working directory after creating the SLURM log directory:

```bash
mkdir -p /project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss/logs
```

Full 10,000-molecule analysis on chromosomes 1–22, X and Y for all three regions:

```bash
sbatch /project/spott/cshan/fiber-seq/code/clustering_methods/auto_correlation/tss/06_run_tss_autocorrelation.sh
```

One region only, e.g. `downstream_1000_tss` (array task 2); the combined
cross-region job is submitted only when the array covers all three tasks:

```bash
sbatch --array=2 /project/spott/cshan/fiber-seq/code/clustering_methods/auto_correlation/tss/06_run_tss_autocorrelation.sh
```

When all three tasks run, the lowest task submits `08_plot_sliding_combined.R`
as a dependent job (`afterok` on the whole array). It can also be run on its own:

```bash
Rscript 08_plot_sliding_combined.R [--root /project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss]
```

The `.sh` runner contains `#SBATCH` directives and requires a SLURM allocation.
It requests `pi-spott` / `bigmem`, one node, two CPUs, 300 GB RAM and 30 hours;
the previous full run used about 1 GB and 11 minutes per task before the sliding
windows were added. All numerical-library and Numba thread counts are capped at
one. The disk-backed pool avoids keeping all reads or their m6A signals in memory.

The full run writes to
`/project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss/<region>/`.
Runner options: `--out-dir`
(replaces the region folder; disables the combined job), `--bins-tsv`,
`--canonical-bed`, `--ft-root`, `--ref` (reference FASTA for the A/T control),
`--per-bin`, `--seed`, `--max-lag`, `--n-pcs`, `--n-neighbors`, `--resolution`,
`--chrom chr21,chr22`, the NRL/sliding options `--win-width`, `--win-step`,
`--nrl-min`, `--nrl-max`, `--min-lag`, `--min-prominence`,
`--prominence-quantile`, `--n-null`, `--decay-max-lag`, `--flat-threshold`, and
`--no-combined`. Their values are logged at the start of each task. Omit
`--chrom` for the full genome. `--max-lag` is **inclusive** and defaults to all
lags of the window (width - 1); the joint clustering pipeline requires at least
three lag features. SLURM logs are shared by all regions as
`logs/slurm_<job-id>_<task>.out`; the combined job logs to
`logs/slurm_<job-id>_combined.out`.

## Inputs and expression-table recovery

Defaults:

- Canonical BED:
  `/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed`.
- Expression table:
  `/project/spott/cshan/fiber-seq/macrophage_project/expr_access/tables/tss_expression_bins.tsv`.
- Indexed Fiber-seq BED12 files under
  `macrophage_project/FiberHMM/extract/ft_result_dir/{LPS_0,LPS_5,LPS_10,LPS_15}/extracted_results/m6a_by_chr/`.
- Reference FASTA for the A/T sequence control:
  `/project/spott/reference/human/GRCh38/hg38.fa` (as in
  `code/accessibility/expr_access/02_tss_m6a_profiles.py`).

When the default expression table is missing, `00_prepare_expression_bins.R`
evaluates the unchanged `code/accessibility/expr_access/01_expression_bins.R`,
redirecting only its `OUT_DIR` assignment to a working directory. It uses the
original STAR counts, GTF, TPM calculation and quartile definitions. It never
writes to the existing expression-analysis directory. If an explicitly supplied
`--bins-tsv` is missing, the runner stops instead. A successfully reconstructed
table is published as `<out-dir>/tables/tss_expression_bins.tsv`. If the default
source remains missing, each new submission rebuilds the table; supply that saved
table with `--bins-tsv` to reuse it explicitly.

The expression code's mean TPM pools the twelve RNA samples at LPS 0, 5, 10 and
15 minutes. These fixed bins apply to all four Fiber-seq samples. No quartiles
are recalculated from selected molecules or separately by timepoint.

## Molecule selection and signal definition

The sampler verifies each retained gene against the canonical BED by gene ID,
gene name, chromosome, TSS and strand. BED annotations must say `protein_coding`
and contain `Ensembl_canonical`. As in the existing expression code, the BED's
20 bp interval defines `tss = start + 10` in zero-based coordinates.
`not_expressed` and other bins are excluded.

Windows are defined relative to the gene strand. For the TSS-relative interval
`[s, e)` (`--window-start`/`--window-end` of `01_sample_tss_molecules.py`,
negative = upstream, end exclusive):

- `+` strand: genomic interval `[tss + s, tss + e)`;
- `-` strand: genomic interval `[tss - e + 1, tss - s + 1)`.

`canonical_genes()` computes this genomic interval per gene and stores it as
`window_start`/`window_end` next to the offsets `window_offset_start`/
`window_offset_end`. The tabix scan, interval merging and full-span
eligibility (`read_start <= window_start` and `read_end >= window_end`) run on
genomic coordinates; the pool scan anchors on the sorted genomic window starts,
which is the same test on both strands because every window shares one width.
Span coverage follows the existing BED extraction representation; BED12 does
not provide a CIGAR with which to assess internal alignment gaps.

`m6a_path()` and `merge_intervals()` are reused from
`code/accessibility/expr_access/02_tss_m6a_profiles.py`. Tabix scans merged
windows to build a SQLite pool of fully spanning reads, without parsing the
m6A arrays or calculating ACFs. One global `read_id` primary key deduplicates
across chromosomes and samples. For a read with multiple eligible alignments,
the longest eligible alignment is retained. Seeded hash ordering resolves equal
length alignments and multiple eligible genes. Each physical read is consequently
assigned exactly one gene and one expression bin before sampling. The original
BED record digest identifies its exact representation for subsequent extraction.

A separate, independently keyed hash of `(seed, read_id)` supplies a reproducible
random ordering within each bin. Taking its first `--per-bin` records samples
without replacement. Default seed is 0. Pool counts refer to these **disjoint,
uniquely assigned** molecules, not duplicate read–gene memberships; multi-TSS reads
are not counted in multiple bins. If any assigned pool is too small, the sampler
reports all four available counts and stops without oversampling or calculating
ACFs. There are no gene-level or timepoint-level sampling quotas.

Only selected molecules are fetched again to build the binary matrix. Interior
BED12 m6A blocks produce ones; first/last sentinel blocks are dropped, and all
other covered bases are zero. Each row is built in genomic order and then
**reversed for `-` strand genes**, so column `j` always equals TSS-relative
offset `s + j` (`orientation = "transcriptional"` in the metadata). A check
confirms on both strands that the column for offset 0 maps back to genomic
base `tss`. No reference A/T
filter, smoothing, tapering, detrending, binning or FFT is introduced.
`m6a_count` is the number of ones **inside the window**. Zero-call and all-one
signals remain in the balanced sample and are flagged downstream.

The ACF is invariant to sequence reversal, so per-read ACF values change only
through which reads and bases fall in each strand-oriented window. Heatmaps,
the mean oriented binary profile, sliding windows and any other position-based
output do depend on the orientation.

Metadata includes `read_id`, `sample`, `timepoint`, numeric `lps_minutes`, gene ID
and name, `chrom`, `tss`, gene `strand`, `read_strand`, `gene_type`, canonical flag,
`mean_tpm`, `expr_bin`, `m6a_count`, read coordinates, genomic
`window_start`/`window_end`, `window_offset_start`/`window_offset_end`,
`orientation` and the source-record hash. `row_index` is zero-based and records
the original matrix order in the saved tables.

## Unchanged scientific functions

`02_compute_tss_autocorrelations.py` imports `autocorrelations()` from
`../03_compute_autocorrelations.py`:

```text
ACF(k) = sum((x[t]-mean(x)) * (x[t+k]-mean(x))) / (N * var(x))
```

Every ACF in this folder (fixed windows, per-molecule NRL, sliding windows and
the A/T control) is computed by this parent function on the unbinned per-base
0/1 signal at every integer lag from 0 to `max_lag` (1 bp resolution), with no
binning, smoothing, tapering, detrending, downsampling or FFT truncation. Step
02 retains all nonnegative lags of the window by default (0–1999 for `2000_tss`,
0–899 for the one-sided regions) and checks that each saved ACF row has exactly
`max_lag + 1` columns. Valid reads must have lag 0 equal to one within `1e-10`.
Zero-variance rows have all-NaN ACFs and a false validity flag.

`03_leiden_cluster_autocorrelations.py` imports `cluster_profiles()` from
`../04_cluster_autocorrelations.py` and calls it once for the entire sample:

```text
ACF → PCA (50 PCs) → correlation-distance kNN (10 neighbors)
    → Leiden (resolution 0.4, seed 0) → UMAP for visualization
```

The parent's directed, weighted `leidenalg` workflow, PCA centering, feature
handling and small-data safeguards are retained. Zero-variance or otherwise
unclusterable rows keep their status and are displayed as `Unclustered`. Thus the
balanced sample always contains 10,000 molecules, while the number entering
Leiden can be smaller. UMAP is never used as clustering input.

`04_prepare_plot_tables.py` summarizes clusters and expression bins without
changing the profiles. It reuses `repeat_peak()` from
`../05_summarize_autocorrelations.py` for the optional positive local 140–250 bp
peak. This is descriptive, not a periodicity significance test; peaks are missing
if the configured lag range cannot cover the complete peak band.

## Single-molecule footprints

`04b_extract_tss_footprints.py` runs after step 04, before step 05, alongside
the existing NRL steps. It reads the temporary oriented `binary_m6a.npy`,
`row_ids.npy` and `plot_metadata.tsv`; it never resamples or clusters.
FiberHMM inputs are tabix-indexed BED12 files at
`macrophage_project/FiberHMM/extract/firehmm_{footprint,tf}/LPS0/`
`LPS0_hmm_extracted_{footprint,tf}_chr22.bed.gz` (analogously for each
sample/chromosome; remove underscores from `LPS_0` to obtain `LPS0`).
Every BED block is real, including first/last; zero-block records contribute no
calls. A unique read-ID/sample/chromosome record is accepted, including records
trimmed to the outermost calls. If several distinct records match, the exact
sampled read span must resolve them to one; otherwise they are reported as
ambiguous and excluded.

Nucleosome sizes **>90 bp** and TF sizes **<60 bp** are retained; 60–90 bp
is excluded. Clip in genomic coordinates, map the included interval endpoints
with the existing sampler coordinate helper, then restore half-open column
intervals. This preserves transcriptional orientation on both strands.
The original-sized category matrix uses 0 = white/no call, 1 = dark grey
`#4d4d4d`/nucleosome, 2 = orange `#f16913`/TF, 3 = purple `#800080`/m6A,
in increasing precedence. m6A is copied directly from the oriented binary matrix.
Checks enforce shape, codes, per-row counts, **category == 3 iff binary == 1
at every base**, endpoint round trips and mirrored minus-strand masks.

Published tables: `footprint_records.tsv.gz` (read_id, row_index, track,
start, end, size; clipped genomic half-open coordinates and original call size),
`footprint_matches.tsv` (per-source match/ambiguity status),
`footprint_validation.tsv` (match counts, fraction with any retained footprint,
and validation results), and `footprint_minus_strand_checks.tsv` (example
mapped intervals). The row-major uint8 `intermediate/footprint_categories.bin`
and its shape TSV are temporary and are included in the R plotting signature.

`08_single_fiber_footprints.pdf` is the 16th figure: one raster at one pixel per
clustered molecule per base, in exactly the ACF heatmap order. Cluster,
expression-bin and timepoint strips, axes, boundaries and labels are vector.
The x axis is position relative to TSS in the transcription direction
(negative = upstream); the dashed TSS line appears only inside the window.

## Per-molecule NRL and regularity

`tss_nrl.py` holds pure functions (no I/O) that take one ACF row and its window
width `N`. `04b_molecule_nrl.py` applies them to the fixed-window ACFs of step
02; `04c_sliding_windows.py` applies them to every sliding window.

**NRL** follows the SAMOSA secondary-peak scan (Abdulhay et al. 2020, eLife):

- lags below `--min-lag` (default 60) are ignored to skip the lag-0 shoulder;
- `scipy.signal.find_peaks` with `prominence = min_prominence` finds local maxima;
- the NRL is the lag of the first positive local maximum in
  [`--nrl-min`, `--nrl-max`] (default 120–300 bp);
- the second peak is the first positive local maximum within ±50 bp of 2 × NRL
  that lies inside the available lags;
- the parent `repeat_peak()` 140–250 bp result is reported for comparison
  (`repeat_peak_lag`, `repeat_peak_value`).

`nrl_status` is `ok`, or `no_peak` (no in-band local maximum passes the
prominence), `peak_negative` (in-band maxima exist but none is positive),
`window_too_short` (fewer available lags than `nrl_min`) or `zero_variance`
(undefined ACF); NRL and the second peak are NaN unless `ok`.

**`min_prominence` is calibrated from the data, not fixed.** A per-base
single-molecule ACF has sampling noise of about 1/√N (0.02 for N = 2,000, 0.03
for N = 500), so a small fixed value such as 0.01 makes every structureless row
yield a "peak" and places the first local maximum of real reads on a noise bump
at the lower band edge. `tss_nrl.calibrate_prominence()` therefore:

1. draws up to `--n-null` (default 2,000) rows of the matrix being analysed
   (seeded);
2. shuffles the 0/1 positions within each row (same m6A count, no spatial
   structure) and computes the parent ACF of the shuffled rows;
3. records, per row, the largest prominence of any positive local maximum in
   [`nrl_min`, `nrl_max`] at lags ≥ `min_lag` (0 when there is none);
4. sets `min_prominence` to the `--prominence-quantile` (default 0.95) of these
   null maxima, so at most about 5% of structureless rows of that window
   length pass.

The calibration is done at the ACF length actually analysed: on the fixed
window in each region's `04b` (N = 2,000 or 900) and on random (read, window)
slices of width `--win-width` in `04c` (N = 500, shared by the m6A windows and
the A/T control). The value and the null quantiles are written to
`tables/nrl_prominence_calibration.tsv` and
`tables/sliding_prominence_calibration.tsv`, repeated as a `min_prominence`
column in the summary tables and logged. `--min-prominence` sets a fixed
override instead. On the macrophage 2 kb windows the calibrated value is
about 0.13; at that level about half of the molecules have a detected peak
with a median NRL near 187 bp, whereas 0.01 returns a 124 bp "NRL" for almost
every molecule.

**Regularity from the ACF decay** (Baldi et al. 2018, Mol Cell: the longer the
ACF takes to dampen, the more regular the array). The parent ACF divides by `N`,
so it decays as `(N − k)/N` even for a perfectly periodic signal. For the decay
metrics only, the ACF is rescaled by `N/(N − k)` (the unbiased form) and
restricted to lags ≤ `--decay-max-lag` (default `floor(N/2)`); the raw ACF is
unchanged everywhere else (saved matrices, heatmaps, mean curves, clustering,
NRL scan). On the rescaled ACF:

1. `peak1_height`/`peak2_height` are the rescaled heights at the NRL and the
   second peak, and `peak_ratio = peak2_height / peak1_height`
   (`nrl_peak_height` in the tables is the raw ACF at the NRL lag; the second
   peak may lie beyond `decay_max_lag` in short windows, where its rescaled
   height is noisy);
2. a damped cosine `A·exp(−k/λ)·cos(2πk/P + φ) + c` is fitted over
   [`min_lag`, `decay_max_lag`] with `scipy.optimize.curve_fit`, bounds
   `P ∈ [nrl_min, nrl_max]` and `0 < λ ≤ 10 × decay_max_lag`, initialising `P`
   from the NRL peak when one exists; reported as `decay_length_bp` (λ),
   `fit_period_bp` (P), `fit_r2` and `fit_converged`, with NaN on failure. A
   decay length at the upper bound (2,500 bp for 500 bp windows, 10,000 bp for
   the 2 kb window) is censored: no decay is resolvable within the fitted
   lags. In 500 bp windows most reads sit at this bound, so compare medians
   and the fraction below the bound rather than interpreting single values;
3. `damping_lag_bp` is the smallest lag after which the envelope stays below
   `--flat-threshold` (default 0.05) for all remaining lags up to
   `decay_max_lag`; the envelope is a centred rolling maximum of |ACF| over one
   period (the NRL if found, else 190 bp) on lags 1..`decay_max_lag` (lag 0 is
   identically one and excluded); NaN when the envelope never flattens. Note
   that the rescaled ACF noise (≈ 1/√N × N/(N−k), i.e. 0.045–0.09 at
   N = 500) rarely stays below an absolute 0.05 within 250 lags, so in the
   sliding windows the damping lag is NaN for most reads and the fixed-window
   values cluster near `decay_max_lag`; a threshold relative to 1/√N would be
   needed for short windows and is left to `--flat-threshold`.

`04b_molecule_nrl.py` writes `tables/nrl_per_molecule.tsv` (one row per sampled
molecule: `read_id, sample, timepoint, lps_minutes, gene_id, gene_name, chrom,
tss, strand, expr_bin, mean_tpm, m6a_count, window, window_offset_start,
window_offset_end, cluster, nrl_bp, nrl_status, nrl_peak_height, peak2_lag_bp,
peak2_height, peak_ratio, repeat_peak_lag, repeat_peak_value, decay_length_bp,
fit_period_bp, fit_r2, fit_converged, damping_lag_bp`), and
`tables/nrl_summary_by_bin.tsv` / `tables/nrl_summary_by_cluster.tsv` with n,
medians, IQRs and the fraction of molecules with a detected NRL peak.

## Sliding windows

`04c_sliding_windows.py` runs in each array task on that task's own sampled
molecules and oriented binary matrix. Nothing is resampled, so every molecule
contributes to every window of its region. The span is always the task's own
window and never exceeds 2,000 bp. Window starts run from the span start to
`span_end − win_width` inclusive (`--win-width` default 500, `--win-step`
default 100; 200 is allowed), and every window must lie inside the span:

| Region (array task) | Span | Window starts (width 500, step 100) | Windows |
|---|---|---|---|
| `2000_tss` (`span_2kb`) | `[-1000, +1000)` | −1000 … +500 | 16 |
| `upstream_1000_tss` (`upstream`) | `[-1000, -100)` | −1000 … −600 | 5 |
| `downstream_1000_tss` (`downstream`) | `[+100, +1000)` | +100 … +500 | 5 |

The one-sided regions therefore never touch `[-100, +100)`. For each window
the oriented row is sliced, the parent ACF is computed with
`max_lag = win_width − 1` (lags 0–499 at 1 bp) and the `tss_nrl.py` metrics are
applied with decay lags ≤ `floor(N/2)`. Computation is chunked; no reads ×
windows × lags array is stored.

Windows of `2000_tss` that overlap `[-100, +100)` carry `overlaps_ndr_core =
True` (six windows, centres −250 … +250 at the defaults). Their "NRL" partly
reflects the edge of the nucleosome-depleted region rather than nucleosome
spacing (Clarkson et al. 2019, NAR), so they are shaded in the plots and left
out of the paired mirror-window tests.

**Paired versus unpaired.** Within `2000_tss`, the same molecules appear in
every window, so window-to-window comparisons are paired by read (mirror
windows at centres ∓750, ∓650, ∓550, ∓450, ∓350: Wilcoxon signed-rank).
`upstream_1000_tss` and `downstream_1000_tss` are different reads, so upstream
versus downstream is an unpaired Wilcoxon rank-sum test on gene-level medians
per expression bin. BH correction is applied within each test family across
bins and metrics (`nrl`, `peak_fraction`, `decay`, `damping`, `m6a`).

**A/T control.** The same windows of the same molecules are analysed on the
indicator of the strand-oriented reference base being A or T (fetched from the
reference FASTA with the same reversal for `-` genes; the indicator is
unchanged by complementation). Its median NRL, fraction with a peak and median
decay/damping per window are in `tables/sliding_at_control.tsv` and drawn as
dashed grey lines.

Limitations: 500 bp holds only about 2–3 nucleosome periods, so single-molecule
NRL and decay estimates are noisy; compare medians between groups. At lag `k`
only `N − k` base pairs contribute to the ACF.

Output tables per region: `tables/sliding_nrl_long.tsv.gz` (one row per read ×
window with the read metadata plus `region, win_start, win_end, win_center,
overlaps_ndr_core, m6a_count_win, nrl_bp, nrl_status, nrl_peak_height,
peak2_height, peak_ratio, decay_length_bp, damping_lag_bp, fit_r2`),
`tables/sliding_summary.tsv` (per expression bin × timepoint × window, plus
`timepoint = all`: n, median NRL and IQR, fraction with a detected peak, median
decay length, median damping lag, mean m6A density),
`tables/sliding_mean_acf.tsv.gz` (mean ACF per bin × window × lag for the m6A
signal and the A/T control), `tables/sliding_at_control.tsv` and
`tables/sliding_prominence_calibration.tsv`.

## Results and plots

| Directory | Contents |
|---|---|
| `<region>/tables/` | Eligible genes and pool counts, sampled/clustered metadata, full ACF plot matrix and summaries, per-molecule NRL tables, prominence calibrations, sliding-window tables; reconstructed expression table when needed |
| `<region>/plots/` | Sixteen final R-generated PDF figures |
| `combined/tables/` | `sliding_tests.tsv`, `sliding_bootstrap_ci.tsv`, `sliding_combined_summary.tsv` |
| `combined/plots/` | Ten cross-region PDF figures |
| `logs/` | SLURM `slurm_<job-id>_<task>.out/.err` and `slurm_<job-id>_combined.out/.err` |

The runner publishes tables and plots only after all stages and their checks
succeed, replacing each final file atomically. SLURM logs use the fixed paths in
the runner's `#SBATCH` directives, including for `--out-dir` runs;
override them with `sbatch --output ... --error ...` if needed. No persistent
`intermediate/`, `validation/`, `inputs/`, plot manifest or lock file is
created by the runner; existing artifacts from older runs are not deleted.

All matrices and molecule metadata retain their original sampled order.
`tables/plot_metadata.tsv` adds a separate one-based `heatmap_rank` ordered by
numeric Leiden label, descending m6A count within that label, and original row
index to break ties. `Unclustered` appears last. The wide
`tables/acf_heatmap.tsv.gz` retains every sampled molecule and every lag in
original matrix order; sorting is for display only.

Per region, `05_plot_tss_autocorrelation.R` creates sixteen PDFs:

- Strand-oriented single-molecule footprint raster (`08_single_fiber_footprints`).
- All-molecule ACF heatmap with cluster annotations and within-cluster m6A ranks.
- Four expression-bin heatmaps in a shared-label panel.
- Cluster mean ACF curves and expression-bin mean ACF curves (from lag 25 bp).
- Cluster fractions within each expression bin and reciprocal bin composition.
- m6A-count distributions by cluster and expression bin, plus pooled cluster distributions.
- UMAP colored by cluster, expression bin and m6A count.
- Violin/box plots of NRL, decay length (log10) and damping lag by expression
  bin (`08_nrl_metrics_by_expression`) and by cluster (`08b_nrl_metrics_by_cluster`).
- Per-cluster NRL histograms in the SAMOSA Fig. 3D style (`09_nrl_histogram_by_cluster`).
- Fraction of molecules with a detected NRL peak per bin and per cluster (`10_nrl_peak_fraction`).

`08_plot_sliding_combined.R` creates, in `combined/plots/`, one figure per
sliding metric (median NRL, fraction with a detected peak, median decay length
with median damping lag, m6A density), each with the 2 kb-span panel (NDR
overlapping windows shaded), the one-sided panel (upstream and downstream on
one axis with the gap at `[-100, +100)`), the A/T control as a dashed grey line,
one line per expression bin and bootstrap 95% CIs computed by resampling genes
(not reads); a `_by_timepoint` version of each; and mean-ACF heatmaps (lags
0–499 × window centre) per expression bin for the 2 kb and one-sided runs. The
x axis is the window centre relative to the TSS with upstream on the left.

Heatmaps embed full-resolution rasters with no interpolation. They include lag 0;
their symmetric color range uses the 99th percentile of absolute nonzero-lag ACFs
to keep nonzero-lag structure visible. Values outside that display range saturate
and the legend says so; saved ACF values are unchanged.

**Figures show clustered molecules only** where clusters are involved.
`Unclustered` molecules (zero variance or otherwise unclusterable) stay in
`tables/` with that label, but are dropped from the heatmaps, boxplots, UMAPs,
composition bars and per-cluster NRL plots; the composition bars are
renormalized over clustered reads within each expression bin. In the TSVs,
`cluster_composition.tsv` keeps the original denominator of **all sampled reads
within each expression bin** including the explicit `Unclustered` category.
Per-expression-bin NRL plots and all sliding-window outputs use every sampled
molecule.

Final plotting uses base R; Python produces matrices and tables only.
The individual worker scripts retain their working-file behavior when invoked
directly. For example, `05_plot_tss_autocorrelation.R --out-dir PATH --png` can
generate PNGs while the working directory still contains the footprint binary
and shape TSV, but also writes validation and manifest files there. Published
tables alone are insufficient to rerun the full plotting step after cleanup.

## Runtime

Runtime defaults are the existing
`/project/spott/cshan/envs/Jupyter-notebook/bin/python` and
`/software/R-4.4.1-el8-x86_64/bin/Rscript`; override via `TSS_PYTHON` and
`TSS_RSCRIPT`. Python requires NumPy, pandas, pysam, SciPy, Scanpy, anndata,
scikit-learn, umap-learn, igraph and leidenalg. Only expression-table recovery
additionally requires R `data.table`. No environment installation or parent-code
modification is required.

Checks are embedded in the numbered analysis/plot scripts. They enforce
canonical protein-coding metadata, expressed-only bins, full spans, unique read
IDs, balanced sampling, window-width binary rows, m6A counts, matrix/metadata
order, transcriptional orientation (offset 0 ↔ genomic `tss` on both strands),
lag-zero normalization, `max_lag + 1` ACF columns, invalid row exclusions,
successful joint Leiden/UMAP, global labels, heatmap ordering, composition
denominators, NRL status consistency, sliding-window counts inside the span,
read-ID identity between the sliding table and the sampled molecules, and
complete R plot generation.
