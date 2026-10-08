# Macrophage canonical-TSS autocorrelation

This pipeline samples molecules **before** constructing binary signals or
calculating autocorrelations. The full run selects exactly 10,000 distinct raw
read IDs: 2,500 each from `Q1_low`, `Q2`, `Q3`, and `Q4_high`. All expression bins
are clustered together. Parent scripts are imported and remain unchanged: the
ACF, clustering and group summaries use exactly the method of the parent
analysis in `..` (unsmoothed per-base 0/1 signal, parent `autocorrelations()`,
parent `cluster_profiles()`, parent `repeat_peak()`).

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

Files in this folder:

| File | Role |
|---|---|
| `00_prepare_expression_bins.R` | Rebuilds the expression table from the unchanged expression code when it is absent |
| `01_process_tss_molecules.py` | Every Python processing stage, run in order in one process: `sample`, `acf`, `cluster`, `tables`, `footprints`, `nrl` (plus the shared helpers and the NRL/null functions) |
| `02_plot_tss_autocorrelation.R` | The seventeen per-region figures and their checks |
| `03_run_tss_autocorrelation.sh` | SLURM array runner: one task per region, atomic publishing |
| `tss_autocorrelation.Rmd` | Report embedding the published figures and tables |

## Submit

Run from any working directory after creating the SLURM log directory:

```bash
mkdir -p /project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss/logs
```

Full 10,000-molecule analysis on chromosomes 1–22, X and Y for all three regions:

```bash
sbatch /project/spott/cshan/fiber-seq/code/clustering_methods/auto_correlation/tss/03_run_tss_autocorrelation.sh
```

One region only, e.g. `downstream_1000_tss` (array task 2):

```bash
sbatch --array=2 /project/spott/cshan/fiber-seq/code/clustering_methods/auto_correlation/tss/03_run_tss_autocorrelation.sh
```

The `.sh` runner contains `#SBATCH` directives and requires a SLURM allocation.
It requests `pi-spott` / `bigmem`, one node, two CPUs, 300 GB RAM and 30 hours;
an earlier full run used about 1 GB and 11 minutes per task, and the group null
adds about 0.8 GB for the shuffled ACF pool of `2000_tss`. All numerical-library
and Numba thread counts are capped at one. The disk-backed pool avoids keeping
all reads or their m6A signals in memory.

The full run writes to
`/project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss/<region>/`.
Runner options: `--out-dir`
(replaces the region folder), `--bins-tsv`,
`--canonical-bed`, `--ft-root`, `--ref` (reference FASTA for the A/T-restricted
null shuffle), `--per-bin`, `--seed`, `--max-lag`, `--n-pcs`, `--n-neighbors`,
`--resolution`, `--chrom chr21,chr22`, and the NRL options `--nrl-min`,
`--nrl-max`, `--min-lag`, `--min-prominence`, `--prominence-quantile`,
`--n-null`, `--null-shuffle at|all`, `--n-shuffles-per-molecule` and
`--n-null-groups`. Their values are logged at the start of each task. Omit
`--chrom` for the full genome. `--max-lag` is **inclusive** and defaults to all
lags of the window (width - 1); the joint clustering pipeline requires at least
three lag features. SLURM logs are shared by all regions as
`logs/slurm_<job-id>_<task>.out`.

## Inputs and expression-table recovery

Defaults:

- Canonical BED (1 bp rows, start = 0-based TSS):
  `/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed`.
- Expression table:
  `/project/spott/cshan/fiber-seq/macrophage_project/expr_access/tables/tss_expression_bins.tsv`.
- Indexed Fiber-seq BED12 files under
  `macrophage_project/FiberHMM/extract/ft_result_dir/{LPS_0,LPS_5,LPS_10,LPS_15}/extracted_results/m6a_by_chr/`.
- Reference FASTA for the A/T-restricted null shuffle:
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
1 bp interval defines `tss = start` in zero-based coordinates.
`not_expressed` and other bins are excluded.

Windows are defined relative to the gene strand. For the TSS-relative interval
`[s, e)` (`--window-start`/`--window-end` of `01_process_tss_molecules.py`,
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
the mean oriented binary profile, the A/T null mask and any other
position-based output do depend on the orientation.

Metadata includes `read_id`, `sample`, `timepoint`, numeric `lps_minutes`, gene ID
and name, `chrom`, `tss`, gene `strand`, `read_strand`, `gene_type`, canonical flag,
`mean_tpm`, `expr_bin`, `m6a_count`, read coordinates, genomic
`window_start`/`window_end`, `window_offset_start`/`window_offset_end`,
`orientation` and the source-record hash. `row_index` is zero-based and records
the original matrix order in the saved tables.

## Unchanged scientific functions

Stage `acf` of `01_process_tss_molecules.py` imports `autocorrelations()` from
`../03_compute_autocorrelations.py`:

```text
ACF(k) = sum((x[t]-mean(x)) * (x[t+k]-mean(x))) / (N * var(x))
```

Every ACF in this folder (fixed windows, per-molecule NRL and the shuffled
nulls) is computed by this parent function on the unbinned per-base
0/1 signal at every integer lag from 0 to `max_lag` (1 bp resolution), with no
binning, smoothing, tapering, detrending, downsampling or FFT truncation. Step
02 retains all nonnegative lags of the window by default (0–1999 for `2000_tss`,
0–899 for the one-sided regions) and checks that each saved ACF row has exactly
`max_lag + 1` columns. Valid reads must have lag 0 equal to one within `1e-10`.
Zero-variance rows have all-NaN ACFs and a false validity flag.

Stage `cluster` imports `cluster_profiles()` from
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

Stage `tables` summarizes clusters and expression bins without
changing the profiles, following the parent `summarize()`: mean (and median)
ACF per group, number of molecules, `fraction_clustered`, `mean_m6a_call_fraction`
(m6A calls over window width), and `repeat_peak()` from
`../05_summarize_autocorrelations.py` for the strongest positive local
140–250 bp peak of the group mean ACF. This is descriptive, not a periodicity
significance test (see the group null below); peaks are missing if the
configured lag range cannot cover the complete peak band. `cluster_summary.tsv`
and `cluster_composition.tsv` carry the zero-variance and unclustered counts
(`n_zero_variance`, `n_unclustered`, `n_bin_zero_variance`, `n_bin_unclustered`),
and `unclustered_counts.tsv` breaks them down per expression bin and timepoint.

## Single-molecule footprints

Stage `footprints` runs after the tables stage and before the NRL stage. It
reads the temporary oriented `binary_m6a.npy`,
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

`08_single_fiber_footprints.pdf` is one raster at one pixel per
clustered molecule per base, in exactly the ACF heatmap order. Cluster,
expression-bin and timepoint strips, axes, boundaries and labels are vector.
The x axis is position relative to TSS in the transcription direction
(negative = upstream); the dashed TSS line appears only inside the window.

## Per-molecule NRL

Stage `nrl` of `01_process_tss_molecules.py` applies the pure NRL functions
(`find_nrl`, `molecule_metrics`) and the shuffled-null functions (`at_mask`,
`shuffle_rows`, `calibrate_prominence`, `shuffle_pool`, `group_null`) of the same
script to the fixed-window ACFs of stage `acf`.

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
for N = 900), so a small fixed value such as 0.01 makes every structureless row
yield a "peak" and places the first local maximum of real reads on a noise bump
at the lower band edge. `calibrate_prominence()` therefore, separately in each
region:

1. draws up to `--n-null` (default 2,000) of the region's sampled molecules
   (seeded);
2. shuffles each molecule's m6A calls **only among the A/T positions of its own
   strand-oriented reference window** (`--null-shuffle at`, the default; the
   reference sequence is fetched from hg38.fa with the same reversal for `-`
   genes, G/C positions stay 0 and the number of calls is unchanged).
   `--null-shuffle all` restores a shuffle over every position;
3. computes the parent ACF of the shuffled molecules and records, per molecule,
   the largest prominence of any positive local maximum in [`nrl_min`,
   `nrl_max`] at lags ≥ `min_lag` (0 when there is none);
4. sets `min_prominence` to the `--prominence-quantile` (default 0.95) of these
   null maxima, so at most about 5% of structureless molecules of that window
   length pass.

The value, the null median and 90th/95th/99th percentiles, the shuffle mode,
window length and number of null molecules are written to
`tables/nrl_prominence_calibration.tsv`, the per-molecule null maxima to
`tables/nrl_null_prominences.tsv` (plotted in `11_null_prominence_distribution`),
repeated as a `min_prominence` column in the summary tables and logged.
`--min-prominence` sets a fixed override instead and is the only way to use a
fixed value. Checks confirm that every shuffled vector keeps its source
molecule's m6A count and, in `at` mode, calls only A/T positions.

**Group-level null for the mean-ACF peak** (per expression bin and per
cluster). Every sampled molecule is shuffled `--n-shuffles-per-molecule`
(default 5) times with the same A/T-restricted shuffle and parent ACF (seeded;
the pool is kept in memory, about 0.8 GB for `2000_tss`). For each group,
`--n-null-groups` (default 200) null means are built by drawing, with
replacement, a group-sized set of shuffled ACFs from that group's own valid
molecules' shuffles and averaging them; a check confirms that no draw comes
from outside the group. The observed statistic is the prominence of the most
prominent positive local maximum of the group's mean ACF inside the NRL band;
its empirical p-value is `(1 + #null ≥ observed) / (1 + n_null)`.
`tables/group_peak_tests.tsv` holds, per group, `n_valid`, the observed peak
lag, height and prominence, the number of null means, how many reach the
observed value, the p-value and the null median and 95th percentile.

**Removed metrics.** The damped-cosine decay length, fit period, fit R² and
damping lag (Baldi et al. 2018), together with the unbiased `N/(N − k)`
rescaling and the `--decay-max-lag` / `--flat-threshold` options, are no longer
computed. `peak_ratio` is now the raw second-peak height over the raw NRL-peak
height.

Stage `nrl` writes `tables/nrl_per_molecule.tsv` (one row per sampled
molecule: `read_id, sample, timepoint, lps_minutes, gene_id, gene_name, chrom,
tss, strand, expr_bin, mean_tpm, m6a_count, window, window_offset_start,
window_offset_end, cluster, nrl_bp, nrl_status, nrl_peak_height, peak2_lag_bp,
peak2_height, peak_ratio, repeat_peak_lag, repeat_peak_value`), and
`tables/nrl_summary_by_bin.tsv` / `tables/nrl_summary_by_cluster.tsv` with
`n`, `n_acf_valid`, `n_zero_variance`, `n_unclustered`, `n_with_nrl_peak`,
`fraction_with_nrl_peak`, median NRL, IQR, `mad_nrl_bp`
(`scipy.stats.median_abs_deviation`, `nan_policy='omit'`), `median_peak_ratio`,
the group mean-ACF peak lag, height, prominence and p-value, and
`min_prominence`.

## Archived sliding-window results

The 500 bp sliding-window analysis (per-window NRL, decay and damping, A/T
sequence control, cross-region combined figures and tests) has been removed
from the pipeline. The first run of this version moves each region's
`sliding_nrl_long.tsv.gz`, `sliding_summary.tsv`, `sliding_mean_acf.tsv.gz`,
`sliding_at_control.tsv` and `sliding_prominence_calibration.tsv` from
`<region>/tables/` to `<root>/_archive_sliding_<date>/<region>/tables/` and
logs every move; nothing is deleted. The `combined/` folder was removed
earlier, so it has no archive copy.

## Results and plots

| Directory | Contents |
|---|---|
| `<region>/tables/` | Eligible genes and pool counts, sampled/clustered metadata, full ACF plot matrix and summaries, unclustered counts, per-molecule NRL tables, prominence calibration and null prominences, group peak tests; reconstructed expression table when needed |
| `<region>/plots/` | Seventeen final R-generated PDF figures |
| `logs/` | SLURM `slurm_<job-id>_<task>.out/.err` |
| `_archive_sliding_<date>/<region>/tables/` | The five `sliding_*` tables of the removed sliding-window analysis, moved there by the first run of this version (never deleted) |

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

Per region, `02_plot_tss_autocorrelation.R` creates seventeen PDFs:

- Strand-oriented single-molecule footprint raster (`08_single_fiber_footprints`).
- All-molecule ACF heatmap with cluster annotations and within-cluster m6A ranks.
- Four expression-bin heatmaps in a shared-label panel.
- Cluster mean ACF curves and expression-bin mean ACF curves (from lag 25 bp).
- Cluster fractions within each expression bin and reciprocal bin composition.
- m6A-count distributions by cluster and expression bin, plus pooled cluster distributions.
- UMAP colored by cluster, expression bin and m6A count.
- Violin/box plots of NRL by expression bin (`08_nrl_metrics_by_expression`) and
  by cluster (`08b_nrl_metrics_by_cluster`).
- Per-cluster NRL histograms in the SAMOSA Fig. 3D style (`09_nrl_histogram_by_cluster`).
- Fraction of molecules with a detected NRL peak per bin and per cluster (`10_nrl_peak_fraction`).
- Null distribution of the largest in-band peak prominence of the A/T-shuffled
  molecules with the calibrated threshold marked (`11_null_prominence_distribution`).

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
Per-expression-bin NRL plots and tables use every sampled molecule.

Final plotting uses base R; Python produces matrices and tables only.
The individual worker scripts retain their working-file behavior when invoked
directly. For example, `02_plot_tss_autocorrelation.R --out-dir PATH` writes
seventeen PDFs while the working directory still contains the footprint binary
and shape TSV, plus validation and manifest files there. The wrapper publishes
only the PDFs from `plots/`. Published tables alone are insufficient to rerun
the full plotting step after cleanup.

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
denominators, NRL status consistency, null shuffles that keep each molecule's
m6A count and (in `at` mode) call only A/T positions, group nulls drawn only
from the group's own molecules, complete R plot generation, and, after
publishing, the absence of any sliding-window, decay or damping file name,
table column or report reference.
