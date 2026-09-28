# Macrophage canonical-TSS autocorrelation

This pipeline samples molecules **before** constructing binary signals or
calculating autocorrelations. The full run selects exactly 10,000 distinct raw
read IDs: 2,500 each from `Q1_low`, `Q2`, `Q3`, and `Q4_high`. All expression bins
are clustered together. Parent scripts are imported and remain unchanged.

## Submit

Run from any working directory after creating the SLURM log directory:

```bash
mkdir -p /project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss/logs
```

Small real-data test (100 per bin, chromosome 22, all four LPS timepoints):

```bash
sbatch /project/spott/cshan/fiber-seq/code/clustering_methods/auto_correlation/tss/06_run_tss_autocorrelation.sh --test
```

Full 10,000-molecule analysis on chromosomes 1–22, X and Y:

```bash
sbatch /project/spott/cshan/fiber-seq/code/clustering_methods/auto_correlation/tss/06_run_tss_autocorrelation.sh
```

The analysis was validated on 2026-09-28: SLURM job `59660099` completed the chromosome-22 test
with exactly 100 unique molecules per bin, 400 × 2,000 binary/ACF matrices,
400 valid profiles, four joint Leiden clusters, and all eleven R PDFs. All
32 R validation checks passed; maximum lag-zero error was `4.30e-14`.
The assigned eligible pools contained 8,084 / 7,492 / 11,267 / 13,158 reads
for Q1–Q4. These results describe the earlier runner; its scientific functions
and checks remain unchanged by the temporary-workspace update.

The `.sh` runner contains `#SBATCH` directives and requires a SLURM allocation.
It requests `pi-spott` / `bigmem`, one node, two CPUs, 300 GB RAM and 30 hours.
All numerical-library and Numba thread counts are capped at one. The disk-backed
pool avoids keeping all reads or their m6A signals in memory. The 10,000 × 2,000
binary and double-precision ACF matrices occupy approximately 20 MB and 160 MB;
working memory also includes the exact neighbor graph calculation, R matrix
copies, and full-resolution plot rasters. Resource requests reflect the current
runner settings.

The full run writes to
`/project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss/`;
`--test` writes to its `test/` subdirectory. Put `--test` before explicit overrides.
Supported runner options include `--out-dir`, `--bins-tsv`, `--canonical-bed`,
`--ft-root`, `--per-bin`, `--seed`, `--max-lag`, `--n-pcs`, `--n-neighbors`,
`--resolution`, and `--chrom chr21,chr22`. Omit `--chrom` for the full genome.
`--max-lag` is **inclusive** and defaults to 1999; the joint clustering pipeline
requires at least three lag features, so use a maximum of at least 2.

## Inputs and expression-table recovery

Defaults:

- Canonical BED:
  `/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed`.
- Expression table:
  `/project/spott/cshan/fiber-seq/macrophage_project/expr_access/tables/tss_expression_bins.tsv`.
- Indexed Fiber-seq BED12 files under
  `macrophage_project/FiberHMM/extract/ft_result_dir/{LPS_0,LPS_5,LPS_10,LPS_15}/extracted_results/m6a_by_chr/`.

The specified expression table was absent when this analysis was created.
When that default file is missing, `00_prepare_expression_bins.R` evaluates the
unchanged `code/accessibility/expr_access/01_expression_bins.R`, redirecting only
its `OUT_DIR` assignment to `inputs/` inside the temporary workspace. It uses the original STAR counts,
GTF, TPM calculation and quartile definitions. It never writes to the existing
expression-analysis directory. If an explicitly supplied `--bins-tsv` is missing,
the runner stops instead. A successfully reconstructed table is published as
`<out-dir>/tables/tss_expression_bins.tsv`. If the default source remains missing,
each new submission rebuilds the table; supply that saved table with `--bins-tsv`
to reuse it explicitly.

The expression code's mean TPM pools the twelve RNA samples at LPS 0, 5, 10 and
15 minutes. These fixed bins apply to all four Fiber-seq samples. No quartiles
are recalculated from selected molecules or separately by timepoint.

## Molecule selection and signal definition

The sampler verifies each retained gene against the canonical BED by gene ID,
gene name, chromosome, TSS and strand. BED annotations must say `protein_coding`
and contain `Ensembl_canonical`. As in the existing expression code, the BED's
20 bp interval defines `tss = start + 10` in zero-based coordinates.
`not_expressed` and other bins are excluded.

Each window is exactly `[tss - 1000, tss + 1000)`. A BED alignment must satisfy
`read_start <= window_start` and `read_end >= window_end`. Span coverage follows
the existing BED extraction representation; BED12 does not provide a CIGAR with
which to assess internal alignment gaps.

`m6a_path()` and `merge_intervals()` are reused from
`code/accessibility/expr_access/02_tss_m6a_profiles.py`. Tabix scans merged TSS
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

Only selected molecules are fetched again to build the binary matrix. As in the
parent autocorrelation representation, positions run in ascending genomic order
on both gene strands. Interior BED12 m6A blocks produce ones; first/last sentinel
blocks are dropped, and all other covered bases are zero. No reference A/T filter,
strand reversal, smoothing, tapering, detrending, binning or FFT is introduced.
`m6a_count` is the number of ones **inside the 2 kb window**. Zero-call and all-one
signals remain in the balanced sample and are flagged downstream.

Metadata includes `read_id`, `sample`, `timepoint`, numeric `lps_minutes`, gene ID
and name, `chrom`, `tss`, gene `strand`, `read_strand`, `gene_type`, canonical flag,
`mean_tpm`, `expr_bin`, `m6a_count`, read/window coordinates and source-record hash.
`row_index` is zero-based and records the original matrix order in the saved tables.

## Unchanged scientific functions

`02_compute_tss_autocorrelations.py` imports `autocorrelations()` from
`../03_compute_autocorrelations.py`:

```text
ACF(k) = sum((x[t]-mean(x)) * (x[t+k]-mean(x))) / (N * var(x))
```

It retains all nonnegative lags 0–1999 by default, matching the parent's
full-window default. Valid reads must have lag 0 equal to one within `1e-10`.
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

## Results and plots

| Directory | Contents |
|---|---|
| `tables/` | Ten analysis tables: eligible genes and pool counts, sampled/clustered metadata, full ACF plot matrix and summaries; reconstructed expression table when needed |
| `plots/` | Eleven final R-generated PDF figures |
| `logs/` | Standard SLURM `slurm_<job-id>.out` and `slurm_<job-id>.err` logs |

The runner publishes tables and plots only after all stages and their checks
succeed, replacing each final file atomically. SLURM logs use the fixed paths in
the runner's `#SBATCH` directives, including for `--test` or `--out-dir` runs;
override them with `sbatch --output ... --error ...` if needed.

No persistent `intermediate/`, `validation/`, `inputs/`, plot manifest,
`pipeline_<job-id>.log` or `.pipeline.lock` is created by the runner. Existing
artifacts from older runs are not deleted. Jobs submitted before this update,
including job `59661870`, retain the runner script captured at submission.

All matrices and molecule metadata retain their original sampled order.
`tables/plot_metadata.tsv` adds a separate one-based `heatmap_rank` ordered by
numeric Leiden label, descending m6A count within that label, and original row
index to break ties. `Unclustered` appears last. The wide
`tables/acf_heatmap.tsv.gz` retains every sampled molecule and every lag in
original matrix order; sorting is for display only.

R creates eleven PDFs covering the requested views:

- All-molecule ACF heatmap with cluster annotations and within-cluster m6A ranks.
- Four expression-bin heatmaps in a shared-label panel.
- Cluster mean ACF curves and expression-bin mean ACF curves.
- Cluster fractions within each expression bin and reciprocal bin composition.
- m6A-count distributions by cluster and expression bin, plus pooled cluster distributions.
- UMAP colored by cluster, expression bin and m6A count.

Heatmaps embed full-resolution rasters with no interpolation. They include lag 0;
their symmetric color range uses the 99th percentile of absolute nonzero-lag ACFs
to keep nonzero-lag structure visible. Values outside that display range saturate
and the legend says so; saved ACF values are unchanged. Undefined ACF rows are
gray. Aggregate plots show unsmoothed curves, with separate panels including
and omitting lag 0 so that the nonzero-lag structure is visible.

Cluster fractions use **all sampled reads within each expression bin** as their
denominator, including an explicit `Unclustered` category. Cluster mean curves
use their cluster's members; expression-bin curves use all valid ACFs in that
bin, including any valid reads excluded from clustering. Means and medians are
both available in the TSVs.

Final plotting uses base R; Python produces matrices and tables only.
The individual worker scripts retain their working-file behavior when invoked
directly. For example, `05_plot_tss_autocorrelation.R --out-dir PATH --png` can
generate PNGs from the published tables, but also writes validation and manifest
files there. Use a temporary copy of the tables for standalone plotting if you
want to keep only its resulting figures.

## Temporary workspace and validation

Steps 00–05 run in a unique temporary directory under `SLURM_TMPDIR`, falling back
to `TMPDIR` and then `/tmp`. Eligible-read databases, binary/ACF arrays, package
caches and validation files stay there. The runner removes the workspace on
normal exit, errors, `TERM` and `INT`; cleanup cannot be guaranteed after `SIGKILL`
or node failure.

Rerunning the same `sbatch` command starts the analysis again; there is no cached
resume across submissions. A lock on the output directory itself prevents
concurrent updated runners from writing there without creating a lock file.
If an older `.pipeline.lock` exists, the runner opens it read-only and checks its
lock as well, refusing to proceed while an older runner holds it.

Checks are embedded in the numbered analysis/plot scripts; there is no separate
`07_validate_tss_run.py`. They enforce canonical protein-coding metadata,
expressed-only bins, full spans, unique read IDs, balanced sampling, 2,000-base
binary rows, m6A counts, matrix/metadata order, lag-zero normalization, invalid
row exclusions, successful joint Leiden/UMAP, global labels, heatmap ordering,
composition denominators and complete R plot generation.

Runtime defaults are the existing
`/project/spott/cshan/envs/Jupyter-notebook/bin/python` and
`/software/R-4.4.1-el8-x86_64/bin/Rscript`; override via `TSS_PYTHON` and
`TSS_RSCRIPT`. Python requires NumPy, pandas, pysam, SciPy, Scanpy, anndata,
scikit-learn, umap-learn, igraph and leidenalg. Only expression-table recovery
additionally requires R `data.table`. No environment installation or parent-code
modification is required.
