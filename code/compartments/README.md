# Compartment covariation analysis

Implements `PLAN.md` with unchanged shared parsing and plotting helpers. All
coordinates in tables are BED-style, except read call positions inside the
parser. Bins are 20 kb; chromosome tails shorter than 20 kb are dropped.

## Run

Use the provided R environment wrapper from this directory:

```bash
cd /project/spott/cshan/fiber-seq/code/compartments
bash run_R.sh 01_m6a_per_read.R --list
```

First run chr21–22 in a separate output directory (selected automatically by
`COMPARTMENTS_CHROMS`). The two arrays may run concurrently:

```bash
export COMPARTMENTS_CHROMS=chr21,chr22
bash run_R.sh 01_m6a_per_read.R --list  # 21 chunks
sbatch --array=1-21 01_m6a_per_read.sh
sbatch 02_nucleosome_fraction.sh
```

After both arrays succeed, run stages 03, 04, and 05 in order, preferably in a
SLURM allocation (caslake, pi-spott, 1 CPU, 24 GB, 3 hours). Each step stops on
missing required input. `run_R.sh` accepts any R script and trailing arguments:

```bash
bash run_R.sh 03_bin_tables.R
bash run_R.sh 04_bin_similarity.R
bash run_R.sh 05_plots.R
```

For the full genome, unset `COMPARTMENTS_CHROMS` and repeat the two arrays
using their defaults (587 chunks and 31 samples), then stages 03–05. The output
root is `LCL_project/compartments/`; the pilot uses its `chr21_chr22/` subdirectory.
`COMPARTMENTS_OUT` can explicitly select another output directory, including
scratch. Keep the same environment in all stages. Do not mix chromosome sets in
one output directory. Array task indices can also be passed directly to R for
single-task debugging. No jobs are submitted automatically by these scripts.

## Outputs and conventions

- Per-read chunk tables retain zero-call reads in successfully assembled
  sample-windows, including reads with zero observed denominator. Only the
  parser's `Insufficient ...` errors become empty windows; other errors stop.
- Denominators use the union of observed sites across all successfully assembled
  samples in each bin. A sample with no assembled window contributes no sites or
  reads. Read spans use inclusive `start <= position <= end`.
- Pooled rates use summed numerators and denominators. Distribution tables use
  only reads overlapping at least 5 kb with a positive denominator. Their
  `n_reads` is the number of qualifying reads. Missing sample-bins remain in
  aggregate tables, with zero counts and NA rates/distribution summaries.
- Raw m6A tables include all overlapping reads in `n_reads`. Normalization masks
  fewer than 10 reads, nucleosome depth below 5, nonpositive rates, and missing
  values, then applies log and `limma::normalizeQuantiles` across samples.
- `normalization_quantiles_*.tsv` records empirical percentiles for QC. Missing
  values, unequal sample sizes and tied values can prevent exact equality of
  empirical quantiles after limma's interpolation/tie handling.
- Pair tables retain every lagged pair, including NA correlations and excluded
  TAD classes. Constant profiles yield NA. Summary counts include only finite
  correlations with at least 20 shared samples; means, medians and IQRs are
  exact genome-wide summaries, not averages of chromosome summaries.
- Domain membership is half-open; edges must lie strictly between midpoints.
  Any intervening edge takes precedence over shared domain membership,
  including nested-domain edges. Pairs outside domains are excluded unless
  they cross an edge, as specified by the plan's across-boundary rule.
- Plots use only saved tables, including the saved Rao domains. Distance curves
  require both classes to have at least 1,000 valid pairs and shade lags 1–2.
  Blank distance panels mean no lag passed this count threshold. Correlation
  maps include pairs regardless of TAD class and rasterize the rotated lattice.
  Normalization heatmaps share a colour scale across chromosomes per signal.

## Validation

Validation scripts have moved to `/scratch/midway3/cshan/tmp/compartments-validation-20261105/`.
Use their absolute paths with `run_R.sh`; background logs and results are saved there.

```bash
bash run_R.sh /scratch/midway3/cshan/tmp/compartments-validation-20261105/validate.R              # synthetic numerics and syntax
bash run_R.sh /scratch/midway3/cshan/tmp/compartments-validation-20261105/validate.R --real       # also two real windows and AL8 pileup
bash run_R.sh /scratch/midway3/cshan/tmp/compartments-validation-20261105/validate_plots.R        # all five PDFs from a synthetic fixture
bash run_R.sh /scratch/midway3/cshan/tmp/compartments-validation-20261105/validate_pipeline.R     # synthetic array outputs through stages 03-05
```

Real checks compare direct per-bin parser calls with window-derived per-read
counts and denominators for all 31 samples at chr21:20,000,000–20,100,000 and
21,000,000–21,100,000. They also compare exact bigWig integrals to AL8 pileup
coverage over three bins. The plot fixture writes only to
`/tmp/compartments-plot-validation/`. The pipeline fixture uses
`/tmp/compartments-pipeline-validation/`.

Run validation in a background SLURM job when desired:

```bash
sbatch --account=pi-spott --partition=caslake --cpus-per-task=1 --mem=24G \
  --time=01:00:00 --output=compartments-validation-%j.log \
  run_R.sh /scratch/midway3/cshan/tmp/compartments-validation-20261105/validate.R --real
```

The full chr21–22 pilot remains a separate integration check before running
chr1–22. Inspect its normalization percentiles, boundary class counts, and PDFs;
large lags can legitimately have no within-domain pairs. No significance tests
are performed.
