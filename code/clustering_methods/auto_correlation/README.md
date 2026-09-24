# LCL Fiber-seq single-molecule autocorrelation

The selected top-ten FIRE-SNP regions are clustered separately from their raw
binary m6A autocorrelations. Only fully spanning, resolved phased heterozygous
reads are retained. Default Leiden resolution is **0.4**, configurable with
`--resolution VALUE`. The separate `tss_test` workflow is unchanged.

## Run

```bash
cd /project/spott/cshan/fiber-seq/code/clustering_methods/auto_correlation
/project/spott/cshan/envs/Jupyter-notebook/bin/python -B 01_run_autocorrelation.py
```

The default run root is
`/project/spott/cshan/fiber-seq/LCL_project/auto_correlation/resolution04`.
The driver refuses to overwrite a nonempty run. For another run, choose a new,
undated directory inside `LCL_project/auto_correlation`:

```bash
/project/spott/cshan/envs/Jupyter-notebook/bin/python -B 01_run_autocorrelation.py \
  --resolution 1.0 \
  --output-dir /project/spott/cshan/fiber-seq/LCL_project/auto_correlation/resolution1_repeat
```

Run **01** once; it calls 02–05, 07 and 09. No notebook is knitted by the driver.

- `--region REGION_ID`: select one exact saved region ID; repeat for several.
  Omit to run all selected regions.
- `--window-size 2000`: SNP-centered window width, at least three bases.
- `--n-features N`: retain the first N ACF lags including lag zero. Omit for all
  window-size lags (0–1999 for the default 2-kb window).
- `--n-pcs 50`, `--n-neighbors 10`, `--resolution 0.4`, `--seed 0`: clustering parameters.
- `--selection-signature PATH`: existing RDS containing selected regions and
  sample metadata. No accompanying TSV is required.
- `--sample-table PATH`: source CSV mapping sample identifiers to cell lines.

The original selection RDS is used if present. Otherwise the default is the
preserved `resolution1/inputs/selection_signature.rds`. This input was recovered
from the original Leiden regional RDS files and matched against the previous
FIRE-SNP selection; it remains necessary for reruns. New runs read it directly
and do not write another selection RDS or copy input tables.

## Output policy and layout

```text
resolution04/
  inputs/
    rank01_<SNP>_<chr>_<start>_<end>/
      extraction.log
  outputs/
    rank01_<SNP>_<chr>_<start>_<end>/
      allele_average_acf.pdf
      ... other final PDFs ...
      plotting.log
```

All regions use this layout. PDFs are saved directly in their region folder;
no nested `plots/` folder remains. Existing source input files and diagnostic
logs are preserved. Future runs produce **PDFs and diagnostic logs only**:
no JSON manifests, TSV/CSV result tables, assignments, summaries, RDS caches,
or methylation/ACF/PCA/graph matrices are written. Region metadata, read audits,
assignments, allele summaries and occupancy tables stay in memory. Python
passes plotting metadata and arrays to R through standard input; the serialized
pipe payload is never saved as a JSON file. R returns its completion count
through standard output.

There is no saved-assignment refresh mode because assignments are no longer
persisted. To regenerate plots, rerun the selected region(s) in a new output
directory. The optional notebook now discovers PDFs directly, without reading
deleted JSON/TSV files. Notebook parameters apply only when launching a new
analysis; existing plots are not relabeled with those parameters.

## Input and method

`02_build_m6a_input.r` reuses `lcl_build_phased_m6a_input()` and the original
BED extraction and phasing functions. Every genomic base is retained, including
non-A positions and zero-call columns. A covered call is 1 and a covered non-call
is 0; only completely covered windows are retained. The existing extractor
omits spanning molecules with no in-window calls. The workflow additionally
requires HP1/HP2 and a resolved heterozygous focal genotype (`0|1` or `1|0`).
No smoothing, binning, interpolation or sample balancing is applied.

For window length N, the linear per-read ACF is
`sum((x[t]-mean(x)) * (x[t+k]-mean(x))) / (N * population_variance(x))`.
All retained lags enter centered, unscaled PCA, followed by correlation-distance
kNN and weighted Leiden (`leidenalg`, directed, fixed seed). UMAP is for display.
Constant signals are excluded; small or degenerate datasets are reported without
artificial clusters. Cluster numbers are local to each region.

Pooled REF/ALT plots use unique molecules mapped through their phased focal
genotypes. Cluster fractions use clustered reads; ACF means and peak fractions
use valid ACF reads. These are descriptive molecule summaries, not independent
biological-replicate tests. The 140–250-bp positive-local-peak feature is not a
validated nucleosome repeat-length estimate.

`allele_average_acf.pdf` contains three views of the same pooled REF/ALT curves:
all retained lags, 0–500 bp, and 500–2000 bp. Vertical scales are independent.
A 2-kb window supplies lags only through 1999; no value at lag 2000 is invented.
The other allele PDFs show coverage, cluster proportions and peak features.
The region PDFs also include cluster ACFs, heatmaps, UMAP, kNN graphs,
occupancy and annotated Fiber-seq views.

## Haplotype repeat-length plots

`09_haplotype_repeat_plots.py` adds four PDF plots per region:

- `haplotype_median_acf.pdf`: separate REF/ALT median raw ACF curves, with all
  retained lags and a 100–500-bp zoom. A star marks the primary 120–250-bp peak;
  circles mark qualifying observed peaks near integer multiples of that lag.
- `haplotype_fiber_acf_heatmap.pdf`: every unique valid fiber, sorted first by
  phased REF/ALT and then primary peak lag; no-peak fibers sort last. Both full
  positive lags and the first 500 bp are shown. Lag 0 is omitted from the heatmap;
  symmetric color limits use the 99th percentile of absolute positive-lag ACF.
- `locus_haplotype_median_acf.pdf`: two rows (locus × REF/ALT), with candidate
  lags 120–250 bp as columns. Median ACF is color-coded and the primary peak is
  marked. A combined version in `outputs/` compares all selected loci with a
  common color scale. Gray means unavailable, not zero.
- `haplotype_repeat_length_distribution.pdf`: REF/ALT per-fiber primary peak
  lags. Boxes summarize peak-bearing fibers; histogram bins are normalized by
  all valid fibers. Labels give peak-bearing/valid counts.

Every unique retained fiber with a valid ACF contributes, including unclustered
fibers. Physical-molecule deduplication and genotype-to-REF/ALT mapping reuse the
existing allele layer; HP1/HP2 tags are sample-local and are never pooled directly.
Undefined constant ACFs are excluded and counted. Median ACF peaks and medians of
per-fiber peak lags are different summaries and are not substituted for each other.

Peak detection uses `scipy.signal.find_peaks` on the **unsmoothed raw ACF**, with
prominence ≥0.01 ACF units and minimum distance 20 bp. The strongest positive
qualifying peak in 120–250 bp is the candidate repeat length (ties choose the
smaller lag). At least lags 0–251 must be available to evaluate the complete
band. Higher-order peaks require a positive qualifying observed maximum within
±20 bp of each integer multiple, up to the available lag limit. No qualifying
peak means no estimate; peaks are never forced or extrapolated. These fixed
thresholds are display criteria, not a calibrated reliability test. The estimated
chromatin/nucleosome repeat length is a descriptive ACF-based candidate.

These features are annotation-only; none enter PCA, kNN or Leiden. Existing
mean-ACF and 140–250-bp legacy peak plots retain their definitions. The new
median/120–250-bp views are additional plots. All outputs remain PDFs and logs;
there are no JSON/TSV files or persisted per-fiber matrices.

## Validation

```bash
/project/spott/cshan/envs/Jupyter-notebook/bin/python -B 06_test_allele_summaries.py
/project/spott/cshan/envs/Jupyter-notebook/bin/python -B 10_test_haplotype_repeats.py
```

This validates pooled allele summaries and the three-panel ACF figures using
temporary outputs. It does not launch a full FIRE-SNP or all-TSS analysis.
