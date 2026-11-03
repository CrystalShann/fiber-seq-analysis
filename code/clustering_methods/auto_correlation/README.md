# LCL Fiber-seq single-molecule autocorrelation in R

Files follow the analysis steps:

1. `01_prepare_m6a_input.R`: select SNP-centered windows and extract fully spanning,
   phased heterozygous m6A reads using the shared LCL input builder.
2. `02_compute_autocorrelations.R`: calculate mean-centered per-read ACFs.
3. `03_cluster_autocorrelations.R`: centered PCA, exact correlation-distance
   neighbors, a uwot fuzzy graph, igraph Leiden clustering, and display UMAP.
4. `04_summarize_autocorrelations.R`: cluster means/SD/SEM, sample composition,
   repeat peaks, REF/ALT annotation, allele summaries and duplicate-molecule audit.
5. `05_run_autocorrelation.R`: run steps 1–4, cache inputs, collect nucleosome and
   FiberHMM footprints, and save portable `report_data.rds` files.
6. `06_autocorrelation.Rmd`: load saved results and produce figures and the report.

All workflow code is R. No Python or retained test scripts are needed.

## Run

Use an R environment with `dplyr`, `Matrix`, `GenomicRanges`, `Rsamtools`,
`IRanges`, `data.table`, `igraph`, `irlba`, `uwot`, and the existing shared
Fiber-seq plotting dependencies installed. Run from the project root:

```bash
Rscript code/clustering_methods/auto_correlation/05_run_autocorrelation.R
```

For one region:

```bash
Rscript code/clustering_methods/auto_correlation/05_run_autocorrelation.R \
  --region rank01_rs12470189_chr2_119290945_119291178
```

Defaults: 2,000-bp windows, all 2,000 lags, 50 PCs, 10 neighbors (including self,
following the previous graph convention), resolution 0.4, seed 0. Options include
`--window-size`, `--n-features`, `--n-pcs`, `--n-neighbors`, `--resolution`,
`--seed`, repeated `--region`, `--selection-signature`, `--sample-table`,
`--output-dir`, `--project`, `--ft-root`, and `--tracks-root`; see `--help`.

The default output remains `LCL_project/auto_correlation/resolution04`. Existing
region outputs are replaced with native R results when their provenance differs
or they came from the former runner. Matching native R results are skipped.
Changed input sources rebuild the corresponding R input cache. Files are replaced
after successful computation; no backup run is created. `--region` limits which
regions are recomputed, while the default processes every selected region.

The ACF calculation reproduces the previous normalization. On the first locus,
the final R implementation also reproduced all four previous clusters exactly;
UMAP coordinates differed. Native R PCA/graph/Leiden can change partitions at
other loci. Each new report bundle records the method and package versions.

The selection signature is reused from the requested output or the existing
Leiden/ACF runs, or recovered from original saved Leiden input signatures.
Binary m6A matrices and read metadata are cached in
`inputs/<region_id>/m6a_input.rds`. Existing NPZ caches are left untouched;
new R runs extract from the original indexed BED and phasing sources. A failed
report-preparation step can resume from the R input cache. Completed matching
report bundles are skipped on rerun.

Knit `06_autocorrelation.Rmd` after preparation. Its default `output_dir` uses
the same `resolution04` path. Figures and HTML are refreshed by knitting, using
the existing R plotting code; the runner itself updates numerical report inputs.

## Computation and compatibility

```r
source("code/clustering_methods/auto_correlation/02_compute_autocorrelations.R")
result <- autocorrelations(binary_calls)
# result$profiles: reads x lags (starting at lag zero)
# result$valid: nonconstant reads
```

The input is a fully covered binary matrix with one column per consecutive
base. `stats::acf()` uses mean-centering and the fixed sum-of-squares denominator
without lag-dependent overlap correction. Constant rows return all `NA` and
`valid = FALSE`. Optional `positions` and `window` select inclusive bounds;
`n_features` selects the number of lags returned.

PCA retains lag zero, centers features without variance-scaling, and uses a
fixed loading-sign convention. Neighbors use correlation of PCA score rows;
UMAP does not determine cluster membership. Leiden optimizes weighted modularity
on the symmetric fuzzy graph. R cluster labels are zero-based, ordered by size;
the report displays them starting at one. Results record the clustering method
and package versions. Invalid or insufficient profiles retain explicit statuses.

The Pol2 report directly sources `02_compute_autocorrelations.R` and
`04_summarize_autocorrelations.R`. Per-molecule arrays are saved as
`autocorrelation_per_molecule_{upstream,downstream}.rds`; summary TSV names and
columns are unchanged. Re-run the Pol2 `--stage autocorrelation` to compute new
RDS arrays. Existing NPZ files are preserved, and knitting can still use existing
summary TSVs.

Methylation display uses calls divided by the number of pooled observed
m6A-call positions (`n_observed_sites`), with NA when none are observed.
Allele summaries validate HP/genotype-to-REF/ALT mapping, deduplicate physical
molecules only when their annotations and ACFs agree, and weight retained
molecules equally.

## Numerical comparison already completed

The R ACF was compared against NumPy on the report's three-read example and 32
cached LCL reads: maximum absolute differences were `2.22e-16` and `2.42e-14`
(tolerance `1e-12`). The Pol2 computation chunk, truncated lags, both windows,
zero-variance rows, cluster means/SD/SEM, sample composition and plateau/tie peak
behavior also passed comparison. The comparison scripts were removed as requested.
These checks establish ACF/summary equivalence, not identical R/Scanpy clustering.
