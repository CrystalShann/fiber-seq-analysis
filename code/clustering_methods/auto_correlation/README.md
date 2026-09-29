# LCL Fiber-seq single-molecule autocorrelation

The scripts are numbered consecutively:

1. `01_run_autocorrelation.py`: prepare reads, calculate ACFs, cluster, and save report inputs.
2. `02_build_m6a_input.r`: extract phased, fully spanning m6A reads for the runner.
3. `03_compute_autocorrelations.py`: calculate per-read ACFs.
4. `04_cluster_autocorrelations.py`: cluster valid ACF profiles and calculate UMAP coordinates.
5. `05_summarize_autocorrelations.py`: calculate cluster and allele summaries.
6. `06_autocorrelation.Rmd`: define the R plotting functions, calculate haplotype repeat candidates, and knit the report.

Run the analysis on the cluster with the original read sources available. The
default command fills missing report inputs in the existing `resolution04` run:

```bash
cd /project/spott/cshan/fiber-seq/code/clustering_methods/auto_correlation
/project/spott/cshan/envs/Jupyter-notebook/bin/python -B 01_run_autocorrelation.py
```

The runner checks R packages before starting. It tries `Rscript` on `PATH`,
then the cluster's R 4.4.1 installation with its matching user library if needed.
It supplies `R_LIBS_USER` explicitly because `--vanilla` skips `.Renviron`.
The selected R executable and version are printed once. An explicit
`AUTOCOR_RSCRIPT=/path/to/Rscript` override is respected and reports missing
packages instead of falling back to another R.

Extracted binary m6A matrices and read metadata are saved together in
`inputs/<region_id>/m6a_input.npz` before clustering. Reruns reuse this compressed
cache when the selection, region, extraction code, and source-file metadata
match. Report preparation can therefore be retried without rereading the BEDs.
The cache is written atomically; it contains NumPy arrays and JSON metadata
and can be opened with `numpy.load(path, allow_pickle=False)`.

```bash
env RSTUDIO_PANDOC=/software/pandoc-2.17.1.1-el8-x86_64/bin \
  R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4" \
  LD_LIBRARY_PATH=/software/gcc-12.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib \
  /software/R-4.4.1-el8-x86_64/bin/Rscript --vanilla -e \
  'rmarkdown::render("06_autocorrelation.Rmd", output_file="autocorrelation.html")'
```
