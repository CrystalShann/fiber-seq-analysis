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

Extracted binary m6A matrices and read metadata are saved together in
`inputs/<region_id>/m6a_input.npz` before clustering