# m6A vs Pol II Pausing Index

This folder compares m6A levels around TSSs for genes with high and low Pol II pausing index, using the pausing-index outputs in:

`/project/spott/cshan/fiber-seq/results/PolII/annotations`

The default script behavior matches the high/low PI setup in `PolII_FP.ipynb`:

- High PI: genes in the top quartile of PI.
- Low PI: genes in the bottom quartile of PI.
- TSS anchor: CAGE TSS when available, otherwise GENCODE TSS.
- m6A calls: `modbase_code = "a"` from the phased modBAM with `MOD_PROB_THRESHOLD=0.9`.

## Run

```bash
module load R
Rscript /project/spott/cshan/fiber-seq/code/footprintR_modbam_code/m6a_pausing_index/run_m6a_pausing_index.R
```

Useful overrides:

```bash
SAMPLE=AL10_bc2178_19130 \
CHROMS=chr1,chr2 \
WINDOW_BP=10000 \
BIN_SIZE=100 \
TSS_SOURCE_FILTER=ALL \
LOW_Q=0.25 \
HIGH_Q=0.75 \
Rscript /project/spott/cshan/fiber-seq/code/footprintR_modbam_code/m6a_pausing_index/run_m6a_pausing_index.R
```

Outputs are written under:

`/project/spott/cshan/fiber-seq/results/PolII/m6a_pausing_index`

Main outputs:

- `tables/*_high_low_pausing_gene_groups.tsv`
- `tables/*_m6A_tss_aligned_positions.tsv.gz`
- `tables/*_gene_m6A_summary_high_low_PI.tsv`
- `tables/*_gene_bin_m6A_summary_high_low_PI.tsv.gz`
- `tables/*_gene_bin_m6A_modified_call_histogram_high_low_PI.tsv`
- `tables/*_m6A_metaprofile_high_low_PI.tsv`
- `tables/*_m6A_group_summary_high_low_PI.tsv`
- `plots/*_m6A_metaprofile_high_low_PI.pdf`
- `plots/*_gene_m6A_fraction_boxplot_high_low_PI.pdf`
- `plots/*_gene_bin_m6A_modified_call_histogram_high_low_PI.pdf`

## Plot Existing Outputs

If the m6A calls have already been generated and you only want to make the extra boxplot/histogram:

```bash
module load R

OUTPUT_DIR=/project/spott/cshan/fiber-seq/results/PolII/m6a_pausing_index/AL10_bc2178_19130_1kb_bin10_modthresh0.9_all \
SAMPLE=AL10_bc2178_19130 \
WINDOW_BP=1000 \
BIN_SIZE=10 \
Rscript /project/spott/cshan/fiber-seq/code/footprintR_modbam_code/m6a_pausing_index/plot_m6a_pausing_index_outputs.R
```

`MAX_MODIFIED_HIST_BIN` controls the top category in the histogram plot. The default is `10`, so genes with at least 10 modified m6A calls in a bin are grouped as `10+`.
