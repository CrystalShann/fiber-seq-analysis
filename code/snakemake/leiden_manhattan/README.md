# Leiden + Manhattan workflow

Run commands from `code/snakemake/leiden_manhattan`

Submit jobs using `snakemake-executor-plugin-slurm`

Build the LCL sample sheet once:

```bash
module load R/4.4.1
Rscript --vanilla config/make_lcl_sample_sheet.R \
  /project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv \
  /project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples \
  config/samples_lcl.tsv
```

All samples require unique `sample_name`, `sample_label`, and `fire_dir`
  - phasing also needs `phasing_dir` and `cell_line`

Dry-run one dataset, or submit it to SLURM:

```bash
snakemake -n --config datasets_to_run=macrophage
snakemake --profile profiles/slurm --config datasets_to_run=macrophage
```

All settings are in `config/config.yaml`

`parsing_functions` names the shared
`code/parsing_functions/parsing_footprints_functions.r`, relative to `project_root`.
The R scripts source this configured file through `scripts/common.R`. It contains
the BED12/tabix readers, read metadata and methylation matrix helpers, and
`read_sample_region_reads()` for per-sample m6A/CpG loading. The topic-model
assembly and the default Leiden assembly share that loader; the full-span Leiden
branch keeps its existing read selection.

All BED12 block expansion goes through `convert_ft_bed12_to_bed6()`. The shared
`footprint_format_columns()` table in `parsing_footprints_functions.r` maps explicit
format names to column counts; column counts never select the format. The
converter handles sentinel policy, longest-alignment selection, and optional block
validation. Existing strict checks apply only to `bed13_fiberhmm`; fibertools
paths retain their previous acceptance behavior. BED6/BED4 readers keep their
column checks and region filtering. Configured format values are unchanged.

`plotting_functions` names the single shared file
`code/parsing_functions/plotting_functions.r`, relative to `project_root`.
`common.R` sources it after the parsing file and before the Leiden plot library.
It defines palettes, `timepoint_palette()`, `theme_fiberseq()`, `save_figure()`,
`plot_stacked_proportion()`, `plot_interval_track()`, `plot_group_profile()` and
`plot_read_heatmap()`
without loading packages or drawing anything when sourced. Timepoint colors are
LPS_0 `#bdbdbd`, LPS_5 `#6baed6`, LPS_10 `#2171b5`, and LPS_15 `#08306b`;
LCL sample and footprint-class palettes remain separate.

Both files are included in `CORE`, including the plot and summarize rule inputs,
so changes make dependent rules eligible to rerun. Plot filenames, dimensions,
PDF devices, `panels.tsv` statuses and manifests retain their existing behavior.
Raw-read proportions use count normalization; precomputed fractions retain their
original denominators. Profile drawing uses the shared renderer while preparation
keeps each existing denominator, smoothing window and promoter orientation. Ordered
layer specifications preserve ribbons, outlines, columns and annotation placement.
Read heatmaps share one renderer with genomic, clustering-feature and FIRE-window
layouts. Each layout retains its row ordering, annotations and raster settings;
PDF filenames, dimensions and manifest panel identifiers are unchanged.
The macrophage composition figure keeps its stacked first
panel and dodged second panel. Direct heatmap PDF exports inherit `pdf.options()`
for the background; ggplot exports retain their configured background.

Dry-run all configured datasets without executing analysis jobs:

```bash
snakemake -n
```

Or dry-run each dataset with `--config datasets_to_run=macrophage`,
`--config datasets_to_run=lcl_asfire`, or `--config datasets_to_run=lcl_promoters`.

```yaml
clustering:
  window_size: [0, 10]
  k_neighbors: [10, 20]
  resolution: [0.5, 1]
  seed: [1, 2]
  kernel_sigma: null
```

Output path are
`<results_dir>/<dataset>/<region_id>/bin<B>/k<K>/res<R>/seed<S>/{tables,plots}`.
The method remains exact `dist(..., method="manhattan")`


Custom regions can also be supplied as a BED file with `regions.bed`. BED column 4
is used as the unique `region_id`; when it is absent or `.`, IDs are generated as
`chr_start_end` using the workflow's 1-based inclusive coordinates. BED coordinates are converted from 0-based,
half-open to the workflow's 1-based, inclusive analysis coordinates. For example:

```yaml
regions:
  bed: config/my_regions.bed
```

The macrophage configuration uses this adapter for its two configured custom
regions in `config/macrophage_custom_regions.bed`.

To change nucleosome size, edit the dataset's `footprints.nucleosome`:

```yaml
min_size: 101
max_size: null
```
