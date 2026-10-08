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

