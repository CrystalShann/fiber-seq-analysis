# parsing_functions

Shared R functions used across the project:

| File | Contents |
|---|---|
| `parsing_footprints_functions.r` | Readers for fibertools / FiberHMM BED files, and `assemble_region_m6a()`, which builds the read x position m6A (or CpG) matrix used for clustering |
| `plotting_functions.r` | Shared colour palettes, ggplot theme, figure saving and plot builders (bar charts, interval tracks, profiles, read heatmaps) |

Every function has a comment block above it in the source file listing its
inputs and output. This README gives the overview, the input files the
functions read, and the output files they write.

## Loading

```r
library(dplyr); library(GenomicRanges); library(Rsamtools); library(Matrix)
source("/project/spott/cshan/fiber-seq/code/parsing_functions/parsing_footprints_functions.r")
source("/project/spott/cshan/fiber-seq/code/parsing_functions/plotting_functions.r")
```

Sourcing only defines functions and palettes; nothing is read, drawn or
written. The parsing functions also use `data.table`, `IRanges`; the plotting
functions use `ggplot2`, `cowplot`, `ComplexHeatmap`, `grid` and `scales`.
`plot_group_profile(coordinates = "relative")` and `plot_read_heatmap()` also
need helpers from `code/clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r`.

Sourced by the Snakemake workflow (`code/snakemake/leiden_manhattan/scripts/common.R`,
configured in `config/config.yaml`), the Leiden notebooks
(`leiden_LCL.Rmd`, `marcophage_leiden_manhattan.Rmd`, `testing_clustering_parameters.Rmd`),
the autocorrelation pipeline (`02_build_m6a_input.r`, `06_autocorrelation.Rmd`),
`fire_frequency/05_vis_fire_freq.Rmd` and `TF_enrichment.ipynb`, the topic model
(`topic_modelling_functions.r`, `read_region_extraction/region_data_utils.R`),
`nucleosome_positioning/process_nuc_pos.R` and `co-accessibility/coaccess_plot_functions.R`.
`plotting_functions.r` is also sourced by the plotting scripts of those analyses.

## Input files

### fibertools `ft extract` BED12 (m6A, CpG, nucleosomes)

One row per read alignment: `chr, start, end, RID, score, strand, thickStart,
thickEnd, rgb, blockCount, blockSizes, blockStarts`. Each block is one call
(1 bp for m6A / CpG) or one nucleosome. The first and last block of every read
are sentinels and are dropped. A call's 1-based position is
`BED start + blockStart + blockSize`. Files are bgzipped and tabix-indexed
(`.tbi` next to each file).

**Macrophage** (samples `LPS_0`, `LPS_5`, `LPS_10`, `LPS_15`; chr1-22, X, Y, M)

```
/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/ft_result_dir/<sample>/extracted_results/<modality>_by_chr/<sample>.ft_extracted_<modality>.<chr>.bed.gz
```
`<modality>` is `m6a`, `cpg` or `nuc`. Example:
`/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/ft_result_dir/LPS_0/extracted_results/m6a_by_chr/LPS_0.ft_extracted_m6a.chr1.bed.gz`

**LCL** (the 31 samples in `code/snakemake/leiden_manhattan/config/samples_lcl.tsv`, e.g. `AL10_bc2178_19130`)

```
/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results/<sample>/extracted_results/<modality>_by_chr/<sample>.ft_extracted_<modality>.<chr>.bed.gz
```
Examples:
`/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results/AL10_bc2178_19130/extracted_results/m6a_by_chr/AL10_bc2178_19130.ft_extracted_m6a.chr1.bed.gz`,
`/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results/AL10_bc2178_19130/extracted_results/nuc_by_chr/AL10_bc2178_19130.ft_extracted_nuc.chr2.bed.gz`
(nucleosomes, read by `extract_nucleosomes()` in `leiden_manhattan_plots.r`).

### FiberHMM footprints

The format names below are the ones accepted by `footprint_format_columns()` and
`convert_ft_bed12_to_bed6(format = ...)`.

| Format | Columns | Files | Read by |
|---|---|---|---|
| `bed13_fiberhmm` | BED12 + per-block scores, no sentinels | `/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_footprint/<LPS0\|LPS5\|LPS10\|LPS15>/<label>_hmm_extracted_footprint_<chr>.bed.gz` (+ `.tbi`); per-region slices `/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/region_plots/<region>/<LPS_x>/parsed/region.fiberhmm_footprint.bed` | Snakemake macrophage nucleosomes (`scripts/footprint_io.R`), `nucleosome_positioning/process_nuc_pos.R`, `co-accessibility/coaccess_plot_functions.R`, `read_region_extraction/region_data_utils.R` |
| `bed15_fiberhmm_tf` | BED12 + per-block tq, left edge, right edge | `/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/region_plots/<region>/<LPS_x>/parsed/region.fiberhmm_tf.bed`, cut from `/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_tf/<label>/<label>_hmm_extracted_tf_<chr>.bed.gz` | `read_region_extraction/region_data_utils.R` (via `05_vis_fire_freq.Rmd`) |
| `bed6_per_sample` | one footprint per row, one file per sample | `/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_tf/ft_by_size/<size10-30\|size40-60\|size60-80>/<label>/<label>_tf_size<size>_<chr>.bed.gz` (+ `.tbi`) | Snakemake macrophage TF footprints (`scripts/footprint_io.R`), `nucleosome_positioning/process_nuc_pos.R` |
| `bed4_pooled` | one footprint per row, samples pooled | `/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged/combined/joint_trained_tracks/<chr>/combined_<chr>_<size>bp_fps.bed.gz` (chr1-22, no `.tbi`, so it is streamed rather than queried with tabix) | Snakemake `lcl_asfire` / `lcl_promoters` TF footprints, `leiden_LCL.Rmd` |

## `parsing_footprints_functions.r`

| Function | What it does | Returns |
|---|---|---|
| `read_tabix_region()` | Queries a tabix-indexed `.bed.gz` (any file above) for the records overlapping a region | data.frame with unnamed columns `V1, V2, ...` |
| `footprint_format_columns()` | Column count expected for a format name (table above) | integer |
| `footprint_column_error()` | Stops with a message naming the file, its format and the wrong column count | (error) |
| `keep_longest_alignment()` | Keeps one alignment per read (the longest) when a read has supplementary alignments | data.frame, one row per RID |
| `read_ft_bed12()` | Reads a fibertools BED12, the whole file or one region with tabix, and names its 12 columns | data.frame, one row per alignment |
| `convert_ft_bed12_to_bed6()` | Expands BED12 blocks into one row per call / footprint, handles sentinels, per-block scores and format checks | data.frame `chr, start, end, RID, score, strand` (+ optional columns) |
| `assemble_region_m6a()` | Builds the pooled read x position matrix for one region across samples (details below) | list `met_mat, rids_df, qc, signature` |
| `extract_ft_region_reads()` | Calls (1-based positions) of the reads overlapping a region, from one sample's BED12 | data.frame, one row per call |
| `extract_ft_read_info()` | One row per read (RID, chr, start, end, strand) from a calls table or straight from a BED12 | data.frame, one row per read |
| `read_sample_region_reads()` | `extract_ft_region_reads()` for each sample under an `ft_result_dir` | list, one calls table per sample |

`extract_ft_region_reads()`, `extract_ft_read_info()` and
`read_sample_region_reads()` are kept for interactive use; no script currently
calls them (matrix building goes through `assemble_region_m6a()`).

### `assemble_region_m6a()`

For each sample it reads the region from `<fire_dir>/extracted_results/<modality>_by_chr/<sample>.ft_extracted_<modality>.<chr>.bed.gz`
(or the paths given in `m6a_paths`), keeps each read's longest alignment and
turns the calls inside the analysis window into a sparse read x position matrix.
Read IDs are `<sample_name>::<read name>`; `rids_df$original_RID` keeps the read name.

| Option | Values |
|---|---|
| `modality` | `"m6a"` (default) or `"cpg"` |
| `full_span` | `TRUE`: only reads whose alignment covers the whole window; entries are 1 (call) / 0 (no call). `FALSE` (default): every overlapping read; NA where the read does not reach |
| `positions` | `NULL` (default): positions with at least one call; or a vector, e.g. every base of the window |
| `matrix_dir`, `reuse` | optional cache folder (see outputs); the cache is reused only when its region, samples, paths and options are identical |

Callers: the Leiden workflows (`leiden_manhattan_cluster_regions()`, `run_lcl_clustering()`,
`leiden_LCL.Rmd`, `marcophage_leiden_manhattan.Rmd`, `05_vis_fire_freq.Rmd`, `TF_enrichment.ipynb`),
the Snakemake step `scripts/assemble_raw_reads.R`, the LCL autocorrelation input
(`code/haplotype_phasing/LCL_phased_m6a_input.r`, every base as columns) and the
topic model (`assemble_region_met_data()` in `code/topic_model/topic_modelling_functions.r`,
the only CpG caller).

## `plotting_functions.r`

| Function / object | What it does | Returns |
|---|---|---|
| `LEIDEN_CLUSTER_COLORS` | 25 cluster colours | character vector |
| `LEIDEN_TIMEPOINT_COLORS` | Colours for `LPS_0`, `LPS_5`, `LPS_10`, `LPS_15` | named character vector |
| `cluster_palette()` | Colours clusters by their position in the level order | named colours |
| `cluster_id_palette()` | Colours clusters by the number in "cluster<N>", so a cluster keeps its colour when others are absent | named colours |
| `timepoint_palette()` | Timepoint colours for labels written `LPS_5`, `5 min` or `5` | named colours |
| `theme_fiberseq()` | Shared ggplot theme (classic, bw or cowplot base) with optional overrides | ggplot theme |
| `save_figure()` | Saves a ggplot, or a grid / ComplexHeatmap drawing through a callback, creating the folder | the path (invisibly); writes the file |
| `plot_stacked_proportion()` | Stacked / filled / dodged bars, e.g. cluster composition per timepoint | ggplot |
| `plot_interval_track()` | Genomic intervals (cCREs, peaks, tested regions) as rectangles or segments | ggplot |
| `plot_group_profile()` | Value along the genome per group (e.g. m6A fraction per cluster) as lines, areas, columns or ribbons | ggplot |
| `plot_read_heatmap()` | Read x position m6A heatmap split by cluster, in `genomic`, `features` or `fire` layout | ComplexHeatmap object (or drawn) |

## Output files

Only `assemble_region_m6a()` (parsing) and `save_figure()` (plotting) write
files; every other function returns its result in memory.

### `assemble_region_m6a()` caches (when `matrix_dir` is set)

```
<matrix_dir>/<region_id>/<region_id>_m6a_matrix.rds     # the returned list
<matrix_dir>/<region_id>/<region_id>_read_qc.tsv        # per sample: overlapping_reads, full_span_reads
```
CpG caches are named `<region_id>_cpg_matrix.rds` and `<region_id>_cpg_read_qc.tsv`
(none exist yet).

Existing: the 10 LCL top AS-FIRE regions written by `leiden_LCL.Rmd`, e.g.
`/project/spott/cshan/fiber-seq/LCL_project/Leiden_manhattan/m6a_summary/rank01_rs12470189_chr2_119290945_119291178/rank01_rs12470189_chr2_119290945_119291178_m6a_matrix.rds`
and `.../rank01_rs12470189_chr2_119290945_119291178_read_qc.tsv`.

The Snakemake step `assemble_raw_reads.R` saves the returned list (plus region
and sample metadata) to
`/project/spott/cshan/fiber-seq/snakemake_results/leiden_man_clustering/<dataset>/cache/raw/<region_id>.rds`.
These exist for the 10 macrophage regions, e.g.
`/project/spott/cshan/fiber-seq/snakemake_results/leiden_man_clustering/macrophage/cache/raw/IL1B.rds`;
`lcl_asfire` and `lcl_promoters` have not been run yet (only `regions.tsv`).

### `save_figure()` outputs (path chosen by the caller)

| Caller | Output | Existing example |
|---|---|---|
| Snakemake `scripts/plot.R` | `/project/spott/cshan/fiber-seq/snakemake_results/leiden_man_clustering/<dataset>/<region_id>/bin<w>/k<k>/res<r>/seed<s>/plots/<panel>.pdf` | `/project/spott/cshan/fiber-seq/snakemake_results/leiden_man_clustering/macrophage/IL1B/bin0/k10/res1/seed1/plots/heatmap_IL1B.pdf` |
| `save_fiberseq_plots()` (`leiden_manhattan_plots.r`, from `leiden_LCL.Rmd`) | `/project/spott/cshan/fiber-seq/LCL_project/Leiden_manhattan/<region_id>/bin0/k10/resolution1/plots/*.pdf` | `/project/spott/cshan/fiber-seq/LCL_project/Leiden_manhattan/rank01_rs12470189_chr2_119290945_119291178/bin0/k10/resolution1/plots/fig1_read_footprints.pdf` |
| `fire_frequency/05_vis_fire_freq.Rmd` | `/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/top10_example_plots/rank<NN>_<region>/rank<NN>_panel_*.pdf` | `/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/top10_example_plots/rank01_chr11_128582973_128583191/rank01_panel_A_aggregate_m6a.pdf` |
| `fire_frequency/06_footprints_in_FIRE.Rmd` | `/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/plots/footprints_in_FIRE/*.pdf` | `/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/plots/footprints_in_FIRE/fire_freq_vs_footprints_groups.pdf` |
| `co-accessibility/coaccess_examples.Rmd` | `/project/spott/cshan/fiber-seq/macrophage_project/co-accessibility/coaccess/plots/<pair_id>_{reads,config_bars}.pdf` | `/project/spott/cshan/fiber-seq/macrophage_project/co-accessibility/coaccess/plots/pair1_TIFAB_chr5_135456305_reads.pdf` |
