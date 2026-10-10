# parsing_functions



| File | Contents |
|---|---|
| `parsing_footprints_functions.r` | Readers for fibertools / FiberHMM BED files, and `assemble_region_m6a()`, which builds the read x position m6A (or CpG) matrix used for clustering; per-region extraction and loading for read-level region plots; fibers, display tracks and 2x2 tables for cCRE-pair co-accessibility |
| `plotting_functions.r` | Shared colour palettes, ggplot theme, figure saving and plot builders (bar charts, interval tracks, profiles, stacked read-level region panels, read heatmaps); one-call region plots and cCRE-pair co-accessibility figures |
| `extract_region_result_macrophage.sh` | Cuts one region for one macrophage timepoint into `<out_root>/<outname>/<sample>/parsed/` |
| `extract_region_result_lcl.sh` | The same for one LCL sample, from tabix slices only |

Every function has a comment block above it in the source file listing its
inputs and output. This README gives the overview, the input files the
functions read, and the output files they write.

## Loading

```r
library(dplyr); library(GenomicRanges); library(Rsamtools); library(Matrix)
source("/project/spott/cshan/fiber-seq/code/parsing_functions/parsing_footprints_functions.r")
source("/project/spott/cshan/fiber-seq/code/parsing_functions/plotting_functions.r")
```

Source `parsing_footprints_functions.r` first: `plot_region_example()` and
`plot_coaccess_pair()` in `plotting_functions.r`

`plot_group_profile(coordinates = "relative")` and `plot_read_heatmap()` also
need helpers from `code/clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r`.



Sourced by the Snakemake workflow (`code/snakemake/leiden_manhattan/scripts/common.R`,
configured in `config/config.yaml`), the Leiden notebooks
(`leiden_LCL.Rmd`, `marcophage_leiden_manhattan.Rmd`, `testing_clustering_parameters.Rmd`),
the autocorrelation pipeline (`02_build_m6a_input.r`, `06_autocorrelation.Rmd`),
`fire_frequency/05_vis_fire_freq.Rmd`, `fire_frequency/06_footprints_in_FIRE.Rmd` and
`TF_enrichment.ipynb`, the topic model (`topic_modelling_functions.r`,
`topic_model/1_plot_promoter_regions.Rmd`), `nucleosome_positioning/process_nuc_pos.R`
and `co-accessibility/macrophage/coaccess_macrophage.Rmd`.
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
(nucleosomes, read by `extract_nucleosomes()` in `leiden_manhattan_plots.r`, by
`load_ft_tracks(nuc_source = "ft")`, and, sliced by `extract_region_result_lcl.sh`,
by `read_ft_nuc_region()`).

The m6A files of both datasets are also read by `load_ft_tracks()` (the m6A
raster of the co-accessibility figures) and, sliced into `parsed/`, by
`read_ft_mod_region()`.

### FIRE outputs

| File | Columns | Files | Read by |
|---|---|---|---|
| FIRE elements | per-read FIRE elements (FDR <= 0.05) | `<FIRE root>/<sample>/additional-outputs-v0.1/fire-peaks/<sample>-v0.1-fire-elements.bed.gz` | `load_region()` (the accessibility call of the co-accessibility figures); sliced to `parsed/fire_elements.bed` by both scripts |
| FIRE peaks | 29 (v0.1) | `<FIRE root>/<sample>/<sample>-fire-v0.1-peaks.bed.gz` | sliced to `parsed/fire_peaks.bed` by both scripts, read by `read_fire_peaks_region()` |
| `fire_all` (`ft fire --extract --all`) | 11 (`FIRE_BED_COLS`): per-read segments coloured nucleosome (169,169,169), linker (147,112,219) or FIRE, with the HP tag | LCL only: `/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results/<sample>/extracted_results/<sample>.fire_all.bed.gz` | `load_region(span_source = "fire_all")`; sliced to `parsed/fire.bed` by `extract_region_result_lcl.sh` |
| read spans | one primary alignment per row | macrophage only: `/project/spott/cshan/fiber-seq/macrophage_project/co-accessibility/<sample>/<sample>.read_spans.bed.gz` (`co-accessibility/macrophage/02_read_spans.sh`) | `load_region(span_source = "read_spans")` |
| aligned blocks | primary BED12, split at CIGAR D/N; no sentinels | LCL: `LCL_project/co-accessibility/<sample>/<sample>.aligned_blocks.bed.gz` | `load_region(span_source = "aligned_blocks")` |

`<FIRE root>` is `/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE`
(macrophage) or `/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results` (LCL).
Example: `/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE/LPS_0/additional-outputs-v0.1/fire-peaks/LPS_0-v0.1-fire-elements.bed.gz`.
The macrophage `parsed/fire.bed` is not a slice: `extract_region_result_macrophage.sh`
recomputes it from the region BAM (`ft add-nucleosomes | ft fire | ft fire --extract`),
and `read_fire_region()` reads it.

### FiberHMM footprints

The format names below are the ones accepted by `footprint_format_columns()` and
`convert_ft_bed12_to_bed6(format = ...)`.

| Format | Columns | Files | Read by |
|---|---|---|---|
| `bed13_fiberhmm` | BED12 + per-block scores, no sentinels | `/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_footprint/<LPS0\|LPS5\|LPS10\|LPS15>/<label>_hmm_extracted_footprint_<chr>.bed.gz` (+ `.tbi`); per-region slices `/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/region_plots/<region>/<LPS_x>/parsed/region.fiberhmm_footprint.bed` | Snakemake macrophage nucleosomes (`scripts/footprint_io.R`), `nucleosome_positioning/process_nuc_pos.R`, `load_ft_tracks(nuc_source = "fiberhmm")`, `read_fiberhmm_region()` (slices, via `load_region_results(nucleosome_source = "fiberhmm")`) |
| `bed15_fiberhmm_tf` | BED12 + per-block tq, left edge, right edge | macrophage: `/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/region_plots/<region>/<LPS_x>/parsed/region.fiberhmm_tf.bed`, cut from `/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_tf/<label>/<label>_hmm_extracted_tf_<chr>.bed.gz`; LCL (FiberHMM v2): `/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/FiberHMM_v2/results/results_v1_trained_model/<sample>/<sample>.recalled_tf.bed.gz` (+ `.tbi`) | `read_fiberhmm_region()` (slices, via `load_region_results()`), `load_ft_tracks(tf_source = "recalled_tf")` (LCL, score >= 50) |
| `bed6_per_sample` | one footprint per row, one file per sample | `/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_tf/ft_by_size/<size10-30\|size40-60\|size60-80>/<label>/<label>_tf_size<size>_<chr>.bed.gz` (+ `.tbi`) | Snakemake macrophage TF footprints (`scripts/footprint_io.R`), `nucleosome_positioning/process_nuc_pos.R`, `load_ft_tracks(tf_source = "by_size")` |
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

### Read-level region plots: extraction and loading

These read the `parsed/` folders written by the two extraction scripts (see
[Extraction scripts](#extraction-scripts)). Figures are drawn by
`plot_region_panels()` / `plot_region_example()` in `plotting_functions.r`.

| Function / object | What it does | Returns |
|---|---|---|
| `FT_BED12_COLS`, `FIRE_BED_COLS`, `RGB_NUCLEOSOME`, `RGB_LINKER` | Column names of an `ft extract` BED12 and an `ft fire --extract` BED; FIRE segment colours | character |
| `read_bed()` | Reads a headerless BED-like file, naming its leading columns | data.frame (empty for a 0-byte file) |
| `read_ft_mod_region()` | `parsed/extracted.{m6a,cpg}.bed.gz` -> one row per modified base inside the window (longest alignment, sentinels dropped) | data.frame `chr, RID, strand, read_start, read_end, pos` |
| `read_fire_region()` | `parsed/fire.bed` -> segments labelled nucleosome / linker / FIRE by colour; FIRE scores capped at 1 | data.frame `FIRE_BED_COLS` + `class` |
| `read_fire_peaks_region()` | `parsed/fire_peaks.bed` (29 columns) -> peak interval and FDR | data.frame `chr, start, end, FDR, logFDR`, or NULL |
| `assign_size_class()` | Bins footprint sizes into closed `[lo, hi]` classes; NA outside the breaks | factor `"<lo>-<hi> bp"` |
| `read_fiberhmm_region()` | `parsed/region.fiberhmm_{tf,footprint}.bed` -> one row per footprint, score-filtered and size-classed | data.frame `chr, start, end, RID, size, score, read_start, read_end, class` |
| `read_ft_nuc_region()` | `parsed/extracted.nuc.bed.gz` (LCL) -> one row per nucleosome, default window 130-160 bp | same columns as `read_fiberhmm_region()` (score NA) |
| `subset_footprints_in_fire()` | Keeps footprints lying entirely inside a FIRE element of the same read | subset of the footprints |
| `build_rids_df()` | One row per read from its FIRE segments: span, strand, HP tag, arrow direction | data.frame |
| `build_met_mat()` | Read x position matrix: 1 modified, 0 covered but unmodified, NA not covered | matrix |
| `smooth_mean_slidewindow()` | Sliding-window mean over positions | numeric |
| `pileup_reads()` | Modification fraction per position, optionally per read group (denominator = reads in the group) | data.frame `pos, base, group, cov, met, frac, smooth_cov, smooth_frac` |
| `calculate_dist()` | Euclidean distance over the positions where both reads are non-NA | dist |
| `hclust_reads_by_footprints()` | Clusters reads on their size-classed footprints inside the region | hclust or NULL |
| `order_reads_by_footprints()` | Plot order: reads without footprints first, then hclust order | character vector of RIDs |
| `extract_region_results()` | Runs an extraction script once per sample (skipped when `parsed/fire.bed` exists for the same region, unless `regenerate = TRUE`) | character vector of `parsed/` folders |
| `load_region_results()` | Loads one region across samples from the `parsed/` folders; `nucleosome_source = "fiberhmm"` (macrophage) or `"ft"` (LCL) | list `region, sample_names, reads, rids_df, fire, fire_peaks, fps, fps_infire, nucs, pileup, ...` |
| `add_read_groups()` | Attaches a per-read grouping (timepoint, cluster, ...) and recomputes the pileup per group | the result with `group_col` set |
| `add_topic_clusters()` | `add_read_groups()` from a topic-model `read_topic_assignments_<outname>.tsv` | the grouped result |
| `add_haplotype_groups()` | `add_read_groups()` by the HP tag in `fire.bed` | the grouped result |

Callers: `05_vis_fire_freq.Rmd` and `06_footprints_in_FIRE.Rmd` (through
`plot_region_example()`; 05 also calls `load_region_results()` and
`add_read_groups()` for its top-10 example panels), `TF_enrichment.ipynb` (with
`plot_tf_cluster_profiles.R`) and `topic_model/1_plot_promoter_regions.Rmd`
(`add_topic_clusters()`).

### cCRE-pair co-accessibility

Shared by the macrophage timecourse and the LCL samples; every path is an
argument. The accessibility call is always the FIRE elements; m6A, nucleosome
and TF footprint tracks are display only. Fibers are keyed `"<sample> <read name>"`.

| Function / object | What it does | Returns |
|---|---|---|
| `FP_SIZE_BINS` | TF footprint size bins `size10-30`, `size40-60`, `size60-80` (half-open `[lo, hi)`) | character |
| `CONFIG_LEVELS` | Fiber configurations: both accessible, CRE1 only, CRE2 only, neither | character |
| `CCRE_CLASSES` | ENCODE SCREEN v4 cCRE classes: `PLS`, `pELS`, `dELS`, `CA-H3K4me3`, `CA-CTCF`, `CA-TF`, `CA`, `TF` | character |
| `tabix_region()` | Tabix query of a BED region with the command line tool (1-based query) | data.table |
| `load_region()` | Fiber spans and FIRE elements; `span_source = "aligned_blocks"` for LCL co-accessibility, `"read_spans"` for macrophage, or legacy `"fire_all"` for HP displays | list `region, samples, spans, elements, blocks`; aligned mode clips elements while preserving their IDs |
| `load_ft_tracks()` | Adds the display tracks: `ft extract` m6A; nucleosomes (`nuc_source = "fiberhmm"` or `"ft"`); TF footprints (`tf_source = "by_size"`, `"recalled_tf"` or `"none"`) | the result with `m6a, nuc, size_fps` |
| `label_reads()` | Per-fiber configuration at the pair: accessible when one of its FIRE elements overlaps the cCRE by >= 1 bp; `read_rule = "any"` or `"contain"` | data.table `key, sample_name, shared, acc1, acc2, config` |
| `order_reads()` | Raster order: sample, configuration, start | character vector of keys |
| `coaccess_m6a_fraction()` | Unbinned m6A fraction per position and sample over covering fibers | data.table `pos, sample_name, frac` |
| `coaccess_2x2()` | The pair's 2x2 from the fibers covering both cCREs, with `fisher.test(table + 1)` | one-row data.frame `both, cre1_only, cre2_only, neither, n_shared, fisher_or, fisher_p` |
| `cre_class_from_id()` | cCRE class from a `CRE_ID` (`accession1.accession2.class`) | character |
| `select_cre_pairs()` | Keeps the pairs whose classes match `cre_types` (NULL = any); `cre_match = "either"`, `"both"` or `"pair"` (with `cre_types2`); pairs are unordered | rows of the pair table |

Caller: `co-accessibility/macrophage/coaccess_macrophage.Rmd` (through `plot_coaccess_pair()`
and `select_cre_pairs()`).

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
| `plot_region_panels()` | Stacked read-level region figure (pileup, peaks, FIRE peaks, m6A reads, FIRE / nucleosome / TF footprint reads) from a `load_region_results()`-style list | list `panels, combined` (patchwork) |
| `plot_read_heatmap()` | Read x position m6A heatmap split by cluster, in `genomic`, `features` or `fire` layout | ComplexHeatmap object (or drawn) |
| `M6A_COL`, `BACKBONE_COL`, `CONFIG_COLS` | m6A tile, fiber backbone and configuration colours of the co-accessibility figures | character |
| `coaccess_region_result()` | Turns a cCRE pair's fibers into the `plot_region_panels()` input | list |
| `configuration_proportion_inputs()` | Configuration proportions per bar | list `rows, totals` |
| `plot_coaccess_pair()` | Read-level figure of one cCRE pair plus configuration bars, for either dataset (see below) | list `res, labels, counts, panels, bars, region, height` |
| `plot_region_example()` | One-call region plot: `extract_region_results()` (when `script` is given), `load_region_results()`, optional `group` function, `plot_region_panels()` | list `res, plot, result_dirs` |

### `plot_coaccess_pair()`

| Dataset | `span_source` | `nuc_source` | `tf_source` | `facet_by` |
|---|---|---|---|---|
| Macrophage | `"read_spans"` | `"fiberhmm"` | `"by_size"` | `"sample"` (one facet per timepoint) |
| LCL co-accessibility | `"aligned_blocks"` | `"ft"` | `"recalled_tf"` | `"pooled"` or `"sample"` |
| Legacy LCL HP display | `"fire_all"` | `"ft"` | `"recalled_tf"` | `"hp"` (H1 / H2 / UNK); not the co-accessibility denominator |


## Extraction scripts

Both take `<sample_name> <chr:start-end> <outname> <out_root>` and write
`<out_root>/<outname>/<sample_name>/parsed/`. They are normally run from R by
`extract_region_results()` / `plot_region_example()`, which pass the window as
`chr:start-end` and log the script output; they can also be run directly with bash.

| `parsed/` file | `extract_region_result_macrophage.sh` | `extract_region_result_lcl.sh` |
|---|---|---|
| `extracted.m6a.bed.gz`, `extracted.cpg.bed.gz` (+ `.tbi`) | tabix slice of the macrophage `ft extract` files | tabix slice of the LCL `ft extract` files |
| `extracted.nuc.bed.gz` (+ `.tbi`) | - | tabix slice of `nuc_by_chr` |
| `region.fiberhmm_tf.bed` | slice of `firehmm_tf/<label>/..._tf_<chr>.bed.gz` | slice of the FiberHMM v2 `<sample>.recalled_tf.bed.gz` |
| `region.fiberhmm_footprint.bed` | slice of `firehmm_footprint/<label>/..._footprint_<chr>.bed.gz` (nucleosomes) | - |
| `fire_elements.bed`, `fire_peaks.bed` | slices of the FIRE elements and peaks (`fire_peaks.bed` may be empty) | same |
| `fire.bed` | `ft add-nucleosomes \| ft fire \| ft fire --extract` over `<sample>/region.bam` (samtools slice of the merged BAM) | slice of `<sample>.fire_all.bed.gz` |

Both scripts stop when a required output is empty or has the wrong column
count. The macrophage script needs samtools, ft (`/project/spott/cshan/envs/fire-env/bin/ft`),
tabix and bgzip; the LCL script reads no BAM and runs no ft command.

`extract_region_results()` also writes, next to `parsed/`:

```
<out_root>/<outname>/<sample>/region.txt     # chr:start-end of the last extraction
<out_root>/<outname>/<sample>/extract.log    # script stdout and stderr
```

## Output files

Only `assemble_region_m6a()` (parsing), `extract_region_results()` (parsing,
through the extraction scripts) and `save_figure()` (plotting) write files;
every other function returns its result in memory.

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


### Region extractions (`extract_region_results()`)

| Caller | `out_root` | Existing example |
|---|---|---|
| `fire_frequency/05_vis_fire_freq.Rmd`, `06_footprints_in_FIRE.Rmd` | `/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/region_plots` | `/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/region_plots/chr11_128582973_128583191/LPS_0/parsed/fire.bed` |
| `topic_model/1_plot_promoter_regions.Rmd` | `/project/spott/cshan/fiber-seq/macrophage_project/topic_model/region_plots/early_repsonse_genes/promoters` | none yet |

The existing `region_plots/` folders were extracted before `extract_region_results()`
existed, so they have no `region.txt` / `extract.log`; they are reused as long as
`parsed/fire.bed` is there.

`TF_enrichment.ipynb` runs `extract_region_result_macrophage.sh` from Python
instead, into `/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/tf_motif/0_vs_15_min/example_regions/<TF>_<region>/<sample>/`,
with its own stamp `extraction_complete.json` (region, sample and the script's
modification time).


The pooled all31 LCL notebook is `code/co-accessibility/LCL/LCL_co-access.Rmd`.
Its aligned blocks are produced by `code/co-accessibility/LCL/LCL_read_spans.sh`.
The LCL notebook contains separate chunks that call `load_region()`,
`label_reads()`, `coaccess_2x2()`, and `load_ft_tracks()` directly, then use
`coaccess_region_result()`, `plot_region_panels()`,
`configuration_proportion_inputs()`, and `plot_stacked_proportion()` to draw the
figures. The macrophage notebook uses the `plot_coaccess_pair()` wrapper.
Both notebooks explicitly pass `fire_overlap_fraction = 0.5` to their respective
entry points: one FIRE element on that fiber must overlap at least 50% of the
cCRE length. The shared helper default remains any positive overlap for
backwards compatibility with other callers.
In LCL aligned-block mode, only aligned bases of the same FIRE element count
toward that threshold. D/N gaps are unobserved. `clip_to_aligned_blocks()` also
breaks display annotations at those gaps, and the m6A profile denominator counts
only fibers aligned at each position. Macrophage behavior is unchanged.
