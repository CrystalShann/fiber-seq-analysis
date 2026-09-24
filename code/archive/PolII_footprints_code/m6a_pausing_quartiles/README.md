# m6A Pol II pausing quartiles

This folder contains the workflow that summarizes m6A signal around transcription start sites after splitting genes into four Pol II pausing index (PI) quartiles. The runner reads the pausing annotation table, assigns Q1-Q4 groups, queries modBAM signal around each TSS, and writes per-gene, per-bin, per-read, and quartile-level summary tables plus plots.

## Workflow

The main entry point is [run_m6a_pausing_quartiles.R](run_m6a_pausing_quartiles.R). It directly implements the quartile workflow and also sources the shared utility file [modbam_footprintR_functions.R](../../footprintR_modbam_code/modbam_footprintR_functions.R), which provides the modBAM readers, long-table conversion helpers, and shared plot/table utilities used elsewhere in the project.

The runner performs the following steps:

1. Reads the pausing index table and filters to standard chromosomes.
2. Chooses the TSS anchor from CAGE or GENCODE coordinates, depending on `TSS_SOURCE_FILTER`.
3. Splits genes into PI quartiles using the 25th, 50th, and 75th percentiles.
4. Builds promoter windows around each TSS.
5. Reads modBAM calls at summary level for the full promoter window and at quickread level for the TSS accessibility window.
6. Aligns both position-level and read-level calls back to the TSS coordinates.
7. Aggregates the aligned calls into the output tables and plots.

## Inputs

- Pausing index table in TSV format with at least:
  - `PI`, `gene_id`, `gene_name`, `chrom`, `strand`, `tss`
  - optional `cage_tss`, `tss_gencode`
- modBAM file with m6A calls and its `.bai` index

## Output Directory

All outputs are written beneath:

`${OUTPUT_ROOT}/${SAMPLE}_${WINDOW_BP/1000}kb_bin${BIN_SIZE}_modthresh${MOD_PROB_THRESHOLD}_${TSS_SOURCE_FILTER}/`

The directory contains:

- `tables/`
- `plots/`

## Output Tables

The table names below match the exact files written by the runner.

### 1. Gene quartile assignment table

File: [tables/${SAMPLE}_quartile_gene_groups.tsv](tables/${SAMPLE}_quartile_gene_groups.tsv)

Columns:

- `gene_name`
- `gene_id`
- `gene_id_base`
- `chrom`
- `strand`
- `tss_anchor`
- `tss_source`
- `cage_tss`
- `tss_gencode`
- `PI`
- `pausing_group`
- `q1_cut`
- `q2_cut`
- `q3_cut`

This table is the gene-level input list after filtering and quartile assignment.

### 2. TSS-aligned m6A positions

File: [tables/${SAMPLE}_m6A_tss_aligned_positions.tsv.gz](tables/${SAMPLE}_m6A_tss_aligned_positions.tsv.gz)

Columns written by the alignment step:

- `call_chromosome`
- `position`
- `call_strand`
- `sample_name`
- `modbase_code`
- `total_calls`
- `modified_calls`
- `unmodified_calls`
- `frac_modified`
- `tss_uid`
- `chromosome`
- `tss_coordinate`
- `tss_strand`
- `gene_name_or_tss_id`
- `tss_id`
- `gene_id`
- `gene_id_base`
- `gene_name`
- `PI`
- `pausing_group`
- `tss_source`
- `promoter_start`
- `promoter_end`
- `rel_pos`

This is the summary-level per-position m6A table aligned to each TSS promoter window, used for the metaprofile and boxplot

### 3. Per-read TSS accessibility

File: [tables/${SAMPLE}_tss_read_accessibility.tsv.gz](tables/${SAMPLE}_tss_read_accessibility.tsv.gz)

Columns:

- `read_id`
- `tss_uid`
- `gene_id`
- `gene_id_base`
- `gene_name`
- `chromosome`
- `tss_coordinate`
- `tss_strand`
- `tss_source`
- `PI`
- `pausing_group`
- `n_positions_covered`
- `n_m6a_modified`
- `n_m6a_unmodified`
- `tss_accessible`
- `sample_name`
- `modbase_code`
- `modbase_label`
- `call_chromosome`
- `position`
- `call_strand`
- `mod_status`
- `call_state`
- `promoter_start`
- `promoter_end`
- `rel_pos`

The read-level function uses quickread modBAM calls within `TSS_ACCESS_WINDOW` bp of each TSS and marks `tss_accessible = TRUE` when a read carries at least one modified m6A call in that window.

### 4. Gene-level m6A summary

File: [tables/${SAMPLE}_gene_m6A_summary_quartiles.tsv](tables/${SAMPLE}_gene_m6A_summary_quartiles.tsv)

Columns:

- `tss_uid`
- `gene_id`
- `gene_id_base`
- `gene_name`
- `chromosome`
- `tss_coordinate`
- `tss_strand`
- `tss_source`
- `PI`
- `pausing_group`
- `m6A_modified_calls`
- `m6A_unmodified_calls`
- `m6A_total_calls`
- `m6A_positions_with_coverage`
- `m6A_fraction_modified`

This is the per-TSS, full-window m6A summary used for the boxplot and quartile statistics.

### 5. Gene-by-bin m6A summary

File: [tables/${SAMPLE}_gene_bin_m6A_summary_quartiles.tsv.gz](tables/${SAMPLE}_gene_bin_m6A_summary_quartiles.tsv.gz)

Columns:

- `tss_uid`
- `gene_id`
- `gene_id_base`
- `gene_name`
- `chromosome`
- `tss_coordinate`
- `tss_strand`
- `tss_source`
- `PI`
- `pausing_group`
- `bin`
- `m6A_modified_calls`
- `m6A_unmodified_calls`
- `m6A_total_calls`
- `m6A_fraction_modified`

This table keeps every TSS across every bin in the promoter window, with missing bins filled as zero counts and `NA` fractions when there is no coverage.

### 6. Histogram of modified calls per gene per bin

File: [tables/${SAMPLE}_gene_bin_m6A_histogram_quartiles.tsv](tables/${SAMPLE}_gene_bin_m6A_histogram_quartiles.tsv)

Columns:

- `pausing_group`
- `bin`
- `modified_call_bin`
- `n_genes`
- `n_tss`

The histogram bins the number of modified calls observed per gene per bin, capped at `MAX_MODIFIED_HIST_BIN`.

### 7. Quartile metaprofile

File: [tables/${SAMPLE}_m6A_metaprofile_quartiles.tsv](tables/${SAMPLE}_m6A_metaprofile_quartiles.tsv)

Columns:

- `pausing_group`
- `bin`
- `modified_calls`
- `total_calls`
- `n_tss`
- `fraction_modified`
- `modified_calls_per_tss`

This is the quartile-level metaprofile used for the main line plot.

### 8. Quartile summary statistics

File: [tables/${SAMPLE}_m6A_group_summary_quartiles.tsv](tables/${SAMPLE}_m6A_group_summary_quartiles.tsv)

Columns:

- `pausing_group`
- `n_genes`
- `n_tss`
- `median_PI`
- `mean_m6A_fraction`
- `median_m6A_fraction`
- `mean_total_calls`
- `median_total_calls`

This table gives one row per PI quartile.

## Output Plots

- [plots/${SAMPLE}_m6A_metaprofile_quartiles.pdf](plots/${SAMPLE}_m6A_metaprofile_quartiles.pdf)
- [plots/${SAMPLE}_gene_m6A_fraction_boxplot_quartiles.pdf](plots/${SAMPLE}_gene_m6A_fraction_boxplot_quartiles.pdf)
- [plots/${SAMPLE}_gene_bin_m6A_histogram_quartiles.pdf](plots/${SAMPLE}_gene_bin_m6A_histogram_quartiles.pdf)

## Running

Option 1: SLURM

```bash
sbatch run_m6a_pausing_quartiles.sh
```

Option 2: direct Rscript with environment variables

```bash
export SAMPLE=AL10_bc2178_19130
export CHROMS=AUTO
export WINDOW_BP=1000
export BIN_SIZE=10
export MOD_PROB_THRESHOLD=0.9
export MODBAM_TSS_CHUNK_SIZE=200
export TSS_ACCESS_WINDOW=100
export MAX_MODIFIED_HIST_BIN=10
export TSS_SOURCE_FILTER=ALL
export POLII_ROOT=/project/spott/cshan/fiber-seq/results/PolII
export PAUSING_PATH=${POLII_ROOT}/annotations/pausing_index_principal_with_CAGE_TSS_all_genes.tsv
export OUTPUT_ROOT=${POLII_ROOT}/m6a_pausing_quartiles
export BAM=/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples/${SAMPLE}.5mC.6mA.aligned.phased.bam

Rscript run_m6a_pausing_quartiles.R
```

## Environment Variables

- `SAMPLE`: sample name used in file naming, default `AL10_bc2178_19130`
- `CHROMS`: comma-separated list of chromosomes or `AUTO` for standard chr1-22,X,Y
- `WINDOW_BP`: TSS window size in bp, default `10000`
- `BIN_SIZE`: bin size in bp, default `100`
- `MOD_PROB_THRESHOLD`: probability threshold for m6A calls, default `0.9`
- `MODBAM_TSS_CHUNK_SIZE`: number of TSS per chunk, default `250`
- `TSS_ACCESS_WINDOW`: half-window used for per-read accessibility, default `100`
- `MAX_MODIFIED_HIST_BIN`: cap for histogram bin labels, default `10`
- `TSS_SOURCE_FILTER`: `ALL`, `CAGE`, or `GENCODE`, default `ALL`
- `POLII_ROOT`: base output folder, default `/project/spott/cshan/fiber-seq/results/PolII`
- `PAUSING_PATH`: pausing index annotation table
- `OUTPUT_ROOT`: output folder for quartile results
- `BAM`: path to the modBAM file with m6A calls

## Function Map

### Functions defined in [run_m6a_pausing_quartiles.R](run_m6a_pausing_quartiles.R)

- `get_env_chr(name, default)`: read a string environment variable with fallback.
- `get_env_int(name, default)`: read an integer environment variable with fallback.
- `get_env_num(name, default)`: read a numeric environment variable with fallback.
- `standard_chromosomes(x)`: keep chr1-22, X, and Y only.
- `strip_gene_version(x)`: remove Ensembl version suffixes.
- `make_pausing_quartiles(pausing_path, tss_source_filter)`: load PI annotations, choose TSS anchors, and assign quartiles.
- `make_tss_windows(pausing_groups, chrom_name, window_bp, seqinfo)`: build promoter windows and metadata for one chromosome.
- `append_tsv(dt, path, wrote_once)`: append a table chunk to a TSV file.
- `summarize_group_meta(aligned_dt, groups_dt, window_bp, bin_size)`: build quartile metaprofiles.
- `summarize_gene_m6a(aligned_dt, groups_dt)`: summarize m6A at the TSS level.
- `summarize_gene_bin_m6a(aligned_dt, groups_dt, window_bp, bin_size)`: summarize m6A by TSS and bin.
- `summarize_gene_bin_histogram(gene_bin_dt, max_modified_bin)`: count modified-call histograms.
- `label_read_tss_accessibility(...)`: compute per-read accessibility around each TSS.
- `plot_m6a_quartile_meta(...)`: write the quartile metaprofile PDF.
- `plot_gene_m6a_boxplot_quartiles(...)`: write the gene-level boxplot PDF.
- `plot_gene_bin_histogram_quartiles(...)`: write the histogram heatmap PDF.

### Shared helpers sourced from [modbam_footprintR_functions.R](../../footprintR_modbam_code/modbam_footprintR_functions.R)

- `read_bam_seqinfo(bam_path)`: read chromosome names and lengths from the BAM header.
- `coerce_iranges_coord(x, label)`: validate and coerce coordinates for GenomicRanges.
- `read_modbam_summary_positions(...)`: call `footprintR::readModBam` at summary level and return per-position m6A counts.
- `getBinarizedModificationsLongTable(...)`: convert quickread mod probability matrices into a long binary-call table.
- `write_gz_tsv(x, path)`: write a gzipped TSV.
- `make_chapter5_example_plot(...)`: produce the example TSS footprintR plot when requested.

## Notes

- The quartile boundaries are based on the PI distribution after filtering to genes with valid PI values and standard chromosomes.
- `tss_accessible` is defined at the read/TSS level, not at the position level.
- The read-level accessibility table is generated from quickread modBAM data, while the TSS-aligned summary table is generated from summary-level modBAM data.
