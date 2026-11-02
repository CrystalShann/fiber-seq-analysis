# Pooled Pol II footprint-positive clustering

One global Leiden–Manhattan analysis of all 31 LCL samples. Only molecules with a
qualifying footprint enter clustering. There are no QC reports or excluded-molecule
tables. Necessary input validation and a memory check remain.

## Selection and m6A proportion

Use the dominant FANTOM5 TSSs (chr1–22, X). A TSS anchor must
overlap an ENCSR000EAD GRCh38 Pol II IDR peak (`ENCFF886PSD`). Each TSS also gets
`pol2_tss_signal`, the mean `ENCFF328MMS` fold-enrichment over its covered bases at
−50..+299 (Day's TSS region, strand-oriented); missing coverage stays NA.
A qualifying footprint is a FiberHMM v2 TF block 30–60 bp long (inclusive) that
overlaps the 2 kb TSS window, taken from the molecule's own sample
(`<footprint_dir>/<S>/<S>.recalled_tf.bed.gz`, BED15, one block per footprint, tq < 50
already dropped upstream). A molecule may carry any number of footprints. There is no
downstream-position requirement. Identical interval/read-name records are counted once.
Because footprints are read per sample, a read name repeated in another sample never
qualifies the wrong molecule.

The molecule must span the full window, relative positions -1000 through +999.
All full-window molecules define the observed-site denominator; only
footprint-positive full-span molecules enter clustering:

```
whole-window proportion = molecule's distinct m6A call count / pooled observed sites
50-bp feature           = molecule's distinct calls in bin / pooled observed sites in bin
```

The site set is the union of called positions across all full-window molecules of all
31 samples **at that TSS**, using the shared BED12 parser and the existing
`bin_read_matrix()`. It is not all A/T reference positions. Fully covered positions
without a call are zero; eligible zero-call reads are retained if other full-window
reads define all bins. TSSs with an empty observed-site bin are omitted rather than
imputing zero. The global matrix has 40 aligned 50-bp columns. No per-read scaling
is applied. Minus-strand positions are reversed before binning. All eligible
read/TSS pairs define features before sampling; their per-bp calls are kept for the
molecule heatmap.

Feature chunks are extracted by chromosome (100 neighbouring TSSs per tabix query),
then joined for **one** graph and partition. Draw up to 5,000 unique sample/read
molecules without replacement, uniformly from all 31 samples combined, with seed 1.
For a molecule qualifying at multiple TSSs, choose one eligible TSS at random. Only
these sampled molecules enter the single global run; other eligible molecules receive nearest-neighbour assignments in the `project` stage. Sample/cell-line metadata remain attached; no per-sample clustering
or balancing is performed. Existing `manhattan_knn_graph()`,
`leiden_partition()`, cluster relabeling and profile functions are called directly:
k=10, resolution=1, seed=1, automatic kernel width. This preserves the exact dense
distance calculation. A conservative memory estimate stops a run that cannot fit;
saved features remain available. The graph is exact on the requested random sample.
Bulk Pol II fold-enrichment (`ENCFF328MMS`) is measured over the footprints of the
sampled molecules only; missing bigWig signal remains missing.

## Cluster order, PI and Pol II signal at the TSS

PI comes from `LCL_project/pausing/pause_index.tsv` (Day et al. 2016, RefSeq isoforms)
and is joined by Ensembl gene ID (`ensg`, versions ignored). Its RefSeq TSS can differ
from the FANTOM5 TSS used here; genes without a finite PI stay NA. A cluster's PI is
the median over its TSSs, each TSS counted once. The enrichment composite displays clusters lowest median
PI first; clusters without a PI follow in Leiden size order. All other figures
display clusters in numeric order (cluster1, cluster2, …). Cluster names stay the
Leiden size labels (`cluster1` = largest) and keep their colours in every figure.

Because every TSS is in a Pol II peak, the Fig. 1c bars use Pol II ChIP signal instead
of peak enrichment: `pol2_tss_signal` of each molecule's TSS, log2 of the cluster mean
over the mean of all clustered molecules, Wilcoxon test of the cluster against all
others, BH over clusters (`adjust_bh()`), with stars (*** < 0.001, ** < 0.01,
* < 0.05). Molecules without coverage are left out. Molecules from the same TSS share
its signal, so the test treats them as independent.

## Run

From the project root:

```bash
# Runs prepare, extraction of all chromosomes, global clustering and reporting.
sbatch code/pausing/clustering_Pol2/00_run.sh
```

`pol2_chip_cluster_utils/config.yaml` contains input paths, sample size and the memory budget. The Slurm default is
32 GiB with a 24 GiB estimated clustering budget. The full distance matrix alone
costs 8*N^2 bytes; the check budgets additional workspaces. Large cohorts may need
more memory and runtime than this queue allocation.

Stages can be resumed independently; use a compute allocation for extraction and clustering:

```bash
bash code/pausing/clustering_Pol2/00_run.sh --stage prepare
sbatch --cpus-per-task=4 code/pausing/clustering_Pol2/00_run.sh --stage count
bash code/pausing/clustering_Pol2/00_run.sh --stage features --chrom chr1
bash code/pausing/clustering_Pol2/00_run.sh --stage features
bash code/pausing/clustering_Pol2/00_run.sh --stage cluster
sbatch --cpus-per-task=8 code/pausing/clustering_Pol2/00_run.sh --stage project
sbatch --cpus-per-task=8 code/pausing/clustering_Pol2/00_run.sh --stage tracks
bash code/pausing/clustering_Pol2/00_run.sh --stage report
```

The last two stages require all configured chromosome chunks. Caches check input
paths, sizes, modification times and source hashes. Shared parsers do not write
their usual QC files because this workflow filters candidate reads before matrix
construction and never invokes the shared matrix cache writer.

`count` only checks the specified footprint overlap, sample/read identity
and full-window coverage. It writes `tables/eligible_molecules.rds` and prints the total;
it is independent of the m6A denominator and bin-completeness requirement. For
parallel feature extraction after `prepare`, submit `sbatch --array=1-23%4
code/pausing/clustering_Pol2/00_run.sh`; array indices 1–22 and 23 select chr1–22 and
chrX. Run `cluster` once after every extraction task has finished.

## Outputs and interpretation

`LCL_project/pausing/clustering_Pol2/` holds `pol2_clustering_report.html`,
`run_config.yaml` and `session_info.txt` at the top level, and two subfolders.

`tables/` holds every RDS and TSV:

- Stage caches: `prepared.rds`, `features/<chrom>.rds`, `clustering.rds`,
  `annotations.rds`, and `eligible_molecules.rds` from `count`.
- `molecule_assignments.tsv.gz`: cluster assignments, with PI and `pol2_tss_signal`.
- `candidate_footprints.tsv.gz`: candidate footprints, with Pol II signal and coverage.
- `bin_denominators.tsv.gz`: observed-site denominators.
- `tss_cluster_fractions.tsv.gz` and `sample_cluster_fractions.tsv`: per-TSS and
  per-sample cluster fractions.
- `cluster_profiles.tsv`: per-cluster 50-bp-bin profiles.
- `cluster_summary.tsv`: size, median PI, TSS signal ratio, Wilcoxon FDR and display order.

`plots/` holds the PDF figures. They are drawn by
`pol2_clustering_plots.R` with the shared Leiden–Manhattan helpers
(`plot_read_heatmap()`, `plot_group_profile()`, `plot_stacked_proportion()`,
`cluster_id_palette()`, `save_figure()`):

- `molecule_heatmap.pdf`: per-bp m6A calls of every clustered molecule, rows split by
  cluster, cluster and sample bars, per-cluster call-fraction lines, TSS dashed.
- `footprint_positions.pdf`: candidate footprint midpoints per cluster.
- `footprint_molecules.pdf`: each clustered molecule's footprints, one row per molecule
  in the heatmap's row order (`stack_rows()`, `plot_molecules()`), coloured by
  `MOLECULE_TRACK_COLORS` class. The data come from `tables/tracks/molecule_footprints.tsv.gz`,
  written by the `tracks` stage (`pol2_molecule_tables.R`) from the molecule's own sample:
  - fibertools nucleosome calls of 130-160 bp (`nuc_by_chr`, longest alignment);
  - FiberHMM v2 `recalled_tf` blocks binned as 10-30, 30-60 (inclusive, the clustering
    rule) and 60-80 bp.
  Other sizes are left out. `report` needs this table.
- `pol2_clusters_and_enrichment.pdf`: heatmap, per-cluster profiles with median PI, a pausing
  index wedge (lowest PI at the top) and Pol II ChIP signal at the TSS, in the layout
  of Fig. 1c of the nano-NOMe-seq paper.
- `loci/<gene>_TSS_molecules.pdf`: one figure per gene in `locus_genes`. It shows every
  eligible molecule at the gene's TSS, not only the sampled ones. Clustered molecules
  keep their Leiden cluster. The others take the majority cluster of their `k_neighbors`
  nearest clustered molecules (Manhattan distance, 40 bins). The figure has genes over
  ±250 kb (`gene_models`), the gene's in-window GTEx eVariants (`eqtl_pairs`), and the
  cluster proportions at the TSS versus all clustered molecules. Below those are the
  per-cluster per-base m6A fraction and one row per molecule. It is drawn by the shared
  `plot_locus_molecules()` in `parsing_functions/plotting_functions.r`, which also draws
  the HSCA domain/CTCF figures; the CTCF, loop and domain tracks are left out here.

Clusters are not given state names. The footprint criterion is candidate evidence for
Pol II association, not direct protein identification. Selecting footprint-positive
molecules means the reported fractions cannot estimate the abundance of states among
all TSS molecules.

## Imprinting and cluster features

`project` writes `tables/tracks/molecule_projection.tsv.gz` for every eligible read/TSS
observation. It preserves sampled Leiden assignments and projects other observations
by k-nearest-neighbour voting in batches, parallel across chromosomes.

`plots/imprinting_by_cluster.pdf` shows the fraction of unique physical molecules
at DMR-overlapping TSS windows per cluster, and cluster composition at imprinted
versus other TSSs. A molecule seen in both groups belongs only to the imprinted
composition bar. Counts are descriptive, with no statistical test.
`tables/imprinting_loci.tsv` includes DMR names, molecule/sample counts, cluster4
fraction and fold over its overall fraction, dominant cluster, TPM, PI and Pol II.
The panel data are saved as `imprinting_by_cluster.tsv` and `imprinting_composition.tsv`.

`plots/cluster_features.pdf` summarizes only the 5,000 clustered molecules:

- m6A: stored per-molecule fraction and pooled sum(calls) / sum(positions).
- CpG: pb-CpG-tools count mode, reference CpGs, minimum coverage 1. Cluster BAMs
  contain one longest matching alignment per molecule, retaining modification tags.
  The first sample header is restricted to selected chromosomes and reference IDs
  are remapped, since unused patch contigs are absent from `hg38.fa`.
  Sites are restricted to that cluster's TSS windows. Boxplots use the cluster/TSS
  fraction assigned to each molecule; the pooled fraction counts each genomic site
  once per cluster, avoiding double-counting overlapping windows.
- GRO-cap: full-window sense-strand read sum, with absent reads set to zero.
- DNase: full-window mean. Both signal panels use log10(x + 1); diamonds show
  molecule-weighted means, while methylation diamonds show pooled fractions.

CpG extraction is embedded in plotting; temporary BAMs and raw CpG outputs are
removed on exit. Result tables `cluster_features.tsv`, `cluster_cpg_by_tss.tsv`,
and `cluster_feature_molecules.tsv.gz` support re-plotting without extraction.
Locus examples are GSDME, STING1 and CASP4. The shared profile label is “m6A fraction”.

`merged_bam_overrides` corrects the obsolete NA18916 BAM filename locally without
changing the shared sample sheet or invalidating clustering caches.

## Per-cluster autocorrelation from saved assignments

```bash
bash code/pausing/clustering_Pol2/00_run.sh --stage autocorrelation
```

The report computes this analysis if its summary table is absent; the standalone
stage explicitly recomputes it using the same report chunk. It uses the saved, strand-oriented
1-bp m6A matrix in `clustering.rds$sites`, matching sample/read/TSS identity to
`molecule_assignments.tsv.gz`. Clustering assignments, colors and the original
−1000..999 clustering window remain unchanged.

Two windows are analyzed independently: **upstream −500..−25 bp** and
**downstream +25..+500 bp**, inclusive (476 bases each). On the minus strand these
map to genomic TSS+25..+500 and TSS−500..−25, respectively. They are sliced from
the already oriented matrix, without reversing twice or concatenating across the
gap. Every molecule's mean and variance are recomputed separately for each window.

The input is binary m6A accessibility, not nucleosome occupancy or the 40-bin
features. Covered positions without a call are zero. The shared
`autocorrelations` function computes mean-centered linear ACF divided by
N × population variance, without lag-dependent overlap correction. Constant
molecules are excluded separately in each window; the other valid molecules
receive equal weight in their cluster means. No uncertainty estimates or counts
are plotted and no resampling is performed.

`tables/` contains:

- `autocorrelation_per_molecule_upstream.npz` and
  `autocorrelation_per_molecule_downstream.npz`: molecule × lag ACF, row identity,
  all lags 0..475, strand-relative input positions and window name. Undefined
  ACF rows remain NaN.
- `autocorrelation_molecules.tsv.gz`:
  molecule metadata, window boundaries, inclusion status and per-molecule peaks.
- `autocorrelation_cluster_summary.tsv.gz` and `autocorrelation_cluster_peaks.tsv`:
  window-specific mean curves and strongest positive local peaks in 140..250 bp,
  using the shared `repeat_peak` helper. Missing peaks remain missing. These are
  descriptive peaks, not proof of nucleosome spacing.
- `autocorrelation_provenance.rds`: input hashes, windows, normalization and policy.

`plots/autocorrelation/` contains six PDFs: overlay curves, separate cluster
panels and cluster-mean heatmaps, each with `_upstream_-500_to_-25_bp` or
`_downstream_+25_to_+500_bp` suffix.
All use lags 15..475 bp and numeric cluster order. The complete lag range remains
in the tables; no full-range figures are generated. Large-lag values attenuate
under the fixed normalization. The enrichment composite alone retains PI order.


## Script layout and shared analysis functions

- `00_run.sh`, `run_pol2_clustering.R`: stable command-line entry points.
- `pol2_chip_cluster_utils/`: configuration, original clustering/molecule helpers,
  and the BAM helper for cluster feature extraction.
- `select_tss/`: TSS preparation and eligible-molecule counting entry points.
- `pol2_clustering_report.Rmd`: reporting plus the autocorrelation,
  nucleosome-positioning and expression analysis chunks.
- `pol2_clustering_plots.R`: existing Pol II figure functions.
- The reference paper PDF remains at the root.

The root `pol2_clustering_functions.R` is a compatibility symlink to the grouped
original. Its path and unchanged contents remain part of saved cache signatures.
No clustering is rerun merely to reorganize the files.

There are no separate autocorrelation/nucleosome adapters or expression scripts.
The Rmd imports the original Python `autocorrelations`, `summarize`, and
`repeat_peak` functions, and sources the original named ACF plotting chunk from
`06_autocorrelation.Rmd` without running that report. Nucleosome analysis sources
`code/nucleosome_positioning/process_nuc_pos.R` and `plot_nuc_pos.R` directly.
No copied analysis-function definitions are maintained in the Pol II directory.

The standalone stages execute the named Rmd chunks, making the report the only
source of analysis orchestration:

```bash
bash code/pausing/clustering_Pol2/00_run.sh --stage autocorrelation
bash code/pausing/clustering_Pol2/00_run.sh --stage nucleosome
bash code/pausing/clustering_Pol2/00_run.sh --stage expression
```

The report uses existing autocorrelation/nucleosome summary tables when present.
Run the corresponding standalone stage after changing its settings or inputs.
Expression is recalculated from the configured TPM tables when the report renders.
Test and verification scripts have been removed as requested. Validation used
throwaway checks outside the project for exact ACF window equivalence, strand
coordinates, sample-distinct molecule keys, occupancy deduplication, zero calls,
midpoint selection, optional gap omission, and unchanged saved clustering signatures.

## Nucleosome positioning

Uses each molecule's saved FANTOM5 TSS, sample/read/observation identity and exact
strand-oriented −1000..+999 window. Per-sample `recalled_footprint.bed.gz` files
under the configured FiberHMM directory contain BED13 nucleosome blocks; all
blocks are used without a size or block-score filter. Full footprint boundaries
are retained. The original no-returned-call convention is preserved; there is no
separate missing-data analysis, and missing/unreadable inputs stop execution.

The original `occupancy_profile()` computes unique-molecule occupancy at positions
on a 10-bp grid, not bin averages. Profiles pool molecules within each cluster.
Its optional `metric="centers"` mode supplies raw inferred footprint-center
counts in 25-bp bins. Center-density and promoter-balanced outputs are omitted.
Centers are inferred from footprints, not directly measured histone dyads.

The original `read_position_summary(include_gap=FALSE)` chooses −1/+1 by
midpoint: nearest negative and nearest nonnegative center. A TSS-covering footprint
may be selected. Missing flanks remain NA. No inter-nucleosome gap is calculated.
TSS occupancy is reported separately as occupied/spanning molecules. Its 95%
intervals use 500 crossed promoter × biological-sample bootstrap draws, seed 1;
molecules are not resampled independently. Positional summaries include valid
−1/+1 counts and are also saved per promoter/sample stratum.

`tables/nuc/nucleosome_*.tsv.gz` includes molecules, full footprints, inferred centers,
occupancy profiles and raw center counts, cluster summaries, and promoter/sample
strata. `nucleosome_provenance.rds` records sources and definitions.
`plots/nuc/nucleosome_*.pdf` contains molecule occupancy, raw center counts,
midpoint distributions and TSS occupancy. Numeric cluster order and colors are
preserved. These are descriptions conditional on the existing Pol II accessibility
cluster selection, not independent validation of those clusters.


## Table organization

Core saved clustering inputs and assignments remain directly under `tables/`;
the chromosome feature cache remains in `tables/features/`.
Analysis outputs are grouped into `tables/autocorrelation/`, `tables/nuc/`,
`tables/cluster_features/` (including CpG-derived results), and
`tables/tracks/` (projected molecules and extracted footprints).
All corresponding readers and writers use these locations.

Imprinting summaries and expression summaries are computed in memory when needed.
Their former `imprinting_*.tsv` and `cluster_expression.tsv` exports are no longer
written. Expensive extracted data and per-molecule analysis results remain saved.


## Current source layout

Numbered runnable scripts and the report are at the workflow root:
`00_run.sh`, `01_run_pol2_clustering.R`, `02_select_tss.R`,
`03_count_pol2_eligible.R`, `04_pol2_cluster_bams.py`, and
`05_pol2_clustering_report.Rmd`.
Shared R libraries and `config.yaml` are in `pol2_chip_cluster_utils/`.
The README and reference paper remain at the root.

The root `pol2_clustering_functions.R` is only a compatibility symlink:
the existing clustering cache hashes this exact path. The actual library lives
in the utilities directory and is unchanged. The unnumbered Rmd symlink likewise
preserves the existing core report function's entry path.
