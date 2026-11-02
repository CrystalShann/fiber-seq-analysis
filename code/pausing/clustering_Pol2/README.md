# Pooled Pol II footprint-positive clustering

One global Leiden–Manhattan analysis of all 31 LCL samples. Only molecules with a
qualifying footprint enter clustering. There are no QC reports or excluded-molecule
tables. Necessary input validation and a memory check remain.

## Selection and m6A proportion

Use the dominant FANTOM5 TSSs (chr1–22, X) without a TPM filter. A TSS anchor must
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
these sampled molecules enter the single global run; other eligible molecules are not
assigned clusters. Sample/cell-line metadata remain attached; no per-sample clustering
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
the median over its TSSs, each TSS counted once. Clusters are displayed lowest median
PI first; clusters without a PI follow in Leiden size order. Cluster names stay the
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
sbatch code/pausing/clustering_Pol2/run.sh
```

`config.yaml` contains input paths, sample size and the memory budget. The Slurm default is
32 GiB with a 24 GiB estimated clustering budget. The full distance matrix alone
costs 8*N^2 bytes; the check budgets additional workspaces. Large cohorts may need
more memory and runtime than this queue allocation.

Stages can be resumed independently; use a compute allocation for extraction and clustering:

```bash
bash code/pausing/clustering_Pol2/run.sh --stage prepare
sbatch --cpus-per-task=4 code/pausing/clustering_Pol2/run.sh --stage count
bash code/pausing/clustering_Pol2/run.sh --stage features --chrom chr1
bash code/pausing/clustering_Pol2/run.sh --stage features
bash code/pausing/clustering_Pol2/run.sh --stage cluster
bash code/pausing/clustering_Pol2/run.sh --stage report
```

The last two stages require all configured chromosome chunks. Caches check input
paths, sizes, modification times and source hashes. Shared parsers do not write
their usual QC files because this workflow filters candidate reads before matrix
construction and never invokes the shared matrix cache writer.

`count` only checks the specified footprint overlap, sample/read identity
and full-window coverage. It writes `tables/eligible_molecules.rds` and prints the total;
it is independent of the m6A denominator and bin-completeness requirement. For
parallel feature extraction after `prepare`, submit `sbatch --array=1-23%4
code/pausing/clustering_Pol2/run.sh`; array indices 1–22 and 23 select chr1–22 and
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
- `cluster_profiles.pdf`: mean 50-bp-bin m6A proportion per cluster.
- `footprint_positions.pdf`: candidate footprint midpoints per cluster.
- `footprint_molecules.pdf`: each clustered molecule's 30-60 bp footprints, one row per
  molecule in the heatmap's row order (`stack_rows()`, `plot_molecules()`).
- `fig1c_pol2_clusters.pdf`: heatmap, per-cluster profiles with median PI, a pausing
  index wedge (lowest PI at the top) and Pol II ChIP signal at the TSS, in the layout
  of Fig. 1c of the nano-NOMe-seq paper.

Clusters are not given state names. The footprint criterion is candidate evidence for
Pol II association, not direct protein identification. Selecting footprint-positive
molecules means the reported fractions cannot estimate the abundance of states among
all TSS molecules.

## Verification

`test_pol2_clustering.R` covers all-full-window denominators and eligible-only rows, zero-call retention,
partial spans, strand reversal, per-bp calls, footprint size/window requirements,
sample-specific footprints, random sampling without duplicate molecules, bigWig missing
coverage, the TSS-region coordinates, the Ensembl PI join, TSS-signal ratio direction
and PI ordering. Run using the R
environment exported in `run.sh`. No production inputs are modified by the fixtures.
