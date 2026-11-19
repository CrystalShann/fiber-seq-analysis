# Nucleosome-positioning entropy in LCL Fiber-seq

This analysis implements the Shannon-entropy calculation in Li et al. (2026),
*Comprehensive profiling of chromatin occupancy dynamics through the cell cycle*,
[Methods: Nucleosome disorganization quantification](https://doi.org/10.1093/nar/gkaf1385).
Fiber-seq nucleosome coverage fractions replace the paper's two-dimensional
MNase-seq kernel scores. This is an adaptation, not an exact reproduction of the
published occupancy-scoring procedure.

## Run

From `/project/spott/cshan/fiber-seq`:

```bash
# All 150 example genes on Slurm (4 CPUs, 8 GB, 1 hour)
mkdir -p results/logs
sbatch code/entropy/nucleosome/run_entropy.sh

# Regenerate PDF figures and supporting tables from existing entropy results
bash code/entropy/nucleosome/run_entropy_features.sh

# Extract promoter footprint fractions against existing gene entropy, then plot
sbatch code/entropy/nucleosome/run_entropy.sh --footprints-only

# Recalculate a single gene and each of its saved Leiden clusters in a separate folder
ENTROPY_OUT="$PWD/LCL_project/entropy/CD74" \
  bash code/entropy/nucleosome/run_entropy.sh --by-cluster CD74

# Gene-only entropy, skipping original cluster aggregation
ENTROPY_OUT="$PWD/LCL_project/entropy/CD74_gene_only" \
  bash code/entropy/nucleosome/run_entropy.sh --gene-only CD74
```

`ENTROPY_OUT` overrides the default output directory, `LCL_project/entropy`.
`ENTROPY_WORKERS` overrides the allocated Slurm CPU count; local execution defaults
to one worker. Optional trailing gene names select a subset. The runner uses the
same R installation and libraries as the existing PolII workflow.
Successful entropy runs also execute the plotting notebook. Plotting can run
independently on saved partial results: invalid genes are excluded. The current full run has 148 valid genes; SLC44A4 has
no saved cohort and CLEC12B failed the nucleosome alignment consistency check.
`--by-cluster` explicitly enables per-gene Leiden-cluster entropy (also the default,
preserving prior behavior); `--gene-only` disables cluster entropy.
These options do not recluster fibers. Labels come from the saved recurrent-pattern
read table, falling back to the original `<group>/<gene>_read_clusters.tsv`.
Every original cluster is retained. Programmatically use `calculate_gene(..., by_cluster=TRUE)`.
For selected genes, use a separate `ENTROPY_OUT` as above to preserve the full run.
`--footprints-only` reads existing `tables/gene_entropy.tsv` and updates the
footprint fractions without recalculating entropy; omit gene names to refresh all
saved genes. Normal successful entropy runs also calculate these fractions.

## Inputs and shared functions

The gene list is `code/pausing/PolII/example_gene_list.tsv`. TSSs and fiber identities
come from `LCL_project/pausing/PolII/read_features/tables/{gene_features.tsv,read_features.tsv.gz}`.
These retain the existing FANTOM5 LCL dominant TSS, or GENCODE v46 canonical TSS
fallback, and the fibers spanning TSS +/-1,000 bp. Thus entropy, cluster assignments,
and saved accessibility features refer to the same cohort. This does not recruit
additional fibers covering only the 500 bp entropy window. In the current inputs,
149 genes have 75,276 gene/fiber observations; SLC44A4 has no usable saved cohort
and is retained with `NA` entropy and `status=no_fibers`. A fiber overlapping two
genes can contribute separately to each gene.

Original nucleosome calls are queried from each sample's indexed fibertools BED12,
using the 31-sample LCL metatable. The exact alignment coordinates must match the
saved cohort. Fibers with present BED records but no retained nucleosome calls
remain in the denominator; missing records cause an explicit error. Nucleosome
lengths are 130–160 bp inclusive, matching the existing PolII analysis. This size
selection is an input convention of this project, not a requirement of Li et al.

Existing functions are sourced and never changed:

- `parsing_functions/parsing_footprints_functions.r`: `read_ft_bed12()` (and
  `read_tabix_region()`), `convert_ft_bed12_to_bed6()` for sentinel removal, and
  `assign_size_class()` for nucleosome selection.
- `nucleosome_positioning/process_nuc_pos.R`: `build_region_anchors()`,
  `build_read_table()`, `rel_interval()`, `occupancy_profile()` with
  `bin_size=1, metric="occupancy"` (internally `covered_positions()`), and
  `read_position_summary()`.

The existing `extract_region_footprints()` loader expects FiberHMM BED13 in its
explicit-path interface, so a small adapter queries the LCL fibertools BED12
through the shared parser instead. No files in either shared directory are edited.

## Calculation

The TSS base is offset 0. Exactly offsets **+1 through +500**, inclusive, are used:
genomic `[TSS+1, TSS+500]` on `+`, and `[TSS-500, TSS-1]` on `-`, expressed as
1-based closed coordinates. BED calls remain 0-based half-open until converted by
the shared coordinate helper. Read alignment strand does not determine gene orientation.

At each nucleotide, `occupancy = nucleosome_fibers / covering_fibers`. The shared
occupancy function counts the union of a fiber's overlapping nucleosome calls;
it does not count midpoints or multiply-count overlapping calls. Coverage uses
the saved alignment spans, consistent with the existing pipeline, rather than
CIGAR-level alignment blocks. The current full-span cohort has constant depth
within each group, but the estimator also supports varying span coverage.

```text
total_nucleosome_occupancy_signal = sum(occupancy[1:500])
p[i] = occupancy[i] / total_nucleosome_occupancy_signal
nucleosome_entropy_bits = -sum(p[i] * log2(p[i]), over p[i] > 0)
```

No pseudocounts, smoothing, midpoint counting, binning, or entropy normalization
are applied. `total_nucleosome_occupancy_signal` is a sum of fractions, with maximum
500, not the number of nucleosomes or total occupied fiber-bases. If all occupancy
is zero, entropy is `NA` and the total is zero. If any of the 500 bases has no
coverage, both entropy and total are `NA` (`incomplete_coverage` or `no_fibers`),
because an unobserved base is not an observed zero. Coverage and status columns
make these cases explicit. The maximum entropy is `log2(500) = 8.965784` bits.

Each gene is analyzed over all its retained fibers. Saved recurrent-cluster labels
are joined by sample-qualified RID from
`LCL_project/pausing/PolII/recurrent_patterns/tables/<gene>/read_features.tsv.gz`.
Every original cluster is analyzed separately, including ineligible clusters;
their `eligible` flag is retained.
An absent label file leaves gene entropy available and sets
`cluster_labels_available=FALSE`. An inconsistent label cohort is an error.

## Outputs

All tabular outputs are in `LCL_project/entropy/tables/`; all figures are PDF files
in `LCL_project/entropy/nucleosome_entropy_plots/`.

- `gene_entropy.tsv`: exactly one row per requested gene.
- `nucleosome_entropy.tsv`: gene and original cluster rows, identified
  by `level`; includes fiber counts, total occupancy signal, entropy in bits,
  coverage diagnostics, TSS metadata, expression, and supporting features.
- `occupancy_per_base.tsv.gz`: 500 rows per analyzed group, including genomic
  position, transcription-oriented offset, denominator, numerator, occupancy and
  normalized probability. Zero-occupancy groups have undefined probabilities.
- `gene_entropy_features.tsv`: valid genes joined to Puffin motif classes.
- `gene_group_summary.tsv`: group counts, medians, quartiles and ranges.
- `puffin_motif_entropy_summary.tsv`: motif-class summaries pooled across gene groups.
- `highest_lowest_5_gene_entropy.tsv`: the ten genes displayed in the Rmd table.
- `short_footprint_fractions.tsv`: one row per gene and short-footprint size class,
  including all statuses, denominator, numerator, fraction, and exact window.
- `entropy_vs_plus1_start.tsv`, `entropy_vs_short_footprints.tsv`: plotted observations
  for the two gene-level scatter plots.

`code/entropy/entropy_features.rmd` defines nucleosome figures and molecule-cluster figures:

1. `01_gene_group_entropy.pdf`: gene dots and boxplots, with median diamonds,
   printed median entropy, valid gene counts and the theoretical maximum.
2. `02_entropy_by_puffin_motif.pdf`: one violin per motif, all gene groups pooled,
   ordered by median entropy. Each dot is one gene; labels show gene counts,
   boxes show the IQR and black diamonds show medians.
3. `03_entropy_vs_plus1_start.pdf`: gene entropy on x and mean +1 nucleosome start
   position on y, colored by gene group. The mean uses only fibers with a start
   at +1..+299 and is a position rather than a variability statistic.
4. `04_entropy_vs_short_footprints.pdf`: gene entropy on x and fraction of fibers
   with a short footprint overlapping **-100..+100 bp** on y, colored by gene
   group and faceted by short-footprint size class. Entropy remains +1..+500 bp.

Puffin annotations are taken from `04a_puffin_motif_scores.tsv` in the existing
PolII promoter-annotation results. Chromosome, strand and TSS must match. Motif
labels belong to the gene sequence, not individual clusters. The motif violin pools gene groups, so differences can reflect gene-group
composition. Motif comparisons remain descriptive.
The notebook uses shared `theme_fiberseq()` and `save_figure()`
without editing shared functions. The runner executes
chunks with knitr and removes its temporary Markdown output, producing no HTML,
PNG or compiled notebook report in the output folders.

Runs do not generate session-info or run-metadata files.

## Short-footprint fraction

Helpers in `nucleosome_entropy.R` query each sample's indexed FiberHMM v1 file:
`FiberHMM/merged/<sample>/joint_trained_results/<sample>_<chrom>_fp.bed.gz`.
The legacy BED10 records are adapted to BED12 for the existing
`convert_ft_bed12_to_bed6()` with `drop_sentinels=FALSE`, then classified with the
shared `assign_size_class()`. First and last FiberHMM blocks are real calls.
The per-sample calls reproduced the existing CD74 and HLA-DRA per-fiber core
footprint occupancy measurements to floating-point precision in a pilot check.

The footprint window includes both endpoints and the TSS: 201 bases from -100
through +100. A footprint qualifies with at least 1 bp of overlap; its full size,
before clipping, determines its class. The shared class boundaries are
10 <= length < 30, 30 <= length < 60, and 60 <= length <= 80 bp. Each fiber is
counted at most once per class even if it has several footprints; a fiber can
have qualifying calls in multiple classes. The denominator includes every retained
fiber spanning the full footprint window, including fibers with no short calls.
Missing or ambiguous FiberHMM records produce an explicit error, not a zero.
Genes without valid entropy are retained with `footprint_status=entropy_unavailable`.
The helper `short_footprint_fraction()` also handles partial-span fibers by
excluding them from the denominator; current retained promoter fibers all span
the full window. These fractions are presence/absence per fiber, not occupancy
per base. The raw entropy nucleosome selection remains 130-160 bp.

`plus1_start_sd_bp` uses the saved PolII definition: earliest nucleosome starting
at offsets +1..+299, excluding a nucleosome covering the TSS. Its count is
`n_plus1_fibers`. `plus1_nearest_midpoint_sd_bp` uses the shared positioning
helper's definition: the nearest midpoint at or downstream of the TSS among
calls overlapping the +/-1 kb promoter window; this can include a TSS-crossing
nucleosome. Its separate count is `n_plus1_nearest_midpoint_fibers`. Both standard
deviations require at least two observed positions; fibers without a qualifying
nucleosome do not supply a position.

`promoter_m6a_core` is the mean saved per-fiber m6A fraction at observed sites over
-200..+200. `tss_nucleosome_free_fraction` is the fraction of fibers without a
retained nucleosome covering the TSS; it is distinct from m6A accessibility.
TPM values come from the supplied example gene list.

## Interpretation

Higher entropy means a flatter aggregate occupancy profile across the 500 bases;
lower entropy means occupancy is concentrated at fewer bases. These are spatial
occupancy entropies, not entropies of cluster membership or initiation-site usage.
A fixed 147 bp footprint alone already gives `log2(147)` bits, and a completely
flat profile has maximal entropy regardless of its absolute occupancy level.
Therefore this adaptation is a proxy for positioning disorganization, not a
direct measurement of fiber-to-fiber positional variance. Interpret it alongside
total occupancy, depth, +1 variability and accessibility. Pooling fibers weights
samples by their fiber counts; it is not equal weighting of donors. Small clusters
can produce less stable estimates. Gene-level expression is repeated for its
cluster rows; those repeated values are not independent gene observations.


Selected-cluster examples (06-11) are single-figure PDFs generated by the shared Rmd,
two per cluster: `*_m6a.pdf` is an m6A heatmap with methylation track and
`*_footprints.pdf` shows nucleosome and short-footprint intervals with the methylation track. Only the saved IL21R
cluster5 (36 fibers), CTCF cluster4 (70), or TNFRSF13B cluster1 (188) is shown.
Views span TSS +/-1000 bp; the entropy labels still refer to +1..+500 bp.
The notebook reuses `plot_example_m6a()` and `plot_example_footprints()` wrappers
and the unchanged shared `plot_read_heatmap()` and `plot_molecules()` functions.
Indexed sample-specific FiberHMM v1 files supply short footprints; exact saved
alignments supply nucleosomes. The nucleosome counts are checked against the
saved per-base entropy inputs before plotting. A summary table is displayed in
the Rmd.
