# Dominant-TSS fiber features

`03_fiber_features.R` implements Step 3 in R and defines no new functions. It
sources `code/parsing_functions/parsing_footprints_functions.r`, using its
tabix reader, BED12 reader, sentinel-removing block converter, and FIRE constants.

The gene input is a TSV with `gene_id`, `chrom`, `strand`, `anchor_tss`,
`anchor_tes`, and `pause_peak_rel`. Anchors and TES are 1-based genomic positions;
missing pause peaks use relative +20..+80. The default input is
`LCL_project/pausing/dominant_tss/pause_index_dominant.tsv`.

## Generate dominant TSSs and pausing measurements together

From the repository root:

```bash
mkdir -p results/logs
sbatch code/pausing/run_pause_index.sh
```

This produces both `LCL_project/pausing/pause_index.tsv` (canonical anchors) and
`LCL_project/pausing/dominant_tss/pause_index_dominant.tsv` (CAGE anchors).
The default initiation input is the local hg38 FANTOM5 GM12878 consensus,
`annotations/fantom5_LCL/fantom5.hg38.LCL.consensus.TSS.bed.gz` under
`/project/spott/cshan`. Each single-base site is supported in all three
replicates; its score is the summed tag count.

Sites are matched to GENCODE v46 transcript TSSs on the same strand within
250 bp. Each site is assigned to its nearest gene; equal-distance ties between
genes are excluded. For each gene, the highest-count assigned site is the
dominant TSS. Ties prefer smaller annotation distance, a canonical transcript,
lower genomic coordinate, then transcript ID. The TES comes from the matched
transcript; matches with the TES upstream of the proposed anchor are excluded.
Genes with no eligible CAGE site are omitted from the dominant output and counted
in the log. This is the strongest supported **single base**, not a sum over a
promoter cluster. No peak-level TPM filter is applied.

GRO-seq is lifted from hg19 to hg38 and quantified around both anchor sets in
one pass, using the existing -50..+300 promoter and +300..TES+3000 body windows.
Overlapping signal is weighted by overlap length when summing bigWig intervals.
`pause_peak_rel` is the strongest GRO-seq base at transcription-relative +20..+80;
ties choose the earliest relative position and absent signal yields `NA`.
All `anchor_tss`/`anchor_tes` coordinates are 1-based. The dominant table also
records the matched transcript, annotated TSS, distance, CAGE count and source path.
Both outputs cover chr1–22 and chrX; the GRO-seq input has no chrY.

Optional settings (the CAGE override must contain hg38 single-base BED6 or BED9
records with numeric tag counts in column 5):

```bash
sbatch code/pausing/run_pause_index.sh --tss-distance 250 \
  --cage /path/to/hg38_CTSS.bed.gz --out-dir /path/to/pausing
```

After the job succeeds, run `sbatch code/pausing/dominant_tss/run_fiber_features.sh`.
Clustering must use these same anchors and generated features; select chr1–22,X
for clustering to match the upstream chromosome set.

Selection and pause-peak tests: `Rscript code/pausing/dominant_tss/test_dominant_tss.R`
(run from the repository root with the wrapper's R library environment).

The authoritative sample manifest defaults to
`/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv`.
CSV or TSV manifests must contain `sample_name`, `fire_dir`, and `bam_file`.

Eligibility uses the indexed BAM: mapped primary alignments only, excluding
secondary and supplementary flags, with no MAPQ threshold. The reference span
must cover the whole inclusive gene window. CIGAR deletion/skipped-reference gaps
do not break span coverage. The m6A BED record must match the primary read name,
start, and end exactly. A longer supplementary alignment never replaces it.
Zero-call, sentinel-only m6A records remain eligible.

The default window is -1000..+2000 (3001 bases). Custom windows must include at
least -500..+500 for the fixed methylation and footprint features. Downstream
architecture is missing if the window does not include +2000; elongation is
missing when its full defined interval is unavailable. Chromosome-edge genes
whose full windows cannot be represented are skipped with a count.

Outputs are `fiber_features.{chr}.tsv.gz` and `fiber_qc.{chr}.tsv`. QC counts
unique primary reads overlapping or spanning at least one gene window; spanning
counts precede the m6A match. Missing FIRE counts are unique emitted reads with
no overlapping segmentation. `rows_written` counts gene/read pairs. Missing
FIRE data remain `other` with missing segmentation measurements; an absent
matching CpG record gives missing `mcpg_n`, whereas a sentinel-only record gives
zero.

Methylation site definitions follow the observed-position matrices used by
`assemble_region_m6a()` in the shared parsing code, `assemble_region_met_data()`
in `code/topic_model/topic_modelling_functions.r`, and `methylation_by_cluster()`
in `code/clustering_methods/Leiden_Manhattan/leiden_manhattan_plots.r`:

- `m6a_site_n`: number of distinct m6A call coordinates across all eligible
  fibers in all samples, within the full gene window.
- `cpg_site_n`: number of distinct CpG call coordinates across those fibers,
  within genomic TSS +/-500.
- `m6a_n` and `mcpg_n`: number of distinct called coordinates on each fiber in
  the corresponding window. Adjacent CpG call coordinates remain separate,
  matching the existing matrices; there is no reference-based G-to-C conversion.

The denominators are shared by every fiber for a gene, including fibers with
zero calls. They exclude positions found only on rejected alignments. A gene
with no calls has a zero denominator; a derived proportion is then undefined.
These count observed modification sites, not all possible A/T or CpG sites;
they depend on which fibers/samples are included. `m6a_site_n` and `cpg_site_n`
replace the former reference-based `at_n` and `cpg_n` columns.

No reference FASTA is read or required, and there is no `--ref` option. Chromosome
lengths come from BAM headers. The shared high-level matrix assembler itself
selects longest alignments, so Step 3 reuses its lower-level BED parsers and
observed-coordinate definition while preserving primary-only eligibility.

There is no pooled-footprint audit or sample repair stage. The local
`fp_ambiguous` flag detects a name in multiple manifest samples' overlapping
m6A alignments (including supplementary and non-spanning alignments). Such
rows retain non-footprint features but have missing footprint measurements.
This check cannot detect contributors absent from the manifest. Identical
pooled footprint rows are deduplicated. Unindexed pooled gzip files are streamed
once per chromosome, retaining only records overlapping the requested windows;
indexed files use tabix. Segments and pooled footprints do not
carry BAM flags; they are associated by name and coordinates within the primary
span, so overlapping supplementary-derived features cannot be distinguished.

Run one chromosome, optionally restricting genes for a smoke test:

```bash
bash code/pausing/dominant_tss/run_fiber_features.sh --chrom chr22 \
  --genes /path/to/pause_index_dominant.tsv --gene-ids ENSG00000100320
```

Use `--help` to list all path/window options. Submit all 23 chromosomes with:

```bash
sbatch code/pausing/dominant_tss/run_fiber_features.sh
```

The wrapper follows the repository's R 4.4.1 library settings. Required packages
are data.table, GenomicRanges, Rsamtools, GenomicAlignments, and dplyr
(used by the existing BED12 converter). No genome-wide jobs are submitted by
the script itself.

Integration tests (synthetic indexed BAM/BED inputs, both gene strands):

```bash
Rscript code/pausing/dominant_tss/test_fiber_features.R
```
