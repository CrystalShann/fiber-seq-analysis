# Dominant TSS

`select_dominant_tss.R` defines the CAGE dominant-TSS selection and GRO-seq
pause-peak functions sourced by `code/pausing/pause_index.R`.

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

Selection and pause-peak tests: `Rscript code/pausing/dominant_tss/test_dominant_tss.R`
(run from the repository root with the wrapper's R library environment).
