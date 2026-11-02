# Dominant TSS

CAGE-supported TSS per gene, written to
`/project/spott/cshan/annotations/fantom5_LCL/fantom5.hg38.LCL.dominant_TSS.bed.gz` (+ `.tbi`).

| File | Role |
| --- | --- |
| `build_ctss_consensus.R` | FANTOM5 GM12878 replicate CTSS files -> `fantom5.hg38.LCL.consensus.TSS.bed.gz` |
| `write_dominant_tss.R` | `select_dominant_tss()` (the rules below), run on GENCODE v46; writes the BED |
| `run_write_dominant_tss.sh` | wrapper with the R environment (bash on a login node is fine) |

the replicates are in `/project/spott/cshan/annotations/fantom5_LCL/fantom5.hg38.LCL.rep1.CTSS.bed.gz`
  - 5th column is `score` = raw tag count = the # of CAGE reads whose 5' end falls on this base


```bash
bash code/pausing/dominant_tss/run_write_dominant_tss.sh
```

The input is the hg38 FANTOM5 GM12878 consensus: single-base sites detected in all three
replicates, score = summed tag count. 

Sites are matched to GENCODE v46 transcript TSSs (chr1-22, chrX) on the same strand within 250 bp. Each site is assigned to its nearest gene; equal-distance ties between genes are excluded. 

For each gene, the highest-count assigned site is the dominant TSS. Ties prefer smaller annotation distance, a canonical transcript, lower genomic coordinate, then transcript ID. The TES comes from the matched transcript; matches with the TES upstream of the proposed anchor are excluded. Genes with no eligible CAGE site are omitted (counted in the log). 

