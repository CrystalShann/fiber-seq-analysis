# LCL TSS autocorrelation vs RNA Pol II pausing

Does single-molecule nucleosome patterning at the TSS differ between high- and
low-pausing promoters, at matched expression? Pilot on 6 expression-matched gene pairs
from the 31-sample LCL Fiber-seq cohort.

**Nothing scientific is reimplemented.** The ACF, the repeat-peak finder and the read
extraction are sourced unchanged from the existing code:

| Reused as-is | From |
|---|---|
| `autocorrelations()` | `../03_compute_autocorrelations.py` |
| `repeat_peak()` | `../05_summarize_autocorrelations.py` |
| `assemble_region_m6a()` | `code/parsing_functions/parsing_footprints_functions.r` |

The Python functions are loaded into R with `reticulate::source_python()` and called
directly; verified against `python -B` on a real region matrix at 6.3e-14 max absolute
difference. New code here is only the LCL-specific part: pausing index, gene selection,
and plots.

## Run

```bash
./run_lcl_tss.sh            # stages 00 -> 03 in order
./run_lcl_tss.sh 02 03      # only the named stages
```

The wrapper exists to set the library environment and must be used:

- `LD_LIBRARY_PATH` needs `/project/spott/cshan/envs/dimelo/lib` **prepended**, or
  rtracklayer fails to load with
  `libcurl.so: undefined symbol: curl_easy_option_next`.
- `RETICULATE_PYTHON` must point at an interpreter with numpy/scipy
  (`/project/spott/cshan/envs/Jupyter-notebook/bin/python`).
- Override the interpreter with `LCL_TSS_RSCRIPT` if needed.

Runtime: stage 00 about 1 minute, stage 02 about 1 minute, the others seconds. No SLURM
needed at this scale.

## Stages

| Script | Does | Writes |
|---|---|---|
| `00_pause_index.R` | lifts GM12878 GRO-seq hg19 to hg38, scores the Day et al. PI windows on canonical TSS/TES | `tables/pause_index.tsv`, `tables/liftover_qc.tsv` |
| `01_select_tss_genes.R` | filters, takes the PI tails, matches 1:1 on expression | `tables/tss_regions.tsv` |
| `02_tss_autocorrelation.R` | extracts spanning reads, orients, ACF, paired tests | `outputs/<gene>/acf.rds`, `tables/*.tsv` |
| `03_plot_tss_autocorrelation.R` | 6 figures | `results/plots/*.pdf` |
| `tss_autocorrelation.Rmd` | report, incl. full pausing-index provenance | — |

Outputs go to `/project/spott/cshan/fiber-seq/LCL_project/auto_correlation/lcl_tss/`.

### Figures come from stage 03, not from the Rmd

`03_plot_tss_autocorrelation.R` is what writes the PDFs. The Rmd only *embeds* PDFs that
already exist and performs no computation, so knitting it would not produce them — and
it should not be knitted at all. Verify the Rmd with `knitr::purl()` then `Rscript` on
the extracted chunks.

To refresh figures after a code change: `./run_lcl_tss.sh 03`
(or `./run_lcl_tss.sh 02 03` if the tables also need rebuilding).

## Where the pausing index comes from

The definition and the data are from two different papers. This trips people up, so it
is spelled out in the Rmd as well.

- **Definition** (TSSR = −50..+300, body = +300..TES+3kb, PI = density ratio):
  Day et al. 2016, *Genome Biology* 17:120, PMID 27259512, PMC4893286.
  Computed there from **RNAP2 ChIP-seq** on hg19/RefSeq, and it publishes **no per-gene
  PI table** (Additional file 2 = datasets, 3 = GO of top/bottom 25%, 4 = primers), so
  PI is recomputed here.
- **Nascent data**: GEO `GSM1480326` (GSE60454/GSE60456), Core et al. 2014,
  *Nature Genetics* 46:1311-20, PMID 25383968. GM12878 **GRO-seq**, hg19, one replicate.
  Not from Day et al.
- **Expression**: ENCODE `ENCSR820PHH`, GM12878 total RNA-seq, GRCh38, files
  `ENCFF345SHY` and `ENCFF873VWU`. Older GENCODE IDs, so the TSS join strips the
  `.version` suffix.

Two things that must not be changed casually:

1. **PI needs the canonical TSS, not the Ensembl gene record.** For ACTB the gene record
   sits 33 kb from the real promoter. The TES comes from the *same* canonical transcript.
2. **PI is undefined for silent genes** — both windows are empty, giving 0/Inf/NA
   (measured: RHO 0/3, GFAP 0/16, NEUROD1 0/0, MYOD1 3/0). Non-expressed genes are
   excluded; that is why this is a pausing analysis and not an expression analysis.

## Design notes

- **The gene is the unit of analysis.** Molecules within a gene share a promoter and are
  not independent, so `tables/paired_tests.tsv` (paired across matched pairs) is the
  primary statistic. The pooled group curves in figure 01 are molecule-weighted and
  descriptive only; per-gene molecule counts range 283–669, so pooling weights genes
  unevenly.
- **Expression matching is load-bearing.** PI is weakly *anti*-correlated with
  expression here (Spearman log10 PI vs log10 TPM about −0.27), so unmatched PI tails
  would differ in transcription rate too.
- **The 5–60 TPM band** keeps PI measurable (below ~5 TPM both windows are sparse) and
  avoids a pool dominated by ribosomal-protein promoters.
- **Per-molecule NRL is uncalibrated.** `repeat_peak()` has no prominence constraint and
  single-molecule ACF noise is about 1/sqrt(2000), so structureless molecules still
  yield a peak. An A/T-restricted shuffled null belongs to the genome-wide stage.

## Result of the pilot

No significant paired difference in phasing amplitude (p = 0.625) or repeat length
(p = 0.125) across 6 pairs. Repeat length is the more promising direction: shorter at
high-PI promoters in 4 of 5 evaluable pairs, median −14 bp. The three-fold amplitude gap
visible between the pooled curves is **not** reproduced pair by pair and should not be
reported as an effect of pausing.

## Known shared-data issue, worked around locally

One shared m6A extract is broken: `AL6_bc2162_19093.ft_extracted_m6a.chr1.bed.gz` is
**plain gzip, not BGZF**. It is the only non-BGZF file among all 775 (31 samples x 25
chromosomes); every other one is BGZF. Because it is not BGZF it cannot be
tabix-indexed at all — `tabix -f -p bed` refuses with "the compression is not BGZF" —
and htslib therefore returns nothing for any chr1 query while reporting the failure on
stderr only. The accompanying `.tbi` is the original and is fine; it simply indexes a
layout the file no longer has. The symptom is that the sample silently contributes
**zero** reads to every chr1 region.

The data itself is intact (221,378 records, chr1:13,423-248,945,697), so the repair is a
lossless recompression, not a re-index.

The shared tree is read-only for this analysis, so the corrected file lives under
`LCL_project/fixed_inputs/<sample>/extracted_results/m6a_by_chr/` and stage 02
substitutes it through `assemble_region_m6a()`'s own `m6a_paths` argument. The
decompressed content is byte-identical to the shared original (same md5). Any
sample/chromosome placed in `fixed_inputs/` with a `.tbi` is picked up automatically.

To regenerate it:

```bash
zcat <shared>/AL6_bc2162_19093.ft_extracted_m6a.chr1.bed.gz \
  | bgzip -c > <fixed_inputs>/.../AL6_bc2162_19093.ft_extracted_m6a.chr1.bed.gz
tabix -p bed <fixed_inputs>/.../AL6_bc2162_19093.ft_extracted_m6a.chr1.bed.gz
```

The proper fix is to recompress the shared file in place with `bgzip` and re-index it;
that rewrites shared data, so it is left to whoever owns that tree. Stage 02 keeps the
dropout detector (`tables/sample_contribution.tsv` plus a warning) either way, so a
future instance cannot pass silently.

## Caveats

- Expression and GRO-seq are both GM12878, but the Fiber-seq is 31 other LCL donors.
  A genotype-matched GM12878 Fiber-seq CRAM exists at
  `FIRE/results/GM12878_PM00001_PacBio/` but has no `extracted_results/`, so it would
  need an `ft extract --m6a` pass first.
- Six pairs: a plumbing-and-direction pilot, not a powered test.
- One GRO-seq replicate, so PI has no replicate-level error.
- Day et al.'s windows applied to GRO-seq; a narrower TSS window is common for GRO-seq
  and is worth a sensitivity check.
