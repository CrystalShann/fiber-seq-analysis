# Promoter shape from LCL CAGE

Shape index (SI) and peaked / broad classes for the promoters of the 100 expressed PolII example
genes (housekeeping and cell type specific in `code/pausing/PolII/example_gene_list.tsv`), from nucleotide-resolution FANTOM5 GM12878 CAGE.
Method after Hoskins et al. 2011 (Genome Res 21:182, doi 10.1101/gr.112466.110) and Schor et
al. 2017 (Nat Genet 49:550, doi 10.1038/ng.3791). Published steps and the adaptations made here
are listed in [Published methods vs adaptations](#published-methods-vs-adaptations).

Scope: only the example genes. The 50 not_expressed genes were removed on request (2026-11-13;
`GROUPS` in `01_promoter_shape.R`): 46 had no CAGE dominant TSS and almost none had tags. Fiber-seq clusters, archetypes and recurrent families are not used.

## Run

```bash
bash code/promoter_shape/run_promoter_shape.sh        # 01: all tables (~1 min, login node is fine)
/project/spott/cshan/envs/puffin/bin/python code/promoter_shape/03_puffin_dominant_motif_site.py   # ~2 min, after 01
```

Then open `02_promoter_shape_figures.Rmd` in RStudio and run its chunks (or render it): each chunk
shows its figure and writes the same PDF to `figures/`.

The wrapper runs `/software/R-4.4.1-el8-x86_64/bin/Rscript` with the project `R_LIBS_USER`
(data.table, GenomicRanges, ggplot2, patchwork). Avoid the `amd` partition (SIGILL).

| File | Role |
| --- | --- |
| `01_promoter_shape.R` | paths, parameters and helpers (input readers, promoter matching, SI / width / shift / bootstrap), then `run_promoter_shape()` for Steps 1-8, 10, 11. `source()` gives the helpers only; Rscript runs the pipeline |
| `02_promoter_shape_figures.Rmd` | Step 9 figures 01-07, the Puffin motif section (tables and figures 08-09) and the representative genes (section 9) |
| `03_puffin_dominant_motif_site.py` | reruns Puffin at dominant_TSS to place each gene's dominant motif (only TATA / YY1 positions are in 04a), for section 9 |
| `run_promoter_shape.sh` | R environment wrapper |

## Inputs

All inputs already existed in the project; nothing was downloaded for this analysis. Inputs are
read only.

| Input | Path | Source |
| --- | --- | --- |
| Example genes | `code/pausing/PolII/example_gene_list.tsv` | PolII curated list; the 50 housekeeping and 50 cell type specific genes are used |
| dominant_TSS | `annotations/fantom5_LCL/fantom5.hg38.LCL.dominant_TSS.bed.gz` | `code/pausing/dominant_tss/write_dominant_tss.R` |
| dominant_TSS fallback | `annotations/gencodev46_Ensembl_canonical_TSS.bed` | GENCODE v46 canonical TSS (1-bp BED) |
| CTSS, 3 replicates | `annotations/fantom5_LCL/fantom5.hg38.LCL.rep{1,2,3}.CTSS.bed.gz` | FANTOM5 hg38 reprocessed GM12878 hCAGE CNhs12331 / CNhs12332 / CNhs12333; download links in `annotations/fantom5_LCL/README.md` |
| Promoter intervals | `annotations/fantom5_LCL/fantom5.hg38.CAGE_peaks.bed.gz` | FANTOM5 `hg38_fair+new_CAGE_peaks_phase1and2.bed.gz` |
| Peak-to-gene | `annotations/fantom5_LCL/hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt.gz` + `hgnc_complete_set.txt` | FANTOM5 hg38 annotation (HGNC IDs) -> Ensembl gene IDs via HGNC |
| Peak names | `annotations/fantom5_LCL/human_phase1and2_CAGE_Peak_name.txt.gz` | FANTOM5 `p<n>@<gene>` names (display only) |

`annotations/` = `/project/spott/cshan/annotations`.

Data notes:

- **dominant_TSS is not recalculated.** It is looked up exactly as `resolve_example_tss()` in
  `code/pausing/PolII/pol2_pausing_functions.r` does: the FANTOM5 LCL dominant TSS by gene name
  (99 of the 100 analysed genes), else the GENCODE v46 canonical TSS (FCGR2B only, which has no
  promoter). TSS = BED end (1-based). `dominant_TSS_source` records which; the test checks every
  value against `LCL_project/pausing/PolII/recurrent_patterns/tables/example_genes.tsv`.
- **hg38.** The FANTOM5 hg38 readme (fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/00_readme.txt)
  states "target genome: hg38 (the raw set of hg38 sequences, except for alt, random, Un.)". Spot
  check: the top GAPDH CTSS is chr12:6,534,517 (+), 5 bp from the GENCODE hg38 TSS.
- **Read filtering is upstream and cannot be redone here** (the CTSS files are counts, not reads).
  Same readme: "5'-end of the CAGE read alignments with mapping quality above 20 (which is
  equivalent to single mapped reads) and percent identity 85% are counted at 1bp resolution".
  Multi-mapped reads are therefore already excluded, which also removes tags at multi-copy loci
  (FCGR2B has no FANTOM5 peak; HLA genes keep only uniquely mapped tags).
- **Samples.** The three GM12878 libraries are the only "B lymphoblastoid" libraries in the FANTOM5
  hg38 cell-line set. FANTOM5 also has one library of RPMI 1788 (CNhs10744, labelled "b cell line";
  Cellosaurus CVCL_2710, EBV-transformed, male). It is not used: a single library from another
  donor, while dominant_TSS was built from GM12878.

## Method

Coordinates are 1-based inclusive throughout (`start` = FANTOM5 BED start + 1, `end` = BED end;
TSS and CTSS positions = BED end). Offsets and shifts are taken in the direction of
transcription: negative = upstream, positive = downstream, on both strands.

### Step 1. Promoter regions

Candidates for each gene:

1. same-strand FANTOM5 peaks that the hg38 FANTOM5 annotation assigns to the gene (HGNC ID ->
   Ensembl gene ID);
2. the peak that contains the gene's dominant_TSS on its strand, whoever FANTOM5 assigns it to
   (same-strand FANTOM5 peaks never overlap, so there is at most one).

Kept as promoters (one row each, user choice):

- the dominant_TSS peak always, also without tags (`contains_dominant_TSS = TRUE`; shape then
  unavailable);
- every other candidate with >= 1 tag in each of the three replicates (the consensus rule of
  `build_ctss_consensus.R`, applied per peak).

Candidates on another chromosome or strand are dropped (one: the chr19 ERF peak that FANTOM5 also
assigns to ETF1). Ambiguous matches seen in the 2026-11-13 run (the qc tables that listed them
were removed on request the same day):

- dominant_TSS peak not assigned to the gene by FANTOM5: CFLAR (`p7@CFLAR`, no HGNC link in the
  hg38 annotation), ENSG00000285218 (in `p1@CLDN11`);
- peaks shared with other genes: HLA-DRB1 (with HLA-DRB3/4/5), TNFRSF13B (with TBC1D27P);
- genes with several promoters: all kept, the dominant_TSS peak flagged;
- genes with no promoter (dominant_TSS outside every FANTOM5 peak and no LCL-active annotated
  peak): FCGR2B. It has no row in the main table.

Result of the 2026-11-13 run: 518 promoters for 99 genes (99 dominant_TSS promoters, all passing
QC, and 419 alternative promoters).

### Steps 2-3. CTSS profiles and aggregation

- Each replicate's tag count is taken at every nucleotide of every promoter on the promoter's
  strand; positions without a CTSS are 0. Replicate profiles are written before aggregation.
- Depth: 10.33 M, 9.45 M and 10.00 M tags (largest / smallest = 1.09). Below `MAX_DEPTH_RATIO`
  (1.5) the raw counts are summed, so the aggregate keeps integer tag counts for the 20-tag rule
  and the multinomial bootstrap. Above it, 01 stops and asks for normalisation.
- Comparability (checked on 2026-11-13, no longer recomputed): pairwise Pearson r of
  log2(CPM + 1) was 0.994-0.995 for promoter totals and 0.978-0.979 per nucleotide; each
  replicate's own class matched the pooled class for 98.9% of promoters.
- Profiles are not smoothed and no positional offset is applied (single assay).

### Steps 4-8. Features

With c_i the aggregated count at nucleotide i of the promoter and p_i = c_i / sum_j c_j:

| Feature | Definition |
| --- | --- |
| `total_CAGE_tags` | sum of c_i |
| `n_active_CTSS` | positions with c_i > 0 |
| `shape_index` | SI = 2 + sum over positions with c_i > 0 of p_i log2 p_i (2 = one nucleotide; -1 = eight equal nucleotides); NA without tags |
| `strongest_CAGE_TSS` | position of max c_i; ties go to the tied position nearest dominant_TSS, then the 5'-most |
| `dominant_TSS_shift` | strongest_CAGE_TSS - dominant_TSS, sign flipped on the minus strand (> 0 = downstream) |
| `dominant_TSS_fraction` | max c_i / total (fraction at strongest_CAGE_TSS) |
| `TSS_width_10_90` | positions ordered 5' -> 3'; q10 / q90 = first position where the cumulative count reaches 10% / 90%; width = abs(q90 - q10) bp (0 when both fall on one nucleotide) |
| `Hoskins_shape_class` | peaked if SI > -1, broad if SI <= -1; unclassified unless `QC_status` = pass |
| `Schor_shape_class` | narrow if SI > -1, broad if SI <= -1; unclassified unless `QC_status` = pass |
| `classification_stability` | fraction of bootstrap resamples on the same side of SI = -1 as the observed SI; NA below 20 tags |
| `QC_status` | `no_CAGE_tags` (shape unavailable), `low_coverage` (1-19 tags), `unstable_classification` (stability < 0.9), else `pass` |

Both class columns use the published SI = -1 cut (user choice, 2026-11-13; the task text had SI > 1
for Hoskins and -1.5 for Schor, which do not match the papers). They differ only in their labels,
so their proportions and stability are identical. SI is reported for low-coverage promoters but
they stay unclassified.

### Step 7. Classification stability

For each promoter with >= 20 tags: 500 multinomial resamples, each drawing the observed tag total
from the observed per-nucleotide proportions (nonparametric bootstrap); SI of each resample;
stability = fraction on the same side of -1 as the observed SI; stable if >= 0.9. The seed of a
promoter is 20261113 + its FANTOM5 peak number (hg_<number>.1), so a promoter's result does not
depend on which other promoters are analysed. RNG: R 4.4.1 default (Mersenne-Twister / Inversion /
Rejection).

## Parameters

| Parameter | Value | Published or adapted |
| --- | --- | --- |
| `SI_CUT` | -1 | published (Hoskins Methods; Schor Supplementary Table 7) |
| `MIN_TAGS` | 20 | follows Hoskins: "20 tags are required, on average, to confidently infer promoter class" |
| `N_RESAMPLE` | 500 | adaptation (Hoskins' subsampling settings are in their Supplemental Methods, not accessible) |
| resampling depth | observed total tags | adaptation |
| `SEED` | 20261113 + peak number | adaptation |
| `STABILITY_MIN` | 0.9 | adaptation |
| `MAX_DEPTH_RATIO` | 1.5 | adaptation (observed 1.09: raw counts summed) |

## Puffin motifs and promoter shape

Section 8 of the Rmd joins the Puffin scores of the example genes
(`LCL_project/pausing/PolII/promoter_annotation/tables/04a_puffin_motif_scores.tsv`, scored at
dominant_TSS - 1, checked) to the dominant_TSS promoters passing QC. Motif contributions are summed
over strands per family. Each feature is tested against SI by Spearman rho and by
rank(SI) ~ rank(feature) + gene group + log2 TPM + log10 CAGE tags, on all genes and on genes whose
dominant Puffin motif is not TATA; BH over features. Classes are compared by Kruskal-Wallis
(classes with >= 5 genes, with and without TATA) and TATA class vs peaked by Fisher's test. Puffin
is trained on FANTOM CAGE, so its scores are not independent of CAGE-derived shape.

## Representative genes

Section 9 of the Rmd shows 11 genes picked by agreement between the Puffin dominant motif and CAGE
shape (peaked: PPIA, POU2AF1, ACTG1, NFKBIA; broad: SEC61A1, PSMB1, RAN, CD9; discordant: SRSF3,
DYNLL1, RPS6), with their table first. Per gene:

- `<gene>_CAGE_peaks.pdf`: dominant_TSS -150..+150 bp, per-replicate CTSS (scaled to each
  replicate's maximum), FANTOM5 decomposed peaks as grey bands and a numbered track coloured by
  Hoskins class with SI, and the dominant Puffin motif (largest contribution at the TSS, U1 snRNP
  excluded) as a yellow centre +/- 4 bp band. Its centre is the peak activation of that motif from
  `03_puffin_dominant_motif_site.py`: TATA+ in -60..-10 and YY1+ in 0..+50 (04a's ranges, so equal
  to `tata_center` / `yy1_center`), all other motifs in -100..+50.
- `<gene>_read_level.pdf`: `plot_read_level()` of `code/pausing/PolII/pol2_pausing_functions.r`,
  as in the TATA section of `summarize_recurrent_patterns.Rmd`, with the Puffin TATA box. Each gene
  reads all 31 samples (minutes per gene); existing PDFs are reused unless `REDRAW_MOLECULES` is TRUE.

## Published methods vs adaptations

| Item | Published | Here |
| --- | --- | --- |
| SI formula | Hoskins Methods: SI = 2 + sum_i^L p_i log2 p_i, "L is the set of base positions that have at least one TSS tag" | same |
| Classes | Hoskins Methods: "Promoter regions with shape index score >-1 were classified as peaked (P); all others were classified as broad (B)". Schor Supplementary Table 7: "shape: whether the peak is broad or narrow (threshold is -1)" | same cut, both columns |
| Unclassified | Hoskins: "unclassified due to low tag count (2487) or class-instability (982) as determined by subsampling from the existing TSS distribution"; test details in their Supplemental Methods | adapted: < 20 tags, or bootstrap stability < 0.9 (500 resamples at observed depth) |
| Promoter units | Hoskins: own iterative hierarchical clustering of tags into promoter regions (CAGE peaks < 50 bp apart merged), consensus over CAGE, RLM-RACE and cap-trapped ESTs | adapted: FANTOM5 hg38 DPI peaks, kept separate (user instruction) |
| Assays | Hoskins: three assays, each smoothed to absorb their 1-2 bp relative offsets; Schor: CAGE across 81 *Drosophila* lines | adapted: CAGE only, 3 GM12878 replicates summed, no smoothing or offset correction |
| Read filtering | FANTOM5: MAPQ > 20 (single-mapped), identity >= 85% | inherited from FANTOM5, nothing further |
| strongest_CAGE_TSS, dominant_TSS_shift, dominant_TSS_fraction, n_active_CTSS | not in either paper | added for this project |
| TSS_width_10_90 | not in either paper; follows the CAGE 10-90% interquantile width (CAGEr reports q90 - q10 + 1) | distance q90 - q10, as requested |

## Outputs

`/project/spott/cshan/fiber-seq/LCL_project/promoter_shape/`

| File | Content |
| --- | --- |
| `promoter_shape_annotations.tsv` | one row per gene x promoter, columns as in Steps 4-8 plus `gene_id`, `gene_name`, `gene_group`, `promoter_id` (FANTOM5 hg38 peak ID), `chr`, `start`, `end`, `strand`, `dominant_TSS`, `dominant_TSS_source`, `contains_dominant_TSS` |
| `CTSS_profiles/sample_CTSS_profiles.tsv.gz` | every nucleotide x replicate: `pos`, `offset_from_dominant_TSS`, raw `count`, `cpm` (per genome-wide library size) |
| `CTSS_profiles/aggregated_CTSS_profiles.tsv.gz` | every nucleotide: summed `count`, `fraction` (p_i) |
| `figures/01-09_*.pdf` | figures (PDF only); 08-09 = Puffin motifs vs SI; `05` picks, per class, the 4 passing dominant_TSS promoters closest to the class median SI among those with at least the class median of tags |
| `puffin_dominant_motif_sites.tsv` | per gene: `dominant_motif`, `dominant_contribution`, `motif_center` (bp from dominant_TSS), `motif_activation` |
| `figures/representative_genes/<gene>/` | `<gene>_CAGE_peaks.pdf` and `<gene>_read_level.pdf` (section 9) |

## Caveats

- **DPI peaks split broad initiation regions.** 63 of the 99 dominant_TSS promoters have a kept
  alternative promoter of the same gene within 50 bp (70 within 300 bp). For these 99 promoters, the dominant peak holds a median 95% (IQR 85-100%) of the tags within
  300 bp, but SI per peak leans toward peaked, and `TSS_width_10_90` cannot exceed the peak
  (dominant_TSS peaks are 7-112 bp, median 43).
- **Alternative promoters are mixed.** They include gene-body CAGE peaks of highly expressed genes
  (e.g. 18 for ACTG1, 31 for CFLAR) and distant internal TSSs (up to 381 kb from dominant_TSS).
  For one value per gene use `contains_dominant_TSS == TRUE`.
- **SI is biased upward at low depth**, and the bootstrap inherits the bias: the resampled median
  SI exceeds the observed SI for 99% of promoters (median +0.21 at 20-50 tags, +0.01 above 1,000).
  Near the cut, promoters just below -1 are flagged unstable more often than those just above
  (within 0.25 of -1: 24 of 32 below vs 5 of 42 above in the 2026-11-13 run).
