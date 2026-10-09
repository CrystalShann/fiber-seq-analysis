# co-accessibility: same-molecule co-accessibility of ENCODE cCRE pairs around gene TSSs

Measures whether pairs of ENCODE cCREs inside a gene window are accessible on the
**same fiber** more often than chance, separately for each LPS timepoint
(`LPS_0`, `LPS_5`, `LPS_10`, `LPS_15` — minutes of stimulation, merged replicates).

## Pipeline

```
make_gencode_v46_all_tss.sh  ->  01_make_fire_universe.sh  ->  02_read_spans.sh
                                       ->  03_coaccess_cres.py  ->  coaccess_examples.Rmd
```

### 0. `make_gencode_v46_all_tss.sh` — unchanged

The only script carried over. Builds 20 bp TSS intervals from the GENCODE v46 GTF,
the `Ensembl_canonical` subset, and the 1 bp canonical TSS bed used as gene anchors:

- `/project/spott/cshan/annotations/gencode.v46.annotation_all_tss.bed` (254,070 transcript TSSs, 20 bp)
- `/project/spott/cshan/annotations/TSS_interval_gencodev46_Ensembl_canonical.bed` (63,086 canonical TSSs, one per gene, 20 bp)
- `/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed` (same 63,086 TSSs as 1 bp rows, start = 0-based TSS)

### 1. `01_make_fire_universe.sh` — the shared element universe

Run once. A few minutes; no `sbatch` needed.

1. **FIRE peak union** — peaks from all four timepoints pooled and merged. Columns
   1–3 of the peaks bed (`peak_start`/`peak_end`); columns 4–5 are a *narrower core*
   interval and are not what Kevin selects. No `pass_coverage` filter, as he applies
   none when screening cCREs.
2. **cCRE universe** — cCREs overlapping a merged peak by ≥1 bp, **cCRE coordinates
   kept**, `CRE_ID = accession1.accession2.CRE_label`.
3. **Gene windows** — TSS ± 10 kb around each canonical TSS, keyed on `gene_id`.
4. **Membership** — cCRE → gene window on any overlap, plus per-timepoint peak
   membership flags so a per-timepoint-peak-restricted view stays a filter rather
   than a different universe.

Outputs in `.../co-accessibility/universe/`:

| file | contents |
|---|---|
| `fire_peaks_union.bed` | 143,070 merged peak intervals |
| `cre_universe.bed` | 146,925 cCREs (86,361 dELS / 29,161 pELS / 20,474 PLS / …) |
| `gene_windows.bed` | 63,049 gene windows |
| `cre_gene_map.tsv.gz` | 142,632 (cCRE, gene) memberships across 40,240 genes; 30,140 carry >1 cCRE |
| `cre_in_timepoint_peaks.tsv.gz` | per-timepoint peak membership per cCRE |

Note the universe is **promoter-enriched by construction** — 39% of PLS survive the
FIRE-peak filter against 4.8% of dELS. That is inherent to Kevin's design, not a bug,
but it shapes what the pair set can contain.

### 2. `02_read_spans.sh` (SLURM array, one task per timepoint)

Genome-wide aligned fiber spans from the FIRE CRAM, `chr1–22, X, Y`:

```
samtools view -T hg38.fa -F 0x900 <cram> <chrom> | bedtools bamtobed
```

`-F 0x900` drops secondary and supplementary alignments — without it the same read
name appears at several loci and the name-keyed join to FIRE elements silently fuses
distinct alignments. Kevin has no equivalent filter. Output is bgzipped and
tabix-indexed, and reused when present.

### 3. `03_coaccess_cres.py`

Pairs are enumerated **within each gene window**: all cCRE pairs sharing a window
with `500 < gap < 20000`. Gap is the GenomicRanges-style distance between intervals
(0 if they overlap or are bookended), matching his `distance()` — not a midpoint
distance.

Per timepoint, a fiber's state at a cCRE:

- **covered** (denominator) — the fiber overlaps the cCRE (`--read-rule any`,
  Kevin's rule) or spans it entirely (`--read-rule contain`). From the fiber spans
  built in step 2 out of the `-filtered` FIRE CRAM.
- **accessible** (numerator) — covered AND one of that fiber's **FIRE elements**
  overlaps the cCRE. From
  `lizarraga_FIRE/<s>/additional-outputs-v0.1/fire-peaks/<s>-v0.1-fire-elements.bed.gz`,
  the FIRE pipeline's own per-read element call at FDR ≤ 0.05 — the same
  information Kevin reads out of `fire.bed`'s FIRE-class segments, already computed
  genome-wide and tabix-indexed.

`FiberHMM/extract/ft_result_dir` is **not** used here. It holds `ft extract`
m6A/CpG/nucleosome calls with no FIRE scoring, from the unfiltered BAM; it feeds the
figures only.

2×2, matching `table(fire_region1, fire_region2)` with levels forced to `c(FALSE, TRUE)`:

```
             elem2 FALSE      elem2 TRUE
elem1 FALSE  co_closed        CRE2_access
elem1 TRUE   CRE1_access      co_access
```

- `pval` — two-sided Fisher exact on the table **+1 in every cell**, exactly
  `fisher.test(contingency_table + pseudocount)`
- `fisher_estimate` — the **conditional MLE** odds ratio `fisher.test` reports
- `OR` — his separate cross-product `((co_access+1)(co_closed+1)) / ((CRE1_access+1)(CRE2_access+1))`.
  This is **not** the same number as `fisher_estimate`; both are his and both are kept.
- `pval_raw`, `or_haldane` — added. The `+1` pseudocount perturbs the null in an
  uncontrolled direction, so **use `pval_raw` for inference**.

Flags: `--read-rule {any,contain}`, `--min-dist`, `--max-dist`, `--pseudocount`,
`--kevin-compat`, `--chrom` for a single-chromosome run.

Outputs in `.../co-accessibility/coaccess/`:

| file | contents |
|---|---|
| `<s>_coaccess_stat.tsv.gz` | **one row per (gene, cCRE pair)** — Kevin's `coaccess_stat_df` columns first, in his order and names |
| `<s>_coaccess_pairs.tsv.gz` | **one row per distinct cCRE pair**, with `n_genes`, `gene_ids`, and BH `fdr` / `fdr_raw` |

`co_closed` is also emitted as `co_inaccess`: his builder uses the former, his
`example_figures.Rmd` filters on the latter. Both names, one quantity — so his
snippets run unchanged:

```r
coaccess_stat_df <- fread("LPS_15_coaccess_stat.tsv.gz")
sig_coaccess_stat_df <- coaccess_stat_df[pval < 0.01 & dist > 200 &
                                         co_access > 10 & co_inaccess > 5]
sig_coaccess_stat_df[gene_name == "NTRK1" & dist == 2436]
```

Use `*_coaccess_stat.tsv.gz` for gene-anchored lookup and `*_coaccess_pairs.tsv.gz`
for anything statistical — the gene-anchored table repeats a pair once per
containing gene.

### 4. `coaccess_examples.Rmd` + `coaccess_plot_functions.R`

Read-level figures for selected pairs, one facet per timepoint, fibers sorted by
configuration at the pair (the counterpart of his `cluster = "fire_configs"`), with
the pair shaded on **every** panel. Four panels:

| panel | draws | from |
|---|---|---|
| `accessibility` | fraction of covering fibers accessible across the window | spans + FIRE elements |
| `cres` | the cCRE track, tested pair marked red | `cre_universe.bed` |
| `m6a` | one row per fiber: backbone + methylated adenines | `ft_result_dir/*/m6a_by_chr` |
| `reads` | one row per fiber: backbone + **nucleosomes** (navy) + **FIRE elements** (orange) | `ft_result_dir/*/nuc_by_chr` + `lizarraga_FIRE` |

m6A and nucleosomes come from the **existing** `ft extract` run under
`macrophage_project/FiberHMM/extract/ft_result_dir` — no `ft` invocation is needed
anywhere in this pipeline. Those files are BED12 with one row per fiber and the
features as blocks; note every row carries a leading **size-0 sentinel block** at
offset 0, which the shared `convert_ft_bed12_to_bed6()` drops.

They are **display only**. The accessibility call the 2×2 is built from is always
the FIRE elements from `lizarraga_FIRE`. `ft_result_dir` has no FIRE scoring at all
(it is `ft extract --m6a/--cpg/--nuc`, not `ft fire`), and it was extracted from the
unfiltered BAM, so it covers ~8% more fibers than the `-filtered` FIRE CRAM —
using it as a denominator would count fibers FIRE itself excluded.

`coaccess_plot_functions.R` is an independent implementation. Kevin's `plots.R`
cannot be reused here regardless: `rids_df` requires a `sample_name` with exactly
three underscore-delimited fields (`LPS_0` has two), `cluster = "fire_configs"`
hard-`stop()`s unless exactly two `cluster_regions` are passed, and the panels his
notebooks call (`pileup_haps`, `dimelo_reads`) no longer exist in the checked-out
function.

| function | purpose |
|---|---|
| `load_region()` | tabix fiber spans + FIRE elements for a window, across timepoints |
| `load_ft_tracks()` | attach m6A + nucleosome blocks from `ft_result_dir` (display only) |
| `read_tabix_region()` + `convert_ft_bed12_to_bed6()` | read and expand a BED12 slice using the shared LCL helpers |
| `label_reads()` | per-fiber configuration at the pair, under the same rules as the table |
| `order_reads()` | sort fibers by timepoint, then configuration, then position |
| `coaccess_region_result()` | adapt the retained, ordered pair reads for the shared region builder |
| `coaccess_m6a_fraction()` | per-position fraction of covering fibers methylated, by timepoint |
| `plot_region_panels()` (shared `parsing_functions/plotting_functions.r`) | the four stacked panels, pair shaded on each |
| `plot_stacked_proportion()` | configuration proportions per timepoint |

## Running it

```bash
cd /project/spott/cshan/fiber-seq/code/co-accessibility

# 0. TSS annotation (once; skip if the canonical TSS bed exists)
bash make_gencode_v46_all_tss.sh

# 1. shared cCRE universe (once)
bash 01_make_fire_universe.sh

# 2. genome-wide fiber spans, four timepoints in parallel (~25 min each)
sbatch 02_read_spans.sh

# 3. co-accessibility + Fisher tests
sbatch run_coaccess.sh
#    single-chromosome smoke test:
#    python3 03_coaccess_cres.py --chrom chr21 --out-dir /tmp/chr21

# 4. read-level figures
module load R/4.4.1 pandoc/2.17.1.1
Rscript -e 'rmarkdown::render("coaccess_examples.Rmd")'
