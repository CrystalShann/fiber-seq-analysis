# Same-fiber cCRE co-accessibility

This analysis looks for ENCODE cCRE pairs that are accessible on the same fiber more often than expected from each cCRE's accessibility. LCL pools counts across 31 samples. Macrophage calculates co-accessibility separately for each LPS timepoint.

## LCL

### Analysis overview

1. Merge FIRE peaks from all 31 samples -> union FIRE peaks on chr1–22, X and Y.
2. Define the cCRE universe -> keep cCREs overlapping the union by at least 1 bp, keeping the cCRE coordinates.
3. Link cCREs to genes -> build TSS ± 10 kb windows around GENCODE v46 `Ensembl_canonical` TSSs and assign cCREs by overlap.
4. Extract each sample's primary aligned blocks from the filtered FIRE CRAM -> BED12 files. Deletions and skipped regions are excluded from aligned coverage.
5. Build candidate cCRE pairs -> both cCREs must share a gene window, with 500 < interval gap < 20,000 bp.
6. Keep fibers whose aligned blocks fully cover both cCREs.
7. Call accessibility at each cCRE -> one FIRE element on the same fiber must cover at least 50% of the cCRE. Only the part of that element overlapping the fiber's aligned blocks contributes; separate FIRE elements are not combined.
8. Count four configurations for each pair: neither accessible, cCRE1 only, cCRE2 only, both accessible. Sum the raw counts across all 31 samples, keeping read identities separate by sample.
9. Add one pseudocount to each cell, calculate Fisher's exact p-value and the odds ratio, and apply BH correction across unique pairs genome-wide. Raw-count p-values and their FDR are also saved.
10. Use [LCL_co-access.Rmd](LCL/CRE_co-accessibility_FIRE/LCL_co-access.Rmd) to select up to ten pairs with FDR < 0.05, OR > 1, and observed co-accessibility above the expectation under independence, then draw the single-fiber figures.

### Input

Sample metadata:

`/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv`

All 31 samples are included, including samples without genotypes. The CSV's `fire_dir` gives each sample's FIRE directory, normally under `/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results/`.

| Input | Used for |
|---|---|
| `<fire_dir>/<sample>-fire-v0.1-peaks.bed.gz` | Union FIRE peaks, using columns 1–3 |
| `<fire_dir>/<sample>-fire-v0.1-filtered.cram` and `.crai` | Primary aligned blocks |
| `<fire_dir>/additional-outputs-v0.1/fire-peaks/<sample>-v0.1-fire-elements.bed.gz` and `.tbi` | Accessibility calls on each fiber |
| `/project/spott/cshan/annotations/GRCh38-cCREs.bed` | cCRE coordinates and classes |
| `/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed` | Gene TSS anchors |
| `/project/spott/reference/human/GRCh38.p14/hg38.fa` | Reference for reading CRAM files |

The notebook also reads display tracks from each sample's `extracted_results/` and the FiberHMM results under `/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/FiberHMM_v2/results/results_v1_trained_model`. These tracks provide m6A, nucleosome and footprint annotations for the figures. Accessibility counts come from the FIRE elements and aligned blocks.

### Output

Output root: `/project/spott/cshan/fiber-seq/LCL_project/co-accessibility/CRE_co-accessibility_FIRE/`.

| Step | Output relative to the output root | Contents |
|---|---|---|
| `LCL_fire_universe.sh` | `universe/sample_manifest.tsv`, `universe/sample_inputs.tsv` | Sample identities, donor labels and input paths |
| `LCL_fire_universe.sh` | `universe/fire_peaks_union.bed`, `universe/cre_universe.bed` | Merged FIRE peaks and selected cCREs |
| `LCL_fire_universe.sh` | `universe/gene_windows.bed`, `universe/cre_gene_map.tsv.gz` | TSS windows and cCRE-to-gene membership |
| `LCL_fire_universe.sh` | `universe/cre_in_sample_peaks.tsv.gz`, `universe/universe.complete` | Per-sample peak membership and completion marker |
| `LCL_read_spans.sh` | `<sample>/<sample>.aligned_blocks.bed.gz`, `.tbi`, `.source.tsv` | Indexed primary aligned blocks and source information |
| `LCL_run_coaccess.sh` | `coaccess/LCL_coaccess_pairs.tsv.gz` | One row per unique pair, with pooled counts, odds ratios, p-values and FDR |
| `LCL_run_coaccess.sh` | `coaccess/LCL_coaccess_stat.tsv.gz` | Gene-anchored results; pairs can appear under multiple genes |
| `LCL_run_coaccess.sh` | `coaccess/LCL_pair_sample_counts.tsv.gz` | Four configuration counts per sample and pair |
| `LCL_run_coaccess.sh` | `coaccess/LCL_count_qc.tsv.gz`, `coaccess/sample_manifest.tsv` | Counting summary and input provenance |
| `LCL_co-access.Rmd` | `coaccess/all_chr/tables/`, `coaccess/all_chr/plots/` | Selected pairs, fiber labels, count summaries and figures |

### Run the analysis

Run the universe script first, then submit the aligned-block array. The counting job waits for all 31 array tasks to finish successfully.

```bash
cd /project/spott/cshan/fiber-seq/code/co-accessibility/LCL/CRE_co-accessibility_FIRE
mkdir -p /project/spott/cshan/fiber-seq/results/logs

# 1. Build the cCRE universe and sample manifest
bash LCL_fire_universe.sh

# 2. Extract primary aligned blocks for all 31 samples
span_job=$(sbatch --parsable LCL_read_spans.sh)

# 3. Pool counts and calculate co-accessibility
sbatch --dependency="afterok:${span_job}" LCL_run_coaccess.sh \
    --fire-overlap-fraction 0.5
```

`LCL_run_coaccess.sh` runs [LCL_coaccess_cres.py](LCL/CRE_co-accessibility_FIRE/LCL_coaccess_cres.py). After it finishes, run the chunks in `LCL_co-access.Rmd` in order or knit the notebook.

To change the output root or sample metadata, set `LCL_COACCESS_ROOT` or `LCL_SAMPLE_METATABLE` before running the shell scripts. Use the same settings for all stages and update the notebook parameters to match. The manifest must contain 31 unique samples.

### m6A permutation test on the top ten pairs

[LCL_co-access_v2_perm_methy_prop.Rmd](LCL/CRE_co-accessibility_FIRE/LCL_co-access_v2_perm_methy_prop.Rmd) extends the analysis in `LCL_co-access.Rmd`. I used Fisher's exact test to test whether cCRE pairs are co-accessible more often than expected, but there is high donor heterogeneity: a pair accessible in one donor may not be accessible in another, and mixing donors can create an association on its own.

This notebook asks a within-sample question instead: within a sample, does a fiber with more m6A at CRE1 also have more m6A at CRE2? It uses continuous m6A proportions rather than a binary open/closed call, and the permutation test shuffles only within each sample.

Two sets of pairs are tested:

- `pooled`: the top ten pairs from the 31-sample Fisher test, using all 31 samples
- `noXCI`: the top ten pairs after removing female chrX. chrX pairs use male samples only, and all other chromosomes use all 31 samples

#### Analysis overview

1. Load the top ten pairs (ranked by FDR) for each set, and each pair's per-sample counts from `LCL_pair_sample_counts.tsv.gz`.
2. Reload each pair's fibers with `load_region()` and label them with `label_reads()`, using the same settings as the counting (`read_rule = "contain"`, `fire_overlap_fraction = 0.5`). Keep the fibers covering both cCREs. The per-sample counts of the four states must match the saved counts exactly, so the m6A test uses the same fibers as the Fisher test.
3. Calculate the m6A proportion of every fiber at each cCRE with `fiber_m6a_proportions()` from `parsing_footprints_functions.r`:
   - for each sample, query the sample's m6A calls over the cCRE, keeping the longest alignment of each read
   - a read is kept only if it covers the whole cCRE
   - `assemble_region_m6a()` builds one reads × positions matrix from the full-span reads of all included samples: 1 = m6A call, 0 = no call. The columns are the positions with at least one m6A call on any of those reads
   - proportion = the fiber's m6A calls ÷ the number of columns. Every fiber at a cCRE has the same denominator (`n_positions1`, `n_positions2`), and it is not the cCRE width or the A/T count
4. Keep the `qualifying` fibers: full-span at both cCREs, with finite proportions. Each fiber now has `prop1` (CRE1) and `prop2` (CRE2).
5. Run the permutation test on each pair (below), then apply BH correction across the ten pairs within each set.
6. Summarize the mean proportions per pair, per donor and per FIRE state, and draw the figures.

#### How the permutation test works

The test is `permute_pooled_spearman()` in [stats_testing_functions.r](../parsing_functions/stats_testing_functions.r). The test statistic is the Spearman correlation (rho) between `prop1` and `prop2` across the fibers of a pair. It measures how consistently the two values rise and fall together, using their ranks:

- rho = +1: fibers with higher CRE1 m6A always have higher CRE2 m6A
- rho = 0: no relationship
- rho = −1: higher CRE1 always goes with lower CRE2

Observed rho

1. Rank `prop1` and `prop2` across all fibers of the pair; ties get the average of their ranks.
2. Center each rank: c = rank − mean rank.
3. rho = Σ(c1 × c2) / √(Σc1² × Σc2²). A fiber above average at both cCREs (+ × +) or below average at both (− × −) adds a positive product. A fiber high at one and low at the other adds a negative product, so a large positive sum means the two tend to agree.

Two versions of rho are reported:

- `spearman_rho` (pooled): ranks are centered on the mean rank of all fibers. This includes differences between samples, for example a sample that is more methylated at both cCREs
- `within_sample_adjusted_rho`: ranks are centered on each sample's own mean rank, so only the variation among fibers within a sample counts

Null rho

1. Within each sample, randomly reassign the CRE2 ranks among that sample's fibers. CRE1 and the sample labels stay fixed, so each sample keeps its own CRE2 values and only which fiber they belong to changes. Samples whose CRE2 proportion does not vary are not shuffled.
2. Recompute both rhos, and repeat 10,000 times to get 10,000 null rhos.
3. Calculate the one-sided p-value for a positive association with `empirical_pvalue()`:
   `p = (1 + #{null rho ≥ observed rho}) / (1 + 10,000)`. The smallest possible p-value is 1/10,001 ≈ 1e-4.

Notes

- The shuffle keeps each sample's mean rank, so the pooled null still contains the between-sample part of the correlation. A pooled null centered well above 0 means samples differ systematically, with some high at both cCREs and some low at both. The adjusted null is centered near 0.
- For the same reason, `p_positive` and `adjusted_p_positive` are always identical: under every shuffle, the pooled statistic equals the adjusted statistic plus a constant set by the sample mean ranks, rescaled. There is one test, and the two rhos are effect sizes. If they are close, the correlation comes from fibers within samples. If the adjusted rho is much lower, part of the pooled trend comes from samples that are more methylated at both cCREs.
- Each pair has its own seed, `seed + pair_rank` (+1000 for `noXCI`), so a pair's result does not depend on which pairs are run together.
- BH correction (`adjust_bh()`) uses the ten pairs within each set as the family.
- The shuffled rhos are also converted to nominal one-sided Spearman p-values (`spearman_nominal_logp()`), with t = r·√((n−2)/(1−r²)) and n − 2 degrees of freedom. The histograms and QQ plots of these p-values are diagnostics only; the reported p-values are the permutation p-values.

#### Input

The `coaccess/`, `universe/` and `<sample>/` paths are under the LCL output root `/project/spott/cshan/fiber-seq/LCL_project/co-accessibility/CRE_co-accessibility_FIRE/`.

| Input | Used for |
|---|---|
| `coaccess/LCL_top10_pairs.tsv`, `coaccess/LCL_top10_pairs_noXCI.tsv` | The two top-ten sets, saved by `LCL_co-access.Rmd` |
| `coaccess/LCL_pair_sample_counts.tsv.gz` | Per-sample configuration counts, used to check the reloaded fibers |
| `coaccess/noXCI/LCL_sample_sex_noXCI.tsv` | Sex of each sample, written by `LCL_coaccess_noXCI.py`; selects the male samples for noXCI chrX pairs |
| `universe/sample_manifest.tsv`, `coaccess/sample_manifest.tsv` | Sample identities, donors and input paths; the two manifests must agree |
| `<sample>/<sample>.aligned_blocks.bed.gz` and `.tbi` | Fiber coverage of each cCRE |
| `<fire_dir>/additional-outputs-v0.1/fire-peaks/<sample>-v0.1-fire-elements.bed.gz` and `.tbi` | FIRE accessibility of each fiber |
| `<fire_dir>/extracted_results/m6a_by_chr/<sample>.ft_extracted_m6a.<chr>.bed.gz` and `.tbi` | m6A calls on each fiber (`ft extract` BED12) |
| `/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv` | Donor of each sample |

`AL6_bc2162_19093`'s chr1 m6A file is plain gzip rather than BGZF, so tabix returns no reads from it. None of the current top-ten pairs are on chr1; for chr1 pairs those fibers would appear as `missing_ft_record` in the exclusion tables.

#### Output

Output root: `/project/spott/cshan/fiber-seq/LCL_project/co-accessibility/CRE_co-accessibility_FIRE/co-access_v2_perm_methy_prop/`. Files already on disk are kept and only missing files are written, so delete a file to regenerate it. `<set>` is `pooled` or `noXCI`.

| Output relative to the output root | Contents |
|---|---|
| `tables/<set>_top10_original_fisher_statistics.tsv` | The selected pairs with their original Fisher statistics |
| `tables/per_fiber_m6a_proportions.tsv` | One row per fiber per pair: donor, FIRE state, m6A calls, denominators and proportions at CRE1 and CRE2, and `qualifying` |
| `tables/pooled_m6a_means.tsv` | Mean `prop1` and `prop2` for each pair, with the number of fibers and the denominators |
| `tables/donor_m6a_means.tsv`, `tables/fire_state_m6a_means.tsv` | The same means by donor and by FIRE state |
| `tables/donor_accessibility_observed_expected.tsv` | Observed and expected co-accessibility for each donor |
| `tables/m6a_exclusions_by_sample.tsv`, `tables/m6a_exclusions_by_pair.tsv` | Fibers lost between the Fisher counts and the m6A test, and why |
| `tables/pooled_m6a_permutation_statistics.tsv` | Both rhos, permutation p-values, BH p-values, null mean and SD, and status for each pair |
| `plots/<set>_donor_accessibility.pdf` | FIRE configuration proportions by donor |
| `plots/<set>_per_fiber_m6a_scatter.pdf` | `prop1` against `prop2` for every fiber, coloured by FIRE state |
| `plots/<set>_permutation_null.pdf` | The 10,000 null rhos with the observed rho (red line), pooled and adjusted |
| `plots/<set>_fire_state_m6a_means.pdf` | Mean m6A proportion at each cCRE by FIRE state |
| `plots/<set>_shuffled_p_value_histograms.pdf`, `plots/<set>_shuffled_p_value_qq.pdf` | Nominal Spearman p-values of the shuffled datasets |

#### Run the analysis

Run this after `LCL_co-access.Rmd` has saved both top-ten tables. Run the chunks in `LCL_co-access_v2_perm_methy_prop.Rmd` in order or knit the notebook; the twenty pairs take about five minutes.

## Macrophage

### Analysis overview

1. Merge FIRE peaks from `LPS_0`, `LPS_5`, `LPS_10` and `LPS_15` -> one shared peak union on chr1–22, X and Y. Each timepoint contains merged replicates.
2. Define the cCRE universe -> keep cCREs overlapping the union by at least 1 bp, keeping the cCRE coordinates. Record peak membership for each timepoint.
3. Link cCREs to genes -> build TSS ± 10 kb windows around canonical TSSs and assign cCREs by overlap.
4. Extract primary read spans from each timepoint's filtered FIRE CRAM -> indexed BED files. Secondary and supplementary alignments are excluded.
5. Build candidate cCRE pairs -> both cCREs must share a gene window, with 500 < interval gap < 20,000 bp.
6. Keep fibers whose read span covers the full length of both cCREs. Macrophage uses continuous read spans for coverage.
7. Call accessibility at each cCRE -> one FIRE element on the same fiber must cover at least 50% of the cCRE. Separate FIRE elements are not combined.
8. Count the four accessibility configurations separately at each timepoint. Add one pseudocount per cell, calculate Fisher's exact p-value and the odds ratio, and apply BH correction across unique pairs within each timepoint. Raw-count p-values and their FDR are also saved.
9. Use [coaccess_macrophage.Rmd](macrophage/coaccess_macrophage.Rmd) to select pairs with a CA-CTCF partner and at least one protein-coding gene anchor. The notebook looks for weak coupling at baseline followed by an upward trend in the odds ratio, then draws up to five pairs across all four timepoints.

### Input

FIRE root: `/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE/`.

Here, `<sample>` is `LPS_0`, `LPS_5`, `LPS_10` or `LPS_15`.

| Input | Used for |
|---|---|
| `<FIRE root>/<sample>/<sample>-fire-v0.1-peaks.bed.gz` | Shared FIRE peak union, using columns 1–3 |
| `<FIRE root>/<sample>/<sample>-fire-v0.1-filtered.cram` and `.crai` | Primary read spans |
| `<FIRE root>/<sample>/additional-outputs-v0.1/fire-peaks/<sample>-v0.1-fire-elements.bed.gz` and `.tbi` | Accessibility calls on each fiber |
| `/project/spott/cshan/annotations/GRCh38-cCREs.bed` | cCRE coordinates and classes |
| `/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed` | Gene TSS anchors |
| `/project/spott/reference/human/GRCh38/hg38.fa` | Reference for reading CRAM files |

The notebook also reads m6A from `macrophage_project/FiberHMM/extract/ft_result_dir/`, nucleosomes from `macrophage_project/FiberHMM/extract/firehmm_footprint/`, and size-binned TF footprints from `macrophage_project/FiberHMM/extract/firehmm_tf/ft_by_size/`, all under `/project/spott/cshan/fiber-seq/`. These are display tracks; the configuration counts use the filtered FIRE read spans and FIRE elements.

### Output

Output root: `/project/spott/cshan/fiber-seq/macrophage_project/co-accessibility/`.

| Step | Output relative to the output root | Contents |
|---|---|---|
| `01_make_fire_universe.sh` | `universe/fire_peaks_union.bed`, `universe/cre_universe.bed` | Merged FIRE peaks and selected cCREs |
| `01_make_fire_universe.sh` | `universe/gene_windows.bed`, `universe/cre_gene_map.tsv.gz` | TSS windows and cCRE-to-gene membership |
| `01_make_fire_universe.sh` | `universe/cre_in_timepoint_peaks.tsv.gz` | Peak membership at each timepoint |
| `02_read_spans.sh` | `<sample>/<sample>.read_spans.bed.gz` and `.tbi` | Indexed primary read spans |
| `run_coaccess.sh` | `coaccess/<sample>_coaccess_pairs.tsv.gz` | One row per unique pair at each timepoint, with counts, odds ratios, p-values and FDR |
| `run_coaccess.sh` | `coaccess/<sample>_coaccess_stat.tsv.gz` | Gene-anchored results; pairs can appear under multiple genes |
| `coaccess_macrophage.Rmd` | `coaccess/plots/<pair_id>_reads.pdf` | Single-fiber panels across the four timepoints |
| `coaccess_macrophage.Rmd` | `coaccess/plots/<pair_id>_config_bars.pdf` | Accessibility configuration proportions at each timepoint |

### Run the analysis

Run the universe script once, then submit the four-task read-span array. The counting job waits for all four tasks to finish successfully.

```bash
cd /project/spott/cshan/fiber-seq/code/co-accessibility/macrophage
mkdir -p /project/spott/cshan/fiber-seq/results/logs

# 1. Build the shared cCRE universe
bash 01_make_fire_universe.sh

# 2. Extract primary read spans for all four timepoints
span_job=$(sbatch --parsable 02_read_spans.sh)

# 3. Calculate co-accessibility separately at each timepoint
sbatch --dependency="afterok:${span_job}" run_coaccess.sh \
    --fire-overlap-fraction 0.5
```

`run_coaccess.sh` runs [03_coaccess_cres.py](macrophage/03_coaccess_cres.py). After it finishes, run the chunks in `coaccess_macrophage.Rmd` in order or knit the notebook. Keep `fire_overlap_fraction` in the notebook matched to the counting setting.

## Shared TSS annotation

Both analyses use the same canonical TSS BED. If it is missing, build it before running either universe script:

```bash
cd /project/spott/cshan/fiber-seq/code/co-accessibility
bash make_gencode_v46_all_tss.sh
```

Input: `/project/spott/cshan/annotations/gencode.v46.annotation.gtf.gz`. The script downloads the GENCODE v46 GTF if it is missing.

Outputs under `/project/spott/cshan/annotations/`:

| File | Contents |
|---|---|
| `gencode.v46.annotation_all_tss.bed` | 20 bp intervals around transcript TSSs |
| `TSS_interval_gencodev46_Ensembl_canonical.bed` | 20 bp intervals for `Ensembl_canonical` transcripts |
| `gencodev46_Ensembl_canonical_TSS.bed` | 1 bp canonical TSS anchors used by both analyses |

The shell scripts use the project tool paths already set in each script. Array and counting jobs run through SLURM; their logs are saved under `/project/spott/cshan/fiber-seq/results/logs/`.
