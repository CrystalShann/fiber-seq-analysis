# Same-fiber cCRE co-accessibility

This analysis looks for ENCODE cCRE pairs that are accessible on the same fiber more often than expected from each cCRE's accessibility. LCL pools counts across 31 samples. Macrophage calculates co-accessibility separately for each LPS timepoint.

## LCL

### Analysis overview

1. Merge FIRE peaks from all 31 samples -> union FIRE peaks on chr1–22, X and Y.
2. Define the cCRE universe -> keep cCREs overlapping the union by at least 1 bp, keeping the cCRE coordinates.
3. Link cCREs to genes -> build TSS ± 10 kb windows around GENCODE v46 `Ensembl_canonical` TSSs and assign cCREs by overlap.
4. Extract each sample's primary aligned blocks from the filtered FIRE CRAM -> BED12 files. Deletions and skipped regions are excluded from aligned coverage.
5. Build candidate cCRE pairs -> both cCREs must share a gene window, with 500 < interval gap < 20,000 bp.
6. Keep fibers whose aligned blocks fully cover both cCREs (`--read-rule contain`).
7. Call accessibility at each cCRE -> one FIRE element on the same fiber must cover at least 50% of the cCRE. Only the part of that element overlapping the fiber's aligned blocks contributes; separate FIRE elements are not combined.
8. Count four configurations for each pair: neither accessible, cCRE1 only, cCRE2 only, both accessible. Sum the raw counts across all 31 samples, keeping read identities separate by sample.
9. Add one pseudocount to each cell, calculate Fisher's exact p-value and the odds ratio, and apply BH correction across unique pairs genome-wide. Raw-count p-values and their FDR are also saved.
10. Use [LCL_co-access.Rmd](LCL/LCL_co-access.Rmd) to select up to ten pairs with FDR < 0.05, OR > 1, and observed co-accessibility above the expectation under independence, then draw the single-fiber figures.

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

Output root: `/project/spott/cshan/fiber-seq/LCL_project/co-accessibility/`.

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
cd /project/spott/cshan/fiber-seq/code/co-accessibility/LCL
mkdir -p /project/spott/cshan/fiber-seq/results/logs

# 1. Build the cCRE universe and sample manifest
bash LCL_fire_universe.sh

# 2. Extract primary aligned blocks for all 31 samples
span_job=$(sbatch --parsable LCL_read_spans.sh)

# 3. Pool counts and calculate co-accessibility
sbatch --dependency="afterok:${span_job}" LCL_run_coaccess.sh \
    --read-rule contain --fire-overlap-fraction 0.5
```

`LCL_run_coaccess.sh` runs [LCL_coaccess_cres.py](LCL/LCL_coaccess_cres.py). After it finishes, run the chunks in `LCL_co-access.Rmd` in order or knit the notebook.

To change the output root or sample metadata, set `LCL_COACCESS_ROOT` or `LCL_SAMPLE_METATABLE` before running the shell scripts. Use the same settings for all stages and update the notebook parameters to match. The manifest must contain 31 unique samples.

## Macrophage

### Analysis overview

1. Merge FIRE peaks from `LPS_0`, `LPS_5`, `LPS_10` and `LPS_15` -> one shared peak union on chr1–22, X and Y. Each timepoint contains merged replicates.
2. Define the cCRE universe -> keep cCREs overlapping the union by at least 1 bp, keeping the cCRE coordinates. Record peak membership for each timepoint.
3. Link cCREs to genes -> build TSS ± 10 kb windows around canonical TSSs and assign cCREs by overlap.
4. Extract primary read spans from each timepoint's filtered FIRE CRAM -> indexed BED files. Secondary and supplementary alignments are excluded.
5. Build candidate cCRE pairs -> both cCREs must share a gene window, with 500 < interval gap < 20,000 bp.
6. Keep fibers whose read spans overlap both cCREs (`--read-rule any`). Macrophage uses continuous read spans for coverage.
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
    --read-rule any --fire-overlap-fraction 0.5
```

`run_coaccess.sh` runs [03_coaccess_cres.py](macrophage/03_coaccess_cres.py). After it finishes, run the chunks in `coaccess_macrophage.Rmd` in order or knit the notebook. Keep `read_rule` and `fire_overlap_fraction` in the notebook matched to the counting settings.

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
