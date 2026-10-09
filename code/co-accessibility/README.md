# Same-fiber cCRE co-accessibility

The two datasets have separate scripts. Both retain cCRE pairs sharing canonical
TSS +/-10 kb windows, with strict `500 < interval gap < 20000` on chr1-22/X/Y.
A same-read FIRE element must cover **at least 50% of the cCRE length**. The
denominator uses primary alignment spans overlapping both cCREs (`any` by default).
Different FIRE elements are not combined to reach 50%; the overlap is not reciprocal.

FIRE peaks select the candidate cCRE universe by >=1 bp overlap, retaining cCRE
coordinates. Per-read FIRE elements provide accessibility calls. These are different
steps: the 50% rule belongs to the per-read accessibility call.

| Dataset | Code | Results | Testing |
|---|---|---|---|
| LCL | [LCL/](LCL/) | `/project/spott/cshan/fiber-seq/LCL_project/co-accessibility/` | One pooled analysis of all 31 samples; donor metadata only |
| Macrophage | [macrophage/](macrophage/) | `/project/spott/cshan/fiber-seq/macrophage_project/co-accessibility/` | Separately for LPS_0/5/10/15 |

The shared annotation builder remains [make_gencode_v46_all_tss.sh](make_gencode_v46_all_tss.sh).
The macrophage scripts and detailed documentation moved to [macrophage/README.md](macrophage/README.md).
Old flat script/notebook paths must be replaced with these dataset subfolders.

## LCL inputs

The source of sample identities and donor metadata is:

`/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv`

All 31 samples are required, including four without genotypes. The 27-sample
phased subset and genotype VCF are not used. Inputs are read from each CSV
`fire_dir`, normally `/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results/<sample>/`:

- `<sample>-fire-v0.1-peaks.bed.gz`: full peak coordinates in columns 1-3.
- `<sample>-fire-v0.1-filtered.cram` plus `.crai`: primary alignment spans.
- `additional-outputs-v0.1/fire-peaks/<sample>-v0.1-fire-elements.bed.gz` plus `.tbi`: per-read FIRE calls.
- `extracted_results/{m6a,nuc}_by_chr/`: display-only tracks, restricted to tested read IDs.

Reference: `/project/spott/reference/human/GRCh38.p14/hg38.fa`.
Annotations: `/project/spott/cshan/annotations/GRCh38-cCREs.bed` and
`/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed`.

## LCL workflow

| Script | Action |
|---|---|
| [LCL_fire_universe.sh](LCL/LCL_fire_universe.sh) | Validate the 31-sample manifest; build genome-wide FIRE peak union, cCRE/gene universe, and sample peak metadata |
| [LCL_read_spans.sh](LCL/LCL_read_spans.sh) | SLURM array, one task per sample, at most six concurrent; extract/index primary alignment BED6 spans |
| [LCL_coaccess_cres.py](LCL/LCL_coaccess_cres.py) | Count per sample/chromosome, pool raw counts across all31, test once per unique pair, correct genome-wide |
| [LCL_run_coaccess.sh](LCL/LCL_run_coaccess.sh) | SLURM launcher for the standalone LCL Python script |
| [LCL_plot_coaccess.sh](LCL/LCL_plot_coaccess.sh) | Render the top-five notebook after the analysis job succeeds |
| [LCL_co-access.Rmd](LCL/LCL_co-access.Rmd) | Read completed results, select up to five significant positive pairs, and verify plotted read states against saved counts |

```bash
cd /project/spott/cshan/fiber-seq/code/co-accessibility/LCL
mkdir -p /project/spott/cshan/fiber-seq/results/logs
bash LCL_fire_universe.sh
span_job=$(sbatch --parsable LCL_read_spans.sh)
analysis_job=$(sbatch --parsable --dependency="afterok:${span_job}" LCL_run_coaccess.sh)
sbatch --dependency="afterok:${analysis_job}" LCL_plot_coaccess.sh
```

`LCL_COACCESS_ROOT` can set a separate output root and `LCL_SAMPLE_METATABLE`
can supply another manifest with exactly 31 unique sample names. Use the same
values across stages. There is no timepoint field and no per-donor coverage floor,
phasing restriction, stratification, weighting, or per-donor hypothesis test.
Read names are resolved within each sample before pooling counts, so names cannot
accidentally join fibers between samples.

## Statistical outputs

The four raw cells are `co_closed`, `CRE1_access`, `CRE2_access`, and `co_access`.
They are summed across all samples **before** adding one pseudocount per cell.
As in the macrophage script, `pval` is two-sided Fisher on the table+1;
`fisher_estimate` is the conditional odds estimate and `OR` is the corrected
cross-product odds ratio. `pval_raw` tests unmodified counts; `fdr_raw` and
`or_haldane` are also retained. The corrected and raw tests are distinct outputs.

BH correction is applied once across all unique tested pairs, not separately by
chromosome or repeated gene membership. Expected pooled co-accessibility is
`(co_access + CRE1_access) * (co_access + CRE2_access) / n_shared_reads`.
Donor counts are bookkeeping only, not the margins of a stratified null model.

Outputs under `LCL_project/co-accessibility/`:

- `universe/sample_manifest.tsv`, `sample_inputs.tsv`: all31 input identities and paths.
- `universe/{fire_peaks_union.bed,cre_universe.bed,gene_windows.bed,cre_gene_map.tsv.gz,cre_in_sample_peaks.tsv.gz}`.
- `universe/universe.complete`: published last; interrupted universe builds cannot be used for counting.
- `<sample>/<sample>.read_spans.bed.gz`, `.tbi`, `.source.tsv`.
- `coaccess/LCL_coaccess_pairs.tsv.gz`: one row per unique pair with both FDRs, effects and rule metadata.
- `coaccess/LCL_coaccess_stat.tsv.gz`: gene-oriented results; a pair may occur for multiple genes.
- `coaccess/LCL_pair_sample_counts.tsv.gz`: raw per-sample cells, including zeros, with donor labels.
- `coaccess/LCL_count_qc.tsv.gz`, `sample_manifest.tsv`, `run_info.json`: provenance and completion marker.
- `coaccess/plots/top5/`: selected pairs, fiber metadata, read plots and configuration bars.

The notebook selects `fdr < 0.05`, `OR > 1`, and `co_access > expected_co_access`,
then orders by FDR, descending OR, descending shared coverage and pair ID. It
shows fewer than five if necessary. Display downsampling never changes counts
or tests. Recompute older outputs made with >=1 bp FIRE overlap; the notebook
requires matching `fire_overlap_fraction=0.5` and `read_rule=any` metadata.

The plotting steps are directly in `LCL/LCL_co-access.Rmd`, in separate chunks
for loading fibers, verifying counts, exporting metadata, preparing display
tracks, drawing fiber panels and configuration bars, and saving figures. It
sources `parsing_footprints_functions.r` and `plotting_functions.r` from
`code/parsing_functions/` and calls their reusable functions at each step.
Run the chunks in order or knit the notebook; `LCL_plot_coaccess.sh` remains an
optional SLURM launcher. The m6A profile and rasters use the displayed subset
(at most 120 fibers per pair by default); counts and bars use all shared fibers.

`--chrom chr21` is an explicit smoke-test mode, written by default under
`smoke/chr21/`. Its FDR is chromosome-scoped and it is not a genome-wide result.
The notebook rejects partial analyses unless `allow_partial_results=TRUE` is
explicitly supplied. Full analysis has no gene or chromosome restriction.
