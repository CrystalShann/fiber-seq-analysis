# Pooled Pol II footprint-positive clustering

One global Leiden–Manhattan analysis of all 31 LCL samples at fantom5 consensus TSS with that have Pol II signal and a FiberHMM footprint 30-60 bp anywhere long the read

## Selection criteria

Use the dominant FANTOM5 TSSs (chr1–22, X). A TSS anchor must overlap an ENCSR000EAD GRCh38 Pol II IDR peak (`ENCFF886PSD`). Each TSS also gets
`pol2_tss_signal`, the mean `ENCFF328MMS` fold-enrichment over its covered bases at −50..+299 (Day's TSS region, strand-oriented); missing coverage stays NA

- IDR stands for Irreproducible Discovery Rate. It's ENCODE's standard way of picking ChIP-seq peaks that reproduce between replicates (Li et al. 2011). An "IDR peak" is a peak that passed that reproducibility filter

A qualifying footprint is a FiberHMM v2 TF block 30–60 bp long that overlaps the 2 kb TSS window, taken from the molecule's own sample
(`<footprint_dir>/<S>/<S>.recalled_tf.bed.gz`)

For a molecule qualifying at multiple TSSs, choose one eligible TSS at random. Bulk Pol II fold-enrichment (`ENCFF328MMS`) is measured over the footprints of the sampled molecules only; missing bigWig signal remains missing.

The molecule must span the full window, relative positions -1000 through +999. All full-window molecules define the observed-site denominator; only
footprint-positive full-span molecules enter clustering:

```
whole-window proportion = molecule's distinct m6A call count / pooled observed sites
50-bp feature           = molecule's distinct calls in bin / pooled observed sites in bin
```

`bin_read_matrix()`: methylation fraction is calculated in 40 bins of each 50 bp

## How to run

All commands are run from the project root (`/project/spott/cshan/fiber-seq`).
Settings and input paths are in `pol2_chip_cluster_utils/config.yaml`; outputs go to
its `output` folder (`LCL_project/pausing/clustering_Pol2`). Slurm logs are written to
`results/logs/pol2_clusters_<jobid>.out` / `.err`.

`00_run.sh` is the single entry point (caslake, 32 GiB, 24 h). It calls
`01_run_pol2_clustering.R` with:

```
--stage prepare|count|features|cluster|project|tracks|autocorrelation|nucleosome|expression|promoter-motifs|motif-bars|report|all
--config YAML     (default: pol2_chip_cluster_utils/config.yaml)
```

`bash code/pausing/clustering_Pol2/00_run.sh ...` runs the same command on the current
node instead of submitting it.

## All together

```bash
sbatch code/pausing/clustering_Pol2/00_run.sh
```

Runs prepare -> features (all chromosomes) -> cluster -> project -> tracks -> report in one
job. It does not run the optional `count` stage. The report computes the
autocorrelation and nucleosome results itself if their tables do not exist yet.

## Step by step

Run the stages in this order; each stage reads the outputs of the earlier ones
(all paths below are under the output folder).

1. **prepare**: select dominant FANTOM5 TSSs that overlap a Pol II peak.
   Writes `tables/prepared.rds`.

   ```bash
   sbatch code/pausing/clustering_Pol2/00_run.sh --stage prepare
   ```

2. **count** (optional, not used by later stages): count eligible molecules per
   chromosome. Uses `--cpus-per-task` cores. Writes `tables/eligible_molecules.rds`.

   ```bash
   sbatch --cpus-per-task=4 code/pausing/clustering_Pol2/00_run.sh --stage count
   ```

3. **features**: extract the 50-bp m6A features per chromosome. Writes
   `tables/features/<chrom>.rds`. Either all chromosomes in one job, one chromosome,
   or one array task per chromosome (tasks 1–22 = chr1–22, task 23 = chrX):

   ```bash
   sbatch code/pausing/clustering_Pol2/00_run.sh --stage features
   sbatch code/pausing/clustering_Pol2/00_run.sh --stage features --chrom chr1
   sbatch --array=1-23 code/pausing/clustering_Pol2/00_run.sh
   ```

4. **cluster**: one pooled Leiden clustering of the sampled molecules. Needs all 23
   feature files. Writes `tables/clustering.rds`.

   ```bash
   sbatch code/pausing/clustering_Pol2/00_run.sh --stage cluster
   ```

5. **project**: assign the remaining eligible molecules to their nearest cluster.
   Writes `tables/tracks/molecule_projection.tsv.gz`.

   ```bash
   sbatch --cpus-per-task=8 code/pausing/clustering_Pol2/00_run.sh --stage project
   ```

6. **tracks**: collect footprints and nucleosomes of the clustered molecules.
   Writes `tables/tracks/molecule_footprints.tsv.gz`.

   ```bash
   sbatch --cpus-per-task=8 code/pausing/clustering_Pol2/00_run.sh --stage tracks
   ```

7. **report**: write the assignment/summary tables and figures (`plots/`), then render
   `pol2_clustering_report.html` into the output folder.

   ```bash
   sbatch code/pausing/clustering_Pol2/00_run.sh --stage report
   ```

8. **autocorrelation / nucleosome / expression** (optional): recompute one report
   section on its own, e.g. after changing its settings or inputs. Needs a finished
   report (reads `tables/molecule_assignments.tsv.gz`). Run `report` again afterwards
   to show the new results in the HTML.

   ```bash
   sbatch code/pausing/clustering_Pol2/00_run.sh --stage autocorrelation
   sbatch code/pausing/clustering_Pol2/00_run.sh --stage nucleosome
   sbatch code/pausing/clustering_Pol2/00_run.sh --stage expression
   ```

## Rerunning

Saved tables record the inputs and code they were built from. If a stage stops with a
"stale" or "inputs changed" error, rerun from the earliest affected stage onwards.

## Promoter motif examples

`05_motif_scan.py` scans explicit promoter windows using local JASPAR models and
FIMO. `06_pol2_clustering_report.Rmd` produces one PDF per selected promoter in
`plots/promoter_motif_cluster_examples`. Its single motif panel shows the highest
FIMO-score hit per family (TATA, YY1, NFY, SP, ETS, ZNF143, NRF1, CREB), all black.
Yellow shading marks every displayed motif; families without hits are labeled.
The accessibility profiles and individual reads reuse the existing plotting
functions. Nucleosome blocks of 130–160 bp remain included. U1 is omitted.
Support tables count 10–30-bp blocks covering at least half the displayed motif.

```bash
bash code/pausing/clustering_Pol2/00_run.sh --stage promoter-motifs
bash code/pausing/clustering_Pol2/00_run.sh --stage motif-bars
```

The first command scans the seven selected windows and regenerates their PDFs
without rendering the full report. The second runs `07_motif_cluster_bars.R`,
preserving the original 5,000-read cluster bars and their existing data/output
folders named `motif_puffin_scan`. The former code folder of that name is removed.
All raw FIMO hits remain in the tables; selection affects the promoter display
and its support table, not the original cohort cluster bars.

Methylation profiles now use observed sites without smoothing. The display-only
`pol2_observed_profiles.R` adapters recover the original full-window-read site
sets at each TSS and cache them in `tables/observed_site_profiles`. Pooled per-bp
profiles divide calls by reads whose TSS has an observed site at that position;
undefined positions are omitted. The 50-bp clustering features retain their
existing observed-site denominators. Promoter profiles retain motif highlighting
and use filled areas. ACF input stays on the consecutive bp grid; only its
reported methylation fractions use observed sites.

## Prepare observed sites before rendering

Run from the project root. The observed-site stage only prepares missing or stale
chromosome masks; it does not rerun clustering. Four Slurm CPUs match its workers.

```bash
sbatch --cpus-per-task=4 code/pausing/clustering_Pol2/00_run.sh --stage observed-sites
# After that job completes successfully:
sbatch code/pausing/clustering_Pol2/00_run.sh --stage render
```

`render` knits the saved report directly; `report` is the original analysis and
report stage. Setup loads masks only and never starts the full BED extraction.
The cache remains under `tables/observed_site_profiles`; completed chromosomes
are reused. Existing analyses must already be prepared for the other sections.

For one selected promoter, independent of the full mask job:

```bash
bash code/pausing/clustering_Pol2/00_run.sh --stage promoter-motifs --gene ARF4
```

This prepares only ARF4's mask if needed, reuses the saved promoter read cache,
scans its window, and writes its PDF to `plots/promoter_motif_cluster_examples`.
Its supporting tables go into `tables/promoter_motif_cluster_examples/ARF4`,
leaving the seven-promoter tables intact. Replace ARF4 with MCOLN1, PSMC1, GAS2,
FNTB, TPK1 or EIF2S1. Without `--gene`, `promoter-motifs` uses the complete cache.

### Preparation jobs for knitting

FIMO and histone bigWig extraction no longer run during knitting. Submit these
once before rendering (and rerun them after changing their inputs):

```bash
sbatch --job-name=pol2_motif_inputs code/pausing/clustering_Pol2/00_run.sh --stage motif-inputs
sbatch --job-name=pol2_histone_inputs code/pausing/clustering_Pol2/00_run.sh --stage histone-inputs
```

Both stages use `08_prepare_report_inputs.R`. `motif-inputs` calls the standalone
`05_motif_scan.py` and prepares/reuses the promoter read cache. `histone-inputs`
uses the original histone signal function and saves its profile table.
Missing caches now stop rendering with a preparation command, rather than
triggering BED extraction, FIMO, or bigWig queries. For a single example, the
plotting command now also requires a prepared mask: run `observed-sites --gene
ARF4` first if the full mask job has not completed. The shared seven-promoter
motif/read cache can serve an ARF4-only figure.
