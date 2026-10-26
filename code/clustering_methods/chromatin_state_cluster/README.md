# Shared promoter chromatin subclusters

New standalone R workflow. Existing clustering and parsing code is sourced without modification. No new helper functions are defined. Code has been written but not executed or submitted.

The unit is one `(gene_id, sample_name, rid)` row from the original promoter table. All samples and genes contribute to one shared clustering within each fixed broad state: `closed`, `open_fire`, `open_linker`, and `other`. A read overlapping two gene windows contributes one row to each gene. Highly covered genes/samples therefore contribute more observations; there is no balancing or discovery subsampling.

## Inputs

Defaults in `chromatin_state_cluster.R`:

- `--feature-dir`: `LCL_project/pausing/dominant_tss`, containing `fiber_features.{chr}.tsv.gz` from Step 3.
- `--genes`: `LCL_project/pausing/dominant_tss/pause_index_dominant.tsv`, with `gene_id`, `chrom`, `strand`, `anchor_tss` (1-based).
- `--samples`: the shared 31-sample manifest; `sample_name` and `fire_dir` locate `{fire_dir}/extracted_results/{sample_name}.fire_all.bed.gz` and its tabix index.
- `--out-dir`: `LCL_project/pausing/dominant_tss/chromatin_state_cluster`.

The promoter table is authoritative for primary-alignment eligibility and full-window coverage. This workflow checks its BED alignment spans but does not reread BAM flags or select a longest alignment. FIRE records are joined by sample, genomic window, and read ID; the FIRE format cannot distinguish multiple alignments of the same read at the same locus.

The default gene table is an upstream Step 2 output, not a bundled file. Supply
`--genes` with the same anchor table used to generate the promoter features.
The two-sample smoke test used
`/scratch/midway3/cshan/tmp/dominant-fibers-smoke-gwvgubpk/genes.tsv`, together with
that directory's `samples.csv` and `observed_sites/` feature directory. This
two-gene fixture is not a substitute for the full gene-anchor table.

## Stages

1. **prepare** runs per chromosome. Existing `read_tabix_region()` and `occupancy_matrix()` build four occupancy channels: nucleosomes at least 90 bp long, shorter nucleosomes, FIRE, and linker. Original segment lengths are classified before clipping. Profiles are oriented by gene strand over inclusive −1000 through +2000 bp and averaged into 25-bp bins (120 full bins plus one 1-bp bin). No reference genome is needed.
2. **cluster** pools all selected chromosomes, genes and samples within each broad state. Manhattan distances use bin-width weights, `width / (2 * 3001)` per channel, without z-scaling. An inline RcppAnnoy Manhattan index constructs a union KNN graph (`k=10`, 50 trees, `search_k=5000`); exponential edge affinities feed the existing `leiden_partition()`. Closest mean profiles are merged, or heterogeneous groups split using deterministic distant seeds, to obtain exactly three subclusters where at least three distinct profiles exist. These are adjusted Leiden partitions. Original Leiden assignments are retained separately.
3. **report** joins assignments back to all original rows, preserving original columns and broad states. It writes composition and feature summaries and plots for interpreting the subclusters. This is a clustering report, not the pausing regression/Step 4.

The broad state is checked using the original priority within −10 to +10 bp: qualifying nucleosome → `closed`; otherwise FIRE → `open_fire`; otherwise linker → `open_linker`; otherwise `other`.

All eligible rows remain in the output. Clustering requires full, unambiguous segmentation: missing segments, uncovered bases, overlaps between different channels, or a disagreement with the stored broad state produce a QC status and no subcluster. Same-channel overlaps are unioned. Missing segmentation is never treated as a zero-valued closed profile. A broad state with fewer than three distinct usable profiles receives `insufficient_variation`, rather than arbitrary groups.

Subclusters are named by descending size, with a stable read-key tie break, e.g. `closed_S1`–`closed_S3`. Labels are specific to a fitted run. Biological names should follow inspection of their occupancy profiles and feature summaries. m6A, CpG, pausing, expression, sample identity and gene identity are not distance features.

## Invocation templates (not executed)

The default submission prepares one chromosome per array task: tasks 1–22
are chr1–chr22, task 23 is chrX, and task 24 is chrY. Each task requests
12 hours (previously 48 hours) and 128 GB. These limits are allocations, not
measured runtime estimates.

From this directory, submit preparation and then the pooled clustering/report
job, replacing PREPARE_JOB_ID with the ID returned by the first command:

```bash
sbatch run_chromatin_state_cluster.sh
sbatch --dependency=afterok:PREPARE_JOB_ID run_pooled_chromatin_state_cluster.sh
```

The pooled job requests 24 hours and 128 GB. It fits each broad state using all
24 chromosomes together, then writes the report. It has no array directive.
Neither script submits jobs automatically. Do not submit `report` directly
with the preparation script's default array. The explicit `all` mode remains
available through `bash` without an array, but is not the default.

Both jobs use default input/output paths. The full gene-anchor table and all
24 chromosome feature files must already exist; directory creation does not
generate those inputs. Pass `--genes` and `--feature-dir` to preparation for
custom input paths, and pass the same `--out-dir` and profile settings to both
jobs when overriding them. Full-data jobs have not been submitted.

A restricted chromosome run must use the same selection at every stage, for example `--chrom chr22`. Pass any custom window/bin/nucleosome settings consistently at all stages. `--chroms chr1,chr2` selects multiple chromosomes for cluster/report. Inputs are checked by path, size and modification time before reporting; changing prepared inputs requires reclustering.

Required packages include data.table, GenomicRanges, Rsamtools, igraph, RcppAnnoy, and dependencies loaded by the existing shared plotting sources. The wrapper uses the existing R 4.4.1 installation and project library paths. No packages are installed automatically.

## Outputs

- `profiles/`: prepared matrices, original rows plus internal identifiers and QC statuses, and per-chromosome QC counts.
- `states/`: per-state assignments, graph, model/settings, mean occupancy profiles, merge/split log, and optional neighbor-validation statistics.
- `tables/fiber_features.clustered.{chr}.tsv.gz`: original promoter table plus `leiden_community`, `promoter_subcluster`, and `clustering_status`.
- `tables/composition_by_{sample,gene,gene_sample}.tsv`: counts and fractions. Denominators include retained QC failures; `fraction_within_broad` uses all rows of the same broad state in that grouping.
- `tables/subcluster_counts.tsv`, `subcluster_feature_summary.tsv`, `clustering_qc.tsv`, and `session_info.txt`.
- `plots/`: mean occupancy profiles, a bounded display sample of fiber profiles, and representative fibers nearest each cluster's mean profile. Representatives are not claimed to be medoids.

The optional `--validate-size 2000` diagnostic compares approximate and exact neighbors on a separate bounded subset, including tie-tolerant distance-radius agreement and partition ARI. It does not measure full-population neighbor recall or change fitting membership. Set it to `0` to disable. `--plot-per-cluster 100` controls display sampling only. Summary methylation fractions use the promoter table's pooled observed-site denominators, not reference-derived A/T or CpG counts. No pooled-footprint audit or pausing regression is included.
