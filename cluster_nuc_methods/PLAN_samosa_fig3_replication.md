# SAMOSA Figure-3 replication: oligonucleosome-pattern clustering on macrophage Fiber-seq (LPS_0)

## Context

Replicate the core analysis of Abdulhay & McNally et al. (eLife 2020, SAMOSA): per-molecule autocorrelograms of methylation signal → unsupervised Leiden clustering into "oligonucleosome patterns" (regular vs irregular, graded by NRL) → per-molecule nucleosome-repeat-length (NRL) estimation via autocorrelogram peak calling → Figure-3-style characterization panels. Substrate: THP-1 macrophage Fiber-seq LPS_0 `ft extract` outputs under `macrophage_project/FiberHMM/extract/ft_result_dir/LPS_0/extracted_results/`.

**User decisions (confirmed):** use the **m6A signal** (`m6a_by_chr`) as the per-molecule substrate — closest analog to SAMOSA's modification probability and avoids circularity with `ft add-nucleosomes`' HMM-enforced nucleosome-scale footprints (the nuc-occupancy version stays available as a `--signal nuc` variant for a QC concordance check). Scope: **LPS_0, Figure 3 only**; the build step is parameterized by sample so LPS_5/10/15 are trivial reruns later.

## Key data facts (verified)

- `m6a_by_chr/LPS_0.ft_extracted_m6a.<chr>.bed.gz`: BED12, **one row per aligned read** (col4 = read ID, col5 = MAPQ), interior blocks = 1 bp m6A calls (~29/kb). **First block (0-length) and last block (1 bp) are sentinels — strip them** (same convention as `nuc_by_chr`; documented in `code/topic_model/read_region_extraction/region_data_utils.R:44-53`).
- ~18M reads total for LPS_0; median span 8.2 kb, p10 2.5 kb; ~0.5% of reads have multiple rows per chrom (supplementary alignments) → dedupe keep-longest-span per RID.
- `replicate_by_chr/LPS_0.rid_replicate.<chr>.tsv.gz`: RID → R1/R2/R3, for the per-replicate reproducibility panel.
- Independent ground truth: median start-to-start nucleosome spacing in `nuc_by_chr` is **188 bp** — the mean ACF secondary peak must land near this.
- No haplotype info in these files (fine for Fig 3).

## Method adaptation (paper → our data)

- Paper smooths modification probability with a 33 bp rolling mean and autocorrelates the first 1000 bp from the 5' MNase cut. Fiber-seq reads have **no MNase anchor** and are ~8 kb, so: compute the **full-read ACF evaluated at lags 0–1000 bp** (one autocorrelogram per read, span ≥ 2000 bp). Full-read averaging (~7k products per lag) rescues the sparsity of binary m6A calls.
- ACF spec per molecule: binary vector (1 at m6A bases, sentinels stripped) → `np.convolve(ones(33)/33, mode="valid")` → mean-center → autocorrelation via one zero-padded `rfft/irfft` → unbiased normalize by `(N−k)·var`. Drop and count zero-variance (no-m6A) reads. Store lags 0–1000 as float32; cluster on lags 1–1000.
- Filters: MAPQ ≥ 20, span ≥ 2000 bp, drop chrM at clustering (nucleosome-free mtDNA is a trivial cluster).
- Scale: deterministic subsample keep read iff `zlib.crc32(rid) % 60 == 0` → **~300k molecules** genome-wide (reproducible, no RNG; even a 0.05% pattern still yields ~150 molecules > the 100-molecule floor). ACF matrix 300k × 1000 float32 ≈ 1.2 GB.

## Files to create

Code in `/project/spott/cshan/fiber-seq/code/cluster_nuc_methods/` (empty dir, intended home). Results in `/project/spott/cshan/fiber-seq/macrophage_project/cluster_nuc_methods/results/{acf_by_chr,tables,plots}`.

1. **`build_acf_matrix.py`** — per-chrom worker. CLI: `--sample LPS_0 --chrom chr1 --signal m6a|nuc --out-dir`. Streams the bed.gz in pandas chunks (`usecols=[0,1,2,3,4,5,10,11]`, chunksize 50k), applies MAPQ + crc32 hash filter **before** splitting block strings (only ~1/60 of rows get parsed — keeps chr1 tractable), strips sentinel blocks, vectorizes, computes ACF, dedupes keep-longest per RID, joins replicate labels. Writes one `.npz` per chrom: `acf` (n×1001 f32), `sig2k` (n×2000 f16 smoothed signal over first 2 kb, NaN-padded — feeds heatmap panels without refetch), metadata arrays (rid, chrom, start, end, strand, mapq, n_events, span, replicate). Parsing idioms ported from `code/enhancer_accessibility/fiberseq_footprints.py` and `process_nuc_pos.R:109-133`.
2. **`run_build_acf.sbatch`** — array driver, `--array=1-25` → chr1..22,X,Y,M; `SAMPLE=${1:-LPS_0}`. Header cloned from `macrophage_project/FiberHMM/extract/code/ft_extract_nuc.sh` (`--account=pi-spott --partition=caslake`, logs `/project/spott/cshan/fiber-seq/results/logs/build_acf_%A_%a.out`). `--cpus-per-task=2 --mem=16G --time=4:00:00`. Interpreter: `/project/spott/cshan/envs/Jupyter-notebook/bin/python3` (scanpy 1.11.5, leidenalg, scipy verified installed).
3. **`cluster_acf.py`** — single job: concat 25 npz, drop chrM, global dedupe; `AnnData(acf[:,1:1001])` → `sc.tl.pca(n_comps=50)` → `sc.pp.neighbors(n_neighbors=15, metric="correlation")` → `sc.tl.leiden(resolution=0.5, random_state=0)` (pattern from `code/archive/PolII_footprints_code/cluster_footprints_notebooks/cluster_by_footprint_class.ipynb`; retune resolution 0.3–0.8 until every cluster ≥ 100 molecules and 4–10 clusters total). Per-molecule NRL: `scipy.signal.find_peaks(acf[100:401], distance=100, prominence=0.005)`, NRL = 100 + first peak, NaN + per-cluster failure tally on no peak; calibrate prominence once against cluster-mean ACFs (peaks should fall 170–200 bp — paper thresholds can't be copied since signal amplitudes differ). Outputs to `results/tables/`: `acf_clusters.LPS_0.tsv.gz` (per-read assignments + NRL), `cluster_summary.LPS_0.tsv`, `cluster_profiles.LPS_0.npz` (mean ACF/sig2k per cluster + 100 sampled read indices per cluster). `--mem=32G --cpus-per-task=8`.
4. **`plot_clusters.ipynb`** — login-node-safe (reads only small tables/npz). Panels: (3A) cluster abundance bars per replicate R1–R3; (3B) heatmap of per-cluster mean ACF over lags 0–1000 **and** per-cluster mean smoothed m6A over first 1000 bp (caption: anchor = alignment start, not MNase cut); (3C/D) violin + histogram of per-molecule NRL per cluster with failure rates; (3E) single-molecule heatmaps, 100 molecules/cluster over first 2 kb.

Constants at the top of each file (SMOOTH_BP=33, MIN_SPAN=2000, MAX_LAG=1000, MAPQ_MIN=20, HASH_MOD=60, PCA=50, KNN=15, LEIDEN_RES=0.5, SEED=0, NRL_RANGE=100–400, MIN_CLUSTER=100) — no config layer, per CLAUDE.md simplicity.

## Execution order & verification

1. **Smoke test**: run `build_acf_matrix.py` on chr21 only (~4k sampled reads) and push it through `cluster_acf.py` end-to-end before submitting the array — validates sentinel stripping, ACF shape, memory. Verify: mean ACF of chr21 shows a secondary peak near 185–190 bp.
2. Submit `sbatch run_build_acf.sbatch LPS_0`; confirm 25 npz files, row counts ≈ chrom read count / 60.
3. Run `cluster_acf.py` via sbatch. Gates: (a) dataset-wide mean ACF secondary peak at 180–190 bp (must agree with the independently measured 188 bp start-to-start median — strongest correctness check); (b) median per-molecule NRL 185–195 bp; (c) no cluster < 100 molecules, 4–10 clusters; (d) per-replicate cluster fractions agree (print max deviation).
4. Make plots; visually confirm regular (banded) vs irregular (unstructured) single-molecule heatmaps, and that NRL-failure rate is higher in irregular clusters (the paper's own signature).
5. QC concordance: rerun chr21 with `--signal nuc`, report ARI between m6a- and nuc-based cluster labels once (flags HMM-artifact clusters).

## Known caveats (state in outputs, don't engineer around)

- Cluster identities won't map 1:1 to the paper's seven (different cell state, chemistry, caller); the deliverable is the method + the regular/irregular, NRL-graded taxonomy.
- Full-read ACF averages within-read heterogeneity; tiling reads into pseudo-molecules is the named v2 if needed.
- Read-start-anchored positional panels are texture summaries, not boundary-anchored biology; quantitative claims ride on the ACF panels.

## Future (named only)

Fig-6-style chromatin-state enrichment (assignments table already carries coordinates → bedtools intersect); cross-timepoint comparison via `sbatch run_build_acf.sbatch LPS_{5,10,15}`.
