# FIRE_MSP: global same-molecule codependency of FIRE peaks across the LPS timecourse

Measures whether pairs of FIRE peaks tend to be accessible **on the same
fiber** more often than expected if they were independent, genome-wide
(autosomes, all genes), separately for each LPS timepoint (`LPS_0`, `LPS_5`,
`LPS_10`, `LPS_15` — minutes of LPS stimulation, merged replicates). Each
timepoint is tested by itself; timepoints are shown side by side but never
tested against a baseline timepoint.

## Pipeline

```
intersect_msp_fire.sh   ->  fire_codependency.py (run_fire_codependency.sh)  ->  plot_codependency.R
```

### 1. `intersect_msp_fire.sh` (SLURM array, one task per timepoint)

Inputs: the merged FIRE results per timepoint under
`/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE/<s>/`
(`<s>-fire-v0.1-filtered.cram`, `<s>-fire-v0.1-peaks.bed.gz`).

Steps, writing to `/project/spott/cshan/fiber-seq/macrophage_project/FIRE_MSP/<s>/`:

1. **Peak filter** — FIRE peaks with `pass_coverage == true`, autosomes only,
   score = fire_coverage/coverage → `<s>.FIRE_peaks_pass_coverage.autosomes.bed.gz`.
   The **LPS_0** peak set is the shared peak universe every timepoint is
   intersected against, so pair identities match across the timecourse (this is
   a shared coordinate set, not a baseline comparison).
2. **MSP extraction** — `ft extract -x "len(msp)>150"` from the FIRE CRAM,
   BED12 → one row per MSP (BED6) → `<s>.MSP_BED6_single_150.bed.gz`.
3. **MSP × peak intersect** — MSPs overlapping an LPS_0 peak by ≥50% (of either
   feature) → `<s>.MSP_150_FIRE_intersect.vsLPS0peaks.bed.gz` (numerator).
4. **Read × peak intersect** — primary alignments **fully containing** an LPS_0
   peak (`bedtools intersect -F 1`, secondary/supplementary dropped)
   → `<s>.FIRE_read_intersect.vsLPS0peaks.bed.gz` (denominator).

Run LPS_0 first (`bash intersect_msp_fire.sh LPS_0`, builds the universe), then
`sbatch intersect_msp_fire.sh` for the array.

### 2. `fire_codependency.py` (via `sbatch run_fire_codependency.sh`)

Per timepoint, a fiber's state at a peak:

- **covered** — the read fully contains the peak
- **actuated** — covered AND an MSP ≥50%-overlaps the peak

All peaks are included by default (`--min-act` 0); setting `--min-act` higher
drops peaks whose maximum actuation across timepoints is below the threshold.
For every peak pair co-spanned by one fiber (same chromosome), the four
molecule states are counted: `[neither, first, second, both]`. There is no
coverage minimum by default (`--min-cov` 0); a pair is scored as long as each
peak is actuated on at least one of the co-spanning fibers (pairs where either
peak is never actuated among the shared fibers are dropped):

- `score = (obs_both − p1·p2) × 4`, with p1/p2/obs from the co-spanning fibers
  (>0 means the peaks open together more than expected; ×4 scales the max to ~1)
- `fs_score = (obs_both − act1·act2) × 4`, expectation from each peak's overall
  actuation at that timepoint
- `odds_ratio = ((both+1)(neither+1)) / ((first+1)(second+1))` —
  Haldane–Anscombe; robust to marginal actuation shifts between timepoints, so
  pair strength stays comparable as overall accessibility changes
- `fisher_p` — two-sided Fisher exact test of the pair's 2×2 table
  (`[[neither, second], [first, both]]`) with a +1 pseudocount in every cell.
  The null is that the two peaks' open/closed states are independent across
  the co-spanning fibers, given each peak's marginal rate.

Distances are peak-midpoint distances, binned by `--bin-len` (100 bp; bins
below 2×bin_len merged into bin 1).

Outputs in `.../FIRE_MSP/codependency/`:

| file | contents |
|---|---|
| `peak_actuation_by_timepoint.tsv` | per peak: nCov / nMSP / actuation per timepoint, max_act, pass flag |
| `codep_<s>.csv.gz` | per pair: peak coords/actuations, 4 counts, score, fs_score, odds_ratio, dist, dist_bin, fisher_p |
| `bin_codep_<s>_avg.csv`, `_median.csv` | distance-binned score summaries |

### 3. `plot_codependency.R` (`module load R/4.4.1 && Rscript plot_codependency.R`)

Re-bins per-pair scores by `floor(log2(dist))` (<512 bp merged, ≥16 kb pooled)
and writes to `codependency/figures/` (each figure with its `.tbl.gz` table):

- `codep_distance_decay_global.pdf` — mean binned score, one line per timepoint
- `codep_fisher_sig_by_bin.csv` — per timepoint and distance bin: number and
  fraction of pairs with `fisher_p < 0.05`

## `archive/`

TSS-/early-response-gene-related code kept for later, currently unused:

- `make_gencode_v46_all_tss.sh` — GENCODE v46 transcript-TSS BED builder
  (20 bp windows; all TSSs + Ensembl_canonical subset)

cd /project/spott/cshan/fiber-seq/code/FIRE_MSP


# 1a. LPS_0 intersect first (builds the peak universe)
sbatch --array=1 intersect_msp_fire.sh

# 1b. once it completes, the remaining three timepoints
sbatch --array=2-4 intersect_msp_fire.sh
cd /project/spott/cshan/fiber-seq/code/FIRE_MSP

# 1a. LPS_0 intersect first (builds the peak universe)
sbatch --array=1 intersect_msp_fire.sh

# 1b. once it completes, the remaining three timepoints
sbatch --array=2-4 intersect_msp_fire.sh

# 2. codependency + Fisher tests, after all four intersects are done
sbatch run_fire_codependency.sh

# 3. figures, after the codependency job finishes
module load R/4.4.1
Rscript plot_codependency.R

