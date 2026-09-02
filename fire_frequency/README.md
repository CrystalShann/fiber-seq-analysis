# FIRE frequency across the LPS timecourse

Genome-wide, per-region fraction of reads carrying a FIRE element, per timepoint
(LPS_0, LPS_5, LPS_10, LPS_15). Groundwork for comparing FIRE frequency across
timepoints (differential accessibility); this step only computes and saves the
frequencies — no cross-timepoint testing yet.

## Definition

For each region of the shared peak universe and each timepoint:

- `n_reads` — distinct reads whose aligned span overlaps the region (>= 1 bp)
- `n_fire`  — those reads with >= 1 FIRE element (FDR <= 0.05) overlapping the region
- `freq`    — `n_fire / n_reads`, NA when `n_reads == 0`

Read-level, any-overlap: the same definitions as `03_coaccess_cres.py
--read-rule any` in `code/co-accessibility/`, applied to single regions instead of
cCRE pairs.

## Regions

Union of the four timepoints' FIRE peaks (`peak_start`/`peak_end`, merged):
`macrophage_project/co-accessibility/universe/fire_peaks_union.bed`
(143,070 intervals, chr1-22/X/Y, built by `01_make_cre_universe.sh`). A shared
universe means every region has a value at every timepoint, so frequencies are
directly comparable; per-timepoint peak sets differ (~115-120k each) and would not
line up.

## Inputs (per timepoint s)

- Denominator: `macrophage_project/co-accessibility/<s>/<s>.read_spans.bed.gz`
  (from `02_read_spans.sh`; aligned spans, `-F 0x900` filtered, col4 = read name)
- Numerator: `<FIRE>/<s>/additional-outputs-v0.1/fire-peaks/<s>-v0.1-fire-elements.bed.gz`
  with `FIRE = /project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE`
  (the FIRE pipeline's own per-read element calls, FDR <= 0.05, col4 = read name).
  Both derive from the same `-fire-v0.1-filtered.cram`, so read names match. The
  script counts "orphan" element hits on uncovered reads as a sanity check: the CRAM
  keeps supplementary alignments (`-F 260`) while the spans dropped them (`-F 0x900`),
  so a small orphan count (~3.5-4.3k per timepoint genome-wide, ~0.04% of element
  hits) is expected — those elements are excluded, keeping the denominator
  primary-alignment-only like `03_coaccess_cres.py`. A large count would flag a
  real mismatch.

## Outputs

`macrophage_project/fire_frequency/`

- `fire_frequency_long.tsv.gz` — one row per (region, timepoint):
  `chrom start end region_id timepoint n_reads n_fire freq` (572,280 rows).
  The shape for downstream stats: per-region binomial GLM / prop tests consume
  `(n_fire, n_reads, timepoint)` directly.
- `fire_frequency_wide.tsv.gz` — one row per region, columns
  `<s>_n_reads <s>_n_fire <s>_freq` per timepoint, genomic order.

## Run

```bash
sbatch run_fire_frequency.sh
```

## 02: pairwise differential FIRE frequency

`02_pairwise_differential.R` tests, for every union region and each of the 6
timepoint pairs (0v5, 0v10, 0v15, 5v10, 5v15, 10v15), whether FIRE frequency
differs: two-sided Fisher exact on the 2x2 of `(n_fire, n_reads - n_fire)` for
timepoint A vs B — the same test the FIRE pipeline's `hap-diffs.R` uses for
H1 vs H2, with its coverage window (per-timepoint median +/- 3*sqrt(median),
floor 10 reads; a region must be inside the window at both timepoints to pass).
No pseudocount.

```bash
module load R/4.4.2+gcc-13.2.0
Rscript 02_pairwise_differential.R      # a few minutes, foreground
```

Output: `macrophage_project/fire_frequency/pairwise_fisher.tsv.gz`, one row per
(region, comparison) = 858,420 rows:

```
chrom start end region_id comparison
n_reads_A n_fire_A freq_A n_reads_B n_fire_B freq_B
delta_freq log2_or p_value rank_p pass_coverage
```

