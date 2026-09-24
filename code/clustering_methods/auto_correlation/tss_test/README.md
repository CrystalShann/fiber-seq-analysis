# protein-coding TSS Fiber-seq autocorrelation workflow

## Files and execution

- `run_tss.py`: canonical windows, R pipe, original ACF/clustering imports,
  annotation and expression orchestration; returns all results in memory.
- `tss_input.R`: sources the original expression script with local output
  interception; extracts spanning reads with the original shared R functions.
- `nrl_annotation.py`: separate 33-bp smoothing, annotation ACF, extrema and NRL.
- `expression_summary.py`: sample/TSS/gene summaries, gene ACFs and composition.
- `plots.py`: final PDF figures; smoke figures are rendered into `BytesIO` only.
- `run_tss.slurm.sh`: SLURM wrapper for the complete long-running Python workflow.


```bash
cd /project/spott/cshan/fiber-seq/code/clustering_methods/auto_correlation/tss_test
/project/spott/cshan/envs/Jupyter-notebook/bin/python -B run_tss.py --help
```

```bash
TSS_OUTPUT=/project/spott/cshan/fiber-seq/LCL_project/auto_correlation/tss_test
mkdir -p "$TSS_OUTPUT"
tss_job_id=$(sbatch --parsable run_tss.slurm.sh --all-tss)
echo "Submitted job: $tss_job_id"
```

The wrapper runs `run_tss.py`, which executes these steps sequentially in one
process with an R subprocess for input:

1. Define exact 2-kb canonical TSS windows, load the original expression bins,
   and retain only protein-coding, expression-matched TSSs within reference bounds.
2. Extract fully spanning m6A reads through `tss_input.R`.
3. Run the imported raw ACF and PCA/correlation-kNN/weighted Leiden functions.
4. Run `nrl_annotation.py` for the separate 33-bp NRL/regularity annotations.
5. Run `expression_summary.py` for the gene/TSS-level expression summaries.
6. Run `plots.py` to save only the final PDF figures.



## Input

Input: `/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed`.

Default Fiber-seq input is the macrophage
`macrophage_project/FiberHMM/extract/ft_result_dir`, matching the expression
script's LPS timepoints. All inputs must have existing tabix indexes


`code/topic_model/topic_modelling_functions.r`:

- `extract_ft_region_reads` (internally `read_tabix_region` and
  `convert_ft_bed12_to_bed6`): indexed retrieval, longest alignment per RID,
  removal of both fibertools sentinel blocks, one-based m6A positions.
- `extract_ft_read_info`: original read span and RID metadata.
- `get_sparse_met_mat`: consecutive 2,000-base binary 0/1 matrix, keeping all
  zero columns and requiring full coverage.



## Methods and parameters

1. Start from each retained **same 2,000-base binary signal**. Apply a uniform
   centered rolling mean of **33 bp**, with 16 neighbors on either side.
   Use complete windows only (`np.convolve(..., mode="valid")`): 1,968 values
   at original centers 16–1,983. The 16 undefined values at each edge are
   omitted, with no zero padding, reflection or partial windows.
2. Compute a separate linear, mean-centered ACF divided by the full centered
   sum of squares. This matches statsmodels `acf(adjusted=False, fft=False)`.
   Annotation lags are 0–1,967; statsmodels is not a runtime dependency.
3. Positive maxima: `scipy.signal.find_peaks(acf, prominence=0.02,
   distance=100, height=0.02)`. Negative minima: the identical call on `-acf`.
   Prominence uses the full signal (`wlen=None`); no width/plateau constraint.
   The height requirement enforces genuinely positive maxima and negative
   minima with absolute amplitude at least 0.02. Spacing is at least **100 bp
   between candidates of the same sign**, not between a maximum and minimum.
4. Keep extrema at lags **50–1,468 inclusive**, requiring at least **500 paired
   smoothed positions** at each lag. Peaks are detected before the lag filter,
   so boundary candidates use their actual neighboring values. This excludes
   the smoothing-related near-zero lobe and low-support tail. Keep **every**
   qualifying extremum in this range, with no cap of eight or any other count.
   No hand-picked cycle subset or significance claim is made.
5. Nucleosome-scale peak lag is the strongest accepted positive maximum in
   **140–250 bp** (inclusive), or missing if none exists. Count accepted
   positive peaks and all accepted extrema. With at least two positive peaks,
   fit `peak_lag = intercept + period * peak_number`, where peak numbers are
   consecutive `1..K`; report slope (NRL), free intercept and regression R².
   Do not force a zero intercept or silently impute missing cycles.
6. `extrema_abs_sum = sum(abs(acf[accepted_extrema]))`;
   `regularity = extrema_abs_sum / n_extrema`. Lag zero is never included.
   With no extrema, sum is 0 but regularity is missing; with constant smoothed
   signal, ACF and both amplitude metrics are missing and status is recorded.
7. A **reliable** NRL requires at least **3** positive peaks, an accepted
   140–250-bp peak, regression slope **140–250 bp**, **R² >= 0.90**, positive
   peak gap population CV **<= 0.25**, and maximum absolute regression residual
   **<= 0.20 × fitted period**. Missing cycles can yield an unreliable estimate;
   their original fit and all reasons are retained. A two-peak fit is always
   unreliable despite its necessarily perfect R².



Reference concepts (read-only, no downloaded code executed):
[statsmodels ACF example](https://github.com/YifanLab/ipdTrimming/blob/main/autocorr_withstatmodel.py),
[positive/negative find_peaks example](https://github.com/YifanLab/ipdTrimming/blob/main/m6A_axdistribution.py),
and [rolling-mean/TSS display examples](https://github.com/YifanLab/ipdTrimming/blob/main/enhancerPlot.R).
Those scripts do not define the exact 33-bp window and reliability rules here;
these are declared above rather than attributed to the reference scripts.

