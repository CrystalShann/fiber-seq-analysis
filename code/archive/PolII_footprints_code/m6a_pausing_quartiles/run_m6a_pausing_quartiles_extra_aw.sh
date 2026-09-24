#!/bin/bash

# Re-run run_m6a_pausing_quartiles.R for an additional TSS_ACCESS_WINDOW
# (250 / 500 / 1000 bp) WITHOUT recomputing the position-level outputs.
# Submit one job per window via submit_extra_aw_runs.sh.
#
# Required env: TSS_ACCESS_WINDOW (must be passed by the caller, e.g. via
# `sbatch --export=...` in submit_extra_aw_runs.sh). The accessibility
# table is written to:
#   ${OUTPUT_ROOT}/<run_dir>/tables/${SAMPLE}_tss_read_accessibility_aw${TSS_ACCESS_WINDOW}bp.tsv.gz
# leaving the existing 100 bp output (no aw suffix) untouched.

#SBATCH --job-name=m6a_quartiles_aw
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G
#SBATCH --time=24:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/m6a_quartiles_aw%x_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/m6a_quartiles_aw%x_%j.err

set -euo pipefail

if [[ -z "${TSS_ACCESS_WINDOW:-}" ]]; then
  echo "TSS_ACCESS_WINDOW must be set (e.g. 250, 500, 1000)" >&2
  exit 1
fi

module load R

export SAMPLE=AL10_bc2178_19130
export CHROMS=AUTO
export WINDOW_BP=1000
export BIN_SIZE=10
export MOD_PROB_THRESHOLD=0.9
export MODBAM_TSS_CHUNK_SIZE=200
export MAX_MODIFIED_HIST_BIN=10
export TSS_SOURCE_FILTER=ALL
export DO_POSITION_SCAN=0

export POLII_ROOT=/project/spott/cshan/fiber-seq/results/PolII
export PAUSING_PATH=${POLII_ROOT}/annotations/pausing_index_principal_with_CAGE_TSS_all_genes.tsv
export OUTPUT_ROOT=${POLII_ROOT}/m6a_pausing_quartiles
export BAM=/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples/${SAMPLE}.5mC.6mA.aligned.phased.bam

echo "Sample           : ${SAMPLE}"
echo "Window (bp)      : ${WINDOW_BP}"
echo "Bin (bp)         : ${BIN_SIZE}"
echo "AccessWin (bp)   : ${TSS_ACCESS_WINDOW}"
echo "DoPositionScan   : ${DO_POSITION_SCAN}"
echo "Output root      : ${OUTPUT_ROOT}"
echo "Started          : $(date)"

Rscript /project/spott/cshan/fiber-seq/code/PolII_footprints_code/m6a_pausing_quartiles/run_m6a_pausing_quartiles.R

echo "Finished         : $(date)"
