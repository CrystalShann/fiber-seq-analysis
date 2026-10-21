#!/bin/bash
#SBATCH --job-name=LCL_cre_footprints
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --qos=bigmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=400G
#SBATCH --time=24:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/LCL_tf_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/LCL_tf_%A_%a.err
set -euo pipefail
SCRIPT_DIR=${LCL_TF_CODE:-/project/spott/cshan/fiber-seq/code/co-accessibility/LCL/cre_summary}
ROOT=${LCL_COACCESS_ROOT:-/project/spott/cshan/fiber-seq/LCL_project/co-accessibility}
OUT=${LCL_TF_OUT:-$ROOT}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
stage=${1:-help}
shift || true
case "$stage" in
  footprints)
    # all FiberHMM footprints per CRE, one sample per task (LCL_coaccess_cre_summary.Rmd):
    #   sbatch --partition=caslake --qos=caslake --mem=48G --array=1-31%6 ... footprints
    module load R/4.4.1
    LCL_COACCESS_ROOT="$OUT" Rscript "$SCRIPT_DIR/count_cre_footprints.R" count "${SLURM_ARRAY_TASK_ID:-1}" ;;
  footprints_combine)
    module load R/4.4.1
    LCL_COACCESS_ROOT="$OUT" Rscript "$SCRIPT_DIR/count_cre_footprints.R" combine ;;
  *) echo 'Usage: bash LCL_tf_cooccupancy.sh {footprints|footprints_combine}' >&2; exit 2 ;;
esac
